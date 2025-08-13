import pandas as pd
import glob
import os
import pyranges as pr
from tqdm import tqdm
import pysam
import numpy as np
import gc 
from joblib import Parallel, delayed 

## Data type and column specifications for efficiency
DTYPE_SPEC = {
    'chrom': 'category', 'start_anchor': 'int32', 'end_anchor': 'int32',
    'name': 'object', 'score': 'int32', 'strand': 'category',
    'item_rgb_orig': 'object', 'block_count_orig': 'int16',
    'block_sizes_orig': 'object', 'block_starts_orig': 'object'
}
COLS_TO_USE = ['chrom', 'start_anchor', 'end_anchor', 'name', 'score', 'strand', 'item_rgb_orig', 'block_sizes_orig']

def parseJunctionFile(file_path):
    # This function is unchanged from the previous optimization
    sample_id = os.path.basename(file_path).split('.')[0]
    df = pd.read_csv(
        file_path, sep='\t', header=None, names=DTYPE_SPEC.keys(),
        usecols=COLS_TO_USE, dtype={'chrom': str}
    )
    df['sample_id_source'] = sample_id
    std_chroms = {f'chr{i}' for i in range(1, 23)} | {'chrX', 'chrY', 'chrM'}
    parsed_blocks = df['block_sizes_orig'].str.strip(',').str.split(',', expand=True)
    overhang_left = pd.to_numeric(parsed_blocks[0], errors='coerce')
    overhang_right = pd.to_numeric(parsed_blocks.iloc[:, -1], errors='coerce')
    transformed_df = df.assign(
        chrom=df['chrom'].astype('category'),
        chromStart=(df['start_anchor'] + overhang_left).astype('int32'),
        chromEnd=(df['end_anchor'] - overhang_right).astype('int32'),
    ).dropna(subset=['chromStart', 'chromEnd'])
    transformed_df = transformed_df[
        transformed_df['chrom'].isin(std_chroms) &
        (transformed_df['chromStart'] < transformed_df['chromEnd'])
    ].copy()
    final_df = pd.DataFrame({
        'chrom': transformed_df['chrom'], 'chromStart': transformed_df['chromStart'],
        'chromEnd': transformed_df['chromEnd'], 'name': transformed_df['name'],
        'score': transformed_df['score'].astype('int32'), 'strand': transformed_df['strand'].astype('category'),
        'thickStart': transformed_df['chromStart'], 'thickEnd': transformed_df['chromEnd'],
        'itemRgb': transformed_df['item_rgb_orig'], 'blockCount': 1,
        'blockSizes': (transformed_df['chromEnd'] - transformed_df['chromStart']).astype('int32'),
        'blockStarts': 0, 'sample_id_source': transformed_df['sample_id_source']
    })
    return final_df

def findExitrons(transformed_df, feature_pr, feature_name):
    """
    Finds both 'contained' and 'spanning' exitrons using a memory-efficient approach.
    """
    if transformed_df.empty:
        return pr.PyRanges()

    # Create the initial PyRanges object for junctions
    junction_df_for_pr = pd.DataFrame({
        'Chromosome': transformed_df['chrom'], 'Start': transformed_df['chromStart'],
        'End': transformed_df['chromEnd'], 'Strand': transformed_df['strand'],
        'title': (transformed_df['chrom'].astype(str) + ':' +
                  transformed_df['chromStart'].astype(str) + ':' +
                  transformed_df['chromEnd'].astype(str) + ':' +
                  transformed_df['strand'].astype(str)),
        'reads': transformed_df['score'], 'sourceID': transformed_df['sample_id_source']
    })
    junction_pr = pr.PyRanges(junction_df_for_pr)

    # 1. FIND CONTAINED EXITRONS
    contained_junctions = junction_pr.overlap(feature_pr, how="containment", strandedness="same")
    if not contained_junctions.empty:
        contained_df = contained_junctions.df
        contained_df['exitron_type'] = 'contained'
    else:
        contained_df = pd.DataFrame()
    print(f"Found {len(contained_df)} 'contained' exitrons in {feature_name} regions.")

    # 2. FIND SPANNING EXITRONS (New Memory-Efficient Logic)
    spanning_df = pd.DataFrame()
    joined_pr = junction_pr.join(feature_pr, strandedness="same")

    if not joined_pr.empty:
        joined_df = joined_pr.df
        joined_df['exon_number'] = pd.to_numeric(joined_df['exon_number'], errors='coerce')
        joined_df.dropna(subset=['exon_number'], inplace=True)
        joined_df['exon_number'] = joined_df['exon_number'].astype(int)
        
        group_keys = ['title', 'transcript_id']
        group_sizes = joined_df.groupby(group_keys)['exon_number'].transform('size')
        two_exon_overlaps_df = joined_df[group_sizes == 2].copy()

        if not two_exon_overlaps_df.empty:
            grouped = two_exon_overlaps_df.groupby(group_keys)
            exon_diff = grouped['exon_number'].max() - grouped['exon_number'].min()
            adjacent_junction_groups = exon_diff[exon_diff == 1]
            if not adjacent_junction_groups.empty:
                spanning_junction_titles = adjacent_junction_groups.index.get_level_values('title').unique()
                spanning_pr = junction_pr[junction_pr.title.isin(spanning_junction_titles)]
                if not spanning_pr.empty:
                    spanning_df = spanning_pr.df
                    spanning_df['exitron_type'] = 'spanning'
    print(f"Found {len(spanning_df)} 'spanning' exitrons in {feature_name} regions.")

    # 3. COMBINE RESULTS
    combined_df = pd.concat([contained_df, spanning_df], ignore_index=True)

    if combined_df.empty:
        return pr.PyRanges()
    return pr.PyRanges(combined_df)

def process_file_and_save(file_path, features, output_dir):
    """
    ## OPTIMIZATION: This is the 'map' step.
    Processes a single file and saves its results to a temporary parquet file.
    This function is designed to be run in a loop or in parallel.
    """
    try:
        sample_id = os.path.basename(file_path).split('.')[0]
        print(f"\nProcessing file: {sample_id}")
        
        output_filepath = os.path.join(output_dir, f"{sample_id}_exitrons.parquet")

        # Skip if already processed
        if os.path.exists(output_filepath):
            print(f"Skipping {sample_id}, output already exists.")
            return

        # 1. Parse and transform data
        transformed_df = parseJunctionFile(file_path)

        # 2. Find exitrons for each feature type
        all_exitron_info = []
        for feature_name, feature_pr in features.items():
            exitron_pr = findExitrons(transformed_df, feature_pr, feature_name)
            if not exitron_pr.empty:
                exitron_df = exitron_pr.df
                exitron_df['FeatureType'] = feature_name
                all_exitron_info.append(exitron_df)
        
        # 3. Concatenate and save results for THIS FILE
        if all_exitron_info:
            final_df = pd.concat(all_exitron_info, ignore_index=True)
            final_df.to_parquet(output_filepath, index=False)
            print(f"-> Successfully saved results for {sample_id} to {output_filepath}")

        # 4. Free up memory
        del transformed_df, all_exitron_info
        gc.collect()

    except Exception as e:
        print(f"!! Error processing file {os.path.basename(file_path)}: {e}")

def combine_exitron_data(temp_dir, final_output_path):
    """
    ## OPTIMIZATION: This is the 'reduce' step.
    Combines all the temporary parquet files into a single final file.
    """
    print("\nCombining all temporary exitron files...")
    temp_files = glob.glob(os.path.join(temp_dir, "*.parquet"))
    if not temp_files:
        print("No temporary files found to combine.")
        return
    
    df_list = [pd.read_parquet(f) for f in tqdm(temp_files)]
    combined_df = pd.concat(df_list, ignore_index=True)
    
    print(f"Saving combined data with {len(combined_df)} rows to {final_output_path}")
    combined_df.to_parquet(final_output_path, index=False)
    return combined_df

def filterExitronData(exitron_data_filepath, output_filepath, min_count=30, min_peeps=10):
    """
    ## OPTIMIZATION: Pivot can still use a lot of memory. This logic is kept,
    ## but if it fails, the data may be too large for a wide format on this machine.
    """
    long_df = pd.read_parquet(exitron_data_filepath)
    long_df.drop_duplicates(subset=['title', 'sourceID'], inplace=True)

    try:
        wide_df = long_df.pivot(index='title', columns='sourceID', values='reads').fillna(0)
    except MemoryError:
        print("Error: Pivoting the data frame caused a MemoryError. The resulting matrix is too large.")
        print("Consider increasing memory allocation or using sparse matrices/Dask for this step.")
        return None

    keep_juncs = wide_df[wide_df > min_count].sum(axis=1) >= min_peeps
    
    filtered_exitron_data = wide_df[keep_juncs].copy()
    filtered_exitron_data.to_parquet(output_filepath)
    print(f"Found {len(filtered_exitron_data)} junctions passing filter. Saved to {output_filepath}")
    return filtered_exitron_data

def normalize_one_sample(sample_id, filtered_exitrons_for_sample, annotations, bam_dir, extend):
    """
    ## OPTIMIZATION: Helper function to normalize one sample (column).
    This can be run in parallel.
    """
    bam_filepath = os.path.join(bam_dir, f'{sample_id}.bam')
    normalized_counts = np.zeros(len(filtered_exitrons_for_sample))

    try:
        with pysam.AlignmentFile(bam_filepath, 'rb') as bam:
            for i, (junc_title, numerator) in enumerate(filtered_exitrons_for_sample.items()):
                if numerator > 0:
                    ann = annotations.loc[junc_title]
                    depth = bam.count_coverage(
                        ann['chr'], ann['start'] - extend, ann['end'] + extend, quality_threshold=0
                    )
                    totals = np.sum(depth, axis=0)
                    denom = np.median(np.concatenate([totals[:extend], totals[-extend:]]))
                    
                    if denom > 0:
                        normalized_counts[i] = numerator / denom
    except FileNotFoundError:
        print(f"Warning: BAM file not found for {sample_id}. Skipping.")
    except Exception as e:
        print(f"Error processing {sample_id}: {e}")
        
    return normalized_counts

def normalizeExitronData(filtered_exitron_data, output_filepath, bam_dir, n_jobs=-1):
    """
    ## OPTIMIZATION: Parallelized normalization using joblib.
    This dramatically speeds up the I/O-bound normalization step.
    """
    annots = pd.DataFrame({"exitron": filtered_exitron_data.index})
    annots[['chr', 'start', 'end', 'strand']] = annots['exitron'].str.split(':', expand=True)
    annots['start'] = pd.to_numeric(annots['start'])
    annots['end'] = pd.to_numeric(annots['end'])
    annots.set_index('exitron', inplace=True)
    
    EXTEND = 6

    # Use joblib to parallelize the loop over samples
    results = Parallel(n_jobs=n_jobs)(
        delayed(normalize_one_sample)(
            col, filtered_exitron_data[col], annots, bam_dir, EXTEND
        )
        for col in tqdm(filtered_exitron_data.columns, desc="Normalizing Samples")
    )
    
    # Transpose results to match original shape and save
    normalized_df = np.array(results).T
    np.save(output_filepath, normalized_df)
    print(f"Successfully saved normalized data to {output_filepath}")
    return normalized_df


def main():
    # --- Step 1: Process individual junction files ---
    # Setup directories
    junctions_dir = "/gpfs/commons/groups/knowles_lab/atokolyi/als/juncs_min10bp/"
    temp_output_dir = "/gpfs/commons/home/ncui/project/temp_exitron_results"
    final_unfiltered_path = "/gpfs/commons/home/ncui/project/unfil_exitrons.parquet"
    os.makedirs(temp_output_dir, exist_ok=True)

    # Load genomic annotations once
    print("Reading GFF annotations...")
    gff = pr.read_gff3("/gpfs/commons/home/ncui/project/reference_data/gencode.v48.annotation.gff3.gz")
    features = {
        "CDS": gff[gff.Feature == "CDS"],
        "exon": gff[gff.Feature == "exon"]
        # Add other features if needed, e.g., "five_prime_UTR", "three_prime_UTR"
    }
    
    # Process each file
    file_paths = glob.glob(os.path.join(junctions_dir, "*.junc"))
    for file_path in file_paths:
        process_file_and_save(file_path, features, temp_output_dir)

    # --- Step 2: Combine processed files ---
    combine_exitron_data(temp_output_dir, final_unfiltered_path)
    
    # --- Step 3: Filter the combined data ---
    filtered_path = "/gpfs/commons/home/ncui/project/fil_exitrons.parquet"
    filtered_data = filterExitronData(final_unfiltered_path, filtered_path)

    # --- Step 4: Normalize the filtered data ---
    if filtered_data is not None:
        bam_dir = "/gpfs/commons/projects/ALS_Consortium_resource/bam_files/"
        normalized_path = 'norm_exitrons.npy'
        normalizeExitronData(filtered_data, normalized_path, bam_dir, n_jobs=8) # Adjust n_jobs as needed

if __name__ == "__main__":
    main()
