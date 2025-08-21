import pandas as pd
import glob
import os
import pyranges as pr
from tqdm import tqdm
import pysam
import numpy as np
import gc
from joblib import Parallel, delayed

# PARSING AND TRANSFORMATION
def parse_and_transform_junctions(file_path):
    """
    Parses a 12-column junction file, validates it, and transforms it into a
    standard BED-like format with corrected intron coordinates.
    """
    # Column names for RegTools junction files
    regtools_column_names = [
        'chrom', 'start_anchor', 'end_anchor', 'name', 'score', 'strand',
        'thick_start_orig', 'thick_end_orig', 'item_rgb_orig',
        'block_count_orig', 'block_sizes_orig', 'block_starts_orig'
    ]
    sample_id = os.path.basename(file_path).split('.')[0]

    if os.path.getsize(file_path) == 0:
        print(f"Warning: Skipping empty file {os.path.basename(file_path)}")
        return pd.DataFrame()

    # Read the file and validate column count
    try:
        df = pd.read_csv(file_path, sep='\t', header=None)
        if df.shape[1] != 12:
            print(f"Warning: Skipping {os.path.basename(file_path)} due to unexpected column count.")
            return pd.DataFrame()
        df.columns = regtools_column_names
    except Exception as e:
        print(f"Error reading or validating {os.path.basename(file_path)}: {e}")
        return pd.DataFrame()

    # Data cleaning and transformation
    df['sample_id_source'] = sample_id
    df.dropna(subset=['start_anchor', 'end_anchor', 'score', 'block_sizes_orig'], inplace=True)
    for col in ['start_anchor', 'end_anchor', 'score']:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    df.dropna(inplace=True)
    for col in ['start_anchor', 'end_anchor', 'score']:
        df[col] = df[col].astype(int)

    # Coord correction
    parsed_blocks = df['block_sizes_orig'].str.strip(',').str.split(',', expand=True)

    # Filter rows for valid blocks
    has_sufficient_blocks = parsed_blocks[1].notna()
    df = df[has_sufficient_blocks].copy()
    parsed_blocks = parsed_blocks[has_sufficient_blocks]
    
    #Recalculate junction coords
    overhang_left = pd.to_numeric(parsed_blocks[0])
    overhang_right = pd.to_numeric(parsed_blocks[1])

    junc_start = df['start_anchor'] + overhang_left
    junc_end = df['end_anchor'] - overhang_right

    # Create final transformed df
    transformed_df = pd.DataFrame({
        'chrom': df['chrom'],
        'chromStart': junc_start,
        'chromEnd': junc_end,
        'name': df['name'],
        'score': df['score'],
        'strand': df['strand'],
        'thickStart': junc_start,
        'thickEnd': junc_end,
        'itemRgb': df['item_rgb_orig'],
        'blockCount': 1,
        'blockSizes': (junc_end - junc_start).astype(str),
        'blockStarts': "0",
        'sample_id_source': df['sample_id_source']
    })

    # Filter out invalid junctions where start >= end
    transformed_df = transformed_df[transformed_df['chromStart'] < transformed_df['chromEnd']].copy()
    return transformed_df


# EXITRON FINDING 
def findExitrons(transformed_df, feature_pr, feature_name):
    """
    Finds both 'contained' and 'spanning' exitrons.
    """
    if transformed_df.empty:
        return pr.PyRanges()

    # Create initial pr object
    junction_df_for_pr = pd.DataFrame({
        'Chromosome': transformed_df['chrom'],
        'Start': transformed_df['chromStart'],
        'End': transformed_df['chromEnd'],
        'Strand': transformed_df['strand'],
        'title': (transformed_df['chrom'].astype(str) + ':' +
                  transformed_df['chromStart'].astype(str) + ':' +
                  transformed_df['chromEnd'].astype(str) + ':' +
                  transformed_df['strand'].astype(str)),
        'reads': transformed_df['score'],
        'sourceID': transformed_df['sample_id_source']
    })
    junction_pr = pr.PyRanges(junction_df_for_pr)

    # 1. FIND CONTAINED EXITRONS
    contained_junctions = junction_pr.overlap(feature_pr, how="containment", strandedness="same")
    contained_df = contained_junctions.df
    contained_df['exitron_type'] = 'contained'
    print(f"Found {len(contained_df)} 'contained' exitrons in {feature_name} regions.")

    # 2. FIND SPANNING EXITRONS
    spanning_df = pd.DataFrame()
    joined_pr = junction_pr.join(feature_pr)

    if not joined_pr.empty:
        joined_df = joined_pr.df
        joined_df['exon_number'] = pd.to_numeric(joined_df['exon_number'])

        # group by junction title and transcript_id
        grouped = joined_df.groupby(['title', 'transcript_id'])

        # define filter to find junctions spanning two adjacent exons
        def is_adjacent(g):

            # junction overlaps with 2 exons
            if len(g) != 2:
                return False
            
            # consecutive exons
            exon_numbers = g['exon_number'].unique()
            return len(exon_numbers) == 2 and abs(exon_numbers[0] - exon_numbers[1]) == 1

        adjacent_groups = grouped.filter(is_adjacent)

        if not adjacent_groups.empty:
            spanning_junctions = adjacent_groups['title'].unique()

            # filter original junctions to get spanning exitrons
            spanning_junctions_pr = junction_pr[junction_pr.title.isin(spanning_junctions)]
            spanning_df = spanning_junctions_pr.df
            spanning_df['exitron_type'] = 'spanning'

    print(f"Found {len(spanning_df)} 'spanning' exitrons in {feature_name} regions.")

    # Combine results
    combined_df = pd.concat([contained_df, spanning_df], ignore_index=True)
    if combined_df.empty:
        return pr.PyRanges()
    return pr.PyRanges(combined_df)


# PARALLEL PROCESSING AND NORMALIZATION
def process_file_and_save(file_path, features, output_dir):
    """
    Processes one file and saves its results to disk.
    """
    try:
        sample_id = os.path.basename(file_path).split('.')[0]
        output_filepath = os.path.join(output_dir, f"{sample_id}_exitrons.parquet")

        if os.path.exists(output_filepath):
            print(f"Skipping {sample_id}, output already exists.")
            return

        # 1. Parse and transform data
        transformed_df = parse_and_transform_junctions(file_path)
        if transformed_df.empty:
            print(f"No valid junctions found in {sample_id} after parsing.")
            return

        # 2. Find exitrons for each feature type
        all_exitron_info = []
        for feature_name, feature_pr in features.items():
            exitron_pr = findExitrons(transformed_df, feature_pr, feature_name)
            if not exitron_pr.empty:
                exitron_df = exitron_pr.df
                exitron_df['FeatureType'] = feature_name
                all_exitron_info.append(exitron_df)

        # 3. Concatenate and save results for this single file
        if all_exitron_info:
            final_df = pd.concat(all_exitron_info, ignore_index=True)
            final_df.to_parquet(output_filepath, index=False)
            print(f"-> Successfully saved results for {sample_id} to {output_filepath}")
        
        del transformed_df, all_exitron_info
        gc.collect()

    except Exception as e:
        print(f"!! Unhandled error processing file {os.path.basename(file_path)}: {e}")

def normalize_one_sample(sample_id, read_counts, annotations, bam_dir, extend):
    """Helper function to normalize one sample (a column in the df). Can be run in parallel."""
    bam_filepath = os.path.join(bam_dir, f'{sample_id}.bam')
    normalized_counts = np.zeros(len(read_counts))
    
    try:
        with pysam.AlignmentFile(bam_filepath, 'rb') as bam:
            for i, (junc_title, numerator) in enumerate(read_counts.items()):
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
        print(f"Warning: BAM file not found for {sample_id}, skipping normalization.")
    except Exception as e:
        print(f"Error processing BAM for {sample_id}: {e}")
        
    return normalized_counts

def normalizeExitronData_parallel(filtered_exitron_data, output_filepath, bam_dir, n_jobs=-1):
    """Parallelized version of the normalization function."""
    # Prepare annotations for all junctions
    annots = pd.DataFrame({"exitron": filtered_exitron_data.index})
    annots[['chr', 'start', 'end', 'strand']] = annots['exitron'].str.split(':', expand=True)
    annots['start'] = pd.to_numeric(annots['start'])
    annots['end'] = pd.to_numeric(annots['end'])
    annots.set_index('exitron', inplace=True)
    
    EXTEND = 6

    # Use joblib to run normalization for each sample in parallel
    results = Parallel(n_jobs=n_jobs)(
        delayed(normalize_one_sample)(
            col, filtered_exitron_data[col], annots, bam_dir, EXTEND
        )
        for col in tqdm(filtered_exitron_data.columns, desc="Normalizing Samples")
    )
    
    # Combine results and save
    normalized_matrix = np.array(results).T
    np.save(output_filepath, normalized_matrix)
    print(f"Successfully saved normalized data to {output_filepath}")
    return normalized_matrix


#  MAIN EXECUTION
def main():
    # 0. Setup Paths and Parameters
    junctions_dir = "/gpfs/commons/groups/knowles_lab/atokolyi/als/juncs_min10bp/"
    temp_output_dir = "/gpfs/commons/home/ncui/project/temp_exitron_results" # For intermediate files
    final_unfiltered_path = "/gpfs/commons/home/ncui/project/unfil_exitrons.parquet"
    filtered_path = "/gpfs/commons/home/ncui/project/fil_exitrons.parquet"
    normalized_path = 'norm_exitrons.npy'
    bam_dir = "/gpfs/commons/projects/ALS_Consortium_resource/bam_files/"
    N_JOBS = 8 # Adjust based on your machine's core count

    os.makedirs(temp_output_dir, exist_ok=True)

    # 1. Load Genomic Annotations
    print("Reading GFF annotations...")
    gff = pr.read_gff3("/gpfs/commons/home/ncui/project/reference_data/gencode.v48.annotation.gff3.gz")
    features = {
        "CDS": gff[gff.Feature == "CDS"],
        "exon": gff[gff.Feature == "exon"],
        "five_prime_UTR": gff[gff.Feature == "five_prime_UTR"],
        "three_prime_UTR": gff[gff.Feature == "three_prime_UTR"]
    }
    print("Annotations loaded.")

    # 2. Process Junction Files in Parallel
    file_paths = glob.glob(os.path.join(junctions_dir, "*.junc")) # [:10] try first 10 
    print(f"\nFound {len(file_paths)} junction files to process.")
    Parallel(n_jobs=N_JOBS)(
        delayed(process_file_and_save)(path, features, temp_output_dir)
        for path in tqdm(file_paths, desc="Processing Junction Files")
    )
    
    # 3. Combine Temporary Files
    print("\nCombining all temporary exitron files...")
    temp_files = glob.glob(os.path.join(temp_output_dir, "*.parquet"))
    if not temp_files:
        print("No temporary files were created. Exiting.")
        return
        
    combined_df = pd.concat(
        [pd.read_parquet(f) for f in tqdm(temp_files, desc="Combining Files")],
        ignore_index=True
    )
    combined_df.to_parquet(final_unfiltered_path, index=False)
    print(f"Combined data has {len(combined_df)} rows. Saved to {final_unfiltered_path}")

    # 4. Filter Combined Data
    print("\nFiltering combined exitron data...")
    long_df = pd.read_parquet(final_unfiltered_path)
    long_df.drop_duplicates(subset=['title', 'sourceID'], inplace=True)
    try:
        wide_df = long_df.pivot(index='title', columns='sourceID', values='reads').fillna(0)
    except MemoryError:
        print("Error: Pivoting caused a MemoryError")
        return
    
    min_count = 30
    min_peeps = 10
    keep_juncs = wide_df[wide_df > min_count].sum(axis=1) >= min_peeps
    filtered_data = wide_df[keep_juncs].copy()
    filtered_data.to_parquet(filtered_path)
    print(f"Found {len(filtered_data)} junctions passing filter. Saved to {filtered_path}")

    # 5. Normalize Filtered Data in Parallel
    if not filtered_data.empty:
        print("\nNormalizing filtered exitron data...")
        normalizeExitronData_parallel(filtered_data, normalized_path, bam_dir, n_jobs=N_JOBS)

if __name__ == "__main__":
    main()
