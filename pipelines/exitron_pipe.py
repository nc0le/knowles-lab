import pandas as pd
import glob
import os
from tqdm import tqdm # progress tracker
import pyranges as pr # parsing gff
import numpy as np
import pysam

def parseJunctionFile(file_path):
    # column names for RegTools junction files
    regtools_column_names = [
        'chrom', 'start_anchor', 'end_anchor', 'name', 'score', 'strand',
        'thick_start_orig', 'thick_end_orig', 'item_rgb_orig',
        'block_count_orig', 'block_sizes_orig', 'block_starts_orig'
    ]
    
    # extract sample ID from the filename
    sample_id = os.path.basename(file_path).split('.')[0]
    
    # read the file into a pandas DataFrame
    df = pd.read_csv(
        file_path, sep='\t', header=None, names=regtools_column_names,
        dtype={'chrom': str, 'block_sizes_orig': str, 'block_starts_orig': str}
    )
        
    df['sample_id_source'] = sample_id

    # convert relevant columns to numeric types, coercing errors
    for col in ['start_anchor', 'end_anchor', 'score']:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
    
    # drop rows if info is missing
    df.dropna(subset=['start_anchor', 'end_anchor', 'score', 'block_sizes_orig'], inplace=True)
    
    # ensure int types
    for col in ['start_anchor', 'end_anchor', 'score']:
        df[col] = df[col].astype(int)

    return df

def transformJunctionData(raw_df):
    
    # CHROMOSOME FILTERING
    original_row_count = len(raw_df)
    
    # allowed chromosomes
    allowed_chrom_numbers = [str(i) for i in range(1, 23)]
    allowed_sex_chroms_upper = ['X', 'Y'] 
    allowed_chromosomes = set()
    for num_chrom in allowed_chrom_numbers:
        allowed_chromosomes.add(num_chrom)
        allowed_chromosomes.add(f"chr{num_chrom}")
    for sex_chrom in allowed_sex_chroms_upper:
        allowed_chromosomes.add(sex_chrom)
        allowed_chromosomes.add(sex_chrom.lower())
        allowed_chromosomes.add(f"chr{sex_chrom}")
        allowed_chromosomes.add(f"chr{sex_chrom.lower()}")
    
    raw_df_filtered = raw_df[raw_df['chrom'].isin(allowed_chromosomes)].copy()
    filtered_row_count = len(raw_df_filtered)
    print(f"Removed {original_row_count - filtered_row_count} rows with non-standard chromosomes.")


    # JUNCTION COORD CORRECTION
    # filter rows for valid blocks
    parsed_blocks_list = raw_df_filtered['block_sizes_orig'].str.strip(',').str.split(',')
    has_sufficient_blocks = parsed_blocks_list.str.len() >= 2
    raw_df_filtered = raw_df_filtered[has_sufficient_blocks].copy()
    parsed_blocks_list = parsed_blocks_list[has_sufficient_blocks]
    
    # recalculating junction coordinates
    raw_df_filtered.loc[:, 'overhang_left'] = parsed_blocks_list.str[0].astype(int)
    raw_df_filtered.loc[:, 'overhang_right'] = parsed_blocks_list.str[1].astype(int)

    junc_start = raw_df_filtered['start_anchor'] + raw_df_filtered['overhang_left']
    junc_end = raw_df_filtered['end_anchor'] - raw_df_filtered['overhang_right']

    # filter out invalid junctions
    valid_junction = junc_start < junc_end
    raw_df_filtered = raw_df_filtered[valid_junction].copy()
    junc_start = junc_start[valid_junction]
    junc_end = junc_end[valid_junction]


    junc_length = junc_end - junc_start

    # create df
    transformed_df = pd.DataFrame()
    transformed_df['chrom'] = raw_df_filtered['chrom']
    transformed_df['chromStart'] = junc_start
    transformed_df['chromEnd'] = junc_end
    transformed_df['name'] = raw_df_filtered['name']
    transformed_df['score'] = raw_df_filtered['score']
    transformed_df['strand'] = raw_df_filtered['strand']
    transformed_df['thickStart'] = junc_start
    transformed_df['thickEnd'] = junc_end
    transformed_df['itemRgb'] = raw_df_filtered['item_rgb_orig']
    transformed_df['blockCount'] = 1
    transformed_df['blockSizes'] = junc_length.astype(str)
    transformed_df['blockStarts'] = "0"
    transformed_df['sample_id_source'] = raw_df_filtered['sample_id_source']

    print(f"Transformed {len(transformed_df)} junction records.")
    
    return transformed_df

gff = pr.read_gff3("gencode.v48.annotation.gff3.gz")
cds = gff[gff.Feature == "CDS"]
print(f"Found {len(cds)} CDS regions.")

def findExitrons(transformed_df):
    transformed_df = transformed_df[transformed_df['strand'].isin(['+', '-'])]

    # generate a unique ID for each junction (chrom:start:end:strand
    unique_id = transformed_df['chrom'].astype(str) + ':' + \
                transformed_df['chromStart'].astype(str) + ':' + \
                transformed_df['chromEnd'].astype(str) + ':' + \
                transformed_df['strand'].astype(str)

    # convert junction data to PyRanges object
    junction_pr = pr.PyRanges({'Chromosome': transformed_df['chrom'],
                    'Start': transformed_df['chromStart'],
                    'End': transformed_df['chromEnd'],
                    'Strand': transformed_df['strand'],
                    'title': unique_id,
                    'reads': transformed_df['score'],
                    'sourceID': transformed_df['sample_id_source']}) 

    # find overlapping junctions
    contained_junctions = junction_pr.overlap(cds, contained_intervals_only=True, strand_behavior='same')
    print(f"Found {len(contained_junctions)} junctions contained within CDS regions.")
            
    return contained_junctions

def compileExitronData(directory_path, output_filepath, file_pattern="*.bam.junc"):

    all_exitron_info = []
    file_paths = glob.glob(os.path.join(directory_path, file_pattern))
    print(f"Found {len(file_paths)} files to process.")

    # testing first 5 out of 100
    '''
    files_to_process = file_paths[:5]
    print(f"Processing the first {len(files_to_process)} files.")
    '''

    for file_path in tqdm(file_paths):
        print("Parsing new file...")
        file_name_only = os.path.basename(file_path)
        try:
            # 1.
            parsed_data = parseJunctionFile(file_path)
            # 2.
            transformed_df = transformJunctionData(parsed_data)
            # 3.
            gr_file = findExitrons(transformed_df)
            #4.
            all_exitron_info.append(gr_file)

        # skip to the next file if an error occurs
        except Exception as e:
            print(f"An error occurred while processing file {file_name_only}: {e}")
            import traceback
            traceback.print_exc()
            continue 

    # concatenate all individual data into matrix 
    final_gr = pr.concat(all_exitron_info)
    print(f"\nSuccessfully compiled exitron data from {len(all_exitron_info)} files.")
    final_gr.to_parquet(output_filepath, index=False)
    print(f"Successfully saved data to {output_filepath}")
    return final_gr

def filterExitronData(exitron_data_filepath, min_count=30, min_peeps=10):
    long_df= pd.read_parquet(exitron_data_filepath)
    # drop duplicates created by pyranges .overlap method
    long_df.drop_duplicates(subset=['title', 'sourceID'], inplace=True)

    # create filtered wide df
    wide_df = long_df.pivot(index='title', columns='sourceID', values='reads').fillna(0)

    wide_is_gt_count = wide_df > min_count
    wide_is_gt_count_sums = wide_is_gt_count.sum(axis=1)
    keep_juncs = wide_is_gt_count_sums[wide_is_gt_count_sums >= min_peeps].index
    print(f"Found {len(keep_juncs)} junctions that passed the filter.")

    # construct final matrix
    filtered_exitron_data = wide_df.loc[keep_juncs].copy()
    filtered_exitron_data = filtered_exitron_data.where(filtered_exitron_data >= min_count, 0)

    return filtered_exitron_data

def normalizeExitronData(filtered_exitron_data, output_filepath):
    annots = pd.DataFrame({"exitron": list(filtered_exitron_data.index)})
    annots = pd.concat([annots, annots['exitron'].str.split(':', expand=True)], axis=1)
    annots.columns = ["exitron", "chr", "start", "end", "strand"]
    annots['start'] = annots['start'].astype(int)
    annots['end'] = annots['end'].astype(int)

    EXTEND = 6
    normalized_df = np.zeros(filtered_exitron_data.shape)
    ##tqdm(range(len(filtered_exitron_data.columns))):
    for i in range(len(filtered_exitron_data.columns)):
        print("Doing",i)
        source_id = filtered_exitron_data.columns[i]
        bam_filepath = f'/gpfs/commons/projects/ALS_Consortium_resource/bam_files/{source_id}.bam'

        with pysam.AlignmentFile(bam_filepath, 'rb') as bam:
            for j in range(len(filtered_exitron_data.index)):
                numerator = filtered_exitron_data.iat[j, i]

                if numerator > 0:
                    junc_chr = annots['chr'][j]
                    junc_start = annots['start'][j]
                    junc_end = annots['end'][j]

                    # calculate denominator
                    depth = bam.count_coverage(junc_chr, junc_start - EXTEND, junc_end + EXTEND, quality_threshold=0)
                    totals = [sum(values) for values in zip(*depth)]
                    
                    denom = np.median(totals[:EXTEND] + totals[-EXTEND:])
                    if denom > 0:
                        normalized_df[j, i] = numerator / denom
    np.save(output_filepath, normalized_df)
    return normalized_df

# change first argument to directory containing junction data 
# change second argument to parquet file for exitron data output 
# change third argument to correct junction file identifier
compileExitronData("/gpfs/commons/groups/knowles_lab/atokolyi/als/juncs_min6bp/", "/gpfs/commons/home/ncui/project/exitron_data.parquet", file_pattern="*.bam.junc")

# change argument to name of exitron data parquet file 
filtered_exitron_data = filterExitronData('final_exitron_data.parquet')

# change argument to name of normalized exitron expression data parquet file 
normalizeExitronData(filtered_exitron_data, 'normalized_data.npy')
