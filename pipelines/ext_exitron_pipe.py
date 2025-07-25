import pandas as pd
import glob
import os
import pyranges as pr # parsing gff

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

print("reading gff")
# load gff3 and extract features
gff = pr.read_gff3("/gpfs/commons/home/ncui/project/reference_data/gencode.v48.annotation.gff3.gz")

features = {
    "CDS": gff[gff.Feature == "CDS"],
    "exon": gff[gff.Feature == "exon"],
    "five_prime_UTR": gff[gff.Feature == "five_prime_UTR"],
    "three_prime_UTR": gff[gff.Feature == "three_prime_UTR"]
}

for name, data in features.items():
    print(f"Found {len(data)} {name} regions.")

def findExitrons(transformed_df, feature, feature_name):
    transformed_df = transformed_df[transformed_df['strand'].isin(['+', '-'])]

    # create unique junction ID
    '''
    junction_id = transformed_df['chrom'].astype(str) + ':' + \
                  transformed_df['chromStart'].astype(str) + ':' + \
                  transformed_df['chromEnd'].astype(str) + ':' + \
                  transformed_df['strand'].astype(str)
    '''

    junction_df_for_pr = pd.DataFrame({
        'Chromosome': transformed_df['chrom'],
        'Start': transformed_df['chromStart'],
        'End': transformed_df['chromEnd'],
        'Strand': transformed_df['strand'],
        'title': transformed_df['chrom'].astype(str) + ':' + \
                 transformed_df['chromStart'].astype(str) + ':' + \
                 transformed_df['chromEnd'].astype(str) + ':' + \
                 transformed_df['strand'].astype(str),
        'reads': transformed_df['score'],
        'sourceID': transformed_df['sample_id_source']
    })

    # Convert the prepared DataFrame to a PyRanges object
    junction_pr = pr.PyRanges(junction_df_for_pr)

    # FIND CONTAINED EXITRONS
    contained_junctions = junction_pr.overlap(feature, how="containment", strandedness="same")
    contained_df = contained_junctions.df
    if not contained_df.empty:
        contained_df['exitron_type'] = 'contained'
    print(f"Found {len(contained_df)} 'contained' exitrons in {feature_name} regions.")

    # FIND SPANNING EXITRONS
    bridging_df = pd.DataFrame()  # Initialize empty dataframe

    required_cols = ['transcript_id', 'exon_number']
    if not all(col in feature.columns for col in required_cols):
        print(f"Warning: Required columns {required_cols} not in {feature_name} features. Skipping bridging exitron search.")
    else:
        # Join junctions with features to find all overlaps
        joined_pr = junction_pr.join(feature)

        if not joined_pr.empty:
            joined_df = joined_pr.df
            joined_df['exon_number'] = pd.to_numeric(joined_df['exon_number'])

            # Group by junction ('title') and transcript_id
            grouped = joined_df.groupby(['title', 'transcript_id'])

            # Define a filter to find junctions spanning two adjacent exons
            def is_adjacent_bridge(g):
                if len(g) != 2:
                    return False
                exon_numbers = g['exon_number'].unique()
                return len(exon_numbers) == 2 and abs(exon_numbers[0] - exon_numbers[1]) == 1

            adjacent_groups = grouped.filter(is_adjacent_bridge)

            if not adjacent_groups.empty:
                # Get the unique titles of the bridging junctions
                bridging_junction_titles = adjacent_groups['title'].unique()

                # Filter original junctions to get the bridging exitrons
                bridging_junctions_pr = junction_pr[junction_pr.title.isin(bridging_junction_titles)]
                bridging_df = bridging_junctions_pr.df
                bridging_df['exitron_type'] = 'bridging'

    print(f"Found {len(bridging_df)} 'bridging' exitrons in {feature_name} regions.")

    # combine and return results
    combined_df = pd.concat([contained_df, bridging_df], ignore_index=True)

    if combined_df.empty:
        return pr.PyRanges()  # Return an empty PyRanges object if no exitrons were found

    return pr.PyRanges(combined_df)


def compileExitronData(directory_path, output_filepath, file_pattern):
    all_exitron_info = []
    file_paths = glob.glob(os.path.join(directory_path, file_pattern))
    print(f"Found {len(file_paths)} files to process.")
    print("starting parse...")

    for file_path in file_paths:
        print(f"\n--- Processing file: {os.path.basename(file_path)} ---")
        try:
            # 1. Parse and transform data
            parsed_data = parseJunctionFile(file_path)
            transformed_df = transformJunctionData(parsed_data)

            # 2. Find exitrons for each feature type
            print("Finding exitrons...")
            for feature_name, feature_pr in features.items():
                print(f"-> Analyzing '{feature_name}' regions...")
                # Pass the feature_name to the function
                exitron_pr = findExitrons(transformed_df, feature_pr, feature_name)

                # Only append if exitrons were found
                if not exitron_pr.empty:
                    exitron_pr.FeatureType = feature_name
                    all_exitron_info.append(exitron_pr)

        except Exception as e:
            print(f"An error occurred while processing file {os.path.basename(file_path)}: {e}")
            import traceback
            traceback.print_exc()
            continue

    # 3. Concatenate and save results
    if not all_exitron_info:
        print("\nNo exitrons found across all files. No output file will be generated.")
        return None

    final_gr = pr.concat(all_exitron_info)
    print(f"\nSuccessfully compiled exitron data from {len(file_paths)} processed files.")
    final_gr.df.to_parquet(output_filepath, index=False)
    print(f"Successfully saved data to {output_filepath}")
    return final_gr

cDNA_unfil_exitrons = compileExitronData(
    "/gpfs/commons/home/ncui/project/falsitron_pipeline/",
    "/gpfs/commons/home/ncui/project/cDNA_unfil_exitrons.parquet",
    file_pattern="*cDNA.fixed.stranded.junc"
)

dRNA_unfil_exitrons = compileExitronData(
    "/gpfs/commons/home/ncui/project/falsitron_pipeline/",
    "/gpfs/commons/home/ncui/project/dRNA_unfil_exitrons.parquet",
    file_pattern="*dRNA.fixed.stranded.junc"
)
