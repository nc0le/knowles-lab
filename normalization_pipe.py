import pandas as pd
from tqdm import tqdm # progress tracker
import pyranges as pr # parsing gff
import numpy as np
import pysam

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

filtered_exitron_data = filterExitronData('final_exitron_data.parquet')

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
    normalized_df.to_parquet(output_filepath, index=False)
    return normalized_df

normalizeExitronData(filtered_exitron_data, 'normalized_data.parquet')
