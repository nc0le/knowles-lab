import pandas as pd
from pyfaidx import Fasta
import sys

def find_falsitrons(
    cDNA_path: str,
    dRNA_path: str,
    genome_fasta_path: str,
    min_cDNA_reads: int = 5,
    filter_splice_sites: bool = True
) -> pd.DataFrame:

    cDNA_df = pd.read_parquet(cDNA_path)
    dRNA_df = pd.read_parquet(dRNA_path)

    cDNA_counts = cDNA_df.groupby('title').size().reset_index(name='cDNA_junc_count')
    dRNA_counts = dRNA_df.groupby('title').size().reset_index(name='dRNA_junc_count')
    
    merged_df = pd.merge(cDNA_counts, dRNA_counts, on='title', how='outer')
    merged_df.fillna(0, inplace=True)
    merged_df['cDNA_junc_count'] = merged_df['cDNA_junc_count'].astype(int)
    merged_df['dRNA_junc_count'] = merged_df['dRNA_junc_count'].astype(int)
    
    # Filter for potential falsitrons based on read counts
    print(f"Filtering for candidates with cDNA reads > {min_cDNA_reads} and dRNA reads == 0...")
    potential_falsitrons = merged_df[
        (merged_df['cDNA_junc_count'] > min_cDNA_reads) &
        (merged_df['dRNA_junc_count'] == 0)
    ].copy()
    
    print(f"Found {len(potential_falsitrons)} candidates after read count filter.")

    genome = Fasta(genome_fasta_path)

    # Parse coordinates from the 'title' column
    coords = potential_falsitrons['title'].str.split(':', expand=True)
    potential_falsitrons['chrom'] = coords[0]
    potential_falsitrons['start'] = pd.to_numeric(coords[1])
    potential_falsitrons['end'] = pd.to_numeric(coords[2])
    potential_falsitrons['strand'] = coords[3]
    
    def get_splice_sites(row):
        try:
            donor = genome[row['chrom']][row['start']-1 : row['start']+1].seq.upper()
            acceptor = genome[row['chrom']][row['end']-2 : row['end']].seq.upper()
            if row['strand'] == '-':
                rc = str.maketrans('ATCG', 'TAGC')
                rev_donor = acceptor.translate(rc)[::-1]
                rev_acceptor = donor.translate(rc)[::-1]
                return f"{rev_donor}{rev_acceptor}"
            return f"{donor}{acceptor}"
        except (KeyError, ValueError):
            return "ERROR"

    potential_falsitrons['splice_site'] = potential_falsitrons.apply(get_splice_sites, axis=1)

    canonical_sites = ['GTAG', 'GCAG']
    df_valid = potential_falsitrons[potential_falsitrons['splice_site'] != 'ERROR']
    high_confidence_falsitrons = df_valid[~df_valid['splice_site'].isin(canonical_sites)].copy()

    # Clean up intermediate columns before returning
    high_confidence_falsitrons.drop(columns=['chrom', 'start', 'end', 'strand'], inplace=True)
    df_sorted = high_confidence_falsitrons.sort_values(by='cDNA_junc_count', ascending=False)

    return df_sorted


if __name__ == "__main__":
    
    cDNA_file = '/gpfs/commons/home/ncui/project/falsitron_pipeline/cDNA_unfil_exitrons.parquet'
    dRNA_file = '/gpfs/commons/home/ncui/project/falsitron_pipeline/dRNA_unfil_exitrons.parquet'
    GENOME_FASTA = '/gpfs/commons/home/ncui/project/reference_data/GRCh38.primary_assembly.genome.fa'
    output_file = '/gpfs/commons/home/ncui/project/falsitron_pipeline/final_falsitrons.csv'

    final_falsitrons = find_falsitrons(
        cDNA_path=cDNA_file,
        dRNA_path=dRNA_file,
        genome_fasta_path=GENOME_FASTA,
        min_cDNA_reads=5,
        filter_splice_sites=True
    )

    print(f"Found {len(final_falsitrons)} high-confidence falsitron candidates.")
    
    print("\nSample of final results:")
    print(final_falsitrons.head())

    final_falsitrons.to_csv(output_file, index=False)
    print(f"\nFinal results saved to {output_file}")
