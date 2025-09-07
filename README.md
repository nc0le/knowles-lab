## 🌱 Falsitron prediction with machine learning (currently)
<img width="677" height="186" alt="Screenshot 2025-07-22 at 11 00 39 AM" src="https://github.com/user-attachments/assets/9538b16f-2efa-4e21-8964-a01f86076c23" />

## Extended pipeline for detecting exitrons and falsitrons

### Usage
Clone ext_exitron_pipe.py and falsitron_finder.py
```
python ext_exitron_pipe.py
```
```
python falsitron_finder.py
```
### Pipeline overview:
1. Force stradedness in regtools data to accomodate nanopore data ('?' -> '-').
2. Parse raw junction data from regtools files.
3. Transform data by recalculating junction coordinates. Outputs data in bed12 format.
4. Extract following features from gencode v48: CDS, exon, five_prime_UTR, three_prime_UTR.
5. Find exitrons *contained* within these regions.
6. Find exitrons spanning adjacent exons.

FINDING FALSITRONS:

7. Use pipeline above to find exitrons for cDNA and dRNA data.
8. Generate matrix of exitrons vs cDNA junction counts (reads) and dRNA junction counts.
9. Filter for exitrons with cDNA counts > min_cDNA_reads (default = 5) and zero dRNA counts.
10. Filter out canonical splice sites.

## 🌱 Pipeline for detecting exitrons and calculating expression levels (for ALS analysis)
<img width="733" height="320" alt="Screenshot 2025-07-22 at 10 23 01 AM" src="https://github.com/user-attachments/assets/a61ae84e-570c-4deb-a2c3-03fab4140af6" />

### Usage
Clone exitron_pipe.py and adjust file locations found at end of script
```
# change first argument to directory containing junction data 
# change second argument to parquet file for exitron data output 
# change third argument to correct junction file identifier
compileExitronData("/gpfs/commons/groups/knowles_lab/atokolyi/als/juncs_min6bp/", "/gpfs/commons/home/ncui/project/exitron_data.parquet", file_pattern="*.bam.junc")

# change argument to name of filtered exitron data parquet file 
filtered_exitron_data = filterExitronData('final_exitron_data.parquet')

# change argument to name of normalized exitron expression data parquet file 
normalizeExitronData(filtered_exitron_data, 'normalized_data.npy')
```
Run script
```
python exitron_pipe.py
```
exitron_pipe.py generates three files:
* Parquet file of all exitrons compiled from raw junction data, contained within CDS regions of gencode v38
* Parquet file of unique exitrons filtered for significance (default: present in at least 10 sources, with at least 30 reads each)
* .npy file of all unique exitrons and their normalized expression levels
  
### Pipeline overview:
1. Parse raw junction data from regtools files.
2. Transform data by recalculating junction coordinates. Outputs data in bed12 format.
3. Extract annotated CDS regions from gencode v48.
4. Find exitrons *contained* within CDS regions. Reference ext_exitron_pipe.py docs for intron spanning exitrons. 
6. Filter exitron data for exitrons present in >= 10 sources, with >= 30 reads each.
7. Normalize exitron data by finding exitron expression proportionate to expression of its surrounding regions.

## 🌱 ALS analysis
### Modeling Exitron Usage in ALS vs Healthy Individuals 
ALSvsHealthy_model.ipynb
#### Key observations:
* Strong tissue specificity for frontal cortex
* Negative coefficients indicate these exitrons are **less used** in individuals with ALS
#### Top results:
```
  Coef.  Std.Err.         z         P>|z|     [0.025    0.975]          tissue                     exitron     p_fdr
-2.660590  0.519358 -5.122847  3.009560e-07  -3.678512 -1.642668  Cortex_Frontal   chr19:35787559:35787643:+  0.000130
-1.527758  0.299147 -5.107042  3.272404e-07  -2.114077 -0.941440  Cortex_Frontal   chr20:63678183:63678267:+  0.000130
-7.592406  1.486782 -5.106603  3.280017e-07 -10.506445 -4.678367  Cortex_Frontal    chr7:83155079:83155169:-  0.000130
-4.972107  1.012073 -4.912797  8.978622e-07  -6.955733 -2.988481  Cortex_Frontal  chr1:240207870:240208101:+  0.000267
-4.962993  1.087790 -4.562456  5.055875e-06  -7.095022 -2.830964  Cortex_Frontal  chr1:240207639:240208134:+  0.000893
```
### Modeling Effect of Exitrons on ALS Severity 
ALS_severity_model.ipynb

#### Key observations:
* Did not find significant exitrons in their correlation with c9orf72 status
