# Summer 2025 Research @ Knowles Lab

## Falsitron prediction with machine learning (currently)
<img width="677" height="186" alt="Screenshot 2025-07-22 at 11 00 39 AM" src="https://github.com/user-attachments/assets/9538b16f-2efa-4e21-8964-a01f86076c23" />

## Extended pipeline for detecting exitrons and falsitrons

### Overview of ext_exitron_pipe.py:
1. Parse raw junction data from regtools files.
2. Transform data by recalculating junction coordinates. Outputs data in bed12 format.
   
## Pipeline for detecting exitrons and calculating expression levels (for ALS analysis)
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
  
### Overview of exitron_pipe.py:
1. Parse raw junction data from regtools files.
2. Transform data by recalculating junction coordinates. Outputs data in bed12 format.
3. Extract annotated CDS regions from gencode v48.
4. Find exitrons *contained* within CDS regions. Reference ext_exitron_pipe.py docs for intron spanning exitrons. 
6. Filter exitron data for exitrons present in >= 10 sources, with >= 30 reads each.
7. Normalize exitron data by finding exitron expression proportionate to expression of its surrounding regions.

## ALS analysis
