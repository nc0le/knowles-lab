## Summer 2025 Research @ Knowles Lab

### RNA-seq pipeline for detecting exitrons (exonic introns) 
pipeline_v1 overview: 

0. Install dependencies.
1. Parse raw junction data from regtools files.
2. Transform data by recalculating junction coordinates. Outputs data in bed12 format.
3. Function to find exitrons within data.
4. Compile all exitron data from all files.
5. Filter exitron data for exitrons present in >= 10 sources, with >= 30 reads each.
6.  

### ML to predict falsitrons (fake exitrons) 
