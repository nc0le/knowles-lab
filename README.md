## Summer 2025 Research @ Knowles Lab

### Detecting exitrons 
<img width="733" height="320" alt="Screenshot 2025-07-22 at 10 23 01 AM" src="https://github.com/user-attachments/assets/a61ae84e-570c-4deb-a2c3-03fab4140af6" />

#### Usage
```python exitron_pipeline.py```
#### Overview
1. Parse raw junction data from regtools files.
2. Transform data by recalculating junction coordinates. Outputs data in bed12 format.
3. Extract annotated CDS regions from gencode v48.
4. Find exitrons *contained* within CDS regions. Reference ____ docs for intron spanning exitrons. 
6. Filter exitron data for exitrons present in >= 10 sources, with >= 30 reads each.
7. Normalize exitron data by finding exitron expression proportionate to expression of its surrounding regions.
### Extended pipeline for detecting exitrons and classifying falsitrons 
