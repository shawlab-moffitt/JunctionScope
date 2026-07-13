# 2026.06.24
- [update] updated junctionDB table - AO

# 2026.07.08
- [feature] QC step readout, is optional, is in config, will add extra time - AO
- [feature] velocyto function added, is optional, is in config, will add extra time - AO
- [feature] juncScope setup will now write out a slurm array batch job file - AO

- [update] config file not required, can enter required command on comman line during project setup - AO

- [bug] Fixed NA barcode lines being output from only raw IDs available - AO

- [to-do] Organize summarized result output - AO
- [to-do] ~~Allow velocyto loom files to be made externally and included in input - AO~~
- [to-do] ~~Derive intron/exon ratio per sample and per barcode from loom file - AO~~
- [to-do] Generate easy to load in R Seurat output - AO
- [to-do] Test scanpy single cell python library for clustering use - AO 

# 2026.07.10
- [update] Updated skeleton json config file with added parameters - AO
- [update] Updated README.md tutorial setup and text - AO
- [update] Conda environment now includes requirments for velocyto and scvelo steps - AO
- [update] Allow velocyto loom files to be made externally and included in input - AO
- [update] Derive intron/exon ratio per sample and per barcode from loom file - AO
- [update] Adjusted barcode summary merge of all samples to include intron/exon ratio columns if available - AO

- [bug] Edited the slurm script output to be an independent environment. You can now deploy the batch script while inside you conda environment - AO



