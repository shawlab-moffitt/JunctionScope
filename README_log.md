# 2026.06.24
- [update] updated junctionDB table - AO

# 2026.07.08
- [feature] QC step readout, is optional, is in config, will add extra time - AO
- [feature] velocyto function added, is optional, is in config, will add extra time - AO
- [feature] juncScope setup will now write out a slurm array batch job file - AO

- [update] config file not required, can enter required command on comman line during project setup - AO

- [bug-fix] Fixed NA barcode lines being output from only raw IDs available - AO

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

- [bug-fix] Edited the slurm script output to be an independent environment. You can now deploy the batch script while inside you conda environment - AO

# 2026.07.15
- [bug] ~~Running main script wthout `--setup-only` argument wrongly assumes running only one sample~~
	~~Error Text: TypeError: main.<locals>._run_one() missing 1 required positional argument: 'loom_file' - AO~~

- [to-do] ~~Make output after setup more clear on run options - AO~~
- [to-do] Check lst file format, see if full direct path required - AO
- [to-do] Simplify setup commands, maybe make default settings for easier running - AO

# 2026.07.22
- [bug-fix] Running main script without `--setup-only` argument wrongly assumes running only one sample
	Fixed. User can now run a list of bam files via config from initial command and bypass setup - AO
- [bug-fix] Found and fixed issue when chromosome name may not include 'chr' prefix - AO
- [bug-fix] Fixed 'bulk' mode summary functionality - AO

- [update] Added spacing in setup output on commands for next steps - AO
- [update] Removed loom stat output checking from verbose, was only temporary to debug - AO
- [update] read_loom duplicate variable names warning is now hidden and reformatted to annoatte the issue has been address in the code - AO

- [idea] activate 'debug' mode

