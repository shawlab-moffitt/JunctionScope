# junctionScope

# Requirments
## Version Control

## Input Files

### Sample table
- Tab-delim text file of sample name and cram/bam file.
- Optional inclusion of loom file as third column


|  |  |  |
| --- | --- | --- |
| Sample1 | path/to/sample1.bam | path/to/sample1.loom |
| Sample2 | path/to/sample2.bam | path/to/sample2.loom |


### Junction Table
- Tab-delim text file of junctions to quantify.
- Header is optional
- Nucleotide sequence column is optional
- Gene column should be formatted as the HUGO gene symbol
- junction name ID should be unique to junction


| junction_name | gene | junction | nt_seq |
| --- | --- | --- | --- |
| junctionName1 | gene1 | chr1:1000-2000 | ...ATCGATCG... | 
| junctionName2 | gene2 | chr2:1000-2000 | ...ATCGATCG... |


## Environment Setup
### Create conda environment

```bash
ml Anaconda3/2024.02-1
# Updated conda environment velocyto compatible
conda env create -f juncScope_velocyto_scvelo.yml
conda activate juncScope_velocyto_scvelo
```

## Run Single Sample
- Arguments in square brackets are optional
```bash
# With config file
python junctionScope.py -c <txt|json config file> --sample <sample name> --bam <cram|bam file>

# OR #

# Without config file
python junctionScope.py --sample <sample name> --bam <cram|bam file> -j <junction table> --mode <singlecell|bulk> [--loom <loom file>] [--qc] [--velocyto]
```

## Multi-Sample Run
### Setup project
Including the argument `--setup-only` will prepare the project folders for execution.
If the `--setup-only` argument is not included the script with setup and run immediately.
```bash
python junctionScope.py -c <txt|json config file> --setup-only
```
The setup will generate the project folder and scripts that need to be run.
It will output the suggested commands to run next.

### Run
Users will be provided a command that can perform interactive job execution which can be run via `sh <juncScope_exc.sh script>`

Or the setup will generate a ready-made slurm array job script that can be submitted to a computing cluster with `sbatch <juncScope_slurm_array.sh script>`

### Post-Process
#### Summarize
The setup function will also provide the command to perform the project summarization. This can be run after all of the samples have completed processing.
```bash
python path/to/output/project/junctScopeSummarize.py
```



