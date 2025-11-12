junctionScope

# Requirments
Built with SAMtools (v1.9) and RegTools (v0.5.2).

# Basic Inputs
![junctionScope_basicInput](https://github.com/user-attachments/assets/1f11269e-5d58-4076-b9e9-c19a3dbca042)

## Config File Setup

```
# junctionSope Config
INPUT="alignment_sample_files.lst"
REGION="chr1:1000:2000"
NTSEQ="ATCGGACATT"
FASTA="ref.fasta"
GTF="ref.gtf"
OUTPUT="junctionScope_Output_Name"
BUFFER=200
STRAND=0
THREADS=8
```

```
sh junctionScope.sh -c junctionScope.conf
sh junctionScope_Output_Name_exc.sh
```

## Command Line Run

```
sh junctionScope.sh -i alignment_sample_files.lst -r chr1:1000:2000 -n ATCGGACATT -f ref.fasta -g ref.gtf -o junctionScope_Output_Name -b 200 -s 0 -t 8
sh junctionScope_Output_Name_exc.sh
```


# Workflows
## samtools method
![junctionScope_samtoolsWorkflow](https://github.com/user-attachments/assets/38ac9345-56b6-47a2-9e1b-c0445cbcc7c4)

## regtools method
![junctionScope_regtoolsWorkflow](https://github.com/user-attachments/assets/7ea9bca2-d4e4-4169-9bbe-783515780a95)
