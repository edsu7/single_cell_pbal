# single_cell_pbal

## Example for running 
```bash 
singularity exec -B \
<FASTQ DIRECTORY>:/fq,<OUT DIRECTORY>/:/out_dir,<REFERENCE DIRECTORY>:/ref \
<CONTAINER> \
python <SCRIPT> \
/out_dir/<INDEX.TSV> <JOB NAME> <REFERENCE>
```
```bash
fq="/home/usr/projects/PX0740/fastqs"
ref="/home/usr/references/"
out_dir="/home/usr/projects/PX0740/"
master_file="/home/usr/projects/PX0740/sample_list.tsv"
file=$(basename $master_file)

singularity exec -B \
$fq:/fq,$out_dir/:/out_dir,$ref:/ref \
/home/usr/software/PBAL_container/scPBAL_container.sif \
python run_example.py \
/out_dir/$file PX0740 hg38_no_alt
```
### Example of sample_list.tsv
#### File should be tab deliminated with no header
SC_IDENTIFIER|Read1|Read2
|------------ | -------------|-----
|AAACATA|/fq/PX0740_AAACATA/HFJMGCCXY_7_1_AAACATA_150bp.concat.fastq.gz|/fq/PX0740_AAACATA/HFJMGCCXY_7_2_AAACATA_150bp.concat.fastq.gz|
|AAAGCAA|/fq/PX0740_AAAGCAA/HFJMGCCXY_7_1_AAAGCAA_150bp.concat.fastq.gz|/fq/PX0740_AAAGCAA/HFJMGCCXY_7_2_AAAGCAA_150bp.concat.fastq.gz|
|AAATGCA|/fq/PX0740_AAATGCA/HFJMGCCXY_7_1_AAATGCA_150bp.concat.fastq.gz|/fq/PX0740_AAATGCA/HFJMGCCXY_7_2_AAATGCA_150bp.concat.fastq.gz|
|AACCCCA|/fq/PX0740_AACCCCA/HFJMGCCXY_7_1_AACCCCA_150bp.concat.fastq.gz|/fq/PX0740_AACCCCA/HFJMGCCXY_7_2_AACCCCA_150bp.concat.fastq.gz|
|AACTTGA|/fq/PX0740_AACTTGA/HFJMGCCXY_7_1_AACTTGA_150bp.concat.fastq.gz|/fq/PX0740_AACTTGA/HFJMGCCXY_7_2_AACTTGA_150bp.concat.fastq.gz|
|AAGACTA|/fq/PX0740_AAGACTA/HFJMGCCXY_7_1_AAGACTA_150bp.concat.fastq.gz|/fq/PX0740_AAGACTA/HFJMGCCXY_7_2_AAGACTA_150bp.concat.fastq.gz|
