# viral-variant

## Description

Calling variants with peptide impact and comparison between samples (via transfer annotation): From assembled genome, reference genome with its annotation and mapping of reads on assembled genome.


## Test and check results 


```

python3 data/test/check_results.py data/test/input/expected_results/batch1/P0.indel.csv data/test/input/expected_results/batch1/P0.snp.csv data/test/out_dir/var_batch_filtered/batch1_batchFiltered/b1_P0_corrected_batchFiltered.tsv data/test/input/ref/random.fa
python3 data/test/check_results.py data/test/input/expected_results/batch1/P1.indel.csv data/test/input/expected_results/batch1/P1.snp.csv data/test/out_dir/var_batch_filtered/batch1_batchFiltered/b1_P1_corrected_batchFiltered.tsv data/test/input/ref/random.fa

python3 data/test/check_results.py data/test/input/expected_results/batch2/P0.indel.csv data/test/input/expected_results/batch2/P0.snp.csv data/test/out_dir/var_batch_filtered/batch2_batchFiltered/b2_P0_corrected_batchFiltered.tsv data/test/input/ref/refseq.fa
python3 data/test/check_results.py data/test/input/expected_results/batch2/P1.indel.csv data/test/input/expected_results/batch2/P1.snp.csv data/test/out_dir/var_batch_filtered/batch2_batchFiltered/b2_P1_corrected_batchFiltered.tsv data/test/input/ref/refseq.fa
python3 data/test/check_results.py data/test/input/expected_results/batch2/P2.indel.csv data/test/input/expected_results/batch2/P2.snp.csv data/test/out_dir/var_batch_filtered/batch2_batchFiltered/b2_P2_corrected_batchFiltered.tsv data/test/input/ref/refseq.fa
python3 data/test/check_results.py data/test/input/expected_results/batch2/P3.indel.csv data/test/input/expected_results/batch2/P3.snp.csv data/test/out_dir/var_batch_filtered/batch2_batchFiltered/b2_P3_corrected_batchFiltered.tsv data/test/input/ref/refseq.fa

batch=1
passage=0
pos=434
ref=data/test/input/ref/random.fa
region=OZ035168.1
samtools tview -p $region:$pos -d T -w 30 data/test/out_dir/aln/b${batch}_P${passage}_abra2.bam $ref > raw.txt ; head -n 2 raw.txt | uniq -c ; cat raw.txt | sed '1,2d' | sort | uniq -c | sort -nr

```
