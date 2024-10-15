#!/usr/bin/bash

# usage from root repo: './data/test/get_data_test.sh'

cd data/test/
mkdir input/
cd input/

## Import genome
# Random hRSV_A fasta (assembled at chromosome level whithout annotation) => 'mapping_ref' for 'batch1'
mkdir ref/
wget -O ref/random.fa.gz "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/964/032/155/GCA_964032155.1_DW-RAT-268_NFW_4C_72h_162145/GCA_964032155.1_DW-RAT-268_NFW_4C_72h_162145_genomic.fna.gz"
gunzip ref/random.fa.gz
## RefSeq hRSV_A (fasta + gff) => 'mapping_ref' for 'batch2' and 'annot_ref' for 'batch1'
wget -O "ref/refseq.fa.gz" "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/815/375/GCF_002815375.1_ASM281537v1/GCF_002815375.1_ASM281537v1_genomic.fna.gz"
wget -O "ref/refseq.gff.gz" "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/815/375/GCF_002815375.1_ASM281537v1/GCF_002815375.1_ASM281537v1_genomic.gff.gz"
gunzip ref/refseq.fa.gz ref/refseq.gff.gz

## Reads and exepted results
git clone https://github.com/nexomis/genvar.git
# Generated
cd genvar/batch1/
python3 ../genvar.py ../../ref/random.fa config.yml 
cd ../batch2/
python3 ../genvar.py ../../ref/refseq.fa config.yml 
# Saved fq
cd ../../
mkdir -p reads/batch1/ reads/batch2/
mv genvar/batch1/*.fq.gz reads/batch1/
mv genvar/batch2/*.fq.gz reads/batch2/
# Saved exepted results
mkdir -p expected_results/batch1/ expected_results/batch2/
mv genvar/batch1/*.csv expected_results/batch1/
mv genvar/batch2/*.csv expected_results/batch2/

rm -rf genvar/

wget https://raw.githubusercontent.com/nexomis/db-kraken2-custom/refs/heads/main/build_custom_db.sh
bash build_custom_db.sh -i GCF_002815375.1 -n ASM281537v1 -j 162145 -o k2test

## RUN
#cd ../../../
#NXF_VER='24.04.4' nextflow run main.nf -profile docker --ref_ratio_threshold_variant 0.95
