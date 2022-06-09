#!usr/bin/bash

# download DNA.fa the file
wget https://raw.githubusercontent.com/HackBio-Internship/wale-home-tasks/main/DNA.fa

#count the number of sequences
grep -c "^>" DNA.fa

#Get the total A, T, G & C counts for all the sequences in the file above
grep -Eo 'A|T|G|C' DNA.fa | sort | uniq -c | awk '{print $2": "$1}'

#Set up a conda environment
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

# download required sofwares(fastqc, fastp, bbtools, bwa)
conda install \
fastqc \
fastp \
bbtools \
bwa

# create a working directory
mkdir workdir && cd $_

# creata a dataset directory and download the required datasets.
mkdir dataset && cd $_
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/ACBarrie_R1.fastq.gz?raw=true -O ACBarrie_R1.fastq.gz
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/ACBarrie_R2.fastq.gz?raw=true -O ACBarrie_R2.fastq.gz
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/Alsen_R1.fastq.gz?raw=true -O Alsen_R1.fastq.gz
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/Alsen_R2.fastq.gz?raw=true -O Alsen_R2.fastq.gz
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/Baxter_R1.fastq.gz?raw=true -O Baxter_R1.fastq.gz
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/Baxter_R2.fastq.gz?raw=true -O Baxter_R2.fastq.gz
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/Chara_R1.fastq.gz?raw=true -O Chara_R1.fastq.gz
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/Chara_R2.fastq.gz?raw=true -O Chara_R2.fastq.gz
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/Drysdale_R1.fastq.gz?raw=true -O Drysdale_R1.fastq.gz
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/Drysdale_R2.fastq.gz?raw=true -O Drysdale_R2.fastq.gz


# create an output folder: all the results will be sent to that folder
cd ~/workdir
mkdir output

# Perform quality control with fastqc

## Create a QC_report folder under outout and sent the output to it.
QC_DIR=~/workdir/output/QC_Report
mkdir -p $QC_DIR

fastqc ~/workdir/dataset/*.fastq.gz -o $QC_DIR


# Trim the low quality reads
TRIMMED_DIR=~/workdir/output/trimmed_data # send all the results under trimmed reads under trimmed_data folder.
DATA_DIR=~/workdir/dataset
mkdir -p $TRIMMED_DIR # make  directory trimmed

samples="ACBarrie Alsen Baxter Chara Drysdale"

for sample in $samples;
do
    fastp \
    -i $DATA_DIR/${sample}_R1.fastq.gz \
    -I $DATA_DIR/${sample}_R2.fastq.gz \
    -o $TRIMMED_DIR/${sample}_trimmed_R1.fastq.gz \
    -O $TRIMMED_DIR/${sample}_trimmed_R2.fastq.gz \
    --html $TRIMMED_DIR/${sample}_fastp.html
done


# Perform the alignment with bwa

## download the reference genome
REFERENCE_DIR=~/workdir/ref_genome
mkdir $REFERENCE_DIR
cd $REFERENCE_DIR
wget wget https://github.com/josoga2/yt-dataset/raw/main/dataset/references/reference.fasta


## index the reference genome
bwa index $REFERENCE_DIR/reference.fasta


## repairing desordered reads and performing alignment
cd ~/workdir
ALIGNEMENT_DIR=~/workdir/output/alignment
REPAIRED_DIR=~/workdir/output/repaired
mkdir -p $REPAIRED_DIR
mkdir -p $ALIGNEMENT_DIR

for sample in $samples;
do
    repair.sh in1=$TRIMMED_DIR/${sample}_trimmed_R1.fastq.gz in2=$TRIMMED_DIR/${sample}_trimmed_R2.fastq.gz out1=$REPAIRED_DIR/${sample}_R1_rep.fastq.gz out2=$REPAIRED_DIR/${sample}_R2_rep.fastq.gz outsingle=$REPAIRED_DIR/${sample}_single.fastq
    bwa mem -r 1 \
    $REFERENCE_DIR/reference.fasta \
    $REPAIRED_DIR/${sample}_R1_rep.fastq.gz $REPAIRED_DIR/${sample}_R2_rep.fastq.gz > $ALIGNEMENT_DIR/${sample}.sam
done