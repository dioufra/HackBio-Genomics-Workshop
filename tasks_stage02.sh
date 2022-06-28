#!usr/bin/bash


# software installation
conda install -y -c bioconda \
bwa \
fastqc \
trimmomatic \
samtools \
bamtools\
snpeff \


# download the dataset	
mkdir -p raw_data 

wget -P raw_data/ https://zenodo.org/record/2582555/files/SLGFSK-N_231335_r1_chr5_12_17.fastq.gz
wget -P raw_data/ https://zenodo.org/record/2582555/files/SLGFSK-N_231335_r2_chr5_12_17.fastq.gz
wget -P raw_data/ https://zenodo.org/record/2582555/files/SLGFSK-T_231336_r1_chr5_12_17.fastq.gz
wget -P raw_data/ https://zenodo.org/record/2582555/files/SLGFSK-T_231336_r2_chr5_12_17.fastq.gz


# Quality control
mkdir Fastqc_Reports

fastqc raw_data/*.fastq.gz -o Fastqc_Reports
multiqc Fastqc_Reports -o Fastqc_Reports


# Trimming of low quality reads
mkdir -p trimmed_reads
mkdir -p trimmed_reads/Fastqc_Results

for sample in ${samples[@]}
do
    trimmomatic PE -threads 8 raw_data/${sample}_r1_chr5_12_17.fastq.gz raw_data/${sample}_r2_chr5_12_17.fastq.gz \
             trimmed_reads/${sample}_r1_paired.fq.gz trimmed_reads/${sample}_r1_unpaired.fq.gz \
             trimmed_reads/${sample}_r2_paired.fq.gz trimmed_reads/${sample}_r2_unpaired.fq.gz \
               ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:8:keepBothReads \
               LEADING:3 TRAILING:10 MINLEN:25
       
    fastqc  trimmed_reads/${sample}_r1_paired.fq.gz  trimmed_reads/${sample}_r2_paired.fq.gz \
                 -o trimmed_reads/Fastqc_Results
done

multiqc  trimmed_reads/Fastqc_results  -o trimmed_reads/Fastqc_results



#READS MAPPING
mkdir reference
mkdir Mapping

# download reference genome 
wget -P reference/ https://zenodo.org/record/2582555/files/hg19.chr5_12_17.fa.gz

# unzip the reference genome
gunzip reference/hg19.chr5_12_17.fa.gz

# index the reference genome
bwa index reference/hg19.chr5_12_17.fa

# perform alignment
bwa mem -R '@RG\tID:231335\tSM:Normal' reference/hg19.chr5_12_17.fa \
        trimmed_reads/SLGFSK-N_231335_r1_paired.fq.gz trimmed_reads/SLGFSK-N_231335_r2_paired.fq.gz > Mapping/SLGFSK-N_231335.sam

bwa mem -R '@RG\tID:231336\tSM:Tumor' reference/hg19.chr5_12_17.fa \
        trimmed_reads/SLGFSK-T_231336_r1_paired.fq.gz trimmed_reads/SLGFSK-T_231336_r2_paired.fq.gz > Mapping/SLGFSK-T_231336.sam



# Conversion of the SAM file to BAM file, sorting and indexing
for sample in ${samples[@]}
do
    #convert sam to bam and sort it by name
    samtools view -@ 20 -S -b Mapping/${sample}.sam | samtools sort -n -@ 32 > Mapping/${sample}.namesorted.bam

    # index bam files
    samtools index Mapping/${sample}.namesorted.bam
done


# Filtering mapped reads
for sample in ${samples[@]}
do
    # filter bam files
    samtools view -q 1 -f 0x2 -F 0x8 -b Mapping/${sample}.sorted.bam > Mapping/${sample}.filtered1.bam
done


# Remove duplicates
for sample in ${samples[@]}
do
    samtools collate Mapping/${sample}.filtered1.bam Mapping/${sample}.namecollate
    samtools fixmate -m Mapping/${sample}.namecollate.bam Mapping/${sample}.fixmate.bam
    samtools sort -@ 32 -o Mapping/${sample}.positionsort.bam Mapping/${sample}.fixmate.bam
    samtools markdup -@ 32 -r Mapping/${sample}.positionsort.bam Mapping/${sample}.clean.bam
done


# letf-align bam files.
for sample in ${samples[@]}
do
    cat Mapping/${sample}.clean.bam | bamleftalign -f hg19.chr5_12_17.fa -m 5 -c > Mapping/${sample}.leftAlign.bam
done


# recalibrate reads mapping qualities.
for sample in ${samples[@]}
do
    samtools calmd -@ 32 -b Mapping/${sample}.leftAlign.bam hg19.chr5_12_17.fa > Mapping/${sample}.recalibrate.bam
done

# refiltter reads mapping qualities
for sample in ${samples[@]}
do
    bamtools filter -in Mapping/${sample}.recalibrate.bam -mapQuality "<=254" > Mapping/${sample}.refilter.bam
done


# VARIANT CALLING AND CLASSIFICATION

mkdir Variants

# convert data to pileup
for sample in ${samples[@]}
do
	samtools mpileup -f hg19.chr5_12_17.fa Mapping/${sample}.refilter.bam --min-MQ 1 --min-BQ 28 \
                > Variants/${sample}.pileup
done

#Call variants
varscan somatic Variants/SLGFSK-N_231335.pileup \
        Variants/SLGFSK-T_231336.pileup Variants/SLGFSK \
        --normal-purity 1  --tumor-purity 0.5 --output-vcf 1 

# merge output vcf files.
bgzip Variants/SLGFSK.snp.vcf > Variants/SLGFSK.snp.vcf.gz
bgzip Variants/SLGFSK.indel.vcf > Variants/SLGFSK.indel.vcf.gz
tabix Variants/SLGFSK.snp.vcf.gz
tabix Variants/SLGFSK.indel.vcf.gz
bcftools merge Variants/SLGFSK.snp.vcf.gz Variants/SLGFSK.indel.vcf.gz > Variants/SLGFSK.vcf


# VARIANT ANNOTATION

#download jar file if snpeff is not installed.
#wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip

# Unzip file
#unzip snpEff_latest_core.zip

#download snpEff database
#java -jar snpEff.jar download hg19

# functionnal annotation with snpeff
snpEff hg19 Variants/SLGFSK.vcf > Variants/SLGFSK.ann.vcf

# clinical annotation with gemini

# install gemini
wget https://raw.github.com/arq5x/gemini/master/gemini/scripts/gemini_install.py
python gemini_install.py /usr/local /usr/local/share/gemini


# load the variant.
gemini load -v Variants/SLGFSK.ann.vcf -t snpEff Annotation/gemini.db