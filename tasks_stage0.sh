#!/usr/bin/bash
firstname="François Pape"
lastname="DIOUF"

# Print the firtname and lastname on the same line.
echo $firstname $lastname


# Print the first anme and the second name on different lines.
echo $firstname
echo $lastname

# Bash story one

# step1 create a directory named team_newton
mkdir francois

# step2 create another directory named biocomputing 
mkdir biocomputing && cd $_

# step3 downloas the three files
wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.fna
wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.gbk
wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.gbk

# step4 move the .fna file to the folder titled team_newton.
mv wildtype.fna ../francois

# step5 delete the duplicated file
rm wildtype.gbk.1

# step6 confirm if the file is mutant or wildtype
grep "tatatata" ../francois/wildtype.fna

# step7 print all the line
grep "tatatata" ../francois/wildtype.fna > mutant.txt

# step8 clear terminal and print all the command used today
clear && history

# step9 list the files in the two folders
ls ~/francois ~/biocomputing 


# bash story two

# step1 graphical representation

#dwnload the figlet command
sudo apt-get install figlet
figlet "François Pape Diouf"

#step3 create a folder named compare in our home directory
cd ~ 
mkdir compare

#step 3-a download the file
cd compare
wget https://www.bioinformatics.babraham.ac.uk/training/Introduction%20to%20Unix/unix_intro_data.tar.gz

# b unzip the file 
gunzip unix_intro_data.tar.gz

# c untar the file 
tar -xvf unix_intro_data.tar

# d  indentify the rRNAs present in Mito.dat
cd "seqmonk_genomes/Saccharomyces cerevisiae/EF4"
grep "rRNA" Mito.dat

# copy Mito.dat to compare 
cp Mito.dat ~/compare/

# f-i change Mito to Mitochondrion in the ID and AC header lines
nano Mito.dat

#  f-ii save the file using ctrl + s and exit with ctrl + x

# f-iii rename Mito.dat to Mitichondrion
cp Mito.dat Mitochondrion.txt

# step4 move to compare and cd to FastQ_Data directory
cd ~/compare/
cd FastQ_Data 

# calculate the total number of lines in lane8_DD_P4_TTAGGC_L008_R1.fastq.gz
cat lane8_DD_P4_TTAGGC_L008_R1.fastq.gz | wc -l

# Print the total number of lines in all fastq.gz files and save it as a new file.
cat *.fastq.gz | wc -l > number_of_lines.txt
