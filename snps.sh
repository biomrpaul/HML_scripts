#!/bzin/bash

FILE1="$1"

FILE2="$2"
sample_name="$3"


mkdir inter_files/$sample_name.dir

echo  "Trimming Files"
java -jar /mnt/gv1/software/tool-deps/Trimmomatic-0.30/trimmomatic-0.30.jar PE -threads 4 \
-phred33 -trimlog inter_files/$sample_name.dir/trimming.log $FILE1 $FILE2 inter_files/$sample_name.dir/$sample_name.trimmed_f inter_files/$sample_name.dir/uo1 \
inter_files/$sample_name.dir/$sample_name.trimmed_r inter_files/$sample_name.dir/uo2 \
ILLUMINACLIP:/mnt/gv1/software/tool-deps/Trimmomatic-0.30/adapters/TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 \
SLIDINGWINDOW:4:20 MINLEN:36


#aln fastq files 
echo "Creating index for forward file"
sudo bwa aln -n 0.04 -o 1 -e -1 -i 5 -d 10 -l 32 -k 2 -m 2000000 -t 1 -M 3 -O 11 -E 4 \
-R 30 -q 0 annot_files/SNPiR.fa inter_files/$sample_name.dir/$sample_name.trimmed_f > inter_files/$sample_name.dir/$sample_name.trimmed_f.sai

echo "Creating index for reverse file"
sudo bwa aln -n 0.04 -o 1 -e -1 -i 5 -d 10 -l 32 -k 2 -m 2000000 -t 1 -M 3 -O 11 -E 4 \
-R 30 -q 0 annot_files/SNPiR.fa inter_files/$sample_name.dir/$sample_name.trimmed_r > inter_files/$sample_name.dir/$sample_name.trimmed_r.sai

#bwa forward and reverse
echo "Aligning trimmed files to SNPiR.fa"
sudo bwa mem annot_files/SNPiR.fa inter_files/$sample_name.dir/$sample_name.trimmed_f inter_files/$sample_name.dir/$sample_name.trimmed_r \
> inter_files/$sample_name.dir/$sample_name.first.sam

#Convert Coordinates
echo "Converting coordinates of aligned file to SNPiR format"
sudo java -Xmx2g SNPiR/convertCoordinates < inter_files/$sample_name.dir/$sample_name.first.sam > inter_files/$sample_name.dir/$sample_name.converted.sam

#Convert sam files to bam files
echo

#Remove duplicates



