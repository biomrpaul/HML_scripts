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

#Add read group info (required for GATK)
java -Xmx2g -jar software/picard-tools-1.110/AddOrReplaceReadGroups.jar INPUT=inter_files/$sample_name.dir/$sample_name.filtered.bam OUTPUT=inter_files/$sample_name.dir/$sample_name.filtered_w_RG.bam RGLB=$sample_name RGPL=Illumina RGPU=Lane1 RGSM=$sample_name VALIDATION_STRINGENCY=SILENT

#reordering bam
java -jar software/picard-tools-1.110/ReorderSam.jar INPUT=inter_files/$sample_name.dir/$sample_name.filtered_w_RG.bam output=inter_files/$sample_name.dir/$sample_name.sorted.bam REFERENCE=annot_files/SNPiR.fa VALIDATION_STRINGENCY=SILENT

#filtering unmapped reads
sudo samtools view -b -F 4  inter_files/$sample_name.dir/$sample_name.sorted.bam > inter_files/$sample_name.dir/$sample_name.filtered.bam

#create index for bam
samtools index inter_files/$sample_name.dir/$sample_name.filtered.bam

#create realigner target
java -jar software/gatk/GenomeAnalysisTK.jar -T RealignerTargetCreator  -nct 24 -R annot_files/SNPiR.fa  -I inter_files/$sample_name.dir/$sample_name.filtered.bam  -o inter_files/$sample_name.dir/$sample_name.intervals -U ALLOW_N_CIGAR_READS

#indel realigner
java -jar software/gatk/GenomeAnalysisTK.jar -T IndelRealigner  -nct 24 -R annot_files/SNPiR.fa  -I inter_files/$sample_name.dir/$sample_name.filtered.bam -targetIntervals inter_files/$sample_name.dir/$sample_name.intervals -o inter_files/$sample_name.dir/$sample_name.indel.aligned.bam -U ALLOW_N_CIGAR_READS --maxReadsForRealignment 100000

#BaseRecalibrator(creates grp file)
java -Xmx4g -jar software/gatk/GenomeAnalysisTK.jar  -T BaseRecalibrator -nct 24 -I inter_files/$sample_name.dir/$sample_name.indel.aligned.bam  -R annot_files/SNPiR.fa  -U ALLOW_N_CIGAR_READS   -knownSites annot_files/dbsnp_hg19.vcf -o inter_files/$sample_name.dir/$sample_name.grp

#Print Reads (creates recalibrated bam from grp file)
java -jar software/gatk/GenomeAnalysisTK.jar -T PrintReads -nct 24 -R annot_files/SNPiR.fa -I inter_files/$sample_name.dir/$sample_name.indel.aligned.bam -U ALLOW_N_CIGAR_READS -BQSR inter_files/$sample_name.dir/$sample_name.grp -o inter_files/$sample_name.dir/$sample_name.final.bam

#UnifiedGenotyper
java -jar software/gatk/GenomeAnalysisTK.jar -T UnifiedGenotyper  -nct 24 -R annot_files/genome.fa -I inter_files/$sample_name.dir/$sample_name.final.bam --dbsnp annot_files/dbsnp_hg19.vcf -o inter_files/$sample_name.dir/$sample_name.snps.raw.vcf -stand_call_conf 0 -stand_emit_conf 0  -U ALLOW_N_CIGAR_READS --output_mode EMIT_VARIANTS_ONLY -glm SNP

#Convert vcf to SNPiR format (using SNPiR tool)
sh SNPiR/convertVCF.sh inter_files/$sample_name.dir/$sample_name.snps.raw.vcf inter_files/$sample_name.dir/$sample_name.snpir.txt 20

#filter mismatches at read end
perl SNPiR/filter_mismatch_first6bp.pl -illumina -infile inter_files/$sample_name.dir/$sample_name.snpir.txt -outfile inter_files/$sample_name.dir/$sample_name.snpir.trimmed.txt -bamfile inter_files/$sample_name.dir/$sample_name.final.bam

#Remove variants that lie in repetitive regions
awk '{OFS="\t";$2=$2-1"\t"$2;print $0}' inter_files/$sample_name.dir/$sample_name.snpir.trimmed.txt | intersectBed -a stdin -b annot_files/hg19_rmsk -v | cut -f1,3-7 > inter_files/$sample_name.dir/$sample_name.snpir.rmsk.txt

#Filter variants in intronic regions (on this step)
perl SNPiR/filter_intron_near_splicejuncts.pl -infile inter_files/$sample_name.dir/$sample_name.snpir.rmsk.txt -outfile inter_files/$sample_name.dir/$sample_name.snpir.rmsk.rmintron.txt -genefile annot_files/gene_file_sorted

#filter variants in homopolymers
perl SNPiR/filter_homopolymer_nucleotides.pl -infile inter_files/$sample_name.dir/$sample_name.snpir.rmsk.rmintron.txt -outfile inter_files/$sample_name.dir/$sample_name.snpir.rmsk.rmintron.rmhom.txt -refgenome annot_files/genome.fa

#filter variants that were caused by mismapped reads
perl SNPiR/BLAT_candidates.pl -infile inter_files/$sample_name.dir/$sample_name.snpir.rmsk.rmintron.rmhom.txt -outfile inter_files/$sample_name.dir/$sample_name.snpir.rmsk.rmintron.rmhom.rmblat.txt -bamfile inter_files/$sample_name.dir/$sample_name.final.bam -refgenome annot_files/SNPiR.fa

mkdir variants

#remove known RNA editing sites
awk '{OFS="\t";$2=$2-1"\t"$2;print $0}' inter_files/$sample_name.dir/$sample_name.snpir.rmsk.rmintron.rmhom.rmblat.txt | intersectBed -a stdin -b annot_files/Human_AG_all_hg19.bed -v > variants/$sample_name.variants.bed 
