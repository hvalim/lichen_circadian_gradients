---
title: 'Part 1: pool-seq with popoolation2'
author: "Francesco Dal Grande, Henrique Valim"
date: "2022-12-16"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

Populations used:

Umbilicaria pustulata:
- Upust_IT_Pool1-6: Mt. Limbara in Sardinia, Italy 
- Upust_ESi_Pool1-3: Sierra de Gredos, Spain
- Upust_ESii_Pool1-6: Puerto de Pico, Spain

U. phaea:
- Uph16-19: Sierra Nevada
- Uph22-28: Mt. Jacinto 

# Introduction: running Pool-seq data analysis with popoolation2

Tools needed to be installed for this process are:

bwa: conda install -c bioconda bwa
popoolation2: download and unzip popoolation2 from here: https://sourceforge.net/projects/popoolation2/
picard tools: conda install -c bioconda picard
samtools: conda install -c bioconda samtools


# 1. Initial cleanup and filtering of raw reads

```{bash Example: cleanup and filtering}

nohup perl ~/Codes/Supplementary_FIle1-FastQFS.pl.pl -filtering Yes -fw P2_cleaned_R1_paired.fq -rw P2_cleaned_R2_paired.fq -fwo ESii_Pool2_Q26_L80_BQ3_R1_paired.fastq -rwo ESii_Pool2_Q26_L80_BQ3_R2_paired.fastq -sngl ESii_Pool2_Q26_L80_BQ3_singletons.fastq -mq 3 -q 26 -l 80 -plotting No -gsize 150 &

```

# 2. Create genome index

The only step needed is to copy the genome fasta file and index it. 

We choose here Umbilicaria_phaea_Uphaea_ref_cold.scaffolds.fa as our genome file, as U_phaea_TBG_1112 (high elevation/cold temperate) is the currently published U. phaea genome:

```{bash Genome index with bwa}

bwa index Umbilicaria_phaea_Uphaea_ref_cold.scaffolds.fa 
bwa index Lpus_4dec_AN3.masked.fasta

```

Therefore, we can skip the usual next step, which is generating the interleaved reads. An example of the code is included here for completion:

```{bash Example: creating interleaved reads}

perl /opt/omega/scripts/shuffleSequences_fastq.pl Uph16_paired.R1.fastq Uph16_paired.R2.fastq Uph16_paired_interleaved.fq &

perl /opt/omega/scripts/shuffleSequences_fastq.pl ESii_Pool1_Q26_L80_BQ3_paired.tr_1.fastq ESii_Pool1_Q26_L80_BQ3_paired.tr_2.fastq ESii_Pool1_Q26_L80_BQ3_paired_interleaved.fq 

```


# 3. Add group reads to bam headers to be used in GATK

We can then proceed to the next step, which is to add the interleaved reads to bam headers. 

```{bash Create sam headers from interlieved reads}

# U. phaea

bwa mem -M -t 12 -R '@RG\tID:group1\tSM:Uph16\tPL:illumina\tLB:lib1\tPU:Uph16' -p Umbilicaria_phaea_Uphaea_ref_cold.scaffolds.fa Uph16_FDSW202335553-1a_H33VHDSXY_L4_paired_interleaved.fq > Uph16_aligned_reads.sam 
bwa mem -M -t 12 -R '@RG\tID:group1\tSM:Uph17\tPL:illumina\tLB:lib1\tPU:Uph17' -p Umbilicaria_phaea_Uphaea_ref_cold.scaffolds.fa Uph17_FDSW202335554-1a_H33VHDSXY_L4_paired_interleaved.fq > Uph17_aligned_reads.sam 
bwa mem -M -t 12 -R '@RG\tID:group1\tSM:Uph18\tPL:illumina\tLB:lib1\tPU:Uph18' -p Umbilicaria_phaea_Uphaea_ref_cold.scaffolds.fa Uph18_FDSW202335555-1a_H33VHDSXY_L4_paired_interleaved.fq > Uph18_aligned_reads.sam 
bwa mem -M -t 12 -R '@RG\tID:group1\tSM:Uph19\tPL:illumina\tLB:lib1\tPU:Uph19' -p Umbilicaria_phaea_Uphaea_ref_cold.scaffolds.fa Uph19_FDSW202335556-1a_H33VHDSXY_L4_paired_interleaved.fq > Uph19_aligned_reads.sam 

bwa mem -M -t 12 -R '@RG\tID:group1\tSM:Uph22\tPL:illumina\tLB:lib1\tPU:Uph22' -p Umbilicaria_phaea_Uphaea_ref_cold.scaffolds.fa Uph22_FDSW202335557-1a_H33VHDSXY_L4_paired_interleaved.fq > Uph22_aligned_reads.sam 
bwa mem -M -t 12 -R '@RG\tID:group1\tSM:Uph23\tPL:illumina\tLB:lib1\tPU:Uph23' -p Umbilicaria_phaea_Uphaea_ref_cold.scaffolds.fa Uph23_FDSW202335558-1a_H33VHDSXY_L4_paired_interleaved.fq > Uph23_aligned_reads.sam 
bwa mem -M -t 12 -R '@RG\tID:group1\tSM:Uph24\tPL:illumina\tLB:lib1\tPU:Uph24' -p Umbilicaria_phaea_Uphaea_ref_cold.scaffolds.fa Uph24_FDSW202335559-1a_H33VHDSXY_L4_paired_interleaved.fq > Uph24_aligned_reads.sam 
bwa mem -M -t 12 -R '@RG\tID:group1\tSM:Uph25\tPL:illumina\tLB:lib1\tPU:Uph25' -p Umbilicaria_phaea_Uphaea_ref_cold.scaffolds.fa Uph25_FDSW202335560-1a_H33VHDSXY_L4_paired_interleaved.fq > Uph25_aligned_reads.sam 
bwa mem -M -t 12 -R '@RG\tID:group1\tSM:Uph26\tPL:illumina\tLB:lib1\tPU:Uph26' -p Umbilicaria_phaea_Uphaea_ref_cold.scaffolds.fa Uph26_FDSW202335561-1a_H33VHDSXY_L4_paired_interleaved.fq > Uph26_aligned_reads.sam 
bwa mem -M -t 12 -R '@RG\tID:group1\tSM:Uph27\tPL:illumina\tLB:lib1\tPU:Uph27' -p Umbilicaria_phaea_Uphaea_ref_cold.scaffolds.fa Uph27_FDSW202335562-1a_H33VHDSXY_L4_paired_interleaved.fq > Uph27_aligned_reads.sam 
bwa mem -M -t 12 -R '@RG\tID:group1\tSM:Uph28\tPL:illumina\tLB:lib1\tPU:Uph28' -p Umbilicaria_phaea_Uphaea_ref_cold.scaffolds.fa Uph28_FDSW202335563-1a_H33VHDSXY_L4_paired_interleaved.fq > Uph28_aligned_reads.sam 

# U. pustulata

bwa mem -M -t 12 -R '@RG\tID:group1\tSM:IT1\tPL:illumina\tLB:lib1\tPU:IT1' -p Lpus_4dec_AN3.masked.fasta IT_Pool1_paired_interleaved.fq > IT_Pool1_aligned_reads.sam
bwa mem -M -t 12 -R '@RG\tID:group1\tSM:IT2\tPL:illumina\tLB:lib1\tPU:IT2' -p Lpus_4dec_AN3.masked.fasta IT_Pool1_paired_interleaved.fq > IT_Pool2_aligned_reads.sam
bwa mem -M -t 12 -R '@RG\tID:group1\tSM:IT3\tPL:illumina\tLB:lib1\tPU:IT3' -p Lpus_4dec_AN3.masked.fasta IT_Pool1_paired_interleaved.fq > IT_Pool3_aligned_reads.sam
bwa mem -M -t 12 -R '@RG\tID:group1\tSM:IT4\tPL:illumina\tLB:lib1\tPU:IT4' -p Lpus_4dec_AN3.masked.fasta IT_Pool1_paired_interleaved.fq > IT_Pool4_aligned_reads.sam
bwa mem -M -t 12 -R '@RG\tID:group1\tSM:IT5\tPL:illumina\tLB:lib1\tPU:IT5' -p Lpus_4dec_AN3.masked.fasta IT_Pool1_paired_interleaved.fq > IT_Pool5_aligned_reads.sam
bwa mem -M -t 12 -R '@RG\tID:group1\tSM:IT6\tPL:illumina\tLB:lib1\tPU:IT6' -p Lpus_4dec_AN3.masked.fasta IT_Pool1_paired_interleaved.fq > IT_Pool6_aligned_reads.sam

bwa mem -M -t 12 -R '@RG\tID:group1\tSM:ESi1\tPL:illumina\tLB:lib1\tPU:ESi1' -p Lpus_4dec_AN3.masked.fasta ESi_Pool1_paired_interleaved.fq > ESi_Pool1_aligned_reads.sam
bwa mem -M -t 12 -R '@RG\tID:group1\tSM:ESi2\tPL:illumina\tLB:lib1\tPU:ESi2' -p Lpus_4dec_AN3.masked.fasta ESi_Pool2_paired_interleaved.fq > ESi_Pool2_aligned_reads.sam
bwa mem -M -t 12 -R '@RG\tID:group1\tSM:ESi3\tPL:illumina\tLB:lib1\tPU:ESi3' -p Lpus_4dec_AN3.masked.fasta ESi_Pool3_paired_interleaved.fq > ESi_Pool3_aligned_reads.sam

bwa mem -M -t 12 -R '@RG\tID:group1\tSM:ESii1\tPL:illumina\tLB:lib1\tPU:ESii1' -p Lpus_4dec_AN3.masked.fasta ESii_Pool1_paired_interleaved.fq > ESii_Pool1_aligned_reads.sam
bwa mem -M -t 12 -R '@RG\tID:group1\tSM:ESii2\tPL:illumina\tLB:lib1\tPU:ESii2' -p Lpus_4dec_AN3.masked.fasta ESii_Pool2_paired_interleaved.fq > ESii_Pool2_aligned_reads.sam
bwa mem -M -t 12 -R '@RG\tID:group1\tSM:ESii3\tPL:illumina\tLB:lib1\tPU:ESii3' -p Lpus_4dec_AN3.masked.fasta ESii_Pool3_paired_interleaved.fq > ESii_Pool3_aligned_reads.sam
bwa mem -M -t 12 -R '@RG\tID:group1\tSM:ESii4\tPL:illumina\tLB:lib1\tPU:ESii4' -p Lpus_4dec_AN3.masked.fasta ESii_Pool4_paired_interleaved.fq > ESii_Pool4_aligned_reads.sam
bwa mem -M -t 12 -R '@RG\tID:group1\tSM:ESii5\tPL:illumina\tLB:lib1\tPU:ESii5' -p Lpus_4dec_AN3.masked.fasta ESii_Pool5_paired_interleaved.fq > ESii_Pool5_aligned_reads.sam
bwa mem -M -t 12 -R '@RG\tID:group1\tSM:ESii6\tPL:illumina\tLB:lib1\tPU:ESii6' -p Lpus_4dec_AN3.masked.fasta ESii_Pool6_paired_interleaved.fq > ESii_Pool6_aligned_reads.sam

```

# 4. Prepare pooled populations for combining

Next, we convert the sam files to bam files, and extract only the mapped reads using a shell file (04_sam_to_bam.sh). 

```{bash 04_sam_to_bam.sh}

#!/bin/sh

echo "name of the population:"

read pop

nohup samtools view -Sb ${pop}_aligned_reads.sam > ${pop}_aligned_reads.bam
samtools view -h -F 4 -b ${pop}_aligned_reads.bam > ${pop}_aligned_reads_only_mapped.bam
java -Xmx8g -jar ~/tools/picard.jar SortSam -I ${pop}_aligned_reads_only_mapped.bam -O ${pop}_aligned_reads.sort.bam -VALIDATION_STRINGENCY SILENT -SO coordinate
java -Xmx8g -jar ~/tools/picard.jar MarkDuplicates -I ${pop}_aligned_reads.sort.bam -O ${pop}_aligned_reads.sort.rmd.bam -M ${pop}_dupstat.txt -VALIDATION_STRINGENCY SILENT -REMOVE_DUPLICATES true
samtools view -q 20 -f 0x0002 -F 0x0004 -F 0x0008 -b ${pop}_aligned_reads.sort.rmd.bam > ${pop}_aligned_reads.sort.rmd.q20.bam


```

# 5. Generating the sync file as a primary popoolation2 input

The next step is to to combine the sorted and mapped reads into an mpileup file, which is used to speed up the process of creating the sync file, which is the main input format for popoolation2.

```{bash Pools are combined using mpileup}

samtools mpileup -B -Q 0 -f Umbilicaria_phaea_Uphaea_ref_cold.scaffolds.fa Uph16_aligned_reads.sort.rmd.q20.bam Uph17_aligned_reads.sort.rmd.q20.bam Uph18_aligned_reads.sort.rmd.q20.bam Uph19_aligned_reads.sort.rmd.q20.bam Uph22_aligned_reads.sort.rmd.q20.bam Uph23_aligned_reads.sort.rmd.q20.bam Uph24_aligned_reads.sort.rmd.q20.bam Uph25_aligned_reads.sort.rmd.q20.bam Uph26_aligned_reads.sort.rmd.q20.bam Uph27_aligned_reads.sort.rmd.q20.bam Uph28_aligned_reads.sort.rmd.q20.bam> Uph_both_combined.mpileup 

# U. pustulata

samtools mpileup -B -Q 0 -f Lpus_4dec_AN3.masked.fasta ESii_Pool1_aligned_reads.sort.rmd.q20.bam ESii_Pool2_aligned_reads.sort.rmd.q20.bam ESii_Pool3_aligned_reads.sort.rmd.q20.bam ESii_Pool4_aligned_reads.sort.rmd.q20.bam ESii_Pool5_aligned_reads.sort.rmd.q20.bam ESii_Pool6_aligned_reads.sort.rmd.q20.bam ESi_Pool1_aligned_reads.sort.rmd.q20.bam ESi_Pool2_aligned_reads.sort.rmd.q20.bam ESi_Pool3_aligned_reads.sort.rmd.q20.bam IT_Pool1_aligned_reads.sort.rmd.q20.bam IT_Pool2_aligned_reads.sort.rmd.q20.bam IT_Pool3_aligned_reads.sort.rmd.q20.bam IT_Pool4_aligned_reads.sort.rmd.q20.bam IT_Pool5_aligned_reads.sort.rmd.q20.bam IT_Pool6_aligned_reads.sort.rmd.q20.bam > IT123456_ESii123456_ESi123.mpileup

```

Next, we identify genomic indel regions and remove them.

```{bash Genomic indel regions are identified}

# U. phaea

perl popoolation2_1201/indel_filtering/identify-indel-regions.pl --indel-window 5 --min-count 2 --input Uph_combined.mpileup --output Uph_combined.indel-regions.gtf &

# U. pustulata

perl popoolation2_1201/indel_filtering/identify-indel-regions.pl --indel-window 5 --min-count 2 --input IT123456_ESii123456_ESi123.mpileup --output IT123456_ESii123456_ESi123.indel-regions.gtf &


```


Then, we utilize these indel regions to filter the sync file. In popoolation1, this was done by filtering the mpileup file, and then creating the sync file; it isn't clear to me what difference it should make to filter the sync file after, rather than filtering the mpileup before generating the sync file, other than presumably a longer sync file generation time. Since mpileup->sync file conversion is apparently a bottle neck, it's not clear why this was changed in the new version of popoolation (ppopoolation2). Nonetheless, we will follow the popoolation2 pipeline and next generate the sync file:


```{bash mpileup file converted to sync file}

# U. phaea

java -ea -Xmx400g -jar popoolation2_1201/mpileup2sync.jar --input Uph_both_combined.mpileup --output Uph_both_combined.sync --fastq-type sanger --min-qual 20 --threads 70

# U. pustulata

java -ea -Xmx400g -jar popoolation2_1201/mpileup2sync.jar --input IT123456_ESii123456_ESi123.mpileup --output IT123456_ESii123456_ESi123.sync --fastq-type sanger --min-qual 20 --threads 70

```

Now, we can use the *indel-regions.gtf file to filter the sync file:

```{bash Indel regions are used to filter the mpileup file}

# U. phaea

perl popoolation2_1201/indel_filtering/filter-sync-by-gtf.pl --gtf Uph_both_combined.indel-regions.gtf --input Uph_both_combined.sync --output filtered.Uph_both_combined.sync &

# U. pustulata

perl popoolation2_1201/indel_filtering/filter-sync-by-gtf.pl --gtf IT123456_ESii123456_ESi123.indel-regions.gtf --input IT123456_ESii123456_ESi123.sync --output filtered.IT123456_ESii123456_ESi123.sync &

```

Next, we subsample the sync file:

```{bash sync files are subsampled}

# U. phaea

perl popoolation2_1201/subsample-synchronized.pl --input filtered.Uph_both_combined.sync --output subsampled_30_filtered.Uph_both_combined.sync --target-coverage 30 --max-coverage 2% --method withoutreplace &

# U. pustulata

perl popoolation2_1201/subsample-synchronized.pl --input filtered.IT123456_ESii123456_ESi123.sync --output subsampled_30_filtered.IT123456_ESii123456_ESi123.sync --target-coverage 30 --max-coverage 2% --method withoutreplace &

```

# 6. Calculating SNP frequency differences and FST sliding windows using the subsampled file

Now we can proceed to creating the outputs of Popoolation2, namely SNP frequency differences and Fst scores.

The first step is to calculate SNP frequencies, which yields two files with the extensions: "_rc" and "_pwc".

    _rc: this file contains the major and minor alleles for every SNP in a concise format
    _pwc: this file contains the differences in allele frequencies for every pairwise comparision of the populations present in the synchronized file
    For details see the man pages of the script

```{bash SNP frequencies}

# U. phaea

nohup perl popoolation2_1201/snp-frequency-diff.pl --input subsampled_30_filtered.Uph_both_combined.sync --output-prefix subsampled_30_filtered.Uph_both_combined --max-coverage 100 --min-count 3 &

# U. pustulata
nohup perl popoolation2_1201/snp-frequency-diff.pl --input subsampled_30_filtered.IT123456_ESii123456_ESi123.sync --output-prefix subsampled_30_filtered.IT123456_ESii123456_ESi123 --max-coverage 100 --min-count 3 &

```

Next, we proceed to calculating Fst values using a sliding window approach. Because we want to get the Fst value for each individual SNP, we use a window size of 1; this automatically suppresses any output for bases that have no SNPs (otherwise, you would presumably have windows in your output with no differences between any of the populations).

```{bash Fst sliding windows}
# U. phaea

nohup perl popoolation2_1201/fst-sliding.pl --input subsampled_30_filtered.Uph_both_combined.sync --output subsampled_30_filtered.Uph_both_combined.fst --suppress-noninformative --min-count 3 --max-coverage 100 --min-covered-fraction 1 --window-size 1 --step-size 1 --pool-size 100 &

# U. pustulata
nohup perl popoolation2_1201/fst-sliding.pl --input subsampled_30_filtered.IT123456_ESii123456_ESi123.sync --output subsampled_30_filtered.IT123456_ESii123456_ESi123.fst --suppress-noninformative --min-count 3 --max-coverage 100 --min-covered-fraction 1 --window-size 1 --step-size 1 --pool-size 100 &

```

Finally, we use a subsetted gtf file (e.g. with only genes, exons, or some other set of features) to create a gene-wise sync file. We then calculate SNP-wise Fst values for only genic SNPs, and whole-gene Fst values.

```{bash make gene-wise sync files and calculate their Fst values}

# U. phaea:

## prepare gene-wise sync file
nohup perl popoolation2_1201/create-genewise-sync.pl --input subsampled_30_filtered.Uph_both_combined.sync --gtf Uphaea_TBG_1112_all_genes.gff3 --output GENES_subsampled_30_filtered.Uph_both_combined.sync

## calculate Fst values for only genic SNPs
perl popoolation2_1201/fst-sliding.pl --suppress-noninformative --min-count 3 --max-coverage 100 --pool-size 100 --min-covered-fraction 1 --window-size 1 --step-size 1 --input GENES_subsampled_30_filtered.Uph_both_combined.sync --output GENES_subsampled_30_filtered.Uph_both_combined.fst

## calculate whole-gene Fst values
perl popoolation2_1201/fst-sliding.pl --min-count 3 --max-coverage 100 --pool-size 100 --min-covered-fraction 1 --window-size 1000000 --step-size 1000000 --input GENES_subsampled_30_filtered.Uph_both_combined.sync --output GENES_subsampled_30_filtered.Uph_both_combined_whole_gene.fst

# U. pustulata:

## prepare gene-wise sync file
nohup perl popoolation2_1201/create-genewise-sync.pl --input subsampled_30_filtered.IT123456_ESii123456_ESi123.sync --gtf all_genes.gtf --output GENES_subsampled_30_filtered.IT123456_ESii123456_ESi123.sync

## calculate Fst values for only genic SNPs
perl popoolation2_1201/fst-sliding.pl --suppress-noninformative --input GENES_subsampled_30_filtered.IT123456_ESii123456_ESi123.sync --output GENES_subsampled_30_filtered.IT123456_ESii123456_ESi123.fst --min-count 3 --max-coverage 100 --pool-size 100 --min-covered-fraction 1 --window-size 1 --step-size 1

## calculate whole-gene Fst values
perl popoolation2_1201/fst-sliding.pl --min-count 3 --max-coverage 100 --pool-size 100 --min-covered-fraction 1 --window-size 1000000 --step-size 1000000 --input GENES_subsampled_30_filtered.IT123456_ESii123456_ESi123.sync --output GENES_subsampled_30_filtered.IT123456_ESii123456_ESi123_whole_gene.fst

```


## CMH tests

Another potential test is the Cochran-Mantel-Haenzel test, which looks for consistent allele frequency changes in biological "replicates" (e.g. similar populations from distinct sites, like those from our multiple gradients).

The comparisons need to be specified in advance; they should take the form of independent contrasts. For example, with two populations of high/low elevation pops (low: 1, 3; high: 2, 4) the contrasts should be 1-2, 3-4

For U. pustulata: IT123456_ESii123456_ESi123
Med: 1, 2, 3, 4, 7, 8, 13,  
Int: 5, 9, 
Cold: 6, 10, 11, 12, 14, 15


For U. phaea: 
  Med.: 1, 2
  Int.: 3
  Cold: 4
  
  Med.: 5, 6, 7 
  Int.: 8
  Cold: 9, 10, 11
  
We compare the distribution of the -log10p values of the SNPs in the circadian/temp genes with 
  i) a random choice of the same number of other genes <- 
  ii) a random choice of non-genic stretches of the same length as the circadian/temp genes 
  iii) a random set of the same number of genome-wide SNPs as in your target genes. <- done
  
If the circadian/temp genes are higher differentiated than expected by distance/population structure/-history/chance they should show a significantly higher mean -log10p value than all the three others.

```{bash CMH test}

# run CMH test on all SNPs, specifying the extremes
perl popoolation2_1201/cmh-test.pl --input subsampled_30_filtered.Uph_both_combined.sync --output subsampled_30_filtered.Uph_both_combined.cmh --min-count 4 --min-coverage 4 --max-coverage 400 --population 1-4,5-11

perl popoolation2_1201/cmh-test.pl --input subsampled_30_filtered.IT123456_ESii123456_ESi123.sync --output subsampled_30_filtered.IT123456_ESii123456_ESi123.cmh --min-count 4 --min-coverage 4 --max-coverage 400 --population 1-6,7-12,13-15

#######
###  get a random distribution of 50 genes 
shuf -n 50 GENES_subsampled_30_filtered.IT123456_ESii123456_ESi123_whole_gene.fst > Upust_random_50_genes.fst

shuf -n 50 GENES_subsampled_30_filtered.Uph_both_combined.fst > Uph_random_50_genes.fst

# extract first column with gene names, then use this as a key to extract genes in gtf
awk '{print $1}' Upust_random_50_genes.fst > Upust_random_50_genes.txt
grep -f Upust_random_50_genes.txt all_genes.gtf > Upust_random_50_genes.gtf
awk -v OFS='\t' '{print $1, $4, $5}' Upust_random_50_genes.gtf > Upust_random_50_genes_range.txt

awk '{print $1}' Uph_random_50_genes.fst > Uph_random_50_genes.txt
grep -f Uph_random_50_genes.txt Uphaea_TBG_1112_all_genes.gtf > Uph_random_50_genes.gtf
awk -v OFS='\t' '{print $1, $4, $5}' Uph_random_50_genes.gtf > Uph_random_50_genes_range.txt

# convert to a bed file
mv Upust_random_50_genes_range.txt > Upust_random_50_genes_range.bed
mv Uph_random_50_genes_range.txt > Uph_random_50_genes_range.bed

# filter bam using the bed file for gene ranges
for pop in IT_Pool1 IT_Pool2 IT_Pool3 IT_Pool4 IT_Pool5 IT_Pool6 ESii_Pool1 ESii_Pool2 ESii_Pool3 ESii_Pool4 ESii_Pool5 ESii_Pool6 ESi_Pool1 ESi_Pool2 ESi_Pool3; do bedtools intersect -abam ${pop}_aligned_reads.sort.rmd.q20.bam -b Upust_random_50_genes_range.bed > ${pop}_random_50_genes_range.bam; done

# bam to mpileup
samtools mpileup -B -Q 0 -f Lpus_4dec_AN3.masked.fasta ESii_Pool1_random_50_genes_range.bam ESii_Pool2_random_50_genes_range.bam ESii_Pool3_random_50_genes_range.bam ESii_Pool4_random_50_genes_range.bam ESii_Pool5_random_50_genes_range.bam ESii_Pool6_random_50_genes_range.bam ESi_Pool1_random_50_genes_range.bam ESi_Pool2_random_50_genes_range.bam ESi_Pool3_random_50_genes_range.bam IT_Pool1_random_50_genes_range.bam IT_Pool2_random_50_genes_range.bam IT_Pool3_random_50_genes_range.bam IT_Pool4_random_50_genes_range.bam IT_Pool5_random_50_genes_range.bam IT_Pool6_random_50_genes_range.bam > Upust_random_50_genes_range.mpileup

# mpileup to sync
java -ea -Xmx400g -jar popoolation2_1201/mpileup2sync.jar --input Upust_random_50_genes_range.mpileup --output Upust_random_50_genes_range.sync --fastq-type sanger --min-qual 20 --threads 70

# run CMH
#perl popoolation2_1201/cmh-test.pl --input Upust_random_50_genes_range.sync --output Upust_random_50_genes_range.cmh --min-count 6 --min-coverage 10 --max-coverage 100 --population 1-6,7-12,13-15

#perl popoolation2_1201/cmh-test.pl --input Upust_random_50_genes_range.sync --output Upust_random_50_genes_range_alt_order.cmh --min-count 6 --min-coverage 10 --max-coverage 100 --population 1-6,7-9,10-15
```


# Checking CMH test results 

Next we want to compare the CMH test results to the subset files we have for the circadian- and temperature-associated SNPs.


```{r CMH test analysis: Uph}

# load data
cmh <- read.table('~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/subsampled_30_filtered.Uph_both_combined.cmh',header = F, sep = "\t")
random_50 <- read.table('~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/Uph_random_50_genes_range.sync',header = F, sep = "\t")
circ <- read.table('~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/subsampled_30_filtered.Uph16_circadian.sync',header = F, sep = "\t")
temp <- read.table('~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/subsampled_30_filtered.Uph16_temperature.sync',header = F, sep = "\t")

# subset down
cmh1 <- cmh[,c(1,2,15)]
random_50.1 <- random_50[,c(1,2)]
circ1 <- circ[,c(1,2)]
temp1 <- temp[,c(1,2)]

cmh1$label <- paste(cmh1$V1,cmh1$V2)
random_50.1$label <- paste(random_50.1$V1,random_50.1$V2)
circ1$label <- paste(circ1$V1,circ1$V2)
temp1$label <- paste(temp1$V1,temp1$V2)

cmh1$fdr <- p.adjust(cmh1$V15, method="BH")

cmh_random <- merge(cmh1,random_50.1, by = c("label","V1","V2"))
cmh_circ <- merge(cmh1,circ1, by = c("label","V1","V2"))
cmh_temp <- merge(cmh1,temp1, by = c("label","V1","V2"))

cmh_random1 <- na.omit(cmh_random)
cmh_circ1 <- na.omit(cmh_circ)
cmh_temp1 <- na.omit(cmh_temp)

cmh_random1$random <- -log10(cmh_random1$V15)
cmh_circ1$circ <- -log10(cmh_circ1$V15)
cmh_temp1$temp <- -log10(cmh_temp1$V15)

# FDR calc
cmh_random1$random_fdr <- p.adjust(cmh_random1$V15, method="BH")
cmh_circ1$circ_fdr <- p.adjust(cmh_circ1$V15, method="BH")
cmh_temp1$temp_fdr <- p.adjust(cmh_temp1$V15, method="BH")

cmh_random1$random_fdr2 <- -log10(cmh_random1$random_fdr)
cmh_circ1$circ_fdr2 <- -log10(cmh_circ1$circ_fdr)
cmh_temp1$temp_fdr2 <- -log10(cmh_temp1$temp_fdr)

cmh_random1$random_fdr3 <- -log10(cmh_random1$fdr)
cmh_circ1$circ_fdr3 <- -log10(cmh_circ1$fdr)
cmh_temp1$temp_fdr3 <- -log10(cmh_temp1$fdr)

random <- cmh_random1[,c(2,3,6)]
circadian <- cmh_circ1[,c(2,3,6)]
temperature <- cmh_temp1[,c(2,3,6)]

library(dplyr)

random_anti_circ <- anti_join(random,circadian, by = c("V1","V2"))
random_anti_temp <- anti_join(random_anti_circ,temperature, by = c("V1","V2"))

mean(random_anti_temp$random) # 7.439486
mean(circadian$circ) # 7.382804
mean(temperature$temp) # 6.870254
wilcox.test(random_anti_temp$random, circadian$circ, alternative = "less") # W = 80624205, p-value = 0.8016
wilcox.test(random_anti_temp$random, temperature$temp, alternative = "less") # W = 93481336, p-value = 1

mean(cmh_random1$random_fdr2) # 7.011796
mean(cmh_circ1$circ_fdr2) # 6.954924
mean(cmh_temp1$temp_fdr2) # 6.442636
wilcox.test(cmh_random1$random_fdr2, cmh_circ1$circ_fdr2, alternative = "less") # W = 81786952, p-value = 0.9978
wilcox.test(cmh_random1$random_fdr2, cmh_temp1$temp_fdr2, alternative = "less") # W = 95226832, p-value = 1

mean(cmh_random1$random_fdr3) # 7.011796
mean(cmh_circ1$circ_fdr3) # 6.954924
mean(cmh_temp1$temp_fdr3) # 6.442636
wilcox.test(cmh_random1$random_fdr3, cmh_circ1$circ_fdr3, alternative = "less") # W = 80618634, p-value = 0.7989
wilcox.test(cmh_random1$random_fdr3, cmh_temp1$temp_fdr3, alternative = "less") # W = 93477434, p-value = 1

library(ggplot2)

ggplot() +
  #geom_density(data = nongenic_anti_temp, aes(x = nongenic, fill = "nongenic_50_3.5kb"), adjust=1.5, alpha=.4) +
  geom_density(data = random_anti_temp, aes(x = random, fill = "random genes (50)"), adjust=1.5, alpha=.4) +
  geom_density(data = circadian, aes(x = circ, fill = "circadian genes (50)"), adjust=1.5, alpha=.4) +
  geom_density(data = temperature, aes(x = temp, fill = "temperature genes (37)"), adjust=1.5, alpha=.4) +
  scale_x_continuous(name = "-log10(p-val)") +
  theme_classic(base_size = 12)

pdf(file.path("Uph_CMH_boxplots.pdf"), width = 6, height = 3)
ggplot() +
  #geom_boxplot(data = nongenic_anti_temp, aes(x = 1, y = nongenic, fill = "nongenic_50_3.5kb")) +
  geom_boxplot(data = random_anti_temp, aes(x = 3, y = random, fill = "Randomly selected genes (50)")) +
  geom_boxplot(data = circadian, aes(x = 4, y = circ, fill = "Circadian-associated homologs (50)")) +
  geom_boxplot(data = temperature, aes(x = 5, y = temp, fill = "Temperature-associated homologs (37)")) +
  scale_fill_manual(values=c("salmon", "gray", "cadetblue")) +
  scale_x_discrete(name = "Gene set") +
  scale_y_continuous(name = "-log(p-value)") +
  theme_classic(base_size = 12)
dev.off()

pdf(file.path("Uph_CMH_boxplots_FDR.pdf"), width = 6, height = 3)
ggplot() +
  #geom_boxplot(data = nongenic_anti_temp, aes(x = 1, y = nongenic, fill = "nongenic_50_3.5kb")) +
  geom_boxplot(data = cmh_random1, aes(x = 3, y = random_fdr3, fill = "Randomly selected genes (50)")) +
  geom_boxplot(data = cmh_circ1, aes(x = 4, y = circ_fdr3, fill = "Circadian-associated homologs (50)")) +
  geom_boxplot(data = cmh_temp1, aes(x = 5, y = temp_fdr3, fill = "Temperature-associated homologs (37)")) +
  scale_fill_manual(values=c("salmon", "gray", "cadetblue")) +
  scale_x_discrete(name = "Gene set") +
  scale_y_continuous(name = "-log(p-value)") +
  theme_classic(base_size = 12)
dev.off()

ggplot() +
  #geom_point(data = nongenic_anti_temp, aes(x = 1, y = nongenic, fill = "nongenic_50_3.5kb"), position = "jitter") +
  geom_point(data = random_anti_temp, aes(x = 2, y = random, col = "random genes (50)"), position = "jitter") +
  geom_point(data = circadian, aes(x = 3, y = circ, col = "circadian genes (50)"), position = "jitter") +
  geom_point(data = temperature, aes(x = 4, y = temp, col = "temperature genes (37)"), position = "jitter") +
  scale_x_discrete(name = "Gene set") +
  theme_classic(base_size = 12)

```

```{r CMH analysis: Upust}

# load data
cmh <- read.table('~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/subsampled_30_filtered.IT123456_ESii123456_ESi123.cmh',header = F, sep = "\t")
random_50_genes <- read.table('~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/Upust_random_50_genes_range.sync',header = F, sep = "\t")
circ <- read.table('~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/subsampled_30_filtered.Upust_IT6_circadian.sync',header = F, sep = "\t")
temp <- read.table('~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/subsampled_30_filtered.Upust_IT6_temperature.sync',header = F, sep = "\t")

# subset down
cmh1 <- cmh[,c(1,2,19)]
random_50.1 <- random_50_genes[,c(1,2)]
circ1 <- circ[,c(1,2)]
temp1 <- temp[,c(1,2)]

cmh1$fdr <- p.adjust(cmh1$V19, method="BH")

cmh1$label <- paste(cmh1$V1,cmh1$V2)
random_50.1$label <- paste(random_50.1$V1,random_50.1$V2)
circ1$label <- paste(circ1$V1,circ1$V2)
temp1$label <- paste(temp1$V1,temp1$V2)

cmh_50 <- merge(cmh1,random_50.1, by = c("label","V1","V2"))
cmh_circ <- merge(cmh1,circ1, by = c("label","V1","V2"))
cmh_temp <- merge(cmh1,temp1, by = c("label","V1","V2"))

random_50.2 <- na.omit(cmh_50)
cmh_circ1 <- na.omit(cmh_circ)
cmh_temp1 <- na.omit(cmh_temp)

random_50.2$random2 <- -log10(random_50.2$V19)
cmh_circ1$circ <- -log10(cmh_circ1$V19)
cmh_temp1$temp <- -log10(cmh_temp1$V19)

# FDR calculation
random_50.2$random2_fdr <- p.adjust(random_50.2$V19, method="BH")
cmh_circ1$circ_fdr <- p.adjust(cmh_circ1$V19, method="BH")
cmh_temp1$temp_fdr <- p.adjust(cmh_temp1$V19, method="BH")

random_50.2$random2_fdr2 <- -log10(random_50.2$random2_fdr)
cmh_circ1$circ_fdr2 <- -log10(cmh_circ1$circ_fdr)
cmh_temp1$temp_fdr2 <- -log10(cmh_temp1$temp_fdr)

random_50.2$random2_fdr3 <- -log10(random_50.2$fdr)
cmh_circ1$circ_fdr3 <- -log10(cmh_circ1$fdr)
cmh_temp1$temp_fdr3 <- -log10(cmh_temp1$fdr)

random2 <- random_50.2[,c(2,3,6)]
circadian <- cmh_circ1[,c(2,3,6)]
temperature <- cmh_temp1[,c(2,3,6)]

library(dplyr)

mean(random_50.2$random2) # 9.41509
mean(cmh_circ1$circ) # 12.6261
mean(cmh_temp1$temp) # 10.06436
wilcox.test(random_50.2$random2, cmh_circ1$circ, alternative = "less") # W = 612880, p-value = 6.267e-07
wilcox.test(random_50.2$random2, cmh_temp1$temp, alternative = "less") # W = 596751, p-value = 0.7492

mean(random_50.2$random2_fdr2) # 8.986929
mean(cmh_circ1$circ_fdr2) # 12.2009
mean(cmh_temp1$temp_fdr2) # 9.640087
wilcox.test(random_50.2$random2_fdr2, cmh_circ1$circ_fdr2, alternative = "less") # W = 606778, p-value = 9.183e-08
wilcox.test(random_50.2$random2_fdr2, cmh_temp1$temp_fdr2, alternative = "less") # W = 599730, p-value = 0.809

mean(random_50.2$random2_fdr3) # 9.041975
mean(cmh_circ1$circ_fdr3) # 12.15129
mean(cmh_temp1$temp_fdr3) # 9.671536
wilcox.test(random_50.2$random2_fdr3, cmh_circ1$circ_fdr3, alternative = "less") # W = 612933, p-value = 6.368e-07
wilcox.test(random_50.2$random2_fdr3, cmh_temp1$temp_fdr3, alternative = "less") # W = 596889, p-value = 0.7521

library(ggplot2)

ggplot() +
  geom_density(data = random_50.2, aes(x = random2, fill = "random_50"), adjust=1.5, alpha=.4) +
  geom_density(data = circadian, aes(x = circ, fill = "circ"), adjust=1.5, alpha=.4) +
  geom_density(data = temperature, aes(x = temp, fill = "temp"), adjust=1.5, alpha=.4) +
  theme_classic(base_size = 12)

pdf(file.path("Upust_CMH_boxplots.pdf"), width = 6, height = 3)
ggplot() +
  geom_boxplot(data = random_50.2, aes(x = 3, y = random2, fill = "Randomly selected genes (50)")) +
  geom_boxplot(data = circadian, aes(x = 4, y = circ, fill = "Circadian-associated homologs (51)")) +
  geom_boxplot(data = temperature, aes(x = 5, y = temp, fill = "Temperature-associated homologs (37)")) +
  scale_fill_manual(values=c("salmon", "gray", "cadetblue")) +
  scale_x_discrete(name = "Gene set") +
  scale_y_continuous(name = "-log(p-value)") +
  theme_classic(base_size = 12)
dev.off()

pdf(file.path("Upust_CMH_boxplots_FDR.pdf"), width = 6, height = 3)
ggplot() +
  geom_boxplot(data = random_50.2, aes(x = 3, y = random2_fdr3, fill = "Randomly selected genes (50)")) +
  geom_boxplot(data = cmh_circ1, aes(x = 4, y = circ_fdr3, fill = "Circadian-associated homologs (51)")) +
  geom_boxplot(data = cmh_temp1, aes(x = 5, y = temp_fdr3, fill = "Temperature-associated homologs (37)")) +
  scale_fill_manual(values=c("salmon", "gray", "cadetblue")) +
  scale_x_discrete(name = "Gene set") +
  scale_y_continuous(name = "-log(p-value)") +
  theme_classic(base_size = 12)
dev.off()

ggplot() +
  geom_point(data = random2_anti_temp, aes(x = 3, y = random2, col = "random_50"), position = "jitter") +
  geom_point(data = circadian, aes(x = 4, y = circ, col = "circ"), position = "jitter") +
  geom_point(data = temperature, aes(x = 5, y = temp, col = "temp"), position = "jitter") +
  theme_classic(base_size = 12)

```

# Check variation in diversity between U. phaea and U. pustulata

```{r}

library(dplyr)
library(plyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(R.devices)

# load data
Fst_Upust <- read.table('GENES_subsampled_30_filtered.IT123456_ESii123456_ESi123.fst',header = F, sep = "\t")
Fst_Uph <- read.table('GENES_subsampled_30_filtered.Uph_both_combined.fst',header = F, sep = "\t")

# Uph

# subset down
Fst_Uph_grad1 <- Fst_Uph[,c(1,6,7,8,16,17,25)]
Fst_Uph_grad2 <- Fst_Uph[,c(1,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60)]

## SN

#make column names
colnames(Fst_Uph_grad1) <- c("gene","1-2","1-3","1-4","2-3","2-4","3-4")
Fst_Uph_grad1_cleaned <- as.data.frame(lapply(Fst_Uph_grad1, function(x) gsub("^.*=","",x)))

# melt the data frame
melted <- melt(Fst_Uph_grad1_cleaned, id = c("gene"), variable.name = "Pairwise", value.name = "Fst")

melted$Fst1 <- as.numeric(melted$Fst)

plot1  <- ggplot(melted, aes(x = Pairwise, y = Fst1)) +
    geom_boxplot() +
    scale_x_discrete("Pairwise") +
    scale_y_continuous("Fst") +
    theme_classic(base_size = 12)

suppressGraphics(ggsave('Uph_SN_genic_SNPS_Fst_all_pairwise_boxplots.png', plot1))

melted %>%
  group_by(Pairwise) %>%
  summarise_at(vars(Fst1), list(name = mean))

#   Pairwise   Fst    Elev_diff   
#   <fct>     <dbl>
# 1 X1.2     0.0208
# 2 X1.3     0.115 
# 3 X1.4     0.130 
# 4 X2.3     0.108 
# 5 X2.4     0.124 
# 6 X3.4     0.0564

## MJ

#make column names
colnames(Fst_Uph_grad2) <- c("gene","1-2","1-3","1-4","1-5","1-6","1-7","2-3","2-4","2-5","2-6","2-7","3-4","3-5","3-6","3-7","4-5","4-6","4-7","5-6","5-7","6-7")
Fst_Uph_grad2_cleaned <- as.data.frame(lapply(Fst_Uph_grad2, function(x) gsub("^.*=","",x)))

# melt the data frame
melted2 <- melt(Fst_Uph_grad2_cleaned, id = c("gene"), variable.name = "Pairwise", value.name = "Fst")

melted2$Fst1 <- as.numeric(melted2$Fst)

plot2  <- ggplot(melted2, aes(x = Pairwise, y = Fst1)) +
    geom_boxplot() +
    scale_x_discrete("Pairwise") +
    scale_y_continuous("Fst") +
    theme_classic(base_size = 12)

suppressGraphics(ggsave('Uph_MJ_genic_SNPS_Fst_all_pairwise_boxplots.png', plot2))

melted2 %>%
  group_by(Pairwise) %>%
  summarise_at(vars(Fst1), list(name = mean)) %>% print(n=21)

#    Pairwise   name
#    <fct>     <dbl>
#  1 X1.2     0.0177
#  2 X1.3     0.0325
#  3 X1.4     0.0896
#  4 X1.5     0.152 
#  5 X1.6     0.132 
#  6 X1.7     0.119 
#  7 X2.3     0.0330
#  8 X2.4     0.0819
#  9 X2.5     0.144 
# 10 X2.6     0.121 
# 11 X2.7     0.107 
# 12 X3.4     0.0857
# 13 X3.5     0.143 
# 14 X3.6     0.121 
# 15 X3.7     0.107 
# 16 X4.5     0.0947
# 17 X4.6     0.0958
# 18 X4.7     0.0865
# 19 X5.6     0.117 
# 20 X5.7     0.115 
# 21 X6.7     0.0771

# Upust

# subset down
Fst_Upust_gradIT <- Fst_Upust[,c(1,6,7,8,9,10,20,21,22,23,33,34,35,45,46,56)]
Fst_Upust_gradESii <- Fst_Upust[,c(1,75,76,77,78,79,83,84,85,86,90,91,92,96,97,101)]
Fst_Upust_gradESi <- Fst_Upust[,c(1,108,109,110)]

## IT

#make column names
colnames(Fst_Upust_gradIT) <- c("gene","1-2","1-3","1-4","1-5","1-6","2-3","2-4","2-5","2-6 ","3-4","3-5","3-6","4-5","4-6","5-6")
Fst_Upust_gradIT_cleaned <- as.data.frame(lapply(Fst_Upust_gradIT, function(x) gsub("^.*=","",x)))

# melt the data frame
melted3 <- melt(Fst_Upust_gradIT_cleaned, id = c("gene"), variable.name = "Pairwise", value.name = "Fst")

melted3$Fst1 <- as.numeric(melted3$Fst)

plot3  <- ggplot(melted3, aes(x = Pairwise, y = Fst1)) +
    geom_boxplot() +
    scale_x_discrete("Pairwise") +
    scale_y_continuous("Fst") +
    theme_classic(base_size = 12)

suppressGraphics(ggsave('Upust_IT_genic_SNPS_Fst_all_pairwise_boxplots.png', plot3))

melted3 %>%
  group_by(Pairwise) %>%
  summarise_at(vars(Fst1), list(name = mean))

#    Pairwise   name
#    <fct>     <dbl>
#  1 X1.2     0.0289
#  2 X1.3     0.0250
#  3 X1.4     0.0214
#  4 X1.5     0.0703
#  5 X1.6     0.116 
#  6 X2.3     0.0384
#  7 X2.4     0.0346
#  8 X2.5     0.0783
#  9 X2.6     0.124 
# 10 X3.4     0.0246
# 11 X3.5     0.0709
# 12 X3.6     0.112 
# 13 X4.5     0.0694
# 14 X4.6     0.116 
# 15 X5.6     0.0385

## ESii

#make column names
colnames(Fst_Upust_gradESii) <- c("gene","1-2","1-3","1-4","1-5","1-6","2-3","2-4","2-5","2-6 ","3-4","3-5","3-6","4-5","4-6","5-6")
Fst_Upust_gradESii_cleaned <- as.data.frame(lapply(Fst_Upust_gradESii, function(x) gsub("^.*=","",x)))

# melt the data frame
melted4 <- melt(Fst_Upust_gradESii_cleaned, id = c("gene"), variable.name = "Pairwise", value.name = "Fst")

melted4$Fst1 <- as.numeric(melted4$Fst)

plot4  <- ggplot(melted4, aes(x = Pairwise, y = Fst1)) +
    geom_boxplot() +
    scale_x_discrete("Pairwise") +
    scale_y_continuous("Fst") +
    theme_classic(base_size = 12)

suppressGraphics(ggsave('Upust_ESii_genic_SNPS_Fst_all_pairwise_boxplots.png', plot4))

melted4 %>%
  group_by(Pairwise) %>%
  summarise_at(vars(Fst1), list(name = mean))

#    Pairwise    name
#    <fct>      <dbl>
#  1 X1.2     0.0264 
#  2 X1.3     0.104  
#  3 X1.4     0.135  
#  4 X1.5     0.133  
#  5 X1.6     0.137  
#  6 X2.3     0.0630 
#  7 X2.4     0.0909 
#  8 X2.5     0.0890 
#  9 X2.6.    0.0927 
# 10 X3.4     0.0187 
# 11 X3.5     0.0177 
# 12 X3.6     0.0207 
# 13 X4.5     0.00944
# 14 X4.6     0.0117 
# 15 X5.6     0.0108 

## ESi

#make column names
colnames(Fst_Upust_gradESi) <- c("gene","1-2","1-3","2-3")
Fst_Upust_gradESi_cleaned <- as.data.frame(lapply(Fst_Upust_gradESi, function(x) gsub("^.*=","",x)))

# melt the data frame
melted5 <- melt(Fst_Upust_gradESi_cleaned, id = c("gene"), variable.name = "Pairwise", value.name = "Fst")

melted5$Fst1 <- as.numeric(melted5$Fst)

plot5  <- ggplot(melted5, aes(x = Pairwise, y = Fst1)) +
    geom_boxplot() +
    scale_x_discrete("Pairwise") +
    scale_y_continuous("Fst") +
    theme_classic(base_size = 12)

suppressGraphics(ggsave('Upust_ESi_genic_SNPS_Fst_all_pairwise_boxplots.png', plot5))

melted5 %>%
  group_by(Pairwise) %>%
  summarise_at(vars(Fst1), list(name = mean))

#   Pairwise   name
#   <fct>     <dbl>
# 1 X1.2     0.141 
# 2 X1.3     0.142 
# 3 X2.3     0.0113


```


Isolation by distance analysis

```{r}

library(dplyr)
library(plyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(R.devices)
library(ggpubr)
library(tidyverse)

# load data
isolation <- read.csv('isolation_by_distance.csv', header = TRUE)

pdf(file.path("Upust_Uph_genic_Fst_isolation_by_distance.pdf"), width = 5, height = 8)
(plot  <- ggplot(isolation, aes(x = Elev_diff, y = Fst_mean, col = Gradient)) +
    geom_point() +
    geom_smooth(method='lm', se=FALSE, color="black", formula = y ~ x) +
    stat_cor(
     # aes(label = after_stat(rr.label)), 
      method = "pearson",
      cor.coef.name = c("R"),
      color = "black", 
      #geom = "label"
      ) +
    stat_regline_equation(color = "black", label.y = 0.13) +
    scale_x_continuous("Pairwise distances (elevation, meters)") +
    scale_color_manual(values=c("skyblue", "skyblue2", "skyblue4","gray20","gray40"), label = c("ESi","ESii","IT","MJ","SN")) +
    scale_y_continuous("FST") +
    #facet_grid(Species~., ) +
    facet_wrap(~fct_rev(Species), dir = "v") +
    theme_classic(base_size = 12))
dev.off()

model <- lm(data = isolation, formula = Fst_mean ~ Elev_diff:Species)
summary(model)

suppressGraphics(ggsave('Upust_Uph_genic_Fst_isolation_by_distance.png', plot, width = 4, height = 8))

```

Next we need to compare the genome-wide diversity metrics for U pustulata and U phaea to determine whether there are differences in how differentiated the two species are along the gradients.

To do this we will look at three metrics:
1. SNPs/base for each species (to see the density of SNPs in the genome)
2. Generate genewise nucleotide diversity (pi) for each population along the gradients and plot them
3. Generate genewise Tajima's D for each population along the gradients and plot them

# SNPs/base

for U pustulata:
SNPs: 798029
bases: 35.7 Mb
ratio: 0.022

For U phaea:
SNPs: 2770905
bases: 35.1 Mb
ratio: 0.079

# Generate pi and Tajima's D for all genes

```{bash}

## U pustulata

# make individual mpileup files
for pop in IT_Pool1 IT_Pool2 IT_Pool3 IT_Pool4 IT_Pool5 IT_Pool6 ESii_Pool1 ESii_Pool2 ESii_Pool3 ESii_Pool4 ESii_Pool5 ESii_Pool6 ESi_Pool1 ESi_Pool2 ESi_Pool3; do samtools mpileup -B -Q 0 -f Lpus_4dec_AN3.masked.fasta ${pop}_aligned_reads.sort.rmd.q20.bam > Upust.${pop}.mpileup; done

#filter mpileup files
for pop in IT_Pool1 IT_Pool2 IT_Pool3 IT_Pool4 IT_Pool5 IT_Pool6 ESii_Pool1 ESii_Pool2 ESii_Pool3 ESii_Pool4 ESii_Pool5 ESii_Pool6 ESi_Pool1 ESi_Pool2 ESi_Pool3; do perl popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl --gtf IT123456_ESii123456_ESi123.indel-regions.gtf --input Upust.${pop}.mpileup --output filtered.Upust.${pop}.mpileup; done

# then we need to subsample the mpileup files ideally

### for U. pustulata circadian loci:
# IT123456_ESii123456_ESi123: 112,102,94,125,116,123,143,114,137,147,136,144,128,166,168

# subsample_mpileup_genes.sh
#!/bin/sh
pop="IT_Pool1 IT_Pool2 IT_Pool3 IT_Pool4 IT_Pool5 IT_Pool6 ESii_Pool1 ESii_Pool2 ESii_Pool3 ESii_Pool4 ESii_Pool5 ESii_Pool6 ESi_Pool1 ESi_Pool2 ESi_Pool3"
cov="112 102 94 125 116 123 143 114 137 147 136 144 128 166 168"
set -- $cov
for i in $pop
do
    echo "pop:" $i "2% max cov:" $1
    perl popoolation_1.2.2/basic-pipeline/subsample-pileup.pl --input filtered.Upust_${i}.mpileup --output subsampled_30_filtered.Upust.${i}.mpileup --target-coverage 30 --max-coverage $cov --method withoutreplace
    shift
done

# get pi for genes
for pop in IT_Pool1 IT_Pool2 IT_Pool3 IT_Pool4 IT_Pool5 IT_Pool6 ESii_Pool1 ESii_Pool2 ESii_Pool3 ESii_Pool4 ESii_Pool5 ESii_Pool6 ESi_Pool1 ESi_Pool2 ESi_Pool3; do perl popoolation_1.2.2/Variance-at-position.pl --pool-size 100 --min-qual 20 --min-coverage 4 --min-count 2 --max-coverage 400  --measure pi --fastq-type sanger  --min-covered-fraction 0.5 --pileup subsampled_30_filtered.Upust.${pop}.mpileup --gtf all_genes.gtf --output subsampled_30_filtered.Upust.${pop}.genes.pi; done

# get tajima's D for genes
for pop in IT_Pool1 IT_Pool2 IT_Pool3 IT_Pool4 IT_Pool5 IT_Pool6 ESii_Pool1 ESii_Pool2 ESii_Pool3 ESii_Pool4 ESii_Pool5 ESii_Pool6 ESi_Pool1 ESi_Pool2 ESi_Pool3; do perl popoolation_1.2.2/Variance-at-position.pl --pool-size 100 --min-qual 20 --min-coverage 4 --min-count 2 --max-coverage 400  --measure D --fastq-type sanger  --min-covered-fraction 0.5 --pileup subsampled_30_filtered.Upust.${pop}.mpileup --gtf all_genes.gtf --output subsampled_30_filtered.Upust.${pop}.genes.tajD; done

## U phaea

# make individual mpileup files
for i in 16 17 18 19; do samtools mpileup -B -Q 0 -f Umbilicaria_phaea_Uphaea_ref_cold.scaffolds.fa Uph${i}_aligned_reads.sort.rmd.q20.bam > Uph${i}.mpileup; done

for i in 22 23 24 25 26 27 28; do samtools mpileup -B -Q 0 -f Umbilicaria_phaea_Uphaea_ref_cold.scaffolds.fa Uph${i}_aligned_reads.sort.rmd.q20.bam > Uph${i}.mpileup; done

#filter mpileup files
for i in 16 17 18 19 22 23 24 25 26 27 28; do perl popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl --gtf Uph_both_combined.indel-regions.gtf --input Uph${i}.mpileup --output filtered.Uph${i}.mpileup; done


# then we need to subsample the mpileup files ideally

### for U. phaea circadian loci:
# Result: '--max-coverage 2%' is equivalent to '--max-coverage 196,147,161,124,202,164,190,199,157,222,186'

# subsample_mpileup_genes.sh
#!/bin/sh
pop="16 17 18 19 22 23 24 25 26 27 28"
cov="196 147 161 124 202 164 190 199 157 222 186"
set -- $cov
for i in $pop
do
    echo "pop:" $i "2% max cov:" $1
    perl popoolation_1.2.2/basic-pipeline/subsample-pileup.pl --input filtered.Uph${i}.mpileup --output subsampled_30_filtered.Uph${i}.mpileup --target-coverage 30 --max-coverage $cov --method withoutreplace
    shift
done

# get pi for genes
for i in 16 17 18 19 22 23 24 25 26 27 28; do perl popoolation_1.2.2/Variance-at-position.pl --pool-size 100 --min-qual 20 --min-coverage 4 --min-count 2 --max-coverage 400  --measure pi --fastq-type sanger  --min-covered-fraction 0.5  --pileup subsampled_30_filtered.Uph${i}.mpileup --gtf Uphaea_TBG_1112_all_genes.gtf --output subsampled_30_filtered.Uph${i}.genes.pi; done

# get tajima's D for genes
for i in 16 17 18 19 22 23 24 25 26 27 28; do perl popoolation_1.2.2/Variance-at-position.pl --pool-size 100 --min-qual 20 --min-coverage 4 --min-count 2 --max-coverage 400  --measure D --fastq-type sanger  --min-covered-fraction 0.5  --pileup subsampled_30_filtered.Uph${i}.mpileup --gtf Uphaea_TBG_1112_all_genes.gtf --output subsampled_30_filtered.Uph${i}.genes.tajD; done



```

Combine the nucleotide diversity figures and the tajima's D for each gradient

```{r Upust, genes, nucleotide diversity}

library(tidyr)
library(ape)
library(vegan)
library(dplyr)
library(stringr)
library(reshape2)
library(ggplot2)
library(R.devices)

# import each pi file 
pi_1 <- read.delim("subsampled_30_filtered.Upust.IT_Pool1.genes.pi", header = F)
pi_2 <- read.delim("subsampled_30_filtered.Upust.IT_Pool2.genes.pi", header = F)
pi_3 <- read.delim("subsampled_30_filtered.Upust.IT_Pool3.genes.pi", header = F)
pi_4 <- read.delim("subsampled_30_filtered.Upust.IT_Pool4.genes.pi", header = F)
pi_5 <- read.delim("subsampled_30_filtered.Upust.IT_Pool5.genes.pi", header = F)
pi_6 <- read.delim("subsampled_30_filtered.Upust.IT_Pool6.genes.pi", header = F)

pi_7 <- read.delim("subsampled_30_filtered.Upust.ESi_Pool1.genes.pi", header = F)
pi_8 <- read.delim("subsampled_30_filtered.Upust.ESi_Pool2.genes.pi", header = F)
pi_9 <- read.delim("subsampled_30_filtered.Upust.ESi_Pool3.genes.pi", header = F)

pi_10 <- read.delim("subsampled_30_filtered.Upust.ESii_Pool1.genes.pi", header = F)
pi_11 <- read.delim("subsampled_30_filtered.Upust.ESii_Pool2.genes.pi", header = F)
pi_12 <- read.delim("subsampled_30_filtered.Upust.ESii_Pool3.genes.pi", header = F)
pi_13 <- read.delim("subsampled_30_filtered.Upust.ESii_Pool4.genes.pi", header = F)
pi_14 <- read.delim("subsampled_30_filtered.Upust.ESii_Pool5.genes.pi", header = F)
pi_15 <- read.delim("subsampled_30_filtered.Upust.ESii_Pool6.genes.pi", header = F)

# add column names
colnames(pi_1) <- c("ID", "num.snps", "frac", "pi_1")
colnames(pi_2) <- c("ID", "num.snps", "frac", "pi_2")
colnames(pi_3) <- c("ID", "num.snps", "frac", "pi_3")
colnames(pi_4) <- c("ID", "num.snps", "frac", "pi_4")
colnames(pi_5) <- c("ID", "num.snps", "frac", "pi_5")
colnames(pi_6) <- c("ID", "num.snps", "frac", "pi_6")
colnames(pi_7) <- c("ID", "num.snps", "frac", "pi_7")
colnames(pi_8) <- c("ID", "num.snps", "frac", "pi_8")
colnames(pi_9) <- c("ID", "num.snps", "frac", "pi_9")
colnames(pi_10) <- c("ID", "num.snps", "frac", "pi_10")
colnames(pi_11) <- c("ID", "num.snps", "frac", "pi_11")
colnames(pi_12) <- c("ID", "num.snps", "frac", "pi_12")
colnames(pi_13) <- c("ID", "num.snps", "frac", "pi_13")
colnames(pi_14) <- c("ID", "num.snps", "frac", "pi_14")
colnames(pi_15) <- c("ID", "num.snps", "frac", "pi_15")

#remove rows with na
pi_1 <- pi_1[!pi_1$pi_1 == "na",]
pi_2 <- pi_2[!pi_2$pi_2 == "na",]
pi_3 <- pi_3[!pi_3$pi_3 == "na",]
pi_4 <- pi_4[!pi_4$pi_4 == "na",]
pi_5 <- pi_5[!pi_5$pi_5 == "na",]
pi_6 <- pi_6[!pi_6$pi_6 == "na",]
pi_7 <- pi_7[!pi_7$pi_7 == "na",]
pi_8 <- pi_8[!pi_8$pi_8 == "na",]
pi_9 <- pi_9[!pi_9$pi_9 == "na",]
pi_10 <- pi_10[!pi_10$pi_10 == "na",]
pi_11 <- pi_11[!pi_11$pi_11 == "na",]
pi_12 <- pi_12[!pi_12$pi_12 == "na",]
pi_13 <- pi_13[!pi_13$pi_13 == "na",]
pi_14 <- pi_14[!pi_14$pi_14 == "na",]
pi_15 <- pi_15[!pi_15$pi_15 == "na",]

#subset dataset
pi_1 <- pi_1[,c(1,4)]
pi_2 <- pi_2[,c(1,4)]
pi_3 <- pi_3[,c(1,4)]
pi_4 <- pi_4[,c(1,4)]
pi_5 <- pi_5[,c(1,4)]
pi_6 <- pi_6[,c(1,4)]
pi_7 <- pi_7[,c(1,4)]
pi_8 <- pi_8[,c(1,4)]
pi_9 <- pi_9[,c(1,4)]
pi_10 <- pi_10[,c(1,4)]
pi_11 <- pi_11[,c(1,4)]
pi_12 <- pi_12[,c(1,4)]
pi_13 <- pi_13[,c(1,4)]
pi_14 <- pi_14[,c(1,4)]
pi_15 <- pi_15[,c(1,4)]

# merge only by gradient
pi_all1 <- merge(pi_1, pi_2, by="ID")
pi_all2 <- merge(pi_all1, pi_3, by="ID")
pi_all3 <- merge(pi_all2, pi_4, by="ID")
pi_all4 <- merge(pi_all3, pi_5, by="ID")
pi_IT <- merge(pi_all4, pi_6, by="ID")

pi_all1 <- merge(pi_7, pi_8, by="ID")
pi_ESi <- merge(pi_all1, pi_9, by="ID")

pi_all1 <- merge(pi_10, pi_11, by="ID")
pi_all2 <- merge(pi_all1, pi_12, by="ID")
pi_all3 <- merge(pi_all2, pi_13, by="ID")
pi_all4 <- merge(pi_all3, pi_14, by="ID")
pi_ESii <- merge(pi_all4, pi_15, by="ID")

# convert to numeric
pi_IT$pi_1 <- as.numeric(as.character(pi_IT$pi_1))
pi_IT$pi_2 <- as.numeric(as.character(pi_IT$pi_2))
pi_IT$pi_3 <- as.numeric(as.character(pi_IT$pi_3))
pi_IT$pi_4 <- as.numeric(as.character(pi_IT$pi_4))
pi_IT$pi_5 <- as.numeric(as.character(pi_IT$pi_5))
pi_IT$pi_6 <- as.numeric(as.character(pi_IT$pi_6))
pi_ESi$pi_7 <- as.numeric(as.character(pi_ESi$pi_7))
pi_ESi$pi_8 <- as.numeric(as.character(pi_ESi$pi_8))
pi_ESi$pi_9 <- as.numeric(as.character(pi_ESi$pi_9))
pi_ESii$pi_10 <- as.numeric(as.character(pi_ESii$pi_10))
pi_ESii$pi_11 <- as.numeric(as.character(pi_ESii$pi_11))
pi_ESii$pi_12 <- as.numeric(as.character(pi_ESii$pi_12))
pi_ESii$pi_13 <- as.numeric(as.character(pi_ESii$pi_13))
pi_ESii$pi_14 <- as.numeric(as.character(pi_ESii$pi_14))
pi_ESii$pi_15 <- as.numeric(as.character(pi_ESii$pi_15))

# prepare boxplots
colnames(pi_IT) <- c("ID","IT-1","IT-2","IT-3","IT-4","IT-5","IT-6")
colnames(pi_ESi) <- c("ID","ESi-1","ESi-2","ESi-3")
colnames(pi_ESii) <- c("ID","ESii-1","ESii-2","ESii-3","ESii-4","ESii-5","ESii-6")

melted1 <- melt(pi_IT, id = c("ID"), variable.name = "Pop", value.name = "pi")
melted2 <- melt(pi_ESi, id = c("ID"), variable.name = "Pop", value.name = "pi")
melted3 <- melt(pi_ESii, id = c("ID"), variable.name = "Pop", value.name = "pi")

melted4 <- rbind(melted1,melted2)
melted <- rbind(melted4,melted3)

melted[c('gradient', 'pop')] <- str_split_fixed(melted$Pop, '-', 2)

# get the means 
melted %>%
  group_by(gradient,pop) %>%
  summarise_at(vars(pi), list(name = mean))

#    gradient pop       name
#    <chr>    <chr>    <dbl>
#  1 ESi      1     0.00147 
#  2 ESi      2     0.000720
#  3 ESi      3     0.000649
#  4 ESii     1     0.00127 
#  5 ESii     2     0.00275 
#  6 ESii     3     0.00219 
#  7 ESii     4     0.000866
#  8 ESii     5     0.000956
#  9 ESii     6     0.000886
# 10 IT       1     0.00178 
# 11 IT       2     0.00282 
# 12 IT       3     0.00171 
# 13 IT       4     0.00167 
# 14 IT       5     0.00293 
# 15 IT       6     0.00173 

plot1  <- ggplot(melted, aes(x = pop, y = pi)) +
    geom_boxplot() +
    stat_summary(aes(x = as.numeric(as.character(pop))), fun=mean, colour="blue", geom="line") +
    scale_x_discrete("Population") +
    scale_y_continuous("Nucleotide diversity, all genes") +
    coord_cartesian(ylim = c(0, 0.1)) +
    facet_grid(.~gradient, scales = "free") +
    theme_classic(base_size = 12)

# save file and plot
suppressGraphics(ggsave('Upust_genes_pi_boxplots.png', plot1, width = 8, height = 4))
write.table(pi_all, file = "Upust_genes_all.pi", sep = "\t", na = "", quote = F, row.names = F, col.names = T)

```


```{r Upust, genes, Tajima's D}

library(tidyr)
library(ape)
library(vegan)
library(dplyr)
library(stringr)
library(reshape2)
library(ggplot2)
library(R.devices)

# import each tajD file 
tajD_1 <- read.delim("subsampled_30_filtered.Upust.IT_Pool1.genes.tajD", header = F)
tajD_2 <- read.delim("subsampled_30_filtered.Upust.IT_Pool2.genes.tajD", header = F)
tajD_3 <- read.delim("subsampled_30_filtered.Upust.IT_Pool3.genes.tajD", header = F)
tajD_4 <- read.delim("subsampled_30_filtered.Upust.IT_Pool4.genes.tajD", header = F)
tajD_5 <- read.delim("subsampled_30_filtered.Upust.IT_Pool5.genes.tajD", header = F)
tajD_6 <- read.delim("subsampled_30_filtered.Upust.IT_Pool6.genes.tajD", header = F)

tajD_7 <- read.delim("subsampled_30_filtered.Upust.ESi_Pool1.genes.tajD", header = F)
tajD_8 <- read.delim("subsampled_30_filtered.Upust.ESi_Pool2.genes.tajD", header = F)
tajD_9 <- read.delim("subsampled_30_filtered.Upust.ESi_Pool3.genes.tajD", header = F)

tajD_10 <- read.delim("subsampled_30_filtered.Upust.ESii_Pool1.genes.tajD", header = F)
tajD_11 <- read.delim("subsampled_30_filtered.Upust.ESii_Pool2.genes.tajD", header = F)
tajD_12 <- read.delim("subsampled_30_filtered.Upust.ESii_Pool3.genes.tajD", header = F)
tajD_13 <- read.delim("subsampled_30_filtered.Upust.ESii_Pool4.genes.tajD", header = F)
tajD_14 <- read.delim("subsampled_30_filtered.Upust.ESii_Pool5.genes.tajD", header = F)
tajD_15 <- read.delim("subsampled_30_filtered.Upust.ESii_Pool6.genes.tajD", header = F)

# add column names
colnames(tajD_1) <- c("ID", "num.snps", "frac", "tajD_1")
colnames(tajD_2) <- c("ID", "num.snps", "frac", "tajD_2")
colnames(tajD_3) <- c("ID", "num.snps", "frac", "tajD_3")
colnames(tajD_4) <- c("ID", "num.snps", "frac", "tajD_4")
colnames(tajD_5) <- c("ID", "num.snps", "frac", "tajD_5")
colnames(tajD_6) <- c("ID", "num.snps", "frac", "tajD_6")
colnames(tajD_7) <- c("ID", "num.snps", "frac", "tajD_7")
colnames(tajD_8) <- c("ID", "num.snps", "frac", "tajD_8")
colnames(tajD_9) <- c("ID", "num.snps", "frac", "tajD_9")
colnames(tajD_10) <- c("ID", "num.snps", "frac", "tajD_10")
colnames(tajD_11) <- c("ID", "num.snps", "frac", "tajD_11")
colnames(tajD_12) <- c("ID", "num.snps", "frac", "tajD_12")
colnames(tajD_13) <- c("ID", "num.snps", "frac", "tajD_13")
colnames(tajD_14) <- c("ID", "num.snps", "frac", "tajD_14")
colnames(tajD_15) <- c("ID", "num.snps", "frac", "tajD_15")

#remove rows with na
tajD_1 <- tajD_1[!tajD_1$tajD_1 == "na",]
tajD_2 <- tajD_2[!tajD_2$tajD_2 == "na",]
tajD_3 <- tajD_3[!tajD_3$tajD_3 == "na",]
tajD_4 <- tajD_4[!tajD_4$tajD_4 == "na",]
tajD_5 <- tajD_5[!tajD_5$tajD_5 == "na",]
tajD_6 <- tajD_6[!tajD_6$tajD_6 == "na",]
tajD_7 <- tajD_7[!tajD_7$tajD_7 == "na",]
tajD_8 <- tajD_8[!tajD_8$tajD_8 == "na",]
tajD_9 <- tajD_9[!tajD_9$tajD_9 == "na",]
tajD_10 <- tajD_10[!tajD_10$tajD_10 == "na",]
tajD_11 <- tajD_11[!tajD_11$tajD_11 == "na",]
tajD_12 <- tajD_12[!tajD_12$tajD_12 == "na",]
tajD_13 <- tajD_13[!tajD_13$tajD_13 == "na",]
tajD_14 <- tajD_14[!tajD_14$tajD_14 == "na",]
tajD_15 <- tajD_15[!tajD_15$tajD_15 == "na",]

#subset dataset
tajD_1 <- tajD_1[,c(1,4)]
tajD_2 <- tajD_2[,c(1,4)]
tajD_3 <- tajD_3[,c(1,4)]
tajD_4 <- tajD_4[,c(1,4)]
tajD_5 <- tajD_5[,c(1,4)]
tajD_6 <- tajD_6[,c(1,4)]
tajD_7 <- tajD_7[,c(1,4)]
tajD_8 <- tajD_8[,c(1,4)]
tajD_9 <- tajD_9[,c(1,4)]
tajD_10 <- tajD_10[,c(1,4)]
tajD_11 <- tajD_11[,c(1,4)]
tajD_12 <- tajD_12[,c(1,4)]
tajD_13 <- tajD_13[,c(1,4)]
tajD_14 <- tajD_14[,c(1,4)]
tajD_15 <- tajD_15[,c(1,4)]

# merge only by gradient
tajD_all1 <- merge(tajD_1, tajD_2, by="ID")
tajD_all2 <- merge(tajD_all1, tajD_3, by="ID")
tajD_all3 <- merge(tajD_all2, tajD_4, by="ID")
tajD_all4 <- merge(tajD_all3, tajD_5, by="ID")
tajD_IT <- merge(tajD_all4, tajD_6, by="ID")

tajD_all1 <- merge(tajD_7, tajD_8, by="ID")
tajD_ESi <- merge(tajD_all1, tajD_9, by="ID")

tajD_all1 <- merge(tajD_10, tajD_11, by="ID")
tajD_all2 <- merge(tajD_all1, tajD_12, by="ID")
tajD_all3 <- merge(tajD_all2, tajD_13, by="ID")
tajD_all4 <- merge(tajD_all3, tajD_14, by="ID")
tajD_ESii <- merge(tajD_all4, tajD_15, by="ID")

# convert to numeric
tajD_IT$tajD_1 <- as.numeric(as.character(tajD_IT$tajD_1))
tajD_IT$tajD_2 <- as.numeric(as.character(tajD_IT$tajD_2))
tajD_IT$tajD_3 <- as.numeric(as.character(tajD_IT$tajD_3))
tajD_IT$tajD_4 <- as.numeric(as.character(tajD_IT$tajD_4))
tajD_IT$tajD_5 <- as.numeric(as.character(tajD_IT$tajD_5))
tajD_IT$tajD_6 <- as.numeric(as.character(tajD_IT$tajD_6))
tajD_ESi$tajD_7 <- as.numeric(as.character(tajD_ESi$tajD_7))
tajD_ESi$tajD_8 <- as.numeric(as.character(tajD_ESi$tajD_8))
tajD_ESi$tajD_9 <- as.numeric(as.character(tajD_ESi$tajD_9))
tajD_ESii$tajD_10 <- as.numeric(as.character(tajD_ESii$tajD_10))
tajD_ESii$tajD_11 <- as.numeric(as.character(tajD_ESii$tajD_11))
tajD_ESii$tajD_12 <- as.numeric(as.character(tajD_ESii$tajD_12))
tajD_ESii$tajD_13 <- as.numeric(as.character(tajD_ESii$tajD_13))
tajD_ESii$tajD_14 <- as.numeric(as.character(tajD_ESii$tajD_14))
tajD_ESii$tajD_15 <- as.numeric(as.character(tajD_ESii$tajD_15))

# prepare boxplots
colnames(tajD_IT) <- c("ID","IT-1","IT-2","IT-3","IT-4","IT-5","IT-6")
colnames(tajD_ESi) <- c("ID","ESi-1","ESi-2","ESi-3")
colnames(tajD_ESii) <- c("ID","ESii-1","ESii-2","ESii-3","ESii-4","ESii-5","ESii-6")

melted1 <- melt(tajD_IT, id = c("ID"), variable.name = "Pop", value.name = "tajD")
melted2 <- melt(tajD_ESi, id = c("ID"), variable.name = "Pop", value.name = "tajD")
melted3 <- melt(tajD_ESii, id = c("ID"), variable.name = "Pop", value.name = "tajD")

melted4 <- rbind(melted1,melted2)
melted <- rbind(melted4,melted3)

melted[c('gradient', 'pop')] <- str_split_fixed(melted$Pop, '-', 2)

# get the means 
melted %>%
  group_by(gradient,pop) %>%
  summarise_at(vars(tajD), list(name = mean))

#    gradient pop      name
#    <chr>    <chr>   <dbl>
#  1 ESi      1     -0.158 
#  2 ESi      2     -0.0787
#  3 ESi      3      0.0501
#  4 ESii     1     -0.435 
#  5 ESii     2     -0.0229
#  6 ESii     3     -0.385 
#  7 ESii     4     -0.328 
#  8 ESii     5     -0.459 
#  9 ESii     6     -0.306 
# 10 IT       1      0.203 
# 11 IT       2     -0.869 
# 12 IT       3      0.126 
# 13 IT       4      0.254 
# 14 IT       5      1.10  
# 15 IT       6     -0.256 

plot1  <- ggplot(melted, aes(x = pop, y = tajD)) +
    geom_boxplot() +
    stat_summary(aes(x = as.numeric(as.character(pop))), fun=mean, colour="red", geom="line") +
    scale_x_discrete("Population") +
    scale_y_continuous("Tajima's D, all genes") +
    coord_cartesian(ylim = c(-3, 3)) +
    facet_grid(.~gradient, scales = "free") +
    theme_classic(base_size = 12)

# save file and plot
suppressGraphics(ggsave('Upust_genes_tajD_boxplots.png', plot1, width = 8, height = 4))
write.table(tajD_all, file = "Upust_genes_all.tajD", sep = "\t", na = "", quote = F, row.names = F, col.names = T)

```


```{r Uph, genes, nucleotide diversity}

library(tidyr)
library(ape)
library(vegan)
library(dplyr)
library(stringr)
library(reshape2)
library(ggplot2)
library(R.devices)

# import each pi file 
pi_16 <- read.delim("subsampled_30_filtered.Uph16.genes.pi", header = F)
pi_17 <- read.delim("subsampled_30_filtered.Uph17.genes.pi", header = F)
pi_18 <- read.delim("subsampled_30_filtered.Uph18.genes.pi", header = F)
pi_19 <- read.delim("subsampled_30_filtered.Uph19.genes.pi", header = F)

pi_22 <- read.delim("subsampled_30_filtered.Uph22.genes.pi", header = F)
pi_23 <- read.delim("subsampled_30_filtered.Uph23.genes.pi", header = F)
pi_24 <- read.delim("subsampled_30_filtered.Uph24.genes.pi", header = F)
pi_25 <- read.delim("subsampled_30_filtered.Uph25.genes.pi", header = F)
pi_26 <- read.delim("subsampled_30_filtered.Uph26.genes.pi", header = F)
pi_27 <- read.delim("subsampled_30_filtered.Uph27.genes.pi", header = F)
pi_28 <- read.delim("subsampled_30_filtered.Uph28.genes.pi", header = F)

# add column names
colnames(pi_16) <- c("ID", "num.snps", "frac", "pi_16")
colnames(pi_17) <- c("ID", "num.snps", "frac", "pi_17")
colnames(pi_18) <- c("ID", "num.snps", "frac", "pi_18")
colnames(pi_19) <- c("ID", "num.snps", "frac", "pi_19")

colnames(pi_22) <- c("ID", "num.snps", "frac", "pi_22")
colnames(pi_23) <- c("ID", "num.snps", "frac", "pi_23")
colnames(pi_24) <- c("ID", "num.snps", "frac", "pi_24")
colnames(pi_25) <- c("ID", "num.snps", "frac", "pi_25")
colnames(pi_26) <- c("ID", "num.snps", "frac", "pi_26")
colnames(pi_27) <- c("ID", "num.snps", "frac", "pi_27")
colnames(pi_28) <- c("ID", "num.snps", "frac", "pi_28")

#remove rows with na
pi_16 <- pi_16[!pi_16$pi_16 == "na",]
pi_17 <- pi_17[!pi_17$pi_17 == "na",]
pi_18 <- pi_18[!pi_18$pi_18 == "na",]
pi_19 <- pi_19[!pi_19$pi_19 == "na",]

pi_22 <- pi_22[!pi_22$pi_22 == "na",]
pi_23 <- pi_23[!pi_23$pi_23 == "na",]
pi_24 <- pi_24[!pi_24$pi_24 == "na",]
pi_25 <- pi_25[!pi_25$pi_25 == "na",]
pi_26 <- pi_26[!pi_26$pi_26 == "na",]
pi_27 <- pi_27[!pi_27$pi_27 == "na",]
pi_28 <- pi_28[!pi_28$pi_28 == "na",]

#subset dataset
pi_16 <- pi_16[,c(1,4)]
pi_17 <- pi_17[,c(1,4)]
pi_18 <- pi_18[,c(1,4)]
pi_19 <- pi_19[,c(1,4)]

pi_22 <- pi_22[,c(1,4)]
pi_23 <- pi_23[,c(1,4)]
pi_24 <- pi_24[,c(1,4)]
pi_25 <- pi_25[,c(1,4)]
pi_26 <- pi_26[,c(1,4)]
pi_27 <- pi_27[,c(1,4)]
pi_28 <- pi_28[,c(1,4)]

#merge by gradient
pi_all1 <- merge(pi_16, pi_17, by="ID")
pi_all2 <- merge(pi_all1, pi_18, by="ID")
pi_SN <- merge(pi_all2, pi_19, by="ID")

pi_all5 <- merge(pi_22, pi_23, by="ID")
pi_all6 <- merge(pi_all5, pi_24, by="ID")
pi_all7 <- merge(pi_all6, pi_25, by="ID")
pi_all8 <- merge(pi_all7, pi_26, by="ID")
pi_all9 <- merge(pi_all8, pi_27, by="ID")
pi_MJ<- merge(pi_all9, pi_28, by="ID")

# convert to numeric
pi_SN$pi_16 <- as.numeric(as.character(pi_SN$pi_16))
pi_SN$pi_17 <- as.numeric(as.character(pi_SN$pi_17))
pi_SN$pi_18 <- as.numeric(as.character(pi_SN$pi_18))
pi_SN$pi_19 <- as.numeric(as.character(pi_SN$pi_19))

pi_MJ$pi_22 <- as.numeric(as.character(pi_MJ$pi_22))
pi_MJ$pi_23 <- as.numeric(as.character(pi_MJ$pi_23))
pi_MJ$pi_24 <- as.numeric(as.character(pi_MJ$pi_24))
pi_MJ$pi_25 <- as.numeric(as.character(pi_MJ$pi_25))
pi_MJ$pi_26 <- as.numeric(as.character(pi_MJ$pi_26))
pi_MJ$pi_27 <- as.numeric(as.character(pi_MJ$pi_27))
pi_MJ$pi_28 <- as.numeric(as.character(pi_MJ$pi_28))

# prepare boxplots
colnames(pi_SN) <- c("ID","SN-1", "SN-2", "SN-3", "SN-4")
colnames(pi_MJ) <- c("ID","MJ-1","MJ-2","MJ-3","MJ-4","MJ-5","MJ-6","MJ-7")

melted1 <- melt(pi_SN, id = c("ID"), variable.name = "Pop", value.name = "pi")
melted2 <- melt(pi_MJ, id = c("ID"), variable.name = "Pop", value.name = "pi")

melted <- rbind(melted1,melted2)

melted[c('gradient', 'pop')] <- str_split_fixed(melted$Pop, '-', 2)

# get the means 
melted %>%
  group_by(gradient, pop) %>%
  summarise_at(vars(pi), list(name = mean))

#    gradient pop      name
#    <chr>    <chr>   <dbl>
#  1 MJ       1     0.0117 
#  2 MJ       2     0.0167 
#  3 MJ       3     0.0184 
#  4 MJ       4     0.0227 
#  5 MJ       5     0.00650
#  6 MJ       6     0.0108 
#  7 MJ       7     0.0143 
#  8 SN       1     0.0129 
#  9 SN       2     0.0152 
# 10 SN       3     0.0186 
# 11 SN       4     0.0158 

plot1  <- ggplot(melted, aes(x = pop, y = pi)) +
    geom_boxplot() +
    stat_summary(aes(x = as.numeric(as.character(pop))), fun=mean, colour="blue", geom="line") +
    scale_x_discrete("Population") +
    scale_y_continuous("Nucleotide diversity, all genes") +
    coord_cartesian(ylim = c(0, 0.1)) +
    facet_grid(.~gradient, scales = "free") +
    theme_classic(base_size = 12)

# save file and plot
suppressGraphics(ggsave('Uph_genes_pi_boxplots.png', plot1, width = 8, height = 4))
write.table(melted, file = "Uph_genes_pi_all.pi", sep = "\t", na = "", quote = F, row.names = F, col.names = T)

```

```{r Uph, genes, tajima's D}

library(tidyr)
library(ape)
library(vegan)
library(dplyr)
library(stringr)
library(reshape2)
library(ggplot2)
library(R.devices)

# import each tajD file 
tajD_16 <- read.delim("subsampled_30_filtered.Uph16.genes.tajD", header = F)
tajD_17 <- read.delim("subsampled_30_filtered.Uph17.genes.tajD", header = F)
tajD_18 <- read.delim("subsampled_30_filtered.Uph18.genes.tajD", header = F)
tajD_19 <- read.delim("subsampled_30_filtered.Uph19.genes.tajD", header = F)

tajD_22 <- read.delim("subsampled_30_filtered.Uph22.genes.tajD", header = F)
tajD_23 <- read.delim("subsampled_30_filtered.Uph23.genes.tajD", header = F)
tajD_24 <- read.delim("subsampled_30_filtered.Uph24.genes.tajD", header = F)
tajD_25 <- read.delim("subsampled_30_filtered.Uph25.genes.tajD", header = F)
tajD_26 <- read.delim("subsampled_30_filtered.Uph26.genes.tajD", header = F)
tajD_27 <- read.delim("subsampled_30_filtered.Uph27.genes.tajD", header = F)
tajD_28 <- read.delim("subsampled_30_filtered.Uph28.genes.tajD", header = F)

# add column names
colnames(tajD_16) <- c("ID", "num.snps", "frac", "tajD_16")
colnames(tajD_17) <- c("ID", "num.snps", "frac", "tajD_17")
colnames(tajD_18) <- c("ID", "num.snps", "frac", "tajD_18")
colnames(tajD_19) <- c("ID", "num.snps", "frac", "tajD_19")

colnames(tajD_22) <- c("ID", "num.snps", "frac", "tajD_22")
colnames(tajD_23) <- c("ID", "num.snps", "frac", "tajD_23")
colnames(tajD_24) <- c("ID", "num.snps", "frac", "tajD_24")
colnames(tajD_25) <- c("ID", "num.snps", "frac", "tajD_25")
colnames(tajD_26) <- c("ID", "num.snps", "frac", "tajD_26")
colnames(tajD_27) <- c("ID", "num.snps", "frac", "tajD_27")
colnames(tajD_28) <- c("ID", "num.snps", "frac", "tajD_28")

#remove rows with na
tajD_16 <- tajD_16[!tajD_16$tajD_16 == "na",]
tajD_17 <- tajD_17[!tajD_17$tajD_17 == "na",]
tajD_18 <- tajD_18[!tajD_18$tajD_18 == "na",]
tajD_19 <- tajD_19[!tajD_19$tajD_19 == "na",]

tajD_22 <- tajD_22[!tajD_22$tajD_22 == "na",]
tajD_23 <- tajD_23[!tajD_23$tajD_23 == "na",]
tajD_24 <- tajD_24[!tajD_24$tajD_24 == "na",]
tajD_25 <- tajD_25[!tajD_25$tajD_25 == "na",]
tajD_26 <- tajD_26[!tajD_26$tajD_26 == "na",]
tajD_27 <- tajD_27[!tajD_27$tajD_27 == "na",]
tajD_28 <- tajD_28[!tajD_28$tajD_28 == "na",]

#subset dataset
tajD_16 <- tajD_16[,c(1,4)]
tajD_17 <- tajD_17[,c(1,4)]
tajD_18 <- tajD_18[,c(1,4)]
tajD_19 <- tajD_19[,c(1,4)]

tajD_22 <- tajD_22[,c(1,4)]
tajD_23 <- tajD_23[,c(1,4)]
tajD_24 <- tajD_24[,c(1,4)]
tajD_25 <- tajD_25[,c(1,4)]
tajD_26 <- tajD_26[,c(1,4)]
tajD_27 <- tajD_27[,c(1,4)]
tajD_28 <- tajD_28[,c(1,4)]

#merge by gradient
tajD_all1 <- merge(tajD_16, tajD_17, by="ID")
tajD_all2 <- merge(tajD_all1, tajD_18, by="ID")
tajD_SN <- merge(tajD_all2, tajD_19, by="ID")

tajD_all5 <- merge(tajD_22, tajD_23, by="ID")
tajD_all6 <- merge(tajD_all5, tajD_24, by="ID")
tajD_all7 <- merge(tajD_all6, tajD_25, by="ID")
tajD_all8 <- merge(tajD_all7, tajD_26, by="ID")
tajD_all9 <- merge(tajD_all8, tajD_27, by="ID")
tajD_MJ<- merge(tajD_all9, tajD_28, by="ID")

# convert to numeric
tajD_SN$tajD_16 <- as.numeric(as.character(tajD_SN$tajD_16))
tajD_SN$tajD_17 <- as.numeric(as.character(tajD_SN$tajD_17))
tajD_SN$tajD_18 <- as.numeric(as.character(tajD_SN$tajD_18))
tajD_SN$tajD_19 <- as.numeric(as.character(tajD_SN$tajD_19))

tajD_MJ$tajD_22 <- as.numeric(as.character(tajD_MJ$tajD_22))
tajD_MJ$tajD_23 <- as.numeric(as.character(tajD_MJ$tajD_23))
tajD_MJ$tajD_24 <- as.numeric(as.character(tajD_MJ$tajD_24))
tajD_MJ$tajD_25 <- as.numeric(as.character(tajD_MJ$tajD_25))
tajD_MJ$tajD_26 <- as.numeric(as.character(tajD_MJ$tajD_26))
tajD_MJ$tajD_27 <- as.numeric(as.character(tajD_MJ$tajD_27))
tajD_MJ$tajD_28 <- as.numeric(as.character(tajD_MJ$tajD_28))

# prepare boxplots
colnames(tajD_SN) <- c("ID","SN-1", "SN-2", "SN-3", "SN-4")
colnames(tajD_MJ) <- c("ID","MJ-1","MJ-2","MJ-3","MJ-4","MJ-5","MJ-6","MJ-7")

melted1 <- melt(tajD_SN, id = c("ID"), variable.name = "Pop", value.name = "tajD")
melted2 <- melt(tajD_MJ, id = c("ID"), variable.name = "Pop", value.name = "tajD")

melted <- rbind(melted1,melted2)

melted[c('gradient', 'pop')] <- str_split_fixed(melted$Pop, '-', 2)

# get the means 
melted %>%
  group_by(gradient,pop) %>%
  summarise_at(vars(tajD), list(name = mean))

#    gradient pop      name
#    <chr>    <chr>   <dbl>
#  1 MJ       1     -0.778 
#  2 MJ       2     -0.344 
#  3 MJ       3      0.676 
#  4 MJ       4     -0.0558
#  5 MJ       5     -0.905 
#  6 MJ       6     -0.771 
#  7 MJ       7     -0.433 
#  8 SN       1     -0.431 
#  9 SN       2     -0.144 
# 10 SN       3     -0.0811
# 11 SN       4     -0.394 

plot1  <- ggplot(melted, aes(x = pop, y = tajD)) +
    geom_boxplot() +
    stat_summary(aes(x = as.numeric(as.character(pop))), fun=mean, colour="red", geom="line") +
    scale_x_discrete("Population") +
    scale_y_continuous("Tajima's D, all genes") +
    coord_cartesian(ylim = c(-3, 3)) +
    facet_grid(.~gradient, scales = "free") +
    theme_classic(base_size = 12)

# save file and plot
suppressGraphics(ggsave('Uph_genes_tajD_boxplots.png', plot1, width = 8, height = 4))
write.table(melted, file = "Uph_genes_tajD_all.tajD", sep = "\t", na = "", quote = F, row.names = F, col.names = T)

```
