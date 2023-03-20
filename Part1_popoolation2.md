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
- Upust_IT_Pool1-6: Mt. Limbara in Sardinia, Italy (IT1-4: Medit., IT5: mid, IT6: cold temperate)
- Upust_ESi_Pool1-3:
- Upust_ESii_Pool1-6: 

U. phaea:
- Uph16-19: Sierra Nevada (16&17: Medit., 18: mid, 19: cold temperate)
- Uph22-28: Mt. Jacinto (22-23-24: Medit., 25: mid, 26-28: cold temperate) 

U. hispanica:
- H1-6: bottom to H6:high elevation  


# Introduction: Re-running Pool-seq data analysis with popoolation2

Because The files are already trimmed, cleaned, and interleaved, steps 1-2 can be skipped. These files are located in /gendata_is/fdalgrande/Poolseq_Uphaea*/X201SC20030623-Z01-F001/raw_data/Uph*/ as of this writing (May 2022).

Tools needed to be installed for this process are:

bwa: conda install -c bioconda bwa
popoolation2: download and unzip popoolation2 from here: https://sourceforge.net/projects/popoolation2/
picard tools: conda install -c bioconda picard
samtools: conda install -c bioconda samtools


# 1. Initial cleanup and filtering of raw reads

Usually, we would begin by cleaning (trimming) and filtering the raw paired reads. However, this step has already been performed, so we include an example code for the sake of completion:

```{bash Example: cleanup and filtering}

nohup perl ~/Codes/Supplementary_FIle1-FastQFS.pl.pl -filtering Yes -fw P2_cleaned_R1_paired.fq -rw P2_cleaned_R2_paired.fq -fwo ESii_Pool2_Q26_L80_BQ3_R1_paired.fastq -rwo ESii_Pool2_Q26_L80_BQ3_R2_paired.fastq -sngl ESii_Pool2_Q26_L80_BQ3_singletons.fastq -mq 3 -q 26 -l 80 -plotting No -gsize 150 &

```


# 2. Create genome index

The only step needed is to copy the genome fasta file and index it. 

We choose here Umbilicaria_phaea_Uphaea_ref_cold.scaffolds.fa as our genome file, as U_phaea_TBG_1112 (high elevation/cold temperate) is the currently published U. phaea genome:

```{bash Genome index with bwa}

bwa index Umbilicaria_phaea_Uphaea_ref_cold.scaffolds.fa # file taken from /home/shared/Umbilicaria_2021_annotation/U_phaea/Uphaea_ref_cold-no-train/

```

Therefore, we can skip the usual next step, which is generating the interleaved reads. An example of the code is included here for completion:

```{bash Example: creating interleaved reads}

perl /opt/omega/scripts/shuffleSequences_fastq.pl Uph16_paired.R1.fastq Uph16_paired.R2.fastq Uph16_paired_interleaved.fq &

perl /opt/omega/scripts/shuffleSequences_fastq.pl ESii_Pool1_Q26_L80_BQ3_paired.tr_1.fastq ESii_Pool1_Q26_L80_BQ3_paired.tr_2.fastq ESii_Pool1_Q26_L80_BQ3_paired_interleaved.fq 

```


# 3. Add group reads to bam headers to be used in GATK

We can then proceed to the next step, which is to add the interleaved reads to bam headers. 

```{bash Create sam headers from interlieved reads}

bwa mem -M -t 12 -R '@RG\tID:group1\tSM:Uph16\tPL:illumina\tLB:lib1\tPU:Uph16' -p Umbilicaria_phaea_Uphaea_ref_cold.scaffolds.fa /gendata_is/fdalgrande/Poolseq_Uphaea_X201SC20030623-Z01-F001_raw_data/raw_data/Uph16/Uph16_FDSW202335553-1a_H33VHDSXY_L4_paired_interleaved.fq > Uph16_aligned_reads.sam &
bwa mem -M -t 12 -R '@RG\tID:group1\tSM:Uph17\tPL:illumina\tLB:lib1\tPU:Uph17' -p Umbilicaria_phaea_Uphaea_ref_cold.scaffolds.fa /gendata_is/fdalgrande/Poolseq_Uphaea_X201SC20030623-Z01-F001_raw_data/raw_data/Uph17/Uph17_FDSW202335554-1a_H33VHDSXY_L4_paired_interleaved.fq > Uph17_aligned_reads.sam &
bwa mem -M -t 12 -R '@RG\tID:group1\tSM:Uph18\tPL:illumina\tLB:lib1\tPU:Uph18' -p Umbilicaria_phaea_Uphaea_ref_cold.scaffolds.fa /gendata_is/fdalgrande/Poolseq_Uphaea_X201SC20030623-Z01-F001_raw_data/raw_data/Uph18/Uph18_FDSW202335555-1a_H33VHDSXY_L4_paired_interleaved.fq > Uph18_aligned_reads.sam &
bwa mem -M -t 12 -R '@RG\tID:group1\tSM:Uph19\tPL:illumina\tLB:lib1\tPU:Uph19' -p Umbilicaria_phaea_Uphaea_ref_cold.scaffolds.fa /gendata_is/fdalgrande/Poolseq_Uphaea_X201SC20030623-Z01-F001_raw_data/raw_data/Uph19/Uph19_FDSW202335556-1a_H33VHDSXY_L4_paired_interleaved.fq > Uph19_aligned_reads.sam &
bwa mem -M -t 12 -R '@RG\tID:group1\tSM:Uph22\tPL:illumina\tLB:lib1\tPU:Uph22' -p Umbilicaria_phaea_Uphaea_ref_cold.scaffolds.fa /gendata_is/fdalgrande/Poolseq_Uphaea_X201SC20030623-Z01-F001_raw_data/raw_data/Uph22/Uph22_FDSW202335557-1a_H33VHDSXY_L4_paired_interleaved.fq > Uph22_aligned_reads.sam &
bwa mem -M -t 12 -R '@RG\tID:group1\tSM:Uph23\tPL:illumina\tLB:lib1\tPU:Uph23' -p Umbilicaria_phaea_Uphaea_ref_cold.scaffolds.fa /gendata_is/fdalgrande/Poolseq_Uphaea_X201SC20030623-Z01-F001_raw_data/raw_data/Uph23/Uph23_FDSW202335558-1a_H33VHDSXY_L4_paired_interleaved.fq > Uph23_aligned_reads.sam &
bwa mem -M -t 12 -R '@RG\tID:group1\tSM:Uph24\tPL:illumina\tLB:lib1\tPU:Uph24' -p Umbilicaria_phaea_Uphaea_ref_cold.scaffolds.fa /gendata_is/fdalgrande/Poolseq_Uphaea_X201SC20030623-Z01-F001_raw_data/raw_data/Uph24/Uph24_FDSW202335559-1a_H33VHDSXY_L4_paired_interleaved.fq > Uph24_aligned_reads.sam &
bwa mem -M -t 12 -R '@RG\tID:group1\tSM:Uph25\tPL:illumina\tLB:lib1\tPU:Uph25' -p Umbilicaria_phaea_Uphaea_ref_cold.scaffolds.fa /gendata_is/fdalgrande/Poolseq_Uphaea_X201SC20030623-Z01-F001_raw_data/raw_data/Uph25/Uph25_FDSW202335560-1a_H33VHDSXY_L4_paired_interleaved.fq > Uph25_aligned_reads.sam &
bwa mem -M -t 12 -R '@RG\tID:group1\tSM:Uph26\tPL:illumina\tLB:lib1\tPU:Uph26' -p Umbilicaria_phaea_Uphaea_ref_cold.scaffolds.fa /gendata_is/fdalgrande/Poolseq_Uphaea_X201SC20030623-Z01-F001_raw_data/raw_data/Uph26/Uph26_FDSW202335561-1a_H33VHDSXY_L4_paired_interleaved.fq > Uph26_aligned_reads.sam &
bwa mem -M -t 12 -R '@RG\tID:group1\tSM:Uph27\tPL:illumina\tLB:lib1\tPU:Uph27' -p Umbilicaria_phaea_Uphaea_ref_cold.scaffolds.fa /gendata_is/fdalgrande/Poolseq_Uphaea_X201SC20030623-Z01-F001_raw_data/raw_data/Uph27/Uph27_FDSW202335562-1a_H33VHDSXY_L4_paired_interleaved.fq > Uph27_aligned_reads.sam &
bwa mem -M -t 12 -R '@RG\tID:group1\tSM:Uph28\tPL:illumina\tLB:lib1\tPU:Uph28' -p Umbilicaria_phaea_Uphaea_ref_cold.scaffolds.fa /gendata_is/fdalgrande/Poolseq_Uphaea_X201SC20030623-Z01-F001_raw_data/raw_data/Uph28/Uph28_FDSW202335563-1a_H33VHDSXY_L4_paired_interleaved.fq > Uph28_aligned_reads.sam &

# hispanica
bwa mem -M -t 12 -R '@RG\tID:group1\tSM:H1\tPL:illumina\tLB:lib1\tPU:H1' -p Lasallia_hispanica_U_hispanica_TBG_2337.scaffolds.fa H1_Q26_L80_BQ3_paired_interleaved.fq > H1_aligned_reads.sam

```

# 4. Prepare pooled populations for combining

Next, we convert the sam files to bam files, and extract only the mapped reads using a shell file (04_sam_to_bam.sh). 

In the following code, remember that a double ampersand (&&) tells the bash to run the code sequentially, and only if the previous step ran successfully. 

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

This way, when running the shell file (with the command "bash 04_sam_to_bam.sh") you will be prompted to provide a population name; you can give the prefix of the file (e.g. "Uph16") and the file will cycle through each command in turn. To create a shell file, simply open up the text editor (e.g. "vim 04_sam_to_bam.sh") and paste the above text inside the file.

An alternate version of this code is to simply run it with the following code:

```{bash Convert sam to bam and then extract only mapped reads (-F 4) (i.e. removing ambiguously mapped reads)}

### a. Uph16-19: Sierra Nevada (16&17: Medit., 18: mid, 19: cold temperate)
samtools view -Sb Uph19_aligned_reads.sam > Uph19_aligned_reads.bam &&
samtools view -h -F 4 -b Uph19_aligned_reads.bam > Uph19_aligned_reads_only_mapped.bam &&
java -Xmx8g -jar ~/tools/picard.jar SortSam -I Uph19_aligned_reads_only_mapped.bam -O Uph19_aligned_reads.sort.bam -VALIDATION_STRINGENCY SILENT -SO coordinate &&
java -Xmx8g -jar ~/tools/picard.jar MarkDuplicates -I Uph19_aligned_reads.sort.bam -O Uph19_aligned_reads.sort.rmd.bam -M Uph19_dupstat.txt -VALIDATION_STRINGENCY SILENT -REMOVE_DUPLICATES true &&
samtools view -q 20 -f 0x0002 -F 0x0004 -F 0x0008 -b Uph19_aligned_reads.sort.rmd.bam > Uph19_aligned_reads.sort.rmd.q20.bam &&

# Etc....

### Uph22-28: Mt. Jacinto (22-23-24: Medit., 25: mid, 26-28: cold temperate) 

samtools view -Sb Uph28_aligned_reads.sam > Uph28_aligned_reads.bam &&
samtools view -h -F 4 -b Uph28_aligned_reads.bam > Uph28_aligned_reads_only_mapped.bam &&
java -Xmx8g -jar ~/tools/picard.jar SortSam -I Uph28_aligned_reads_only_mapped.bam -O Uph28_aligned_reads.sort.bam -VALIDATION_STRINGENCY SILENT -SO coordinate &&
java -Xmx8g -jar ~/tools/picard.jar MarkDuplicates -I Uph28_aligned_reads.sort.bam -O Uph28_aligned_reads.sort.rmd.bam -M Uph28_dupstat.txt -VALIDATION_STRINGENCY SILENT -REMOVE_DUPLICATES true &&
samtools view -q 20 -f 0x0002 -F 0x0004 -F 0x0008 -b Uph28_aligned_reads.sort.rmd.bam > Uph28_aligned_reads.sort.rmd.q20.bam 

# Etc....


```


# 5. Generating the sync file as a primary popoolation2 input

The next step is to to combine the sorted and mapped reads into an mpileup file, which is used to speed up the process of creating the sync file, which is the main input format for popoolation2.

```{bash Pools are combined using mpileup}

samtools mpileup -B -Q 0 -f ../Umbilicaria_phaea_Uphaea_ref_cold.scaffolds.fa ../SierraNevada/Uph16_aligned_reads.sort.rmd.q20.bam ../SierraNevada/Uph17_aligned_reads.sort.rmd.q20.bam ../SierraNevada/Uph18_aligned_reads.sort.rmd.q20.bam ../SierraNevada/Uph19_aligned_reads.sort.rmd.q20.bam ../MtJacinto/Uph22_aligned_reads.sort.rmd.q20.bam ../MtJacinto/Uph23_aligned_reads.sort.rmd.q20.bam ../MtJacinto/Uph24_aligned_reads.sort.rmd.q20.bam ../MtJacinto/Uph25_aligned_reads.sort.rmd.q20.bam ../MtJacinto/Uph26_aligned_reads.sort.rmd.q20.bam ../MtJacinto/Uph27_aligned_reads.sort.rmd.q20.bam ../MtJacinto/Uph28_aligned_reads.sort.rmd.q20.bam> Uph_both_combined.mpileup &

samtools mpileup -B -Q 0 -f ../Umbilicaria_phaea_Uphaea_ref_cold.scaffolds.fa ../MtJacinto/Uph22_aligned_reads.sort.rmd.q20.bam ../MtJacinto/Uph23_aligned_reads.sort.rmd.q20.bam ../MtJacinto/Uph24_aligned_reads.sort.rmd.q20.bam ../MtJacinto/Uph25_aligned_reads.sort.rmd.q20.bam ../MtJacinto/Uph26_aligned_reads.sort.rmd.q20.bam ../MtJacinto/Uph27_aligned_reads.sort.rmd.q20.bam ../MtJacinto/Uph28_aligned_reads.sort.rmd.q20.bam > Uph_MJ_combined.mpileup &

# hispanica:
samtools mpileup -B -Q 0 -f Lasallia_hispanica_U_hispanica_TBG_2337.scaffolds.fa H1_aligned_reads.sort.rmd.q20.bam H2_aligned_reads.sort.rmd.q20.bam H3_aligned_reads.sort.rmd.q20.bam H4_aligned_reads.sort.rmd.q20.bam H5_aligned_reads.sort.rmd.q20.bam H6_aligned_reads.sort.rmd.q20.bam > Uhi_combined.mpileup &

```

Next, we identify genomic indel regions and remove them.

```{bash Genomic indel regions are identified}

perl /home/hvalim/tools/popoolation2_1201/indel_filtering/identify-indel-regions.pl --indel-window 5 --min-count 2 --input Uph_combined.mpileup --output Uph_combined.indel-regions.gtf &

# hispanica:
perl /home/hvalim/tools/popoolation2_1201/indel_filtering/identify-indel-regions.pl --indel-window 5 --min-count 2 --input Uhi_combined.mpileup --output Uhi_combined.indel-regions.gtf &


```


Then, we utilize these indel regions to filter the sync file. In popoolation1, this was done by filtering the mpileup file, and then creating the sync file; it isn't clear to me what difference it should make to filter the sync file after, rather than filtering the mpileup before generating the sync file, other than presumably a longer sync file generation time. Since mpileup->sync file conversion is apparently a bottle neck, it's not clear why this was changed in the new version of popoolation (ppopoolation2). Nonetheless, we will follow the popoolation2 pipeline and next generate the sync file:


```{bash mpileup file converted to sync file}

java -ea -Xmx400g -jar /home/hvalim/tools/popoolation2_1201/mpileup2sync.jar --input Uph_both_combined.mpileup --output Uph_both_combined.sync --fastq-type sanger --min-qual 20 --threads 70 &

# hispanica:
java -ea -Xmx400g -jar /home/hvalim/tools/popoolation2_1201/mpileup2sync.jar --input Uhi_combined.mpileup --output Uhi_combined.sync --fastq-type sanger --min-qual 20 --threads 70 &

```

Now, we can use the *indel-regions.gtf file to filter the sync file:

```{bash Indel regions are used to filter the mpileup file}

perl /home/hvalim/tools/popoolation2_1201/indel_filtering/filter-sync-by-gtf.pl --gtf Uph_both_combined.indel-regions.gtf --input Uph_both_combined.sync --output filtered.Uph_both_combined.sync &

# hispanica
perl /home/hvalim/tools/popoolation2_1201/indel_filtering/filter-sync-by-gtf.pl --gtf Uhi_combined.indel-regions.gtf --input Uhi_combined.sync --output filtered.Uhi_combined.sync &

```

Next, we subsample the sync file:

```{bash sync files are subsampled}

perl /home/hvalim/tools/popoolation2_1201/subsample-synchronized.pl --input filtered.Uph_both_combined.sync --output subsampled_30_filtered.Uph_both_combined.sync --target-coverage 30 --max-coverage 2% --method withoutreplace &

# hispanica
perl /home/hvalim/tools/popoolation2_1201/subsample-synchronized.pl --input filtered.Uhi_combined.sync --output subsampled_30_filtered.Uhi_combined.sync --target-coverage 30 --max-coverage 2% --method withoutreplace &

# better to increase the max coverage (i.e. decrease the number) than decrease coverage
# look at the sync file to find out what the coverage is at those positions
# also: generate non-subsampled Fst file to check for fixation at the extremes of each of those SNPs of interest

```

# 6. Calculating SNP frequency differences and FST sliding windows using the subsampled file

Now we can proceed to creating the outputs of Popoolation2, namely SNP frequency differences and Fst scores.

The first step is to calculate SNP frequencies, which yields two files with the extensions: "_rc" and "_pwc".

    _rc: this file contains the major and minor alleles for every SNP in a concise format
    _pwc: this file contains the differences in allele frequencies for every pairwise comparision of the populations present in the synchronized file
    For details see the man pages of the script

```{bash SNP frequencies}

nohup perl /home/hvalim/tools/popoolation2_1201/snp-frequency-diff.pl --input subsampled_30_filtered.Uph_both_combined.sync --output-prefix subsampled_30_filtered.Uph_both_combined --max-coverage 100 --min-count 3 &

# hispanica
nohup perl /home/hvalim/tools/popoolation2_1201/snp-frequency-diff.pl --input subsampled_30_filtered.Uhi_combined.sync --output-prefix subsampled_30_filtered.Uhi_combined --max-coverage 100 --min-count 3 &

```

Next, we proceed to calculating Fst values using a sliding window approach. Because we want to get the Fst value for each individual SNP, we use a window size of 1; this automatically suppresses any output for bases that have no SNPs (otherwise, you would presumably have windows in your output with no differences between any of the populations).

```{bash Fst sliding windows}

nohup perl /home/hvalim/tools/popoolation2_1201/fst-sliding.pl --input subsampled_30_filtered.Uph_both_combined.sync --output subsampled_30_filtered.Uph_both_combined.fst --suppress-noninformative --min-count 3 --max-coverage 100 --min-covered-fraction 1 --window-size 1 --step-size 1 --pool-size 100 &

# hispanica
nohup perl /home/hvalim/tools/popoolation2_1201/fst-sliding.pl --input subsampled_30_filtered.Uhi_combined.sync --output subsampled_30_filtered.Uhi_combined.fst --suppress-noninformative --min-count 3 --max-coverage 100 --min-covered-fraction 1 --window-size 1 --step-size 1 --pool-size 100 &

```

Finally, we use a subsetted gtf file (e.g. with only genes, exons, or some other set of features) to create a gene-wise sync file. We then calculate SNP-wise Fst values for only genic SNPs, and whole-gene Fst values.

```{bash make gene-wise sync files and calculate their Fst values}

# U. phaea:

## prepare gene-wise sync file
nohup perl /home/hvalim/tools/popoolation2_1201/create-genewise-sync.pl --input subsampled_30_filtered.Uph_both_combined.sync --gtf ../Uphaea_TBG_1112_all_genes.gff3 --output New_Poolseq_assembly/GENES_subsampled_30_filtered.Uph_both_combined.sync

## calculate Fst values for only genic SNPs
perl /home/hvalim/tools/popoolation2_1201/fst-sliding.pl --suppress-noninformative --min-count 3 --max-coverage 100 --pool-size 100 --min-covered-fraction 1 --window-size 1 --step-size 1 --input GENES_subsampled_30_filtered.Uph_both_combined.sync --output GENES_subsampled_30_filtered.Uph_both_combined.fst

## calculate whole-gene Fst values
perl /home/hvalim/tools/popoolation2_1201/fst-sliding.pl --min-count 3 --max-coverage 100 --pool-size 100 --min-covered-fraction 1 --window-size 1000000 --step-size 1000000 --input GENES_subsampled_30_filtered.Uph_both_combined.sync --output GENES_subsampled_30_filtered.Uph_both_combined_whole_gene.fst

# U. hispanica:

## prepare gene-wise sync file
nohup perl /home/hvalim/tools/popoolation2_1201/create-genewise-sync.pl --input subsampled_30_filtered.Uhi_combined.sync --gtf Uhispanica_TBG_2337_all_genes.gtf --output GENES_subsampled_30_filtered.Uhi_combined.sync

## calculate Fst values for only genic SNPs
perl /home/hvalim/tools/popoolation2_1201/fst-sliding.pl --suppress-noninformative --input GENES_subsampled_30_filtered.Uhi_combined.sync --output GENES_subsampled_30_filtered.Uhi_combined.fst --min-count 3 --max-coverage 100 --pool-size 100 --min-covered-fraction 1 --window-size 1 --step-size 1 

## calculate whole-gene Fst values
perl /home/hvalim/tools/popoolation2_1201/fst-sliding.pl --min-count 3 --max-coverage 100 --pool-size 100 --min-covered-fraction 1 --window-size 1000000 --step-size 1000000 --input GENES_subsampled_30_filtered.Uhi_combined.sync --output GENES_subsampled_30_filtered.Uhi_combined_whole_gene.fst

```
