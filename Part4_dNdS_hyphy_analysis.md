---
title: 'Part 4: dN/dS analysis with HyPhy'
author: "Henrique Valim"
date: "2022-12-20"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# dN/dS analyses of gradients

Begin by extracting just the reads for the clock genes from the reference genomes.

Then we make an indexed "genome" of just the clock genes, and then follow the code to map and align these separately using samtools and bcftools

  NOTE: you may need to manually re-name the genes to make sure there are no overlaps (i.e. same chromosome and no gene name); here the names match Table S2 in the manuscript.

## 1. Index set of circadian homologs (in conda env poolseq)

```{bash}

conda activate poolseq

bwa index Uph_circadian.fasta # etc
```
	

## 2. Create SAM files

```{bash}

### for U. pustulata circadian loci:
for i in 1 2 3 4 5 6; do bwa mem -M -t 12 -R '@RG\tID:group1\tSM:Upust_IT${i}\tPL:illumina\tLB:lib1\tPU:Upust_IT${i}' -p Upust_circadian.fasta /phylodata/fdalgrande/TRIMMED_PoolSeq_3gradients/IT_Pool${i}_paired_interleaved.fq > Upust_IT${i}_circadian_aligned_reads.sam; done

for i in 1 2 3 4 5 6; do bwa mem -M -t 12 -R '@RG\tID:group1\tSM:Upust_ESii${i}\tPL:illumina\tLB:lib1\tPU:Upust_ESii${i}' -p Upust_circadian.fasta /phylodata/fdalgrande/TRIMMED_PoolSeq_3gradients/ESii_Pool${i}_paired_interleaved.fq > Upust_ESii${i}_circadian_aligned_reads.sam; done

for i in 1 2 3; do bwa mem -M -t 12 -R '@RG\tID:group1\tSM:Upust_ESi${i}\tPL:illumina\tLB:lib1\tPU:Upust_ESi${i}' -p Upust_circadian.fasta /phylodata/fdalgrande/TRIMMED_PoolSeq_3gradients/ESi_Pool${i}_paired_interleaved.fq > Upust_ESi${i}_circadian_aligned_reads.sam; done


### for U. phaea circadian loci:
for i in 16 17 18 19 22 23 24 25 26 27 28; do bwa mem -M -t 12 -R '@RG\tID:group1\tSM:Uph${i}\tPL:illumina\tLB:lib1\tPU:Uph${i}' -p Uph_circadian.fasta /gendata_is/fdalgrande/Poolseq_Uphaea_X201SC20030623-Z01-F001_raw_data/raw_data/Uph${i}/Uph${i}_*_L4_paired_interleaved.fq > Uph${i}_circadian_aligned_reads.sam; done



### for U. pustulata temperature loci:
for i in 1 2 3 4 5 6; do bwa mem -M -t 12 -R '@RG\tID:group1\tSM:Upust_IT${i}\tPL:illumina\tLB:lib1\tPU:Upust_IT${i}' -p Upust_temperature.fasta /phylodata/fdalgrande/TRIMMED_PoolSeq_3gradients/IT_Pool${i}_paired_interleaved.fq > Upust_IT${i}_temperature_aligned_reads.sam; done

for i in 1 2 3 4 5 6; do bwa mem -M -t 12 -R '@RG\tID:group1\tSM:Upust_ESii${i}\tPL:illumina\tLB:lib1\tPU:Upust_ESii${i}' -p Upust_temperature.fasta /phylodata/fdalgrande/TRIMMED_PoolSeq_3gradients/ESii_Pool${i}_paired_interleaved.fq > Upust_ESii${i}_temperature_aligned_reads.sam; done

for i in 1 2 3; do bwa mem -M -t 12 -R '@RG\tID:group1\tSM:Upust_ESi${i}\tPL:illumina\tLB:lib1\tPU:Upust_ESi${i}' -p Upust_temperature.fasta /phylodata/fdalgrande/TRIMMED_PoolSeq_3gradients/ESi_Pool${i}_paired_interleaved.fq > Upust_ESi${i}_temperature_aligned_reads.sam; done


### for U. phaea temperature loci:
for i in 16 17 18 19 22 23 24 25 26 27 28; do bwa mem -M -t 12 -R '@RG\tID:group1\tSM:Uph${i}\tPL:illumina\tLB:lib1\tPU:Uph${i}' -p Uph_temperature.fasta /gendata_is/fdalgrande/Poolseq_Uphaea_X201SC20030623-Z01-F001_raw_data/raw_data/Uph${i}/Uph${i}_*_L4_paired_interleaved.fq > Uph${i}_temperature_aligned_reads.sam; done

```
	
## 3. Convert to BAM and extract only mapped reads
	
```{bash}
conda activate bcftools

### for U. pustulata circadian loci:
for pop in IT1 IT2 IT3 IT4 IT5 IT6 ESii1 ESii2 ESii3 ESii4 ESii5 ESii6 ESi1 ESi2 ESi3; do samtools view -Sb Upust_${pop}_circadian_aligned_reads.sam > Upust_${pop}_circadian_aligned_reads.bam; done

for pop in IT1 IT2 IT3 IT4 IT5 IT6 ESii1 ESii2 ESii3 ESii4 ESii5 ESii6 ESi1 ESi2 ESi3; do samtools view -h -F 4 -b Upust_${pop}_circadian_aligned_reads.bam > Upust_${pop}_circadian_only_mapped.bam; done

for pop in IT1 IT2 IT3 IT4 IT5 IT6 ESii1 ESii2 ESii3 ESii4 ESii5 ESii6 ESi1 ESi2 ESi3; do java -Xmx8g -jar ~/tools/picard.jar SortSam -I Upust_${pop}_circadian_only_mapped.bam -O Upust_${pop}_circadian_only_mapped.sort.bam -VALIDATION_STRINGENCY SILENT -SO coordinate; done

for pop in IT1 IT2 IT3 IT4 IT5 IT6 ESii1 ESii2 ESii3 ESii4 ESii5 ESii6 ESi1 ESi2 ESi3; do java -Xmx8g -jar ~/tools/picard.jar MarkDuplicates -I Upust_${pop}_circadian_only_mapped.sort.bam  -O Upust_${pop}_circadian_only_mapped.sort.rmd.bam -M Upust_${pop}_circadian_dupstat.txt -VALIDATION_STRINGENCY SILENT -REMOVE_DUPLICATES true; done

for pop in IT1 IT2 IT3 IT4 IT5 IT6 ESii1 ESii2 ESii3 ESii4 ESii5 ESii6 ESi1 ESi2 ESi3; do samtools view -q 20 -f 0x0002 -F 0x0004 -F 0x0008 -b Upust_${pop}_circadian_only_mapped.sort.rmd.bam > Upust_${pop}_circadian_only_mapped.sort.rmd.q20.bam; done



### for U. pustulata temperature-associated loci:
for pop in IT1 IT2 IT3 IT4 IT5 IT6 ESii1 ESii2 ESii3 ESii4 ESii5 ESii6 ESi1 ESi2 ESi3; do samtools view -Sb Upust_${pop}_temperature_aligned_reads.sam > Upust_${pop}_temperature_aligned_reads.bam; done

for pop in IT1 IT2 IT3 IT4 IT5 IT6 ESii1 ESii2 ESii3 ESii4 ESii5 ESii6 ESi1 ESi2 ESi3; do samtools view -h -F 4 -b Upust_${pop}_temperature_aligned_reads.bam > Upust_${pop}_temperature_only_mapped.bam; done

for pop in IT1 IT2 IT3 IT4 IT5 IT6 ESii1 ESii2 ESii3 ESii4 ESii5 ESii6 ESi1 ESi2 ESi3; do java -Xmx8g -jar ~/tools/picard.jar SortSam -I Upust_${pop}_temperature_only_mapped.bam -O Upust_${pop}_temperature_only_mapped.sort.bam -VALIDATION_STRINGENCY SILENT -SO coordinate; done

for pop in IT1 IT2 IT3 IT4 IT5 IT6 ESii1 ESii2 ESii3 ESii4 ESii5 ESii6 ESi1 ESi2 ESi3; do java -Xmx8g -jar ~/tools/picard.jar MarkDuplicates -I Upust_${pop}_temperature_only_mapped.sort.bam  -O Upust_${pop}_temperature_only_mapped.sort.rmd.bam -M Upust_${pop}_temperature_dupstat.txt -VALIDATION_STRINGENCY SILENT -REMOVE_DUPLICATES true; done

for pop in IT1 IT2 IT3 IT4 IT5 IT6 ESii1 ESii2 ESii3 ESii4 ESii5 ESii6 ESi1 ESi2 ESi3; do samtools view -q 20 -f 0x0002 -F 0x0004 -F 0x0008 -b Upust_${pop}_temperature_only_mapped.sort.rmd.bam > Upust_${pop}_temperature_only_mapped.sort.rmd.q20.bam; done




### for U. phaea circadian loci:
for i in 16 17 18 19 22 23 24 25 26 27 28; do samtools view -Sb Uph${i}_circadian_aligned_reads.sam > Uph${i}_circadian_aligned_reads.bam; done

for i in 16 17 18 19 22 23 24 25 26 27 28; do samtools view -h -F 4 -b Uph${i}_circadian_aligned_reads.bam > Uph${i}_circadian_only_mapped.bam; done

for i in 16 17 18 19 22 23 24 25 26 27 28; do java -Xmx8g -jar ~/tools/picard.jar SortSam -I Uph${i}_circadian_only_mapped.bam -O Uph${i}_circadian_only_mapped.sort.bam -VALIDATION_STRINGENCY SILENT -SO coordinate; done

for i in 16 17 18 19 22 23 24 25 26 27 28; do java -Xmx8g -jar ~/tools/picard.jar MarkDuplicates -I Uph${i}_circadian_only_mapped.sort.bam  -O Uph${i}_circadian_only_mapped.sort.rmd.bam -M Uph${i}_circadian_dupstat.txt -VALIDATION_STRINGENCY SILENT -REMOVE_DUPLICATES true; done

for i in 16 17 18 19 22 23 24 25 26 27 28; do samtools view -q 20 -f 0x0002 -F 0x0004 -F 0x0008 -b Uph${i}_circadian_only_mapped.sort.rmd.bam > Uph${i}_circadian_only_mapped.sort.rmd.q20.bam; done





### for U. phaea temperature-associated loci:
for i in 16 17 18 19 22 23 24 25 26 27 28; do samtools view -Sb Uph${i}_temperature_aligned_reads.sam > Uph${i}_temperature_aligned_reads.bam; done

for i in 16 17 18 19 22 23 24 25 26 27 28; do samtools view -h -F 4 -b Uph${i}_temperature_aligned_reads.bam > Uph${i}_temperature_only_mapped.bam; done

for i in 16 17 18 19 22 23 24 25 26 27 28; do java -Xmx8g -jar ~/tools/picard.jar SortSam -I Uph${i}_temperature_only_mapped.bam -O Uph${i}_temperature_only_mapped.sort.bam -VALIDATION_STRINGENCY SILENT -SO coordinate; done

for i in 16 17 18 19 22 23 24 25 26 27 28; do java -Xmx8g -jar ~/tools/picard.jar MarkDuplicates -I Uph${i}_temperature_only_mapped.sort.bam  -O Uph${i}_temperature_only_mapped.sort.rmd.bam -M Uph${i}_temperature_dupstat.txt -VALIDATION_STRINGENCY SILENT -REMOVE_DUPLICATES true; done

for i in 16 17 18 19 22 23 24 25 26 27 28; do samtools view -q 20 -f 0x0002 -F 0x0004 -F 0x0008 -b Uph${i}_temperature_only_mapped.sort.rmd.bam > Uph${i}_temperature_only_mapped.sort.rmd.q20.bam; done

```


## 4. Convert to mpileup file, then to FASTQ, then to FASTA (in conda env sambcfenv)

```{bash}

conda activate sambcfenv

# for U. pustulata circadian genes:

for pop in IT1 IT2 IT3 IT4 IT5 IT6 ESii1 ESii2 ESii3 ESii4 ESii5 ESii6 ESi1 ESi2 ESi3; do bcftools mpileup --fasta-ref Upust_circadian.fasta Upust_${pop}_circadian_only_mapped.sort.rmd.q20.bam | bcftools call -c | vcfutils.pl vcf2fq > Upust_${pop}_circadian.fastq; done

for pop in IT1 IT2 IT3 IT4 IT5 IT6 ESii1 ESii2 ESii3 ESii4 ESii5 ESii6 ESi1 ESi2 ESi3; do seqtk seq -aQ64 -q20 -n N Upust_${pop}_circadian.fastq > Upust_${pop}_circadian.fasta; done


# for U. pustulata temperature genes:

for pop in IT1 IT2 IT3 IT4 IT5 IT6 ESii1 ESii2 ESii3 ESii4 ESii5 ESii6 ESi1 ESi2 ESi3; do bcftools mpileup --fasta-ref Upust_temperature.fasta Upust_${pop}_temperature_only_mapped.sort.rmd.q20.bam | bcftools call -c | vcfutils.pl vcf2fq > Upust_${pop}_temperature.fastq; done

for pop in IT1 IT2 IT3 IT4 IT5 IT6 ESii1 ESii2 ESii3 ESii4 ESii5 ESii6 ESi1 ESi2 ESi3; do seqtk seq -aQ64 -q20 -n N Upust_${pop}_temperature.fastq > Upust_${pop}_temperature.fasta; done

 
 
# for U. phaea circadian genes:
for i in 16 17 18 19 22 23 24 25 26 27 28; do bcftools mpileup --fasta-ref Uph_circadian.fasta Uph${i}_circadian_only_mapped.sort.rmd.q20.bam | bcftools call -c | vcfutils.pl vcf2fq > Uph${i}_circadian.fastq; done

for i in 16 17 18 19 22 23 24 25 26 27 28; do seqtk seq -aQ64 -q20 -n N Uph${i}_circadian.fastq > Uph${i}_circadian.fasta; done



# for U. phaea temperature genes:
for i in 16 17 18 19 22 23 24 25 26 27 28; do bcftools mpileup --fasta-ref Uph_temperature.fasta Uph${i}_temperature_only_mapped.sort.rmd.q20.bam | bcftools call -c | vcfutils.pl vcf2fq > Uph${i}_temperature.fastq; done

for i in 16 17 18 19 22 23 24 25 26 27 28; do seqtk seq -aQ64 -q20 -n N Uph${i}_temperature.fastq > Uph${i}_temperature.fasta; done


```


## 5. Sort these by gene (labeling and combining the copies from each population) and then perform an alignment and generate a tree using RAxML. Export the tree as *.nex (combined alignment and tree)

Add labels for each population to end of each gene per pop with sed:

```{bash}

# for the population chimeric sequences:

for i in 16 17 18 19 22 23 24 25 26 27 28; do sed -i "/^>/ s/$/_Uph${i}/" Uph${i}_circadian.fasta; done

for i in 16 17 18 19 22 23 24 25 26 27 28; do sed -i "/^>/ s/$/_Uph${i}/" Uph${i}_temperature.fasta; done

for pop in IT1 IT2 IT3 IT4 IT5 IT6 ESii1 ESii2 ESii3 ESii4 ESii5 ESii6 ESi1 ESi2 ESi3; do sed -i "/^>/ s/$/_Upust${pop}/" Upust_${pop}_circadian.fasta; done

for pop in IT1 IT2 IT3 IT4 IT5 IT6 ESii1 ESii2 ESii3 ESii4 ESii5 ESii6 ESi1 ESi2 ESi3; do sed -i "/^>/ s/$/_Upust${pop}/" Upust_${pop}_temperature.fasta; done

# for the reference sequences:

sed -i "/^>/ s/$/_Uph-cold/" Uph_circadian.fasta

sed -i "/^>/ s/$/_Uph-cold/" Uph_temperature.fasta

sed -i "/^>/ s/$/_Upust-cold/" Upust_circadian.fasta

sed -i "/^>/ s/$/_Upust-cold/" Upust_temperature.fasta
```

We use these gene names to create a nested for loop that will extract all genes for all populations separately and then combine them with the cold reference copy.

First for U. phaea:

```{bash}

# for circadian genes:

for gene in ACR2 ccg-1 ccg-9 CDC55 CHO2 CKB1 CKB2 con-6 con-10 csn-1 csn-2 csn-4 csn-6 csn-7a csp-2 DUN1 ecm33 FUN30 fwd-1 GAT1 GDP1-2 GLC7 hrp3 HRR25-1 HRR25-3 ISW2 KIN28 MgSsk1 MTR4 OPI1 PPG1 PPH3 PPT1 pzh1 ras1 ras2 RPN11 RRI1 RSR1 SEC4 SKI2-2 SNF2 SOD1 SSN3 TRR1 UME6 upf-1 VAS1 wc-1 wc-2 YCK2
  do for i in 16 17 18 19 22 23 24 25 26 27 28
    do seqkit grep -r -p "${gene}" Uph${i}_circadian.fasta > Uph${i}_${gene}.fasta
    cat Uph*_${gene}.fasta > Uph_${gene}_combined.fasta
    seqkit grep -r -p "${gene}" Uph_circadian.fasta > Uph_cold_${gene}.fasta | cat Uph_cold_${gene}.fasta Uph_${gene}_combined.fasta > Uph_${gene}_combined1.fasta
done
done

mv Uph*combined1.fasta 01b_input_data/


# for temperature genes:

for gene in CCP1 DED1 hsp30 KIN28 LSP1 PRP28 PYC1 cut-1b cat-4 GET3 MKC1 NTH1_1 PRF1 DBP3 LSP1b cr-1 KAR2 PAN2 BUR1 gcy-3 HSP82 PHO85 SSB1 cat1b xyl1 DBP4 DBP10 wos2 gcy-1a HSP78_1 PRP5 CAT1 DBP2 SFL1 CTA8 gcy-1b SSA2 FBP1 DES1 SSC1 GSY1 LSK1 MRH4 dak1 HSP98 cut-1a GPA1 GPA1_2 GPH1 SSE1 NTH1_1
  do for i in 16 17 18 19 22 23 24 25 26 27 28
    do seqkit grep -r -p "${gene}" Uph${i}_temperature.fasta > Uph${i}_${gene}.fasta
    cat Uph*_${gene}.fasta > Uph_${gene}_combined.fasta
    seqkit grep -r -p "${gene}" Uph_temperature.fasta > Uph_cold_${gene}.fasta | cat Uph_cold_${gene}.fasta Uph_${gene}_combined.fasta > Uph_${gene}_combined1.fasta
done
done

mv Uph*combined1.fasta 01b_input_data/

```

Then for U. pustulata:

```{bash}

# for circadian genes:

for gene in acr-2 adv-1 ccg-1 ccg-7 ccg-9 Chd3 chol-1 ck-1a cka-1 cka-2 ckb-1-1 ckb-1-2 con-6 con-10 crf6-1 crf10-1 csn-1 csn-2 csn-4 csn-5-1 csn-5-2 csn-6 csn-7a csp-2 cys-9 ecm33 frh frh-2 frq fwd-1 isw-1 nit-2 pph-3-1 pph-3-2 pph-3-3 pph-3-4 prd-4-1 prd-4-2 ras-1-1 ras-1-2 ras-1-3 ras-2 rent-1 rgb-1 rrg-1 sod-1 vvd-1 vvd-2 wc-1 wc-2
  do for pop in IT1 IT2 IT3 IT4 IT5 IT6 ESii1 ESii2 ESii3 ESii4 ESii5 ESii6 ESi1 ESi2 ESi3
    do seqkit grep -r -p "${gene}" Upust_${pop}_circadian.fasta > Upust_${pop}_${gene}.fasta
    cat Upust*_${gene}.fasta > Upust_${gene}_combined.fasta
    seqkit grep -r -p "${gene}" Upust_circadian.fasta > Upust_cold_${gene}.fasta | cat Upust_cold_${gene}.fasta Upust_${gene}_combined.fasta > Upust_${gene}_combined1.fasta
done
done

mv Upust*combined1.fasta 01_input_data/

# for temperature genes:

for gene in hog1 cut-1b xyl1 hsp90 GNAQ hsp30 SSB1 fbp-1 hsp70b treB bipA acu-6 cr-1 hsf-2 hsp70 cat-4 hsp98 gcy-3 mrh-4 Pfd get3 gcy-1 Q1K6J1 LSP1 hsf-1 GNA3 GlyP cat-4b dak-1 AGA1 suc katG hsp88 cut-1 GNA2 hsp90a gsy-1 
  do for pop in IT1 IT2 IT3 IT4 IT5 IT6 ESii1 ESii2 ESii3 ESii4 ESii5 ESii6 ESi1 ESi2 ESi3
    do seqkit grep -r -p "${gene}" Upust_${pop}_temperature.fasta > Upust_${pop}_${gene}.fasta
    cat Upust*_${gene}.fasta > Upust_${gene}_combined.fasta
    seqkit grep -r -p "${gene}" Upust_temperature.fasta > Upust_cold_${gene}.fasta | cat Upust_cold_${gene}.fasta Upust_${gene}_combined.fasta > Upust_${gene}_combined1.fasta
done
done

mv Upust*combined1.fasta 01_input_data/

```


Next we perform the alignment and tree building

For U. phaea:

```{bash}

# for circadian genes

for gene in ACR2 ccg-1 ccg-9 CDC55 CHO2 CKB1 CKB2 con-6 con-10 csn-1 csn-2 csn-4 csn-6 csn-7a csp-2 DUN1 ecm33 FUN30 fwd-1 GAT1 GDP1-2 GLC7 hrp3 HRR25-1 HRR25-3 ISW2 KIN28 MgSsk1 MTR4 OPI1 PPG1 PPH3 PPT1 pzh1 ras1 ras2 RPN11 RRI1 RSR1 SEC4 SKI2-2 SNF2 SOD1 SSN3 TRR1 UME6 upf-1 VAS1 wc-1 wc-2 YCK2
  do mafft Uph_${gene}_combined1.fasta > Uph_${gene}_combined1.msa
  raxmlHPC -m GTRGAMMA  -s Uph_${gene}_combined1.msa -n Uph_${gene}_combined1_tree
done

# for temperature genes

for gene in CCP1 DED1 hsp30 KIN28 LSP1 PRP28 PYC1 cut-1b cat-4 GET3 MKC1 NTH1_1 PRF1 DBP3 LSP1b cr-1 KAR2 PAN2 BUR1 gcy-3 HSP82 PHO85 SSB1 cat1b xyl1 DBP4 DBP10 wos2 gcy-1a HSP78_1 PRP5 CAT1 DBP2 SFL1 CTA8 gcy-1b SSA2 FBP1 DES1 SSC1 GSY1 LSK1 MRH4 dak1 HSP98 cut-1a GPA1 GPA1_2 GPH1 SSE1 NTH1_1
  do mafft Uph_${gene}_combined1.fasta > Uph_${gene}_combined1.msa
  raxmlHPC -m GTRGAMMA  -s Uph_${gene}_combined1.msa -n Uph_${gene}_combined1_tree
done

```

For U. pustulata:

```{bash}

# for circadian genes

for gene in acr-2 adv-1 ccg-1 ccg-7 ccg-9 Chd3 chol-1 ck-1a cka-1 cka-2 ckb-1-1 ckb-1-2 con-6 con-10 crf6-1 crf10-1 csn-1 csn-2 csn-4 csn-5-1 csn-5-2 csn-6 csn-7a csp-2 cys-9 ecm33 frh frh-2 frq fwd-1 isw-1 nit-2 pph-3-1 pph-3-2 pph-3-3 pph-3-4 prd-4-1 prd-4-2 ras-1-1 ras-1-2 ras-1-3 ras-2 rent-1 rgb-1 rrg-1 sod-1 vvd-1 vvd-2 wc-1 wc-2
  do mafft Upust_${gene}_combined1.fasta > Upust_${gene}_combined1.msa
  raxmlHPC -m GTRGAMMA  -s Upust_${gene}_combined1.msa -n Upust_${gene}_combined1_tree
done

# for temperature genes

for gene in hog1 cut-1b xyl1 hsp90 GNAQ hsp30 SSB1 fbp-1 hsp70b treB bipA acu-6 cr-1 hsf-2 hsp70 cat-4 hsp98 gcy-3 mrh-4 Pfd get3 gcy-1 Q1K6J1 LSP1 hsf-1 GNA3 GlyP cat-4b dak-1 AGA1 suc katG hsp88 cut-1 GNA2 hsp90a gsy-1 
  do mafft Upust_${gene}_combined1.fasta > Upust_${gene}_combined1.msa
  raxmlHPC -m GTRGAMMA  -s Upust_${gene}_combined1.msa -n Upust_${gene}_combined1_tree
done

```


## 6. Run HyPhy analyses (using conda)

Generate the tree by loading aligned (.msa) files into Geneious and running RAxML (100 bootstraps) before downloading and saving files as .nex trees.

Then, use the CleanStopCodons.bf batch file that is a part of the HyPhy pipeline to convert all stop codons in the alignment into ambiguous loci:

```{bash}

hyphy ~/Dropbox\ \(Senckenberg\)/Valim/04_Bioinformatic_work/05_Population_genetics/10_hyphy_analysis/hyphy-develop/res/TemplateBatchFiles/CleanStopCodons.bf

```

Just follow the prompts, inputting your file name. Use the basic Universal codon code and option 1 (keep all sequences and sites). 
Save the file as an *.msa file with the same name as the input file, plus "_cleaned" or something similar. 

Next, perform the test itself by using hyphy busted in this case (or whatever other test may be of interest, check http://hyphy.org/getting-started/#typical-uses-of-hyphy). 
As input files, use the _cleaned.msa file generated by the CleanStopCodons.bf batch file above, and the .nex file generated by RAxML in Geneious. 
May want to test MEME as it is the standard test of positive selection at specific loci, as well as BUSTED

For U. phaea:

```{bash}

hyphy busted --alignment Upust_ccg-9_combined1_cleaned.msa --tree Upust_ccg-9_combined1_RAxML_Tree.nex

# for circadian genes:

for gene in ACR2 ccg-1 ccg-9 CDC55 CHO2 CKB1 CKB2 con-6 con-10 csn-1 csn-2 csn-4 csn-6 csn-7a csp-2 DUN1 ecm33 FUN30 fwd-1 GAT1 GDP1-2 GLC7 hrp3 HRR25-1 HRR25-3 ISW2 KIN28 MgSsk1 MTR4 OPI1 PPG1 PPH3 PPT1 pzh1 ras1 ras2 RPN11 RRI1 RSR1 SEC4 SKI2-2 SNF2 SOD1 SSN3 TRR1 UME6 upf-1 VAS1 wc-1 wc-2 YCK2
  do echo ${gene}
  hyphy busted --alignment Uph_${gene}_combined1_cleaned.msa --tree Uph_${gene}_combined1_RAxML_Tree.nex | tail -n 1
done

# for temperatuere genes:

for gene in CCP1 DED1 hsp30 KIN28 LSP1 PRP28 PYC1 cut-1b cat-4 GET3 MKC1 NTH1_1 PRF1 DBP3 LSP1b cr-1 KAR2 PAN2 BUR1 gcy-3 HSP82 PHO85 SSB1 cat1b xyl1 DBP4 DBP10 wos2 gcy-1a HSP78_1 PRP5 CAT1 DBP2 SFL1 CTA8 gcy-1b SSA2 FBP1 DES1 SSC1 GSY1 LSK1 MRH4 dak1 HSP98 cut-1a GPA1 GPA1_2 GPH1 SSE1 NTH1_1
  do echo ${gene}
  hyphy busted --alignment Uph_${gene}_combined1_cleaned.msa --tree Uph_${gene}_combined1_cleaned_RAxML_Tree.nex | tail -n 1
done


```

For U. pustulata:

```{bash}

hyphy busted --alignment Upust_ccg-9_combined1_cleaned.msa --tree Upust_ccg-9_combined1_RAxML_Tree.nex

# for circadian genes

for gene in acr-2 adv-1 ccg-1 ccg-7 ccg-9 Chd3 chol-1 ck-1a cka-1 cka-2 ckb-1-1 ckb-1-2 con-6 con-10 crf6-1 crf10-1 csn-1 csn-2 csn-4 csn-5-1 csn-5-2 csn-6 csn-7a csp-2 cys-9 ecm33 frh frh-2 frq fwd-1 isw-1 nit-2 pph-3-1 pph-3-2 pph-3-3 pph-3-4 prd-4-1 prd-4-2 ras-1-1 ras-1-2 ras-1-3 ras-2 rent-1 rgb-1 rrg-1 sod-1 vvd-1 vvd-2 wc-1 wc-2
  do echo ${gene}
  hyphy busted --alignment Upust_${gene}_combined1_cleaned.msa --tree Upust_${gene}_combined1_RAxML_Tree.nex | tail -n 1
done

# for temperature genes

for gene in hog1 cut-1b xyl1 hsp90 GNAQ hsp30 SSB1 fbp-1 hsp70b treB bipA acu-6 cr-1 hsf-2 hsp70 cat-4 hsp98 gcy-3 mrh-4 Pfd get3 gcy-1 Q1K6J1 LSP1 hsf-1 GNA3 GlyP cat-4b dak-1 AGA1 suc katG hsp88 cut-1 GNA2 hsp90a gsy-1 
  do echo ${gene}
  hyphy busted --alignment Upust_${gene}_combined1_cleaned.msa --tree Upust_${gene}_combined1_cleaned_RAxML_Tree.nex | tail -n 1
done



```




