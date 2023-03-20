---
title: "Part 2: diversity analysis with popoolation1 and vegan"
author: "Henrique Valim"
date: "2022-12-20"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```


# 1. Calculating nucleotide diversity using popoolation1

First extract just the chromosome, start and end columns (without the header) of the homologs of N. crassa circadian clock and temperature-associated genes.

Save this as a tab-delimited file, and rename the extension as *.bed and then proceed.

## a. extract full gene list (scaffold, start and stop locations in BED file) using bedtools

```{bash bed to bam}
### for U. pustulata circadian loci:
for i in 1 2 3 4 5 6; do bedtools intersect -a /phylodata/fdalgrande/TRIMMED_PoolSeq_3gradients/IT_Pool${i}_aligned_reads.sort.rmd.q20.bam -b Upust_circadian_GO_labeled.bed > Upust_IT${i}_circadian.bam; done

for i in 1 2 3 4 5 6; do bedtools intersect -a /phylodata/fdalgrande/TRIMMED_PoolSeq_3gradients/ESii_Pool${i}_aligned_reads.sort.rmd.q20.bam -b Upust_circadian_GO_labeled.bed > Upust_ESii${i}_circadian.bam; done

for i in 1 2 3; do bedtools intersect -a /phylodata/fdalgrande/TRIMMED_PoolSeq_3gradients/ESi_Pool${i}_aligned_reads.sort.rmd.q20.bam -b Upust_circadian_GO_labeled.bed > Upust_ESi${i}_circadian.bam; done

### for U. phaea circadian loci:
for i in 16 17 18 19; do bedtools intersect -a ../SierraNevada/Uph${i}_aligned_reads.sort.rmd.q20.bam -b Uph_circadian_GO_labeled.bed > Uph${i}_circadian.bam; done

for i in 22 23 24 25 26 27 28; do bedtools intersect -a ../MtJacinto/Uph${i}_aligned_reads.sort.rmd.q20.bam -b Uph_circadian_GO_labeled.bed > Uph${i}_circadian.bam; done

### For U. hispanica circadian loci:
for i in 1 2 3 4 5 6; do bedtools intersect -a H${i}_aligned_reads.sort.rmd.q20.bam -b Uhis_circadian_GO_labeled.bed > Uhis${i}_circadian.bam; done

#######

### for U. pustulata temperature loci:
for i in 1 2 3 4 5 6; do bedtools intersect -a /phylodata/fdalgrande/TRIMMED_PoolSeq_3gradients/IT_Pool${i}_aligned_reads.sort.rmd.q20.bam -b Upust_temperature_GO_labeled.bed > Upust_IT${i}_temperature.bam; done

for i in 1 2 3 4 5 6; do bedtools intersect -a /phylodata/fdalgrande/TRIMMED_PoolSeq_3gradients/ESii_Pool${i}_aligned_reads.sort.rmd.q20.bam -b Upust_temperature_GO_labeled.bed > Upust_ESii${i}_temperature.bam; done

for i in 1 2 3; do bedtools intersect -a /phylodata/fdalgrande/TRIMMED_PoolSeq_3gradients/ESi_Pool${i}_aligned_reads.sort.rmd.q20.bam -b Upust_temperature_GO_labeled.bed > Upust_ESi${i}_temperature.bam; done

### for U. phaea temperature loci:
for i in 16 17 18 19; do bedtools intersect -a ../SierraNevada/Uph${i}_aligned_reads.sort.rmd.q20.bam -b Uph_temperature_GO_labeled.bed > Uph${i}_temperature.bam; done

for i in 22 23 24 25 26 27 28; do bedtools intersect -a ../MtJacinto/Uph${i}_aligned_reads.sort.rmd.q20.bam -b Uph_temperature_GO_labeled.bed > Uph${i}_temperature.bam; done
```

you can compare the BAM files' number of mapped reads to check if it worked:

```{bash bam check}
samtools view -c -F 260 ../SierraNevada/Uph16_aligned_reads.sort.rmd.q20.bam 
# 16218568
samtools view -c -F 260 Uph16_circadian.bam 
# 65380
```

## b. next, mpileup files must be prepared for each sample individually

```{bash bam to mpileup}
### for U. pustulata circadian loci:
for pop in IT1 IT2 IT3 IT4 IT5 IT6 ESii1 ESii2 ESii3 ESii4 ESii5 ESii6 ESi1 ESi2 ESi3; do samtools mpileup -B -Q 0 -f Lpus_4dec_AN3.masked.fasta Upust_${pop}_circadian.bam > Upust_${pop}_circadian.mpileup; done

### for U. phaea circadian loci:
for i in 16 17 18 19 22 23 24 25 26 27 28; do samtools mpileup -B -Q 0 -f ../Umbilicaria_phaea_Uphaea_ref_cold.scaffolds.fa Uph${i}_circadian.bam > Uph${i}_circadian.mpileup; done

### For U. hispanica circadian loci:
for i in 1 2 3 4 5 6; do samtools mpileup -B -Q 0 -f Lasallia_hispanica_U_hispanica_TBG_2337.scaffolds.fa Uhis${i}_circadian.bam > Uhis${i}_circadian.mpileup; done

##########

### for U. pustulata temperature loci:
for pop in IT1 IT2 IT3 IT4 IT5 IT6 ESii1 ESii2 ESii3 ESii4 ESii5 ESii6 ESi1 ESi2 ESi3; do samtools mpileup -B -Q 0 -f Lpus_4dec_AN3.masked.fasta Upust_${pop}_temperature.bam > Upust_${pop}_temperature.mpileup; done

### for U. phaea temperature loci:
for i in 16 17 18 19 22 23 24 25 26 27 28; do samtools mpileup -B -Q 0 -f ../Umbilicaria_phaea_Uphaea_ref_cold.scaffolds.fa Uph${i}_temperature.bam > Uph${i}_temperature.mpileup; done
```


## c. get nucleotide diversity directly from mpileup files

```{bash mpileup to pi}
### for U. pustulata circadian loci:
for pop in IT1 IT2 IT3 IT4 IT5 IT6 ESii1 ESii2 ESii3 ESii4 ESii5 ESii6 ESi1 ESi2 ESi3; do perl /home/hvalim/tools/popoolation_1.2.2/Variance-sliding.pl --measure pi --input Upust_${pop}_circadian.mpileup --output Upust_${pop}_circadian.pi --min-count 2 --min-coverage 10 --max-coverage 100 --window-size 1 --step-size 1 --pool-size 100 --fastq-type sanger; done

### for U. phaea circadian loci:
for i in 16 17 18 19 22 23 24 25 26 27 28; do perl /home/hvalim/tools/popoolation_1.2.2/Variance-sliding.pl  --measure pi --input Uph${i}_circadian.mpileup --output Uph${i}_circadian.pi --min-count 2 --min-coverage 10 --max-coverage 100 --window-size 1 --step-size 1 --pool-size 100 --fastq-type sanger; done

############

### for U. pustulata temperature loci:
for pop in IT1 IT2 IT3 IT4 IT5 IT6 ESii1 ESii2 ESii3 ESii4 ESii5 ESii6 ESi1 ESi2 ESi3; do perl /home/hvalim/tools/popoolation_1.2.2/Variance-sliding.pl --measure pi --input Upust_${pop}_temperature.mpileup --output Upust_${pop}_temperature.pi --min-count 2 --min-coverage 10 --max-coverage 100 --window-size 1 --step-size 1 --pool-size 100 --fastq-type sanger; done

### for U. phaea temperature loci:
for i in 16 17 18 19 22 23 24 25 26 27 28; do perl /home/hvalim/tools/popoolation_1.2.2/Variance-sliding.pl  --measure pi --input Uph${i}_temperature.mpileup --output Uph${i}_temperature.pi --min-count 2 --min-coverage 10 --max-coverage 100 --window-size 1 --step-size 1 --pool-size 100 --fastq-type sanger; done
```

# 2. Combining PI files and initial analysis

Next, we move over to R and import the PI files and process them.

## a. U. pustulata

```{r circadian loci PI PCA, U. pustulata }

library(tidyr)
library(ape)
library(vegan)
library(dplyr)

# import each pi file 
pi_1 <- read.delim("Upust_IT1_circadian.pi", header = F)
pi_2 <- read.delim("Upust_IT2_circadian.pi", header = F)
pi_3 <- read.delim("Upust_IT3_circadian.pi", header = F)
pi_4 <- read.delim("Upust_IT4_circadian.pi", header = F)
pi_5 <- read.delim("Upust_IT5_circadian.pi", header = F)
pi_6 <- read.delim("Upust_IT6_circadian.pi", header = F)

pi_7 <- read.delim("Upust_ESi1_circadian.pi", header = F)
pi_8 <- read.delim("Upust_ESi2_circadian.pi", header = F)
pi_9 <- read.delim("Upust_ESi3_circadian.pi", header = F)

pi_10 <- read.delim("Upust_ESii1_circadian.pi", header = F)
pi_11 <- read.delim("Upust_ESii2_circadian.pi", header = F)
pi_12 <- read.delim("Upust_ESii3_circadian.pi", header = F)
pi_13 <- read.delim("Upust_ESii4_circadian.pi", header = F)
pi_14 <- read.delim("Upust_ESii5_circadian.pi", header = F)
pi_15 <- read.delim("Upust_ESii6_circadian.pi", header = F)


# add column names
colnames(pi_1) <- c("Chromosome", "window", "num.snps", "frac", "pi_1")
colnames(pi_2) <- c("Chromosome", "window", "num.snps", "frac", "pi_2")
colnames(pi_3) <- c("Chromosome", "window", "num.snps", "frac", "pi_3")
colnames(pi_4) <- c("Chromosome", "window", "num.snps", "frac", "pi_4")
colnames(pi_5) <- c("Chromosome", "window", "num.snps", "frac", "pi_5")
colnames(pi_6) <- c("Chromosome", "window", "num.snps", "frac", "pi_6")
colnames(pi_7) <- c("Chromosome", "window", "num.snps", "frac", "pi_7")
colnames(pi_8) <- c("Chromosome", "window", "num.snps", "frac", "pi_8")
colnames(pi_9) <- c("Chromosome", "window", "num.snps", "frac", "pi_9")
colnames(pi_10) <- c("Chromosome", "window", "num.snps", "frac", "pi_10")
colnames(pi_11) <- c("Chromosome", "window", "num.snps", "frac", "pi_11")
colnames(pi_12) <- c("Chromosome", "window", "num.snps", "frac", "pi_12")
colnames(pi_13) <- c("Chromosome", "window", "num.snps", "frac", "pi_13")
colnames(pi_14) <- c("Chromosome", "window", "num.snps", "frac", "pi_14")
colnames(pi_15) <- c("Chromosome", "window", "num.snps", "frac", "pi_15")

# make ID column
pi_1$ID <- paste(pi_1$Chromosome, pi_1$window, sep = '_')
pi_2$ID <- paste(pi_2$Chromosome, pi_2$window, sep = '_')
pi_3$ID <- paste(pi_3$Chromosome, pi_3$window, sep = '_')
pi_4$ID <- paste(pi_4$Chromosome, pi_4$window, sep = '_')
pi_5$ID <- paste(pi_5$Chromosome, pi_5$window, sep = '_')
pi_6$ID <- paste(pi_6$Chromosome, pi_6$window, sep = '_')
pi_7$ID <- paste(pi_7$Chromosome, pi_7$window, sep = '_')
pi_8$ID <- paste(pi_8$Chromosome, pi_8$window, sep = '_')
pi_9$ID <- paste(pi_9$Chromosome, pi_9$window, sep = '_')
pi_10$ID <- paste(pi_10$Chromosome, pi_10$window, sep = '_')
pi_11$ID <- paste(pi_11$Chromosome, pi_11$window, sep = '_')
pi_12$ID <- paste(pi_12$Chromosome, pi_12$window, sep = '_')
pi_13$ID <- paste(pi_13$Chromosome, pi_13$window, sep = '_')
pi_14$ID <- paste(pi_14$Chromosome, pi_14$window, sep = '_')
pi_15$ID <- paste(pi_15$Chromosome, pi_15$window, sep = '_')

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
pi_1 <- pi_1[,c(6,5)]
pi_2 <- pi_2[,c(6,5)]
pi_3 <- pi_3[,c(6,5)]
pi_4 <- pi_4[,c(6,5)]
pi_5 <- pi_5[,c(6,5)]
pi_6 <- pi_6[,c(6,5)]
pi_7 <- pi_7[,c(6,5)]
pi_8 <- pi_8[,c(6,5)]
pi_9 <- pi_9[,c(6,5)]
pi_10 <- pi_10[,c(6,5)]
pi_11 <- pi_11[,c(6,5)]
pi_12 <- pi_12[,c(6,5)]
pi_13 <- pi_13[,c(6,5)]
pi_14 <- pi_14[,c(6,5)]
pi_15 <- pi_15[,c(6,5)]

#merge all pi results into one file
pi_all1 <- merge(pi_1, pi_2, by="ID")
pi_all2 <- merge(pi_all1, pi_3, by="ID")
pi_all3 <- merge(pi_all2, pi_4, by="ID")
pi_all4 <- merge(pi_all3, pi_5, by="ID")
pi_all5 <- merge(pi_all4, pi_6, by="ID")
pi_all6 <- merge(pi_all5, pi_7, by="ID")
pi_all7 <- merge(pi_all6, pi_8, by="ID")
pi_all8 <- merge(pi_all7, pi_9, by="ID")
pi_all9 <- merge(pi_all8, pi_10, by="ID")
pi_all10 <- merge(pi_all9, pi_11, by="ID")
pi_all11 <- merge(pi_all10, pi_12, by="ID")
pi_all12 <- merge(pi_all11, pi_13, by="ID")
pi_all13 <- merge(pi_all12, pi_14, by="ID")
pi_all <- merge(pi_all13, pi_15, by="ID")

pi_all$pi_1 <- as.numeric(pi_all$pi_1)
pi_all$pi_2 <- as.numeric(pi_all$pi_2)
pi_all$pi_3 <- as.numeric(pi_all$pi_3)
pi_all$pi_4 <- as.numeric(pi_all$pi_4)
pi_all$pi_5 <- as.numeric(pi_all$pi_5)
pi_all$pi_6 <- as.numeric(pi_all$pi_6)
pi_all$pi_7 <- as.numeric(pi_all$pi_7)
pi_all$pi_8 <- as.numeric(pi_all$pi_8)
pi_all$pi_9 <- as.numeric(pi_all$pi_9)
pi_all$pi_10 <- as.numeric(pi_all$pi_10)
pi_all$pi_11 <- as.numeric(pi_all$pi_11)
pi_all$pi_12 <- as.numeric(pi_all$pi_12)
pi_all$pi_13 <- as.numeric(pi_all$pi_13)
pi_all$pi_14 <- as.numeric(pi_all$pi_14)
pi_all$pi_15 <- as.numeric(pi_all$pi_15)

dds.pcoa=pcoa(vegdist(t((pi_all[,2:16])), method="euclidean")/100)
scores=dds.pcoa$vectors
percent <- dds.pcoa$values$Eigenvalues
percent / sum(percent) #percent for each axes; change in caption below

labels1 <- c("IT-1","IT-2","IT-3","IT-4","IT-5","IT-6","ESi-1","ESi-2","ESi-3","ESii-1","ESii-2","ESii-3","ESii-4","ESii-5","ESii-6")

pdf(file.path("PCA_circadian_nucl_div_Upust_15pops_PC1_PC2_NEW.pdf"))
plot(scores[,1], scores[,2],
     col=c("red","red","red","red","purple","blue", "red","blue","blue", "red","red","purple","blue","blue","blue"), 
     pch = c(19,19,19,19,19,19,19,19,19,19,19,19,19,19,19),
     xlab = "PC1 (34.26%)", ylab = "PC2 (28.40%)",main="Nucleotide diversity, circadian loci",sub="Umbilicaria pustulata, Sardinian and Spanish gradients",)
text(scores[,1], scores[,2], labels=labels1, cex= 0.7, pos=4)
dev.off()

pdf(file.path("PCA_circadian_nucl_div_Upust_15pops_PC2_PC3_NEW.pdf"))
plot(scores[,2], scores[,3],
     col=c("red","red","red","red","purple","blue", "red","blue","blue", "red","red","purple","blue","blue","blue"), 
     pch = c(19,19,19,19,19,19,19,19,19,19,19,19,19,19,19),
     xlab = "PC2 (28.40%)", ylab = "PC3 (9.28%)",main="Nucleotide diversity, circadian loci",sub="Umbilicaria pustulata, Sardinian and Spanish gradients",)
text(scores[,2], scores[,3], labels=labels1, cex= 0.7, pos=4)
dev.off()

#test site (gradient) and climate (location on gradient) differences
site_df <- data.frame(
  ID = c("pi_1", "pi_2", "pi_3", "pi_4", "pi_5", "pi_6", "pi_7", "pi_8", "pi_9", "pi_10", "pi_11", "pi_12", "pi_13", "pi_14", "pi_15"),
  line = c("IT","IT","IT","IT","IT","IT","ESi","ESi","ESi","ESii","ESii","ESii","ESii","ESii","ESii")
)

clim_df <- data.frame(
  ID = c("pi_1", "pi_2", "pi_3", "pi_4", "pi_5", "pi_6", "pi_7", "pi_8", "pi_9", "pi_10", "pi_11", "pi_12", "pi_13", "pi_14", "pi_15"),
  line = c("red","red","red","red","purple","blue", "red","blue","blue", "red","red","purple","blue","blue","blue")
)


line_site <- site_df$line
line_clim <- clim_df$line

adonis2_site <- adonis2(t(pi_all[,2:16])~line, data = site_df, permutations = 1000000, method = "manhattan") # 0.1299
adonis2_clim <- adonis2(t(pi_all[,2:16])~line, data = clim_df, permutations = 1000000, method = "manhattan") # 0.006145 **

capture.output(adonis2_site, file = "Upust_adonis2_circadian_site_NEW.txt")
capture.output(adonis2_clim, file = "Upust_adonis2_circadian_clim_NEW.txt")

# climate zone is highly significant, but site is not

write.table(pi_all, file = "Upust_circadian_GO_all.pi", sep = "\t", na = "", quote = F, row.names = F, col.names = T)
```

a
```{r temperature loci PI PCA, U. pustulata }

library(tidyr)
library(ape)
library(vegan)
library(dplyr)

# import each pi file 
pi_1 <- read.delim("Upust_IT1_temperature.pi", header = F)
pi_2 <- read.delim("Upust_IT2_temperature.pi", header = F)
pi_3 <- read.delim("Upust_IT3_temperature.pi", header = F)
pi_4 <- read.delim("Upust_IT4_temperature.pi", header = F)
pi_5 <- read.delim("Upust_IT5_temperature.pi", header = F)
pi_6 <- read.delim("Upust_IT6_temperature.pi", header = F)

pi_7 <- read.delim("Upust_ESi1_temperature.pi", header = F)
pi_8 <- read.delim("Upust_ESi2_temperature.pi", header = F)
pi_9 <- read.delim("Upust_ESi3_temperature.pi", header = F)

pi_10 <- read.delim("Upust_ESii1_temperature.pi", header = F)
pi_11 <- read.delim("Upust_ESii2_temperature.pi", header = F)
pi_12 <- read.delim("Upust_ESii3_temperature.pi", header = F)
pi_13 <- read.delim("Upust_ESii4_temperature.pi", header = F)
pi_14 <- read.delim("Upust_ESii5_temperature.pi", header = F)
pi_15 <- read.delim("Upust_ESii6_temperature.pi", header = F)


# add column names
colnames(pi_1) <- c("Chromosome", "window", "num.snps", "frac", "pi_1")
colnames(pi_2) <- c("Chromosome", "window", "num.snps", "frac", "pi_2")
colnames(pi_3) <- c("Chromosome", "window", "num.snps", "frac", "pi_3")
colnames(pi_4) <- c("Chromosome", "window", "num.snps", "frac", "pi_4")
colnames(pi_5) <- c("Chromosome", "window", "num.snps", "frac", "pi_5")
colnames(pi_6) <- c("Chromosome", "window", "num.snps", "frac", "pi_6")
colnames(pi_7) <- c("Chromosome", "window", "num.snps", "frac", "pi_7")
colnames(pi_8) <- c("Chromosome", "window", "num.snps", "frac", "pi_8")
colnames(pi_9) <- c("Chromosome", "window", "num.snps", "frac", "pi_9")
colnames(pi_10) <- c("Chromosome", "window", "num.snps", "frac", "pi_10")
colnames(pi_11) <- c("Chromosome", "window", "num.snps", "frac", "pi_11")
colnames(pi_12) <- c("Chromosome", "window", "num.snps", "frac", "pi_12")
colnames(pi_13) <- c("Chromosome", "window", "num.snps", "frac", "pi_13")
colnames(pi_14) <- c("Chromosome", "window", "num.snps", "frac", "pi_14")
colnames(pi_15) <- c("Chromosome", "window", "num.snps", "frac", "pi_15")

# make ID column
pi_1$ID <- paste(pi_1$Chromosome, pi_1$window, sep = '_')
pi_2$ID <- paste(pi_2$Chromosome, pi_2$window, sep = '_')
pi_3$ID <- paste(pi_3$Chromosome, pi_3$window, sep = '_')
pi_4$ID <- paste(pi_4$Chromosome, pi_4$window, sep = '_')
pi_5$ID <- paste(pi_5$Chromosome, pi_5$window, sep = '_')
pi_6$ID <- paste(pi_6$Chromosome, pi_6$window, sep = '_')
pi_7$ID <- paste(pi_7$Chromosome, pi_7$window, sep = '_')
pi_8$ID <- paste(pi_8$Chromosome, pi_8$window, sep = '_')
pi_9$ID <- paste(pi_9$Chromosome, pi_9$window, sep = '_')
pi_10$ID <- paste(pi_10$Chromosome, pi_10$window, sep = '_')
pi_11$ID <- paste(pi_11$Chromosome, pi_11$window, sep = '_')
pi_12$ID <- paste(pi_12$Chromosome, pi_12$window, sep = '_')
pi_13$ID <- paste(pi_13$Chromosome, pi_13$window, sep = '_')
pi_14$ID <- paste(pi_14$Chromosome, pi_14$window, sep = '_')
pi_15$ID <- paste(pi_15$Chromosome, pi_15$window, sep = '_')

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
pi_1 <- pi_1[,c(6,5)]
pi_2 <- pi_2[,c(6,5)]
pi_3 <- pi_3[,c(6,5)]
pi_4 <- pi_4[,c(6,5)]
pi_5 <- pi_5[,c(6,5)]
pi_6 <- pi_6[,c(6,5)]
pi_7 <- pi_7[,c(6,5)]
pi_8 <- pi_8[,c(6,5)]
pi_9 <- pi_9[,c(6,5)]
pi_10 <- pi_10[,c(6,5)]
pi_11 <- pi_11[,c(6,5)]
pi_12 <- pi_12[,c(6,5)]
pi_13 <- pi_13[,c(6,5)]
pi_14 <- pi_14[,c(6,5)]
pi_15 <- pi_15[,c(6,5)]

#merge all pi results into one file
pi_all1 <- merge(pi_1, pi_2, by="ID")
pi_all2 <- merge(pi_all1, pi_3, by="ID")
pi_all3 <- merge(pi_all2, pi_4, by="ID")
pi_all4 <- merge(pi_all3, pi_5, by="ID")
pi_all5 <- merge(pi_all4, pi_6, by="ID")
pi_all6 <- merge(pi_all5, pi_7, by="ID")
pi_all7 <- merge(pi_all6, pi_8, by="ID")
pi_all8 <- merge(pi_all7, pi_9, by="ID")
pi_all9 <- merge(pi_all8, pi_10, by="ID")
pi_all10 <- merge(pi_all9, pi_11, by="ID")
pi_all11 <- merge(pi_all10, pi_12, by="ID")
pi_all12 <- merge(pi_all11, pi_13, by="ID")
pi_all13 <- merge(pi_all12, pi_14, by="ID")
pi_all <- merge(pi_all13, pi_15, by="ID")

pi_all$pi_1 <- as.numeric(pi_all$pi_1)
pi_all$pi_2 <- as.numeric(pi_all$pi_2)
pi_all$pi_3 <- as.numeric(pi_all$pi_3)
pi_all$pi_4 <- as.numeric(pi_all$pi_4)
pi_all$pi_5 <- as.numeric(pi_all$pi_5)
pi_all$pi_6 <- as.numeric(pi_all$pi_6)
pi_all$pi_7 <- as.numeric(pi_all$pi_7)
pi_all$pi_8 <- as.numeric(pi_all$pi_8)
pi_all$pi_9 <- as.numeric(pi_all$pi_9)
pi_all$pi_10 <- as.numeric(pi_all$pi_10)
pi_all$pi_11 <- as.numeric(pi_all$pi_11)
pi_all$pi_12 <- as.numeric(pi_all$pi_12)
pi_all$pi_13 <- as.numeric(pi_all$pi_13)
pi_all$pi_14 <- as.numeric(pi_all$pi_14)
pi_all$pi_15 <- as.numeric(pi_all$pi_15)

dds.pcoa=pcoa(vegdist(t((pi_all[,2:16])), method="euclidean")/100)
scores=dds.pcoa$vectors
percent <- dds.pcoa$values$Eigenvalues
percent / sum(percent) #percent for each axes; change in caption below

labels1 <- c("IT-1","IT-2","IT-3","IT-4","IT-5","IT-6","ESi-1","ESi-2","ESi-3","ESii-1","ESii-2","ESii-3","ESii-4","ESii-5","ESii-6")

pdf(file.path("PCA_temperature_nucl_div_Upust_15pops_PC1_PC2_NEW.pdf"))
plot(scores[,1], scores[,2],
     col=c("red","red","red","red","purple","blue", "red","blue","blue", "red","red","purple","blue","blue","blue"), 
     pch = c(19,19,19,19,19,19,19,19,19,19,19,19,19,19,19),
     xlab = "PC1 (30.59%)", ylab = "PC2 (20.19%)",main="Nucleotide diversity, temperature loci",sub="Umbilicaria pustulata, Sardinian and Spanish gradients",)
text(scores[,1], scores[,2], labels=labels1, cex= 0.7, pos=4)
dev.off()

pdf(file.path("PCA_temperature_nucl_div_Upust_15pops_PC2_PC3_NEW.pdf"))
plot(scores[,2], scores[,3],
     col=c("red","red","red","red","purple","blue", "red","blue","blue", "red","red","purple","blue","blue","blue"), 
     pch = c(19,19,19,19,19,19,19,19,19,19,19,19,19,19,19),
     xlab = "PC2 (20.19%)", ylab = "PC3 (9.22%)",main="Nucleotide diversity, temperature loci",sub="Umbilicaria pustulata, Sardinian and Spanish gradients",)
text(scores[,2], scores[,3], labels=labels1, cex= 0.7, pos=4)
dev.off()

#test site (gradient) and climate (location on gradient) differences
site_df <- data.frame(
  ID = c("pi_1", "pi_2", "pi_3", "pi_4", "pi_5", "pi_6", "pi_7", "pi_8", "pi_9", "pi_10", "pi_11", "pi_12", "pi_13", "pi_14", "pi_15"),
  line = c("IT","IT","IT","IT","IT","IT","ESi","ESi","ESi","ESii","ESii","ESii","ESii","ESii","ESii")
)

clim_df <- data.frame(
  ID = c("pi_1", "pi_2", "pi_3", "pi_4", "pi_5", "pi_6", "pi_7", "pi_8", "pi_9", "pi_10", "pi_11", "pi_12", "pi_13", "pi_14", "pi_15"),
  line = c("red","red","red","red","purple","blue", "red","blue","blue", "red","red","purple","blue","blue","blue")
)


line_site <- site_df$line
line_clim <- clim_df$line

adonis2_site <- adonis2(t(pi_all[,2:16])~line, data = site_df, permutations = 1000000, method = "manhattan") # 
adonis2_clim <- adonis2(t(pi_all[,2:16])~line, data = clim_df, permutations = 1000000, method = "manhattan") # 

capture.output(adonis2_site, file = "Upust_adonis2_temperature_site_NEW.txt")
capture.output(adonis2_clim, file = "Upust_adonis2_temperature_clim_NEW.txt")

write.table(pi_all, file = "Upust_temperature_GO_all.pi", sep = "\t", na = "", quote = F, row.names = F, col.names = T)

```


## b. U. phaea

```{r circadian loci PI PCA, U. phaea }

library(tidyr)
library(ape)
library(vegan)
library(dplyr)

#import each pi file 

# circadian
pi_16 <- read.delim("Uph16_circadian.pi", header = F)
pi_17 <- read.delim("Uph17_circadian.pi", header = F)
pi_18 <- read.delim("Uph18_circadian.pi", header = F)
pi_19 <- read.delim("Uph19_circadian.pi", header = F)

pi_22 <- read.delim("Uph22_circadian.pi", header = F)
pi_23 <- read.delim("Uph23_circadian.pi", header = F)
pi_24 <- read.delim("Uph24_circadian.pi", header = F)
pi_25 <- read.delim("Uph25_circadian.pi", header = F)
pi_26 <- read.delim("Uph26_circadian.pi", header = F)
pi_27 <- read.delim("Uph27_circadian.pi", header = F)
pi_28 <- read.delim("Uph28_circadian.pi", header = F)


# add column names
colnames(pi_16) <- c("Chromosome", "window", "num.snps", "frac", "pi_16")
colnames(pi_17) <- c("Chromosome", "window", "num.snps", "frac", "pi_17")
colnames(pi_18) <- c("Chromosome", "window", "num.snps", "frac", "pi_18")
colnames(pi_19) <- c("Chromosome", "window", "num.snps", "frac", "pi_19")

colnames(pi_22) <- c("Chromosome", "window", "num.snps", "frac", "pi_22")
colnames(pi_23) <- c("Chromosome", "window", "num.snps", "frac", "pi_23")
colnames(pi_24) <- c("Chromosome", "window", "num.snps", "frac", "pi_24")
colnames(pi_25) <- c("Chromosome", "window", "num.snps", "frac", "pi_25")
colnames(pi_26) <- c("Chromosome", "window", "num.snps", "frac", "pi_26")
colnames(pi_27) <- c("Chromosome", "window", "num.snps", "frac", "pi_27")
colnames(pi_28) <- c("Chromosome", "window", "num.snps", "frac", "pi_28")

# make ID column
pi_16$ID <- paste(pi_16$Chromosome, pi_16$window, sep = '_')
pi_17$ID <- paste(pi_17$Chromosome, pi_17$window, sep = '_')
pi_18$ID <- paste(pi_18$Chromosome, pi_18$window, sep = '_')
pi_19$ID <- paste(pi_19$Chromosome, pi_19$window, sep = '_')

pi_22$ID <- paste(pi_22$Chromosome, pi_22$window, sep = '_')
pi_23$ID <- paste(pi_23$Chromosome, pi_23$window, sep = '_')
pi_24$ID <- paste(pi_24$Chromosome, pi_24$window, sep = '_')
pi_25$ID <- paste(pi_25$Chromosome, pi_25$window, sep = '_')
pi_26$ID <- paste(pi_26$Chromosome, pi_26$window, sep = '_')
pi_27$ID <- paste(pi_27$Chromosome, pi_27$window, sep = '_')
pi_28$ID <- paste(pi_28$Chromosome, pi_28$window, sep = '_')

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
pi_16 <- pi_16[,c(6,5)]
pi_17 <- pi_17[,c(6,5)]
pi_18 <- pi_18[,c(6,5)]
pi_19 <- pi_19[,c(6,5)]

pi_22 <- pi_22[,c(6,5)]
pi_23 <- pi_23[,c(6,5)]
pi_24 <- pi_24[,c(6,5)]
pi_25 <- pi_25[,c(6,5)]
pi_26 <- pi_26[,c(6,5)]
pi_27 <- pi_27[,c(6,5)]
pi_28 <- pi_28[,c(6,5)]

#merge all pi results into one file
pi_all1 <- merge(pi_16, pi_17, by="ID")
pi_all2 <- merge(pi_all1, pi_18, by="ID")
pi_all3 <- merge(pi_all2, pi_19, by="ID")
pi_all4 <- merge(pi_all3, pi_22, by="ID")
pi_all5 <- merge(pi_all4, pi_23, by="ID")
pi_all6 <- merge(pi_all5, pi_24, by="ID")
pi_all7 <- merge(pi_all6, pi_25, by="ID")
pi_all8 <- merge(pi_all7, pi_26, by="ID")
pi_all9 <- merge(pi_all8, pi_27, by="ID")
pi_all <- merge(pi_all9, pi_28, by="ID")

pi_all$pi_16 <- as.numeric(pi_all$pi_16)
pi_all$pi_17 <- as.numeric(pi_all$pi_17)
pi_all$pi_18 <- as.numeric(pi_all$pi_18)
pi_all$pi_19 <- as.numeric(pi_all$pi_19)

pi_all$pi_22 <- as.numeric(pi_all$pi_22)
pi_all$pi_23 <- as.numeric(pi_all$pi_23)
pi_all$pi_24 <- as.numeric(pi_all$pi_24)
pi_all$pi_25 <- as.numeric(pi_all$pi_25)
pi_all$pi_26 <- as.numeric(pi_all$pi_26)
pi_all$pi_27 <- as.numeric(pi_all$pi_27)
pi_all$pi_28 <- as.numeric(pi_all$pi_28)

dds.pcoa=pcoa(vegdist(t((pi_all[,2:12])), method="euclidean")/100)
scores=dds.pcoa$vectors
percent <- dds.pcoa$values$Eigenvalues
percent / sum(percent) #percent for each axes; change in caption below

labels1 <- c("SN-1", "SN-2", "SN-3", "SN-4","MJ-1","MJ-2","MJ-3","MJ-4","MJ-5","MJ-6","MJ-7")

pdf(file.path("PCA_circadian_nucl_div_Uph_11pops_PC1_PC2_NEW.pdf"))
plot(scores[,1], scores[,2],
     col=c("red","red","purple","blue", "red","red","red","purple","blue","blue","blue"),
     pch = c(19,19,19,19,19,19,19,19,19,19,19),
     xlab = "PC1 (23.89%)", ylab = "PC2 (21.73%)",main="Nucleotide diversity, temperature loci",sub="Umbilicaria phaea, Sierra Nevada and Mt Jacinto Gradients",)
text(scores[,1], scores[,2], labels=labels1, cex= 0.7, pos=4)
dev.off()

pdf(file.path("PCA_circadian_nucl_div_Uph_11pops_PC2_PC3_NEW.pdf"))
plot(scores[,2], scores[,3],
     col=c("red","red","purple","blue", "red","red","red","purple","blue","blue","blue"),
     pch = c(19,19,19,19,19,19,19,19,19,19,19),
     xlab = "PC2 (21.73%)", ylab = "PC3 (13.86%)",main="Nucleotide diversity, temperature loci",sub="Umbilicaria phaea, Sierra Nevada and Mt Jacinto Gradients",)
text(scores[,2], scores[,3], labels=labels1, cex= 0.7, pos=4)
dev.off()

#test site (gradient) and climate (location on gradient) differences

site_df <- data.frame(
  ID = c("pi_16", "pi_17", "pi_18", "pi_19", "pi_22", "pi_23", "pi_24", "pi_25", "pi_26", "pi_27", "pi_28"),
  line = c("SN", "SN", "SN", "SN","MJ", "MJ", "MJ", "MJ","MJ", "MJ", "MJ")
)

clim_df <- data.frame(
  ID = c("pi_16", "pi_17", "pi_18", "pi_19", "pi_22", "pi_23", "pi_24", "pi_25", "pi_26", "pi_27", "pi_28"),
  line = c("red","red","purple","blue", "red","red","red","purple","blue","blue","blue")
)

line_site <- site_df$line
line_clim <- clim_df$line

adonis2_site <- adonis2(t(pi_all[,2:12])~line, data = site_df, permutations = 1000000, method = "manhattan") # 0.1212
adonis2_clim <- adonis2(t(pi_all[,2:12])~line, data = clim_df, permutations = 1000000, method = "manhattan") # 0.000155 ***

capture.output(adonis2_site, file = "Uph_adonis2_circadian_site_NEW.txt")
capture.output(adonis2_clim, file = "Uph_adonis2_circadian_clim_NEW.txt")
# climate zone is highly significant, but site is not

write.table(pi_all, file = "Uph_circadian_GO_all.pi", sep = "\t", na = "", quote = F, row.names = F, col.names = T)

```


```{r temperature loci PI PCA, U. phaea }

library(tidyr)
library(ape)
library(vegan)
library(dplyr)

#import each pi file 

# temperature
pi_16 <- read.delim("Uph16_temperature.pi", header = F)
pi_17 <- read.delim("Uph17_temperature.pi", header = F)
pi_18 <- read.delim("Uph18_temperature.pi", header = F)
pi_19 <- read.delim("Uph19_temperature.pi", header = F)

pi_22 <- read.delim("Uph22_temperature.pi", header = F)
pi_23 <- read.delim("Uph23_temperature.pi", header = F)
pi_24 <- read.delim("Uph24_temperature.pi", header = F)
pi_25 <- read.delim("Uph25_temperature.pi", header = F)
pi_26 <- read.delim("Uph26_temperature.pi", header = F)
pi_27 <- read.delim("Uph27_temperature.pi", header = F)
pi_28 <- read.delim("Uph28_temperature.pi", header = F)

# add column names
colnames(pi_16) <- c("Chromosome", "window", "num.snps", "frac", "pi_16")
colnames(pi_17) <- c("Chromosome", "window", "num.snps", "frac", "pi_17")
colnames(pi_18) <- c("Chromosome", "window", "num.snps", "frac", "pi_18")
colnames(pi_19) <- c("Chromosome", "window", "num.snps", "frac", "pi_19")

colnames(pi_22) <- c("Chromosome", "window", "num.snps", "frac", "pi_22")
colnames(pi_23) <- c("Chromosome", "window", "num.snps", "frac", "pi_23")
colnames(pi_24) <- c("Chromosome", "window", "num.snps", "frac", "pi_24")
colnames(pi_25) <- c("Chromosome", "window", "num.snps", "frac", "pi_25")
colnames(pi_26) <- c("Chromosome", "window", "num.snps", "frac", "pi_26")
colnames(pi_27) <- c("Chromosome", "window", "num.snps", "frac", "pi_27")
colnames(pi_28) <- c("Chromosome", "window", "num.snps", "frac", "pi_28")

# make ID column
pi_16$ID <- paste(pi_16$Chromosome, pi_16$window, sep = '_')
pi_17$ID <- paste(pi_17$Chromosome, pi_17$window, sep = '_')
pi_18$ID <- paste(pi_18$Chromosome, pi_18$window, sep = '_')
pi_19$ID <- paste(pi_19$Chromosome, pi_19$window, sep = '_')

pi_22$ID <- paste(pi_22$Chromosome, pi_22$window, sep = '_')
pi_23$ID <- paste(pi_23$Chromosome, pi_23$window, sep = '_')
pi_24$ID <- paste(pi_24$Chromosome, pi_24$window, sep = '_')
pi_25$ID <- paste(pi_25$Chromosome, pi_25$window, sep = '_')
pi_26$ID <- paste(pi_26$Chromosome, pi_26$window, sep = '_')
pi_27$ID <- paste(pi_27$Chromosome, pi_27$window, sep = '_')
pi_28$ID <- paste(pi_28$Chromosome, pi_28$window, sep = '_')

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
pi_16 <- pi_16[,c(6,5)]
pi_17 <- pi_17[,c(6,5)]
pi_18 <- pi_18[,c(6,5)]
pi_19 <- pi_19[,c(6,5)]

pi_22 <- pi_22[,c(6,5)]
pi_23 <- pi_23[,c(6,5)]
pi_24 <- pi_24[,c(6,5)]
pi_25 <- pi_25[,c(6,5)]
pi_26 <- pi_26[,c(6,5)]
pi_27 <- pi_27[,c(6,5)]
pi_28 <- pi_28[,c(6,5)]

#merge all pi results into one file
pi_all1 <- merge(pi_16, pi_17, by="ID")
pi_all2 <- merge(pi_all1, pi_18, by="ID")
pi_all3 <- merge(pi_all2, pi_19, by="ID")
pi_all4 <- merge(pi_all3, pi_22, by="ID")
pi_all5 <- merge(pi_all4, pi_23, by="ID")
pi_all6 <- merge(pi_all5, pi_24, by="ID")
pi_all7 <- merge(pi_all6, pi_25, by="ID")
pi_all8 <- merge(pi_all7, pi_26, by="ID")
pi_all9 <- merge(pi_all8, pi_27, by="ID")
pi_all <- merge(pi_all9, pi_28, by="ID")

pi_all$pi_16 <- as.numeric(pi_all$pi_16)
pi_all$pi_17 <- as.numeric(pi_all$pi_17)
pi_all$pi_18 <- as.numeric(pi_all$pi_18)
pi_all$pi_19 <- as.numeric(pi_all$pi_19)

pi_all$pi_22 <- as.numeric(pi_all$pi_22)
pi_all$pi_23 <- as.numeric(pi_all$pi_23)
pi_all$pi_24 <- as.numeric(pi_all$pi_24)
pi_all$pi_25 <- as.numeric(pi_all$pi_25)
pi_all$pi_26 <- as.numeric(pi_all$pi_26)
pi_all$pi_27 <- as.numeric(pi_all$pi_27)
pi_all$pi_28 <- as.numeric(pi_all$pi_28)

dds.pcoa=pcoa(vegdist(t((pi_all[,2:12])), method="euclidean")/100)
scores=dds.pcoa$vectors
percent <- dds.pcoa$values$Eigenvalues
percent / sum(percent) #percent for each axes; change in caption below

labels1 <- c("SN-1", "SN-2", "SN-3", "SN-4","MJ-1","MJ-2","MJ-3","MJ-4","MJ-5","MJ-6","MJ-7")

pdf(file.path("PCA_temperature_nucl_div_Uph_11pops_PC1_PC2_NEW.pdf"))
plot(scores[,1], scores[,2],
     col=c("red","red","purple","blue", "red","red","red","purple","blue","blue","blue"),
     pch = c(19,19,19,19,19,19,19,19,19,19,19),
     xlab = "PC1 (22.53%)", ylab = "PC2 (18.14%)",main="Nucleotide diversity, temperature loci",sub="Umbilicaria phaea, Sierra Nevada and Mt Jacinto Gradients",)
text(scores[,1], scores[,2], labels=labels1, cex= 0.7, pos=4)
dev.off()

pdf(file.path("PCA_temperature_nucl_div_Uph_11pops_PC2_PC3_NEW.pdf"))
plot(scores[,2], scores[,3],
     col=c("red","red","purple","blue", "red","red","red","purple","blue","blue","blue"),
     pch = c(19,19,19,19,19,19,19,19,19,19,19),
     xlab = "PC2 (18.14%)", ylab = "PC3 (13.15%)",main="Nucleotide diversity, temperature loci",sub="Umbilicaria phaea, Sierra Nevada and Mt Jacinto Gradients",)
text(scores[,2], scores[,3], labels=labels1, cex= 0.7, pos=4)
dev.off()

#test site (gradient) and climate (location on gradient) differences

site_df <- data.frame(
  ID = c("pi_16", "pi_17", "pi_18", "pi_19", "pi_22", "pi_23", "pi_24", "pi_25", "pi_26", "pi_27", "pi_28"),
  line = c("SN", "SN", "SN", "SN","MJ", "MJ", "MJ", "MJ","MJ", "MJ", "MJ")
)

clim_df <- data.frame(
  ID = c("pi_16", "pi_17", "pi_18", "pi_19", "pi_22", "pi_23", "pi_24", "pi_25", "pi_26", "pi_27", "pi_28"),
  line = c("red","red","purple","blue", "red","red","red","purple","blue","blue","blue")
)

line_site <- site_df$line
line_clim <- clim_df$line

adonis2_site <- adonis2(t(pi_all[,2:12])~line, data = site_df, permutations = 1000000, method = "manhattan")
adonis2_clim <- adonis2(t(pi_all[,2:12])~line, data = clim_df, permutations = 1000000, method = "manhattan") 

capture.output(adonis2_site, file = "Uph_adonis2_temperature_site_NEW.txt")
capture.output(adonis2_clim, file = "Uph_adonis2_temperature_clim_NEW.txt")
# climate zone is highly significant, but site is not

write.table(pi_all, file = "Uph_temperature_GO_all.pi", sep = "\t", na = "", quote = F, row.names = F, col.names = T)

```

# 3. Analyzing the PI data using NMDS

Finally, we analyze the nucleotide diversity data via NMDS using the vegan package in R.

## a. U. pustulata

```{r circadian loci PI NMDS, U. pustulata}

library(vegan)

# load and prepare files
pi_all <- read.delim("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/Upust_circadian_GO_all.pi", header = T)
pi_all1 <- pi_all[,-1]
rownames(pi_all1) <- pi_all[,1]
pi_all2 <- t(pi_all1)

bioclim <- read.csv("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/bioclim_15pops_raw.csv", header = T)
bioclim1 <- bioclim[,-1] 
rownames(bioclim1) <- bioclim[,1]
bioclim2 <- t(bioclim1)

elevation <- bioclim2[,20]
bio1 <- bioclim2[,1]

treat=c("IT","IT","IT","IT","IT","IT","ESi","ESi","ESi","ESii","ESii","ESii","ESii","ESii","ESii")
labels1 <- c("IT-1","IT-2","IT-3","IT-4","IT-5","IT-6","ESi-1","ESi-2","ESi-3","ESii-1","ESii-2","ESii-3","ESii-4","ESii-5","ESii-6")

# run NMDS
nmds1 = metaMDS(pi_all2, # our data matrix (may need to be transformed with t() in order to get a matrix?)
                k=2) # The number of reduced dimensions

# fit our environmental data and save envfit output
en = envfit(nmds1, bioclim2, permutations = 999, na.rm = TRUE)
en
capture.output(en, file = "~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/Upust_circadian_envfit.txt")

# run ANOSIM (bio1 as categorical variable) and save output
treat=c("red","red","red","red","purple","blue", "red","blue","blue", "red","red","purple","blue","blue","blue")
ano = anosim(pi_all2, treat, distance = "bray", permutations = 9999)
ano
capture.output(ano, file = "~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/Upust_circadian_ANOSIM.txt")


# run Mantel test (bio1 as continuous variable) and save output
dist.pi = vegdist(pi_all2, method = "bray") #abundance data frame - bray curtis dissimilarity
dist.bio1 = dist(bio1, method = "euclidean") #environmental vector - euclidean distance
pi_bio1 = mantel(dist.pi, dist.bio1, method = "spearman", permutations = 9999, na.rm = TRUE) #pi vs bio1 
pi_bio1
capture.output(pi_bio1, file = "~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/Upust_circadian_Mantel.txt")

# check the stress plot
stressplot(nmds1)
plot(nmds1)
plot(nmds1,display=c("sites"))

# plot elevation/bio1 using ordisurf and label the plot with labels1
pdf(file.path("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/NMDS_circadian_nucl_div_Upust_15pops_elevation.pdf"))
ordisurf(nmds1,elevation,main="",col="forestgreen")
orditorp(nmds1,display="sites",col=c("red","red","red","red","purple","blue", "red","blue","blue", "red","red","purple","blue","blue","blue"),
         labels = labels1, air=0.01,cex=1.25)
dev.off()

pdf(file.path("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/NMDS_circadian_nucl_div_Upust_15pops_BIO1.pdf"))
ordisurf(nmds1,bio1,main="",col="black")
#orditorp(nmds1,display="species",col="gray",air=0.01)
orditorp(nmds1,display="sites",col=c("red","red","red","red","purple","blue", "red","blue","blue", "red","red","purple","blue","blue","blue"),
         labels = labels1, air=0.01,cex=1.25)
dev.off()

# save a supplemental figure with all bioclim vectors
pdf(file.path("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/NMDS_circadian_nucl_div_Upust_15pops_w_bioclim_vectors.pdf"))
ordiplot(nmds1,type="n")
orditorp(nmds1,display="sites",col=c("red","red","red","red","purple","blue", "red","blue","blue", "red","red","purple","blue","blue","blue"),
         labels = labels1, air=0.01,cex=1.25)
plot(en, col = "gray40")
dev.off()

```


```{r temperature loci PI NMDS, U. pustulata}

library(vegan)

# load and prepare files
pi_all <- read.delim("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/Upust_temperature_GO_all.pi", header = T)
pi_all1 <- pi_all[,-1]
rownames(pi_all1) <- pi_all[,1]
pi_all2 <- t(pi_all1)

bioclim <- read.csv("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/bioclim_15pops_raw.csv", header = T)
bioclim1 <- bioclim[,-1] 
rownames(bioclim1) <- bioclim[,1]
bioclim2 <- t(bioclim1)

elevation <- bioclim2[,20]
bio1 <- bioclim2[,1]

treat=c("IT","IT","IT","IT","IT","IT","ESi","ESi","ESi","ESii","ESii","ESii","ESii","ESii","ESii")
labels1 <- c("IT-1","IT-2","IT-3","IT-4","IT-5","IT-6","ESi-1","ESi-2","ESi-3","ESii-1","ESii-2","ESii-3","ESii-4","ESii-5","ESii-6")

# run NMDS
nmds1 = metaMDS(pi_all2, # our data matrix (may need to be transformed with t() in order to get a matrix?)
                k=2) # The number of reduced dimensions

# fit our environmental data and save envfit output
en = envfit(nmds1, bioclim2, permutations = 999, na.rm = TRUE)
en
capture.output(en, file = "~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/Upust_temperature_envfit.txt")

# run ANOSIM (bio1 as categorical variable) and save output
treat=c("red","red","red","red","purple","blue", "red","blue","blue", "red","red","purple","blue","blue","blue")
ano = anosim(pi_all2, treat, distance = "bray", permutations = 9999)
ano
capture.output(ano, file = "~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/Upust_temperature_ANOSIM.txt")


# run Mantel test (bio1 as continuous variable) and save output
dist.pi = vegdist(pi_all2, method = "bray") #abundance data frame - bray curtis dissimilarity
dist.bio1 = dist(bio1, method = "euclidean") #environmental vector - euclidean distance
pi_bio1 = mantel(dist.pi, dist.bio1, method = "spearman", permutations = 9999, na.rm = TRUE) #pi vs bio1 
pi_bio1
capture.output(pi_bio1, file = "~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/Upust_temperature_Mantel.txt")

# check the stress plot
stressplot(nmds1)
plot(nmds1)
plot(nmds1,display=c("sites"))

# plot elevation/bio1 using ordisurf and label the plot with labels1
pdf(file.path("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/NMDS_temperature_nucl_div_Upust_15pops_elevation.pdf"))
ordisurf(nmds1,elevation,main="",col="forestgreen")
orditorp(nmds1,display="sites",col=c("red","red","red","red","purple","blue", "red","blue","blue", "red","red","purple","blue","blue","blue"),
         labels = labels1, air=0.01,cex=1.25)
dev.off()

pdf(file.path("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/NMDS_temperature_nucl_div_Upust_15pops_BIO1.pdf"))
ordisurf(nmds1,bio1,main="",col="black")
#orditorp(nmds1,display="species",col="gray",air=0.01)
orditorp(nmds1,display="sites",col=c("red","red","red","red","purple","blue", "red","blue","blue", "red","red","purple","blue","blue","blue"),
         labels = labels1, air=0.01,cex=1.25)
dev.off()

# save a supplemental figure with all bioclim vectors
pdf(file.path("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/NMDS_temperature_nucl_div_Upust_15pops_w_bioclim_vectors.pdf"))
ordiplot(nmds1,type="n")
orditorp(nmds1,display="sites",col=c("red","red","red","red","purple","blue", "red","blue","blue", "red","red","purple","blue","blue","blue"),
         labels = labels1, air=0.01,cex=1.25)
plot(en, col = "gray40")
dev.off()

```

## b. U. phaea 

```{r circadian loci PI NMDS, U. phaea}

library(vegan)

# load and prepare files
pi_all <- read.delim("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/Uph_circadian_GO_all.pi", header = T)
pi_all1 <- pi_all[,-1]
rownames(pi_all1) <- pi_all[,1]
pi_all2 <- t(pi_all1)

bioclim <- read.csv("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/bioclim_11pops_raw.csv", header = T)
bioclim$ID <- paste(bioclim$site, bioclim$pop, sep = '_') # make id column
bioclim1 <- bioclim[,-c(1:2,15)] 
rownames(bioclim1) <- bioclim[,15]
bioclim2 <- as.matrix(bioclim1)

elevation <- bioclim2[,12]
bio1 <- bioclim2[,1]

labels1 <- c("SN-1", "SN-2", "SN-3", "SN-4","MJ-1","MJ-2","MJ-3","MJ-4","MJ-5","MJ-6","MJ-7")

# run NMDS
nmds1 = metaMDS(pi_all2, # our data matrix (may need to be transformed with t() in order to get a matrix?)
                k=2) # The number of reduced dimensions

# fit our environmental data and save envfit output
en = envfit(nmds1, bioclim2, permutations = 999, na.rm = TRUE)
en
capture.output(en, file = "~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/Uph_circadian_envfit.txt")

# run ANOSIM and save output
treat=c("red","red","purple","blue", "red","red","red","purple","blue","blue","blue")
ano = anosim(pi_all2, treat, distance = "bray", permutations = 9999)
ano
capture.output(ano, file = "~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/Uph_circadian_ANOSIM.txt")


# run Mantel test (bio1 as continuous variable) and save output
dist.pi = vegdist(pi_all2, method = "bray") #abundance data frame - bray curtis dissimilarity
dist.bio1 = dist(bio1, method = "euclidean") #environmental vector - euclidean distance
pi_bio1 = mantel(dist.pi, dist.bio1, method = "spearman", permutations = 9999, na.rm = TRUE) #pi vs bio1 
pi_bio1
capture.output(pi_bio1, file = "~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/Uph_circadian_Mantel.txt")

# check the stress plot
stressplot(nmds1)
plot(nmds1)
plot(nmds1,display=c("sites"))

# plot elevation/bio1 using ordisurf and label the plot with labels1
pdf(file.path("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/NMDS_circadian_nucl_div_Uph_11pops_elevation.pdf"))
ordisurf(nmds1,elevation,main="",col="forestgreen")
orditorp(nmds1,display="sites",col=c("red","red","purple","blue", "red","red","red","purple","blue","blue","blue"),
         labels = labels1, air=0.01,cex=1.25)
dev.off()

pdf(file.path("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/NMDS_circadian_nucl_div_Uph_11pops_BIO1.pdf"))
ordisurf(nmds1,bio1,main="",col="black")
#orditorp(nmds1,display="species",col="gray",air=0.01)
orditorp(nmds1,display="sites",col=c("red","red","purple","blue", "red","red","red","purple","blue","blue","blue"),
         labels = labels1, air=0.01,cex=1.25)
dev.off()

# save a supplemental figure with all bioclim vectors
pdf(file.path("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/NMDS_circadian_nucl_div_Uph_11pops_w_bioclim_vectors.pdf"))
ordiplot(nmds1,type="n")
orditorp(nmds1,display="sites",col=c("red","red","purple","blue", "red","red","red","purple","blue","blue","blue"),
         labels = labels1, air=0.01,cex=1.25)
plot(en, col = "gray40")
dev.off()

```


```{r temperature loci PI NMDS, U. phaea}

library(vegan)

# load and prepare files
pi_all <- read.delim("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/Uph_temperature_GO_all.pi", header = T)
pi_all1 <- pi_all[,-1]
rownames(pi_all1) <- pi_all[,1]
pi_all2 <- t(pi_all1)

bioclim <- read.csv("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/bioclim_11pops_raw.csv", header = T)
bioclim$ID <- paste(bioclim$site, bioclim$pop, sep = '_') # make id column
bioclim1 <- bioclim[,-c(1:2,15)] 
rownames(bioclim1) <- bioclim[,15]
bioclim2 <- as.matrix(bioclim1)

elevation <- bioclim2[,12]
bio1 <- bioclim2[,1]

labels1 <- c("SN-1", "SN-2", "SN-3", "SN-4","MJ-1","MJ-2","MJ-3","MJ-4","MJ-5","MJ-6","MJ-7")

# run NMDS
nmds1 = metaMDS(pi_all2, # our data matrix (may need to be transformed with t() in order to get a matrix?)
                k=2) # The number of reduced dimensions

# fit our environmental data and save envfit output
en = envfit(nmds1, bioclim2, permutations = 999, na.rm = TRUE)
en
capture.output(en, file = "~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/Uph_temperature_envfit.txt")

# run ANOSIM and save output
treat=c("red","red","purple","blue", "red","red","red","purple","blue","blue","blue")
ano = anosim(pi_all2, treat, distance = "bray", permutations = 9999)
ano
capture.output(ano, file = "~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/Uph_temperature_ANOSIM.txt")


# run Mantel test (bio1 as continuous variable) and save output
dist.pi = vegdist(pi_all2, method = "bray") #abundance data frame - bray curtis dissimilarity
dist.bio1 = dist(bio1, method = "euclidean") #environmental vector - euclidean distance
pi_bio1 = mantel(dist.pi, dist.bio1, method = "spearman", permutations = 9999, na.rm = TRUE) #pi vs bio1 
pi_bio1
capture.output(pi_bio1, file = "~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/Uph_temperature_Mantel.txt")

# check the stress plot
stressplot(nmds1)
plot(nmds1)
plot(nmds1,display=c("sites"))

# plot elevation/bio1 using ordisurf and label the plot with labels1
pdf(file.path("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/NMDS_temperature_nucl_div_Uph_11pops_elevation.pdf"))
ordisurf(nmds1,elevation,main="",col="forestgreen")
orditorp(nmds1,display="sites",col=c("red","red","purple","blue", "red","red","red","purple","blue","blue","blue"),
         labels = labels1, air=0.01,cex=1.25)
dev.off()

pdf(file.path("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/NMDS_temperature_nucl_div_Uph_11pops_BIO1.pdf"))
ordisurf(nmds1,bio1,main="",col="black")
#orditorp(nmds1,display="species",col="gray",air=0.01)
orditorp(nmds1,display="sites",col=c("red","red","purple","blue", "red","red","red","purple","blue","blue","blue"),
         labels = labels1, air=0.01,cex=1.25)
dev.off()

# save a supplemental figure with all bioclim vectors
pdf(file.path("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/NMDS_temperature_nucl_div_Uph_11pops_w_bioclim_vectors.pdf"))
ordiplot(nmds1,type="n")
orditorp(nmds1,display="sites",col=c("red","red","purple","blue", "red","red","red","purple","blue","blue","blue"),
         labels = labels1, air=0.01,cex=1.25)
plot(en, col = "gray40")
dev.off()

```

