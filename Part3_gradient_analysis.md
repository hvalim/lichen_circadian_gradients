---
title: "Part 3: analysis of circadian alleles along elevation gradients"
author: "Henrique Valim"
date: "4/20/2022"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Setup

First loading packages and data for visualization of (major) allele frequencies across the three gradients:

1. IT: Italy (Limbara, Sardinia)
2. ES1: Spain 1 (Gredos)
3. ES2: Spain 2 (Puerto de Pico)

To get the data for the number of each specific base pair, we turn to the sync file, the main input files for popoolation2. They contain the allele frequencies for every population and every base in the reference genome. We can then extract just the alleles that we identified from the *_rc file that demonstrated a similar turnover in the identity of the major allele (by looking at the "major_allele.maa" column).

## Full pipeline: from significant Fisher's exact test scores across extremes to gradient analysis

Now that the above method is validated (extracting allele frequencies from the sync file), we will combine this approach with an Fst analysis to pull out significant loci (rather than by eye, as done above).

To do this, we will first extract the circadian loci of interest and then use a Fisher's exact test analysis performed between and within the two climatic extremes (i.e. ignoring intermediate populations).

```{bash extracting circadian genes from sync file}

# U. pustulata

awk -F'\t' '$1~/scaffold_12/ && $2 >= 838899 && $2 <= 842573' /phylodata/fdalgrande/TRIMMED_PoolSeq_3gradients/subsampled_30_filtered.IT123456_ESii123456_ESi123.sync > Upust_WC1_alleles.sync
awk -F'\t' '$1~/scaffold_1/ && $2 >= 3239713 && $2 <= 3241276' /phylodata/fdalgrande/TRIMMED_PoolSeq_3gradients/subsampled_30_filtered.IT123456_ESii123456_ESi123.sync > Upust_WC2_alleles.sync
awk -F'\t' '$1~/scaffold_1/ && $2 >= 1991889 && $2 <= 1994328' /phylodata/fdalgrande/TRIMMED_PoolSeq_3gradients/subsampled_30_filtered.IT123456_ESii123456_ESi123.sync > Upust_WC_like1_alleles.sync
awk -F'\t' '$1~/scaffold_6/ && $2 >= 3773 && $2 <= 6045' /phylodata/fdalgrande/TRIMMED_PoolSeq_3gradients/subsampled_30_filtered.IT123456_ESii123456_ESi123.sync > Upust_WC_like2_alleles.sync
awk -F'\t' '$1~/scaffold_9/ && $2 >= 515221 && $2 <= 518690' /phylodata/fdalgrande/TRIMMED_PoolSeq_3gradients/subsampled_30_filtered.IT123456_ESii123456_ESi123.sync > Upust_FRHa_alleles.sync
awk -F'\t' '$1~/scaffold_5/ && $2 >= 1250325 && $2 <= 1254491' /phylodata/fdalgrande/TRIMMED_PoolSeq_3gradients/subsampled_30_filtered.IT123456_ESii123456_ESi123.sync > Upust_FRHb_alleles.sync
awk -F'\t' '$1~/scaffold_12/ && $2 >= 981549 && $2 <= 984260' /phylodata/fdalgrande/TRIMMED_PoolSeq_3gradients/subsampled_30_filtered.IT123456_ESii123456_ESi123.sync > Upust_FRQ_alleles.sync
awk -F'\t' '$1~/scaffold_14/ && $2 >= 221669 && $2 <= 223362' /phylodata/fdalgrande/TRIMMED_PoolSeq_3gradients/subsampled_30_filtered.IT123456_ESii123456_ESi123.sync > Upust_CK1_alleles.sync
awk -F'\t' '$1~/scaffold_13/ && $2 >= 393821 && $2 <= 394250' /phylodata/fdalgrande/TRIMMED_PoolSeq_3gradients/subsampled_30_filtered.IT123456_ESii123456_ESi123.sync > Upust_CCG1_alleles.sync
awk -F'\t' '$1~/scaffold_24/ && $2 >= 136236 && $2 <= 136960' /phylodata/fdalgrande/TRIMMED_PoolSeq_3gradients/subsampled_30_filtered.IT123456_ESii123456_ESi123.sync > Upust_CCG6_alleles.sync
awk -F'\t' '$1~/scaffold_5/ && $2 >= 293617 && $2 <= 295671' /phylodata/fdalgrande/TRIMMED_PoolSeq_3gradients/subsampled_30_filtered.IT123456_ESii123456_ESi123.sync > Upust_CCG7_alleles.sync
awk -F'\t' '$1~/scaffold_4/ && $2 >= 131206 && $2 <= 132764' /phylodata/fdalgrande/TRIMMED_PoolSeq_3gradients/subsampled_30_filtered.IT123456_ESii123456_ESi123.sync > Upust_CCG8_alleles.sync
awk -F'\t' '$1~/scaffold_7/ && $2 >= 1769080 && $2 <= 1771706' /phylodata/fdalgrande/TRIMMED_PoolSeq_3gradients/subsampled_30_filtered.IT123456_ESii123456_ESi123.sync > Upust_CCG9_alleles.sync

# Finally, combine all files together:
cat Upust*alleles.sync > Upust_circadian_alleles.sync


# U. phaea

awk -F'\t' '$1~/scaffold_3/ && $2 >= 1785861 && $2 <= 1789495' subsampled_30_filtered.Uph_both_combined.sync > Uph_WC1_alleles.sync
awk -F'\t' '$1~/scaffold_9/ && $2 >= 1270142 && $2 <= 1271703' subsampled_30_filtered.Uph_both_combined.sync > Uph_WC2_alleles.sync
awk -F'\t' '$1~/scaffold_15/ && $2 >= 807734 && $2 <= 811227' subsampled_30_filtered.Uph_both_combined.sync > Uph_FRHa_alleles.sync
awk -F'\t' '$1~/scaffold_12/ && $2 >= 448168 && $2 <= 452096' subsampled_30_filtered.Uph_both_combined.sync > Uph_FRHb_alleles.sync
awk -F'\t' '$1~/scaffold_2/ && $2 >= 1364622 && $2 <= 1366037' subsampled_30_filtered.Uph_both_combined.sync > Uph_CK1a_alleles.sync
awk -F'\t' '$1~/scaffold_20/ && $2 >= 110017 && $2 <= 111112' subsampled_30_filtered.Uph_both_combined.sync > Uph_CK1b_alleles.sync
awk -F'\t' '$1~/scaffold_15/ && $2 >= 309342 && $2 <= 310926' subsampled_30_filtered.Uph_both_combined.sync > Uph_CK1c_alleles.sync
awk -F'\t' '$1~/scaffold_19/ && $2 >= 801922 && $2 <= 802328' subsampled_30_filtered.Uph_both_combined.sync > Uph_CCG1_alleles.sync
awk -F'\t' '$1~/scaffold_13/ && $2 >= 351947 && $2 <= 353378' subsampled_30_filtered.Uph_both_combined.sync > Uph_CCG7_alleles.sync
awk -F'\t' '$1~/scaffold_17/ && $2 >= 793512 && $2 <= 795138' subsampled_30_filtered.Uph_both_combined.sync > Uph_CCG8_alleles.sync
awk -F'\t' '$1~/scaffold_18/ && $2 >= 1025947 && $2 <= 1028342' subsampled_30_filtered.Uph_both_combined.sync > Uph_CCG9_alleles.sync


# Finally, combine all files together:
cat *alleles.sync > Uph_circadian_alleles.sync

```

We can now run Fisher's exact test (FET) on just the circadian loci of interest:

```{bash Fisher's Exact Test (FET)}

# U. pustulata

perl /home/hvalim/tools/popoolation2_1201/fisher-test.pl --input Upust_circadian_alleles.sync --output Upust_circadian_alleles3.fet --min-count 3 --min-coverage 10 --max-coverage 100 --suppress-noninformative

# tried this for the not downsampled fisher's exact test:
perl /home/hvalim/tools/popoolation2_1201/fisher-test.pl --input Upust_circadian_alleles2.sync --output Upust_circadian_alleles2.fet --min-count 3 --min-coverage 3 --max-coverage 200 --suppress-noninformative


# U. phaea

perl /home/hvalim/tools/popoolation2_1201/fisher-test.pl --input Uph_circadian_alleles.sync --output Uph_circadian_alleles3.fet --min-count 3 --min-coverage 10 --max-coverage 100 --suppress-noninformative

```


Now, to continue on to extracting significant loci from each circadian clock gene

#### U. pustulata

```{r Part 1: finding significantly variable regions across gradient extremes in circadian clock genes}
# install.packages("ggplot2")
# install.packages("emmeans")
# install.packages("dplyr")
# install.packages("plyr")
# install.packages("tidyr")
# install.packages("reshape2")
library(ggplot2)
library(emmeans)
library(dplyr)
library(plyr)
library(tidyr)
library(reshape2)

# load FET (Fisher's exact test) significances for all circadian SNPs
Fst_scores_snps <- read.table("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/5_Gradient_analysis_Upustulata/Upust_circadian_alleles4.fet", sep = "\t", header = F, check.names = F)

# clean up file for R usage
names(Fst_scores_snps) <- Fst_scores_snps[1,]

# rename columns that will be used
colnames(Fst_scores_snps)[1] <- "chr"
colnames(Fst_scores_snps)[2] <- "pos"
colnames(Fst_scores_snps)[5] <- "cov"

names(Fst_scores_snps) <- gsub("=.*$", "", names(Fst_scores_snps))


# replace colons with dashes to match bioclim data
names(Fst_scores_snps) <- gsub("X", "", names(Fst_scores_snps))
names(Fst_scores_snps) <- gsub("\\:", "-", names(Fst_scores_snps))

# remove extra columns
Fst_scores_snps1 <- Fst_scores_snps[,-c(3:4)]

# save the col.names for later
col_names <- colnames(Fst_scores_snps1)

# clean up data columns
Fst_scores_snps1_cleaned <- as.data.frame(lapply(Fst_scores_snps1, function(x) gsub("^.*=","",x)))
names(Fst_scores_snps1_cleaned) <- col_names

# melt the data frame
melted <- melt(Fst_scores_snps1_cleaned, id = c("chr","pos","cov"), variable.name = "Pairwise", value.name = "FET")

# merge with Pairwise labels from  Fst per-gene data set (we are only using the labels, so it doesn't matter that only FRQ, WC1 and WC2 are in this data set)
Fst_scores__genes_circadian <- read.csv("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/5_Gradient_analysis_Upustulata/FRQ_WC1_WC2_Fst_data.csv", header = T)
Pairwise_labels <- Fst_scores__genes_circadian[,c(1:4)]
merged <- merge(melted, Pairwise_labels, by ="Pairwise")
merged$pos <- as.factor(merged$pos) 
merged$pos_num <- as.numeric(as.character(merged$pos))

merged_unique <- merged %>% group_by(pos,chr) %>% slice(which.max(FET))

# get a list of highest values for each position and gene (should all be mostly in Across) which are above -log(0.05)
merged_unique_sig <- merged %>% group_by(pos,chr) %>% slice(which.max(FET)) %>% slice(which(FET>2.996))

merged_unique_sig1 <- merged_unique_sig[merged_unique_sig$Elevation_Type == "Mixed",]

#merged_unique_sig2 <- merged_unique_sig[merged_unique_sig$Pairwise == "1-6" | merged_unique_sig$Pairwise == "1-12" | merged_unique_sig$Pairwise == "5-15" | 
#                                          merged_unique_sig$Pairwise == "6-7" | merged_unique_sig$Pairwise == "7-12" | merged_unique_sig$Pairwise == "7-15" | 
#                                          merged_unique_sig$Pairwise == "6-13" | merged_unique_sig$Pairwise == "12-13" | merged_unique_sig$Pairwise == "13-15" ,]

# write this file 
write.table(merged, file = "~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/5_Gradient_analysis_Upustulata/SNPs_circadian_unique_all4.fet", sep = "\t", na = "", quote = F, row.names = F, col.names = T)
write.table(merged_unique_sig1, file = "~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/5_Gradient_analysis_Upustulata/SNPs_circadian_significant4.fet", sep = "\t", na = "", quote = F, row.names = F, col.names = T)

```

Before continuing, we can make our FET graphs:

```{r Part 1b: FET figures}

subset_WC1 <- merged[merged$chr == "scaffold_12" & merged$pos_num_num > 838899 & merged$pos_num_num < 842573,]
subset_WC1$region <- "WC-1"
subset_WC_like1 <- merged[merged$chr == "scaffold_1" & merged$pos_num > 1991889 & merged$pos_num < 1994328,]
subset_WC_like1$region <- "WC-like_1"
subset_WC_like2 <- merged[merged$chr == "scaffold_6" & merged$pos_num > 3773 & merged$pos_num < 6045,]
subset_WC_like2$region <- "WC-like_2"
subset_FRHa <- merged[merged$chr == "scaffold_9" & merged$pos_num > 515221 & merged$pos_num < 518690,]
subset_FRHa$region <- "FRHa"
subset_FRHb <- merged[merged$chr == "scaffold_5" & merged$pos_num > 1250325 & merged$pos_num < 1254491,]
subset_FRHb$region <- "FRHb"
subset_FRQ <- merged[merged$chr == "scaffold_12" & merged$pos_num > 981549 & merged$pos_num < 984260,]
subset_FRQ$region <- "FRQ"
subset_CCG1 <- merged[merged$chr == "scaffold_13" & merged$pos_num > 393821 & merged$pos_num < 394250,]
subset_CCG1$region <- "CCG1"
subset_CCG6 <- merged[merged$chr == "scaffold_24" & merged$pos_num > 136236 & merged$pos_num < 136960,]
subset_CCG6$region <- "CCG6"
subset_CCG7 <- merged[merged$chr == "scaffold_5" & merged$pos_num > 293617 & merged$pos_num < 295671,]
subset_CCG7$region <- "CCG7"
subset_CCG8 <- merged[merged$chr == "scaffold_4" & merged$pos_num > 131206 & merged$pos_num < 132764,]
subset_CCG8$region <- "CCG8"
subset_CCG9 <- merged[merged$chr == "scaffold_7" & merged$pos_num > 1769080 & merged$pos_num < 1771706,]
subset_CCG9$region <- "CCG9"
subset_WC2 <- merged[merged$chr == "scaffold_1" & merged$pos_num > 3239713 & merged$pos_num < 3241276,]
subset_WC2$region <- "WC-2"
subset_CK1 <- merged[merged$chr == "scaffold_9" & merged$pos_num > 837740 & merged$pos_num < 839329,]
subset_CK1$region <- "CK1"

subset_all <- rbind(subset_FRHa,subset_FRHb,subset_FRQ,subset_WC2,subset_WC_like2,subset_CCG1,subset_CCG6,subset_CCG7,subset_CCG8,subset_CCG9)
subset_all$region <- as.factor(subset_all$region)

subset_all1 <- subset_all
subset_all1$FET <- as.numeric(subset_all1$FET)

# NOTE: WC1, WC2, and FRQ are in reverse relative to genome positions
(plot1  <- ggplot(subset_all1, aes(x = pos_num, y = FET, col = Elevation_Binary)) +
    #geom_point() +
    stat_summary(geom = "line") +
    scale_x_continuous("SNP position") +
    scale_y_continuous("Fisher's exact test scores (Fst)") +
    scale_color_discrete("Zones") +
    geom_hline(yintercept = -log(0.05), col = "black", linetype = "dashed") +
    theme_minimal() +
    facet_grid(.~region, scales = "free") +
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = T),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black")))

ggsave(plot1, filename = "~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/5_Gradient_analysis_Upustulata/circadian_FET_SNPs.pdf", 
       height = 3, width = 15, 
       device = cairo_pdf)


bioclim <- read.csv("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/5_Gradient_analysis_Upustulata/bioclim_15pops_raw.csv", header = T, check.names = F)

(plot1  <- ggplot(bioclim, aes(x = site, y = Elevation, col = BIO1)) +
    geom_point() +
    stat_summary(geom = "line") +
    scale_x_discrete("Site") +
    scale_y_continuous("Elevation (m)") +
    scale_color_continuous("BIO1") +
    theme_minimal(base_size=12) +
    #facet_grid(.~region, scales = "free") +
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = T),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black")))

library("RColorBrewer")
myPalette <- colorRampPalette(brewer.pal(11, "Spectral"))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(1, 1600))

(plot2  <- ggplot(bioclim, aes(x = site, y = Elevation, col = BIO12)) +
    geom_point(size = 3) +
    stat_summary(geom = "line") +
    scale_x_discrete("Site") +
    scale_y_continuous("Elevation (m)") +
    scale_color_continuous("BIO12 (mm)") +
    theme_minimal(base_size=12) +
    #facet_grid(.~region, scales = "free") +
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = T),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black")) +sc)

ggsave(plot1, filename = "~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/5_Gradient_analysis_Upustulata/Elevation_gradients.pdf", 
       height = 4, width = 4, 
       device = cairo_pdf)

ggsave(plot2, filename = "~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/5_Gradient_analysis_Upustulata/Elevation_gradients_precip.pdf", 
       height = 4, width = 4, 
       device = cairo_pdf)

```

Next, we use the extreme SNPs to pull out the relevant allele frequency data from the sync file (in this case, a subset of the sync file with only the circadian loci of interest).

```{r Part 2: creating allele ratio file from sync file -- Spearman sig. results}
# install.packages("ggplot2")
# install.packages("emmeans")
# install.packages("dplyr")
# install.packages("plyr")
# install.packages("tidyr")
# install.packages("reshape2")
library(ggplot2)
library(emmeans)
library(dplyr)
library(plyr)
library(tidyr)
library(reshape2)

# load circadian sync file and unique, significant circadian SNPs
sync <- read.table("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/5_Gradient_analysis_Upustulata/Upust_circadian_alleles4.sync", sep = "\t", header = F, check.names = F)
#unique <- read.table("SNPs_circadian_significant2.fet", sep = "\t", header = T, check.names = F)
unique <- read.table("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/5_Gradient_analysis_Upustulata/SNPs_circadian_unique_all4.fet", sep = "\t", header = T, check.names = F)

unique_subset <- unique[,c(2:4)]

results <- unique_subset %>% dplyr::count(unique_subset$chr,unique_subset$pos) 

# label columns
colnames(sync) <- c("chr","pos","rc","IT_1","IT_2","IT_3","IT_4","IT_5","IT_6","ES1_1","ES1_2","ES1_3","ES1_4","ES1_5","ES1_6","ES2_1","ES2_2","ES2_3")

# filter based on positions in unique df
sync_significant <- sync %>%
      filter(pos %in% unique$pos)

# now we add in the coverage information
sync_significant1 <- merge(sync_significant, unique_subset, by =c("chr","pos"))

# drop sync now since it takes up a lot of memory
rm(sync)
gc()

# melt data frame
sync_significant2 <- melt(sync_significant1, id = c("chr","pos","cov","rc"))

# separate value and value columns
mdata <- sync_significant2 %>% 
      separate(variable, sep = "_", into = c("site","pop")) %>% 
      separate(value, sep = ":", into = c("A","T","C","G","N","del"))

# re-melt data, and set variable types properly
mdata1 <- melt(mdata, id=c("chr","pos","cov","rc","site","pop")) 
mdata1$pop <- as.numeric(mdata1$pop)
#mdata1$pos <- as.factor(mdata1$pos)
mdata1$cov <- as.numeric(mdata1$cov)
mdata1$value <- as.numeric(mdata1$value)
mdata1$pos <- as.numeric(mdata1$pos)

# create an actual ratio (out of the 30 subsampled populations)
mdata1$value_ratio <- mdata1$value / mdata1$cov

# remove N and del columns, which are not necessary
mdata1.1 <- mdata1[mdata1$variable != "N" & mdata1$variable != "del",]

# subset all individual gene data sets
# WC1
subset_WC1 <- mdata1.1[mdata1.1$chr == "scaffold_12" & mdata1.1$pos > 838899 & mdata1.1$pos < 842573,]
subset_WC1$region <- "WC-1"
# WC2
subset_WC2 <- mdata1.1[mdata1.1$chr == "scaffold_1" & mdata1.1$pos > 3239713 & mdata1.1$pos < 3241276,]
subset_WC2$region <- "WC-2"
# WC-like 1
subset_WC_like1 <- mdata1.1[mdata1.1$chr == "scaffold_1" & mdata1.1$pos > 1991889 & mdata1.1$pos < 1994328,]
subset_WC_like1$region <- "WC-like_1"
# WC-like 1
subset_WC_like2 <- mdata1.1[mdata1.1$chr == "scaffold_6" & mdata1.1$pos > 3773 & mdata1.1$pos < 6045,]
subset_WC_like2$region <- "WC-like_2"

# FRHa
subset_FRHa <- mdata1.1[mdata1.1$chr == "scaffold_9" & mdata1.1$pos > 515221 & mdata1.1$pos < 518690,]
subset_FRHa$region <- "FRHa"
# FRHb
subset_FRHb <- mdata1.1[mdata1.1$chr == "scaffold_5" & mdata1.1$pos > 1250325 & mdata1.1$pos < 1254491,]
subset_FRHb$region <- "FRHb"
# FRQ
subset_FRQ <- mdata1.1[mdata1.1$chr == "scaffold_12" & mdata1.1$pos > 981549 & mdata1.1$pos < 984260,]
subset_FRQ$region <- "FRQ"

# CCG1
subset_CCG1 <- mdata1.1[mdata1.1$chr == "scaffold_13" & mdata1.1$pos > 393821 & mdata1.1$pos < 394250,]
subset_CCG1$region <- "CCG1"
# CCG6
subset_CCG6 <- mdata1.1[mdata1.1$chr == "scaffold_24" & mdata1.1$pos > 136236 & mdata1.1$pos < 136960,]
subset_CCG6$region <- "CCG6"
# CCG7
subset_CCG7 <- mdata1.1[mdata1.1$chr == "scaffold_5" & mdata1.1$pos > 293617 & mdata1.1$pos < 295671,]
subset_CCG7$region <- "CCG7"
# CCG8
subset_CCG8 <- mdata1.1[mdata1.1$chr == "scaffold_4" & mdata1.1$pos > 131206 & mdata1.1$pos < 132764,]
subset_CCG8$region <- "CCG8"
# CCG9
subset_CCG9 <- mdata1.1[mdata1.1$chr == "scaffold_7" & mdata1.1$pos > 1769080 & mdata1.1$pos < 1771706,]
subset_CCG9$region <- "CCG9"

# CK1
subset_CK1 <- mdata1.1[mdata1.1$chr == "scaffold_9" & mdata1.1$pos > 837740 & mdata1.1$pos < 839329,]
subset_CK1$region <- "CK1"

# re-join subsets with their new region labels
subset_all <- rbind(subset_FRHa,subset_FRHb,subset_FRQ,subset_WC2,subset_WC_like2,subset_CCG1,subset_CCG6,subset_CCG7,subset_CCG8,subset_CCG9)
subset_all$region <- as.factor(subset_all$region)


# write this file 
write.table(subset_all, file = "~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/5_Gradient_analysis_Upustulata/circadian_alleles_sync_labeled_melted4.txt", sep = "\t", na = "", quote = F, row.names = F, col.names = T)

```


```{r Part 2b: creating allele ratio file from sync file -- FET sig. results}
# install.packages("ggplot2")
# install.packages("emmeans")
# install.packages("dplyr")
# install.packages("plyr")
# install.packages("tidyr")
# install.packages("reshape2")
library(ggplot2)
library(emmeans)
library(dplyr)
library(plyr)
library(tidyr)
library(reshape2)

# load circadian sync file and unique, significant circadian SNPs
sync <- read.table("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/5_Gradient_analysis_Upustulata/Upust_circadian_alleles4.sync", sep = "\t", header = F, check.names = F)
#unique <- read.table("SNPs_circadian_significant2.fet", sep = "\t", header = T, check.names = F)
unique_sig <- read.table("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/5_Gradient_analysis_Upustulata/SNPs_circadian_significant4.fet", sep = "\t", header = T, check.names = F)
#unique_sig <- read.table("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/5_Gradient_analysis_Upustulata/SNPs_circadian_significant4b.fet", sep = "\t", header = T, check.names = F)

unique_subset <- unique_sig[,c(2:4)]

results <- unique_subset %>% dplyr::count(unique_subset$chr,unique_subset$pos) 

# label columns
colnames(sync) <- c("chr","pos","rc","IT_1","IT_2","IT_3","IT_4","IT_5","IT_6","ES1_1","ES1_2","ES1_3","ES1_4","ES1_5","ES1_6","ES2_1","ES2_2","ES2_3")

# filter based on positions in unique df
sync_significant <- sync %>%
      filter(pos %in% unique_sig$pos)

# now we add in the coverage information
sync_significant1 <- merge(sync_significant, unique_subset, by =c("chr","pos"))

# drop sync now since it takes up a lot of memory
rm(sync)
gc()

# melt data frame
sync_significant2 <- melt(sync_significant1, id = c("chr","pos","cov","rc"))

# separate value and value columns
mdata <- sync_significant2 %>% 
      separate(variable, sep = "_", into = c("site","pop")) %>% 
      separate(value, sep = ":", into = c("A","T","C","G","N","del"))

# re-melt data, and set variable types properly
mdata1 <- melt(mdata, id=c("chr","pos","cov","rc","site","pop")) 
mdata1$pop <- as.numeric(mdata1$pop)
#mdata1$pos <- as.factor(mdata1$pos)
mdata1$cov <- as.numeric(mdata1$cov)
mdata1$value <- as.numeric(mdata1$value)
mdata1$pos <- as.numeric(mdata1$pos)

# create an actual ratio (out of the 30 subsampled populations)
mdata1$value_ratio <- mdata1$value / mdata1$cov

# remove N and del columns, which are not necessary
mdata1.1 <- mdata1[mdata1$variable != "N" & mdata1$variable != "del",]

# subset all individual gene data sets
# WC1
subset_WC1 <- mdata1.1[mdata1.1$chr == "scaffold_12" & mdata1.1$pos > 838899 & mdata1.1$pos < 842573,]
subset_WC1$region <- "WC-1"
# WC2
subset_WC2 <- mdata1.1[mdata1.1$chr == "scaffold_1" & mdata1.1$pos > 3239713 & mdata1.1$pos < 3241276,]
subset_WC2$region <- "WC-2"
# WC-like 1
subset_WC_like1 <- mdata1.1[mdata1.1$chr == "scaffold_1" & mdata1.1$pos > 1991889 & mdata1.1$pos < 1994328,]
subset_WC_like1$region <- "WC-like_1"
# WC-like 1
subset_WC_like2 <- mdata1.1[mdata1.1$chr == "scaffold_6" & mdata1.1$pos > 3773 & mdata1.1$pos < 6045,]
subset_WC_like2$region <- "WC-like_2"

# FRHa
subset_FRHa <- mdata1.1[mdata1.1$chr == "scaffold_9" & mdata1.1$pos > 515221 & mdata1.1$pos < 518690,]
subset_FRHa$region <- "FRHa"
# FRHb
subset_FRHb <- mdata1.1[mdata1.1$chr == "scaffold_5" & mdata1.1$pos > 1250325 & mdata1.1$pos < 1254491,]
subset_FRHb$region <- "FRHb"
# FRQ
subset_FRQ <- mdata1.1[mdata1.1$chr == "scaffold_12" & mdata1.1$pos > 981549 & mdata1.1$pos < 984260,]
subset_FRQ$region <- "FRQ"

# CCG1
subset_CCG1 <- mdata1.1[mdata1.1$chr == "scaffold_13" & mdata1.1$pos > 393821 & mdata1.1$pos < 394250,]
subset_CCG1$region <- "CCG1"
# CCG6
subset_CCG6 <- mdata1.1[mdata1.1$chr == "scaffold_24" & mdata1.1$pos > 136236 & mdata1.1$pos < 136960,]
subset_CCG6$region <- "CCG6"
# CCG7
subset_CCG7 <- mdata1.1[mdata1.1$chr == "scaffold_5" & mdata1.1$pos > 293617 & mdata1.1$pos < 295671,]
subset_CCG7$region <- "CCG7"
# CCG8
subset_CCG8 <- mdata1.1[mdata1.1$chr == "scaffold_4" & mdata1.1$pos > 131206 & mdata1.1$pos < 132764,]
subset_CCG8$region <- "CCG8"
# CCG9
subset_CCG9 <- mdata1.1[mdata1.1$chr == "scaffold_7" & mdata1.1$pos > 1769080 & mdata1.1$pos < 1771706,]
subset_CCG9$region <- "CCG9"

# CK1
subset_CK1 <- mdata1.1[mdata1.1$chr == "scaffold_9" & mdata1.1$pos > 837740 & mdata1.1$pos < 839329,]
subset_CK1$region <- "CK1"

# re-join subsets with their new region labels
subset_all <- rbind(subset_FRHa,subset_FRHb,subset_FRQ,subset_WC2,subset_WC_like2,subset_CCG1,subset_CCG6,subset_CCG7,subset_CCG8,subset_CCG9)
subset_all$region <- as.factor(subset_all$region)


# write this file 
write.table(subset_all, file = "~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/5_Gradient_analysis_Upustulata/circadian_alleles_sync_labeled_melted4_sig.txt", sep = "\t", na = "", quote = F, row.names = F, col.names = T)

```

Next, we can analyze these data: identify which loci have trends along the gradient, merge the data with the environmental/climatic data, and finally graph it.

```{r Part 3: allele ratio gradient analysis -- extracting signifant Spearman corr. loci}
# install.packages("ggplot2")
# install.packages("emmeans")
# install.packages("dplyr")
# install.packages("plyr")
# install.packages("tidyr")
# install.packages("reshape2")
# install.packages("emmeans")
# install.packages("multcompView")
# install.packages("corrr")
library(ggplot2)
library(emmeans)
library(dplyr)
library(plyr)
library(tidyr)
library(reshape2)
library(emmeans)
library(multcompView)

# clear the list! 
rm(list=ls())
gc()

# load the labeled, melted circadian sync file
data <- read.table("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/5_Gradient_analysis_Upustulata/circadian_alleles_sync_labeled_melted4.txt", sep = "\t", header = T, check.names = F)

# set variables properly
data$region <- as.factor(data$region)
data$site <- as.factor(data$site)
data$pos <- as.factor(data$pos)
data$pop <- as.numeric(data$pop)
data$variable <- as.factor(data$variable)

# briefly calculate how many significant SNPs we have:
data1 <- data[,c(2,10)]
unique <- data1 %>% distinct()
sig_count <- unique %>% dplyr::count(unique$region) 


## generating for loop to cycle through all correlations of allele ratios across pop for each site

## subset each nucleotide and run each separately
# A
data_A <- data[data$variable == "A",]
data_A$genom_pos <- paste(data_A$chr,data_A$pos)
data_A$genom_pos <- as.factor(data_A$genom_pos)

results_A <- data.frame(row.names = 1:1)

for(i in 1:length(levels(data_A$genom_pos))){
  for(j in 1:length(levels(data_A$site))){
    subset_focus <- subset(data_A, genom_pos == levels(data_A$genom_pos)[i] & site == levels(data_A$site)[j])
    corr_focus <- suppressWarnings(cor(subset_focus$pop,subset_focus$value_ratio, method = "spearman"))
      if(is.na(corr_focus)){
        results_A[1000*(j-1)+i,1] <- levels(data_A$genom_pos)[i]
        results_A[1000*(j-1)+i,2] <- levels(data_A$site)[j]
        results_A[1000*(j-1)+i,2] <- levels(data_A$region)[j]
        results_A[1000*(j-1)+i,3] <- NA
        results_A[1000*(j-1)+i,4] <- NA
      }
    else{
        p_val_focus <- suppressWarnings(cor.test(subset_focus$pop,subset_focus$value_ratio, method = "spearman"))[3][[1]]
        results_A[1000*(j-1)+i,1] <- levels(data_A$genom_pos)[i]
        results_A[1000*(j-1)+i,2] <- levels(data_A$site)[j]
        results_A[1000*(j-1)+i,2] <- levels(data_A$region)[j]
        results_A[1000*(j-1)+i,3] <- corr_focus
        results_A[1000*(j-1)+i,4] <- p_val_focus
    }
  colnames(results_A) <- c("genom_pos","site","spearman_corr","p_val_corr")
  }
}

# T
data_T <- data[data$variable == "T",]
data_T$genom_pos <- paste(data_T$chr,data_T$pos)
data_T$genom_pos <- as.factor(data_T$genom_pos)

results_T <- data.frame(row.names = 1:1)

for(i in 1:length(levels(data_T$genom_pos))){
  for(j in 1:length(levels(data_T$site))){
    subset_focus <- subset(data_T, genom_pos == levels(data_T$genom_pos)[i] & site == levels(data_T$site)[j])
    corr_focus <- suppressWarnings(cor(subset_focus$pop,subset_focus$value_ratio, method = "spearman"))
      if(is.na(corr_focus)){
        results_T[144*(j-1)+i,1] <- levels(data_T$genom_pos)[i]
        results_T[144*(j-1)+i,2] <- levels(data_T$site)[j]
        results_T[144*(j-1)+i,2] <- levels(data_T$region)[j]
        results_T[144*(j-1)+i,3] <- NA
        results_T[144*(j-1)+i,4] <- NA
      }
    else{
        p_val_focus <- suppressWarnings(cor.test(subset_focus$pop,subset_focus$value_ratio, method = "spearman"))[3][[1]]
        results_T[144*(j-1)+i,1] <- levels(data_T$genom_pos)[i]
        results_T[144*(j-1)+i,2] <- levels(data_T$site)[j]
        results_T[144*(j-1)+i,2] <- levels(data_T$region)[j]
        results_T[144*(j-1)+i,3] <- corr_focus
        results_T[144*(j-1)+i,4] <- p_val_focus
    }
  colnames(results_T) <- c("genom_pos","site","spearman_corr","p_val_corr")
  }
}

# C
data_C <- data[data$variable == "C",]
data_C$genom_pos <- paste(data_C$chr,data_C$pos)
data_C$genom_pos <- as.factor(data_C$genom_pos)

results_C <- data.frame(row.names = 1:1)

for(i in 1:length(levels(data_C$genom_pos))){
  for(j in 1:length(levels(data_C$site))){
    subset_focus <- subset(data_C, genom_pos == levels(data_C$genom_pos)[i] & site == levels(data_C$site)[j])
    corr_focus <- suppressWarnings(cor(subset_focus$pop,subset_focus$value_ratio, method = "spearman"))
      if(is.na(corr_focus)){
        results_C[144*(j-1)+i,1] <- levels(data_C$genom_pos)[i]
        results_C[144*(j-1)+i,2] <- levels(data_C$site)[j]
        results_C[144*(j-1)+i,2] <- levels(data_C$region)[j]
        results_C[144*(j-1)+i,3] <- NA
        results_C[144*(j-1)+i,4] <- NA
      }
    else{
        p_val_focus <- suppressWarnings(cor.test(subset_focus$pop,subset_focus$value_ratio, method = "spearman"))[3][[1]]
        results_C[144*(j-1)+i,1] <- levels(data_C$genom_pos)[i]
        results_C[144*(j-1)+i,2] <- levels(data_C$site)[j]
        results_C[144*(j-1)+i,2] <- levels(data_C$region)[j]
        results_C[144*(j-1)+i,3] <- corr_focus
        results_C[144*(j-1)+i,4] <- p_val_focus
    }
  colnames(results_C) <- c("genom_pos","site","spearman_corr","p_val_corr")
  }
}

# G
data_G <- data[data$variable == "G",]
data_G$genom_pos <- paste(data_G$chr,data_T$pos)
data_G$genom_pos <- as.factor(data_G$genom_pos)

results_G <- data.frame(row.names = 1:1)

for(i in 1:length(levels(data_G$genom_pos))){
  for(j in 1:length(levels(data_G$site))){
    subset_focus <- subset(data_G, genom_pos == levels(data_G$genom_pos)[i] & site == levels(data_G$site)[j])
    corr_focus <- suppressWarnings(cor(subset_focus$pop,subset_focus$value_ratio, method = "spearman"))
      if(is.na(corr_focus)){
        results_G[144*(j-1)+i,1] <- levels(data_G$genom_pos)[i]
        results_G[144*(j-1)+i,2] <- levels(data_G$site)[j]
        results_G[144*(j-1)+i,2] <- levels(data_G$region)[j]
        results_G[144*(j-1)+i,3] <- NA
        results_G[144*(j-1)+i,4] <- NA
      }
    else{
        p_val_focus <- suppressWarnings(cor.test(subset_focus$pop,subset_focus$value_ratio, method = "spearman"))[3][[1]]
        results_G[144*(j-1)+i,1] <- levels(data_G$genom_pos)[i]
        results_G[144*(j-1)+i,2] <- levels(data_G$site)[j]
        results_G[144*(j-1)+i,2] <- levels(data_G$region)[j]
        results_G[144*(j-1)+i,3] <- corr_focus
        results_G[144*(j-1)+i,4] <- p_val_focus
    }
  colnames(results_G) <- c("genom_pos","site","spearman_corr","p_val_corr")
  }
}

# remove all NAs
results_A.1 <- results_A[!is.na(results_A$spearman_corr),]
results_T.1 <- results_T[!is.na(results_T$spearman_corr),]
results_C.1 <- results_C[!is.na(results_C$spearman_corr),]
results_G.1 <- results_G[!is.na(results_G$spearman_corr),]

# add back variable information
results_A.1$variable <- "A"
results_T.1$variable <- "T"
results_C.1$variable <- "C"
results_G.1$variable <- "G"

# re-join subsets with their new variable labels
results_all <- rbind(results_A.1,results_T.1,results_C.1,results_G.1)
results_all$variable <- as.factor(results_all$variable)
results_all$site <- as.factor(results_all$site)

#subset only significant results
results_significant <- results_all[results_all$p_val_corr < 0.05,]

# we can then make a vector of all the significant loci for ES1 and IT (ES2 only has 3 populations in the gradient, so there are very few significant hits, presumably due to the low resolution of the 3-pop gradient)
significant_loci <- results_significant$genom_pos

# join all chr and pos values
data$genom_pos <- paste(data$chr,data$pos)

# join with significant loci
data_significant_subset <- data[data$genom_pos %in% significant_loci, ]

# write this file 
write.table(data_significant_subset, file = "~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/5_Gradient_analysis_Upustulata/circadian_alleles_spearman_sig_0_05_loci4.txt", sep = "\t", na = "", quote = F, row.names = F, col.names = T)

```

Next, it would be interesting to plot these values with their respective bioclim data. 

```{r Part 4: allele ratio gradient analysis -- Spearman sig. results}
# install.packages("ggplot2")
# install.packages("emmeans")
# install.packages("dplyr")
# install.packages("plyr")
# install.packages("tidyr")
# install.packages("reshape2")
# install.packages("emmeans")
# install.packages("multcompView")
# install.packages("corrr")
library(ggplot2)
library(emmeans)
library(dplyr)
library(plyr)
library(tidyr)
library(reshape2)
library(emmeans)
library(multcompView)

# clear the list! 
rm(list=ls())
gc()

# load the file with significant loci and the bioclim file
data <- read.table("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/5_Gradient_analysis_Upustulata/circadian_alleles_spearman_sig_0_05_loci4.txt", sep = "\t", header = T, check.names = F)
bioclim <- read.csv("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/5_Gradient_analysis_Upustulata/bioclim_15pops_raw.csv", header = T, check.names = F)

joined_data <- merge(data,bioclim, by = c("site","pop"))


# looking at loci significant in any of the sites (rather than all three)

# need to manually set the axis range for each of the y-axes for setting an additional x axis:
# ylim.prim <- c(0, 1)   # in this example, allele ratio
# ylim.sec <- c(350, 950)    # in this example, elevation
# b <- diff(ylim.prim)/diff(ylim.sec)
# a <- ylim.prim[1] - b*ylim.sec[1] 
### Then, add this to the plot:
# geom_line(aes(y = a + BIO12*b), color = "cornflowerblue", linetype = "dashed") + # to add the data 
# scale_y_continuous("% base frequency", breaks = c(0,0.5,1), sec.axis = sec_axis(~ (. - a)/b, name = "Mean annual precip. (mm)")) + # modify the scale like this


WC1_joined_data <- joined_data[joined_data$region == "WC-1",]
WC2_joined_data <- joined_data[joined_data$region == "WC-2",]
WC_like1_joined_data <- joined_data[joined_data$region == "WC-like_1",]
WC_like2_joined_data <- joined_data[joined_data$region == "WC-like_2",]
FRQ_joined_data <- joined_data[joined_data$region == "FRQ",]
FRHa_joined_data <- joined_data[joined_data$region == "FRHa",]
FRHb_joined_data <- joined_data[joined_data$region == "FRHb",]
CK1_joined_data <- joined_data[joined_data$region == "CK1",]
CK1c_joined_data <- joined_data[joined_data$region == "CK1c",]
CCG1_joined_data <- joined_data[joined_data$region == "CCG1",]
CCG6_joined_data <- joined_data[joined_data$region == "CCG6",]
CCG7_joined_data <- joined_data[joined_data$region == "CCG7",]
CCG8_joined_data <- joined_data[joined_data$region == "CCG8",]
CCG9_joined_data <- joined_data[joined_data$region == "CCG9",]


(plot1  <- ggplot(WC2_joined_data, aes(x = BIO1, y = value_ratio, color = variable, linetype = site)) +
    geom_point(size = 1) +
    stat_summary(fun.data = mean_se, geom = "line") +
    scale_color_manual("Base", values = c("darkblue","gold","darkred","darkgreen")) +
    scale_x_continuous("Mean annual temp. (Celsius)") +
    scale_y_continuous("% base frequency", breaks = c(0,0.5,1)) +
    facet_grid(genom_pos~., scales = "free") +
    theme_minimal(base_size=10) +
    theme(strip.text.y.right = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(color = "black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black")))

(plot2  <- ggplot(WC_like2_joined_data, aes(x = BIO1, y = value_ratio, color = variable, linetype = site)) +
    geom_point(size = 1) +
    stat_summary(fun.data = mean_se, geom = "line") +
    scale_color_manual("Base", values = c("darkblue","gold","darkred","darkgreen")) +
    scale_x_continuous("Mean annual temp. (Celsius)") +
    scale_y_continuous("% base frequency", breaks = c(0,0.5,1)) +
    facet_grid(genom_pos~., scales = "free") +
    theme_minimal(base_size=10) +
    theme(strip.text.y.right = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(color = "black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black")))

(plot3  <- ggplot(FRHa_joined_data, aes(x = BIO1, y = value_ratio, color = variable, linetype = site)) +
    geom_point(size = 1) +
    stat_summary(fun.data = mean_se, geom = "line") +
    scale_color_manual("Base", values = c("darkblue","gold","darkred","darkgreen")) +
    scale_x_continuous("Mean annual temp. (Celsius)") +
    scale_y_continuous("% base frequency", breaks = c(0,0.5,1)) +
    facet_grid(genom_pos~., scales = "free") +
    theme_minimal(base_size=10) +
    theme(strip.text.y.right = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(color = "black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black")))

FRHa_joined_data_unfixed <- FRHa_joined_data[FRHa_joined_data$genom_pos == "scaffold_9 517371",]
FRHa_joined_data_fixed <- FRHa_joined_data[FRHa_joined_data$genom_pos == "scaffold_9 515299",]
FRHa_joined_data_negative <- FRHa_joined_data[FRHa_joined_data$genom_pos == "scaffold_9 516472",]

FRHa_joined_data_focus <- rbind(FRHa_joined_data_unfixed,FRHa_joined_data_fixed,FRHa_joined_data_negative)


(plot3a  <- ggplot(FRHa_joined_data_focus, aes(x = BIO1, y = value_ratio, color = variable, linetype = site)) +
    geom_point(size = 1) +
    stat_summary(fun.data = mean_se, geom = "line") +
    scale_color_manual("Base", values = c("darkblue","gold","darkred","darkgreen")) +
    scale_x_continuous("Mean annual temp. (Celsius)") +
    scale_y_continuous("% base frequency", breaks = c(0,0.5,1)) +
    facet_grid(.~genom_pos, scales = "free") +
    theme_minimal(base_size=10) +
    theme(strip.text.y.right = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(color = "black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black")))




(plot4  <- ggplot(FRHb_joined_data, aes(x = BIO1, y = value_ratio, color = variable, linetype = site)) +
    geom_point(size = 1) +
    stat_summary(fun.data = mean_se, geom = "line") +
    scale_color_manual("Base", values = c("darkblue","gold","darkred","darkgreen")) +
    scale_x_continuous("Mean annual temp. (Celsius)") +
    scale_y_continuous("% base frequency", breaks = c(0,0.5,1)) +
    facet_grid(genom_pos~., scales = "free") +
    theme_minimal(base_size=10) +
    theme(strip.text.y.right = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(color = "black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black")))

(plot5  <- ggplot(CCG1_joined_data, aes(x = BIO1, y = value_ratio, color = variable, linetype = site)) +
    geom_point(size = 1) +
    stat_summary(fun.data = mean_se, geom = "line") +
    scale_color_manual("Base", values = c("darkblue","gold","darkred","darkgreen")) +
    scale_x_continuous("Mean annual temp. (Celsius)") +
    scale_y_continuous("% base frequency", breaks = c(0,0.5,1)) +
    facet_grid(genom_pos~., scales = "free") +
    theme_minimal(base_size=10) +
    theme(strip.text.y.right = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(color = "black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black")))

(plot6  <- ggplot(CCG6_joined_data, aes(x = BIO1, y = value_ratio, color = variable, linetype = site)) +
    geom_point(size = 1) +
    stat_summary(fun.data = mean_se, geom = "line") +
    scale_color_manual("Base", values = c("darkblue","gold","darkred","darkgreen")) +
    scale_x_continuous("Mean annual temp. (Celsius)") +
    scale_y_continuous("% base frequency", breaks = c(0,0.5,1)) +
    facet_grid(genom_pos~., scales = "free") +
    theme_minimal(base_size=10) +
    theme(strip.text.y.right = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(color = "black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black")))

(plot7  <- ggplot(CCG7_joined_data, aes(x = BIO1, y = value_ratio, color = variable, linetype = site)) +
    geom_point(size = 1) +
    stat_summary(fun.data = mean_se, geom = "line") +
    scale_color_manual("Base", values = c("darkblue","gold","darkred","darkgreen")) +
    scale_x_continuous("Mean annual temp. (Celsius)") +
    scale_y_continuous("% base frequency", breaks = c(0,0.5,1)) +
    facet_grid(genom_pos~., scales = "free") +
    theme_minimal(base_size=10) +
    theme(strip.text.y.right = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(color = "black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black")))

(plot8  <- ggplot(CCG8_joined_data, aes(x = BIO1, y = value_ratio, color = variable, linetype = site)) +
    geom_point(size = 1) +
    stat_summary(fun.data = mean_se, geom = "line") +
    scale_color_manual("Base", values = c("darkblue","gold","darkred","darkgreen")) +
    scale_x_continuous("Mean annual temp. (Celsius)") +
    scale_y_continuous("% base frequency", breaks = c(0,0.5,1)) +
    facet_grid(genom_pos~., scales = "free") +
    theme_minimal(base_size=10) +
    theme(strip.text.y.right = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(color = "black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black")))

(plot9  <- ggplot(CCG9_joined_data, aes(x = BIO1, y = value_ratio, color = variable, linetype = site)) +
    geom_point(size = 1) +
    stat_summary(fun.data = mean_se, geom = "line") +
    scale_color_manual("Base", values = c("darkblue","gold","darkred","darkgreen")) +
    scale_x_continuous("Mean annual temp. (Celsius)") +
    scale_y_continuous("% base frequency", breaks = c(0,0.5,1)) +
    facet_grid(genom_pos~., scales = "free") +
    theme_minimal(base_size=10) +
    theme(strip.text.y.right = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(color = "black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black")))


(plot10  <- ggplot(FRQ_joined_data, aes(x = BIO1, y = value_ratio, color = variable, linetype = site)) +
    geom_point(size = 1) +
    stat_summary(fun.data = mean_se, geom = "line") +
    scale_color_manual("Base", values = c("darkblue","gold","darkred","darkgreen")) +
    scale_x_continuous("Mean annual temp. (Celsius)") +
    scale_y_continuous("% base frequency", breaks = c(0,0.5,1)) +
    facet_grid(genom_pos~., scales = "free") +
    theme_minimal(base_size=10) +
    theme(strip.text.y.right = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(color = "black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black")))


pdf("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/5_Gradient_analysis_Upustulata/WC2_gradient_spearman_sig_loci_x_BIO1_sync.pdf", height = 20, width = 4)
plot1
dev.off() 

pdf("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/5_Gradient_analysis_Upustulata/WC_like2_gradient_spearman_sig_loci_x_BIO1_sync.pdf", height = 40, width = 4)
plot2
dev.off() 

pdf("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/5_Gradient_analysis_Upustulata/FRHa_gradient_spearman_sig_loci_x_BIO1_sync.pdf", height = 40, width = 4)
plot3
dev.off() 

pdf("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/5_Gradient_analysis_Upustulata/FRHa_gradient_spearman_sig_loci_x_BIO1_sync_Fig_focus.pdf", height = 2, width = 6)
plot3a
dev.off() 

pdf("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/5_Gradient_analysis_Upustulata/FRHb_gradient_spearman_sig_loci_x_BIO1_sync.pdf", height = 40, width = 4)
plot4
dev.off() 

pdf("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/5_Gradient_analysis_Upustulata/CCG1_gradient_spearman_sig_loci_x_BIO1_sync.pdf", height = 20, width = 4)
plot5
dev.off() 

pdf("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/5_Gradient_analysis_Upustulata/CCG6_gradient_spearman_sig_loci_x_BIO1_sync.pdf", height = 20, width = 4)
plot6
dev.off() 

pdf("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/5_Gradient_analysis_Upustulata/CCG7_gradient_spearman_sig_loci_x_BIO1_sync.pdf", height = 40, width = 4)
plot7
dev.off() 

pdf("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/5_Gradient_analysis_Upustulata/CCG8_gradient_spearman_sig_loci_x_BIO1_sync.pdf", height = 30, width = 4)
plot8
dev.off() 

pdf("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/5_Gradient_analysis_Upustulata/CCG9_gradient_spearman_sig_loci_x_BIO1_sync.pdf", height = 60, width = 4)
plot9
dev.off() 

pdf("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/5_Gradient_analysis_Upustulata/FRQ_gradient_spearman_sig_loci_x_BIO1_sync.pdf", height = 60, width = 4)
plot10
dev.off() 


```

```{r Part 4b: allele ratio gradient analysis -- FET sig. results}
# install.packages("ggplot2")
# install.packages("emmeans")
# install.packages("dplyr")
# install.packages("plyr")
# install.packages("tidyr")
# install.packages("reshape2")
# install.packages("emmeans")
# install.packages("multcompView")
# install.packages("corrr")
library(ggplot2)
library(emmeans)
library(dplyr)
library(plyr)
library(tidyr)
library(reshape2)
library(emmeans)
library(multcompView)

# clear the list! 
rm(list=ls())
gc()

# load the file with significant loci and the bioclim file
data <- read.table("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/5_Gradient_analysis_Upustulata/circadian_alleles_sync_labeled_melted4_sig.txt", sep = "\t", header = T, check.names = F)
#data <- read.table("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/5_Gradient_analysis_Upustulata/circadian_alleles_sync_labeled_melted4_sig2.txt", sep = "\t", header = T, check.names = F)
bioclim <- read.csv("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/5_Gradient_analysis_Upustulata/bioclim_15pops_raw.csv", header = T, check.names = F)

joined_data <- merge(data,bioclim, by = c("site","pop"))


joined_data$genom_pos <- paste(joined_data$chr,joined_data$pos)
joined_data$genom_pos <- as.factor(joined_data$genom_pos)

# looking at loci significant in any of the sites (rather than all three)

# need to manually set the axis range for each of the y-axes for setting an additional x axis:
# ylim.prim <- c(0, 1)   # in this example, allele ratio
# ylim.sec <- c(350, 950)    # in this example, elevation
# b <- diff(ylim.prim)/diff(ylim.sec)
# a <- ylim.prim[1] - b*ylim.sec[1] 
### Then, add this to the plot:
# geom_line(aes(y = a + BIO12*b), color = "cornflowerblue", linetype = "dashed") + # to add the data 
# scale_y_continuous("% base frequency", breaks = c(0,0.5,1), sec.axis = sec_axis(~ (. - a)/b, name = "Mean annual precip. (mm)")) + # modify the scale like this


WC1_joined_data <- joined_data[joined_data$region == "WC-1",]
WC2_joined_data <- joined_data[joined_data$region == "WC-2",]
WC_like1_joined_data <- joined_data[joined_data$region == "WC-like_1",]
WC_like2_joined_data <- joined_data[joined_data$region == "WC-like_2",]
FRQ_joined_data <- joined_data[joined_data$region == "FRQ",]
FRHa_joined_data <- joined_data[joined_data$region == "FRHa",]
FRHb_joined_data <- joined_data[joined_data$region == "FRHb",]
CK1_joined_data <- joined_data[joined_data$region == "CK1",]
CK1c_joined_data <- joined_data[joined_data$region == "CK1c",]
CCG1_joined_data <- joined_data[joined_data$region == "CCG1",]
CCG6_joined_data <- joined_data[joined_data$region == "CCG6",]
CCG7_joined_data <- joined_data[joined_data$region == "CCG7",]
CCG8_joined_data <- joined_data[joined_data$region == "CCG8",]
CCG9_joined_data <- joined_data[joined_data$region == "CCG9",]


(plot1  <- ggplot(WC2_joined_data, aes(x = BIO1, y = value_ratio, color = variable, linetype = site)) +
    geom_point(size = 1) +
    stat_summary(fun.data = mean_se, geom = "line") +
    scale_color_manual("Base", values = c("darkblue","gold","darkred","darkgreen")) +
    scale_x_continuous("Mean annual temp. (Celsius)") +
    scale_y_continuous("% base frequency", breaks = c(0,0.5,1)) +
    facet_grid(genom_pos~., scales = "free") +
    theme_minimal(base_size=10) +
    theme(strip.text.y.right = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(color = "black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black")))

(plot2  <- ggplot(WC_like2_joined_data, aes(x = BIO1, y = value_ratio, color = variable, linetype = site)) +
    geom_point(size = 1) +
    stat_summary(fun.data = mean_se, geom = "line") +
    scale_color_manual("Base", values = c("darkblue","gold","darkred","darkgreen")) +
    scale_x_continuous("Mean annual temp. (Celsius)") +
    scale_y_continuous("% base frequency", breaks = c(0,0.5,1)) +
    facet_grid(genom_pos~., scales = "free") +
    theme_minimal(base_size=10) +
    theme(strip.text.y.right = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(color = "black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black")))

(plot3  <- ggplot(FRHa_joined_data, aes(x = BIO1, y = value_ratio, color = variable, linetype = site)) +
    geom_point(size = 1) +
    stat_summary(fun.data = mean_se, geom = "line") +
    scale_color_manual("Base", values = c("darkblue","gold","darkred","darkgreen")) +
    scale_x_continuous("Mean annual temp. (Celsius)") +
    scale_y_continuous("% base frequency", breaks = c(0,0.5,1)) +
    facet_grid(genom_pos~., scales = "free") +
    theme_minimal(base_size=10) +
    theme(strip.text.y.right = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(color = "black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black")))

(plot4  <- ggplot(FRHb_joined_data, aes(x = BIO1, y = value_ratio, color = variable, linetype = site)) +
    geom_point(size = 1) +
    stat_summary(fun.data = mean_se, geom = "line") +
    scale_color_manual("Base", values = c("darkblue","gold","darkred","darkgreen")) +
    scale_x_continuous("Mean annual temp. (Celsius)") +
    scale_y_continuous("% base frequency", breaks = c(0,0.5,1)) +
    facet_grid(genom_pos~., scales = "free") +
    theme_minimal(base_size=10) +
    theme(strip.text.y.right = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(color = "black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black")))

(plot5  <- ggplot(CCG1_joined_data, aes(x = BIO1, y = value_ratio, color = variable, linetype = site)) +
    geom_point(size = 1) +
    stat_summary(fun.data = mean_se, geom = "line") +
    scale_color_manual("Base", values = c("darkblue","gold","darkred","darkgreen")) +
    scale_x_continuous("Mean annual temp. (Celsius)") +
    scale_y_continuous("% base frequency", breaks = c(0,0.5,1)) +
    facet_grid(genom_pos~., scales = "free") +
    theme_minimal(base_size=10) +
    theme(strip.text.y.right = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(color = "black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black")))

(plot6  <- ggplot(CCG6_joined_data, aes(x = BIO1, y = value_ratio, color = variable, linetype = site)) +
    geom_point(size = 1) +
    stat_summary(fun.data = mean_se, geom = "line") +
    scale_color_manual("Base", values = c("darkblue","gold","darkred","darkgreen")) +
    scale_x_continuous("Mean annual temp. (Celsius)") +
    scale_y_continuous("% base frequency", breaks = c(0,0.5,1)) +
    facet_grid(genom_pos~., scales = "free") +
    theme_minimal(base_size=10) +
    theme(strip.text.y.right = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(color = "black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black")))

(plot7  <- ggplot(CCG7_joined_data, aes(x = BIO1, y = value_ratio, color = variable, linetype = site)) +
    geom_point(size = 1) +
    stat_summary(fun.data = mean_se, geom = "line") +
    scale_color_manual("Base", values = c("darkblue","gold","darkred","darkgreen")) +
    scale_x_continuous("Mean annual temp. (Celsius)") +
    scale_y_continuous("% base frequency", breaks = c(0,0.5,1)) +
    facet_grid(genom_pos~., scales = "free") +
    theme_minimal(base_size=10) +
    theme(strip.text.y.right = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(color = "black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black")))

(plot8  <- ggplot(CCG8_joined_data, aes(x = BIO1, y = value_ratio, color = variable, linetype = site)) +
    geom_point(size = 1) +
    stat_summary(fun.data = mean_se, geom = "line") +
    scale_color_manual("Base", values = c("darkblue","gold","darkred","darkgreen")) +
    scale_x_continuous("Mean annual temp. (Celsius)") +
    scale_y_continuous("% base frequency", breaks = c(0,0.5,1)) +
    facet_grid(genom_pos~., scales = "free") +
    theme_minimal(base_size=10) +
    theme(strip.text.y.right = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(color = "black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black")))

(plot9  <- ggplot(CCG9_joined_data, aes(x = BIO1, y = value_ratio, color = variable, linetype = site)) +
    geom_point(size = 1) +
    stat_summary(fun.data = mean_se, geom = "line") +
    scale_color_manual("Base", values = c("darkblue","gold","darkred","darkgreen")) +
    scale_x_continuous("Mean annual temp. (Celsius)") +
    scale_y_continuous("% base frequency", breaks = c(0,0.5,1)) +
    facet_grid(genom_pos~., scales = "free") +
    theme_minimal(base_size=10) +
    theme(strip.text.y.right = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(color = "black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black")))


(plot10  <- ggplot(FRQ_joined_data, aes(x = BIO1, y = value_ratio, color = variable, linetype = site)) +
    geom_point(size = 1) +
    stat_summary(fun.data = mean_se, geom = "line") +
    scale_color_manual("Base", values = c("darkblue","gold","darkred","darkgreen")) +
    scale_x_continuous("Mean annual temp. (Celsius)") +
    scale_y_continuous("% base frequency", breaks = c(0,0.5,1)) +
    facet_grid(genom_pos~., scales = "free") +
    theme_minimal(base_size=10) +
    theme(strip.text.y.right = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(color = "black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black")))


pdf("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/5_Gradient_analysis_Upustulata/WC2_gradient_FET_sig_loci_x_BIO1_sync.pdf", height = 3, width = 4)
plot1
dev.off() 

pdf("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/5_Gradient_analysis_Upustulata/WC_like2_gradient_FET_sig_loci_x_BIO1_sync.pdf", height = 6, width = 4)
plot2
dev.off() 

pdf("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/5_Gradient_analysis_Upustulata/FRHa_gradient_FET_sig_loci_x_BIO1_sync.pdf", height = 10, width = 4)
plot3
dev.off() 

pdf("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/5_Gradient_analysis_Upustulata/FRHb_gradient_FET_sig_loci_x_BIO1_sync.pdf", height = 40, width = 4)
plot4
dev.off() 

pdf("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/5_Gradient_analysis_Upustulata/CCG1_gradient_FET_sig_loci_x_BIO1_sync.pdf", height = 20, width = 4)
plot5
dev.off() 

pdf("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/5_Gradient_analysis_Upustulata/CCG6_gradient_FET_sig_loci_x_BIO1_sync.pdf", height = 20, width = 4)
plot6
dev.off() 

pdf("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/5_Gradient_analysis_Upustulata/CCG7_gradient_FET_sig_loci_x_BIO1_sync.pdf", height = 40, width = 4)
plot7
dev.off() 

pdf("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/5_Gradient_analysis_Upustulata/CCG8_gradient_FET_sig_loci_x_BIO1_sync.pdf", height = 30, width = 4)
plot8
dev.off() 

pdf("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/5_Gradient_analysis_Upustulata/CCG9_gradient_FET_sig_loci_x_BIO1_sync.pdf", height = 60, width = 4)
plot9
dev.off() 

pdf("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/5_Gradient_analysis_Upustulata/FRQ_gradient_FET_sig_loci_x_BIO1_sync.pdf", height = 60, width = 4)
plot10
dev.off() 


```


#### U. phaea analysis


```{r Part 1: finding significantly variable regions across gradient extremes in circadian clock genes}
# install.packages("ggplot2")
# install.packages("emmeans")
# install.packages("dplyr")
# install.packages("plyr")
# install.packages("tidyr")
# install.packages("reshape2")
library(ggplot2)
library(emmeans)
library(dplyr)
library(plyr)
library(tidyr)
library(reshape2)

# load FET (Fisher's exact test) significances for all circadian SNPs
Fst_scores_snps <- read.table("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/6_Umbilicaria_phaea_analyses/Uph_circadian_alleles4.fet", sep = "\t", header = F, check.names = F)

# clean up file for R usage
names(Fst_scores_snps) <- Fst_scores_snps[1,]

# rename columns that will be used
colnames(Fst_scores_snps)[1] <- "chr"
colnames(Fst_scores_snps)[2] <- "pos"
colnames(Fst_scores_snps)[5] <- "cov"

names(Fst_scores_snps) <- gsub("=.*$", "", names(Fst_scores_snps))


# replace colons with dashes to match bioclim data
names(Fst_scores_snps) <- gsub("X", "", names(Fst_scores_snps))
names(Fst_scores_snps) <- gsub("\\:", "-", names(Fst_scores_snps))

# remove extra columns
Fst_scores_snps1 <- Fst_scores_snps[,-c(3:4)]

# save the col.names for later
col_names <- colnames(Fst_scores_snps1)

# clean up data columns
Fst_scores_snps1_cleaned <- as.data.frame(lapply(Fst_scores_snps1, function(x) gsub("^.*=","",x)))
names(Fst_scores_snps1_cleaned) <- col_names

# melt the data frame
melted <- melt(Fst_scores_snps1_cleaned, id = c("chr","pos","cov"), variable.name = "Pairwise2", value.name = "FET")

# merge with Pairwise labels from  Fst per-gene data set (we are only using the labels, so it doesn't matter that only FRQ, WC1 and WC2 are in this data set)
Fst_scores_genes_circadian <- read.csv("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/6_Umbilicaria_phaea_analyses/population_contrast_labels.csv", header = T)
Pairwise_labels <- Fst_scores_genes_circadian[,c(1:5)]
merged <- merge(melted, Pairwise_labels, by ="Pairwise2")
merged$position <- as.factor(merged$pos) 
merged$pos_num <- as.numeric(as.character(merged$pos))

merged_unique <- merged %>% group_by(pos,chr) %>% slice(which.max(FET))

# get a list of highest values for each position and gene (should all be mostly in Across) which are above -log(0.05)
merged_unique_sig <- merged %>% group_by(pos,chr) %>% slice(which.max(FET)) %>% slice(which(FET>2.996))

merged_unique_sig1 <- merged_unique_sig[merged_unique_sig$Elevation_Binary == "Across",]

#merged_unique_sig2 <- merged_unique_sig[merged_unique_sig$Pairwise2 == "1-4" | merged_unique_sig$Pairwise2 == "1-11" | merged_unique_sig$Pairwise2 == "5-4" | merged_unique_sig$Pairwise2 == "5-11",]

# write this file 
write.table(merged, file = "~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/6_Umbilicaria_phaea_analyses/SNPs_circadian_unique_all4.fet", sep = "\t", na = "", quote = F, row.names = F, col.names = T)
write.table(merged_unique_sig1, file = "~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/6_Umbilicaria_phaea_analyses/SNPs_circadian_significant4.fet", sep = "\t", na = "", quote = F, row.names = F, col.names = T)

```

Before continuing, we can make our FET graphs:

```{r Part 1b: FET figures}

subset_WC1 <- merged[merged$chr == "scaffold_3" & merged$pos_num_num > 1785861 & merged$pos_num_num < 1789495,]
subset_WC1$region <- "WC-1"
subset_WC2 <- merged[merged$chr == "scaffold_9" & merged$pos_num > 1270142 & merged$pos_num < 1271703,]
subset_WC2$region <- "WC-2"

subset_FRHa <- merged[merged$chr == "scaffold_15" & merged$pos_num > 807734 & merged$pos_num < 811227,]
subset_FRHa$region <- "FRHa"
subset_FRHb <- merged[merged$chr == "scaffold_12" & merged$pos_num > 448168 & merged$pos_num < 452096,]
subset_FRHb$region <- "FRHb"

subset_CK1a <- merged[merged$chr == "scaffold_2" & merged$pos_num > 1364622 & merged$pos_num < 1366037,]
subset_CK1a$region <- "CK1a"
subset_CK1b <- merged[merged$chr == "scaffold_20" & merged$pos_num > 110017 & merged$pos_num < 111112,]
subset_CK1b$region <- "CK1b"
subset_CK1c <- merged[merged$chr == "scaffold_15" & merged$pos_num > 309342 & merged$pos_num < 310926,]
subset_CK1c$region <- "CK1c"

subset_CCG1 <- merged[merged$chr == "scaffold_19" & merged$pos_num > 801922 & merged$pos_num < 802328,]
subset_CCG1$region <- "CCG1"
subset_CCG7 <- merged[merged$chr == "scaffold_13" & merged$pos_num > 351947 & merged$pos_num < 353378,]
subset_CCG7$region <- "CCG7"
subset_CCG8 <- merged[merged$chr == "scaffold_17" & merged$pos_num > 793512 & merged$pos_num < 795138,]
subset_CCG8$region <- "CCG8"
subset_CCG9 <- merged[merged$chr == "scaffold_18" & merged$pos_num > 1025947 & merged$pos_num < 1028342,]
subset_CCG9$region <- "CCG9"

subset_all <- rbind(subset_WC2,subset_FRHa,subset_FRHb,subset_CK1a,subset_CK1b,subset_CK1c,subset_CCG1,subset_CCG7,subset_CCG8,subset_CCG9)
subset_all$region <- as.factor(subset_all$region)

subset_all1 <- subset_all
subset_all1$FET <- as.numeric(subset_all1$FET)

# NOTE: WC1, WC2, and FRQ are in reverse relative to genome positions
(plot1  <- ggplot(subset_all1, aes(x = pos_num, y = FET, col = Elevation_Binary)) +
    #geom_point() +
    stat_summary(geom = "line") +
    scale_x_continuous("SNP position") +
    scale_y_continuous("Fisher's exact test scores (Fst)") +
    scale_color_discrete("Zones") +
    geom_hline(yintercept = -log(0.05), col = "black", linetype = "dashed") +
    theme_minimal() +
    facet_grid(.~region, scales = "free") +
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = T),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black")))

ggsave(plot1, filename = "~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/6_Umbilicaria_phaea_analyses/circadian_FET_SNPs.pdf", 
       height = 3, width = 15, 
       device = cairo_pdf)


bioclim <- read.csv("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/6_Umbilicaria_phaea_analyses/bioclim_11pops_raw.csv", header = T, check.names = F)

(plot1  <- ggplot(bioclim, aes(x = site, y = Elevation, col = BIO1)) +
    geom_point(size = 3) +
    stat_summary(geom = "line") +
    scale_x_discrete("Site") +
    scale_y_continuous("Elevation (m)") +
    scale_color_continuous("BIO1") +
    theme_minimal(base_size=12) +
    #facet_grid(.~region, scales = "free") +
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = T),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black")))

library("RColorBrewer")
myPalette <- colorRampPalette(brewer.pal(11, "Spectral"))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(1, 1600))

(plot2  <- ggplot(bioclim, aes(x = site, y = Elevation, col = BIO12)) +
    geom_point(size = 3) +
    stat_summary(geom = "line") +
    scale_x_discrete("Site") +
    scale_y_continuous("Elevation (m)") +
    scale_color_continuous("BIO12 (mm)", ) +
    theme_minimal(base_size=12) +
    #facet_grid(.~region, scales = "free") +
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = T),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black")) +sc)

ggsave(plot1, filename = "~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/6_Umbilicaria_phaea_analyses/Elevation_gradients.pdf", 
       height = 4, width = 4, 
       device = cairo_pdf)

ggsave(plot2, filename = "~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/6_Umbilicaria_phaea_analyses/Elevation_gradients_precip.pdf", 
       height = 4, width = 4, 
       device = cairo_pdf)

```

Next, we use the extreme SNPs to pull out the relevant allele frequency data from the sync file (in this case, a subset of the sync file with only the circadian loci of interest).

```{r Part 2: creating allele ratio file from sync file}
# install.packages("ggplot2")
# install.packages("emmeans")
# install.packages("dplyr")
# install.packages("plyr")
# install.packages("tidyr")
# install.packages("reshape2")
library(ggplot2)
library(emmeans)
library(dplyr)
library(plyr)
library(tidyr)
library(reshape2)

# load circadian sync file and unique, significant circadian SNPs
sync <- read.table("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/6_Umbilicaria_phaea_analyses/Uph_circadian_alleles4.sync", sep = "\t", header = F, check.names = F)
#unique <- read.table("SNPs_circadian_significant2.fet", sep = "\t", header = T, check.names = F)
unique <- read.table("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/6_Umbilicaria_phaea_analyses/SNPs_circadian_unique_all4.fet", sep = "\t", header = T, check.names = F)

unique_subset <- unique[,c(2:4)]

# label columns
colnames(sync) <- c("chr","pos","rc","MJ_1","MJ_2","MJ_3","MJ_4","MJ_5","MJ_6","MJ_7","SN_1","SN_2","SN_3","SN_4")

# filter based on positions in unique df
sync_significant <- sync %>%
      filter(pos %in% unique$pos)

# now we add in the coverage information
sync_significant1 <- merge(sync_significant, unique_subset, by =c("chr","pos"))

# drop sync now since it takes up a lot of memory
rm(sync)
gc()

# melt data frame
sync_significant2 <- melt(sync_significant1, id = c("chr","pos","cov","rc"))

# separate value and value columns
mdata <- sync_significant2 %>% 
      separate(variable, sep = "_", into = c("site","pop")) %>% 
      separate(value, sep = ":", into = c("A","T","C","G","N","del"))

# re-melt data, and set variable types properly
mdata1 <- melt(mdata, id=c("chr","pos","cov","rc","site","pop")) 
mdata1$pop <- as.numeric(mdata1$pop)
#mdata1$pos <- as.factor(mdata1$pos)
mdata1$cov <- as.numeric(mdata1$cov)
mdata1$value <- as.numeric(mdata1$value)
mdata1$pos <- as.numeric(mdata1$pos)

# create an actual ratio (out of the 30 subsampled populations)
mdata1$value_ratio <- mdata1$value / mdata1$cov

# remove N and del columns, which are not necessary
mdata1.1 <- mdata1[mdata1$variable != "N" & mdata1$variable != "del",]

# subset all individual gene data sets
subset_WC1 <- mdata1.1[mdata1.1$chr == "scaffold_3" & mdata1.1$pos > 1785861 & mdata1.1$pos < 1789495,]
subset_WC1$region <- "WC-1"
subset_WC2 <- mdata1.1[mdata1.1$chr == "scaffold_9" & mdata1.1$pos > 1270142 & mdata1.1$pos < 1271703,]
subset_WC2$region <- "WC-2"

subset_FRHa <- mdata1.1[mdata1.1$chr == "scaffold_15" & mdata1.1$pos > 807734 & mdata1.1$pos < 811227,]
subset_FRHa$region <- "FRHa"
subset_FRHb <- mdata1.1[mdata1.1$chr == "scaffold_12" & mdata1.1$pos > 448168 & mdata1.1$pos < 452096,]
subset_FRHb$region <- "FRHb"

subset_CK1a <- mdata1.1[mdata1.1$chr == "scaffold_2" & mdata1.1$pos > 1364622 & mdata1.1$pos < 1366037,]
subset_CK1a$region <- "CK1a"
subset_CK1b <- mdata1.1[mdata1.1$chr == "scaffold_20" & mdata1.1$pos > 110017 & mdata1.1$pos < 111112,]
subset_CK1b$region <- "CK1b"
subset_CK1c <- mdata1.1[mdata1.1$chr == "scaffold_15" & mdata1.1$pos > 309342 & mdata1.1$pos < 310926,]
subset_CK1c$region <- "CK1c"

subset_CCG1 <- mdata1.1[mdata1.1$chr == "scaffold_19" & mdata1.1$pos > 801922 & mdata1.1$pos < 802328,]
subset_CCG1$region <- "CCG1"
subset_CCG7 <- mdata1.1[mdata1.1$chr == "scaffold_13" & mdata1.1$pos > 351947 & mdata1.1$pos < 353378,]
subset_CCG7$region <- "CCG7"
subset_CCG8 <- mdata1.1[mdata1.1$chr == "scaffold_17" & mdata1.1$pos > 793512 & mdata1.1$pos < 795138,]
subset_CCG8$region <- "CCG8"
subset_CCG9 <- mdata1.1[mdata1.1$chr == "scaffold_18" & mdata1.1$pos > 1025947 & mdata1.1$pos < 1028342,]
subset_CCG9$region <- "CCG9"

# re-join subsets with their new region labels
subset_all <- rbind(subset_WC2,subset_FRHa,subset_FRHb,subset_CK1a,subset_CK1b,subset_CK1c,subset_CCG1,subset_CCG7,subset_CCG8,subset_CCG9)
subset_all$region <- as.factor(subset_all$region)

# write this file 
write.table(subset_all, file = "~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/6_Umbilicaria_phaea_analyses/circadian_alleles_sync_labeled_melted4.txt", sep = "\t", na = "", quote = F, row.names = F, col.names = T)

```


```{r Part 2b: creating allele ratio file from sync file -- FET sig. results}
# install.packages("ggplot2")
# install.packages("emmeans")
# install.packages("dplyr")
# install.packages("plyr")
# install.packages("tidyr")
# install.packages("reshape2")
library(ggplot2)
library(emmeans)
library(dplyr)
library(plyr)
library(tidyr)
library(reshape2)

# load circadian sync file and unique, significant circadian SNPs
sync <- read.table("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/6_Umbilicaria_phaea_analyses/Uph_circadian_alleles4.sync", sep = "\t", header = F, check.names = F)
#unique <- read.table("SNPs_circadian_significant2.fet", sep = "\t", header = T, check.names = F)
unique_sig <- read.table("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/6_Umbilicaria_phaea_analyses/SNPs_circadian_significant4.fet", sep = "\t", header = T, check.names = F)
#unique_sig <- read.table("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/6_Umbilicaria_phaea_analyses/SNPs_circadian_significant4b.fet", sep = "\t", header = T, check.names = F)

unique_subset <- unique_sig[,c(2:4)]

# label columns
colnames(sync) <- c("chr","pos","rc","MJ_1","MJ_2","MJ_3","MJ_4","MJ_5","MJ_6","MJ_7","SN_1","SN_2","SN_3","SN_4")

# filter based on positions in unique df
sync_significant <- sync %>%
      filter(pos %in% unique_sig$pos)

# now we add in the coverage information
sync_significant1 <- merge(sync_significant, unique_subset, by =c("chr","pos"))

# drop sync now since it takes up a lot of memory
rm(sync)
gc()

# melt data frame
sync_significant2 <- melt(sync_significant1, id = c("chr","pos","cov","rc"))

# separate value and value columns
mdata <- sync_significant2 %>% 
      separate(variable, sep = "_", into = c("site","pop")) %>% 
      separate(value, sep = ":", into = c("A","T","C","G","N","del"))

# re-melt data, and set variable types properly
mdata1 <- melt(mdata, id=c("chr","pos","cov","rc","site","pop")) 
mdata1$pop <- as.numeric(mdata1$pop)
#mdata1$pos <- as.factor(mdata1$pos)
mdata1$cov <- as.numeric(mdata1$cov)
mdata1$value <- as.numeric(mdata1$value)
mdata1$pos <- as.numeric(mdata1$pos)

# create an actual ratio (out of the 30 subsampled populations)
mdata1$value_ratio <- mdata1$value / mdata1$cov

# remove N and del columns, which are not necessary
mdata1.1 <- mdata1[mdata1$variable != "N" & mdata1$variable != "del",]

# subset all individual gene data sets
subset_WC1 <- mdata1.1[mdata1.1$chr == "scaffold_3" & mdata1.1$pos > 1785861 & mdata1.1$pos < 1789495,]
subset_WC1$region <- "WC-1"
subset_WC2 <- mdata1.1[mdata1.1$chr == "scaffold_9" & mdata1.1$pos > 1270142 & mdata1.1$pos < 1271703,]
subset_WC2$region <- "WC-2"

subset_FRHa <- mdata1.1[mdata1.1$chr == "scaffold_15" & mdata1.1$pos > 807734 & mdata1.1$pos < 811227,]
subset_FRHa$region <- "FRHa"
subset_FRHb <- mdata1.1[mdata1.1$chr == "scaffold_12" & mdata1.1$pos > 448168 & mdata1.1$pos < 452096,]
subset_FRHb$region <- "FRHb"

subset_CK1a <- mdata1.1[mdata1.1$chr == "scaffold_2" & mdata1.1$pos > 1364622 & mdata1.1$pos < 1366037,]
subset_CK1a$region <- "CK1a"
subset_CK1b <- mdata1.1[mdata1.1$chr == "scaffold_20" & mdata1.1$pos > 110017 & mdata1.1$pos < 111112,]
subset_CK1b$region <- "CK1b"
subset_CK1c <- mdata1.1[mdata1.1$chr == "scaffold_15" & mdata1.1$pos > 309342 & mdata1.1$pos < 310926,]
subset_CK1c$region <- "CK1c"

subset_CCG1 <- mdata1.1[mdata1.1$chr == "scaffold_19" & mdata1.1$pos > 801922 & mdata1.1$pos < 802328,]
subset_CCG1$region <- "CCG1"
subset_CCG7 <- mdata1.1[mdata1.1$chr == "scaffold_13" & mdata1.1$pos > 351947 & mdata1.1$pos < 353378,]
subset_CCG7$region <- "CCG7"
subset_CCG8 <- mdata1.1[mdata1.1$chr == "scaffold_17" & mdata1.1$pos > 793512 & mdata1.1$pos < 795138,]
subset_CCG8$region <- "CCG8"
subset_CCG9 <- mdata1.1[mdata1.1$chr == "scaffold_18" & mdata1.1$pos > 1025947 & mdata1.1$pos < 1028342,]
subset_CCG9$region <- "CCG9"

# re-join subsets with their new region labels
subset_all <- rbind(subset_WC2,subset_FRHa,subset_FRHb,subset_CK1a,subset_CK1b,subset_CK1c,subset_CCG1,subset_CCG7,subset_CCG8,subset_CCG9)
subset_all$region <- as.factor(subset_all$region)

# write this file 
write.table(subset_all, file = "~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/6_Umbilicaria_phaea_analyses/circadian_alleles_sync_labeled_melted4_sig.txt", sep = "\t", na = "", quote = F, row.names = F, col.names = T)

```

Next, we can analyze these data: identify which loci have trends along the gradient, merge the data with the environmental/climatic data, and finally graph it.

```{r Part 3: allele ratio gradient analysis -- extracting signifant Spearman corr. loci}
# install.packages("ggplot2")
# install.packages("emmeans")
# install.packages("dplyr")
# install.packages("plyr")
# install.packages("tidyr")
# install.packages("reshape2")
# install.packages("emmeans")
# install.packages("multcompView")
# install.packages("corrr")
library(ggplot2)
library(emmeans)
library(dplyr)
library(plyr)
library(tidyr)
library(reshape2)
library(emmeans)
library(multcompView)

# clear the list! 
rm(list=ls())
gc()

# load the labeled, melted circadian sync file
data <- read.table("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/6_Umbilicaria_phaea_analyses/circadian_alleles_sync_labeled_melted4.txt", sep = "\t", header = T, check.names = F)

# set variables properly
data$region <- as.factor(data$region)
data$site <- as.factor(data$site)
data$pos <- as.factor(data$pos)
data$pop <- as.numeric(data$pop)
data$variable <- as.factor(data$variable)

# briefly calculate how many significant SNPs we have:
data1 <- data[,c(2,10)]
unique <- data1 %>% distinct()
sig_count <- unique %>% dplyr::count(unique$region) 

## generating for loop to cycle through all correlations of allele ratios across pop for each site

## subset each nucleotide and run each separately
# A
data_A <- data[data$variable == "A",]
data_A$genom_pos <- paste(data_A$chr,data_A$pos)
data_A$genom_pos <- as.factor(data_A$genom_pos)

results_A <- data.frame(row.names = 1:12720)

for(i in 1:length(levels(data_A$genom_pos))){
  for(j in 1:length(levels(data_A$site))){
    subset_focus <- subset(data_A, genom_pos == levels(data_A$genom_pos)[i] & site == levels(data_A$site)[j])
    corr_focus <- suppressWarnings(cor(subset_focus$pop,subset_focus$value_ratio, method = "spearman"))
      if(is.na(corr_focus)){
        results_A[1000*(j-1)+i,1] <- levels(data_A$genom_pos)[i]
        results_A[1000*(j-1)+i,2] <- levels(data_A$site)[j]
        results_A[1000*(j-1)+i,2] <- levels(data_A$region)[j]
        results_A[1000*(j-1)+i,3] <- NA
        results_A[1000*(j-1)+i,4] <- NA
      }
    else{
        p_val_focus <- suppressWarnings(cor.test(subset_focus$pop,subset_focus$value_ratio, method = "spearman"))[3][[1]]
        results_A[1000*(j-1)+i,1] <- levels(data_A$genom_pos)[i]
        results_A[1000*(j-1)+i,2] <- levels(data_A$site)[j]
        results_A[1000*(j-1)+i,2] <- levels(data_A$region)[j]
        results_A[1000*(j-1)+i,3] <- corr_focus
        results_A[1000*(j-1)+i,4] <- p_val_focus
    }
  colnames(results_A) <- c("genom_pos","site","spearman_corr","p_val_corr")
  }
}

# T
data_T <- data[data$variable == "T",]
data_T$genom_pos <- paste(data_T$chr,data_T$pos)
data_T$genom_pos <- as.factor(data_T$genom_pos)

results_T <- data.frame(row.names = 1:432)

for(i in 1:length(levels(data_T$genom_pos))){
  for(j in 1:length(levels(data_T$site))){
    subset_focus <- subset(data_T, genom_pos == levels(data_T$genom_pos)[i] & site == levels(data_T$site)[j])
    corr_focus <- suppressWarnings(cor(subset_focus$pop,subset_focus$value_ratio, method = "spearman"))
      if(is.na(corr_focus)){
        results_T[144*(j-1)+i,1] <- levels(data_T$genom_pos)[i]
        results_T[144*(j-1)+i,2] <- levels(data_T$site)[j]
        results_T[144*(j-1)+i,2] <- levels(data_T$region)[j]
        results_T[144*(j-1)+i,3] <- NA
        results_T[144*(j-1)+i,4] <- NA
      }
    else{
        p_val_focus <- suppressWarnings(cor.test(subset_focus$pop,subset_focus$value_ratio, method = "spearman"))[3][[1]]
        results_T[144*(j-1)+i,1] <- levels(data_T$genom_pos)[i]
        results_T[144*(j-1)+i,2] <- levels(data_T$site)[j]
        results_T[144*(j-1)+i,2] <- levels(data_T$region)[j]
        results_T[144*(j-1)+i,3] <- corr_focus
        results_T[144*(j-1)+i,4] <- p_val_focus
    }
  colnames(results_T) <- c("genom_pos","site","spearman_corr","p_val_corr")
  }
}

# C
data_C <- data[data$variable == "C",]
data_C$genom_pos <- paste(data_C$chr,data_C$pos)
data_C$genom_pos <- as.factor(data_C$genom_pos)

results_C <- data.frame(row.names = 1:432)

for(i in 1:length(levels(data_C$genom_pos))){
  for(j in 1:length(levels(data_C$site))){
    subset_focus <- subset(data_C, genom_pos == levels(data_C$genom_pos)[i] & site == levels(data_C$site)[j])
    corr_focus <- suppressWarnings(cor(subset_focus$pop,subset_focus$value_ratio, method = "spearman"))
      if(is.na(corr_focus)){
        results_C[144*(j-1)+i,1] <- levels(data_C$genom_pos)[i]
        results_C[144*(j-1)+i,2] <- levels(data_C$site)[j]
        results_C[144*(j-1)+i,2] <- levels(data_C$region)[j]
        results_C[144*(j-1)+i,3] <- NA
        results_C[144*(j-1)+i,4] <- NA
      }
    else{
        p_val_focus <- suppressWarnings(cor.test(subset_focus$pop,subset_focus$value_ratio, method = "spearman"))[3][[1]]
        results_C[144*(j-1)+i,1] <- levels(data_C$genom_pos)[i]
        results_C[144*(j-1)+i,2] <- levels(data_C$site)[j]
        results_C[144*(j-1)+i,2] <- levels(data_C$region)[j]
        results_C[144*(j-1)+i,3] <- corr_focus
        results_C[144*(j-1)+i,4] <- p_val_focus
    }
  colnames(results_C) <- c("genom_pos","site","spearman_corr","p_val_corr")
  }
}

# G
data_G <- data[data$variable == "G",]
data_G$genom_pos <- paste(data_G$chr,data_T$pos)
data_G$genom_pos <- as.factor(data_G$genom_pos)

results_G <- data.frame(row.names = 1:432)

for(i in 1:length(levels(data_G$genom_pos))){
  for(j in 1:length(levels(data_G$site))){
    subset_focus <- subset(data_G, genom_pos == levels(data_G$genom_pos)[i] & site == levels(data_G$site)[j])
    corr_focus <- suppressWarnings(cor(subset_focus$pop,subset_focus$value_ratio, method = "spearman"))
      if(is.na(corr_focus)){
        results_G[144*(j-1)+i,1] <- levels(data_G$genom_pos)[i]
        results_G[144*(j-1)+i,2] <- levels(data_G$site)[j]
        results_G[144*(j-1)+i,2] <- levels(data_G$region)[j]
        results_G[144*(j-1)+i,3] <- NA
        results_G[144*(j-1)+i,4] <- NA
      }
    else{
        p_val_focus <- suppressWarnings(cor.test(subset_focus$pop,subset_focus$value_ratio, method = "spearman"))[3][[1]]
        results_G[144*(j-1)+i,1] <- levels(data_G$genom_pos)[i]
        results_G[144*(j-1)+i,2] <- levels(data_G$site)[j]
        results_G[144*(j-1)+i,2] <- levels(data_G$region)[j]
        results_G[144*(j-1)+i,3] <- corr_focus
        results_G[144*(j-1)+i,4] <- p_val_focus
    }
  colnames(results_G) <- c("genom_pos","site","spearman_corr","p_val_corr")
  }
}

# remove all NAs
results_A.1 <- results_A[!is.na(results_A$spearman_corr),]
results_T.1 <- results_T[!is.na(results_T$spearman_corr),]
results_C.1 <- results_C[!is.na(results_C$spearman_corr),]
results_G.1 <- results_G[!is.na(results_G$spearman_corr),]

# add back variable information
results_A.1$variable <- "A"
results_T.1$variable <- "T"
results_C.1$variable <- "C"
results_G.1$variable <- "G"

# re-join subsets with their new variable labels
results_all <- rbind(results_A.1,results_T.1,results_C.1,results_G.1)
results_all$variable <- as.factor(results_all$variable)
results_all$site <- as.factor(results_all$site)

#subset only significant results
results_significant <- results_all[results_all$p_val_corr < 0.05,]

# we can then make a vector of all the significant loci for ES1 and IT (ES2 only has 3 populations in the gradient, so there are very few significant hits, presumably due to the low resolution of the 3-pop gradient)
significant_loci <- results_significant$genom_pos

# join all chr and pos values
data$genom_pos <- paste(data$chr,data$pos)

# join with significant loci
data_significant_subset <- data[data$genom_pos %in% significant_loci, ]

# write this file 
write.table(data_significant_subset, file = "~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/6_Umbilicaria_phaea_analyses/circadian_alleles_spearman_sig_0_05_loci4.txt", sep = "\t", na = "", quote = F, row.names = F, col.names = T)

```

Next, we plot the results with their respective bioclim data. 

```{r Part 4: allele ratio gradient analysis -- Spearman sig. results}
# install.packages("ggplot2")
# install.packages("emmeans")
# install.packages("dplyr")
# install.packages("plyr")
# install.packages("tidyr")
# install.packages("reshape2")
# install.packages("emmeans")
# install.packages("multcompView")
# install.packages("corrr")
library(ggplot2)
library(emmeans)
library(dplyr)
library(plyr)
library(tidyr)
library(reshape2)
library(emmeans)
library(multcompView)

# clear the list! 
rm(list=ls())
gc()

# load the file with significant loci and the bioclim file
data <- read.table("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/6_Umbilicaria_phaea_analyses/circadian_alleles_spearman_sig_0_05_loci4.txt", sep = "\t", header = T, check.names = F)
bioclim <- read.csv("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/6_Umbilicaria_phaea_analyses/bioclim_11pops_raw.csv", header = T, check.names = F)

joined_data <- merge(data,bioclim, by = c("site","pop"))

joined_data$genom_pos <- paste(joined_data$chr,joined_data$pos)
joined_data$genom_pos <- as.factor(joined_data$genom_pos)

# looking at loci significant in any of the sites (rather than all three)

# need to manually set the axis range for each of the y-axes for setting an additional x axis:
# ylim.prim <- c(0, 1)   # in this example, allele ratio
# ylim.sec <- c(350, 950)    # in this example, elevation
# b <- diff(ylim.prim)/diff(ylim.sec)
# a <- ylim.prim[1] - b*ylim.sec[1] 
### Then, add this to the plot:
# geom_line(aes(y = a + BIO12*b), color = "cornflowerblue", linetype = "dashed") + # to add the data 
# scale_y_continuous("% base frequency", breaks = c(0,0.5,1), sec.axis = sec_axis(~ (. - a)/b, name = "Mean annual precip. (mm)")) + # modify the scale like this

WC1_joined_data <- joined_data[joined_data$region == "WC-1",]
WC2_joined_data <- joined_data[joined_data$region == "WC-2",]
FRHa_joined_data <- joined_data[joined_data$region == "FRHa",]
FRHb_joined_data <- joined_data[joined_data$region == "FRHb",]
CK1a_joined_data <- joined_data[joined_data$region == "CK1a",]
CK1b_joined_data <- joined_data[joined_data$region == "CK1b",]
CK1c_joined_data <- joined_data[joined_data$region == "CK1c",]
CCG1_joined_data <- joined_data[joined_data$region == "CCG1",]
CCG7_joined_data <- joined_data[joined_data$region == "CCG7",]
CCG8_joined_data <- joined_data[joined_data$region == "CCG8",]
CCG9_joined_data <- joined_data[joined_data$region == "CCG9",]


(plot2  <- ggplot(WC2_joined_data, aes(x = BIO1, y = value_ratio, color = variable, linetype = site)) +
    geom_point(size = 1) +
    stat_summary(fun.data = mean_se, geom = "line") +
    scale_color_manual("Base", values = c("darkblue","gold","darkred","darkgreen")) +
    scale_x_continuous("Mean annual temp. (Celsius)") +
    scale_y_continuous("% base frequency", breaks = c(0,0.5,1)) +
    facet_grid(genom_pos~., scales = "free") +
    theme_minimal(base_size=10) +
    theme(strip.text.y.right = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(color = "black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black")))

(plot3  <- ggplot(FRHa_joined_data, aes(x = BIO1, y = value_ratio, color = variable, linetype = site)) +
    geom_point(size = 1) +
    stat_summary(fun.data = mean_se, geom = "line") +
    scale_color_manual("Base", values = c("darkblue","gold","darkred","darkgreen")) +
    scale_x_continuous("Mean annual temp. (Celsius)") +
    scale_y_continuous("% base frequency", breaks = c(0,0.5,1)) +
    facet_grid(genom_pos~., scales = "free") +
    theme_minimal(base_size=10) +
    theme(strip.text.y.right = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(color = "black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black")))

(plot4  <- ggplot(FRHb_joined_data, aes(x = BIO1, y = value_ratio, color = variable, linetype = site)) +
    geom_point(size = 1) +
    stat_summary(fun.data = mean_se, geom = "line") +
    scale_color_manual("Base", values = c("darkblue","gold","darkred","darkgreen")) +
    scale_x_continuous("Mean annual temp. (Celsius)") +
    scale_y_continuous("% base frequency", breaks = c(0,0.5,1)) +
    facet_grid(genom_pos~., scales = "free") +
    theme_minimal(base_size=10) +
    theme(strip.text.y.right = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(color = "black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black")))

FRHb_joined_data_linear <- FRHb_joined_data[FRHb_joined_data$genom_pos == "scaffold_12 451039",]
FRHb_joined_data_bell <- FRHb_joined_data[FRHb_joined_data$genom_pos == "scaffold_12 449105",]
FRHb_joined_data_negative <- FRHb_joined_data[FRHb_joined_data$genom_pos == "scaffold_12 448180",]

FRHb_joined_data_focus <- rbind(FRHb_joined_data_linear,FRHb_joined_data_bell,FRHb_joined_data_negative)


(plot4a  <- ggplot(FRHb_joined_data_focus, aes(x = BIO1, y = value_ratio, color = variable, linetype = site)) +
    geom_point(size = 1) +
    stat_summary(fun.data = mean_se, geom = "line") +
    scale_color_manual("Base", values = c("darkblue","gold","darkred","darkgreen")) +
    scale_x_continuous("Mean annual temp. (Celsius)") +
    scale_y_continuous("% base frequency", breaks = c(0,0.5,1)) +
    facet_grid(.~genom_pos, scales = "free") +
    theme_minimal(base_size=10) +
    theme(strip.text.y.right = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(color = "black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black")))





(plot5  <- ggplot(CK1a_joined_data, aes(x = BIO1, y = value_ratio, color = variable, linetype = site)) +
    geom_point(size = 1) +
    stat_summary(fun.data = mean_se, geom = "line") +
    scale_color_manual("Base", values = c("darkblue","gold","darkred","darkgreen")) +
    scale_x_continuous("Mean annual temp. (Celsius)") +
    scale_y_continuous("% base frequency", breaks = c(0,0.5,1)) +
    facet_grid(genom_pos~., scales = "free") +
    theme_minimal(base_size=10) +
    theme(strip.text.y.right = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(color = "black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black")))

(plot6  <- ggplot(CK1b_joined_data, aes(x = BIO1, y = value_ratio, color = variable, linetype = site)) +
    geom_point(size = 1) +
    stat_summary(fun.data = mean_se, geom = "line") +
    scale_color_manual("Base", values = c("darkblue","gold","darkred","darkgreen")) +
    scale_x_continuous("Mean annual temp. (Celsius)") +
    scale_y_continuous("% base frequency", breaks = c(0,0.5,1)) +
    facet_grid(genom_pos~., scales = "free") +
    theme_minimal(base_size=10) +
    theme(strip.text.y.right = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(color = "black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black")))

(plot7  <- ggplot(CK1c_joined_data, aes(x = BIO1, y = value_ratio, color = variable, linetype = site)) +
    geom_point(size = 1) +
    stat_summary(fun.data = mean_se, geom = "line") +
    scale_color_manual("Base", values = c("darkblue","gold","darkred","darkgreen")) +
    scale_x_continuous("Mean annual temp. (Celsius)") +
    scale_y_continuous("% base frequency", breaks = c(0,0.5,1)) +
    facet_grid(genom_pos~., scales = "free") +
    theme_minimal(base_size=10) +
    theme(strip.text.y.right = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(color = "black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black")))

(plot8  <- ggplot(CCG1_joined_data, aes(x = BIO1, y = value_ratio, color = variable, linetype = site)) +
    geom_point(size = 1) +
    stat_summary(fun.data = mean_se, geom = "line") +
    scale_color_manual("Base", values = c("darkblue","gold","darkred","darkgreen")) +
    scale_x_continuous("Mean annual temp. (Celsius)") +
    scale_y_continuous("% base frequency", breaks = c(0,0.5,1)) +
    facet_grid(genom_pos~., scales = "free") +
    theme_minimal(base_size=10) +
    theme(strip.text.y.right = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(color = "black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black")))

(plot9  <- ggplot(CCG7_joined_data, aes(x = BIO1, y = value_ratio, color = variable, linetype = site)) +
    geom_point(size = 1) +
    stat_summary(fun.data = mean_se, geom = "line") +
    scale_color_manual("Base", values = c("darkblue","gold","darkred","darkgreen")) +
    scale_x_continuous("Mean annual temp. (Celsius)") +
    scale_y_continuous("% base frequency", breaks = c(0,0.5,1)) +
    facet_grid(genom_pos~., scales = "free") +
    theme_minimal(base_size=10) +
    theme(strip.text.y.right = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(color = "black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black")))

(plot10  <- ggplot(CCG8_joined_data, aes(x = BIO1, y = value_ratio, color = variable, linetype = site)) +
    geom_point(size = 1) +
    stat_summary(fun.data = mean_se, geom = "line") +
    scale_color_manual("Base", values = c("darkblue","gold","darkred","darkgreen")) +
    scale_x_continuous("Mean annual temp. (Celsius)") +
    scale_y_continuous("% base frequency", breaks = c(0,0.5,1)) +
    facet_grid(genom_pos~., scales = "free") +
    theme_minimal(base_size=10) +
    theme(strip.text.y.right = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(color = "black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black")))

(plot11  <- ggplot(CCG9_joined_data, aes(x = BIO1, y = value_ratio, color = variable, linetype = site)) +
    geom_point(size = 1) +
    stat_summary(fun.data = mean_se, geom = "line") +
    scale_color_manual("Base", values = c("darkblue","gold","darkred","darkgreen")) +
    scale_x_continuous("Mean annual temp. (Celsius)") +
    scale_y_continuous("% base frequency", breaks = c(0,0.5,1)) +
    facet_grid(genom_pos~., scales = "free") +
    theme_minimal(base_size=10) +
    theme(strip.text.y.right = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(color = "black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black")))

pdf("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/6_Umbilicaria_phaea_analyses/WC2_gradient_spearman_sig_loci_x_BIO1_sync.pdf", height = 120, width = 4)
plot2
dev.off() 

pdf("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/6_Umbilicaria_phaea_analyses/FRHa_gradient_spearman_sig_loci_x_BIO1_sync.pdf", height = 180, width = 4)
plot3
dev.off() 

pdf("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/6_Umbilicaria_phaea_analyses/FRHb_gradient_spearman_sig_loci_x_BIO1_sync.pdf", height = 180, width = 4)
plot4
dev.off() 


pdf("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/6_Umbilicaria_phaea_analyses/FRHb_gradient_spearman_sig_loci_x_BIO1_sync_Fig_focus.pdf", height = 2, width = 6)
plot4a
dev.off() 


pdf("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/6_Umbilicaria_phaea_analyses/CK1a_gradient_spearman_sig_loci_x_BIO1_sync.pdf", height = 100, width = 4)
plot5
dev.off() 

pdf("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/6_Umbilicaria_phaea_analyses/CK1b_gradient_spearman_sig_loci_x_BIO1_sync.pdf", height = 150, width = 4)
plot6
dev.off() 

pdf("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/6_Umbilicaria_phaea_analyses/CK1c_gradient_spearman_sig_loci_x_BIO1_sync.pdf", height = 50, width = 4)
plot7
dev.off() 

pdf("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/6_Umbilicaria_phaea_analyses/CCG1_gradient_spearman_sig_loci_x_BIO1_sync.pdf", height = 40, width = 4)
plot8
dev.off() 

pdf("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/6_Umbilicaria_phaea_analyses/CCG7_gradient_spearman_sig_loci_x_BIO1_sync.pdf", height = 80, width = 4)
plot9
dev.off() 

pdf("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/6_Umbilicaria_phaea_analyses/CCG8_gradient_spearman_sig_loci_x_BIO1_sync.pdf", height = 100, width = 4)
plot10
dev.off() 

pdf("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/6_Umbilicaria_phaea_analyses/CCG9_gradient_spearman_sig_loci_x_BIO1_sync.pdf", height = 150, width = 4)
plot11
dev.off() 

```




```{r Part 4b: allele ratio gradient analysis -- FET sig. results}
# install.packages("ggplot2")
# install.packages("emmeans")
# install.packages("dplyr")
# install.packages("plyr")
# install.packages("tidyr")
# install.packages("reshape2")
# install.packages("emmeans")
# install.packages("multcompView")
# install.packages("corrr")
library(ggplot2)
library(emmeans)
library(dplyr)
library(plyr)
library(tidyr)
library(reshape2)
library(emmeans)
library(multcompView)

# clear the list! 
rm(list=ls())
gc()

# load the file with significant loci and the bioclim file
data <- read.table("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/6_Umbilicaria_phaea_analyses/circadian_alleles_sync_labeled_melted4_sig.txt", sep = "\t", header = T, check.names = F)
#data <- read.table("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/6_Umbilicaria_phaea_analyses/circadian_alleles_sync_labeled_melted4_sig2.txt", sep = "\t", header = T, check.names = F)
bioclim <- read.csv("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/6_Umbilicaria_phaea_analyses/bioclim_11pops_raw.csv", header = T, check.names = F)

joined_data <- merge(data,bioclim, by = c("site","pop"))

joined_data$genom_pos <- paste(joined_data$chr,joined_data$pos)
joined_data$genom_pos <- as.factor(joined_data$genom_pos)

# looking at loci significant in any of the sites (rather than all three)

# need to manually set the axis range for each of the y-axes for setting an additional x axis:
# ylim.prim <- c(0, 1)   # in this example, allele ratio
# ylim.sec <- c(350, 950)    # in this example, elevation
# b <- diff(ylim.prim)/diff(ylim.sec)
# a <- ylim.prim[1] - b*ylim.sec[1] 
### Then, add this to the plot:
# geom_line(aes(y = a + BIO12*b), color = "cornflowerblue", linetype = "dashed") + # to add the data 
# scale_y_continuous("% base frequency", breaks = c(0,0.5,1), sec.axis = sec_axis(~ (. - a)/b, name = "Mean annual precip. (mm)")) + # modify the scale like this

WC1_joined_data <- joined_data[joined_data$region == "WC-1",]
WC2_joined_data <- joined_data[joined_data$region == "WC-2",]
FRHa_joined_data <- joined_data[joined_data$region == "FRHa",]
FRHb_joined_data <- joined_data[joined_data$region == "FRHb",]
CK1a_joined_data <- joined_data[joined_data$region == "CK1a",]
CK1b_joined_data <- joined_data[joined_data$region == "CK1b",]
CK1c_joined_data <- joined_data[joined_data$region == "CK1c",]
CCG1_joined_data <- joined_data[joined_data$region == "CCG1",]
CCG7_joined_data <- joined_data[joined_data$region == "CCG7",]
CCG8_joined_data <- joined_data[joined_data$region == "CCG8",]
CCG9_joined_data <- joined_data[joined_data$region == "CCG9",]


(plot2  <- ggplot(WC2_joined_data, aes(x = BIO1, y = value_ratio, color = variable, linetype = site)) +
    geom_point(size = 1) +
    stat_summary(fun.data = mean_se, geom = "line") +
    scale_color_manual("Base", values = c("darkblue","gold","darkred","darkgreen")) +
    scale_x_continuous("Mean annual temp. (Celsius)") +
    scale_y_continuous("% base frequency", breaks = c(0,0.5,1)) +
    facet_grid(genom_pos~., scales = "free") +
    theme_minimal(base_size=10) +
    theme(strip.text.y.right = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(color = "black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black")))

(plot3  <- ggplot(FRHa_joined_data, aes(x = BIO1, y = value_ratio, color = variable, linetype = site)) +
    geom_point(size = 1) +
    stat_summary(fun.data = mean_se, geom = "line") +
    scale_color_manual("Base", values = c("darkblue","gold","darkred","darkgreen")) +
    scale_x_continuous("Mean annual temp. (Celsius)") +
    scale_y_continuous("% base frequency", breaks = c(0,0.5,1)) +
    facet_grid(genom_pos~., scales = "free") +
    theme_minimal(base_size=10) +
    theme(strip.text.y.right = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(color = "black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black")))

(plot4  <- ggplot(FRHb_joined_data, aes(x = BIO1, y = value_ratio, color = variable, linetype = site)) +
    geom_point(size = 1) +
    stat_summary(fun.data = mean_se, geom = "line") +
    scale_color_manual("Base", values = c("darkblue","gold","darkred","darkgreen")) +
    scale_x_continuous("Mean annual temp. (Celsius)") +
    scale_y_continuous("% base frequency", breaks = c(0,0.5,1)) +
    facet_grid(genom_pos~., scales = "free") +
    theme_minimal(base_size=10) +
    theme(strip.text.y.right = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(color = "black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black")))

(plot5  <- ggplot(CK1a_joined_data, aes(x = BIO1, y = value_ratio, color = variable, linetype = site)) +
    geom_point(size = 1) +
    stat_summary(fun.data = mean_se, geom = "line") +
    scale_color_manual("Base", values = c("darkblue","gold","darkred","darkgreen")) +
    scale_x_continuous("Mean annual temp. (Celsius)") +
    scale_y_continuous("% base frequency", breaks = c(0,0.5,1)) +
    facet_grid(genom_pos~., scales = "free") +
    theme_minimal(base_size=10) +
    theme(strip.text.y.right = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(color = "black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black")))

(plot6  <- ggplot(CK1b_joined_data, aes(x = BIO1, y = value_ratio, color = variable, linetype = site)) +
    geom_point(size = 1) +
    stat_summary(fun.data = mean_se, geom = "line") +
    scale_color_manual("Base", values = c("darkblue","gold","darkred","darkgreen")) +
    scale_x_continuous("Mean annual temp. (Celsius)") +
    scale_y_continuous("% base frequency", breaks = c(0,0.5,1)) +
    facet_grid(genom_pos~., scales = "free") +
    theme_minimal(base_size=10) +
    theme(strip.text.y.right = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(color = "black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black")))

(plot7  <- ggplot(CK1c_joined_data, aes(x = BIO1, y = value_ratio, color = variable, linetype = site)) +
    geom_point(size = 1) +
    stat_summary(fun.data = mean_se, geom = "line") +
    scale_color_manual("Base", values = c("darkblue","gold","darkred","darkgreen")) +
    scale_x_continuous("Mean annual temp. (Celsius)") +
    scale_y_continuous("% base frequency", breaks = c(0,0.5,1)) +
    facet_grid(genom_pos~., scales = "free") +
    theme_minimal(base_size=10) +
    theme(strip.text.y.right = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(color = "black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black")))

(plot8  <- ggplot(CCG1_joined_data, aes(x = BIO1, y = value_ratio, color = variable, linetype = site)) +
    geom_point(size = 1) +
    stat_summary(fun.data = mean_se, geom = "line") +
    scale_color_manual("Base", values = c("darkblue","gold","darkred","darkgreen")) +
    scale_x_continuous("Mean annual temp. (Celsius)") +
    scale_y_continuous("% base frequency", breaks = c(0,0.5,1)) +
    facet_grid(genom_pos~., scales = "free") +
    theme_minimal(base_size=10) +
    theme(strip.text.y.right = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(color = "black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black")))

(plot9  <- ggplot(CCG7_joined_data, aes(x = BIO1, y = value_ratio, color = variable, linetype = site)) +
    geom_point(size = 1) +
    stat_summary(fun.data = mean_se, geom = "line") +
    scale_color_manual("Base", values = c("darkblue","gold","darkred","darkgreen")) +
    scale_x_continuous("Mean annual temp. (Celsius)") +
    scale_y_continuous("% base frequency", breaks = c(0,0.5,1)) +
    facet_grid(genom_pos~., scales = "free") +
    theme_minimal(base_size=10) +
    theme(strip.text.y.right = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(color = "black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black")))

(plot10  <- ggplot(CCG8_joined_data, aes(x = BIO1, y = value_ratio, color = variable, linetype = site)) +
    geom_point(size = 1) +
    stat_summary(fun.data = mean_se, geom = "line") +
    scale_color_manual("Base", values = c("darkblue","gold","darkred","darkgreen")) +
    scale_x_continuous("Mean annual temp. (Celsius)") +
    scale_y_continuous("% base frequency", breaks = c(0,0.5,1)) +
    facet_grid(genom_pos~., scales = "free") +
    theme_minimal(base_size=10) +
    theme(strip.text.y.right = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(color = "black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black")))

(plot11  <- ggplot(CCG9_joined_data, aes(x = BIO1, y = value_ratio, color = variable, linetype = site)) +
    geom_point(size = 1) +
    stat_summary(fun.data = mean_se, geom = "line") +
    scale_color_manual("Base", values = c("darkblue","gold","darkred","darkgreen")) +
    scale_x_continuous("Mean annual temp. (Celsius)") +
    scale_y_continuous("% base frequency", breaks = c(0,0.5,1)) +
    facet_grid(genom_pos~., scales = "free") +
    theme_minimal(base_size=10) +
    theme(strip.text.y.right = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(color = "black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black")))

pdf("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/6_Umbilicaria_phaea_analyses/WC2_gradient_FET_sig_loci_x_BIO1_sync.pdf", height = 100, width = 4)
plot2
dev.off() 

pdf("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/6_Umbilicaria_phaea_analyses/FRHa_gradient_FET_sig_loci_x_BIO1_sync.pdf", height = 150, width = 4)
plot3
dev.off() 

pdf("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/6_Umbilicaria_phaea_analyses/FRHb_gradient_FET_sig_loci_x_BIO1_sync.pdf", height = 150, width = 4)
plot4
dev.off() 

pdf("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/6_Umbilicaria_phaea_analyses/CK1a_gradient_FET_sig_loci_x_BIO1_sync.pdf", height = 100, width = 4)
plot5
dev.off() 

pdf("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/6_Umbilicaria_phaea_analyses/CK1b_gradient_FET_sig_loci_x_BIO1_sync.pdf", height = 150, width = 4)
plot6
dev.off() 

pdf("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/6_Umbilicaria_phaea_analyses/CK1c_gradient_FET_sig_loci_x_BIO1_sync.pdf", height = 50, width = 4)
plot7
dev.off() 

pdf("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/6_Umbilicaria_phaea_analyses/CCG1_gradient_FET_sig_loci_x_BIO1_sync.pdf", height = 40, width = 4)
plot8
dev.off() 

pdf("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/6_Umbilicaria_phaea_analyses/CCG7_gradient_FET_sig_loci_x_BIO1_sync.pdf", height = 100, width = 4)
plot9
dev.off() 

pdf("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/6_Umbilicaria_phaea_analyses/CCG8_gradient_FET_sig_loci_x_BIO1_sync.pdf", height = 100, width = 4)
plot10
dev.off() 

pdf("~/Dropbox (Senckenberg)/Valim/04_Bioinformatic_work/05_Population_genetics/6_Umbilicaria_phaea_analyses/CCG9_gradient_FET_sig_loci_x_BIO1_sync.pdf", height = 150, width = 4)
plot11
dev.off() 

```


#### Trying a different approach: searching in qbGLM for most differentiated, and checking their GO terms




