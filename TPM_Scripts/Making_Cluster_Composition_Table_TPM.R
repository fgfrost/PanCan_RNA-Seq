rm(list = ls())
library(data.table)
library(magrittr)
library(dplyr)
library(tidyr)

setwd("/Users/a703402454/Desktop/UCSC_Project_Testing/Test_Data_Set_70%/TPM_Analysis/")

SamplePrimaries <-
  fread("Tumor_Meta_Common_Primaries.csv",
        sep = ",",
        header = TRUE,
        data.table = FALSE
        )
SamplePrimaries <- SamplePrimaries[,c("sample", "_primary_site")]
colnames(SamplePrimaries) <- c("sample", "PrimarySite")
SamplePrimaries %<>% arrange(sample)

# Uncorrected data (aka raw)
rawClusters <-
  fread("WGCNA_Plots_Uncorrected/HC_Group_Assignments_Uncorrected.csv") %>%
  as_tibble
rawClusters %<>% arrange(Sample)
rawClusters %<>% rename(Cluster = Group)
rawClusters$PrimarySite <- SamplePrimaries$PrimarySite

# Tissue-corrected data (aka grandmean)
gmClusters <-
  fread("WGCNA_Plots_Tissue-Corrected/HC_Group_Assignments_Tissue-Corrected.csv") %>%
  as_tibble
gmClusters %<>% arrange(Sample)
gmClusters %<>% rename(Cluster = Group)
gmClusters$PrimarySite <- SamplePrimaries$PrimarySite

# Grand mean-corrected data (confusingly, aka praveen mean)
pmClusters <-
  fread("WGCNA_Plots_Praveen-Mean/HC_Group_Assignments_Uncorrected.csv") %>%
  as_tibble
pmClusters %<>% arrange(Sample)
pmClusters %<>% rename(Cluster = Group)
pmClusters$PrimarySite <- SamplePrimaries$PrimarySite

rm(SamplePrimaries)

# Make everything even

rawCounts <-
  rawClusters %>%
  group_by(PrimarySite, Cluster) %>%
  count %>%
  spread(PrimarySite, n) %>%
  as.data.frame
rm(rawClusters)
rawCounts %<>% data.frame("Dataset" = "Uncorrected", .)

gmCounts <-
  gmClusters %>%
  group_by(PrimarySite, Cluster) %>%
  count %>%
  spread(PrimarySite, n) %>%
  as.data.frame
rm(gmClusters)
gmCounts %<>% data.frame("Dataset" = "Tissue-corrected", .)


pmCounts <- 
  pmClusters %>%
  group_by(PrimarySite, Cluster) %>%
  count %>%
  spread(PrimarySite, n) %>%
  as.data.frame
rm(pmClusters)
pmCounts %<>% data.frame("Dataset" = "Grand mean-corrected", .)

totalCounts <- 
  rawCounts %>%
  list(gmCounts, pmCounts) %>%
  rbindlist
totalCounts[is.na(totalCounts)] <- 0

fwrite(totalCounts,
       "Clusterwise_Primary_Site_Counts_All_Datasets.csv",
       sep = ",",
       col.names = TRUE,
       row.names = FALSE
       )

totalCounts %>%
  filter(Dataset == "Uncorrected") %>%
  select(-Dataset, 
         -Cluster
         ) %>%
  as.matrix %>%
  colSums %>%
  matrix(nrow = 1, dimnames = list(NULL, names(.))) %>%
  as.data.frame %>%
  fwrite("Tissue_Counts_TCGA-TARGET.csv",
         sep = ",",
         col.names = TRUE,
         row.names = FALSE
         )


