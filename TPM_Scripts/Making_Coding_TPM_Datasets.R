rm(list = ls())
library(data.table)
library(dplyr)
library(magrittr)
library(stringr)

setwd("/Users/a703402454/Desktop/UCSC_Project_Testing/Xena_Hub_Data_Files")

TTGtpmData <- 
  fread("TcgaTargetGtex_rsem_gene_tpm",
        sep = "\t",
        header = TRUE,
        showProgress = TRUE
        )
geneIDs <- TTGtpmData$sample
geneIDs %<>% str_remove( "\\..*")
TTGtpmData$sample <- geneIDs

setwd("/Users/a703402454/Desktop/UCSC_Project_Testing/Test_Data_Set_70%/")
ENSGidsCoding <- fread("ENSGids_Coding.csv")
ENSGidsCoding %<>% pull(ENSGids)
TTGtpmData %<>% filter(sample %in% ENSGidsCoding)

setwd("/Users/a703402454/Desktop/UCSC_Project_Testing/Xena_Hub_Data_Files/")
Meta <- fread("TcgaTargetGTEX_phenotype.txt")

TumorMeta <-  
  Meta %>%
  filter(grepl("^T", sample))
NormalMeta <-  
  Meta %>%
  filter(grepl("^G", sample))
  
TumorPrimaries <-
  TumorMeta %>%
  pull(`_primary_site`) %>%
  unique
NormalPrimaries <-
  NormalMeta %>%
  pull(`_primary_site`) %>%
  unique

CommonPrimaries <- intersect(NormalPrimaries, TumorPrimaries)
CommonPrimaries <- CommonPrimaries[CommonPrimaries != ""]

TumorMeta %<>% filter(`_primary_site` %in% CommonPrimaries)
NormalMeta %<>% filter(`_primary_site` %in% CommonPrimaries)

Tumor <- TTGtpmData[,TumorMeta$sample] 
Normal <-TTGtpmData[,NormalMeta$sample]
rm(TTGtpmData)

# Convert expression units to TPM + 1
Tumor %<>% apply(2, function(x) {(2^x) + 0.999})
Normal %<>% apply(2, function(x) {(2^x) + 0.999})

setwd("/Users/a703402454/Desktop/UCSC_Project_Testing/Test_Data_Set_70%/TPM_Analysis")
Tumor %>%
  as.data.frame %>%
  fwrite("Test_Set_TCGA_TARGET_TPM+1_coding.csv",
         sep = ",",
         col.names = TRUE,
         row.names = FALSE
         )
Normal %>%
  as.data.frame %>%
  fwrite("Test_Set_GTEX_TPM+1_coding.csv",
         sep = ",",
         col.names = TRUE,
         row.names = FALSE
         )
TumorMeta %>%
  as.data.frame %>%
  fwrite("Tumor_Meta_Common_Primaries.csv",
         sep = ",",
         col.names = TRUE,
         row.names = FALSE
         )
NormalMeta %>%
  as.data.frame %>%
  fwrite("Normal_Meta_Common_Primaries.csv",
         sep = ",",
         col.names = TRUE,
         row.names = FALSE
         )



