rm(list = ls())
library(data.table)
library(dplyr)
library(magrittr)
library(stringr)
library(doParallel)

nCores <- 
  detectCores() %>% 
  subtract(2)
cl <- makeCluster(nCores)
registerDoParallel(cl)

setwd("/Users/a703402454/Desktop/UCSC_Project_Testing/Test_Data_Set_70%/TPM_Analysis")
Tumor <- 
  fread("Test_Set_TCGA_TARGET_TPM+1_coding.csv",
        sep = ",",
        header = TRUE,
        data.table = FALSE
        )

Normal <- 
  fread("Test_Set_GTEX_TPM+1_coding.csv",
        sep = ",",
        header = TRUE
        )
TumorMeta <- 
  fread("Tumor_Meta_Common_Primaries.csv",
        sep = ",",
        header = TRUE
        )
NormalMeta <- 
  fread("Normal_Meta_Common_Primaries.csv",
        sep = ",",
        header = TRUE
        )

allSites <- unique(NormalMeta$`_primary_site`)

normalMeans <- 
  foreach(site = allSites) %dopar% {
    tmpNames <- NormalMeta$sample[NormalMeta$`_primary_site` == site]
    means <- rowMeans(Normal[,colnames(Normal) %in% tmpNames])
    return(means)
  }
rm(Normal)
names(normalMeans) <- allSites

for (site in allSites) {
  tmpNames <- TumorMeta$sample[TumorMeta$`_primary_site` == site]
  Tumor[,colnames(Tumor) %in% tmpNames] <- log(Tumor[,colnames(Tumor) %in% tmpNames]/normalMeans[[site]])
}

Tumor %>%
  fwrite("Test_Set_TCGA_TARGET_TPM_Tissue-Corrected_coding.csv",
         sep = ",",
         col.names = TRUE,
         row.names = FALSE
         )
stopCluster(cl)



