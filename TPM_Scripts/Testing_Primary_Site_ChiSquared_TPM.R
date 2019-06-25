rm(list = ls())
library(data.table)
library(magrittr)
library(dplyr)

setwd("/Users/a703402454/Desktop/UCSC_Project_Testing/Test_Data_Set_70%/TPM_Analysis/")
clusterCounts <- 
  fread("Clusterwise_Primary_Site_Counts_All_Datasets.csv") %>%
  as_tibble


nullFreq <-
  clusterCounts %>%
  filter(Dataset == "Uncorrected") %>%
  select(-Dataset, -Cluster) %>%
  colSums %>%
  divide_by(sum(.))

tissueCounts <- list()
allData <-
  clusterCounts %>%
  pull(Dataset) %>%
  unique
  
for (data in allData) {
  tissueCounts[[data]] <- 
    clusterCounts %>%
    filter(Dataset == data) %>%
    select(-Dataset, -Cluster)
}

allSites <- colnames(tissueCounts[[1]])

pHyperResults <- tissueCounts
pHyperResults %<>%
  lapply(function(x){
    nc <- ncol(x)
    nr <- nrow(x)
    x <- matrix(0, nrow = nr, ncol = nc )
    x <- as.data.frame(x)
    colnames(x) <- allSites
    return(x)
    }
    )

for (data in allData) {
  j <- 
    pHyperResults[[data]] %>%
    nrow
  tmp <- 
    tissueCounts[[data]] %>%
    as.matrix
  for (site in allSites) {
    for (cluster in 1:j) {
      pHyperResults[[data]][cluster,site] <-
        phyper(q = tmp[cluster,site],
               m = sum(tmp[,site]),
               n = sum(tmp)-sum(tmp[,site]),
               k = sum(tmp[cluster,]),
               lower.tail = FALSE
               )
    }
  }
  pHyperResults[[data]] %<>% data.frame("Dataset" = data, .)
}

pHyperResults %<>% rbindlist

pHyperResults %>%
  fwrite("Hypergeometric_Test_By_Primary_Site_All_Datasets.csv",
         sep = ",",
         col.names = TRUE,
         row.names = FALSE
         )


