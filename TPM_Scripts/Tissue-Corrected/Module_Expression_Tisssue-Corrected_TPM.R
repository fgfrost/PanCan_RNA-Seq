rm(list = ls())
library(magrittr)
library(dplyr)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(car)
library(DescTools)
library(stringr)
library(wordspace)
library(dendextend)
library(flashClust)
library(pheatmap)

setwd("/Users/a703402454/Desktop/UCSC_Project_Testing/Test_Data_Set_70%/TPM_Analysis/WGCNA_Plots_Tissue-Corrected/Module_Genes/")
ModuleColors <- list.files(pattern = ".txt")
ModuleColors <- sub("module_",
                    "",
                    ModuleColors
                    )
ModuleColors <- sub("\\.txt",
                    "",
                    ModuleColors
                    )

ModuleGenes.lst <- list()

for (color in ModuleColors){
  tmp <-
    paste("/Users/a703402454/Desktop/UCSC_Project_Testing/Test_Data_Set_70%/TPM_Analysis/WGCNA_Plots_Tissue-Corrected/Module_Genes/module_",
          color,
          ".txt",
          sep = ""
    ) %>%
    fread(header = FALSE,
          data.table = FALSE
          )
  ModuleGenes.lst[[color]] <- tmp[[1]]
}

EnrichmentDF <-
  fread("/Users/a703402454/Desktop/UCSC_Project_Testing/Test_Data_Set_70%/TPM_Analysis/WGCNA_Plots_Tissue-Corrected/Module_Enrichment_Tissue-Corrected.csv",
        sep = ",",
        header = TRUE,
        data.table = FALSE
        )
EnrichmentDF <- EnrichmentDF[!is.na(EnrichmentDF$Module_File) &  !is.na(EnrichmentDF$Geneset),]

EnrichedColors <-
  unique(EnrichmentDF$Module_File) %>%
  str_remove_all("module_") %>%
  str_remove_all("\\.txt") %>%
  str_remove_all("/")

ModuleGenes.lst <- ModuleGenes.lst[EnrichedColors]

Tumor <-
  fread("/Users/a703402454/Desktop/UCSC_Project_Testing/Test_Data_Set_70%/TPM_Analysis/Test_Set_TCGA_TARGET_TPM_Tissue-Corrected_coding.csv",
        sep = ",",
        header = TRUE,
        data.table = FALSE
        )

ENSGids <- 
  fread("/Users/a703402454/Desktop/UCSC_Project_Testing/Test_Data_Set_70%/ENSGids_Coding.csv",
        sep = ",",
        header = TRUE,
        data.table = FALSE
        )
rownames(Tumor) <- ENSGids <- ENSGids[,"ENSGids"]

SamplePrimaries <-
  fread("/Users/a703402454/Desktop/UCSC_Project_Testing/Test_Data_Set_70%/TPM_Analysis/Tumor_Meta_Common_Primaries.csv",
        sep = ",",
        header = TRUE,
        data.table = FALSE
        )
SamplePrimaries <- SamplePrimaries[,c("sample", "_primary_site")]
colnames(SamplePrimaries) <- c("sample", "PrimarySite")
tmpDF <-
  fread("/Users/a703402454/Desktop/UCSC_Project_Testing/Test_Data_Set_70%/TPM_Analysis/WGCNA_Plots_Tissue-Corrected/HC_Group_Assignments_Tissue-Corrected.csv",
        sep = ",",
        header = TRUE,
        data.table = FALSE
        )
Clusters <- tmpDF$Group
names(Clusters) <- tmpDF$Sample
rm(tmpDF)
Clusters <- Clusters[colnames(Tumor)]

ClustNum <- 
  Clusters %>%
  max

rownames(SamplePrimaries) <- SamplePrimaries$sample
SamplePrimaries <- SamplePrimaries[,-1,drop = FALSE]
colnames(SamplePrimaries) <- "Tissue"
SamplePrimaries <- SamplePrimaries[colnames(Tumor),,drop = FALSE]
GroupingDF <- 
  data.frame("PrimarySite" = SamplePrimaries$Tissue,
             "Cluster" = Clusters,
             row.names = colnames(Tumor)
             )
GroupingDF <- GroupingDF[!is.na(GroupingDF$PrimarySite),]
Tumor <- Tumor[,rownames(GroupingDF)]
GroupingDF$Sample <- rownames(GroupingDF)
rownames(GroupingDF) <- NULL

ModuleExpression.lst <- list()
for (color in EnrichedColors) {
  ModuleExpression.lst[[color]] <- Tumor[ModuleGenes.lst[[color]],]
}
AllModuleGenes <- 
  ModuleGenes.lst %>%
  unlist
names(AllModuleGenes) <- NULL
NonModuleGenes <- ENSGids[!(ENSGids %in% AllModuleGenes)]

GroupingDF <- GroupingDF[!is.na(GroupingDF$PrimarySite),]

ModulePlotDFs.lst <- list()
for (color in EnrichedColors) {
  ModulePlotDFs.lst[[color]] <-
    data.frame("TissueGrp" = GroupingDF$PrimarySite,
               "MCgrp" = GroupingDF$Cluster %>%
                 paste("Cluster",.) %>%
                 factor(levels = paste("Cluster",c(1:ClustNum)
                                       )
                        ),
               "Module" = color
               )
}

ModulePlotDF <- do.call(rbind, ModulePlotDFs.lst)
rm(ModulePlotDFs.lst)


tmp <- 
  ModuleExpression.lst %>%
  lapply(nrow)
for (col in names(tmp)) {
  tmp[[col]] <- rep(col, tmp[[col]])
}

ModuleColorLabels <- unlist(tmp) 
names(ModuleColorLabels) <- NULL

ModuleExpression <- do.call(rbind, ModuleExpression.lst)
rm(ModuleExpression.lst)

palMaker <-
  brewer.pal(8,"Set2") %>%
  colorRampPalette
ColPal <-
  ModulePlotDF$TissueGrp %>%
  unique %>% 
  length %>%
  palMaker

setwd("/Users/a703402454/Desktop/UCSC_Project_Testing/Test_Data_Set_70%/TPM_Analysis/WGCNA_Plots_Tissue-Corrected/")

clusterTable <- table(ModulePlotDF$MCgrp,ModulePlotDF$TissueGrp)
clusterTable %<>%
  rowSums %>%
  divide_by(clusterTable,.) %>%
  multiply_by(100) %>%
  round(digits = 2) %>% 
  as.matrix %>%
  as.data.frame.matrix
rm(ModulePlotDF)

getwd() %>%
  paste0("/Cluster_Composition_Tissue-Corrected.csv") %>%
  fwrite(clusterTable,
         .,
         sep = ",",
         row.names = TRUE,
         col.names = TRUE
         )

sampleDist <- 
  ModuleExpression %>%
  as.matrix %>%
  dist.matrix(byrow = FALSE, method = "cosine") %>%
  as.dist

heatDendro <-
  sampleDist %>%
  flashClust(method = "ward") 
k <- 
  GroupingDF$Cluster %>%
  max
tmp <- 
  brewer.pal(8, "Dark2") %>%
  colorRampPalette
dendColPal <- 
  GroupingDF$Cluster %>%
  max %>%
  tmp
rm(tmp)

rownames(ModuleExpression) <- NULL

ColPalTissue <- ColPal
names(ColPalTissue) <- unique(GroupingDF$PrimarySite)

HeatPal <- brewer.pal(11,"RdBu")

ModuleExpression <- as.matrix(ModuleExpression)


primarySiteColsRaw <- c()
for (sample in GroupingDF$Sample) {
  primarySiteColsRaw[sample] <- as.character(GroupingDF$PrimarySite[GroupingDF$Sample == sample])
}
names(primarySiteColsRaw) <- NULL

primarySiteColsRaw <- 
  ColPalTissue[primarySiteColsRaw] %>%
  names
primarySiteColsDendro <- primarySiteColsRaw[order.dendrogram(as.dendrogram(heatDendro))] 

ModuleColorLabels %<>%
  data.frame("Module" = ., 
             stringsAsFactors = FALSE,
             check.names = FALSE
             )
primarySiteColsDendro <-
  data.frame("Cluster" = as.character(GroupingDF$Cluster[order.dendrogram(as.dendrogram(heatDendro))]),
             stringsAsFactors = FALSE,
             check.names = FALSE
             )
primarySiteColsRaw %<>%
  data.frame("Primary Site" = .,
             stringsAsFactors = FALSE,
             check.names = FALSE
             ) %>%
  arrange(`Primary Site`)

rownames(ModuleExpression) <- rownames(ModuleColorLabels) <-  unlist(ModuleGenes.lst)
rownames(primarySiteColsDendro) <- rownames(primarySiteColsRaw) <- colnames(ModuleExpression)

annoColors <- 
  list("Module" = unique(ModuleColorLabels$Module),
       "Primary Site" = ColPalTissue,
       "Cluster" = dendColPal
       )
names(annoColors[["Module"]]) <- annoColors[["Module"]]
names(annoColors[["Primary Site"]])  <-
  GroupingDF$PrimarySite %>%
  as.character %>%
  unique
names(annoColors[["Cluster"]])  <-
  GroupingDF$Cluster %>%
  max %>%
  seq(1, ., 1) %>%
  as.character


primarySiteColsDendro$tmp <- rownames(primarySiteColsDendro)
primarySiteColsDendro %<>% arrange(Cluster)
rownames(primarySiteColsDendro) <- primarySiteColsDendro$tmp
primarySiteColsDendro %<>% select(-tmp)


breaksClust <- 
  primarySiteColsDendro$Cluster %>%
  table %>%
  as.vector
breaksClust <- c(0, breaksClust)
for (i in 2:(length(breaksClust))) {
  breaksClust[i] <- breaksClust[i] + breaksClust[i-1]
}
breaksClust <- breaksClust[-length(breaksClust)]

breaksTissue <- 
  primarySiteColsRaw$`Primary Site` %>%
  table
breaksTissue <- breaksTissue[unique(primarySiteColsRaw$`Primary Site`)]
breaksTissue %<>% as.vector
breaksTissue <- c(0, breaksTissue)

for (i in 2:(length(breaksTissue))) {
  breaksTissue[i] <- breaksTissue[i] + breaksTissue[i-1]
}
breaksTissue <- breaksTissue[-length(breaksTissue)]


# Plotting Cluster-Wise Expression
ModuleExpression[,rownames(primarySiteColsDendro)] %>%
  pheatmap(col = HeatPal,
           kmeans_k = NA,
           border_color = NA,
           scale = "none",
           
           annotation_row = ModuleColorLabels,
           annotation_names_row = TRUE,
           annotation_col = primarySiteColsDendro,
           annotation_names_col = TRUE,
           annotation_legend = TRUE,
           annotation_colors = annoColors,
           
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           legend = TRUE,
           show_rownames = FALSE,
           show_colnames = FALSE,
           gaps_col = breaksClust,
           
           filename = "Cluster-Wise_Expression_of_Modules_Tissue-Corrected.png",
           height = 11,
           width = 12
           )



# Plotting Tissue-Wise Expression
ModuleExpression[,rownames(primarySiteColsRaw)] %>%
  pheatmap(col = HeatPal,
           kmeans_k = NA,
           border_color = NA,
           scale = "none",
           
           annotation_row = ModuleColorLabels,
           annotation_names_row = TRUE,
           
           annotation_col = primarySiteColsRaw,
           annotation_names_col = TRUE,
           
           annotation_legend = TRUE,
           annotation_colors = annoColors,
           
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           legend = TRUE,
           show_rownames = FALSE,
           show_colnames = FALSE,
           gaps_col = breaksTissue,
           
           filename = "Tissue-Wise_Expression_of_Modules_Tissue-Corrected.png",
           height = 11,
           width = 12
           )
