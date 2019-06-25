rm(list = ls())
library(magrittr)
library(dplyr)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(car)
library(DescTools)
library(stringr)

setwd("/Users/a703402454/Desktop/UCSC_Project_Testing/Test_Data_Set_70%/TPM_Analysis/WGCNA_Plots_Praveen-Mean/Module_Genes/")
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
    paste("/Users/a703402454/Desktop/UCSC_Project_Testing/Test_Data_Set_70%/TPM_Analysis/WGCNA_Plots_Praveen-Mean/Module_Genes/module_",
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
  fread("/Users/a703402454/Desktop/UCSC_Project_Testing/Test_Data_Set_70%/TPM_Analysis/WGCNA_Plots_Praveen-Mean/Module_Enrichment_Praveen-Mean.csv",
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
  fread("/Users/a703402454/Desktop/UCSC_Project_Testing/Test_Data_Set_70%/TPM_Analysis/Test_Set_TCGA_TARGET_TPM_Praveen_Mean_coding.csv",
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
  fread("/Users/a703402454/Desktop/UCSC_Project_Testing/Test_Data_Set_70%/TPM_Analysis/WGCNA_Plots_Praveen-Mean/HC_Group_Assignments_Uncorrected.csv",
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
rm(SamplePrimaries)
rm(Clusters)



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
               "avg" = colMeans(ModuleExpression.lst[[color]], na.rm = TRUE),
               "Module" = color
               )
}
ModulePlotDF <- do.call(rbind,ModulePlotDFs.lst)
palMaker <-
  brewer.pal(8,"Set2") %>%
  colorRampPalette
ColPal <-
  ModulePlotDF$TissueGrp %>%
  unique %>% 
  length %>%
  palMaker

setwd("/Users/a703402454/Desktop/UCSC_Project_Testing/Test_Data_Set_70%/TPM_Analysis/WGCNA_Plots_Praveen-Mean/")

clusterTable <- table(ModulePlotDF$MCgrp,ModulePlotDF$TissueGrp)
clusterTable %<>%
  rowSums %>%
  divide_by(clusterTable,.) %>%
  multiply_by(100) %>%
  round(digits = 2) %>% 
  as.matrix %>%
  as.data.frame.matrix


getwd() %>%
  paste0("/Cluster_Composition_Praveen-Mean.csv") %>%
  fwrite(clusterTable,
         .,
         sep = ",",
         row.names = TRUE,
         col.names = TRUE
         )

png("Tissue-Wise_Expression_of_Modules_Praveen-Mean.png",
    height = 27,
    width = 95,
    units = "cm",
    res = 400
    )

ggplot(ModulePlotDF,
       aes(y = avg,
           fill = TissueGrp
           )
       ) + 
  geom_boxplot(aes(color = TissueGrp)) +
  geom_boxplot(aes(fill = TissueGrp), 
               outlier.color = NA
               ) +
  geom_hline(yintercept = 0,
             lty = "dashed",
             color = "red"
             ) +
  scale_color_manual(values = ColPal) +
  scale_fill_manual(values = ColPal) +
  theme_bw() +
  theme(strip.text.x = element_text(size = 25,
                                    face = "bold",
                                    margin = margin(0.35, 0, 0.35, 0, unit = "cm")
                                    ),
        axis.title = element_text(size = 27, face = "bold"),
        legend.title = element_text(size = 27, face = "bold"),
        legend.title.align = 0.5,
        axis.text = element_text(size = 25),
        legend.text = element_text(size = 25),
        legend.spacing.x = unit(0.75, "cm"),
        legend.spacing.y = unit(0.75, "cm"),
        legend.key.height = unit(1, "cm"),
        plot.margin = margin(1, 1, 1, 1, unit = "cm"),
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank()
        ) +
  labs(y = "Tumor Expression (Ln(Tumor/Normal))",
       fill = "Tissue",
       color = "Tissue") +
  facet_wrap(~Module,
             nrow = 1,
             strip.position = "top"
             )  
dev.off()


png("Cluster-Wise_Expression_of_Modules_Praveen-Mean.png",
    height = 27,
    width = 95,
    units = "cm",
    res = 400
    )
ggplot(ModulePlotDF,
       aes(y = avg,
           fill = Module
           )
       ) + 
  geom_boxplot(aes(color = Module)) +
  geom_boxplot(aes(fill = Module),
               outlier.color = NA) +
  geom_hline(yintercept = 0,
             lty = "dashed",
             color = "red"
             ) +
  scale_fill_manual(values = EnrichedColors,
                    labels = EnrichedColors
                    ) +
  scale_color_manual(values = EnrichedColors,
                     labels = EnrichedColors
                     ) +
  theme_bw() +
  theme(strip.text.x = element_text(size = 25,
                                    face = "bold",
                                    margin = margin(0.35, 0, 0.35, 0, unit = "cm")
                                    ),
        axis.title = element_text(size = 27, face = "bold"),
        legend.title = element_text(size = 27, face = "bold"),
        legend.title.align = 0.5,
        axis.text = element_text(size = 25),
        legend.text = element_text(size = 25),
        legend.spacing.x = unit(0.75, "cm"),
        legend.spacing.y = unit(0.75, "cm"),
        legend.key.height = unit(1, "cm"),
        plot.margin = margin(1, 1, 1, 1, unit = "cm"),
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank()
        ) +
  facet_wrap(~MCgrp,
             ncol = ClustNum,
             strip.position = "top"
             ) +
  labs(y = "Tumor Expression (Ln(Tumor/Normal))")

dev.off()

# Stats Time!


clusterWiseModuleAOV.lst <- list()
for (color in EnrichedColors) {
  varTest <- 
    leveneTest(avg ~ MCgrp,
               data = ModulePlotDF[ModulePlotDF$Module == color,]
               )
  if (varTest$`Pr(>F)`[1] > 0.05) {
    clusterWiseModuleAOV.lst[[color]]$model <- 
      aov(avg ~ MCgrp,
          data = ModulePlotDF[ModulePlotDF$Module == color,]
      )
    clusterWiseModuleAOV.lst[[color]]$test <- "anova"
  } else {
    clusterWiseModuleAOV.lst[[color]]$model <-
      kruskal.test(avg ~ MCgrp,
                   data = ModulePlotDF[ModulePlotDF$Module == color,]
                   )
    clusterWiseModuleAOV.lst[[color]]$test <- "kruskal"
    }
  }

clusterWiseModulePostHoc.lst <- list()
for (color in EnrichedColors) {
  if (clusterWiseModuleAOV.lst[[color]]$test == "kruskal") {
    eval <- clusterWiseModuleAOV.lst[[color]]$model$p.value 
  } else {
    tmp <- 
      clusterWiseModuleAOV.lst[[color]]$model %>%
      summary %>%
      unlist
    eval <- tmp["Pr(>F)1"]
    rm(tmp)
  }
  eval %>%
    paste(color,.,sep="-") %>%
    print
  tmp <- 
    DunnTest(avg ~ MCgrp,
             data = ModulePlotDF[ModulePlotDF$Module == color,]
             )
  tmpDF <-
    tmp[[1]] %>%
    as.data.frame %>%
    select(.,pval) %>%
    data.frame("P-Value" = .)
  rm(tmp)
  colnames(tmpDF) %<>% paste0("-",color)
  clusterWiseModulePostHoc.lst[[color]] <- tmpDF
}

clusterWiseModulePostHoc <-
  clusterWiseModulePostHoc.lst %>%
  do.call(cbind,.)

getwd() %>%
  paste0("/Cluster-Cluster_Module_Expression_PostHoc_Uncorrected.csv") %>%
  fwrite(clusterWiseModulePostHoc,
         .,
         sep = ",",
         row.names = TRUE,
         col.names = TRUE
  )



aov.lst <- list()
for (color in EnrichedColors) {
  aov.lst[[color]] <-
    aov(avg ~ TissueGrp,
        data = ModulePlotDF[ModulePlotDF$Module == color,]
        )
  }
aovResults.lst <- lapply(aov.lst, summary)
aovResults.lst <- lapply(aovResults.lst, unlist)
aovResultsDF <- 
  aovResults.lst %>%
  do.call(rbind,.) %>%
  data.frame("Module" = rownames(.),.)
rm(aovResults.lst)
rownames(aovResultsDF) <- NULL

PostHoc.lst <- lapply(aov.lst, TukeyHSD)

for (color in EnrichedColors) {
  PostHoc.lst[[color]] <-
    PostHoc.lst[[color]]$TissueGrp %>%
    data.frame("Comparison" = rownames(.),
               "Module" = color,
               .)
}

PrimarySites <- 
  GroupingDF$PrimarySite %>%
  unique
comparisonsDF <-
  PostHoc.lst[[1]][1] %>%
  as.matrix %>%
  strsplit("-") %>%
  unlist %>%
  matrix(ncol = 2, byrow = TRUE) %>%
  as.data.frame
colnames(comparisonsDF) <- c("Tissue1","Tissue2")
comparisonsDF$Col1 <- match(comparisonsDF$Tissue1,PrimarySites)
comparisonsDF$Col2 <- match(comparisonsDF$Tissue2,PrimarySites)

pValueDFs.lst <- list()
for (color in EnrichedColors) {
  tempDF <-
    PrimarySites %>%
    length %>%
    matrix(1,
           ncol = .,
           nrow = .,
           dimnames = list(PrimarySites,PrimarySites)
    ) %>%
    as.data.frame
  for (i in 1:nrow(comparisonsDF)) {
    x <- comparisonsDF$Col1[i]
    y <- comparisonsDF$Col2[i]
    tempDF[x,y] <- PostHoc.lst[[color]]$p.adj[i]
    tempDF[y,x] <- PostHoc.lst[[color]]$p.adj[i]
    rm(x)
    rm(y)
  }
  pValueDFs.lst[[color]] <- tempDF
  rm(tempDF)
}

for (color in EnrichedColors) {
  pValueDFs.lst[[color]] %<>%
    data.frame("Module" = color,
               "P-Value_Matrix" = rownames(.),
               .)
  rownames(pValueDFs.lst[[color]]) <- NULL
}

pValueDFsAll <-
  pValueDFs.lst %>%
  do.call(rbind,.)
rownames(pValueDFsAll) <- NULL

getwd() %>%
  paste0("/Tissue-Wise_Module_Expression_P_Values_Uncorrected.csv") %>%
  fwrite(pValueDFsAll,
         .,
         sep = ",",
         row.names = FALSE,
         col.names = TRUE
  )
