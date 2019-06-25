rm(list = ls())
library(mclust)
library(RColorBrewer)
library(data.table)
library(dplyr)
library(magrittr)
library(stringr)
library(umap)
library(RColorBrewer)
library(ggplot2)
library(gridExtra)




setwd("/Users/a703402454/Desktop/UCSC_Project_Testing/Test_Data_Set_70%/TPM_Analysis/WGCNA_Plots_Uncorrected/Module_Genes/")
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
  ModuleGenes.lst[[color]] <-
    color %>%
    paste("/Users/a703402454/Desktop/UCSC_Project_Testing/Test_Data_Set_70%/TPM_Analysis/WGCNA_Plots_Uncorrected/Module_Genes/module_",
          .,
          ".txt",
          sep = ""
          ) %>%
    fread(
      header = FALSE,
      data.table = FALSE
      )
}

AllModuleGenes <- unlist(ModuleGenes.lst)
names(AllModuleGenes) <- NULL

TcgaTarget <- 
  fread("/Users/a703402454/Desktop/UCSC_Project_Testing/Test_Data_Set_70%/TPM_Analysis/Test_Set_TCGA_TARGET_TPM+1_coding.csv",
        sep = ",",
        header = TRUE,
        data.table = FALSE
        )

ENSGids <- 
  fread("/Users/a703402454/Desktop/UCSC_Project_Testing/Test_Data_Set_70%/ENSGids_Coding.csv",
        header = TRUE,
        data.table = FALSE
        )
ENSGids$index <- c(1:nrow(ENSGids))
ENSGidsInModules <- ENSGids[ENSGids$ENSGids %in% AllModuleGenes,]

TT_Modules <- TcgaTarget[ENSGidsInModules$index,]
rm(TcgaTarget)

SampleNames <- colnames(TT_Modules)
TT_Modules <- as.matrix(TT_Modules)
TT_Modules.t <- t(TT_Modules)
gene.names <- ENSGidsInModules$ENSGids

# BIC <- mclustBIC(TT_Modules.t,
#                  G = c(1:30),
#                  verbose = TRUE
#                  )
# setwd("/Users/a703402454/Desktop/UCSC_Project_Testing/Test_Data_Set_70%/Praveen_Mean_Data/Mclust_Plots/")
# png("BIC_Model_Fitting_Praveen-Mean.png",
#     width = 800,
#     height = 800
#     )

# plot.mclustBIC(BIC,
#                G = c(1:30),
#                legendArgs = list(x = "bottomright")
#                )
#   title(main = "Bayesian Information Criterion (BIC) Model Fitting")
# dev.off()

# summary(BIC)
# mod1 <- Mclust(TT_Modules.t,
#                x = BIC,
#                verbose = TRUE
#                )
# summary(mod1,
#         parameters = TRUE
#         )

# clustersMod1 <- mod1$classification

tmpDF <-
  fread("/Users/a703402454/Desktop/UCSC_Project_Testing/Test_Data_Set_70%/TPM_Analysis/WGCNA_Plots_Uncorrected/HC_Group_Assignments_Uncorrected.csv",
        sep = ",",
        header = TRUE,
        data.table = FALSE
        )
clustersMod1 <- tmpDF$Group
names(clustersMod1) <- tmpDF$Sample
rm(tmpDF)
clustersMod1 <- clustersMod1[SampleNames]


ClustNumMod1 <- max(clustersMod1)
MclustClustersMod1 <- paste("Cluster",
                            c(1:ClustNumMod1)
                            )

UMAPtest <- 
  umap(TT_Modules.t,
       verbose = TRUE
       )

colnames(UMAPtest$layout) <- c("UMAP 1","UMAP 2")
ClustNum <- max(clustersMod1)
if(ClustNum > 8){
  Set2 <- brewer.pal(8,
                     "Set2"
                     )
  ColorFun <- colorRampPalette(Set2)
  ColorPal <- ColorFun(ClustNum)
} else {
  ColorPal <- brewer.pal(ClustNum,
                         "Set2"
                         )
  ColorFun <- 
    brewer.pal(8,
               "Set2"
               ) %>%
    colorRampPalette

}

ClusterColorsMod1 <- ColorPal[clustersMod1]

######### Tissues Assignment Nonsense #########

SamplePrimaries <-
  fread("/Users/a703402454/Desktop/UCSC_Project_Testing/Test_Data_Set_70%/TPM_Analysis/Tumor_Meta_Common_Primaries.csv",
        sep = ",",
        header = TRUE,
        data.table = FALSE
        )
SamplePrimaries <- SamplePrimaries[,c("sample", "_primary_site")]
colnames(SamplePrimaries) <- c("sample", "PrimarySite")
PrimarySites <- SamplePrimaries
SampleNames %<>% str_remove("\\..")
rownames(PrimarySites) <- SampleNames <- PrimarySites$sample
PrimarySites <- PrimarySites[SampleNames,]
PrimarySites <- PrimarySites[complete.cases(PrimarySites),]
MetaMasterPrimaries <- PrimarySites$PrimarySite
rm(PrimarySites)
allSites <- unique(MetaMasterPrimaries)
allSites <- allSites[order(allSites)]

palMaker <-
  brewer.pal(8,"Set2") %>%
  colorRampPalette
TissuePal <- 
  MetaMasterPrimaries %>%
  unique %>%
  length %>%
  palMaker
names(TissuePal) <- allSites

ShapePal <-
  allSites %>%
  length %>%
  divide_by(3) %>%
  rep(c(15,16,17,18),.)
ShapePal <- ShapePal[1:length(allSites)]
names(ShapePal) <- allSites


UMAPplot <- as.data.frame(UMAPtest$layout)
UMAPplot$grpMC <- 
  paste("Cluster",
        clustersMod1
        ) %>%
  factor(levels = paste("Cluster",c(1:ClustNum)))
UMAPplot <- UMAPplot[SampleNames,]
UMAPplot$grpTissues <- SamplePrimaries$PrimarySite



setwd("/Users/a703402454/Desktop/UCSC_Project_Testing/Test_Data_Set_70%/TPM_Analysis/Mclust_Plots_Uncorrected/")
png("mClust_VS_Tissue_UMAP.png",
    width = 1600,
    height = 800
    )


mcPlot <-
  ggplot(UMAPplot,
         aes(x = `UMAP 1`,
             y = `UMAP 2`,
             color = grpMC,
             shape = grpMC
             )
         ) +
  geom_point(size = 2.2) +
  guides(color = guide_legend(ncol = 1),
         shape = guide_legend(override.aes = list(size = 5))
         ) +
  scale_color_manual(values = ColorPal) +
  scale_shape_manual(values = rep(16, ClustNum) ) +
  theme_bw() +
  theme(plot.title = element_text(size = 24,
                                  face = "bold"
                                   ),
        axis.title = element_text(size = 20,
                                  face = "bold"
                                  ),
        legend.title = element_text(size = 18,
                                    face = "bold"
                                    ),
        axis.text = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.key.height = unit(1, "cm"), 
        plot.margin = margin(t = 1,
                             r = 1,
                             b = 1,
                             l = 1,    
                             unit = "cm"
                             )
        ) +
  labs(x = "UMAP 1",
       y = "UMAP 2",
       title = "Expression Profile-Grouped Tumor Transcriptomes",
       color = "Cluster",
       shape = "Cluster"
       )

tissuePlot <- 
  ggplot(UMAPplot,
         aes(x = `UMAP 1`,
             y = `UMAP 2`,
             color = grpTissues,
             shape = grpTissues
             )
         ) +
  geom_point(size = 2.2) +
  guides(color = guide_legend(ncol = 1),
         shape = guide_legend(override.aes = list(size = 5))
         ) +
  scale_color_manual(values = TissuePal) +
  scale_shape_manual(values = ShapePal) +
  theme_bw() +
  theme(plot.title = element_text(size = 24,
                                  face = "bold"
                                  ),
        axis.title = element_text(size = 20,
                                  face = "bold"
                                  ),
        legend.title = element_text(size = 18,
                                    face = "bold"
                                    ),
        axis.text = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.key.height = unit(1, "cm"),
        plot.margin = margin(t = 1,
                             r = 1,
                             b = 1,
                             l = 1,
                             unit = "cm"
                             )
        ) +
  labs(x = "UMAP 1",
       y = "UMAP 2",
       title = "Tissue-Grouped Tumor Transcriptomes",
       color = "Primary Site",
       shape = "Primary Site"
       )
grid.arrange(mcPlot,
             tissuePlot,
             ncol = 2
             )

dev.off()

# ClustOut <-
#   clustersMod1 %>%
#   as.data.frame
# colnames(ClustOut) <- "Cluster"
# ClustOut$Sample <-
#   rownames(ClustOut) %>%
#   gsub("\\.","-",.)
# rownames(ClustOut) <- NULL

# fwrite(ClustOut,
#        "/Users/a703402454/Desktop/UCSC_Project_Testing/Test_Data_Set_70%/Praveen_Mean_Data/Mclust_Plots/Cluster_Assignment_Praveen-Mean.csv",
#        sep = ",",
#        col.names = TRUE,
#        row.names = FALSE
#        )
