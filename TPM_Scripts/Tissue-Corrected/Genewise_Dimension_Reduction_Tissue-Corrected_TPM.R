rm(list = ls())
library(WGCNA)
library(data.table)
library(RColorBrewer)
library(flashClust)
library(magrittr)
library(wordspace)
library(dendextend)
library(dplyr)

enableWGCNAThreads()

Tumor <- 
  fread("/Users/a703402454/Desktop/UCSC_Project_Testing/Test_Data_Set_70%/TPM_Analysis/Test_Set_TCGA_TARGET_TPM_Tissue-Corrected_coding.csv",
        sep = ",",
        header = TRUE,
        data.table = FALSE
        )
TumorNames <- colnames(Tumor)


nrow(Tumor)
ncol(Tumor)
ENSGids <- 
  fread("/Users/a703402454/Desktop/UCSC_Project_Testing/Test_Data_Set_70%/ENSGids_Coding.csv",
        sep = ",",
        header = TRUE,
        data.table = FALSE
        )
gene.names <-ENSGids$ENSGids

Tumor.t <- t(Tumor)
rm(Tumor)

sftDF <- 
  fread("/Users/a703402454/Desktop/UCSC_Project_Testing/Test_Data_Set_70%/TPM_Analysis/WGCNA_Plots_Tissue-Corrected/sftPower_Tissue-Corrected.csv",
        sep = ",",
        header = TRUE,
        data.table = FALSE
        )

softPower <- 
  sftDF$x[sftDF$y >= 0.8] %>%
  min

softPower %>%
  paste("The Soft Power is",.) %>%
  print

adj <- 
  adjacency(Tumor.t,
            type = "unsigned",
            power = softPower
            )
TOM <- 
  TOMsimilarity(adj,
                TOMType = "signed"
                )

colnames(TOM) <- rownames(TOM) <- gene.names
dissTOM <- 1-TOM
geneTree <- 
  dissTOM %>%
  as.dist %>%
  flashClust(method = "average")

setwd("/Users/a703402454/Desktop/UCSC_Project_Testing/Test_Data_Set_70%/TPM_Analysis/WGCNA_Plots_Tissue-Corrected")

minModuleSize <- 20

dynamicMods <- 
  cutreeDynamic(dendro = geneTree,
                method = "tree",
                minClusterSize = minModuleSize
                )
table(dynamicMods)
dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)

par(mfrow = c(1,1))
png("All_Gene_Dendrogram_With_Colors_Test_Set_Tissue-Corrected.png",
    height = 800,
    width = 800
    )

plotDendroAndColors(geneTree,
                    dynamicColors,
                    "Dynamic Tree Cut",
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05,
                    main = "Gene dendrogram and module colors"
                    )
dev.off()

restGenes <-  (dynamicColors != "grey")
diss1 <- 1-TOMsimilarityFromExpr(Tumor.t[,restGenes],
                                 power = softPower
                                 )
colnames(diss1) <- rownames(diss1) <- gene.names[restGenes]
hier1 <- flashClust(as.dist(diss1),
                    method="average"
                    )

png("Module_Gene_Dendrogram_With_Colors_Test_Set_Tissue-Corrected.png",
    height = 800,
    width = 800
    )

plotDendroAndColors(hier1,
                    dynamicColors[restGenes],
                    "Dynamic Tree Cut",
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05,
                    main = "Gene dendrogram and module colors"
                    )
dev.off()

diag(diss1) <- NA

png("Test_Set_TOMplot_Tissue-Corrected.png",
    height = 800,
    width = 800
    )

TOMplot(diss1,
        hier1,
        as.character(dynamicColors[restGenes]
                     )
        )

dev.off()

module_colors <- setdiff(unique(dynamicColors),
                         "grey"
                         )

for (color in module_colors){
  module <- gene.names[which(dynamicColors == color)]
  paste("/Users/a703402454/Desktop/UCSC_Project_Testing/Test_Data_Set_70%/TPM_Analysis/WGCNA_Plots_Tissue-Corrected/Module_Genes/module_",
        color,
        ".txt",
        sep = ""
        ) %>%
    write.table(module,
                .,
                sep = "\t",
                row.names = FALSE,
                col.names = FALSE,
                quote = FALSE
                )
  }
moduleGenes.lst <- list()
for (color in module_colors){
  moduleGenes.lst[[color]]  <- gene.names[which(dynamicColors == color)]
}

module.order <-
  dynamicColors %>%
  as.factor() %>%
  tapply(1:ncol(Tumor.t),
         .,
         I
         ) %>%
  unlist


m <- t(t(Tumor.t[,module.order])/apply(Tumor.t[,module.order],
                                       2,
                                       max
                                       )
       )
SamplePrimaries <-
  fread("/Users/a703402454/Desktop/UCSC_Project_Testing/Test_Data_Set_70%/TPM_Analysis/Tumor_Meta_Common_Primaries.csv",
        sep = ",",
        header = TRUE,
        data.table = FALSE
        )
SamplePrimaries <- SamplePrimaries[,c("sample", "_primary_site")]
colnames(SamplePrimaries) <- c("sample", "PrimarySite")
rownames(SamplePrimaries) <- SamplePrimaries$sample
SamplePrimaries <- SamplePrimaries["PrimarySite"]
SamplePrimaries <- SamplePrimaries[rownames(Tumor.t),]
TissueVec <- unique(SamplePrimaries)
palMaker <-
  brewer.pal(8,"Set2") %>%
  colorRampPalette()
cvec <- 
  TissueVec %>%
  length %>%
  palMaker
rm(palMaker)
names(cvec) <- TissueVec
heatmapCols <- cvec[SamplePrimaries]
names(heatmapCols) <- NULL


sampleDist <- 
  Tumor.t[,!grepl("grey",names(module.order))] %>%
  dist.matrix(method = "cosine") %>%
  as.dist


heatDendro <-
  sampleDist %>%
  flashClust(method = "ward") %>%
  as.dendrogram
SampleOrder <- order.dendrogram(heatDendro)
k <- 
  heatDendro %>%
  find_k
colorpalDend <-
  brewer.pal(k$k,"Set2")
heatDendro %<>% 
  color_branches(k = k$k,
                 col = colorpalDend,
                 groupLabels = TRUE
                 )
groups <-
  heatDendro %>%
  cutree(k = k$k)
groups %>%
  data.frame("Sample" = names(.),
             "Group" = .
             ) %>%
  fwrite("/Users/a703402454/Desktop/UCSC_Project_Testing/Test_Data_Set_70%/TPM_Analysis/WGCNA_Plots_Tissue-Corrected/HC_Group_Assignments_Tissue-Corrected.csv",
         sep = ",",
         row.names = FALSE,
         col.names = TRUE
         )


png("All_Genes_Heatmap_Raw.png",
    height = 800,
    width = 800
    )

m %>%
  t %>%
  heatmap(zlim = c(0,1),
          col = gray.colors(100),
          Rowv = NA,
          Colv = heatDendro,
          labRow = NA,
          labCol = FALSE,
          scale = "none",
          RowSideColors = dynamicColors[module.order],
          ColSideColors = heatmapCols[SampleOrder]
          )
legend("right", 
       legend = names(cvec), 
       fill = cvec,
       bty = "n",
       xpd = TRUE
       )

dev.off()

png("Module_Genes_Heatmap_Raw.png",
    height = 800,
    width = 800
    )

heatmap(t(m[,!grepl("grey",names(module.order))]),
        zlim = c(0,1),
        col = gray.colors(100),
        Rowv = NA,
        Colv = heatDendro,
        labRow = NA,
        labCol = FALSE,
        scale = "none",
        RowSideColors = dynamicColors[module.order[!grepl("grey",names(module.order))]],
        ColSideColors = heatmapCols[SampleOrder]
        )
legend("right", 
       legend = names(cvec), 
       fill = cvec,
       bty = "n",
       xpd = TRUE
       )

dev.off()

MEList <- moduleEigengenes(Tumor.t,
                           colors = dynamicColors
                           )
MEs <- MEList$eigengenes

png("Modules_Heatmap_Raw.png",
    height = 800,
    width = 800
    )

plotEigengeneNetworks(MEs,
                      "",
                      marDendro = c(0,4,1,2),
                      marHeatmap = c(3,4,1,2)
                      )
dev.off()
