rm(list = ls())
library(WGCNA)
library(data.table)
library(RColorBrewer)
library(flashClust)
library(magrittr)
library(dplyr)

enableWGCNAThreads()

Tumor <- 
  fread("/Users/a703402454/Desktop/UCSC_Project_Testing/Test_Data_Set_70%/TPM_Analysis/Test_Set_TCGA_TARGET_TPM_Praveen_Mean_coding.csv",
        sep = ",",
        header = TRUE,
        data.table = FALSE
        )

nrow(Tumor)
ncol(Tumor)
ENSGids <-
  fread("/Users/a703402454/Desktop/UCSC_Project_Testing/Test_Data_Set_70%/ENSGids_Coding.csv",
        sep = ",",
        header = TRUE,
        data.table = FALSE
        )



Tumor.t <- t(Tumor)
rm(Tumor)

powers <- c(1:25)

sft <- 
  pickSoftThreshold(Tumor.t,
                    dataIsExpr = TRUE,
                    powerVector = powers,
                    corFnc = cor,
                    corOptions = list(method = "spearman"),
                    networkType = "signed",
                    RsquaredCut = 0.8,
                    verbose = 10
                    )


sizeGrWindow(9,5)
par(mfrow = c(1,2))
cex1 <- 0.9

setwd("/Users/a703402454/Desktop/UCSC_Project_Testing/Test_Data_Set_70%/TPM_Analysis/WGCNA_Plots_Praveen-Mean")
png("Soft_Power_Plot1.png",
    width = 800,
    height = 800
    )

plot(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit, signed R^2",
     type="n",
     main = paste("Scale independence")
     ) +
text(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels = powers,
     cex = cex1,
     col = "red"
     ) +
abline(h = 0.80,
       col = "red"
       )

dev.off()

png("Soft_Power_Plot2.png",
    width = 800,
    height = 800
    )

plot(sft$fitIndices[,1],
     sft$fitIndices[,5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
     ) +
text(sft$fitIndices[,1],
     sft$fitIndices[,5], labels = powers,
     cex = cex1,
     col = "red"
     )

dev.off()

power.df <- 
  data.frame(x = sft$fitIndices[,1],
             y = -sign(sft$fitIndices[,3])*sft$fitIndices[,2]
             )

fwrite(power.df,
       "sftPower_Praveen-Mean.csv",
       sep = ",",
       col.names = TRUE,
       row.names = FALSE
       )
