rm(list = ls())
library(data.table)
library(magrittr)
library(dplyr)
library(WebGestaltR)

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


ModuleFolder <- "/Users/a703402454/Desktop/UCSC_Project_Testing/Test_Data_Set_70%/TPM_Analysis/WGCNA_Plots_Praveen-Mean/Module_Genes/"
Module_Enrichment <-
  WebGestaltRBatch(enrichMethod = "ORA",
                   organism = "hsapiens",
                   enrichDatabase = "pathway_Wikipathway",
                   enrichDatabaseType = "ensembl_gen",
                   interestGeneFolder = ModuleFolder,
                   interestGeneType = "ensembl_gene_id",
                   collapseMethod = "mean",
                   referenceSet = "genome_protein-coding",
                   fdrMethod = "BH",
                   isOutput = FALSE,
                   nThreads = 7
                   )

#listGeneSet()

names(Module_Enrichment) <- ModuleColors[order(ModuleColors)]
Lengths <- c()
for(i in 1:length(ModuleColors)){
  Lengths[i] <- length(unlist(Module_Enrichment[[i]]))
  }

NamesCol <- c("Geneset","Description","Link","C","O","E","R","P_Value",
              "FDR","Overlap_Genes_Entrez", "Overlap_Genes_ENSG")
namesLen <- length(NamesCol) 
Ncol <- length(NamesCol)
for (i in 1:length(Module_Enrichment)) {
  Temp <- unlist(Module_Enrichment[[i]])
  TempFile <- Temp[1]
  TempFile <- 
    gsub("/Users/a703402454/Desktop/UCSC_Project_Testing/Test_Data_Set_70%/TPM_Analysis/WGCNA_Plots_Praveen-Mean/Module_Genes/",
         "",
         TempFile
         )
  Temp <- Temp[-1]
  if(length(Temp) != 0) {
    Nrow <- 
      Temp %>%
      length %>%
      divide_by(namesLen)
    names(Temp) <- NULL
    TempMat <- matrix(Temp,
                      nrow =  Nrow
                      )
    Module_Enrichment[[i]] <- 
      data.frame("Module_File" = rep(TempFile,
                                     Nrow
                                     ),
                 TempMat
                 )
    } else {
      TempMat <- matrix(c(TempFile,rep("NA",
                                       Ncol
                                       )
                          ),
                        ncol = Ncol + 1

                        )

      Module_Enrichment[[i]] <- data.frame(TempMat)
  }
}

for (i in 1:length(Module_Enrichment)) {
  if(nrow(Module_Enrichment[[i]]) == 0) {
    Module_Enrichment[[i]] <- 
      names(Module_Enrichment)[i] %>%
      paste("module_",.,".txt",sep = "") %>%
      c(.,rep("NA",(Ncol-1))) %>%
      matrix %>%
      t
  } else {}
}

for (i in 1:length(Module_Enrichment)) {
  colnames(Module_Enrichment[[i]]) <- c("Module_File",NamesCol)
  rownames(Module_Enrichment[[i]]) <- NULL
}
Module_Enrichment.df <- do.call(rbind,Module_Enrichment)
enrichLengths <- c()
for (i in 1:nrow(Module_Enrichment.df)) {
  enrichLengths[i] <-
    Module_Enrichment.df$Overlap_Genes_Entrez[i] %>%
    as.character %>%
    strsplit(";") %>%
    unlist %>%
    length
  
}
Module_Enrichment.df$Number_of_Enriched_Genes <- enrichLengths

fwrite(Module_Enrichment.df,
       "/Users/a703402454/Desktop/UCSC_Project_Testing/Test_Data_Set_70%/TPM_Analysis/WGCNA_Plots_Praveen-Mean/Module_Enrichment_Praveen-Mean.csv",
       sep = ",",
       col.names = TRUE,
       row.names = FALSE
       )
