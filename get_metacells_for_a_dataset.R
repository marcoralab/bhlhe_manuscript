library(dplyr)
library(data.table)
library(Seurat)
library(scater)
source('pisces_2019_functions/functions/process-utils.R')
source('pisces_2019_functions/functions/cluster-functions.R')
source('pisces_2019_functions/functions/viper-utils.R')
library(ggplot2)
library(ggpubr)
library(viper)
library(pheatmap)
library(RColorBrewer)
library(MUDAN)
library(umap)
library(data.table)
library(Matrix)

### Functions below were adapted from  the Single Cell Analysis Boot Camp at Columbia University (instructors Lukas Vlahos and Pasquale Laise) in 2019 and are also provided in the pisces_2019_functions folder. Updated version of the PISCES package can be found here: 
### https://github.com/califano-lab/PISCES 

QCTransform_count <- function(raw.mat, minCount = minCount, maxCount = 100000, minGeneReads = 1) {
  filt.mat <- raw.mat[, colSums(raw.mat) > minCount & colSums(raw.mat) < maxCount]
  filt.mat <- filt.mat[ rowSums(raw.mat) >= minGeneReads ,]
  rem.genes <- nrow(raw.mat) - nrow(filt.mat); rem.cells <- ncol(raw.mat) - ncol(filt.mat)
  print(paste('Removed ', rem.genes, ' genes and ', rem.cells, ' cells.', sep =''))
  return(filt.mat)
}

get_meta_cells_matrix <- function(data, ### matrix with rows as gene names and columns as cell IDs/column names of choice
                                  numNeighbors, ### number of neighbors for MetaCells
                                  subSize,     ### number of MetaCells
                                  file_name,   ### output file name
                                  gene="symbol", ### either "symbol" OR "ensemble" depending on the dataset
                                  species) ### either mouse or human
  {
  mt.table <- read.table("single-cell-pipeline-master/mt-geneList.upd.csv", sep = '\t', header = TRUE, stringsAsFactors = FALSE)
  if (gene=="ensemble" & species=="human") {
    print(paste("Gene setting is in Ensemble and species = human"))
    mt.ensg <- mt.table$Gene.stable.ID
  } else if (gene=="ensemble" & species=="mouse") {
    print(paste("ene setting is in Ensemble and species = mouse"))
    mt.ensg <- mt.table$Mouse.gene.stable.ID
  } else if (gene=="symbol" & species=="human") {
    print(paste("Gene setting is in Symbol and species = human"))
    mt.ensg <- mt.table$HGNC.symbol
  } else if (gene=="symbol" & species=="mouse") {
    print(paste("Gene setting is in Symbol and species = mouse"))
    mt.ensg <- mt.table$MGI.symbol
  }
  mt.perc <- MTPercent(data, mt.ensg)
  thresh.cells <- names(mt.perc)[which(mt.perc < 0.1)]
  print(paste(length(thresh.cells), 'out of',ncol(data),'cells survived after cleaning due to mitochodrial contamination', sep =' '))
  data <- data[, thresh.cells ]
  filt.mat <- QCTransform_count(data,200)
  print(paste('The final number of cells is ', ncol(filt.mat), sep=""))
  print(paste('CPM normalizing...'))
  cpm.mat <- CPMTransform(filt.mat)
  print(paste('Computing a distance matrix...'))
  pearson <- cor(cpm.mat, method = c("pearson"))
  dist.mat <- (1 - pearson)/2
  print(paste('Running MetaCells...'))
  meta_cells_out <- MetaCells(filt.mat, dist.mat, numNeighbors = numNeighbors, subSize=subSize)
  cpm.meta.mat <- CPMTransform(meta_cells_out)
  print(paste('Writing out a file for ARACNe...'))
  ARACNeTable(cpm.meta.mat,file_name)
}

library(Matrix)


set.seed(0)
raw.mat_init <- readRDS("Example_dataset.rds")
mat <- as.matrix(raw.mat_init)
get_meta_cells_matrix(mat,15,200,"Example_dataset.neighb_15_subsize_200_MetaCells.for_ARACNe", gene="symbol",species="mouse")
