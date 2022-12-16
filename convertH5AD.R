library(Seurat)
library(sceasy)
library(reticulate)
data <- readRDS("pbmc-1C-dm.rds")
obj <- CreateSeuratObject(counts = data[[1]])
obj$label <- data[[2]]
convertFormat(obj, from="seurat", to="anndata", outFile='pbmc-1C-dm.h5ad')