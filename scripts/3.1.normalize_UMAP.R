
# Libraries ---------------------------------------------------------------

library(Seurat)
library(tidyverse)

object = readRDS("data/HFC_fresh.rds")

#lets check falsecodes and negprb in RNA Assay
counts = GetAssayData(object, assay = "RNA")

#check negprb first
dim(counts)
head(row.names(counts))
row.names(counts) %>% grep("NegPrb", ., value = TRUE)

#FalseCodes
row.names(counts) %>% grep("FalseCode", ., value = TRUE)


#transform with SCTransform
object = SCTransform(object, assay = "RNA", new.assay.name = "SCT")

#PCA
object = RunPCA(object, assay = "SCT", reduction.name = "PCA", npcs = 50)
#view PCA
DimPlot(object, reduction = "PCA", raster = FALSE)


#run UMAP
object = RunUMAP(object, reduction = "PCA", reduction.name = "UMAP", dims = 1:30, repulsion.strength = 10)
#view UMAP
DimPlot(object, reduction = "UMAP", raster = FALSE)


#finding clusters
object = FindNeighbors(object, reduction = "PCA", dims = 1:30)
object = FindClusters(object, algorithm = 1, resolution = 0.3)

#look at clusters on UMAP
DimPlot(object, reduction = "UMAP", group.by = "SCT_snn_res.0.3", raster = FALSE)
-