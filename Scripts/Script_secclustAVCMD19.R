## D19 datasets AVCM, secondary clustering

library(Seurat)
library(URD)
library(Matrix)
library(fossil) 
library(dplyr)
library(plyr)

#import dataset with all 4 subtypes
setwd("/path/to/directory")
AVCM<-readRDS("00_AVCMsubsetofD19.rds")
AVCM
#An object of class Seurat 
#38611 features across 442 samples within 3 assays 
#Active assay: RNA (18444 features, 0 variable features)
#2 other assays present: SCT, integrated
#2 dimensional reductions calculated: pca, umap
##SCTransform + batch correction

data.list <- SplitObject(AVCM, split.by = "batch")

for (i in 1:length(data.list)) { data.list[[i]] <- SCTransform(data.list[[i]], verbose = TRUE)}

#select features for downstream integration, and run PrepSCTIntegration, which 
#ensures that all necessary Pearson residuals have been calculated

data.features <- SelectIntegrationFeatures(object.list = data.list, nfeatures = 3000)
data.list <- PrepSCTIntegration(object.list = data.list, anchor.features = data.features, 
                                verbose = TRUE)

#identify anchors and integrate the datasets
#make sure to set normalization.method = 'SCT'

data.anchors <- FindIntegrationAnchors(object.list = data.list, normalization.method = "SCT", 
                                       anchor.features = data.features, verbose = TRUE)
data.integrated <- IntegrateData(anchorset = data.anchors, normalization.method = "SCT",verbose = TRUE)

#proceed with downstream analysis (i.e. visualization, clustering) on the integrated dataset
#do not scale integrated data

data.integrated <- RunPCA(data.integrated, verbose = FALSE)
DimPlot(data.integrated, reduction = "pca", dims = c(1,2))
DimHeatmap(data.integrated, dims = 1, cells = 500, balanced = TRUE)

ElbowPlot(data.integrated, ndims = 60)
data.integrated <- RunUMAP(data.integrated, dims = 1:20)
DimPlot(data.integrated, reduction = "umap")
data.integrated <- FindNeighbors(data.integrated, dims = 1:20)
data.integrated <- FindClusters(object = data.integrated, reduction = "pca",dims = 1:20, resolution = 0.6,random.seed = 2020)
head(Idents(data.integrated), 5)
data.integrated <- RunUMAP(data.integrated, dims = 1:20)
DimPlot(data.integrated, reduction = "umap", label = TRUE, repel = TRUE)

DimPlot(data.integrated, group.by = "orig.ident", label = TRUE, repel = TRUE)

DefaultAssay(data.integrated) <- "RNA"

# Normalize RNA data for visualization purposes

data.integrated <- NormalizeData(data.integrated, verbose = FALSE)
data.integrated <- ScaleData(data.integrated, verbose = FALSE)

#UMAP Fig.1E
cols<-c("#9C00E7", "#ce1256", "#8c96c6", "#fbb4b9")
DimPlot(data.integrated, reduction = "umap", label = TRUE, repel = TRUE, cols = cols, pt.size = 1.5)

#FeaturePlots Fig. 1F
FeaturePlot(data.integrated, features = c("TNNT2", "ACTN2"), cols = rev(brewer.pal(n = 11, name = "RdGy")))

#ViolinPlots Fig. 1G
VlnPlot(object = data.integrated, features =c("TNNT2", "TBX2", "TBX3", "BMP2", "NPPA", "HEY1"), cols=cols)

#create Heatmap Fig. 1H
all.markers <- FindAllMarkers(object = data.integrated, only.pos = TRUE, min.pct = 0.25,  thresh.use = 0.25)
top10 <- all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(object = data.integrated, features  = top10$gene, label = TRUE, group.colors = cols)+ scale_fill_gradientn(colors =rev(brewer.pal(n = 11, name = "RdBu")))


#annotate AVCM in D19 dataset as shown in Fig. S4F
#import dataset containing all 4 subtypes and annotate clustering identified in this analysis
#D19<-readRDS("subset.rds")
D19@meta.data[rownames(data.integrated@meta.data), "AVCM.clusters"] <- data.integrated@meta.data$seurat_clusters

DimPlot(D19, group.by = "AVCM.clusters", label = TRUE, repel = TRUE, cols = cols, pt.size = 1)


#obtain DE list for all 4 clusters (Supplement Table 2)
#markers 0
markers_0 <- FindMarkers(AVCM, ident.1 = "0", ident.2 = NULL, only.pos = TRUE) 
markers_0[
  with(markers_0, order( avg_log2FC, decreasing = TRUE)),
]

write.csv(markers_0, file="/path/to/directory/DE_list_AVCMD19/markers_0.csv", row.names = TRUE)

## markers 1
markers_1 <- FindMarkers(AVCM, ident.1 = "1", ident.2 = NULL, only.pos = TRUE) 
markers_1[
  with(markers_1, order( avg_log2FC, decreasing = TRUE)),
]

write.csv(markers_1, file="/path/to/directory/DE_list_AVCMD19/markers_1.csv", row.names = TRUE)


#marker 2
markers_2 <- FindMarkers(AVCM, ident.1 = "2", ident.2 = NULL, only.pos = TRUE) 
markers_2[
  with(markers_2, order( avg_log2FC, decreasing = TRUE)),
]

write.csv(markers_2, file="/path/to/directory/DE_list_AVCMD19/markers_2.csv", row.names = TRUE)


#markers 3
markers_3 <- FindMarkers(AVCM, ident.1 = "3", ident.2 = NULL, only.pos = TRUE) 
markers_3[
  with(markers_3, order( avg_log2FC, decreasing = TRUE)),
]

write.csv(markers_3, file="/path/to/directory/DE_list_AVCMD19/markers_3.csv", row.names = TRUE)


#### SessionInfo ####
writeVersions <- function(sessionDir="/path/to/directory/Documentation"){
  write(paste0("Bioconductor version ", capture.output(tools:::.BioC_version_associated_with_R_version()),"\n"), 
        paste0(sessionDir,"/sessionInfo.txt"))
  write(capture.output(sessionInfo()), paste0(sessionDir,"/sessionInfo.txt"), append=TRUE)
}

#writeLines(capture.output(sessionInfo()), "../Documentation/sessionInfo_secclust_AVCMD19.txt")
writeVersions()
