#projection of AVCM cells onto in vivo data using integration anchors
library(Seurat)
library(Matrix)
library(fossil) 
library(dplyr)
library(plyr)
library(ggplot2)

#import datasets (Hill et al., 2019; De Soysa et al., 2019) which were additionally annotated with CellType and timepoint
setwd("/path/to/directory/AVC_invivo/allCM/E9.25")
E9.25<-readRDS("E9.25_annotated_timepoint_CellType.rds")
setwd("/path/to/directory/AVC_invivo/allCM/E10.5")
E10.5<-readRDS("E10.5_annotated_timepoint_CellType.rds")
setwd("/path/to/directory/AVC_invivo/allCM/E13.5")
E13.5<-readRDS("E13.5_annotated_timepoint_CellType.rds")

#import hiPSC-AVCM dataset with converted gene names
setwd("L:/basic/Personal Archive/A/awiesinger/scRNAseq_analyses/AVNCMs/analysis/mouse_Orthologs")
hAVCM<-readRDS("hiPSC_mouseOrthologue_GeneNames.rds") #this dataset was generated as described in "Script_convertgenenames.R"


hAVCM@meta.data$timepoint<- "hiPSC"
hAVCM@meta.data$CellType <- "hiPSC"
head(hAVCM[[]])

#merge all datasets
data<-merge(x= E9.25, y= c(E10.5, E13.5, hAVCM), project = "AVC-CM")

data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^mt-")
head(data[[]])

# Visualize QC metrics as a violin plot
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "timepoint")

plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1+plot2

##SCTranform normalization + batch correction
data.list <- SplitObject(data, split.by = "timepoint")

##Prior to finding anchors, we perform standard preprocessing, 
##and identify variable features individually for each. Variable feature selection is based on a 
##variance stabilizing transformation ("vst")

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




DefaultAssay(data.integrated) <- "integrated"

#proceed with downstream analysis (i.e. visualization, clustering) on the integrated dataset
#do not scale integrated data

data.integrated <- RunPCA(data.integrated, verbose = FALSE)
DimPlot(data.integrated, reduction = "pca", dims = c(1,2))
DimHeatmap(data.integrated, dims = 1, cells = 500, balanced = TRUE)

ElbowPlot(data.integrated, ndims = 50)
data.integrated <- RunUMAP(data.integrated, dims = 1:20)
DimPlot(data.integrated, reduction = "umap")
data.integrated <- FindNeighbors(data.integrated, dims = 1:20)
data.integrated <- FindClusters(object = data.integrated, reduction = "pca",dims = 1:20, resolution = 0.6,random.seed = 2020)
head(Idents(data.integrated), 5)
data.integrated <- RunUMAP(data.integrated, dims = 1:20)
DimPlot(data.integrated, reduction = "umap", label = TRUE, repel = TRUE)

DimPlot(data.integrated, group.by = "timepoint", label = TRUE, repel = TRUE)
DimPlot(data.integrated, group.by = "AVC", label = TRUE, repel = TRUE)

DefaultAssay(data.integrated) <- "RNA"

#continue # Normalize RNA data for visualization purposes

data <- NormalizeData(data.integrated, verbose = FALSE)
data <- ScaleData(data, verbose = FALSE)

#import hiPSC-AVCMs with mouse Orthologue genes separately to be able to annotate hiPSC-AVCMs in dataset
setwd("L:/basic/Personal Archive/A/awiesinger/scRNAseq_analyses/AVNCMs/analysis/mouse_Orthologs")
hAVCM<-readRDS("hiPSC_mouseOrthologue_GeneNames.rds")
data@meta.data[rownames(hAVCM@meta.data), "AVC"] <- hAVCM@meta.data$seurat_clusters


#annotate clusters identified in hiPSC-AVCMs in entire dataset
data@meta.data[rownames(hAVCM@meta.data), "AVC"] <- hAVCM@meta.data$seurat_clusters

DimPlot(data, group.by = "AVC", label = TRUE, repel = TRUE)
DimPlot(data, group.by = "CellType", label = TRUE, repel = TRUE)

FeaturePlot(data, c("Tnnt2", "Actn2"))
Idents(data)<-"seurat_clusters"
subset<-subset(data, idents = c("5", "8","13", "14", "15", "16", "17", "18", "19", "20"), invert=TRUE)
DimPlot(subset, label = TRUE, repel = TRUE)

all.markers <- FindAllMarkers(object = subset, only.pos = TRUE, min.pct = 0.25,  thresh.use = 0.25)
top10 <- all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
library(RColorBrewer)
DoHeatmap(object = subset, features  = top10$gene, label = T) + scale_fill_gradientn(colors =rev(brewer.pal(n = 11, name = "RdBu")))

#annotate clusters
current.cluster.ids <- c(0,1,2,3,4,6,7,9,10,11,12)
new.cluster.ids <- c("VCM", "VCM", "VCM", "ACM", "AVCM", "OFT", "SV_SANCM", "VCM", "VCM", "ACM", "VCM")
subset@active.ident <- plyr::mapvalues(x = subset@active.ident, from = current.cluster.ids, to = new.cluster.ids)
#and ordered
my_levels <- c("SV_SANCM", "ACM", "AVCM", "VCM", "OFT")
subset@active.ident <- factor(x = subset@active.ident, levels = my_levels)

DimPlot(subset, reduction = "umap", label = TRUE, repel = TRUE)


#get DE lists (Supplement Table 3)
#AVCM
markers_AVCM <- FindMarkers(subset, ident.1 = "AVCM", ident.2 = NULL , only.pos = TRUE) 
markers_AVCM[
  with(markers_AVCM, order(avg_log2FC, decreasing = TRUE)),
]

write.csv(markers_AVCM, file="/path/to/directory/markers_AVCM.csv", row.names = TRUE)

#SANCM
markers_SV_SANCM <- FindMarkers(subset, ident.1 = "SV_SANCM", ident.2 = NULL , only.pos = TRUE) 
markers_SV_SANCM[
  with(markers_SV_SANCM, order(avg_log2FC, decreasing = TRUE)),
]

write.csv(markers_SV_SANCM, file="/path/to/directory/markers_SV_SANCM.csv", row.names = TRUE)

#ACM
markers_ACM <- FindMarkers(subset, ident.1 = "ACM", ident.2 = NULL , only.pos = TRUE) 
markers_ACM[
  with(markers_ACM, order(avg_log2FC, decreasing = TRUE)),
]

write.csv(markers_ACM, file="/path/to/directory/markers_ACM.csv", row.names = TRUE)

#VCM
markers_VCM <- FindMarkers(subset, ident.1 = "VCM", ident.2 = NULL , only.pos = TRUE) 
markers_VCM[
  with(markers_VCM, order(avg_log2FC, decreasing = TRUE)),
]

write.csv(markers_VCM, file="/path/to/directory/markers_VCM.csv", row.names = TRUE)

#OFT
markers_OFT <- FindMarkers(subset, ident.1 = "OFT", ident.2 = NULL , only.pos = TRUE) 
markers_OFT[
  with(markers_OFT, order(avg_log2FC, decreasing = TRUE)),
]

write.csv(markers_OFT, file="/path/to/directory/markers_OFT.csv", row.names = TRUE)



#figures for paper
#Fig. 2A
cols=c( "#abd9e9", "#d7191c","#fdae61", "#2c7bb6", "#dfc27d")
DimPlot(subset, reduction = "umap", label = TRUE, repel = TRUE, cols = cols, pt.size = 0.7)

#Fig. 2B
#AVCM cluster
FeaturePlot(subset, c("Tbx2", "Tbx3", "Bmp2", "Hcn4"), cols= c("lightgrey", "#e66101"), pt.size = 0.7) 

#Fig. S4B
#SV_SANCM cluster
FeaturePlot(subset, c("Shox2", "Vsnl1", "Tbx18"), cols= c("lightgrey", "#2b8cbe"), pt.size = 0.7) 

#Atrial cluster
FeaturePlot(subset, c("Nppa", "Nr2f2", "Hey1"), cols= c("lightgrey", "#d7191c"), pt.size = 0.7) 

#ventricular clusters
FeaturePlot(subset, c("Myl2", "Myh7","Hey2"), cols= c("lightgrey","#2c7bb6"), pt.size = 0.7) 

#OFT cluster
FeaturePlot(subset, c("Rspo3", "Isl1", "Bmp4"), cols= c("lightgrey","#fe9929"), pt.size = 0.7) 

#Fig. 2C
#identify AVC clusters
Idents(data)<-"timepoint"
hiPSC <- WhichCells(data, idents = c("hiPSC"))
DimPlot(subset, label=T, repel=TRUE, cells.highlight= hiPSC, pt.size = 0.9)

#extract cells
Idents(subset)<-"seurat_clusters"
AVC<-WhichCells(subset, ident="4")
VCM <- WhichCells(subset, idents = c("0", "1", "2", "9", "10", "12"))
OFT <- WhichCells(subset, idents = "6")
ACM <- WhichCells(subset, idents = c("3", "11"))
SV_SANCM <- WhichCells(subset, idents = c("7"))

VCM<-subset(subset, idents= c("0", "1", "2", "9", "10", "12"))
OFT<-subset(subset, idents= "6")
ACM<-subset(subset, idents= c("3", "11"))
SV_SANCM<-subset(subset, idents= "7")
AVC<-subset(subset, idents= "4")

VCM@meta.data$annotation<-"VCM"
OFT@meta.data$annotation<-"OFT"
ACM@meta.data$annotation<-"ACM"
SV_SANCM@meta.data$annotation<-"SV_SANCM"
AVC@meta.data$annotation<-"AVC"

#annotate those cells in entire dataset

subset@meta.data[rownames(AVC@meta.data), "annotation"] <- AVC@meta.data$annotation
subset@meta.data[rownames(VCM@meta.data), "annotation"] <- VCM@meta.data$annotation
subset@meta.data[rownames(OFT@meta.data), "annotation"] <- OFT@meta.data$annotation
subset@meta.data[rownames(ACM@meta.data), "annotation"] <- ACM@meta.data$annotation
subset@meta.data[rownames(SV_SANCM@meta.data), "annotation"] <- SV_SANCM@meta.data$annotation

Idents(subset)<-"annotation"

# Heatmap Fig. S4A
all.markers <- FindAllMarkers(object = subset, only.pos = TRUE, min.pct = 0.25,  thresh.use = 0.25)
top10 <- all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(object = subset, features  = top10$gene, label = T) + scale_fill_gradientn(colors =rev(brewer.pal(n = 11, name = "RdBu")))


#### SessionInfo ####
writeVersions <- function(sessionDir="/path/to/directory/Documentation"){
  write(paste0("Bioconductor version ", capture.output(tools:::.BioC_version_associated_with_R_version()),"\n"), 
        paste0(sessionDir,"/sessionInfo_integration_invivo_invitro.txt"))
  write(capture.output(sessionInfo()), paste0(sessionDir,"/sessionInfo.txt"), append=TRUE)
}

#writeLines(capture.output(sessionInfo()), "../Documentation/sessionInfo_integration_invivo_invitro.txt")
writeVersions()
