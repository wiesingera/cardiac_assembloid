#projection of AVCM cells onto in vivo data using integration anchors
library(Seurat)
library(Matrix)
library(fossil) 
library(dplyr)
library(plyr)
library(ggplot2)

#import datasets (Hill et al., 2019; De Soysa et al., 2019) which were additionally annotated with CellType and timepoint
setwd("L:/basic/Personal Archive/A/awiesinger/scRNAseq_analyses/AVNCMs/AVC_invivo/allCM/E9.25")
E9.25<-readRDS("E9.25_annotated_timepoint_CellType.rds")
setwd("L:/basic/Personal Archive/A/awiesinger/scRNAseq_analyses/AVNCMs/AVC_invivo/allCM/E10.5")
E10.5<-readRDS("E10.5_annotated_timepoint_CellType.rds")
setwd("L:/basic/Personal Archive/A/awiesinger/scRNAseq_analyses/AVNCMs/AVC_invivo/allCM/E13.5")
E13.5<-readRDS("E13.5_annotated_timepoint_CellType.rds")

setwd("L:/basic/Personal Archive/A/awiesinger/scRNAseq_analyses/AVNCMs/analysis/mouse_Orthologs")
hAVCM<-readRDS("data_for_projection.rds") #this dataset was generated as described in "Script_convertgenenames.R"


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

#annotate clusters identified in hiPSC-AVCMs in entire dataset
data@meta.data[rownames(hAVCM@meta.data), "AVC"] <- hAVCM@meta.data$seurat_clusters

DimPlot(data, group.by = "AVC", label = TRUE, repel = TRUE)
DimPlot(data, group.by = "CellType", label = TRUE, repel = TRUE)

FeaturePlot(data, c("Tnnt2", "Actn2"))
Idents(data)<-"seurat_clusters"
subset<-subset(data, idents = c("8", "17", "18", "19", "20"), invert=TRUE)
DimPlot(subset, label = TRUE, repel = TRUE)

all.markers <- FindAllMarkers(object = subset, only.pos = TRUE, min.pct = 0.25,  thresh.use = 0.25)
top10 <- all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
library(RColorBrewer)
DoHeatmap(object = subset, features  = top10$gene, label = T) + scale_fill_gradientn(colors =rev(brewer.pal(n = 11, name = "RdBu")))

#remove more unwanted clusters
subset<-subset(subset, idents = c("13", "14", "15", "16"), invert=TRUE)
DimPlot(subset, label = TRUE, repel = TRUE)

all.markers <- FindAllMarkers(object = subset, only.pos = TRUE, min.pct = 0.25,  thresh.use = 0.25)
top10 <- all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
library(RColorBrewer)
DoHeatmap(object = subset, features  = top10$gene, label = T) + scale_fill_gradientn(colors =rev(brewer.pal(n = 11, name = "RdBu")))


#identify AVC clusters
Idents(subset)<-"timepoint"
hiPSC <- WhichCells(subset, idents = c("hiPSC"))
DimPlot(subset, label=T, cells.highlight= hiPSC)

Idents(subset)<-"AVC"
AVC_E10.5 <- WhichCells(subset, idents = "11")
AVC_13.5 <- WhichCells(subset, idents = "9")
AVC_E9.25 <- WhichCells(subset, idents = c("8"))
AVC_3 <- WhichCells(subset, idents = c("3"))
AVC_4 <- WhichCells(subset, idents = c("4"))

DimPlot(subset, label=T, cells.highlight= AVC_E10.5)
DimPlot(subset, label=T, cells.highlight= AVC_13.5)
DimPlot(subset, label=T, cells.highlight= AVC_E9.25)
DimPlot(subset, label=T, cells.highlight= AVC_3)
DimPlot(subset, label=T, cells.highlight= AVC_4)

#remove non-CM from hiPSC dataset

subset<-subset(subset, cells= AVC_3, invert=TRUE)
subset<-subset(subset, cells= AVC_4, invert=TRUE)

##cluster 5 appears as a smooth muscle cell cluster and will be removed
Idents(subset)<- "seurat_clusters"
subset<-subset(subset, idents="5", invert=TRUE)
DimPlot(subset, reduction = "umap", label = TRUE, repel = TRUE)

#annotate clusters
current.cluster.ids <- c(0,1,2,3,4,6,7,9,10,11,12)
new.cluster.ids <- c("VCM", "VCM", "VCM", "ACM", "AVCM", "OFT", "SV_SANCM", "VCM", "VCM", "ACM", "VCM")
subset@active.ident <- plyr::mapvalues(x = subset@active.ident, from = current.cluster.ids, to = new.cluster.ids)
#and ordered
my_levels <- c("SV_SANCM", "ACM", "AVCM", "VCM", "OFT")
subset@active.ident <- factor(x = subset@active.ident, levels = my_levels)

DimPlot(subset, reduction = "umap", label = TRUE, repel = TRUE)

#get DE lists
#AVCM
markers_AVCM <- FindMarkers(subset, ident.1 = "AVCM", ident.2 = NULL , only.pos = TRUE) 
markers_AVCM[
  with(markers_AVCM, order(avg_log2FC, decreasing = TRUE)),
]

write.csv(markers_AVCM, file="L:/basic/Personal Archive/A/awiesinger/scRNAseq_analyses/AVNCMs/AVC_invivo/allCM/integration/DElist_annotation/markers_AVCM.csv", row.names = TRUE)

#SANCM
markers_SV_SANCM <- FindMarkers(subset, ident.1 = "SV_SANCM", ident.2 = NULL , only.pos = TRUE) 
markers_SV_SANCM[
  with(markers_SV_SANCM, order(avg_log2FC, decreasing = TRUE)),
]

write.csv(markers_SV_SANCM, file="L:/basic/Personal Archive/A/awiesinger/scRNAseq_analyses/AVNCMs/AVC_invivo/allCM/integration/DElist_annotation/markers_SV_SANCM.csv", row.names = TRUE)

#ACM
markers_ACM <- FindMarkers(subset, ident.1 = "ACM", ident.2 = NULL , only.pos = TRUE) 
markers_ACM[
  with(markers_ACM, order(avg_log2FC, decreasing = TRUE)),
]

write.csv(markers_ACM, file="L:/basic/Personal Archive/A/awiesinger/scRNAseq_analyses/AVNCMs/AVC_invivo/allCM/integration/DElist_annotation/markers_ACM.csv", row.names = TRUE)

#VCM
markers_VCM <- FindMarkers(subset, ident.1 = "VCM", ident.2 = NULL , only.pos = TRUE) 
markers_VCM[
  with(markers_VCM, order(avg_log2FC, decreasing = TRUE)),
]

write.csv(markers_VCM, file="L:/basic/Personal Archive/A/awiesinger/scRNAseq_analyses/AVNCMs/AVC_invivo/allCM/integration/DElist_annotation/markers_VCM.csv", row.names = TRUE)

#OFT
markers_OFT <- FindMarkers(subset, ident.1 = "OFT", ident.2 = NULL , only.pos = TRUE) 
markers_OFT[
  with(markers_OFT, order(avg_log2FC, decreasing = TRUE)),
]

write.csv(markers_OFT, file="L:/basic/Personal Archive/A/awiesinger/scRNAseq_analyses/AVNCMs/AVC_invivo/allCM/integration/DElist_annotation/markers_OFT.csv", row.names = TRUE)



#figures for paper
#Fig. 2a
cols=c( "#abd9e9", "#d7191c","#fdae61", "#2c7bb6", "#dfc27d")
DimPlot(subset, reduction = "umap", label = TRUE, repel = TRUE, cols = cols, pt.size = 0.7)

#Fig. 2b
#AVCM cluster
FeaturePlot(subset, c("Tbx2", "Tbx3", "Bmp2", "Hcn4"), cols= c("lightgrey", "#e66101"), pt.size = 0.7) 

#Fig. S4b
#SV_SANCM cluster
FeaturePlot(subset, c("Shox2", "Vsnl1", "Tbx18"), cols= c("lightgrey", "#2b8cbe"), pt.size = 0.7) 

#Atrial cluster
FeaturePlot(subset, c("Nppa", "Nr2f2", "Hey1"), cols= c("lightgrey", "#d7191c"), pt.size = 0.7) 

#ventricular clusters
FeaturePlot(subset, c("Myl2", "Myh7","Hey2"), cols= c("lightgrey","#2c7bb6"), pt.size = 0.7) 

#OFT cluster
FeaturePlot(subset, c("Rspo3", "Isl1", "Bmp4"), cols= c("lightgrey","#fe9929"), pt.size = 0.7) 

#Fig. 2c
DimPlot(subset, label=T, repel=TRUE, cells.highlight= hiPSC, pt.size = 0.9)

#extract cells
Idents(subset)<-"seurat_clusters"
AVC<-WhichCells(subset, ident="4")
VCM <- WhichCells(subset, idents = c("0", "1", "2", "9", "10", "12"))
OFT <- WhichCells(subset, idents = "6")
ACM <- WhichCells(subset, idents = c("3", "11"))
SV_SANCM <- WhichCells(subset, idents = c("7"))

#check whether correct cells were selected
DimPlot(subset, label=T, cells.highlight= VCM)
DimPlot(subset, label=T, cells.highlight= OFT)
DimPlot(subset, label=T, cells.highlight= ACM)
DimPlot(subset, label=T, cells.highlight= SV_SANCM)
DimPlot(subset, label=T, cells.highlight= AVC)

VCM<-subset(subset, idents= c("0", "1", "2", "9", "10", "12"))
OFT<-subset(subset, idents= "6")
ACM<-subset(subset, idents= c("3", "11"))
SV_SANCM<-subset(subset, idents= "7")
AVC<-subset(subset, idents= "4")

DimPlot(VCM, label = TRUE, repel = TRUE)
DimPlot(ACM, label = TRUE, repel = TRUE)
DimPlot(OFT, label = TRUE, repel = TRUE)
DimPlot(SV_SANCM, label = TRUE, repel = TRUE)
DimPlot(AVC, label = TRUE, repel = TRUE)

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

head(subset[[]])
tail(subset[[]])

Idents(subset)<-"annotation"

# Fig. S2a
all.markers <- FindAllMarkers(object = subset, only.pos = TRUE, min.pct = 0.25,  thresh.use = 0.25)
top10 <- all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(object = subset, features  = top10$gene, label = T) + scale_fill_gradientn(colors =rev(brewer.pal(n = 11, name = "RdBu")))


# Fig. 2f
Idents(subset)<- "CellType"
subset_invivo <- subset(subset, idents="hiPSC", invert=TRUE)

Idents(subset_invivo)<- "annotation"
features = c("Syne2", "Msi2",	"Rcsd1",	"Ltbp1",	"Bambi",	"Palld",	"Mef2a",	"Unc5b",	"Mical2",
             "Atp1a1",	"Ppp1r14c",	"Atp1b1",	"Fbn2",	"Nebl",	"Rspo3",	"Cpne5",	"Msx2",	"Bmp2",	"Tbx3",	"Tbx2")

DotPlot(subset_invivo, features = features, cols = c("lightgrey", "darkorange")) + coord_flip()

#### SessionInfo ####
writeVersions <- function(sessionDir="L:/basic/Personal Archive/A/awiesinger/scripts and data for github/Documentation"){
  write(paste0("Bioconductor version ", capture.output(tools:::.BioC_version_associated_with_R_version()),"\n"), 
        paste0(sessionDir,"/sessionInfo.txt"))
  write(capture.output(sessionInfo()), paste0(sessionDir,"/sessionInfo.txt"), append=TRUE)
}

#writeLines(capture.output(sessionInfo()), "../Documentation/sessionInfo_integration_invivo_invitro.txt")
writeVersions()
