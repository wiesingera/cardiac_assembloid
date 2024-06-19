## cluster analysis of hiPSC-SANCM, -VCM, -ACM and W+R (AVCM)

library(Seurat)
library(URD)
library(Matrix)
library(fossil) 
library(dplyr)
library(plyr)
library(RColorBrewer)
library(babelgene)

setwd("/path/to/directory/")

data<-readRDS("00_merged_D19reps_SAN_AM_VM_AVCM.rds")

## vizualize QC
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")

head(data@meta.data)

VlnPlot(data, "nFeature_RNA")

data <-subset(data, subset = nFeature_RNA >1000 & nFeature_RNA < 8000 & percent.mt < 50)

rownames<- rownames(data)
new.rownames<- gsub("--chr.*","", rownames)
head(new.rownames)

rownames(data@assays$RNA@counts)<- new.rownames
rownames(data@assays$RNA@data)<- new.rownames

data
#An object of class Seurat 
#18447 features across 2412 samples within 1 assay 
#Active assay: RNA (18447 features, 0 variable features)

##SCTransform + batch correction

data.list <- SplitObject(data, split.by = "batch")

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
DimHeatmap(data.integrated, dims = 10:25, cells = 500, balanced = TRUE)

ElbowPlot(data.integrated, ndims = 60)
data.integrated <- RunUMAP(data.integrated, dims = 1:25)
DimPlot(data.integrated, reduction = "umap")
data.integrated <- FindNeighbors(data.integrated, dims = 1:25)
data.integrated <- FindClusters(object = data.integrated, reduction = "pca",dims = 1:25, resolution = 0.6,random.seed = 2020)
head(Idents(data.integrated), 5)
data.integrated <- RunUMAP(data.integrated, dims = 1:25)
DimPlot(data.integrated, reduction = "umap", label = TRUE, repel = TRUE)
DimPlot(data.integrated, group.by = "orig.ident", label = TRUE, repel = TRUE)

DefaultAssay(data.integrated) <- "RNA"

# Normalize RNA data for visualization purposes

data.integrated <- NormalizeData(data.integrated, verbose = FALSE)
data.integrated <- ScaleData(data.integrated, verbose = FALSE)

# explore clusters
all.markers <- FindAllMarkers(object = data.integrated, only.pos = TRUE, min.pct = 0.25,  thresh.use = 0.25)

top10 <- all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(object = data.integrated, features  = top10$gene, label = TRUE)

# cluster 8 and cluster 13 show very high spike-in counts
VlnPlot(data.integrated, "nFeature_RNA") #cluster 8 and 13 also show very low feature counts
plot1 <- FeatureScatter(data.integrated, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data.integrated, feature1 = "nFeature_RNA", feature2 = "nCount_RNA")
plot1
plot2

#remove cluster 8 and 13 (this dataset was used to annotate AVCM subclustering named D19 as shown in Fig. S3e)
subset <- subset(data.integrated, ident = c("0", "1", "2", "3", "4", "5", "6", "7", "9", "10", "11", "12"))

#visualize subset
DimPlot(subset, group.by = "orig.ident", label = TRUE, repel = TRUE)
DimPlot(subset, reduction = "umap", label = TRUE, repel = TRUE)

#Fig. 1C 
cols<-c("#F77DD9", "#33AFF2", "#9C00E7", "#FFAF6F", "#31a354", "#8856a7", "#fec44f", "#7fcdbb", "#43a2ca", "#c994c7", "#de2d26", "#dd1c77")
DimPlot(subset, reduction = "umap", label = TRUE, repel = TRUE, cols = cols, pt.size = 0.7)

#export DE genes for each cluster (Supplement Table 1)
##markers 0
markers_0 <- FindMarkers(subset, ident.1 = "0", ident.2 = NULL, only.pos = TRUE) 
markers_0[
  with(markers_0, order( avg_log2FC, decreasing = TRUE)),
]

write.csv(markers_0, file="/path/to/directory/DE_list/markers_0.csv", row.names = TRUE)


## markers 1
markers_1 <- FindMarkers(subset, ident.1 = "1", ident.2 = NULL, only.pos = TRUE) 
markers_1[
  with(markers_1, order( avg_log2FC, decreasing = TRUE)),
]

write.csv(markers_1, file="/path/to/directory/DE_list/markers_1.csv", row.names = TRUE)


#marker 2
markers_2 <- FindMarkers(subset, ident.1 = "2", ident.2 = NULL, only.pos = TRUE) 
markers_2[
  with(markers_2, order( avg_log2FC, decreasing = TRUE)),
]

write.csv(markers_2, file="/path/to/directory/DE_list/markers_2.csv", row.names = TRUE)


#markers 3
markers_3 <- FindMarkers(subset, ident.1 = "3", ident.2 = NULL, only.pos = TRUE) 
markers_3[
  with(markers_3, order( avg_log2FC, decreasing = TRUE)),
]

write.csv(markers_3, file="/path/to/directory/DE_list/markers_3.csv", row.names = TRUE)


#markers 4
markers_4 <- FindMarkers(subset, ident.1 = "4", ident.2 = NULL, only.pos = TRUE) 
markers_4[
  with(markers_4, order( avg_log2FC, decreasing = TRUE)),
]

write.csv(markers_4, file="/path/to/directory/DE_list/markers_4.csv", row.names = TRUE)


#markers 5
markers_5 <- FindMarkers(subset, ident.1 = "5", ident.2 = NULL, only.pos = TRUE) 
markers_5[
  with(markers_5, order( avg_log2FC, decreasing = TRUE)),
]

write.csv(markers_5, file="/path/to/directory/DE_list/markers_5.csv", row.names = TRUE)


#markers 6
markers_6 <- FindMarkers(subset, ident.1 = "6", ident.2 = NULL, only.pos = TRUE) 
markers_6[
  with(markers_6, order( avg_log2FC, decreasing = TRUE)),
]

write.csv(markers_6, file="/path/to/directory/DE_list/markers_6.csv", row.names = TRUE)


#markers 7
markers_7 <- FindMarkers(subset, ident.1 = "7", ident.2 = NULL, only.pos = TRUE) 
markers_7[
  with(markers_7, order( avg_log2FC, decreasing = TRUE)),
]

write.csv(markers_7, file="/path/to/directory/DE_list/markers_7.csv", row.names = TRUE)


#markers 9
markers_9 <- FindMarkers(subset, ident.1 = "9", ident.2 = NULL, only.pos = TRUE) 
markers_9[
  with(markers_9, order( avg_log2FC, decreasing = TRUE)),
]

write.csv(markers_9, file="/path/to/directory/DE_list/markers_9.csv", row.names = TRUE)


#cluster 10 
markers_10 <- FindMarkers(subset, ident.1 = "10", ident.2 = NULL, only.pos = TRUE) 
markers_10[
  with(markers_10, order( avg_log2FC, decreasing = TRUE)),
]

write.csv(markers_10, file="/path/to/directory/DE_list/markers_10.csv", row.names = TRUE)


#cluster 11 
markers_11 <- FindMarkers(subset, ident.1 = "11", ident.2 = NULL, only.pos = TRUE) 
markers_11[
  with(markers_11, order( avg_log2FC, decreasing = TRUE)),
]

write.csv(markers_11, file="/path/to/directory/DE_list/markers_11.csv", row.names = TRUE)


#cluster 12 
markers_12 <- FindMarkers(subset, ident.1 = "12", ident.2 = NULL, only.pos = TRUE) 
markers_12[
  with(markers_12, order( avg_log2FC, decreasing = TRUE)),
]


##extract CMs only
#according to TNNT2 and ACTN2 expresseion, cluster 5, 7, 9, 10 and 12 are non-CMs
#Fig. 1D
FeaturePlot(subset, features = c("TNNT2", "ACTN2"), cols = rev(brewer.pal(n = 11, name = "RdGy")))
VlnPlot(subset, features=(c("TNNT2", "ACTN2")), cols= cols)

#subset for CMs only
subset<-subset(subset, idents = c("0", "1", "2", "3", "4", "6"))

#Fig. S3A
cols<-c("#F77DD9", "#33AFF2", "#9C00E7", "#FFAF6F", "#31a354", "#fec44f")
DimPlot(subset, reduction = "umap", label = TRUE, repel = TRUE, cols = cols, pt.size = 0.9)

#Fig. S3B
FeaturePlot(subset, features = c("TNNT2", "ACTN2"), cols = rev(brewer.pal(n = 11, name = "RdGy")))

#Fig. S3C-H
FeaturePlot(subset, features = c("TBX2", "TBX3", "MSX2", "RSPO3")) & scale_colour_gradient(low = "lightgrey", high = "purple")
FeaturePlot(subset, features = c("SHOX2", "ISL1"))& scale_colour_gradient(low = "lightgrey", high = "steelblue3")
FeaturePlot(subset, features = c("TCF21", "COL3A1")) & scale_colour_gradient(low = "lightgrey", high = "darkgreen")
FeaturePlot(subset, features = c("NPPA", "NR2F2")) & scale_colour_gradient(low = "lightgrey", high = "maroon2")
FeaturePlot(subset, features = c("MYL2", "HOPX")) & scale_colour_gradient(low = "lightgrey", high = "orangered2")
FeaturePlot(subset, features = c("HAPLN1", "PITX2")) & scale_colour_gradient(low = "lightgrey", high = "orange")

##########################################################3
##### Module Scores Fig.2

#Get top 50 genes enriched in mouse AVC cells Fig. 2E
#marker_AVCM retrieval is described in script: integration_invivo_invitro.R
top50 <- markers_AVCM %>% top_n(n = 50, wt = avg_log2FC)

head(top50)

##convert to human orthologues
library(babelgene)

genes<-rownames(top50)
orthologs <-orthologs(genes = genes, species = "mouse", human = FALSE)
TOP50<-orthologs$human_symbol

#add module score of TOP 50 mouse AVC genes to our dataset
subset <- AddModuleScore(subset,
                         features = list(TOP50),
                         name="mouseAVC_enriched")
# Plot scores
FeaturePlot(subset,
            features = "mouseAVC_enriched1", label = TRUE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))


#Get top 50 genes enriched in human fetal AVN (Farah 2024) Fig. 2F
AVN_enriched <-c("RSPO3",	"BMP2",	"GNG11",	"IGFBP5",	"TBX3",	"MSX2",	"MB",	"HAMP",	"FHL2",	"ACTA2",
                 "PRSS35",	"SFRP1",	"ID2",	"ID3",	"CXCL12",	"DOK4",	"TSC22D1",	"CRHBP",	"SPEG",	"TSTD1",
                 "RRAD",	"EMC10",	"MYH6",	"ID1",	"COL18A1",	"EPHA4",	"VCAM1",	"LMO4",	"BAMBI",	"IGFBP6",
                 "PRSS23",	"CPNE5",	"TIMP3",	"SLC2A3",	"FAM189A2",	"CKB",	"CRABP2",	"ADAMTS19",	"PKIG",	"JUNB",
                 "NME4",	"ADM",	"ID4",	"SAA1",	"MARCKS",	"ZFP36L1",	"JUN",	"ITGA10",	"SMAD6",	"ATP1A1")
#add module score
subset <- AddModuleScore(subset,
                         features = list(AVN_enriched),
                         name="hAVNfetal_enriched")

# Plot scores
FeaturePlot(subset,
            features = "hAVNfetal_enriched1", label = TRUE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

#Get top 50 genes enriched in human AVN (Kanemaru 2023) Fig. 2F
AVN_enriched <-c("AL358335.2",	"RSPO3",	"SLITRK2",	"PTGER3",	"BMP2",	"SEMA3D",	"MOXD1",	"PDZRN4",	"MYL2",	"GNG11",
                 "HMGCLL1",	"GPR158",	"KCND2",	"NPNT",	"KEL",	"DPYSL5",	"ARHGAP6",	"AL356234.2",	"MPPED2",	"SPOCK3",
                 "ACTA1",	"MGARP",	"MSX2",	"GADL1",	"AC026765.3",	"CGNL1",	"PDE1A",	"KDR",	"NLGN1",	"LRRC4C",
                 "PLXNA2",	"ATP1B4",	"CRNDE",	"TBX3",	"SNAP91",	"APOA1",	"ANO4",	"KLHL4",	"LRFN5",	"SEMA3E",
                 "LINC02284",	"NEGR1",	"MMP23B",	"EDIL3",	"TNFRSF11B",	"MAP2",	"LRRC38",	"ADAMTS5",	"UBASH3B",	"PCDH20")
#add module score
subset <- AddModuleScore(subset,
                         features = list(AVN_enriched),
                         name="AVN_enriched")

#Plot scores
FeaturePlot(subset,
            features = "AVN_enriched1", label = TRUE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))


#Get top 50 genes enriched in PSC-AVC cells (Schmidt 2023) Fig. S4E
#Schmidt data can be retrieved from the following link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE239890

markers_AVC <- FindMarkers(Schmidtdata_seurat_object, ident.1 = "AVC", ident.2 = NULL, only.pos = TRUE) 
markers_AVC[
  with(markers_AVC, order( avg_log2FC, decreasing = TRUE)),
]
top50 <- markers_SchmidtAVCM %>% top_n(n = 50, wt = avg_log2FC)

head(top50)

#add module score of TOP 50 Schmidt AVC genes to our dataset
subset <- AddModuleScore(subset,
                         features = list(rownames(top50)),
                         name="SchmidtAVC_enriched")
# Plot scores
FeaturePlot(subset,
            features = "SchmidtAVC_enriched1", label = TRUE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

#### SessionInfo ####
writeVersions <- function(sessionDir="/path/to/directory/Documentation"){
  write(paste0("Bioconductor version ", capture.output(tools:::.BioC_version_associated_with_R_version()),"\n"), 
        paste0(sessionDir,"/sessionInfo.txt"))
  write(capture.output(sessionInfo()), paste0(sessionDir,"/sessionInfo.txt"), append=TRUE)
}

#writeLines(capture.output(sessionInfo()), "../Documentation/sessionInfo_D19_all4subtypes.txt")
writeVersions()
