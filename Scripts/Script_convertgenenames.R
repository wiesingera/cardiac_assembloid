library(babelgene) #https://cran.r-project.org/web/packages/babelgene/vignettes/babelgene-intro.html

AVCM1<-read.delim(file="L:/basic/Personal Archive/A/awiesinger/scRNAseq_analyses/AVNCMs/counts/AMC-HD-s071_HH2T2BGXH_S5_R2.TranscriptCounts.tsv",header=TRUE,row.names = 1, as.is = TRUE,sep="\t")
AVCM2<-read.delim(file="L:/basic/Personal Archive/A/awiesinger/scRNAseq_analyses/AVNCMs/counts/AMC-HD-s072_HH2T2BGXH_S6_R2.TranscriptCounts.tsv",header=TRUE,row.names = 1, as.is = TRUE,sep="\t")


##seems to work?
AVCM1<-CreateSeuratObject(counts = AVCM1, min.cells=3,min.features = 200, project = "AVCM1", assay = "RNA")
AVCM2<-CreateSeuratObject(counts = AVCM2, min.cells=3,min.features = 200, project = "AVCM2", assay = "RNA")

AVN<- merge(AVCM1, y = AVCM2, add.cell.ids = c("1", "2"), project = "scRNAseq")


rownames<- rownames(AVN)
new.rownames<- gsub("--chr.*","", rownames)
head(new.rownames)

rownames(AVN@assays$RNA@counts)<- new.rownames
rownames(AVN@assays$RNA@data)<- new.rownames

genes<-rownames(AVN@assays$RNA@counts)
musOrthologs<-orthologs(genes = genes, species = "mouse")
head(musOrthologs)

#convert back to human for sorting
hgOrthologs<-musOrthologs$human_symbol
head(hgOrthologs)

AVN2<-subset(AVN, features = hgOrthologs)
dim(AVN2)
#[1] 12835   612
dim(AVN)
#[1] 15921   612

genes<-rownames(AVN2)
musOrthologs<-orthologs(genes = genes, species = "mouse")
head(musOrthologs)

# Transfer clusterings to main object (in new column D19)
musOrthologs2<-subset(musOrthologs,!duplicated(musOrthologs$human_symbol))
#rownames(musOrthologs) <- musOrthologs$human_symbol

#AVN2@assays$RNA@counts[rownames(musOrthologs$human_symbol)] <-musOrthologs$symbol


AVN2@assays$RNA@counts@Dimnames[[1]] <- musOrthologs2$symbol
AVN2@assays$RNA@data@Dimnames[[1]] <- musOrthologs2$symbol

head(AVN2@assays$RNA@counts)
head(AVN2@assays$RNA@data)

setwd("/path/to/directory/mouse_orthologues")
saveRDS(AVN2, "hiPSC_mouseOrthologue_GeneNames.rds")

##continue analysis as in Script_integration_invivo_invitro
