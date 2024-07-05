### To include in revision
### To include in revision
### To include in revision



if (!requireNamespace("BSgenome.Mmusculus.UCSC.mm10", quietly = TRUE))
  BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")

if (!requireNamespace("EnsDb.Mmusculus.v79", quietly = TRUE))
  BiocManager::install("EnsDb.Mmusculus.v79")
library(EnsDb.Mmusculus.v79)


library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(EnsDb.Mmusculus.v79)

library(stringr)
library(Seurat)
library(Signac)
library(ggplot2)
library(dplyr)

#86k cells
adata <- readRDS("/corgi/cellbuster/argelaguet2023/bab/ftp1.babraham.ac.uk/data/processed/rna/seurat.rds")
#adata$sample <- str_split_fixed(colnames(adata),"#",2)[,1]
#adata$stage <- str_split_fixed(adata$sample,"_",2)[,1]

#Add metadata
cellmeta <- read.csv("/corgi/cellbuster/argelaguet2023/preprocessed/GSE205117_cell_metadata.txt.gz",sep="\t")
rownames(cellmeta) <- cellmeta$cell
cellmeta <- cellmeta[colnames(adata),]
adata@meta.data <- cbind(adata@meta.data, cellmeta[,-(1:2)])




########################################################################
################# RNAseq dim red & clustering ##########################
########################################################################


#Random subset
adata <- adata[,sample(colnames(adata),10000)]

## Normalization
adata<- NormalizeData(adata,normalization.method = 'LogNormalize')

#cell cycle scoring
if(organism=="human"){
  s.genes <- cc.genes.updated.2019$s.genes
  g2m.genes <- cc.genes.updated.2019$g2m.genes
}
if(organism=="mouse"){
  ccgenes <- read.csv("/corgi/cellbuster/py/cc_mouse.csv")
  s.genes <- ccgenes$symbol[ccgenes$ccstage=="s"]
  g2m.genes <- ccgenes$symbol[ccgenes$ccstage=="g2m"]
}
adata <- CellCycleScoring(adata, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
ccgenes <- rbind(data.frame(gene=s.genes, cc="s"),data.frame(gene=g2m.genes, cc="g2m"))

## HIGHLY VARIABLE FEATURE ANALYSIS. gets 850 with previous automatic method
adata <- FindVariableFeatures(adata, selection.method = 'vst',binning.method = 'equal_frequency', verbose = T, nfeatures = 1000)
adata <- ScaleData(adata, features = rownames(adata)) #default is to only use highly variable features






########################################################################
################# RNAseq dim red & clustering ##########################
########################################################################


DefaultAssay(adata) <- "RNA"

## Dimensional reduction
adata <- RunPCA(adata,reduction.name = 'PCA_RNA',ndims.print = c(1,3,5), verbose = T)
adata <- RunUMAP(object = adata, dims = 1:30, reduction = "PCA_RNA", reduction.name="UMAP_RNA")
adata <- RunUMAP(adata, dims = 1:30,reduction = 'PCA_RNA',n.components = 3L, reduction.name ='UMAP_RNA_3D' )
adata <- FindNeighbors(adata, dims = 1:30,reduction = 'PCA_RNA')





################################################################################
################### DE analysis ################################################
################################################################################


################### Compare RNA-based clusters - fine
DefaultAssay(adata) <- "RNA"
Idents(adata) <- adata$rna_clusters1
cluster_markers_rna1 <- FindAllMarkers(adata, only.pos = F, test.use = 'wilcox')
cluster_markers_rna1 %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) -> top10_rna





DimPlot(object = adata, reduction = 'UMAP_RNA', group.by = "stage", label = TRUE) #8.75 is suspicious
DimPlot(object = adata, reduction = 'UMAP_RNA', group.by = "sample", label = TRUE)
DimPlot(object = adata, reduction = 'UMAP_RNA', group.by = "Phase", label = TRUE)




################################################################################
################### add telo info ##############################################
################################################################################


## Barcode info for multiome 10x
bc10x <- data.frame(
  atac=read.table("/data/henlab/software/cellranger-arc-2.0.0/lib/python/atac/barcodes/737K-arc-v1.txt.gz")[,1],
  rna=read.table("/data/henlab/software/cellranger-arc-2.0.0/lib/python/cellranger/barcodes/737K-arc-v1.txt.gz")[,1]
)

mapping <- read.csv("/corgi/cellbuster/argelaguet2023/samplemeta_atac.csv",sep="\t")


sampleids <- unique(adata$sample)
all_telocnt <- list()
for(cur_sampleid in sampleids){
  print(cur_sampleid)

  #Read the corresponding telo file  
  dat <- read.csv(paste("/corgi/cellbuster/argelaguet2023/fromaws/",cur_sampleid,"_ATAC_R1.fastq.gz_telo.csv",sep=""))
  colnames(dat) <- c("c1","c2","c6","c7","atac")
  
  telocnt <- sqldf::sqldf("select atac, count(atac) as cnt from dat group by atac")
  telocnt$atac <- str_to_upper(sapply(telocnt$atac, spgs::reverseComplement))
  telocnt <- merge(telocnt, bc10x)
  telocnt$bc <- sprintf("%s#%s-1", cur_sampleid, telocnt$rna)
  
  all_telocnt[[cur_sampleid]] <- telocnt
}
telocnt <- do.call(rbind, all_telocnt)


#### Add to adata
rownames(telocnt) <- telocnt$bc
thecnt <- telocnt[colnames(adata),]$cnt
thecnt[is.na(thecnt)] <- 0
adata$telocnt <- thecnt
adata$rel_telocnt <- adata$telocnt / adata$nFrags_atac
adata$log_rel_telocnt <- log10((1+adata$telocnt) / adata$nFrags_atac)

adata$rank_rel_telocnt <- NA
for(cursample in unique(adata$sample)){
  adata$rank_rel_telocnt[adata$sample==cursample] <- rank(adata$rel_telocnt[adata$sample==cursample])/sum(adata$sample==cursample) #try again
}





########################################################################
################# Figures for paper ####################################
########################################################################


unique(adata$stage)

adata <-  adata[,!is.na(adata$stage)]
#adata <- adata[,adata$sample!="E7.75_rep1"] #oddly few telomeres

adata <- adata[,!(adata$sample=="E8.5_CRISPR_T_KO" | adata$sample=="E8.5_CRISPR_T_WT")]

unique(adata$sample)

p1 <- DimPlot(object = adata, reduction = 'UMAP_RNA', group.by = "stage", label = TRUE) + ggtitle("Emb.Stage")
p2 <- DimPlot(object = adata, reduction = 'UMAP_RNA', group.by = "Phase", label = TRUE) + ggtitle("CC Phase")

p3 <- FeaturePlot(object = adata, reduction = 'UMAP_RNA', features = c("rank_rel_telocnt")) + 
  scale_colour_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 3, name = "RdBu"))) + 
  ggtitle("rank nTA")

p4 <- FeaturePlot(object = adata, reduction = 'UMAP_RNA', features = c("nCount_RNA")) + ggtitle("nCount RNA")
p5 <- FeaturePlot(object = adata, reduction = 'UMAP_RNA', features = c("nFrags_atac")) + ggtitle("nFrags ATAC")
p6 <- FeaturePlot(object = adata, reduction = 'UMAP_RNA', features = c("NucleosomeRatio_atac")) + 
  scale_colour_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 3, name = "RdBu"))) +
  ggtitle("N.S.")
p7 <- FeaturePlot(object = adata, reduction = 'UMAP_RNA', features = c("mitochondrial_percent_RNA")) + ggtitle("RNA MT%")

ptot <- egg::ggarrange(p1,p2,p3,p4,p5,p6,p7)
ptot
#ggsave(plot=ptot, "arg_all_qc.svg", width = 10, height = 10)
ggsave(plot=ptot, "arg_all_qc.png", width = 10, height = 10)


DimPlot(object = adata, reduction = 'UMAP_RNA', group.by = "celltype.mapped", label = TRUE)
ggsave("arg_all_ct.svg", width = 20, height = 10)
ggsave("arg_all_ct.png", width = 20, height = 10)


