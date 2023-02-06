if(FALSE){

  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
  
  devtools::install_github("stuart-lab/signac", ref = "develop")
  
  BiocManager::install("hdf5r", force = T)
  BiocManager::install("biovizBase", force = T)
  BiocManager::install("EnsDb.Hsapiens.v86", force = T)
  BiocManager::install("JASPAR2020", force = T)
  BiocManager::install("BSgenome.Hsapiens.UCSC.hg38", force = T)
  BiocManager::install("celldex", force = T)
  BiocManager::install('SingleR', force = T)
  BiocManager::install('SingleCellExperiment', force=T)
  BiocManager::install("chromVAR")
  
  devtools::install_github("GreenleafLab/chromVARmotifs")
  
  devtools::install_github("cellgeni/sceasy")  
  
  devtools::install_github('quadbiolab/Pando')
  
  install.packages('grr')
  install.packages('JASPAR2020')
  install.packages('rhdf5')
  install.packages('plotly')
  #install.packages('enrichR')
  #install.packages("scran")
  install.packages('rsvd')
  install.packages('SeuratWrappers')
  install.packages('Seurat')
  install.packages('celldex')
  install.packages('SingleR')
  install.packages("Signac")
  install.packages("umap")
  install.packages("strawr")
    
  #library(scran)
  library(enrichR)

  install.packages('pander')
  library(pander)
  library(strawr)
}

set.seed(1234)


library(Signac)
library(Seurat)
library(chromVAR) 
library(motifmatchr)
library(Matrix)
library(SummarizedExperiment)
library(BiocParallel)
library(hdf5r)
library(BSgenome.Hsapiens.UCSC.hg38)
library(EnsDb.Hsapiens.v86)
library(JASPAR2020)
library(TFBSTools)
library(plotly)
library(dplyr)
library(stringr)
library(ggplot2)
library(patchwork)
library(SingleR)
library(celldex)
library(umap)
library(devtools)
library(reticulate)
library(cowplot)
library(sceasy)
library(multiHiCcompare)
library(sqldf)
library(reshape2)
library(motifmatchr)


`%!in%` <- Negate(`%in%`)

#https://github.com/GreenleafLab/chromVARmotifs
library(chromVARmotifs)


################################################################################
################### Load data and metadata #####################################
################################################################################



# load the RNA and ATAC data
counts <- Read10X_h5("/corgi/cellbuster/10xpbmc/aggr/outs/filtered_feature_bc_matrix.h5")
fragpath <- "/corgi/cellbuster/10xpbmc/aggr/outs/atac_fragments.tsv.gz"

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"

# create a Seurat object containing the RNA adata
adata <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# create ATAC assay and add it to the object
adata[["peaks"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)

adata$lib <- str_split_fixed(colnames(adata),"-",2)[,2]

remove(counts)


## telomemore information --- for aggregated data. list below is from libs.csv ... /corgi/cellbuster/10xpbmc/libs.csv 
map_ds2i <- data.frame(row.names=c(
  "10k_PBMC_Multiome_nextgem_Chromium_Controller",
  "10k_PBMC_Multiome_nextgem_Chromium_X",
  "pbmc_granulocyte_sorted_10k",  #got more cells than 1,2
  "pbmc_granulocyte_unsorted_10k"),num=c(1,2,3,4))
teloinfo <- read.csv("/corgi/cellbuster/10xpbmc/summary_kmer.java.csv")
rownames(teloinfo) <- sprintf("%s-%s",str_sub(teloinfo$barcode,1,16),map_ds2i[teloinfo$dataset,,drop=FALSE]$num)
sum(rownames(teloinfo) %in% names(adata$orig.ident))
teloinfo <- teloinfo[names(adata$orig.ident),]
adata$rawcnt_telo <- teloinfo$totalcnt_CCCTAA
adata$dedupcnt_telo <- teloinfo$dedupcnt_CCCTAA
adata$totalcnt_telo <- teloinfo$total
adata$norm_telo <- log10(adata$dedupcnt_telo/(adata$totalcnt_telo+1)+1e-3)
hist(adata$norm_telo)

#should compute rank within each batch!
adata$rank_norm_telo<-NA
for(i in c(1:4)){
  #adata$rank_norm_telo[adata$lib==i] <- rank(adata$norm_telo[adata$lib==i])
  allrank <- rank(adata$norm_telo[adata$lib==i])
  adata$rank_norm_telo[adata$lib==i] <- allrank/max(allrank)
}


################################################################################
################### RNA processing with batch correction #######################
################################################################################


DefaultAssay(adata) <- "RNA"

adata[["percent.mt"]] <- PercentageFeatureSet(adata, pattern = "^MT")
adata[['percent.ribo']]<- PercentageFeatureSet(adata, "^RP[SL]")

## Visualize QC metrics as a violin plot
VlnPlot(adata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 2,pt.size = 0)
plot1 <- FeatureScatter(adata, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(adata, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

## Remove bad cells
adata<- subset(adata, 
               subset = percent.mt < 40 & 
                 nFeature_RNA > 500 & nFeature_RNA < 10000 & 
                 nCount_RNA > 200 & nCount_RNA<30000)

###############

adata.list <- SplitObject(adata, split.by = "lib")

adata.list <- lapply(X = adata.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = adata.list)
immune.anchors <- FindIntegrationAnchors(object.list = adata.list, anchor.features = features)
remove(adata.list)
adata <- IntegrateData(anchorset = immune.anchors, new.assay.name = "integrated")
remove(immune.anchors)


# Run the standard workflow for visualization and clustering
DefaultAssay(adata) <- "integrated"
adata <- ScaleData(adata, verbose = FALSE)
adata <- RunPCA(adata, npcs = 30, verbose = FALSE)
adata <- RunUMAP(adata, reduction = "pca", dims = 1:30)


#cell cycle scoring
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
adata <- CellCycleScoring(adata, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
if(FALSE){
  #Possible CC removal
  adata <- ScaleData(adata, vars.to.regress = c('S.Score','G2M.Score'))
}




################################################################################
################### ATAC basic processing ######################################
################################################################################

DefaultAssay(adata) <- "peaks"
adata <- NucleosomeSignal(adata)
adata <- TSSEnrichment(adata, fast=FALSE) 

## QC plots
VlnPlot(
  object = adata,
  features = c("nCount_RNA", "nCount_peaks", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)

adata$high.tss <- ifelse(adata$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(adata, group.by = 'high.tss') + NoLegend()

adata$nucleosome_group <- ifelse(adata$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = adata, group.by = 'nucleosome_group')

## Filter out low quality cells
adata <- subset(
  x = adata,
  subset = nCount_peaks < 100e03 &  #up to 200k in some. below 50k likely fine
    nCount_peaks > 10e03 & #below this value, nTA cannot be too trusted 
    nucleosome_signal < 2 &
    TSS.enrichment > 1 #very loose. likely fine
)




################################################################################
################### Load motif data ############################################
################################################################################


DefaultAssay(adata) <- "peaks"

if(TRUE){
  
  ## Subset to only use standard chromosomes
  main.chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
  keep.peaks <- which(as.character(seqnames(granges(adata))) %in% main.chroms)
  adata[["peaks"]] <- subset(adata[["peaks"]], features = rownames(adata[["peaks"]])[keep.peaks])
  
  ## Get and add motifs from JASPAR
  pfm <- getMatrixSet(
    x = JASPAR2020,
    opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
  )
  adata <- AddMotifs(adata, genome = BSgenome.Hsapiens.UCSC.hg38,pfm = pfm,verbose = T)
} else {
  pfm <- human_pwms_v1
  adata <- AddMotifs(adata, genome = BSgenome.Hsapiens.UCSC.hg38,pfm = pfm,verbose = T)
}


## Load motif data
motif2names <- data.frame(
  row.names=names(pfm),
  id=names(pfm),
  name=sapply(1:length(pfm), function(x) name(pfm[[x]]))
)
names2motifs <- motif2names
rownames(names2motifs) <- names2motifs$name



################################################################################
################### automatic cell annotation by SingleR #######################
################################################################################

#HPCA has a lot of B_cell: subtypes
#Monaco Has "Non-switched memory B cells", "Naive B cells","Exhausted B cells","Plasmablasts"
#DICE is different T cells

DefaultAssay(adata) <- "integrated"

dice <-   celldex::DatabaseImmuneCellExpressionData()
pred.dice <- SingleR(test = as.SingleCellExperiment(adata), ref = dice, assay.type.test=1,
                     labels = dice$label.fine)
adata$pred.dice <- pred.dice$labels
#10k_PBMC_Multiome_nextgem_Chromium_X  : CD16+ mono has hi telo. and NK. and some T cells toward NK


hpca <-   celldex::HumanPrimaryCellAtlasData()
pred.HPCA <- SingleR(test = as.SingleCellExperiment(adata), ref = hpca, assay.type.test=1,
                     labels = hpca$label.fine)
adata$pred.hpca <- pred.HPCA$labels

monaco <- celldex::MonacoImmuneData()
pred.monaco <- SingleR(test = as.SingleCellExperiment(adata), ref = monaco, assay.type.test=1,
                       labels = monaco$label.fine)
adata$pred.monaco <- pred.monaco$labels

## Filtering by HPCA
#adata <- adata[,str_starts(adata$pred.hpca,"B_cell")]

## Filtering by Monaco. cluster is detected as doublets to great extent as well; likely technical reasons
#adata <- adata[,adata$pred.monaco %!in% c("MAIT cells", "Follicular helper T cells", "Non-Vd2 gd T cells","Plasmacytoid dendritic cells")]
#unique(adata$pred.monaco)



########################################################################
################# RNAseq dim red & clustering ##########################
########################################################################


DefaultAssay(adata) <- "integrated"

## Dimensional reduction
adata <- RunPCA(adata,reduction.name = 'PCA_RNA',ndims.print = c(1,3,5), verbose = T)
adata <- RunUMAP(object = adata, dims = 1:30, reduction = "PCA_RNA", reduction.name="UMAP_RNA")
adata <- RunUMAP(adata, dims = 1:30,reduction = 'PCA_RNA',n.components = 3L, reduction.name ='UMAP_RNA_3D' )
adata <- FindNeighbors(adata, dims = 1:30,reduction = 'PCA_RNA')


## Split everything 
adata <- FindClusters(adata, resolution = 0.5)
adata$rna_clusters <- factor(sprintf("r%s",adata$seurat_clusters))

## Split everything 
adata <- FindClusters(adata, resolution = 0.05)
adata$rna_clusters2 <- factor(sprintf("r2_%s",adata$seurat_clusters))

## Split everything 
adata <- FindClusters(adata, resolution = 1)
adata$rna_clusters3 <- factor(sprintf("r3_%s",adata$seurat_clusters))

## Split everything 
adata <- FindClusters(adata, resolution = 0.1)
adata$rna_clusters4 <- factor(sprintf("r4_%s",adata$seurat_clusters))


DimPlot(object = adata, reduction = 'UMAP_RNA', group.by = "rna_clusters3", label = TRUE)


################################################################################
################# ATACseq dim red & clustering #################################
################################################################################


DefaultAssay(adata) <- "peaks"

## Dimensional reduction
adata <- RunTFIDF(adata)
adata <- FindTopFeatures(adata, min.cutoff = 3, verbose = T )
adata <- RunSVD(adata, verbose = T, reduction.name = 'SVD_ATAC')
DepthCor(adata, reduction = 'SVD_ATAC',n = 20)  

#Removing first component due to strong correlation of depth.
#PC2-3 still seems to hold biology as well
use_dim <- 2:30 

adata <- RunUMAP(object = adata, reduction = 'SVD_ATAC', dims = use_dim, reduction.name = 'UMAP_ATAC')
adata <- RunUMAP(object = adata, reduction = 'SVD_ATAC', dims = use_dim, n.components = 3L, reduction.name = 'UMAP_ATAC_3D')

## Split everything
adata <- FindNeighbors(object = adata, reduction = 'SVD_ATAC', dims = use_dim)

adata <- FindClusters(object = adata, verbose = TRUE, algorithm = 3)
adata$atac_clusters <- factor(sprintf("a%s",adata$seurat_clusters))

adata <- FindClusters(object = adata, verbose = TRUE, algorithm = 3, resolution = 1.3)
adata$atac_clusters3 <- factor(sprintf("a3_%s",adata$seurat_clusters))



################################################################################
################### DE analysis ################################################
################################################################################


################### Compare ATAC-based clusters - 1
DefaultAssay(adata) <- "RNA" #yes!
Idents(adata) <- adata$atac_clusters
cluster_markers_atac1 <- FindAllMarkers(adata, only.pos = F, test.use = 'wilcox', )
cluster_markers_atac1 %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) -> top10_atac

################### Compare RNA-based clusters - 3
DefaultAssay(adata) <- "RNA"
Idents(adata) <- adata$rna_clusters3
cluster_markers_rna3 <- FindAllMarkers(adata, only.pos = F, test.use = 'wilcox')
cluster_markers_rna3 %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) -> top10_rna3

################### Compare ATAC-based clusters - a3
DefaultAssay(adata) <- "RNA" #yes!
Idents(adata) <- adata$atac_clusters3
cluster_markers_atac3 <- FindAllMarkers(adata, only.pos = F, test.use = 'wilcox', )
cluster_markers_atac3 %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) -> top10_atac3



################### What happens in the T cell tip, with nTA progression?
DefaultAssay(adata) <- "RNA"
Idents(adata) <- adata$rna_clusters3
cluster_markers_rna3_t <- FindAllMarkers(adata[,adata$rna_clusters3 %in% c("r3_8","r3_4","r3_0")], only.pos = F, test.use = 'wilcox')
cluster_markers_rna3_t %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) -> top10_rna3_t


################################################################################
################### Plots of metadata ##########################################
################################################################################

#nucleosome signal as a mitosis signal?

#nucleosome_signal = inv(rank_norm_telo), for monocytes, not for T cells

# plot(adata$nucleosome_signal_ctcf,adata$nucleosome_signal,pch=19,cex=0.1,xlim=c(0,1),ylim=c(0,1))
# adata$rank_nucleosome_signal_ctcf <- rank(adata$nucleosome_signal_ctcf)

adata$rank_nucleosome_signal <- rank(adata$nucleosome_signal)
plot(adata$rank_norm_telo,adata$nucleosome_signal,pch=19, cex=0.1)
plot(adata$norm_telo,adata$nucleosome_signal,pch=19, cex=0.1) #complex plot. some correlation inverse correlation at places. different in Bcell


FeaturePlot(object = adata, reduction = 'UMAP_RNA', features = c("nCount_RNA","percent.mt","nFeature_RNA"))
FeaturePlot(object = adata, reduction = 'UMAP_RNA', features = c("nCount_peaks", "TSS.enrichment", "nucleosome_signal"))
FeaturePlot(object = adata, reduction = 'UMAP_RNA', features = c("norm_telo"))
FeaturePlot(object = adata, reduction = 'UMAP_RNA', features = c("rank_norm_telo"))
FeaturePlot(object = adata, reduction = 'UMAP_RNA', features = c("rank_nucleosome_signal"))
#FeaturePlot(object = adata, reduction = 'UMAP_RNA', features = c("rank_nucleosome_signal_ctcf"))
#FeaturePlot(object = adata, reduction = 'UMAP_RNA', features = c("nucleosome_signal_ratio"))
#FeaturePlot(object = adata, reduction = 'UMAP_RNA', features = c("rank_norm_telo","nucleosome_signal"))
#FeaturePlot(object = adata, reduction = 'UMAP_RNA', features = c("rank_norm_telo"), split.by = "lib")




plot(adata$rank_norm_telo, adata$nucleosome_signal)
plot(adata$norm_telo, adata$nucleosome_signal) 
#for pbmc: linear negative correlation. possibly 2 groups
#for b cells: positive correlation


#FeaturePlot(object = adata, reduction = 'UMAP_RNA', features = c("prob_doublet"))
DimPlot(object = adata, reduction = 'UMAP_RNA', group.by = "Phase", label = TRUE)
DimPlot(object = adata, reduction = 'UMAP_RNA', group.by = "lib")
DimPlot(object = adata, reduction = 'UMAP_RNA', pt.size = 0.3, group.by = "trust4", cols=c("black","green","white","red"))

DimPlot(object = adata, reduction = 'UMAP_RNA', group.by = "rna_clusters", label = TRUE)
DimPlot(object = adata, reduction = 'UMAP_RNA', group.by = "rna_clusters2", label = TRUE)
DimPlot(object = adata, reduction = 'UMAP_RNA', group.by = "rna_clusters3", label = TRUE)
DimPlot(object = adata, reduction = 'UMAP_RNA', group.by = "rna_clusters4", label = TRUE)
#DimPlot(object = adata, reduction = 'UMAP_RNA', group.by = "ct", label = TRUE)

######### Predicted cell types
DimPlot(object = adata, reduction = 'UMAP_RNA', group.by = "pred.dice", label = TRUE)
DimPlot(object = adata, reduction = 'PCA_RNA', group.by = "pred.dice", label = TRUE)  
#new CD14 monocyte cluster in lib3,4

DimPlot(object = adata, reduction = 'UMAP_RNA', group.by = "pred.hpca", label = TRUE)
DimPlot(object = adata, reduction = 'PCA_RNA', group.by = "pred.hpca", label = TRUE)  


DimPlot(object = adata, reduction = 'UMAP_RNA', group.by = "pred.monaco", label = TRUE)  #cleaner
DimPlot(object = adata, reduction = 'PCA_RNA', group.by = "pred.monaco", label = TRUE)  


#Down in G1 phase here; confounded by cell type!
Idents(adata) <- adata$Phase
VlnPlot(adata, c("rank_norm_telo"))

#Highly correlated
plot(adata$rawcnt_telo, adata$nCount_ATAC)
plot(adata$dedupcnt_telo, adata$nCount_ATAC)
#less correlated! highest rank have the least #ATAC. below #ATAC, 1e04, be careful!
plot(adata$rank_norm_telo, adata$nCount_ATAC)


if(FALSE){
  ########################## For paper!
  #NK cells up
  #Monocyte CD16- down
  df <- data.frame(
    id=adata$pred.hpca,
    telo=adata$rank_norm_telo)
  df <- df[df$id %in% names(table(df$id))[table(df$id)>30],]
  df$id <- factor(df$id, levels=sqldf("select id, avg(telo) as m from df group by id order by m desc")$id)
  ggplot(df, aes(id,telo))+ geom_violin(scale="width", fill="gray") + geom_point(size=0.5) + coord_flip()
  
  
  #NK cells up
  #naive cells down; mature ones up
  df <- data.frame(
    id=adata$pred.dice,
    telo=adata$rank_norm_telo)
  df <- df[df$id %in% names(table(df$id))[table(df$id)>30],]
  df$id <- factor(df$id, levels=sqldf("select id, avg(telo) as m from df group by id order by m desc")$id)
  ggplot(df, aes(id,telo))+ geom_violin(scale="width", fill="gray") + geom_point(size=0.5) + coord_flip()
  
  #mature up
  #monocytes down
  df <- data.frame(
    id=adata$pred.monaco,
    telo=adata$rank_norm_telo)
  df <- df[df$id %in% names(table(df$id))[table(df$id)>30],]
  df$id <- factor(df$id, levels=sqldf("select id, avg(telo) as m from df group by id order by m desc")$id)
  ggplot(df, aes(id,telo))+ geom_violin(scale="width", fill="gray") + geom_point(size=0.5) + coord_flip()
  
  
  #also holds on the absolute scale
  df <- data.frame(
    id=adata$pred.monaco,
    telo=adata$rawcnt_telo)
  df <- df[df$id %in% names(table(df$id))[table(df$id)>30],]
  df$id <- factor(df$id, levels=sqldf("select id, avg(telo) as m from df group by id order by m desc")$id)
  ggplot(df, aes(id,telo))+ geom_violin(scale="width", fill="gray") + geom_point(size=0.5) + coord_flip()
  
}



######
FeaturePlot(object = adata, reduction = 'UMAP_ATAC', features = c("nCount_RNA","percent.mt","nFeature_RNA"))
FeaturePlot(object = adata, reduction = 'UMAP_ATAC', features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"))
FeaturePlot(object = adata, reduction = 'UMAP_ATAC', features = c("norm_telo"))
FeaturePlot(object = adata, reduction = 'UMAP_ATAC', features = c("rank_norm_telo"))
FeaturePlot(object = adata, reduction = 'UMAP_ATAC', features = c("prob_doublet"))
DimPlot(object = adata, reduction = 'UMAP_ATAC', group.by = "Phase", label = TRUE, cols=c("black","green","white","red"))
DimPlot(object = adata, reduction = 'UMAP_ATAC', group.by = "donor_id", label = TRUE)
DimPlot(object = adata, reduction = 'UMAP_ATAC', group.by = "lib")
DimPlot(object = adata, reduction = 'UMAP_ATAC', pt.size = 0.3, group.by = "trust4", cols=c("black","green","white","red"))

DimPlot(object = adata, reduction = 'UMAP_ATAC', group.by = "atac_clusters", label = TRUE)
DimPlot(object = adata, reduction = 'UMAP_ATAC', group.by = "ct", label = TRUE)
#atac#6 or 9 are gdT. mainly 6. what is the extra 9?


################################################################################
################### The PCA view ###############################################
################################################################################

##### RNA

DimPlot(object = adata, reduction = 'PCA_RNA', group.by = "pred.hpca", label = TRUE)  
DimPlot(object = adata, reduction = 'PCA_RNA', group.by = "pred.monaco", label = TRUE)  
DimPlot(object = adata, reduction = 'PCA_RNA', group.by = "Phase", label = TRUE, cols=c("black","green","red"))

DimPlot(object = adata, reduction = 'PCA_RNA', group.by = "rna_clusters", label = TRUE)  


FeaturePlot(object = adata, reduction = 'PCA_RNA', features=c("CXCR4","CD83"), label = TRUE)  
FeaturePlot(object = adata, reduction = 'UMAP_RNA', features=c("FCRL4"), label = TRUE)  
FeaturePlot(object = adata, reduction = 'UMAP_RNA', features=c("FCRL3"), label = TRUE)  #in new r2 if anything




#### ATAC

use_atac_dims <- c(2,4)
DimPlot(object = adata, reduction = 'SVD_ATAC', group.by = "pred.hpca", label = TRUE, dims = use_atac_dims)  
DimPlot(object = adata, reduction = 'SVD_ATAC', group.by = "pred.monaco", label = TRUE, dims = use_atac_dims)  
DimPlot(object = adata, reduction = 'SVD_ATAC', group.by = "Phase", label = TRUE, cols=c("black","green","red"), dims = use_atac_dims)

DimPlot(object = adata, reduction = 'SVD_ATAC', group.by = "rna_clusters", label = TRUE, dims = use_atac_dims)

FeaturePlot(object = adata, reduction = 'SVD_ATAC', dims = use_atac_dims, features=c("CXCR4","CD83"), label = TRUE)  



################################################################################
################### Comparison of ATAC and RNA clusters ########################  UMAP
################################################################################

DimPlot(object = adata, reduction = "UMAP_ATAC", pt.size = 1,group.by = 'Phase')+ ggtitle('scATAC seq')


p1 <- DimPlot(object = adata, reduction = 'UMAP_RNA',  pt.size = 1,group.by = "rna_clusters", label = TRUE)+ ggtitle('rna clusters on scRNA umap')
p2 <- DimPlot(object = adata, reduction = "UMAP_ATAC", pt.size = 1,group.by = 'atac_clusters', label = TRUE)+ ggtitle('atac clusters on scATAC umap')
p3 <- DimPlot(object = adata, reduction = "UMAP_ATAC", pt.size = 1,group.by = 'rna_clusters', label = TRUE)  + ggtitle('rna clusters on ATAC umap')
p4 <- DimPlot(object = adata, reduction = 'UMAP_RNA',  pt.size = 1,group.by = "atac_clusters", label = TRUE) + ggtitle("atac clusters on rna umap")

p5 <- DimPlot(object = adata, reduction = 'UMAP_RNA',  pt.size = 1,group.by = "Phase", label = TRUE)+ ggtitle("phase on rna umap")
p6 <- DimPlot(object = adata, reduction = 'UMAP_ATAC', pt.size = 1,group.by = "Phase", label = TRUE)+ ggtitle("phase on atac umap")

p1/p2|p3/p4
p1/p2|p3/p4|p5/p6

############################

p1 <- DimPlot(object = adata, reduction = 'UMAP_RNA',  pt.size = 1,group.by = "rna_clusters3", label = TRUE)+ ggtitle('rna clusters on scRNA umap')
p2 <- DimPlot(object = adata, reduction = "UMAP_ATAC", pt.size = 1,group.by = 'atac_clusters3', label = TRUE)+ ggtitle('atac clusters on scATAC umap')
p3 <- DimPlot(object = adata, reduction = "UMAP_ATAC", pt.size = 1,group.by = 'rna_clusters3', label = TRUE)  + ggtitle('rna clusters on ATAC umap')
p4 <- DimPlot(object = adata, reduction = 'UMAP_RNA',  pt.size = 1,group.by = "atac_clusters3", label = TRUE) + ggtitle("atac clusters on rna umap")

p5 <- DimPlot(object = adata, reduction = 'UMAP_RNA',  pt.size = 1,group.by = "Phase", label = TRUE)+ ggtitle("phase on rna umap")
p6 <- DimPlot(object = adata, reduction = 'UMAP_ATAC', pt.size = 1,group.by = "Phase", label = TRUE)+ ggtitle("phase on atac umap")

p1/p2|p3/p4
p1/p2|p3/p4|p5/p6






################################################################################
################### Plot genes 2D ##############################################
################################################################################


DefaultAssay(adata) <- "RNA"

showgenes <- c(
  "SOX5"   
)
FeaturePlot(adata, features = showgenes, reduction = "UMAP_RNA") 
FeaturePlot(adata, features = showgenes, reduction = "UMAP_ATAC") 
FeaturePlot(adata, features = showgenes, reduction = "PCA_RNA") 
FeaturePlot(adata, features = showgenes, reduction = "SVD_ATAC", dims = use_atac_dims) 





########################################################################
################ 3d umap ###############################################
########################################################################

plot_umap_groupby <- function(adata, usered, colorby, textby=colorby, 
                              colors=c(
                                "darkgreen",
                                "red",
                                "black",
                                "yellow4",
                                "royalblue1",
                                "darkmagenta")){
  axname <- names(adata@reductions[[usered]])
  plot.data <- FetchData(object = adata, vars = c(axname,colorby,textby))
  colnames(plot.data) <- c("x","y","z","colorby","textby")
  if(textby!=colorby){
    plot.data$textby <- sprintf("%s %s",plot.data$colorby,plot.data$textby)
  }
  fig_atac <- plot_ly(data = plot.data, 
                      x = ~x, y = ~y, z = ~z, 
                      color = ~colorby, 
                      colors = colors,
                      type = "scatter3d", 
                      mode = "markers", 
                      marker = list(size = 1.5, width=2), 
                      text=~textby, 
                      hoverinfo="text") 
  fig_atac
}

plot_umap_gene <- function(adata, usered, colorby, textby=colorby, dorank=FALSE){
  axname <- names(adata@reductions[[usered]])
  plot.data <- FetchData(object = adata, vars = c(axname,colorby,textby))
  colnames(plot.data) <- c("x","y","z","colorby","textby")
  if(dorank){
    plot.data$colorby <- rank(plot.data$colorby)
  }
  fig_atac <- plot_ly(data = plot.data, 
                      x = ~x, y = ~y, z = ~z, 
                      color = ~colorby, 
                      type = "scatter3d", 
                      mode = "markers", 
                      marker = list(size = 1.5, width=2), 
                      text=~textby, 
                      hoverinfo="text") 
  fig_atac
}


plot_umap_groupby(adata,"UMAP_RNA_3D","Phase", colors=c("red","green","blue"))
plot_umap_groupby(adata,"UMAP_ATAC_3D","Phase", colors=c("red","green","blue"))

plot_umap_groupby(adata,"UMAP_RNA_3D","rna_clusters","Phase")
plot_umap_groupby(adata,"UMAP_RNA_3D","rna_clusters2","Phase")
plot_umap_groupby(adata,"UMAP_RNA_3D","rna_clusters3","Phase") #Good level to follow transitions!  
plot_umap_groupby(adata,"UMAP_RNA_3D","rna_clusters4","Phase")
plot_umap_groupby(adata,"UMAP_RNA_3D","ct","Phase")
plot_umap_groupby(adata,"UMAP_RNA_3D","Phase","ct")

plot_umap_groupby(adata,"UMAP_ATAC_3D","rna_clusters","Phase")
plot_umap_groupby(adata,"UMAP_ATAC_3D","rna_clusters3","Phase")
plot_umap_groupby(adata,"UMAP_ATAC_3D","rna_clusters4","Phase")
plot_umap_groupby(adata,"UMAP_ATAC_3D","atac_clusters3","ct")   #high nTA T cells are r3_8 and r3_4; r3_11 is random
plot_umap_groupby(adata,"UMAP_ATAC_3D","ct","atac_clusters3")



plot_umap_groupby(adata,"UMAP_RNA_3D","pred.dice","Phase")  
plot_umap_groupby(adata,"UMAP_RNA_3D","pred.hpca","Phase")  
plot_umap_groupby(adata,"UMAP_RNA_3D","pred.monaco","Phase")  
plot_umap_groupby(adata,"UMAP_RNA_3D","Phase","pred.hpca")


plot_umap_gene(adata,"UMAP_RNA_3D","rank_norm_telo","pred.hpca")
plot_umap_gene(adata,"UMAP_RNA_3D","totalcnt_telo","pred.hpca")  #not much affected!
plot_umap_gene(adata,"UMAP_RNA_3D","nCount_RNA","pred.hpca")  #up in CD16+ monocytes
plot_umap_gene(adata,"UMAP_RNA_3D","nCount_ATAC","pred.hpca")  #varies; complex pattern
plot_umap_groupby(adata,"UMAP_RNA_3D","pred.hpca")
plot_umap_groupby(adata,"UMAP_RNA_3D","pred.monaco")


plot_umap_gene(adata,"UMAP_ATAC_3D","rank_norm_telo","pred.hpca")
plot_umap_gene(adata,"UMAP_ATAC_3D","totalcnt_telo","pred.hpca")  #not much affected!
plot_umap_gene(adata,"UMAP_ATAC_3D","nCount_RNA","pred.hpca")  #up in CD16+ monocytes
plot_umap_gene(adata,"UMAP_ATAC_3D","nCount_ATAC","pred.hpca")  #varies; complex pattern
plot_umap_groupby(adata,"UMAP_ATAC_3D","pred.dice")
plot_umap_groupby(adata,"UMAP_ATAC_3D","pred.hpca")
plot_umap_groupby(adata,"UMAP_ATAC_3D","pred.monaco")

adata$nCount_ATAC



plot_umap_gene(adata,"UMAP_ATAC_3D","CD8A","rna_clusters3")    #distinct from the rest
plot_umap_gene(adata,"UMAP_ATAC_3D","GZMK","rna_clusters3")    #QUITE UP
plot_umap_gene(adata,"UMAP_ATAC_3D","GZMH","rna_clusters3")    #MOST UP
plot_umap_gene(adata,"UMAP_ATAC_3D","GZMA","rna_clusters3")    #well correlating all over 

plot_umap_gene(adata,"UMAP_ATAC_3D","BACH2","rna_clusters3")   #really well anticorrelating; for paper
plot_umap_gene(adata,"UMAP_ATAC_3D","LEF1","rna_clusters3")    #really well anticorrelating; for paper

plot_umap_gene(adata,"UMAP_ATAC_3D","CCL5","rna_clusters3")
plot_umap_gene(adata,"UMAP_ATAC_3D","NKG7","rna_clusters3")
plot_umap_gene(adata,"UMAP_ATAC_3D","GZMH","rna_clusters3")

plot_umap_gene(adata,"UMAP_ATAC_3D","ZEB2","rna_clusters3")

############ Make sense of clusters
plot_umap_gene(adata,"UMAP_RNA_3D","GZMK","rna_clusters3")  
plot_umap_gene(adata,"UMAP_RNA_3D","CD8A","rna_clusters3")  
plot_umap_gene(adata,"UMAP_RNA_3D","CCL5","rna_clusters3")  
plot_umap_gene(adata,"UMAP_RNA_3D","NKG7","rna_clusters3")  
plot_umap_gene(adata,"UMAP_RNA_3D","CD4","rna_clusters3")   
plot_umap_gene(adata,"UMAP_RNA_3D","GZMH","rna_clusters3")  


plot_umap_gene(adata,"UMAP_ATAC_3D","CD14","rna_clusters3")
plot_umap_gene(adata,"UMAP_ATAC_3D","FCGR3A","rna_clusters3") #CD16


plot_umap_gene(adata,"UMAP_RNA_3D","CDKN1C","rna_clusters3") 
plot_umap_gene(adata,"UMAP_RNA_3D","CD14","rna_clusters3")
plot_umap_gene(adata,"UMAP_RNA_3D","FCGR3A","rna_clusters3") #CD16
plot_umap_gene(adata,"UMAP_RNA_3D","FUT4","rna_clusters3") #CD15
plot_umap_gene(adata,"UMAP_RNA_3D","ITGAM","rna_clusters3") #CD11b


 #CD11b

########################################################################
################ comparison to other annotations, DICE #################
########################################################################


plot_umap_groupby(adata,"UMAP_RNA_3D","pred.hpca")
plot_umap_groupby(adata,"UMAP_ATAC_3D","pred.hpca")

plot_umap_groupby(adata,"UMAP_RNA_3D","pred.monaco","rna_clusters")
plot_umap_groupby(adata,"UMAP_ATAC_3D","pred.monaco")

count_overlaps_clusterings <- function(the_rna_clusters,the_atac_clusters){
  list_rna_clusters <- unique(the_rna_clusters)
  list_atac_clusters <- unique(the_atac_clusters)
  confusion_matrix <- matrix(nrow=length(list_atac_clusters), ncol=length(list_rna_clusters))
  for(i in 1:length(list_atac_clusters)){
    for(j in 1:length(list_rna_clusters)){
      confusion_matrix[i,j] <- sum(the_rna_clusters==list_rna_clusters[j] & the_atac_clusters==list_atac_clusters[i])
    }
  }
  colnames(confusion_matrix) <- list_rna_clusters
  rownames(confusion_matrix) <- list_atac_clusters
  # confusion_matrix <- confusion_matrix[order(rowSums(confusion_matrix)),]
  # confusion_matrix <- confusion_matrix[,order(colSums(confusion_matrix))]
  confusion_matrix <- confusion_matrix[order(apply(confusion_matrix,1,max)),]
  confusion_matrix <- confusion_matrix[,order(apply(confusion_matrix,2,max))]
  confusion_matrix
}

count_overlaps_clusterings(adata$rna_clusters, adata$pred.monaco)
count_overlaps_clusterings(adata$rna_clusters, adata$pred.hpca)

########################################################################
################ Celltype assignment ###################################
########################################################################

assignCelltypes <- function(the_clusters, ct_assignment){
  ct_assignment <- str_split_fixed(ct_assignment,"\\:",2)
  ct_assignment <- data.frame(row.names = ct_assignment[,1], newname=factor(ct_assignment[,2]))
  missing_ct <- unique(the_clusters[!(the_clusters %in% rownames(ct_assignment))])
  if(length(missing_ct)>0){
    print("Missing ct")
    print(missing_ct)
  }
  ct_assignment[as.character(the_clusters),,drop=FALSE]$newname  
}

# #todo
# ct_assignment <- c(
#   "r4_12:Plasmablast",  #no doubt high in xbp1
#   "r4_13:nonGC naive odd"  #up in IGHG1 oddly. somewhat toward plasma
# )
# adata$ct <- assignCelltypes(adata$rna_clusters4, ct_assignment)
#
# Idents(adata) <- adata$ct
# DotPlot(adata, features = unique(features), assay="RNA") + RotatedAxis()



############################
#Might not make sense with a dotplot for motif activity, since it is a dense matrix!
#Here, heatmap
listct <- levels(adata$ct)
meanact <- matrix(nrow=length(listct), ncol=nrow(adata@assays$chromvar@data))
for(i in 1:length(listct)){
  meanact[i,] <- rowMeans(adata@assays$chromvar@data[,adata$ct==listct[i]])
}
rownames(meanact) <- listct
colnames(meanact) <- rownames(adata@assays$chromvar@data)

#meanact <- meanact[rownames(meanact)!="Plasmablast",]
for(i in 1:ncol(meanact)){
  meanact[,i] <- meanact[,i]-min(meanact[,i])
  meanact[,i] <- meanact[,i]/max(meanact[,i])
}

m <- merge(melt(meanact),data.frame(Var2=motif2names$id, tf=motif2names))
m <- m[m$tf.name %in% features,]
m$tf.name <- factor(m$tf.name, levels = unique(features))

ggplot(m, aes(tf.name, Var1, fill=value)) +   #    geom_tile() + 
  geom_point(aes(size=value, color=value)) + 
  scale_size(range = c(0, 5)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


########################################################################
################ Chromvar analysis #####################################
########################################################################

library(BiocParallel)
register(MulticoreParam(2)) #avoid swap trashing

DefaultAssay(adata) <- "peaks"

## Compute motif activity; this is deviations in counts for each motif across peaks
## as compared to total count. it does not care /which/ peaks the motif is in, just the total sum
##https://www.nature.com/articles/nmeth.4401
adata<- RunChromVAR(
  object = adata,
  genome = BSgenome.Hsapiens.UCSC.hg38, assay = 'peaks'
)

### possible that this should be run on each separately TODO

########################################################################
################ Find DE motif activity ################################
########################################################################


#Based on RNA
DefaultAssay(adata) <- 'chromvar'
Idents(adata) <- adata$rna_clusters3
cluster_motifs_rna3 <- FindAllMarkers(adata, only.pos = F, test.use = 'wilcox')
cluster_motifs_rna3 %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) -> topmotifs_rna3
topmotifs_rna3$name <- motif2names[topmotifs_rna3$gene,]$name


DefaultAssay(adata) <- 'chromvar'
Idents(adata) <- adata$rna_clusters4
cluster_motifs_rna4 <- FindAllMarkers(adata, only.pos = F, test.use = 'wilcox')
cluster_motifs_rna4 %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) -> topmotifs_rna4
topmotifs_rna4$name <- motif2names[topmotifs_rna4$gene,]$name


DefaultAssay(adata) <- 'chromvar'
Idents(adata) <- adata$rna_clusters
cluster_motifs_rna <- FindAllMarkers(adata, only.pos = F, test.use = 'wilcox')
cluster_motifs_rna %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) -> topmotifs_rna
topmotifs_rna$name <- motif2names[topmotifs_rna$gene,]$name

#specifically between 2 naive groups; PDE4D vs MAST4
DefaultAssay(adata) <- 'chromvar'
Idents(adata) <- adata$rna_clusters4
cluster_motifs_twomem <- FindAllMarkers(adata[,adata$rna_clusters4 %in% c("r4_1","r4_2")], only.pos = F, test.use = 'wilcox')
cluster_motifs_twomem %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) -> topmotifs_twomem
topmotifs_twomem$name <- motif2names[topmotifs_twomem$gene,]$name


#Based on atac
DefaultAssay(adata) <- 'chromvar'
Idents(adata) <- adata$atac_clusters3
cluster_motifs_atac3 <- FindAllMarkers(adata, only.pos = F, test.use = 'wilcox')
cluster_motifs_atac3 %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) -> topmotifs_atac3
topmotifs_atac3$name <- motif2names[topmotifs_atac3$gene,]$name


DefaultAssay(adata) <- 'chromvar'
Idents(adata) <- adata$atac_clusters
cluster_motifs_atac <- FindAllMarkers(adata, only.pos = F, test.use = 'wilcox')
cluster_motifs_atac %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) -> topmotifs_atac
topmotifs_atac$name <- motif2names[topmotifs_atac$gene,]$name




#Based on phase
DefaultAssay(adata) <- 'chromvar'
Idents(adata) <- adata$Phase
cluster_motifs_phase <- FindAllMarkers(adata, only.pos = F, test.use = 'wilcox')
cluster_motifs_phase %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) -> topmotifs_phase
topmotifs_phase$name <- motif2names[topmotifs_phase$gene,]$name



########################################################################
################ Motif compared to axes ################################   JASPAR motifs
########################################################################


compare_gene_with_axis <- function(adata, the_axis){
  use_genes <- adata@assays$integrated@var.features
  allcor <- NULL
  for(i in use_genes){
    print(i)
    allcor <- c(allcor, cor(as.double(adata@assays$integrated[i,]),the_axis,method = "spearman"))
  }
  df <- data.frame(
    sumexp=as.double(rowSums(adata@assays$integrated[adata@assays$integrated@var.features,])),
    cor=allcor,
    gene=use_genes)
#    motif=rownames(adata@assays$RNA))
  df <- df[order(df$cor),]
  df <- df[!is.na(df$cor),]
  df  
}


compare_motif_with_axis <- function(adata, the_axis){
  allcor <- NULL
  for(i in 1:nrow(adata@assays$chromvar)){
    allcor <- c(allcor, cor(as.double(adata@assays$chromvar[i,]),the_axis,method = "spearman"))
  }
  df <- data.frame(
    cor=allcor,
    motif=rownames(adata@assays$chromvar))
  df$name <- motif2names[df$motif,]$name
  df <- df[order(df$cor),]
  df <- df[!is.na(df$cor),]
  df  
}




################################### nTelo axis #######################


###### nTelo axis  -- within CD4 T
axis_telo <- adata$rank_norm_telo
use_sub <- adata$rna_clusters2 %in% c("r2_0")
genecor_telo <- compare_gene_with_axis(adata[,use_sub], axis_telo[use_sub]) #BACH2 still negative corr. pos corr, no clue
rbind(head(genecor_telo,n=20),tail(genecor_telo, n=20))

motifcor_telo <- compare_motif_with_axis(adata[,use_sub], axis_telo[use_sub])
rbind(head(motifcor_telo),tail(motifcor_telo))
### new
#TCF7L2 / LEF1 / Ahr::Arnt negative corr (-15%)
#BATF::JUN positive corr (21%)
write.csv(motifcor_telo, "/home/mahogny/jupyter/telomere_paper/pbmc/corr/cor_motif_cd4.csv")

###### nTelo axis  -- within CD8 T 
axis_telo <- adata$rank_norm_telo
use_sub <- adata$rna_clusters %in% c("r3","r5","r6","r13")
genecor_telo <- compare_gene_with_axis(adata[,use_sub], axis_telo[use_sub])
rbind(head(genecor_telo,n=20),tail(genecor_telo, n=20))
###new. still holds
# cor       gene
# 634  -0.3556189       LEF1 * 
# 421  -0.3522123      BACH2 *
# 1404 -0.3178975      PDE3B
# 1654 -0.3138543        TXK
# 408  -0.2869429      NELL2
# 1971 -0.2723586 AL008638.1
# 1488 -0.2706669  LINC01865
# 456  -0.2636072      MARK1
# 303  -0.2612516      NRCAM
# 491  -0.2499166      ROBO1
# 643  -0.2482325 AL162414.1
# 749  -0.2471796      PTPRK
# 1074 -0.2447889      CAMK4
# 659  -0.2414941       SPIC
# 1772 -0.2404082      IGF1R
# 170  -0.2394886        ERG
# 1749 -0.2369103     ABLIM1
# 341  -0.2368333       FHIT
# 645  -0.2345046       CDH9
# 600  -0.2320711   LY86-AS1
# 1641  0.2529153      MS4A2
# 383   0.2530362      BAALC
# 17    0.2535954       CCL4
# 797   0.2542701       PLEK
# 1761  0.2561637       SRPX
# 1301  0.2638406     S100A4
# 207   0.2686289       ZEB2
# 285   0.2695755      KLRD1
# 305   0.2721459       CST7
# 1361  0.2762273      KLRG1
# 9     0.2784065       GNLY
# 601   0.2793171        PZP
# 98    0.2875601       PRF1
# 115   0.2885893       GZMH
# 376   0.2912276        A2M
# 1265  0.3057225     PYHIN1
# 192   0.3165828       GZMA
# 535   0.3306708     TGFBR3
# 50    0.3514048       NKG7
# 72    0.3702894       CCL5

motifcor_telo <- compare_motif_with_axis(adata[,use_sub], axis_telo[use_sub])
rbind(head(motifcor_telo),tail(motifcor_telo))
write.csv(motifcor_telo, "/home/mahogny/jupyter/telomere_paper/pbmc/corr/cor_motif_cd8.csv")
#new
#TCF7L2, LEF1 negative corr. (-40%)
#BATF3 most positive (41%)



###### nTelo axis  -- within T and NK
axis_telo <- adata$rank_norm_telo
use_sub <- adata$rna_clusters2 %in% c("r2_0","r2_2","r2_3")
genecor_telo <- compare_gene_with_axis(adata[,use_sub], axis_telo[use_sub])
rbind(head(genecor_telo,n=20),tail(genecor_telo, n=20))
# cor       gene
# 634  -0.2636125       LEF1
# 421  -0.2608387      BACH2
# 1404 -0.2237882      PDE3B
# 303  -0.1909175      NRCAM
# 408  -0.1904496      NELL2
# 1074 -0.1869140      CAMK4
# 1971 -0.1823106 AL008638.1
# 1477 -0.1820758  LINC02160
# 1747 -0.1778979      HMHB1
# 341  -0.1765482       FHIT
# 1772 -0.1742620      IGF1R
# 1150 -0.1740258        HBD
# 313  -0.1708721 AL163932.1
# 1330 -0.1706591      CASC9
# 1616 -0.1691921    RASGRF2
# 1654 -0.1650074        TXK
# 1643 -0.1631966      NR3C2
# 1749 -0.1620801     ABLIM1
# 1660 -0.1614174     CACHD1
# 1079 -0.1612326 AC005580.1
# 1706  0.1390440      KCNA5
# 98    0.1400705       PRF1
# 207   0.1422277       ZEB2
# 1761  0.1443296       SRPX
# 457   0.1455621    PPP2R2B
# 1638  0.1455854       MSLN
# 305   0.1541805       CST7
# 50    0.1606647       NKG7
# 1484  0.1641191 AL158817.1
# 1367  0.1661920      SYNE2
# 285   0.1671994      KLRD1
# 1615  0.1672195     VPREB1
# 1039  0.1752974     ARPP21
# 192   0.1808239       GZMA
# 1301  0.1841026     S100A4
# 72    0.1888682       CCL5
# 535   0.1906018     TGFBR3
# 1265  0.1937972     PYHIN1
# 1063  0.1993802     TENT5B
# 525   0.2079478       KRT5

motifcor_telo <- compare_motif_with_axis(adata[,use_sub], axis_telo[use_sub])
rbind(head(motifcor_telo),tail(motifcor_telo))
#most neg LEF1 (27%)
#most pos BATF3 (32%)



###### nTelo axis  -- within just T
axis_telo <- adata$rank_norm_telo
use_sub <- adata$rna_clusters2 %in% c("r2_0") #not updated

genecor_telo <- compare_gene_with_axis(adata[,use_sub], axis_telo[use_sub])
rbind(head(genecor_telo,n=20),tail(genecor_telo, n=20))

motifcor_telo <- compare_motif_with_axis(adata[,use_sub], axis_telo[use_sub])
rbind(head(motifcor_telo),tail(motifcor_telo))


###### nTelo axis  -- within monocytes
axis_telo <- adata$rank_norm_telo
use_sub <- adata$rna_clusters3 %in% c("r3_19", "r3_13", "r3_6", "r3_9", "r3_4", "r3_10", "r3_5")
#use_sub <- grep("Monocyte",adata$pred.hpca)   # %in% c("r2_1")
#use_sub <- grep("monocytes",adata$pred.monaco)   # %in% c("r2_1")
#use_sub <- adata$rna_clusters2 %in% c("r2_1")
genecor_telo <- compare_gene_with_axis(adata[,use_sub], axis_telo[use_sub])
rbind(head(genecor_telo,n=20),tail(genecor_telo, n=20))

head(genecor_telo[genecor_telo$sumexp>100,],n=20)
tail(genecor_telo[genecor_telo$sumexp>100,],n=20)
genecor_telo[genecor_telo$gene=="",]


FeaturePlot(object = adata[,use_sub], reduction = 'UMAP_RNA', 
            features = c("rank_norm_telo","SLCO4C1","PRICKLE1","PSAP","HLA-DRA","CD74","CST3",
                         "CDKN1C"))
#HLA-DRA, HLA-DRB5, HLA-DPA1, CD74,   etc. 
plot_umap_gene(adata,"UMAP_RNA_3D","norm_telo")
plot_umap_gene(adata,"UMAP_RNA_3D","rank_norm_telo")
plot_umap_gene(adata,"UMAP_RNA_3D","nCount_RNA")
plot_umap_gene(adata,"UMAP_RNA_3D","nucleosome_signal")
plot_umap_gene(adata,"UMAP_RNA_3D","SLCO4C1")
plot_umap_gene(adata,"UMAP_RNA_3D","PRICKLE1")
plot_umap_gene(adata,"UMAP_RNA_3D","IL31RA")

plot_umap_groupby(adata,"UMAP_RNA_3D","Phase","rna_clusters3") 
plot_umap_groupby(adata,"UMAP_RNA_3D","rna_clusters3","Phase") 


FeaturePlot(object = adata[,use_sub], reduction = 'UMAP_RNA', features = c("IGHV1-2","CDC45","IGKV1-39","SCN3A","TBX10","IL1A"))
FeaturePlot(object = adata[,use_sub], reduction = 'UMAP_RNA', features = c("rank_norm_telo","CDKN1C"))
#adata$rank_norm_telo
genecor_telo[genecor_telo$gene=="FCGR3A",]



## old
# cor       gene
# 85  -0.1890237       VCAN
# 340 -0.1876373      CSF3R
# 463 -0.1854591       DYSF
# 725 -0.1810160       NCF1
# 433 -0.1731174   ARHGAP26
# 991 -0.1699722       DPYD
# 176 -0.1662971      PLCB1
# 289 -0.1535904      RBM47
# 238 -0.1493702   ARHGAP24
# 188 -0.1459236     PLXDC2
# 451 -0.1432703      PADI4
# 175 -0.1420926       CD36
# 241 -0.1408704      LRMDA
# 788 -0.1376078     JARID2
# 300 -0.1316663       GAB2
# 612 -0.1314273      MEGF9
# 71  -0.1314020     S100A8
# 206 -0.1274630      ACSL1
# 121 -0.1260128      NAMPT
# 517 -0.1152899    SLC11A1
# 191  0.1697942     FCER1G
# 165  0.1703035       PSAP
# 13   0.1722380   HLA-DPB1
# 366  0.1741609 AC104809.2
# 404  0.1749029      CSF1R
# 84   0.1772764   HLA-DRB5
# 580  0.1776707       RHOC
# 38   0.1806031   HLA-DQA1
# 213  0.1810147       AIF1
# 691  0.1827384       PFN1
# 218  0.1833545       HES4
# 18   0.1880104    HLA-DRA
# 232  0.1923366      IFI30
# 246  0.1977650       CTSS
# 25   0.2145828   HLA-DPA1
# 201  0.2222134       LST1
# 145  0.2251262      MTSS1
# 21   0.2320239       CST3
# 34   0.2429956     CDKN1C  ### Cyclin Dependent Kinase Inhibitor 1C
# 31   0.2449951     FCGR3A  **CD16

motifcor_telo <- compare_motif_with_axis(adata[,use_sub], axis_telo[use_sub])
rbind(head(motifcor_telo),tail(motifcor_telo))
#Zic1::Zic2 most negative (-27%)
#Klf12 most positive (12%)


#part of CD4
axis_telo <- adata$rank_norm_telo
use_sub <- adata$rna_clusters4 %in% c("r4_9","r4_0")
genecor_telo <- compare_gene_with_axis(adata[,use_sub], axis_telo[use_sub])
rbind(head(genecor_telo,n=20),tail(genecor_telo, n=20))
#IL32 positive. CDK6 6% negative. 

motifcor_telo <- compare_motif_with_axis(adata[,use_sub], axis_telo[use_sub])
rbind(head(motifcor_telo),tail(motifcor_telo))


########################################################################
################ Plot motif activity ###################################   JASPAR motifs
########################################################################


DefaultAssay(adata) <- 'chromvar'

plot_umap_motif <- function(usered, motifname, textby, dorank=FALSE){
  themotif <- motif2names$id[motif2names$name==motifname][1]
  print(themotif)
  DefaultAssay(adata) <- 'chromvar'
  plot_umap_gene(adata,usered, themotif, textby, dorank=dorank)
}  

plot_umap_motif("UMAP_RNA_3D", "CTCF","rna_clusters3", dorank=TRUE)
plot_umap_motif("UMAP_ATAC_3D", "CTCF","rna_clusters3", dorank=TRUE)


DefaultAssay(adata) <- 'chromvar'
FeaturePlot(adata, reduction = "UMAP_RNA", features = motif2names$id[motif2names$name=="CTCF"][1])



########################################################################
########### UMAP comparing TF activities ###############################
########################################################################

DefaultAssay(adata) <- 'chromvar'

#Use seurat to compare motif activities
DefaultAssay(adata) <- 'chromvar'
tmp <- as.matrix(adata@assays$chromvar@data)
tmp[is.na(tmp)]<-0 #or remove these entries?
#rownames(tmp) <- motif2names[rownames(tmp),]$name  ### later
motif_adata <- CreateSeuratObject(counts = t(tmp))
remove(tmp)

motif_adata <- ScaleData(motif_adata, features = rownames(motif_adata))
motif_adata <- RunPCA(motif_adata,features = rownames(motif_adata), verbose = FALSE)
motif_adata <- FindNeighbors(motif_adata, dims = 1:30)
motif_adata <- FindClusters(motif_adata, resolution = 1)
motif_adata <- RunUMAP(object = motif_adata, dims = 1:30, n.neighbors = 20)

DimPlot(object = motif_adata, pt.size = 1, label = TRUE) + ggtitle('Motif activity similarity')


grouped_motifs <- data.frame(
  x=Embeddings(motif_adata, reduction = "umap")[,1],
  y=Embeddings(motif_adata, reduction = "umap")[,2],
  motif=names(motif_adata$seurat_clusters),
  cluster=motif_adata$seurat_clusters)
grouped_motifs <- grouped_motifs[order(grouped_motifs$cluster),]
grouped_motifs$name <- motif2names[grouped_motifs$motif,]$name

grouped_motifs$newname <- grouped_motifs$name
grouped_motifs$newname[grouped_motifs$name %!in% c(
  "XBP1",
  "PPARG"
  )] <- ""

ggplot(grouped_motifs, aes(x=x, y=y, color=cluster, label=newname))+
  geom_point()+
  #geom_point(color="gray")+
  geom_text(color="black")

#Label DE genes; top markers... or use highly variable score?
grouped_motifs$isDE <- grouped_motifs$name %in% cluster_markers_rna3$gene[cluster_markers_rna3$p_val<1e-100]
#grouped_motifs$isDE <- grouped_motifs$name %in% top10_rna3$gene



########### Find which motifs belong to a highly variable gene
map_jaspar_uniprot <- read.csv('/home/mahogny/jupyter/tcellpaper/jaspar_uniprot.csv', header = T, sep = ',')[c(3,5,6)]
colnames(map_jaspar_uniprot) <- c("jaspar_baseid","jasparname","uniprot")
map_jaspar_uniprot <- map_jaspar_uniprot[map_jaspar_uniprot$uniprot!="",]
map_uniprot_symbol <- read.csv('/home/mahogny/jupyter/tcellpaper/uniprot_symbol_ensid.csv', header = T, sep = '\t')
colnames(map_uniprot_symbol) <- c('ensid','uniprot','symbol')

DefaultAssay(adata) <- 'RNA'
map_uniprot_symbol <- map_uniprot_symbol[map_uniprot_symbol$symbol %in% rownames(adata),]
merged_jaspar <- unique(merge(map_uniprot_symbol, map_jaspar_uniprot))
merged_jaspar <- merge(
  merged_jaspar,
  data.frame(
    jaspar_baseid=str_split_fixed(colnames(motif_adata),"\\.",2)[,1],
    jaspar_id=colnames(motif_adata))
)



DefaultAssay(adata) <- 'RNA'
grouped_motifs$is_variable <- grouped_motifs$motif %in% unique(merged_jaspar[merged_jaspar$symbol %in% VariableFeatures(adata),]$jaspar_id)
grouped_motifs[grouped_motifs$is_variable,]
grouped_motifs$varname <- grouped_motifs$name
grouped_motifs$varname[!grouped_motifs$is_variable] <- ""

#ggplot(grouped_motifs, aes(x=x, y=y, label=varname))+
#  geom_point(color="gray")+geom_text()

ggplot(grouped_motifs, aes(x=x, y=y, color=cluster, label=varname))+
  geom_point()+
  #geom_point(color="gray")+
  geom_text(color="black")


########################################################################
############# What is close to CTCF? ###################################
########################################################################


## ugly way to get an overview
grouped_motifs$varname[grouped_motifs$name=="CTCF"] <- "CTCF"
grouped_motifs$varname[grouped_motifs$name=="RBPJ"] <- "RBPJ"

ggplot(grouped_motifs, aes(x=x, y=y, color=cluster, label=varname))+
  geom_point()+
  #geom_point(color="gray")+
  geom_text(color="black")



#################### Which motifs are in peaks with CTCF?



################## Which motifs co-vary with CTCF activity?


DefaultAssay(adata) <- 'chromvar'
motifid_ctcf <- motif2names$id[motif2names$name=="CTCF"]
cortf <- NULL
tfvals <- as.double(adata@assays$chromvar[str_replace_all(motifid_ctcf,"_","-"),])
#tfvals <- as.double(adata@assays$chromvar[motifid_ctcf,])
for(i in 1:nrow(adata)){
  print(i)
  cortf <- c(cortf, cor(tfvals,as.double(adata@assays$chromvar[i,])))
}
cortf <- data.frame(
  #id=rownames(adata),
  id=str_replace_all(rownames(adata),"-","_"),
  cor=cortf
)
cortf <- merge(cortf, motif2names)
cortf <- cortf[!is.na(cortf$cor),]
cortf <- cortf[order(cortf$cor),]
tail(cortf,n=30)
head(cortf,n=30)

######### JASPAR ####
#Most negative: 
#Most positive: 

plot(cortf$cor)
cortf[cortf$name=="RBPJ",] #no corr



########################################################################
################ SCVI analysis #########################################
########################################################################



Sys.setenv(RETICULATE_PYTHON = here::here("/home/mahogny/miniconda3/bin/python"))
sc <- import("scanpy", convert = FALSE)
scvi <- import("scvi", convert = FALSE)

##SCVI need raw rounts; so reload the data and subset to the same size
counts <- Read10X_h5("/corgi/sebas/tcell_multi/aggregated_tlibs/AGG6789_good/outs/filtered_feature_bc_matrix.h5")
raw_adata <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)
remove(counts)
raw_adata <- raw_adata[VariableFeatures(adata),colnames(adata)] #Save time by only using variable features

py_adata <- convertFormat(raw_adata, from="seurat", to="anndata", main_layer="counts", drop_single_values=FALSE)
print(py_adata) # Note generally in Python, dataset conventions are obs x var


#Transfer relevant metadata
themeta <- py_to_r(py_adata$obs)
themeta$S_score <- adata$S.Score
themeta$G2M_score <- adata$G2M.Score
py_adata$obs <- r_to_py(themeta)

# run setup_anndata
scvi$model$SCVI$setup_anndata(
  py_adata,
  continuous_covariate_keys = c("S_score","G2M_score"))

#TODO Training will be faster when sparse matrix is formatted as CSR  ---- find a way to convert?

# categorical_covariate_keys=["cell_source", "donor"],
# continuous_covariate_keys=["percent_mito", "percent_ribo"]

### TODO can feed cell cycle phase scores as covariates

# create the model and train it
model <- scvi$model$SCVI(py_adata)
model$train()

# to specify the number of epochs when training:
# model$train(max_epochs = as.integer(400))

# get the latent represenation and put it back in our original Seurat object
latent <- as.matrix(model$get_latent_representation())
rownames(latent) <- colnames(adata)
adata[["scvi"]] <- CreateDimReducObject(embeddings = latent, key = "scvi_", assay = DefaultAssay(adata))


# Find clusters, then run UMAP, and visualize
adata <- FindNeighbors(adata, dims = 1:10, reduction = "scvi")
adata <- FindClusters(adata, resolution =2.1) 
adata$scvi_rna_clusters <- adata$seurat_clusters
adata$scvi_rna_clusters <- factor(sprintf("sr%s",adata$scvi_rna_clusters))
adata <- RunUMAP(adata, dims = 1:10, reduction = "scvi", n.components = 2, reduction.name="UMAP_SCVI_RNA")
adata <- RunUMAP(adata, dims = 1:10, reduction = "scvi", n.components = 3, reduction.name="UMAP_SCVI_RNA_3D")

#DimPlot(adata, reduction = "UMAP_SCVI_RNA", pt.size = 1)

plot_umap_groupby(adata,"UMAP_SCVI_RNA_3D","Phase","rna_clusters3", colors = c("red","green","blue"))

plot_umap_groupby(adata,"UMAP_SCVI_RNA_3D","rna_clusters","Phase")
plot_umap_groupby(adata,"UMAP_SCVI_RNA_3D","rna_clusters3","Phase")
plot_umap_groupby(adata,"UMAP_SCVI_RNA_3D","scvi_rna_clusters","Phase")

plot_umap_groupby(adata,"UMAP_SCVI_RNA_3D","pred.dice","Phase")

plot_umap_groupby(adata,"UMAP_ATAC_3D","scvi_rna_clusters","Phase")


################### Compare RNA-based clusters - fine
DefaultAssay(adata) <- "RNA"
Idents(adata) <- adata$scvi_rna_clusters
cluster_markers_scvirna1 <- FindAllMarkers(adata, only.pos = F, test.use = 'wilcox')
cluster_markers_scvirna1 %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) -> top10_scvirna


#This defines an alternative cluster next to SELL. it is all over in the ATACseq data
plot_umap_gene(adata,"UMAP_SCVI_RNA_3D","ESYT2","scvi_rna_clusters")  #everewhere but peaks here. FGF linked. only t cell paper: https://pubmed.ncbi.nlm.nih.gov/32879390/  CRAC channel
plot_umap_gene(adata,"UMAP_SCVI_RNA_3D","LZTFL1","scvi_rna_clusters") #highest but broad. splits it all up in another direction
plot_umap_gene(adata,"UMAP_SCVI_RNA_3D","LDLRAD4","scvi_rna_clusters") #like above
plot_umap_gene(adata,"UMAP_SCVI_RNA_3D","BTBD11","scvi_rna_clusters") #in yet a different way
#plot_umap_gene(adata,"UMAP_SCVI_RNA_3D","ESYT2","scvi_rna_clusters") #everewhere but peaks here


plot_umap_gene(adata,"UMAP_SCVI_RNA_3D","ARID5B","scvi_rna_clusters") #CC-origin?

plot_umap_gene(adata,"UMAP_SCVI_RNA_3D","MALAT1","scvi_rna_clusters") #annoyingly high in one naive like cluster. but not only

#### naivish cells are similar:
plot_umap_gene(adata,"UMAP_SCVI_RNA_3D","TCF7","scvi_rna_clusters") #early in the naive cascade
plot_umap_gene(adata,"UMAP_SCVI_RNA_3D","CCR7","scvi_rna_clusters") #in half a huge blob
plot_umap_gene(adata,"UMAP_SCVI_RNA_3D","SELL","scvi_rna_clusters") 
plot_umap_gene(adata,"UMAP_SCVI_RNA_3D","IL7R","scvi_rna_clusters") 
plot_umap_gene(adata,"UMAP_SCVI_RNA_3D","LEF1","scvi_rna_clusters")
plot_umap_gene(adata,"UMAP_SCVI_RNA_3D","RBPJ","scvi_rna_clusters")  #in all but naive! some in big blob

#read https://www.nature.com/articles/s41577-021-00563-6 !! RBPJ and TCF7


### previous scATAC of human T cell development; should cite
# https://www.cell.com/immunity/fulltext/S1074-7613(20)30465-9?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS1074761320304659%3Fshowall%3Dtrue



plot_umap_gene(adata,"UMAP_RNA_3D","CCR7","rna_clusters3") 



if(FALSE){
  plot_umap_gene(adata,"UMAP_SCVI_RNA_3D","JUNB","rna_clusters")
}



########################################################################
################ Final UMAP plots ###################################### For T cells
########################################################################



# bdata <- adata[,adata$rna_clusters2 %in% c("r2_0","r2_1") &
#                  adata$pred.dice %!in% c("B cells, naive") &
#                  adata$pred.monaco %!in% c("Progenitor cells")]

#r2 2 0 3

bdata <- adata[,adata$rna_clusters2 %in% c("r2_0","r2_2","r2_3") &
                 #adata$pred.dice %!in% c("B cells, naive") &
                 adata$pred.monaco %!in% c("Progenitor cells","Intermediate monocytes","MAIT cells","Switched memory B cells","Naive B cells") &
                 (str_starts(adata$pred.dice,"T cells") | str_starts(adata$pred.dice, "NK cells"))]



#bdata <- FindVariableFeatures(bdata, selection.method = 'vst',binning.method = 'equal_frequency', verbose = T, nfeatures = 1000)
#bdata <- ScaleData(bdata, features = rownames(adata)) #default is to only use highly variable features
#bdata <- RunPCA(bdata,reduction.name = 'PCA_RNA',ndims.print = c(1,3,5), verbose = T, features = )
#DefaultAssay(adata) <- "RNA"
#HVFInfo(adata)

bdata <- RunUMAP(object = bdata, dims = 1:30, reduction = "PCA_RNA", reduction.name="UMAP_RNA")

# FeaturePlot(object = bdata, reduction = 'UMAP_RNA', features = c("rank_norm_telo"))
# FeaturePlot(object = bdata, reduction = 'UMAP_RNA', features = c("CD8A","CD4","NKG7"))
# DimPlot(object = bdata, reduction = 'UMAP_RNA', group.by = "Phase", label = TRUE)
# DimPlot(object = bdata, reduction = 'UMAP_RNA', group.by = "pred.monaco", label = TRUE)
# DimPlot(object = bdata, reduction = 'UMAP_RNA', group.by = "pred.dice", label = TRUE)

#better to show ATAC UMAP
styleGG <- function(p){
  p+#xlim(-5,12)+ylim(-10,10)+
    xlab("UMAP1")+ylab("UMAP2")+ 
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
}
p1 <- styleGG(FeaturePlot(raster=TRUE, object = bdata, reduction = 'UMAP_RNA', features = c("rank_norm_telo"))+ggtitle("rank(nTA)"))
p2 <- styleGG(DimPlot(raster=TRUE, object = bdata, reduction = 'UMAP_RNA', group.by = "Phase"))
p3 <- styleGG(DimPlot(raster=TRUE, object = bdata, reduction = 'UMAP_RNA', group.by = "pred.monaco")+ggtitle("Cell type"))#, label = TRUE))  #Monaco is better than DICE
p4 <- styleGG(FeaturePlot(raster=TRUE, object = bdata, reduction = 'UMAP_RNA', features = c("BACH2")))
p5 <- styleGG(FeaturePlot(raster=TRUE, object = bdata, reduction = 'UMAP_RNA', features = c("LEF1")))
p6 <- styleGG(FeaturePlot(raster=TRUE, object = bdata, reduction = 'UMAP_RNA', features = c("CCL5")))
p7 <- styleGG(FeaturePlot(raster=TRUE, object = bdata, reduction = 'UMAP_RNA', features = c("GZMA")))
p8 <- styleGG(FeaturePlot(raster=TRUE, object = bdata, reduction = 'UMAP_RNA', features = c("GZMH")))
p4/p5|p6/p7|p8/p1|p2/p3
ggsave("~/tcell.pdf", width = 15)




#nTA also up in T cells, CD4+, memory TREG (DICE agrees)
#What is up in DICE, T cells, CD4+, memory TREG?
DefaultAssay(bdata) <- "RNA" 
Idents(bdata) <- bdata$pred.dice
cluster_markers_dice <- FindAllMarkers(bdata, only.pos = F, test.use = 'wilcox')
cluster_markers_dice %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) -> top10_dice


styleGG(DimPlot(object = bdata, reduction = 'UMAP_RNA', group.by = "pred.dice", label = TRUE))
DimPlot(object = bdata, reduction = 'UMAP_RNA', group.by = "pred.dice")#, label = TRUE)
(DimPlot(object = bdata, reduction = 'UMAP_RNA', group.by = "pred.monaco", label = TRUE))#))
(DimPlot(object = bdata, reduction = 'UMAP_RNA', group.by = "pred.dice", label = TRUE))#, label = TRUE))



########################################################################
################ Final UMAP plots ###################################### For monocytes
########################################################################

bdata <- adata[,adata$rna_clusters2 %in% c("r2_0","r2_1") &  #should be 1 and 5, or just 1
                 adata$pred.dice %!in% c("B cells, naive") &
                 adata$pred.monaco %in% c("Classical monocytes","Intermediate monocytes","Non classical monocytes") &
                 str_starts(adata$pred.dice,"Monocytes")]


styleGG <- function(p){
  p+xlim(-11,-4)+ylim(-7,5)+
    xlab("UMAP1")+ylab("UMAP2")+ 
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
}
p1 <- styleGG(FeaturePlot(raster=TRUE, object = bdata, reduction = 'UMAP_RNA', features = c("rank_norm_telo"))) 
p2 <- styleGG(FeaturePlot(raster=TRUE, object = bdata, reduction = 'UMAP_RNA', features = c("FCGR3A")))  #CD16
p3 <- styleGG(FeaturePlot(raster=TRUE, object = bdata, reduction = 'UMAP_RNA', features = c("CD14")))  #CD16
p4 <- styleGG(FeaturePlot(raster=TRUE, object = bdata, reduction = 'UMAP_RNA', features = c("CDKN1C")))
p5 <- styleGG(DimPlot(raster=TRUE, object = bdata, reduction = 'UMAP_RNA', group.by = "Phase"))
p6 <- styleGG(DimPlot(raster=TRUE, object = bdata, reduction = 'UMAP_RNA', group.by = "pred.monaco")+ggtitle("Cell type"))
p1/p2|p3/p4|p5/p6
ggsave("~/monocyte.pdf")




########################################################################
############ Possibility of telomere length measurement? ###############
########################################################################


#in newly integrated pbmc, these correlate better
plot(
  adata$nucleosome_signal, 
  adata$norm_telo,pch=19,cex=0.1)  

plot(
  adata$nucleosome_signal, 
  adata$norm_telo/adata$nucleosome_signal,pch=19,cex=0.1)  

plot(
  adata$nucleosome_signal,
  adata$dedupcnt_telo/(adata$totalcnt_telo+1)/adata$nucleosome_signal,pch=19,cex=0.1)

plot(
  adata$nucleosome_signal,
  1/(adata$dedupcnt_telo/(adata$totalcnt_telo+1)/adata$nucleosome_signal),pch=19,cex=0.1)


plot(
  adata$nucleosome_signal,
  log(1/(adata$dedupcnt_telo/(adata$totalcnt_telo+1)/adata$nucleosome_signal)),pch=19,cex=0.1)


# x <- adata$nucleosome_signal-mean(adata$nucleosome_signal)
# y <- (adata$norm_telo)/(adata$nucleosome_signal-mean(adata$nucleosome_signal))
# plot(
#   x[x!=0], 
#   y[x!=0],pch=19,cex=0.1)  

plot(
  adata$nucleosome_signal, 
  exp(adata$norm_telo),pch=19,cex=0.1)  #the exp makes little difference



themod <- lm(
  norm_telo ~ nucleosome_signal,
  data.frame(
    nucleosome_signal=adata$nucleosome_signal,
    norm_telo=adata$norm_telo))

themod <- lm(
  norm_telo ~ nucleosome_signal,
  data.frame(
    nucleosome_signal=adata$nucleosome_signal,
    norm_telo=exp(adata$norm_telo)))

adata$red_telo <- themod$residuals
adata$rank_red_telo <- rank(adata$red_telo)

plot(
  adata$nucleosome_signal,
  adata$red_telo, 
  pch=19,cex=0.1)  #the exp makes little difference
plot(
  adata$nucleosome_signal, 
  adata$norm_telo,pch=19,cex=0.1)  #the exp makes little difference

FeaturePlot(object = adata, reduction = 'UMAP_RNA', features = c("norm_telo","red_telo"))
FeaturePlot(object = adata, reduction = 'UMAP_RNA', features = c("rank_norm_telo","rank_red_telo"))
#FeaturePlot(object = adata, reduction = 'UMAP_RNA', features = c("red_telo"))
  


########################################################################
################ Deduplication analysis for paper ######################
########################################################################



#How many reads per cell?
adata$rawcnt_telo
mean(adata$nCount_RNA) #3804
mean(adata$nCount_peaks) #25769


#Depuplication stats
map_lib2ds <- data.frame(
  row.names=c("1","2","3","4"),
  name=c(
  "10k_PBMC_Multiome_nextgem_Chromium_Controller",
  "10k_PBMC_Multiome_nextgem_Chromium_X",
  "pbmc_granulocyte_sorted_10k",  #got more cells than 1,2
  "pbmc_granulocyte_unsorted_10k"))
adf <- data.frame(
  raw=log10(1+adata$rawcnt_telo),
  dedup=log10(1+adata$dedupcnt_telo),
  color=map_lib2ds[adata$lib,,drop=FALSE]$name
)
ggplot(
  adf[sample(1:nrow(adf), size=20000),],
  aes(raw,dedup,col=color)
) + geom_point(size=0.1) +xlab("Log10 1+Raw count") +ylab("Log10 1+Deduplicate count") +xlim(0.5,3)+ylim(0.5,3)
ggsave("~/pbmc_telo_dedup.pdf", width = 7, height = 3)

sum(adata$dedupcnt_telo, na.rm = TRUE) / sum(adata$rawcnt_telo, na.rm = TRUE)
# 0.4043948






