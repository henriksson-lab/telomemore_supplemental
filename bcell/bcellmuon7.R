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
counts <- Read10X_h5("/corgi/cellbuster/bigb/aggr/outs/filtered_feature_bc_matrix.h5")
fragpath <- "/corgi/cellbuster/bigb/aggr/outs/atac_fragments.tsv.gz"

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"

# create a Seurat object containing the RNA adata
adata <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# create ATAC assay and add it to the object
adata[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)

remove(counts)


## Barcode info for multiome 10x
bc10x <- data.frame(
  atac=read.table("/data/henlab/software/cellranger-arc-2.0.0/lib/python/atac/barcodes/737K-arc-v1.txt.gz")[,1],
  rna=read.table("/data/henlab/software/cellranger-arc-2.0.0/lib/python/cellranger/barcodes/737K-arc-v1.txt.gz")[,1]
)

## Donor information
donorinfo <- read.csv("/corgi/cellbuster/bigb/all_matched_donor_ids.tsv")
rownames(donorinfo) <- sprintf("%s-%s",str_sub(donorinfo$cell,1,16),as.integer(str_sub(donorinfo$batch,4)))  
sum(rownames(donorinfo) %in% names(adata$orig.ident))
donorinfo <- donorinfo[names(adata$orig.ident),]
adata$prob_doublet <- donorinfo$prob_doublet
adata$donor_id <- donorinfo$donor_id
adata$lib <- donorinfo$batch

## telomemore information
teloinfo <- read.csv("/corgi/cellbuster/bigb/summary_kmer.java.csv")
rownames(teloinfo) <- sprintf("%s-%s",str_sub(teloinfo$barcode,1,16),as.integer(str_sub(teloinfo$dataset,4))) 
sum(rownames(teloinfo) %in% names(adata$orig.ident))
teloinfo <- teloinfo[names(adata$orig.ident),]
adata$rawcnt_telo <- teloinfo$totalcnt_CCCTAA
adata$dedupcnt_telo <- teloinfo$dedupcnt_CCCTAA
adata$totalcnt_telo <- teloinfo$total
adata$norm_telo <- log10(adata$dedupcnt_telo/(adata$totalcnt_telo+1)+1e-3)
hist(adata$norm_telo)


adata$rank_norm_telo <- NA
for(i in unique(adata$lib)){
  theranks <- rank(adata$norm_telo[adata$lib==i])
  adata$rank_norm_telo[adata$lib==i] <- theranks/max(theranks)
}

#rank telo per dataset


## TRUST4 info
adata$trust4 <- ""
trust4 <- NULL
for(onelib in 1:4){ 
  dat <- read.csv(sprintf("/corgi/cellbuster/bigb/trust4/lib%s/TRUST_gex_possorted_bam_barcode_report.tsv",onelib),sep="\t")
  trust4 <- rbind(trust4, data.frame(
    bc=sprintf("%s-%s",str_split_fixed(dat$X.barcode,"-",2)[,1],onelib),
    celltype=dat$cell_type
  ))
}
rownames(trust4) <- trust4$bc
adata$trust4 <- trust4[names(adata$orig.ident),]$celltype
adata$trust4[is.na(adata$trust4)] <- ""


## Gather annotation previously transferred from King et al
king_annot <- read.csv("/corgi/cellbuster/bigb/aggr/annot_from_king.csv")
rownames(king_annot) <- king_annot$X
king_annot <- king_annot[,-1]
adata@meta.data <- cbind(adata@meta.data,king_annot[colnames(adata),])

## Only keep B cells
adata <- adata[,adata$king_lineage=="B Cells"]
  
## Filtering by vireo
adata <- adata[,adata$donor_id %!in% c("doublet","unassigned")]


################################################################################
################### RNA basic processing #######################################
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
               subset = percent.mt < 25 & 
                 nFeature_RNA > 1000 & nFeature_RNA < 7000 & 
                 nCount_RNA > 500 & nCount_RNA<20000)

## Normalization
adata<- NormalizeData(adata,normalization.method = 'LogNormalize')

#cell cycle scoring
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
adata <- CellCycleScoring(adata, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
if(FALSE){
  #Possible CC removal
  adata <- ScaleData(adata, vars.to.regress = c('S.Score','G2M.Score'))
}


## HIGHLY VARIABLE FEATURE ANALYSIS. gets 850 with previous automatic method
adata <- FindVariableFeatures(adata, selection.method = 'vst',binning.method = 'equal_frequency', verbose = T, nfeatures = 1000)
adata <- ScaleData(adata, features = rownames(adata)) #default is to only use highly variable features

#Ensure AICDA is here
length(adata@assays$RNA@var.features)
"AICDA" %in% adata@assays$RNA@var.features

#Summary scores for Ig genes
adata$IGHD_sum <- colSums(adata@assays$RNA@scale.data[setdiff(rownames(adata)[str_starts(rownames(adata),"IGHD")],"IGHD"),])
adata$IGHV_sum <- colSums(adata@assays$RNA@scale.data[setdiff(rownames(adata)[str_starts(rownames(adata),"IGHV")],"IGHV"),])
adata$IGHG_sum <- colSums(adata@assays$RNA@scale.data[setdiff(rownames(adata)[str_starts(rownames(adata),"IGHG")],"IGHG"),])
adata$IGHA_sum <- colSums(adata@assays$RNA@scale.data[setdiff(rownames(adata)[str_starts(rownames(adata),"IGHA")],"IGHA"),])


################################################################################
################### ATAC basic processing ######################################
################################################################################


DefaultAssay(adata) <- "ATAC"

adata <- NucleosomeSignal(adata)
adata <- TSSEnrichment(adata, fast=FALSE) 

## QC plots
VlnPlot(
  object = adata,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
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
  subset = nCount_ATAC < 100e03 & #very loose. half is fine
    nucleosome_signal < 2 &
    TSS.enrichment > 1 #very loose
)

## peak calling
peaks <- CallPeaks(adata, macs2.path = '/home/mahogny/miniconda3/bin/macs2')

## Remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

## quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(adata),
  features = peaks,
  cells = colnames(adata), 
  verbose = T
)

## create a new assay using the MACS2 peak set and add it to the Seurat object
adata[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = annotation
)




################################################################################
################### Load motif data ############################################
################################################################################


DefaultAssay(adata) <- "peaks"

if(TRUE){
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


DefaultAssay(adata) <- "RNA"

hpca <-   celldex::HumanPrimaryCellAtlasData()
pred.HPCA <- SingleR(test = as.SingleCellExperiment(adata), ref = hpca, assay.type.test=1,
                     labels = hpca$label.fine)
adata$pred.hpca <- pred.HPCA$labels

monaco <- celldex::MonacoImmuneData()
pred.monaco <- SingleR(test = as.SingleCellExperiment(adata), ref = monaco, assay.type.test=1,
                       labels = monaco$label.fine)
adata$pred.monaco <- pred.monaco$labels

## Filtering by HPCA
adata <- adata[,str_starts(adata$pred.hpca,"B_cell")]

## Filtering by Monaco. cluster is detected as doublets to great extent as well; likely technical reasons
adata <- adata[,adata$pred.monaco %!in% c("MAIT cells", "Follicular helper T cells", "Non-Vd2 gd T cells","Plasmacytoid dendritic cells")]
#unique(adata$pred.monaco)



########################################################################
################# RNAseq dim red & clustering ##########################
########################################################################


DefaultAssay(adata) <- "RNA"

## Dimensional reduction
adata <- RunPCA(adata,reduction.name = 'PCA_RNA',ndims.print = c(1,3,5), verbose = T)
adata <- RunUMAP(object = adata, dims = 1:30, reduction = "PCA_RNA", reduction.name="UMAP_RNA")
adata <- RunUMAP(adata, dims = 1:30,reduction = 'PCA_RNA',n.components = 3L, reduction.name ='UMAP_RNA_3D' )
adata <- FindNeighbors(adata, dims = 1:30,reduction = 'PCA_RNA')

#DimPlot(adata,dims=c(2,3))
#DimHeatmap(adata,dims = c(1,2,3,4), reduction="PCA_RNA")

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
adata <- FindClusters(adata, resolution = 1.3)
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
#The first 3 components correlate with depth


#HVFInfo(adata)

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


################### Compare RNA-based clusters - fine
DefaultAssay(adata) <- "RNA"
Idents(adata) <- adata$rna_clusters
cluster_markers_rna1 <- FindAllMarkers(adata, only.pos = F, test.use = 'wilcox')
cluster_markers_rna1 %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) -> top10_rna

#DoHeatmap(adata, features = top10_rna$gene, size=2) #+ NoLegend()

################### Compare RNA-based clusters - 2 clusters
DefaultAssay(adata) <- "RNA"
Idents(adata) <- adata$rna_clusters2
cluster_markers_rna2 <- FindAllMarkers(adata, only.pos = F, test.use = 'wilcox')
cluster_markers_rna2 %>% group_by(cluster) %>% top_n(n = 40, wt = avg_log2FC) -> top10_rna2

#DoHeatmap(adata, features = top10_rna$gene, size=2) #+ NoLegend()

################### Compare ATAC-based clusters
DefaultAssay(adata) <- "RNA" #yes!
Idents(adata) <- adata$atac_clusters
cluster_markers_atac1 <- FindAllMarkers(adata, only.pos = F, test.use = 'wilcox', )

cluster_markers_atac1 %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) -> top10_atac

#DoHeatmap(adata, features = top10_atac$gene) #+ NoLegend()




################### Compare RNA-based clusters - fine
DefaultAssay(adata) <- "RNA"
Idents(adata) <- adata$rna_clusters3
cluster_markers_rna3 <- FindAllMarkers(adata, only.pos = F, test.use = 'wilcox')
cluster_markers_rna3 %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) -> top10_rna3

################### Compare RNA-based clusters - fine
DefaultAssay(adata) <- "RNA"
Idents(adata) <- adata$rna_clusters4
cluster_markers_rna4 <- FindAllMarkers(adata, only.pos = F, test.use = 'wilcox')
cluster_markers_rna4 %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) -> top10_rna4



################### Compare two memory RNA-based clusters 
DefaultAssay(adata) <- "RNA"
Idents(adata) <- adata$rna_clusters4
cluster_markers_rna_twomem <- FindAllMarkers(adata[,adata$rna_clusters4 %in% c("r4_1","r4_2")], only.pos = F, test.use = 'wilcox')
cluster_markers_rna_twomem %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) -> top10_rna_twomem


################### Compare three naive RNA-based clusters 
DefaultAssay(adata) <- "RNA"
Idents(adata) <- adata$rna_clusters4
cluster_markers_rna_naive <- FindAllMarkers(adata[,adata$rna_clusters4 %in% c("r4_0","r4_3","r4_5")], only.pos = F, test.use = 'wilcox')
cluster_markers_rna_naive %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) -> top10_rna_naive



################### Compare GC RNA-based clusters 
DefaultAssay(adata) <- "RNA"
Idents(adata) <- adata$rna_clusters4
cluster_markers_rna_gc <- FindAllMarkers(adata[,adata$rna_clusters4 %in% c("r4_8","r4_10","r4_11","r4_9","r4_6")], only.pos = F, test.use = 'wilcox')
cluster_markers_rna_gc %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) -> top10_rna_gc


################### Compare ATAC-based clusters - fine
DefaultAssay(adata) <- "RNA" #yes!
Idents(adata) <- adata$atac_clusters3
cluster_markers_atac3 <- FindAllMarkers(adata, only.pos = F, test.use = 'wilcox', )
cluster_markers_atac3 %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) -> top10_atac3



################### Compare Phase-based clusters - fine
DefaultAssay(adata) <- "RNA"
Idents(adata) <- adata$Phase
cluster_markers_phase <- FindAllMarkers(adata, only.pos = F, test.use = 'wilcox')
cluster_markers_phase %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) -> top10_phase




################################################################################
################### Plots of metadata ##########################################
################################################################################


FeaturePlot(object = adata, reduction = 'UMAP_RNA', features = c("nCount_RNA","percent.mt","nFeature_RNA"))
FeaturePlot(object = adata, reduction = 'UMAP_RNA', features = c("nCount_ATAC", "TSS.enrichment", "nucleosome_signal"))
FeaturePlot(object = adata, reduction = 'UMAP_RNA', features = c("norm_telo"))
FeaturePlot(object = adata, reduction = 'UMAP_RNA', features = c("rank_norm_telo"))
FeaturePlot(object = adata, reduction = 'UMAP_RNA', features = c("prob_doublet"))
DimPlot(object = adata, reduction = 'UMAP_RNA', group.by = "Phase", label = TRUE)
DimPlot(object = adata, reduction = 'UMAP_RNA', group.by = "donor_id", label = TRUE)
DimPlot(object = adata, reduction = 'UMAP_RNA', group.by = "lib")
DimPlot(object = adata, reduction = 'UMAP_RNA', pt.size = 0.3, group.by = "trust4", cols=c("black","green","white","red"))

DimPlot(object = adata, reduction = 'UMAP_RNA', group.by = "rna_clusters4", label = TRUE)
DimPlot(object = adata, reduction = 'UMAP_RNA', group.by = "ct", label = TRUE)

DimPlot(object = adata, reduction = 'UMAP_RNA', group.by = "pred.hpca", label = TRUE)
DimPlot(object = adata, reduction = 'PCA_RNA', group.by = "pred.hpca", label = TRUE)  

DimPlot(object = adata, reduction = 'UMAP_RNA', group.by = "pred.monaco", label = TRUE)  
DimPlot(object = adata, reduction = 'PCA_RNA', group.by = "pred.monaco", label = TRUE)  



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
################### Comparison of ATAC and RNA clusters ########################  by UMAP
################################################################################

#todo we got the tool to compare 2 umaps! 

DimPlot(object = adata, reduction = "UMAP_ATAC", group.by = 'Phase')+ ggtitle('scATAC seq')


p1 <- DimPlot(object = adata, reduction = 'UMAP_RNA',  group.by = "rna_clusters", label = TRUE)+ ggtitle('rna clusters on scRNA umap')
p2 <- DimPlot(object = adata, reduction = "UMAP_ATAC", group.by = 'atac_clusters', label = TRUE)+ ggtitle('atac clusters on scATAC umap')
p3 <- DimPlot(object = adata, reduction = "UMAP_ATAC", group.by = 'rna_clusters', label = TRUE)  + ggtitle('rna clusters on ATAC umap')
p4 <- DimPlot(object = adata, reduction = 'UMAP_RNA',  group.by = "atac_clusters", label = TRUE) + ggtitle("atac clusters on rna umap")

p5 <- DimPlot(object = adata, reduction = 'UMAP_RNA',  group.by = "Phase", label = TRUE)+ ggtitle("phase on rna umap")
p6 <- DimPlot(object = adata, reduction = 'UMAP_ATAC', group.by = "Phase", label = TRUE)+ ggtitle("phase on atac umap")

p1/p2|p3/p4
p1/p2|p3/p4|p5/p6

############################

p1 <- DimPlot(object = adata, reduction = 'UMAP_RNA',  group.by = "rna_clusters3", label = TRUE)+ ggtitle('rna clusters on scRNA umap')
p2 <- DimPlot(object = adata, reduction = "UMAP_ATAC", group.by = 'atac_clusters3', label = TRUE)+ ggtitle('atac clusters on scATAC umap')
p3 <- DimPlot(object = adata, reduction = "UMAP_ATAC", group.by = 'rna_clusters3', label = TRUE)  + ggtitle('rna clusters on ATAC umap')
p4 <- DimPlot(object = adata, reduction = 'UMAP_RNA',  group.by = "atac_clusters3", label = TRUE) + ggtitle("atac clusters on rna umap")

p5 <- DimPlot(object = adata, reduction = 'UMAP_RNA',  group.by = "Phase", label = TRUE)+ ggtitle("phase on rna umap")
p6 <- DimPlot(object = adata, reduction = 'UMAP_ATAC', group.by = "Phase", label = TRUE)+ ggtitle("phase on atac umap")

p1/p2|p3/p4
p1/p2|p3/p4|p5/p6






################################################################################
################### Plot genes 2D ##############################################
################################################################################

DefaultAssay(adata) <- "RNA"

#### For paper -- supplementary
666
if(TRUE){
  DefaultAssay(adata) <- "RNA"
  p1 <- FeaturePlot(object = adata, reduction = 'UMAP_RNA', features = c("AICDA")) #+ ggtitle("RNA")
  p2 <- FeaturePlot(object = adata, reduction = 'UMAP_ATAC', features = c("AICDA"))#+ ggtitle("ATAC")
  p3 <- FeaturePlot(object = adata, reduction = 'UMAP_RNA', features = c("rank_norm_telo")) #+ ggtitle("RNA")
  p4 <- FeaturePlot(object = adata, reduction = 'UMAP_ATAC', features = c("rank_norm_telo"))#+ ggtitle("ATAC")
  p5 <- DimPlot(object = adata, reduction = 'UMAP_RNA',  group.by = "Phase", label = TRUE)+ ggtitle("Phase | RNA")
  p6 <- DimPlot(object = adata, reduction = 'UMAP_ATAC', group.by = "Phase", label = TRUE)+ ggtitle("Phase | ATAC")
  p1/p2|p3/p4|p5/p6
} else {
  DefaultAssay(adata) <- "RNA"
  #bdata <- adata[,str_starts(adata$ct,"GC")]
  #bdata <- adata[,adata$rna_clusters2 %in% c("r2","r3")]  #misses a small tip
  #DimPlot(object = adata, reduction = 'UMAP_RNA', group.by = "rna_clusters2", label = TRUE)
  
  bdata <- adata[,FetchData(object = adata, vars = names(adata@reductions[["UMAP_RNA"]]))[,1] > 8]
  
  p1 <- FeaturePlot(object = bdata, reduction = 'UMAP_RNA', features = c("AICDA")) #+ ggtitle("RNA")
  p3 <- FeaturePlot(object = bdata, reduction = 'UMAP_RNA', features = c("rank_norm_telo")) #+ ggtitle("RNA")
  p5 <- DimPlot(object = bdata, reduction = 'UMAP_RNA',  group.by = "Phase", label = TRUE)+ ggtitle("phase | RNA")
  
  bdata <- adata[,FetchData(object = adata, vars = names(adata@reductions[["UMAP_ATAC"]]))[,2] < -1]
  
  p2 <- FeaturePlot(object = bdata, reduction = 'UMAP_ATAC', features = c("AICDA"))#+ ggtitle("ATAC")
  p4 <- FeaturePlot(object = bdata, reduction = 'UMAP_ATAC', features = c("rank_norm_telo"))#+ ggtitle("ATAC")
  p6 <- DimPlot(object = bdata, reduction = 'UMAP_ATAC', group.by = "Phase", label = TRUE)+ ggtitle("phase | ATAC")
  p1/p2|p3/p4|p5/p6
  
}



showgenes <- c(
  "IGHD",  #immature, then entry to DZ. PC2 is really a IGHD axis
  
  "FOXP1",  #toward DZ/immature B. not so specific. so up in the most mature only
  
  "PDE4D",  #only subset is low in nonGC. what are those?
  
  "ATXN1",  #mature nonGC only. DNA binding. linked to EAE https://pubmed.ncbi.nlm.nih.gov/?term=ATXN1+b+cell&size=200
  "MAST4",  #links GC/nonGC. poorly studied  https://pubmed.ncbi.nlm.nih.gov/?term=mast4&size=200
  
  "MAML3",  #immature B, then GC
  
  "RAPGEF5", #beta catenin. https://pubmed.ncbi.nlm.nih.gov/?term=RAPGEF5&size=200
  
  "TOX", #all GC, then very specifically mature B
  
   ## Cell cycle cascade
  "TOP2A",  #DZ specifically
  "MKI67",  #DZ specifically
  "HMGB2",  #less DZ specific
  "SMC4",   #DZ mainly, some nonGC
  "CCND2",  #nonGZ only!
  
  "IGHM", #naive and DZ
  
  
  "CD83",   #up in naive and DZ
  "IRF4",   #up mainly in LZ and nonGC
  "NFKB1",  #up everywhere
  "FCRL4",  #up in mature only
  "ITGAX",   #up in mature only
  "SOX5"   #up in mature and most GC -- sits together in ATAC
)
FeaturePlot(adata, features = showgenes, reduction = "PCA_RNA") 
FeaturePlot(adata, features = showgenes, reduction = "SVD_ATAC", dims = use_atac_dims) 




showgenes <- c(
  "TNFRSF13B", #

  "ZEB2",  #up in all nonGC, some DZ -- 
  
  #These don't group well. another axis? r7
  "BCL11B",
  "THEMIS",
  "LEF1",
  "IL7R",
  
  #Plasmablast
  "PRDM1",
  "XBP1",
  "MZB1",
  "IGHG4",  

  #Certain subset of nonGC
  "IFI44L",
  "HERC5",
  "STAT1",  #Mainly in nonGC
  "CD55"  #MAYBE elevated in nonGC; and LZ?
  
  #Clustering thinks these are extensions of IFI44L + , into nonGC. it stretches it
  #"MX2",
  #"ISG15",
)
FeaturePlot(adata, features = showgenes, reduction = "PCA_RNA") 
FeaturePlot(adata, features = showgenes, reduction = "SVD_ATAC", dims = use_atac_dims) 


### From r2
showgenes <- c(
  "MAST4",    #microtubule-associated serine/threonine protein kinases; MAST2 regulates IL12 in macrophages
  "ADAMTS6",
  "PTPRJ",
  "COL4A3",
  "ARHGAP24",
  "GALNTL6",  #should be a new cluster on its own TODO
  "CCSER1",
  "BCL2"  # All nonGC
)
FeaturePlot(adata, features = showgenes, reduction = "PCA_RNA") 
FeaturePlot(adata, features = showgenes, reduction = "SVD_ATAC", dims = use_atac_dims) 






showgenes <- c(
  "CD83",
  "CXCR4",
  "CCND2"
)
FeaturePlot(adata, features = showgenes, reduction = "UMAP_RNA") 
FeaturePlot(adata, features = showgenes, reduction = "UMAP_ATAC") 
FeaturePlot(adata, features = showgenes, reduction = "PCA_RNA") 
FeaturePlot(adata, features = showgenes, reduction = "SVD_ATAC", dims = use_atac_dims) 




#Markers from https://www.frontiersin.org/files/Articles/602539/fimmu-12-602539-HTML/image_m/fimmu-12-602539-g001.jpg
showgenes <- c(
  #DN4
  "HOPX",
  "PDE4D",
  
  #DN2
  "EMP3",
  "ZEB2",
  
  #DN1
  "JCHAIN",
  
  #C-mem2 (poor markers)
  "MT-ATP8",
  "CRIP1",
  #C-mem1 similar to C-mem2
  
  #Naive
  "FCER2",
  "PLPP5",
  
  #Trans
  "SOX4",
  
  
  #### from elsewhere. B mem marker.
  "CD27"
)
FeaturePlot(adata, features = showgenes, reduction = "UMAP_RNA") 
FeaturePlot(adata, features = showgenes, reduction = "UMAP_ATAC") 
FeaturePlot(adata, features = showgenes, reduction = "PCA_RNA") 
FeaturePlot(adata, features = showgenes, reduction = "SVD_ATAC", dims = use_atac_dims) 








#Markers for ATAC clustering, different from rnaseq
showgenes <- c(
  
  #a3
  "COL19A1",
  "FCER2",
  "CD69",   #messy
  "ZBTB16",
  "FOS",  #bit messy. check
  
  #a4 -- definitely a subtype
  "IL4R",   
  "FCRL1",
  
  #a5 
  "DAAM1",
  "PELI1",
  "AFF2",
  "IGLC3",
  "MYBL1",
  #a5 - suggesting link to GC; toward r6 (memB)
  "TOX",
  "SOX5",

  #a6 enriched
  "CCND3",
  
  #a6 - def a subtype
  "ZBTB16",
  "GAB1",
  "SKAP1",
  "CD200",
  "FCER2"
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



#DNMT3A-lo

plot_umap_groupby(adata,"UMAP_RNA_3D","Phase", colors=c("red","green","blue"))
plot_umap_groupby(adata,"UMAP_ATAC_3D","Phase")

plot_umap_groupby(adata,"UMAP_RNA_3D","rna_clusters","Phase")
plot_umap_groupby(adata,"UMAP_RNA_3D","rna_clusters2","Phase")
plot_umap_groupby(adata,"UMAP_RNA_3D","rna_clusters3","Phase")  #brings out exit from r0 to plasma more clearly. GC split 3-way
plot_umap_groupby(adata,"UMAP_RNA_3D","rna_clusters4","Phase")
plot_umap_groupby(adata,"UMAP_RNA_3D","ct","Phase")
plot_umap_groupby(adata,"UMAP_RNA_3D","Phase","ct")

plot_umap_groupby(adata,"UMAP_ATAC_3D","rna_clusters","Phase")
plot_umap_groupby(adata,"UMAP_ATAC_3D","rna_clusters3","Phase")
plot_umap_groupby(adata,"UMAP_ATAC_3D","rna_clusters4","Phase")
plot_umap_groupby(adata,"UMAP_ATAC_3D","atac_clusters3","ct")
plot_umap_groupby(adata,"UMAP_ATAC_3D","ct","atac_clusters3")
#MLLT3 enrich in ... atac_a6. goes into FCRL4hi=atac3_a10

plot_umap_gene(adata,"UMAP_RNA_3D","FCRL5","rna_clusters")
plot_umap_gene(adata,"UMAP_RNA_3D","rank_norm_telo","rna_clusters")

plot_umap_gene(adata,"UMAP_ATAC_3D","AICDA","rna_clusters") 
plot_umap_gene(adata,"UMAP_ATAC_3D","rank_norm_telo","rna_clusters") 

#B1 vs B2. B1 has no memory. cannot say anything
plot_umap_gene(adata,"UMAP_RNA_3D","MS4A1","ct") #up in LZ tip, FCRL4 hi
plot_umap_gene(adata,"UMAP_RNA_3D","CD27","ct") #up in LZ early, mem
plot_umap_gene(adata,"UMAP_RNA_3D","CD5","ct") #



#Easiest to use for filtering. note that the T cell cluster has some B too! monocyte cluster got some B too, all types
plot_umap_groupby(adata,"UMAP_RNA_3D","pred.hpca","Phase")  



plot_umap_gene(adata,"UMAP_RNA_3D","IGHD","rna_clusters")  #immature
plot_umap_gene(adata,"UMAP_RNA_3D","FOXP1","rna_clusters")  #end point is the jump to GC

plot_umap_gene(adata,"UMAP_RNA_3D","PDE4D","rna_clusters")  #weird pattern. start, end, mid...
plot_umap_gene(adata,"UMAP_RNA_3D","ATXN1","rna_clusters")  #toward switched, plasmablast
plot_umap_gene(adata,"UMAP_RNA_3D","MAST4","rna_clusters")  #from jump to r2
plot_umap_gene(adata,"UMAP_RNA_3D","CD27","rna_clusters")   #Bmem marker. in r1,r2, GC, plasma
plot_umap_gene(adata,"UMAP_RNA_3D","CD83","rna_clusters")   #r5, then jumping to tip of r3 which is DZ. also up in 
plot_umap_gene(adata,"UMAP_RNA_3D","CCND2","rna_clusters")  #r5, CC checkpoint into GC
plot_umap_gene(adata,"UMAP_RNA_3D","","rna_clusters")
plot_umap_gene(adata,"UMAP_RNA_3D","AICDA","rna_clusters")
plot_umap_gene(adata,"UMAP_RNA_3D","rank_norm_telo","ct")
plot_umap_gene(adata,"UMAP_RNA_3D","FCRL4","rna_clusters")
plot_umap_gene(adata,"UMAP_RNA_3D","XBP1","rna_clusters4")
plot_umap_gene(adata,"UMAP_RNA_3D","JUN","rna_clusters4")


plot_umap_gene(adata,"UMAP_RNA_3D","FCER2","rna_clusters")   ### Most clear separator, naive vs switched
plot_umap_gene(adata,"UMAP_RNA_3D","FOS","rna_clusters")    #all over. maybe more naive
plot_umap_gene(adata,"UMAP_RNA_3D","IL4R","rna_clusters")   #distinct in part of r0. this cluster need be split
plot_umap_gene(adata,"UMAP_RNA_3D","FCRL1","rna_clusters")  #gradient from naive. gone in plasma
plot_umap_gene(adata,"UMAP_RNA_3D","CCND3","rna_clusters")  #gradient from naive. gone in plasma !!
plot_umap_gene(adata,"UMAP_RNA_3D","CD200","rna_clusters")  #r0 subset. need to split
plot_umap_gene(adata,"UMAP_RNA_3D","ZBTB16","rna_clusters") #r0 subset. need to split
plot_umap_gene(adata,"UMAP_RNA_3D","","rna_clusters")
plot_umap_gene(adata,"UMAP_RNA_3D","","rna_clusters")
plot_umap_gene(adata,"UMAP_RNA_3D","","rna_clusters")



######## Comparison to king
plot_umap_gene(adata,"UMAP_RNA_3D","FCRL4","rna_clusters")  #specific. r6 and r4_7. specific. they defines as part of mbc, I agree

#"activated". not sure what cluster this is. spread
plot_umap_gene(adata,"UMAP_RNA_3D","CD69","rna_clusters") 
plot_umap_gene(adata,"UMAP_RNA_3D","JUN","rna_clusters") 
plot_umap_gene(adata,"UMAP_RNA_3D","EGR1","rna_clusters") 
plot_umap_gene(adata,"UMAP_RNA_3D","DUSP2","rna_clusters")  #mainly r5 actually
plot_umap_gene(adata,"UMAP_RNA_3D","IL6","rna_clusters4")   #maybe more r2, r4_1

#naive - we agree
plot_umap_gene(adata,"UMAP_RNA_3D","FCER2","rna_clusters") 

#pre-gC. ccnd2. we agree
plot_umap_gene(adata,"UMAP_RNA_3D","CCND2","rna_clusters") 
plot_umap_gene(adata,"UMAP_RNA_3D","PSME2","rna_clusters") 

#LZ GC. yes, but small tip!
plot_umap_gene(adata,"UMAP_RNA_3D","EBI3","rna_clusters")  #also in pre-gC!

#DZ GC. 
plot_umap_gene(adata,"UMAP_RNA_3D","SUGCT","rna_clusters")  #interesting, this is the mid point of our GC
plot_umap_gene(adata,"UMAP_RNA_3D","ISG20","rna_clusters")  #interesting, this is the mid point of our GC
plot_umap_gene(adata,"UMAP_RNA_3D","AICDA","rna_clusters")  #this is mid point, but also separate cells!!!! rather specific for them. but then a bit in "cycling" linked to GC

#FCRL3+
plot_umap_gene(adata,"UMAP_RNA_3D","FCRL3","rna_clusters4")  #separate cluster in king. not for us. might correspond to r2 or r4_1; but their is close to GC

#Pre-PB
plot_umap_gene(adata,"UMAP_RNA_3D","FRZB","rna_clusters4")   #barely any
plot_umap_gene(adata,"UMAP_RNA_3D","HOPX","rna_clusters4")   #in LZ tip
plot_umap_gene(adata,"UMAP_RNA_3D","BTNL9","rna_clusters4")  #in LZ tip
plot_umap_gene(adata,"UMAP_RNA_3D","FGFR1","rna_clusters4")  #broad in GC

#cycling. matches pretty well
plot_umap_gene(adata,"UMAP_RNA_3D","HMGB2","rna_clusters4")
plot_umap_gene(adata,"UMAP_RNA_3D","TUBA1B","rna_clusters4")
plot_umap_gene(adata,"UMAP_RNA_3D","MKI67","rna_clusters4")
plot_umap_gene(adata,"UMAP_RNA_3D","","rna_clusters4")



plot_umap_gene(adata,"UMAP_RNA_3D","CCSER1","rna_clusters")




########################################################################
################ comparison to King annotation #########################
########################################################################

#We disagree a bit on the preB (their is T); a single myeloid
plot_umap_groupby(adata,"UMAP_RNA_3D","king_lineage")
plot_umap_groupby(adata,"UMAP_ATAC_3D","king_lineage")
table(adata$king_lineage)

#It is quite a mess! overall agreement
plot_umap_groupby(adata,"UMAP_RNA_3D","king_celltype")
plot_umap_groupby(adata,"UMAP_ATAC_3D","king_celltype")
table(adata$king_celltype)

#interesting but not sure how useful. bit all over
plot_umap_groupby(adata,"UMAP_RNA_3D","king_MBC")
plot_umap_groupby(adata,"UMAP_ATAC_3D","king_MBC")

#not sure how to use
plot_umap_groupby(adata,"UMAP_RNA_3D","king_subset")
plot_umap_groupby(adata,"UMAP_ATAC_3D","king_subset")


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
##################
#r0 naive B
#r1 mostly switched, some exhausted
#r2 mostly nonswitched
#r3 exhausted, confused annotation
#r4 exhausted
#r5 mostly naive
#r6 mostly exhausted
#r7 total spread (proB)
#r8 total spread (plasmablast)
#r9 mostly naive

count_overlaps_clusterings(adata$rna_clusters, adata$pred.hpca)
##################
#r0 naive
#r1 memory, some naive
#r2 naive, then memory
#r3 GC
#r4 GC; CXCR4+
#r5 naive
#r6 naive/memory 50:50
#r7 naive/memory
#r8 plasmablast
#r9 naive


count_overlaps_clusterings(adata$rna_clusters3, adata$king_celltype)
# Quite useless!


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


ct_assignment <- c(
  "r4_12:Plasmablast",  #no doubt high in xbp1
  
  "r4_8:GC DZ S",
  "r4_10:GC DZ G2M",
  "r4_11:GC LZ1",
  "r4_9:GC LZ2",
  "r4_6:GC LZ tip",
  "r4_4:Follicular to GC",

  "r4_7:Follicular FCRL4_hi",
  "r4_1:Follicular mem TOX-hi FOXP1-lo",  #was towardmem PDE4D_maybe_from_LZ
  "r4_2:Follicular mem FOXP1-lo TOX-hi",  #was MAST4
  "r4_14:Follicular interferon",  #quite spread in memory. and naive?
  "r4_5:Follicular naive MLLT3", #MLLT3 is in polymerase super elongation complex. furthest from GC; closest to FOXP1-lo.
  "r4_0:Follicular naive DNMT3A-lo", #mixes in ATACseq; DNMT3A differs. bit up in ZBTB16
  "r4_3:Follicular naive DNMT3A-hi", #DNMT3A does DNA methylation. these up in IKZF2
  
  "r4_13:Follicular intermediate"  #up in IGHG1 oddly. somewhat toward plasma   ---- was naive odd
)
adata$ct <- assignCelltypes(adata$rna_clusters4, ct_assignment)

Idents(adata) <- adata$ct
DotPlot(adata, features = unique(features), assay="RNA") + RotatedAxis()


if(FALSE){
  ## King markers
  features <- c(
    "FCRL4",  #specific. r6 and r4_7. specific. they defines as part of mbc, I agree
    
    #"activated". not sure what cluster this is. spread
    "CD69", 
    "JUN",
    "EGR1", 
    "DUSP2",
    "IL6",
    
    #naive - we agree
    "FCER2",
    
    #pre-gC. ccnd2. we agree
    "CCND2",
    "PSME2",
    
    #LZ GC. yes, but small tip!
    "EBI3",  #also in pre-gC!
    
    #DZ GC. 
    "SUGCT",  #interesting, this is the mid point of our GC
    "ISG20",  #interesting, this is the mid point of our GC
    "AICDA",  #this is mid point, but also separate cells!!!! rather specific for them. but then a bit in "cycling" linked to GC
    
    #FCRL3+
    "FCRL3",  #separate cluster in king. not for us. might correspond to r2 or r4_1; but their is close to GC
    
    #Pre-PB
    "FRZB",   #barely any
    "HOPX",   #in LZ tip
    "BTNL9",  #in LZ tip
    "FGFR1",  #broad in GC
    
    #cycling. matches pretty well
    "HMGB2",
    "TUBA1B",
    "MKI67"
  )
  Idents(adata) <- adata$rna_clusters4
  DotPlot(adata, features = features, assay="RNA") + RotatedAxis()
}







#towardmem PDE4D & MAST4: also TNFRSF13B, ANK3, 
#PDE4D: TBXAS1, SSPN, bit more TOX, bit less FOXP1 & ZEB2
#note nonGC FCRL4 hi, seems like a more developed version of above

#AKAP6: also SCN3A, SYN3, PCDH9
#
#all naivish: ZBTB16


## try split nonGC
features <- c(
  "AKAP6"
)
Idents(adata) <- adata$ct
DotPlot(adata, features = unique(features), assay="RNA") + RotatedAxis()


# split naive
features <- c(
  "SH3RF3",
  "ACTB",
  "CCSER1",
  "YBX3",
  "IL7",
  "IKZF2",
  "DNMT3A",
  "IGLC1",
  "IGLC2",
  "IGLC3",
  "BACH2",
  "SNTG2",
  "MLLT3",
  "FCRL5",
  "IL4R"
)
Idents(adata) <- adata$ct
DotPlot(adata, features = unique(features), assay="RNA") + RotatedAxis()


features <- c(

  "Sox5",
  
  #DZ G2M: 
  "POU2F3",
  "POU3F4", 
  "TCF4", 
  "SNAI1",
  "ZEB1",
  
  #DZ S
  "BATF::JUN", 
  "BATF3",
  "Smad2::Smad3", 
  "BACH2",
  
  #LZ tip
  "Klf1", 
  "NR3C2", 
  "Ar",
  "CTCF",
  "GATA3",
  "FOXP1",  #also for memory
  
  #FCR4+
  "YY1",
  "ATF7",
  "JUNB",

  #nonGC to GC
  "NFKB2",
  
  #PB
  "XBP1",
  "KLF15", #also near CTCF/GATA3
  "E2F6",  #also near CTCF/GATA3
  "ETV6",
  
  #Interferon B cells
  "IRF9",
  "IRF4",
  
  #Memory
  #"TOX",  #HMG-box  -- but no specific motif
  #"FOXP1",
  
  #GC LZ2, but largely nonGC
  "SPIC",
  "IKZF1",
  
  
  "ZBTB7A",
  "NFKB2",
  "GABPA",
  "EWSR1-FLI1",
  "ELF5",
  "STAT1::STAT2",
  "ETV1"
)
Idents(adata) <- adata$ct
DotPlot(adata[,adata$ct!="Plasmablast"], features = unique(names2motifs[features,]$id), assay="chromvar") + RotatedAxis()
#DotPlot(adata, features = unique(names2motifs[features,]$id), assay="chromvar") + RotatedAxis()
#DotPlot(adata, features = motif2names$id[motif2names$name %in% features], assay="chromvar") + RotatedAxis()

if(FALSE){
  adata_for_motif_display <- adata[,adata$ct!="Plasmablast"]
  rownames(adata_for_motif_display@assays$chromvar@data) <- motif2names[rownames(adata_for_motif_display),]$name
  Idents(adata_for_motif_display) <- adata_for_motif_display$ct
  DefaultAssay(adata_for_motif_display) <- "chromvar"
  
  DotPlot(adata_for_motif_display, features = unique(features), assay="chromvar") + RotatedAxis()
666
    
  remove(adata_for_motif_display)  
} else {
  
  #Might not make sense with a dotplot for motif activity, since it is a dense matrix!
  
  listct <- levels(adata$ct)
  meanact <- matrix(nrow=length(listct), ncol=nrow(adata@assays$chromvar@data))
  for(i in 1:length(listct)){
    meanact[i,] <- rowMeans(adata@assays$chromvar@data[,adata$ct==listct[i]])
  }
  rownames(meanact) <- listct
  colnames(meanact) <- rownames(adata@assays$chromvar@data)
  
  meanact <- meanact[rownames(meanact)!="Plasmablast",]
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
}



#rownames(adata_for_motif_display) 
#aw.data, data, and scale.data slots directly, 

#interferon
motif2names[c("MA0653.1","MA1419.1"),] #IRF9 and IRF4

motif2names[c("MA0095.2","MA0834.1"),] #FCRL4: YY1, maybe ATF7
motif2names[c("MA0627.2","MA0830.2","MA1558.1","MA0789.1"),] #DZ G2M: POU2F3/POU3F4, TCF4, SNAI1, 
motif2names[c("MA0493.1","MA1622.1","MA0462.2","MA0835.2"),] #DZ S: Klf1, Smad2::Smad3, BATF::JUN, BATF3
motif2names[c("MA0727.1","MA0778.1","MA0007.3"),] #LZ tip: especially Nc3C2, NFKB2. lesser Ar
motif2names[c("MA0844.1",""),] #PB: XBP1,    but many more
#CTCF up in PB, and DZ S. bit in G2M



########################################################################
################ Chromvar analysis #####################################
########################################################################

DefaultAssay(adata) <- "peaks"

## Compute motif activity; this is deviations in counts for each motif across peaks
## as compared to total count. it does not care /which/ peaks the motif is in, just the total sum
##https://www.nature.com/articles/nmeth.4401
adata<- RunChromVAR(
  object = adata,
  genome = BSgenome.Hsapiens.UCSC.hg38, assay = 'peaks'
)



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



#Compare specifically LZ and DZ   
DefaultAssay(adata) <- 'chromvar'
Idents(adata) <- adata$rna_clusters3
cluster_motifs_gc <- FindAllMarkers(adata[,adata$rna_clusters3 %in% c("r3_5", "r3_8", "r3_9", "r3_10", "r3_7")], only.pos = F, test.use = 'wilcox')
cluster_motifs_gc %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) -> topmotifs_gc_rna3
topmotifs_gc_rna3$name <- motif2names[topmotifs_gc_rna3$gene,]$name






########################################################################
################ Motif compared to axes ################################   JASPAR motifs
########################################################################


compare_gene_with_axis <- function(adata, the_axis){
  use_genes <- adata@assays$RNA@var.features
  allcor <- NULL
  for(i in use_genes){
    print(i)
    allcor <- c(allcor, cor(as.double(adata@assays$RNA[i,]),the_axis,method = "spearman"))
  }
  df <- data.frame(
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


###### nTelo axis  -- within GC
axis_telo <- adata$norm_telo
use_sub <- adata$rna_clusters2 %in% c("r2","r3")
motifcor_telo <- compare_motif_with_axis(adata[,use_sub], axis_telo[use_sub])

genecor_telo <- compare_gene_with_axis(adata[,use_sub], axis_telo[use_sub])
rbind(head(genecor_telo,n=20),tail(genecor_telo, n=20))
# cor       gene
# 657 -0.2320954      KCNQ5 ## channel. little known
# 428 -0.2228774  LINC01572
# 643 -0.2103337       GMDS
# 109 -0.2066552     DIAPH3 # very anticorr in GC. CC gene for sure https://elifesciences.org/articles/61974
# 292 -0.2045190       CDK1 ## cell cycle!
# 174 -0.2037290     KIF18B

# neg: CENPI, HISH1H1B

# 89   0.1602015 AC104170.1
# 76   0.1635084      BACH2 ## TF!
# 401  0.1797303   IFNG-AS1
# 816  0.1959182    RNF144B
# 578  0.2191169      SUGCT
# 150  0.2204567      KANK1

genecor_telo[genecor_telo$gene=="AICDA",] 
#10% corr

rbind(head(motifcor_telo),tail(motifcor_telo))
#most neg: NFKB2, EBF1, JUNB, JUND, RELA, FOSL2::JUNB 
#most pos: TCF3, FOXK2, TCF*, FOXP3 !!! 

plot_umap_motif("UMAP_RNA_3D", "BACH2","rna_clusters3", dorank=TRUE) #up in LZ and cycling!
plot_umap_motif("UMAP_RNA_3D", "EBF1","rna_clusters3", dorank=TRUE) #up in the tip in LZ primarily!

plot_umap_gene(adata,"UMAP_RNA_3D","BACH2","rna_clusters3") #gene exp up in dark zone
plot_umap_gene(adata,"UMAP_RNA_3D","AICDA","rna_clusters3") #gene exp up in dark zone
plot_umap_gene(adata,"UMAP_RNA_3D","CXCR4","rna_clusters3") #
plot_umap_gene(adata,"UMAP_RNA_3D","EBF1","rna_clusters3") # up in the LZ tip
plot_umap_gene(adata,"UMAP_RNA_3D","CDK1","rna_clusters3") #
plot_umap_gene(adata,"UMAP_RNA_3D","rank_norm_telo","rna_clusters3") #gene exp up in dark zone
plot_umap_groupby(adata,"UMAP_RNA_3D","Phase","rna_clusters3", colors=c("red","green","blue")) #


plot_umap_gene(adata,"UMAP_RNA_3D","MZB1","rna_clusters3") #up where CDK1 up. wut?


plot_umap_motif("UMAP_RNA_3D", "CTCF","rna_clusters3", dorank=TRUE) #up in LZ and cycling!


DefaultAssay(adata) <- 'RNA'
FeaturePlot(adata, features = c("rank_norm_telo",genecor_telo$gene[1:5],tail(genecor_telo$gene,n=5)), reduction = "PCA_RNA")   #positive corr genes are in mature and vice versa
FeaturePlot(adata, features = c("rank_norm_telo",genecor_telo$gene[1:5],tail(genecor_telo$gene,n=5)), reduction = "UMAP_RNA") 

head(motifcor_maturity[order(motifcor_maturity$cor, decreasing = TRUE),])

# DefaultAssay(adata) <- 'chromvar'
# FeaturePlot(adata, features = c("rank_norm_telo",motif2names[motif2names$name=="Sox5",]$id), reduction = "PCA_RNA") 
# FeaturePlot(adata, features = c("rank_norm_telo",motif2names[motif2names$name=="Sox5",]$id), reduction = "UMAP_RNA") 

#FeaturePlot(adata, features = c("rank_norm_telo",motif2names[str_starts(str_to_lower(motif2names$name),"sox"),]$id), reduction = "PCA_RNA") 
#FeaturePlot(adata, features = c("rank_norm_telo",motif2names[str_starts(str_to_lower(motif2names$name),"sox"),]$id), reduction = "UMAP_RNA") 


adata$tempact <- rank(as.double(adata@assays$chromvar[motif2names[motif2names$name=="Sox5",]$id,]))
FeaturePlot(adata, features = c("rank_norm_telo","tempact"), reduction = "PCA_RNA") 
FeaturePlot(adata, features = c("rank_norm_telo","tempact"), reduction = "UMAP_RNA") 
plot(adata$rank_norm_telo,adata$tempact,pch=19,cex=0.5)



###### PC2 is maturity axis; or SVD4?
use_sub <- adata$rna_clusters2=="r0"
axis_maturity <- FetchData(object = adata, vars = names(adata@reductions[["PCA_RNA"]]))[,2]
motifcor_maturity <- compare_motif_with_axis(adata[,use_sub], axis_maturity[use_sub])

genecor_maturity <- compare_gene_with_axis(adata[,use_sub], axis_maturity[use_sub])


###### SVD2 is GC axis, in ATAC?
use_sub <- adata$rna_clusters2 %in% c("r0","r1")
axis_gc <- FetchData(object = adata, vars = names(adata@reductions[["SVD_ATAC"]]))[,2]
motifcor_gc <- compare_motif_with_axis(adata[,use_sub], axis_gc[use_sub])

genecor_gc <- compare_gene_with_axis(adata[,use_sub], axis_gc[use_sub])
#Or use limma voom

########################

#Maturity - Mainly pulls GC motifs for some reason...
FeaturePlot(adata, features = rev(motifcor_maturity$motif)[1:5], reduction = "PCA_RNA") 

##### GC - This works well
FeaturePlot(adata, features = c(motifcor_gc$motif[1:5],tail(motifcor_gc$motif,n=5)), reduction = "PCA_RNA") 
head(motifcor_gc,n=10)
tail(motifcor_gc,n=20)
#ELF5/BCL6 is the most nonGC. 
#POU2F3 and family are the most gC; less certain about LIN54; BATF::JUN
FeaturePlot(adata, features = c("MA0627.2",), reduction = "PCA_RNA") 



#Telo - I think this works well!
use_sub <- adata$rna_clusters2=="r1"
FeaturePlot(adata[,use_sub], features = c("rank_norm_telo",rev(motifcor_telo$motif)[1:5]), reduction = "PCA_RNA") #maybe positive cor with FOXL1/FOXO/TCF4 -- same as in large motif analysis!
FeaturePlot(adata[,use_sub], features = c("rank_norm_telo",rev(motifcor_telo$motif)[1:5]), reduction = "UMAP_RNA") 

FeaturePlot(adata[,use_sub], features = c("rank_norm_telo",(motifcor_telo$motif)[1:5]), reduction = "PCA_RNA") ### this might work. NFKB2, EBF1, JUNB anticorrelate
FeaturePlot(adata[,use_sub], features = c("rank_norm_telo",(motifcor_telo$motif)[1:5]), reduction = "UMAP_RNA") 


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

#Clearly up in cycling
plot_umap_motif("UMAP_RNA_3D", "SP1","rna_clusters", dorank=TRUE)
plot_umap_motif("UMAP_RNA_3D", "Klf12","rna_clusters", dorank=TRUE) 

#Up in tip of cycling, and LZ, and exit from nonGC. this suggests that REL is related to the first entry into GC.
#also likely possibly to split dividing cells somehow?
plot_umap_motif("UMAP_RNA_3D", "REL","rna_clusters", dorank=TRUE)
plot_umap_motif("UMAP_RNA_3D", "NFKB2","rna_clusters", dorank=TRUE)
plot_umap_motif("UMAP_RNA_3D", "JUN","rna_clusters", dorank=TRUE)

#Complex pattern but suggesting related to entry into GC/LZ; but more diffuse
plot_umap_motif("UMAP_RNA_3D", "ZKSCAN5","rna_clusters", dorank=TRUE)


plot_umap_motif("UMAP_RNA_3D", "NR3C2","rna_clusters", dorank=TRUE)
plot_umap_motif("UMAP_RNA_3D", "Ar","rna_clusters", dorank=TRUE)


#The most extreme tip of cycling, and plasmablast. This is an exit TF
#More in switched B cells. later division?
plot_umap_motif("UMAP_RNA_3D", "IRF9","rna_clusters", dorank=TRUE)

#DZ and plasmablast
plot_umap_motif("UMAP_RNA_3D", "FOXI1","rna_clusters", dorank=TRUE)

#Complex pattern; plasma, but broad
plot_umap_motif("UMAP_RNA_3D", "STAT1::STAT2","rna_clusters", dorank=TRUE)


#Very confusing pattern. broad, more immature, but also in LZ and cycling
plot_umap_motif("UMAP_RNA_3D", "IKZF1","rna_clusters", dorank=TRUE)

#nonGC
plot_umap_motif("UMAP_RNA_3D", "BCL6","rna_clusters", dorank=TRUE)


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

DimPlot(object = motif_adata,  label = TRUE) + ggtitle('Motif activity similarity')


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
  "PPARG",
  "RORA",
  "STAT1"
  # "SREBF2",
  # "NFAT",
  # "C3",
  # "RORC",
  # "MAF",
  # "Stat5a",
  # "BCL6",
  # "YY2",
  # "SOX9",
  # "GATA3",
  # #HOXA2 cluster 2 unclear. find!
  # "CTCF",
  # "Arid3a",
  # "FOXP3",
  # 
  # "JUN",
  # "TBX21",
  # "FLI1",
  # "IRF4",
  # "BHLHE40",
  # "CLOCK",
  # "HIF1A", #and ARNT::HIF1A
  # "KLF3",
  # "BATF",
  # "BACH2",
  # "ZEB1",
  # "Rbpjl"
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
############ Motif to gene expression correlation ######################
########################################################################

correlateMotifActivityExpression <- function(adata){
  
  ### Load mapping motif-gene
  map_jaspar_uniprot <- read.csv('/home/mahogny/jupyter/tcellpaper/jaspar_uniprot.csv', header = T, sep = ',')[c(3,5,6)]
  colnames(map_jaspar_uniprot) <- c("jaspar_baseid","jasparname","uniprot")
  map_jaspar_uniprot <- map_jaspar_uniprot[map_jaspar_uniprot$uniprot!="",]
  map_uniprot_symbol <- read.csv('/home/mahogny/jupyter/tcellpaper/uniprot_symbol_ensid.csv', header = T, sep = '\t')
  colnames(map_uniprot_symbol) <- c('ensid','uniprot','symbol')
  
  ## Only keep genes which are in adata
  #DefaultAssay(adata) <- 'RNA'
  map_uniprot_symbol <- map_uniprot_symbol[map_uniprot_symbol$symbol %in% rownames(adata@assays$RNA),]
  merged_jaspar <- unique(merge(map_uniprot_symbol, map_jaspar_uniprot))

  merged_jaspar <- merge(
    merged_jaspar,
    data.frame(
      jaspar_baseid=str_split_fixed(rownames(adata@assays$chromvar),"\\.",2)[,1],
      jaspar_id=rownames(adata@assays$chromvar))
  )

  #Filter: keep the first item having symbol:jaspar_id  -- seems duplicated now
  merged_jaspar <- merged_jaspar[isUnique(sprintf("%s--%s",merged_jaspar$symbol, merged_jaspar$jaspar_id)),]
  merged_jaspar
  
  ## Only keep highly variable genes ... ??
  
  merged_jaspar$pval <- NA
  merged_jaspar$cor <- NA
  for(i in 1:nrow(merged_jaspar)){
    ctres <- cor.test(
      as.double(adata@assays$RNA[merged_jaspar$symbol[i],]),
      as.double(adata@assays$chromvar[merged_jaspar$jaspar_id[i],]))
    merged_jaspar$cor[i] <- ctres$estimate
    merged_jaspar$pval[i] <- ctres$p.value  
  }
  merged_jaspar$log_pval <- -log10(merged_jaspar$pval)
  merged_jaspar <- merged_jaspar[!is.na(merged_jaspar$cor),]
  merged_jaspar
}

corMAE.adata <- correlateMotifActivityExpression(adata)
corMAE.adata <- corMAE.adata[order(corMAE.adata$pval),]
corMAE.adata <- corMAE.adata[!is.infinite(corMAE.adata$log_pval),]

ggplot(corMAE.adata, aes(cor, log_pval, label=symbol)) + geom_point() + geom_label()

666
#possible figure for paper


########################################################################
################ Zoom in on GC #########################################
########################################################################

bdata <- adata[,str_starts(adata$ct,"GC")]

#note: must do this BEFORE changing IDs in bdata
corMAE.bdata <- correlateMotifActivityExpression(bdata)
corMAE.bdata <- corMAE.bdata[order(corMAE.bdata$pval),]
ggplot(corMAE.bdata, aes(cor, log_pval, label=symbol)) + geom_point() + geom_label()
unique(corMAE.bdata$symbol)[1:20]
#"BACH2"  "NFKB2"  "TCF3"   "NFKB1"  "BATF"   "FOXP1"  "BATF3"  "EBF1"   "E2F7"   "E2F8"   
#"ELK3"   "TFDP1"  "KLF12"  "MEF2B"  "RELB"   "CUX1"   "E2F1"   "POU3F1" "SPIB"   "JUND"  
top_features_GC_rna <- c(
  "BACH2","NFKB2","TCF3","NFKB1","BATF","FOXP1","EBF1","E2F7","E2F8",
  "ELK3","KLF12","MEF2B","RELB","CUX1","E2F1"
  #removed BATF3, JUND, SPIB, TFDP1
)

unique(corMAE.bdata$jasparname)[1:20]
# [1] "BACH2"      "NFKB2"      "TCF3"       "NFKB1"      "BATF::JUN"  "BATF"       "FOXP1"      "BATF3"      "EBF1"       "TAL1::TCF3"
# [11] "E2F7"       "E2F8"       "ELK3"       "TFDP1"      "KLF12"      "MEF2B"      "RELB"       "CUX1"       "E2F1"       "POU3F1"    
top_features_GC_motif <- c(
  "BACH2","NFKB2","TCF3","NFKB1","BATF::JUN","BATF","FOXP1","BATF3","EBF1","TAL1::TCF3",
  "E2F7","E2F8","ELK3","TFDP1","KLF12","MEF2B","RELB","CUX1","E2F1","POU3F1"
)


############### Redo dimensional reduction
bdata <- RunPCA(bdata,reduction.name = 'PCA_RNA',ndims.print = c(1,3,5), verbose = T)
bdata <- RunUMAP(object = bdata, dims = 1:30, reduction = "PCA_RNA", reduction.name="UMAP_RNA", spread=0.5)
bdata <- RunUMAP(bdata, dims = 1:30,reduction = 'PCA_RNA',n.components = 3L, reduction.name ='UMAP_RNA_3D' )
bdata <- RunSVD(bdata, verbose = T, reduction.name = 'SVD_ATAC')
use_dim_2 <- 2:30 
bdata <- RunUMAP(object = bdata, reduction = 'SVD_ATAC', dims = use_dim_2, reduction.name = 'UMAP_ATAC')
#bdata <- RunUMAP(object = bdata, reduction = 'SVD_ATAC', dims = use_dim_2, n.components = 3L, reduction.name = 'UMAP_ATAC_3D')



## New clustering at higher resolution
# DefaultAssay(bdata) <- "RNA"
# bdata <- FindNeighbors(bdata, dims = 1:30,reduction = 'PCA_RNA')
# bdata <- FindClusters(bdata, resolution = 1.3)
# bdata$rna_clusters5 <- factor(sprintf("r5_%s",bdata$seurat_clusters))
# cluster_markers_rna5 <- FindAllMarkers(bdata, only.pos = F, test.use = 'wilcox')
# cluster_markers_rna5 %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) -> top10_rna5
# 
# DefaultAssay(bdata) <- 'chromvar'
# Idents(bdata) <- bdata$rna_clusters5
# cluster_motifs_rna5 <- FindAllMarkers(bdata, only.pos = F, test.use = 'wilcox')
# cluster_motifs_rna5 %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) -> topmotifs_rna5
# topmotifs_rna5$name <- motif2names[topmotifs_rna5$gene,]$name
# DimPlot(object = bdata, reduction = 'UMAP_RNA',  group.by = c("ct","rna_clusters5"), label = TRUE)


## New clustering at higher resolution -- rna6
DefaultAssay(bdata) <- "RNA"
bdata <- FindNeighbors(bdata, dims = 1:30,reduction = 'PCA_RNA')
bdata <- FindClusters(bdata, resolution = 1.5)
bdata$rna_clusters6 <- factor(sprintf("r6_%s",bdata$seurat_clusters))
cluster_markers_rna6 <- FindAllMarkers(bdata, only.pos = F, test.use = 'wilcox')
cluster_markers_rna6 %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) -> top10_rna6

DefaultAssay(bdata) <- 'chromvar'
Idents(bdata) <- bdata$rna_clusters6
cluster_motifs_rna6 <- FindAllMarkers(bdata, only.pos = F, test.use = 'wilcox')
cluster_motifs_rna6 %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) -> topmotifs_rna6
topmotifs_rna6$name <- motif2names[topmotifs_rna6$gene,]$name

#DimPlot(object = bdata, reduction = 'UMAP_RNA',  group.by = c("ct","rna_clusters6","Phase"), label = TRUE)



####### Update cell type assignment for the GC, also in the adata object
ct_assignment_gc <- c(
  "r6_1:GC DZ S",       #RRM, POLQ?
  "r6_5:GC DZ S",       #DTL  MCM4,  other MCM*   repair program
  "r6_6:GC DZ G2M",     #PLK1   not so specific?
  "r6_8:GC DZ G2M",     #AURKB  etc  not specific    merge with 6?
  "r6_10:GC DZ IGHD S/G2M",    #hmm. CD44
  "r6_2:GC DZ AICDA KMT2A-hi",   #AFF3  FOXP1  CDK13 
  "r6_9:GC DZ AICDA BMP7-hi",  #also RIMS2+  ZNF385B+  #### BMP7 FOXP1     PAX5-low! 
  "r6_0:GC LZ MAML3+",   #MAML3, DNER
  "r6_3:GC LZ intermediate",  #CD22 FOXN3   !!   check in 3d!! 
  "r6_4:GC LZ exiting IGHM+",   #BANK1   IGHM  CD83+
  "r6_7:GC LZ exiting EBI3+NFKB1+"   #EBI3   NFKB1   CD83
)
bdata$ct <- assignCelltypes(bdata$rna_clusters6, ct_assignment_gc)

## Transfer celltype to adata as well
newct <- data.frame(
  row.names = colnames(bdata), 
  ct = bdata$ct
)
newct <- as.character(newct[colnames(adata),,drop=FALSE]$ct)
newct[is.na(newct)] <- as.character(adata$ct)[is.na(newct)]
adata$ct <- newct





if(FALSE){
  plot_umap_groupby(bdata,"UMAP_RNA_3D","Phase", colors=c("red","green","blue"))
  plot_umap_gene(bdata,"UMAP_RNA_3D","AICDA","rna_clusters6")
  plot_umap_gene(bdata,"UMAP_RNA_3D","CD83","rna_clusters6")
  plot_umap_gene(bdata,"UMAP_RNA_3D","MAML3","rna_clusters6")  ### ideal!
  plot_umap_gene(bdata,"UMAP_RNA_3D","rank_norm_telo","rna_clusters6")
  plot_umap_groupby(bdata,"UMAP_RNA_3D","rna_clusters6")
  plot_umap_gene(bdata,"UMAP_RNA_3D","KANK1","rna_clusters6")  #top correlating!
  plot_umap_gene(bdata,"UMAP_RNA_3D","CDK13","rna_clusters6") 
  plot_umap_gene(bdata,"UMAP_RNA_3D","POLQ","rna_clusters6")   #not always overlapping?
  plot_umap_gene(bdata,"UMAP_RNA_3D","CD44","rna_clusters6")   #for IGHM!
  plot_umap_gene(bdata,"UMAP_RNA_3D","EBI3","rna_clusters6")   #specific for r6_7!
  plot_umap_gene(bdata,"UMAP_RNA_3D","IGHD","rna_clusters6") 
  plot_umap_gene(bdata,"UMAP_RNA_3D","IGHA_sum","rna_clusters6") 
  plot_umap_gene(bdata,"UMAP_RNA_3D","IGHG_sum","rna_clusters6") 
  plot_umap_gene(bdata,"UMAP_RNA_3D","IGHE","rna_clusters6") 
  plot_umap_gene(bdata,"UMAP_RNA_3D","BACH2","rna_clusters6") 
  plot_umap_gene(bdata,"UMAP_RNA_3D","NFKB1","rna_clusters6") 
  plot_umap_gene(bdata,"UMAP_RNA_3D","FOXN3","rna_clusters6")   #down near AICDA... especially around BMP7?
  plot_umap_gene(bdata,"UMAP_RNA_3D","E2F1","rna_clusters6") #E2F1 and E2F8 probably represent slightly different cells. but for another paper 
  plot_umap_gene(bdata,"UMAP_RNA_3D","E2F8","rna_clusters6") 
  plot_umap_gene(bdata,"UMAP_RNA_3D","CD22","rna_clusters6") 
  plot_umap_gene(bdata,"UMAP_RNA_3D","KMT2A","rna_clusters6") 
  plot_umap_gene(bdata,"UMAP_RNA_3D","MAML3","rna_clusters6") 
  plot_umap_gene(bdata,"UMAP_RNA_3D","PAX5","rna_clusters6") 
  plot_umap_gene(bdata,"UMAP_RNA_3D","MEF2B","rna_clusters6") #mainly in late IGH*! 
  plot_umap_gene(bdata,"UMAP_RNA_3D","CD83","rna_clusters6") #maybe
  plot_umap_gene(bdata,"UMAP_RNA_3D","SNX29","rna_clusters6") 
  plot_umap_gene(bdata,"UMAP_RNA_3D","HLA-DQA1","rna_clusters6") #maybe
  plot_umap_gene(bdata,"UMAP_RNA_3D","FOXN3","rna_clusters6") 
  plot_umap_gene(bdata,"UMAP_RNA_3D","PIP4K2A","rna_clusters6") 
  plot_umap_gene(bdata,"UMAP_RNA_3D","IGHA_sum","rna_clusters6") 
  plot_umap_gene(bdata,"UMAP_RNA_3D","IGHG_sum","rna_clusters6") 
  plot_umap_gene(bdata,"UMAP_RNA_3D","CUX1","rna_clusters6") 
  plot_umap_gene(bdata,"UMAP_RNA_3D","KLF12","rna_clusters6") 
  plot_umap_gene(bdata,"UMAP_RNA_3D","PAX5","rna_clusters6")   #down in the mid, cluster r6_9
  plot_umap_gene(bdata,"UMAP_RNA_3D","CD22","rna_clusters6") #Looks like it forms a gap in LZ. Inhibits BCR. down in DZ, esp AICDA
  plot_umap_gene(bdata,"UMAP_RNA_3D","BMP7","rna_clusters6")   #up in the mid, cluster r6_9.... but not all of them!!
  plot_umap_gene(bdata,"UMAP_RNA_3D","FOXP1","rna_clusters6")   #up in the mid, cluster r6_9  and r6_2
}





#################################### Good markers for GC
features <- c(
  "MAML3","DNER",  #
  "RRM2","POLQ",
  "AFF3","FOXP1","CDK13",  #agree
  "CD22","FOXN3",
  "BANK1","IGHM",
  "DTL","MCM4",
  "PLK1",
  "EBI3","NFKB1","CD83",
  "AURKB",
  "RIMS2","BMP7",

  "IGHA_sum",
  "IGHG_sum",
  #"IGHA1", "IGHA2",
  #"IGHG1", "IGHG2", "IGHG3", "IGHG4",
  "IGHE"  #Furthest away
  
  #"IGHC"
)
DefaultAssay(bdata) <- "RNA"
FeaturePlot(object = bdata, reduction = 'UMAP_RNA', features = c("rank_norm_telo",unique(features)),ncol=5, pt.size = 0.5) #+ ggtitle("RNA")

features <- c(
  "nCount_RNA",
  "LTB",
  "CD52", #up in all but tip and part of CC. biomarker for cancer https://www.frontiersin.org/articles/10.3389/fgene.2020.578002/full
  "BANK1", #up in LZ2 = early BC. hmm.
  "IL7", #limited to LZ2 = early BC
  "CD22", #mature B. prevents overactivation https://en.wikipedia.org/wiki/CD22 inhibitory for BCR
  "SNX29", #up in LZ. we do not understand it
  "PIP4K2A",  #additional marker for late-stage GC B
  "FOXN3", #possible marker for late-stage GC B. possibly doing DNA damage CC checkpointing. expression inhibits prolif https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5190042/
  "EBI3", #up in tip. part of IL27 and IL35. IL35 blocks T cells!!! NEW! (also Bregs wtf) 
  "SGPP2", #up in tip; up in inflammatory responses. we know nothing
  "DUSP2", #up in tip; upstream of STAT3 https://www.nature.com/articles/ni.3278  ; STAT3 important for B cell activity and GC maintenance https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4875824/ 
  "CD83", #maturation marker, in tip. rename based on it. oddly, also up in LZ2=IGHM/IGHD
  "BCL2A1" #anti-apoptotic https://www.nature.com/articles/cdd2011158
)
DefaultAssay(bdata) <- "RNA"
FeaturePlot(object = bdata, reduction = 'UMAP_RNA', features = c("rank_norm_telo",unique(features))) #+ ggtitle("RNA")
FeaturePlot(object = bdata, reduction = 'UMAP_ATAC', features = c("rank_norm_telo",unique(features))) #+ ggtitle("RNA")






features <- c(
  #### vs telomere
  #From gene exp
  "BACH2",
  
  #Most negative
  "NFKB2",
  "EBF1",
  "JUNB",
  "RELA",
  
  #Positive
  "FOXP3",
  "FOXK2",
  "TCF3"
)
DefaultAssay(bdata) <- "chromvar"
FeaturePlot(object = bdata, reduction = 'UMAP_RNA', features = c("rank_norm_telo",unique(features))) #+ ggtitle("RNA")
FeaturePlot(object = bdata, reduction = 'UMAP_ATAC', features = c("rank_norm_telo",unique(features))) #+ ggtitle("RNA")


features <- c(
  "KCNQ5", #no doubt
  "GMDS",  #no doubt
  "DIAPH3",
  "CDK1",
  "MME",
  
  "FCER2",
  
  #signature weak/strong here https://www.nature.com/articles/ni.3460/figures/2
  "BATF",
  
  "BACH2",  #strong correlation
  "RNF144B", #perfect
  "SUGCT", #perfect
  "KANK1" #perfect
)
DefaultAssay(bdata) <- "RNA"
FeaturePlot(object = bdata, reduction = 'UMAP_RNA', features = c("rank_norm_telo",unique(features))) #+ ggtitle("RNA")
FeaturePlot(object = bdata, reduction = 'UMAP_ATAC', features = c("rank_norm_telo",unique(features))) #+ ggtitle("RNA")



#### CORRELATING GENES
#### Final supplemental figure? genes to be sold
features <- c(
  "BACH2",  #strong correlation
  "KANK1",
  "BATF"
)
DefaultAssay(bdata) <- "RNA"
p1 <- FeaturePlot(object = bdata, reduction = 'UMAP_RNA', features = c("rank_norm_telo",unique(features)), ncol=4) #+ ggtitle("RNA")
p2 <- FeaturePlot(object = bdata, reduction = 'UMAP_ATAC', features = c("rank_norm_telo",unique(features)), ncol=4) #+ ggtitle("RNA")
p1/p2




#### SHM genes -- mismatch DNA repair (MMR)
list_hypermutation <- c(
  #mismatch DNA repair (MMR)
  "MLH1",
  "MLH3",
  "MSH5",
  "RPA1",
  "RPA2",
  "ABL1",
  #"RNASEH2A",
  "MCM9",
  "RPA3",
  "RNASEH2C",
  "SETD2",
  "TDG",
  #"TREX1",
  "EXO1",
  "PMS1",
  "PMS2",
  "TP73",
  "XPC",
  "MSH3",
  
  #If expressed, not near AICDA
  "RNASEH2B",
  "LIG1",
  "MSH2",
  "MCM8",
  "POLD3",
  "PCNA"
)
DefaultAssay(bdata) <- "RNA"
FeaturePlot(object = bdata, reduction = 'UMAP_RNA', features = c("rank_norm_telo",unique(list_hypermutation))) #+ ggtitle("RNA")
FeaturePlot(object = bdata, reduction = 'UMAP_ATAC', features = c("rank_norm_telo",unique(list_hypermutation))) #+ ggtitle("RNA")




#### IGH genes
#https://en.wikipedia.org/wiki/IGH@
features <- c(
  
  #Symbols for variable (V) immunoglobulin gene segments start with IGHV and include two or three numbers separated by dashes. Examples:
  #IGHV1-2, IGHV1-3, …, IGHV1-69-2, IGHV2-5, …, IGHV7-4-1
  "IGHV_sum",
  
  #Symbols for diversity (D) immunoglobulin gene segments start with IGHD and include two numbers separated by dashes. Examples:
  #IGHD1-1, IGHD1-7, …, IGHD7-27
  "IGHD_sum",
  
  #"IGHJ1", "IGHJ2", "IGHJ3", "IGHJ4", "IGHJ5", "IGHJ6",  #Nothing to see
  
  #Symbols for constant region (C) immunoglobulin genes: 
  
  "IGHA1", "IGHA2",
  "IGHG1", "IGHG2", "IGHG3", "IGHG4",
  "IGHE",  #Furthest away
  "IGHD", #up in LZ2
  "IGHM", #up in LZ2. note that LZ2 is a bit down in nTA; how to explain?
  
  "AICDA"
)
DefaultAssay(bdata) <- "RNA"
FeaturePlot(object = bdata, reduction = 'UMAP_RNA', features = c("rank_norm_telo",unique(features))) #+ ggtitle("RNA")
FeaturePlot(object = bdata, reduction = 'UMAP_ATAC', features = c("rank_norm_telo",unique(features))) #+ ggtitle("RNA")





########################################################################
################ Do ATACseq pileups over GC ############################
########################################################################

library(Rsamtools) #tabix file
library(stringi)
library(fastmatch)


GetReadsInRegion <- function(
    cellmap,
    region,
    tabix.file,
    cells = NULL,
    verbose = TRUE,
    ...
) {
  file.to.object <- names(x = cellmap)
  names(x = file.to.object) <- cellmap
  
  if (verbose) {
    message("Extracting reads in requested region")
  }
  if (!is(object = region, class2 = "GRanges")) {
    region <- StringToGRanges(regions = region, ...)
  }
  # remove regions that aren't in the fragment file
  common.seqlevels <- intersect(
    x = seqlevels(x = region),
    y = seqnamesTabix(file = tabix.file)
  )
  region <- keepSeqlevels(
    x = region,
    value = common.seqlevels,
    pruning.mode = "coarse"
  )
  reads <- scanTabix(file = tabix.file, param = region)
  reads <- TabixOutputToDataFrame(reads = reads)
  reads <- reads[
    fmatch(x = reads$cell, table = cellmap, nomatch = 0L) > 0,
  ]
  # convert cell names to match names in object
  reads$cell <- file.to.object[reads$cell]
  if (!is.null(x = cells)) {
    reads <- reads[reads$cell %in% cells, ]
  }
  if (nrow(reads) == 0) {
    return(reads)
  }
  reads$length <- reads$end - reads$start
  return(reads)
}


SetIfNull <- function(x, y) {
  if (is.null(x = x)) {
    return(y)
  } else {
    return(x)
  }
}

MultiGetReadsInRegion <- function(
    object,
    region,
    fragment.list = NULL,
    assay = NULL,
    ...
) {
  if (inherits(x = object, what = "Seurat")) {
    # pull the assay
    assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
    object <- object[[assay]]
  }
  fragment.list <- SetIfNull(
    x = fragment.list,
    y = Fragments(object = object)
  )
  if (length(x = fragment.list) == 0) {
    # no fragments set
    stop("No fragment files found")
  }
  res <- data.frame()
  for (i in seq_along(along.with = fragment.list)) {
    tbx.path <- GetFragmentData(object = fragment.list[[i]], slot = "path")
    cellmap <- GetFragmentData(object = fragment.list[[i]], slot = "cells")
    tabix.file <- TabixFile(file = tbx.path)
    open(con = tabix.file)
    reads <- GetReadsInRegion(
      cellmap = cellmap,
      region = region,
      tabix.file = tabix.file,
      ...
    )
    res <- rbind(res, reads)
    close(con = tabix.file)
  }
  return(res)
}

TabixOutputToDataFrame <- function(reads, record.ident = TRUE) {
  if (record.ident) {
    nrep <- elementNROWS(x = reads)
  }
  reads <- unlist(x = reads, use.names = FALSE)
  if (length(x = reads) == 0) {
    df <- data.frame(
      "chr" = "",
      "start" = "",
      "end" = "",
      "cell" = "",
      "count" = ""
    )
    df <- df[-1, ]
    return(df)
  }
  reads <- stri_split_fixed(str = reads, pattern = "\t")
  n <- length(x = reads[[1]])
  unlisted <- unlist(x = reads)
  e1 <- unlisted[n * (seq_along(along.with = reads)) - (n - 1)]
  e2 <- as.numeric(x = unlisted[n * (seq_along(along.with = reads)) - (n - 2)])
  e3 <- as.numeric(x = unlisted[n * (seq_along(along.with = reads)) - (n - 3)])
  e4 <- unlisted[n * (seq_along(along.with = reads)) - (n - 4)]
  e5 <- as.numeric(x = unlisted[n * (seq_along(along.with = reads)) - (n - 5)])
  df <- data.frame(
    "chr" = e1,
    "start" = e2,
    "end" = e3,
    "cell" = e4,
    "count" = e5,
    stringsAsFactors = FALSE,
    check.rows = FALSE,
    check.names = FALSE
  )
  if (record.ident) {
    df$ident <- rep(x = seq_along(along.with = nrep), nrep)
  }
  return(df)
}


makePileUp <- function(region, ct){
  sumCountAtac <- sum(adata$nCount_ATAC[adata$ct == ct])
  cells <- colnames(adata)[adata$ct == ct]
  reads <- MultiGetReadsInRegion(
    object = adata, assay = "ATAC", 
    region = region, cells = cells, verbose = FALSE)
  pileupStart <- min(reads$start)
  pileupEnd <- max(reads$end)
  pileup <- rep(0,pileupEnd-pileupStart+1)
  for(i in 1:nrow(reads)){
    ind <- (reads$start[i]-pileupStart+1):(reads$end[i]-pileupStart+1)
    pileup[ind] <- pileup[ind]+1  
  }
  d <- data.frame(pos=pileupStart:pileupEnd, height=pileup)
#  plot((pileupStart:pileupEnd)[seq(1, length(pileup), 100)], pileup[seq(1, length(pileup), 100)],type="l")
  d$cs <- cumsum(d$height)
  d <- d[seq(1, nrow(d), 100),]
  d$sumh <- c(d$cs[2:nrow(d)]-d$cs[1:(nrow(d)-1)],0)
  plot(d$pos, d$sumh/sumCountAtac, type="l", ylim=c(0,1e-3))
}


#1e-03
#old thereg <- "chr12-113234845-116010625"
thereg <- str_replace_all("chr14-105,571,527-107,043,718",",","")

makePileUp(thereg, "GC DZ S")
makePileUp(thereg, "GC DZ G2M")

makePileUp(thereg, "GC LZ1")     #
makePileUp(thereg, "GC LZ2")
makePileUp(thereg, "GC LZ tip")  #one peak up
makePileUp(thereg, "Follicular to GC")

makePileUp(thereg, "Plasmablast")
makePileUp(thereg, "Follicular mem TOXhi FOXP1lo")
makePileUp(thereg, "Follicular mem FOXP1lo TOXhi")


unique(adata$ct)
length(pileup)

# "r4_8:GC DZ S",
# "r4_10:GC DZ G2M",
# "r4_11:GC LZ1",
# "r4_9:GC LZ2",
# "r4_6:GC LZ tip",
# "r4_4:Follicular to GC",






########################################################################
################ Do RNAseq pileups over GC #############################
########################################################################

library(Rsamtools) 
library(GenomicAlignments)


makePileUpRna <- function(region, ct){
  sumCountRna <- sum(adata$nCount_RNA[adata$ct == ct])
  
  #region <- "chr12:113234845-116010625"
  #ct <- "GC DZ S"
  
  reads <- NULL
  for(curlib in c("lib1","lib2","lib3","lib4")){
    bv <- BamViews(sprintf("/corgi/cellbuster/bigb/cellranger-arc/%s/outs/gex_possorted_bam.bam",curlib),
                   bamRanges=as(region, "GRanges"))
    foo <- readGAlignments(bv, param = ScanBamParam(tag=c("CB")))
    gex_cb <- str_split_fixed(foo$gex_possorted_bam.bam@elementMetadata$CB,"-",2)[,1]
    use_cb <- str_split_fixed(colnames(adata)[adata$lib == curlib & adata$ct == ct],"-",2)[,1]
    reads <- rbind(
      reads,
      unique(data.frame(
        start=start(foo$gex_possorted_bam.bam)[gex_cb %in% use_cb],
        end=end(foo$gex_possorted_bam.bam)[gex_cb %in% use_cb]
      ))
    )
  }
  
  pileupStart <- min(reads$start)
  pileupEnd <- max(reads$end)
  pileup <- rep(0,pileupEnd-pileupStart+1)
  for(i in 1:nrow(reads)){
    ind <- (reads$start[i]-pileupStart+1):(reads$end[i]-pileupStart+1)
    pileup[ind] <- pileup[ind]+1  
  }
  d <- data.frame(pos=pileupStart:pileupEnd, height=pileup)
  #  plot((pileupStart:pileupEnd)[seq(1, length(pileup), 100)], pileup[seq(1, length(pileup), 100)],type="l")
  d$cs <- cumsum(d$height)
  d <- d[seq(1, nrow(d), 100),]
  d$sumh <- c(d$cs[2:nrow(d)]-d$cs[1:(nrow(d)-1)],0)
  d$sumh <- d$sumh/sumCountRna
  d$ct <- ct
  plot(d$pos, d$sumh, type="l", ylim=c(0,0.004), main=ct)
  d
}


thereg <- str_replace_all("chr14:105,571,527-107,043,718",",","")
thereg

# makePileUpRna(thereg, "GC DZ S")  #0.0015
# makePileUpRna(thereg, "GC DZ G2M")  #0.0020
# makePileUpRna(thereg, "GC LZ1")     #0.004
# makePileUpRna(thereg, "GC LZ2")   #0.004
# makePileUpRna(thereg, "GC LZ tip")  #0.0030
# 
# makePileUpRna(thereg, "Follicular to GC") #0.0025
# 
# makePileUpRna(thereg, "Plasmablast") #0.004
# makePileUpRna(thereg, "Follicular mem TOXhi FOXP1lo") #0.003
# makePileUpRna(thereg, "Follicular mem FOXP1lo TOXhi") #0.003



allpileup <- rbind(
  makePileUpRna(thereg, "GC DZ S"),  #0.0015
  makePileUpRna(thereg, "GC DZ G2M"),  #0.0020
  makePileUpRna(thereg, "GC LZ1"),     #0.004
  makePileUpRna(thereg, "GC LZ2"),   #0.004
  makePileUpRna(thereg, "GC LZ tip"),
  makePileUpRna(thereg, "Plasmablast"))  #0.0030
  
  #makePileUpRna(thereg, "Follicular to GC") #0.0025

ggplot(allpileup[allpileup$ct!="Plasmablast",], aes(pos, sumh))+geom_line()+ facet_wrap(~ct, nrow=1)
ggplot(allpileup, aes(pos, sumh))+geom_line()+ facet_wrap(~ct, ncol=1)
ggplot(allpileup[allpileup$ct!="Plasmablast",], aes(pos, sumh, color=ct))+geom_line()#+ylim(0,0.005)
#114e6
#chr12:114000000-114000000










########################################################################
################ Final plots for the paper #############################
########################################################################

#Simplified cell types
ct_no_gc <- as.character(adata$ct)
ct_no_gc[str_starts(ct_no_gc,"GC DZ")] <- "GC DZ"
ct_no_gc[str_starts(ct_no_gc,"GC LZ")] <- "GC LZ"
adata$ct_nogc <- factor(ct_no_gc, levels = sort(unique(ct_no_gc)))

#Ensure alphabetic order
adata$ct <- factor(as.character(adata$ct),levels = sort(unique(as.character(adata$ct))))



#################### all overview
DefaultAssay(bdata) <- "RNA"
#p1 <- FeaturePlot(raster=TRUE, pt.size = 3, object = adata, reduction = 'UMAP_RNA', features = c("AICDA")) #+ ggtitle("RNA")
p3 <- FeaturePlot(raster=TRUE, pt.size = 3, object = adata, reduction = 'UMAP_RNA', features = c("rank_norm_telo")) + ggtitle("rank(nTA)")
p5 <- DimPlot(raster=TRUE, pt.size = 3, object = adata, reduction = 'UMAP_RNA',  group.by = "Phase", label = TRUE)+ ggtitle("Phase")
p7 <- DimPlot(raster=TRUE, pt.size = 3, object = adata, reduction = 'UMAP_RNA',  group.by = "ct_nogc")+ ggtitle("Cell type")#, label = TRUE)

#p2 <- FeaturePlot(object = adata, reduction = 'UMAP_ATAC', features = c("AICDA"))#+ ggtitle("ATAC")
p4 <- FeaturePlot(raster=TRUE, pt.size = 3, object = adata, reduction = 'UMAP_ATAC', features = c("rank_norm_telo")) + ggtitle("rank(nTA)")
p6 <- DimPlot(raster=TRUE, pt.size = 3, object = adata, reduction = 'UMAP_ATAC', group.by = "Phase", label = TRUE)+ ggtitle("Phase")
p8 <- DimPlot(raster=TRUE, pt.size = 3, object = adata, reduction = 'UMAP_ATAC', group.by = "ct_nogc")+ ggtitle("Cell type")#, label = TRUE)
#p1/p2|
p3/p4|p5/p6|p7/p8
ggsave("~/bcell_all_plot1.pdf", width = 14)


#################### GC overview
DefaultAssay(bdata) <- "RNA"
p1 <- FeaturePlot(raster=TRUE, pt.size = 3, object = bdata, reduction = 'UMAP_RNA', features = c("AICDA")) #+ ggtitle("RNA")
p3 <- FeaturePlot(raster=TRUE, pt.size = 3, object = bdata, reduction = 'UMAP_RNA', features = c("rank_norm_telo")) + ggtitle("rank(nTA)")
p5 <- DimPlot(raster=TRUE, pt.size = 3, object = bdata, reduction = 'UMAP_RNA',  group.by = "Phase", label = TRUE)+ ggtitle("Phase")
p7 <- DimPlot(raster=TRUE, pt.size = 3, object = bdata, reduction = 'UMAP_RNA',  group.by = "ct")+ ggtitle("Cell type")#, label = TRUE)

p2 <- FeaturePlot(raster=TRUE, pt.size = 3, object = bdata, reduction = 'UMAP_ATAC', features = c("AICDA"))#+ ggtitle("ATAC")
p4 <- FeaturePlot(raster=TRUE, pt.size = 3, object = bdata, reduction = 'UMAP_ATAC', features = c("rank_norm_telo")) + ggtitle("rank(nTA)")
p6 <- DimPlot(raster=TRUE, pt.size = 3, object = bdata, reduction = 'UMAP_ATAC', group.by = "Phase", label = TRUE)+ ggtitle("Phase")
p8 <- DimPlot(raster=TRUE, pt.size = 3, object = bdata, reduction = 'UMAP_ATAC', group.by = "ct")+ ggtitle("Cell type")#, label = TRUE)
p1/p2|p3/p4|p5/p6|p7/p8
ggsave("~/bcell_gc_plot1.pdf", width = 16)




#head(corMAE.adata)
corMAE.adata$symbol2 <- corMAE.adata$symbol
corMAE.adata$symbol2[corMAE.adata$log_pval<10] <- NA
ggplot(corMAE.adata, aes(cor, log_pval, label=symbol2)) + geom_point() + geom_label() + 
  ylim(0,100) + xlim(-0.3,0.3) + xlab("Correlation") + ylab("-Log10 pval")
ggsave("~/bcell_tfcor_all.pdf", width = 4)

corMAE.bdata$symbol2 <- corMAE.bdata$symbol
corMAE.bdata$symbol2[corMAE.bdata$log_pval<4] <- NA
ggplot(corMAE.bdata, aes(cor, log_pval, label=symbol2)) + geom_point() + geom_label() +
  ylim(0,60) + xlim(-0.5,0.3) + xlab("Correlation") + ylab("-Log10 pval")
ggsave("~/bcell_tfcor_gc.pdf", width = 4)










############################# Marker genes at the highest level
features <- c(
  
  #Entry into GC
  "CCND2",
  
  #TFs all up in GC
  "TCL1A",
  "BCL6",
  "MEF2B",
  
  ##### GC subsets. to work out!
  #"MKI67", #DZ*
  #"TOP2A", #DZ*
  #"BRIP1", #DZ*
  #"DIAPH3", #DZ*
  
  "SMC4",
  "MME",
  "MAML3",
  "SOX5",
  "EZH2",
  "AFF2",
  "KLHL6",
  
  #GC DZ S
  "E2F8",
  "E2F1",
  
  #GC tip
  "EBF1",
  
  #GZ LZ1
  "BCL7A",
  "POU4F1",
  "BACH2",
  
  
  "CXCR4", #Should be gradient through CC in GC
  "CD83",  #Supposed LZ marker
  "AICDA",
  
  #Memory subsets. Need to discuss FOXP1 and TOX.
  "TOX",     #TOX up in GC; then later down. 
  "FOXP1",   #being FOXP1 is the normal thing. so FOXP1 down at the very end?
  "PDE4D",
  "PTPRJ",
  "MAST4",
  
  #Naive B. IgM but not IgD is downregulated on autoreactive B cells. 
  "IGHD",
  "IGHM",
  
  #All naive; spread out in RNA
  "ZBTB16",
  "CCND3",
  "FCER2",
  
  #naive GAB1
  "GAB1",
  "IKZF2",
  "DNMT3A",
  
  #Naive AKAP6
  "AKAP6",  # MTOCs https://pubmed.ncbi.nlm.nih.gov/33295871/
  "SNTG2",
  "MLLT3",
  "IL4R",
  
  #FCRL4 hi
  "FCRL4",
  "TNFRSF13B",
  "ANK3",
  "ZEB2",
  
  #all nonGC
  "CD69",  
  "JUN",
  "FCER2",
  
  #Plasmablast  
  "IRF4",
  "PRDM1",
  "XBP1",
  "MZB1",
  "IGHG1",
  
  #Memory B cell marker; down in naive. up also in GC
  "CD27",
  
  #Interferon B
  "IFI44L"
  
)
DefaultAssay(adata) <- "RNA"
Idents(adata) <- adata$ct_nogc
DotPlot(adata, features = unique(features), assay="RNA") + RotatedAxis()
ggsave("~/bcell_all_plot2.pdf", width = 18, height = 4)





############################# Marker genes at GC level
features <- c(
  "BMP7","ZNF385B",
  
  "POLQ",
  "UNG",
  
  "AICDA",
  "CXCR4",
  "KANK1",
  "CDK13",
  
  "FOXP1",
  "KLF12",
  "ELK3",
  "BACH2",
  "FOXO1",
  
  "NFKB1",
  "CD83",
  
  "CUX1",
  "MAML3",
  "CD22",
  "PAX5", #complex pattern
  
  "FOXP1",
  
  "IGHM",
  #"IGHD_sum",
  "IGHG_sum",
  
  "DNER"
)
DefaultAssay(bdata) <- "RNA"
Idents(bdata) <- bdata$ct
DotPlot(bdata, features = unique(features), assay="RNA") + RotatedAxis()
ggsave("~/bcell_gc_plot2.pdf", width = 10, height = 3.3)






#################### GC specific genes
DefaultAssay(bdata) <- "RNA"
p1 <- FeaturePlot(raster=TRUE, pt.size = 3, object = bdata, reduction = 'UMAP_RNA', features = c("MAML3")) #+ ggtitle("RNA")
p2 <- FeaturePlot(raster=TRUE, pt.size = 3, object = bdata, reduction = 'UMAP_RNA', features = c("BACH2"))#+ ggtitle("ATAC")
p3 <- FeaturePlot(raster=TRUE, pt.size = 3, object = bdata, reduction = 'UMAP_RNA', features = c("CD83")) #+ ggtitle("RNA")
p4 <- FeaturePlot(raster=TRUE, pt.size = 3, object = bdata, reduction = 'UMAP_RNA', features = c("CDK13"))#+ ggtitle("ATAC")
p5 <- FeaturePlot(raster=TRUE, pt.size = 3, object = bdata, reduction = 'UMAP_RNA', features = c("CXCR4")) #+ ggtitle("RNA")
p6 <- FeaturePlot(raster=TRUE, pt.size = 3, object = bdata, reduction = 'UMAP_RNA', features = c("FOXP1"))#+ ggtitle("ATAC")
p1/p2|p3/p4|p5/p6#|p7/p8
ggsave("~/bcell_gc_plot3.pdf", height = 6, width = 10)







dev.off()



####################### TF activities
features <- c(
  #"Sox5",
  #DZ G2M: 
  #"POU2F3",
  #"POU3F4", 
  "TCF3", 
  "SNAI1",
  "ZEB1",
  
  #DZ S
  #"BATF::JUN", 
  #"BATF3",
  #"Smad2::Smad3", 
  "BACH2",
  
  #LZ tip
  #"Klf1", 
  "NR3C2", 
  #"Ar",
  #"CTCF",
  "GATA3",
  "FOXP1",  #also for memory
  
  #FCR4+
  "YY1",
  "ATF7",
  "JUNB",
  
  #nonGC to GC
  "NFKB2",
  
  #PB
  "XBP1",
  #"KLF15", #also near CTCF/GATA3
  #"E2F6",  #also near CTCF/GATA3
  #"ETV6",
  
  #Interferon B cells
  #"IRF9",
  #"IRF4",
  "IRF1", #possibly more likely
  
  #Memory
  #"TOX",  #HMG-box  -- but no specific motif
  #"FOXP1",
  
  #GC LZ2, but largely nonGC
  "SPIC",
  "IKZF1",
  
  "ELK3",
  "TFDP1",
  "SPIB",
  
  
  
  #Correlates well with transcription levels; so probably what to present
  "BACH2",
  "BATF",
  "NFKB2",
  "TCF3",
  "FOXP1",
  "EBF1",
  "E2F7",
  "KLF12",
  "MEF2B",
  #"CUX1",
  "POU3F1",
  
  "JUND",
  "HES6",
  "JUNB",
  #"BACH1",
  "IRF1",
  "PAX5"
)


listct <- levels(adata$ct)
meanact <- matrix(nrow=length(listct), ncol=nrow(adata@assays$chromvar@data))
for(i in 1:length(listct)){
  meanact[i,] <- rowMeans(adata@assays$chromvar@data[,adata$ct==listct[i]])
}
rownames(meanact) <- listct
colnames(meanact) <- rownames(adata@assays$chromvar@data)

meanact <- meanact[rownames(meanact)!="Plasmablast",]
for(i in 1:ncol(meanact)){
  meanact[,i] <- meanact[,i]-min(meanact[,i])
  meanact[,i] <- meanact[,i]/max(meanact[,i])
}

m <- merge(reshape2::melt(meanact),data.frame(Var2=motif2names$id, tf=motif2names))
m <- m[m$tf.name %in% features,]
m$tf.name <- factor(m$tf.name, levels = unique(features))

ggplot(m, aes(tf.name, Var1, fill=value)) +   #    geom_tile() + 
  geom_point(aes(size=value, color=value)) + 
  scale_size(range = c(0, 5)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("~/bcell_tfactivity.pdf", width = 8, height = 4.5)



########################################################################
######################### export for VAE analysis ######################
########################################################################


if(FALSE){
  # if (!requireNamespace("remotes", quietly = TRUE)) {
  #   install.packages("remotes")
  # }
  # remotes::install_github("mojaveazure/seurat-disk")
  # 
  # library(SeuratDisk)
  # SaveH5Seurat(CreateSeuratObject(bdata[["RNA"]]), "/corgi/johanlehto/bdata_gc.h5", TRUE)
  # 
  
  write.csv(bdata@meta.data, "/corgi/johanlehto/bdata_gc.meta.csv")
  
}



########################################################################
######################### potential viewer #############################
########################################################################


data.frame(
  ct=adata$ct,
  as.double(adata@assays$RNA["XBP1",]))


axname <- names(adata@reductions[["UMAP_RNA_3D"]])
plot.data <- FetchData(object = adata, vars = c(axname,"XBP1","JUNB","IGHD","ct"))
colnames(plot.data)[1:3] <- c("x","y","z")#,"colorby","textby")
plot.data

write.csv(plot.data,"~/umap3d_bcell.csv")



########################################################################
######################### cellxgene viewer #############################
########################################################################


library(Seurat)
library(SeuratData)
library(SeuratDisk)


#new object with just key metadata

# SaveH5Seurat(adata, filename = "/corgi/websites/bcellatlas/allb.h5Seurat")
# Convert("/corgi/websites/bcellatlas/allb.h5Seurat", dest = "h5ad")

store_h5ad <- function(bdata, fname){
  stupidob <- CreateSeuratObject(bdata[[c("RNA")]])
  stupidob[["chromvar"]] <- bdata[[c("chromvar")]]
  stupidob@meta.data <- bdata@meta.data[,c("pred.hpca","pred.monaco","ct","Phase","donor_id","rawcnt_telo","rank_norm_telo")]
  stupidob$ct <- as.character(stupidob$ct)
  
  theax <- FetchData(object = bdata, vars = c(names(bdata@reductions[["UMAP_RNA"]])))
  colnames(theax) <- c("UMAP1","UMAP2")
  
  stupidob[["UMAP"]] <- CreateDimReducObject(embeddings = as.matrix(theax), key = "UMAP_", assay = "RNA")
  
  SaveH5Seurat(stupidob, filename = fname, overwrite = TRUE)
  Convert(fname, dest = "h5ad",overwrite = TRUE)
  
  remove(stupidob)
}

store_h5ad(bdata, "/corgi/websites/bcellatlas/gcb.h5Seurat")
store_h5ad(adata, "/corgi/websites/bcellatlas/allb.h5Seurat")  #TODO update annotation from bdata!!

