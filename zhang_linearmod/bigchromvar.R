

#devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())
library(stringr)
library(limma)
library(ArchR)
library(stringr)


############### read data ################################
combined_dataset <- readRDS("/corgi/henriksson/telomere/zhang2021_all/chromvar_proj.rds")
ArrowFiles <- combined_dataset$ArrowFiles
proj <- combined_dataset$proj
datadir <- combined_dataset$datadir
archr_genome <- combined_dataset$genome


set.seed(1)
addArchRThreads(threads = 16)
addArchRGenome(archr_genome)   #hg38 or mm10



# combined_dataset <- readRDS("chromvar_proj.rds")




################################################################
############# Prepare chromvar #################################
################################################################

print("============= activity scores ============")

# pseudobulk
proj@cellColData$Clusters <- proj@cellColData$zhang_ct
proj <- addGroupCoverages(ArchRProj = proj, groupBy = "Clusters")

# identify peaks, using cell groupings to aid this
pathToMacs2 <- findMacs2()
proj <- addReproduciblePeakSet(
    ArchRProj = proj, 
    groupBy = "Clusters", 
    pathToMacs2 = pathToMacs2
)

proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif")

if("Motif" %ni% names(proj@peakAnnotation)){
    proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif")
}
proj <- addPeakMatrix(ArchRProj = proj)
proj <- addBgdPeaks(proj)
proj <- addDeviationsMatrix(
  ArchRProj = proj, 
  peakAnnotation = "Motif",
  force = TRUE
)

#could subset to look at the most varying ones. hm. #especially after fitting a linear model
#plotVarDev <- getVarDeviations(proj, name = "MotifMatrix", plot = TRUE)

mm <- getMatrixFromProject(
  ArchRProj = proj,
  useMatrix = "MotifMatrix")
#deviations <- as.data.frame(assay(mm))
deviations <- as.data.frame(as.matrix(assay(mm)))


x <- apply(deviations,1,sd)
dfsd <- data.frame(
	gene=names(x),
	sd=x
)


#should compute rank within each batch!
proj@cellColData$log_telo_norm <- log10(1+proj@cellColData$dedup_telo/proj@cellColData$nFrags*1e8)
proj$rank_norm_telo<-NA
for(i in unique(proj$Sample)){
  allrank <- rank(proj$log_telo_norm[proj$Sample==i])
  proj$rank_norm_telo[proj$Sample==i] <- allrank/max(allrank)
}



################################################################
############# Calculate corr for TFs ###########################
################################################################

### correlate to one variable
corr_to_var <- function(proj,deviations,use_var) {
	print(use_var)

	metadata <- data.frame(
		telo=proj@cellColData[,use_var],
		cluster=factor(proj$Clusters)
	)

	if(length(unique(metadata$cluster))==1) {
		design <- model.matrix(data=metadata, ~telo)
	} else {
		design <- model.matrix(data=metadata, ~telo + cluster)
	}

	fit <- lmFit(deviations, design)
	fit <- eBayes(fit)
	tt <- topTable(fit,n=ncol(deviations),coef="telo")

	df <- data.frame(
		gene=rownames(tt),
		cor_telo=tt$logFC,
		p_telo=tt$P.Value
	)
	df
}

### correlate to all interesting vars
corr_to_allvar <- function(proj, deviations){
	df1 <- corr_to_var(proj, deviations, "rank_norm_telo")
	colnames(df1) <- c("gene","cor_cnt_telo","p_rank_telo")
	df1
}

print("======== running correlations ========")

####### correlate by cell type
proj@cellColData$Clusters <- proj@cellColData$zhang_ct
all_df <- corr_to_allvar(proj, deviations)
all_df <- all_df[order(all_df$cor_cnt_telo),]

#neg: MEF; AR; 
#pos: GRHL1, TP63, HNF1B,   FOSL up there too
tail(all_df,n=20)

write.csv(all_df,"/home/mahogny/jupyter/zhang_for_telo/allcorr_ct.csv")



####### correlate by sample					error. set to ct before
proj@cellColData$Clusters <- proj@cellColData$Sample
all_df <- corr_to_allvar(proj, deviations)
write.csv(all_df,"/home/mahogny/jupyter/zhang_for_telo/allcorr_sample.csv")



####### correlate by sample and cell type
proj@cellColData$Clusters <- sprintf("%s__%s",proj@cellColData$zhang_ct,proj@cellColData$Sample)
all_df <- corr_to_allvar(proj, deviations)
write.csv(all_df,"/home/mahogny/jupyter/zhang_for_telo/allcorr_ct_sample.csv")


################################################################
############# Paper plotting - nTA matrix #######################
################################################################


#in other zhang.R file


################################################################
############# Paper plotting - volcanoes #######################
################################################################


####### Uncorrected

all_df <- read.csv("/home/mahogny/jupyter/zhang_for_telo/allcorr_ct.csv")
all_df$log_p <- -log10(all_df$p_rank_telo)
all_df <- all_df[order(all_df$p_rank_telo),]

all_df$gene2 <- str_split_fixed(all_df$gene,"_",2)[,1]

all_df$gene2[all_df$log_p < 5] <- NA
all_df$gene2[(all_df$gene2 %in% c("MEF2B","MEF2A","SNAI2","HNF1A","AC0021266"))] <- NA

ggplot(all_df, aes(cor_cnt_telo, log_p, label=gene2))+geom_point()+geom_label(size=2.5)+
  xlab("Motif ~ nTA")+ylab("-Log10 p-value")+
  coord_flip()
ggsave("~/zhang_motif_uncorrected.pdf", width = 4, height = 6)

head(all_df)
all_df

########  by sample
# AR high here

all_df <- read.csv("/home/mahogny/jupyter/zhang_for_telo/allcorr_sample.csv")
all_df$log_p <- -log10(all_df$p_rank_telo)
all_df <- all_df[order(all_df$p_rank_telo),]

all_df$gene2 <- str_split_fixed(all_df$gene,"_",2)[,1]

all_df$gene2[all_df$log_p < 2] <- NA
#all_df$gene2[(all_df$gene2 %in% c("MEF2B","MEF2A","SNAI2","HNF1A","AC0021266))] <- NA

ggplot(all_df, aes(cor_cnt_telo, log_p, label=gene2))+geom_point()+geom_label(size=2.5)+
  coord_flip()

head(all_df)




########  by sample and ct
# AR high here

all_df <- read.csv("/home/mahogny/jupyter/zhang_for_telo/allcorr_ct_sample.csv")
all_df$log_p <- -log10(all_df$p_rank_telo)
all_df <- all_df[order(all_df$p_rank_telo),]

all_df$gene2 <- str_split_fixed(all_df$gene,"_",2)[,1]

all_df$gene2[all_df$log_p < 2] <- NA
all_df$gene2[(all_df$gene2 %in% c("RELA","MEF2A","MEF2B","AC0021266","PBX1","HOXD11","TFCP2L1"))] <- NA

ggplot(all_df, aes(cor_cnt_telo, log_p, label=gene2))+geom_point()+geom_label(size=2.5)+
  xlab("Motif ~ nTA + cell type:sample")+ylab("-Log10 p-value")+
  coord_flip()
ggsave("~/zhang_motif_corrected.pdf", width = 4, height = 6)

head(all_df,n=40)


# ################################################################
# ############# Calculate corr for TFs, ct specific ##############
# ################################################################
# 
# 
# print("=============== for each ct ================")
# all_ct <- proj@cellColData$zhang_ct
# proj@cellColData$Clusters <- proj@cellColData$zhang_ct
# for(current_ct in all_ct){
# 	print(current_ct)
# 
# 	sub_deviations <- deviations[proj@cellColData$Clusters==current_ct,]
# 	sub_proj <- proj[proj@cellColData$Clusters==current_ct,]
# 	#proj@cellColData$zhang_ct <- "one_ct"
# 	print(dim(sub_deviations))
# 
# 	all_df <- corr_to_allvar(sub_proj, sub_deviations)
# 	write.csv(all_df,sprintf("corrforct_%s.csv",current_ct))
# }
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ############## post-analysis
# if(FALSE){
# 
# 
# 	dat <- read.csv("allcorr_ct.csv")
# 	dat$cor_cnt_telo <- dat$cor_cnt_telo*dat$sd
# 	dat$cor_log_telo <- dat$cor_log_telo*dat$sd
# 	dat$cor_cnt_log_telo_norm <- dat$cor_cnt_log_telo_norm*dat$sd
# 	dat <- dat[order(dat$cor_cnt_telo),]
# 
# 	head(dat[,c("gene","cor_cnt_telo","cor_log_telo","cor_cnt_log_telo_norm")],n=20)
# 	tail(dat[,c("gene","cor_cnt_telo","cor_log_telo","cor_cnt_log_telo_norm")],n=20)
# 
# 	dat <- dat[order(dat$cor_cnt_log_telo_norm),]
# 	head(dat[,c("gene","cor_cnt_telo","cor_log_telo","cor_cnt_log_telo_norm")],n=20)
# 	tail(dat[,c("gene","cor_cnt_telo","cor_log_telo","cor_cnt_log_telo_norm")],n=20)
# 
# 
# library(chromVARmotifs)
# library(TFBSTools)
# data("human_pwms_v2")
# pfm <- human_pwms_v2
# id2sym <- data.frame(
#         id=sapply(1:length(pfm),  function(x) ID(pfm[[x]]) ),
#         symbol=sapply(1:length(pfm),  function(x) name(pfm[[x]]) )
# )
# 
# id2sym[grepl("FOSL2",id2sym$id),]
# id2sym[grepl("CGBP",id2sym$id),]
# pfm[["BAD931201_LINE1717_CGBP_D_N3"]]
# 
# }
# 
# 
# # ..(A/T)ACGT(A/T)   is ENSG00000250542_156   AACGTT (T box)   not really
