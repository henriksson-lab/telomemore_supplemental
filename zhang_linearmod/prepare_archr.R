library(stringr)

datadir <- readLines("setting.dataset")

################################################################
############### Prepare archR ##################################
################################################################
archr_genome <- readLines(sprintf("%s/setting.archr_genome",datadir))

library(ArchR)
set.seed(1)
addArchRThreads(threads = 16) 
addArchRGenome(archr_genome)   #hg38 or mm10


################################################################
################ set up data ###################################
################################################################


if(file.exists("list_samples.csv")){
	#Use selected subset
	ArrowFiles <- read.table("list_samples.csv",header=FALSE,stringsAsFactors=FALSE)[,1]
	ArrowFiles <- sprintf("%s.arrow",ArrowFiles)
	print("using subset")
	print(ArrowFiles)
} else {
	#Pull all from the directory
	ArrowFiles <- list.files(sprintf("%s/arrows",datadir))
	ArrowFiles <- ArrowFiles[str_ends(ArrowFiles,"arrow")]
}


proj <- ArchRProject(
  ArrowFiles = sprintf("%s/arrows/%s",datadir,ArrowFiles), 
  outputDirectory = "archr_output",
  copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)

################################################################
############# Add telomere info ################################
################################################################


if(TRUE){
	####################### new java counter ################################

	if(grepl("zhang", datadir)){
		### hack only for zhang!
		all_kmers <- read.csv(sprintf("%s/summary_kmer.withct.csv.gz", datadir))
	        all_kmers$dataset <- str_replace(all_kmers$dataset,"_1","_rep1")
		#"ENC-1JKYN-012-SM-JF1O6_snATAC_body_of_pancreas_rep1"    to become
		#adipose_omentum_SM-CSSD4_rep1
	} else {
		###Normal way
		all_kmers <- read.csv(sprintf("%s/summary_kmer.java.csv", datadir))
	}
	rownames(all_kmers) <- sprintf("%s#%s",  all_kmers$dataset, all_kmers$barcode)
	#Only keep cells that have a kmer count (sufficient quality according to cellranger)
	proj <- proj[rownames(proj@cellColData) %in% rownames(all_kmers),]

	#Transfer telocount
	proj@cellColData$total_telo <- all_kmers[rownames(proj@cellColData),]$totalcnt_CCCTAA
	proj@cellColData$dedup_telo <- all_kmers[rownames(proj@cellColData),]$dedupcnt_CCCTAA

	if(grepl("zhang", datadir)){
		#Transfer celltype
		proj@cellColData$zhang_ct <- all_kmers[rownames(proj@cellColData),"cell.type"]		
	}

	if(nrow(proj@cellColData)==0){
	        printf("============666 no cells left")
	}



} else {
	####################### original counter by william ################################

	all_kmers <- read.csv(sprintf("%s/summary_kmer.csv", datadir))
	rownames(all_kmers) <- all_kmers[,1]
	### hack only for zhang!
	if(grepl("zhang", datadir)){
	        rownames(all_kmers) <- str_replace(rownames(all_kmers),"_1","_rep1")
	}
	#Only keep cells that have a kmer count (sufficient quality according to cellranger)
	proj <- proj[rownames(proj@cellColData) %in% rownames(all_kmers),]

	#Transfer telocount
	proj@cellColData$cnt_telo <- all_kmers[rownames(proj@cellColData),]$cnt / proj@cellColData$nFrags
	proj@cellColData$log_telo <- log10(0.001+proj@cellColData$cnt_telo)

	if(nrow(proj@cellColData)==0){
	        printf("============666 no cells left")
	}
}

################################################################
############# Basic analysis ###################################
################################################################

## cluster ... if not zhang_all, too large
if(!grepl("zhang_all", datadir)){

	proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI")

	if(TRUE){
	        proj <- addHarmony(
	            ArchRProj = proj,
	            reducedDims = "IterativeLSI",
	            name = "Harmony",
	            groupBy = "Sample"
	        )
	        proj <- addClusters(input = proj, reducedDims = "Harmony",force = TRUE)
	        proj <- addUMAP(ArchRProj = proj, reducedDims = "Harmony",force = TRUE)
	} else {
	        proj <- addClusters(input = proj, reducedDims = "IterativeLSI")
	        proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI")
	}
	proj <- addImputeWeights(proj)

}

################################################################
#################### save to file ##############################
################################################################

combined_dataset <- list(
        ArrowFiles=ArrowFiles,
        datadir=datadir,
        genome=archr_genome,
        proj=proj)
saveRDS(combined_dataset, "archr_proj.rds")




############### plotting #########################


## do not attempt to plot all of zhang
if(!grepl("zhang_all", datadir)){
	source("/corgi/henriksson/telomere/plot_one_archr.R")
}
