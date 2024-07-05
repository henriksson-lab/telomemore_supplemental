### To include in revision
### To include in revision
### To include in revision


################################################################################
######################### nTA per chromosome ###################################
################################################################################


library(matrixStats) 

samplemeta <- read.csv("/corgi/otherdataset/xingxi_atacsee/cellcycle/samplemeta.csv",sep="\t")

rels <- read.csv("/home/mahogny/jupyter/xingxi/relsize.csv",sep="\t")  #file sizes??? get proper bam!
rownames(rels) <- rels$sra

allfiles <- list.files("/home/mahogny/jupyter/xingxi/")
allfiles <- allfiles[str_ends(allfiles, "telo.bam")]
cnt <- matrix(ncol=length(allfiles), nrow=26)
colnames(cnt) <- allfiles

for(f in allfiles){
  print(f)
  dat <- read.table(pipe(paste0("samtools idxstat /home/mahogny/jupyter/xingxi/",f)))
  cnt[,f] <- dat[,3]#/sum(dat[,3])
  print(sum(dat[,3]))
}
rownames(cnt) <- dat[,1]
colnames(cnt) <- str_replace_all(colnames(cnt),"_1.fastq.gz.unfiltered.bam.telo.bam","")

#normalize by reads in file
for(i in 1:ncol(cnt)){
  cnt[,i] <- cnt[,i]/rels[colnames(cnt)[i],]$numb * 1e8
}


toplot <- reshape2::melt(as.matrix(as.data.frame(cnt)))
colnames(toplot) <- c("chr","Run","nTA")
toplot <- merge(toplot, samplemeta)

ggplot(toplot,aes(paste(cc_phase,Run),chr, fill=nTA)) + geom_tile() + 
  coord_flip() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  xlab("") + ylab("") + 
  theme(panel.background = element_blank())
ggsave("xingxi_nta_per_chrom.svg",width = 6, height = 2)

# Run	Cell_Line	Cell_type	geo
# SRR3336945	GM12878	FACS sorted G1 high group	GSM2108582
# SRR3336946	GM12878	FACS sorted G1 high group	GSM2108583
# SRR3336947	GM12878	FACS sorted G1 low group	GSM2108584
# SRR3336948	GM12878	FACS sorted G1 low group	GSM2108585
# SRR3336949	GM12878	FACS sorted S group	GSM2108586
# SRR3336950	GM12878	FACS sorted S group	GSM2108587
# SRR3336951	GM12878	FACS sorted G2 group	GSM2108588
# SRR3336952	GM12878	FACS sorted G2 group	GSM2108589


##### TODO why no G1 early? one missing file  -- 45


################################################################################
####################### histogram of tandem repeats ############################
################################################################################

samplemeta <- read.csv("/corgi/otherdataset/xingxi_atacsee/cellcycle/samplemeta.csv",sep="\t")


alldat <- list()
for(f in list.files("/home/mahogny/jupyter/xingxi")){  #not a good place! TODO
  if(str_ends(f,"telo.bam.hist")){
    dat <- read.csv(paste0("/home/mahogny/jupyter/xingxi/",f),sep="\t")

    onet <- table(rowMaxs(as.matrix(dat)))
    onet <- data.frame(
      reps=names(onet),
      cnt=as.integer(onet)
    )
    onet$cnt <- onet$cnt/sum(onet$cnt)
    onet$Run <- f
    alldat[[f]] <- onet
  }
}
alldat <- do.call(rbind, alldat)  
alldat$reps <- as.integer(alldat$reps)
alldat$Run <- str_split_fixed(alldat$Run,"_",2)[,1]
alldat <- merge(alldat, samplemeta)

ggplot(alldat, aes(reps, cnt, fill=cc_phase, group=paste(cc_phase,Run))) + 
  geom_bar(stat="identity", position = "dodge") +
  xlab("Number of tandem repeats") +
  ylab("Frequency") +
  scale_x_continuous(breaks=seq(3,10,by=1))
ggsave("xingxi_tandem_rep.svg",width = 6, height = 2)










################################################################################
######################### pileup for paper #####################################
################################################################################



bw1 <- rtracklayer::import.bw("/corgi/otherdataset/xingxi_atacsee/cellcycle/t2t_bam/SRR3336946_1.fastq.gz.unfiltered.bam.telo.bw")


numpos <- round(300e6/50) #maximum size
aq <- makeGRangesFromDataFrame(data.frame(seqname=rep("chr1",numpos),start=((1:numpos)-1)*50+1,end=((1:numpos)-1)*50+1))

tol <- findOverlaps(aq, bw1)
all_score <- data.frame(
  pos=start(ranges(aq))[as.data.frame(tol)$queryHits],
  score1=bw1$score[as.data.frame(tol)$subjectHits]
)


#### Full size
binsize <- max(red_score$pos_red)/10000 #if too small, cannot see when zoomed out
all_score$pos_red <- round(all_score$pos/binsize)*binsize
red_score <- sqldf::sqldf("select pos_red, avg(score1) as score1 from all_score group by pos_red")
red_score$score1 <- red_score$score1/sum(red_score$score1)

p0 <- ggplot(red_score, aes(pos_red, score1)) + geom_line() + ylab("") + xlab("Chr1 (bp)")
ggsave(plot = p0, "out_xingxi_pileup_all.svg",width = 8,height = 1.5)

#### Zoomed in
binsize <- 100
all_score$pos_red <- round(all_score$pos/binsize)*binsize
red_score <- sqldf::sqldf("select pos_red, avg(score1) as score1 from all_score group by pos_red")
red_score$score1 <- red_score$score1/sum(red_score$score1)

dx <- 20000
p1 <- ggplot(red_score, aes(pos_red, score1)) + geom_line() + xlim(c(0,dx)) + ylab("") + xlab("Chr1 (bp)")
p2 <- ggplot(red_score, aes(pos_red, score1)) + geom_line() + xlim(c(max(red_score$pos_red)-dx,max(red_score$pos_red))) + ylab("") + xlab("Chr1 (bp)")
ptot <- p1|p2
ptot
ggsave(plot = ptot, "out_xingxi_pileup_telo.svg",width = 8,height = 1.5)




################################################################################
######################### MASH comparison of libs ##############################
################################################################################

mash_telo <- read.table("/corgi/otherdataset/xingxi_atacsee/cellcycle/mash.telo.csv")[,1:3]
colnames(mash_telo) <- c("lib1","lib2","dist")
mash_telo$lib1 <- str_replace_all(mash_telo$lib1, ".telo.fq", "")
mash_telo$lib2 <- str_replace_all(mash_telo$lib2, ".telo.fq", "")

mash_telo <- mash_telo[mash_telo$lib1!=mash_telo$lib2,]

ggplot(mash_telo, aes(lib1,lib2, fill=dist)) + 
  geom_tile() + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  xlab("Telo reads") + ylab("") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



mash_all <- read.table("/corgi/otherdataset/xingxi_atacsee/cellcycle/mash.all.csv")[,1:3]
colnames(mash_all) <- c("lib1","lib2","dist")
mash_all$lib1 <- str_replace_all(mash_all$lib1, ".all.fq", "")
mash_all$lib2 <- str_replace_all(mash_all$lib2, ".all.fq", "")

mash_all <- mash_all[mash_all$lib1!=mash_all$lib2,]

ggplot(mash_all, aes(lib1,lib2, fill=dist)) + 
  geom_tile() + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  xlab("All reads") + ylab("") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

