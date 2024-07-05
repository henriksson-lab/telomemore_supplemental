if(FALSE){
  devtools::install_github("omarwagih/ggseqlogo")
  
}


### To include in revision
### To include in revision
### To include in revision



################################################################################
########################### Motifs near ends ###################################
################################################################################

 
library(stringr)

# just color all reads
# grep --color -E '^|CTGTCTCTTATACACATCT|AGATGTGTATAAGAGACAG' reads_manyinsert.txt

# GREP_COLOR='01;36' grep --color -E '^|CTGTCTCTTATACACATCT|AGATGTGTATAAGAGACAG' reads_manyinsert.txt --color=always | grep --color -E '^|TTAGGG|CCCTAA'


compute_motif <- function(fname){
  dat <- read.table(fname)
  #dat <- dat[str_sub(dat$V1,20,20)=="G",,drop==FALSE]
  #dat <- dat[str_sub(dat$V1,20,20)=="C",,drop=FALSE]
  motif_cnt <- matrix(ncol=str_length(dat$V1[1]), nrow=4)
  for(i in 1:ncol(motif_cnt)){
    motif_cnt[,i] <- as.integer(table(factor(str_sub(dat$V1,i,i), levels=c("A","T","C","G"))))
    motif_cnt[,i] <- motif_cnt[,i]/sum(motif_cnt[,i])
  }
  rownames(motif_cnt) <- c("A","T","C","G")
  #round(motif_cnt, digit=2)-0.25
  motif_cnt
}

compute_motif("/husky/fromsequencer/ira20240508/out/motif1s.txt")
compute_motif("/husky/fromsequencer/ira20240508/out/motif1e.txt")
compute_motif("/husky/fromsequencer/ira20240508/out/motif2s.txt")
compute_motif("/husky/fromsequencer/ira20240508/out/motif2e.txt")

#round(motif_cnt, digit=3)


#library(motifStack)
#motif <- new("pcm", mat=compute_motif("/husky/fromsequencer/ira20240508/out/motif1s.txt"), name="foo")
#plot(motif)




#the two are largely equivalent, as we should have expected
library(ggseqlogo)
ptot <- ggseqlogo(compute_motif("/husky/fromsequencer/ira20240508/out/motif1s.txt")) | ggseqlogo(compute_motif("/husky/fromsequencer/ira20240508/out/motif1e.txt"))
ptot
ggsave(plot = ptot, "./motif_at_end.svg")
ggseqlogo(compute_motif("/husky/fromsequencer/ira20240508/out/motif2s.txt")) | ggseqlogo(compute_motif("/husky/fromsequencer/ira20240508/out/motif2e.txt"))



################################################################################
################## Lengths of fragments ########################################
################################################################################

if(FALSE){
  
  telolen <- read.table("/husky/fromsequencer/ira20240508/out/telo_vs_len.txt")
  telolen
  
  colnames(telolen) <- c("p1s","p1e","p2s","p2e","t1","t2","len")
  hist(telolen$len)
  
  hist(telolen$len[telolen$p1s>1])
  hist(telolen$len[telolen$p1s>0])
  ggplot(telolen, aes(len, fill=p1s>0)) + geom_histogram(aes(y=after_stat(density)),bins=100)
  #NOTE: there are some really short reads. extract < 5kb
  
  ggplot(telolen, aes(len, fill=p1e>0)) + geom_histogram(aes(y=after_stat(density)),bins=100)
  ggplot(telolen, aes(len, fill=p2s>0)) + geom_histogram(aes(y=after_stat(density)),bins=100)
  ggplot(telolen, aes(len, fill=p2e>0)) + geom_histogram(aes(y=after_stat(density)),bins=100)
  
  
  ggplot(telolen, aes(len, fill=t2>1)) + geom_histogram(aes(y=after_stat(density)),bins=100)
  
  
  #ggplot(telolen, aes(p1s, ))
  
  #table(telolen$p1s) 
  table(telolen$p1e)
  table(telolen$p2s)
  table(telolen$p2e)
  table(telolen$t1)
  table(telolen$t2)
  
  
  mean(telolen$p1s) 
  mean(telolen$p1e)
  mean(telolen$p2s)
  mean(telolen$p2e)
  mean(telolen$t1)
  mean(telolen$t2)
  
  
  # samtools view pr_080_001.bam   | grep -E "CCCTAACCCTAACCCTAACCCTAA|TTAGGGTTAGGGTTAGGGTTAGGG" >  4telo.pr_080_001.txt
  
}





################################################################################
################## Targeted analysis of telomere motif #########################
################################################################################


library(Biostrings)
#dna = DNAStringSet(c("ATCTCGGCGCGCATCGCGTACGCTACTAGC", "ACCGCTA"))
#complement(dna)


read_motifs_len <- function(len){
  dat <- c(
    read.table("/husky/fromsequencer/ira20240508/out/motif1s.txt")$V1,
    read.table("/husky/fromsequencer/ira20240508/out/motif2s.txt")$V1,
    as.character(reverseComplement(DNAStringSet(read.table("/husky/fromsequencer/ira20240508/out/motif1e.txt")$V1))),
    as.character(reverseComplement(DNAStringSet(read.table("/husky/fromsequencer/ira20240508/out/motif2e.txt")$V1)))
  )
  dat_s <- c(
    read.table("/husky/fromsequencer/ira20240508/out/motif1s.txt")$V1,
    read.table("/husky/fromsequencer/ira20240508/out/motif2s.txt")$V1
  )
  dat_s <- str_sub(dat_s,str_length(dat_s[1])-len+1,str_length(dat_s[1]))
  dat_e <- c(
    read.table("/husky/fromsequencer/ira20240508/out/motif1e.txt")$V1,
    read.table("/husky/fromsequencer/ira20240508/out/motif2e.txt")$V1
  )
  dat_e <- str_sub(dat_e,1,len)
  df <- as.data.frame(table(c(dat_s,dat_e)))
  colnames(df) <- c("kmer","cnt")
  rownames(df) <- df$kmer
  df <- df[order(df$cnt, decreasing = TRUE),]
  df
}
 



read_motifs_special <- function(){
  len <- 12
  dat <- c(
    read.table("/husky/fromsequencer/ira20240508/out/motif1s.txt")$V1,
    read.table("/husky/fromsequencer/ira20240508/out/motif2s.txt")$V1,
    as.character(reverseComplement(DNAStringSet(read.table("/husky/fromsequencer/ira20240508/out/motif1e.txt")$V1))),
    as.character(reverseComplement(DNAStringSet(read.table("/husky/fromsequencer/ira20240508/out/motif2e.txt")$V1)))
  )
  dat_s <- c(
    read.table("/husky/fromsequencer/ira20240508/out/motif1s.txt")$V1,
    read.table("/husky/fromsequencer/ira20240508/out/motif2s.txt")$V1
  )
  dat_s <- str_sub(dat_s,str_length(dat_s[1])-len+1,str_length(dat_s[1]))
  dat_e <- c(
    read.table("/husky/fromsequencer/ira20240508/out/motif1e.txt")$V1,
    read.table("/husky/fromsequencer/ira20240508/out/motif2e.txt")$V1
  )
  dat_e <- str_sub(dat_e,1,len)
  dat <- c(dat_s,dat_e)
  
  dat <- paste0(
    str_sub(dat,12-8,12-8),
    str_sub(dat,12-5,12-5),
    str_sub(dat,12-4,12-4),
    str_sub(dat,12-3,12-3),
    str_sub(dat,12,12)
  )
  
  df <- as.data.frame(table(dat))
  colnames(df) <- c("kmer","cnt")
  rownames(df) <- df$kmer
  df <- df[order(df$cnt, decreasing = TRUE),]
  #G**TAG**G
  df
} 


plot_complex_motif <- function(df, kmers){
  df <- df[kmers,]
  df$kmer <- kmers
  df$cnt[is.na(df$cnt)] <- 0
  df$frac <- df$cnt/sum(df$cnt)
  ggplot(df, aes(kmer, frac)) + geom_bar(stat = "identity") + coord_flip() + ylab("Fraction") + theme(axis.text.x = element_text(family="courier"))
}


df <- read_motifs_len(6)
plot_complex_motif(df, c("TTAGGG","TAGGGT","AGGGTT","GGGTTA","GGTTAG","GTTAGG")) + ylim(0,1) 
ggsave("./targeted_6bp.svg", height = 2, width = 3)


df <- read_motifs_len(12)
plot_complex_motif(df,c("TTAGGGTTAGGG","TAGGGTTAGGGT","AGGGTTAGGGTT","GGGTTAGGGTTA","GGTTAGGGTTAG","GTTAGGGTTAGG")) + ylim(0,1)
ggsave("./targeted_12bp.svg", height = 2, width = 4)



################################################################################
########################### Motifs around middle inserts #######################
################################################################################


compute_motif <- function(fname){
  dat <- read.table(fname)
  dat <- dat$V1
  dat <- c(dat, as.character(reverseComplement(DNAStringSet(dat))))
  
  motif_cnt <- matrix(ncol=str_length(dat[1]), nrow=4)
  for(i in 1:ncol(motif_cnt)){
    motif_cnt[,i] <- as.integer(table(factor(str_sub(dat,i,i), levels=c("A","T","C","G"))))
    motif_cnt[,i] <- motif_cnt[,i]/sum(motif_cnt[,i])
  }
  rownames(motif_cnt) <- c("A","T","C","G")
  #round(motif_cnt, digit=2)-0.25
  motif_cnt
}


#the two are largely equivalent, as we should have expected
library(ggseqlogo)
ggseqlogo(compute_motif("/husky/fromsequencer/ira20240508/out/motif2.noc.txt"))
ggsave("./motif_mid.svg", width = 8, height = 9.3/2)


if(FALSE){
  ## Is CCCTAA uniformly distributed?
  dat <- read.table("/husky/fromsequencer/ira20240508/out/motif2.noc.txt")
  dat <- dat$V1
  dat <- c(dat, as.character(reverseComplement(DNAStringSet(dat))))
  
  library(stringr)
  str_locate_all(dat,"TTAGGG")
  
  pos_motif <- c(
    do.call(rbind,str_locate_all(dat,"TTAGGG"))[,1],
    do.call(rbind,str_locate_all(dat,"CCCTAA"))[,1]
  )
  
  df <- data.frame(pos=pos_motif)
  pos_hist <- sqldf::sqldf("select pos, count(pos) as cnt from df group by pos")
  ggplot(pos_hist, aes(pos+3-1-str_length(dat[1])/2, cnt)) + geom_bar(stat = "identity")
}

################################################################################
########################### Insertion pos vs read ##############################
################################################################################


if(FALSE){
  
  dat <- read.table("/husky/fromsequencer/ira20240508/out/pos_t1.txt")
  colnames(dat) <- c("at","len","chr","read_pos")
  p1 <- ggplot(dat,aes(at,len)) + geom_point(cex=0.1) + geom_abline(slope=1,color="red") + xlab("pos, TTAGGG x 8")
  dat <- read.table("/husky/fromsequencer/ira20240508/out/pos_t2.txt")
  colnames(dat) <- c("at","len","chr","read_pos")
  p2 <- ggplot(dat,aes(at,len)) + geom_point(cex=0.1) + geom_abline(slope=1,color="red") + xlab("pos, CCCTAA x 8")
  dat <- read.table("/husky/fromsequencer/ira20240508/out/pos_tn.txt")
  colnames(dat) <- c("at","len","chr","read_pos")
  p3 <- ggplot(dat,aes(at,len)) + geom_point(cex=0.1) + geom_abline(slope=1,color="red") + xlab("pos, CTGTCTCTTATACACATCTA (tn5 insert)")
  egg::ggarrange(p1,p2,p3)
}


################################################################################
########################### Does tn5 truncate telomeric reads? #################
################################################################################

if(FALSE){
  # read length if there is telomeres "to the right", with and wihout an insertion to the end
  #
  #(base) [mahogny@r33 align]$ grep -c CTGTCTCTTATACACATCTA  xteloreads.txt 
  #29
  #(base) [mahogny@r33 align]$ wc -l xteloreads.txt 
  #288 xteloreads.txt
  
  
  #compare_insert_pos.py
  
  if(FALSE){
    dat <- read.table("/husky/fromsequencer/ira20240508/out/pos_trans_vs_ins.txt")
    colnames(dat) <- c("type","totlen","avg_pos_tn","totlen2","num_t1","num_t2", "avgpos_t1","avgpos_t2")
    #4x hexamer
    dat <- dat[dat$num_t1>10 | dat$num_t2>10,]
    
    #dat[dat$avg_pos_tn<dat$totlen*0.2,]
    dat[dat$avg_pos_tn<200,]
    
    dat$has_tn_left <- dat$avg_pos_tn<200
  }
  
  
  dat <- read.table("/husky/fromsequencer/ira20240508/out/pos_trans_vs_ins2.txt")
  colnames(dat) <- c("type","totlen","num_tn","avg_pos_tn","num_t1","num_t2", "avgpos_t1","avgpos_t2")
  
  #dat <- dat[dat$num_t1>10 | dat$num_t2>10,]
  dat$has_tn_left  <- dat$num_tn>0 & dat$avg_pos_tn<200
  dat$has_tn_right <- dat$num_tn>0 & dat$avg_pos_tn>dat$totlen-200 
  
  ggplot(dat, aes(totlen, avgpos_t2, color=has_tn_left)) + geom_point()
  ggplot(dat, aes(totlen, avgpos_t1, color=has_tn_left)) + geom_point() #if cut, then mostly on one side
  
  ggplot(dat, aes(totlen, num_t1, color=has_tn_left)) + geom_point() 
  ggplot(dat, aes(totlen, num_t2, color=has_tn_left)) + geom_point() #this is a control for above! or can be made one
  
  ggplot(dat, aes(totlen, num_t2*4*6, color=has_tn_right)) + geom_point() + 
    geom_smooth(method = "lm",formula=y~0+x) +
    xlab("Read length") + ylab("Hexamer cumulative length")  ### not ideal still
  
  ggplot(dat, aes(totlen, num_t1*4*6, color=has_tn_left)) + geom_point() + 
    geom_smooth(method = "lm",formula=y~0+x) +
    xlab("Read length") + ylab("Hexamer cumulative length")
  
  ##### can also look at fraction of read length! but not a pretty plot
  ggplot(dat, aes(has_tn_right, num_t2*4*6/totlen)) + geom_point()
  ggplot(dat, aes(has_tn_left, num_t1*4*6/totlen)) + geom_point()
}


################################################################################
#################### telomere length vs cell line ##############################    ratio: all of the data
################################################################################ 


dat <- read.csv("/husky/fromsequencer/ira20240508/out/bulk_count.txt")
dat$freq <- dat$cnt/dat$tot
p_telocnt_all <- ggplot(dat, aes(cell, freq, fill=motif)) + 
  geom_bar(stat = "identity", position = "dodge") +
  xlab("")+
  ylab("Frequency")+
  theme_bw() +
  theme(
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+  
  theme(panel.grid.major.y = element_line(colour = "#808080")) 
p_telocnt_all

df <- merge(
  data.frame(motif=dat$motif[dat$cell=="CCRF"],freq_ccrf=dat$freq[dat$cell=="CCRF"]),
  data.frame(motif=dat$motif[dat$cell=="Jurkat"],freq_jurkat=dat$freq[dat$cell=="Jurkat"])
)
df$ratio <- df$freq_jurkat/df$freq_ccrf
df
p_all <- ggplot(df, aes(motif, ratio)) + geom_bar(stat = "identity", position = "dodge") + 
  theme_bw() +
  theme(
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+  
  theme(panel.grid.major.y = element_line(colour = "#808080")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



################################################################################
#################### telomere length vs cell line ##############################   ratio: reads with inserts
################################################################################

dat <- read.csv("/husky/fromsequencer/ira20240508/out/wi.telocnt")
dat <- dat[!str_ends(dat$motif,"6"),]
dat$freq <- dat$cnt/dat$tot
p_telocnt_ins <- ggplot(dat, aes(cell, freq, fill=motif)) + 
  geom_bar(stat = "identity", position = "dodge")+
  xlab("")+
  ylab("Frequency")+
  theme_bw() +
  theme(
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+  
  theme(panel.grid.major.y = element_line(colour = "#808080"))

p_telocnt_all / p_telocnt_ins

df <- merge(
  data.frame(motif=dat$motif[dat$cell=="CCRF"],freq_ccrf=dat$freq[dat$cell=="CCRF"]),
  data.frame(motif=dat$motif[dat$cell=="Jurkat"],freq_jurkat=dat$freq[dat$cell=="Jurkat"])
)
df$ratio <- df$freq_jurkat/df$freq_ccrf
df
p_ins <- ggplot(df, aes(motif, ratio)) + geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  theme(
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+  
  theme(panel.grid.major.y = element_line(colour = "#808080")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#ptot <- (p_all + ylab("Ratio: all reads")) / (p_ins + ylab("Ratio: reads with inserts"))
ptot <- (p_all + ylab("Ratio") + xlab("") + ylim(0,1)) / (p_ins + ylab("Ratio") + xlab("") + ylim(0,1))
ptot
ggsave(plot=ptot, "./ratio.svg", height = 4, width = 2)


ptot <- p_telocnt_all / p_telocnt_ins
ptot
ggsave(plot=ptot, "./abundance_wgs.svg", height = 4, width = 4)

################################################################################
#################### telomere coverage vs the rest #############################
################################################################################ 

dat <- read.csv("telocov.csv")
dat$freq <- dat$ins/dat$tot
ggplot(dat, aes(reg, freq, fill=celltype)) + 
  geom_bar(stat="identity", position="dodge") + 
  xlab("") +
  ylab("Fraction of reads with Tn5 inserts") + 
  theme_bw() +
  theme(
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+  
  theme(panel.grid.major.y = element_line(colour = "#808080")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("./fraction_tn5.svg", height = 4, width = 2)


################################################################################
#################### Normalized number of inserts per bin / kb  ################
################################################################################ 

if(FALSE){
  library(rtracklayer)
  
  bw <- rtracklayer::import.bw("/husky/fromsequencer/ira20240508/bw/bc2085.bw") #'http://people.binf.ku.dk/~bwp982/RNAseq_Claudia/01_Ars2-tot-1-repplus.bw')
  curchr <- "chr1"
  
  df <- data.frame(
    chr=seqnames(bw),
    posleft=floor(start(ranges(bw))/5000000),
    score=bw$score
  )
  mean_depth <- sqldf::sqldf("select chr, posleft, avg(score) as depth from df group by posleft, chr order by chr, posleft")
  mean_depth <- merge(mean_depth, sqldf::sqldf("select chr, avg(depth) as mean_depth, max(posleft) as len from mean_depth group by chr"))
  #mean_depth$depth <- mean_depth$depth - mean_depth$mean_depth
  mean_depth$norm_pos <- mean_depth$posleft/mean_depth$len
  
  ggplot(mean_depth[mean_depth$chr=="chr1",], aes(posleft, depth)) + geom_line()
  ggplot(mean_depth[mean_depth$chr=="chr3",], aes(posleft, depth)) + geom_line() #why so weird?
  ggplot(mean_depth[mean_depth$chr=="chr4",], aes(posleft, depth)) + geom_line()
  ggplot(mean_depth[mean_depth$chr=="chr5",], aes(posleft, depth)) + geom_line()
  ggplot(mean_depth[mean_depth$chr=="chr6",], aes(posleft, depth)) + geom_line() #why so weird?
  ggplot(mean_depth[mean_depth$chr=="chr7",], aes(posleft, depth)) + geom_line() #why so weird?
  ggplot(mean_depth[mean_depth$chr=="chr8",], aes(posleft, depth)) + geom_line() #why so weird?
  ggplot(mean_depth[mean_depth$chr!="chrM" & mean_depth$chr!="chrY" & mean_depth$chr!="chr6",], aes(posleft, depth, color=chr)) + geom_line()
  ggplot(mean_depth[mean_depth$chr!="chrM" & mean_depth$chr!="chrY" & mean_depth$chr!="chr6",], aes(norm_pos, depth, color=chr)) + geom_line()
  
  ggplot(mean_depth[mean_depth$chr!="chrM" & mean_depth$chr!="chrY" & mean_depth$chr!="chr6",], aes(norm_pos, depth)) + 
    geom_point() + 
    geom_smooth(method="lm", formula=  y ~ splines::bs(x, knots = c(25, 40, 60), degree = 3))
  
  
  mean_depth[mean_depth$depth>100,]
  
  
  ######## compare ins vs depth
  #chr14 is odd. coverage up at ends. not tn5? 
  allins <- read.table("/husky/fromsequencer/ira20240508/bw/bc2085.absinspos.bed")
  df <- data.frame(
    chr=allins$V1,
    posleft=floor(allins$V2)/5000000
  )
  sum_ins <- sqldf::sqldf("select chr, posleft, count(*) as cnt_ins from df group by posleft, chr order by chr, posleft")
  forchr <- "chr22"
  ggplot(sum_ins[sum_ins$chr==forchr,], aes(posleft, cnt_ins)) + 
    geom_line() +
    geom_line(data = mean_depth[mean_depth$chr==forchr,], mapping = aes(posleft, depth),color="blue")  
  
  
  
  cumpos <- allins[allins$V1=="chr1",]
  cumpos$rank <- 1:nrow(cumpos)
  cumpos$frac_len <- cumpos$V2/max(cumpos$V2)
  plot(cumpos$frac_len, (cumpos$rank - cumpos$frac_len*nrow(cumpos))/nrow(cumpos), type="l")
  for(i in 2:20){
    cumpos <- allins[allins$V1==paste0("chr",i),]
    cumpos$rank <- 1:nrow(cumpos)
    cumpos$frac_len <- cumpos$V2/max(cumpos$V2)
    lines(cumpos$frac_len, (cumpos$rank - cumpos$frac_len*nrow(cumpos))/nrow(cumpos), type="l")
  }
  
}

