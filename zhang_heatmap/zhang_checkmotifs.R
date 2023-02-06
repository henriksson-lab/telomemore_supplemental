if(FALSE){
  
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("TFBSTools")
  
  devtools::install_github("GreenleafLab/chromVARmotifs")
} 

library(TFBSTools)
library(stringr)
library(seqLogo)
library(chromVARmotifs) #https://github.com/GreenleafLab/chromVARmotifs
library(ggplot2)

library(tidyr)
library(dplyr)
library(tibble)
library(S4Vectors)

######################################
# load motifs
data("human_pwms_v2")

pfm <- human_pwms_v2
id2sym <- data.frame(
  id=sapply(1:length(pfm),  function(x) ID(pfm[[x]]) ),
  symbol=sapply(1:length(pfm),  function(x) name(pfm[[x]]) )
)




.summarizeChromVARMotifs <- function(motifs = NULL){
  
  motifNames <- lapply(seq_along(motifs), function(x){
    namex <- make.names(motifs[[x]]@name)
    if(grepl("LINE", namex)){
      splitNamex <- stringr::str_split(motifs[[x]]@ID, pattern="\\_", simplify = TRUE)
      namex <- splitNamex[1, grep("LINE",splitNamex[1,]) + 1]
    }
    if(substr(namex,nchar(namex),nchar(namex))=="."){
      namex <- substr(namex,1,nchar(namex)-1)
    }
    namex <- paste0(namex, "_", x)
    namex
  }) %>% unlist(.)
  
  motifNames2 <- lapply(seq_along(motifs), function(x){
    namex <- make.names(motifs[[x]]@name)
    if(grepl("LINE", namex)){
      splitNamex <- stringr::str_split(motifs[[x]]@ID, pattern="\\_", simplify = TRUE)
      namex <- splitNamex[1, grep("LINE",splitNamex[1,]) + 1]
    }
    if(substr(namex,nchar(namex),nchar(namex))=="."){
      namex <- substr(namex,1,nchar(namex)-1)
    }
    namex
  }) %>% unlist(.)
  
  motifDF <- lapply(seq_along(motifs), function(x){
    df <- data.frame(
      row.names = motifNames[x],
      name = motifNames2[[x]],
      ID = motifs[[x]]@ID,
      strand = motifs[[x]]@strand,
      stringsAsFactors = FALSE
    )
  }) %>% Reduce("rbind", .) %>% DataFrame
  
  names(motifs) <- motifNames
  
  out <- list(motifs = motifs, motifSummary = motifDF)
  
  return(out)
  
}
chromvarlist <- .summarizeChromVARMotifs(human_pwms_v2)


#### Read ranked motifs
dat <- read.csv("zhang/corr/allcorr_ct_sample.csv")  #most conservative
dat <- read.csv("zhang/corr/allcorr_sample.csv")
dat <- read.csv("zhang/corr/allcorr_ct.csv")

#### Make linear fits comparable as not all motifs have the same SD (for better or worse)
dat$cor_cnt_telo <- dat$cor_cnt_telo*dat$sd
dat$cor_log_telo <- dat$cor_log_telo*dat$sd
dat$cor_cnt_log_telo_norm <- dat$cor_cnt_log_telo_norm*dat$sd
# 
# dat <- dat[order(dat$cor_cnt_telo),]
# head(dat[,c("gene","cor_cnt_telo","cor_log_telo","cor_cnt_log_telo_norm")],n=20)
# tail(dat[,c("gene","cor_cnt_telo","cor_log_telo","cor_cnt_log_telo_norm")],n=20)
# 
# dat <- dat[order(dat$cor_log_telo),]
# head(dat[,c("gene","cor_cnt_telo","cor_log_telo","cor_cnt_log_telo_norm")],n=20)
# tail(dat[,c("gene","cor_cnt_telo","cor_log_telo","cor_cnt_log_telo_norm")],n=20)
# 
# dat <- dat[order(dat$cor_cnt_log_telo_norm),]
# head(dat[,c("gene","cor_cnt_telo","cor_log_telo","cor_cnt_log_telo_norm")],n=20)
# tail(dat[,c("gene","cor_cnt_telo","cor_log_telo","cor_cnt_log_telo_norm")],n=20)
# 
# ######## compare on p-values
# 
# dat <- dat[order(dat$p_cnt_telo),]
# head(dat[,c("gene","cor_cnt_telo","cor_log_telo","cor_cnt_log_telo_norm","p_cnt_telo")],n=20) #most related
# tail(dat[,c("gene","cor_cnt_telo","cor_log_telo","cor_cnt_log_telo_norm","p_cnt_telo")],n=20) #least related
# ### brings out KLF1_179, CTCFL_198, CTCF_177, MGA_103 .. TP63_704    FOSL2_105  
# 
# dat <- dat[order(dat$p_log_telo),]
# head(dat[,c("gene","cor_cnt_telo","cor_log_telo","cor_cnt_log_telo_norm","p_log_telo")],n=20) #most related
# tail(dat[,c("gene","cor_cnt_telo","cor_log_telo","cor_cnt_log_telo_norm","p_log_telo")],n=20) #least related
# ### TP63_704  NR5A2_667   .... AR_689        least CREBs 
# 
# dat <- dat[order(dat$p_log_telo_norm),]
# head(dat[,c("gene","cor_cnt_telo","cor_log_telo","cor_cnt_log_telo_norm","p_log_telo_norm")],n=20) #most related
# tail(dat[,c("gene","cor_cnt_telo","cor_log_telo","cor_cnt_log_telo_norm","p_log_telo_norm")],n=20) #least related
# ### GRHL1_391   TP63_704  TP63_704

### TODO Maybe plot as a scatter, and point to interesting motifs with a pic?

#plot_motif(which(id2sym$id=="HNF4G_688"))

plotall <- function(usecor, use_p){
  
  nr <- nrow(dat)
  dat$one_cor <- dat[,usecor]
  dat$one_logp <- -log10(dat[,use_p])
  dat$similar_cor <- c(FALSE,abs(dat$one_cor[1:(nr-1)]-dat$one_cor[2:nr])/abs(dat$one_cor[2:nr])<0.1)
  dat$similar_p <- c(FALSE,abs(dat$one_logp[1:(nr-1)]-dat$one_logp[2:nr])/abs(dat$one_logp[2:nr])<0.1)
  dat$similar_both <- dat$similar_cor & dat$similar_p
  dat <- dat[!dat$similar_both,]

  dat$gene_part <- str_split_fixed(dat$gene,"_",2)[,1]
  plot(dat$one_cor, dat$one_logp,pch=19,cex=0, xlab="Corr.factor",ylab="-Log 10 p-value")
  text(dat$one_cor, dat$one_logp,labels = dat$gene_part,cex=0.5)
}
plotall("cor_cnt_telo","p_cnt_telo")
#MEF2D
#KLF1
#NAIF1
#ZHX1


plotall("cor_log_telo","p_log_telo")
pdf("zhang/corr/plots/volcano_log_telo.pdf",width = 5, height = 5)
plotall("cor_log_telo","p_log_telo")
dev.off()
#MEF2D*
#KLF1
#TP63
#GRHL1
#CEBPA/B/D

plotall("cor_cnt_log_telo_norm","p_log_telo_norm")
#AC0021266
#MEF2A  --- always relevant family
#NR4A1 vs AR?
#CEBP*



####################################################################
############### cool motifs ########################################
####################################################################




#########################
# Plot one motif
plot_motif <- function(i){
  m <- 2**as.matrix(pfm[[i]])
  for(i in 1:ncol(m)){
    m[,i] <- m[,i]/sum(m[,i])
  }
  seqLogo(m,xaxis = FALSE, yaxis = FALSE)
}

plot_motif_sym <- function(sym){
  print(sym)
  ID <-chromvarlist$motifSummary[sym,]$ID
  #print(ID)
  plot_motif(which(id2sym$id==ID))
  pdf(sprintf("zhang/corr/plots/motif_%s.pdf",sym), width = 3, height = 3)
  plot_motif(which(id2sym$id==ID))
  dev.off()
}


#id2sym$id


plot_motif_sym("ENSG00000250542_156")
plot_motif_sym("KLF9_192")
plot_motif_sym("FOSL2_105")
plot_motif_sym("JUNB_139")
plot_motif_sym("CGBP_298")    #Just CG

plot_motif_sym("OTX1_437")
plot_motif_sym("BACH2_113")
plot_motif_sym("TP63_704")  #distinct
plot_motif_sym("KDM2B_299")
plot_motif_sym("E2F1_312")
plot_motif_sym("CENPB_289")

plot_motif_sym("AR_689")
plot_motif_sym("HOXD11_464")
plot_motif_sym("GRHL1_391")
plot_motif_sym("ZBED1_816")
plot_motif_sym("IRF6_628")
plot_motif_sym("ENSG00000250542_156")
#plot_motif_sym("")
#plot_motif_sym("")
plot_motif_sym("CTCFL_198")
plot_motif_sym("KLF1_179")
plot_motif_sym("CTCF_177")
plot_motif_sym("SMARCC1_651")



plot_motif_sym("NR5A1_677")
plot_motif_sym("NR3C1_666")
plot_motif_sym("MEF2D_642")
#plot_motif_sym("")
#plot_motif_sym("")

#dat[str_starts(dat$gene,"NR5A1"),]
dat[str_starts(dat$gene,"MEF2D"),]


#oneid <- chromvarlist$motifSummary["",]$ID
#plot_motif(which(id2sym$id==oneid))
#pfm[[which(id2sym$id==oneid)]]



##########################
# Plot all motifs
for(i in 1:2){
  pdf(sprintf("motif_%s_%s",i,id2sym$symbol[i]))
  plot_motif(i)
  dev.off()
}


