library(ggplot2)


### To include in revision
### To include in revision
### To include in revision


################################################################################
######################## ML scores #############################################
################################################################################

subtelo_score <- read.table("/corgi/otherdataset/xingxi_atacsee/cellcycle/ml_pred.csv", header = FALSE, row.names = 1,sep=",")


subtelo_score2 <- data.frame(
  Run = str_split_fixed(rownames(subtelo_score),"_",2)[,1],
  subtelo_score1 = rowSums(subtelo_score[,1:25])/rowSums(subtelo_score),
  subtelo_score2 = rowSums(subtelo_score[,40:75])/rowSums(subtelo_score),
  subtelo_score3 = rowSums(subtelo_score[,76:100])/rowSums(subtelo_score)
)



meta <- read.csv("/corgi/otherdataset/xingxi_atacsee/cellcycle/samplemeta.csv",sep="\t")


toplot <- merge(meta, subtelo_score2)


toplot$renamed_ct <- str_remove_all(str_remove_all(toplot$Cell_type, "FACS sorted ")," group")
toplot$renamed_ct <- str_replace_all(toplot$renamed_ct, " low"," early")
toplot$renamed_ct <- str_replace_all(toplot$renamed_ct, " high"," late")

p1 <- ggplot(toplot, aes(renamed_ct, subtelo_score1)) + geom_point() + coord_flip() + xlab("Telo + subtelo") + ylab("Fraction reads")
p2 <- ggplot(toplot, aes(renamed_ct, subtelo_score2)) + geom_point() + coord_flip() + xlab("Chrom. proper") + ylab("Fraction reads")
p3 <- ggplot(toplot, aes(renamed_ct, subtelo_score3)) + geom_point() + coord_flip() + xlab("Centromere") + ylab("Fraction reads")
ptot <- egg::ggarrange(p1,p2,p3,ncol=1)
ggsave("xingxi_cc_ml.svg",ptot, width = 5, height = 4)


#subtelo_score2$subtelo_score1_n <-  subtelo_score2$subtelo_score1/subtelo_score2$subtelo_score2
#subtelo_score2$subtelo_score3_n <-  subtelo_score2$subtelo_score3/subtelo_score2$subtelo_score2
#p1 <- ggplot(toplot, aes(Cell_type, subtelo_score1_n)) + geom_point() + coord_flip()
#p3 <- ggplot(toplot, aes(Cell_type, subtelo_score3_n)) + geom_point() + coord_flip()
#egg::ggarrange(p1,p3,ncol=1)



################################################################################
#################### T cell time course ########################################
################################################################################


############ human 

library(ggplot2)

subtelo_score <- read.table("/corgi/otherdataset/tcelltimecourse/human_atac/ml_pred.csv", header = FALSE, row.names = 1,sep=",")

subtelo_score2 <- data.frame(
  ena = str_split_fixed(rownames(subtelo_score),"_",2)[,1],
  subtelo_score1 = rowSums(subtelo_score[,1:25])/rowSums(subtelo_score),
  subtelo_score2 = rowSums(subtelo_score[,40:75])/rowSums(subtelo_score),
  subtelo_score3 = rowSums(subtelo_score[,76:100])/rowSums(subtelo_score)
)

meta <- read.csv("/corgi/otherdataset/tcelltimecourse/samplemeta.csv",sep="\t")
meta$time <- str_split_fixed(meta$sampleid,"_",4)[,3]
meta$time[grep("Naive",meta$sampleid)] <- 0
meta$time <- as.integer(str_replace_all(meta$time,"h",""))

toplot <- merge(meta, subtelo_score2)

use_knots <- unique(meta$time)
use_knots <- c(0,4,24,72)
p1 <- ggplot(toplot, aes(time, subtelo_score1*1e6)) + geom_point() + 
  #geom_smooth(method="lm", formula=  y ~ splines::bs(x, knots = use_knots, degree = 2))  + 
  xlab("Human Th2,time (h)") + ylab("ML telo+subtelo (CPM)")
#p2 <- ggplot(toplot, aes(time, subtelo_score2)) + geom_point() + geom_smooth(method="lm", formula=  y ~ splines::bs(x, knots = use_knots, degree = 2)) #+ coord_flip()
#p3 <- ggplot(toplot, aes(time, subtelo_score3)) + geom_point() + geom_smooth(method="lm", formula=  y ~ splines::bs(x, knots = use_knots, degree = 2)) #+ coord_flip()
#egg::ggarrange(p1,p2,p3,ncol=1)



############ mouse 

library(ggplot2)

subtelo_score <- read.table("/corgi/otherdataset/tcelltimecourse/mouse_atac/ml_pred.csv", header = FALSE, row.names = 1,sep=",")

subtelo_score2 <- data.frame(
  ena = str_split_fixed(rownames(subtelo_score),"_",2)[,1],
  subtelo_score1 = rowSums(subtelo_score[,1:25])/rowSums(subtelo_score),
  subtelo_score2 = rowSums(subtelo_score[,40:75])/rowSums(subtelo_score),
  subtelo_score3 = rowSums(subtelo_score[,76:100])/rowSums(subtelo_score)
)

meta <- read.csv("/corgi/otherdataset/tcelltimecourse/samplemeta.csv",sep="\t")
meta$time <- str_split_fixed(meta$sampleid,"_",4)[,3]
meta$time[grep("Naive",meta$sampleid)] <- 0
meta$time <- as.integer(str_replace_all(meta$time,"h",""))

toplot <- merge(meta, subtelo_score2)

use_knots <- unique(meta$time)
use_knots <- c(0,4,24,72)
p2 <- ggplot(toplot, aes(time, subtelo_score1*1e6)) + geom_point() + 
  #geom_smooth(method="lm", formula=  y ~ splines::bs(x, knots = use_knots, degree = 2)) + 
  xlab("Mouse Th2,time (h)") + ylab("ML telo+subtelo (CPM)")



ptot <- p1|p2
ptot
ggsave(plot=ptot, "ml_tcell_tc.svg",width = 8, height = 3)

