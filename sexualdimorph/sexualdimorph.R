library(ggplot2)
library(stringr)

# 54 men, 66 women
# about 2.5 reps/sample
# (54+66)*2.5

telocnt <- read.table("telo.csv",sep=",")
colnames(telocnt) <- c("srr","teloc","totc")
telocnt$srr <- str_replace_all(telocnt$srr,"_1.fastq","")

runtable <- read.csv("SraRunTable.csv")
runtable <- data.frame(srr=runtable$Run, SUBJECT_ID=runtable$submitted_subject_id, libname=runtable$Library.Name)

length(unique(runtable$SUBJECT_ID))  #54+66=120. in here, 121

dbgap <- read.csv("samplemeta.csv",sep="\t")
dbgap <- dbgap[,c("SUBJECT_ID","SEX","AGE","BMI")]

dat <- merge(merge(runtable,telocnt),dbgap)
dat$AGE <- str_replace_all(dat$AGE,"\\+","")
dat$AGE <- as.double(dat$AGE)
dat$frac_telo <- dat$teloc / dat$totc
dat$logfrac_telo <- log10((1+dat$teloc) / dat$totc)


################# read quality of genomes! #####################
flist <- list.files("nucleosomesignal/")
fragsize <- matrix(nrow=500, ncol=length(flist))#NULL
colnames(fragsize) <- str_split_fixed(flist,"_",2)[,1]
for(i in seq_along(flist)){
  fcontent <- read.table(sprintf("nucleosomesignal/%s", flist[i]),sep="\t")#, col.names = FALSE)
  head(fcontent)
  fragsize[,i] <- fcontent$V2/sum(fcontent$V2)
}

plot(fragsize[,1],type="l",ylim=c(0,0.02))
for(i in seq_along(flist)){
  lines(fragsize[,i])
}

nucsig <- data.frame(
  srr=colnames(fragsize),
  nucsic=colSums(fragsize[147:294,])/colSums(fragsize[1:146,])
)


dat <- merge(dat,nucsig)



################## plotting ###############################



plot(dat$frac_telo, dat$nucsic)
plot(log(dat$frac_telo), dat$nucsic)

dat2 <- dat[dat$nucsic<1 & dat$frac_telo < 0.5e-4 & dat$logfrac_telo > -6,]
plot(dat2$frac_telo, dat2$nucsic)

### Maybe don't trust those with really low nucsig.
ggplot(dat[dat$nucsic<1 & dat$frac_telo < 0.5e-4 & dat$logfrac_telo > -6,], aes(AGE,logfrac_telo)) + 
#ggplot(dat[dat$nucsic<1 & log(dat$frac_telo) < -9,], aes(AGE,logfrac_telo)) + 
  geom_point() + 
  #  geom_text() +
  geom_smooth(method = "lm", se = FALSE)

ggplot(dat[dat$nucsic>1,], aes(AGE,logfrac_telo)) + 
  geom_point() + 
  #  geom_text() +
  geom_smooth(method = "lm", se = FALSE)




#plot(dat$AGE,dat$frac_telo)
#plot(dat$AGE,log10((1+dat$teloc) / dat$totc))


keep <- dat$nucsic<1 & dat$frac_telo < 0.5e-4 & dat$logfrac_telo > -6
cor.test(as.double(dat$AGE[keep]),log10((1+dat$teloc[keep]) / dat$totc[keep]))

cor.test(as.double(dat$AGE),log10((1+dat$teloc) / dat$totc))

#### Currently in paper
ggplot(dat[dat$logfrac_telo> -6,], aes(AGE,logfrac_telo)) + 
  geom_point() + 
  xlab("Age") +
  ylab("Log10 1+nTA") +
  geom_smooth(method = "lm", se = FALSE)
ggsave("nTA vs age filtered.pdf", height = 2)

cor.test(dat[dat$logfrac_telo> -6,]$logfrac_telo, dat[dat$logfrac_telo> -6,]$AGE)

# line.lm <- lm(logfrac_telo ~ AGE, dat)
# ggplot(dat, aes(AGE,logfrac_telo)) + geom_point(color="blue") + geom_abline(slope = coef(line.lm)[["AGE"]], 
#             intercept = coef(line.lm)[["(Intercept)"]])



########### telomere abundance ##############

sum(telocnt$teloc)/sum(telocnt$totc)









