###############################################################################
#
# A single-cell atlas of chromatin accessibility in the human genome
# https://www.cell.com/cell/fulltext/S0092-8674(21)01279-4?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867421012794%3Fshowall%3Dtrue
# https://www.cell.com/cell/pdf/S0092-8674(21)01279-4.pdf

# 1.3M cells
# 155 sci-atac datasets

###############################################################################
#
# https://data.mendeley.com/datasets/yv4fzv6cnm/1 
#  get 1B_Cell_metadata.tsv.gz
# from https://data.mendeley.com/public-files/datasets/yv4fzv6cnm/files/8aa0d83c-156b-4f30-b85b-160b7792f3b3/file_downloaded
#
###############################################################################


library(stringr)
library(sqldf)
library(ggplot2)

### Read processed telo counts + cell types (see zhang_summarizetelo.R)
info <- read.csv("zhang/summary_kmer.withct.csv.gz")
info$tissue <- str_replace(str_split_fixed(info$dataset,"-",2)[,1],"_SM","")

## For each sample, figure out average seq depth
avg_total <- sqldf("select dataset, avg(total) as avg_total from info group by dataset")
info <- merge(info, avg_total)

#info$norm_telo <- log(1+info$cnt_telo)
#info$norm_telo <- log(1+info$cnt_telo/info$avg_total*10000)
#info$norm_telo <- log(1+info$dedupcnt_CCCTAA)-log(info$avg_total)  #this even shifts 0!
info$cnt_telo <- info$dedupcnt_CCCTAA
info$log_telo <- log10(1+info$dedupcnt_CCCTAA)


#should compute rank within each batch!
#info$log_telo_norm <- log10(1+proj@cellColData$dedup_telo/proj@cellColData$nFrags*1e8)
info$rank_norm_telo<-NA
for(i in unique(info$dataset)){
  allrank <- rank(info$log_telo[info$dataset==i])
  info$rank_norm_telo[info$dataset==i] <- allrank/max(allrank)
}




######### Merge some types to clean up the plots and increase cell numbers
info$cell.type <- as.character(info$cell.type)
info$cell.type[str_starts(info$cell.type,"Colon Epithelial Cell")] <- "Colon Epithelial Cell"
info$cell.type[str_starts(info$cell.type,"Pericyte \\(General\\)")] <- "Pericyte (General)"
info$cell.type[str_starts(info$cell.type,"Pancreatic Beta Cell")] <- "Pancreatic Beta Cell"
info$cell.type[str_starts(info$cell.type,"Pancreatic Alpha Cell")] <- "Pancreatic Alpha Cell"
info$cell.type[str_starts(info$cell.type,"Mammary Luminal Epithelial Cell")] <- "Mammary Luminal Epithelial Cell"
info$cell.type[str_starts(info$cell.type,"Endothelial Cell \\(General\\)")] <- "Endothelial Cell (General)"
info$cell.type[str_starts(info$cell.type,"Cardiac Pericyte")] <- "Cardiac Pericyte"
info$cell.type[str_starts(info$cell.type,"Colon Epithelial Cell")] <- "Colon Epithelial Cell"
info$cell.type[str_starts(info$cell.type,"Keratinocyte")] <- "Keratinocyte"
info$cell.type[str_starts(info$cell.type,"Endothelial Cell \\(General\\)")] <- "Endothelial Cell (General)"
info$cell.type[str_starts(info$cell.type,"Smooth Muscle \\(Esophageal Muscularis\\)")] <- "Smooth Muscle (Esophageal Muscularis)"
info$cell.type[str_starts(info$cell.type,"Smooth Muscle \\(Colon\\)")] <- "Smooth Muscle (Colon)"
info$cell.type[str_starts(info$cell.type,"Vascular Smooth Muscle")] <- "Vascular Smooth Muscle"
# info$cell.type[str_starts(info$cell.type,"")] <- ""
#Alverolar Type 2,Immune   is this Alveolar Type 2 (AT2) Cell??


########### telomeres vs sex ########################
ggplot(info, aes(x=sex, y=log_telo)) + geom_violin()+ coord_flip()

table(info$sex)

#telo_reps <- sqldf("select dataset, avg(log_telo) as log_telo, count(log_telo) as cnt, `cell.type`, sex, tissue from info group by dataset,`cell.type`, sex order by log_telo")
telo_reps <- sqldf("select avg(log_telo) as log_telo, count(log_telo) as cnt, `cell.type`, sex, tissue from info group by tissue,`cell.type`, sex order by log_telo")
telo_reps <- telo_reps[telo_reps$cnt>100,]
telo_reps_male <- telo_reps[telo_reps$sex=="male",]
telo_reps_female <- telo_reps[telo_reps$sex=="female",]
telo_reps_male <- data.frame(
  log_telo_male=telo_reps_male$log_telo, 
  cell_type=telo_reps_male$cell.type, 
  tissue=telo_reps_male$tissue, 
  cnt_male=telo_reps_male$cnt)
telo_reps_female <- data.frame(
  log_telo_female=telo_reps_female$log_telo, 
  cell_type=telo_reps_female$cell.type, 
  tissue=telo_reps_female$tissue, 
  cnt_female=telo_reps_female$cnt)

telo_comp <- merge(telo_reps_female,telo_reps_male)
telo_comp$diff <- telo_comp$log_telo_female - telo_comp$log_telo_male
telo_comp <- telo_comp[order(telo_comp$diff),]
## more dimorphisms????

telo_reps_male[telo_reps_male$cell_type=="Chief Cell",]
telo_comp[telo_comp$cell_type=="Chief Cell",]


################################################################################################
########### telomeres vs cell type heatmap ######################################
################################################################################################

######### Some cells are not very abundant. Remove these
ct_count <- sqldf("select count(`cell.type`) as cnt, `cell.type` as ct from info group by `cell.type` order by cnt")
keep_ct <- ct_count$ct[ct_count$cnt>100]
info <- info[info$cell.type %in% keep_ct,]
ct_count

########### telomeres vs cell type ###########
avg_telo_ct <- sqldf("select `cell.type` as ct, avg(rank_norm_telo) as avg_telo from info group by `cell.type` order by avg_telo")

#Order cell types by averages
info$cell.type <- factor(info$cell.type, levels = avg_telo_ct$ct)
info

#ggplot(info, aes(x=cell.type, y=norm_telo)) + geom_violin()+ coord_flip()
ggplot(info, aes(x=cell.type, y=rank_norm_telo)) + geom_violin()+ coord_flip()

#info$rank_norm_telo
#ggsave("zhang/plot_violin_telo_ct.pdf", height = 20, width = 5)
#ggsave("zhang/plot_violin_telo_ct.png", height = 15, width = 5)




################################################################################################
##################### telomeres vs tissue ###################################################### new
################################################################################################

# 
# ### Comparison to RTL in healthy tissue
# demanelis <- read.csv("zhang/table_demanelis2020.csv",sep = "\t")
# demanelis$have_rtl <- demanelis$rtl_sd!=0
# demanelis <- demanelis[order(demanelis$rtl_mean),]
# demanelis$Tissue <- factor(demanelis$Tissue, levels = demanelis$Tissue)
# 
# p_rtl<- ggplot(demanelis, aes(x=Tissue, y=rtl_mean, fill=have_rtl)) +  #, 
#   geom_bar(stat="identity", color="black", 
#            position=position_dodge()) +
#   geom_errorbar(aes(ymin=rtl_mean-rtl_sd, ymax=rtl_mean+rtl_sd), width=.2,
#                 position=position_dodge(.9)) + coord_flip() 
#   #labs(title="Tooth length per dose", x="Dose (mg)", y = "Length")+
#   #theme_classic() 
#   #scale_fill_manual(values=c('#999999','#E69F00'))
# #print(p_rtl)
# 
# info$tissue <- factor(info$tissue, levels = demanelis$zhang)
# extras <- data.frame(
#   tissue=demanelis$zhang[demanelis$missing_zhang=="y"],
#   log_telo=rep(0,sum(demanelis$missing_zhang=="y"))  
# )
# info2 <- rbind(
#   info[,c("tissue","rank_log_telo")],extras,extras,extras
# )
# p_telo <- ggplot(info2, aes(x=tissue, y=rank_log_telo)) + 
#   geom_violin(width = 1.5)+ coord_flip() 
#   #theme_classic() 
# print(p_rtl)
# 
# library(gridExtra)
# ptot <- grid.arrange(p_rtl, p_telo, ncol=2, nrow =1)
# ggsave("zhang/plots/plot_rtl_vs_telo.pdf",width = 8,height = 5, plot=ptot)
# 

################################################################################################
##################### telomeres vs tissue ###################################################### old 
################################################################################################
# 
# 
# ########### telomeres vs tissue ###########
# avg_telo_tissue <- sqldf("select tissue, avg(log_telo) as avg_telo from info group by tissue order by avg_telo desc")
# info$tissue <- factor(info$tissue, levels = avg_telo_tissue$tissue)
# 
# avg_telo_tissue <- sqldf("select tissue, avg(norm_telo) as avg_telo from info group by tissue order by avg_telo desc")
# info$tissue <- factor(info$tissue, levels = avg_telo_tissue$tissue)
# 
# ggplot(info, aes(x=tissue, y=log_telo)) + 
#   geom_violin()+ coord_flip()
# ggsave("zhang/plot_violin_telo_tissue.pdf", height = 20, width = 5)
# ggsave("zhang/plot_violin_telo_tissue.png", height = 15, width = 5)
# 
# info$log_telo <- log10(1+info$cnt_telo)
# ggplot(info, aes(x=tissue, y=log_telo)) + 
#   geom_violin()+ coord_flip()  # unnormalized
# 
# unique(info$sample)
# unique(info$tissue)


################################################################################################
##################### adrenal glands? ###################################################### 
################################################################################################

infoa <- info[info$tissue=="adrenal_gland",]
infoa <- infoa[infoa$cell.type %in%c(
  "Endothelial (Exocrine Tissues)",
  "Macrophage (General)",
  "Fibroblast (Liver Adrenal)",
  "Cortical Epithelial-like",
  "Zona Glomerulosa Cortical Cell",
  "Zona Fasciculata Cortical Cell",
  "Transitional Zone Cortical Cell"),]
ggplot(infoa, aes(x=cell.type, y=rank_norm_telo)) + 
 geom_violin()+ coord_flip()
sqldf("select `cell.type`, avg(rank_norm_telo) as avg_telo from infoa group by `cell.type` order by avg_telo desc")
#Cortical Epithelial-like

#sort(table(info$cell.type[info$tissue=="adrenal_gland"]))

# 
# 1        Cortical Epithelial-like -6.869029
# 2  Zona Fasciculata Cortical Cell -6.893700
# 3 Transitional Zone Cortical Cell -6.954343

### . It is also a secondary site of androgen synthesis !!!!!




################################################################################################
###################################? ###################################################### old 
################################################################################################



########### telomeres vs cell type and tissue ###########
avg_telo_ct_tissue <- sqldf("select `cell.type` as ct, tissue, avg(rank_norm_telo) as avg_telo from info group by `cell.type`, tissue order by avg_telo")
avg_telo_ct_tissue <- avg_telo_ct_tissue[avg_telo_ct_tissue$avg_telo!=0,] #remove those with too few samples


sqldf("select ct, avg(avg_telo) as a from avg_telo_ct_tissue group by ct order by a")


ggplot(avg_telo_ct_tissue, aes(tissue,ct, fill= avg_telo)) + 
  geom_tile() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("zhang/plots/plot_heatmap_telo_ct_tissue.pdf", height = 15, width = 7)
#ggsave("zhang/plots/plot_heatmap_telo_ct_tissue.png", height = 20, width = 10)




unique(info$sample[info$Life.stage=="Adult"])

