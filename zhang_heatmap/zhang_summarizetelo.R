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


### Load kmer data from java telomemore
all_kmers <- read.csv("zhang/summary_kmer.java.csv.gz")

### Load mapping sample <-> annotation sample name
#map to sample names used in their metadata table
map_fname_sample <- read.csv("zhang/map_fname_sampleid.csv",stringsAsFactors = FALSE)
colnames(map_fname_sample) <- c("dataset","sample")
map_fname_sample$dataset <- str_replace(map_fname_sample$dataset,".demultiplexed.R1.fastq.gz.telo","")   ### Islet cells might not be treated well!

## Merge info
all_kmers <- merge(all_kmers, map_fname_sample)
all_kmers$dataset <- all_kmers$sample
all_kmers$cellID <- sprintf("%s+%s", all_kmers$sample, all_kmers$barcode)
head(all_kmers)

#colnames(sexinfo)
all_kmers$cellID <- sprintf("%s+%s", all_kmers$sample, all_kmers$barcode)


######### Same thing for donor info (infered)

sexinfo <- read.csv("zhang/sexinfo.csv.gz")[,-1]
sexinfo$dataset <- str_split_fixed(sexinfo$cb,"#",2)[,1]
sexinfo$barcode <- str_split_fixed(sexinfo$cb,"#",2)[,2]
sexinfo$sample <- str_replace(sexinfo$dataset,"_rep1","_1")   ### Islet cells might not be treated well!
sexinfo$cellID <- sprintf("%s+%s", sexinfo$sample, sexinfo$barcode)

#sexinfo$cellID %in% all_kmers$cellID



######### Read and integrate into the Zhang cell type annotation
info <- read.csv("zhang/1B_Cell_metadata.tsv.gz",sep="\t")
# only keep those we have kmers for (speeds up). only keeping 500k cells
info <- info[info$sample %in% map_fname_sample$sample,]

info <- merge(all_kmers,info)
info <- merge(info, sexinfo[,c("cellID","sex")])

## cleanup
info <- info[,c("dataset","barcode","dedupcnt_CCCTAA","totalcnt_CCCTAA","total","logUMI","cell.type","sex")]

## Save
write.csv(info,"zhang/summary_kmer.withct.csv")  


## later maybe. or use sample
#info$tissue <- str_replace_all(str_replace(str_split_fixed(info$tissue,"-",2)[,1],"_SM",""),"_"," ")



