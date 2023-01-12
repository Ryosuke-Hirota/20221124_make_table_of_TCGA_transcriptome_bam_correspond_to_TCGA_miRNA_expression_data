# this script is to correspond TCGA transcriptome bam to TCGA miRNA expression data 
# 2022/11/24 made

# activate package to get TCGA ID
library(stringr)
library(TCGAutils)
library(GenomicDataCommons)
library(magrittr)

# make new directory
setwd("C:/Rdata")
dir.create("20221124_correlation_analysis_using_TCGA_colon_transcriptome_bam_and_miRNA_quantification")

# import manifest of TCGA colon transcriptome bam
# this manifest is located at "https://www.dropbox.com/home/Okamura%20Lab%20share%20folder/Hirota/results_and_matterials/20221005_TCGA_salmon_quant_transcriptome"
setwd("C:/Rdata/20221005_TCGA_colon_salmon_quant_transcriptome")
transcript.list <-read.table("gdc_manifest.2022-10-18_only_TCGA_colon_transcriptome.txt",sep="\t",header = T,stringsAsFactors = F)
t.uuid <-transcript.list[,1]
t.uuid <-UUIDtoBarcode(t.uuid,from_type = "file_id")
t.uuid[,3] <-transcript.list[,2]
colnames(t.uuid)[2:3] <-c("aliquot_submitter_id","file_name")

# import manifest of TCGA miRNA quantification
# this manifest is located at "\\fsw-q02\okamura-lab\Files_related_to_M1_Projects\Hirota\TCGA_plot\20230105_TCGA_colon_miRNA_quantification"
miRNA.list <-read.table("gdc_manifest.2023-01-05_TCGA_colon_miRNA_quantification.txt",sep="\t",header = T,stringsAsFactors = F)
m.uuid <-miRNA.list[,1]
m.uuid <-UUIDtoBarcode(m.uuid,from_type = "file_id")
m.uuid[,3] <-miRNA.list[,2]
colnames(m.uuid)[2:3] <-c("aliquot_submitter_id","file_name")

# import table about TCGA transcriptome aliquot
# this table is located at "\\fsw-q02\okamura-lab\Files_related_to_M1_Projects\Hirota\TCGA_plot\20221005_TCGA_colon_salmon_quant_transcriptome"
t.aliquot <-read.table("20221005_TCGA_STARtwopass_transcriptome_aliquot.tsv",sep="\t",header = T,stringsAsFactors = F,fill = T,quote = "")
t.aliquot <-t.aliquot[,c(11,10,2,4)]

# search rows which have same barcode
t <-match(t.uuid[,2],t.aliquot[,1])

# remove unnecessary rows and add file id and file name
t.df <-t.aliquot[t,]
t.df[,c(5,6)] <-t.uuid[,c(1,3)]
colnames(t.df)[1:2] <-paste0("transcript_",colnames(t.df)[1:2])
colnames(t.df)[5:6] <-paste0("transcript_",colnames(t.df)[5:6])

# import table about TCGA miRNA aliquot
# this table is located at "\\fsw-q02\okamura-lab\Files_related_to_M1_Projects\Hirota\TCGA_plot\20221005_TCGA_colon_salmon_quant_transcriptome"
m.aliquot <-read.table("20221006_TCGA_miRNA_aliquot.tsv",sep="\t",header = T,stringsAsFactors = F,fill = T,quote = "")
m.aliquot <-m.aliquot[,c(11,10,2,4)]

# search rows which have same barcode
m <-match(m.uuid[,2],m.aliquot[,1])

# remove unnecessary rows and add file id and file name
m.df <-m.aliquot[m,]
m.df[,c(5,6)] <-m.uuid[,c(1,3)]
colnames(m.df)[1:2] <-paste0("miRNA_",colnames(m.df)[1:2])
colnames(m.df)[5:6] <-paste0("miRNA_",colnames(m.df)[5:6])

# merge aliquot of transcriptome and miRNA
tm.df <-merge(t.df,m.df,by=c("sample_id","case_id"))

cor.table <-tm.df[,c(5,6,3,9,10,7,1,2)]
rownames(cor.table) <-NULL

# extract rows which have duplicated file name
t.dup <-cor.table[duplicated(cor.table[,1])|duplicated(cor.table[,1],fromLast = T),]
m.dup <-cor.table[duplicated(cor.table[,4])|duplicated(cor.table[,4],fromLast = T),]

# search row number of duplicated file name 
t <-as.numeric(rownames(t.dup))
m <-as.numeric(rownames(m.dup))
tm <-sort(unique(append(t,m)))

# divide duplicated and not duplicated
dup.table <-cor.table[tm,]
cor.table <-cor.table[-tm,]

# set function to find simillar barcode 
find.sim.id <-function(x,y,z){
  library(stringr)
  correct.row <-NULL
  for (i in 1:nrow(z)) {
  x.id <-str_sub(x[i],start = -7,end = -5)
  y.id <-str_sub(y[i],start = -7,end = -5)
  if(x.id==y.id){
    correct.row <-append(correct.row,i)
  }}
  return(correct.row)
}

# search simillar barcode from duplicated file name and remain them
sim <-find.sim.id(dup.table[,3],dup.table[,6],dup.table)
no.dup.table <-dup.table[sim,]

# output final result
final.cor.table <-rbind(cor.table,no.dup.table)
setwd("C:/Rdata/20221124_correlation_analysis_using_TCGA_colon_transcriptome_bam_and_miRNA_quantification")
write.table(final.cor.table,"correspondence_table_between_TCGA_colon_transcriptome_bam_and_miRNA_qunatification_file.txt",sep = "\t",row.names = F,quote = F)
