
# Sadanandam et al., A blood transcriptome-based system-level analysis of disease progression, immune regulation and symptoms in coronavirus-infected patients
# ResearchSquare - published preprint - DOI:10.21203/rs.3.rs-30473/v1
# ssGSEA analysis for CoV-Up-gene and CoV-Down-gene signatures 
# One dataset from Xiong et al. has been shown. Similarly other datasets can be applied
# this code applies to immune cell type analysis in the manuscript

dir <- "/Users/asadanandam/Dropbox (SPM)/SysPreMed/Anguraj/COVID/manuscript/script/test/ssGSEA_cov_sig"
gene_set <- "/Users/asadanandam/Dropbox (SPM)/SysPreMed/Anguraj/COVID/manuscript/Data/Xiaong/20200426_GSE1739_patient_healthy_SAM_genes.gmt"
#prepare a .gct file as per GenePattern - https://www.genepattern.org/file-formats-guide
input <- "/Users/asadanandam/Dropbox (SPM)/SysPreMed/Anguraj/COVID/manuscript/Data/Xiaong/20200426_Xiaong_supplementary_PBMC.gct"
output <- "20200426_Xiaong_supplementary_PBMC_Patient_ssGSEA_results.gct"

setwd(dir)
  
#ssGSEA - make sure you install ssGSEA code from GenePattern
source("/Users/asadanandam/Dropbox (SPM)/Anguraj laptop/laptop-backup/asadanandam/software/ssGSEA/common.R")
source("/Users/asadanandam/Dropbox (SPM)/Anguraj laptop/laptop-backup/asadanandam/software/ssGSEA/ssGSEAProjection.Library.R")
#source("/Users/asadanandam/Dropbox (SPM)/Anguraj laptop/laptop-backup/asadanandam/software/ssGSEA/ssGSEAProjection.R")

ssGSEA.project.dataset(input.ds = input, output.ds = output, gene.sets.dbfile.list=gene_set)


#remove first two rows manually before reading the file
ssGSEA_RNAseq <- read.delim("20200426_Xiaong_supplementary_PBMC_Patient_ssGSEA_results.gct", stringsAsFactors = FALSE)
dim(ssGSEA_RNAseq)

ssGSEA_RNAseq_m <- data.matrix(ssGSEA_RNAseq[,3:ncol(ssGSEA_RNAseq)])
dim(ssGSEA_RNAseq_m)
rownames(ssGSEA_RNAseq_m) <- ssGSEA_RNAseq[,1]

cl_xiaong <- rep(c("N","P"),c(3,3))

boxplot(ssGSEA_RNAseq_m[1,]~cl_xiaong)
kruskal.test(ssGSEA_RNAseq_m[1,]~cl_xiaong)

ssGSEA_RNAseq_xiang_m_d <- data.frame(unlist(ssGSEA_RNAseq_m[1,]),cl_xiaong)
names(ssGSEA_RNAseq_xiang_m_d) <- c("data","class")

library(ggplot2)

pdf(file= "20200426_Xiaong_supplementary_PBMC_Patient_ssGSEA_results_violin.pdf")

ggplot(data = ssGSEA_RNAseq_xiang_m_d, aes(x=reorder(as.factor(class), -as.numeric(as.character(data))),y=as.numeric(as.character(data)))) + geom_violin(trim=FALSE,fill="gray" )+ geom_jitter(shape=16, position=position_jitter(0.2)) + geom_boxplot(width=0.1, fill = c("#999999", "#E69F00")) 
dev.off()

kruskal.test(ssGSEA_RNAseq_m[1,]~as.factor(cl_xiaong))

ssGSEA_RNAseq_xiang_m_d_n <- data.frame(unlist(ssGSEA_RNAseq_m[2,]),cl_xiaong)
names(ssGSEA_RNAseq_xiang_m_d_n) <- c("data","class")
pdf("20200426_Xiaong_supplementary_PBMC_Patient_ssGSEA_results_violin_normal.pdf", useDingbats = FALSE)

ggplot(data = ssGSEA_RNAseq_xiang_m_d_n, aes(x=as.factor(class),y=as.numeric(as.character(data)))) + geom_violin(trim=FALSE,fill="gray" )+ geom_jitter(shape=16, position=position_jitter(0.2)) + geom_boxplot(width=0.1, fill = c("#999999", "#E69F00")) 
dev.off()

