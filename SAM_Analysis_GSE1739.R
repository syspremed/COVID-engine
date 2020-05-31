setwd("/Users/asadanandam/Dropbox (SPM)/SysPreMed/Anguraj/COVID/GSE1739/again_new/another_new")
setwd("/Users/asadanandam/Dropbox (SPM)/SysPreMed/Anguraj/COVID/GSE1739/again_new/another_new")



sam_COVID <- function(dir)
{

setwd(dir)
  
# Reading the CEL files from Affymetrix for GSE1739  
library(affy)

data <- ReadAffy()
data_rma <- rma(data)
data_exp <- exprs(data_rma)
dim(data_exp)
#head(data_exp[,1:5])


# Reading phenotypic information for GSE1739
library(GEOquery)

gse <- getGEO("GSE1739")

length(gse)

gse_1 <- gse[[1]]

gse_1_fdata <- fData(gse_1)
dim(gse_1_fdata)

#head(gse_1_fdata[,10:16])

#head(gse_1_fdata$`Gene Symbol`)

m <- match(rownames(data_exp), gse_1_fdata[,1] )
w <- which(!is.na(m))
length(w)

data_exp_1 <- cbind(gse_1_fdata$`Gene Symbol`[m[w]],data_exp)
dim(data_exp_1)
colnames(data_exp_1)[1] <- "Genes"
head(data_exp_1[,1:5])


write.table(data_exp_1, "20200426_GSE1739_data_genes_only.txt", row.names = FALSE, sep = "\t", quote = FALSE)

# parsing sample IDs
str <- NULL

for(i in 1:nrow(data_exp_1)){
  str[i] <- unlist(strsplit(data_exp_1[i,1]," "))[1]
}

length(str)

head(str)

data_exp_2 <- data_exp_1

data_exp_2[,1] <- str
dim(data_exp_2)  

head(data_exp_2[,1:5]) 

write.table(data_exp_2, "20200426_GSE1739_data_genes_only_corrected.txt", row.names = FALSE, sep = "\t", quote = FALSE)

# selecting variable genes with SD > 0.8
source("/Users/asadanandam/Dropbox (SPM)/Anguraj laptop/laptop-backup/asadanandam/software/screenExpr.r")
screenExpr("20200426_GSE1739_data_genes_only_corrected.txt", 0.8)

data_sd0.8 <- read.delim("20200426_GSE1739_data_genes_only_corrected_sd0.8.txt", stringsAsFactors = FALSE, row.names = 1)
dim(data_sd0.8)

cl <- rep(c(0,1),c(10,4))

#SAM analysis
library(siggenes)

sam.out <- sam(data_sd0.8,cl,rand=123,gene.names = rownames(data_sd0.8))
sam.out

pdf("sam_delta.pdf")
plot(sam.out)
dev.off()
pdf("sam_fold.pdf")
plot(sam.out,1.1)
dev.off()

genes <- summary(sam.out, 1.1)@row.sig.genes

write.table(summary(sam.out,1.1)@mat.sig,"20200426_samout.1.1_290genes.txt",  sep = "\t", quote = FALSE)

signi <- summary(sam.out,1.1)@mat.sig
dim(signi)

patient <- which(signi$d.value < 0)
length(patient)

healthy <- which(signi$d.value >= 0)
length(healthy)

data_sd0.8_sam <- data_sd0.8[genes,]
dim(data_sd0.8_sam)

write.table(data_sd0.8_sam, "20200426_GSE1739_data_genes_only_corrected_SAM_genes.txt", sep = "\t", quote = FALSE)
}
