# Sadanandam et al., A blood transcriptome-based system-level analysis of disease progression, immune regulation and symptoms in coronavirus-infected patients
# ResearchSquare - published preprint - DOI:10.21203/rs.3.rs-30473/v1
# Enrichment analysis for CoV gene signatures

#names(caf) <- "Down genes"
library(hypeR)

library(qusage)

KEGG <- msigdb_download("Homo sapiens", "C2", "CP:KEGG")
REACTOME <- msigdb_download("Homo sapiens", "C2", "CP:REACTOME")
H <-msigdb_download("Homo sapiens", "H", "")
C7 <-msigdb_download("Homo sapiens", "C7", "")

# gmt files were downloaded from Enrichr to local drive
subcell <- read.gmt("/Users/asadanandam/Dropbox (SPM)/SysPreMed/Anguraj/Analysis/Internal/gene_sets/Jensen_COMPARTMENTS.gmt")
biogps <- read.gmt("/Users/asadanandam/Dropbox (SPM)/SysPreMed/Anguraj/COVID/genesets/BioGPS_Human_cell_type/biogps_gene_set_library_up_crisp.gmt")
biogps_down <- read.gmt("/Users/asadanandam/Dropbox (SPM)/SysPreMed/Anguraj/COVID/genesets/BioGPS_Human_cell_type/biogps_gene_set_library_dn_crisp.gmt")

#Hallmarks
hyp_H_d <- hypeR(caf$`Down genes`,H, test="hypergeometric")
hyp_show(hyp_H_d)
write.table(hyp_H_d$data,"20200531_Hallmarks_hypeR_data.txt", row.names = FALSE, sep = "\t", quote = FALSE)
pdf(file="Hallmarks_plot.pdf", useDingbats = FALSE)
hyp_dots(hyp_H_d, fdr = 0.2)
dev.off()
hyp_emap(hyp_H_d, top=30, val="fdr", fdr = 1)


#MsigDB C7
hyp_C7_d <- hypeR(caf$`Down genes`, C7, test="hypergeometric")
hyp_show(hyp_C7_d)
write.table(hyp_C7_d$data,"20200531_C7_hypeR_data.txt", row.names = FALSE, sep = "\t", quote = FALSE)

pdf(file="C7_immune_genes_plot.pdf")
hyp_dots(hyp_C7_d, fdr = 0.2)
dev.off()
pdf(file="C7_immune_genes_plot_network.pdf")
hyp_emap(hyp_C7_d, top=30, val="fdr")
dev.off()

#REACTOME
hyp_REACTOME_d <- hypeR(caf$`Down genes`, REACTOME, test="hypergeometric")
hyp_show(hyp_REACTOME_d)
write.table(hyp_REACTOME_d$data,"20200531_REACTOME_hypeR_data.txt", row.names = FALSE, sep = "\t", quote = FALSE)
pdf(file="REACTOME_plot.pdf", useDingbats = FALSE)
hyp_dots(hyp_REACTOME_d, fdr = 0.2)
dev.off()
hyp_emap(hyp_REACTOME_d, top=15, val="fdr", similarity_metric = "overlap_similarity", similarity_cutoff = 0.5)

#KEGG
hyp_KEGG_d <- hypeR(caf$`Down genes`,KEGG, test="hypergeometric")
hyp_show(hyp_KEGG_d)
write.table(hyp_KEGG_d$data,"20200531_KEGG_hypeR_data.txt", row.names = FALSE, sep = "\t", quote = FALSE)
pdf(file="KEGG_plot.pdf")
hyp_dots(hyp_KEGG_d, fdr = 0.2)
dev.off()
hyp_emap(hyp_KEGG_d, top=40, val="fdr", fdr=0.2, similarity_metric =  "overlap_similarity", similarity_cutoff = 0.25)

#bioGPS
hyp_biogps_d <- hypeR(caf$`Patients Genes`,biogps, test="hypergeometric")
hyp_show(hyp_biogps_d)
write.table(hyp_biogps_d$data,"20200531_biogps_hypeR_data.txt", row.names = FALSE, sep = "\t", quote = FALSE)

pdf(file="biogps_up_plot_p_0.00001.pdf", useDingbats = FALSE)
hyp_dots(hyp_biogps_d, fdr = 0.00001,top=50)
dev.off()
hyp_emap(hyp_biogps_d, fdr=0.2, val="fdr")

hyp_biogps_u <- hypeR(caf$`Healthy Genes`,biogps, test="hypergeometric")
hyp_show(hyp_biogps_u)
write.table(hyp_biogps_u$data,"20200531_biogps_DOWN_hypeR_data.txt", row.names = FALSE, sep = "\t", quote = FALSE)

pdf(file="biogps_up_healthy_plot_0.0000001.pdf", useDingbats = FALSE)
hyp_dots(hyp_biogps_u, fdr = 0.0000001, top=50)
dev.off()
hyp_emap(hyp_biogps_d, top=100, val="fdr")

