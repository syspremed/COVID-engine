
# Sadanandam et al., A blood transcriptome-based system-level analysis of disease progression, immune regulation and symptoms in coronavirus-infected patients
# ResearchSquare - published preprint - DOI:10.21203/rs.3.rs-30473/v1
# NTP analysis for acute vs. recovering patients from Lee et al. 

setwd("/Users/asadanandam/Dropbox (SPM)/SysPreMed/Anguraj/COVID/manuscript/script/test/NTP_covid")
source("/Users/asadanandam/Dropbox (SPM)/Anguraj laptop/laptop-backup/asadanandam/software/NTPez.R")
NTPez("20200415_GSE1739_RMA_cell_genes_only_corrected_names_sd0.gct","/Users/asadanandam/Dropbox (SPM)/SysPreMed/Anguraj/COVID/GSE1739/again_new/additional_dataset/Lee_NTP_marker_only_39_of_52.txt")

#reading the ntp result file
ntp <- read.delim("NTP_prediction_result.xls", stringsAsFactors = FALSE)

dim(ntp)

ntp_d <- data.frame(ntp[1:10,])
dim(ntp_d)

ntp_w <- which(ntp_d$predict.label == 1)
ntp_d$dist.to.template[ntp_w] = -1*ntp_d$dist.to.template[ntp_w] 

library(ggplot2)

#pdf("2020417_acute_vs_recovering_NTP_lolipop_plot_GSE1739.pdf", useDingbats = FALSE)
theme_set(theme_bw())

# this plot was modified as in Figure 1G in Adobe Illustrator
ggplot(ntp_d, aes(x=reorder(sample.names, dist.to.template), y=dist.to.template, label=predict.label)) + 
  geom_point(stat='identity', color="darkgreen", size=8)  +
  geom_segment(aes(y = 0, 
                   x = sample.names, 
                   yend = dist.to.template, 
                   xend = sample.names), 
               color = "gray", size = 2) +
  geom_text(color="white", size=2,lineend = "butt") +
  labs(title="Diverging Lollipop Chart", 
       subtitle="Distance to acute vs recovering patients: Lollipop") + 
  ylim(-1.3, 1.3) +
  coord_flip()

#dev.off()
