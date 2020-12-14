# @Plots for MSc thesis cheetah MBD-seq  
# @Author: Tania P. Guerrero
# @October, 2019
# @IZW, Berlin. 

#This script contains the code used to plot most figures shown on my research outputs
#with some exceptions, such as the PCA that is found in the script medipReport.Rmd 

# All data will be loaded from resultsDirectory :) , so:
source("dataNames_AciJub.R")
setwd(resultsDirectory)

####packages 
library(ggplot2)
library(dplyr)
library(reshape2)
library(tidyr)
library(viridis)
library(clusterProfiler)
library(enrichplot)

# Plots from QC and data prep. --------------------------------------------
load(paste(resultsDirectory, "/qcLists.RData", sep=""))

## plotting saturation analysis
Sat_Plots <- for (i in satList[1:10]) {
  (MEDIPS.plotSaturation(i, main = paste("Saturation analysis WS 200", seqFiles$sampleName, sep="")))
}
###still need to fix the naming and how to store the output

# plot pattern coverage with pie chart
Cov_Plots_pie <- for (i in coverList[1:10]) {
  MEDIPS.plotSeqCoverage(seqCoverageObj= i, 
                         type = "pie", cov.level = c(0, 1, 2, 3, 4, 5))
}


# plot pattern coverage with histogram 
Cov_Plots_hist <- for (i in coverList[1:10]) {
  MEDIPS.plotSeqCoverage(seqCoverageObj= i, 
                         type = "hist", t= 30, main="Sequence pattern coverage")
}


#Calibration plot using the biggest scaffold and sample 1
png("./calibrationPlot.png", width = 1300, height = 900, pointsize = 20)
MEDIPS.plotCalibrationPlot(CSet = CS, main = "Calibration Plot", MSet =mset1[[1]], 
                           plot_chr = "NW_020834736.1", rpkm = T, xrange = T)
dev.off()

### Correlation matrix plots
#devtools::install_github("kassambara/ggcorrplot")
library(ggcorrplot)
library(ggplot2)
ggcorrplot(corMatrix)

##Plot corrMatrix
mCorMatrix <- melt(corMatrix, na.rm = T)

gc <- ggplot(data=mCorMatrix, aes(Var2, Var1, fill=value))+ 
  geom_tile(color="white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  coord_fixed()
corMat_plot <- gc + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, 
                                   hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        legend.justification = c(1, 0),
        legend.position = c(0.5, 0.8),
        legend.direction = "horizontal", legend.text = element_text(size = 12))+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))

ggsave(filename = "CorrMatrix.png",
       plot = corMat_plot,
       bg = "transparent",
       width = 25, height = 15, units = "cm", path = resultsDirectory)

ggsave(filename = "CorrMatrix_white.png",
       plot = corMat_plot,
       width = 25, height = 15, units = "cm", path = resultsDirectory)


# Plots from Differential methylation (methResults, etc. tables) ------------------------------------------------------

# Data wrangling to tailor the MA-plot  ---------------------------------------------
#This section contains basically a lot of manual data wrangling to add the top 10 ranked promoters
#in the MA-plot. But this code was actually quite useful for manually ranking tables in general
#which can be easily exported. Although, there're some sanity checks missing to confirm
#whether the coordinates (genomic ranges are correct)

load("ranked_promoters.Rdata") #which contain ranked promoters: ranked_p_g, ranked_p_l
#be careful if ncbi column might interfere ... so just drop it if necessary 

### Filtering top 10 promoters & creating the legend and labels
# The legends and labels were tailored very especifically to the coordinates. Therefore, this would need to 
#be manually adjusted for further samples/results (there's a ggplot package that simplifies this, see below!)

  #Hypermethylated promoters 
ranked_p_g %>% 
  filter(edgeR.logFC > 0, #I think this argument it's now kinda redundant but works as double-checkpoint
         gene != "") %>% 
  top_n(10, desc(edgeR.adj.p.value)) %>% 
  mutate(code = LETTERS[1:10]) -> tbl_pg_filtered1

tbl_pg_filtered1 %>%
  mutate(code = forcats::fct_collapse(code, "(CDEF)" = c("C", "D", "E", "F")
                               )) %>% 
  group_by(code) %>% slice(1) %>% data.frame() -> tbl_pg_filtered

tbl_pg_filtered1 %>% 
  dplyr::select(code, gene) %>%
  unite(code, gene, col = "legend", sep = " = ") %>% 
  pull(legend) -> gene_names_pg

  #Hypomethylated promoters 
ranked_p_l %>% 
  filter(edgeR.logFC < 0,
         gene != "") %>% 
  top_n(10, desc(edgeR.adj.p.value)) %>% 
  mutate(code = letters[1:10]) -> tbl_pl_filtered1

tbl_pl_filtered1 %>%
  mutate(code = forcats::fct_collapse(code, "(fg)" = c("f", "g")
  )) %>% 
  group_by(code) %>% slice(1) %>% data.frame() -> tbl_pl_filtered

tbl_pl_filtered1 %>% 
  dplyr::select(code, gene) %>%
  unite(code, gene, col = "legend", sep = " = ") %>% 
  pull(legend) -> gene_names_pl

# Behold the actual plot:
png("../results/plot-MA_Free-ranging_vs_Zoo_cheetahs_DER.png", width = 1300, height = 900, pointsize = 20)
par(las = 1, mar = c(5, 5, 1, 1))
smoothScatter(x=methResults$edgeR.logCPM,
              y=methResults$edgeR.logFC,
              xlim = c(-5, 1), ylim = c(-10, 6.5),
              xlab="log2(average counts per million)", 
              ylab="log2(Fold-change:Free-ranging/Zoo)")
legend("topright", c("p-val < 0.005", "adjusted p-val < 0.05"), pch = c(20,18), col = c("orange", "red"), bty = "n", title = "DMRs:")
points(x=p_val_relaxed$edgeR.logCPM, y= p_val_relaxed$edgeR.logFC, pch=".", col=ggplot2::alpha("orange", 0.7))
points(x=p_val_adj$edgeR.logCPM, y= p_val_adj$edgeR.logFC, pch=18, cex= 0.9, col=ggplot2::alpha("red", 0.5))
# points(x=p_val_adj_rel$edgeR.logCPM, y= p_val_adj_rel$edgeR.logFC, pch=".", cex= 0.5, col="blue")
abline(h=0, col="purple", lty = 2, lwd = 2)
text(x=tbl_pg_filtered$edgeR.logCPM, y= tbl_pg_filtered$edgeR.logFC, labels = tbl_pg_filtered$code, cex = 0.8, font = 2)
text(x=tbl_pl_filtered$edgeR.logCPM, y= tbl_pl_filtered$edgeR.logFC, labels = tbl_pl_filtered$code, cex = 1, font = 2)
legend("bottomleft", gene_names_pg, title = "Top 10 hypermethylated promoters:",  bty = "n")
legend("bottomright", gene_names_pl, title = "Top 10 hypomethylated promoters:",  bty = "n") # fill = to set coloured boxes next to legend text
dev.off()

#This version apparently does not work ... but the png looks good enough 
bitmap("../results/plot-MA_Free-ranging_vs_Zoo_cheetahs_DER.tiff", width = 1300, height = 900, type = "tifflzw", pointsize = 20, res = 300)
par(las = 1, mar = c(5, 5, 1, 1))
smoothScatter(x=methResults$edgeR.logCPM,
              y=methResults$edgeR.logFC,
              xlim = c(-5, 1), ylim = c(-10, 6.5),
              xlab="log2(average counts per million)", 
              ylab="log2(Fold-change:Free-ranging/Zoo)")
legend("topright", c("p-val < 0.005", "adjusted p-val < 0.05"), pch = c(20,18), col = c("orange", "red"), bty = "n", title = "Differentially methylated regions:")
points(x=p_val_relaxed$edgeR.logCPM, y= p_val_relaxed$edgeR.logFC, pch=".", col=ggplot2::alpha("orange", 0.7))
points(x=p_val_adj$edgeR.logCPM, y= p_val_adj$edgeR.logFC, pch=18, cex= 0.9, col=ggplot2::alpha("red", 0.5))
# points(x=p_val_adj_rel$edgeR.logCPM, y= p_val_adj_rel$edgeR.logFC, pch=".", cex= 0.5, col="blue")
abline(h=0, col="purple", lty = 2, lwd = 2)
text(x=tbl_pg_filtered$edgeR.logCPM, y= tbl_pg_filtered$edgeR.logFC, labels = tbl_pg_filtered$code, cex = 0.8, font = 2)
text(x=tbl_pl_filtered$edgeR.logCPM, y= tbl_pl_filtered$edgeR.logFC, labels = tbl_pl_filtered$code, cex = 1, font = 2)
legend("bottomleft", gene_names_pg, title = "Top 10 hypermethylated genes:",  bty = "n")
legend("bottomright", gene_names_pl, title = "Top 10 hypomethylated genes:",  bty = "n") # fill = to set coloured boxes next to legend text
dev.off()

### This is an alternative to MA base R plot, that might need less data wrangling to add legends, but I didn't have time to try it out
##This one needs the whole annotated table  
#library(ggpubr)
#ggmaplot()


# BarPlot with annotated genomic regions with MTC p-val 0.05  ---------------
load("/FuncAnn_tables.RData") #Which contain gain_ann & loss_ann tables 

gain05p <- table(gain_ann$LOCATION)
barplot(gain05p, main= "Fraction of differentially hypermethylated regions")

g05_splice <- unique(gain_ann$GENEID[gain_ann$LOCATION=="spliceSite"])
g05_intron <- unique(gain_ann$GENEID[gain_ann$LOCATION=="intron"])
g05_5UTR <- unique(gain_ann$GENEID[gain_ann$LOCATION=="fiveUTR"])
g05_3UTR <- unique(gain_ann$GENEID[gain_ann$LOCATION=="threeUTR"])
g05_coding <- unique(gain_ann$GENEID[gain_ann$LOCATION=="coding"])
g05_intergen <- unique(gain_ann$GENEID[gain_ann$LOCATION=="intergenic"])  ##mm I need R to actually count the NAs = 991 
g05_prom <- unique(gain_ann$GENEID[gain_ann$LOCATION=="promoter"])
length(g05_intergen)

loss05p <- table(loss_ann$LOCATION)
barplot(loss05p, main= "Fraction of differentially hypomethylated regions")

l05_splice <- unique(loss_ann$GENEID[loss_ann$LOCATION=="spliceSite"])
l05_intron <- unique(loss_ann$GENEID[loss_ann$LOCATION=="intron"])
l05_5UTR <- unique(loss_ann$GENEID[loss_ann$LOCATION=="fiveUTR"])
l05_3UTR <- unique(loss_ann$GENEID[loss_ann$LOCATION=="threeUTR"])
l05_coding <- unique(loss_ann$GENEID[loss_ann$LOCATION=="coding"])
l05_intergen <- unique(loss_ann$GENEID[loss_ann$LOCATION=="intergenic"]) ##mm I need R to actually count the NAs = 86
l05_prom <- unique(loss_ann$GENEID[loss_ann$LOCATION=="promoter"])
length(g05_intron)

regions <- c("spliceSite", "intron", "fiveUTR", "threeUTR", "coding", "intergenic", "promoter")
DMR_rel05=matrix(NA, 7,2, dimnames=list(regions, c("hypermethylated", "hypomethylated")))

DMR_rel05[1,1]= length(g05_splice) 
DMR_rel05[1,2]= length(l05_splice)

DMR_rel05[2,1]=length(g05_intron)
DMR_rel05[2,2]=length(l05_intron)

DMR_rel05[3,1]=length(g05_5UTR)
DMR_rel05[3,2]=length(l05_5UTR)

DMR_rel05[4,1]=length(g05_3UTR)
DMR_rel05[4,2]=length(l05_3UTR)

DMR_rel05[5,1]=length(g05_coding)
DMR_rel05[5,2]=length(l05_coding)

DMR_rel05[6,1]=991
DMR_rel05[6,2]=86

DMR_rel05[7,1]=length(g05_prom)
DMR_rel05[7,2]=length(l05_prom)

DMR_rel05

DMR_df05 <- as.data.frame(DMR_rel05)
rownames(DMR_df05) <- NULL
DMR_df05$Genomic_Region <- regions

write.csv(DMR_df05, row.names = F, quote = F, file="MeDIPS/results/DiffM_Genomic_Regions.csv")

# Code to the barplot
#base R 
png("./MeDIPS/results/DMRs_genomic_regions_05.png", width =  600, height= 600)
barplot(t(DMR_rel05),  main="Fraction of differentially methylated regions", beside=T, ylab="Fraction of 200 pb windows")
legend("topright", c("hypermethylated", "hypomethylated"), fill=grey.colors(2))
dev.off()

#Long format for ggplot
DMR_long05 <- gather(DMR_df05, key = "condition", value = count, hypomethylated:hypermethylated, factor_key = T)

#ggplot
gb <- ggplot(DMR_long05, aes(x=DMR_long05$Genomic_Region, y=DMR_long05$count, fill=factor(DMR_long05$condition))) +
  geom_bar(stat = "identity", width = 0.8, position = "dodge") + 
  xlab("Genomic region")+ylab("Fraction of annotated DMRs") + 
  labs(title= "Differentially methylated genomic regions [DMRs adjusted p-val <0.05]") +
  scale_fill_manual("", values = c("light blue", "blue"), labels= c("hypomethylated", "hypermethylated"))+ 
  scale_y_continuous(breaks = c(0, 100, 250, 500, 900, 1300))
#+ coord_cartesian(ylim = c(10:5000))
plotfile_gb <- gb + theme(
  # Remove panel border
  panel.border = element_blank(),  
  # Remove panel grid lines
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  # Remove panel background
  panel.background = element_blank(),
  # Add axis line
  axis.line = element_line(colour = "grey"),
  #Customize ticks text
  axis.text = element_text(face = "bold", size = 14),
  plot.title = element_text(face = "bold", size = 16), 
  axis.title = element_text(size = 14, color = "gray31"),
  legend.text = element_text(size = 14), 
  legend.position = "top")

ggsave(filename = "DMRs_genomic_regions05_a.png",
       plot = plotfile_gb,
       bg = "transparent",
       width = 25, height = 15, units = "cm", 
       path = resultsDirectory)

#transparent.plot <- gb + ggpubr::theme_transparent()

#ggsave(filename = "DMRs_genomic_regions05.png",
#  plot = transparent.plot,
#  bg = "transparent",
#   width = 10, height = 12, units = "cm", path = resultsDirectory)

#Diverging barchart
g3 <- ggplot(DMR_long05, aes(x=DMR_long05$Genomic_Region, y=DMR_long05$count)) + 
  geom_bar(stat='identity', aes(fill=DMR_long05$condition), width=.5)  +
  scale_fill_manual(name="", 
                    labels = c("Hypomethylated", "Hypermethylated"), 
                    values = c("hypomethylated"="light blue", "hypermethylated"="blue")) + 
  labs(title= "Fraction of differentially methylated regions [DMRs adjusted p-val <0.05]") + 
  xlab("Genomic region")+ylab("Fraction of 200 pb windows mapped to genes") +
  coord_flip() + 
  scale_y_continuous(breaks = c(0, 100, 250, 500, 1000, 3000))

plotfile_b <- g3 + theme(
  # Remove panel border
  panel.border = element_blank(),  
  # Remove panel grid lines
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  # Remove panel background
  panel.background = element_blank(),
  # Add axis line
  axis.line = element_line(colour = "grey"),
  #Customize ticks text
  axis.text = element_text(face = "bold", size = 14),
  plot.title = element_text(face = "bold", size = 16), 
  axis.title = element_text(size = 14, color = "gray31"),
  legend.text = element_text(size = 14), 
  legend.position = "top"
)
plotfile_b
ggsave(filename = "DMRs_genomic_regions_b.png",
       plot = plotfile_b,
       bg = "transparent",
       width = 25, height = 15, units = "cm", 
       path = resultsDirectory)


# Heatmap by p-val and meth status  ---------------------------------------
#This should be improved I don't really like it thaaat much 
load(paste(resultsDirectory, "/KeggEnrichResults.RData", sep="")) #which contains k_gp, k_lp, *prom_ranked_pval* & k_Allprom
#Here we only need the prom_ranked_pval 

# 1. Selecting the necessary info. 
prom_ranked_pval %>% dplyr::select(FC1 = edgeR.logFC, 
                                  edgeR.adj.p.value, gene) %>%  
                                  mutate(FC2=ifelse(FC1 <0, 0,1)) -> hm_prom

hm_prom <- as.data.frame(hm_prom, na.rm=T, stringsAsFactors=F)
hm_prom <- hm_prom[c(-127, -18),]
row.names(hm_prom) <- hm_prom$gene

#Just some checks 
#hm_prom1 <- scale(hm_prom[, c(-3, -4)])
#heatmap(hm_prom1, scale = "row")

class(hm_prom)

###ggplot tile code: 
#ggplot2::ggplot(hm_prom, ggplot2::aes(x=FC2, y = gene)) + ggplot2::geom_tile(ggplot2::aes(fill=edgeR.adj.p.value))
prom_methstatus <-  ggplot2::ggplot(hm_prom, ggplot2::aes(x = FC2, y=gene, fill=edgeR.adj.p.value)) + 
                    ggplot2::geom_tile() + 
                    scale_fill_viridis(name = "Adjusted p-value", option="C") + 
                    ggplot2::ylab("Putative gene promoters") +
                  ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                  panel.border = ggplot2::element_blank(),
                  panel.background = ggplot2::element_blank(),
                  axis.ticks = ggplot2::element_blank(), 
                  axis.title.x =ggplot2::element_blank(), 
                  axis.text.x= ggplot2::element_blank(), 
                  axis.text.y = ggplot2::element_text(vjust = 1, 
                                                     hjust = 1, size = 10, face ="bold", 
                                                     margin =ggplot2::margin(t = 1, b = 1)) 
                  )

ggplot2::ggsave(filename = "prom_methstatus.png",
                plot = prom_methstatus,
                bg = "transparent",
                width = 15, height = 50, units = "cm", path = resultsDirectory)

ggplot2::ggsave(filename = "prom_methstatus.tiff",
                plot = prom_methstatus,
                bg = "transparent",
                width = 15, height = 50, units = "cm", path = resultsDirectory)


# Plots from KEGG pathway enrichment (ClusterProfiler) --------------------
#If "/KeggEnrichResults.RData" has not been loaded (previous plot), please do so: Contains k_gp, k_lp, prom_ranked_pval & k_Allprom 
#Note that these plots are merely a guide to the possibilities since the best visualization will depend on the future results

## 1. Plots enrichment for significant/interesting results from separate lists of promoters 

#Quick way to explore which might work better 
dotplot(k_lp)
barplot(k_lp)
emapplot(k_lp) + ggtitle("KEGG Pathway enrichment map of hypomethylated promoters")

#clusterProfiler & enrichplot have a lot of cool options but it depends on the input we have
#e.g. cnetplot(), upsetplot(), etc. 

# Formatting and saving plots

tiff("./dotplotKEGGlp.tiff", width=3800, height=3000, res=400)
dotpp <- dotplot(k_lp, showCategory=7) + 
  ggtitle("KEGG Pathway enrichment for hypomethylated promoters [p-val<0.05]")
dotpp + ggplot2::theme(plot.title = ggplot2::element_text(face="bold", 
                                                          hjust = 0.5, color = "gray31"))
dev.off()

png("./barplotKEGGgp.png",width = 800, height =1000, pointsize = 30)
barplot(k_gp, color ="pvalue", showCategory = 5) + 
  ggtitle("KEGG Pathway enrichment for hypermethylated promoters [p-val<0.05]")
dev.off()

png("./enrichmapKEGGlp.png",width = 1200, height = 1000, pointsize = 30)
emapplot(k_lp, color = "pvalue") + 
  ggtitle("KEGG Pathway enrichment map of hypomethylated promoters")
dev.off()

png("./enrichmapKEGGlp.tiff", width=3600, height=3300, res=400)
emap_lp <- emapplot(k_lp, color = "pvalue") + 
  ggtitle("KEGG Pathway enrichment map of hypomethylated promoters")
emap_lp + ggplot2::theme(plot.title = ggplot2::element_text(face ="bold", 
                                                            hjust = 0.5, color = "gray31"))
dev.off()

png("./enrichmapKEGGgp.tiff", width = 4500, height= 4000, res=400)
emap_gp <- emapplot(k_gp, color = "pvalue") + 
  ggtitle("KEGG Pathway enrichment map of hypomethylated promoters")
emap_gp + ggplot2::theme(plot.title = ggplot2::element_text(face ="bold", 
                                                            hjust = 0.5, color = "gray31"))
dev.off()

## 2. Plots for combined methylation status (hypo & hyper)

#Creating an object with the logFC to fil in the foldChange argument of the heatplot 
geneList_prom <- prom_ranked_pval$edgeR.logFC
names(geneList_prom) <- prom_ranked_pval$ncbi

###thesis
tiff("./heat_enrich.tiff", width = 5000, height = 2500, res = 400)
heat_enrich <- heatplot(k_Allprom, foldChange = geneList_prom) #thesis 
heat_enrich + ggplot2::theme(axis.text.y = ggplot2::element_text(vjust= 1, hjust = 1,  size = 12, face ="bold"), 
                             axis.text.x = ggplot2::element_text(size = 12, face ="bold") ) + 
  ggplot2::scale_x_discrete(labels=c("ACTN1", "PIK3R1", "LOC106981787 IL-8 like", "ITGB2", "CLDN11", "CCR9", "CCL1")) 
dev.off()                            

heat_enrich <- heatplot(enrichK_prom, foldChange = geneList_prom) #thesis 
heat_enrich + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 15, face ="bold"), 
                             axis.text.x = ggplot2::element_text(size = 12, face ="bold") ) + 
  ggplot2::scale_x_discrete(labels=c("ACTN1", "PIK3R1", "LOC106981787 IL-8 like", "ITGB2", "CLDN11", "CCR9", "CCL1"))

# You've reached the end of this script :P ------
