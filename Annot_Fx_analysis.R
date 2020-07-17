# @Base and functional annotation for MSc thesis cheetah MBD-seq  
# @Author: Tania P. Guerrero
# @September, 2019
# @IZW, Berlin. 

# The first 3 chunks of this script prepare for and perform the basic and functional annotation 
#of the DMR - putative promotes obtained via MEDIPS package for differential methylation (edgeR). 
# The last chunk (d) runs one of the pathway enrichment tests: 
  #KEGG overrepresentation (hypergeometric test) enrichKEGG()
# NOTE: The BSgenome and org.db were also used for another steps not shown here.

# a) Forging BSgenome Ajubatus  ----------------------------------------------
# MANUAL https://www.bioconductor.org/packages/devel/bioc/vignettes/BSgenome/inst/doc/BSgenomeForge.pdf
### NOTE: If individual fasta files (individual scaffolds) are not yet available (as advised by J. Paijmans) 
#execute the following steps in the terminal:

  #First, generate a fai file for your reference (it contains information on scaffold names, lengths, etc, 
  #and is very handy to have with your reference)
  #samtools faidx AciJub_2_GCF_003709585.1.fasta

  #Then write a simple loop to extract the individual scaffolds from the reference, using your fai file:
  #for i in $(awk '{print $1}' AciJub_2_GCF_003709585.1.fasta.fai); do samtools faidx AciJub_2_GCF_003709585.1.fasta $i > AciJub_2_"$i".fasta; done;

  #mseqnames <-  paste("NW_0",seq(from="15130604",to="15144985",by=1),".1", sep = "")
  #mseqnames <-  paste("NW_0", c("015130604":"015144985"),".1", sep = "")

## 1. Using BSgenome forge according to the manual as follows: 
#BiocManager::install("BSgenome")
library(BSgenome)

#1.1 Create a seed :The format of this file is DCF (Debian Control File), 
#which is also the format used for the DESCRIPTION file of any R package.
#It's better to use a pre-existing BSgenome seed and just tailor it as needed 
#Seeds can be accessed as follows:
seed_files <- system.file("extdata",
                          "GentlemanLab", 
                          package="BSgenome")
tail(list.files(seed_files, pattern="-seed$"))

## E.g. Display seed file for musFur1 (or any other organism from the previous list):
musFur1_seed <- list.files(seed_files, 
                           pattern="\\.musFur1-seed$", 
                           full.names=TRUE) 
cat(readLines(musFur1_seed), sep="\n")

# The seed I used can be found inside the built package documentation 
# @ /data/fg2/guerrero/cheetah_methylation/MeDIPS/seed_Aci_Jub2

# 1.2 Run the forge package code 
  # Any prefix needed to regenerate filename from seqnames
seqfiles_prefix <- "AciJub_2_"
  # Any suffix needed to regenerate filename from seqnames
seqfiles_suffix <- ".fasta"
  #getSeqSrcpaths <- function(seqnames, prefix="AciJub_2_", suffix=".fasta", seqs_srcdir=seqs_srcdir)
seqs_srcdir = "/data/fg2/guerrero/refGenome"
  #trace("forgeBSgenomeDataPkg", browser, exit=browser) # for debugging only Alex. C. 
forgeBSgenomeDataPkg("BSgenome.AciJubGCF.NCBI.Acijub2-seed", 
                     seqs_srcdir = "/data/fg2/guerrero/refGenome", 
                     destdir = "/data/fg2/guerrero/cheetah_methylation/MeDIPS_test/seed_Aci_Jub2")
  #If for some reason this step fails, check the folder containing the individual scaffold files, 
#maybe the presence of the .gff file interferes (though it shouldn't)

## 2. After successfully forging it, go to terminal to the folder where seed file and forged is located and type
# R CMD build package_name
# R CMD check package_namecd .. 
# R CMD INSTALL package_name

# Sanity check 
library(BSgenome.AciJubGCF.NCBI.Acijub2)
seed_files <- "/data/fg2/guerrero/cheetah_methylation/MeDIPS/seed_Aci_Jub2/BSgenome.AciJubGCF.NCBI.Acijub2.Rcheck/BSgenome.AciJubGCF.NCBI.Acijub2"
tail(list.files(seed_files, pattern = "-seed$"))
# BSgenome should be up and running!!!  

# b) Forging non-model organism's annotation org.db packages: A. jubatus ----------------
#BiocManager::install("AnnotationForge")
library(AnnotationForge)

# Using makeOrgPackageFromNCBI #This step takes forever ... if it could be run in the shell
# in background mode would be best! But getting this package it's very useful to avoid 
#further complications on the downstream analysis 
#NOTE: getting this package does not substitute the genomic regions functional annotation  
# obtained via VariantAnnotaiton (next step "c") 

#https://www.ncbi.nlm.nih.gov/genome/?term=txid32536[Organism:noexp]

makeOrgPackageFromNCBI(version = "0.1",
                       author = "NCBI",
                       maintainer = "Tania Guerrero",
                       outputDir = ".",
                       tax_id = "32536",
                       genus = "Acinonyx",
                       species = "jubatus")

install.packages("./org.Ajubatus.eg.db", repos = NULL)

# c) Annotation: Functional genomic regions & genes   ---------------------------------------------
## 1. Packages 
# "https://www.bioconductor.org/packages/devel/workflows/vignettes/annotation/inst/doc/Annotating_Genomic_Ranges.html"
#BiocManager::install("MEDIPS")
library(MEDIPS)
#BiocManager::install("GenomicFeatures")
library(GenomicFeatures)
#BiocManager::install("VariantAnnotation")
library(VariantAnnotation)
#library(GenomicRanges) careful with overlapping fxs 
source("dataNames_AciJub.R")

## 2. gff_acijub needs to be TxDb or GRanges List to be used as argument in locateVariants() 
  #NOTE: Alternatively, could also use forged org.eg.db Ajubatus genome

  #So let's create a TxDb with GenomicRanges from gff file 
TxDb <- makeTxDbFromGFF(annotationGFF, 
                        dataSource = "gff annotation from cheetah assembly AciJub_2_Rico", 
                         organism = "Acinonyx jubatus", taxonomyId = 32536)

#  2.1. (optional) Importing annotation file as data frame. This one will not be used for annotating with VariantAnnotation 
  #BUT this object is might be useful if further steps (esp. assigning identifiers and other metadata) are required to be done manually. 
  # Which is why I decided to keep it here because this scripts contains all annotation-related wrangling.
#rawGFF <- readGFF(annotationGFF) 

## 3. Now, let's import results & filter out hypo and hypermethylated sites (0.05 MTC p-val)
      #***These annotation outputs are the ones used for further analysis: enrichment, plots, etc. 

#Getting the target subset from methResults  
  #loading raw results (if necessary) from differentially methylated regions (DMRs) obtained via MEDIPS 
  # (scripts medipAnalysisAdd_AciJub.R, medipProcessing_AciJub.R)

methResults <- load(paste(resultsDirectory, comparisonNames[an], "/methResults.RData", sep = ""))

#Selecting the MTC p-val we are interested in
p_val_0.05 <- MEDIPS.selectSig(methResults, p.value = 0.05, adj = T)

#Splitting tables by their methylation status 
  # gain = hypermethylation 
gain = p_val_0.05[which(p_val_0.05$edgeR.logFC > 0),]
  # loss = hypomethylation
loss = p_val_0.05[which(p_val_0.05$edgeR.logFC < 0),]

## 4. Formatting the data for the annotation
#Merging using MEDIPS.mergeFrames
gain_m = MEDIPS.mergeFrames(gain)
loss_m = MEDIPS.mergeFrames(loss)

head(gain_m)
head(loss_m)
tail(gain_m)
tail(loss_m)

gain_mg <- makeGRangesFromDataFrame(gain_m, keep.extra.columns = T)
loss_mg <- makeGRangesFromDataFrame(loss_m, keep.extra.columns = T)

TxDb2 <- TxDb #back up to intrtersect the seqlevels in next step 
              # (this needs to be adjusted per run, 
              # would be worth finding a way to avoid this step) 

# 5. Performing annotation: 
#Hypermethylated tables p-val 0.05
seqlevels(TxDb2) <- intersect(seqlevels(TxDb), seqlevels(gain_mg))
gain_ann <- locateVariants(gain_mg, TxDb2, AllVariants())
table(gain_ann$LOCATION)

#Hypomethylated tables p-val 0.05
seqlevels(TxDb2) <- intersect(seqlevels(TxDb), seqlevels(loss_mg))
loss_ann <- locateVariants(loss_mg, TxDb2, AllVariants())
table(loss_ann$LOCATION)

save(gain_ann, loss_ann, file = paste(resultsDirectory, "/FuncAnn_tables.RData", sep=""))

# d) Pathway enrichment analysis with ClusterProfiler ---------------------------
#Although this step can be performed with other genomic regions though I constrained it to promoters for the 
# better-studied functional outcome. I would not recommend this analysis for remaining enriched regions (UTRs, introns, intergenic, etc.)
#Those regions need another type of analysis if they are to be looked at. 

# Additional packages
#BiocManager::install("KEGGREST")
library(KEGGREST)
#BiocManager::install("clusterProfiler")
library(clusterProfiler)
library(readr)
library(dplyr)
library("org.Ajubatus.eg.db")

# 1. Data wrangling & sanitation to rank promoters by p-val
      ##--First for the hypermethylated promoters MTC p-val <0.05
p_gnames <- gain_ann$GENEID[gain_ann$LOCATION == "promoter"]
p_gchr <- gain_ann@seqinfo@seqnames[gain_ann$LOCATION == "promoter"] ##careful cuz this is a factor...

rangesp <- as.data.frame(gain_ann@ranges[gain_ann$LOCATION == "promoter"])

class(rangesp)
length(p_gnames)
length(p_gchr)
unique(p_gnames)

tbl_pg <- bind_cols(chr = p_gchr, 
                    start = rangesp$start, 
                    stop = rangesp$end, 
                    gene = p_gnames)
#removing duplicates
tbl_pg <- distinct(tbl_pg, gene, .keep_all = T)

#getting the rest of the info for ranking and further analysis
tbl_pg <- gain %>% left_join(tbl_pg, by = "start") ##removing the "stop" binding column so that they all ranges are preserved
unique(tbl_pg$gene)

gProm <- distinct(tbl_pg, gene, .keep_all = T) 
gProm %>% as_tibble() %>% 
  dplyr::select(chr = chr.x, start, 
                stop = stop.x, edgeR.logFC, 
                edgeR.logCPM, edgeR.p.value, 
                edgeR.adj.p.value, gene) %>% 
  arrange(edgeR.adj.p.value) -> ranked_p_g

        ##--Now for the hypomethylated promoters :) adj p-val <0.05
p_lnames <- loss_ann$GENEID[loss_ann$LOCATION == "promoter"]
p_lchr <- loss_ann@seqinfo@seqnames[loss_ann$LOCATION == "promoter"] 

rangesl <- as.data.frame(loss_ann@ranges[loss_ann$LOCATION == "promoter"])

class(rangesl)
length(p_lnames)
length(p_lchr)
dim(rangesl)

tbl_pl <- bind_cols(chr = p_lchr, 
                    start = rangesl$start, 
                    stop = rangesl$end,
                    gene = p_lnames)

#removing duplicates
tbl_pl <- distinct(tbl_pl, gene, .keep_all = T)

#getting the rest of the info for ranking and further analysis
tbl_pl <- loss %>% left_join(tbl_pl, by = "start")

lProm <- distinct(tbl_pl, gene, .keep_all = T)

lProm %>% as_tibble() %>% 
  dplyr::select(chr = chr.x, start, stop,
                edgeR.logFC, edgeR.logCPM, edgeR.p.value, 
                edgeR.adj.p.value, gene) %>% 
  arrange(edgeR.adj.p.value) -> ranked_p_l

# 2. (more data wrangling ¬¬) Adding NCBI ids to promoter lists 
#NOTE 1 !!!: The following steps use the org.eg.db package, if this package was not available
# for any reason, the wrangling in order to assign NCBI ids or any other identifier would
#be less straightforward in order to maintain the ranking. From my experience, I'd recommend biomaRt but since cheetah
#it's not available on biomart it'd need to use cat (which is what I used to get GO terms lists)
# but comes at the cost that some genes are dropped for lack of homologous b/w cat & cheetah.
#NOTE 2: apparently each time the org.eg.db is built, some categories change their names, so locating the 
# entrez ids might need some exploration, but the general steps should remain the same! :)
#THIS package should eventually be officially submitted to Bioconductor to avoid having to repeat
#the step every time! So if someone wants to help me! we should do it! 

#Locating the available attributes (e.g. identifiers) 
keys(org.Ajubatus.eg.db)
head(keys(org.Ajubatus.eg.db))

#Using the required identifier table to overlap it with our results 
z <- as.data.frame(org.Ajubatus.egALIAS2EG)
# Removing identifiers that do not map to any entrez gene id
z$alias_symbol <- z$alias_symbol[!is.na(z$alias_symbol)]

  ## For hypermethylated promoters:
#Mapping the NCBI id (gene id) to the gene symbols 
ranked_p_g$ncbi <- z$gene_id[match(ranked_p_g$gene, 
                                           z$alias_symbol) ]

  ## For hypomethylated promoters:
#Mapping the NCBI id (gene id) to the gene symbols 
ranked_p_l$ncbi <- z$gene_id[match(ranked_p_l$gene, 
                                   z$alias_symbol) ]

#Extracting NCBI ids for enrichKEGG()  
NCBI_gp <- c(as.character(ranked_p_g$ncbi))
NCBI_gp <- NCBI_gp[!is.na(NCBI_gp)]

NCBI_lp <- c(as.character(ranked_p_l$ncbi))
NCBI_lp <- NCBI_lp[!is.na(NCBI_lp)]

#Save ranked lists for downstream analysis and UniProt, GO, etc. 
save(ranked_p_g, ranked_p_l, file = paste(resultsDirectory, "/ranked_promoters.RData", sep=""))

# 3. Performing pathway enrichment analysis with enrichKEGG()
#This code shows the Over Representation Analysis (ORA) (Boyle et al. 2004), a widely used approach to 
# determine whether known biological functions or processes are over-represented (= enriched) in an experimentally-derived gene list, 
# e.g. a list of differentially expressed genes (DEGs).
# The p-value is calculated hypergeometric distribution.

#The clusterProfiler package offers more elaborated/robust options like proper GSEA, as well as combinations with GO, & other databases
#but with this type of data I think the ORA was the best suited, however we could explore more options.

#Some limitations of this enrichment tests still are: 1) depends on annotated modules and the ones available for cheetah are based on cat, 
#which in turn are based on human, dog, etc.; 2) it needs a solution to reduce the "noise" caused by non-tissue
# specific signals. For my thesis I manually filtered out all genes not recorded to be expressed (based on UniProt) in blood but 
#this wasn't enough and it's also not very feasible for long gene lists. 

# 3.1 Verifying supported organism 
search_kegg_organism("aju", by = "kegg_code")

#There are 173 modules for Cheetah 
modules_aju <- keggLink("module", "aju")
modules_aju <- unique(modules_aju)

#Input ID type can be kegg, ncbi-geneid, ncbi-proteinid or UniProt. I used ncbi-geneid
#There's an optional argument called "universe" which I think might help reducing the "noise" by providing 
#background genes. These background genes would probably need to be manually mined 
# to account for the fact that these are epigenetic and not expression data ... But not really sure where from, though. 

#NOTE FOR PLOTTING: pvalueCutoff and qvalueCutoff can be tailored ONLY for visualization purposes
#(e.g. setting a less astringent threshold so that the plots show more results) 

# 3.2 Performing the analysis for separate lists of promoters 

#Hypermethylated putative promoters
k_gp <- enrichKEGG(NCBI_gp, 
                      organism = "aju", 
                      keyType = "kegg", 
                      pvalueCutoff = 0.05, 
                      pAdjustMethod = "BH", 
                      minGSSize = 10, 
                      maxGSSize = 100,
                      qvalueCutoff = 0.2, #This was just the default & I'd like to read more about these setting in general
                      use_internal_data = F)

#Hypomethylated putative promoters 
k_lp <- enrichKEGG(NCBI_lp, 
                      organism = "aju", 
                      keyType = "kegg", 
                      pvalueCutoff = 0.05, 
                      pAdjustMethod = "BH",
                      minGSSize = 10, 
                      maxGSSize = 100,
                      qvalueCutoff = 0.2, 
                      use_internal_data = F)


# 3.3 For all promoters combining methylation status (hypo & hyper)
#These steps need to be adjusted according to results but the general idea is merging hypo & hypermethylated results 
#using the ranked tables and the identifier tables (ncbi id or other identifier needed) to provide the necessary information 
#for the analysis and for plotting 

prom <- rbind(ranked_p_l, ranked_p_g)

prom %>% dplyr::select(edgeR.logFC, 
                       edgeR.adj.p.value, gene, ncbi) %>% 
  dplyr::arrange(edgeR.adj.p.value) -> prom_ranked_pval

gl_prom <- prom_ranked_pval$ncbi[!is.na(prom_ranked_pval$ncbi)]

k_Allprom <- enrichKEGG(gl_prom, organism = "aju", 
                           keyType = "kegg", 
                           pvalueCutoff = 0.05, 
                           pAdjustMethod = "fdr",
                           minGSSize = 10, 
                           maxGSSize = 120,
                           qvalueCutoff = 0.2, 
                           use_internal_data = F)

save(k_gp, k_lp, prom_ranked_pval, k_Allprom, file = paste(resultsDirectory, "/KeggEnrichResults.RData", sep=""))

# You've reached the end of this script :P ------
