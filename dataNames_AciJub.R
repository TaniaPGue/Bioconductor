## Created 6/29/2015 by Daniel Beck
## Last modified 5/31/2016
#Modified by Johana P. 20/03/2019
#Modified by Tania G. 25/03/2019
#Modified by Tania G. Oct/2019 (these settings were the ones used for my MSc data analysis)
#Please note that at some stages is necessary to go back and "switch-on/off" 
#some options or modified directories (particularly for the file generateReports.R and medipReport.Rmd)
#maybe that wouldn't be necessary if directories hadn't become a mess but just a general fyi: FILE NAMES MATTER, NAMES GIVEN IN THIS 
# FILE, MATTER! a lot of trouble-shooting has to do with just fixing some paths and variable names!
#very importantly: THESE SETTINGS REQUIRE LOGIC WHEN ASSIGNING THE VECTORS AND VARIABLES. E.G. IF THERE ARE 4 SAMPLES INSTEAD OF 5, ADJUST 
#ACCORDINGLY!

## This is the configuration script for the MeDIP-seq analysis pipeline. It holds
## most analysis parameters and options. It also holds sample/filename information
## and defines which analyses are performed. It also performs some preliminary 
## checks to ensure the analysis can procede. 

################################
## General project attributes ##
################################
## This section sets the folder structure. The default is to have a single project
## folder with three main subfolders for code,f data, and results. An additional 
## genome directory holds the reference genome and is located within the data folder.
## It should be fine to modify this structure, however, extensive testing of any
## alternatives has not been done.

projectName <- "cheetah_methylation"

# Project directory (main folder)
projectDirectory <- "/data/fg2/guerrero/cheetah_methylation/MeDIPS/"
# Data directory (holds data files)
dataDirectory <- paste(projectDirectory, "data/", sep="")
# Code directory (holds all code used for the MeDIP-seq analysis)
codeDirectory <- paste(projectDirectory, "master_code/", sep="")
# Results directory (folder for all results)
resultsDirectory <- paste(projectDirectory, "results/", sep="") 
# Genome directory (holds reference genome and BSgenome package)
genomeDirectory <- paste(dataDirectory, "genome/refGenome/", sep="")

# Trimmed sequences ##Tania G.
#trimmedDirectory <- paste(dataDirectory, "trimmed_seq/", sep="")

####################
## Analysis flags ##
####################
## These are TRUE/FALSE flags for whether a particular part of the analysis should be
## performed. If these are set to false, the resulting output files should already be 
## present in the appropriate directories.

# Download the reference genome from ftp/http site?
downloadGenome <- FALSE
# Generate FastQC reports for raw datafiles?
generateFastQC <- TRUE #09 Aug #Some fastq reports were run from here, others from the terminal (see my_MBD_seq_cheetah.txt and related txt files)
# Clean raw reads using default parameters?
cleanReads <- TRUE 
# Use cleaned reads for mapping? (This allows the read cleaning step to be separated
# from the mapping step.)
useCleanReads <- TRUE
# Build Bowtie2 index?
buildBTindex <- TRUE #I am not sure is it really necessary to build the index again...if I have it from BWA..
# Build BSgenome library?
buildBSgenome <- TRUE
# Map reads with Bowtie2?
mapReads <- TRUE
# Convert SAM files to sorted BAM files?
convertToBam <- TRUE
# Generate index files for BAM files?
generateIndex <- TRUE
# Download annotation files?
downloadAnnotation <- TRUE


######################
## Reference genome ## TANIA'S NOTE: At the time I did my analysis I couldn't make the script annotateDMRs_AciJub.R work
# but I corrected a few things that maybe should make it work now if people wanna give it a try. Though I STILL think 
#the issue lies at the conflict of chromosome-scaffold level. 
######################
## These variables store the source of the reference genome files and associated annotation.
## Even if they aren't used, these variables should be specified to ensure the source of 
## the files is recorded and easily determinable for all projects.

# Reference genome source
genomeSourceFiles <- "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/709/585/GCF_003709585.1_Aci_jub_2"
# Annotation source (if an annotation file is used)
annotationSource <- "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/709/585/GCF_003709585.1_Aci_jub_2/GCF_003709585.1_Aci_jub_2_genomic.gff.gz"
# version gff 2: "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/709/585/GCF_003709585.1_Aci_jub_2/GCF_003709585.1_Aci_jub_2_genomic.gtf.gz"

################################
## BSgenome package variables ##
################################
## These options are used for the BSgenome package creation. The seqnames parameter
## is often the most tedious (especially for large numbers of scaffolds). This will
## need to be improved in the future. 
#Tania G. It can be easier to create the BSgenome using the forgeBSGenome2_Cheetah.R script
# in which case the output will be and work as a more 
#"strict" BSgenome package (look for notes along scripts for differences on its usage)
## BSgenome package. 08/07/2019

# The seed file name for BSgenome package creation
seedName <- "BSgenome.AciJubGCF.NCBI.Acijub2-seed"
# The package name (NO UNDERSCORES OR DASHES!!!)
bsgenomePackageName <- "BSgenome.AciJub.NCBI.Acijub2"
# Name for the reference genome
referenceName <- "AciJub"
# Species name
species <- "Acinonyx jubatus"
# Source of reference genome
provider <- "NCBI_RefSeq"
# Reference genome version
version <- "2_Rico"
# Release date
release_date <- "2018/11/13"
# Names of all scaffolds or chromosomes
seqnames <-  "paste('NW_0',seq(from=020834726,to=020837944,by=1),'.1', sep='')"
# Any prefix needed to regenerate filename from seqnames
seqfiles_prefix <- "AciJub_2_"
# Any suffix needed to regenerate filename from seqnames
seqfiles_suffix <- ".fasta"

###############################
## Raw file / sample pairing ##
###############################
## This section associates files with samples. The files are assumed to be located in the
## dataDirectory folder. The full path will be based on that directory.

# The seqFiles data frame holds all sample information. The ctFlag is an identifier for the 
# treatment group, but any attribute of the sample can be used. 
seqFiles <- data.frame(
  sampleName = c("FRA135", "FRA138", "FRA171", "FRP096", "FRP114", "Z192", "Z194", "Z195", "Z196", "Z197"), 
  p1FileName = c("FRA135_R1.fastq.gz", "FRA138_R1.fastq.gz", "FRA171_R1.fastq.gz", "FRP096_R1.fastq.gz", "FRP114_R1.fastq.gz", "Z192_R1.fastq.gz", "Z194_R1.fastq.gz", "Z195_R1.fastq.gz", "Z196_R1.fastq.gz", "Z197_R1.fastq.gz"),
  p2FileName = c("FRA135_R2.fastq.gz", "FRA138_R2.fastq.gz", "FRA171_R2.fastq.gz", "FRP096_R2.fastq.gz", "FRP114_R2.fastq.gz", "Z192_R2.fastq.gz", "Z194_R2.fastq.gz", "Z195_R2.fastq.gz", "Z196_R2.fastq.gz", "Z197_R2.fastq.gz"),
  ctFlag = c("FR","FR","FR","FR","FR", "Zoo","Zoo","Zoo","Zoo","Zoo"),
  stringsAsFactors = F
)


#seqFiles <- data.frame(
 # sampleName = c("FRA135", "FRA138", "FRA171", "FRP096", "FRP114", "Z192", "Z194", "Z195", "Z196", "Z197"), 
  #p1FileName = c("FRA135_R1_clipped.fastq.gz", "FRA138_R1_clipped.fastq.gz", "FRA171_R1_clipped.fastq.gz", "FRP096_R1_clipped.fastq.gz", "FRP114_R1_clipped.fastq.gz", "Z192_R1_clipped.fastq.gz", "Z194_R1_clipped.fastq.gz", "Z195_R1_clipped.fastq.gz", "Z196_R1_clipped.fastq.gz", "Z197_R1_clipped.fastq.gz"),
  #p2FileName = c("FRA135_R2_clipped.fastq.gz", "FRA138_R2_clipped.fastq.gz", "FRA171_R2_clipped.fastq.gz", "FRP096_R2_clipped.fastq.gz", "FRP114_R2_clipped.fastq.gz", "Z192_R2_clipped.fastq.gz", "Z194_R2_clipped.fastq.gz", "Z195_R2_clipped.fastq.gz", "Z196_R2_clipped.fastq.gz", "Z197_R2_clipped.fastq.gz"),
  #ctFlag = c("Z", "FR"),
  #stringsAsFactors = F
#)

# The Trimmomatic cleaned files include four files for each sample, two paired files and 
# two unpaired files. By default, the file names are generated from the raw file names in 
# seqFiles (above). If the cleaning was done manually, or if the clean data filenames are
# different than the default below, cleanFileNames needs to be modified.
cleanFileNames <- data.frame(
  cp1FileName = paste(seqFiles$sampleName, "_R1_paired.fastq.gz", sep = ""),
  cp2FileName = paste(seqFiles$sampleName, "_R2_paired.fastq.gz", sep = ""),
  cs1FileName = paste(seqFiles$sampleName, "_R1_unpaired.fastq.gz", sep = ""),
  cs2FileName = paste(seqFiles$sampleName, "_R2_unpaired.fastq.gz", sep = ""),
  stringsAsFactors = F
)

# Files are automatically given names when converted to SAM and sorted BAM formats. All
# files are kept. These can be changed if necessary.
samFileName <- paste(seqFiles$sampleName, ".sam", sep="")
bamFileName <- paste(seqFiles$sampleName, ".bam", sep="")
sbamFileName <- paste(seqFiles$sampleName, ".bam", sep="")
#sbamFileName <- paste("s", bamFileName, sep="")


##################################
## Bowtie2 and SAMtools options ##
##################################
## These options determine how Bowtie2 maps the sample reads to the reference. Any
## Bowtie2 options can be specified using the otherMappingOptions variable.

# This defines the number of processing threads used by Bowtie2 for the mapping. It
# is also used for Trimmomatic read cleanning and trimming if necessary.
numThreads <- 7 
# All other parameters can be added here
otherMappingOptions <- "--no-unal"
# The index name can be specified, if different from referenceName due to manual creation.
indexName <- "AciJub"
# Multithreaded samtools may use excessive memory. Used for view and sort (conversion of
# SAM to BAM and sorted BAM.
numSamThreads <- 7

#######################
## MEDIPS parameters ##
#######################
## These parameters modify the MEDIPS analysis. The defaults shown here have not been fully
## explored. It isn't clear what the ideal values might be. Most of these are inherited from
## Haque's code. See MEDIPS documentation for details.

##Tania G. These parameters are the ****CORE**** of the analysis in many ways. 
#An effort for understanding them better and selecting them properly should be done!!!
#Reading the docs from Leinhard et al. helps a lot: Lienhard M., Chavez L. (2016). Quantitative Comparison of Large-Scale DNA Enrichment Sequencing Data chapter 10. 
#In: Ewy Mathe ́ and Sean Davis (Eds.), Statistical Genomics: Methods and Protocols, Methods in Molecular Biology, vol. 1418, DOI 10.1007/978-1-4939-3578-9_10, Springer Science-Business Media New York.
#Lienhard M, Grimm C, Morkel M, Herwig R, Chavez L (2014). 
#“MEDIPS: genome-wide differential coverage analysis of sequencing data derived from DNA enrichment experiments.” Bioinf, 30, 284-286.

# uniq = max allowed number of stacked reads per genomic position
#uniq = 1 would mean absolutely no dups are allowed
uniq <- 1
# extend = all reads will be extended to a length of xx nt according to strand information
extend <- 150
# shift is alterntaive to extend: reads are not extended but shifted by the specified number of nucleotides
shift <- 0
# window size
ws = 200
# subselection of chromosomes  
chr.select <- "NW_020834736.1" 
#TANIA'S NOTE:  chr.select THIS SETTING IS ONLY USED FOR RUNNING SOME QC tests 
# OR TRIALS BE SURE THAT IT'S NOT USED WHEN RUNNING THE ACTUAL ANALYSIS!!!

# correct p.values for multiple testing. Available are: holm, hochberg, hommel, bonferroni, BH, BY, fdr, none
p.adj <- "fdr"
# methods for calculating differential coverage (options: ttest and edgeR. ttest will be calculated only in case there are at least three replicates per group. edgeR can be applied when there are less than three replicates per group
diff.method <- "edgeR"

## Tania G. the following (particularly the following 3) parameters need more exploration. 
#Performing a trial with and without is still missing and should be attempted using maybe one scaffold to spot any differences.
# This parameter determines, if a CpG density dependent relative methylation scores (rms) will be calculated for the MEDIPS SETs given at the slots MSet1 and MSet2.
MeDIP <- FALSE #TANIA: DUNNO WHY THIS WAS SET TO FALSE IN THE ORIGINAL SCRIPTS FROM SKINNER'S LAB!!! BUT I THINK IT SHOULD BE SET TO TRUE!!! See more comments on 
#this below ... 

#TANIA: For CpG density dependent normalization of MeDIP-seq data, we need to generate a coupling set. 
#The coupling set must be created based on the same reference genome, the same set of chromosomes, and with the same window size used for the MEDIPS SETs. 
#For this, we specify the first MEDIPS SET in the hESCs object as reference, but any of the other MEDIPS SETs would be fine as well, 
#because all of them consist of the same set of chromosomes and have been generated with the same window size.
CScalc <- FALSE ### Tania G. In original script this was set to FALSE BUT again would be BETTER to re-run (using chr.select for trial) 
#WITH THIS SET TO TRUE and see if something changes. This it's recommended to compare to bisulfite data but since we don't have input-data to correct for 
#genetic background this might help a little by normilizing coverage based on CpG density
#Using this argument allows to use ttest instead of edgeR for differential coverage testing. 
#The advantage of ttest is that it can be applied to CpG density normalized rms (relative methylation score) values.  
#This would need the following arguments (which also need to be included in the medipAnalysis.R code): 
#CScalc <- TRUE 
#diff.method <- "ttest"
#MeDIP <- TRUE
#type <- "rms" 

#  copy number variation will be tested by applying the package DNAcopy to the window-wise log-ratios calculated based on the the means per group. 
#TANIA: THIS IS AN IMPORTANT ONE! BUT APPARENTLY IT NEEDS INPUT-SEQ DATA, WHICH I think would be having technical replicates w/o the enrichment for each condition
CNV <- FALSE

# threshold for a minimum sum of counts across all samples per window (default=10). Windows with lower coverage will not be tested for differential coverage.
minRowSum <- 10

# This vector holds all raw p-value thresholds to use for the analyses. 
pValues <- c(1e-03, 1e-04, 1e-05, 1e-06, 1e-07)
# This vector holds all multiple testing adjusted p-value thresholds to use for the analyses. 
MTCpValues <- c(0.3, 0.2, 0.1, 0.05)

####################
## DMR parameters ##  = Differentially methylated regions
####################
## These parameters define how DMR edges are defined. These values are completely arbitrarily
## chosen. It isn't clear how these values should be chosen to accuratly reflect a biologically
## significant feature.

# This p-value threshold defines DMR boundaries TANIA: AGAIN these settings need to be thorughly studied and chosen! 
dmrBoundPvalue <- 0.1
# Adjacency distance (=1 when windows must be exactly adjacent). This determines how far apart
# significant windows can be and remain in the same DMR.
adjDist <- 1000

# The maxDMRnum variable gives a maximum number of DMRs on which to calculate CpG density and other
# information. This speeds up the pipeline. However, this will need to be increased if the p-value 
# of interest has more DMR than this number.
maxDMRnum <- 5000

###########################
## Annotation parameters ##
###########################
## These variables hold information on the annotation. Three types of annotation are possible.
## The annotation is still being incorporated into the pipeline.

# The annotation can be from Biomart ("biomart"), from a GFF file ("gff"), or from BLAST results
# ("blast").
#Tania G: Annotation from gff did not work...it runs but does not add annotation...
#maybe a second look without the crazy rush might turn up into discovering that it was smth quite simple.
annotationType <- "gff"
# If annotation is from a GFF file, include the filename. TANIA: I WAS NEVER SURE WHETHER IT NEEDED THE PATH OR THE FILE NAME! BUT I TRIED BOTH
# AND STILL DIDN'T WORK LOL 
annotationGFF <- "/data/fg2/guerrero/cheetah_methylation/MeDIPS/data/genome/refGenome/Acijub_2_GCF_003709585.1_genomic.gff"
chrPrefix <- "NW_"

# If annotation is from Biomart, include biomaRt dataset information.
biomartHost <- ""
biomartDataset <- ""

############################
## Comparisons / analyses ##
############################
## This section defines which comparisons should be made. All analyses will iterate through
## each of these.

# Names for the comparisons ### TANIA: Not sure how this affected the analysis, this was the "default" but I got the feeling that it should be smth else
comparisonNames <- c("all")

# Which samples are being compared. The pairs flag was included to allow for a pairwise
# type analysis. It isn't currently functional, but some code requires it.
comparison <- list()
comparison[[1]] <- data.frame(mset1=c(1:5), mset2=c(6:10), pairs=F)

# Description added directly to report
comparisonDescription <- list()
comparisonDescription[[1]] <- "All free-ranging vs zoo cheetahs"

#################
## Data checks ##
#################
## These are very basic checks to detect problems in this configuration file. This section
## should be extended each time an unexpected error occurs in this script.

if (length(comparison)!=length(unique(comparisonNames))) { 
  warning("The number of comparison names does not match the number of comparisons.\n",
          "This could be caused by duplicate comparison names.")
}
if (CScalc & !MeDIP) {
  warning("The CSset coupling set is being calculated (CScalc = TRUE) 
          but not used (MeDIP = FALSE).")
}
if (!CScalc & MeDIP) {
  warning("The CSset coupling set is not calculated (CScalc = FALSE) 
          but is needed (MeDIP = TRUE).")
}

