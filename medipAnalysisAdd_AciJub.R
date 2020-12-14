## Created 5/27/2015 by Daniel Beck
## Last modified 2/22/2017
# Last modified Tania G. September, 2019 

## This code performs the MeDIP analysis. It is the second step in the analysis
## pipeline. The prepareData.R script should typically be run first. The dataNames.R
## configuration file is also used for this script.

## This latest version only loads the samples needed for the comparisons selected. It 
## does not perform any QC. It can be used when additional comparisons are added
## to the project without the addition of new samples.


# Load relevant libraries, custom functions, and the configuration script. Note, the
# order is important, as dataNames.R holds information about the bsgenomePackageName.
library(MEDIPS)
source("dataNames_AciJub.R")
source("customFunctions.R") 
library(bsgenomePackageName, character.only = T, lib.loc = genomeDirectory) #using the M.Skinner manually crafted genome forging

## TANIA: alternatively: 
#library(BSgenome.AciJubGCF.NCBI.Acijub2) # using forged genome according to manual
#BSgenome = "BSgenome.AciJubGCF.NCBI.Acijub2"  # using forged genome according to manual

#####################
## Quality control ## (TANIA: I moved this section to the beginning because some results from this section are actually necessary to assign analysis settings 
#e.g. window size depends on saturation analysis performance)
#####################
## This section calculates various quality control metrics and visualizations. While 
## several are included here, many of them take considerable time. I've commented out
## the more computationally intensive ones. These can be uncommented if more extensive
## QC is required.

## Next calculate saturation for each of the bam files (for quality control)
# nit= defines the number of subsets created frm the full sets of available regions (default=10)
# nrit=methods which randomly select data entries may be processed several times inorder to obtain more stable results. 
#By specifying the nrit parameter (default=1)it is possible to run the saturation analysis several times. 
#The final results re-turned to the saturation results object are the averaged results of each randomiteration step.
satList<-lapply(paste(dataDirectory, sbamFileName, sep = ""), 
               function(i) {
                  MEDIPS.saturation(file = i, BSgenome = bsgenomePackageName, uniq = uniq, 
                                    extend = extend, shift = shift, window_size = ws, 
                                    chr.select = chr.select, nit = 10, nrit = 1, 
                                    empty_bins = TRUE, rank = FALSE)
               })


## Next calculate CpG enrichment
enrichList<-lapply(paste(dataDirectory, sbamFileName, sep = ""),
                   function(i) {
                     MEDIPS.CpGenrich(file = i, BSgenome = bsgenomePackageName, uniq = uniq,
                                      extend = extend, shift = shift, chr.select = chr.select)
                   })

## This looks at coverage levels for the reference genome CpGs.
coverList <- lapply(paste(dataDirectory, sbamFileName, sep=""),
                    function(i) {
                      MEDIPS.seqCoverage(file =i, pattern = "CG", 
                                         BSgenome = bsgenomePackageName, 
                                         chr.select = chr.select, extend = extend, 
                                         shift = shift, uniq = uniq)
                    })


# Tania G. This measures the Pearson correlation in read depth between samples (to run this is necessary to first create processedBamFiles object, next section)
corMatrix = MEDIPS.correlation(MSets = processedBamFiles, plot = T, method = "pearson")
# pdf output manually moved to Results directory & ggplot formatted matrix it's on plots.R file 


# The QC results are saved to a RData file in the results directory. This line should be changed
# if the commented analyses above are used (to include satList and/or enrichList).
save(satList, coverList, corMatrix, enrichList, file = paste(resultsDirectory, "/qcLists.RData", sep=""))


#############################
## Read files and reformat ##
#############################
## This step reads in the sorted BAM sample files and converts them to a matrix of
## genomic windows with associated coverage (read counts). These options are defined
## in the dataNames.R configuration file.
#Tania G. ATTENTION HERE!!!!!!: 
## run with chr.select (dataNames.R) parameter for qc analysis
##run without chr.select parameter for actual data analysis

processedBamFiles <- lapply(X = paste(dataDirectory, sbamFileName, sep = ""), 
                            FUN = MEDIPS.createSet, 
                            BSgenome = bsgenomePackageName, 
                            extend = extend, 
                            shift = shift, 
                            uniq = uniq, 
                            window_size = ws 
                          )


# 27/03/2019
######## ###########
## Normalization ##
###################
## This step normalizes the data by calculating local CpG density. This is set to false
## by default in our pipline. This was taken from Haque's analysis and hasn't been 
## sufficiently explored. TANIA: THESE ARE THE SETTINGS I ALSO FLAGGED FOR MORE CONSIDERATION IN dataNames_AciJub.R!!!

if (CScalc) {
  CS <- MEDIPS.couplingVector(pattern = "CG", refObj = processedBamFiles[[1]])
} else {
  CS <- NULL
}
#CS <- MEDIPS.couplingVector(pattern = "CG", refObj = processedBamFiles[[1]])
######*Attention here: For my MSc I did not change the default (which was not running the CS) I RAN THE CS AFTER HAVING PERFORMED THE WHOLE ANALYSIS 
#JUST TO SEE HOW IT LOOKED
#SO IF RUN AGAIN IT might CHANGE, SINCE NOW THE CS IS NOT EMPTY----##### 27 SEP 2019


###################
## Identify DMRs ##
###################
## This step performs the actual DMR analysis. Each genomic window is assigned a 
## probability of being a DMR. All samples in mset1 are compared with all samples 
## in mset2. This code loops over all comparisons specified in the dataNames.R 
## configuration script.

for (analysis in 1:length(comparison)) {  
  mset1<-processedBamFiles[comparison[[analysis]]$mset1]
  mset2<-processedBamFiles[comparison[[analysis]]$mset2]
  # I'm not sure when this would occur. I think it was an early test.
  if (length(mset1)==0) {
    mset1 <- NULL
  }
  if (length(mset2)==0) {
    mset2 <- NULL
  }
  # Perform analysis
  methResults <- MEDIPS.meth(MSet1 = mset1, 
                             MSet2 = mset2,
                             p.adj = p.adj,
                             diff.method = diff.method,
                             MeDIP = MeDIP,
                             CNV = CNV,
                             CSet = CS,
                             minRowSum = minRowSum)
  # Save results to a comparison specific folder in the results directory
  system(paste("mkdir ", resultsDirectory, comparisonNames[analysis], sep=""))
  save(methResults, file=paste(resultsDirectory, 
                               comparisonNames[analysis], 
                               "/methResults.RData", sep=""))
  # Clean up unnecessary objects
  rm(methResults)
}

