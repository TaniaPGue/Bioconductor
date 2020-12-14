## Created 8/5/2015 by Daniel Beck
## Last modified 4/5/2016
##Last modified Tania G. September, 2019 

## This code processes the medipAnalysis results. It identifies DMR, defines DMR 
## edges, and generates the preliminary DMR tables.

# Load relevant libraries, custom functions, and the configuration script. Note, the
# order is important, as dataNames.R holds information about the bsgenomePackageName.
library(MEDIPS)
source("dataNames_AciJub.R")
source("customFunctions.R")
#library(BSgenome)
#library(BSgenome.AciJubGCF.NCBI.Acijub1)
#BSgenome= "BSgenome.AciJubGCF.NCBI.Acijub1"

library(bsgenomePackageName, character.only = T, lib.loc = genomeDirectory)
 

# This code loops over all comparisons specified in the dataNames.R configuration script.
for (analysis in 1:length(comparison)) {  
  # Results of the medipAnalysis.R script are loaded. 
  load(paste(resultsDirectory, comparisonNames[analysis], "/methResults.RData", sep=""))
  
  
  ############################
  ## Identify all DMR (raw) ##
  ############################
  ## The DMR are then identified using all p-value thresholds specified in the dataNames.R 
  ## configuration script. DMR edges are restricted to chromosome boundaries and CpG density 
  ## is calculated. This first analysis does not adjust p-values for multiple testing.
  
  methList <- list()
  methListEtc <- list()
  MTCmethList <- list()
  MTCmethListEtc <- list()
  temp <- extractDMR(x = methResults, distance = adjDist, 
                     pValueIdentify = pValues, pValueExtend = dmrBoundPvalue, 
                     mtcIdentify = FALSE, mtcExtend = FALSE)
  methList <- temp[[1]]
  methListEtc <- temp[[2]]
  rm(temp)
  # Loop over all p-value thresholds
  for (pV in 1:length(methList)) {
    if (!is.null(methList[[pV]])) {
      oneUse <- strsplit(as.character(rbind(methListEtc[[pV]])[, "edgeR.p.value"]), 
                         split = ";")
      methList[[pV]]$minP <- sapply(oneUse, function(i) min(as.numeric(i)))
      # Correct DMR edges that fall outside of chromosome boundaries
      methList[[pV]] <- modifyStop(dmrList = methList[[pV]], 
                                   refGenome = eval(parse(text = referenceName)), 
                                   maxDMR = maxDMRnum)
      # Add CpG density
      methList[[pV]] <- calcCpGdensity(methList[[pV]], maxDMRnum)
    }
  }
  

  ############################
  ## Identify all DMR (MTC) ##
  ############################
  ## This step identifies DMR in the same manner as the previous section. However, the p-value
  ## threshold used here is adjusted to account for multiple testing. The type of adjustment
  ## is defined in dataNames.R (but should probably be FDR).
  
  temp <- extractDMR(x = methResults, distance = adjDist, pValueIdentify = MTCpValues, 
                     pValueExtend = dmrBoundPvalue, mtcIdentify = TRUE, mtcExtend = FALSE)
  MTCmethList <- temp[[1]]
  MTCmethListEtc <- temp[[2]]
  rm(temp)
  # Loop over all p-value thresholds
  for (pV in 1:length(MTCmethList)) {
    if (!is.null(MTCmethList[[pV]])) {
      oneUse <- strsplit(as.character(rbind(MTCmethListEtc[[pV]])[, "edgeR.p.value"]), 
                         split = ";")
      MTCmethList[[pV]]$minP <- sapply(oneUse, function(i) min(as.numeric(i)))
      # Correct DMR edges that fall outside of chromosome boundaries
      MTCmethList[[pV]] <- modifyStop(dmrList = MTCmethList[[pV]], 
                                     refGenome = eval(parse(text = referenceName)),
                                      maxDMR = maxDMRnum)
      # Add CpG density
      MTCmethList[[pV]] <- calcCpGdensity(MTCmethList[[pV]], maxDMRnum)
    }
  }

  
  #################################
  ## Extract multiple window DMR ##
  #################################
  ## This section extracts the DMR that include multiple significant windows. These DMR tables 
  ## are named with a trailing "2p". This nomenclature is a holdover from Haque's code.
  
  methList2p <- lapply(methList, function(i) {
    if(!is.null(i)) {
      if (!is.na(i)) {
        i <- i[which(i$numSigWin >= 2), ]
      } else {
        i <- NA
      }
    } else {
      i <- NULL
    }
  })
  MTCmethList2p <- lapply(MTCmethList, function(i) {
    if (!is.null(i)) {
      if (!is.na(i)) {
        i <- i[which(i$numSigWin >= 2),]
      } else {
        i <- NA
      }
    } else {
      i <- NULL
    }
  })

  
  ################################
  ## Summarize and save results ##
  ################################
  ## This section creates tables of DMR numbers and saves the results to file. All files are
  ## saved in the results directory for each specific comparison.
  
  # Generate DMR number table summaries and save them in CSV format
  dmrNumberTable<-cbind("p-value" = pValues,
                        "allWindow" = sapply(methList, nrow),
                        "twoWindow" = sapply(methList2p, nrow))
  MTCdmrNumberTable<-cbind("p-value" = MTCpValues,
                           "allWindow" = sapply(MTCmethList, nrow),
                           "twoWindow" = sapply(MTCmethList2p, nrow))
  write.csv(file = paste(resultsDirectory, comparisonNames[analysis], 
                         "/dmrNumber.csv", sep = ""),
            x = dmrNumberTable, quote = F, row.names = FALSE) 
  write.csv(file = paste(resultsDirectory, comparisonNames[analysis], 
                         "/MTCdmrNumber.csv", sep = ""), 
            x = MTCdmrNumberTable, quote = F, row.names = FALSE) 
  
  # Save DMR tables and expanded DMR results (methListETC). The expanded results are saved
  # in a seperate file to speed loading of the more commonly used basic results.
  save(methList, methList2p, MTCmethList, MTCmethList2p, dmrNumberTable, MTCdmrNumberTable, 
       file = paste(resultsDirectory, comparisonNames[analysis], "/methLists.RData", sep = ""))
  save(methListEtc, MTCmethListEtc, 
       file = paste(resultsDirectory, comparisonNames[analysis], "/methListEtc.RData", sep = ""))
  
}


###Tania G. At this step some tables can be subsetted using the MEDIPS package functions (these tables were also used for plotting on plots.R)
# Note that these options are similar to the subsetting used by the Skinner's MEDIPS script, but that one apparently uses chromosome positions to merge.
#I mainly used the subsetting functions from MEDIPS package and not the ones from these scripts. 
####p-values of interest in a more simplified way 
#THIS CAN BE DONE AT ANY OTHER STEP AS NEEDED! 
p_val_relaxed <- MEDIPS.selectSig(methResults, p.value = 0.005, adj=F) 
p_val_adj <- MEDIPS.selectSig(methResults, p.value = 0.05, adj=T)
p_val_adj_rel <- MEDIPS.selectSig(methResults, p.value = 0.1, adj=T)




