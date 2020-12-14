## Daniel Beck
## Created 10/5/2015
## Modified 10/8/2015
## Modified by Tania G. September, 2019

## This code is intended to add annotation to DMRs identified by the medipProcessing script.

library(MEDIPS)
source("dataNames_AciJub.R")
source("customFunctions.R")
#library(BSgenome.AciJubGCF.NCBI.Acijub2) # using forged genome according to manual
#BSgenome= "BSgenome.AciJubGCF.NCBI.Acijub2"  # using forged genome according to manual

library(bsgenomePackageName, character.only = T, lib.loc = genomeDirectory)
#library(biomaRt)

# Moving this here so it is only called once. It fails occationally for unknown reasons. (TANIA: OCCASIONALLY?? ¬¬ ) 

#Tania G. So yes indeed fails, not occasionally but every time I tried...
# I think it has to do with the absence of actual chromosomes, bc this code seems to be written for chromosomes, and the cheetah genome is still in scaffolds...dunno
# maybe the problem might is at fx (folder functions) addAnnotationGFF.R in the findOveralps() 
#but might be something else or multiple things...
#it happens that it runs but the annotation columns are empty (NA's)
if (annotationType == "biomart") {
 annotationObject <- useMart(biomart = "ensembl",
                             dataset = biomartDataset, 
                             host = biomartHost)
}

for (analysis in 1:length(comparisonNames)) { 
  if (!(comparison[[analysis]]$pairs[1])) {  # Don't calculate annotation for pair analysis
    load(paste(resultsDirectory, comparisonNames[analysis], "/methLists.RData", sep = ""))
       
    ####################
    ## Add Annotation ##
    ####################
    
    annMat <- list()
    MTCannMat <- list()
    # add Annotation information for each DMR. This is separate from the loop so that the 
    # import.gff and getAnnotation functions are only run once
    if (annotationType == "gff") {
      gff <- import.gff(paste(genomeDirectory, annotationGFF, sep = "")) ### Tania G. maybe the issue is here...17 dic/19 although it does fetch the right file...
      #changing gff names so that the findOverlaps fx works... maybe...
      
      methList <- lapply(methList, function(i) {
                           addAnnotationGFF(dmrList = i, gff = gff, ### Tania G. I think the problem is that in the fx addAnnotationGFF the findOveralps fails to find any overlap
                                            maxDMR = maxDMRnum, 
                                            chrPrefix = chrPrefix)
          #02 Sept Tomas: tried using the gff@seqnames@values as chrPrefix 
                           })
      MTCmethList <- lapply(MTCmethList, function(i) {
                              addAnnotationGFF(dmrList = i, gff = gff, 
                                               maxDMR = maxDMRnum, 
                                               chrPrefix = chrPrefix)
                              })
    }
    
   if (annotationType == "biomart") {
      annotationObject <- useMart(biomart = "ensembl",
                                   dataset = biomartDataset, 
                                host = biomartHost)
            
      for (i in 1:length(methList)) {
       a <- addAnnotationBiomart(dmrList = methList[[i]], extension = 10000,
                               annotationObject = annotationObject, 
                               maxDMR = maxDMRnum, chrPrefix = chrPrefix)
       methList[[i]] <- a$dmrList
       annMat[[i]] <- a$annMat
    
     }
      for (i in 1:length(MTCmethList)) {
      b <- addAnnotationBiomart(dmrList = MTCmethList[[i]], extension = 10000,
                                  annotationObject = annotationObject, 
                                  maxDMR = maxDMRnum, chrPrefix = chrPrefix)
        MTCmethList[[i]] <- b$dmrList
        MTCannMat[[i]] <- b$annMat
      }
    }
       
    ## Re-calculate twoWindow DMRs to include annotation information
       
    ##################################
    ## Multiple significant windows ##
    ##################################
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
       
    save(methList, methList2p, MTCmethList, MTCmethList2p, 
         dmrNumberTable, MTCdmrNumberTable, annMat, MTCannMat, 
         file = paste(resultsDirectory, comparisonNames[analysis], "/methLists.RData", sep = ""))
  }
  print(paste("Progress: ", analysis / length(comparison), "%", sep=""))
   
}


#Tania G. export MTCmethlist for testing annotation outside bioconductor (THIS IS JUST AN EXAMPLE TABLES CAN BE SAVED AT ANY CONVENIENT MOMENT)
write.table(as.data.frame(MTCmethList2p[[4]]),file="MTCmethlistresults_p-val_0.05.csv", quote=F,sep=",",row.names=F)
write.table(as.data.frame(methList2p[[5]]),file="methlistresults_p-val_1e-07.csv", quote=F,sep=",",row.names=F)

