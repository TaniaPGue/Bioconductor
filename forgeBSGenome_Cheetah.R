### R code from vignette source 'BSgenomeForge.Rnw'
### Encoding: UTF-8

##################################################



###################################################
### code chunk number 2: BSgenomeForge.Rnw:429-440
###################################################


###################################################
### code chunk number 3: BSgenomeForge.Rnw:453-455 (eval = FALSE)
###################################################
## library(BSgenome)
## forgeBSgenomeDataPkg("path/to/my/seed") 

###THIS IS AN EXAMPLE...
###################################################
### code chunk number 4: BSgenomeForge.Rnw:678-683
###################################################
library(BSgenome)
seed_files <- system.file("extdata", "GentlemanLab", package="BSgenome")
tail(list.files(seed_files, pattern="\\.masked-seed$"))
rn4_masked_seed <- list.files(seed_files, pattern="\\.rn4\\.masked-seed$", full.names=TRUE)
cat(readLines(rn4_masked_seed), sep="\n")


###################################################
### Tania Guerrero: THIS IS THE ACTUAL GENOME FORGE I DID THE VERY FRIST ATTEMPT WITH THE FIRST GENOME RELEASE ACIJUB_1 
#package library and documentation found in folder called "seed"  /data/fg2/guerrero/cheetah_methylation/MeDIPS/seed
###################################################
library(BSgenome)

#mseqnames <-  paste("NW_0",seq(from="15130604",to="15144985",by=1),".1", sep="")
# mseqnames <-  paste("NW_0", c("015130604":"015144985"),".1", sep="")

# Any prefix needed to regenerate filename from seqnames
seqfiles_prefix <- "AciJub_2_"
# Any suffix needed to regenerate filename from seqnames
seqfiles_suffix <- ".fasta"
seqs_srcdir = "/data/fg2/guerrero/cheetah_methylation/MeDIPS_test/genome/scaffolds"
#trace("forgeBSgenomeDataPkg", browser, exit=browser) # for debugging, Alex C. only
#forgeBSgenomeDataPkg("BSgenome.AciJubGCF.NCBI.Acijub1-seed", seqs_srcdir = "/data/fg2/guerrero/cheetah_methylation/MeDIPS_test/genome/scaffolds", destdir = "/data/fg2/guerrero/cheetah_methylation/MeDIPS_test/seed")

#after succesfully forging it go to terminal where seed file and forged is located and type
# R CMD build package_name
# R CMD check package_name
# R CMD INSTALL package_name
#
###################################################

#
###################################################
library(BSgenome.AciJubGCF.NCBI.Acijub1)
seed_files <- "/data/fg2/guerrero/cheetah_methylation/MeDIPS_test/seed/BSgenome.AciJubGCF.NCBI.Acijub1.Rcheck/BSgenome.AciJubGCF.NCBI.Acijub1"
tail(list.files(seed_files, pattern="-seed$"))

###Another Examples
## Display seed file for musFur1:
musFur1_seed <- list.files(seed_files, pattern="\\.musFur1-seed$", full.names=TRUE)
cat(readLines(musFur1_seed), sep="\n")

## Display seed file for rn4:
rn4_seed <- list.files(seed_files, pattern="\\.rn4-seed$", full.names=TRUE)
cat(readLines(rn4_seed), sep="\n")
### code chunk number 6: BSgenomeForge.Rnw:742-743
###################################################
sessionInfo()

