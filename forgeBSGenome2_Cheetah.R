### R code from vignette source 'BSgenomeForge.Rnw'

library(BSgenome)

###To run on the terminal if individual fasta files are not yet available (script from Johanna Paijmans)
#To generate the BSgenome yourself, you just need the individual scaffolds and then run Skinner's script. Here is how I did it:

#First, generate a fai file for your reference (it contains information on scaffold names, lengths, etc, and is very handy to have with your reference)
#samtools faidx AciJub_2_GCF_003709585.1.fasta

#Then I wrote a simple loop to extract the individual scaffolds from the reference, using your fai file:
#for i in $(awk '{print $1}' AciJub_2_GCF_003709585.1.fasta.fai); do samtools faidx AciJub_2_GCF_003709585.1.fasta $i > AciJub_2_"$i".fasta; done;

#mseqnames <-  paste("NW_0",seq(from="15130604",to="15144985",by=1),".1", sep="")
# mseqnames <-  paste("NW_0", c("015130604":"015144985"),".1", sep="")

## From there either run Skinner's script or the forgegenome accordign to the manual as follows: :)
# Any prefix needed to regenerate filename from seqnames
seqfiles_prefix <- "AciJub_2_"
# Any suffix needed to regenerate filename from seqnames
seqfiles_suffix <- ".fasta"
#getSeqSrcpaths <- function(seqnames, prefix="AciJub_2_", suffix=".fasta", seqs_srcdir=seqs_srcdir)
seqs_srcdir = "/data/fg2/guerrero/refGenome"
#trace("forgeBSgenomeDataPkg", browser, exit=browser) # for debugging, Alex C. only
forgeBSgenomeDataPkg("BSgenome.AciJubGCF.NCBI.Acijub2-seed", seqs_srcdir = "/data/fg2/guerrero/refGenome", destdir = "/data/fg2/guerrero/cheetah_methylation/MeDIPS_test/seed_Aci_Jub2")

#If for some reason this would not work, check the folder containing the individual scaffold files, maybe the presence of the .gff file interferes, it shouldn't but it might...

#after succesfully forging it go to terminal, to the folder where seed file and forged is located and type
# R CMD build package_name
# R CMD check package_namecd .. 
# R CMD INSTALL package_name
#
###################################################

# SANITY CHECK
###################################################
library(BSgenome.AciJubGCF.NCBI.Acijub2)
seed_files <- "/data/fg2/guerrero/cheetah_methylation/MeDIPS/seed_Aci_Jub2/BSgenome.AciJubGCF.NCBI.Acijub2.Rcheck/BSgenome.AciJubGCF.NCBI.Acijub2"
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

