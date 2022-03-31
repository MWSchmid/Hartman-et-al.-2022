#!/usr/bin/env Rscript

rm(list=ls())
rDir <- "/media/mwschmid/myData/MWSchmid/MarcelVanDerHeijden_rhizospheres/GitIgnore_results_16S_semiRarefied/OTU"
rDir <- "/media/mwschmid/myData/MWSchmid/MarcelVanDerHeijden_rhizospheres/GitIgnore_results_16S_semiRarefied_noRhizobia/OTU"

rDir <- "/media/mwschmid/myData/MWSchmid/MarcelVanDerHeijden_rhizospheres/GitIgnore_results_ITS_semiRarefied/OTU"
rDir <- "/media/mwschmid/myData/MWSchmid/MarcelVanDerHeijden_rhizospheres/GitIgnore_results_ITS_semiRarefied_noGlomeromycota/OTU"
rm(list=ls())

##########################################################################################
### arguments from commandline
myarg <- commandArgs()
argPos <- grep("--args", myarg, fixed = TRUE)
rDir <- as.character(myarg[argPos+1])

##########################################################################################
### libraries
suppressPackageStartupMessages({
  library("DESeq2")
  library("RColorBrewer")
  library("Rtsne")
  library("gplots")
  source("/media/mwschmid/myData/MWSchmid/Development/R/RNAseqWrapper.R")
  source("/media/mwschmid/myData/MWSchmid/Development/R/RNAseqWrapperAds.R")
  source("/media/mwschmid/myData/MWSchmid/Development/R/OTU_functions.R")
  source("/media/mwschmid/myData/MWSchmid/Development/R/lm_wrapper.R")
})

f.print.message <- function(...) { cat("=== ", format(Sys.time(), "%Y %b %d %X"), paste("===", ...,"\n")) }

pcname <- system('uname -n',intern=T)
numCoresPerMachine <- c("nuke" = 4, "styx" = 15, "marc-IEU" = 7, "piftel" = 30)
numCoresAvailable <- numCoresPerMachine[pcname]

##########################################################################################
### load data
load(file.path(rDir, "forDownstreamAnalyses.Rdata")) # contains: seqDat (filtered), seqDatRarefied (filtered and rarefied), normData (based on seqDat), sampleTab (has color and Pch), taxMod, taxModAll, otusToRemove, removeRhizobia, removeGlomeromycota
load(file.path(rDir, "allDataUsedInTests.Rdata")) # contains: forTest
load(file.path(rDir, "diffAbundance.Rdata")) # contains: diffAbunResultsDesign, diffAbunResultsDesignYesNo, diffAbunResultsDesignSpecies, noOutlierRefitting
load(file.path(rDir, "diffAbundanceOneToOne.Rdata")) # contains: diffAbunResultsOneToOne, specificOtus, speciesLFCmatrices, speciesFDRmatrices

otuDat <- seqDat[,rownames(forTest)]
normData <- normData[,rownames(sampleTab)]

sigThreshold <- 0.05
LFCthreshold <- 1

if (length(unique(sampleTab$Type_sample)) > 1) {
  f.print.message("THIS SCRIPT ASSUMES THAT YOU USE ONLY THE ROOT SAMPLES!")
  quit("no", 1)
}

sigOTUs <- lapply(c(diffAbunResultsDesignYesNo, diffAbunResultsDesignSpecies), function(x) x$get_significant_entries(sigThreshold, sigThreshold, LFCthreshold))
sigOTUsCoreMicrobiome <- list()
sigOTUsCoreMicrobiome$nn <- intersect(sigOTUs$AMFwithinNonRhizo$down, sigOTUs$bothVsNone$down)
sigOTUsCoreMicrobiome$An <- intersect(sigOTUs$AMFwithinNonRhizo$up, sigOTUs$RhizobiaWithinAMF$down)
sigOTUsCoreMicrobiome$AR <- intersect(sigOTUs$RhizobiaWithinAMF$up, sigOTUs$bothVsNone$up)

toUseInNetwork <- Reduce(union, sigOTUsCoreMicrobiome)
temp <- normData[toUseInNetwork, ]
minBinOccForMine <- ceiling(ncol(temp)*0.1)
binSum <- rowSums(seqDat[toUseInNetwork,]>0)
forMINE <- binSum>minBinOccForMine
write.csv(t(temp[forMINE,]), file.path(rDir, "coreMicrobiomeForNetwork.csv"), quote = FALSE, row.names = FALSE)