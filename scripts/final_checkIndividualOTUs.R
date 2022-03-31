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

if (grepl("GitIgnore_results_ITS", rDir)) {
  isITS <- TRUE
  dataDir <- "/media/mwschmid/myData/MWSchmid/MarcelVanDerHeijden_rhizospheres/GitIgnore_usearchOutput_ITS"
  fgAnno <- read.table(file.path(dataDir, "otus_funGuild_matched_forR.txt"), sep = '\t', stringsAsFactors = FALSE, header = TRUE, row.names = 1)
  colnames(fgAnno) <- gsub("\\.", "", colnames(fgAnno))
  funguildColumnsToAdd <- c("TrophicMode", "Guild", "Trait", "ConfidenceRanking", "GrowthMorphology")
} else {
  isITS <- FALSE
}

otuDat <- seqDat[,rownames(forTest)]
normData <- normData[,rownames(sampleTab)]
meanData <- f.summarize.columns(normData, data.frame(sample = rownames(sampleTab), group = sampleTab$Species_name), mean)

sigThreshold <- 0.05
LFCthresholdForOneToOne <- 0
minSigCompsForOneToOne <- 10

if (length(unique(sampleTab$Type_sample)) > 1) {
  f.print.message("THIS SCRIPT ASSUMES THAT YOU USE ONLY THE ROOT SAMPLES!")
  quit("no", 1)
}

##########################################################################################
### Annotate the specific OTUs
allSpecies <- unique(forTest$SpeciesChar)
sigOtusForAnnot <- list()
for (curSpecies in allSpecies) {
  temp <- rowSums((speciesFDRmatrices[[curSpecies]] < sigThreshold) & (speciesLFCmatrices[[curSpecies]] > LFCthresholdForOneToOne), na.rm = TRUE)
  mask <- temp >= minSigCompsForOneToOne
  curSigOtus <- rownames(speciesFDRmatrices[[curSpecies]])[mask]
  curSigOtuCounts <- temp[mask]
  if (length(curSigOtus) > 0) {
    sigOtusForAnnot[[curSpecies]] <- data.frame(otu = curSigOtus, plantSpecies = rep(curSpecies, length(curSigOtus)), numSigComps = curSigOtuCounts, row.names = NULL, stringsAsFactors = FALSE)
  }
}
sigOtusForAnnot <- do.call("rbind", sigOtusForAnnot); rownames(sigOtusForAnnot) <- NULL
for (toAdd in colnames(taxMod)) {
  sigOtusForAnnot[[toAdd]] <- taxMod[sigOtusForAnnot$otu, toAdd]
}
if (isITS) {
  for (toAdd in funguildColumnsToAdd) {
    mask <- sigOtusForAnnot$otu %in% rownames(fgAnno)
    sigOtusForAnnot[[toAdd]] <- NA
    sigOtusForAnnot[mask, toAdd] <- fgAnno[sigOtusForAnnot$otu[mask], toAdd]
  }
}

write.csv(sigOtusForAnnot, file.path(rDir, "plantSpeciesSpecificOtusAnnotated.csv"), row.names = FALSE)

##########################################################################################
### Do enrichment for all of them (i.e., are there specific OTUs that are more likely to be plant species specific)
universe <- rownames(normData)
testSet <- unique(sigOtusForAnnot$otu)

for (curTaxum in setdiff(colnames(taxMod), "domain")) {
  annoMap <- taxMod[,curTaxum]; names(annoMap) <- rownames(taxMod)
  allTaxum <- annoMap[universe]; length(table(allTaxum))
  enriTest <- f.test.annotation.enrichment(testSet, universe, annoMap)
  numSigGroups <- sum(enriTest$FDR < 0.05)
  if (numSigGroups > 0) {
    f.print.message("number significant groups:", numSigGroups, "in", curTaxum)
    print(subset(enriTest, FDR < 0.05))
  }
}

if (isITS) {
  for (curThing in funguildColumnsToAdd) {
    annoMap <- fgAnno[,curThing]; names(annoMap) <- rownames(fgAnno)
    allTaxum <- annoMap[universe]; length(table(curThing))
    enriTest <- f.test.annotation.enrichment(testSet, universe, annoMap)
    numSigGroups <- sum(enriTest$FDR < 0.05)
    if (numSigGroups > 0) {
      f.print.message("number significant groups:", numSigGroups, "in", curThing)
      print(subset(enriTest, FDR < 0.05))
    }
  }
}

##########################################################################################
### do this also for the core microbiomes
sigThreshold <- 0.05
LFCthreshold <- 1
sigOTUs <- lapply(c(diffAbunResultsDesignYesNo, diffAbunResultsDesignSpecies), function(x) x$get_significant_entries(sigThreshold, sigThreshold, LFCthreshold))
sigOTUsCoreMicrobiome <- list()
sigOTUsCoreMicrobiome$nn <- intersect(sigOTUs$AMFwithinNonRhizo$down, sigOTUs$bothVsNone$down)
sigOTUsCoreMicrobiome$An <- intersect(sigOTUs$AMFwithinNonRhizo$up, sigOTUs$RhizobiaWithinAMF$down)
sigOTUsCoreMicrobiome$AR <- intersect(sigOTUs$RhizobiaWithinAMF$up, sigOTUs$bothVsNone$up)

universe <- rownames(normData)
for (curSet in names(sigOTUsCoreMicrobiome)) {
  testSet <- unique(sigOTUsCoreMicrobiome[[curSet]])
  for (curTaxum in setdiff(colnames(taxMod), "domain")) {
    annoMap <- taxMod[,curTaxum]; names(annoMap) <- rownames(taxMod)
    allTaxum <- annoMap[universe]; length(table(allTaxum))
    enriTest <- f.test.annotation.enrichment(testSet, universe, annoMap)
    numSigGroups <- sum(enriTest$FDR < 0.05)
    if (numSigGroups > 0) {
      f.print.message("number significant groups:", numSigGroups, "in", curTaxum, "(", curSet, ")")
      print(subset(enriTest, FDR < 0.05))
      write.csv(enriTest, file.path(rDir, paste0("coreMicrobiome_tax_", curTaxum, '_', curSet,"_enrichmentAny.csv")))
    }
  }
  if (isITS) {
    for (curThing in funguildColumnsToAdd) {
      annoMap <- fgAnno[,curThing]; names(annoMap) <- rownames(fgAnno)
      allTaxum <- annoMap[universe]; length(table(curThing))
      enriTest <- f.test.annotation.enrichment(testSet, universe, annoMap)
      numSigGroups <- sum(enriTest$FDR < 0.05)
      if (numSigGroups > 0) {
        f.print.message("number significant groups:", numSigGroups, "in", curThing, "(", curSet, ")")
        print(subset(enriTest, FDR < 0.05))
        write.csv(enriTest, file.path(rDir, paste0("coreMicrobiome_tax_", curTaxum, '_', curSet,"_enrichmentAny.csv")))
      }
    }
  }
}






















