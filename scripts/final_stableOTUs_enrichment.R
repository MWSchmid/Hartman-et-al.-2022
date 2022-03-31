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
numCoresPerMachine <- c("nuke" = 4, "styx" = 15, "marc-IEU" = 7, "mininuke" = 2, "piftel" = 30)
numCoresAvailable <- numCoresPerMachine[pcname]

##########################################################################################
### load data
load(file.path(rDir, "forDownstreamAnalyses.Rdata")) # contains: seqDat (filtered), seqDatRarefied (filtered and rarefied), normData (based on seqDat), sampleTab (has color and Pch), taxMod, taxModAll, otusToRemove, removeRhizobia, removeGlomeromycota
load(file.path(rDir, "allDataUsedInTests.Rdata")) # contains: forTest
load(file.path(rDir, "diffAbundance.Rdata")) # contains: diffAbunResultsDesign, diffAbunResultsDesignYesNo, diffAbunResultsDesignSpecies, noOutlierRefitting
load(file.path(rDir, "diffAbundanceOneToOne.Rdata")) # contains: diffAbunResultsOneToOne, specificOtus, speciesLFCmatrices, speciesFDRmatrices

otuDat <- seqDat[,rownames(forTest)]
normData <- normData[,rownames(sampleTab)]
meanData <- f.summarize.columns(normData, data.frame(sample = rownames(sampleTab), group = sampleTab$Species_name), mean)

if (grepl("ITS", rDir)) {
  dataDir <- "/media/mwschmid/myData/MWSchmid/MarcelVanDerHeijden_rhizospheres/GitIgnore_usearchOutput_ITS"
  fgAnno <- read.table(file.path(dataDir, "otus_funGuild_matched_forR.txt"), sep = '\t', stringsAsFactors = FALSE, header = TRUE, row.names = 1)
  colnames(fgAnno) <- gsub("\\.", "", colnames(fgAnno))
  colsToTest <- c("TrophicMode", "Guild", "Trait", "ConfidenceRanking", "GrowthMorphology")
}

if (length(unique(sampleTab$Type_sample)) > 1) {
  f.print.message("THIS SCRIPT ASSUMES THAT YOU USE ONLY THE ROOT SAMPLES!")
  quit("no", 1)
}

##########################################################################################
### get an invariant set of OTUs as the OTUs that are not significant in any of the one-to-one comparisons
sigOTUsOneToOne <- lapply(diffAbunResultsOneToOne, function(y) lapply(y, function(x) x$get_significant_entries(0.05, 0.05, 0)$any))
sigOTUsOneToOne <- lapply(sigOTUsOneToOne, function(x) Reduce(union, x))
sigOTUsOneToOne <- Reduce(union, sigOTUsOneToOne)
sigOTUsOthers <- lapply(c(diffAbunResultsDesignYesNo, diffAbunResultsDesignSpecies), function(x) x$get_significant_entries(0.05, 0.05, 0)$any)
sigOTUsOthers <- Reduce(union, sigOTUsOthers)
anySig <- union(sigOTUsOthers, sigOTUsOneToOne)
stableOTUs <- setdiff(rownames(normData), anySig)
f.print.message("found", length(stableOTUs), "stable OTUs.")
write.csv(normData[stableOTUs,], file.path(rDir, "stable_OTUs.csv"))

##########################################################################################
### annotation enrichment, spider plot and heatmaps
universe <- rownames(normData)
for (curTaxum in setdiff(colnames(taxMod), "domain")) {
  annoMap <- taxMod[,curTaxum]; names(annoMap) <- rownames(taxMod)
  allTaxum <- annoMap[universe]; length(table(allTaxum))
  allTaxumCounts <- sort(table(allTaxum), decreasing = TRUE)
  # enrichment test
  enriTest <- f.test.annotation.enrichment(stableOTUs, universe, annoMap)
  if (sum(enriTest$FDR < 0.1) > 0) { f.print.message("number significant groups (", curTaxum, "):", sum(enriTest$FDR < 0.1)) }
  #f.print.message("number significant groups (", curTaxum, "):", sum(enriTest$pValue < 0.05))
  if (curTaxum == "family") {
    print(subset(enriTest, pValue < 0.05))
  }
  write.csv(enriTest, file.path(rDir, paste0("stableOTUs_tax_", curTaxum, "_enrichment.csv")))
  annoMap <- taxModAll[,curTaxum]; names(annoMap) <- rownames(taxModAll)
  allTaxum <- annoMap[universe]; length(table(allTaxum))
  allTaxumCounts <- sort(table(allTaxum), decreasing = TRUE)
  # enrichment test
  enriTest <- f.test.annotation.enrichment(stableOTUs, universe, annoMap)
  if (sum(enriTest$FDR < 0.1) > 0) { f.print.message("number significant groups, entire annotation (", curTaxum, "):", sum(enriTest$FDR < 0.1)) }
  #f.print.message("number significant groups (", curTaxum, "):", sum(enriTest$pValue < 0.05))
  if (curTaxum == "family") {
    print(subset(enriTest, pValue < 0.05))
  }
  write.csv(enriTest, file.path(rDir, paste0("stableOTUs_tax_entireAnnotation_", curTaxum, "_enrichment.csv")))
}

