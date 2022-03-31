#!/usr/bin/env Rscript

rDir <- "/media/mwschmid/myData/MWSchmid/MarcelVanDerHeijden_rhizospheres/GitIgnore_results_16S/OTU"

rDir <- "/media/mwschmid/myData/MWSchmid/MarcelVanDerHeijden_rhizospheres/GitIgnore_results_ITS/OTU"

##########################################################################################
### arguments from commandline
myarg <- commandArgs()
argPos <- grep("--args", myarg, fixed = TRUE)
rDir <- as.character(myarg[argPos+1])
noOutlierRefitting <- TRUE
if ("--refitOutliers" %in% myarg) {
  noOutlierRefitting <- FALSE
}

##########################################################################################
### libraries
suppressPackageStartupMessages({
  options(java.parameters = "-Xmx12g")
  library("parallel")
  library("DESeq2")
  source("/media/mwschmid/myData/MWSchmid/Development/R/RNAseqWrapper.R")
  source("/media/mwschmid/myData/MWSchmid/Development/R/RNAseqWrapperAds.R")
  source("/media/mwschmid/myData/MWSchmid/Development/R/lm_wrapper.R")
  source("/media/mwschmid/myData/MWSchmid/Development/R/OTU_functions.R")
  source("/media/mwschmid/myData/MWSchmid/MarcelVanDerHeijden_rhizospheres/scripts/createAnalysisSets.R")
})

f.print.message <- function(...) { cat("=== ", format(Sys.time(), "%Y %b %d %X"), paste("===", ...,"\n")) }

pcname <- system('uname -n',intern=T)
numCoresPerMachine <- c("nuke" = 4, "styx" = 15, "marc-IEU" = 7, "piftel" = 30)
numCoresAvailable <- numCoresPerMachine[pcname]

##########################################################################################
### load data
load(file.path(rDir, "forDownstreamAnalyses.Rdata")) # contains: seqDat (filtered and eventually rarefied), normData, sampleTab (has color and Pch), taxMod, taxModAll, keepOnlyRootSamples, removeRhizobia, useRarefaction
load(file.path(rDir, "allDataUsedInTests.Rdata")) # contains: forTest
otuDat <- seqDat[,rownames(forTest)]

if (length(unique(sampleTab$Type_sample)) > 1) {
  f.print.message("THIS SCRIPT ASSUMES THAT YOU USE ONLY THE ROOT SAMPLES!")
  quit("no", 1)
}

##########################################################################################
### two things, one is the design, the other the different contrasts with the species
### Though, for the design, use Block and Pot, then test AMF/Rhizobia at the end, so all kind of the same
forDiffAbundanceDesign <- f.get.simple.analysis.sets.for.da()
diffAbunResultsDesign <- f.run.DESeq.models.test.only.last.term.in.parallel.likelihood.test(forDiffAbundanceDesign, otuDat, forTest, numCores = min(c(length(forDiffAbundanceDesign), numCoresPerMachine)), noOutlierRefitting = TRUE)

##########################################################################################
### "high" vs "low" (yes vs no) to get senseful fold-changes, 
entryToTranslator <- c("nonRhizo" = "low", "nonAMF" = "low", "Rhizo" = "high", "AMF" = "high")
forTest$RhizobiaBinary <- factor(entryToTranslator[as.character(forTest$Rhizobia)], levels = c("low", "high"))
forTest$AMFbinary <- factor(entryToTranslator[as.character(forTest$AMF)], levels = c("low", "high"))
forDiffAbundanceDesign <- f.get.simple.analysis.sets.for.da.high.vs.low()
diffAbunResultsDesignYesNo <- f.run.DESeq.models.test.only.last.term.in.parallel.highest.vs.lowest(forDiffAbundanceDesign, otuDat, forTest, numCores = min(c(length(forDiffAbundanceDesign), numCoresPerMachine)), noOutlierRefitting = TRUE)
#forDiffAbundanceDesign <- f.get.simple.analysis.sets.for.da.high.vs.low.species.last()
#diffAbunResultsDesignYesNoSpeciesLast <- f.run.DESeq.models.test.only.second.last.term.in.parallel.highest.vs.lowest(forDiffAbundanceDesign, otuDat, forTest, numCores = min(c(1, numCoresPerMachine)), noOutlierRefitting = TRUE)

##########################################################################################
### "high" vs "low" (yes vs no) to get senseful fold-changes, one species vs rest. 
detailedSpeciesContrastsNames <- grep("_vs_others$", colnames(forTest), value = TRUE)
#detailedSpeciesContrastsNames <- detailedSpeciesContrastsNames[detailedSpeciesContrastsNames!="Lupin_vs_others"]
for (curCont in detailedSpeciesContrastsNames) {
  forTest[[curCont]] <- as.character(forTest[[curCont]])
  forTest[[curCont]][forTest[[curCont]] == "others"] <- "low"
  forTest[[curCont]][forTest[[curCont]] != "low"] <- "high"
  forTest[[curCont]] <- factor(forTest[[curCont]], levels = c("low", "high"))
}
detailedSpeciesContrastsList <- lapply(detailedSpeciesContrastsNames, function(x) f.get.species.contrasts.da(x)); names(detailedSpeciesContrastsList) <- detailedSpeciesContrastsNames
diffAbunResultsDesignSpecies <- f.run.DESeq.models.test.only.last.term.in.parallel.highest.vs.lowest(detailedSpeciesContrastsList, otuDat, forTest, numCores = min(c(length(detailedSpeciesContrastsList), numCoresPerMachine)), noOutlierRefitting = TRUE)

save(diffAbunResultsDesign, diffAbunResultsDesignYesNo, diffAbunResultsDesignSpecies, noOutlierRefitting, file = file.path(rDir, "diffAbundance.Rdata"))

if (file.exists(file.path(rDir, "diffAbun_factor_LRT_tests_DE_results.xlsx"))) { file.remove(file.path(rDir, "diffAbun_factor_LRT_tests_DE_results.xlsx")) }
if (file.exists(file.path(rDir, "diffAbun_species_vs_others_DE_results.xlsx"))) { file.remove(file.path(rDir, "diffAbun_species_vs_others_DE_results.xlsx")) }
if (file.exists(file.path(rDir, "diffAbun_treat_yes_vs_no_DE_results.xlsx"))) { file.remove(file.path(rDir, "diffAbun_treat_yes_vs_no_DE_results.xlsx")) }

f.write.DEGtabs.to.workbook(diffAbunResultsDesign, rDir, "diffAbun_factor_LRT_tests")
f.write.DEGtabs.to.workbook(diffAbunResultsDesignYesNo, rDir, "diffAbun_treat_yes_vs_no")
f.write.DEGtabs.to.workbook(diffAbunResultsDesignSpecies, rDir, "diffAbun_species_vs_others")

