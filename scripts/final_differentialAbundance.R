#!/usr/bin/env Rscript

rDir <- "/media/mwschmid/myData/MWSchmid/MarcelVanDerHeijden_rhizospheres/GitIgnore_results_16S_semiRarefied/OTU"
rDir <- "/media/mwschmid/myData/MWSchmid/MarcelVanDerHeijden_rhizospheres/GitIgnore_results_ITS_semiRarefied/OTU"

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
  source("/media/mwschmid/myData/MWSchmid/MarcelVanDerHeijden_rhizospheres/scripts/final_createAnalysisSets.R")
})

f.print.message <- function(...) { cat("=== ", format(Sys.time(), "%Y %b %d %X"), paste("===", ...,"\n")) }

pcname <- system('uname -n',intern=T)
numCoresPerMachine <- c("nuke" = 4, "styx" = 15, "marc-IEU" = 7, "piftel" = 30)
numCoresAvailable <- numCoresPerMachine[pcname]

##########################################################################################
### load data
load(file.path(rDir, "forDownstreamAnalyses.Rdata")) # contains: seqDat (filtered), seqDatRarefied (filtered and rarefied), normData (based on seqDat), sampleTab (has color and Pch), taxMod, taxModAll, otusToRemove, removeRhizobia, removeGlomeromycota
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
diffAbunResultsDesign <- f.run.DESeq.models.test.only.last.term.in.parallel.likelihood.test(forDiffAbundanceDesign, otuDat, forTest, numCores = min(c(length(forDiffAbundanceDesign), numCoresAvailable)), noOutlierRefitting = TRUE)

##########################################################################################
### "high" vs "low" (yes vs no) to get senseful fold-changes, 
entryToTranslator <- c("nonRhizo" = "low", "nonAMF" = "low", "Rhizo" = "high", "AMF" = "high")
forTest$RhizobiaBinary <- factor(entryToTranslator[as.character(forTest$Rhizobia)], levels = c("low", "high"))
forTest$AMFbinary <- factor(entryToTranslator[as.character(forTest$AMF)], levels = c("low", "high"))
entryToTranslatorRhizobiaWithinAMF <- c("AMF_nonRhizo" = "low", "nonAMF_Rhizo" = "notOfInterestA", "AMF_Rhizo" = "high", "nonAMF_nonRhizo" = "notOfInterestB")
forTest$RhizobiaWithinAMF <- factor(entryToTranslatorRhizobiaWithinAMF[as.character(forTest$SymboGroup)], levels = c("low", "notOfInterestA", "notOfInterestB", "high"))
entryToTranslatorAMFwithinNonRhizo <- c("AMF_nonRhizo" = "high", "nonAMF_Rhizo" = "notOfInterestA", "AMF_Rhizo" = "notOfInterestB", "nonAMF_nonRhizo" = "low")
forTest$AMFwithinNonRhizo <- factor(entryToTranslatorAMFwithinNonRhizo[as.character(forTest$SymboGroup)], levels = c("low", "notOfInterestA", "notOfInterestB", "high"))
entryToTranslator <- c("AMF_nonRhizo" = "notOfInterestA", "nonAMF_Rhizo" = "notOfInterestB", "AMF_Rhizo" = "high", "nonAMF_nonRhizo" = "low")
forTest$bothVsNone <- factor(entryToTranslator[as.character(forTest$SymboGroup)], levels = c("low", "notOfInterestA", "notOfInterestB", "high"))
forDiffAbundanceDesign <- f.get.simple.analysis.sets.for.da.high.vs.low()
diffAbunResultsDesignYesNo <- f.run.DESeq.models.test.only.last.term.in.parallel.highest.vs.lowest(forDiffAbundanceDesign, otuDat, forTest, numCores = min(c(length(forDiffAbundanceDesign), numCoresAvailable)), noOutlierRefitting = TRUE)
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
diffAbunResultsDesignSpecies <- f.run.DESeq.models.test.only.last.term.in.parallel.highest.vs.lowest(detailedSpeciesContrastsList, otuDat, forTest, numCores = min(c(length(detailedSpeciesContrastsList), numCoresAvailable)), noOutlierRefitting = TRUE)

save(diffAbunResultsDesign, diffAbunResultsDesignYesNo, diffAbunResultsDesignSpecies, noOutlierRefitting, file = file.path(rDir, "diffAbundance.Rdata"))

if (file.exists(file.path(rDir, "diffAbun_factor_LRT_tests_DE_results.xlsx"))) { file.remove(file.path(rDir, "diffAbun_factor_LRT_tests_DE_results.xlsx")) }
if (file.exists(file.path(rDir, "diffAbun_species_vs_others_DE_results.xlsx"))) { file.remove(file.path(rDir, "diffAbun_species_vs_others_DE_results.xlsx")) }
if (file.exists(file.path(rDir, "diffAbun_treat_yes_vs_no_DE_results.xlsx"))) { file.remove(file.path(rDir, "diffAbun_treat_yes_vs_no_DE_results.xlsx")) }

f.write.DEGtabs.to.workbook(diffAbunResultsDesign, rDir, "diffAbun_factor_LRT_tests")
f.write.DEGtabs.to.workbook(diffAbunResultsDesignYesNo, rDir, "diffAbun_treat_yes_vs_no")
f.write.DEGtabs.to.workbook(diffAbunResultsDesignSpecies, rDir, "diffAbun_species_vs_others")

##########################################################################################
# do one vs one species to get very specific OTUs
f.do.species.contrast <- function(dds, speciesA, speciesB) {
  mainRes <- results(dds, contrast=c("SpeciesChar", speciesA, speciesB))
  mainResDf <- data.frame(generic = mainRes$baseMean, logFC = mainRes$log2FoldChange, pVal = mainRes$pvalue, adjP = mainRes$padj, row.names = rownames(mainRes))
  mainResDEGtab <- c.DEGtab(tool = "DESeq", method = "DESeq_default", table = mainResDf, isPairwise = FALSE, pairOrCont = paste0(speciesA, "_vs_", speciesB))
  return(mainResDEGtab)
}

formulaString = "~ Block + Pot + SpeciesChar"
dds <- DESeqDataSetFromMatrix(countData = otuDat, colData = forTest, design = as.formula(formulaString))
dds <- DESeq(dds, minReplicatesForReplace=Inf)
allSpecies <- unique(forTest$SpeciesChar)
diffAbunResultsOneToOne <- list()
speciesLFCmatrices <- list()
speciesFDRmatrices <- list()
specificOtus <- list()
for (speciesA in allSpecies) {
  allOtherSpecies <- setdiff(allSpecies, speciesA)
  temp <- mclapply(allOtherSpecies, function(x) f.do.species.contrast(dds, speciesA, x), mc.cores = numCoresAvailable)
  names(temp) <- allOtherSpecies
  specificOtus[[speciesA]] <- Reduce(intersect, lapply(temp, function(x) x$get_significant_entries(0.05, 0.05, 0)$up))
  diffAbunResultsOneToOne[[speciesA]] <- temp
  rowsForMatrix <- Reduce(union, lapply(temp, function(x) x$get_significant_entries(0.05, 0.05, 0)$up))
  lfcMat <- matrix(NA, nrow = length(rowsForMatrix), ncol = length(allOtherSpecies), dimnames = list(rowsForMatrix, allOtherSpecies))
  fdrMat <- matrix(NA, nrow = length(rowsForMatrix), ncol = length(allOtherSpecies), dimnames = list(rowsForMatrix, allOtherSpecies))
  for (speciesB in names(temp)) {
    commonOtus <- intersect(rownames(lfcMat), rownames(temp[[speciesB]]$table))
    lfcMat[commonOtus, speciesB] <- temp[[speciesB]]$table[commonOtus, "logFC"]
    fdrMat[commonOtus, speciesB] <- temp[[speciesB]]$table[commonOtus, "adjP"]
  }
  speciesLFCmatrices[[speciesA]] <- lfcMat
  speciesFDRmatrices[[speciesA]] <- fdrMat
} # only Lupin has one highly specific OTU, zot17

save(diffAbunResultsOneToOne, specificOtus, speciesLFCmatrices, speciesFDRmatrices, file = file.path(rDir, "diffAbundanceOneToOne.Rdata"))

















