#!/usr/bin/env Rscript

rDir <- "/media/mwschmid/myData/MWSchmid/MarcelVanDerHeijden_rhizospheres/GitIgnore_results_16S/OTU"

rDir <- "/media/mwschmid/myData/MWSchmid/MarcelVanDerHeijden_rhizospheres/GitIgnore_results_ITS/OTU"

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
load(file.path(rDir, "forDownstreamAnalyses.Rdata")) # contains: seqDat (filtered and eventually rarefied), normData, sampleTab (has color and Pch), taxMod, taxModAll, keepOnlyRootSamples, removeRhizobia, useRarefaction
load(file.path(rDir, "allDataUsedInTests.Rdata")) # contains: forTest
load(file.path(rDir, "diffAbundance.Rdata")) # contains: diffAbunResultsDesign, diffAbunResultsDesignYesNo, diffAbunResultsDesignSpecies, noOutlierRefitting
otuDat <- seqDat[,rownames(forTest)]
normData <- normData[,rownames(sampleTab)]
meanData <- f.summarize.columns(normData, data.frame(sample = rownames(sampleTab), group = sampleTab$Species_name), mean)

sigThreshold <- 0.05
LFCthreshold <- 1
LFCthresholdLRT <- 0

if (length(unique(sampleTab$Type_sample)) > 1) {
  f.print.message("THIS SCRIPT ASSUMES THAT YOU USE ONLY THE ROOT SAMPLES!")
  quit("no", 1)
}

##########################################################################################
### get all significant OTUs
sigOTUsLRT <- lapply(diffAbunResultsDesign, function(x) x$get_significant_entries(sigThreshold, sigThreshold, LFCthresholdLRT))
sigOTUs <- lapply(c(diffAbunResultsDesignYesNo, diffAbunResultsDesignSpecies), function(x) x$get_significant_entries(sigThreshold, sigThreshold, LFCthreshold))
anySigOtuWithLRT <- Reduce(union, lapply(c(sigOTUs, sigOTUsLRT), function(x) x$any))
anySigOtu <- Reduce(union, lapply(sigOTUs , function(x) x$any))
numSigOtus <- do.call("rbind", lapply(sigOTUs, function(x) unlist(lapply(x, length))))
numSigOtusLTR <- do.call("rbind", lapply(sigOTUsLRT, function(x) unlist(lapply(x, length))))
numSigOtusLTR[,"down"] <- "lfcNotUseful"
numSigOtusLTR[,"up"] <- "lfcNotUseful"
write.csv(numSigOtus, file.path(rDir, "numSigOtus.csv"), quote = FALSE)
write.csv(numSigOtusLTR, file.path(rDir, "numSigOtusLTR.csv"), quote = FALSE)

f.print.message("found", length(anySigOtu), "significant OTUs, with LRT included:", length(anySigOtuWithLRT))

##########################################################################################
##### tSNE with significant OTUs
doTSNE <- TRUE
if (doTSNE) {
  # thanks to: https://github.com/lmweber/Rtsne-example/blob/master/Rtsne_viSNE_example_Marrow1.R
  # 101 random seeds, select fit with lowest error
  # run, select, and plot 2D t-SNE projection - color by species richness and pch by soil type and plant history
  forColoring <- list( default = forTest$color )
  forTSNE <- normData[anySigOtu, rownames(forTest)]
  rtsneResults <- f.run.select.plot.tSNE(unique(t(as.matrix(forTSNE))), forColoring, forTest$Pch, rDir, "tSNE_sigAnyCont.svg")
  forTSNE <- normData[anySigOtuWithLRT, rownames(forTest)]
  rtsneResults <- f.run.select.plot.tSNE(unique(t(as.matrix(forTSNE))), forColoring, forTest$Pch, rDir, "tSNE_sigAnyContWithLRT.svg")
}

##########################################################################################
##### Get VST data to account for Block and Pot
# from the DESeq-manual:
# It is possible to visualize the transformed data with batch variation removed, using the removeBatchEffect function from limma.
# This simply removes any shifts in the log2-scale expression data that can be explained by batch. The paradigm for this operation would be:
doVSTplot <- FALSE # it did not look so convincing
if (doVSTplot) {
  formulaString = "~ Block + Pot + Species_name"
  dds <- DESeqDataSetFromMatrix(countData = otuDat, colData = sampleTab, design = formula(formulaString))
  if (noOutlierRefitting) {
    dds <- DESeq(dds, minReplicatesForReplace=Inf)
  } else {
    dds <- DESeq(dds)
  }
  vstData <- vst(dds, blind=FALSE)
  # correction
  mat <- assay(vstData)
  mat <- limma::removeBatchEffect(mat, vstData$Block, vstData$Pot)
  assay(vstData) <- mat
  write.csv(mat, file.path(rDir, "vstCorrectedData.csv"))
  forColoring <- list( default = forTest$color )
  forTSNE <- mat[anySigOtu, rownames(forTest)]
  rtsneResults <- f.run.select.plot.tSNE(unique(t(as.matrix(forTSNE))), forColoring, forTest$Pch, rDir, "tSNE_vstCorrected_sigAnyCont.svg")
  forTSNE <- mat[anySigOtuWithLRT, rownames(forTest)]
  rtsneResults <- f.run.select.plot.tSNE(unique(t(as.matrix(forTSNE))), forColoring, forTest$Pch, rDir, "tSNE_vstCorrected_sigAnyContWithLRT.svg")
}

##########################################################################################
##### heatmap with the species LFCs
speciesToColor <- unique(sampleTab[,c("Species_name", "color")]); rownames(speciesToColor) <- speciesToColor$Species_name

sigOTUsSpecies <- lapply(diffAbunResultsDesignSpecies, function(x) x$get_significant_entries(sigThreshold, sigThreshold, LFCthreshold))
anySigOtuSpecies <- Reduce(union, lapply(sigOTUsSpecies, function(x) x$any))
lfcMatrix <- matrix(0, nrow = length(anySigOtuSpecies), ncol = nrow(speciesToColor), dimnames = list(anySigOtuSpecies, speciesToColor$Species_name))
for (curCont in names(diffAbunResultsDesignSpecies)) {
  curSpecies <- gsub("_vs_others", "", curCont)
  curSigOtus <- sigOTUsSpecies[[curCont]]$any
  lfcMatrix[curSigOtus, curSpecies] <- diffAbunResultsDesignSpecies[[curCont]]$table[curSigOtus, "logFC"]
}

if (max(abs(range(lfcMatrix))) > 5) {
  f.print.message("limiting to a LFC of 5")
  lfcMatrix[lfcMatrix < (-5)] <- (-5)
  lfcMatrix[lfcMatrix > 5] <- 5
}

pdf(file.path(rDir, paste0("speciesContrasts_LFCmatrix_", LFCthreshold, "_maxFDR_", sigThreshold, ".pdf")), width = 10, height = 12) 
heatmap.2(lfcMatrix, col = colorRampPalette(c("orangered3", "white", "slateblue3"))(31), trace="none", scale = "none", margins = c(15,15), ColSideColors = speciesToColor[colnames(lfcMatrix), "color"])
dev.off()

lfcMatrixOnlyPos <- lfcMatrix[rowSums(lfcMatrix > 0) > 0, ]
pdf(file.path(rDir, paste0("speciesContrasts_LFCmatrixOnlyPos_", LFCthreshold, "_maxFDR_", sigThreshold, ".pdf")), width = 10, height = 12) 
heatmap.2(lfcMatrixOnlyPos, col = colorRampPalette(c("orangered3", "white", "slateblue3"))(31), trace="none", scale = "none", margins = c(15,15), ColSideColors = speciesToColor[colnames(lfcMatrixOnlyPos), "color"])
dev.off()

# scaled normalized data
pdf(file.path(rDir, paste0("speciesContrasts_scaledNormData.pdf")), width = 10, height = 12) 
heatmap.2(meanData[anySigOtuSpecies,], col = f.blueblackyellow(31), trace="none", scale = "row", margins = c(15,15), ColSideColors = speciesToColor[colnames(meanData), "color"])
dev.off()


##########################################################################################
##### heatmap with within group - they are actually not as useful as one species vs the rest
withinGroups <- grep("^within_", names(diffAbunResultsDesign), value = TRUE)
sigOTUsWithinGroups <- lapply(diffAbunResultsDesign[withinGroups], function(x) x$get_significant_entries(sigThreshold, sigThreshold, LFCthresholdLRT))
anySigOtuWithinGroups <- Reduce(union, lapply(sigOTUsWithinGroups, function(x) x$any))

# scaled normalized data
pdf(file.path(rDir, paste0("withinGroupsContrasts_scaledNormData.pdf")), width = 10, height = 12) 
heatmap.2(meanData[anySigOtuWithinGroups,], col = f.blueblackyellow(31), trace="none", scale = "row", margins = c(15,15), ColSideColors = speciesToColor[colnames(meanData), "color"])
dev.off()

##########################################################################################
### The rhizobia and AMF comparisons
sigOTUsMainCont <- lapply(diffAbunResultsDesignYesNo, function(x) x$get_significant_entries(sigThreshold, sigThreshold, LFCthreshold))
anySigOtuMainCont <- Reduce(union, lapply(sigOTUsMainCont, function(x) x$any))

# scaled normalized data
pdf(file.path(rDir, paste0("mainContrasts_scaledNormData.pdf")), width = 10, height = 12) 
heatmap.2(meanData[anySigOtuMainCont,], col = f.blueblackyellow(31), trace="none", scale = "row", margins = c(15,15), ColSideColors = speciesToColor[colnames(meanData), "color"])
dev.off()

##########################################################################################
### More TSNEs
if (doTSNE) {
  forColoring <- list( default = forTest$color )
  forTSNE <- normData[anySigOtuMainCont, rownames(forTest)]
  rtsneResults <- f.run.select.plot.tSNE(unique(t(as.matrix(forTSNE))), forColoring, forTest$Pch, rDir, "tSNE_sigMainCont.svg")
  forTSNE <- normData[anySigOtuSpecies, rownames(forTest)]
  rtsneResults <- f.run.select.plot.tSNE(unique(t(as.matrix(forTSNE))), forColoring, forTest$Pch, rDir, "tSNE_sigSpeciesCont.svg")
}






