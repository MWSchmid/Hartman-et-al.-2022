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

# check for Kyle
#temp <- as.data.frame(t(meanData[c("zot918", "zot6280", "zot67", "zot499", "zot428", "zot1494"),]))
#narf <- unique(sampleTab[,c("Species_name", "AMF", "Rhizobia")])
#rownames(narf) <- narf$Species_name
#temp$AMFRHIZO <- paste0(narf[rownames(temp),"AMF"], '_', narf[rownames(temp),"Rhizobia"])
#hurz <- aggregate(temp[,c("zot918", "zot6280", "zot67", "zot499", "zot428", "zot1494")], by = list(type=temp$AMFRHIZO), mean)
#write.csv(temp, "/home/marc/Downloads/forKyle_normAndAve.csv")
#write.csv(hurz, "/home/marc/Downloads/forKyle_normAndAve_aggr.csv", row.names = FALSE)

if (grepl("ITS", rDir)) {
  dataDir <- "/media/mwschmid/myData/MWSchmid/MarcelVanDerHeijden_rhizospheres/GitIgnore_usearchOutput_ITS"
  fgAnno <- read.table(file.path(dataDir, "otus_funGuild_matched_forR.txt"), sep = '\t', stringsAsFactors = FALSE, header = TRUE, row.names = 1)
  colnames(fgAnno) <- gsub("\\.", "", colnames(fgAnno))
  colsToTest <- c("TrophicMode", "Guild", "Trait", "ConfidenceRanking", "GrowthMorphology")
}

sigThreshold <- 0.05
LFCthreshold <- 1
LFCthresholdLRT <- 0
minSigCompsForOneToOnePlot <- 9

if (length(unique(sampleTab$Type_sample)) > 1) {
  f.print.message("THIS SCRIPT ASSUMES THAT YOU USE ONLY THE ROOT SAMPLES!")
  quit("no", 1)
}

#apply(seqDatRarefied, 2, max)
#colSums(seqDatRarefied)
#narf <- rep(0, ncol(seqDatRarefied)); names(narf) <- colnames(seqDatRarefied)
#for (curCol in colnames(seqDatRarefied)) {
#  temp <- seqDatRarefied[,curCol]
#  names(temp) <- rownames(seqDatRarefied)
#  temp <- sort(temp, decreasing = TRUE)
#  narf[curCol] <- round(sum(temp[1:10])/sum(temp)*100, 2)
#}
#mean(narf)
#range(narf)
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

### get "core microbiomes"
sigOTUsCoreMicrobiome <- list()
sigOTUsCoreMicrobiome$nn <- intersect(sigOTUs$AMFwithinNonRhizo$down, sigOTUs$bothVsNone$down)
sigOTUsCoreMicrobiome$An <- intersect(sigOTUs$AMFwithinNonRhizo$up, sigOTUs$RhizobiaWithinAMF$down)
sigOTUsCoreMicrobiome$AR <- intersect(sigOTUs$RhizobiaWithinAMF$up, sigOTUs$bothVsNone$up)
#sigOTUsCoreMicrobiome$AMF <- intersect(sigOTUs$AMFwithinNonRhizo$up, sigOTUs$bothVsNone$up)
#checkAn <- intersect(c("zot918", "zot6280", "zot67", "zot499", "zot428", "zot1494"), sigOTUsCoreMicrobiome$An)
#diffAbunResultsDesignYesNo$AMFwithinNonRhizo$table[checkAn, "logFC"]
#diffAbunResultsDesignYesNo$RhizobiaWithinAMF$table[checkAn, "logFC"]
#checkAR <- intersect(c("zot918", "zot6280", "zot67", "zot499", "zot428", "zot1494"), sigOTUsCoreMicrobiome$AR)
#diffAbunResultsDesignYesNo$RhizobiaWithinAMF$table[checkAR, "logFC"]
#diffAbunResultsDesignYesNo$bothVsNone$table[checkAR, "logFC"]
#intersect(c("zot918", "zot6280", "zot67", "zot499", "zot428", "zot1494"), sigOTUsCoreMicrobiome$AR)
numSigOtusCoreMicrobiome <- cbind(names(sigOTUsCoreMicrobiome), unlist(lapply(sigOTUsCoreMicrobiome, length)))
colnames(numSigOtusCoreMicrobiome) <- c("symboGroup", "numberOfOtus")
write.csv(numSigOtusCoreMicrobiome, file.path(rDir, "numSigOtusCoreMicrobiome.csv"), quote = FALSE, row.names = FALSE)
for (curGroup in names(sigOTUsCoreMicrobiome)) {
  out <- taxMod[sigOTUsCoreMicrobiome[[curGroup]],]
  if (grepl("ITS", rDir)) {
    for (cc in colsToTest) {
      out[,cc] <- fgAnno[sigOTUsCoreMicrobiome[[curGroup]],cc]
    }
  }
  write.csv(out, file.path(rDir, paste0("coreMicrobiome_", curGroup, "_confidentAnnotation.csv")), quote = FALSE)
  out <- taxModAll[sigOTUsCoreMicrobiome[[curGroup]],]
  if (grepl("ITS", rDir)) {
    for (cc in colsToTest) {
      out[,cc] <- fgAnno[sigOTUsCoreMicrobiome[[curGroup]],cc]
    }
  }
  write.csv(out, file.path(rDir, paste0("coreMicrobiome_", curGroup, "_entireAnnotation.csv")), quote = FALSE)
}

#anySigCore <- unique(unlist(sigOTUsCoreMicrobiome))
#sigOtuCounts <- rowSums(seqDat[anySigCore,])
#nonSigOtuCounts <- rowSums(seqDat[setdiff(rownames(seqDat), anySigCore),])
#summary(sigOtuCounts)
#summary(nonSigOtuCounts)

### similar for the 1-to-1 comparisons
allSpecies <- names(speciesFDRmatrices)
numberSigComparisons <- length(allSpecies)-1
highlySpecificOtus <- matrix(0, nrow = length(allSpecies), ncol = numberSigComparisons, dimnames = list(allSpecies, paste0("minSigComps_", 1:numberSigComparisons)))
for (curSpecies in allSpecies) {
  temp <- rowSums((speciesFDRmatrices[[curSpecies]] < 0.05) & (speciesLFCmatrices[[curSpecies]] > 0), na.rm = TRUE)
  highlySpecificOtus[curSpecies, ] <- sapply(1:numberSigComparisons, function(x) sum(temp >= x))
}
write.csv(highlySpecificOtus, file.path(rDir, "numSigOtusOneToOne.csv"), quote = FALSE)

anySigOtuWithLRTandOneToOne <- union(anySigOtuWithLRT, Reduce(union, lapply(speciesFDRmatrices, rownames)))

f.print.message("found", length(anySigOtu), "significant OTUs, with LRT included:", length(anySigOtuWithLRT), " and with one-to-one:", length(anySigOtuWithLRTandOneToOne))

##########################################################################################
##### heatmap/image with the numSigOtusOneToOne.csv
orderedSpecies <- c("Lupin", "Brome", "Petunia", "Tobacco", "Tomato", "Wheat", "Maize", "Alfalfa", "Lotus", "Pea", "Medicago", "Trifolium", "Arabidopsis", "Brassica", "Cardamine", "Spinach", "Sugarbeet")
speciesToColor <- unique(sampleTab[,c("Species_name", "color")]); rownames(speciesToColor) <- speciesToColor$Species_name

f.print.message("at least three significant comparisons")
highlySpecificOtus <- highlySpecificOtus[,3:ncol(highlySpecificOtus)]
highlySpecificOtus <- highlySpecificOtus[orderedSpecies,]
f.open.figure(rDir, "numSigOtusOneToOne_image.svg", TRUE, height = 4, width = 4)
par(mar = c(5,7,0,0))
xChars <- gsub("minSigComps_", "", colnames(highlySpecificOtus))
yChars <- rownames(highlySpecificOtus)
xPos <- 1:ncol(highlySpecificOtus)
yPos <- 1:nrow(highlySpecificOtus)
xLabs <- "minimal number of significant comparisons"
yLabs <- ""#"plant species"
toPlot <- t(highlySpecificOtus)
image(xPos, yPos, toPlot, xlab = xLabs, ylab = yLabs, main = "", yaxt = "n", xaxt = "n")
axis(1, at = xPos, labels = xChars, outer = FALSE)
axis(2, at = yPos, labels = yChars, outer = FALSE, las = 1)
text(rep(xPos, each = length(yPos)), yPos, t(toPlot), cex = 0.7)
f.close.figure()

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
  forTSNE <- normData[anySigOtuWithLRTandOneToOne, rownames(forTest)]
  rtsneResults <- f.run.select.plot.tSNE(unique(t(as.matrix(forTSNE))), forColoring, forTest$Pch, rDir, "tSNE_sigAnyContWithLRTandOneToOne.svg")
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
  forTSNE <- mat[anySigOtuWithLRTandOneToOne, rownames(forTest)]
  rtsneResults <- f.run.select.plot.tSNE(unique(t(as.matrix(forTSNE))), forColoring, forTest$Pch, rDir, "tSNE_vstCorrected_sigAnyContWithLRTandOneToOne.svg")
}

##########################################################################################
##### heatmap with the species LFCs
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
##### heatmaps one to one
for (curSpecies in allSpecies) {
  lfcMatrix <- speciesLFCmatrices[[curSpecies]]
  lfcMatrix <- lfcMatrix[, orderedSpecies[orderedSpecies %in% colnames(lfcMatrix)]]
  pdf(file.path(rDir, paste0("oneToOne_LFCmatrix_", curSpecies, ".pdf")), width = 10, height = 12) 
  heatmap.2(lfcMatrix, col = colorRampPalette(c("orangered3", "white", "slateblue3"))(31), trace="none", scale = "none", margins = c(15,15), dendrogram = "row", Colv = FALSE,ColSideColors = speciesToColor[colnames(lfcMatrix), "color"])
  #heatmap.2(lfcMatrix, col = colorRampPalette(c("orangered3", "white", "slateblue3"))(31), trace="none", scale = "none", margins = c(15,15), ColSideColors = speciesToColor[colnames(lfcMatrix), "color"])
  dev.off()
}

##########################################################################################
### Plot the species-specific OTUs
sigOtusForOneToOnePlot <- c()
for (curSpecies in allSpecies) {
  curSigOtus <- rownames(speciesFDRmatrices[[curSpecies]])[rowSums((speciesFDRmatrices[[curSpecies]] < 0.05) & (speciesLFCmatrices[[curSpecies]] > 0), na.rm = TRUE) >= minSigCompsForOneToOnePlot]
  sigOtusForOneToOnePlot <- union(sigOtusForOneToOnePlot, curSigOtus)
}
pdf(file.path(rDir, paste0("speciesContrasts_minNumSig_", minSigCompsForOneToOnePlot, "_scaledNormData.pdf")), width = 10, height = 12) 
heatmap.2(meanData[sigOtusForOneToOnePlot,orderedSpecies], col = f.blueblackyellow(31), trace="none", scale = "row", margins = c(15,15), ColSideColors = speciesToColor[orderedSpecies, "color"], Colv = FALSE, dendrogram = "row")
dev.off()

##########################################################################################
### Core microbiomes
anyCoreMicrobiome <- Reduce(union, sigOTUsCoreMicrobiome)

# scaled normalized data
forOrder <- rev(c("Lupin", "Brome", "Petunia", "Tobacco", "Tomato", "Wheat", "Maize", "Alfalfa", "Lotus", "Pea", "Medicago", "Trifolium", "Arabidopsis", "Brassica", "Cardamine", "Spinach", "Sugarbeet")) # because of the color
pdf(file.path(rDir, paste0("coreMicrobiome_scaledNormData.pdf")), width = 10, height = 12) 
heatmap.2(meanData[anyCoreMicrobiome,forOrder], col = f.blueblackyellow(31), trace="none", scale = "row", margins = c(15,15), ColSideColors = speciesToColor[forOrder, "color"], Colv = FALSE, dendrogram = "row")
dev.off()

for (curCore in names(sigOTUsCoreMicrobiome)) {
  curSig <- sigOTUsCoreMicrobiome[[curCore]]
  pdf(file.path(rDir, paste0("coreMicrobiome_", curCore, "_scaledNormData.pdf")), width = 10, height = 12) 
  heatmap.2(meanData[curSig, forOrder], col = f.blueblackyellow(31), trace="none", scale = "row", margins = c(15,15), ColSideColors = speciesToColor[forOrder, "color"], Colv = FALSE, dendrogram = "row")
  dev.off()
}

##########################################################################################
### The rhizobia and AMF comparisons
sigOTUsMainCont <- lapply(diffAbunResultsDesignYesNo, function(x) x$get_significant_entries(sigThreshold, sigThreshold, LFCthreshold))
anySigOtuMainCont <- Reduce(union, lapply(sigOTUsMainCont, function(x) x$any))
write.table(anySigOtuMainCont, file.path(rDir, "mainContrastsSigOtus.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)

# scaled normalized data
pdf(file.path(rDir, paste0("mainContrasts_scaledNormData.pdf")), width = 10, height = 12) 
heatmap.2(meanData[anySigOtuMainCont,forOrder], col = f.blueblackyellow(31), trace="none", scale = "row", margins = c(15,15), ColSideColors = speciesToColor[forOrder, "color"], Colv = FALSE, dendrogram = "row")
dev.off()

for (curCont in names(sigOTUsMainCont)){
  curSig <- sigOTUsMainCont[[curCont]]$any
  pdf(file.path(rDir, paste0("mainContrasts_", curCont, "_scaledNormData.pdf")), width = 10, height = 12) 
  heatmap.2(meanData[curSig,forOrder], col = f.blueblackyellow(31), trace="none", scale = "row", margins = c(15,15), ColSideColors = speciesToColor[forOrder, "color"], Colv = FALSE, dendrogram = "row")
  dev.off()
  f.print.message(curCont, "has", length(curSig), "entries.")
}

##########################################################################################
### The rhizobia and AMF comparisons, only pairwise contrasts
sigOTUsPairwiseCont <- lapply(diffAbunResultsDesignYesNo[c("RhizobiaWithinAMF", "AMFwithinNonRhizo", "bothVsNone")], function(x) x$get_significant_entries(sigThreshold, sigThreshold, LFCthreshold))
anySigOtuPairwiseCont <- Reduce(union, lapply(sigOTUsPairwiseCont, function(x) x$any))
f.print.message("pairwise contrasts have:", length(anySigOtuPairwiseCont))

# scaled normalized data
pdf(file.path(rDir, paste0("pairwiseContrasts_scaledNormData.pdf")), width = 10, height = 12) 
heatmap.2(meanData[anySigOtuMainCont,forOrder], col = f.blueblackyellow(31), trace="none", scale = "row", margins = c(15,15), ColSideColors = speciesToColor[forOrder, "color"], Colv = FALSE, dendrogram = "row")
dev.off()

##########################################################################################
### The rhizobia and AMF comparisons, scale for pot
f.get.residuals.after.pot <- function(x, sampleTab) {
  out <- resid(lm(x ~ sampleTab$Pot))
  return(out)
}
sampleTab <- sampleTab[colnames(normData),]
normRes <- t(apply(normData, 1, function(x) f.get.residuals.after.pot(x, sampleTab)))
normResMean <- f.summarize.columns(normRes, data.frame(sample = rownames(sampleTab), group = sampleTab$Species_name, stringsAsFactors = FALSE), mean)

if (max(abs(range(normResMean))) > 5) {
  f.print.message("limiting to a residual of 5")
  normResMean[normResMean < (-5)] <- (-5)
  normResMean[normResMean > 5] <- 5
}

# residuals
pdf(file.path(rDir, paste0("mainContrasts_residuals.pdf")), width = 10, height = 12) 
#heatmap.2(meanData[anySigOtuMainCont,], col = f.blueblackyellow(31), trace="none", scale = "row", margins = c(15,15), ColSideColors = speciesToColor[colnames(meanData), "color"])
heatmap.2(normResMean[anySigOtuMainCont,forOrder], col = f.blueblackyellow(31), trace="none", scale = "none", margins = c(15,15), ColSideColors = speciesToColor[forOrder, "color"], Colv = FALSE, dendrogram = "row")
dev.off()

f.print.message(length(anySigOtuMainCont), "OTUs in the residual heatmap.")

##########################################################################################
### More TSNEs
if (doTSNE) {
  forColoring <- list( default = forTest$color )
  forTSNE <- normData[anySigOtuMainCont, rownames(forTest)]
  rtsneResults <- f.run.select.plot.tSNE(unique(t(as.matrix(forTSNE))), forColoring, forTest$Pch, rDir, "tSNE_sigMainCont.svg")
  forTSNE <- normData[anySigOtuSpecies, rownames(forTest)]
  rtsneResults <- f.run.select.plot.tSNE(unique(t(as.matrix(forTSNE))), forColoring, forTest$Pch, rDir, "tSNE_sigSpeciesCont.svg")
}

#names(diffAbunResultsDesign)
#names(diffAbunResultsDesignYesNo)
#names(diffAbunResultsDesignSpecies)
#names(diffAbunResultsOneToOne)
