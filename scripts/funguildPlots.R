#!/usr/bin/env Rscript

dataDir <- "/media/mwschmid/myData/MWSchmid/MarcelVanDerHeijden_rhizospheres/GitIgnore_usearchOutput_ITS"
rDir <- "/media/mwschmid/myData/MWSchmid/MarcelVanDerHeijden_rhizospheres/GitIgnore_results_ITS/OTU"

##########################################################################################
### arguments from commandline
myarg <- commandArgs()
argPos <- grep("--args", myarg, fixed = TRUE)
dataDir <- as.character(myarg[argPos+1])
rDir <- as.character(myarg[argPos+2])

##########################################################################################
### libraries
suppressPackageStartupMessages({
  library("DESeq2")
  library("RColorBrewer")
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

eDir <- file.path(rDir, "enrichment")
dir.create(eDir, showWarnings = FALSE)

##########################################################################################
### load data
load(file.path(rDir, "forDownstreamAnalyses.Rdata")) # contains: seqDat (filtered and eventually rarefied), normData, sampleTab (has color and Pch), taxMod, taxModAll, keepOnlyRootSamples, removeRhizobia, useRarefaction
load(file.path(rDir, "allDataUsedInTests.Rdata")) # contains: forTest
load(file.path(rDir, "diffAbundance.Rdata")) # contains: diffAbunResultsDesign, diffAbunResultsDesignYesNo, diffAbunResultsDesignSpecies, noOutlierRefitting
otuDat <- seqDat[,rownames(forTest)]
normData <- normData[,rownames(sampleTab)]
meanData <- f.summarize.columns(normData, data.frame(sample = rownames(sampleTab), group = sampleTab$Species_name), mean)

fgAnno <- read.table(file.path(dataDir, "otus_funGuild_matched_forR.txt"), sep = '\t', stringsAsFactors = FALSE, header = TRUE, row.names = 1)
colnames(fgAnno) <- gsub("\\.", "", colnames(fgAnno))
colsToTest <- c("TrophicMode", "Guild", "Trait", "ConfidenceRanking", "GrowthMorphology")

sigThreshold <- 0.05
LFCthreshold <- 1
LFCthresholdLRT <- 0

if (length(unique(sampleTab$Type_sample)) > 1) {
  f.print.message("THIS SCRIPT ASSUMES THAT YOU USE ONLY THE ROOT SAMPLES!")
  quit("no", 1)
}

##########################################################################################
### annotation enrichment, spider plot and heatmaps
imageColors <- colorRampPalette(c("orangered3", "white", "slateblue3"))(31)
allResults <- c(diffAbunResultsDesignYesNo, diffAbunResultsDesignSpecies)
sigOTUs <- lapply(allResults, function(x) x$get_significant_entries(sigThreshold, sigThreshold, LFCthreshold))
universe <- intersect(rownames(normData), rownames(fgAnno))
for (cc in names(sigOTUs)) {
  #curMetric <- gsub("both_", "upper_", "lower", "", cc)
  curSigOtus <- sigOTUs[[cc]]
  if (length(curSigOtus$any) < 10) { next }
  curDataSig <- allResults[[cc]]$table[curSigOtus$any,]
  curDataSig <- curDataSig[intersect(rownames(curDataSig), universe),]
  #f.topLevel.plot.annotation.spider(curSigOtus, taxMod, "phylum", universe, eDir, paste0(cc, "_phylum"))
  for (curThing in colsToTest) {
    annoMap <- fgAnno[,curThing]; names(annoMap) <- rownames(fgAnno)
    allTaxum <- annoMap[universe]; length(table(allTaxum))
    allTaxumCounts <- sort(table(allTaxum), decreasing = TRUE)
    lfcByTaxum <- aggregate(curDataSig$logFC, by = list(taxum = annoMap[rownames(curDataSig)]), mean)
    rownames(lfcByTaxum) <- lfcByTaxum$taxum
    colnames(lfcByTaxum) <- c("Taxum", "averageLFC")
    write.csv(lfcByTaxum, file.path(eDir, paste0("lfcAveragedBy_", curThing, '_', cc, ".csv")), row.names = FALSE)
    # test first for significant differences in the annotation - one set compared to everything
    enriTest <- f.test.annotation.enrichment(curSigOtus$any, universe, annoMap)
    enriTestUp <- f.test.annotation.enrichment(curSigOtus$up, universe, annoMap)
    enriTestDown <- f.test.annotation.enrichment(curSigOtus$down, universe, annoMap)
    f.print.message("number significant groups:", sum(enriTest$FDR < 0.05))
    write.csv(enriTest, file.path(eDir, paste0("tax_", curThing, '_', cc,"_enrichmentAny.csv")))
    write.csv(enriTestUp, file.path(eDir, paste0("tax_", curThing, '_', cc,"_enrichmentUp.csv")))
    write.csv(enriTestDown, file.path(eDir, paste0("tax_", curThing, '_', cc,"_enrichmentDown.csv")))
    if (curThing %in% colsToTest){ # just plot all of them anyway...
      UpEnri <- rownames(enriTestUp)[(enriTestUp$FDR < 0.05) & (enriTestUp$observed > enriTestUp$expected)]
      UpDepl <- rownames(enriTestUp)[(enriTestUp$FDR < 0.05) & (enriTestUp$observed < enriTestUp$expected)]
      DownEnri <- rownames(enriTestDown)[(enriTestDown$FDR < 0.05) & (enriTestDown$observed > enriTestDown$expected)]
      DownDepl <- rownames(enriTestDown)[(enriTestDown$FDR < 0.05) & (enriTestDown$observed < enriTestDown$expected)]
      toPlot <- as.matrix(cbind(lfcByTaxum$averageLFC, lfcByTaxum$averageLFC)); colnames(toPlot) <- c("same", "label")
      rownames(toPlot) <- lfcByTaxum$Taxum
      rownames(toPlot)[rownames(toPlot) %in% UpEnri] <- paste0(rownames(toPlot)[rownames(toPlot) %in% UpEnri], "_upEnr")
      rownames(toPlot)[rownames(toPlot) %in% UpDepl] <- paste0(rownames(toPlot)[rownames(toPlot) %in% UpDepl], "_upDep")
      rownames(toPlot)[rownames(toPlot) %in% DownEnri] <- paste0(rownames(toPlot)[rownames(toPlot) %in% DownEnri], "_downEnr")
      rownames(toPlot)[rownames(toPlot) %in% DownDepl] <- paste0(rownames(toPlot)[rownames(toPlot) %in% DownDepl], "_downDep")
      toPlot <- toPlot[rowSums(is.na(toPlot)) < ncol(toPlot), ]
      if (is.null(nrow(toPlot))) {
        f.print.message("skipping", cc, "(", length(sigOTUs$any), ")", curThing)
        next
      }
      if (nrow(toPlot) < 2 | ncol(toPlot) < 2) {
        f.print.message("skipping", cc, "(", length(sigOTUs$any), ")", curThing)
        next
      }
      f.open.figure(eDir, paste0("heatmap_", curThing, '_', cc, "_LFCs.svg"), TRUE, width = 8, height = 7)
      par(oma = c(5, 0, 0, 15))
      heatmap.2(toPlot, col = imageColors, scale = "none", trace = "none", Rowv=FALSE, Colv=FALSE, dendrogram="none") # keep the order of the individual contrasts
      f.close.figure()
    }
  }
}

##########################################################################################
### collect results
allFiles <- list.files(eDir, "^tax_")
#curFileName <- allFiles[1]
allRes <- list()
resCounter <- 0
for (curFileName in allFiles) {
  curParts <- gsub("^tax_|\\.csv", "", curFileName)
  curParts <- unlist(strsplit(curParts, "_"))
  curThing <- curParts[1]
  contrast <- paste0(curParts[2:(length(curParts)-1)], collapse = '_')
  enrichType <- curParts[length(curParts)]
  temp <- read.csv(file.path(eDir, curFileName), header = TRUE, stringsAsFactors = FALSE, row.names = 1)
  toKeep <- rownames(temp)[temp$pValue < 0.05]
  if (length(toKeep) == 0) {next}
  for (tk in toKeep) {
    resCounter <- resCounter+1
    allRes[[resCounter]] <- c(curThing, contrast, enrichType, tk, unlist(temp[tk, ]))
  }
}
allRes <- as.data.frame(do.call("rbind", allRes), stringsAsFactors = FALSE)
colnames(allRes) <- c("taxLevel", "contrast", "enrichType", "taxon", "observed", "expected", "pValue", "FDR_separate")
isDepleted <- as.numeric(allRes$observed) < as.numeric(allRes$expected)
allRes$enrichType[isDepleted] <- gsub("enrichment", "depletion", allRes$enrichType[isDepleted])
write.csv(allRes, file.path(eDir, "summary_onlySignificant.csv"), row.names = FALSE)












