#!/usr/bin/env Rscript

rDir <- "/media/mwschmid/myData/MWSchmid/MarcelVanDerHeijden_rhizospheres/GitIgnore_results_16S_semiRarefied_noRhizobia/OTU"
rDir <- "/media/mwschmid/myData/MWSchmid/MarcelVanDerHeijden_rhizospheres/GitIgnore_results_ITS_semiRarefied_noGlomeromycota/OTU"

##########################################################################################
### arguments from commandline
myarg <- commandArgs()
argPos <- grep("--args", myarg, fixed = TRUE)
rDir <- as.character(myarg[argPos+1])

##########################################################################################
### libraries
suppressPackageStartupMessages({
  options(java.parameters = "-Xmx12g")
  library("vegan")
  library("phangorn")
  library("RColorBrewer")
  library("GUniFrac")
  library("parallel")
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
load(file.path(rDir, "allDataUsedInTests.Rdata")) # forTest

if (length(unique(sampleTab$Type_sample)) > 1) {
  f.print.message("THIS SCRIPT ASSUMES THAT YOU USE ONLY THE ROOT SAMPLES!")
  quit("no", 1)
}

##########################################################################################
### plot MSR/MER/MEV
forOrder <- rev(c("Lupin", "Alfalfa", "Lotus", "Pea", "Medicago", "Trifolium", "Brome", "Petunia", "Tobacco", "Tomato", "Wheat", "Maize", "Arabidopsis", "Brassica", "Cardamine", "Spinach", "Sugarbeet")) # because of the color
speciesColors <- unique(forTest[,c("SpeciesChar", "color")]); rownames(speciesColors) <- speciesColors$SpeciesChar
toSummarize <- c("MSR", "MER", "MEV")
groupsToUse <- c("SymboGroup", "between_AMF_nonRhizo_families", "between_nonAMF_nonRhizo_families", "within_nonAMF_nonRhizo", "within_AMF_nonRhizo", "within_AMF_Rhizo")
dontShow <- c("others", "notOfInterest")
for (curGroup in groupsToUse) { forTest[[curGroup]] <- as.character(forTest[[curGroup]]) }
forWidth <- unlist(lapply(groupsToUse, function(x) length(setdiff(unique(forTest[[x]]), dontShow)))); names(forWidth) <- groupsToUse
for (toSum in toSummarize) {
  temp <- lapply(groupsToUse, function(x) aggregate(forTest[[toSum]], by = list(group = forTest[[x]]), mean)); names(temp) <- groupsToUse
  groupsToPlot <- setdiff(Reduce(union, lapply(temp, function(x) x$group)), dontShow)
  forYaxisLabels <- pretty(forTest[[toSum]])
  f.open.figure(rDir, paste0("boxplots_", toSum, ".svg"), useSVG = TRUE, height = 8, width = 2+sum(forWidth))
  par(oma = c(5, 5, 0, 0), mar = c(0, 0, 0, 0))
  layout(matrix(1:length(groupsToUse), nrow = 1), widths = forWidth)
  for (curGroup in groupsToUse) {
    plot(NA, main = "", bty = "n", xaxs = "r", yaxs = "r", xlab = "", ylab = "", las = 1, cex = 0.8, tck = 0.01, xlim = c(0.5, forWidth[curGroup]+0.5), ylim = range(forYaxisLabels), xaxt = "n", yaxt = "n")
    curPos <- 1
    subGroups <- setdiff(unique(forTest[[curGroup]]), dontShow)
    if (sum(subGroups %in% forOrder) == length(subGroups)) {
      subGroups <- forOrder[forOrder %in% subGroups]
      curCols <- speciesColors[subGroups, "color"]
    } else {
      if (curGroup == "SymboGroup") {
        subGroups <- c("nonAMF_nonRhizo", "AMF_nonRhizo", "AMF_Rhizo", "nonAMF_Rhizo")
      }
      curCols <- rep("white", length(subGroups))
    }
    names(curCols) <- subGroups
    for (curSubGroup in subGroups) {
      toPlot <- forTest[forTest[[curGroup]] == curSubGroup,toSum]
      boxplot(toPlot, names = c(curGroup), ylim = range(forYaxisLabels), drawRect = TRUE, add = TRUE, at = curPos, xaxt = "n", yaxt = "n", col = curCols[curSubGroup])
      points(curPos+(1-jitter(rep(1, length(toPlot)), 10)), toPlot, cex = 3.3, col = "black", pch = 16)#f.add.alpha("#FF00A9", 0.3)
      #curMean <- mean(toPlot)
      #lines(c(curPos-0.3,curPos+0.3), c(curMean, curMean), col = "black", lwd = 4, lty = 1)
      curPos <- curPos + 1
    }
    if (curGroup == "SymboGroup") { axis(2, at = forYaxisLabels, labels = forYaxisLabels, outer = TRUE, las = 1, line=2, lwd=2, cex.axis=3) }
    axis(1, at = 1:length(subGroups), labels = subGroups, outer = TRUE, las = 2, line=2, lwd=2, cex.axis=3)
  }
  f.close.figure()
}
