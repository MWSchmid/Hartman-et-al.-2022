#!/usr/bin/env Rscript

rm(list=ls())
rDir <- "/media/mwschmid/myData/MWSchmid/MarcelVanDerHeijden_rhizospheres/GitIgnore_results_16S_semiRarefied/OTU"
rDir <- "/media/mwschmid/myData/MWSchmid/MarcelVanDerHeijden_rhizospheres/GitIgnore_results_ITS_semiRarefied/OTU"
rm(list=ls())

##########################################################################################
### arguments from commandline
myarg <- commandArgs()
argPos <- grep("--args", myarg, fixed = TRUE)
rDir <- as.character(myarg[argPos+1])

##########################################################################################
### libraries
suppressPackageStartupMessages({
  options(java.parameters = "-Xmx12g")
  library("RColorBrewer")
  library("grid")
  library("Ternary")
  source("/media/mwschmid/myData/MWSchmid/Development/R/RNAseqWrapper.R")
  source("/media/mwschmid/myData/MWSchmid/Development/R/RNAseqWrapperAds.R")
  #source("/media/mwschmid/myData/MWSchmid/MarcelVanDerHeijden_rhizospheres/scripts/function_ternary_plot.R")
  #source("/media/mwschmid/myData/MWSchmid/Development/R/lm_wrapper.R")
  #source("/media/mwschmid/myData/MWSchmid/Development/R/OTU_functions.R")
  #source("/media/mwschmid/myData/MWSchmid/MarcelVanDerHeijden_rhizospheres/scripts/final_createAnalysisSets.R")
})

f.print.message <- function(...) { cat("=== ", format(Sys.time(), "%Y %b %d %X"), paste("===", ...,"\n")) }

pcname <- system('uname -n',intern=T)
numCoresPerMachine <- c("nuke" = 4, "styx" = 15, "marc-IEU" = 7, "piftel" = 30)
numCoresAvailable <- numCoresPerMachine[pcname]

##########################################################################################
### load data
load(file.path(rDir, "forDownstreamAnalyses.Rdata")) # contains: seqDat (filtered), seqDatRarefied (filtered and rarefied), normData (based on seqDat), sampleTab (has color and Pch), taxMod, taxModAll, otusToRemove, removeRhizobia, removeGlomeromycota
load(file.path(rDir, "allDataUsedInTests.Rdata")) # forTest

coreMicrobiomesFullTables <- list(
  An = read.csv(file.path(rDir, "coreMicrobiome_An_entireAnnotation.csv"), stringsAsFactors = FALSE, row.names = 1, header = TRUE),
  AR = read.csv(file.path(rDir, "coreMicrobiome_AR_entireAnnotation.csv"), stringsAsFactors = FALSE, row.names = 1, header = TRUE),
  nn = read.csv(file.path(rDir, "coreMicrobiome_nn_entireAnnotation.csv"), stringsAsFactors = FALSE, row.names = 1, header = TRUE)
  #An = read.csv(file.path(rDir, "coreMicrobiome_An_confidentAnnotation.csv"), stringsAsFactors = FALSE, row.names = 1, header = TRUE),
  #AR = read.csv(file.path(rDir, "coreMicrobiome_AR_confidentAnnotation.csv"), stringsAsFactors = FALSE, row.names = 1, header = TRUE),
  #nn = read.csv(file.path(rDir, "coreMicrobiome_nn_confidentAnnotation.csv"), stringsAsFactors = FALSE, row.names = 1, header = TRUE)
)
coreMicrobiomes <- lapply(coreMicrobiomesFullTables, function(x) rownames(x))
coreMic <- do.call("rbind", coreMicrobiomesFullTables)
if (grepl("results_ITS_", rDir)) {
  levelToUse <- "class"
  sortedClasses <- sort(unique(coreMic[[levelToUse]]))
} else {
  levelToUse <- "class"
  coreMic$class[grepl("^Acidobacteria", coreMic$class)] <- "Acidobacteria"
  coreMic$class[grepl("^undef_undef_undef_Parcubacteria", coreMic$class)] <- "Parcubacteria"
  coreMic$class[grepl("^undef_undef_undef_Saccharibacteria", coreMic$class)] <- "Saccharibacteria"
  sortedClasses <- sort(unique(coreMic[[levelToUse]]))
  proteoBac <- c("Alphaproteobacteria", "Betaproteobacteria", "Gammaproteobacteria", "Deltaproteobacteria")
  sortedClasses <- c(proteoBac[proteoBac %in% sortedClasses], setdiff(sortedClasses, proteoBac))
}
rownames(coreMic) <- gsub("[[:alpha:]]{2}\\.", "", rownames(coreMic))
phylumColors <- f.blackblueyellowredpinkNICE(length(sortedClasses))
names(phylumColors) <- sortedClasses
useCustomColors <- FALSE
if (useCustomColors){
  if (grepl("results_ITS_", rDir)) {
    phylumColors <- c(
      "ukn" = "black",
      "Ascomycota" = "cyan",
      "Basidiomycota" = "magenta",
      "Olpidiomycota" = "yellow",
      "Glomeromycota" = "darkgreen"
    )
  } else {
    phylumColors <- c(
      "ukn" = "black",
      "Proteobacteria" = "darkgreen",
      "Actinobacteria" = "black",
      "Bacteroidetes" = "black",
      "Firmicutes" = "black",
      "Acidobacteria" = "black",
      "Candidatus_Saccharibacteria" = "black",
      "Armatimonadetes" = "black",
      "Parcubacteria" = "black",
      "Planctomycetes" = "black"
    )
  }
}
otuColorPhylum <- phylumColors[coreMic[[levelToUse]]]; names(otuColorPhylum) <- rownames(coreMic)

if (length(unique(sampleTab$Type_sample)) > 1) {
  f.print.message("THIS SCRIPT ASSUMES THAT YOU USE ONLY THE ROOT SAMPLES!")
  quit("no", 1)
}

##########################################################################################
### keep only the three main conditions
translator <- c("AMF_nonRhizo" = "An", "nonAMF_nonRhizo" = "nn", "nonAMF_Rhizo" = "nR", "AMF_Rhizo" = "AR")
sampleTab$simpleGroup <- translator[paste0(sampleTab$AMF, '_', sampleTab$Rhizobia)]
sampleTab <- subset(sampleTab, simpleGroup != "nR")
normData <- normData[,rownames(sampleTab)]

# average by species, filter for a minimal abundance!, then summarize by type
aveBySpecies <- f.summarize.columns(normData, data.frame(sample = rownames(sampleTab), group = sampleTab$Species_name), mean)
#aveBySpecies <- aveBySpecies[rowSums(aveBySpecies>=4)>0,]
forMoreAve <- unique(sampleTab[,c("Species_name", "simpleGroup")])
aveByGroupNoFilt <- f.summarize.columns(aveBySpecies, data.frame(sample = forMoreAve$Species_name, group = forMoreAve$simpleGroup), mean)
aveByGroup <- aveByGroupNoFilt[rowSums(aveByGroupNoFilt>=1)>0,]

# select colors
otuColor <- rep("grey", nrow(aveByGroup)); names(otuColor) <- rownames(aveByGroup)
otuColor[names(otuColor) %in% coreMicrobiomes$AR] <- "#1b9e77" # green
otuColor[names(otuColor) %in% coreMicrobiomes$nn] <- "#d95f02" # red
otuColor[names(otuColor) %in% coreMicrobiomes$An] <- "#7570b3" # blue

# select scaling
otuProportionPercent <- rowSums(aveByGroup)/sum(aveByGroup)*100
minCex <- 1
otuCex <- minCex + otuProportionPercent
names(otuCex) <- rownames(aveByGroup)

pdf(file.path(rDir, "coreMicrobiome_ternaryPlotWithAll.pdf"), width = 10, height = 10)
# empty plot
TernaryPlot(atip = colnames(aveByGroup)[1],
            btip = colnames(aveByGroup)[2],
            ctip = colnames(aveByGroup)[3],
            # no minor grid lines
            grid.minor.lines = 0
            )

# add grey points
mask <- grepl("grey|gray", otuColor)
TernaryPoints(aveByGroup[mask,],
              cex = otuCex[mask],
              col = otuColor[mask],
              pch = 16
)

# add colored points
TernaryPoints(aveByGroup[!mask,],
              cex = otuCex[!mask],
              col = otuColor[!mask],
              pch = 16
)
invisible(dev.off())


pdf(file.path(rDir, "coreMicrobiome_ternaryPlot.pdf"), width = 10, height = 10)
# empty plot
TernaryPlot(atip = colnames(aveByGroup)[1],
            btip = colnames(aveByGroup)[2],
            ctip = colnames(aveByGroup)[3],
            # no minor grid lines
            grid.minor.lines = 0
)

# add colored points
mask <- grepl("grey|gray", otuColor)
TernaryPoints(aveByGroup[!mask,],
              cex = otuCex[!mask],
              col = otuColor[!mask],
              pch = 16
)
invisible(dev.off())


isInPlot <- intersect(rownames(aveByGroupNoFilt), rownames(coreMic))
otuProportionPercentNoFilt <- rowSums(aveByGroupNoFilt)/sum(aveByGroupNoFilt)*100
otuCexNoFilt <- minCex + otuProportionPercentNoFilt
names(otuCexNoFilt) <- rownames(aveByGroupNoFilt)
pdf(file.path(rDir, "coreMicrobiome_ternaryPlotColByPhylum.pdf"), width = 10, height = 10)
TernaryPlot(atip = colnames(aveByGroup)[1],
            btip = colnames(aveByGroup)[2],
            ctip = colnames(aveByGroup)[3],
            # no minor grid lines
            grid.minor.lines = 0
)
TernaryPoints(aveByGroupNoFilt[isInPlot,],
              cex = otuCexNoFilt[isInPlot],
              col = otuColorPhylum[isInPlot],
              pch = 16
)
legend('topright', 
       legend = names(phylumColors),
       cex = 0.8, bty = 'n', pch = 21, pt.cex = 1.8,
       pt.bg = phylumColors)
invisible(dev.off())

isInPlot <- intersect(rownames(aveByGroupNoFilt), rownames(coreMic))
otuProportionPercentNoFilt <- rowSums(aveByGroupNoFilt)/sum(aveByGroupNoFilt)*100
otuCexNoFilt <- minCex + otuProportionPercentNoFilt
names(otuCexNoFilt) <- rownames(aveByGroupNoFilt)
pdf(file.path(rDir, "coreMicrobiome_ternaryPlotColByPhylum_withGrey.pdf"), width = 10, height = 10)
TernaryPlot(atip = colnames(aveByGroup)[1],
            btip = colnames(aveByGroup)[2],
            ctip = colnames(aveByGroup)[3],
            # no minor grid lines
            grid.minor.lines = 0
)
# add grey points
mask <- names(otuColor)[grepl("grey|gray", otuColor)]
TernaryPoints(aveByGroup[mask,],
              cex = otuCex[mask],
              col = otuColor[mask],
              pch = 16
)

# add colored points
TernaryPoints(aveByGroupNoFilt[isInPlot,],
              cex = otuCexNoFilt[isInPlot],
              col = otuColorPhylum[isInPlot],
              pch = 16
)

legend('topright', 
       legend = names(phylumColors),
       cex = 0.8, bty = 'n', pch = 21, pt.cex = 1.8,
       pt.bg = phylumColors)
invisible(dev.off())

##########################################################################################
# some checks
quit("no", 0)
otuProps <- aveByGroup/rowSums(aveByGroup)
checkThese <- rownames(otuProps)[otuProps[,"An"] > 0.9]
aveByGroup[checkThese,]





