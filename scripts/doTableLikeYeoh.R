#!/usr/bin/env Rscript

rDir <- "/media/mwschmid/myData/MWSchmid/MarcelVanDerHeijden_rhizospheres/GitIgnore_results_16S_semiRarefied/OTU"
rDir <- "/media/mwschmid/myData/MWSchmid/MarcelVanDerHeijden_rhizospheres/GitIgnore_results_16S_semiRarefied_noRhizobia/OTU"

rDir <- "/media/mwschmid/myData/MWSchmid/MarcelVanDerHeijden_rhizospheres/GitIgnore_results_ITS_semiRarefied/OTU"
rDir <- "/media/mwschmid/myData/MWSchmid/MarcelVanDerHeijden_rhizospheres/GitIgnore_results_ITS_semiRarefied_noGlomeromycota/OTU"

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
meanData <- f.summarize.columns(normData, data.frame(sample = rownames(sampleTab), group = sampleTab$Species_name), mean)

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

if (grepl("ITS", rDir)) {
  taxaToUse <- c("phylum", "class", "order", "family", "genus")
} else {
  taxaToUse <- c("phylum", "class", "order")
}


##########################################################################################
### get all significant OTUs
sigOTUsLRT <- lapply(diffAbunResultsDesign, function(x) x$get_significant_entries(sigThreshold, sigThreshold, LFCthresholdLRT))
sigOTUs <- lapply(c(diffAbunResultsDesignYesNo, diffAbunResultsDesignSpecies), function(x) x$get_significant_entries(sigThreshold, sigThreshold, LFCthreshold))

### get "core microbiomes"
sigOTUsCoreMicrobiome <- list()
sigOTUsCoreMicrobiome$nn <- intersect(sigOTUs$AMFwithinNonRhizo$down, sigOTUs$bothVsNone$down)
sigOTUsCoreMicrobiome$An <- intersect(sigOTUs$AMFwithinNonRhizo$up, sigOTUs$RhizobiaWithinAMF$down)
sigOTUsCoreMicrobiome$AR <- intersect(sigOTUs$RhizobiaWithinAMF$up, sigOTUs$bothVsNone$up)
allOTUs <- Reduce(union, sigOTUsCoreMicrobiome)
#print(lapply(sigOTUsCoreMicrobiome, length))

##########################################################################################
sigOTUsCoreMicrobiomeAnnotationConfident <- list()
sigOTUsCoreMicrobiomeAnnotationEntire <- list()
for (curGroup in names(sigOTUsCoreMicrobiome)) {
  sigOTUsCoreMicrobiomeAnnotationConfident[[curGroup]] <- table(apply(taxMod[sigOTUsCoreMicrobiome[[curGroup]],taxaToUse], 1, function(x) paste(x, collapse = '|')))#
  sigOTUsCoreMicrobiomeAnnotationEntire[[curGroup]] <- table(apply(taxModAll[sigOTUsCoreMicrobiome[[curGroup]],taxaToUse], 1, function(x) paste(x, collapse = '|')))#, "family", "genus"
}

allTaxaStrings <- sort(union(Reduce(union, lapply(sigOTUsCoreMicrobiomeAnnotationConfident, function(x) names(x))),
                        Reduce(union, lapply(sigOTUsCoreMicrobiomeAnnotationEntire, function(x) names(x)))))

##########################################################################################
### table like Yeoh et al 2017
out <- as.data.frame(t(sapply(allTaxaStrings, function(x) unlist(strsplit(x, "|", fixed = TRUE)))), stringsAsFactors = FALSE)
colnames(out) <- taxaToUse
rownames(out) <- allTaxaStrings

# all number of OTUs per core microbiome
for (curGroup in names(sigOTUsCoreMicrobiome)) {
  confCol <- paste0("conf_", curGroup)
  entCol <- paste0("entire_", curGroup)
  out[[confCol]] <- 0
  out[[entCol]] <- 0
  out[names(sigOTUsCoreMicrobiomeAnnotationConfident[[curGroup]]), confCol] <- sigOTUsCoreMicrobiomeAnnotationConfident[[curGroup]]
  out[names(sigOTUsCoreMicrobiomeAnnotationEntire[[curGroup]]), entCol] <- sigOTUsCoreMicrobiomeAnnotationEntire[[curGroup]]
}

# add average relative abundance
#temp <- aggregate(meanData[allOTUs,], by = list(taxonomy = apply(taxMod[allOTUs,taxaToUse], 1, function(x) paste(x, collapse = '|'))), sum)
#averageRelAbundanceConf <- t(round(t(temp[,colnames(meanData)])/colSums(meanData)*100, 2))
temp <- aggregate(t(t(meanData)/colSums(meanData))*100, by = list(taxonomy = apply(taxMod[rownames(meanData),taxaToUse], 1, function(x) paste(x, collapse = '|'))), mean)
averageRelAbundanceConf <- temp[,colnames(meanData)]
rownames(averageRelAbundanceConf) <- temp$taxonomy

#temp <- aggregate(meanData[allOTUs,], by = list(taxonomy = apply(taxModAll[allOTUs,taxaToUse], 1, function(x) paste(x, collapse = '|'))), sum)
#averageRelAbundanceEnt <- t(round(t(temp[,colnames(meanData)])/colSums(meanData)*100, 2))
temp <- aggregate(t(t(meanData)/colSums(meanData))*100, by = list(taxonomy = apply(taxModAll[rownames(meanData),taxaToUse], 1, function(x) paste(x, collapse = '|'))), mean)
averageRelAbundanceEnt <- temp[,colnames(meanData)]
rownames(averageRelAbundanceEnt) <- temp$taxonomy

out$averageAbundanceConf <- 0
out$averageAbundanceEntire <- 0
out[intersect(rownames(out), rownames(averageRelAbundanceConf)), "averageAbundanceConf"] <- round(rowMeans(averageRelAbundanceConf[intersect(rownames(out), rownames(averageRelAbundanceConf)),]),4)
out[intersect(rownames(out), rownames(averageRelAbundanceEnt)), "averageAbundanceEntire"] <- round(rowMeans(averageRelAbundanceEnt[intersect(rownames(out), rownames(averageRelAbundanceEnt)),]),4)

#for (curSpecies in colnames(meanData)) {
#  confCol <- paste0("averageAbundanceConf_", curSpecies)
#  entCol <- paste0("averageAbundanceEntire_", curSpecies)
#  out[[confCol]] <- 0
#  out[[entCol]] <- 0
#  out[rownames(averageRelAbundanceConf), confCol] <- averageRelAbundanceConf[,curSpecies]
#  out[rownames(averageRelAbundanceEnt), entCol] <- averageRelAbundanceEnt[,curSpecies]
#}

write.csv(out, file.path(rDir, "tableLikeYeoh.csv"), quote = FALSE, row.names = FALSE)

