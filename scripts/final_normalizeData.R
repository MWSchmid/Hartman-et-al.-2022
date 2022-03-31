#!/usr/bin/env Rscript

baseDir <- "/media/mwschmid/myData/MWSchmid/MarcelVanDerHeijden_rhizospheres"
dataDir <- "/media/mwschmid/myData/MWSchmid/MarcelVanDerHeijden_rhizospheres/GitIgnore_usearchOutput_ITS"
rDir <- "/media/mwschmid/myData/MWSchmid/MarcelVanDerHeijden_rhizospheres/GitIgnore_results_ITS_semiRarefied_noGlomeromycota"
otuType <- "OTU"
samplesToRemoveStr <- "none"

# arguments from commandline
myarg <- commandArgs()
argPos <- grep("--args", myarg, fixed = TRUE)
baseDir <- as.character(myarg[argPos+1])
dataDir <- as.character(myarg[argPos+2])
rDir <- as.character(myarg[argPos+3])
otuType <- as.character(myarg[argPos+4])
samplesToRemoveStr <- as.character(myarg[argPos+5])
if (samplesToRemoveStr == "none") {
  samplesToRemove <- c()
} else {
  samplesToRemove <- unlist(strsplit(samplesToRemoveStr, ",", TRUE))
}
removeRhizobia <- FALSE # whether to remove primary symbionts, keeps only root samples
if ("--removeRhizobia" %in% myarg) {
  removeRhizobia <- TRUE
}
removeGlomeromycota <- FALSE # whether to remove primary symbionts, keeps only root samples
if ("--removeGlomeromycota" %in% myarg) {
  removeGlomeromycota <- TRUE
}

manualTesting <- FALSE
if ("--manualTesting" %in% myarg) {
  manualTesting <- TRUE
}

if (!dir.exists(rDir)) {
  dir.create(rDir)
}

is16S <- grepl("16S", rDir)
isITS <- grepl("ITS", rDir)

# some settings
useRobustGeoMean <- FALSE
pcname <- system('uname -n',intern=T)
runByMarc <- pcname %in% c("nuke", "styx", "marc-IEU", "piftel")

# functions
f.split.utax.taxo <- function(x) {
  x <- gsub("\\(.*?\\)", "", x)
  out <- unlist(strsplit(x, ",", TRUE))
  if (length(out) > 1) {
    out <- t(sapply(out, function(x) unlist(strsplit(x, ":"))))
  } else {
    out <- FALSE
  }
  return(out)
}

f.modify.utax.taxon.table <- function(taxDat, onlyConfident = TRUE) {
  for (i in 1:ncol(taxDat)) {
    if (length(grep("\\(.*?\\)", taxDat[1:100,i])) > 0) {
      colWithAll <- i
    } else {
      if (sum(taxDat[1:100,i] == "+") + sum(taxDat[1:100,i] == "-")) {
        colWithStrand <- i
      } else {
        colWithConfident <- i
      }
    }
  }
  colToUse <- ifelse(onlyConfident, colWithConfident, colWithAll)
  emptyEntry <- rep("ukn", nrow(taxDat))
  out <- data.frame(
    d = emptyEntry,
    p = emptyEntry,
    c = emptyEntry,
    o = emptyEntry,
    f = emptyEntry,
    g = emptyEntry,
    s = emptyEntry,
    row.names = rownames(taxDat),
    stringsAsFactors = FALSE
  )
  for (rn in rownames(taxDat)) {
    temp <- f.split.utax.taxo(taxDat[rn,colToUse])
    if (is.matrix(temp)) { out[rn,temp[,1]] <- temp[,2] }
  }
  colnames(out) <- c("domain", "phylum", "class", "order", "family", "genus", "species")
  # some cases have the genus but nothing in between - replace ukn with undef
  for (i in ncol(out):2) {
    lowerIsKnown <- out[,i] != "ukn"
    upperIsUnknown <- out[,i-1] == "ukn"
    out[lowerIsKnown&upperIsUnknown,i-1] <- paste0("undef_", out[lowerIsKnown&upperIsUnknown,i])
  }
  return(out)
}

f.geomeanForDESeq2 <- function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

f.print.message <- function(...) { cat("=== ", format(Sys.time(), "%Y %b %d %X"), paste("===", ...,"\n")) }

# libraries
suppressPackageStartupMessages({
  library("DESeq2")
  library("Rtsne")
  library("RColorBrewer")
  library("vegan")
})

# load data
rDir <- file.path(rDir, otuType)
dir.create(rDir, showWarnings = FALSE)

#sum(sampleTab$Sample_Idb != sampleTab$Sample_Idf)
sampleTab16S <- read.table(file.path(baseDir, "report", "16S_design_botG.txt"), header = TRUE, stringsAsFactors = FALSE, row.names = 1, sep = '\t')
sampleTabITS <- read.table(file.path(baseDir, "report", "ITS_design_botG.txt"), header = TRUE, stringsAsFactors = FALSE, row.names = 1, sep = '\t')
seqDat <- read.csv(file.path(dataDir, paste0("r", toupper(otuType), "s.csv")), stringsAsFactors = FALSE, row.names = 1); dim(seqDat)
colnames(seqDat) <- gsub("16S|ITS", "", colnames(seqDat))
colnames(seqDat) <- gsub("Sample67mix", "Sample64lo", colnames(seqDat)) # based on the barcode information and the design table I think that Sample67mix should be Sample64lo
colnames(seqDat) <- gsub("Sample39lo", "Sample49lo", colnames(seqDat)) # based on the  design table I think that Sample39lo should be Sample49lo
fullUtax <- read.table(file.path(dataDir, paste0(tolower(otuType), "s_tax.txt")), stringsAsFactors = FALSE, row.names = 1, sep = '\t', header = FALSE); dim(fullUtax)
taxMod <- f.modify.utax.taxon.table(fullUtax, onlyConfident = TRUE)
taxModAll <- f.modify.utax.taxon.table(fullUtax, onlyConfident = FALSE)
if (isITS) {
  f.print.message("removing", sum(taxMod$domain == "undef_Streptophyta"), "plant ITS sequences (confident assignment)")
  toKeep <- rownames(taxMod)[taxMod$domain != "undef_Streptophyta"]
  taxMod <- taxMod[toKeep,]
  taxModAll <- taxModAll[toKeep,]
  seqDat <- seqDat[toKeep,]
  f.print.message("keeping", length(toKeep), "OTUs")
  #quit("no", 0)
}
write.csv(taxMod, file.path(rDir, "taxModOnlyConfident.csv"))
write.csv(taxModAll, file.path(rDir, "taxModAll.csv"))

if (is16S) {
  sampleTab <- sampleTab16S
} else {
  sampleTab <- sampleTabITS
}

missingInTable <- setdiff(colnames(seqDat), rownames(sampleTab))
if (length(missingInTable) > 0) {
  f.print.message("missing samples:", paste0(missingInTable, collapse = '\n'))
  quit("no", 1)
}
sampleTab <- sampleTab[colnames(seqDat),]
sampleTab$seqDepth <- colSums(seqDat)
sampleTab$pch <- 18
sampleTab$pch[(sampleTab$AMF == "AMF") & (sampleTab$Rhizobia == "nonRhizo")] <- 16
sampleTab$pch[(sampleTab$AMF == "AMF") & (sampleTab$Rhizobia == "Rhizo")] <- 15
sampleTab$pch[(sampleTab$AMF == "nonAMF") & (sampleTab$Rhizobia == "nonRhizo")] <- 17

##########################################################################################
# check the read counts, nodules have the lowest counts by far
checkReadCounts <- colSums(seqDat)
print(sort(checkReadCounts, decreasing = TRUE))
sum(checkReadCounts<100)
range(checkReadCounts[checkReadCounts>50])

##########################################################################################
# remove requested samples - not necessary
if (length(samplesToRemove) > 0) {
  f.print.message("REMOVING SAMPLES:", paste0(samplesToRemove, collapse = ', '))
  toKeep <- setdiff(rownames(sampleTab), samplesToRemove)
  sampleTab <- sampleTab[toKeep,]
  seqDat <- seqDat[,rownames(sampleTab)]
}

##########################################################################################
# OTU preprocessing - filter low-abundant OTUs
if (is16S) {
  minTotCou <- 30 # minimal total read count across all samples
  minBinOcc <- 3 # occurs in at least 5 samples
} else {
  minTotCou <- 10 # minimal total read count across all samples
  minBinOcc <- 3 # occurs in at least 5 samples
}
minBinOccForMine <- ceiling(120*0.2)
totSum <- rowSums(seqDat)
binSum <- rowSums(seqDat>0)
combi <- (binSum>minBinOcc)&(totSum>minTotCou)
seqDat <- seqDat[combi,]
f.print.message(paste0("Filter: Removed ", sum(!combi), " OTUs; ", sum(combi), " OTUs remain."))
binSum <- rowSums(seqDat>0)
forMINE <- binSum>minBinOccForMine

print(sum(forMINE))

selectedOtus <- rownames(seqDat)   # highVarOtu
write.table(selectedOtus, file.path(rDir, "selectedOTUsForUniFrac.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)

##########################################################################################
# Define whout would be removed if it were removed
if (is16S) {
  subSampleTab <- subset(sampleTab, Type_sample == "Nodule")
  subDat <- seqDat[, rownames(subSampleTab)]
  totSum <- rowSums(subDat)
  binSum <- rowSums(subDat>0)
  combi <- (binSum>0)&(totSum>=10) # the Nodule samples are anyway the ones with few reads, rarefaction should not do so much
  subDat <- subDat[combi,]
  f.print.message("removed", nrow(seqDat)-nrow(subDat), "OTUs with too little coverage for Rhizobia identification")
  otusToRemove <- rownames(subDat)[taxModAll[rownames(subDat),"genus"] %in% unique(grep("hizobium", taxModAll$genus, value = TRUE))]
  write.table(otusToRemove, file.path(rDir, "symbionts.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
} else {
  otusToRemove <- rownames(seqDat)[rownames(seqDat) %in% rownames(taxModAll)[taxModAll$phylum == "Glomeromycota"]]
  write.table(otusToRemove, file.path(rDir, "symbionts.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
}

##########################################################################################
# If removal is requested, tell what and remove from the OTUs that should be kept
if ((is16S & removeRhizobia) | (isITS & removeGlomeromycota)) {
  f.print.message("removing", length(otusToRemove), "OTUs.")
  numReadsRemoved <- colSums(seqDat[otusToRemove,])
  percentReadsRemoved <- round(numReadsRemoved/colSums(seqDat)*100, 2)
  f.print.message("The following percentage of reads will be removed:")
  print(percentReadsRemoved)
  out <- data.frame(sampleName = names(percentReadsRemoved), percentRemoved = percentReadsRemoved)
  write.csv(out, file.path(rDir, "percentOTUreadsRemoved.csv"), row.names = FALSE, quote = FALSE)
  OTUsToKeep <- setdiff(rownames(seqDat), otusToRemove)
  write.table(otusToRemove, file.path(rDir, "primarySymbiontOtus.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.csv(seqDat[otusToRemove,], file.path(rDir, "primarySymbiontOtus.csv"), quote = FALSE)
} else {
  OTUsToKeep <- rownames(seqDat)
}

##########################################################################################
# Subset OTUs and samples
print(nrow(sampleTab))
sampleTab <- subset(sampleTab, Type_sample == "Root")
print(nrow(sampleTab))
seqDat <- seqDat[OTUsToKeep, rownames(sampleTab)]
totSum <- rowSums(seqDat)
binSum <- rowSums(seqDat>0)
combi <- (binSum>minBinOcc)&(totSum>minTotCou)
seqDat <- seqDat[combi,]
f.print.message(paste0("Filter after sample and OTU cleanup: Removed ", sum(!combi), " OTUs; ", sum(combi), " OTUs remain."))
binSum <- rowSums(seqDat>0)
forMINE <- binSum>minBinOccForMine

##########################################################################################
# Use rarefaction yes/no
rarefyTo <- min(colSums(seqDat))
seqDatRarefied <- t(rrarefy(t(seqDat), rarefyTo))
write.csv(seqDat, file.path(rDir, "seqDat.csv"))
write.csv(seqDatRarefied, file.path(rDir, "seqDatRarefied.csv"))

f.print.message("rarefied to", rarefyTo, "reads!")

##########################################################################################
# Percent of primary symbiont OTUs present
if (manualTesting){
  if (!(removeRhizobia & removeGlomeromycota)) {
    primSymb <- intersect(rownames(seqDat), otusToRemove)
    #numSeqPrimarySymb <- colSums(seqDat[primSymb, ] > 10)
    #sampleTab$numPrimarySymb <- NA
    #sampleTab[names(numSeqPrimarySymb), "numPrimarySymb"] <- numSeqPrimarySymb
    #sampleTab$percPrimSymb <- round(sampleTab$numPrimarySymb/length(primSymb)*100, 2)
    #numSeqPrimarySymb <- colSums(seqDatRarefied[primSymb, ] > 0)
    #percSeqPrimarySymb <- round(numSeqPrimarySymb/length(primSymb)*100, 2)
    numSeqPrimarySymb <- colSums(seqDatRarefied[primSymb, ])
    percSeqPrimarySymb <- round(numSeqPrimarySymb/colSums(seqDatRarefied)*100, 2)
    sampleTab$numPrimarySymb <- NA
    sampleTab[names(numSeqPrimarySymb), "numPrimarySymb"] <- numSeqPrimarySymb
    sampleTab$percPrimSymb <- NA
    sampleTab[names(percSeqPrimarySymb), "percPrimSymb"] <- percSeqPrimarySymb
    
    orderedSpecies <- rev(c("Lupin", "Alfalfa", "Lotus", "Pea", "Medicago", "Trifolium", "Brome", "Petunia", "Tobacco", "Tomato", "Wheat", "Maize", "Arabidopsis", "Brassica", "Cardamine", "Spinach", "Sugarbeet"))
    samplesOrdered <- c()
    for (curSpec in orderedSpecies) {
      samplesOrdered <- c(samplesOrdered, rownames(sampleTab)[sampleTab$Species_name == curSpec])
    }
    toPlot <- sampleTab[samplesOrdered,c("AMF", "Rhizobia", "Species_name", "numPrimarySymb", "percPrimSymb")]
    toPlot$Species_name <- factor(toPlot$Species_name, levels = orderedSpecies)
    svg(file.path(rDir, "percPrimSymbReads.svg"), height = 4, width = 6)
    boxplot(percPrimSymb ~ Species_name, data = toPlot, drawRect = TRUE, las = 1)
    jitterX <- jitter(rep(0, nrow(toPlot)), 10)
    points(as.numeric(toPlot$Species_name)+jitterX, toPlot$percPrimSymb, col = "black", pch = 16, cex = 1)
    invisible(dev.off())
    
    sampleTab$percPrimSymbReads <- sampleTab$percPrimSymb
    
    # same but present calls
    numSeqPrimarySymb <- colSums(seqDatRarefied[primSymb, ] > 0)
    percSeqPrimarySymb <- round(numSeqPrimarySymb/length(primSymb)*100, 2)
    sampleTab$numPrimarySymb <- NA
    sampleTab[names(numSeqPrimarySymb), "numPrimarySymb"] <- numSeqPrimarySymb
    sampleTab$percPrimSymb <- NA
    sampleTab[names(percSeqPrimarySymb), "percPrimSymb"] <- percSeqPrimarySymb
    
    orderedSpecies <- rev(c("Lupin", "Alfalfa", "Lotus", "Pea", "Medicago", "Trifolium", "Brome", "Petunia", "Tobacco", "Tomato", "Wheat", "Maize", "Arabidopsis", "Brassica", "Cardamine", "Spinach", "Sugarbeet"))
    samplesOrdered <- c()
    for (curSpec in orderedSpecies) {
      samplesOrdered <- c(samplesOrdered, rownames(sampleTab)[sampleTab$Species_name == curSpec])
    }
    toPlot <- sampleTab[samplesOrdered,c("AMF", "Rhizobia", "Species_name", "numPrimarySymb", "percPrimSymb")]
    toPlot$Species_name <- factor(toPlot$Species_name, levels = orderedSpecies)
    svg(file.path(rDir, "percPrimSymbPresent.svg"), height = 4, width = 6)
    boxplot(percPrimSymb ~ Species_name, data = toPlot, drawRect = TRUE, las = 1)
    jitterX <- jitter(rep(0, nrow(toPlot)), 10)
    points(as.numeric(toPlot$Species_name)+jitterX, toPlot$percPrimSymb, col = "black", pch = 16, cex = 1)
    invisible(dev.off())
    
    sampleTab$percPrimSymbPres <- sampleTab$percPrimSymb
    
    
    ##########################################################################################
    ### put together all test data
    speciesToFamily <- c(
      "Maize" = "Poaceae",
      "Wheat" = "Poaceae",
      "Brome" = "Poaceae",
      "Petunia" = "Solanaceae",
      "Tomato" = "Solanaceae",
      "Tobacco" = "Solanaceae",
      "Cardamine" = "Brassicaceae",
      "Arabidopsis" = "Brassicaceae",
      "Brassica" = "Brassicaceae",
      "Spinach" = "Amaranthaceae",
      "Sugarbeet" = "Amaranthaceae",
      "Trifolium" = "Fabaceae",
      "Pea" = "Fabaceae",
      "Medicago" = "Fabaceae",
      "Lotus" = "Fabaceae",
      "Alfalfa" = "Fabaceae",
      "Lupin" = "Fabaceae"
    )
        
    forTest <- sampleTab
    forTest$Block <- factor(forTest$Block, ordered = FALSE)
    forTest$Species <- factor(forTest$Species_name, ordered = FALSE)
    forTest$SpeciesChar <- forTest$Species_name
    forTest$AMF <- factor(forTest$AMF, ordered = FALSE)
    forTest$Rhizobia <- factor(forTest$Rhizobia, ordered = FALSE)
    forTest$isLupin <- ifelse(forTest$SpeciesChar == "Lupin", "yes", "no") # fit this before the symboGroup
    forTest$SymboGroup <- factor(paste0(forTest$AMF, "_", forTest$Rhizobia), ordered = FALSE)
    forTest$family <- speciesToFamily[forTest$SpeciesChar]

    ##########################################################################################
    ### get more variables for species contrasts
    forTest$correct_nn <- "no"
    forTest$correct_nn[forTest$SymboGroup == "nonAMF_nonRhizo"] <- "yes"
    forTest$correct_An <- "no"
    forTest$correct_An[forTest$SymboGroup == "AMF_nonRhizo"] <- "yes"
    forTest$correct_AR <- "no"
    forTest$correct_AR[forTest$SymboGroup == "AMF_Rhizo"] <- "yes"
    forTest$within_AMF_nonRhizo <- forTest$SpeciesChar
    forTest$within_AMF_nonRhizo[!((sampleTab$AMF == "AMF") & (sampleTab$Rhizobia == "nonRhizo"))] <- "others"
    forTest$within_AMF_Rhizo <- forTest$SpeciesChar
    forTest$within_AMF_Rhizo[!((sampleTab$AMF == "AMF") & (sampleTab$Rhizobia == "Rhizo"))] <- "others"
    forTest$within_nonAMF_nonRhizo <- forTest$SpeciesChar
    forTest$within_nonAMF_nonRhizo[!((sampleTab$AMF == "nonAMF") & (sampleTab$Rhizobia == "nonRhizo"))] <- "others"
    for (curFactor in grep("^within", colnames(forTest), value = TRUE)) {
      forTest[[curFactor]] <- factor(forTest[[curFactor]], ordered = FALSE)
    }
    
    forTest$between_nonAMF_nonRhizo_families <- forTest$family
    forTest$between_nonAMF_nonRhizo_families[!((sampleTab$AMF == "nonAMF") & (sampleTab$Rhizobia == "nonRhizo"))] <- "notOfInterest"
    forTest$between_AMF_nonRhizo_families <- forTest$family
    forTest$between_AMF_nonRhizo_families[!((sampleTab$AMF == "AMF") & (sampleTab$Rhizobia == "nonRhizo"))] <- "notOfInterest"
    
    for (curSpecies in unique(forTest$SpeciesChar)) {
      newFactorName <- paste0(curSpecies, "_vs_others")
      forTest[[newFactorName]] <- forTest$SpeciesChar
      forTest[[newFactorName]][forTest[[newFactorName]] != curSpecies] <- "others"
      forTest[[newFactorName]] <- factor(forTest[[newFactorName]], levels = c("others", curSpecies), ordered = FALSE)
    }

    # test percent primary symbionts
    source("/media/mwschmid/myData/MWSchmid/Development/R/lm_wrapper.R")
    source("/media/mwschmid/myData/MWSchmid/MarcelVanDerHeijden_rhizospheres/scripts/final_createAnalysisSets.R")
    analysisSets <- f.get.simple.analysis.sets()
    speciesContrasts <- f.get.species.contrasts(c("within_AMF_nonRhizo", "within_AMF_Rhizo", "within_nonAMF_nonRhizo")) # Species does not work now
    speciesContrastsMore <- f.get.species.contrasts(c("between_AMF_nonRhizo_families", "between_nonAMF_nonRhizo_families", "within_AMF_nonRhizo", "within_AMF_Rhizo", "within_nonAMF_nonRhizo")) # Species does not work now
    detailedSpeciesContrastsNames <- grep("_vs_others$", colnames(forTest), value = TRUE)
    detailedSpeciesContrastsList <- lapply(detailedSpeciesContrastsNames, function(x) f.get.species.contrasts(c(x, "Species"))) # Species does not work now
    names(detailedSpeciesContrastsList) <- detailedSpeciesContrastsNames
    detailedSpeciesContrastsList$Lupin_vs_others <- NULL
    detailedSpeciesContrasts <- list()
    for (curName in names(detailedSpeciesContrastsList)) {
      detailedSpeciesContrasts[[curName]] <- detailedSpeciesContrastsList[[curName]]$noBlockNoPotAgainstSpecies
    }
    for (toChange in names(analysisSets)) {
      analysisSets[[toChange]]$targetVars <- c("percPrimSymbReads", "percPrimSymbPres")
    }
    for (toChange in names(speciesContrasts)) {
      speciesContrasts[[toChange]]$targetVars <- c("percPrimSymbReads", "percPrimSymbPres")
    }
    for (toChange in names(speciesContrastsMore)) {
      speciesContrastsMore[[toChange]]$targetVars <- c("percPrimSymbReads", "percPrimSymbPres")
    }
    for (toChange in names(detailedSpeciesContrasts)) {
      detailedSpeciesContrasts[[toChange]]$targetVars <- c("percPrimSymbReads", "percPrimSymbPres")
    }
    
    modelResults <- f.run.and.save.model(analysisSets, forTest, rDir, "percSymb_models.xlsx", skipCSV = TRUE)
    modelResults <- f.run.and.save.model(speciesContrasts, forTest, rDir, "percSymb_models_speciesContrasts.xlsx", skipCSV = TRUE)
    modelResults <- f.run.and.save.model(speciesContrastsMore, forTest, rDir, "percSymb_models_speciesContrastsMore.xlsx", skipCSV = TRUE)
    modelResults <- f.run.and.save.model(detailedSpeciesContrasts, forTest, rDir, "percSymb_models_detailedSpeciesContrasts.xlsx", skipCSV = TRUE)
  }
  f.print.message("MANUAL TESTING WAS ON, QUITTING")
  quit("no", 0)
}



##########################################################################################
# normalize with the design (it actually doesn't matter which formula one takes)
if (length(unique(sampleTab$Type_sample)) == 1) {
  sampleTab$treat__ <- factor(paste0(sampleTab$Species_name, '_', sampleTab$AMF, '_', sampleTab$Rhizobia, '_', sampleTab$Pot))
} else {
  sampleTab$treat__ <- factor(paste0(sampleTab$Species_name, '_', sampleTab$Type_sample, '_', sampleTab$AMF, '_', sampleTab$Rhizobia, '_', sampleTab$Pot))
}
formulaString <- "~treat__"
design <- model.matrix(formula(formulaString), data = sampleTab, contrasts.arg = NULL)
dds <- DESeqDataSetFromMatrix(countData = round(seqDat[,rownames(sampleTab)]), colData = sampleTab, design = formula(formulaString))
if (useRobustGeoMean) {
  f.print.message("WARNING: using manually calculated geometric means accoring to phyloseq manual")
  geoMeans <-  apply(counts(dds), 1, f.geomeanForDESeq2)
  dds <- estimateSizeFactors(dds, geoMeans = geoMeans)
} else {
  dds <- DESeq(dds) # with outlier refitting
}
normData <- log2(DESeq2::counts(dds, normalized = TRUE) + 1)
write.csv(normData, file.path(rDir, "normData.csv"))
#testRLOG <- rlog(dds, blind=FALSE)

write.csv(t(normData), file.path(rDir, paste0("r", toupper(otuType), "s_forMine.csv")), quote = FALSE, row.names = FALSE)
write.csv(t(normData[forMINE,]), file.path(rDir, paste0("r", toupper(otuType), "s_forMine_20perc.csv")), quote = FALSE, row.names = FALSE)

##########################################################################################
# set colors
#sampleTab$symboType <- paste0(sampleTab$AMF, '_', sampleTab$Rhizobia)
#for (curType in unique(sampleTab$symboType)) {
#  print(curType)
#  print(unique(subset(sampleTab, symboType == curType)$Species_name))
#}
lupinColors <- brewer.pal(3, "Purples"); names(lupinColors) <- c("Lupin_Cluster", "Lupin_Nodule", "Lupin_Root") # nonAMF_Rhizo
AMF_nonRhizo <- brewer.pal(6, "Greens"); names(AMF_nonRhizo) <- c("Brome", "Petunia", "Tobacco", "Tomato", "Wheat", "Maize")
AMF_Rhizo <- brewer.pal(5, "Oranges"); names(AMF_Rhizo) <- c("Alfalfa", "Lotus", "Pea", "Medicago", "Trifolium")
nonAMF_nonRhizo <- brewer.pal(5, "Blues"); names(nonAMF_nonRhizo) <- c("Arabidopsis", "Brassica", "Cardamine", "Spinach", "Sugarbeet")
otherColors <- c(AMF_nonRhizo, AMF_Rhizo, nonAMF_nonRhizo)
speciesColors <- c(lupinColors["Lupin_Root"], otherColors)
names(speciesColors) <- gsub("_Root", "", names(speciesColors))

sampleTab$color <- NA
lupinMask <- sampleTab$Species_name == "Lupin"
sampleTab$color[lupinMask] <- lupinColors[paste0(sampleTab$Species_name, "_", sampleTab$Type_sample)[lupinMask]]
sampleTab$color[!lupinMask] <- otherColors[sampleTab$Species_name[!lupinMask]]
if (sum(is.na(sampleTab$color)) > 0) {
  f.print.message("NA colors:")
  print(subset(sampleTab, is.na(color)))
  f.print.message("ERROR - NA colors!")
  quit("no", 1)
}

forPch <- c("Cluster" = 17, "Left_over" = 15, "Nodule" = 18, "Root" = 16)
sampleTab$pchSampleType <- forPch[sampleTab$Type_sample]

##########################################################################################
# boxplots
tiff(file.path(rDir, "boxplots.tiff"), width = 800, height = 1600, compression = "lzw")
layout(matrix(1:2, ncol = 2))
boxplot(log2(seqDat+1), horizontal = TRUE, col = sampleTab[colnames(seqDat), "color"], las = 1, main = "Raw counts, log2(x+1)-transformed")
boxplot(normData, horizontal = TRUE, col = sampleTab[colnames(normData), "color"], las = 1, main = "Normalized counts, log2(x+1)-transformed")
dev.off()

##########################################################################################
# save an Rdata file with all necessary data for future analysis
save(seqDat, seqDatRarefied, normData, sampleTab, taxMod, taxModAll, otusToRemove, removeRhizobia, removeGlomeromycota, file = file.path(rDir, "forDownstreamAnalyses.Rdata"))

##########################################################################################
# write a table with the annotation counts
forSuppInfo <- taxMod[rownames(normData),]
for (toReport in setdiff(colnames(forSuppInfo), "domain")) {
  temp <- table(forSuppInfo[[toReport]])
  out <- data.frame(taxum = names(temp), occurence = as.vector(temp), stringsAsFactors = FALSE)
  out[out == "ukn"] <- "unknown"
  write.csv(out, file.path(rDir, paste0("taxaOccurences_", toReport, ".csv")), row.names = FALSE, quote = FALSE)
}

#quit("no", 0)
##########################################################################################
# additional work that relies on other libraries
if (!runByMarc) { quit("no", 0) }

suppressPackageStartupMessages({
  library("RColorBrewer")
  library("limma")
  library("edgeR")
  library("DESeq2")
  library("vegan")
  source("/media/mwschmid/myData/MWSchmid/Development/R/RNAseqWrapper.R")
  source("/media/mwschmid/myData/MWSchmid/Development/R/RNAseqWrapperAds.R")
  source("/media/mwschmid/myData/MWSchmid/Development/R/RNAseqWrapper_overviewBandwidthFix.R")
  source("/media/mwschmid/myData/MWSchmid/Development/R/OTU_functions.R")
  source("/media/mwschmid/myData/MWSchmid/MethAn/MethAnR.R")
})

f.open.figure(rDir, "speciesLegend.svg", TRUE)
f.plot.legend(speciesColors, 16)
f.close.figure()

##########################################################################################
### add group var
sampleTab$symboGroup <- paste0(sampleTab$AMF, '_', sampleTab$Rhizobia)

##########################################################################################
# scatterplots with few of them
set.seed(333)

selectedSamples <- sample(rownames(sampleTab), 10)
f.do.some.overview.fix(log2(seqDat[,selectedSamples]+1), rDir, "raw_random")
f.do.some.overview.fix(normData[,selectedSamples], rDir, "norm_random")
f.do.some.overview.fix(log2(seqDatRarefied[,selectedSamples]+1), rDir, "rarefied_random")

##########################################################################################
# sample-wise correlation matrices and boxplots with the normalized data 
f.generic.correlation.matrix(normData, rDir, "normDataPearson", corMethod = "pearson", ColSideColors = sampleTab[colnames(seqDatRarefied), "color"])
f.generic.correlation.matrix(normData, rDir, "normDataPearsonOHV", corMethod = "pearson", useOnlyHighVar = TRUE, ColSideColors = sampleTab[colnames(seqDatRarefied), "color"])
f.generic.correlation.matrix(log2(seqDatRarefied+1), rDir, "rarefiedDataPearson", corMethod = "pearson", ColSideColors = sampleTab[colnames(normData), "color"])
f.generic.correlation.matrix(log2(seqDatRarefied+1), rDir, "rarefiedPearsonOHV", corMethod = "pearson", useOnlyHighVar = TRUE, ColSideColors = sampleTab[colnames(normData), "color"])

####################################
### tSNE on all OTUs
forTSNE <- log2(seqDatRarefied+1)[, rownames(sampleTab)]
rtsneResults <- f.run.select.plot.tSNE(t(as.matrix(forTSNE)), list(species = sampleTab[colnames(seqDatRarefied), "color"]), sampleTab[colnames(seqDatRarefied), "Pch"], rDir, "tSNE.svg")
forTSNEaveragedSpecies <- f.summarize.columns(forTSNE, by = data.frame(sample = rownames(sampleTab), group = sampleTab$Species_name, stringsAsFactors = FALSE), mean)
rtsneResults <- f.run.select.plot.tSNE(t(as.matrix(forTSNEaveragedSpecies)), list(species = speciesColors[colnames(forTSNEaveragedSpecies)]), 16, rDir, "tSNE_means.svg", 101, perplexity = 5)

#######################################################################################
### RDA rarefied
selectedOtus <- rownames(seqDatRarefied)   # highVarOtu
selectedSamples <- rownames(sampleTab)  #[!is.na(expDat$root_disease)]
#subOtus <- t(log2(seqDatRarefied+1)[selectedOtus,selectedSamples]) # no log-transformation
subOtus <- t(seqDatRarefied[selectedOtus,selectedSamples]) # no log-transformation
subSampleTab <- sampleTab[selectedSamples,]

if (length(unique(subSampleTab$Type_sample)) == 1) {
  ###RDA2 <- rda(subOtus ~  Condition(Block) + Condition(Pot) + AMF + Rhizobia+ Species_name, subSampleTab)
  RDA2 <- rda(subOtus ~ Condition(seqDepth) + Condition(Pot) + AMF + Rhizobia+ Species_name, subSampleTab)
  ###RDA2 <- rda(subOtus ~ Condition(seqDepth) + Condition(Pot) + symboGroup + Species_name, subSampleTab)
  pcoa <- capscale(subOtus ~ Condition(seqDepth) + Condition(Pot) + AMF + Rhizobia+ Species_name, subSampleTab, distance = "bray")
  # Adding "by=term" to the anova() function you can analyze the model terms separately (Type I ss) or using "by=axis" gives the significance of each axis.
  print(permutest(pcoa, permutations = how(nperm = 999)) )
} else {
  ###RDA2 <- rda(subOtus ~ Type_sample + Block + Pot + AMF + Rhizobia + Species_name, subSampleTab)
  RDA2 <- rda(subOtus ~ Condition(seqDepth) + Condition(Type_sample) + Condition(Pot) + AMF + Rhizobia + Species_name, subSampleTab)
  ###RDA2 <- rda(subOtus ~ Condition(seqDepth) + Condition(Type_sample) + Condition(Pot) + symboGroup + Species_name, subSampleTab)
  pcoa <- capscale(subOtus ~ Condition(seqDepth) + Condition(Pot) + AMF + Rhizobia+ Species_name, subSampleTab, distance = "bray")
  # Adding "by=term" to the anova() function you can analyze the model terms separately (Type I ss) or using "by=axis" gives the significance of each axis.
  print(permutest(pcoa, permutations = how(nperm = 999)) )
}

curColor <- sampleTab[rownames(subOtus), "color"]
curPch <- sampleTab[rownames(subOtus), "pch"]
axes <- summary(RDA2)$site
temp <- summary(RDA2)$cont; temp
f.open.figure(rDir, "RDA_rarefied.svg", TRUE, height = 10, width = 5)
par(mfrow=c(2,1))
plot(axes[,1:2], col=curColor, bg=curColor, pch = curPch, cex = 2) # mono hist is a dot, mix hist is a triangle
barplot(temp$importance["Eigenvalue",], main="PCA Eigenvalues", col="black", xlab = "", ylab = "")
f.close.figure()

# get colors for species for kyle
forKyle <- unique(sampleTab[, c("Species_name", "color")])
#write.csv(forKyle, "/home/marc/Downloads/forKyleSpeciesColors.csv")

axes <- summary(pcoa)$sites
temp <- summary(pcoa)$cont; temp
f.open.figure(rDir, "CCA_rarefied.svg", TRUE, height = 10, width = 5)
par(mfrow=c(2,1))
plot(axes[,1:2], col=curColor, bg=curColor, pch = curPch, cex = 2) # mono hist is a dot, mix hist is a triangle
limitsInOriginal <- par("usr")
barplot(temp$importance["Eigenvalue",], main="PCA Eigenvalues", col="black", xlab = "", ylab = "")
f.close.figure()

f.open.figure(rDir, "CCA_rarefied_withEllipse.svg", TRUE, height = 10, width = 5)
par(mfrow=c(2,1))
plot(axes[,1:2], col=curColor, bg=curColor, pch = curPch, cex = 2) # mono hist is a dot, mix hist is a triangle
ordiellipse(pcoa, subSampleTab$symboGroup, xlim = limitsInOriginal[1:2], ylim = limitsInOriginal[3:4])
barplot(temp$importance["Eigenvalue",], main="PCA Eigenvalues", col="black", xlab = "", ylab = "")
f.close.figure()

print(pcoa)
print(temp$importance[,1:2])

# with the averages
if (length(unique(subSampleTab$Type_sample)) != 1) { 
  f.print.message("not doing average PCA with rarefied data")
} else {
  subOtus <- t(f.summarize.columns(t(subOtus), by = data.frame(sample = rownames(subSampleTab), group = subSampleTab$Species_name, stringsAsFactors = FALSE), mean))
  subSampleTab <- unique(sampleTab[,c("AMF", "Rhizobia", "symboGroup", "Species_name", "pch")])
  rownames(subSampleTab) <- subSampleTab$Species_name
  subSampleTab <- subSampleTab[rownames(subOtus),]
  subSampleTab$color <- speciesColors[rownames(subSampleTab)]
  RDA2 <- rda(subOtus ~ AMF + Rhizobia, subSampleTab)
  #RDA2 <- rda(subOtus ~ symboGroup, subSampleTab)
  axes <- summary(RDA2)$site
  temp <- summary(RDA2)$cont; temp
  curColor <- subSampleTab[rownames(subOtus), "color"]
  curPch <- subSampleTab[rownames(subOtus), "pch"]
  f.open.figure(rDir, "RDA_rarefied_averages.svg", TRUE, height = 10, width = 5)
  par(mfrow=c(2,1))
  plot(axes[,1:2], col=curColor, pch = curPch, bg=curColor, cex = 2) # mono hist is a dot, mix hist is a triangle
  barplot(temp$importance["Eigenvalue",], main="PCA Eigenvalues", col="black", xlab = "", ylab = "")
  f.close.figure()
  #print("Rarefaction average:")
  #print(RDA2)
  #print(temp$importance[,1:2])
  
  pcoa <- capscale(subOtus ~ AMF + Rhizobia, subSampleTab, distance = "bray")
  # Adding "by=term" to the anova() function you can analyze the model terms separately (Type I ss) or using "by=axis" gives the significance of each axis.
  #RDA2 <- rda(subOtus ~ symboGroup, subSampleTab)
  axes <- summary(pcoa)$sites
  temp <- summary(pcoa)$cont; temp
  curColor <- subSampleTab[rownames(subOtus), "color"]
  curPch <- subSampleTab[rownames(subOtus), "pch"]
  f.open.figure(rDir, "CCA_rarefied_averages.svg", TRUE, height = 10, width = 5)
  par(mfrow=c(2,1))
  plot(axes[,1:2], col=curColor, pch = curPch, bg=curColor, cex = 2) # mono hist is a dot, mix hist is a triangle
  limitsInOriginal <- par("usr")
  barplot(temp$importance["Eigenvalue",], main="PCA Eigenvalues", col="black", xlab = "", ylab = "")
  f.close.figure()
  print("Rarefaction average:")
  print(pcoa)
  print(temp$importance[,1:2])
  print(permutest(pcoa, permutations = how(nperm = 999)) )
  
  f.open.figure(rDir, "CCA_rarefied_averages_withEllipse.svg", TRUE, height = 10, width = 5)
  par(mfrow=c(2,1))
  plot(axes[,1:2], col=curColor, bg=curColor, pch = curPch, cex = 2) # mono hist is a dot, mix hist is a triangle
  ordiellipse(pcoa, subSampleTab$symboGroup, xlim = limitsInOriginal[1:2], ylim = limitsInOriginal[3:4])
  barplot(temp$importance["Eigenvalue",], main="PCA Eigenvalues", col="black", xlab = "", ylab = "")
  f.close.figure()
  
}

quit("no", 0)

#######################################################################################
### RDA DESeq
selectedOtus <- rownames(normData)   # highVarOtu
selectedSamples <- rownames(sampleTab)  #[!is.na(expDat$root_disease)]
subOtus <- t(normData[selectedOtus,selectedSamples])
subSampleTab <- sampleTab[selectedSamples,]

if (length(unique(subSampleTab$Type_sample)) == 1) {
  #RDA2 <- rda(subOtus ~ Condition(Block) + Condition(Pot) + AMF + Rhizobia+ Species_name, subSampleTab)
  RDA2 <- rda(subOtus ~ Condition(seqDepth) + Condition(Pot) + AMF + Rhizobia+ Species_name, subSampleTab)
  #RDA2 <- rda(subOtus ~ Condition(seqDepth) + Condition(Pot) + symboGroup + Species_name, subSampleTab)
} else {
  #RDA2 <- rda(subOtus ~ Type_sample + Block + Pot + AMF + Rhizobia + Species_name, subSampleTab)
  RDA2 <- rda(subOtus ~ Condition(seqDepth) + Condition(Type_sample) + Condition(Pot) + AMF + Rhizobia + Species_name, subSampleTab)
  #RDA2 <- rda(subOtus ~ Condition(seqDepth) + Condition(Type_sample) + Condition(Pot) + symboGroup + Species_name, subSampleTab)
}
axes <- summary(RDA2)$site
temp <- summary(RDA2)$cont; temp
curColor <- sampleTab[rownames(subOtus), "color"]
curPch <- sampleTab[rownames(subOtus), "pch"]
f.open.figure(rDir, "RDA_DESeq2.svg", TRUE, height = 10, width = 5)
par(mfrow=c(2,1))
plot(axes[,1:2], col=curColor, pch = curPch, bg=curColor, cex = 2) # mono hist is a dot, mix hist is a triangle
barplot(temp$importance["Eigenvalue",], main="PCA Eigenvalues", col="black", xlab = "", ylab = "")
f.close.figure()

print(RDA2)
print(temp$importance[,1:2])

if (length(unique(subSampleTab$Type_sample)) != 1) { 
  f.print.message("not doing average DESeq with rarefied data")
} else {
  subOtus <- t(f.summarize.columns(t(subOtus), by = data.frame(sample = rownames(subSampleTab), group = subSampleTab$Species_name, stringsAsFactors = FALSE), mean))
  subSampleTab <- unique(sampleTab[,c("AMF", "Rhizobia", "symboGroup", "Species_name", "pch")])
  rownames(subSampleTab) <- subSampleTab$Species_name
  subSampleTab <- subSampleTab[rownames(subOtus),]
  subSampleTab$color <- speciesColors[rownames(subSampleTab)]
  RDA2 <- rda(subOtus ~ AMF + Rhizobia, subSampleTab)
  #RDA2 <- rda(subOtus ~ symboGroup, subSampleTab)
  axes <- summary(RDA2)$site
  temp <- summary(RDA2)$cont; temp
  curColor <- subSampleTab[rownames(subOtus), "color"]
  curPch <- subSampleTab[rownames(subOtus), "pch"]
  f.open.figure(rDir, "RDA_DESeq2_averages.svg", TRUE, height = 10, width = 5)
  par(mfrow=c(2,1))
  plot(axes[,1:2], col=curColor, pch = curPch, bg=curColor, cex = 2) # mono hist is a dot, mix hist is a triangle
  barplot(temp$importance["Eigenvalue",], main="PCA Eigenvalues", col="black", xlab = "", ylab = "")
  f.close.figure()
  print("DESeq average:")
  print(RDA2)
  print(temp$importance[,1:2])
}









