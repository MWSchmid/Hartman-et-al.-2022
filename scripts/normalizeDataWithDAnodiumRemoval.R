#!/usr/bin/env Rscript

baseDir <- "/media/mwschmid/myData/MWSchmid/MarcelVanDerHeijden_rhizospheres"
dataDir <- "/media/mwschmid/myData/MWSchmid/MarcelVanDerHeijden_rhizospheres/GitIgnore_usearchOutput_16S"
rDir <- "/media/mwschmid/myData/MWSchmid/MarcelVanDerHeijden_rhizospheres/GitIgnore_results_16S"
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
removeRhizobia <- FALSE
if ("--removeRhizobia" %in% myarg) {
  removeRhizobia <- TRUE
}
useRarefaction <- FALSE
if ("--useRarefaction" %in% myarg) {
  useRarefaction <- TRUE
}
if (!dir.exists(rDir)) {
  dir.create(rDir)
}

is16S <- grepl("16S", rDir)

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
##### the botXX numbers seem to be correct
#sampleTab16S$botNum <- paste0("bot", sampleTab16S$Sample_id)
#sampleTabITS$botNum <- paste0("bot", sampleTabITS$Sample_id)
#sampleTab16S$hurz <- "no"
#sampleTab16S$hurz[sampleTab16S$AMF == "AMF"] <- "yes"
#sum(moreInfo[paste0("bot", sampleTab16S$Sample_id), "AMF"] != sampleTab16S$hurz)
#sampleTab16S$hurz <- "no"
#sampleTab16S$hurz[sampleTab16S$Rhizobia == "Rhizo"] <- "yes"
#sum(moreInfo[paste0("bot", sampleTab16S$Sample_id), "Rhizobia"] != sampleTab16S$hurz)
#moreInfo[paste0("bot", sampleTab16S$Sample_id), ][moreInfo[paste0("bot", sampleTab16S$Sample_id), "Rhizobia"] != sampleTab16S$hurz,]
#subset(sampleTab16S, Species_name == "Lupin")
###### moreInfo <- read.csv(file.path(baseDir, "report", "design.csv"), header = TRUE, stringsAsFactors = FALSE, row.names = 1, sep = ',')
seqDat <- read.csv(file.path(dataDir, paste0("r", toupper(otuType), "s.csv")), stringsAsFactors = FALSE, row.names = 1); dim(seqDat)
colnames(seqDat) <- gsub("16S|ITS", "", colnames(seqDat))
colnames(seqDat) <- gsub("Sample67mix", "Sample64lo", colnames(seqDat)) # based on the barcode information and the design table I think that Sample67mix should be Sample64lo
colnames(seqDat) <- gsub("Sample39lo", "Sample49lo", colnames(seqDat)) # based on the  design table I think that Sample39lo should be Sample49lo
fullUtax <- read.table(file.path(dataDir, paste0(tolower(otuType), "s_tax.txt")), stringsAsFactors = FALSE, row.names = 1, sep = '\t', header = FALSE); dim(fullUtax)
taxMod <- f.modify.utax.taxon.table(fullUtax, onlyConfident = TRUE)
write.csv(taxMod, file.path(rDir, "taxModOnlyConfident.csv"))
taxModAll <- f.modify.utax.taxon.table(fullUtax, onlyConfident = FALSE)
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
minTotCou <- 30 # minimal total read count across all samples
minBinOcc <- 3 # occurs in at least 5 samples
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
# Use rarefaction yes/no
if (useRarefaction) {
  rarefyTo <- min(colSums(seqDat))
  seqDat <- t(rrarefy(t(seqDat), rarefyTo))
}

##########################################################################################
# Remove Rhizobia if requested (then also remove the nodule and leftover samples)
if (removeRhizobia) {
  if (!is16S) { f.print.message("Rhizobia can only be removed for 16S data."); quit("no", 1) }
  samplesToUseForRhizobiaFinder <- unique(sampleTab$Sample_id[sampleTab$Type_sample == "Nodule"])
  subSampleTab <- subset(sampleTab, Sample_id %in% samplesToUseForRhizobiaFinder)
  subDat <- seqDat[, rownames(subSampleTab)]
  totSum <- rowSums(subDat)
  binSum <- rowSums(subDat>0)
  combi <- (binSum>minBinOcc)&(totSum>minTotCou)
  subDat <- subDat[combi,]
  f.print.message("removed", nrow(seqDat)-nrow(subDat), "OTUs with too little coverage for Rhizobia identification")
  # differential abundance
  formulaString <- "~Pot+Species_name+Type_sample"
  formulaStringLess <- "~Pot+Species_name"
  dds <- DESeqDataSetFromMatrix(countData = subDat, colData = subSampleTab, design = formula(formulaString))
  ddsOverall <- DESeq(dds, minReplicatesForReplace=Inf, test = "LRT", reduced = as.formula(formulaStringLess))
  ddsDetailed <- DESeq(dds, minReplicatesForReplace=Inf)
  daResults <- list()
  daResults$mainEffect <- results(ddsOverall, cooksCutoff=FALSE)
  daResults$specific_nod_vs_leftovers <- results(ddsDetailed, contrast=c("Type_sample","Nodule","Left_over"), cooksCutoff=FALSE)
  daResults$specific_nod_vs_roots <- results(ddsDetailed, contrast=c("Type_sample","Nodule","Root"), cooksCutoff=FALSE)
  f.getUpregulatedOtus <- function(x) {
    out <- rownames(x)[(x$padj < 0.01) & (x$log2FoldChange > 1)]
    out <- out[!is.na(out)]
    return(out)
  }
  f.getSigOtus <- function(x) {
    out <- rownames(x)[x$padj < 0.01]
    out <- out[!is.na(out)]
    return(out)
  }
  sigOtusOverall <- f.getSigOtus(daResults$mainEffect)
  upRegOtus <- lapply(daResults[c("specific_nod_vs_leftovers", "specific_nod_vs_roots")], f.getUpregulatedOtus)
  print(lapply(upRegOtus, length))
  #candidatesForRemoval <- union(intersect(upRegOtus$specific_nod_vs_leftovers, sigOtusOverall), intersect(upRegOtus$specific_nod_vs_roots, sigOtusOverall))
  candidatesForRemoval <- union(upRegOtus$specific_nod_vs_leftovers, upRegOtus$specific_nod_vs_roots)
  f.print.message("Found", length(candidatesForRemoval), "candidates for removal")
  #doubleCheckDirection <- RNAseqWrapper:::f.summarize.columns(normData, byTab = data.frame(sample = rownames(subSampleTab), group = paste0(subSampleTab$Pot, "_", subSampleTab$Species_name, "_", subSampleTab$Type_sample), stringsAsFactors = FALSE), mean)
  #doubleCheckDirection <- RNAseqWrapper:::f.summarize.columns(normData, byTab = data.frame(sample = rownames(subSampleTab), group = subSampleTab$Type_sample, stringsAsFactors = FALSE), mean)
  #head(doubleCheckDirection[candidatesForRemoval,])
  # only remove rhizobia
  otusToRemove <- candidatesForRemoval[taxModAll[candidatesForRemoval,"genus"] == "Rhizobium"]
  #otusToRemove <- candidatesForRemoval
  f.print.message("Found", length(otusToRemove), "candidates for removal that were annotated (all annotations) as Rhizobium.")
  #sort(table(taxModAll[,"genus"])) # there are quite some more
  #sort(table(taxMod[,"order"]))
  #sort(table(taxModAll[rownames(normData),"order"]))
  # check how many reads are removed by this:
  numReadsRemoved <- colSums(seqDat[otusToRemove,])
  percentReadsRemoved <- round(numReadsRemoved/colSums(seqDat)*100, 2)
  f.print.message("The following percentage of reads will be removed:")
  print(percentReadsRemoved)
  # remove Rhizobia OTUs from the seqdata, remove non-root samples, refilter seq data, and if necessary re-rarefy (I kept the thing above to keep the "selectedOTUsForUniFrac.txt" identical)
  samplesToKeep <- rownames(sampleTab)[sampleTab$Type_sample == "Root"]
  OTUsToKeep <- setdiff(rownames(seqDat), otusToRemove)
  seqDat <- seqDat[OTUsToKeep, samplesToKeep]
  totSum <- rowSums(seqDat)
  binSum <- rowSums(seqDat>0)
  combi <- (binSum>minBinOcc)&(totSum>minTotCou)
  seqDat <- seqDat[combi,]
  f.print.message(paste0("Filter after Rhizobia-removal: Removed ", sum(!combi), " OTUs; ", sum(combi), " OTUs remain."))
  binSum <- rowSums(seqDat>0)
  forMINE <- binSum>minBinOccForMine
  if (useRarefaction) {
    rarefyTo <- min(colSums(seqDat))
    seqDat <- t(rrarefy(t(seqDat), rarefyTo))
  }
}

for (curSample in colnames(seqDat)) {
  percOfLib <- round(seqDat[,curSample]/sum(seqDat[,curSample])*100, 2)
  print(sum(sort(percOfLib, decreasing = TRUE)[1:5]))
}

##########################################################################################
# normalize with the design (it actually doesn't matter which formula one takes)
sampleTab$treat__ <- factor(paste0(sampleTab$Species_name, '_', sampleTab$Type_sample, '_', sampleTab$AMF, '_', sampleTab$Rhizobia, '_', sampleTab$Pot))
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
sampleTab$symboType <- paste0(sampleTab$AMF, '_', sampleTab$Rhizobia)
for (curType in unique(sampleTab$symboType)) {
  print(curType)
  print(unique(subset(sampleTab, symboType == curType)$Species_name))
}
speciesToColor <- c()

lupinColors <- brewer.pal(3, "Purples"); names(lupinColors) <- c("Lupin_Cluster", "Lupin_Nodule", "Lupin_Root") # nonAMF_Rhizo
AMF_nonRhizo <- brewer.pal(6, "Greens"); names(AMF_nonRhizo) <- c("Brome", "Petunia", "Tobacco", "Tomato", "Wheat", "Maize")
AMF_Rhizo <- brewer.pal(5, "Oranges"); names(AMF_Rhizo) <- c("Alfalfa", "Lotus", "Pea", "Medicago", "Trifolium")
nonAMF_nonRhizo <- brewer.pal(5, "Blues"); names(nonAMF_nonRhizo) <- c("Arabidopsis", "Brassica", "Cardamine", "Spinach", "Sugarbeet")
otherColors <- c(AMF_nonRhizo, AMF_Rhizo, nonAMF_nonRhizo)

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
sampleTab$Pch <- forPch[sampleTab$Type_sample]

##########################################################################################
# boxplots
tiff(file.path(rDir, "boxplots.tiff"), width = 800, height = 1600, compression = "lzw")
layout(matrix(1:2, ncol = 2))
boxplot(log2(seqDat+1), horizontal = TRUE, col = sampleTab[colnames(seqDat), "color"], las = 1, main = "Raw counts, log2(x+1)-transformed")
boxplot(normData, horizontal = TRUE, col = sampleTab[colnames(normData), "color"], las = 1, main = "Normalized counts, log2(x+1)-transformed")
dev.off()

##########################################################################################
# additional work that relies on other libraries
if (!runByMarc) { quit("no", 0) }

suppressPackageStartupMessages({
  library("RColorBrewer")
  library("RNAseqWrapper")
  library("limma")
  library("edgeR")
  library("DESeq2")
  library("vegan")
  source("/media/mwschmid/myData/MWSchmid/Development/R/RNAseqWrapperAds.R")
  source("/media/mwschmid/myData/MWSchmid/Development/R/RNAseqWrapper_overviewBandwidthFix.R")
  source("/media/mwschmid/myData/MWSchmid/Development/R/OTU_functions.R")
  source("/media/mwschmid/myData/MWSchmid/MethAn/MethAnR.R")
})

##########################################################################################
# scatterplots with few of them
set.seed(333)

selectedSamples <- sample(rownames(sampleTab), 10)
f.do.some.overview.fix(log2(seqDat[,selectedSamples]+1), rDir, "raw_random")
f.do.some.overview.fix(normData[,selectedSamples], rDir, "norm_random")

##########################################################################################
# sample-wise correlation matrices and boxplots with the normalized data 
f.generic.correlation.matrix(normData, rDir, "normDataPearson", corMethod = "pearson", ColSideColors = sampleTab[colnames(normData), "color"])
f.generic.correlation.matrix(normData, rDir, "normDataPearsonOHV", corMethod = "pearson", useOnlyHighVar = TRUE, ColSideColors = sampleTab[colnames(normData), "color"])

####################################
### tSNE on all OTUs
forTSNE <- normData[, rownames(sampleTab)]
rtsneResults <- f.run.select.plot.tSNE(t(as.matrix(forTSNE)), list(species = sampleTab[colnames(normData), "color"]), sampleTab[colnames(normData), "Pch"], rDir, "tSNE.svg")

#######################################################################################
### RDA
selectedOtus <- rownames(normData)   # highVarOtu
selectedSamples <- rownames(sampleTab)  #[!is.na(expDat$root_disease)]
subOtus <- t(normData[selectedOtus,selectedSamples])
subSampleTab <- sampleTab[selectedSamples,]

RDA2 <- rda(subOtus ~ Type_sample + Pot + Species_name + AMF + Rhizobia, subSampleTab)
axes <- summary(RDA2)$site
temp <- summary(RDA2)$cont; temp
curColor <- sampleTab[rownames(subOtus), "color"]
f.open.figure(rDir, "RDA.svg", TRUE, height = 10, width = 5)
par(mfrow=c(2,1))
plot(axes[,1:2], pch=rep(16, length(curColor)), col=curColor, bg=curColor, cex = 2) # mono hist is a dot, mix hist is a triangle
barplot(temp$importance["Eigenvalue",], main="PCA Eigenvalues", col="black", xlab = "", ylab = "")
f.close.figure()

print(RDA2)
print(temp$importance[,1:2])
