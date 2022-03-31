#!/usr/bin/env Rscript

rDir <- "/media/mwschmid/myData/MWSchmid/MarcelVanDerHeijden_rhizospheres/GitIgnore_results_16S/OTU"
myTreeFile <- "/media/mwschmid/myData/MWSchmid/MarcelVanDerHeijden_rhizospheres/GitIgnore_results_16S/OTU/aligned.tre"

rDir <- "/media/mwschmid/myData/MWSchmid/MarcelVanDerHeijden_rhizospheres/GitIgnore_results_ITS/OTU"
myTreeFile <- "/media/mwschmid/myData/MWSchmid/MarcelVanDerHeijden_rhizospheres/GitIgnore_results_ITS/OTU/aligned.tre"

##########################################################################################
### arguments from commandline
myarg <- commandArgs()
argPos <- grep("--args", myarg, fixed = TRUE)
rDir <- as.character(myarg[argPos+1])
myTreeFile <- as.character(myarg[argPos+2])

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
  source("/media/mwschmid/myData/MWSchmid/MarcelVanDerHeijden_rhizospheres/scripts/createAnalysisSets.R")
})

f.print.message <- function(...) { cat("=== ", format(Sys.time(), "%Y %b %d %X"), paste("===", ...,"\n")) }

pcname <- system('uname -n',intern=T)
numCoresPerMachine <- c("nuke" = 4, "styx" = 15, "marc-IEU" = 7, "piftel" = 30)
numCoresAvailable <- numCoresPerMachine[pcname]

##########################################################################################
### load data
load(file.path(rDir, "forDownstreamAnalyses.Rdata")) # contains: seqDat (filtered and eventually rarefied), normData, sampleTab (has color and Pch), taxMod, taxModAll, keepOnlyRootSamples, removeRhizobia, useRarefaction

if (length(unique(sampleTab$Type_sample)) > 1) {
  f.print.message("THIS SCRIPT ASSUMES THAT YOU USE ONLY THE ROOT SAMPLES!")
  quit("no", 1)
}

##########################################################################################
##### bacterial diversity and evenness
set.seed(333)
if (useRarefaction) {
  f.print.message("Data was rarefied, using log2(x+1) transformation.")
  dataForIndicesAndPermanovas <- log2(seqDat + 1)
} else {
  dataForIndicesAndPermanovas <- normData
}

specRich <- colSums(seqDat[, rownames(sampleTab)] > 0)
bioDiv <- vegan:::diversity(t(dataForIndicesAndPermanovas), index = "shannon", base = 2)
effRich <- 2^bioDiv
pielouEvenness <- vegan:::diversity(t(dataForIndicesAndPermanovas), index = "shannon")/log(specRich)
f.open.figure(rDir, "specRich_bioDiv_pielouEven.svg", TRUE, height = 12, width = 4)
par(mfrow = c(4,1))
f.histogram(specRich, main = "SR"); f.histogram(bioDiv, main = "SD"); f.histogram(effRich, main = "ER"); f.histogram(pielouEvenness, main = "PE")
f.close.figure()

##########################################################################################
### put together all test data
forTest <- sampleTab
forTest$Block <- factor(forTest$Block, ordered = FALSE)
forTest$Species <- factor(forTest$Species_name, ordered = FALSE)
forTest$SpeciesChar <- forTest$Species_name
forTest$AMF <- factor(forTest$AMF, ordered = FALSE)
forTest$Rhizobia <- factor(forTest$Rhizobia, ordered = FALSE)
forTest$MSR <- specRich[rownames(forTest)]
forTest$MBD <- bioDiv[rownames(forTest)]
forTest$MER <- effRich[rownames(forTest)]
forTest$MEV <- pielouEvenness[rownames(forTest)]

##########################################################################################
### get more variables for species contrasts
forTest$within_AMF_nonRhizo <- forTest$SpeciesChar
forTest$within_AMF_nonRhizo[!((sampleTab$AMF == "AMF") & (sampleTab$Rhizobia == "nonRhizo"))] <- "others"
forTest$within_AMF_Rhizo <- forTest$SpeciesChar
forTest$within_AMF_Rhizo[!((sampleTab$AMF == "AMF") & (sampleTab$Rhizobia == "Rhizo"))] <- "others"
forTest$within_nonAMF_nonRhizo <- forTest$SpeciesChar
forTest$within_nonAMF_nonRhizo[!((sampleTab$AMF == "nonAMF") & (sampleTab$Rhizobia == "nonRhizo"))] <- "others"
for (curFactor in grep("^within", colnames(forTest), value = TRUE)) {
  forTest[[curFactor]] <- factor(forTest[[curFactor]], ordered = FALSE)
}

for (curSpecies in unique(forTest$SpeciesChar)) {
  newFactorName <- paste0(curSpecies, "_vs_others")
  forTest[[newFactorName]] <- forTest$SpeciesChar
  forTest[[newFactorName]][forTest[[newFactorName]] != curSpecies] <- "others"
  forTest[[newFactorName]] <- factor(forTest[[newFactorName]], levels = c("others", curSpecies), ordered = FALSE)
}

##########################################################################################
### write out the data used in all tests
save(forTest, file = file.path(rDir, "allDataUsedInTests.Rdata"))
write.csv(forTest, file.path(rDir, "allDataUsedInTests.csv"))

##########################################################################################
### get analysis sets, there aren't so many
#source("/media/mwschmid/myData/MWSchmid/MarcelVanDerHeijden_rhizospheres/scripts/createAnalysisSets.R")
analysisSets <- f.get.simple.analysis.sets()
speciesContrasts <- f.get.species.contrasts(c("within_AMF_nonRhizo", "within_AMF_Rhizo", "within_nonAMF_nonRhizo")) # Species does not work now
detailedSpeciesContrastsNames <- grep("_vs_others$", colnames(forTest), value = TRUE)
detailedSpeciesContrastsList <- lapply(detailedSpeciesContrastsNames, function(x) f.get.species.contrasts(c(x, "Species"))) # Species does not work now
names(detailedSpeciesContrastsList) <- detailedSpeciesContrastsNames

# select only the full model for the detailedSpeciescontrast
detailedSpeciesContrastsFullModelFirstRhizo <- list()
detailedSpeciesContrastsFullModelFirstAMF <- list()
for (curName in names(detailedSpeciesContrastsList)) {
  detailedSpeciesContrastsFullModelFirstRhizo[[curName]] <- detailedSpeciesContrastsList[[curName]]$withBlockAndPotAgainstSpecies
  detailedSpeciesContrastsFullModelFirstAMF[[curName]] <- detailedSpeciesContrastsList[[curName]]$withBlockAndPotAgainstSpeciesFirstAMF
  if (curName == "Lupin_vs_others") { # the interaction kills this one
    #forTest[which(paste0(forTest$Rhizobia, forTest$AMF) == "RhizononAMF"),]
    detailedSpeciesContrastsFullModelFirstRhizo[[curName]]$formulaParts <- gsub("Rhizobia:AMF + ", "", detailedSpeciesContrastsFullModelFirstRhizo[[curName]]$formulaParts, fixed = TRUE)
    detailedSpeciesContrastsFullModelFirstRhizo[[curName]]$alternative_FandP <- detailedSpeciesContrastsFullModelFirstRhizo[[curName]]$alternative_FandP[c("Pot", "Rhizobia", "AMF")]
    detailedSpeciesContrastsFullModelFirstAMF[[curName]]$formulaParts <- gsub("AMF:Rhizobia + ", "", detailedSpeciesContrastsFullModelFirstAMF[[curName]]$formulaParts, fixed = TRUE)
    detailedSpeciesContrastsFullModelFirstAMF[[curName]]$alternative_FandP <- detailedSpeciesContrastsFullModelFirstAMF[[curName]]$alternative_FandP[c("Pot", "Rhizobia", "AMF")]
  }
}

##########################################################################################
### LMs (formulas and alternative F and P are also used for PERMANOVA)
modelResults <- f.run.and.save.model(analysisSets, forTest, rDir, "biodiv_models.xlsx", skipCSV = TRUE)
modelResults <- f.run.and.save.model(speciesContrasts, forTest, rDir, "biodiv_models_speciesContrasts.xlsx", skipCSV = TRUE)
modelResults <- f.run.and.save.model(detailedSpeciesContrastsFullModelFirstRhizo, forTest, rDir, "biodiv_models_detailedSpeciesContrastsRhizoFirst.xlsx", skipCSV = TRUE)
modelResults <- f.run.and.save.model(detailedSpeciesContrastsFullModelFirstAMF, forTest, rDir, "biodiv_models_detailedSpeciesContrastsAmfFirst.xlsx", skipCSV = TRUE)
#modelResults <- mclapply(detailedSpeciesContrastsNames, function(x) f.run.and.save.model(detailedSpeciesContrastsList[[x]], forTest, rDir, paste0("biodiv_models_detCont_", x, ".xlsx"), skipCSV = TRUE), mc.cores = numCoresAvailable)

subOtus <- t(dataForIndicesAndPermanovas[,rownames(sampleTab)]) # log2(x+1) counts
subExpDat <- forTest[rownames(subOtus),]
modelResults <- f.run.and.save.adonis(analysisSets, subOtus, subExpDat, rDir, "adonis_models.xlsx", numCores = numCoresAvailable, skipCSV = TRUE, verbose = FALSE)
modelResults <- f.run.and.save.adonis(speciesContrasts, subOtus, subExpDat, rDir, "adonis_models_speciesContrasts.xlsx", numCores = numCoresAvailable, skipCSV = TRUE, verbose = FALSE)
modelResults <- f.run.and.save.adonis(detailedSpeciesContrastsFullModelFirstRhizo, subOtus, subExpDat, rDir, "adonis_models_detailedSpeciesContrastsRhizoFirst.xlsx", numCores = numCoresAvailable, skipCSV = TRUE, verbose = FALSE)
modelResults <- f.run.and.save.adonis(detailedSpeciesContrastsFullModelFirstAMF, subOtus, subExpDat, rDir, "adonis_models_detailedSpeciesContrastsAmfFirst.xlsx", numCores = numCoresAvailable, skipCSV = TRUE, verbose = FALSE)
modelResults <- f.run.and.save.adonis.unifrac(analysisSets, subOtus, subExpDat, myTreeFile, rDir, "adonis_unifrac.xlsx", numCores = numCoresAvailable, skipCSV = TRUE, verbose = FALSE)
modelResults <- f.run.and.save.adonis.unifrac(speciesContrasts, subOtus, subExpDat, myTreeFile, rDir, "adonis_unifrac_speciesContrasts.xlsx", numCores = numCoresAvailable, skipCSV = TRUE, verbose = FALSE)
modelResults <- f.run.and.save.adonis.unifrac(detailedSpeciesContrastsFullModelFirstRhizo, subOtus, subExpDat, myTreeFile, rDir, "adonis_unifrac_detailedSpeciesContrastsRhizoFirst.xlsx", numCores = numCoresAvailable, skipCSV = TRUE, verbose = FALSE)
modelResults <- f.run.and.save.adonis.unifrac(detailedSpeciesContrastsFullModelFirstAMF, subOtus, subExpDat, myTreeFile, rDir, "adonis_unifrac_detailedSpeciesContrastsAmfFirst.xlsx", numCores = numCoresAvailable, skipCSV = TRUE, verbose = FALSE)

#modelResults <- mclapply(detailedSpeciesContrastsNames, function(x) f.run.and.save.adonis(detailedSpeciesContrastsList[[x]], subOtus, subExpDat, rDir, paste0("adonis_models_detCont_", x,".xlsx"), numCores = 1, skipCSV = TRUE, verbose = FALSE), mc.cores = numCoresAvailable)
#modelResults <- mclapply(detailedSpeciesContrastsNames, function(x) f.run.and.save.adonis.unifrac(detailedSpeciesContrastsList[[x]], subOtus, subExpDat, myTreeFile, rDir, paste0("adonis_unifrac_detCont_", x,".xlsx"), numCores = 1, skipCSV = TRUE, verbose = FALSE), mc.cores = numCoresAvailable)








