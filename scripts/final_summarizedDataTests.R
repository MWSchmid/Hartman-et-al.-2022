#!/usr/bin/env Rscript

rDir <- "/media/mwschmid/myData/MWSchmid/MarcelVanDerHeijden_rhizospheres/GitIgnore_results_16S_semiRarefied/OTU"
myTreeFile <- "/media/mwschmid/myData/MWSchmid/MarcelVanDerHeijden_rhizospheres/GitIgnore_results_16S_semiRarefied/OTU/aligned.tre"

rDir <- "/media/mwschmid/myData/MWSchmid/MarcelVanDerHeijden_rhizospheres/GitIgnore_results_16S_semiRarefied_noRhizobia/OTU"
myTreeFile <- "/media/mwschmid/myData/MWSchmid/MarcelVanDerHeijden_rhizospheres/GitIgnore_results_16S_semiRarefied_noRhizobia/OTU/aligned.tre"

rDir <- "/media/mwschmid/myData/MWSchmid/MarcelVanDerHeijden_rhizospheres/GitIgnore_results_ITS_semiRarefied_noGlomeromycota/OTU"
myTreeFile <- "/media/mwschmid/myData/MWSchmid/MarcelVanDerHeijden_rhizospheres/GitIgnore_results_ITS_semiRarefied_noGlomeromycota/OTU/aligned.tre"
rm(list=ls())

##########################################################################################
### arguments from commandline
myarg <- commandArgs()
argPos <- grep("--args", myarg, fixed = TRUE)
rDir <- as.character(myarg[argPos+1])
myTreeFile <- as.character(myarg[argPos+2])
useRarefaction <- TRUE
if ("--noRarefaction" %in% myarg) {
  useRarefaction <- FALSE
}

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
numCoresAvailable <- 1 # I think that it doesn't help, it's rather getting stuck

##########################################################################################
### species to family mapping
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

##########################################################################################
### load data
load(file.path(rDir, "forDownstreamAnalyses.Rdata")) # contains: seqDat (filtered), seqDatRarefied (filtered and rarefied), normData (based on seqDat), sampleTab (has color and Pch), taxMod, taxModAll, otusToRemove, removeRhizobia, removeGlomeromycota

if (length(unique(sampleTab$Type_sample)) > 1) {
  f.print.message("THIS SCRIPT ASSUMES THAT YOU USE ONLY THE ROOT SAMPLES!")
  quit("no", 1)
}

##########################################################################################
##### bacterial diversity and evenness
set.seed(333)
if (useRarefaction) {
  #f.print.message("Data was rarefied, using log2(x+1) transformation.")
  #dataForIndicesAndPermanovas <- log2(seqDatRarefied+1)
  f.print.message("Data was rarefied, not doing the transformation.")
  dataForIndicesAndPermanovas <- seqDatRarefied
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

#summary(effRichWithoutSymb)
#summary(effRichWithSymb)

##########################################################################################
### put together all test data
forTest <- sampleTab
forTest$Block <- factor(forTest$Block, ordered = FALSE)
forTest$Species <- factor(forTest$Species_name, ordered = FALSE)
forTest$SpeciesChar <- forTest$Species_name
forTest$AMF <- factor(forTest$AMF, ordered = FALSE)
forTest$Rhizobia <- factor(forTest$Rhizobia, ordered = FALSE)
forTest$isLupin <- ifelse(forTest$SpeciesChar == "Lupin", "yes", "no") # fit this before the symboGroup
forTest$SymboGroup <- factor(paste0(forTest$AMF, "_", forTest$Rhizobia), ordered = FALSE)
forTest$family <- speciesToFamily[forTest$SpeciesChar]
forTest$MSR <- specRich[rownames(forTest)]
forTest$MBD <- bioDiv[rownames(forTest)]
forTest$MER <- effRich[rownames(forTest)]
forTest$MEV <- pielouEvenness[rownames(forTest)]

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

##########################################################################################
# some testing for Christina
#anova(lm(MSR ~ AMF:Rhizobia + AMF + Rhizobia, data = forTest))
#anova(lm(MSR ~ AMF + Rhizobia + AMF:Rhizobia, data = forTest))
#car:::Anova(lm(MSR ~ AMF + Rhizobia + AMF:Rhizobia, data = forTest), type = 3)
#car:::Anova(lm(MSR ~ AMF + Rhizobia + AMF:Rhizobia, data = forTest), type = 2)
#car:::Anova(lm(MSR ~ AMF + Rhizobia + AMF:Rhizobia, data = forTest), type = 1)

##########################################################################################
### write out the data used in all tests
save(forTest, file = file.path(rDir, "allDataUsedInTests.Rdata"))
write.csv(forTest, file.path(rDir, "allDataUsedInTests.csv"))

##########################################################################################
### get analysis sets, there aren't so many
#source("/media/mwschmid/myData/MWSchmid/MarcelVanDerHeijden_rhizospheres/scripts/createAnalysisSets.R")
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

##########################################################################################
### LMs (formulas and alternative F and P are also used for PERMANOVA)
modelResults <- f.run.and.save.model(analysisSets, forTest, rDir, "biodiv_models.xlsx", skipCSV = TRUE)
modelResults <- f.run.and.save.model(speciesContrasts, forTest, rDir, "biodiv_models_speciesContrasts.xlsx", skipCSV = TRUE)
modelResults <- f.run.and.save.model(speciesContrastsMore, forTest, rDir, "biodiv_models_speciesContrastsMore.xlsx", skipCSV = TRUE)
modelResults <- f.run.and.save.model(detailedSpeciesContrasts, forTest, rDir, "biodiv_models_detailedSpeciesContrasts.xlsx", skipCSV = TRUE)

subOtus <- t(dataForIndicesAndPermanovas[,rownames(sampleTab)]) # now raw counts log2(x+1) counts
subExpDat <- forTest[rownames(subOtus),]
modelResults <- f.run.and.save.adonis.bray.curtis(analysisSets, subOtus, subExpDat, rDir, "adonis_models.xlsx", numCores = numCoresAvailable, skipCSV = TRUE, verbose = FALSE)
modelResults <- f.run.and.save.adonis.bray.curtis(speciesContrasts, subOtus, subExpDat, rDir, "adonis_models_speciesContrasts.xlsx", numCores = numCoresAvailable, skipCSV = TRUE, verbose = FALSE)
modelResults <- f.run.and.save.adonis.bray.curtis(speciesContrastsMore, subOtus, subExpDat, rDir, "adonis_models_speciesContrastsMore.xlsx", numCores = numCoresAvailable, skipCSV = TRUE, verbose = FALSE)
modelResults <- f.run.and.save.adonis.bray.curtis(detailedSpeciesContrasts, subOtus, subExpDat, rDir, "adonis_models_detailedSpeciesContrasts.xlsx", numCores = numCoresAvailable, skipCSV = TRUE, verbose = FALSE)
modelResults <- f.run.and.save.adonis.unifrac(analysisSets, subOtus, subExpDat, myTreeFile, rDir, "adonis_unifrac.xlsx", numCores = numCoresAvailable, skipCSV = TRUE, verbose = FALSE, dataAreCounts = TRUE)
modelResults <- f.run.and.save.adonis.unifrac(speciesContrasts, subOtus, subExpDat, myTreeFile, rDir, "adonis_unifrac_speciesContrasts.xlsx", numCores = numCoresAvailable, skipCSV = TRUE, verbose = FALSE, dataAreCounts = TRUE)
modelResults <- f.run.and.save.adonis.unifrac(speciesContrastsMore, subOtus, subExpDat, myTreeFile, rDir, "adonis_unifrac_speciesContrastsMore.xlsx", numCores = numCoresAvailable, skipCSV = TRUE, verbose = FALSE, dataAreCounts = TRUE)
modelResults <- f.run.and.save.adonis.unifrac(detailedSpeciesContrasts, subOtus, subExpDat, myTreeFile, rDir, "adonis_unifrac_detailedSpeciesContrasts.xlsx", numCores = numCoresAvailable, skipCSV = TRUE, verbose = FALSE, dataAreCounts = TRUE)

#modelResults <- mclapply(detailedSpeciesContrastsNames, function(x) f.run.and.save.adonis(detailedSpeciesContrastsList[[x]], subOtus, subExpDat, rDir, paste0("adonis_models_detCont_", x,".xlsx"), numCores = 1, skipCSV = TRUE, verbose = FALSE), mc.cores = numCoresAvailable)
#modelResults <- mclapply(detailedSpeciesContrastsNames, function(x) f.run.and.save.adonis.unifrac(detailedSpeciesContrastsList[[x]], subOtus, subExpDat, myTreeFile, rDir, paste0("adonis_unifrac_detCont_", x,".xlsx"), numCores = 1, skipCSV = TRUE, verbose = FALSE), mc.cores = numCoresAvailable)

##########################################################################################
### For the basic characterization
taxModOverall <- taxMod[rownames(dataForIndicesAndPermanovas),]
print(table(taxModOverall$domain))
if (grepl("16S", rDir)) {
  taxModOverall <- subset(taxModOverall, domain == "Bacteria")
} else {
  taxModOverall <- subset(taxModOverall, domain == "Fungi")
}
phylumCounts <- sort(table(taxModOverall$phylum), decreasing = TRUE)
perPhylum <- aggregate(normData[rownames(taxModOverall),], by = list(phylum = taxModOverall$phylum), sum)
rownames(perPhylum) <- perPhylum$phylum
perPhylum <- perPhylum[names(phylumCounts),setdiff(colnames(perPhylum), "phylum")]
perPhylumPerc <- t(round(t(perPhylum)/colSums(perPhylum)*100, 1))
write.csv(perPhylumPerc, file.path(rDir, "perPhylumPercentages.csv"))
topPhyla <- names(phylumCounts)[1:4]
print(round(sum(taxModOverall$phylum %in% topPhyla)/nrow(taxModOverall)*100, 1))
print(round(phylumCounts[topPhyla]/sum(phylumCounts)*100, 1))
if (grepl("16S", rDir)) {
  classCounts <- sort(table(taxModOverall$class), decreasing = TRUE)
  perClassPerc <- round(classCounts/sum(classCounts)*100, 1)
  print(perClassPerc[c("Alphaproteobacteria", "Betaproteobacteria", "Gammaproteobacteria", "Deltaproteobacteria")])
}





