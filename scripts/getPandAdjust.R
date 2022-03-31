#!/usr/bin/env Rscript
library("stats")
library("parallel")

myPath <- "/home/marc/temp/network_all_values.csv"
myRandom <- "/media/mwschmid/myData/MWSchmid/MarcelVanDerHeijden_pesticides/GitIgnore_results_ITS/OTU/allRandomResults.csv"

## arguments from commandline
myarg <- commandArgs()
argPos <- grep("--args", myarg, fixed = TRUE)
myPath <- as.character(myarg[argPos+1])
myRandom <- as.character(myarg[argPos+2])

f.print.message <- function(x) { cat("=== ", format(Sys.time(), "%Y %b %d %X"), paste0("=== ", x,"\n")) }

myData <- read.csv(myPath, header = TRUE, stringsAsFactors = FALSE)
myRandom <- read.csv(myRandom, header = TRUE, stringsAsFactors = FALSE)

# note: MICeNL and LR are two-sided
oneSided <- c("MICe","TICe","MASe","MEVe","MCNe")
twoSided <- c("MICeNL","LR")

getPoneSided <- function(x, background) {
  out <- sapply(x, function(y) mean(background > y, na.rm = TRUE))
  return(out)
}

getPtwoSided <- function(x, background) {
  pHigh <- sapply(x, function(y) mean(background > y, na.rm = TRUE))
  pLow <- sapply(x, function(y) mean(background < y, na.rm = TRUE))
  pMin <- apply(cbind(pHigh, pLow), 1, min)
  return(pMin)
}

numReps <- ceiling(nrow(myData)/15)
numParts <- rep(1:15, each = numReps)[1:nrow(myData)]
f.print.message("Calculating P-values.")
for (cn in oneSided) {
  f.print.message(cn)
  nc <- paste0("P_", cn)
  temp <- split(myData[[cn]], numParts)
  background <- myRandom[[cn]]
  results <- mclapply(temp, function(x) getPoneSided(x, background), mc.cores = 15)
  myData[[nc]] <- unlist(results)
}

for (cn in twoSided) {
  f.print.message(cn)
  nc <- paste0("P_", cn)
  temp <- split(myData[[cn]], numParts)
  background <- myRandom[[cn]]
  results <- mclapply(temp, function(x) getPtwoSided(x, background), mc.cores = 15)
  myData[[nc]] <- unlist(results)
}

f.print.message("Adjusting P-values.")
toAdjust <- grep("^P_", colnames(myData), value = TRUE)
for (tA in toAdjust) {
  tAnew <- gsub("^P", "FDR", tA)
  myData[[tAnew]] <- p.adjust(myData[[tA]], "BH")
}

f.print.message(paste0("writing ", paste0(myPath, ".adj")))
write.table(myData, paste0(myPath, ".adj"), sep = '\t', quote = FALSE)

f.print.message(paste0("writing ", paste0(myPath, ".adj.sig")))
sigData <- subset(myData, FDR_TICe < 0.05)
write.table(sigData, paste0(myPath, ".adj.sig"), sep = '\t', quote = FALSE)

f.print.message("finished.")



