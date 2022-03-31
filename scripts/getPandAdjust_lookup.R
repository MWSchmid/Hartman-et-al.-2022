#!/usr/bin/env Rscript
library("stats")
library("parallel")

myPath <- "/home/marc/temp/networkValTest.csv"
myRandom <- "/home/marc/temp/forTest.csv"

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

####################################
### Generate lookup tables first, then get P-values
precisionToUse <- 4
for (cn in c(oneSided, twoSided)) {
  temp <- myData[[cn]]
  temp <- round(temp, precisionToUse)
  minMax <- range(temp)
  toCheck <- round(seq(minMax[1], minMax[2], 1/(10^precisionToUse)), precisionToUse) # there was some imprecision in one case
  numReps <- ceiling(length(toCheck)/60)
  numParts <- rep(1:60, each = numReps)[1:length(toCheck)]
  forRunning <- split(toCheck, numParts)
  background <- myRandom[[cn]]
  if (cn %in% oneSided) {
    #results <- unlist(mclapply(forRunning, function(x) getPoneSided(x, background), mc.cores = 4))
    results <- unlist(lapply(forRunning, function(x) getPoneSided(x, background)))
  } else {
    #results <- unlist(mclapply(forRunning, function(x) getPtwoSided(x, background), mc.cores = 4))
    results <- unlist(lapply(forRunning, function(x) getPtwoSided(x, background)))
  }
  names(results) <- paste0("v_", format(toCheck, scientific = FALSE))
  # assign P-values
  f.print.message(cn)
  nc <- paste0("P_", cn)
  temp <- paste0("v_", format(temp, scientific = FALSE))
  myData[[nc]] <- results[temp]
  naEntries <- sum(is.na(results[temp]))
  if (naEntries > 0) {
    f.print.message(paste0("Found ", naEntries, " NAs, check:", cn))
    print(head(temp[is.na(results[temp])]))
  }
}

f.print.message(paste0("writing ", paste0(myPath, ".withP")))
write.table(myData, paste0(myPath, ".withP"), sep = '\t', quote = FALSE, row.names = FALSE)

####################################
### Generate FDRs

f.print.message("Adjusting P-values.")
toAdjust <- grep("^P_", colnames(myData), value = TRUE)
for (tA in toAdjust) {
  tAnew <- gsub("^P", "FDR", tA)
  myData[[tAnew]] <- p.adjust(myData[[tA]], "BH")
}

f.print.message(paste0("writing ", paste0(myPath, ".adj")))
write.table(myData, paste0(myPath, ".adj"), sep = '\t', quote = FALSE, row.names = FALSE)

f.print.message(paste0("writing ", paste0(myPath, ".adj.sig")))
sigData <- subset(myData, FDR_TICe < 0.05)
write.table(sigData, paste0(myPath, ".adj.sig"), sep = '\t', quote = FALSE, row.names = FALSE)

f.print.message("finished.")



