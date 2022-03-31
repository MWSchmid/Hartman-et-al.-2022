#!/usr/bin/env Rscript

myPath <- "/media/mwschmid/myData/MWSchmid/MarcelVanDerHeijden_pesticides/GitIgnore_combinedNetwork/mergedForMine.csv"

## arguments from commandline
myarg <- commandArgs()
argPos <- grep("--args", myarg, fixed = TRUE)
myPath <- as.character(myarg[argPos+1])
outDir <- as.character(myarg[argPos+2])

f.print.message <- function(x) { cat("=== ", format(Sys.time(), "%Y %b %d %X"), paste0("=== ", x,"\n")) }

if (!dir.exists(outDir)) {
  f.print.message("Directory does not exist!")
  quit("no", 1)
}

myData <- read.csv(myPath, stringsAsFactors = FALSE, header = TRUE)

for (i in 1:1000) {
  set.seed(i)
  # sample 100 otus and shuffle per OTU all values
  toSample <- 100
  if (toSample > ncol(myData)) {
    toSample <- ncol(myData)
  }
  subData <- myData[, sample(colnames(myData), )]
  subData <- apply(subData, 2, function(x) sample(x, length(x)))
  write.csv(subData, file.path(outDir, paste0("randomData_", i, ".csv")), row.names = FALSE, quote = FALSE)
}


