#!/usr/bin/env Rscript

tempDir <- "/home/marc/temp"
rDir <- tempDir

# get arguments from commandline
myarg <- commandArgs()
argPos <- grep("--args", myarg, fixed = TRUE)
tempDir <- as.character(myarg[argPos+1])
rDir <- as.character(myarg[argPos+2])

# Note that the scripts assumes that there are otu and zotu tables present

# a function for printing
f.print.message <- function(x) { cat("=== ", format(Sys.time(), "%Y %b %d %X"), paste0("=== ", x,"\n")) }

for (tabType in c("otu", "zotu", "otu_ms4", "zotu_ms4")) {
  curEnding <- paste0("\\.", tabType, ".count")
  allFileNames <- list.files(tempDir, curEnding)
  if (length(allFileNames) == 0) {
    f.print.message(paste0("skipping", tabType))
    next
  }
  print(allFileNames)
  if (length(allFileNames) == 1) {
    out <- read.table(file.path(tempDir, allFileNames), sep = '\t', row.names = 1, header = TRUE, comment.char = "")
  } else {
    allTabs <- list()
    allRows <- c()
    allCols <- c()
    for (curFileName in allFileNames) {
      curPrefix <- gsub(curEnding, "", curFileName)
      temp <- read.table(file.path(tempDir, curFileName), sep = '\t', row.names = 1, header = TRUE, comment.char = "")
      allRows <- union(allRows, rownames(temp))
      allCols <- union(allCols, colnames(temp))
      allTabs[[paste0("fn", curPrefix)]] <- temp
    }
    out <- matrix(0, ncol = length(allCols), nrow = length(allRows), dimnames = list(allRows, allCols)); dim(out)
    for (cn in names(allTabs)) {
      temp <- allTabs[[cn]]
      out[rownames(temp), colnames(temp)] <- out[rownames(temp), colnames(temp)] + as.matrix(temp)
    }
  }
  write.csv(out, file.path(rDir, paste0("r", toupper(tabType), "s.csv")))
}
