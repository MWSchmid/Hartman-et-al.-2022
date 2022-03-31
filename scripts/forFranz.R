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

# load data
dataDir <- "/path/to/boGaEx"
dataDir <- "/home/marc/Downloads/boGaEx"

# 16S
sampleTab <- read.table(file.path(dataDir, "16S_design_botG.txt"), header = TRUE, stringsAsFactors = FALSE, row.names = 1, sep = '\t')
seqDat <- read.csv(file.path(dataDir, "16S_rOTUs.csv"), stringsAsFactors = FALSE, row.names = 1); dim(seqDat); colnames(seqDat) <- gsub("ITS|16S", "", colnames(seqDat))
colnames(seqDat) <- gsub("Sample67mix", "Sample64lo", colnames(seqDat)) # based on the barcode information and the design table I think that Sample67mix should be Sample64lo
colnames(seqDat) <- gsub("Sample39lo", "Sample49lo", colnames(seqDat)) # based on the  design table I think that Sample39lo should be Sample49lo
fullUtax <- read.table(file.path(dataDir, "16S_otus_tax.txt"), stringsAsFactors = FALSE, row.names = 1, sep = '\t', header = FALSE)
taxMod <- f.modify.utax.taxon.table(fullUtax, onlyConfident = TRUE)

setdiff(colnames(seqDat), rownames(sampleTab))
setdiff(rownames(sampleTab), colnames(seqDat))

# remove the samples that were skipped
sampleTab <- sampleTab[colnames(seqDat),]
write.csv(sampleTab, file.path(dataDir, "sampleTab_16S_onlySamplesFinallyUsed.csv"))

# ITS
sampleTab <- read.table(file.path(dataDir, "ITS_design_botG.txt"), header = TRUE, stringsAsFactors = FALSE, row.names = 1, sep = '\t')
seqDat <- read.csv(file.path(dataDir, "ITS_rOTUs.csv"), stringsAsFactors = FALSE, row.names = 1); dim(seqDat); colnames(seqDat) <- gsub("ITS|16S", "", colnames(seqDat))
colnames(seqDat) <- gsub("Sample67mix", "Sample64lo", colnames(seqDat)) # based on the barcode information and the design table I think that Sample67mix should be Sample64lo
colnames(seqDat) <- gsub("Sample39lo", "Sample49lo", colnames(seqDat)) # based on the  design table I think that Sample39lo should be Sample49lo
fullUtax <- read.table(file.path(dataDir, "ITS_otus_tax.txt"), stringsAsFactors = FALSE, row.names = 1, sep = '\t', header = FALSE)
taxMod <- f.modify.utax.taxon.table(fullUtax, onlyConfident = TRUE)
funGuildAnno <- read.table(file.path(dataDir, "ITS_otus_funGuild_matched_forR.txt"), header = TRUE, stringsAsFactors = FALSE, row.names = 1, sep = '\t')
# In case you want to use the taxonomy assignments from FunGuild:
# NOTE: this only works if you use onlyConfident = FALSE
funGuildTaxMod <- f.modify.utax.taxon.table(funGuildAnno, onlyConfident = FALSE)

setdiff(colnames(seqDat), rownames(sampleTab))
setdiff(rownames(sampleTab), colnames(seqDat))

# remove the samples that were skipped
sampleTab <- sampleTab[colnames(seqDat),]
write.csv(sampleTab, file.path(dataDir, "sampleTab_ITS_onlySamplesFinallyUsed.csv"))


