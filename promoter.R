## R script for getting promoter regions in whole genome
## usage: Rscript promoter.R 10000
## example: define the promoter region as the 10000bp 
## upstream and downstream of the TSS

# setwd('/Users/bioinfor/Desktop/')

args <- commandArgs()
promoter_region <- as.numeric(args[6])

library(TxDb.Hsapiens.UCSC.hg19.knownGene)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
cols(txdb)
keytypes(txdb)
Keys <- keys(txdb, keytype="GENEID")
cols <- c("TXCHROM", "TXSTART", "TXEND", "TXSTRAND")
Table <- select(txdb, keys=Keys, cols=cols, keytype="GENEID")

startpoint <- rep(NA, nrow(Table))
endpoint <- rep(NA, nrow(Table))


for (i in 1:nrow(Table)) {
  row <- Table[i,]
  Strand = row[3]
  if ( Strand == "+") {
    TSS = row[4]
    Start = TSS - promoter_region
    End = TSS + promoter_region
    startpoint[i] = Start
    endpoint[i] = End
  } else if ( Strand == "-") {
    TSS = row[5]
    Start = TSS + promoter_region
    End = TSS - promoter_region
    startpoint[i] = End
    endpoint[i] = Start
  }  
}

Promoter = cbind(Table[,2], startpoint, endpoint)
colnames(Promoter) <- c('chr','start','end')

write.table(Promoter, "promoter.bed", quote=FALSE, sep="\t",
            col.names=FALSE, row.names=FALSE)