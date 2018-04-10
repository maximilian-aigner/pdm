# annotation packages
library(VariantAnnotation)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)   

annotate_snps_with_genes <- function(input)
{
  # Takes a df with columns (id, chr, pos)
  target <- with(input,
                 GRanges( seqnames = Rle(chr),
                          ranges   = IRanges(pos, end=pos, names=rsid),
                          strand   = Rle(strand("*")) ) )
  loc <- locateVariants(target, TxDb.Hsapiens.UCSC.hg19.knownGene, AllVariants())
  names(loc) <- NULL
  out <- as.data.frame(loc)
  out$names <- names(target)[ out$QUERYID ]
  #out <- out[ , c("names", "seqnames", "start", "end", "LOCATION", "GENEID", "PRECEDEID", "FOLLOWID")]
  out <- out[, c("names", "seqnames", "GENEID")]
  #out <- out[!duplicated(out[,c("names", "seqnames")]),]
  out <- unique(out)
  Symbol2id <- as.list( org.Hs.egSYMBOL2EG )
  id2Symbol <- rep( names(Symbol2id), sapply(Symbol2id, length) )
  names(id2Symbol) <- unlist(Symbol2id)
  # out$GENEID <- as.factor(out$GENEID)
  # out$PRECEDEID <- as.factor(out$PRECEDEID)
  # out$FOLLOWID <- as.factor(out$FOLLOWID)
  # x <- unique( with(out, c(levels(GENEID), levels(PRECEDEID), levels(FOLLOWID))) )
  # table( x %in% names(id2Symbol) ) # good, all found
  
  out$GENESYMBOL <- id2Symbol[ as.character(out$GENEID) ]
  
  # out$PRECEDESYMBOL <- id2Symbol[ as.character(out$PRECEDEID) ]
  # out$FOLLOWSYMBOL <- id2Symbol[ as.character(out$FOLLOWID) ]
  return(out)
}