# annotation packages
library(VariantAnnotation)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)   

annotate.snps_with_genes <- function(input)
{
  # Takes a df with columns (id, chr, pos)
  target <- with(input,
                 GRanges( seqnames = Rle(chr),
                          ranges   = IRanges(pos, end=pos, names=rs),
                          strand   = Rle(strand("*")) ) )
  loc <- locateVariants(target, TxDb.Hsapiens.UCSC.hg19.knownGene, AllVariants())
  names(loc) <- NULL
  out <- as.data.frame(loc)
  out$names <- names(target)[ out$QUERYID ]
  out <- unique(out[, c("names", "seqnames", "GENEID")])
  Symbol2id <- as.list( org.Hs.egSYMBOL2EG )
  id2Symbol <- rep( names(Symbol2id), sapply(Symbol2id, length) )
  names(id2Symbol) <- unlist(Symbol2id)
  out$GENESYMBOL <- id2Symbol[ as.character(out$GENEID) ]
  return(out)
}