
#' getFragmentOverlaps
#' 
#' Count the number of overlapping fragments.
#' 
#' @param x The path to a fragments file, or a GRanges object containing the
#' fragments (with the `name` column containing the barcode, and optionally 
#' the `score` column containing the count).
#' @param barcodes Optional character vector of cell barcodes to consider
#' @param minFrags Minimum number of fragments for a barcode to be 
#' considered. If `uniqueFrags=TRUE`, this is the minimum number of unique 
#' fragments. Ignored if `barcodes` is given.
#' @param regionsToExclude A GRanges of regions to exclude. As per the original
#'   Amulet method, we recommend excluding repeats, as well as sex and 
#'   mitochondrial chromosomes.
#' @param uniqueFrags Logical; whether to use only unique fragments.
#' @param maxFragSize Integer indicating the maximum fragment size to consider
#' @param removeHighOverlapSites Logical; whether to remove sites that have
#'   more than two reads in unexpectedly many cells.
#' @param fullInMemory Logical; whether to process all chromosomes together.
#'   This will speed up the process but at the cost of a very high memory 
#'   consumption (as all fragments will be loaded in memory). This is anyway the
#'   default mode when `x` is not Tabix-indexed.
#' @param verbose Logical; whether to print progress messages.
#' @param BPPARAM A `BiocParallel` parameter object for multithreading. Note 
#'   that multithreading will increase the memory usage.
#'
#' @details When used on normal (or compressed) fragment files, this 
#' implementation is relatively fast (except for reading in the data) but it 
#' has a large memory footprint since the overlaps are performed in memory. It 
#' is therefore recommended to compress the fragment files using bgzip and index
#' them with Tabix; in this case each chromosome will be read and processed 
#' separately, leading to a considerably lower memory footprint.
#'
#' @return A data.frame with counts and overlap statistics for each barcode.
#' 
#' @importFrom GenomicRanges reduce GRanges granges
#' @importFrom BiocGenerics score
#' @importFrom S4Vectors mcols
#' @importFrom IRanges overlapsAny IRanges coverage slice
#' @importFrom Rsamtools TabixFile seqnamesTabix
#' @importFrom rtracklayer import
#' @importFrom stats ppois
#' @importFrom GenomeInfoDb keepSeqlevels
#' @export
getFragmentOverlaps <- function(x, barcodes=NULL, regionsToExclude=GRanges(
                                  c("M","chrM","MT","X","Y","chrX","chrY"),
                                  IRanges(1L,width=10^8)), 
                                minFrags=500L, uniqueFrags=TRUE, 
                                maxFragSize=1000L, removeHighOverlapSites=TRUE,
                                fullInMemory=FALSE, BPPARAM=NULL, verbose=TRUE){
  # prepare inputs
  if(is.null(BPPARAM)) BPPARAM <- BiocParallel::SerialParam()
  stopifnot(is.null(barcodes) || is.character(barcodes))
  if(!is.null(barcodes) && length(barcodes)==1 && file.exists(barcodes))
    barcodes <- readLines(barcodes) # assume barcodes to be a text file
  if(!is.null(regionsToExclude)){
    if(length(regionsToExclude)==1 && file.exists(regionsToExclude)){
      if(verbose) message("Reading regions to exclude...")
      regionsToExclude <- rtracklayer::import(regionsToExclude)
    }else{
      stopifnot(is.null(regionsToExclude) || is(regionsToExclude, "GRanges"))
    }
    regionsToExclude <- sort(regionsToExclude)
  }
  # prepare empty output for returns
  emptyOutput <- data.frame(row.names=character(0), nFrags=integer(0), 
                            uniqFrags=integer(0), nAbove2=integer(0))
  
  if(is.character(x) && length(x)==1){
    if(!file.exists(x)) stop("x should be a fragment file!")
    if(!fullInMemory &&
       is(tf <- tryCatch(TabixFile(x), error=function(e) NULL), "TabixFile")){
      if(verbose) message("Reading Tabix-indexed fragment file and ",
                          "computing overlaps...")
      x <- bplapply(seqnamesTabix(tf), BPPARAM=BPPARAM, FUN=function(x){
        if(verbose) cat(paste0(x,", "))
        getFragmentOverlaps(
          rtracklayer::import(tf, format="bed", 
                              which=GRanges(x, IRanges(1,10^8))), 
          barcodes, regionsToExclude=regionsToExclude, verbose=FALSE, 
          minFrags=0.00001, uniqueFrags=uniqueFrags, maxFragSize=maxFragSize,
          removeHighOverlapSites=removeHighOverlapSites, BPPARAM=NULL
        )
      })
      if(verbose){
        cat("\n")
        message("Merging...")
      }
      x <- x[unlist(lapply(x,nrow))>0]
      if(length(x)==0) return(emptyOutput)

      if(is.null(barcodes)){
        barcodes <- rowsum(
          unlist(lapply(x, FUN=function(x) x$nFrags)), 
                        unlist(lapply(x, row.names)))[,1]
        barcodes <- names(barcodes[barcodes>=minFrags])
      }
      return(as.data.frame(Reduce("+", lapply(x, FUN=function(x){
        x <- as.matrix(x[barcodes,])
        x[is.na(x)] <- 0
        x
      })), row.names=barcodes))
    }else{
      if(!fullInMemory){
        warning("Fragment file is not tabix-indexed, requiring the",
                "whole file to be imported in memory.")
      }else if(verbose){
        message("Reading full fragments...")
      }
      gr <- rtracklayer::import(x, format="bed")
    }
  }else{
    if(!is(x, "GRanges") || !("name" %in% colnames(mcols(x))))
      stop("`x` should either be a path to a fragments file, or a GRanges ",
           "with the 'name' column containing the cell barcode (and optionally
           the 'score' column containing the counts).")
    gr <- x
  }
  gr <- gr[(width(gr)<=maxFragSize),]
  gr$name <- as.factor(gr$name)
  if(!is.null(regionsToExclude)){
    regionsToExclude <- keepSeqlevels(regionsToExclude,
                          intersect(seqlevels(gr),seqlevels(regionsToExclude)),
                          pruning.mode="coarse")
    gr <- gr[!overlapsAny(gr, regionsToExclude)]
  }
  
  if(verbose) message("Splitting and subsetting barcodes...")
  uniqFrags <- table(gr$name)
  if(minFrags<1L & minFrags>0L) minFrags <- round(minFrags*length(gr))
  if(is.null(barcodes)){
    uniqFrags <- uniqFrags[uniqFrags>=minFrags]
  }else{
    if((mis <- length(setdiff(barcodes, names(uniqFrags))))>0)
      if(verbose)
        warning("Some barcodes (", mis, " or ",round(100*mis/length(gr),1),"%)",
                " are missing from the fragments file!")
    uniqFrags <- uniqFrags[intersect(names(uniqFrags), barcodes)]
  }
  gr <- gr[gr$name %in% names(uniqFrags)]
  gr$name <- droplevels(gr$name)
  if(length(gr)==0) return(emptyOutput)
  if(isFALSE(uniqueFrags) && !is.null(score(gr))){
    i <- rep(seq_along(gr), as.integer(score(gr)))
    gr <- GenomicRanges::split(granges(gr)[i], gr$name[i])
    rm(i)
  }else{
    gr <- GenomicRanges::split(granges(gr), gr$name)
  }
  if(verbose) message("Obtaining overlaps...")
  d <- data.frame(row.names=names(gr), nFrags=as.integer(lengths(gr)), 
                  uniqFrags=as.integer(uniqFrags[names(gr)]))
  d$nAbove2 <- 0L
  # obtain loci covered with >2 reads:
  grl <- GRangesList(lapply(gr, FUN=function(x){
    x <- slice(coverage(x), lower=3L, rangesOnly=TRUE)
    GRanges(rep(factor(names(x), seqlevels(gr)), lengths(x)),
            unlist(x, use.names=FALSE))
  }))
  gr2 <- unlist(grl, use.names=FALSE)
  gr2$name <- rep(factor(names(grl), row.names(d)),lengths(grl))
  rm(grl,gr)
  gc(FALSE)
  if(removeHighOverlapSites){
    # remove loci that have >2 reads in too many cells
    ho <- reduce(gr2, min.gapwidth=0L, with.revmap=TRUE)
    hol <- lengths(ho$revmap)
    ho$p <- ppois(hol, mean(hol), lower.tail=FALSE)
    if(length(indices2remove <- unlist(ho$revmap[which(ho$p<0.01)]))>0)
      gr2 <- gr2[-indices2remove]
  }
  tt <- table(gr2$name)
  d[names(tt),"nAbove2"] <- as.integer(tt)
  d
}


.overlapsFromCounts <- function(x, maxWidth=500L, ...){
  if(is(x, "SingleCellExperiment")){
    if(!is.null(ranges(x)) && 
       !all(sum(IRanges::width(ranges(x)))==0L) ){
      y <- counts(x)[which(IRanges::width(ranges(x))<=maxWidth),]
    }else{
      y <- counts(x)
    }
  }else if(is.matrix(x)){
    y <- x
  }else{
    stop("Object should be a count matrix or a SingleCellExperiment")
  }
  libsizes <- colSums(y)
  y <- y>2L
  rs <- rowSums(y)
  rp <- p.adjust(ppois(rs, lambda=mean(rs), lower.tail=FALSE))
  y <- y[rp>0.01,]
  data.frame(row.names=colnames(y), nFrags=libsizes, nAbove2=colSums(y))
}

#' amuletFromCounts
#' 
#' A reimplementation of the Amulet doublet detection method for single-cell 
#' ATACseq (Thibodeau, Eroglu, et al., Genome Biology 2021), based on tile/peak 
#' counts. Note that this is only a fast approximation to the original Amulet 
#' method, and performs considerably worse; for an equivalent implementation, 
#' see \code{\link{amulet}}.
#'
#' @param x A `SingleCellExperiment` object, or a matrix of counts with cells
#' as columns. If the rows represent peaks, it is recommended to limite their
#' width (see details).
#' @param maxWidth the maximum width for a feature to be included. This is 
#' ignored unless `x` is a `SingleCellExperiment` with `rowRanges`.
#'
#' @return If `x` is a `SingleCellExperiment`, returns the object with an 
#' additional `amuletFromCounts.q` colData column. Otherwise returns a vector of
#'  the amulet doublet q-values for each cell.
#' 
#' @details 
#' The rationale for the amulet method is that a single diploid cell should not 
#' have more than two reads covering a single genomic location, and the method 
#' looks for cells enriched with sites covered by more than two reads.
#' If the method is applied on a peak-level count matrix, however, larger peaks
#' can however contain multiple reads even though no single nucleotide is 
#' covered more than once. Therefore, in such case we recommend to limit the 
#' width of the peaks used for this analysis, ideally to maximum twice the upper
#' bound of the fragment size. For example, with a mean fragment size of 250bp 
#' and standard deviation of 125bp, peaks larger than 500bp are very likely to 
#' contain non-overlapping fragments, and should therefore be excluded using the
#' `maxWidth` argument.
#' 
#' @seealso \code{\link{amulet}}
#' 
#' @importFrom IRanges width
#' @importFrom SummarizedExperiment ranges
#' @export
#' @examples
#' x <- mockDoubletSCE()
#' x <- amuletFromCounts(x)
#' table(call=x$amuletFromCounts.q<0.05, truth=x$type)
amuletFromCounts <- function(x, maxWidth=500L){
  d <- .overlapsFromCounts(x, maxWidth)
  q <- p.adjust(ppois(d$nAbove2, lambda=mean(d$nAbove2), lower.tail=FALSE), 
                method="fdr")
  if(is(x,"SummarizedExperiment")){
    x$amuletFromCounts.q <- q
    return(x)
  }
  q
}



#' amulet
#' 
#- A reimplementation of the Amulet doublet detection method for single-cell 
#' ATACseq (Thibodeau, Eroglu, et al., Genome Biology 2021)
#'
#' @param x The path to a fragments file, or a GRanges object containing the
#' fragments (with the `name` column containing the barcode, and the `score`
#' column containing the count).
#' @param ... Any argument to \code{\link{getFragmentOverlaps}}.
#'
#' @details When used on normal (or compressed) fragment files, this 
#' implementation is relatively fast (except for reading in the data) but it 
#' has a large memory footprint since the overlaps are performed in memory. It 
#' is therefore recommended to compress the fragment files using bgzip and index
#' them with Tabix; in this case each chromosome will be read and processed 
#' separately, leading to a considerably lower memory footprint.
#'
#' @return A vector of the number of overlaps above `above`, named with cell 
#' barcodes
#' @importFrom GenomicRanges reduce
#' @importFrom S4Vectors splitAsList mcols
#' @importFrom IRanges overlapsAny
#' @export
amulet <- function(x, ...){
  d <- getFragmentOverlaps(x, ...)
  d$p.value <- ppois(d$nAbove2, mean(d$nAbove2), lower.tail=FALSE)
  d$q.value <- p.adjust(d$p.value, method="BH")
  d
}

