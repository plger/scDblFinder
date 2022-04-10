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
#' @param ret What to return, either barcode 'stats' (default), 'loci', or
#'   'coverages'.
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
#' @importFrom GenomicRanges reduce GRanges granges GRangesList
#' @importFrom BiocGenerics score
#' @importFrom S4Vectors mcols
#' @importFrom IRanges overlapsAny IRanges coverage slice
#' @importFrom Rsamtools TabixFile seqnamesTabix
#' @importFrom rtracklayer import
#' @importFrom stats ppois
#' @importFrom GenomeInfoDb keepSeqlevels seqlevels seqlevelsInUse
#' @importFrom GenomeInfoDb seqlengths seqlengths<- 
#' @export
getFragmentOverlaps <- function(x, barcodes=NULL, regionsToExclude=GRanges(
  c("M","chrM","MT","X","Y","chrX","chrY"),
  IRanges(1L,width=10^8)), 
  minFrags=500L, uniqueFrags=TRUE, 
  maxFragSize=1000L, removeHighOverlapSites=TRUE,
  fullInMemory=FALSE, BPPARAM=NULL, verbose=TRUE,
  ret=c("stats", "loci", "coverages")){
  # prepare inputs
  ret <- match.arg(ret)
  if(ret=="coverages" && !fullInMemory)
    stop("Returning coverages is currently only supported with fullInMemory=TRUE")
  if(is.null(BPPARAM)) BPPARAM <- BiocParallel::SerialParam()
  stopifnot(is.null(barcodes) || is.character(barcodes))
  if(!is.null(barcodes) && length(barcodes)==1 && file.exists(barcodes))
    barcodes <- readLines(barcodes) # assume barcodes to be a text file
  if(!is.null(regionsToExclude)){
    if(length(regionsToExclude)==1 && file.exists(regionsToExclude)){
      if(verbose) message(format(Sys.time(), "%X"),
                          " - Reading regions to exclude")
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
      if(verbose) message(format(Sys.time(), "%X"),
                          " - Reading Tabix-indexed fragment file and ",
                          "computing overlaps")
      x <- bplapply(seqnamesTabix(tf), BPPARAM=BPPARAM, FUN=function(x){
        if(verbose) cat(paste0(x,", "))
        getFragmentOverlaps(
          rtracklayer::import(tf, format="bed", 
                              which=GRanges(x, IRanges(1,10^8))), 
          barcodes, regionsToExclude=regionsToExclude, verbose=FALSE, 
          minFrags=0.00001, uniqueFrags=uniqueFrags, maxFragSize=maxFragSize,
          removeHighOverlapSites=removeHighOverlapSites, BPPARAM=NULL, ret=ret
        )
      })
      if(verbose){
        cat("\n")
        message(format(Sys.time(), "%X"), " - Merging")
      }
      if(ret=="loci") return(unlist(GRangesList(x)))
      if(ret=="coverages"){
        return(x[[1]])
        # x <- lapply(x, as.list)
        # names(x) <- NULL
        # x <- lapply(x, FUN=function(x){
        #   x[!sapply(x, FUN=function(x) length(x@values)==1L && all(x@values==0L))]
        # })
        # x <- x[lengths(x)>0L]
        # return(do.call(RleList, unlist(x, recursive=FALSE)))
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
        message("Fragment file is not tabix-indexed, requiring the",
                "whole file to be imported in memory.")
      }else if(verbose){
        message(format(Sys.time(), "%X"), " - Reading full fragments...")
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
  if(!all(!is.na(seqlengths(gr))))
    seqlengths(gr) <- setNames(sapply(split(end(gr), seqnames(gr)), max)
                               [seqlevels(gr)], seqlevels(gr))
  gr <- gr[(width(gr)<=maxFragSize),]
  gr$name <- as.factor(gr$name)
  if(!is.null(regionsToExclude)){
    regionsToExclude <- regionsToExclude[which(
                    as.factor(seqnames(regionsToExclude)) %in% seqlevels(gr))]
    if(length(regionsToExclude)>0){
      regionsToExclude <- keepSeqlevels(regionsToExclude,
                                        value=seqlevelsInUse(regionsToExclude),
                                        pruning.mode="coarse")
      gr <- gr[!overlapsAny(gr, regionsToExclude)]
    }
  }
  
  if(verbose) message(format(Sys.time(), "%X"),
                      " - Splitting and subsetting barcodes...")
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
  if(ret=="coverages"){
    if(verbose) message(format(Sys.time(), "%X"), " - Computing coverages")
    return(lapply(gr, FUN=coverage))
  }
  if(verbose) message(format(Sys.time(), "%X"), " - Obtaining overlaps...")
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
  if(ret=="loci") return(gr2)
  rm(grl,gr)
  gc(FALSE)
  tt <- table(gr2$name)
  d[names(tt),"total.nAbove2"] <- as.integer(tt)
  if(removeHighOverlapSites) gr2 <- .removeHighOverlapSites(gr2)
  tt <- table(gr2$name)
  d[names(tt),"nAbove2"] <- as.integer(tt)
  d
}

.removeHighOverlapSites <- function(gr, pthres=0.01, retExclusionRanges=FALSE){
  # remove loci that have >2 reads in too many cells
  ho <- reduce(gr, min.gapwidth=0L, with.revmap=TRUE)
  hol <- lengths(ho$revmap)
  ho$p <- ppois(hol, mean(hol), lower.tail=FALSE)
  if(length(indices2remove <- unlist(ho$revmap[which(ho$p<0.01)]))>0)
    gr <- gr[-indices2remove]
  if(retExclusionRanges) return(list(gr=gr, exclusion=granges(ho)))
  gr
}
