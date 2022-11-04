#' amulet
#' 
#- A reimplementation of the Amulet doublet detection method for single-cell 
#' ATACseq (Thibodeau, Eroglu, et al., Genome Biology 2021). The rationale is
#' that cells with unexpectedly many loci covered by more than two reads are 
#' more likely to be doublets.
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
#' separately, leading to a considerably lower memory footprint. See the 
#' underlying \code{\link{getFragmentOverlaps}} for details.
#'
#' @return A data.frame including, for each barcode, the number sites covered by
#'   more than two reads, the number of reads, and p- and q-values (low values
#'   indicative of doublets).
#' 
#' @examples
#' # here we use a dummy fragment file for example:
#' fragfile <- system.file( "extdata", "example_fragments.tsv.gz",
#'                          package="scDblFinder" )
#' res <- amulet(fragfile)
#' 
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


#' clamulet
#' 
#' Classification-powered Amulet-like method
#'
#' @param x The path to a fragment file (see \code{\link{getFragmentOverlaps}}
#'   for performance/memory-related guidelines)
#' @param artificialDoublets The number of artificial doublets to generate
#' @param iter The number of learning iterations (should be 1 to)
#' @param k The number(s) of nearest neighbors at which to gather statistics 
#' @param minCount The minimum number of cells in which a locus is detected to
#'   be considered. If lower than 1, it is interpreted as a fraction of the 
#'   number of cells.
#' @param threshold The score threshold used during iterations
#' @param maxN The maximum number of regions per cell to consider to establish
#'   windows for meta-features
#' @param nfeatures The number of meta-features to consider
#' @param max_depth The maximum tree depth
#' @param returnAll Logical; whether to return data also for artificial doublets
#' @param verbose Logical; whether to print progress information
#' @param ... Arguments passed to \code{\link{getFragmentOverlaps}}
#' 
#' @details 
#' `clamulet` operates similarly to the `scDblFinder` method, but generates 
#' doublets by operating on the fragment coverages. This has the advantage that
#' the number of loci covered by more than two reads can be computed for 
#' artificial doublets, enabling the use of this feature (along with the 
#' kNN-based ones) in a classification scheme. It however has the disadvantage
#' of being rather slow and memory hungry, and appears to be outperformed by a
#' simple p-value combination of the two methods (see vignette).
#'
#' @return A data.frame
#' @export
#' @importFrom IRanges Views viewMaxs slice
clamulet <- function(x, artificialDoublets=NULL, iter=2, k=NULL, minCount=0.001, 
                     maxN=500, nfeatures=25, max_depth=5, threshold=0.75, 
                     returnAll=FALSE, verbose=TRUE, ...){
  m <- .clamulet.buildTable(x, nADbls=artificialDoublets, maxN=maxN,
                            nfeatures=nfeatures, verbose=verbose, ...)
  d <- m$d
  labels <- m$labels
  
  if(verbose) message(format(Sys.time(), "%X"), " - Scoring network")
  k <- .defaultKnnKs(k, sum(labels=="real"))
  kd <- .evaluateKNN(m$m, labels, rep(NA_integer_,length(labels)), k=k)$d
  d <- m <- cbind(d, kd[,grep("weighted|ratio", colnames(kd))])
  m <- as.matrix(m)
  
  if(verbose) message(format(Sys.time(), "%X"), " - Iterative training")
  
  d$include.in.training <- TRUE
  d$type <- labels
  d$score <- predict(.xgbtrain(m, labels, max_depth=3), m)
  
  max.iter <-  iter
  while(iter>0){
    # remove cells with a high chance of being doublets from the training
    w <- which(d$type=="real" & d$score > 0.6)
    if(verbose) message("iter=",max.iter-iter,", ", length(w),
                        " cells excluded from training.")
    d$score <- tryCatch({
      fit <- .xgbtrain(m[-w,], d$type[-w], max_depth=max_depth, nthreads=1L)
      predict(fit, m)
    }, error=function(e) d$score)
    iter <- iter-1
  }
  d$include.in.training[w] <- FALSE
  if(!returnAll){
    d <- d[d$type=="real",]
    d$type <- NULL
  }
  if(verbose) message(format(Sys.time(), "%X"), " Done!")
  d
}


#' @importFrom GenomicRanges seqnames
.clamulet.buildTable <- function(x, nADbls=NULL, minCount=0.001, maxN=500, 
                                 nfeatures=25, verbose=TRUE, ...){
  co <- getFragmentOverlaps(x, ..., ret="coverages", fullInMemory=TRUE,
                            verbose=verbose)
  if(verbose) message(format(Sys.time(), "%X"), " - Obtaining windows")
  seqlvls <- unique(unlist(lapply(head(co,100), names)))
  if(minCount<1) minCount <- max(as.integer(round(minCount*length(co))),2L)
  windows <- coverage(unlist(GRangesList(lapply(co, FUN=function(x){
    x <- .sliceGR(x, seqlvls=seqlvls)
    if(length(x)>maxN) x <- x[sample.int(length(x),maxN)]
    x
  }))))
  windows <- reduce(.sliceGR(windows, lower=minCount, seqlvls=seqlvls),
                    min.gapwidth=100L)

  if(verbose) message(format(Sys.time(), "%X"), " - Obtaining window counts")
  windows <- split(ranges(windows), seqnames(windows))
  counts <- .getCovCounts(co, windows)

  if(verbose) message(format(Sys.time(), "%X"), " - Aggregating features")
  counts <- aggregateFeatures(counts, seq_len(min(ncol(counts)-1,10)),
                              k=min(nrow(counts), nfeatures))
  
  if(verbose) message(format(Sys.time(), "%X"),
                      " - Computing features for artificial doublets")
  # get combination of cells for artificial doublets
  if(is.null(nADbls)) nADbls <- length(co)
  comb <- matrix(sample.int(length(co), size=2*floor(nADbls*2.2), replace=TRUE),
                 ncol=2)
  comb <- comb[comb[,1]!=comb[,2],]
  comb <- head(comb[!duplicated(comb),], nADbls)
  
  # get counts for artificial doublets
  acounts <- counts[,comb[,1]]+counts[,comb[,2]]
  colnames(acounts) <- paste("artificial", seq_len(ncol(acounts)))
  
  if(verbose) message(format(Sys.time(), "%X"),
                      " - Counting overlaps for real cells")
  d <- data.frame(row.names=c(colnames(counts), colnames(acounts)), 
                  type=rep(factor(c("real","doublet"), c("real","doublet")),
                           c(ncol(counts),ncol(acounts))),
                  total=c(colSums(counts),colSums(acounts)))
  d$total.nAbove2 <- d$nAbove2 <- 0L
  
  # get overlaps for real cells
  gr2 <- .getOvsFromCovs(co, seqlvls = seqlvls)
  tt <- table(gr2$name)
  d[names(tt),"total.nAbove2"] <- as.integer(tt)
  gr2 <- .removeHighOverlapSites(gr2, retExclusionRanges=TRUE)
  exclusion <- gr2$exclusion
  gr2 <- gr2$gr
  tt <- table(gr2$name)
  d[names(tt),"nAbove2"] <- as.integer(tt)  
  
  if(verbose) message(format(Sys.time(), "%X"),
                      " - Counting overlaps for artificial doublets")
  # get overlaps for artificial cells
  co <- lapply(setNames(seq_len(nrow(comb)), colnames(acounts)),
               FUN=function(x){
    x <- comb[x,]
    .sumRleLists(co[[x[1]]], co[[x[2]]])
  })
  gr2 <- .getOvsFromCovs(co, seqlvls=seqlvls)
  tt <- table(gr2$name)
  d[names(tt),"total.nAbove2"] <- as.integer(tt)
  gr2 <- gr2[!overlapsAny(gr2, exclusion)]
  tt <- table(gr2$name)
  d[names(tt),"nAbove2"] <- as.integer(tt)
  
  m <- t(scale(t(t(10000*cbind(counts, acounts))/d$total)))
  list(m=m, d=d[,-1], labels=d$type)
}

# adds up to RleLists
#' @importFrom S4Vectors runValue
.sumRleLists <- function(x,y){
  names(useqlvls) <- useqlvls <- union(names(x),names(y))
  isEmptyRle <- function(x) length(runValue(x))==1
  as(lapply(useqlvls, FUN=function(seql){
    if(is.null(x[[seql]]) || isEmptyRle(x[[seql]])) return(y[[seql]])
    if(is.null(y[[seql]]) || isEmptyRle(y[[seql]])) return(x[[seql]])
    suppressWarnings(x[[seql]]+y[[seql]])
  }), "RleList")
}

# obtains loci with >lower coverage from a list of coverage RleList
.getOvsFromCovs <- function(co, seqlvls=NULL, lower=3L){
  if(is.null(seqlvls)) seqlvls <- unique(unlist(lapply(co,names)))
  grl <- GRangesList(lapply(co, lower=lower, seqlvls=seqlvls, FUN=.sliceGR))
  gr2 <- unlist(grl, use.names=FALSE)
  gr2$name <- rep(factor(names(grl), names(co)),lengths(grl))
  gr2
}

# obtains maximum count per window from a list of coverage RleList
#' @importFrom S4Vectors Rle
#' @importFrom GenomicRanges seqnames start end
.getCovCounts <- function(co, windows, fun=IRanges::viewMaxs){
  if(is(windows, "GRanges"))
    windows <- split(ranges(windows), seqnames(windows))
  windows <- windows[lengths(windows)>0]
  def <- as(lapply(windows, FUN=function(x) Rle(0L, max(end(x)))), "RleList")
  as(sapply(co, FUN=function(x){
    if(length( missing <- setdiff(names(def), names(x)) )>0)
      x <- c(x, def[missing])
    x <- x[names(windows)]
    unlist(fun(Views(x, windows)), use.names=FALSE)
  }), "CsparseMatrix")
}

# slide a RleList and returns a GR
.sliceGR <- function(x, lower=1L, seqlvls=NULL){
  if(is.null(seqlvls)) seqlvls <- names(x)
  x <- x[names(x) %in% seqlvls]
  x <- slice(x, lower=as.integer(round(lower)), rangesOnly=TRUE)
  if(sum(lengths(x))==0) return(GRanges(factor(levels=seqlvls)))
  GRanges(rep(factor(names(x), seqlvls), lengths(x)), 
          unlist(x, use.names=FALSE))
}



#' amuletFromCounts
#' 
#' A reimplementation of the Amulet doublet detection method for single-cell 
#' ATACseq (Thibodeau, Eroglu, et al., Genome Biology 2021), based on tile/peak 
#' counts. Note that this is only a fast approximation to the original Amulet 
#' method, and *performs considerably worse*; for an equivalent implementation, 
#' see \code{\link{amulet}}.
#'
#' @param x A `SingleCellExperiment` object, or a matrix of counts with cells
#' as columns. If the rows represent peaks, it is recommended to limite their
#' width (see details).
#' @param maxWidth the maximum width for a feature to be included. This is 
#' ignored unless `x` is a `SingleCellExperiment` with `rowRanges`.
#' @param exclude an optional `GRanges` of regions to be excluded. This is 
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
amuletFromCounts <- function(x, maxWidth=500L, exclude=c("chrM","M","Mt")){
  if(is(x, "SingleCellExperiment")){
    if(!is.null(ranges(x)) && 
       !all(sum(IRanges::width(ranges(x)))==0L) ){
      if(!is.null(exclude)){
        if(is.character(exclude)){
          x <- x[!(seqnames(granges(x)) %in% exclude),]
        }else{
          x <- x[!overlapsAny(granges(x),exclude)]
        }
      } 
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
  d <- data.frame(row.names=colnames(y), nFrags=libsizes, nAbove2=colSums(y))
  q <- p.adjust(ppois(d$nAbove2, lambda=mean(d$nAbove2), lower.tail=FALSE), 
                method="fdr")
  if(is(x,"SummarizedExperiment")){
    x$amuletFromCounts.q <- q
    return(x)
  }
  q
}
