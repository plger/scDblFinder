#' doubletThresholding
#'
#' Sets the doublet scores threshold; typically called by 
#' \code{\link[scDblFinder]{scDblFinder}}.
#'
#' @param scores A vector of the doublet score for each cell (real and 
#' artificial); can be anything ranging from 0 to 1, with higher scores 
#' indicating higher change of being a doublet.
#' @param celltypes A vector of the same length as `scores` indicating, for each
#'  cell, whether it is a 'real' cell or an 'artificial' doublet. Missing values
#'   not allowed.
#' @param clusters Optional vector of cluster assignment for each (real) cell, 
#' used for homotypic doublet correction.
#' @param dbr The expected (mean) doublet rate.
#' @param dbr.sd The standard deviation of the doublet rate, representing the 
#' uncertainty in the estimate.
#' @param prop.fullyRandom The proportion of artificical doublets that are fully
#'  random, used for homotypy correction. Default 0.25 (the default value in 
#' `getArtificialDoublets`). Ignored if `clusters=NULL`
#' @param do.plot Logical; whether to plot the thresholding data (default TRUE).
#'
#' @return A scaler indicating the decided threshold.
#' 
#' @examples
#' # random data
#' m <- t(sapply( seq(from=0, to=5, length.out=50), 
#'                FUN=function(x) rpois(30,x) ) )
#' # generate doublets and merge them with real cells
#' doublets <- getArtificialDoublets(m, 30)
#' celltypes <- rep(c("real","artificial"), c(ncol(m), ncol(doublets)))
#' m <- cbind(m,doublets)
#' # dummy doublet scores:
#' scores <- abs(jitter(1:ncol(m),amount=10))
#' scores <- scores/max(scores)
#' # get threshold
#' doubletThresholding(scores, celltypes, do.plot=FALSE)
#' 
#' @importFrom stats pcauchy optimize ecdf
#' @importFrom graphics abline legend lines plot
#' @export
doubletThresholding <- function(scores, celltypes, clusters=NULL, dbr=0.025, 
                                dbr.sd=0.02, prop.fullyRandom=0, do.plot=TRUE){
  if(!all(sort(unique(celltypes))==c("artificial","real"))){
    stop("`celltypes` should be either 'real' or 'artificial'.")
  }
  if(any(is.na(scores))){
    w <- which(is.na(scores))
    scores <- scores[-w]
    celltypes <- celltypes[-w]
    warning(length(w)," cells with NA scores were discarded.")
  }
  if(!is.null(clusters) && length(unique(clusters))>1){
    homotypic.prop <- sum((table(clusters)/length(clusters))^2)
  }else{
    homotypic.prop <- 0
  }
  # if clusters given, adjust expected dbr for expected homotypic doublets
  dbr <- dbr*(1-homotypic.prop)
  f1 <- ecdf(scores[which(celltypes!="artificial")])
  f2 <- ecdf(scores[which(celltypes=="artificial")])
  # accuracy functions
  accfn <- function(x){ 
      min(1,(1-f1(x))+max(0,f2(x)-prop.fullyRandom*homotypic.prop))
  }
  accfn2 <- function(x){ min(1,(1-f1(x))+f2(x)) }
  # deviation from expected doublet proportion
  dbr.dev <- function(x){
    p <- sum(scores>=x & celltypes=="real")/sum(celltypes=="real")
    abs(pcauchy(abs(dbr-p), scale=dbr.sd)*2-1)^2
  }
  totfn <- function(x){
    accfn(x)+dbr.dev(x)
  }
  th <- optimize(totfn, c(0,1), maximum=FALSE)$minimum
  if(do.plot){
    # we plot the thresholding data
    x <- 1:99/100
    acc <- vapply(x, FUN.VALUE=double(1), FUN=accfn)
    dev <- vapply(x, FUN.VALUE=double(1), FUN=dbr.dev)
    ymax <- min(1,max(c(acc,dev)))
    plot(x, acc, type="l", lty=ifelse(is.null(clusters),1,2), col="blue", 
         ylim=c(0,ymax), lwd=2, main="Thresholding", 
         ylab="Proportion", xlab="Ratio of artificial doublets in neighborhood")
    if(!is.null(clusters))
      lines(x, vapply(x, FUN.VALUE=double(1), FUN=accfn2), col="blue", lwd=2)
    lines(x, dev, col="red", lwd=2)
    lines(x, vapply(x, FUN.VALUE=double(1), FUN=function(x){
            sum(scores>=x & celltypes=="real")/sum(celltypes=="real")
        }), col="darkgrey")
    abline(v=th, lty="dashed")
    leg <- c( "Error in classifying real cells vs artificial doublets",
              "Error in classifying artificial doublets, adjusted for homotypy",
              "Deviation from expected doublet rate",
              "Proportion of real cells considered doublet" )
    if(is.null(clusters)){
      leg <- leg[-2]
      legend("top", col=c("blue", "red","darkgrey"), 
             lty=c(1,1,1), lwd=2, legend=leg)
    }else{
      legend("top", col=c("blue","blue", "red","darkgrey"), 
             lty=c(1,2,1,1), lwd=2, legend=leg)
    }
  }
  th
}
