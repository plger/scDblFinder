
#' doubletThresholding
#'
#' 
#'
#' @param scores A vector of the doublet score for each cell (real and artificial); can be
#' anything ranging from 0 to 1, with higher scores indicating higher change of being a 
#' doublet.
#' @param celltypes A vector of the same length as `scores` indicating, for each cell,
#' whether it is a 'real' cell or an 'artificial' doublet. Missing values not allowed.
#' @param clusters Optional vector of cluster assignment for each (real) cell, used for
#' homotypic doublet correction.
#' @param dbr The expected (mean) doublet rate, defaults to 0.021 on the basis of the 
#'  mixology 10x datasets demuxlet results, where the proportion of demuxlet doublet is
#'  estimated respectively at 0.012 and 0.029 (after accounting adding expected homotypic 
#'  doublets).
#' @param dbr.sd The deviation of the doublet rate, representing the uncertainty in the 
#' estimate. Effectively the scale of a Cauchy distribution. Defaults to 0.01.
#' @param prop.fullyRandom The proportion of artificical doublets that are fully random,
#' used for homotypy correction. Default 0.25 (the default value in 
#' `getArtificialDoublets`). Ignored if `clusters=NULL`
#' @param do.plot Logical; whether to plot the thresholding data (default TRUE).
#'
#' @return A scaler indicating the decided threshold.
#' @export
doubletThresholding <- function(scores, celltypes, clusters=NULL, dbr=0.021, dbr.sd=0.01, 
                                prop.fullyRandom=0.25, do.plot=TRUE){
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
  # if clusters given, adjust expected doublet proportion for expected homotypic doublets
  dbr <- dbr*(1-homotypic.prop)
  f1 <- ecdf(scores[which(celltypes!="artificial")])
  f2 <- ecdf(scores[which(celltypes=="artificial")])
  # accuracy functions
  accfn <- function(x){ (1-f1(x))+max(0,f2(x)-prop.fullyRandom*homotypic.prop) }
  accfn2 <- function(x){ (1-f1(x))+f2(x) }
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
    acc <- sapply(x, accfn)
    plot(x, acc, type="l", lty=ifelse(is.null(clusters),1,2), col="blue", ylim=c(0,max(acc)), lwd=2,
         main="Thresholding", ylab="Proportion", 
         xlab="Ratio of artificial doublets in neighborhood")
    if(!is.null(clusters)) lines(x, sapply(x, accfn2), col="blue", lwd=2)
    lines(x, sapply(x, dbr.dev), col="red", lwd=2)
    lines(x, sapply(x, FUN=function(x) sum(scores>=x & celltypes=="real")/sum(celltypes=="real")), col="darkgrey")
    abline(v=th, lty="dashed")
    leg <- c( "Error in classifying real cells vs artificial doublets",
              "Error in classifying artificial doublets, adjusted for homotypy",
              "Deviation from expected doublet rate",
              "Proportion of real cells considered doublet" )
    if(is.null(clusters)){
      leg <- leg[-2]
      legend("top", col=c("blue", "red","darkgrey"), lty=c(1,1,1), lw=2, legend=leg)
    }else{
      legend("top", col=c("blue","blue", "red","darkgrey"), lty=c(1,2,1,1), lw=2, legend=leg)
    }
    th
  }
}
