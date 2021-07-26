#' clusterStickiness
#'
#' Tests for enrichment of doublets created from each cluster (i.e. cluster's
#' stickiness). Only applicable with >=4 clusters.
#' Note that when applied to an multisample object, this functions assumes that
#' the cluster labels match across samples.
#'
#'
#' @param x A table of double statistics, or a SingleCellExperiment on which
#' \link{scDblFinder} was run using the cluster-based approach.
#' @param type The type of test to use (quasibinomial recommended).
#' @param inclDiff Logical; whether to include the difficulty in the model. If
#' NULL, will be used only if there is a significant trend with the enrichment.
#' @param verbose Logical; whether to print additional running information.
#'
#'
#' @return A table of test results for each cluster.
#' @importFrom stats as.formula coef glm
#' @export
#'
#' @examples
#' sce <- mockDoubletSCE(rep(200,5))
#' sce <- scDblFinder(sce, clusters=TRUE, artificialDoublets=500)
#' clusterStickiness(sce)
clusterStickiness <- function(x, type=c("quasibinomial","nbinom","binomial","poisson"),
                              inclDiff=NULL, verbose=TRUE){
  type <- match.arg(type)
  if(is(x,"SingleCellExperiment")){
    x <- metadata(x)$scDblFinder.stats
    if(is.null(x)) stop("No doublet origin statistics; was scDblFinder run ",
                        "with the cluster-based approach?")
  }
  stopifnot(all(c("combination","observed","expected") %in% colnames(x)))

  if(is.null(inclDiff)) inclDiff <- length(unique(x$combination))>15

  ## build the model matrix of stickiness coefficients
  cls <- t(simplify2array(strsplit(x$combination,"+",fixed=TRUE)))
  if(length(unique(as.character(cls)))<4)
    stop("`clusterStickiness` can only be used with at least 4 clusters.")
  d <- as.data.frame(sapply(unique(as.character(cls)), FUN=function(cl){
    as.integer(apply(cls,1,FUN=function(j) any(j==cl)))
  }))
  celltypes <- colnames(d)
  colnames(d) <- paste0("stickiness.",colnames(d))
  x <- cbind(d,x)
  if(type %in% c("binomial","quasibinomial")){
    x$obs.p <- x$observed/sum(x$observed)
    logit <- function(x) log(x/(1 - x))
    x$exp.p <- logit(x$expected/sum(x$expected))
    x$difficulty <- scale(x$difficulty)
    f <- paste( "obs.p~0+offset(exp.p)+", paste(colnames(d),collapse="+"))
    if(inclDiff) f <- paste0(f,"+difficulty")
    mod <- glm(as.formula(f), data=x, family=type, weights=(x$observed+x$expected)/2)
  }else{
    if(type!="nbinom") x$expected <- log(x$expected)
    x$difficulty <- log(x$difficulty)
    f <- paste( "observed~0+offset(expected)+", paste(colnames(d),collapse="+") )
    if(inclDiff) f <- paste0(f,"+difficulty")
    if(type=="nbinom"){
      type <- .getThetaDist(x$observed, x$expected, verbose=verbose)
      x$expected <- log(x$expected)
    }
    mod <- glm(as.formula(f), data=x, family=type)
  }
  co <- coef(summary(mod))
  co <- as.data.frame(co[grep("stickiness",row.names(co)),])
  row.names(co) <- gsub("^stickiness\\.","",row.names(co))
  colnames(co)[4] <- c("p.value")
  co$FDR <- p.adjust(co$p.value)
  co[order(co$p.value, abs(co$Estimate)-co[,2]),]
}

#' doubletPairwiseEnrichment
#'
#' Calculates enrichment in any type of doublet (i.e. specific combination of
#' clusters) over random expectation.
#' Note that when applied to an multisample object, this functions assumes that
#' the cluster labels match across samples.
#'
#' @param x A table of double statistics, or a SingleCellExperiment on which
#' scDblFinder was run using the cluster-based approach.
#' @param lower.tail Logical; defaults to FALSE to test enrichment (instead of
#' depletion).
#' @param sampleWise Logical; whether to perform tests sample-wise in multi-sample
#' datasets. If FALSE (default), will aggregate counts before testing.
#' @param type Type of test to use.
#' @param inclDiff Logical; whether to regress out any effect of the
#' identification difficulty in calculating expected counts
#' @param verbose Logical; whether to output eventual warnings/notes
#'
#' @return A table of significances for each combination.
#'
#' @export
#' @importFrom stats chisq.test pnbinom pnorm ppois fitted
#' @examples
#' sce <- mockDoubletSCE()
#' sce <- scDblFinder(sce, clusters=TRUE, artificialDoublets=500)
#' doubletPairwiseEnrichment(sce)
doubletPairwiseEnrichment <- function(
  x, lower.tail=FALSE, sampleWise=FALSE,
  type=c("poisson","binomial","nbinom","chisq"),
  inclDiff=TRUE, verbose=TRUE){

  type <- match.arg(type)
  if(is(x,"SingleCellExperiment")){
    x <- metadata(x)$scDblFinder.stats
    if(is.null(x)) stop("No doublet origin statistics; was scDblFinder run ",
                        "with the cluster-based approach?")
  }
  stopifnot(all(c("combination","observed","expected") %in% colnames(x)))

  if("difficulty" %in% colnames(x) && inclDiff){
    theta <- .getThetaDist(x$observed, x$expected, verbose=verbose)
    mod <- glm(x$observed~0+offset(log(x$expected))+log(x$difficulty),
               family=theta)
    x$expected <- fitted(mod)
  }
  if(!sampleWise && "sample" %in% colnames(x))
    x <- aggregate(x[,setdiff(colnames(x),c("sample","combination"))],
                   by=x[,"combination",drop=FALSE], FUN=sum)
  if(type=="binomial"){
    p <- pbinom(x$observed,prob=x$expected/sum(x$expected),size=sum(x$observed),
                lower.tail=FALSE)
  }else{
    if(type=="nbinom"){
      theta <- .getThetaDist(x$observed, x$expected, retValue=TRUE,
                             verbose=verbose)
      if(is.infinite(theta)){
        type <- "poisson"
      }else{
        p <- pnbinom(x$observed, size=theta, mu=x$expected,
                     lower.tail=lower.tail)
      }
    }
    if(type=="poisson"){
      p <- ppois(x$observed, x$expected, lower.tail=lower.tail)
    }else if(type=="chisq"){
      x$other <- sum(x$observed)-x$observed
      x$p <- x$expected/sum(x$expected)
      p <- apply( x[,c("observed","other","p")],1,FUN=function(x){
        chisq.test(x[1:2], p=c(x[3],1-x[3]))$p.value
      })
    }
  }
  ler <- log2((1+x$observed)/(1+x$expected))
  if(lower.tail){
    p[which(ler>0)] <- 1
  }else{
    p[which(ler<0)] <- 1
  }
  d <- data.frame(combination=x$combination, log2enrich=ler, p.value=p,
                  FDR=p.adjust(p))
  if(!is.null(x$sample) && sampleWise) d <- cbind(sample=x$sample, d)
  d[order(d$p.value, 1/abs(d$log2enrich)),]
}

#' @importFrom MASS theta.ml
#' @importFrom stats poisson
#' @importFrom MASS negative.binomial
.getThetaDist <- function(y, mu, maxIter=100, verbose=TRUE, retValue=FALSE){
  theta <- try(MASS::theta.ml(y, mu, limit=maxIter), silent=TRUE)
  if(is(theta,"try-error")){
    if(verbose) warning("Not enough dispersion (theta diverges to infinity) ",
                        "- switching to poisson.")
    if(retValue) return(Inf)
    return(poisson())
  }
  if(verbose) message("theta=", theta)
  if(retValue) return(theta)
  return(negative.binomial(theta=theta))
}

