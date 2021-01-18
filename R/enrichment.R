#' clusterStickiness
#'
#' Tests for enrichment of doublets created from each cluster (i.e. cluster's
#' stickiness). Only applicable with >=4 clusters.
#'
#' @param x A table of double statistics, or a SingleCellExperiment on which
#' \link{scDblFinder} was run.
#' @param type The type of test to use (quasibinomial recommended).
#' @param inclDiff Logical; whether to include the difficulty in the model. If NULL,
#' will be used only if there is a significant trend with the enrichment.
#'
#' @return A table of test results for each cluster.
#' @export
#'
#' @examples
#' sce <- mockDoubletSCE(rep(200,5))
#' sce <- scDblFinder(sce, artificialDoublets=500)
#' clusterStickiness(sce)
clusterStickiness <- function(x, type=c("quasibinomial","nbinom1","binomial","poisson"),
                              inclDiff=FALSE){
  type <- match.arg(type)
  if(is(x,"SingleCellExperiment")){
    x <- metadata(sce)$scDblFinder.stats
  }

  if(!is.null(x$FNR)) x$difficulty <- x$FNR+max(x$FNR)*(x$difficulty/max(x$difficulty))
  if(is.null(inclDiff)){
    mod <- glm(x$observed~offset(log(x$expected))+log(x$difficulty), family="nbinom1")
    co.diff <- coef(summary(mod))["log(x$difficulty)",]
    inclDiff <- co.diff[[4]]<0.25
  }

  ## build the model matrix of stickiness coefficients
  cls <- t(simplify2array(strsplit(x$combination,"+",fixed=TRUE)))
  if(length(unique(as.character(cls)))<4)
    stop("`clusterStickiness` can only be used with at least 4 clusters.")
  d <- as.data.frame(sapply(unique(as.character(cls)), FUN=function(cl){
    as.integer(apply(cls,1,FUN=function(j) any(j==cl)))
  }))
  colnames(d) <- paste0("stickiness.",colnames(d))
  x <- cbind(d,x)
  if(type %in% c("binomial","quasibinomial")){
    x$obs.p <- x$observed/sum(x$observed)
    x$exp.p <- gtools::logit(x$expected/sum(x$expected))
    #x$difficulty <- gtools::logit(abs(x$difficulty)/max(x$difficulty))
    f <- paste( "obs.p~0+offset(exp.p)+", paste(colnames(d),collapse="+"))
    if(inclDiff) f <- paste0(f,"+difficulty")
    mod <- glm(as.formula(f), data=x, family=type, weights=(x$observed+x$expected)/2)
  }else{
    x$expected <- log(x$expected)
    x$difficulty <- log(x$difficulty)
    f <- paste( "observed~0+offset(expected)+", paste(colnames(d),collapse="+") )
    if(inclDiff) f <- paste0(f,"+difficulty")
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
#' @param x A table of double statistics, or a SingleCellExperiment on which scDblFinder
#' was run.
#' @param lower.tail Logical; defaults to FALSE to test enrichment (instead of depletion).
#' @param sampleWise Logical; whether to perform tests sample-wise in multi-sample
#' datasets. If FALSE (default), will aggregate counts before testing.
#' @param type Type of test to use.
#'
#' @return A table of significances for each combination.
#'
#' @examples
#' sce <- mockDoubletSCE()
#' sce <- scDblFinder(sce, artificialDoublets=500)
#' doubletPairwiseEnrichment(sce)
doubletPairwiseEnrichment <- function(x, lower.tail=FALSE, sampleWise=FALSE,
                                      type=c("poisson","binomial","chisq","nbinom1")){
  type <- match.arg(type)
  if(is(x,"SingleCellExperiment")) x <- metadata(sce)$scDblFinder.stats
  if(!sampleWise && "sample" %in% colnames(x))
    x <- aggregate(x[,setdiff(colnames(x),c("sample","combination"))],
                   by=x[,"combination",drop=FALSE], FUN=sum)
  if(type=="binomial"){
    p <- pbinom(x$observed,prob=x$expected/sum(x$expected),size=sum(x$observed),
                lower.tail=FALSE)
  }else if(type=="nbinom1"){
    fit <- suppressWarnings(MASS::fitdistr(x=x$observed, mu=x$expected,
                                           densfun="negative binomial"))
    p <- pnbinom(x$observed, size=fit$estimate, x$expected/sum(x$expected),
                 lower.tail=lower.tail)
  }else if(type=="poisson"){
    p <- ppois(x$observed, x$expected, lower.tail=lower.tail)
  }else if(type=="chisq"){
    x$other <- sum(x$observed)-x$observed
    x$p <- x$expected/sum(x$expected)
    p <- apply( x[,c("observed","other","p")],1,FUN=function(x){
      chisq.test(x[1:2], p=c(x[3],1-x[3]))$p.value
    })
  }
  ler <- log2((1+x$observed)/(1+x$expected))
  if(lower.tail){
    p[which(ler>0)] <- 1
  }else{
    p[which(ler<0)] <- 1
  }
  d <- data.frame(combination=x$combination, log2enrich=ler, p.value=p, FDR=p.adjust(p))
  if(!is.null(x$sample) && sampleWise) d <- cbind(sample=x$sample, d)
  d[order(d$p.value, 1/abs(d$log2enrich)),]
}
