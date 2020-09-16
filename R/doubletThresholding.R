#' doubletThresholding
#'
#' Sets the doublet scores threshold; typically called by 
#' \code{\link[scDblFinder]{scDblFinder}}.
#'
#' @param scores A vector of the doublet score for each cell (real and 
#' artificial); can be anything ranging from 0 to 1, with higher scores 
#' indicating higher change of being a doublet.
#' @param celltypes A vector of the same length as `scores` indicating, for each
#'  cell, whether it is a 'real' cell or a 'doublet'. Missing values
#'   not allowed.
#' @param clusters Optional vector of cluster assignment for each (real) cell, 
#' used for homotypic doublet correction.
#' @param dbr The expected (mean) doublet rate.
#' @param dbr.sd The standard deviation of the doublet rate, representing the 
#' uncertainty in the estimate.
#' @param prop.fullyRandom The proportion of artificical doublets that are fully
#'  random, used for homotypy correction. Default 0.25 (the default value in 
#' `getArtificialDoublets`). Ignored if `clusters=NULL`
#'
#' @return A scaler indicating the decided threshold.
#' 
#' @importFrom stats pcauchy optimize ecdf lm predict dnbinom
#' @export
doubletThresholding <- function( d, dbr=0.025, dbr.sd=0.02, local=TRUE ){
  # check that we have all necessary fields:
  if(!all(c("cluster","type","score","mostLikelyOrigin", "originAmbiguous") 
          %in% colnames(d))) stop("Input misses some columns.")
  
  if(!all(sort(as.character(unique(d$type)))==c("doublet","real"))){
    stop("`type` should be either 'real' or 'doublet'.")
  }
  expected <- getExpectedDoublets(d$cluster[!is.na(d$cluster)], dbr=dbr)
  if(local){
    ths <- .optimThresholds(d$score, d$type, d$mostLikelyOrigin, d$difficulty, 
                          expected)
    d$score.global <- d$score
    logit <- function(p) log(p/(1-p))
    adjp <- ths[d$mostLikelyOrigin,1]
    adjp[is.na(adjp)] <- median(ths)
    d$score <- 1/(1+exp(-(logit(d$score)-logit(sqrt(adjp)))))
  }
  w <- NULL
  if(!is.null(d$include.in.training)) w <- which(!d$include.in.training)
  ee <- sum(expected,na.rm=TRUE)
  ee <- c(ee*(1-dbr.sd), ee*(1+dbr.sd))
  th <- .optimThreshold(d$score, d$type, ee, fdr.include=w)
  o <- d$mostLikelyOrigin[d$type=="real" & d$score>=th]
  stats <- .compareToExpectedDoublets(o, expected, dbr)
  stats$combination <- as.character(stats$combination)
  stats$FNR <- vapply(split(as.data.frame(d), d$mostLikelyOrigin), 
                      FUN.VALUE=numeric(1L), 
                  FUN=function(x) .FNR(x$type, x$score, th) )[stats$combination]
  stats$difficulty <- vapply(split(d$difficulty, d$mostLikelyOrigin), 
                             FUN.VALUE=numeric(1L), na.rm=TRUE, 
                             FUN=median)[stats$combination]
  if(local){
    stats$local.threshold <- ths[stats$combination,2]
    stats$local.moderated <- ths[stats$combination,1]
  }
  list(th=th, stats=stats, finalScores=d$score)
}

.optimThreshold <- function(score, type, expected, fdr.include=NULL){
  type <- type=="real"
  if(!is.null(fdr.include)) fdr.include <- seq_along(score)
  totfn <- function(x){
    .FNR(type, score, x)*2 + .FDR(type[fdr.include], score[fdr.include], x) +
          .prop.dev(type,score,expected,x)^2
  }
  optimize(totfn, c(0,1), maximum=FALSE)$minimum
}

.optimThresholds <- function(score, type, origins, difficulty, expected, 
                             moderate=TRUE){
  ll <- split(seq_along(score), origins)
  names(n) <- n <- names(ll)
  ths <- vapply(n, FUN.VALUE=numeric(1L), FUN=function(comb){
    i <- ll[[comb]]
    type <- type=="real"
    totfn <- function(x){
      .FNR(type[i], score[i], x)*2 + .FDR(type[i], score[i], x) + 
        .prop.dev(type[i],score[i],expected[comb],x)^2
    }
    optimize(totfn, c(0,1), maximum=FALSE)$minimum
  })
  if(moderate){
    diff <- sapply(split(difficulty,origins), na.rm=TRUE, FUN=median)
    diff <- diff[names(ths)]
    mod <- tryCatch({ MASS::rlm(ths~diff) }, error=function(e){
        lm(ths~diff)
    })
    ths2 <- (ths+predict(mod))/2
    if(length(w <- which(ths2<0.01))>0)
      ths2[w] <- apply(cbind(ths[w],rep(0.01,length(w))),1,FUN=max)
    ths2[ths2>1] <- 1
    ths <- cbind(moderated=ths2, local=ths)
  }
  cbind(ths)
}

.FNR <- function(type, score, threshold){
  if(!is.logical(type)) type <- type=="real"
  sum(!type & score<threshold, na.rm=TRUE)/sum(!type)
}
.FDR <- function(type, score, threshold, include=FALSE){
  if(sum(score>=threshold, na.rm=TRUE)==0) return(0)
  if(!is.logical(type)) type <- type=="real"
  sum(type & score>=threshold, na.rm=TRUE)/sum(score>=threshold, na.rm=TRUE)
}

.prop.dev <- function(type, score, expected, threshold){
  if(!is.logical(type)) type <- type=="real"
  x <- 1+sum(score>=threshold & type)
  expected <- expected + 1
  if(length(expected)>1 && x>min(expected) && x<max(expected)) return(0)
  min(abs(x-expected)/expected)
}

.compareToExpectedDoublets <- function(origins, expected, dbr=NULL, od=NULL){
  if(length(origins)==0){
    d <- data.frame(combination=names(expected), 
                    observed=rep(0,length(expected)))
    d$deviation <- d$expected <- as.numeric(expected)
  }else{
    n <- names(expected)
    o <- as.numeric(table(origins)[names(expected)])
    o[is.na(o)] <- 0 
    expected <- as.numeric(expected)
    d <- data.frame( combination=n, observed=o, expected=expected,
                     deviation=abs(expected-o))
  }
  d$prop.deviation <- d$deviation/sum(expected)
  if(!is.null(od)) d$p <- dnbinom(d$observed, mu=d$expected, size=1/od)
  d
}
  
.getThresholdStats <- function(x, expected=NULL, thresholds=(0:100/100), 
                               dbr=NULL, od=0.05){
  if(is(x,"SingleCellExperiment")){
    x <- as.data.frame(colData(x))
    if(is.null(x$scDblFinder.type)) x$scDblFinder.type <- "real"
    x <- x[,grep("^scDblFinder",colnames(x))]
    colnames(x) <- gsub("^scDblFinder\\.","",colnames(x))
  }
  real <- which(x$type=="real")
  if(is.null(expected)) expected <- getExpectedDoublets(x$clusters[w], dbr)
  if(is.null(names(thresholds))) names(thresholds) <- thresholds
  a <- lapply( thresholds, FUN=function(i){
    .compareToExpectedDoublets(x$mostLikelyOrigin[which(x$score>=i & 
                                                            x$type=="real")], 
                               expected, dbr=dbr, od=od)
  })
  a <- cbind(threshold=rep(names(a),vapply(a,nrow,integer(1))),
             do.call(rbind, a))
  if(!all(x$type=="real")){
    b <- lapply( thresholds, FUN=function(i){
      y <- vapply(split(x, x$mostLikelyOrigin), FUN.VALUE=numeric(1L), 
                  FUN=function(x) .FNR(x$type, x$score, i) )
      data.frame(combination=names(y), FNR=y)
    })
    b <- cbind(threshold=rep(names(b),vapply(b,nrow,integer(1))),
               do.call(rbind, b))
    a <- merge(a,b,by=c("combination","threshold"),all.x=TRUE)
  }
  a$threshold <- as.numeric(a$threshold)
  a
}

