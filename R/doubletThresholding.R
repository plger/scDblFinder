#' doubletThresholding
#'
#' Sets the doublet scores threshold; typically called by
#' \code{\link[scDblFinder]{scDblFinder}}.
#'
#' @param d A data.frame of cell properties, with each row representing a cell, as
#' produced by `scDblFinder(..., returnType="table")`, or minimally containing a `score`
#' column.
#' @param dbr The expected (mean) doublet rate. If `d` contains a `cluster` column, the
#' doublet rate will be adjusted for homotypic doublets.
#' @param dbr.sd The standard deviation of the doublet rate, representing the
#' uncertainty in the estimate. Ignored if `method!="optim"`.
#' @param stringency A numeric value >0 and <1 which controls the relative weight of false
#'  positives (i.e. real cells) and false negatives (artificial doublets) in setting the
#'  threshold. A value of 0.5 gives equal weight to both; a higher value (e.g. 0.7) gives
#'  higher weight to the false positives, and a lower to artificial doublets.  Ignored if
#'  `method!="optim"`.
#' @param method The thresholding method to use, either 'auto' (default, automatic
#' selection depending on the available fields), 'optim' (optimization of
#' misclassification rate and deviation from expected doublet rate), 'dbr' (strictly
#' based on the expected doublet rate), or 'griffiths' (cluster-wise number of
#' median absolute deviation in doublet score).
#' @param perSample Logical; whether to perform thresholding individually for each sample.
#' @param p The p-value threshold determining the deviation in doublet score.
#' @param returnType The type of value to return, either doublet calls (`call`) or
#' thresholds (`threshold`).
#'
#' @return A vector of doublet calls if `returnType=="call"`, or a threshold (or vector
#' of thresholds) if `returnType=="threshold"`.
#'
#' @importFrom stats pcauchy optimize ecdf lm predict dnbinom quantile
#'
#' @examples
#' sce <- mockDoubletSCE()
#' d <- scDblFinder(sce, verbose=FALSE, returnType="table")
#' th <- doubletThresholding(d, dbr=0.05)
#' th
#'
#' @importFrom stats mad qnorm setNames
#' @export
doubletThresholding <- function( d, dbr=NULL, dbr.sd=NULL, stringency=0.5, p=0.1,
                                 method=c("auto","optim","dbr","griffiths"),
                                 perSample=TRUE, returnType=c("threshold","call")){
  method <- match.arg(method)
  returnType <- match.arg(returnType)
  if(is.null(d$src)) d$src <- d$type
  if(is.null(dbr.sd)) dbr.sd <- mean(0.4*.gdbr(d,dbr))
  dbr <- .estimateHeterotypicDbRate(d, .checkPropArg(dbr))
  if(!is.data.frame(d) && !is(d,"DFrame"))
    stop("`d` should be a data.frame with minimally the 'score' column.")
  conds <- list("optim"=c("type","score"),
                "dbr"=c("score"),
                "griffiths"=c("score"))
  w <- vapply(conds, FUN.VALUE=logical(1), FUN=function(x) all(x %in% colnames(d)))
  if(method=="auto"){
    if(length(w)==0) stop("`d` misses the necessary columns.")
    method <- names(conds)[which(w)[1]]
  }else{
    if(!w[[method]]) stop("`d` misses the necessary columns.")
  }
  if(method=="optim"){
    if(is.null(d$cluster)) d$cluster <- 1L
    if(!all(sort(as.character(unique(d$type)))==c("doublet","real")))
      stop("`type` should be either 'real' or 'doublet'.")
    if(is.null(d$include.in.training)) d$include.in.training <- TRUE
    if(!is.null(d$sample) && perSample){
      si <- split(seq_len(nrow(d)), d$sample)
      if(!is.null(dbr)){
        if(length(dbr)==1) dbr <- setNames(rep(dbr, length(si)), names(si))
        if(!all(names(si) %in% names(dbr)))
          stop("The names of `dbr` do not correspond to samples of `d`")
      }
      th <- sapply(setNames(names(si),names(si)), FUN=function(s){
        .optimThreshold(d[si[[s]],c("type","src","score","cluster","include.in.training")],
                        dbr=dbr[[s]], dbr.sd=dbr.sd, stringency=stringency)
      })
      ret <- as.factor(d$score > th[d$sample])
    }else{
      th <- .optimThreshold(d, dbr=.gdbr(d,dbr), dbr.sd=dbr.sd, stringency=stringency)
      ret <- as.factor(d$score>th)
    }
    if(returnType=="threshold") return(th)
  }else{
    if(!is.null(d$src)) d <- d[d$src=="real",]
    if(method=="dbr"){
      th <- quantile(d$score, 1-dbr)
      if(returnType=="threshold") return(th)
      ret <- as.factor(d$score>th)
    }else if(method=="griffiths"){
      if(!("sample" %in% colnames(d))) d$sample <- "all"
      i <- split(seq_len(nrow(d)), d$sample)
      meds <- vapply(i, FUN.VALUE=numeric(1), FUN=function(x){
        median(d$score[x],na.rm=TRUE)
      })
      d$dev <- d$score-meds[d$sample]
      mad <- vapply(i, FUN.VALUE=numeric(1), FUN=function(x){
        x <- d$dev[x]
        median(x[x>0],na.rm=TRUE)
      }) * formals(mad)$constant
      if(returnType=="threshold"){
        return(qnorm(p, mean=meds, sd=mad, lower.tail=FALSE))
      }else{
        d$p <- pnorm(d$score, mean=meds[d$sample], sd=mad[d$sample],
                     lower.tail=FALSE)
        ret <- as.factor(d$p < p)
      }
    }else{
      stop("Unknown method '",method,"'")
    }
  }
  levels(ret) <- c("singlet","doublet")
  return(ret)
}

# dbr should be already corrected for homotypy
.optimThreshold <- function(d, dbr=NULL, dbr.sd=NULL, ths=NULL, stringency=0.5){
  if(!(stringency > 0) || !(stringency<1))
    stop("`stringency` should be >0 and <1.")
  if(is.null(dbr)) dbr <- .gdbr(d, dbr=.estimateHeterotypicDbRate(d))
  if(!is.null(dbr.sd)) dbr <- c(max(0,dbr-dbr.sd), min(1,dbr+dbr.sd))
  if(is.null(d$cluster)) d$cluster <- 1L
  wR <- which(d$src=="real")
  expected <- dbr*length(wR)
  if(!is.logical(d$type)) d$type <- d$type=="real"
  fdr.include <- which(d$include.in.training)
  eFN <- sum(grepl("^rDbl\\.",row.names(d))) *
    propHomotypic(d$cluster[d$src=="real"])
  if(length(unique(d$cluster))==1) eFN <- 0
  totfn <- function(x){
    edev <- .prop.dev(d$type,d$score,expected,x)^2
    y <- edev + 2*(1-stringency)*.FNR(d$type, d$score, x, expectedFN=eFN)
    if(!is.null(fdr.include))
      y <- y + .FPR(d$type[fdr.include], d$score[fdr.include], x)*2*stringency
    y
  }
  if(is.null(ths)) return(optimize(totfn, c(0,1), maximum=FALSE)$minimum)
  data.frame( threshold=ths,
              FNR=vapply(ths, FUN.VALUE=numeric(1), FUN=function(x){
                .FNR(d$type, d$score, x, expectedFN=eFN)
              }),
              FDR=vapply(ths, FUN.VALUE=numeric(1), FUN=function(x){
                .FDR(d$type[fdr.include], d$score[fdr.include], x)
              }),
              FPR=vapply(ths, FUN.VALUE=numeric(1), FUN=function(x){
                .FPR(d$type[fdr.include], d$score[fdr.include], x)
              }),
              dev=vapply(ths, FUN.VALUE=numeric(1), FUN=function(x){
                .prop.dev(d$type,d$score,expected,x)^2
              }),
              cost=vapply(ths, FUN.VALUE=numeric(1), FUN=totfn)
  )
}

.FNR <- function(type, score, threshold, expectedFN=0){
  if(!is.logical(type)) type <- type=="real"
  max(c(0, (sum(!type & score<threshold, na.rm=TRUE)-expectedFN)))/sum(!type)
}
.FDR <- function(type, score, threshold){
  if(length(type)==0 || sum(score>=threshold, na.rm=TRUE)==0) return(0)
  if(!is.logical(type)) type <- type=="real"
  sum(type & score>=threshold, na.rm=TRUE)/sum(score>=threshold, na.rm=TRUE)
}
.FPR <- function(type, score, threshold){
  if(length(type)==0) return(0)
  if(!is.logical(type)) type <- type=="real"
  sum(type & score>=threshold, na.rm=TRUE)/sum(type)
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

.getDoubletStats <- function( d, th, dbr=NULL, dbr.sd=0.015 ){
  # check that we have all necessary fields:
  fields <- c("cluster","src","type","score","mostLikelyOrigin", "originAmbiguous","difficulty")
  if(!all(fields %in% colnames(d))) stop("Input misses some columns.")
  if(!is.null(d$sample))
    return(dplyr::bind_rows(lapply(split(seq_len(nrow(d)), d$sample), FUN=function(i){
      .getDoubletStats(d[i,fields], th, dbr=dbr, dbr.sd=dbr.sd)
    }), .id="sample"))
  if(is.null(dbr)) dbr <- 0.01*sum(d$src=="real",na.rm=TRUE)/1000
  o <- d$mostLikelyOrigin[d$type=="real" & d$score>=th]
  expected <- getExpectedDoublets(d$cluster[d$src=="real"], dbr=dbr)
  stats <- .compareToExpectedDoublets(o, dbr=dbr, expected=expected)
  stats$combination <- as.character(stats$combination)
  stats$FNR <- vapply(split(as.data.frame(d), d$mostLikelyOrigin),
                      FUN.VALUE=numeric(1L),
                      FUN=function(x) .FNR(x$type, x$score, th) )[stats$combination]
  stats$difficulty <- vapply(split(d$difficulty, d$mostLikelyOrigin),
                             FUN.VALUE=numeric(1L), na.rm=TRUE,
                             FUN=median)[stats$combination]
  stats
}

#' @importFrom stats quantile
.filterUnrecognizableDoublets <- function( d, minSize=5, minMedDiff=0.1 ){
  if(is.null(d$src)) d$src <- d$type
  da <- d[d$src=="artificial" & grepl("+", d$mostLikelyOrigin, fixed=TRUE),]
  dr <- d[d$src=="real",]
  dr.med <- median(dr$score)
  dr <- split(dr$score, dr$cluster)
  rq <- t(vapply(dr, FUN.VALUE=numeric(2), na.rm=TRUE, probs=c(0.5, 0.9), FUN=quantile))
  da <- split(da$score, droplevels(da$mostLikelyOrigin))
  origs <- strsplit(names(da), "+", fixed=TRUE)
  out <- vapply(names(da), FUN.VALUE=logical(1), FUN=function(x){
    z <- da[[x]]
    if(length(z)<minSize) return(FALSE)
    z <- quantile(z, probs=c(0.1,0.5), na.rm=TRUE)
    x <- origs[[x]]
    any(z[1]<rq[x,2]) || (z[2]-max(dr.med, rq[x,1]))<minMedDiff
  })
  out <- names(out)[which(out)]
  d[d$src!="artificial" | !(d$mostLikelyOrigin %in% out),]
}

.estimateHeterotypicDbRate <- function(d, dbr=NULL){
  d <- as.data.frame(d)
  if(!is.null(d$sample)){
    sd <- split(d[,setdiff(colnames(d),"sample")],d$sample)
    if(length(dbr)==1) dbr <- rep(dbr,length(sd))
    if(is.null(dbr)) dbr <- unlist(lapply(sd, FUN=.gdbr))
    return(unlist(mapply(d=sd, dbr=dbr, FUN=.estimateHeterotypicDbRate)))
  }
  dbr <- .gdbr(d, dbr=dbr)
  if(is.null(d$cluster)){
    th <- sum(d$src=="artificial")/nrow(d)
    prop.homo <- sum(d$src=="artificial" & d$score<th)/sum(d$src=="artificial")
    return(dbr * (1-prop.homo))
  }
  if(is.factor(d$cluster)) d$cluster <- droplevels(d$cluster)
  d <- d[d$src=="real",]
  sum(getExpectedDoublets(d$cluster, dbr=dbr))/nrow(d)
}