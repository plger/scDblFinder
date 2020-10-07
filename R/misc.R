.tablevec <- function(x){
  if(!is(x, "table")) x <- table(x)
  y <- as.numeric(x)
  names(y) <- names(x)
  y
}

.getMostLikelyOrigins <- function(knn, known.origins=NULL){
  origins <- t(vapply( seq_len(nrow(knn$orig)), FUN.VALUE=character(2),
                       FUN=function(i){
           if(all(is.na(knn$orig[i,]))) return(c(NA_character_,NA_character_))
           n <- .tablevec(knn$orig[i,])
           if(length(n)==1) return(c(names(n),FALSE))
           if(sum(n==max(n))==1) return(c(names(n)[which.max(n)],FALSE))
           x <- split(knn$distance[i,], knn$orig[i,])
           x <- sort(vapply(x,min,numeric(1)))
           if( (x[2]-x[1])/x[1] > 0.2 ) return(c(names(x)[1], FALSE))
           return(c(names(x)[1], TRUE))
                       }))
  origins <- as.data.frame(origins)
  colnames(origins) <- c("mostLikelyOrigin", "originAmbiguous")
  if(!is.null(known.origins)){
    w <- which(!is.na(known.origins))
    origins[w,1] <- as.character(known.origins[w])
    origins[w,2] <- "FALSE"
  }
  origins[,1] <- factor(origins[,1])
  origins[,2] <- as.logical(origins[,2])
  origins
}

.adjust.dbr.homotypy <- function(clusters, dbr=NULL){
  if(is.null(dbr)) dbr <- (0.01*length(clusters)/1000)
  # adjust expected dbr for expected homotypic doublets
  homotypic.prop <- sum((table(clusters)/length(clusters))^2)
  dbr*(1-homotypic.prop)
}

#' getExpectedDoublets
#'
#' @param x A vector of cluster labels for each cell
#' @param dbr The expected doublet rate.
#' @param only.heterotypic Logical; whether to return expectations only for 
#' heterotypic doublets
#'
#' @return The expected number of doublets of each combination of clusters
#'
#' @examples
#' # random cluster labels
#' cl <- sample(head(LETTERS,4), size=2000, prob=c(.4,.2,.2,.2), replace=TRUE)
#' getExpectedDoublets(cl)
#' @export
getExpectedDoublets <- function(x, dbr=NULL, only.heterotypic=TRUE){
  if(is(x,"SingleCellExperiment")){
    clusters <- x$scDblFinder.clusters
  }else{
    clusters <- x
  }
  clusters <- factor(clusters)
  ncells <- length(clusters)
  if(is.null(dbr)) dbr <- (0.01*ncells/1000)

  cs <- table(clusters)/ncells
  eg <- expand.grid(seq_along(cs), seq_along(cs))
  eg <- eg[eg[,1]<=eg[,2],]
  expected <- apply(eg,1, FUN=function(x){
    cs[[x[[1]]]]*cs[[x[[2]]]]
  })
  expected <- dbr*ncells*expected/sum(expected)
  names(expected) <- apply(eg,1,FUN=function(x){
    paste(names(cs)[x], collapse="+")
  })
  if(only.heterotypic) expected <- expected[which(eg[,1]!=eg[,2])]
  expected
}

.castorigins <- function(e, val=NULL){
  if(is.table(e) || is.null(dim(e))){
    e <- cbind(do.call(rbind, strsplit(names(e),"+",fixed=TRUE)),
               data.frame(val=as.numeric(e)))
  }
  if(is.null(val)) val <- rev(colnames(e))[1]
  names(n) <- n <- unique(as.character(as.matrix(e[,1:2])))
  vapply(n, FUN.VALUE=numeric(length(n)), FUN=function(x){
    vapply(n, FUN.VALUE=numeric(1), FUN=function(y){
      if(x==y) return(NA)
      w <- which(e[,1] %in% c(x,y) & e[,2] %in% c(x,y))
      if(length(w)==0) return(0)
      sum(e[w,val])
    })
  })
}

.clusterTopG <- function(sce, clusters=NULL, nfeatures=1000){
  if(nfeatures>=nrow(sce)) return(row.names(sce))
  if(is.null(clusters))
    return(row.names(sce)[order(Matrix::rowMeans(counts(sce),na.rm=TRUE),
                                decreasing=TRUE)[seq_len(nfeatures)]])
  # get mean expression across clusters
  cli <- split(seq_len(ncol(sce)), clusters)
  cl.means <- vapply(cli, FUN.VALUE=double(nrow(sce)), FUN=function(x){
    Matrix::rowMeans(counts(sce)[,x,drop=FALSE], na.rm=TRUE)
  })
  # grab the top genes in each cluster
  g <- unique(as.numeric(t(apply(cl.means, 2, FUN=function(x){
    order(x, decreasing=TRUE)[seq_len(nfeatures)]
  }))))[seq_len(nfeatures)]
  g
}

#' mockDoubletSCE
#' 
#' Creates a mock random single-cell experiment object with doublets
#'
#' @param ncells A positive integer vector indicating the number of cells per
#' cluster (min 2 clusters)
#' @param ngenes The number of genes to simulate. Ignored if `mus` is given.
#' @param mus A list of cluster averages.
#' @param dbl.rate The doublet rate
#' @param only.heterotypic Whether to create only heterotypic doublets
#'
#' @return A SingleCellExperiment object, with the colData columns `type` 
#' indicating whether the cell is a singlet or doublet, and `cluster` 
#' indicating from which cluster (or cluster combination) it was simulated.
#' 
#' @export
#' @import SingleCellExperiment
#' @importFrom stats rnorm rpois
#' @examples
#' sce <- mockDoubletSCE()
mockDoubletSCE <- function(ncells=c(200,300), ngenes=200, mus=NULL, 
                           dbl.rate=0.1, only.heterotypic=TRUE){
  if(length(ncells)<2)
    stop("ncells should be a positive integer vector of length >=2")
  if(is.null(names(ncells))) names(ncells)<-paste0("cluster",seq_along(ncells))
  if(is.null(mus)){
    mus <- lapply(ncells, FUN=function(x) 2^rnorm(ngenes))
  }else{
    if(!is.list(mus) || length(mus)!=length(ncells) || 
       length(unique(lengths(mus)))!=1)
      stop("If provided, `mus` should be a list of length equal to that of ",
           "`ncells`, with each slot containing a numeric vector of averages ",
           "for each gene.")
    names(mus) <- names(ncells)
  }

  # non-doublets
  counts <- do.call(cbind, lapply(seq_along(mus), FUN=function(i){
    matrix(rpois(ncells[[i]]*ngenes, mus[[i]]), nrow=ngenes)
  }))
  sce <- SingleCellExperiment(list(counts=counts), 
                              colData=data.frame(type="singlet", 
                                           origin=rep(names(ncells), ncells)))
  
  # doublets
  expected <- getExpectedDoublets(sce$origin, dbl.rate, 
                                  only.heterotypic=only.heterotypic)
  simdbl <- rpois(length(expected), expected)
  mus <- lapply(strsplit(names(expected),"+",fixed=TRUE), FUN=function(x){
    mus[[x[[1]]]]+mus[[x[[2]]]]
  })
  doublets <- do.call(cbind, lapply(seq_along(mus), FUN=function(i){
    matrix(rpois(simdbl[[i]]*ngenes, mus[[i]]), nrow=ngenes)
  }))
  sce2 <- SingleCellExperiment( list(counts=doublets), 
                                colData=data.frame(type="doublet", 
                                        origin=rep(names(expected), simdbl)))
  sce <- cbind(sce, sce2)
  
  # set original cluster for homotypic doublets
  sce$cluster <- factor(sce$origin, c(names(ncells),names(expected)))
  names(n) <- n <- levels(sce$cluster)
  n[paste(names(ncells),names(ncells),sep="+")] <- names(ncells)
  levels(sce$cluster) <- as.character(n)
  
  colnames(sce) <- paste0("cell",seq_len(ncol(sce)))
  row.names(sce) <- paste0("gene",seq_len(ngenes))
  sce$type <- factor(sce$type, c("singlet","doublet"))
  sce
}


.checkColArg <- function(sce, x, acceptNull=TRUE){
  arg <- deparse(substitute(x))
  if(is.null(x)){
    if(!acceptNull) stop("Missing argument `",arg,"`!")
    return(NULL)
  }
  if(is.character(x) && length(x)==1){
    if(!(x %in% colnames(colData(sce)))) 
      stop("Could not find `", arg, "` column in colData!")
    x <- colData(sce)[[x]]
  }else if(length(x)!=ncol(sce)){
    stop("`",arg,"` should have a length equal to the number of columns in ",
         "`sce`.")
  }
  x
}
