# This tests the recoverDoublets function.
# library(scDblFinder); library(testthat); source("test-recoverDoublets.R")

set.seed(99000077)
ngenes <- 100
mu1 <- 2^rexp(ngenes) * 5
mu2 <- 2^rnorm(ngenes) * 5

counts.1 <- matrix(rpois(ngenes*100, mu1), nrow=ngenes)
counts.2 <- matrix(rpois(ngenes*100, mu2), nrow=ngenes)
counts.m <- matrix(rpois(ngenes*20, mu1+mu2), nrow=ngenes)

counts <- cbind(counts.1, counts.2, counts.m)
clusters <- rep(1:3, c(ncol(counts.1), ncol(counts.2), ncol(counts.m)))

library(SingleCellExperiment)
sce <- SingleCellExperiment(list(counts=counts))
sce <- scuttle::logNormCounts(sce)

set.seed(99000007)
test_that("recoverDoublets works as expected", {
    known.doublets <- clusters==3 & rbinom(length(clusters), 1, 0.5)==0
    ref <- recoverDoublets(sce, known.doublets, samples=c(1, 1, 1))

    expect_true(min(ref$proportion[ref$predicted]) >= max(ref$proportion[!ref$predicted & !ref$known]))
    expect_false(any(ref$predicted & ref$known))
    expect_true(sum(ref$predicted) <= metadata(ref)$intra)

    # Responds to 'k'.
    alt <- recoverDoublets(sce, known.doublets, samples=c(1, 1, 1), k=20)
    expect_false(identical(ref, alt))

    # Responds to 'samples'
    alt <- recoverDoublets(sce, known.doublets, samples=c(1, 2, 3))
    expect_false(identical(ref, alt))

    # subset.row has the intended effect
    sub <- recoverDoublets(assay(sce), known.doublets, samples=c(1, 1, 1), subset.row=1:50)
    alt <- recoverDoublets(assay(sce)[1:50,], known.doublets, samples=c(1, 1, 1))
    expect_identical(sub, alt)
})

set.seed(99000008)
test_that("recoverDoublets gives the correct results on the toy example", {
    known.doublets <- clusters==3 & 1:2==1 # alternating doublets.
    ref <- recoverDoublets(sce, known.doublets, samples=c(1, 1), k=10)
    expect_identical(clusters==3, ref$known | ref$predicted)

    expect_true(min(ref$proportion[ref$predicted]) >= max(ref$proportion[!ref$predicted & !ref$known]))
})

set.seed(99000008)
test_that("recoverDoublets works for other inputs", {
    known.doublets <- clusters==3 & rbinom(length(clusters), 1, 0.5)==0
    ref <- recoverDoublets(logcounts(sce), known.doublets, samples=c(1, 1, 1))
    alt <- recoverDoublets(sce, known.doublets, samples=c(1, 1, 1))
    expect_identical(ref, alt)

    # Works for transposition
    alt <- recoverDoublets(t(logcounts(sce)), known.doublets, samples=c(1, 1, 1), transposed=TRUE)
    expect_identical(ref, alt)

    # Works by stuffing values in reduced dims.
    reducedDim(sce, "pretend") <- t(logcounts(sce))
    alt <- recoverDoublets(sce, known.doublets, samples=c(1, 1, 1), use.dimred="pretend")
    expect_identical(ref, alt)
})
