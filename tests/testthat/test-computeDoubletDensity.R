# This tests the doublet density machinery.
# library(scDblFinder); library(testthat); source("test-computeDoubletDensity.R")

set.seed(9900001)
ngenes <- 100
mu1 <- 2^rexp(ngenes) 
mu2 <- 2^rnorm(ngenes)

counts.1 <- matrix(rpois(ngenes*100, mu1), nrow=ngenes)
counts.2 <- matrix(rpois(ngenes*100, mu2), nrow=ngenes)
counts.m <- matrix(rpois(ngenes*20, mu1+mu2), nrow=ngenes)

counts <- cbind(counts.1, counts.2, counts.m)
clusters <- rep(1:3, c(ncol(counts.1), ncol(counts.2), ncol(counts.m)))

set.seed(9900002)
test_that("computeDoubletDensity PC spawning works correctly", {
    sf <- runif(ncol(counts))
    y <- log2(t(t(counts)/sf)+1)
    centers <- rowMeans(y)
    SVD <- svd(t(y - centers), nv=20)

    set.seed(12345)
    sim.pcs <- scDblFinder:::.spawn_doublet_pcs(counts, sf, SVD$v, centers, niters=10000L, block=10000L)

    set.seed(12345)
    L <- sample(ncol(counts), 10000L, replace=TRUE)
    R <- sample(ncol(counts), 10000L, replace=TRUE)
    ref.x <- counts[,L] + counts[,R]
    ref.y <- log2(t(t(ref.x)/(sf[L] + sf[R]))+1)
    ref.pcs <- crossprod(ref.y - centers, SVD$v)

    expect_equal(sim.pcs, ref.pcs)

    # Works with multiple iterations.
    set.seed(23456)
    sim.pcs <- scDblFinder:::.spawn_doublet_pcs(counts, sf, SVD$v, centers, niters=25000L, block=10000L)

    set.seed(23456)
    ref1 <- scDblFinder:::.spawn_doublet_pcs(counts, sf, SVD$v, centers, niters=10000L, block=10000L)
    ref2 <- scDblFinder:::.spawn_doublet_pcs(counts, sf, SVD$v, centers, niters=10000L, block=10000L)
    ref3 <- scDblFinder:::.spawn_doublet_pcs(counts, sf, SVD$v, centers, niters=5000L, block=10000L)

    expect_equal(sim.pcs, rbind(ref1, ref2, ref3))
    expect_identical(dim(sim.pcs), c(25000L, ncol(SVD$v)))
})

set.seed(9900003)
test_that("size factor variations in computeDoubletDensity work correctly", {
    # Library sizes get used.
    set.seed(12345)
    out <- computeDoubletDensity(counts)
    set.seed(12345)
    ref <- computeDoubletDensity(counts, size.factors.norm=scuttle::librarySizeFactors(counts))
    expect_equal(out, ref)

    # Normalization size factors get centered.
    sf1 <- runif(ncol(counts))
    set.seed(23456)
    out <- computeDoubletDensity(counts, size.factors.norm=sf1)
    set.seed(23456)
    ref <- computeDoubletDensity(counts, size.factors.norm=sf1/mean(sf1))
    expect_equal(out, ref)

    # Reacts correctly to size.factors.content.
    sf1 <- sf1/mean(sf1)
    sf2 <- runif(ncol(counts))

    set.seed(23456)
    ref <- computeDoubletDensity(counts, size.factors.norm=sf1)

    set.seed(23456)
    out <- computeDoubletDensity(t(t(counts)/sf1), size.factors.norm=rep(1, ncol(counts)), size.factors.content=1/sf1)
    expect_equal(out, ref)

    # take the product, which gets divided out by 's2' to give back 's1' during the actual normalization.
    set.seed(23456)
    prod <- sf1*sf2 
    scaled <- t(t(counts)*sf2)/mean(prod) 
    out <- computeDoubletDensity(scaled, size.factors.norm=prod, size.factors.content=sf2)
    expect_equal(out, ref)

    # scaling of content size factors don't matter.
    set.seed(23456)
    out <- computeDoubletDensity(scaled, size.factors.norm=prod, size.factors.content=sf2*5) 
    expect_equal(out, ref)
})

set.seed(9900004)
test_that("high-level tests for computeDoubletDensity work correctly", {
    mu1 <- 2^rnorm(ngenes) * 100 # using a really high count to reduce variance.
    mu2 <- 2^rnorm(ngenes) * 100
    ncA <- 100
    ncB <- 100
    ncC <- 51

    counts.A <- matrix(rpois(ngenes*ncA, mu1), ncol=ncA, nrow=ngenes)
    counts.B <- matrix(rpois(ngenes*ncB, mu2), ncol=ncB, nrow=ngenes)
    counts.C <- matrix(rpois(ngenes*ncC, mu1+mu2), ncol=ncC, nrow=ngenes)
    clusters <- rep(1:3, c(ncA, ncB, ncC))

    out <- computeDoubletDensity(cbind(counts.A, counts.B, counts.C))
    expect_true(min(out[clusters==3]) > max(out[clusters!=3]))

    # Now with differences in RNA content.
    counts.A <- matrix(rpois(ngenes*ncA, mu1), ncol=ncA, nrow=ngenes)
    counts.B <- matrix(rpois(ngenes*ncB, mu2), ncol=ncB, nrow=ngenes)
    counts.C <- matrix(rpois(ngenes*ncC, (mu1+2*mu2)/3), ncol=ncC, nrow=ngenes)
    sf.spike <- 1/rep(1:3, c(ncA, ncB, ncC))
    
    X <- cbind(counts.A, counts.B, counts.C) 
    out <- computeDoubletDensity(X, size.factors.content=sf.spike)
    expect_true(min(out[clusters==3]) > max(out[clusters!=3]))

    out <- computeDoubletDensity(X) # fails without size factor info; differences are basically negligible.
    expect_true(max(out[clusters==3]) < min(out[clusters!=3]))
})

set.seed(9900005)
test_that("other settings for computeDoubletDensity work correctly", {
    # Subsetting behaves correctly.
    set.seed(1000)
    sim <- computeDoubletDensity(counts, subset.row=1:50)
    set.seed(1000)
    ref <- computeDoubletDensity(counts[1:50,])
    expect_identical(sim, ref)

    # Warnings raised if too many neighbors are requested.
    expect_warning(computeDoubletDensity(counts, k=1000), "'k' capped")

    # IRLBA works correctly.
    set.seed(2000)
    sim <- computeDoubletDensity(counts, d=5)
    set.seed(2000)
    ref <- computeDoubletDensity(counts, BSPARAM=BiocSingular::IrlbaParam(tol=1e-12, extra.work=50, maxit=20000), d=5)
    expect_true(median( abs(sim-ref)/(sim+ref+1e-6) ) < 0.01)

    # Alternative neighbor search method works correctly.
    expect_error(sim <- computeDoubletDensity(counts, BNPARAM=BiocNeighbors::VptreeParam()), NA)

    # Responds correctly to blocking.
    set.seed(3000)
    ref <- computeDoubletDensity(counts)
    sim1 <- computeDoubletDensity(counts, block=1000)
    expect_equal(log1p(sim1), log1p(ref), tol=0.1)
    sim2 <- computeDoubletDensity(counts, niters=20000)
    expect_equal(log1p(sim2), log1p(ref), tol=0.1)
})

set.seed(9900006)
test_that("computeDoubletDensity works correctly for SCE objects", {
    library(SingleCellExperiment)
    sce <- SingleCellExperiment(list(counts=counts))

    set.seed(1000)
    ref <- computeDoubletDensity(counts)
    set.seed(1000)
    dbl <- computeDoubletDensity(sce)
    expect_identical(ref, dbl)

    # With a different assay.
    assay(sce, "whee") <- counts + rpois(length(counts), lambda=2)
    set.seed(1001)
    ref2 <- computeDoubletDensity(assay(sce, "whee"))
    set.seed(1001)
    dbl2 <- computeDoubletDensity(sce, assay.type="whee")
    expect_identical(ref2, dbl2)

    # With subsetting.
    keep <- sample(nrow(sce), 10)

    set.seed(1003)
    dbl5 <- computeDoubletDensity(sce, subset.row=keep)
    set.seed(1003)
    ref4 <- computeDoubletDensity(sce[keep,])
    expect_identical(ref4, dbl5)
})
