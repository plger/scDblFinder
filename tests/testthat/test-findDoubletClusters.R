# This tests the cluster-based doublet discovery machinery. 
# library(scDblFinder); library(testthat); source("test-findDoubletClusters.R")

set.seed(9900001)
ngenes <- 100
mu1 <- 2^rexp(ngenes) 
mu2 <- 2^rnorm(ngenes)

counts.1 <- matrix(rpois(ngenes*100, mu1), nrow=ngenes)
counts.2 <- matrix(rpois(ngenes*100, mu2), nrow=ngenes)
counts.m <- matrix(rpois(ngenes*20, mu1+mu2), nrow=ngenes)

counts <- cbind(counts.1, counts.2, counts.m)
clusters <- rep(1:3, c(ncol(counts.1), ncol(counts.2), ncol(counts.m)))

RENAMER <- function(val, fields, mapping) 
# A convenience function for remapping internal fields upon subsetting or renaming.
# This is necessary for some equality checks below.
{
    new.pairs <- val$all.pairs
    for (f in fields) {
        val[[f]] <- mapping[as.integer(val[[f]])]
        for (i in seq_along(new.pairs)) {
            new.pairs[[i]][[f]] <- mapping[as.integer(new.pairs[[i]][[f]])]
        }
    }
    val$all.pairs <- new.pairs
    val
}

test_that("findDoubletClusters works correctly with vanilla tests", {
    dbl <- findDoubletClusters(counts, clusters)
    expect_identical(rownames(dbl)[1], "3")
    expect_identical(dbl$source1[1], "2")
    expect_identical(dbl$source2[1], "1")

    # Checking the relative library sizes.
    ls1 <- median(colSums(counts.1))
    ls2 <- median(colSums(counts.2))
    ls3 <- median(colSums(counts.m))

    expect_equal(dbl$lib.size1[1], ls2/ls3)
    expect_equal(dbl$lib.size2[1], ls1/ls3)

    # Checking the proportions.
    expect_equal(dbl$prop, as.integer(table(clusters)[rownames(dbl)])/length(clusters))

    # Checking that p-values are reverse-sorted.
    expect_false(is.unsorted(-dbl$p.value))

    # Checking that we get equivalent results with character cluster input.
    re.clusters <- LETTERS[clusters]
    re.dbl <- findDoubletClusters(counts, re.clusters)

    dbl2  <- RENAMER(dbl, c("source1", "source2"), LETTERS)
    rownames(dbl2) <- LETTERS[as.integer(rownames(dbl2))]
    expect_identical(dbl2, re.dbl)
})

test_that("findDoubletClusters agrees with a reference implementation", {
    mu3 <- 2^rnorm(ngenes)
    counts.3 <- matrix(rpois(ngenes*100, mu3), nrow=ngenes)
    counts <- cbind(counts.1, counts.2, counts.3, counts.m)
    clusters <- rep(1:4, c(ncol(counts.1), ncol(counts.2), ncol(counts.3), ncol(counts.m)))

    dbl <- findDoubletClusters(counts, clusters, get.all.pairs=TRUE)
    ref <- scran::findMarkers(scuttle::normalizeCounts(counts), clusters, full.stats=TRUE)

    for (x in rownames(dbl)) {
        stats <- ref[[x]]
        all.pops <- setdiff(rownames(dbl), x)
        combos <- combn(all.pops, 2)

        # Effectively a re-implentation of the two inner loops.
        collected <- apply(combos, 2, function(chosen) {
            fields <- paste0("stats.", chosen)
            stats1 <- stats[[fields[1]]]
            stats2 <- stats[[fields[2]]]
            p <- pmax(exp(stats1$log.p.value), exp(stats2$log.p.value))
            p[sign(stats1$logFC)!=sign(stats2$logFC)] <- 1
            adj.p <- p.adjust(p, method="BH")
            data.frame(best=rownames(stats)[which.min(p)], p.val=min(adj.p), 
                num.de=sum(adj.p <= 0.05), stringsAsFactors=FALSE)
        })

        collected <- do.call(rbind, collected)
        o <- order(collected$num.de, -collected$p.val)

        obs <- dbl[x,"all.pairs"][[1]]
        expect_identical(obs$source1, pmax(combos[2,], combos[1,])[o])
        expect_identical(obs$source2, pmin(combos[1,], combos[2,])[o])
        expect_identical(obs$num.de, collected$num.de[o])
        expect_identical(obs$best, collected$best[o])
        expect_equal(obs$p.value, collected$p.val[o])

        to.use <- o[1]
        expect_identical(dbl[x,"num.de"], collected[to.use, "num.de"])
        expect_equal(dbl[x,"p.value"], collected[to.use, "p.val"])
        expect_identical(dbl[x,"best"], collected[to.use, "best"])
        expect_identical(sort(c(dbl[x,"source1"],dbl[x,"source2"])), sort(combos[,to.use]))
    }
})

test_that("findDoubletClusters works correctly with row subsets", {
    chosen <- sample(ngenes, 20)
    dbl0 <- findDoubletClusters(counts, clusters, subset.row=chosen)
    ref <- findDoubletClusters(counts[chosen,], clusters)
    ref <- RENAMER(ref, "best", as.character(chosen))
    expect_identical(dbl0, ref)

    # Trying out empty rows.
    out <- findDoubletClusters(counts[0,], clusters)
    expect_identical(nrow(out), nrow(ref))
    expect_true(all(is.na(out$best)))
    expect_true(all(is.na(out$p.value)))
    expect_true(all(out$num.de==0L))

    # While we're here, trying out empty columns.
    expect_error(findDoubletClusters(counts[,0], clusters[0]), "need at least three")
})

test_that("findDoubletClusters works correctly with SE/SCEs", {
    library(SingleCellExperiment)
    sce <- SingleCellExperiment(list(counts=counts))
    ref <- findDoubletClusters(counts, clusters)

    dbl <- findDoubletClusters(sce, clusters)
    expect_identical(ref, dbl)

    # Works with the base class.
    dbl2 <- findDoubletClusters(as(sce, "SummarizedExperiment"), clusters)
    expect_identical(ref, dbl2)

    # Works with column labels.
    colLabels(sce) <- clusters
    dbl3 <- findDoubletClusters(sce)
    expect_identical(ref, dbl3)

    # With a different assay.
    assay(sce, "whee") <- counts + rpois(length(counts), lambda=2)
    ref2 <- findDoubletClusters(assay(sce, "whee"), clusters)
    dbl2 <- findDoubletClusters(sce, clusters, assay.type="whee")
    expect_identical(ref2, dbl2)
})
