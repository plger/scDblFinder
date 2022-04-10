sce <- mockDoubletSCE(ncells=c(100,200,150,100), ngenes=250)
sce$fastcluster <- fastcluster(sce, nfeatures=100, verbose=FALSE)
sce$sample <- sample(LETTERS[1:2], ncol(sce), replace=TRUE)

test_that("fastcluster works as expected",{
  expect_equal(sum(is.na(sce$fastcluster)),0)
  expect_gt(sum(apply(table(sce$cluster, sce$fastcluster),1,max)[1:4])/
                  sum(sce$type=="singlet"), 0.8)
  x <- fastcluster(sce, nfeatures=100, k=3, verbose=FALSE, return="preclusters")
  expect_equal(sum(is.na(x)),0)
  expect_gt(sum(apply(table(sce$cluster, x),1,max)[1:3])/
              sum(sce$type=="singlet"), 0.8)

})

sce <- scDblFinder(sce, clusters="fastcluster", samples="sample",
                   artificialDoublets=250, dbr=0.1, verbose=FALSE)

test_that("scDblFinder works as expected", {
    expect_equal(sum(is.na(sce$scDblFinder.score)),0)
    expect(min(sce$scDblFinder.score)>=0 & max(sce$scDblFinder.score)<=1,
           failure_message="scDblFinder.score not within 0-1")
    expect_gt(sum(sce$type==sce$scDblFinder.class)/ncol(sce), 0.8)
    sce <- scDblFinder(sce, samples="sample", artificialDoublets=250, 
                       dbr=0.1, verbose=FALSE)
    expect_equal(sum(is.na(sce$scDblFinder.score)),0)
    expect(min(sce$scDblFinder.score)>=0 & max(sce$scDblFinder.score)<=1,
           failure_message="scDblFinder.score not within 0-1")
    expect_gt(sum(sce$type==sce$scDblFinder.class)/ncol(sce), 0.8)
})

test_that("feature aggregation works as expected", {
  sce2 <- aggregateFeatures(sce, k=20)
  expect_equal(nrow(sce2),20)
  expect_equal(sum(is.na(counts(sce2)) | is.infinite(counts(sce2))), 0)
  sce2 <- scDblFinder( sce2, clusters="fastcluster", processing="normFeatures",
                      artificialDoublets=250, dbr=0.1, verbose=FALSE)
  expect_equal(sum(is.na(sce2$scDblFinder.score)),0)
  expect_gt(sum(sce2$type==sce2$scDblFinder.class)/ncol(sce2), 0.8)
})

test_that("doublet enrichment works as expected", {
  cs <- clusterStickiness(sce)$FDR
  expect_equal(sum(is.na(cs)),0)
})


test_that("amulet works as expected", {
  fragfile <- system.file("extdata","example_fragments.tsv.gz",
                          package="scDblFinder")
  res <- amulet(fragfile)
  expect_equal(res$nFrags, c(878,2401,2325,1882,1355))
  expect_equal(sum(res$nAbove2<=1), 4)
  expect_equal(res["barcode5","nAbove2"], 6)
  expect_lt(res["barcode5","p.value"], 0.01)
})