set.seed(123)
sce <- mockDoubletSCE(ncells = c(40,40,50), ngenes = 500, dbl.rate = 0.3)
isDoublet <- which(sce$type=="doublet")

clf <- identifyDoubletOrigins(sce, "cluster", nArtificial=100, max_rounds=20)

test_that("identifyDoubletOrigins works", {
  tt <- clf$train_contigency
  expect_gt(sum(diag(tt))/sum(tt), 0.8)
})

test_that("predictDoubletOrigins works", {
  pred <- predictDoubletOrigins(clf, sce[,isDoublet])
  expect_true(all(colnames(sce)[isDoublet]==names(pred)))
  tt <- table(true=sce$origin[isDoublet], pred)
  acc <- sum(diag(tt))/sum(tt)
  print(acc)
  expect_gt(acc, 0.65)
})

test_that("per-sample statistics are returned with multiple samples", {
  set.seed(123)
  sce <- mockDoubletSCE(ncells = c(150,150,200), ngenes = 500, dbl.rate = 0.2)
  sce$sample <- sample(c("A","B"), ncol(sce), replace=TRUE)
  sce <- scDblFinder(sce, samples="sample", clusters="cluster", 
                     artificialDoublets=200, verbose=FALSE)
  expect_true("sample" %in% colnames(metadata(sce)$scDblFinder.stats))
  clf <- identifyDoubletOrigins(sce, "cluster", samples="sample", 
                                nArtificial=100, max_rounds=20, verbose=FALSE)
  ps <- clf$stats.per.sample
  expect_true(all(c("sample","combination","observed","expected") %in% colnames(ps)))
  expect_setequal(ps$sample, c("A","B"))
  nd <- table(sce$sample[sce$scDblFinder.class=="doublet"])
  expect_equal(as.numeric(tapply(ps$observed, ps$sample, sum)[names(nd)]), 
               as.numeric(nd))
  # aggregated statistics are consistent with the per-sample ones
  expect_equal(clf$stats$observed, as.numeric(tapply(ps$observed, ps$combination, sum)[clf$stats$combination]))
  expect_equal(clf$stats$expected, as.numeric(tapply(ps$expected, ps$combination, sum)[clf$stats$combination]))
  # not available without samples
  clf2 <- identifyDoubletOrigins(sce, "cluster", nArtificial=100, 
                                 max_rounds=20, verbose=FALSE)
  expect_null(clf2$stats.per.sample)
})
