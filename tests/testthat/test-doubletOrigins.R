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
