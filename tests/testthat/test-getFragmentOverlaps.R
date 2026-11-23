# Tests for Tabix fragment retrieval beyond 100 Mb, default exclusion width,
# and pthres handling in .removeHighOverlapSites.

test_that("getFragmentOverlaps reads Tabix fragments beyond 100 Mb and matches fullInMemory", {
  skip_if_not_installed("Rsamtools")
  skip_if_not_installed("rtracklayer")

  td <- tempfile("scdbl_tabix_")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  bed_f <- file.path(td, "fragments.tsv")
  lines <- c(
    # Barcode A: 3 overlapping fragments <100 Mb, 3 overlapping >150 Mb, 3 overlapping at ~248.9 Mb
    # Barcode B: fragments exclusively >100 Mb (sorted by chr, start for Tabix)
    # Barcode C: >100 Mb fragments on chrX (must be excluded by default regionsToExclude)
    "chr1\t1000\t1300\tCELL_A\t1",
    "chr1\t1050\t1350\tCELL_A\t1",
    "chr1\t1100\t1400\tCELL_A\t1",
    "chr1\t2000\t2300\tCELL_C\t1",
    "chr1\t150000000\t150000300\tCELL_A\t1",
    "chr1\t150000050\t150000350\tCELL_A\t1",
    "chr1\t150000100\t150000400\tCELL_A\t1",
    "chr1\t155000000\t155000300\tCELL_B\t1",
    "chr1\t155000050\t155000350\tCELL_B\t1",
    "chr1\t155000100\t155000400\tCELL_B\t1",
    "chr1\t248956000\t248956300\tCELL_A\t1",
    "chr1\t248956050\t248956350\tCELL_A\t1",
    "chr1\t248956100\t248956400\tCELL_A\t1",
    "chrX\t150000000\t150000300\tCELL_C\t1",
    "chrX\t150000050\t150000350\tCELL_C\t1",
    "chrX\t150000100\t150000400\tCELL_C\t1"
  )
  writeLines(lines, bed_f)

  bgz_f <- Rsamtools::bgzip(bed_f, dest = file.path(td, "fragments.tsv.gz"), overwrite = TRUE)
  Rsamtools::indexTabix(bgz_f, format = "bed")

  res_tbx <- getFragmentOverlaps(bgz_f, minFrags = 1L, removeHighOverlapSites = FALSE, verbose = FALSE)
  res_mem <- getFragmentOverlaps(bgz_f, minFrags = 1L, removeHighOverlapSites = FALSE, fullInMemory = TRUE, verbose = FALSE)

  expect_equal(res_tbx[c("CELL_A", "CELL_B", "CELL_C"), c("nFrags", "uniqFrags", "nAbove2")],
               res_mem[c("CELL_A", "CELL_B", "CELL_C"), c("nFrags", "uniqFrags", "nAbove2")])
  expect_equal(res_tbx["CELL_A", "nFrags"], 9L)
  expect_equal(res_tbx["CELL_A", "nAbove2"], 3L)
  expect_equal(res_tbx["CELL_B", "nFrags"], 3L)
  expect_equal(res_tbx["CELL_B", "nAbove2"], 1L)
  # chrX beyond 100 Mb is excluded by default regionsToExclude
  expect_equal(res_tbx["CELL_C", "nFrags"], 1L)
  expect_equal(res_tbx["CELL_C", "nAbove2"], 0L)

  loci_tbx <- getFragmentOverlaps(bgz_f, minFrags = 1L, removeHighOverlapSites = FALSE,
                                  ret = "loci", verbose = FALSE)
  expect_equal(length(loci_tbx), 4L)
})

test_that(".removeHighOverlapSites respects the pthres parameter", {
  # 4 cells share locus1 (hol=4) alongside 20 singleton loci (hol=1);
  # mean(hol) = 24/21 = 1.143 -> ppois(4, 1.143, lower.tail=FALSE) ~= 0.0060 < 0.01
  # Whereas for pthres = 1e-6 it must NOT be removed.
  gr <- c(
    GenomicRanges::GRanges("chr1", IRanges::IRanges(rep(100L, 4L), width = 50L)),
    GenomicRanges::GRanges("chr1", IRanges::IRanges(seq(1000L, by = 200L, length.out = 20L), width = 50L))
  )
  kept_strict <- scDblFinder:::.removeHighOverlapSites(gr, pthres = 0.01)
  kept_lenient <- scDblFinder:::.removeHighOverlapSites(gr, pthres = 1e-6)
  expect_equal(length(kept_strict), 20L)
  expect_equal(length(kept_lenient), 24L)
})
