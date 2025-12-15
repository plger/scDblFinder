# getFragmentOverlaps

Count the number of overlapping fragments.

## Usage

``` r
getFragmentOverlaps(
  x,
  barcodes = NULL,
  regionsToExclude = GRanges(c("M", "chrM", "MT", "X", "Y", "chrX", "chrY"), IRanges(1L,
    width = 10^8)),
  minFrags = 500L,
  uniqueFrags = TRUE,
  maxFragSize = 1000L,
  removeHighOverlapSites = TRUE,
  fullInMemory = FALSE,
  BPPARAM = NULL,
  verbose = TRUE,
  ret = c("stats", "loci", "coverages")
)
```

## Arguments

- x:

  The path to a fragments file, or a GRanges object containing the
  fragments (with the \`name\` column containing the barcode, and
  optionally the \`score\` column containing the count).

- barcodes:

  Optional character vector of cell barcodes to consider

- regionsToExclude:

  A GRanges of regions to exclude. As per the original Amulet method, we
  recommend excluding repeats, as well as sex and mitochondrial
  chromosomes. (Note that the end coordinate does not need to be exact
  when excluding entire chromosomes, but greater or equal to the
  chromosome length.)

- minFrags:

  Minimum number of fragments for a barcode to be considered. If
  \`uniqueFrags=TRUE\`, this is the minimum number of unique fragments.
  Ignored if \`barcodes\` is given.

- uniqueFrags:

  Logical; whether to use only unique fragments.

- maxFragSize:

  Integer indicating the maximum fragment size to consider

- removeHighOverlapSites:

  Logical; whether to remove sites that have more than two reads in
  unexpectedly many cells.

- fullInMemory:

  Logical; whether to process all chromosomes together. This will speed
  up the process but at the cost of a very high memory consumption (as
  all fragments will be loaded in memory). This is anyway the default
  mode when \`x\` is not Tabix-indexed.

- BPPARAM:

  A \`BiocParallel\` parameter object for multithreading. Note that
  multithreading will increase the memory usage.

- verbose:

  Logical; whether to print progress messages.

- ret:

  What to return, either barcode 'stats' (default), 'loci', or
  'coverages'.

## Value

A data.frame with counts and overlap statistics for each barcode.

## Details

When used on normal (or compressed) fragment files, this implementation
is relatively fast (except for reading in the data) but it has a large
memory footprint since the overlaps are performed in memory. It is
therefore recommended to compress the fragment files using bgzip and
index them with Tabix; in this case each chromosome will be read and
processed separately, leading to a considerably lower memory footprint.
