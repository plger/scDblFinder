# clamulet

Classification-powered Amulet-like method

## Usage

``` r
clamulet(
  x,
  artificialDoublets = NULL,
  iter = 2,
  k = NULL,
  minCount = 0.001,
  maxN = 500,
  nfeatures = 25,
  max_depth = 5,
  threshold = 0.75,
  returnAll = FALSE,
  verbose = TRUE,
  ...
)
```

## Arguments

- x:

  The path to a fragment file (see
  [`getFragmentOverlaps`](https://plger.github.io/scDblFinder/reference/getFragmentOverlaps.md)
  for performance/memory-related guidelines)

- artificialDoublets:

  The number of artificial doublets to generate

- iter:

  The number of learning iterations (should be 1 to)

- k:

  The number(s) of nearest neighbors at which to gather statistics

- minCount:

  The minimum number of cells in which a locus is detected to be
  considered. If lower than 1, it is interpreted as a fraction of the
  number of cells.

- maxN:

  The maximum number of regions per cell to consider to establish
  windows for meta-features

- nfeatures:

  The number of meta-features to consider

- max_depth:

  The maximum tree depth

- threshold:

  The score threshold used during iterations

- returnAll:

  Logical; whether to return data also for artificial doublets

- verbose:

  Logical; whether to print progress information

- ...:

  Arguments passed to
  [`getFragmentOverlaps`](https://plger.github.io/scDblFinder/reference/getFragmentOverlaps.md)

## Value

A data.frame

## Details

\`clamulet\` operates similarly to the \`scDblFinder\` method, but
generates doublets by operating on the fragment coverages. This has the
advantage that the number of loci covered by more than two reads can be
computed for artificial doublets, enabling the use of this feature
(along with the kNN-based ones) in a classification scheme. It however
has the disadvantage of being rather slow and memory hungry, and appears
to be outperformed by a simple p-value combination of the two methods
(see vignette).
