# amulet

ATACseq (Thibodeau, Eroglu, et al., Genome Biology 2021). The rationale
is that cells with unexpectedly many loci covered by more than two reads
are more likely to be doublets.

## Usage

``` r
amulet(x, ...)
```

## Arguments

- x:

  The path to a fragments file, or a GRanges object containing the
  fragments (with the \`name\` column containing the barcode, and the
  \`score\` column containing the count).

- ...:

  Any argument to
  [`getFragmentOverlaps`](https://plger.github.io/scDblFinder/reference/getFragmentOverlaps.md).

## Value

A data.frame including, for each barcode, the number sites covered by
more than two reads, the number of reads, and p- and q-values (low
values indicative of doublets).

## Details

When used on normal (or compressed) fragment files, this implementation
is relatively fast (except for reading in the data) but it has a large
memory footprint since the overlaps are performed in memory. It is
therefore recommended to compress the fragment files using bgzip and
index them with Tabix; in this case each chromosome will be read and
processed separately, leading to a considerably lower memory footprint.
See the underlying
[`getFragmentOverlaps`](https://plger.github.io/scDblFinder/reference/getFragmentOverlaps.md)
for details.

## Examples

``` r
# here we use a dummy fragment file for example:
fragfile <- system.file( "extdata", "example_fragments.tsv.gz",
                         package="scDblFinder" )
res <- amulet(fragfile)
#> Fragment file is not tabix-indexed, requiring thewhole file to be imported in memory.
#> 08:21:03 PM - Splitting and subsetting barcodes...
#> 08:21:03 PM - Obtaining overlaps...
```
