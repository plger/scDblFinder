# plDoubletFinder

This is an attempt at recreating the rough logic of [DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder), which creates random artificial doublets to detect real doublets, but in a more efficient manner (i.e. simpler, faster, and apparently more accurate).

For more information and a comparison with other tools, see [doublets.md](doublets.md).

## Installation

```{r}
devtools::install_github('plger/plDoubletFinder')
```

### Usage

Given an object `sce` of class `SingleCellExperiment`:
```{r}
library(plDoubletFinder)
res <- plDoubletFinder(sce)
```

`res$ratio` contains the scores, and `res$classification` the final calls.
