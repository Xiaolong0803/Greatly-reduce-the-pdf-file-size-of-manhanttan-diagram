# R_ggplot
Some tips about R drawing

## Millions of SNP data to draw PDF version of Manhattan map (file size 200k,dpi=300)

### First list the advantages of this article:

- More than 21 million SNPS: the original pdf size is 773M, and the pdf size is 203k after the point layer raster.

- 7686395 SNPS: original pdf size 275M, point layer raster pdf size 200k.

- qqman sample data 16470 SNP, original pdf size 853k, point layer raster pdf size 254k

library(tidyverse)
library(ggrastr)
library(qqman)

head(gwasResults,3)
