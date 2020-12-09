# DEvis

DEvis is a powerful, integrated solution for the analysis of differential expression data. 
This package includes an array of tools for manipulating and aggregating data, as well as a wide range
of customizable visualizations, and project management functionality that simplify RNA-Seq analysis 
and provide a variety of ways of exploring and analyzing data.  

## Publication

The DEvis paper can be found [here](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2702-z).  


Price, A., Caciula, A., Guo, C., Lee, B., Morrison, J., Rasmussen, A., Lipkin, W.I., Jain, K. (2019). DEvis: an R package for aggregation and visualization of differential expression data. BMC Bioinformatics, 20(1), 110. https://doi.org/10.1186/s12859-019-2702-z


## Installation:

### DEvis is now available through the CRAN repository.

install.packages("DEVis")

### Manual Installation:
#### If using OS X, please ensure that you are using the most up to date version of xcode. 
Instructions for [how to update xcode](https://stackoverflow.com/questions/15417619/how-do-you-update-xcode-on-osx-to-the-latest-version).
```
#Install DESeq2 dependency from bioconductor.
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")

#Install devtools to allow installation from GitHub
if (!require("devtools")) install.packages("devtools")

#Install DEvis from GitHub repository.
devtools::install_github("price0416/DEvis/DEvis")

#Load the package.
library(DEVis)
```


## Documentation:

For a complete usage tutorial, [see the vignette](https://cran.r-project.org/web/packages/DEVis/vignettes/DEVis_vignette.html).

For technical specifications, [see the DEvis manual](https://github.com/price0416/DEvis/blob/master/DEVis/vignettes/DEVis_vignette.pdf).


