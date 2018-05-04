# DEvis

## Installation:
### If using OS X, please ensure that you are using the most up to date version of xcode. 
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

For a complete usage tutorial, [see the vignette](https://github.com/price0416/DEvis/blob/master/DEVis/vignettes/DEVis_vignette.pdf).


For technical specifications, [see the DEvis manual](https://github.com/price0416/DEvis/blob/master/DEVis/man/DEVis.pdf).


# Note to Reviewers:
Thank you for taking the time to review our package.  Once we receive notification of acceptance we will submit this package to the CRAN repository.  Currently the package meets all standards for submission and acceptance to CRAN.

Once accepted to CRAN, installation will be simplified to:
install.packages(DEVis)

