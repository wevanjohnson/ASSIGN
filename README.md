# ASSIGN: Adaptive Signature Selection and InteGratioN

__Author__: Ying Shen, Andrea H. Bild, David Jenkins, and W. Evan Johnson

__Maintainer__: [David Jenkins](https://github.com/dfjenkins3), [W. Evan Johnson](https://github.com/wevanjohnson/), [Mumtehena Rahman](https://github.com/mumtahena), and Ying Shen <yshen3@bu.edu>

ASSIGN is a computational tool to evaluate the pathway
deregulation/activation status in individual patient samples.
ASSIGN employs a flexible Bayesian factor analysis approach
that adapts predetermined pathway signatures derived either
from knowledge-based literatures or from perturbation
experiments to the cell-/tissue-specific pathway signatures.
The deregulation/activation level of each context-specific
pathway is quantified to a score, which represents the extent
to which a patient sample encompasses the pathway
deregulation/activation signature.

## Installation

ASSIGN is available on Bioconductor:

```
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("ASSIGN")
```

Or install the development version from github:

```
# install.packages("devtools")
devtools::install_github("wevanjohnson/ASSIGN")
```