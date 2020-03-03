# RNDClone
Software package `RNDClone` for the paper "RNDClone: Tumor Subclone Reconstruction Based on Integrating DNA and RNA Sequence Data"

## Installation
The `RNDClone` package does not have any dependencies.

The package can be easily installed with the `devtools` package in R. After `devtools` has been installed, run the following commands in R to install the `RNDClone` package.
```
library(devtools)
install_github("tianjianzhou/RNDClone")
```

## Usage
The `RNDClone` package contains four functions: `RNDClone_RJMCMC`, `DClone_RJMCMC`, `RClone_RJMCMC`, and `RNDClone_PT`.

To understand how to use these functions, it is necessary to introduce some notation: `n`, `N`, `m`, `M`, and `g_fun`.
- `n` is a S * T matrix
- `N` is a S * T
