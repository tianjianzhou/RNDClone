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

To understand how to use these functions, it is necessary to introduce some notation: `n`, `N`, `m`, `M`, and `g_fun`. Let `S` denote the number of loci of the nucleotides that are covered by short reads produced by next-generation sequencing experiments. Let `T` denote the number of tissue samples. Let `G` denote the number of genes in which the `S` nucleotide loci reside.
- `n` is a `S * T` matrix, where `n[s, t]` is the number of variant DNA reads at locus `s` for sample `t`
- `N` is a `S * T` matrix , where `N[s, t]` is the total number of DNA reads at locus `s` for sample `t`
- `m` is a `S * T` matrix, where `m[s, t]` is the number of variant RNA reads at locus `s` for sample `t`
- `M` is a `S * T` matrix , where `M[s, t]` is the total number of RNA reads at locus `s` for sample `t`
- (Optional) `g_fun` is a length `S` vector, where `g_fun[s]` is the index of the gene that locus `s` reside in
