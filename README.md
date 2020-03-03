# RNDClone
Software package `RNDClone` for the paper "RNDClone: Tumor Subclone Reconstruction Based on Integrating DNA and RNA Sequence Data"

## Installation
The `RNDClone` package does not have any dependencies.

The package can be easily installed with the `devtools` package in R. After `devtools` has been installed, run the following commands in R to install the `RNDClone` package.
```
library(devtools)
install_github("tianjianzhou/RNDClone")
```


## Example
Below is a simple example to replicate the simulation results in the paper.

```
library(RNDClone)

data(sim1a_C4_T4)

# Retrieve data
n = sim1a_C4_T4$n
N = sim1a_C4_T4$N
m = sim1a_C4_T4$m
M = sim1a_C4_T4$M
g_fun = sim1a_C4_T4$g_fun

set.seed(345)

# Run the trans-dimensional MCMC as described in the paper (may take a while, ~ 1 hr)
MCMC_spls = RNDClone_RJMCMC(n = n, N = N, m = m, M = M, g_fun = g_fun)

# For testing purpose, use (small number of iterations and burnin)
# MCMC_spls = RNDClone_RJMCMC(n = n, N = N, m = m, M = M, g_fun = g_fun, niter = 50, burnin = 200, thin = 2)

# Retrieve posterior samples of the parameters
C_spls = MCMC_spls$sample_list$C_spls
L_spls = MCMC_spls$sample_list$L_spls
Z_spls = MCMC_spls$sample_list$Z_spls
Lambda_spls = MCMC_spls$sample_list$Lambda_spls
W_spls = MCMC_spls$sample_list$W_spls

# Point estimate of C: posterior mode
C_hat = which.max(tabulate(C_spls))

# Point estimates of L, Z, W and Lambda: Maximum A Posteriori (MAP) conditional on C_hat
# First find which sample has the largest log-posterior
logpost_spls = MCMC_spls$sample_list$logpost_spls
logpost_spls[C_spls != C_hat] = -Inf
index_MAP = which.max(logpost_spls)
L_hat = L_spls[[index_MAP]]
Z_hat = Z_spls[[index_MAP]]
Lambda_hat = Lambda_spls[[index_MAP]]
# The last column of W_hat corresponds to w[t0] in the paper, which is used to capture random noise
W_hat = W_spls[[index_MAP]]
```


## Details
The `RNDClone` package contains four functions: `RNDClone_RJMCMC`, `DClone_RJMCMC`, `RClone_RJMCMC`, and `RNDClone_PT`.

To understand how to use these functions, it is necessary to introduce some notation: `n`, `N`, `m`, `M`, and `g_fun`. Let `S` denote the number of loci of the nucleotides that are covered by short reads produced by next-generation sequencing experiments. Let `T` denote the number of tissue samples. Let `G` denote the number of genes in which the `S` nucleotide loci reside.

- `n` is a `S * T` matrix, where `n[s, t]` is the number of variant DNA reads at locus `s` for sample `t`
- `N` is a `S * T` matrix , where `N[s, t]` is the total number of DNA reads at locus `s` for sample `t`
- `m` is a `S * T` matrix, where `m[s, t]` is the number of variant RNA reads at locus `s` for sample `t`
- `M` is a `S * T` matrix , where `M[s, t]` is the total number of RNA reads at locus `s` for sample `t`
- `g_fun` is a length `S` vector, where `g_fun[s]` is the index of the gene that locus `s` reside in (This argument is optional for the four functions)


`RNDClone_RJMCMC` is a function of the four data matrices `n`, `N`, `m` and `M`. It has the following form:
```
RNDClone_RJMCMC(n = n, N = N, m = m, M = M, g_fun = NULL, ...)
```
Hyperparameters, such as `C_min`, `C_max`, `K_min`, `K_max`, `a_w`, `b_w`, `d`, ..., may be specified, although they have default values. MCMC parameters, such as `niter`, `burnin` and `thin`, may also be specified. The default is `niter = 5000`, `burnin = 20000` and `thin = 2`.

