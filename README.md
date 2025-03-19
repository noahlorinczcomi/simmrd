```simmrd``` is an R package that can be used to generate simulated data for univariable or multivariable Mendelian Randomization (MR) under various conditions of weak instrument bias, GWAS sample overlap, uncorrelated horizontal pleiotropy, and correlated horizontal pleiotropy. These data are intended to be used to evaluate the statistical properties of existing or candidiate MR methods.

# Installing `simmrd`
In R, type
```R
remotes::install_github('noahlorinczcomi/simmrd')
```
or
```R
devtools::install_github('noahlorinczcomi/simmrd')
```
# Paper
Lorincz-Comi, N., Yang, Y., & Zhu, X. simmrd: An open-source tool to perform simulations in Mendelian randomization. Genetic epidemiology, _1-15_. [doi.org/10.1002/gepi.22544](https://doi.org/10.1002/gepi.22544)
# Tutorial
Please follow our [Tutorial](https://github.com/noahlorinczcomi/simmrd/wiki/Tutorial) to learn how to generate GWAS summary data for use in MR simulations.
# Output
It is shown in the [Tutorial](https://github.com/noahlorinczcomi/simmrd/wiki/Tutorial) that the output of either `generate_summary()` or `generate_individual()` is a named list with the following elements:
1. `bx`: $m\times p$ matrix of standardized GWAS-estimated associations between the $m$ IVs and $p$ exposures 
   - The type of standardization is specified by the user with the `MR_standardization` item in their list of parameters (see the [Tutorial](https://github.com/noahlorinczcomi/simmrd/wiki/Tutorial)
2. `by`: $m\times 1$ vector of standardized GWAS-estimated associations between the $m$ IVs and outcome
3. `bxse`: $m\times p$ matrix of estimated standard errors corresponding to `bx`
4. `byse`: $m\times 1$ vector of estimated standard errors corresponding to `by`
5. `RhoME`: $(p+1)\times (p+1)$ matrix of correlations between GWAS estimation errors for the outcome and exposures
   - Consider $\hat{\boldsymbol\beta}_j$ as the $p$-length vector of estimated associations between the *j*th IV and the $p$ exposures, $\boldsymbol\beta_j$ as its estimand, $\hat\alpha_j$ as the estimated association between the *j*th IV and the outcome, and $\alpha_j$ as its estimand. `RhoME` stores $\text{Corr}(\hat{\boldsymbol\beta}_j-\boldsymbol\beta_j,\hat\alpha_j-\alpha_j)$, which is used by MR methods such as [MRBEE](https://doi.org/10.1101/2023.01.10.523480) and [MR-cML](https://doi.org/10.1016/j.ajhg.2023.02.014) to correct for measurement error bias, which includes bias from sample overlap.
6. `LDMatrix`: $m\times m$ matrix of true LD correlations between the IVs
7. `LDhatMatrix`: $m\times m$ matrix of estimated LD correlations between the IVs
   - estimates of the true LD matrix are made by randomly sampling from an *n*-degree of freedom Wishart distribution, where *n* is the user-specified LD reference panel size
8. `theta`: $p\times 1$ vector of true causal effects of the exposures
9. `IVtype`: $m\times 1$ vector of classifications of 'valid', 'UHP', and 'CHP' for each IV to indicate the presence or absence of uncorrelated or correlated horizontal pleiotropy.
10. `bx_unstd`: $m\times p$ matrix of unstandardized GWAS-estimated associations between the $m$ IVs and $p$ exposures
11. `by_unstd`: $m\times 1$ vector of unstandardized GWAS-estimated associations between the $m$ IVs and outcome
12. `bxse_unstd`: $m\times p$ matrix of estimated standard errors corresponding to `bx_unstd`
13. `byse_unstd`: $m\times 1$ vector of estimated standard errors corresponding to `by_unstd`
<!--
15. `true_variance_explained`: scalar value representing the true outcome variance explained by the causal effects of the exposures
    - When the user elects to use `generate_summary()` to generate their data, this value is calculated as `true_variance_explained`$:=\rho^2_{\alpha|\beta}=\boldsymbol\theta^\top(\mathbf{B}^\top\mathbf{R}\mathbf{B}+\boldsymbol\Sigma_\varepsilon)\boldsymbol\theta$, where $\mathbf{x}=\mathbf{B}^\top\mathbf{g}+\boldsymbol\varepsilon$, $\mathbf{g}$ is a vector of causal exposure SNP genotypes, $\mathbf{R}=\text{Cov}(\mathbf{g})$, and $\boldsymbol\Sigma_\varepsilon=\text{Cov}(\boldsymbol\epsilon
-->
