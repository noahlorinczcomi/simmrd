# Summary
The R scripts in this repository can be used to generate summary statistics to use in Mendelian Randomization simulations

There are many different statistical methods available to perform Mendelian Randomization (MR). When these methods are introduced in the literature, simulations are performed to evaluate their performance. There is currently no standard for performing these simulations, hence some methods perform differently in separate simulations intended to mirror the same reality (eg [ref](https://doi.org/10.1101/2021.03.26.437168) and [ref](https://doi.org/10.1214/20-AOS2027
), [ref](https://doi.org/10.1016/j.ajhg.2023.02.014) and [ref](https://doi.org/10.1002/gepi.22295), [ref](https://doi.org/10.1093/ije/dyaa262) and [ref](https://doi.org/10.1016/j.ajhg.2021.05.014)).

The code in this repo is open source, meaning **you** can modify it by initiating a [pull request](https://github.com/noahlorinczcomi/simmr/pulls). We encourage you to make changes to our code with the intention that the community will eventually settle on an accepted standard for how simulations in MR should be performed.

![](https://github.com/noahlorinczcomi/simmr/blob/main/simmr_flowchart.svg)

# Who is this for?
This software is primarily for researchers evaluating the performance of their new or existing MR estimator. This software can also be used by reviewers and others to validate the performance of pre-published manuscripts reporting MR simulations using GWAS summary statistics.

# What is in the software?
This software generates GWAS summary statistics under different scenarios of
- Exposure and outcome GWAS sample overlap [(ref)](https://doi.org/10.1101/2021.06.28.21259622)
- Uncorrelated horizontal pleiotropy [(ref)](https://doi.org/10.1093/ije/dyv080)
- Correlated horizontal pleiotropy [(ref)](https://doi.org/10.1038/s41588-020-0631-4)
- Weak instruments [(ref)](https://doi.org/10.1101/2023.01.10.523480)
- Winner's curse [(ref)](https://doi.org/10.1101/2021.06.28.21259622)
- Linkage disequilibrium between instruments [(ref)](https://doi.org/10.1002/gepi.22506)

that researchers may want to use in the evaluation of univariable or multivariable Mendelian Randomization methods.

# How do I use it?
The main task is downloading three files and making them communicate with each other, which is very easy. 

You can follow these steps:
1) Download the (generate_data.R)[https://github.com/noahlorinczcomi/simmr/blob/main/generate_data.R], [set_params.R](https://github.com/noahlorinczcomi/simmr/blob/main/generate_data.R), and [basic_functions.R](https://github.com/noahlorinczcomi/simmr/blob/main/generate_data.R) files into any folder
2) Start an R session and make this folder your working directory
3) Open the `set_params.R` file and change the parameter settings however you'd like
    - There are descriptions of what each parameter represents within the file itself
4) Run the command on line 35, which is `source('generate_data.R')`
5) All done! Your global R environment will now contain the 10 new objects
    - `bx`: standardized estimates of association between the IVs and the exposure(s)
    - `bxse`: standardized standard errors corresponding to `bx`
    - `by`: standardized estimates of association between the IVs and the outcome
    - `byse`: standardized standard errors correspondings to `by`
    - `mStart`: the number of causal variants you specified for the exposure(s)
    - `mIVs`: the final number of instruments returned
    - `LDMatrix`: the LD matrix for the final set of IVs
    - `RhoME`: a $(p+1)\times(p+1)$ matrix of correlations between errors in GWAS estimates for the outcome and $p$ exposure(s)
        - Some multivariable MR methods such as [MVMR-cML](https://doi.org/10.1016/j.ajhg.2023.02.014) and [MRBEE](https://doi.org/10.1101/2023.01.10.523480) use this matrix to correct for bias from weak instruments
    - `theta`: $p$-length vector of true causal effects
    - `IVtype`: Classification for each IV in data generation: UHP, CHP, Valid. Does not consider weakness!
    - `bxunstd`: unstandardized estimates of association between the final IVs and the exposures 
    - `bxseunstd`: unstandardized standard errors corresponding to `bxunstd`

There are at least three ways to generate the simulated data you want, which are described in sections [Method 1](#method-1:-everything-in-one-environment), [Method 2](#Method-2:-Sourcing-the-`set_params.R`-file), or [Method 3](#Method-3:-Sourcing-multiple-`set_params.R`-files-for-different-simulations).

# Can I add ____ to the simulation?
**YES**

We want this work to be a collaborative effort involving the broader community of MR researchers. Please, if you feel any simulation code contains errors or could be improved in any way, **please initiate a [pull request](https://github.com/noahlorinczcomi/simmr/pulls) and directly modify the code yourself**. 

If you have questions, please contact me: noahlorinczcomi@gmail.com or njl96@case.edu

# Tutorial
Assume we have downloaded the `set_params.R`, `generate_data.R`, and `basicfunctions.R` files into the `dir/` folder on your computer. 
## Method 1: Everything in one environment
Begin and R session and open the `set_params.R` file in (eg) Rstudio. You should see something like this:
```R
rm(list=ls(all=TRUE))
library(mvnfast)
################################################################################
### Exposure(s) (Xs) and Outcome (Y) GWAS
sample_size_Xs=2e4 # scalar
sample_size_Y=2e4 # scalar
prop_gwas_overlap_Xs_and_Y=0.1 # scalar
# prop_gwas_overlap_Xs=1 # fixed in current version
### Phenotypic and Genetic Correlations between Exposure(s) (Xs) and Outcome (Y)
number_of_exposures=3 # scalar
phenotypic_correlation_Xs='ar1(0.2)' # scalar or string (string examples: 'toeplitz','ar1(0.5)')
genetic_correlation_Xs='ar1(0.15)' # scalar or string (string examples: 'toeplitz','ar1(0.5)')
### Variances Explained in Exposure(s) (Xs), Confounder (U), and Outcome (U)
Xs_variance_explained_by_U=0.10 # scalar
Y_variance_explained_by_Xs=c(0.5,0.1,0) # scalar (applied to all exposures) or p-length vector
signs_of_causal_effects=c(1,-1,1) # scalar (applied to all exposures) or p-length vector of mixed 1's and -1's
Y_variance_explained_by_U=0.25 # scalar
### Set of SNPs Causal for Exposure(s) 
number_of_causal_SNPs=100 # scalar
Xs_variance_explained_by_g=0.15 # scalar
number_of_UHP_causal_SNPs=10 # scalar
number_of_CHP_causal_SNPs=10 # scalar
Y_variance_explained_by_UHP=0.05 # scalar
U_variance_explained_by_CHP=0.05 # scalar
mafs_of_causal_SNPs=0.3 # scalar
LD_causal_SNPs='ar(0.5)' # scalar or string (string examples: 'toeplitz','ar1(0.5)', or 'I')
### Standardizing MR data
MR_standardization_type='none' # Qi & Chatterjee MRMix paper, or could be 'z' (Z-score) or 'none'
### Performing IV selection
simtype='weak' # or winners
MVMR_IV_selection_type='joint' # either 'joint' (choose IVs significant in a p-degree of freedom joint test) or 'union' (union set of univariate association tests). Ignore if performing UVMR
IV_Pvalue_threshold=5e-5 # only SNPs with P<this threshold using your choice of MVMR_IV_selection_type test will be considered as IVs
LD_pruning_r2=0.1 # upper boundary of squared LD correlation
fix_Fstatistic_at=30 # average across exposures, not conditional F-statistics
source('generate_data.R')
################################################################################
### The global environment will now contain
# `bx`: standardized associations with the exposure(s)
# `bxse`: standardized standard errors corresponding to `bx`
# `by`: standardized associations with the outcome
# `byse`: standardized standard errors correspondings to `by`
# `bx`, `bxse`, `by` and `byse` are each now in the global environment
# `mStart`: the number of causal variants you specified
# 'mIVs`: the final number of instruments returned
# `LDMatrix`: the LD matrix for the final set of IVs
# `RhoME`: a (p+1)x(p+1) matrix of correlations between measurement errors
# `theta`: vector of true causal effects
# `IVtype`: Classification for each IV in data generation: UHP, CHP, Valid. Does not consider weakness!
# `bxunstd`: unstandardized estimates of association between the final IVs and the exposures 
# `bxseunstd`: unstandardized standard errors corresponding to `bxunstd`
################################################################################
```

Every time you run the `source(generate_data.R)` command, new simulation data will be generated. You can view a plot of your simulated data by (exactly) running `plot_simdata()`, which will produce output like this:
![](https://github.com/noahlorinczcomi/simmr/blob/main/example_plotsimdata.svg)

An example of how you use `simmr` to perform MVMR using IVW and MRBEE is this:
```R
# begin session
# set parameters:
# param1=<> ; param2=<> ; etc. etc. etc.
library(MRBEE)
n_simulations=1000
for(iteration in 1:n_simulations) {
    source('generate_data.R')
    weights=1/byse^2
    ivw=lm(by~bx-1,weights=weights)
    pD=prepData(list(R=RhoME,Ncor=1e5,EstHarm=cbind(by,bx),SEHarm=cbind(byse,bxse)))
    mrbee=MRBEE.IMRP(pD,FDR=TRUE)
    # save results etc.
}
```
## Method 2: Sourcing the `set_params.R` file
Or, you can save the `set_params.R` file with the parameters you want and perform simulations like this:
```R
# begin session
library(MRBEE)
n_simulations=1000
for(iteration in 1:n_simulations) {
    source('set_params.R')
    weights=1/byse^2
    ivw=lm(by~bx-1,weights=weights)
    pD=prepData(list(R=RhoME,Ncor=1e5,EstHarm=cbind(by,bx),SEHarm=cbind(byse,bxse)))
    mrbee=MRBEE.IMRP(pD,FDR=TRUE)
    # save results etc.
}
```
## Method 3: Sourcing multiple `set_params.R` files for different simulations

Or you can save multiple parameter files to produce different settings (e.g., save `set_params_scenario1.R`, `set_params_scenario2.R`) and source them in a nested `for()` loop in R like this:
```R
# begin session
library(MRBEE)
scenarios=c('scenario1','scenario2')
for(scenario in scenarios) { # 2 scenarios in this example
    for(iteration in 1:n_simulations) {
        fp=paste0('set_params_',scenario,'.R')
        source(fp)
        weights=1/byse^2
        ivw=lm(by~bx-1,weights=weights)
        pD=prepData(list(R=RhoME,Ncor=1e5,EstHarm=cbind(by,bx),SEHarm=cbind(byse,bxse)))
        mrbee=MRBEE.IMRP(pD,FDR=TRUE)
        # save results etc.
    }
}
```

