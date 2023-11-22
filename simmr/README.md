# Installing `simmr`
In R, type
```R
remotes::install_github('noahlorinczcomi/simmr')
```
or
```R
devtools::install_github('noahlorinczcomi/simmr')
```
# Generating simulation data
`simmr` generates one simulated data set each time the `generate()` function is called. The `generate()` function takes one argument (`params`), which is a list of parameters used to define the simulation settings.
![](https://github.com/noahlorinczcomi/simmr/blob/main/simmr_flowchart.svg)

## Step 1: Define parameters
Create the object `params` like this:
```R
params=list(
    ### Exposure(s) (Xs) and Outcome (Y) GWAS
    sample_size_Xs=5e4, # scalar
    sample_size_Y=5e4, # scalar
    prop_gwas_overlap_Xs_and_Y=1, # scalar
    # prop_gwas_overlap_Xs=1 # fixed in current version
    ### Phenotypic and Genetic Correlations between Exposure(s) (Xs) and Outcome (Y)
    number_of_exposures=4, # scalar
    phenotypic_correlation_Xs=0.2, # scalar or string (string examples: 'toeplitz','ar1(0.5)')
    genetic_correlation_Xs=0, # scalar or string (string examples: 'toeplitz','ar1(0.5)')
    ### Variances Explained in Exposure(s) (Xs), Confounder (U), and Outcome (U)
    Xs_variance_explained_by_U=1/4-0.12, # scalar
    Y_variance_explained_by_Xs=c(0,0.2,0.2,0.4), # scalar (applied to all exposures) or p-length vector
    signs_of_causal_effects=c(1,1,1,1), # scalar (applied to all exposures) or p-length vector of mixed 1's and -1's
    Y_variance_explained_by_U=0.3, # scalar
    ### Set of SNPs Causal for Exposure(s) 
    number_of_causal_SNPs=100, # scalar
    Xs_variance_explained_by_g=0.12, # scalar
    number_of_UHP_causal_SNPs=30, # scalar
    number_of_CHP_causal_SNPs=10, # scalar
    Y_variance_explained_by_UHP=0.05, # scalar
    U_variance_explained_by_CHP=0.05, # scalar
    mafs_of_causal_SNPs=0.3, # scalar
    LD_causal_SNPs='I', # scalar or string (string examples: 'toeplitz','ar1(0.5)', or 'I')
    ### Standardizing MR data
    MR_standardization='none', # could be 'z' (Z-score) or 'none'. Note that all genotypes and phenotypes are already standardized
    ### Performing IV selection
    simtype='winners', # 'weak' or 'winners'
    MVMR_IV_selection_type='joint', # either 'joint' (choose IVs significant in a p-degree of freedom joint test) or 'union' (union set of univariate association tests). Ignore if performing UVMR
    IV_Pvalue_threshold=5e-5, # only SNPs with P<this threshold using your choice of MVMR_IV_selection_type test will be considered as IVs
    LD_pruning_r2=1, # upper boundary of squared LD correlation
    N_of_LD_ref=Inf, # LD matrix in MR is same as true LD matrix. Otherwise put a size of the ref panel (eg `20000` to represent 20k individuals)
    fix_Fstatistic_at=10 # average across exposures, not conditional F-statistics
)
```

## Step 2: Generate data
Now give `params` to `generate()` like this:
```R
library(simmr)
data=generate(params)
```

## Step 3: Use data
The object `data` is a list containing the following values
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

## Step 4: Plot simulated data
You can use the `plot_simdata(data)` function to view the data you simulated, which can produce the following output.
![](https://github.com/noahlorinczcomi/simmr/blob/main/example_plotsimdata.svg)  