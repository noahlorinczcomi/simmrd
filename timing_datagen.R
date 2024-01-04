rm(list=ls(all=TRUE))
library(mvnfast)
ns=seq(1e4,5e4,1e4)
ms=c(10,50,100,250,500)
res=matrix(nr=length(ms),nc=length(ns))
for(ii in 1:length(ms)) {
  for(kk in 1:length(ns)) {
    sample_size_Xs=ns[kk] # scalar
    sample_size_Y=ns[kk] # scalar
    prop_gwas_overlap_Xs_and_Y=1 # scalar
    # prop_gwas_overlap_Xs=1 # fixed in current version
    ### Phenotypic and Genetic Correlations between Exposure(s) (Xs) and Outcome (Y)
    number_of_exposures=3 # scalar
    phenotypic_correlation_Xs='ar1(0.2)' # scalar or string (string examples: 'toeplitz','ar1(0.5)')
    genetic_correlation_Xs='ar1(0.15)' # scalar or string (string examples: 'toeplitz','ar1(0.5)')
    ### Variances Explained in Exposure(s) (Xs), Confounder (U), and Outcome (U)
    Xs_variance_explained_by_U=0.10 # scalar
    Y_variance_explained_by_Xs=0.50 # scalar
    Y_variance_explained_by_U=0.25 # scalar
    ### Set of SNPs Causal for Exposure(s) 
    number_of_causal_SNPs=ms[ii] # scalar
    Xs_variance_explained_by_g=0.15 # scalar
    number_of_UHP_causal_SNPs=0 # scalar
    number_of_CHP_causal_SNPs=0 # scalar
    Y_variance_explained_by_UHP=0.05 # scalar
    U_variance_explained_by_CHP=0.05 # scalar
    mafs_of_causal_SNPs=0.3 # scalar
    LD_causal_SNPs='ar(0.5)' # scalar or string (string examples: 'toeplitz','ar1(0.5)', or 'I')
    ### Standardizing MR data
    MR_standardization_type='z' # Qi & Chatterjee MRMix paper, or could be 'z' (Z-score) or 'none'
    outcome_type='binary' # or anything else, eg 'wQ#4tB @# TQ' will be interpreted as 'continuous'
    exposure_types='binary' # or anything else
    ### Performing IV selection
    simtype='winners' # or winners
    IV_Pvalue_threshold=5e-5 # in a joint test of H0: beta_j1=betaj2=...=betaj3=0 when there are multiple exposures
    LD_pruning_r2=0.1 # upper boundary of squared LD correlation
    fix_Fstatistic_at=30 # average across exposures, not conditional F-statistics
    to=Sys.time()
    source('generate_data.R')
    res[ii,kk]=as.numeric(Sys.time()-to)
    print(res)
  }
}



