rm(list=ls(all=TRUE))
library(mvnfast)
ns=seq(1e4,5e4,1e4)
ms=c(10,50,100,250,500)
res=matrix(nr=length(ms),nc=length(ns))
for(i in 1:length(ms)) {
  for(k in 1:length(ns)) {
    sample_size_Xs=sample_size_Y=ns[k]
    number_of_causal_SNPs=ms[i]
    proportion_overlapping_in_Xs_and_Y_GWAS=0.1 # scalar
    ### Phenotypic and Genetic Correlations between Exposure(s) (Xs) and Outcome (Y)
    number_of_exposures=3 # scalar
    phenotypic_correlation_Xs='ar1(0.2)' # scalar or string (string examples: 'toeplitz','ar1(0.5)')
    genetic_correlation_Xs='ar1(0.15)' # scalar or string (string examples: 'toeplitz','ar1(0.5)')
    ### Variances Explained in Exposure(s) (Xs), Confounder (U), and Outcome (U)
    variance_in_Xs_explained_by_U=0.1 # scalar or vector
    variance_in_Y_explained_by_U=0.25 # scalar
    variance_in_Y_explained_by_Xs=0.5 # scalar or vector
    ### Set of SNPs Causal for Exposure(s) 
    variance_in_Xs_explained_by_all_causal_SNPs=0.15 # scalar or vector
    number_of_weak_causal_SNPs=0 # scalar (weak for all exposures if there's more than 1)
    variance_in_Xs_explained_by_weak_causal_SNPs=0.001 # scalar or vector
    number_of_UHP_causal_SNPs=0 # scalar
    number_of_CHP_causal_SNPs=0 # scalar
    variance_in_Y_explained_by_UHP_causal_SNPs=0.05 # scalar
    variance_in_Y_explained_by_CHP_causal_SNPs=0.15 # scalar
    variance_in_U_explained_by_CHP_causal_SNPs=0.15 # scalar
    mafs_of_causal_SNPs=0.3 # scalar
    LD_causal_SNPs='ar(0.5)' # scalar or string (string examples: 'toeplitz','ar1(0.5)', or 'I')
    ### Standardizing MR data
    MR_standardization_type='z' # Qi & Chatterjee MRMix paper, or could be 'z' (Z-score) or 'none'
    outcome_type='binary' # or anything else, eg 'wQ#4tB @# TQ' will be interpreted as 'continuous'
    exposure_types='binary' # or anything else
    ### Performing IV selection
    instrument_selection_Pvalue_threshold=5e-5 # in a joint test of H0: beta_j1=betaj2=...=betaj3=0 when there are multiple exposures
    instrument_selection_LD_pruning_r2=0.1 # upper boundary of squared LD correlation
    to=Sys.time()
    source('generate_data.R')
    res[i,k]=as.numeric(Sys.time()-to)
    print(res)
  }
}



