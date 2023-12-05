################################################################################
# using summary data
################################################################################
summary_params=list(
  ### parameters distinguishing this setup from others
  sample_size_Xs=100000, # exposure GWAS n's (scalar or vector)
  sample_size_Y=100000, # outcome GWAS n
  prop_gwas_overlap_Xs_and_Y=0, # proportion of overlap between exposures and outcome GWAS (scalar or vector)
  number_of_exposures=1, # number of exposures
  number_of_causal_SNPs=500, # number of SNPs causing each exposure
  number_of_UHP_causal_SNPs=50, # number of UHP causal SNPs
  number_of_CHP_causal_SNPs=50, # number of CHP causal SNPs
  ratio_of_UHP_variance=0.125, # ratio of UHP variance to valid IV variance
  ratio_of_CHP_variance=0.125, # ratio of CHP variance to valid IV variance
  CHP_correlation=0, # (magnitude of CHP) correlation between CHP and valid IV effect sizes
  simtype='winners', # performs IV selection based on P-value
  fix_Fstatistic_at=10, # (ignore if simtype=='winners') fix the mean F-statistic for the exposures
  ### parameters fixed across all default setups
  prop_gwas_overlap_Xs=0.3, # overlap of exposures GWAS (scalar or numeric matrix)
  phenotypic_correlation_Xs=0.3, # phenotypic correlations between exposures (scalar, string, or matrix; string examples are 'ar1(0.5), 'I', 'toeplitz')
  genetic_correlation_Xs=0.15,  # genetic correlation between exposures (scalar, string, or matrix; string examples are 'ar1(0.5), 'I', 'toeplitz')
  phenotypic_correlations_Xs_and_Y=0.3, # phenotypic correlation between exposures and outcome (scalar or vector)
  true_causal_effects=0.3, # true causal effect sizes (scalar or vector)
  Xs_variance_explained_by_g=0.10, # variance in exposures explained by all causal SNPs (scalar or vector)
  LD_causal_SNPs='I', # scalar, string (examples: 'toeplitz','ar1(0.5)'), or numeric matrix
  number_of_LD_blocks=1, # number of independent LD blocks
  MR_standardization='none', # could be 'z' (Z-score) or 'none'. Note that all genotypes and phenotypes are already assumed to be standardized to mean mean 0 and variance 1
  MVMR_IV_selection_type='union', # (irrelevant for UVMR) 'union': union set of exposure-specific IV sets; 'joint': SNPs significant in a p-degree of freedom joint chi-square test of association with the exposures
  IV_Pvalue_threshold=5e-8, # only SNPs with P<this threshold using your choice of MVMR_IV_selection_type test will be considered as IVs
  LD_pruning_r2=1, # the upper LD r2 pruning threshold for SNPs to be considered as IVs
  N_of_LD_ref=Inf # the sample size of the LD reference panel for estimating LD between IVs
)