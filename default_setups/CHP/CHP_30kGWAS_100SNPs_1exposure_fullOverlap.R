################################################################################
# using summary data
################################################################################
summary_params=list(
  ### parameters distinguishing this setup from others
  sample_size_Xs=30000, # exposure GWAS n's (scalar or vector)
  sample_size_Y=30000, # outcome GWAS n
  prop_gwas_overlap_Xs_and_Y=1, # proportion of overlap between exposures and outcome GWAS (scalar or vector)
  number_of_exposures=3, # number of exposures
  number_of_causal_SNPs=100, # number of SNPs causing each exposure
  number_of_UHP_causal_SNPs=0, # number of UHP causal SNPs
  number_of_CHP_causal_SNPs=20, # number of CHP causal SNPs
  ratio_of_UHP_variance=0.15, # ratio of UHP variance to valid IV variance
  ratio_of_CHP_variance=0.25, # ratio of CHP variance to valid IV variance
  CHP_correlation=-0.5, # (magnitude of CHP) correlation between CHP and valid IV effect sizes
  simtype='winners', # performs IV selection based on P-value
  fix_Fstatistic_at=10, # (ignore if simtype=='winners') fix the mean F-statistic for the exposures
  ### parameters fixed across all default setups
  prop_gwas_overlap_Xs=1, # overlap of exposures GWAS (scalar or numeric matrix)
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


params=list(
  sample_size_Xs=c(2e4,5e4,1e5), # exposure GWAS sample sizes
  sample_size_Y=3e4, # outcome GWAS sample size
  prop_gwas_overlap_Xs_and_Y=1/2, # proportion of exposures' and outcome GWAS overlap
  number_of_exposures=3, # number of exposures
  number_of_causal_SNPs=200, # number of SNPs causing each exposure
  number_of_UHP_causal_SNPs=30, # number of UHP causal SNPs
  number_of_CHP_causal_SNPs=60, # number of CHP causal SNPs
  ratio_of_UHP_variance=2, # ratio of UHP variance to valid IV variance
  ratio_of_CHP_variance=1.5, # ratio of CHP variance to valid IV variance
  CHP_correlation=-0.5, # correlation between CHP and valid IV effect sizes
  simtype='winners', # performs IV selection based on P-value
  fix_Fstatistic_at=10, # ignored because simtype='winners'
  prop_gwas_overlap_Xs=1, # overlap of exposures' GWAS
  phenotypic_correlation_Xs=0.3, # phenotypic correlations between exposures
  genetic_correlation_Xs=0.15,  # genetic correlation between exposures
  phenotypic_correlations_Xs_and_Y=0.3, # phenotypic correlations b/w exposures and outcome
  true_causal_effects=c(0,0.1,0.3), # true causal effect sizes
  Xs_variance_explained_by_g=0.10, # exposure variance explained by SNPs
  LD_causal_SNPs='ar1(0.5)', # LD between causal exposure SNPs
  number_of_LD_blocks=3, # number of independent LD blocks
  MR_standardization='none', # does not standardize GWAS estimates
  MVMR_IV_selection_type='union', # SNPs associated with >0 exposures are candidate IVs
  IV_Pvalue_threshold=5e-8, # only SNPs with P<this threshold are candidate IVs
  LD_pruning_r2=0.3, # the upper LD r2 pruning threshold for candidate IVs
  N_of_LD_ref=500 # the sample size of the LD reference panel
)
t0=Sys.time()
data=simmr::generate_summary(params)
Sys.time()-t0


# IVW and MR-Egger example
set.seed(74667)
library(MendelianRandomization)
library(simmr)
data=generate_summary(params)
obj=mr_mvinput(
  bx=data$bx, by=data$by,
  bxse=data$bxse, byse=data$byse,
  correlation=data$LDhatMatrix)
plot_simdata(data,params)
mr_mvivw(obj)
mr_mvmedian(obj)
RhoME=data$RhoME[c(2:4,1),c(2:4,1)]
mr_mvcML(obj,n=2e4,rho_mat=RhoME)



















