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