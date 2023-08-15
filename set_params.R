rm(list=ls(all=TRUE))
library(mvnfast)
sample_size_Y=1e4 # scalar
sample_size_Xs=1e4 # scalar
number_of_exposures=3 # scalar
# proportion_overlapping_between_Xs_GWAS=1 # (fixed)
proportion_overlapping_in_Xs_and_Y_GWAS=0.1 # scalar or vector
# variance_U=1 # (fixed) scalar
phenotypic_correlation_Xs='ar1(0.2)' # scalar or string (string examples: 'toeplitz','ar1(0.5)')
genetic_correlation_Xs='ar1(0.15)' # scalar or string (string examples: 'toeplitz','ar1(0.5)')
variance_in_Xs_explained_by_U=0.1 # scalar or vector
variance_in_Y_explained_by_U=0.1 # scalar
variance_in_Y_explained_by_Xs=0.5 # scalar or vector
number_of_instruments=500 # scalar
variance_in_Xs_explained_by_all_instruments=0.15 # scalar or vector
number_of_weak_instruments=10 # scalar (no conditionally weak IVs)
variance_in_Xs_explained_by_weak_instruments=0.01 # scalar or vector
number_of_UHP_instruments=100 # scalar
number_of_CHP_instruments=100 # scalar
variance_in_Y_explained_by_UHP_instruments=0.05 # scalar
variance_in_Y_explained_by_CHP_instruments=0.15 # scalar
variance_in_U_explained_by_CHP_instruments=0.15 # scalar
scale_of_CHP_effects=-1 # scalar
mafs_of_instruments=0.3 # scalar
LD_instruments='I' # scalar or string (string examples: 'toeplitz','ar1(0.5)', or 'I')
MR_standardization_type='qi' # Qi & Chatterjee (MRMix paper), or could be 'z'
outcome_type='binary' # or anything else
exposure_types='binary' # or anything else, eg 'wQ#4tB @# TQ' will be interpreted as 'continuous'
source('generate_data.R')
# `bx`, `bxse`, `by` and `byse` are each now in the global environment
