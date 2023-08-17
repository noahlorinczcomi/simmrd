### example plot of LD between SNPs
# y=runif(500,-1,1)
# y=abs(y)
# y[50:150]=rnorm(101)
# y[190:230]=rnorm(41,0,2/2.5)
# y=y^2
# y=y[1:300]
# plot(y,pch=23,bg='gray80',cex=1.5)
# leadix=which.max(y)
# points(leadix,y[leadix],pch=23,bg='indianred',cex=1.5)
# o=cbind(1:length(y),y); o=o[rev(order(o[,2])),]
# o=o[o[,1]<180,]
# fivemore=o[,1][2:6]
# tenmore=o[,1][7:16]
# fiveagain=o[,1][17:22]
# points(fivemore,y[fivemore],pch=23,bg='orange',cex=1.5)
# points(tenmore,y[tenmore],pch=23,bg='cornflowerblue',cex=1.5)
points(fiveagain,y[fiveagain],pch=23,bg='darkgreen',cex=1.5)

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
variance_in_Xs_explained_by_U=0.1 # scalar or vector
variance_in_Y_explained_by_U=0.25 # scalar
variance_in_Y_explained_by_Xs=0.5 # scalar or vector
### Set of SNPs Causal for Exposure(s) 
number_of_causal_SNPs=100 # scalar
variance_in_Xs_explained_by_all_causal_SNPs=0.15 # scalar or vector
number_of_weak_causal_SNPs=5 # scalar (weak for all exposures if there's more than 1)
variance_in_Xs_explained_by_weak_causal_SNPs=0.001 # scalar or vector
number_of_UHP_causal_SNPs=10 # scalar
number_of_CHP_causal_SNPs=10 # scalar
variance_in_Y_explained_by_UHP_causal_SNPs=0.05 # scalar
variance_in_U_explained_by_CHP_causal_SNPs=0.10 # scalar
mafs_of_causal_SNPs=0.3 # scalar
LD_causal_SNPs='ar(0.5)' # scalar or string (string examples: 'toeplitz','ar1(0.5)', or 'I')
### Standardizing MR data
MR_standardization_type='z' # Qi & Chatterjee MRMix paper, or could be 'z' (Z-score) or 'none'
outcome_type='binary' # or anything else, eg 'wQ#4tB @# TQ' will be interpreted as 'continuous'
exposure_types='binary' # or anything else
### Performing IV selection
instrument_selection_Pvalue_threshold=5e-5 # in a joint test of H0: beta_j1=betaj2=...=betaj3=0 when there are multiple exposures
instrument_selection_LD_pruning_r2=0.1 # upper boundary of squared LD correlation
source('generate_data.R')
################################################################################
### plot of UHP,CHP,Valid,Weak IVs
library(dplyr);library(ggplot)
df=data.frame(type=IVtype,lp=bx%*%theta,y=c(by),lpse=1,yse=c(byse))
df %>% filter(type=='valid') %>% select(lp,y) %>% as.matrix() %>% cor()
ao=lm(y~lp,data=df %>% filter(type=='valid'))
df %>%
  ggplot(aes(lp,y,fill=type)) +
  geom_errorbar(aes(ymin=y-2*yse,ymax=y+2*yse),width=1/5,color='gray70') +
  geom_errorbarh(aes(xmin=lp-2*lpse,xmax=lp+2*lpse),height=1/5,color='gray70') +
  geom_point(pch=21,size=3) +
  scale_fill_manual('IV type',values=c('blue','orange','#39d445','red')) +
  theme_classic() +
  theme(legend.position='bottom',legend.text=element_text(size=14)) +
  labs(x='linear predictor for SNP-exposure association estimates',
       y='SNP-outcome association estimates') +
  guides(fill=guide_legend(override.aes=list(size=3))) +
  geom_hline(yintercept=0,linetype='dashed') +
  geom_vline(xintercept=0,linetype='dashed') +
  geom_abline(intercept=coef(ao)[1],slope=coef(ao)[2],col='#31b33b',lwd=1)


