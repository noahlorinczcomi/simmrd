#' Set simulation parameters
#'
#' Constructs a named parameter list for use with \code{generate_summary()} or
#' \code{generate_individual()}.  Every argument has a sensible default so you
#' only need to specify the values you want to change.
#'
#' @param type \code{"summary"} (default) or \code{"individual"}.
#'
#' --- Study design ---
#' @param sample_size_Xs Exposure GWAS sample size(s). Scalar or vector with one value per exposure.
#' @param sample_size_Y Outcome GWAS sample size.
#' @param number_of_exposures Number of exposures.
#' @param number_of_causal_SNPs Number of SNPs with a direct effect on each exposure.
#'
#' --- GWAS overlap ---
#' @param prop_gwas_overlap_Xs_and_Y Proportion of overlap between exposure and outcome GWAS. Scalar or vector.
#' @param prop_gwas_overlap_Xs Proportion of overlap among the exposure GWAS (summary only). Scalar or numeric matrix.
#'
#' --- Pleiotropy (summary data) ---
#' @param number_of_UHP_causal_SNPs Number of uncorrelated horizontal pleiotropy (UHP) SNPs.
#' @param number_of_CHP_causal_SNPs Number of correlated horizontal pleiotropy (CHP) SNPs.
#' @param ratio_of_UHP_variance Ratio of UHP variance to valid-IV variance.
#' @param ratio_of_CHP_variance Ratio of CHP variance to valid-IV variance.
#' @param CHP_correlation Correlation between CHP and valid-IV effect sizes (magnitude of CHP).
#'
#' --- Pleiotropy (individual data) ---
#' @param Y_variance_explained_by_UHP Outcome variance explained by UHP SNPs.
#' @param U_variance_explained_by_CHP Confounder variance explained by CHP SNPs.
#'
#' --- Causal effects ---
#' @param true_causal_effects True causal effect size(s). Scalar or vector (summary only).
#' @param Y_variance_explained_by_Xs Outcome variance explained by each exposure. Scalar or vector (individual only).
#' @param signs_of_causal_effects Signs of causal effects. Scalar or vector (individual only).
#'
#' --- Correlations ---
#' @param phenotypic_correlation_Xs Phenotypic correlations among exposures. Scalar, string (\code{'ar1(0.5)'}, \code{'cs(0.3)'}, \code{'toeplitz'}), or matrix.
#' @param genetic_correlation_Xs Genetic correlations among exposures. Same formats as above.
#' @param phenotypic_correlations_Xs_and_Y Phenotypic correlations between each exposure and the outcome. Scalar or vector (summary only).
#'
#' --- Genetic architecture ---
#' @param Xs_variance_explained_by_g Heritability of each exposure (variance explained by all causal SNPs). Scalar or vector.
#' @param LD_causal_SNPs LD structure among causal SNPs. Scalar, string (\code{'I'}, \code{'toeplitz'}, \code{'ar1(0.5)'}), or numeric matrix.
#' @param number_of_LD_blocks Number of independent LD blocks.
#'
#' --- Confounding (individual only) ---
#' @param Xs_variance_explained_by_U Exposure variance explained by the latent confounder.
#' @param Y_variance_explained_by_U Outcome variance explained by the latent confounder.
#'
#' --- IV selection ---
#' @param simtype \code{"winners"} (P-value-based selection) or \code{"weak"} (fix F-statistic).
#' @param IV_Pvalue_threshold P-value threshold for IV selection (used when \code{simtype = "winners"}).
#' @param fix_Fstatistic_at Target mean F-statistic (used when \code{simtype = "weak"}).
#' @param MVMR_IV_selection_type \code{"union"} or \code{"joint"} (multivariable MR only).
#' @param LD_pruning_r2 Upper r² threshold for LD pruning of IVs.
#'
#' --- Output ---
#' @param MR_standardization Standardization applied to GWAS summary statistics. \code{"none"}, \code{"Z"}, or \code{"QC"}.
#' @param N_of_LD_ref Sample size of the LD reference panel (\code{Inf} = use true LD).
#'
#' @return A named list of parameters ready to pass to \code{generate_summary()} or \code{generate_individual()}.
#' @export
#' @examples
#' # Minimal: one exposure, default settings
#' params <- set_params()
#' data <- generate_summary(params)
#'
#' # Two exposures with CHP, no GWAS overlap
#' params <- set_params(
#'   number_of_exposures        = 2,
#'   true_causal_effects        = c(0.3, 0.1),
#'   prop_gwas_overlap_Xs_and_Y = 0,
#'   number_of_CHP_causal_SNPs  = 20,
#'   ratio_of_CHP_variance      = 0.25,
#'   CHP_correlation            = -0.5
#' )
#' data <- generate_summary(params)
set_params <- function(
  type = "summary",
  # study design
  sample_size_Xs = 1e5,
  sample_size_Y = 1e5,
  number_of_exposures = 1,
  number_of_causal_SNPs = 100,
  # overlap
  prop_gwas_overlap_Xs_and_Y = 1,
  prop_gwas_overlap_Xs = 1,
  # pleiotropy (summary)
  number_of_UHP_causal_SNPs = 0,
  number_of_CHP_causal_SNPs = 0,
  ratio_of_UHP_variance = 0,
  ratio_of_CHP_variance = 0,
  CHP_correlation = 0,
  # pleiotropy (individual)
  Y_variance_explained_by_UHP = 0,
  U_variance_explained_by_CHP = 0,
  # causal effects
  true_causal_effects = 0.3,
  Y_variance_explained_by_Xs = 0.3,
  signs_of_causal_effects = 1,
  # correlations
  phenotypic_correlation_Xs = 0.3,
  genetic_correlation_Xs = 0.15,
  phenotypic_correlations_Xs_and_Y = 0.3,
  # genetic architecture
  Xs_variance_explained_by_g = 0.10,
  LD_causal_SNPs = "I",
  number_of_LD_blocks = 1,
  # confounding (individual)
  Xs_variance_explained_by_U = 0,
  Y_variance_explained_by_U = 0,
  # IV selection
  simtype = "winners",
  IV_Pvalue_threshold = 5e-8,
  fix_Fstatistic_at = 10,
  MVMR_IV_selection_type = "union",
  LD_pruning_r2 = 1,
  # output
  MR_standardization = "none",
  N_of_LD_ref = Inf
) {
  type <- tolower(type)
  if (!type %in% c("summary", "individual")) {
    stop("`type` must be \"summary\" or \"individual\".")
  }

  # ---- validation --------------------------------------------------------
  if (number_of_UHP_causal_SNPs + number_of_CHP_causal_SNPs > number_of_causal_SNPs) {
    stop("`number_of_UHP_causal_SNPs` + `number_of_CHP_causal_SNPs` cannot exceed `number_of_causal_SNPs`.")
  }

  if (any(c(ratio_of_UHP_variance, ratio_of_CHP_variance) < 0)) {
    stop("`ratio_of_UHP_variance` and `ratio_of_CHP_variance` must be >= 0.")
  }

  if (!simtype %in% c("winners", "weak")) {
    stop('`simtype` must be "winners" or "weak".')
  }

  if (!tolower(MR_standardization) %in% c("none", "z", "qc")) {
    stop('`MR_standardization` must be "none", "Z", or "QC".')
  }

  if (!MVMR_IV_selection_type %in% c("union", "joint")) {
    stop('`MVMR_IV_selection_type` must be "union" or "joint".')
  }

  if (IV_Pvalue_threshold <= 0 | IV_Pvalue_threshold > 1) {
    stop("`IV_Pvalue_threshold` must be in (0, 1].")
  }

  if (LD_pruning_r2 <= 0 | LD_pruning_r2 > 1) {
    stop("`LD_pruning_r2` must be in (0, 1].")
  }

  # ---- build output list -------------------------------------------------
  if (type == "summary") {
    params <- list(
      sample_size_Xs = sample_size_Xs,
      sample_size_Y = sample_size_Y,
      number_of_exposures = number_of_exposures,
      number_of_causal_SNPs = number_of_causal_SNPs,
      prop_gwas_overlap_Xs_and_Y = prop_gwas_overlap_Xs_and_Y,
      prop_gwas_overlap_Xs = prop_gwas_overlap_Xs,
      number_of_UHP_causal_SNPs = number_of_UHP_causal_SNPs,
      number_of_CHP_causal_SNPs = number_of_CHP_causal_SNPs,
      ratio_of_UHP_variance = ratio_of_UHP_variance,
      ratio_of_CHP_variance = ratio_of_CHP_variance,
      CHP_correlation = CHP_correlation,
      true_causal_effects = true_causal_effects,
      phenotypic_correlation_Xs = phenotypic_correlation_Xs,
      genetic_correlation_Xs = genetic_correlation_Xs,
      phenotypic_correlations_Xs_and_Y = phenotypic_correlations_Xs_and_Y,
      Xs_variance_explained_by_g = Xs_variance_explained_by_g,
      LD_causal_SNPs = LD_causal_SNPs,
      number_of_LD_blocks = number_of_LD_blocks,
      simtype = simtype,
      IV_Pvalue_threshold = IV_Pvalue_threshold,
      fix_Fstatistic_at = fix_Fstatistic_at,
      MVMR_IV_selection_type = MVMR_IV_selection_type,
      LD_pruning_r2 = LD_pruning_r2,
      MR_standardization = MR_standardization,
      N_of_LD_ref = N_of_LD_ref
    )
  } else {
    params <- list(
      sample_size_Xs              = sample_size_Xs,
      sample_size_Y               = sample_size_Y,
      number_of_exposures         = number_of_exposures,
      number_of_causal_SNPs       = number_of_causal_SNPs,
      prop_gwas_overlap_Xs_and_Y  = prop_gwas_overlap_Xs_and_Y,
      number_of_UHP_causal_SNPs   = number_of_UHP_causal_SNPs,
      number_of_CHP_causal_SNPs   = number_of_CHP_causal_SNPs,
      Y_variance_explained_by_UHP = Y_variance_explained_by_UHP,
      U_variance_explained_by_CHP = U_variance_explained_by_CHP,
      Y_variance_explained_by_Xs  = Y_variance_explained_by_Xs,
      signs_of_causal_effects     = signs_of_causal_effects,
      phenotypic_correlation_Xs   = phenotypic_correlation_Xs,
      genetic_correlation_Xs      = genetic_correlation_Xs,
      Xs_variance_explained_by_g  = Xs_variance_explained_by_g,
      Xs_variance_explained_by_U  = Xs_variance_explained_by_U,
      Y_variance_explained_by_U   = Y_variance_explained_by_U,
      LD_causal_SNPs              = LD_causal_SNPs,
      number_of_LD_blocks         = number_of_LD_blocks,
      simtype                     = simtype,
      IV_Pvalue_threshold         = IV_Pvalue_threshold,
      fix_Fstatistic_at           = fix_Fstatistic_at,
      MVMR_IV_selection_type      = MVMR_IV_selection_type,
      LD_pruning_r2               = LD_pruning_r2,
      MR_standardization          = MR_standardization,
      N_of_LD_ref                 = N_of_LD_ref
    )
  }

  params
}

#' Load a named simulation preset
#'
#' Returns a ready-to-use parameter list corresponding to one of the built-in
#' simulation scenarios.  The list is identical to what \code{set_params()}
#' produces, so every element can be overridden afterwards.
#'
#' @param bias Bias scenario. One of \code{"none"}, \code{"UHP"},
#'   \code{"CHP"}, \code{"UHP_CHP"}, \code{"UHP_CHP_WEAK"}, or \code{"WEAK"}.
#' @param n GWAS sample size for both exposures and outcome. \code{3e4} or \code{1e5}.
#' @param snps Number of causal SNPs per exposure. \code{100} or \code{500}.
#' @param exposures Number of exposures. \code{1} or \code{3}.
#' @param overlap \code{"full"} (complete exposure–outcome overlap) or
#'   \code{"none"} (no exposure–outcome overlap).
#'
#' @return A named parameter list, identical in structure to \code{set_params()} output.
#' @export
#' @seealso \code{\link{set_params}}, \code{\link{list_presets}}
#' @examples
#' # CHP scenario, small GWAS, no overlap
#' params <- load_preset("CHP", n = 3e4, snps = 100, exposures = 1, overlap = "none")
#' data <- generate_summary(params)
#'
#' # Start from a preset, then tweak one thing
#' params <- load_preset("UHP_CHP", n = 1e5, snps = 500, exposures = 3, overlap = "full")
#' params$true_causal_effects <- c(0.1, 0.2, 0.3)
#' data <- generate_summary(params)
load_preset <- function(
  bias = "none",
  n = 1e5,
  snps = 100,
  exposures = 1,
  overlap = "full"
) {
  valid_bias <- c("none", "UHP", "CHP", "UHP_CHP", "UHP_CHP_WEAK", "WEAK")
  if (!bias %in% valid_bias) {
    stop("`bias` must be one of: ", paste(valid_bias, collapse = ", "), ".")
  }
  if (!overlap %in% c("full", "none")) {
    stop('`overlap` must be "full" or "none".')
  }

  # overlap params: noOverlap keeps some exposure-exposure overlap when p > 1
  prop_xy <- if (overlap == "full") 1 else 0
  prop_xs <- if (overlap == "full") 1 else if (exposures > 1) 0.3 else 1

  # pleiotropy SNP counts scale with snps; ratios follow the default setups
  bp <- switch(bias,
    none = list(uhp = 0, chp = 0, r_uhp = 0, r_chp = 0, corr = 0, simtype = "winners", fstat = 10),
    UHP = list(uhp = floor(0.2 * snps), chp = 0, r_uhp = 0.25, r_chp = 0, corr = 0, simtype = "winners", fstat = 10),
    CHP = list(uhp = 0, chp = floor(0.2 * snps), r_uhp = 0, r_chp = 0.25, corr = -0.5, simtype = "winners", fstat = 10),
    UHP_CHP = list(uhp = floor(0.1 * snps), chp = floor(0.1 * snps), r_uhp = 0.125, r_chp = 0.125, corr = 0, simtype = "winners", fstat = 10),
    UHP_CHP_WEAK = list(uhp = floor(0.1 * snps), chp = floor(0.1 * snps), r_uhp = 0.125, r_chp = 0.125, corr = -0.5, simtype = "weak", fstat = 5),
    WEAK = list(uhp = 0, chp = 0, r_uhp = 0, r_chp = 0, corr = 0, simtype = "weak", fstat = 5)
  )

  set_params(
    sample_size_Xs             = n,
    sample_size_Y              = n,
    number_of_exposures        = exposures,
    number_of_causal_SNPs      = snps,
    prop_gwas_overlap_Xs_and_Y = prop_xy,
    prop_gwas_overlap_Xs       = prop_xs,
    number_of_UHP_causal_SNPs  = bp$uhp,
    number_of_CHP_causal_SNPs  = bp$chp,
    ratio_of_UHP_variance      = bp$r_uhp,
    ratio_of_CHP_variance      = bp$r_chp,
    CHP_correlation            = bp$corr,
    simtype                    = bp$simtype,
    fix_Fstatistic_at          = bp$fstat,
    true_causal_effects        = 0.3
  )
}

#' List available simulation presets
#'
#' Prints the valid values for each argument of \code{\link{load_preset}}.
#'
#' @return Invisibly returns a named list of valid values.
#' @export
#' @examples
#' list_presets()
list_presets <- function() {
  opts <- list(
    bias = c("none", "UHP", "CHP", "UHP_CHP", "UHP_CHP_WEAK", "WEAK"),
    n = c(3e4, 1e5),
    snps = c(100, 500),
    exposures = c(1, 3),
    overlap = c("full", "none")
  )
  cat("load_preset() options\n")
  cat("---------------------\n")
  for (nm in names(opts)) {
    cat(sprintf("  %-12s %s\n", nm, paste(opts[[nm]], collapse = ", ")))
  }
  invisible(opts)
}
