> [!NOTE]
> `simmrd` can be run from the command line using `pixi` in the [`cli/`](cli/) repo. Follow the [link](cli)/ to learn more.

# Overview

`simmrd` is an R package for generating simulated GWAS data for univariable or multivariable Mendelian Randomization (MR) under various conditions of weak instrument bias, GWAS sample overlap, uncorrelated horizontal pleiotropy (UHP), and correlated horizontal pleiotropy (CHP).

These simulated data are intended to evaluate the statistical properties of existing or candidate MR methods.

# TL;DR

```R
# one exposure, default params
params <- set_params()
simulated_data   <- generate_summary(params)
```

> [!NOTE]
>  Both `simmrd::generate_individual()` and `simmrd::generate_summary()` generate simulated data for a **single** simulation replicate.

# Installing `simmrd`
```R
remotes::install_github('noahlorinczcomi/simmrd')
```

# Command-line interface

A CLI for running simulations from a YAML parameter file — including multi-iteration Monte Carlo runs — is available in [`cli/`](cli/). It uses [pixi](https://pixi.sh) to manage the R and Python environment.

```bash
cd cli/
pixi install && pixi run setup
pixi run simulate --params params/example_summary.yaml --output results.rds -n 500 --seed 42
```

See [`cli/README.md`](cli/README.md) for full documentation.

# Quick start

`set_params()` builds a parameter list with sensible defaults. You only need to specify the values you want to change.

```R
library(simmrd)

# one exposure, all defaults
params <- set_params()
data   <- generate_summary(params)

# two exposures with CHP, no GWAS overlap
params <- set_params(
  number_of_exposures        = 2,
  true_causal_effects        = c(0.3, 0.1),
  prop_gwas_overlap_Xs_and_Y = 0,
  number_of_CHP_causal_SNPs  = 20,
  ratio_of_CHP_variance      = 0.25,
  CHP_correlation            = -0.5
)
data <- generate_summary(params)

# individual-level data with weak instruments and confounding
params <- set_params(
  type                       = "individual",
  number_of_exposures        = 2,
  Y_variance_explained_by_Xs = c(0, 0.5),
  signs_of_causal_effects    = c(1, 1),
  Xs_variance_explained_by_U = 0.12,
  Y_variance_explained_by_U  = 0.10,
  simtype                    = "weak",
  fix_Fstatistic_at          = 10
)
data <- generate_individual(params)
```

`set_params()` validates your inputs immediately and returns a clear error message if something is misconfigured, rather than failing deep inside the simulation.

## Key parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `type` | `"summary"` | `"summary"` or `"individual"` |
| `sample_size_Xs` | `1e5` | Exposure GWAS sample size(s) |
| `sample_size_Y` | `1e5` | Outcome GWAS sample size |
| `number_of_exposures` | `1` | Number of exposures |
| `number_of_causal_SNPs` | `100` | SNPs with a direct effect on each exposure |
| `prop_gwas_overlap_Xs_and_Y` | `1` | Proportion of overlap between exposure and outcome GWAS |
| `true_causal_effects` | `0.3` | True causal effect size(s) |
| `number_of_UHP_causal_SNPs` | `0` | Number of UHP SNPs |
| `number_of_CHP_causal_SNPs` | `0` | Number of CHP SNPs |
| `simtype` | `"winners"` | `"winners"` (P-value selection) or `"weak"` (fix F-statistic) |
| `IV_Pvalue_threshold` | `5e-8` | P-value threshold for IV selection |

See `?set_params` for the full parameter list.

# Built-in presets

`load_preset()` reproduces any of the named simulation scenarios used in the paper. Use `list_presets()` to see all options.

```R
list_presets()
#> load_preset() options
#> ---------------------
#>   bias         none, UHP, CHP, UHP_CHP, UHP_CHP_WEAK, WEAK
#>   n            3e4, 1e5
#>   snps         100, 500
#>   exposures    1, 3
#>   overlap      full, none
```

```R
# load a named scenario
params <- load_preset("CHP", n = 1e5, snps = 100, exposures = 1, overlap = "none")
data   <- generate_summary(params)

# load a preset and override one value
params <- load_preset("UHP_CHP", n = 1e5, snps = 500, exposures = 3, overlap = "full")
params$true_causal_effects <- c(0.1, 0.2, 0.3)
data <- generate_summary(params)
```

Preset pleiotropy parameters scale automatically with `snps` — 20% of SNPs are pleiotropic for single-bias scenarios (UHP or CHP), and 10% each for combined scenarios (UHP_CHP, UHP_CHP_WEAK).

# Output

Both `generate_summary()` and `generate_individual()` return a named list:

| Element | Description |
|---------|-------------|
| `bx` | $m \times p$ matrix of IV–exposure associations |
| `by` | $m \times 1$ vector of IV–outcome associations |
| `bxse` | Standard errors for `bx` |
| `byse` | Standard errors for `by` |
| `RhoME` | $(p+1) \times (p+1)$ measurement-error correlation matrix |
| `LDMatrix` | True LD correlation matrix among IVs |
| `LDhatMatrix` | Estimated LD correlation matrix among IVs |
| `theta` | True causal effects |
| `IVtype` | Per-IV classification: `"valid"`, `"UHP"`, or `"CHP"` |
| `bx_unstd` / `by_unstd` | Unstandardized versions of `bx` / `by` |
| `bxse_unstd` / `byse_unstd` | Standard errors for unstandardized estimates |

## Paper
Lorincz-Comi, N., Yang, Y., & Zhu, X. simmrd: An open-source tool to perform simulations in Mendelian randomization. Genetic epidemiology, _1-15_. [doi.org/10.1002/gepi.22544](https://doi.org/10.1002/gepi.22544)
