#!/usr/bin/env Rscript
# runner.R – called by simmrd_cli.py; not intended to be run directly.

suppressPackageStartupMessages({
  library(simmrd)
  library(yaml)
  library(optparse)
})

# Argument parsing
option_list <- list(
  make_option("--params",
    type = "character", metavar = "FILE",
    help = "Path to YAML parameter file [required]"
  ),
  make_option("--output",
    type = "character", metavar = "FILE",
    help = "Output RDS file path [required]"
  ),
  make_option("--iterations",
    type = "integer", default = 1L,
    metavar = "INT",
    help = "Number of simulation iterations [default: 1]"
  ),
  make_option("--seed",
    type = "integer", default = NULL, metavar = "INT",
    help = "Master random seed for reproducibility [default: none]"
  )
)

opt <- parse_args(OptionParser(
  usage       = "%prog --params FILE --output FILE [-n INT] [--seed INT]",
  option_list = option_list
))

if (is.null(opt$params)) stop("--params is required")
if (is.null(opt$output)) stop("--output is required")
if (opt$iterations < 1L) stop("--iterations must be >= 1")

# Load and validate the YAML parameter file
raw <- tryCatch(
  yaml::read_yaml(opt$params),
  error = function(e) stop("Could not parse YAML file: ", conditionMessage(e))
)

sim_type <- tolower(if (!is.null(raw$type)) raw$type else "summary")
raw$type <- NULL

known_params <- setdiff(names(formals(set_params)), "type")
unknown <- setdiff(names(raw), known_params)
if (length(unknown) > 0) {
  warning("Unknown parameters in YAML (ignored): ",
    paste(unknown, collapse = ", "),
    call. = FALSE
  )
  raw <- raw[intersect(names(raw), known_params)]
}

params <- tryCatch(
  do.call(set_params, c(list(type = sim_type), raw)),
  error = function(e) stop("Invalid parameter value:\n  ", conditionMessage(e))
)

# ---------------------------------------------------------------------------
# Seeding strategy
#
# Both generate_summary() and generate_individual() accept a `seed` argument.
# We generate a per-iteration seed vector so that:
#   --seed given  →  fully reproducible across all iterations and re-runs
#   --seed absent →  a fresh random sequence each run
#
# This also prevents generate_individual()'s internal set.seed(1) default
# from making every iteration identical when looping.
# ---------------------------------------------------------------------------
n_iter <- opt$iterations

if (!is.null(opt$seed)) set.seed(opt$seed)
iter_seeds <- sample.int(.Machine$integer.max, n_iter)

# Simulation loop
cat(sprintf(
  "Running generate_%s() x %d iteration(s) ...\n", sim_type, n_iter
))

results <- vector("list", n_iter)
pb <- utils::txtProgressBar(min = 0, max = n_iter, style = 3)

for (i in seq_len(n_iter)) {
  results[[i]] <- if (sim_type == "individual") {
    generate_individual(params, seed = iter_seeds[i])
  } else {
    generate_summary(params, seed = iter_seeds[i])
  }
  utils::setTxtProgressBar(pb, i)
}

close(pb)
cat("\n")

# Summary to stdout
iv_counts <- sapply(results, function(r) length(r$by))
cat(sprintf("Iterations   : %d\n", n_iter))
cat(sprintf(
  "IVs selected : min=%d  median=%g  max=%d\n",
  min(iv_counts), median(iv_counts), max(iv_counts)
))

# Output (RDS format)
saveRDS(results, opt$output)
cat(sprintf(
  "Saved to: %s  [format: RDS]\n",
  opt$output
))
