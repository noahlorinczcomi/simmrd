test_that("set_params returns a list with expected elements for summary", {
  p <- set_params()
  expect_type(p, "list")
  expect_true(all(c("sample_size_Xs", "sample_size_Y", "number_of_exposures",
                     "true_causal_effects", "simtype") %in% names(p)))
})

test_that("set_params returns a list with expected elements for individual", {
  p <- set_params(type = "individual")
  expect_type(p, "list")
  expect_true(all(c("Y_variance_explained_by_Xs", "signs_of_causal_effects",
                     "Xs_variance_explained_by_U") %in% names(p)))
})

test_that("set_params validates bad inputs", {
  expect_error(set_params(type = "bad"))
  expect_error(set_params(simtype = "bad"))
  expect_error(set_params(number_of_UHP_causal_SNPs = 60,
                           number_of_CHP_causal_SNPs = 60,
                           number_of_causal_SNPs     = 100))
})

test_that("list_presets returns invisibly and prints output", {
  expect_output(list_presets())
  opts <- list_presets()
  expect_named(opts, c("bias", "n", "snps", "exposures", "overlap"))
})

test_that("load_preset returns a valid parameter list", {
  p <- load_preset("CHP", n = 3e4, snps = 100, exposures = 1, overlap = "none")
  expect_type(p, "list")
  expect_equal(p$sample_size_Xs, 3e4)
  expect_equal(p$number_of_CHP_causal_SNPs, 20)
  expect_error(load_preset("BAD"))
})

test_that("generate_summary returns the expected output structure", {
  skip_on_cran()
  p    <- set_params(number_of_causal_SNPs = 50, sample_size_Xs = 1e4, sample_size_Y = 1e4)
  data <- generate_summary(p, seed = 1)
  expect_type(data, "list")
  expect_true(all(c("bx", "bxse", "by", "byse", "RhoME", "LDMatrix",
                     "LDhatMatrix", "theta", "IVtype",
                     "bx_unstd", "bxse_unstd", "by_unstd", "byse_unstd",
                     "beta_true", "alpha_true", "u", "v", "iv_index") %in% names(data)))
  expect_equal(ncol(data$bx), 1)
  expect_equal(length(data$by), nrow(data$bx))
})

test_that("generate_summary seed produces reproducible results", {
  skip_on_cran()
  p  <- set_params(number_of_causal_SNPs = 50, sample_size_Xs = 1e4, sample_size_Y = 1e4)
  d1 <- generate_summary(p, seed = 42)
  d2 <- generate_summary(p, seed = 42)
  expect_equal(d1$bx, d2$bx)
  expect_equal(d1$by, d2$by)
})

test_that("generate_individual returns the expected output structure", {
  skip_on_cran()
  p    <- set_params(type = "individual", number_of_causal_SNPs = 30,
                     sample_size_Xs = 5e3, sample_size_Y = 5e3,
                     simtype = "weak", fix_Fstatistic_at = 10)
  data <- generate_individual(p, seed = 1)
  expect_type(data, "list")
  expect_true(all(c("bx", "bxse", "by", "byse", "RhoME", "LDMatrix",
                     "LDhatMatrix", "theta", "IVtype") %in% names(data)))
})

test_that("adj_overlap returns a labelled square matrix", {
  m <- adj_overlap(0.2, 0.5, 3)
  expect_equal(dim(m), c(4, 4))
  expect_equal(rownames(m), c("Outcome", "Exposure1", "Exposure2", "Exposure3"))
})
