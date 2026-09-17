root <- normalizePath(".")
source(file.path(root, "scripts", "control.R"))
simulation_settings <- function(task) data.frame(
  task = "S1_SS", n = 120L, M = 10L, SS = c(110, 120), script = "S12_simu.R")
out <- tempfile("s1_ss_test_")
run_project(root, c("--task", "S1_SS", "--seeds", "2", "--workers", "2", "--output", out))
group_dir <- file.path(out, "n120_m10")
fs <- list.files(group_dir, "^S1_SS_n.*[.]Rdata$", full.names = TRUE)
stopifnot(length(fs) == 2L, file.exists(file.path(group_dir, "S1_SS_M10.png")),
          file.exists(file.path(group_dir, "S1_SS_figure_data.Rdata")),
          length(list.files(out, "[.](Rdata|png)$")) == 0L)
for (f in fs) {
  e <- new.env()
  load(f, e)
  stopifnot(identical(names(e$report_table),
    c("method", "failed", "MSE_mean", "MSE_SD", "Acc_mean", "Acc_SD", "ell_mean", "ell_SD")),
    identical(e$report_table$method, c("Oracle", "CICR", "E-FICR", "AA", "SC")),
    identical(e$report_table$failed, e$success_counts$failed[1:5]),
    identical(e$report_table$ell_mean, unname(e$all[,3])),
    identical(e$success_counts$method, c("Oracle", "CICR", "E_FICR", "AA", "SC")),
    !any(e$errors$method == "PPFL"), ncol(e$r) == 5L,
    all(vapply(e$runs, function(x) identical(names(x),
      c("Oracle_out", "CICR_out", "E_FICR_out", "AA_out", "SC_out")), logical(1))))
}
log <- readLines(file.path(out, "output.txt"), warn = FALSE)
stopifnot(!any(grepl("SC2|PPFL|requested|successful", log)))
args <- list(beta0 = e$SS * diag(4), num_init = 3, seed = 1, iters = 30,
             K_hat = 4, n = 120, M = 10, shred = .0001, balance = TRUE)
full <- do.call(S12_main, args)
full_rng <- .Random.seed
reduced <- do.call(S12_main, c(args, list(run_sc2 = FALSE, run_ppfl = FALSE)))
stopifnot(identical(.Random.seed, full_rng), identical(full[names(reduced)], reduced),
          identical(e$runs[[1]], reduced))
ppfl_original <- PPFL_alg
PPFL_alg <- function(...) stop("PPFL must not run in S1_SS")
without_ppfl <- do.call(S12_main, c(args, list(run_sc2 = FALSE, run_ppfl = FALSE)))
PPFL_alg <- ppfl_original
stopifnot(identical(without_ppfl, reduced), identical(.Random.seed, full_rng))
cat("PASS: S1_SS skips SC2 and PPFL; summaries, plots, other estimates and RNG state verified.\n")
