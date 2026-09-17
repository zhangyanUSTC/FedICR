root <- normalizePath(".")
source(file.path(root, "scripts", "control.R"))
simulation_settings <- function(task) data.frame(
  task = "S1_N", n = 120L, M = 10L, SS = 5, script = "S12_simu_ppfl.R")
out <- tempfile("s1_n_test_")
run_project(root, c("--task", "S1_N", "--seeds", "2", "--workers", "2", "--output", out))
e <- new.env()
load(list.files(out, "^S1_N_n.*[.]Rdata$", full.names = TRUE), e)
stopifnot(identical(names(e$report_table),
                    c("method", "failed", "MSE_mean", "MSE_SD", "Acc_mean", "Acc_SD")),
          identical(e$report_table$method, c("Oracle", "CICR", "E-FICR", "AA", "SC", "PPFL")),
          ncol(e$r) == 6L, length(e$runs[[1]]) == 6L,
          !"SC_out2" %in% names(e$runs[[1]]))
log <- readLines(file.path(out, "output.txt"), warn = FALSE)
stopifnot(!any(grepl("SC2|ell|requested|successful", log)))
args <- list(beta0 = 5 * diag(4), num_init = 3, seed = 1, iters = 30,
             K_hat = 4, n = 120, M = 10, shred = .0001, balance = TRUE)
full <- do.call(S12_main, args)
full_rng <- .Random.seed
reduced <- do.call(S12_main, c(args, list(run_sc2 = FALSE)))
stopifnot(identical(.Random.seed, full_rng),
          identical(full[names(reduced)], reduced),
          identical(e$runs[[1]], reduced))
cat("PASS: S1_N omits SC2, reports six methods, and preserves other estimates and RNG state.\n")
