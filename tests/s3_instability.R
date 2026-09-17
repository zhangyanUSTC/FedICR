root <- normalizePath(".")
options(device = function(...) grDevices::pdf(file = NULL))
read_core <- function() {
  e <- new.env(parent = globalenv())
  suppressPackageStartupMessages(sys.source("simulation/S0_basic.R", envir = e))
  e
}
core <- read_core()
original_update <- core$update_Z
original_ri <- core$RI_calculate
for (M in c(20L, 40L, 50L)) for (seed in c(1L, 17L)) {
  X <- array(rep(seq_len(M), 12L), c(M, 6L, 2L))
  Y <- matrix(0, M, 6L)
  z0 <- array(0, c(M, 6L, 2L))
  set.seed(seed)
  fold <- sample(rep_len(1:3, M))
  expected_rng <- .Random.seed
  for (helper in c("cva_M", "cva_M_forAA_alg")) {
    train_ids <- test_ids <- list()
    fit <- function(X, Y, beta_init, ...) {
      train_ids[[length(train_ids) + 1L]] <<- as.integer(X[, 1L, 1L])
      list(beta_hat = beta_init, l = 0)
    }
    core$update_Z <- function(X, Y, beta_hat) {
      test_ids[[length(test_ids) + 1L]] <<- as.integer(X[, 1L, 1L])
      NULL
    }
    core$RI_calculate <- function(z1, z2) 0.75
    score <- core[[helper]](X, Y, diag(2), z0, diag(2), 2L,
                           fit_alg = fit, seed = seed, shred = 1e-4, iters = 2L)
    stopifnot(identical(score, 0.25), identical(.Random.seed, expected_rng),
              identical(train_ids, list(which(fold == 1L), which(fold == 2L))),
              identical(test_ids, rep(list(which(fold == 3L)), 2L)))
  }
}
core$update_Z <- original_update
core$RI_calculate <- original_ri
cat("PASS: balanced client folds, full test sets, normalized s and RNG state.\n")

# Compare s/min with RI/max using exactly the same client splits and fits.
reference <- read_core()
as_ri <- function(x) {
  if (!is.call(x)) return(x)
  if (identical(x[[1L]], as.name("-")) && length(x) == 3L &&
      identical(x[[2L]], 1) && is.call(x[[3L]]) &&
      identical(x[[3L]][[1L]], as.name("RI_calculate"))) return(x[[3L]])
  if (identical(x[[1L]], as.name("which.min"))) x[[1L]] <- as.name("which.max")
  as.call(lapply(as.list(x), as_ri))
}
for (name in c("cva_M", "cva_M_forAA_alg", "cva_n", "S3_main")) {
  f <- reference[[name]]
  body(f) <- as_ri(body(f))
  reference[[name]] <- f
}
expected <- vector("list", 2L)
for (seed in 1:2) {
  arguments <- list(beta0 = 10 * diag(4), num_init = 3, seed = seed,
                    iters = 30, n = 120, M = 10, shred = 1e-4,
                    K_choose = 2:6, evaluate_selected = TRUE)
  old <- do.call(reference$S3_main, arguments)
  old_rng <- .Random.seed
  current <- do.call(core$S3_main, arguments)
  stopifnot(identical(.Random.seed, old_rng),
            identical(old$selected_metrics, current$selected_metrics))
  for (method in names(current$cva)) {
    stopifnot(identical(current$cva[[method]][-1L], 1 - old$cva[[method]][-1L]),
              identical(current$cva[[method]][1L], old$cva[[method]][1L]))
  }
  expected[[seed]] <- current
}
cat("PASS: s/min and RI/max select identical K and refits without extra RNG draws.\n")

source("scripts/control.R")
simulation_settings <- function(task) data.frame(
  task = "S3", n = 120L, M = 10L, SS = 10, script = "S3_simu_cva.R")
out <- tempfile("s3_instability_")
run_project(root, c("--task", "S3", "--seeds", "2", "--workers", "2", "--output", out))
files <- list.files(out, "^S3_n.*[.]Rdata$", full.names = TRUE)
stopifnot(length(files) == 1L)
saved <- new.env()
load(files, envir = saved)
stopifnot(identical(saved$runs, expected), !exists("RI_out", envir = saved, inherits = FALSE),
          all(c("Mean_min_s", "SD_min_s") %in% colnames(saved$s_out)))
for (method in names(saved$result_matrices)) {
  scores <- saved$result_matrices[[method]][, -1L, drop = FALSE]
  complete <- apply(is.finite(scores), 1L, all)
  best <- apply(scores, 1L, function(x) if (all(is.finite(x))) min(x) else NA_real_)
  chosen <- apply(scores, 1L, function(x) {
    if (all(is.finite(x))) saved$K_choose[which.min(x)] else NA_integer_
  })
  stopifnot(identical(saved$K_out[[method]], chosen),
            all(is.finite(scores[complete, , drop = FALSE])),
            all(scores[complete, , drop = FALSE] >= 0 & scores[complete, , drop = FALSE] <= 1),
            isTRUE(all.equal(unname(saved$s_out[method, "Mean_min_s"]), mean(best, na.rm = TRUE))),
            isTRUE(all.equal(unname(saved$s_out[method, "SD_min_s"]), sd(best, na.rm = TRUE))))
}
log <- readLines(file.path(out, "output.txt"))
stopifnot(any(grepl("Mean_min_s", log)), any(grepl("SD_min_s", log)),
          !any(grepl("Mean_max_RI|SD_max_RI", log)),
          file.exists(file.path(out, "S3_selected_K.png")),
          sum(log == "\\begin{table}") == 1L,
          tail(log[nzchar(log)], 1L) == "\\end{table}",
          any(grepl("\\multicolumn{2}{c}{$M=10$}", log, fixed = TRUE)))
source("scripts/report_S3_latex.R")
invisible(capture.output(expected_latex <- print_s3_latex(out)))
table_start <- which(log == "\\begin{table}")
stopifnot(identical(log[seq.int(table_start, length.out = length(expected_latex))], expected_latex))
cat("PASS: parallel runner, saved s scores, selection summaries, Figure 3 and automatic LaTeX output.\n")
cat("Synthetic test output:", out, "\n")
