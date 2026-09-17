# Synthetic integration test, not a paper experiment.
root <- normalizePath(".")
source(file.path(root, "scripts", "control.R"))
stopifnot(parse_options(character())$task == "help",
          parse_options(c("--seeds", "100"))$seeds == 100L)
for (bad in c("0", "-1", "1.5", "Inf", "abc"))
  stopifnot(inherits(try(parse_options(c("--seeds", bad)), silent = TRUE), "try-error"))
paper <- simulation_settings("all")
stopifnot(nrow(paper) == 72L, all(paper$n * paper$M == ifelse(paper$task == "S3",12000,4000)))
# Override settings only in this test.
simulation_settings <- function(task) data.frame(
  task = c("S1_M", "S1_N", "S1_SS", "S2", "S3"), n = 120L, M = 10L,
  SS = c(10,5,110,10,10),
  script = c("S1_simu.R", "S12_simu_ppfl.R", "S12_simu.R", "S2_simu_unba.R", "S3_simu_cva.R"))
out <- tempfile("E_FICR_smoke_")
run_project(root, c("--task", "all", "--seeds", "2", "--workers", "2", "--output", out))
files <- list.files(out, "^(S1_M_RES|S1_N_n|S1_SS_n|S2_n|S3_n).*[.]Rdata$", recursive = TRUE, full.names = TRUE)
task_names <- c("S1_M", "S1_N", "S1_SS", "S2", "S3")
stopifnot(setequal(basename(list.dirs(out, recursive = FALSE)), task_names),
          dir.exists(file.path(out, "S1_SS", "n120_m10")),
          !file.exists(file.path(out, "output.txt")),
          file.exists(file.path(out, "run_info.rds")))
for (task_name in task_names) {
  log <- readLines(file.path(out, task_name, "output.txt"), warn = FALSE)
  banners <- log[grepl("^Running ", log)]
  stopifnot(length(banners) == 1L, startsWith(banners, paste0("Running ", task_name, ":")))
}
stopifnot(length(files) == 5L)
for (f in files) {
  e <- new.env()
  load(f, e)
  stopifnot(length(e$runs) == 2L, !any(vapply(e$runs, inherits, logical(1), what = "error")))
  task <- sub("_(RES_)?n[0-9].*$", "", basename(f))
  fun <- switch(task, S1_M = S1_main, S1_N = S12_main, S1_SS = S12_main,
                S2 = S2_main, S3 = S3_main)
  args <- list(beta0 = e$SS * diag(4), num_init = 3, seed = 1,
               iters = 30, n = 120, M = 10, shred = .0001, balance = task != "S2")
  if (task == "S3") args <- c(args, list(K_choose = 2:6, evaluate_selected = TRUE)) else args$K_hat <- 4
  if (task %in% c("S1_N", "S1_SS")) args$run_sc2 <- FALSE
  if (task == "S1_SS") args$run_ppfl <- FALSE
  reference <- do.call(fun, args)
  stopifnot(isTRUE(all.equal(reference, e$runs[[1L]], tolerance = 0)))
}
cat("PASS: all five runners, saved results and serial/parallel seed-1 agreement.\n")
cat("Synthetic test output:", out, "\n")
