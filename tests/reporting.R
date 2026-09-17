root <- normalizePath(".")
source(file.path(root, "scripts", "control.R"))
simulation_settings <- function(task) data.frame(
  task = "S1_M", n = 120L, M = 10L, SS = 10, script = "S1_simu.R")
out <- tempfile("reporting_test_")
run_project(root, c("--task", "S1_M", "--seeds", "2", "--workers", "2", "--output", out))
files <- list.files(out, "^S1_M_RES_.*[.]Rdata$", full.names = TRUE)
stopifnot(length(files) == 1L)
e <- new.env()
load(files, e)
stopifnot(identical(names(e$report_table), c("method", "failed", "MSE", "Acc", "ell")),
          identical(e$report_table$method, c("Oracle", "CICR", "E-FICR", "AA", "SC")),
          identical(e$report_table$failed, e$success_counts$failed[1:5]),
          identical(unname(as.matrix(e$report_table[, c("MSE", "Acc")])),
                    unname(as.matrix(e$result_table[, c("MSE", "Acc")]))),
          identical(e$report_table$ell, e$result_table$ell / (e$n * e$M)),
          e$iters == 30L)
log <- readLines(file.path(out, "output.txt"), warn = FALSE)
stopifnot(any(grepl("iters=30", log, fixed = TRUE)),
          !any(grepl("requested|successful|错误统计|ACC", log)))
expr <- as.list(parse(file.path(root, "simulation", "S1_simu.R")))
starts <- which(vapply(expr, function(x) is.call(x) && identical(x[[1L]], as.name("<-")) &&
                        identical(x[[2L]], as.name("method_map")), logical(1)))
e$success_counts <- e$success_counts[c(6, 3, 1, 5, 2, 4), ]
e$success_counts$failed <- c(9L, 2L, 0L, 4L, 1L, 3L)
for (i in starts:(starts + 2L)) eval(expr[[i]], e)
stopifnot(identical(e$report_table$failed, 0:4))
cat("PASS: merged S1_M output, name-aligned failure counts, iters and Acc labels.\n")
