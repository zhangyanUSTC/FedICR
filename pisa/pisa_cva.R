current_dir <- getwd()
outcome_dir <- getOption("E_FICR.output")
if (!dir.exists(outcome_dir)) {
  dir.create(outcome_dir)
  message("Directory 'outcome' created.")
} else {
  message("Directory 'outcome' already exists.")
}
time_stamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
out_path <- file.path(outcome_dir, paste0("Pisa_CVA_", time_stamp, ".txt"))
sink(out_path, split = TRUE)

section <- function(title) {
  cat("\n============================================================\n")
  cat(title, "\n")
  cat("============================================================\n")
}

print_amp_table <- function(x, digits = 5L) {
  x <- as.matrix(x)
  cat("Statistic & ", paste(colnames(x), collapse = " & "), "\n", sep = "")
  for (i in seq_len(nrow(x))) {
    cat(
      rownames(x)[i], " & ",
      paste(format(x[i, ], digits = digits, trim = TRUE), collapse = " & "),
      "\n", sep = ""
    )
  }
}

set.seed(1)
iters<-40;shred<-0.0001
cat("PISA CICR cluster-number stability analysis\n")
cat("Started:", format(Sys.time()), "\n")
cat("Iterations:", iters, "\n")
cat("Convergence threshold:", shred, "\n")
cat("Candidate K: 2, 3, 4, 5, 6\n")
cat("Repeats per K: 10\n")
federated_data <- readRDS(
  file.path(getOption("E_FICR.data"), "pisa_federated_data_pv1.rds")
)
suppressPackageStartupMessages(
  suppressMessages(
    source(file.path(getOption("E_FICR.root"), "pisa", "pisa_federated.R"))
  )
)
all_vars <- names(federated_data[[1]])
response_var <- tail(all_vars, 1)
predictor_vars <- setdiff(
  all_vars,
  c("CNT", response_var)
)
M <- length(federated_data)
p <- length(predictor_vars)
n_vec <- vapply(
  federated_data,
  nrow,
  integer(1)
)
X <- lapply(federated_data, function(client_data) {
  as.matrix(
    client_data[, predictor_vars, drop = FALSE]
  )
})
Y <- lapply(federated_data, function(client_data) {
  as.numeric(client_data[[response_var]])
})
names(X) <- names(federated_data)
names(Y) <- names(federated_data)

num_cores <- getOption("E_FICR.workers", 10L)
cl <- start_workers(getOption("E_FICR.workers", 10L))
registerDoParallel(cl)
K_choose<-2:6
cva_score <- matrix(
  NA_real_, nrow = 10L, ncol = length(K_choose),
  dimnames = list(paste0("Repeat", 1:10), paste0("K", K_choose))
)
for (kk in K_choose) {
  t1 <- Sys.time()
  section(paste("FITTING K =", kk))
  r <- foreach(
    i = 1:10,
    .combine = rbind,
    .packages = c(
      'clue',
      "mclust",
      "Matrix"
    ),
    .errorhandling = "stop"
  ) %dopar% {
    tryCatch(
      { set.seed(i)
        beta_init<-INIT_beta_M(X,Y,K_hat=kk,lam=0.001)
        cva_n(X,Y,beta_init=beta_init,K_hat=kk,iters=iters,shred=0.0001,fit_alg = CICR_alg)
      },
      error = function(e) {
        NA_real_
      }
    )
  }
  cva_score[, paste0("K", kk)] <- as.numeric(r)
  t2<- Sys.time()
  cat("Elapsed time:", format(t2 - t1), "\n")
  cat("Successful repeats:", sum(is.finite(r)), "of 10\n")
}
stopCluster(cl)

section("REPEAT-LEVEL CVA SCORES")
print_amp_table(cva_score)

cva_summary <- rbind(
  Mean = colMeans(cva_score, na.rm = TRUE),
  SD = apply(cva_score, 2, sd, na.rm = TRUE),
  Variance = apply(cva_score, 2, var, na.rm = TRUE),
  Min = apply(cva_score, 2, min, na.rm = TRUE),
  Max = apply(cva_score, 2, max, na.rm = TRUE),
  Successful_N = colSums(is.finite(cva_score))
)

section("CVA SUMMARY")
cat("Smaller CVA scores indicate more stable clustering.\n")
print_amp_table(cva_summary)

mean_scores <- cva_summary["Mean", ]
if (any(is.finite(mean_scores))) {
  best_index <- which.min(mean_scores)
  best_k <- K_choose[best_index]
  best_mean <- mean_scores[best_index]
  cat("\nSelected K by minimum mean CVA score:", best_k, "\n")
  cat(
    "Mean CVA score for selected K:",
    format(best_mean, digits = 2, trim = TRUE),
    "\n"
  )
}

save(cva_score, cva_summary, K_choose, iters, shred,
     file = sub("[.]txt$", ".Rdata", out_path))
cat("\nFinished:", format(Sys.time()), "\n")
cat("Results file:", normalizePath(out_path, mustWork = FALSE), "\n")

sink()
