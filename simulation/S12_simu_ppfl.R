num_cores <- min(seed_all, getOption("E_FICR.workers", 10L))
cl <- start_workers(num_cores)
registerDoParallel(cl)
error_counts <- list()

t1 <- Sys.time()
error_count <- 0
runs <- tryCatch(foreach(
  i = seq_len(seed_all),
  .packages = c(
    'clue',
    "MASS",
    "foreach",
    "combinat",
    "mclust",
    "Matrix"
  ),
  .errorhandling = "pass"
) %dopar% {
  S12_main(beta0[1:K, , drop = FALSE], num_init = num_init, seed = i,
           iters = iters, K_hat = K, n = n, M = M, shred = 0.0001,
           balance = TRUE, run_sc2 = FALSE)
}, finally = {
  stopCluster(cl)
  foreach::registerDoSEQ()
})
failed <- vapply(runs, inherits, logical(1), what = "error")
errors <- data.frame(seed = which(failed),
                     method = rep("seed", sum(failed)),
                     message = vapply(runs[failed], conditionMessage, character(1)))
error_count <- sum(failed)
error_counts[[as.character(K)]] <- error_count
error_path <- file.path(file_path, sprintf(
  "S1_N_errors_n%d_m%d_gk%d_p%d_s%d_num%d_SS%d.csv", n, M, K, p, seed_all, num_init, SS))
if (nrow(errors)) write.csv(errors, error_path, row.names = FALSE)
if (all(failed)) {
  save(runs, errors, n, M, K, p, seed_all, num_init, SS, iters,
       file = sub("[.]csv$", ".Rdata", error_path))
  stop("All repetitions failed; raw runs saved beside ", error_path)
}
r <- do.call(rbind, lapply(runs[!failed], function(x) matrix(x, nrow = 1L)))
method_names <- c("Oracle", "CICR", "E_FICR", "AA", "SC", "PPFL")
success_counts <- data.frame(method = method_names, requested = seed_all,
                             successful = 0L)
for (j in seq_along(method_names)) {
  valid <- vapply(r[, j], function(x) {
    metrics <- c("mse", "acc", "mse2", "ARI", "l", "ell")
    is.matrix(x) && nrow(x) > 0L && all(metrics %in% colnames(x)) &&
      all(is.finite(x[, metrics, drop = FALSE]))
  }, logical(1))
  success_counts$successful[j] <- sum(valid)
  if (any(!valid)) {
    errors <- rbind(errors, data.frame(
      seed = which(!failed)[!valid], method = method_names[j],
      message = "Missing or non-finite method result"
    ))
    r[!valid, j] <- rep(list(NULL), sum(!valid))
  }
}
success_counts$failed <- seed_all - success_counts$successful
if (nrow(errors)) write.csv(errors, error_path, row.names = FALSE)
if (any(success_counts$successful == 0L)) {
  save(runs, r, errors, success_counts,
       file = sub("[.]csv$", ".Rdata", error_path))
  stop("No valid S1_N results for: ",
       paste(success_counts$method[success_counts$successful == 0L], collapse = ", "),
       "; raw runs and diagnostics saved beside ", error_path)
}

if (nrow(r) > 0) {
  Oracle_res <- do.call(rbind, r[, 1])
  KR_res <- do.call(rbind, r[, 2])
  E_FICR_res <- do.call(rbind, r[, 3])
  AA_res<-do.call(rbind, r[, 4])
  SMA_res <- do.call(rbind, r[, 5])
  PPFL_res <- do.call(rbind, r[, 6])
}
t2<- Sys.time()
cat('time',t2 - t1, 'k', K, '\n')


method_results <- list(Oracle = Oracle_res, CICR = KR_res, E_FICR = E_FICR_res,
                       AA = AA_res, SC = SMA_res, PPFL = PPFL_res)
all <- t(vapply(method_results, function(x) {
  c(MSE_mean = mean(x[, "mse"]), MSE_SD = sd(x[, "mse"]),
    Acc_mean = mean(x[, "acc"]), Acc_SD = sd(x[, "acc"]))
}, numeric(4L)))
report_table <- data.frame(
  method = sub("^E_FICR$", "E-FICR", rownames(all)),
  failed = success_counts$failed[match(rownames(all), success_counts$method)],
  all,
  row.names = NULL
)
print_simulation_table(report_table, row.names = FALSE)

data_path <- file.path(file_path,sprintf("S1_N_n%d_m%d_gk%d_p%d_s%d_num%d_SS%d.Rdata", n, M, K, p,seed_all,num_init,SS))
save(runs, r, all, report_table, errors, error_counts, success_counts,
     file =data_path, n, M, K, p,seed_all,num_init,SS,iters)
