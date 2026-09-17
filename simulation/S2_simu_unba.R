num_cores <- min(seed_all, getOption("E_FICR.workers", 10L))
cl <- start_workers(num_cores)
registerDoParallel(cl)

Oracle_res <- NULL;KR_res <- NULL;E_FICR_res <- NULL
AA_res<-NULL;SMA_res<-NULL;SMA_res2<-NULL;PPFL_res<-NULL
error_counts <- list()

for (kk in K:K) {
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
    S2_main(beta0[1:K, , drop = FALSE], num_init = num_init, seed = i,
            iters = iters, K_hat = kk, n = n, M = M, shred = 0.0001,
            balance = FALSE)
  }, finally = {
    stopCluster(cl)
    foreach::registerDoSEQ()
  })
  failed <- vapply(runs, inherits, logical(1), what = "error")
  errors <- data.frame(seed = which(failed), method = rep("seed", sum(failed)),
                       init = rep(NA_integer_, sum(failed)),
                       message = vapply(runs[failed], conditionMessage, character(1)))
  method_errors <- do.call(rbind, lapply(runs[!failed], `[[`, "errors"))
  if (!is.null(method_errors)) errors <- rbind(errors, method_errors)
  error_count <- sum(failed)
  error_counts[[as.character(kk)]] <- error_count
  error_path <- file.path(file_path, sprintf(
    "S2_errors_n%d_m%d_gk%d_p%d_s%d_num%d_SS%d.csv", n, M, K, p, seed_all, num_init, SS))
  if (nrow(errors)) write.csv(errors, error_path, row.names = FALSE)
  if (all(failed)) {
  save(runs, errors, n, M, K, p, seed_all, num_init, SS, iters,
       file = sub("[.]csv$", ".Rdata", error_path))
  stop("All repetitions failed; raw runs saved beside ", error_path)
}
  r <- do.call(rbind, lapply(runs[!failed], function(x) matrix(x, nrow = 1L)))
  
  if (nrow(r) > 0) {
    Oracle_res <- rbind(Oracle_res, do.call(rbind, r[, 1]))
    KR_res <- rbind(KR_res, do.call(rbind, r[, 2]))
    E_FICR_res <- rbind(E_FICR_res, do.call(rbind, r[, 3]))
    AA_res<-rbind(AA_res, do.call(rbind, r[, 4]))
    SMA_res <- rbind(SMA_res, do.call(rbind, r[, 5]))
    SMA_res2 <- rbind(SMA_res2, do.call(rbind, r[, 6]))
    PPFL_res <- rbind(PPFL_res, do.call(rbind, r[, 7]))
  }
  t2<- Sys.time()
  cat('time',t2 - t1, 'k', kk, '\n')
}

cat("有效结果数（各方法分别统计）：\n")
success_counts <- vapply(list(Oracle = Oracle_res, CICR = KR_res, E_FICR = E_FICR_res,
                              AA = AA_res, SC = SMA_res, SC2 = SMA_res2, PPFL = PPFL_res),
                         function(x) sum(is.finite(x[, "mse"])), integer(1))
print_simulation_table(success_counts, pad = FALSE)

summarize_groups <- function(x) {
  x <- x[, paste0("mse_g", 1:4), drop = FALSE]
  ans <- c(colMeans(x, na.rm = TRUE), matrixStats::colSds(x, na.rm = TRUE))
  ans[!is.finite(ans)] <- NA_real_
  ans
}
last_KR <- summarize_groups(KR_res)
last_E_FICR <- summarize_groups(E_FICR_res)
last_AA <- summarize_groups(AA_res)
last_SMA <- summarize_groups(SMA_res)
last_SMA2 <- summarize_groups(SMA_res2)
last_oracle <- summarize_groups(Oracle_res)
last_PPFL <- summarize_groups(PPFL_res)
result_table <- rbind( last_oracle, last_KR,last_E_FICR,last_AA, last_SMA, last_SMA2,last_PPFL)

result_table <- as.data.frame(result_table)
result_table<-round(result_table,3)

data_path <- file.path(file_path,sprintf("S2_n%d_m%d_gk%d_p%d_s%d_num%d_SS%d.Rdata", n, M, K, p,seed_all,num_init,SS))
save(runs, r, result_table, errors, error_counts, success_counts,
     file =data_path, n, M, K, p,seed_all,num_init,SS,iters)

print_simulation_table(result_table, pad = FALSE)

fig_path2 <- file.path(file_path,sprintf("S2_n%d_m%d_gk%d_p%d_s%d_num%d_result.png", n, M, K, p,seed_all,num_init))
png(filename = fig_path2, width = 900, height = 300)
formatted_table <- result_table
formatted_table[] <- lapply(result_table, function(x) format(x, digits = 4, nsmall = 4))
grid.table(formatted_table)
dev.off()
