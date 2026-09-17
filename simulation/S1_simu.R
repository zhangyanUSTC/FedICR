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
  S1_main(beta0[1:K, , drop = FALSE], num_init = num_init, seed = i,
          iters = iters, K_hat = K, n = n, M = M, shred = 0.0001,
          balance = TRUE)
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
  "S1_M_errors_n%d_m%d_gk%d_p%d_s%d_num%d_SS%d.csv", n, M, K, p, seed_all, num_init, SS))
if (nrow(errors)) write.csv(errors, error_path, row.names = FALSE)
if (all(failed)) {
  save(runs, errors, n, M, K, p, seed_all, num_init, SS, iters,
       file = sub("[.]csv$", ".Rdata", error_path))
  stop("All repetitions failed; raw runs saved beside ", error_path)
}
r <- do.call(rbind, lapply(runs[!failed], function(x) matrix(x, nrow = 1L)))
method_names <- c("Oracle", "CICR", "E_FICR", "AA", "SC", "SC2")
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
  stop("No valid S1 results for: ",
       paste(success_counts$method[success_counts$successful == 0L], collapse = ", "),
       "; raw runs and diagnostics saved beside ", error_path)
}

if (nrow(r) > 0) {
  Oracle_res <- do.call(rbind, r[, 1])
  KR_res <- do.call(rbind, r[, 2])
  E_FICR_res <- do.call(rbind, r[, 3])
  AA_res<-do.call(rbind, r[, 4])
  SMA_res <- do.call(rbind, r[, 5])
  SMA_res2 <- do.call(rbind, r[, 6])
}
t2<- Sys.time()
cat('time',t2 - t1, 'k', K, '\n')


oracle_index <- colMeans(Oracle_res[, c(3, 4, 13:16), drop = FALSE])

res_list <- list(
  KR   = KR_res, 
  E_FICR  = E_FICR_res, 
  AA   = AA_res, 
  SMA  = SMA_res
)
iter_list <- lapply(res_list, function(res) {
  mean_index(res = res[res[, 2] == K, c(1:6, 15:18), drop = FALSE], iters)
})

# Transform the plotted loss gap only; retain raw losses.
plot_iter_list <- lapply(iter_list, function(it) {
  gap <- it[, 11] - tail(it[, 11], 1L)
  positive <- gap > 0
  it[, 11] <- NA_real_
  it[positive, 11] <- log(gap[positive])
  it
})

plots <- list()
name <- c("MSE","Acc",expression(MSE[2]),"ARI","l",expression(log(ell[t] - ell[t[max]])))
for (i in c(1L, 2L, 6L)) {
  current_iters <- ifelse(i == 6, min(20L, iters), iters)
  plots[[i]] <- plot_mse(
    iter_my  = plot_iter_list$E_FICR[1:current_iters, i+5], 
    iter_non = plot_iter_list$KR[1:current_iters, i+5],
    iter_sin = plot_iter_list$SMA[1:current_iters, i+5], 
    iter_ave = plot_iter_list$AA[1:current_iters, i+5],
    oracle   = if (i == 6L) numeric(0) else oracle_index[i],
    name     = name[i], 
    iters    = current_iters,
    na.rm    = i == 6L
  ) + theme(legend.position = "none")
}

legend <- simulation_legend(
  plot_mse(
    iter_my  = iter_list$E_FICR[, 6], 
    iter_non = iter_list$KR[, 6], 
    iter_sin = iter_list$SMA[, 6], 
    iter_ave = iter_list$AA[, 6],
    oracle   = oracle_index[1],
    name     = "Legend", iters = iters
  ) + theme(legend.position = "right"), position = "right"
)

combined_plot <- plot_grid(
  plot_grid(plots[[1]], plots[[2]], plots[[6]], ncol = 3),
  legend,
  nrow = 2,
  rel_heights = c(0.225, 0.1)
)

if (interactive()) print(combined_plot)

last_rows <- lapply(iter_list, function(x) tail(x, 1)[, 6:11])
result_table <- do.call(rbind, c(list(Oracle = oracle_index), last_rows))
col_append <- c(0, sapply(res_list, function(x) {
  mean(x[!duplicated(x[, "seed"]), "t_max"])
}))
result_table <- cbind(result_table, col_append)
result_table <- as.data.frame(result_table)
colnames(result_table) <- c("MSE", "Acc", "MSE2", "ARI", "l", "ell", "mean_t_max")
method_map <- c(Oracle = "Oracle", KR = "CICR", E_FICR = "E_FICR", AA = "AA", SMA = "SC")
report_methods <- unname(method_map[rownames(result_table)])
report_table <- data.frame(
  method = ifelse(report_methods == "E_FICR", "E-FICR", report_methods),
  failed = success_counts$failed[match(report_methods, success_counts$method)],
  result_table[, c("MSE", "Acc", "ell")],
  row.names = NULL
)
report_table$ell <- report_table$ell / (n * M)
print_simulation_table(report_table, row.names = FALSE)

fig_path <- file.path(file_path,sprintf("S1_M_n%d_m%d_gk%d_p%d_s%d_num%d_SS%d.png", n, M, K, p,seed_all,num_init,SS))
data_path <- file.path(file_path,sprintf("S1_M_RES_n%d_m%d_gk%d_p%d_s%d_num%d_SS%d.Rdata", n, M, K, p,seed_all,num_init,SS))
save(runs, res_list, Oracle_res, iter_list, plot_iter_list, oracle_index, result_table, report_table, errors, error_counts, success_counts,
     file = data_path, n, M, K, p, seed_all, num_init, SS,iters)

ggsave(filename = fig_path, plot = combined_plot, width = 8, height = 3.25, dpi = 300)
