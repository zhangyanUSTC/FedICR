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
           balance = TRUE, run_sc2 = FALSE, run_ppfl = FALSE)
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
  "S1_SS_errors_n%d_m%d_gk%d_p%d_s%d_num%d_SS%d.csv", n, M, K, p, seed_all, num_init, SS))
if (nrow(errors)) write.csv(errors, error_path, row.names = FALSE)
if (all(failed)) {
  save(runs, errors, n, M, K, p, seed_all, num_init, SS, iters,
       file = sub("[.]csv$", ".Rdata", error_path))
  stop("All repetitions failed; raw runs saved beside ", error_path)
}
r <- do.call(rbind, lapply(runs[!failed], function(x) matrix(x, nrow = 1L)))
method_names <- c("Oracle", "CICR", "E_FICR", "AA", "SC")
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
  stop("No valid S1_SS results for: ",
       paste(success_counts$method[success_counts$successful == 0L], collapse = ", "),
       "; raw runs and diagnostics saved beside ", error_path)
}

if (nrow(r) > 0) {
  Oracle_res <- do.call(rbind, r[, 1])
  KR_res <- do.call(rbind, r[, 2])
  E_FICR_res <- do.call(rbind, r[, 3])
  AA_res<-do.call(rbind, r[, 4])
  SMA_res <- do.call(rbind, r[, 5])
}
t2<- Sys.time()
cat('time',t2 - t1, 'k', K, '\n')


my_mean<-c(colMeans(E_FICR_res[,c(5,6,18),drop=FALSE]),colSds(E_FICR_res[,c(5,6,18),drop=FALSE]),colMedians(E_FICR_res[,c(5,6,18),drop=FALSE]))
non_mean<-c(colMeans(KR_res[,c(5,6,18),drop=FALSE]),colSds(KR_res[,c(5,6,18),drop=FALSE]),colMedians(KR_res[,c(5,6,18),drop=FALSE]))
sin_mean<-c(colMeans(SMA_res[,c(5,6,18),drop=FALSE]),colSds(SMA_res[,c(5,6,18),drop=FALSE]),colMedians(SMA_res[,c(5,6,18),drop=FALSE]))
ave_mean<-c(colMeans(AA_res[,c(5,6,18),drop=FALSE]),colSds(AA_res[,c(5,6,18),drop=FALSE]),colMedians(AA_res[,c(5,6,18),drop=FALSE]))

ora_mean<-c(colMeans(Oracle_res[,c(3,4,16),drop=FALSE]),colSds(Oracle_res[,c(3,4,16),drop=FALSE]),colMedians(Oracle_res[,c(3,4,16),drop=FALSE]))

all<-rbind(ora_mean,non_mean,my_mean,ave_mean,sin_mean)
rownames(all) <- c("Oracle", "CICR", "E_FICR", "AA", "SC")
colnames(all) <- c("MSE_mean", "Acc_mean", "ell_mean", "MSE_SD", "Acc_SD", "ell_SD",
                   "MSE_median", "Acc_median", "ell_median")
report_table <- data.frame(
  method = sub("^E_FICR$", "E-FICR", rownames(all)),
  failed = success_counts$failed[match(rownames(all), success_counts$method)],
  all[, c("MSE_mean", "MSE_SD", "Acc_mean", "Acc_SD", "ell_mean", "ell_SD")],
  row.names = NULL
)
print_simulation_table(report_table, row.names = FALSE)
all<-cbind(all,n,M,p,K,SS)

data_path <- file.path(file_path,sprintf("S1_SS_n%d_m%d_gk%d_p%d_s%d_num%d_SS%d.Rdata", n, M, K, p,seed_all,num_init,SS))
save(runs, r, all, report_table, errors, error_counts, success_counts,
     file =data_path, n, M, K, p,seed_all,num_init,SS,iters)
