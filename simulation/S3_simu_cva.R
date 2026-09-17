
.s3_cva_source_files <- Filter(
  function(x) is.character(x) && length(x) == 1L,
  lapply(sys.frames(), function(frame) frame$ofile)
)
if (length(.s3_cva_source_files) == 0L) {
  stop("S3_simu_cva.R must be sourced from the simulation directory")
}
.s3_cva_code_dir <- dirname(
  normalizePath(tail(.s3_cva_source_files, 1L)[[1L]])
)
rm(.s3_cva_source_files)

run_s3_cva_batch <- function(beta0, n, M, K, p, num_init, iters,
                             seed_all, SS, K_choose, file_path,
                             code_dir = NULL) {
  stopifnot(
    is.function(S3_main),
    all(dim(beta0) == c(K, p)),
    length(seed_all) == 1L, is.finite(seed_all), seed_all >= 1L,
    seed_all == as.integer(seed_all),
    num_init >= 1L, num_init == as.integer(num_init),
    iters >= 2L, iters == as.integer(iters),
    length(K_choose) >= 1L,
    all(is.finite(K_choose)), all(K_choose == as.integer(K_choose)),
    !anyDuplicated(K_choose),
    K %in% K_choose
  )
  
  if (is.null(code_dir)) {
    source_files <- Filter(
      function(x) is.character(x) && length(x) == 1L,
      lapply(sys.frames(), function(frame) frame$ofile)
    )
    if (length(source_files) == 0L) {
      stop("Please provide code_dir, the directory containing S0_basic.R")
    }
    code_dir <- dirname(normalizePath(tail(source_files, 1L)[[1L]]))
  }
  
  dir.create(file_path, recursive = TRUE, showWarnings = FALSE)
  num_cores <- min(seed_all, getOption("E_FICR.workers", 10L))
  
  cl <- start_workers(num_cores)
  on.exit({
    parallel::stopCluster(cl)
    foreach::registerDoSEQ()
  }, add = TRUE)
  doParallel::registerDoParallel(cl)
  
  parallel::clusterCall(
    cl,
    function(path) {
      suppressPackageStartupMessages(source(path, local = .GlobalEnv))
      NULL
    },
    file.path(code_dir, "S0_basic.R")
  )
  
  fit_args <- list(
    beta0 = beta0,
    num_init = num_init,
    iters = iters,
    n = n,
    M = M,
    shred = 0.0001,
    balance = TRUE,
    K_choose = K_choose,
    evaluate_selected = TRUE
  )
  
  t1 <- Sys.time()
  runs <- foreach::foreach(
    i = seq_len(seed_all),
    .packages = c("clue", "MASS", "combinat", "mclust", "matrixStats",
                  "Matrix", "fossil"),
    .errorhandling = "pass"
  ) %dopar% {
    do.call(S3_main, c(fit_args, list(seed = i)))
  }
  elapsed <- Sys.time() - t1
  
  is_error <- vapply(runs, inherits, logical(1), what = "error")
  errors <- data.frame(
    seed = which(is_error),
    message = vapply(
      runs[is_error], conditionMessage, character(1)
    ),
    row.names = NULL
  )
  successful_runs <- runs[!is_error]
  stem <- sprintf(
    "S3_n%d_m%d_gk%d_p%d_s%d_num%d_SS%d",
    n, M, K, p, seed_all, num_init, SS
  )
  error_path <- file.path(file_path, paste0("errors_", stem, ".csv"))
  if (nrow(errors)) utils::write.csv(errors, error_path, row.names = FALSE)
  if (length(successful_runs) == 0L) {
    save(runs, errors, beta0, n, M, K, p, seed_all, num_init, iters, SS, K_choose,
         file = file.path(file_path, paste0(stem, ".Rdata")))
    stop("All Simulation 3 repetitions failed; see ", error_path)
  }
  
  method_map <- c(CICR = "KR", E_FICR = "E_FICR", AA = "AA", SC = "SMA")
  result_matrices <- lapply(method_map, function(result_name) {
    do.call(rbind, lapply(successful_runs, function(run) {
      run$cva[[result_name]]
    }))
  })
  selected_metrics <- do.call(
    rbind,
    lapply(successful_runs, `[[`, "selected_metrics")
  )
  
  score_columns <- paste0("K", K_choose)
  select_k <- function(result_matrix) {
    scores <- result_matrix[, score_columns, drop = FALSE]
    apply(scores, 1L, function(x) {
      if (!all(is.finite(x))) return(NA_integer_)
      K_choose[which.min(x)]
    })
  }
  min_s <- function(result_matrix) {
    scores <- result_matrix[, score_columns, drop = FALSE]
    apply(scores, 1L, function(x) {
      if (!all(is.finite(x))) return(NA_real_)
      min(x)
    })
  }
  
  selected_k <- lapply(result_matrices, select_k)
  min_scores <- lapply(result_matrices, min_s)
  successful_seeds <- as.integer(result_matrices[[1L]][, "seed"])
  failed_candidates <- do.call(rbind, lapply(names(result_matrices), function(method) {
    mat <- result_matrices[[method]]
    ids <- which(!is.finite(mat[, score_columns, drop = FALSE]), arr.ind = TRUE)
    data.frame(seed = as.integer(mat[ids[, 1L], "seed"]),
               method = rep(method, nrow(ids)), K = K_choose[ids[, 2L]],
               row.names = NULL)
  }))
  
  K_out <- data.frame(
    seed = successful_seeds,
    as.data.frame(selected_k, check.names = FALSE),
    check.names = FALSE
  )
  s_out <- t(vapply(names(method_map), function(method) {
    selected <- selected_k[[method]]
    best_score <- min_scores[[method]]
    c(
      Correct_K_pct = if (any(is.finite(selected))) mean(selected == K, na.rm = TRUE) * 100 else NA_real_,
      Mean_min_s = if (any(is.finite(best_score))) mean(best_score, na.rm = TRUE) else NA_real_,
      SD_min_s = stats::sd(best_score, na.rm = TRUE),
      Successful_N = sum(is.finite(selected)),
      Requested_N = seed_all,
      Incomplete_N = seed_all - sum(is.finite(selected))
    )
  }, numeric(6L)))
  
  finite_stat <- function(x, fun) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) return(NA_real_)
    fun(x)
  }
  selected_summary <- do.call(rbind, lapply(
    c("Oracle", names(method_map)),
    function(method) {
      dat <- selected_metrics[selected_metrics$method == method, , drop = FALSE]
      data.frame(
        Method = method,
        MSE2_mean = finite_stat(dat$mse2, mean),
        MSE2_SD = finite_stat(dat$mse2, stats::sd),
        MSE2_median = finite_stat(dat$mse2, stats::median),
        ARI_mean = finite_stat(dat$ARI, mean),
        ARI_SD = finite_stat(dat$ARI, stats::sd),
        ARI_median = finite_stat(dat$ARI, stats::median),
        Successful_N = sum(is.finite(dat$mse2) & is.finite(dat$ARI)),
        row.names = NULL
      )
    }
  ))
  
  rdata_path <- file.path(file_path, paste0(stem, ".Rdata"))
  
  save(
    runs, result_matrices, K_out, s_out, selected_metrics,
    selected_summary, errors, failed_candidates,
    beta0, n, M, K, p, seed_all, num_init, iters, SS, K_choose,
    file = rdata_path
  )
  
  cat("Simulation 3 completed.\n")
  cat("Elapsed time:", format(elapsed), "\n")
  cat("Successful repetitions:", length(successful_runs), "of", seed_all, "\n")
  cat("Failed repetitions:", nrow(errors), "\n")
  cat("Failed candidate fits:", nrow(failed_candidates), "\n")
  cat("CVA s is normalized pairwise grouping instability (s_n in the paper); smaller is better.\n")
  cat("Mean_min_s and SD_min_s summarize s at each repetition's selected K.\n")
  cat("K-selection summaries use repetitions with finite scores for every candidate K.\n")
  cat("\nCVA selection summary:\n")
  print_simulation_table(s_out, pad = FALSE)
  cat("\nMSE2 and ARI after refitting at the selected K:\n")
  print_simulation_table(selected_summary, row.names = FALSE, pad = FALSE)
  cat("RData file:", normalizePath(rdata_path, mustWork = FALSE), "\n")
  
  invisible(list(
    runs = runs,
    result_matrices = result_matrices,
    K_out = K_out,
    s_out = s_out,
    selected_metrics = selected_metrics,
    selected_summary = selected_summary,
    errors = errors,
    failed_candidates = failed_candidates,
    paths = c(RData = rdata_path)
  ))
}

.s3_cva_batch <- run_s3_cva_batch(
  beta0 = beta0,
  n = n,
  M = M,
  K = K,
  p = p,
  num_init = num_init,
  iters = iters,
  seed_all = seed_all,
  SS = SS,
  K_choose = K_choose,
  file_path = file_path,
  code_dir = .s3_cva_code_dir
)
list2env(.s3_cva_batch, envir = environment())
rm(.s3_cva_batch, .s3_cva_code_dir)
