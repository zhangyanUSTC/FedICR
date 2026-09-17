
main <- function() {
  data_file <- file.path(getOption("E_FICR.data"), "pisa_federated_data_pv.rds")
  function_file <- file.path(getOption("E_FICR.root"), "pisa", "pisa_federated.R")

  K_hat <- 2L
  iters <- 40L
  shred <- 1e-4
  lam <- 0.001
  init_seeds <- 1001:1005

  output_dir <- getOption("E_FICR.output")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  output_file <- file.path(
    output_dir,
    paste0(
      "Pisa_PV_sensitivity_",
      format(Sys.time(), "%Y-%m-%d_%H-%M-%S"),
      ".txt"
    )
  )

  log_connection <- file(output_file, open = "wt")
  sink(log_connection, split = TRUE)
  sink(log_connection, type = "message")
  on.exit({
    sink(type = "message")
    while (sink.number() > 0L) sink()
    close(log_connection)
  }, add = TRUE)

  cat("PISA CICR plausible-value sensitivity analysis\n")
  cat("Started:", format(Sys.time()), "\n")
  cat("Fixed K:", K_hat, "\n")
  cat("Iterations:", iters, "\n")
  cat("Initialization seeds:", paste(init_seeds, collapse = ", "), "\n\n")

  if (!file.exists(data_file)) stop("Data file not found: ", data_file)
  if (!file.exists(function_file)) stop("Function file not found: ", function_file)

  suppressPackageStartupMessages(
    suppressMessages(
      source(function_file)
    )
  )
  federated_data <- readRDS(data_file)

  required_background <- c("1", "CNT", "ANXMAT", "MISCED", "FISCED")
  required_pv <- unlist(lapply(1:10, function(j) {
    paste0("PV", j, c("MATH", "READ", "SCIE"))
  }))
  required_columns <- c(required_background, required_pv)

  invalid_clients <- names(federated_data)[
    !vapply(
      federated_data,
      function(x) all(required_columns %in% names(x)),
      logical(1)
    )
  ]
  if (length(invalid_clients) > 0L) {
    stop("Required columns are missing for: ", paste(invalid_clients, collapse = ", "))
  }

  client_sizes <- vapply(federated_data, nrow, integer(1))
  cat("Countries/economies:", length(federated_data), "\n")
  cat("Students:", sum(client_sizes), "\n\n")

  results <- vector("list", 10L)
  coefficient_names <- c(
    "Intercept", "ANXMAT", "MISCED", "FISCED", "READ", "SCIE"
  )

  for (j in 1:10) {
    cat("============================================================\n")
    cat("Fitting matched PV set", j, "of 10\n")

    response_var <- paste0("PV", j, "MATH")
    predictor_vars <- c(
      "1", "ANXMAT", "MISCED", "FISCED",
      paste0("PV", j, "READ"),
      paste0("PV", j, "SCIE")
    )

    X <- lapply(federated_data, function(client_data) {
      as.matrix(client_data[, predictor_vars, drop = FALSE])
    })
    Y <- lapply(federated_data, function(client_data) {
      as.numeric(client_data[[response_var]])
    })
    names(X) <- names(federated_data)
    names(Y) <- names(federated_data)

    fits <- vector("list", length(init_seeds))
    losses <- rep(Inf, length(init_seeds))

    for (s in seq_along(init_seeds)) {
      seed <- init_seeds[s]
      cat("  initialization seed", seed, "...")

      fit <- tryCatch({
        set.seed(seed)
        beta_init <- INIT_beta_M(
          X = X,
          Y = Y,
          K_hat = K_hat,
          lam = lam
        )
        CICR_alg(
          X = X,
          Y = Y,
          beta_init = beta_init,
          iters = iters,
          shred = shred
        )
      }, error = function(e) {
        cat(" failed:", conditionMessage(e), "\n")
        NULL
      })

      if (!is.null(fit)) {
        loss <- as.numeric(fit$out[K_hat + 2L])
        if (is.finite(loss)) {
          fits[[s]] <- fit
          losses[s] <- loss
          cat(" loss =", format(loss, digits = 8), "\n")
        } else {
          cat(" failed: non-finite loss\n")
        }
      }
    }

    if (!any(is.finite(losses))) {
      stop("All initializations failed for PV", j, ".")
    }

    best_index <- which.min(losses)
    best_fit <- fits[[best_index]]
    z_hat <- update_Z_list(X = X, Y = Y, beta_hat = best_fit$beta_hat)
    labels <- unlist(
      lapply(z_hat, max.col, ties.method = "first"),
      use.names = FALSE
    )

    results[[j]] <- list(
      pv = j,
      beta_hat = best_fit$beta_hat,
      labels = labels,
      loss = losses[best_index],
      seed = init_seeds[best_index]
    )

    cat("Selected seed:", results[[j]]$seed, "\n")
    cat("Selected loss:", format(results[[j]]$loss, digits = 8), "\n")
    cat("Coefficients before label alignment:\n")
    for (k in 1:K_hat) {
      cat(
        paste(formatC(results[[j]]$beta_hat[k, ], format = "f", digits = 5),
                    collapse = " & "),
        "\n"
      )
    }
    cat("\n")

    rm(X, Y, fits, z_hat, best_fit)
    gc(verbose = FALSE)
  }

  # Group 1 has the higher PV1 mean fitted value.
  pv1_predictors <- c(
    "1", "ANXMAT", "MISCED", "FISCED", "PV1READ", "PV1SCIE"
  )
  X_pv1 <- lapply(federated_data, function(client_data) {
    as.matrix(client_data[, pv1_predictors, drop = FALSE])
  })
  total_n <- sum(vapply(X_pv1, nrow, integer(1)))
  pv1_fitted_means <- vapply(1:K_hat, function(k) {
    sum(vapply(X_pv1, function(X_m) {
      sum(drop(X_m %*% results[[1]]$beta_hat[k, ]))
    }, numeric(1))) / total_n
  }, numeric(1))
  if (pv1_fitted_means[1] < pv1_fitted_means[2]) {
    results[[1]]$labels <- 3L - results[[1]]$labels
    results[[1]]$beta_hat <- results[[1]]$beta_hat[c(2, 1), , drop = FALSE]
  }

  # Align all other PV solutions to the PV1 labels.
  reference_labels <- results[[1]]$labels
  for (j in 2:10) {
    original <- results[[j]]$labels
    swapped <- 3L - original
    if (mean(swapped == reference_labels) > mean(original == reference_labels)) {
      results[[j]]$labels <- swapped
      results[[j]]$beta_hat <- results[[j]]$beta_hat[c(2, 1), , drop = FALSE]
    }
  }

  group_summary <- do.call(rbind, lapply(1:10, function(j) {
    group_n <- tabulate(results[[j]]$labels, nbins = K_hat)
    data.frame(
      PV = j,
      Seed = results[[j]]$seed,
      Loss = results[[j]]$loss,
      Group1_n = group_n[1],
      Group2_n = group_n[2],
      Group1_prop = group_n[1] / sum(group_n),
      Group2_prop = group_n[2] / sum(group_n),
      Acc_with_PV1 = if (j == 1L) 100 else {
        mean(results[[j]]$labels == reference_labels) * 100
      }
    )
  }))

  label_matrix <- do.call(
    cbind,
    lapply(results, function(result) result$labels)
  )
  colnames(label_matrix) <- paste0("PV", 1:10)

  acc_matrix <- diag(100, nrow = 10L, ncol = 10L)
  dimnames(acc_matrix) <- list(paste0("PV", 1:10), paste0("PV", 1:10))

  pv_pairs <- utils::combn(1:10, 2)
  pairwise_summary <- data.frame(
    PV_A = paste0("PV", pv_pairs[1, ]),
    PV_B = paste0("PV", pv_pairs[2, ]),
    Acc_percent = NA_real_
  )

  for (pair_index in seq_len(ncol(pv_pairs))) {
    a <- pv_pairs[1, pair_index]
    b <- pv_pairs[2, pair_index]
    labels_a <- label_matrix[, a]
    labels_b <- label_matrix[, b]

    pair_acc <- max(
      mean(labels_a == labels_b),
      mean(labels_a == (3L - labels_b))
    ) * 100

    acc_matrix[a, b] <- acc_matrix[b, a] <- pair_acc
    pairwise_summary$Acc_percent[pair_index] <- pair_acc
  }

  beta_array <- array(
    NA_real_,
    dim = c(10L, K_hat, length(coefficient_names)),
    dimnames = list(
      PV = paste0("PV", 1:10),
      Group = paste0("Group", 1:K_hat),
      Variable = coefficient_names
    )
  )
  for (j in 1:10) beta_array[j, , ] <- results[[j]]$beta_hat

  coefficient_summary <- do.call(rbind, lapply(1:K_hat, function(k) {
    do.call(rbind, lapply(seq_along(coefficient_names), function(v) {
      values <- beta_array[, k, v]
      data.frame(
        Group = k,
        Variable = coefficient_names[v],
        PV1 = values[1],
        Mean = mean(values),
        Min = min(values),
        Max = max(values),
        SD = stats::sd(values),
        Sign = if (all(values > 0)) "Positive" else if (all(values < 0)) {
          "Negative"
        } else {
          "Mixed"
        }
      )
    }))
  }))

  difference_summary <- do.call(rbind, lapply(seq_along(coefficient_names), function(v) {
    values <- beta_array[, 2, v] - beta_array[, 1, v]
    data.frame(
      Variable = coefficient_names[v],
      PV1_difference = values[1],
      Mean_difference = mean(values),
      Min_difference = min(values),
      Max_difference = max(values),
      SD_difference = stats::sd(values),
      Direction = if (all(values > 0)) "Group2 > Group1" else if (all(values < 0)) {
        "Group2 < Group1"
      } else {
        "Mixed"
      }
    )
  }))

  cat("\n\n================ GROUP SUMMARY ================\n")
  cat("Acc_with_PV1 is the percentage of aligned student labels matching PV1.\n")
  print(group_summary, row.names = FALSE, digits = 6)

  cat("\n================ PAIRWISE Acc MATRIX (%) ================\n")
  cat("Acc uses the optimal label permutation for each PV pair.\n")
  cat("The diagonal is 100%. All 45 unique PV pairs are included.\n")
  cat("PV & ", paste(colnames(acc_matrix), collapse = " & "), "\n", sep = "")
  for (i in 1:10) {
    cat(
      rownames(acc_matrix)[i], " & ",
      paste(formatC(acc_matrix[i, ], format = "f", digits = 4), collapse = " & "),
      "\n", sep = ""
    )
  }

  cat("\n================ ALL 45 PAIRWISE COMPARISONS ================\n")
  print(pairwise_summary, row.names = FALSE, digits = 6)
  cat("Mean pairwise Acc (%):", format(mean(pairwise_summary$Acc_percent), digits = 6), "\n")
  cat("Minimum pairwise Acc (%):", format(min(pairwise_summary$Acc_percent), digits = 6), "\n")

  cat("\n================ ALIGNED COEFFICIENTS ================\n")
  cat("Columns:", paste(coefficient_names, collapse = " & "), "\n")
  for (j in 1:10) {
    cat("PV", j, "\n", sep = "")
    for (k in 1:K_hat) {
      cat(
        "Group", k, " & ",
        paste(formatC(beta_array[j, k, ], format = "f", digits = 5),
              collapse = " & "),
        "\n",
        sep = ""
      )
    }
  }

  cat("\n================ COEFFICIENT SENSITIVITY ================\n")
  print(coefficient_summary, row.names = FALSE, digits = 6)

  cat("\n================ CROSS-GROUP DIFFERENCES ================\n")
  cat("Difference is Group 2 minus Group 1.\n")
  print(difference_summary, row.names = FALSE, digits = 6)

  cat("Maximum pairwise Acc (%):", max(pairwise_summary$Acc_percent), "\n")
  save(results, beta_array, label_matrix, acc_matrix, pairwise_summary,
       group_summary, coefficient_summary, difference_summary, client_sizes,
       K_hat, iters, shred, lam, init_seeds,
       file = sub("[.]txt$", ".Rdata", output_file))
  cat("\nFinished:", format(Sys.time()), "\n")
  cat("Results file:", normalizePath(output_file, mustWork = FALSE), "\n")

  invisible(output_file)
}

main()
