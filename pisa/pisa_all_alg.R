current_dir <- getwd()
outcome_dir <- getOption("E_FICR.output")
if (!dir.exists(outcome_dir)) {
  dir.create(outcome_dir)
  message("Directory 'outcome' created.")
} else {
  message("Directory 'outcome' already exists.")
}
time_stamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
out_path <- file.path(outcome_dir, paste0("Pisa_All_Algorithms_", time_stamp, ".txt"))
sink(out_path, split = TRUE)

section <- function(title) {
  cat("\n============================================================\n")
  cat(title, "\n")
  cat("============================================================\n")
}

print_beta_amp <- function(beta, digits = 5L, significant = FALSE) {
  beta <- as.matrix(beta)
  if (is.null(colnames(beta))) colnames(beta) <- paste0("Beta", seq_len(ncol(beta)))
  cat("Columns:", paste(colnames(beta), collapse = " & "), "\n")
  for (k in seq_len(nrow(beta))) {
    values <- if (significant) {
      format(beta[k, ], digits = digits, trim = TRUE)
    } else {
      formatC(beta[k, ], format = "f", digits = digits)
    }
    cat(
      "Group", k, " & ",
      paste(values, collapse = " & "),
      "\n", sep = ""
    )
  }
}

print_named_values <- function(values, digits = 5L) {
  for (nm in names(values)) {
    value_text <- if (tolower(nm) %in% c("mse", "acc", "accuracy")) {
      formatC(values[[nm]], format = "f", digits = 5L)
    } else {
      format(values[[nm]], digits = digits, trim = TRUE)
    }
    cat(nm, " & ", value_text, "\n", sep = "")
  }
}

# Group 1 has the higher full-sample mean fitted value.
align_groups_by_fitted_mean <- function(beta, X) {
  beta <- as.matrix(beta)
  total_n <- sum(vapply(X, nrow, integer(1)))
  fitted_means <- vapply(seq_len(nrow(beta)), function(k) {
    sum(vapply(
      X,
      function(X_m) sum(drop(X_m %*% beta[k, ])),
      numeric(1)
    )) / total_n
  }, numeric(1))
  group_order <- order(fitted_means, decreasing = TRUE)
  list(
    beta = beta[group_order, , drop = FALSE],
    fitted_means = fitted_means[group_order],
    original_order = group_order
  )
}

print_group_sizes <- function(z, fitted_means) {
  K <- length(fitted_means)
  group_sizes <- rowSums(vapply(
    z,
    function(z_m) colSums(z_m),
    numeric(K)
  ))
  total_n <- sum(group_sizes)

  cat("CICR group sizes aligned with the coefficient rows:\n")
  cat("Group & Students & Proportion & Mean_fitted_value\n")
  for (k in seq_len(K)) {
    cat(
      "Group", k, " & ", group_sizes[k], " & ",
      formatC(group_sizes[k] / total_n, format = "f", digits = 5L), " & ",
      formatC(fitted_means[k], format = "f", digits = 5L),
      "\n", sep = ""
    )
  }
  invisible(group_sizes)
}

iters<-40;shred<-0.0001;K_hat<-2;num_init<-5
cat("PISA comparison of CICR, E-FICR, AA and SC\n")
cat("Started:", format(Sys.time()), "\n")
cat("Fixed K:", K_hat, "\n")
cat("Initializations:", num_init, "\n")
cat("Iterations:", iters, "\n")
cat("Convergence threshold:", shred, "\n")
federated_data <- readRDS(
  file.path(getOption("E_FICR.data"), "pisa_federated_data_pv1.rds")
)
source(file.path(getOption("E_FICR.root"), "pisa", "pisa_federated.R"))
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

section("CICR")
cat("Fitting CICR with multiple initializations:\n")
DR_out<-NULL;ini_choose<-NULL
for (ini in 1:5) {
  set.seed(ini)
  beta_init<-INIT_beta_M(X,Y,K_hat=2,lam=0.001)
  d<-CICR_alg(X,Y,beta_init=beta_init,iters=iters,shred=0.0001)
  cat(ini,"*")
  if(is.null(DR_out)){
    DR_out<-d;ini_choose<-ini
  }
  if(!is.null(DR_out)& !is.null(d)){
    if(DR_out$out[4]>d$out[4])
    {DR_out<-d;ini_choose<-ini}
  }
}
cat("\nSelected initialization:", ini_choose, "\n")

cicr<-DR_out
cat("Loss:", format(cicr$out[K_hat + 2L], digits = 8), "\n")
group_alignment <- align_groups_by_fitted_mean(cicr$beta_hat, X)
cicr$beta_hat <- group_alignment$beta
cat(
  "Group alignment: decreasing full-sample mean fitted value; original rows =",
  paste(group_alignment$original_order, collapse = ", "), "\n"
)
colnames(cicr$beta_hat) <- colnames(X[[1]])
cat("CICR coefficients:\n")
print_beta_amp(cicr$beta_hat, digits = 2L, significant = TRUE)

beta0<-cicr$beta_hat
z0<- update_Z_list(
  X = X,
  Y = Y,
  beta_hat = beta0)
cicr_group_sizes <- print_group_sizes(z0, group_alignment$fitted_means)

beta_init<-NULL
for (ini in 1:num_init) {
  beta_init<-rbind(beta_init,INIT_beta_M(X,Y,K_hat=2,lam=0.001))
}

section("E-FICR")
cat("Fitting E-FICR with multiple initializations:\n")
E_FICR_out<-NULL
for (ini in 1:num_init) {
  a<-E_FICR_alg(X,Y,beta_init=beta_init[((ini-1)*K_hat+1):(ini*K_hat),],z0=z0,beta0=beta0,iters=iters,shred=shred)
  a<-t(as.matrix(c(ini,K_hat,a$t,a$out,ini)))
  if(is.null(E_FICR_out)){
    E_FICR_out<-a
  }
  if(!is.null(E_FICR_out)& ncol(a)>3){
    if(ncol(E_FICR_out)==3){
      E_FICR_out<-a
    } else {
      if(E_FICR_out[nrow(E_FICR_out),7]>a[nrow(a),7])  E_FICR_out<-a
    }
  }
  cat(ini,"*")
}
cat("\nSelected E-FICR result:\n")
E_FICR_values <- as.numeric(E_FICR_out[nrow(E_FICR_out), ])
E_FICR_labels <- c("Initialization", "K", "Iterations", "Reported_iterations",
                "MSE", "Accuracy", "Loss", "Initialization_copy")
names(E_FICR_values) <- E_FICR_labels[seq_along(E_FICR_values)]
print_named_values(E_FICR_values)

section("AA")
cat("Fitting local models and aggregating coefficients:\n")
num_cores <- getOption("E_FICR.workers", 10L)
cl <- start_workers(getOption("E_FICR.workers", 10L))
registerDoParallel(cl)
t1 <- Sys.time()
error_count <- 0
r <- foreach(
  m = 1:M,
  .combine = rbind,
  .packages = c(
    'clue',
    "MASS",
    "foreach",
    "combinat",
    "mclust",
    "Matrix"
  ),
  .errorhandling = "remove"
) %dopar% {
  tryCatch(
    {
      SMA_out<-NULL
      for (ini in 1:num_init) {
        d<-SC_alg(X,Y,beta_init=beta_init[((ini-1)*K_hat+1):(ini*K_hat),],z0=z0,beta0=beta0,iters=iters,shred=shred,single_m = m)
        beta_m<-d$beta_hat
        d<-d$out
        
        if(is.null(SMA_out)& !is.null(d)){
          SMA_out<-d;beta_Mm<-beta_m
        }
        if(!is.null(SMA_out)& !is.null(d)){
          if(SMA_out[4]>d[4]) {
            SMA_out<-d;beta_Mm<-beta_m
          }
        }
      }
      beta_Mm
    },
    error = function(e) {
      NULL
    }
  )
}
t2<- Sys.time()
cat('\ntime',t2 - t1, 'k', K_hat, '\n')
stopCluster(cl)

set.seed(1)
beta_hat<-kmeans(r,centers = K_hat,nstart = 20)$centers
colnames(beta_hat) <- colnames(X[[1]])
z_hat <- update_Z_list(X = X,Y = Y,beta_hat = beta_hat)
aa_score <- c(Iterations = iters, SCORE_hat(X,Y,beta0,beta_hat,z0,z_hat))
cat("AA metrics:\n")
print_named_values(aa_score)
cat("AA coefficients:\n")
print_beta_amp(beta_hat)

section("SC")
cat("Fitting each client separately:\n")
client_ids <- seq_len(min(M, length(X)))
n_cores <- getOption("E_FICR.workers", 10L)
cl <- start_workers(getOption("E_FICR.workers", 10L))
clusterExport(
  cl,
  varlist = c(
    "X", "Y", "beta_init", "z0", "beta0","iters","shred","num_init","K_hat",
    "SC_alg", "SCORE_hat", "l_calculate","update_Z_list"
  ),
  envir = environment()
)
results <- parLapply(cl, client_ids, function(m) {
  SMA_out<-NULL
  for (ini in 1:num_init) {
    d<-SC_alg(X,Y,beta_init=beta_init[((ini-1)*K_hat+1):(ini*K_hat),],z0=z0,beta0=beta0,iters=iters,shred=shred,single_m = m)
    d<-t(as.matrix(c(ini,K_hat,d$t,d$out,ini)))
    if(is.null(SMA_out)){
      SMA_out<-d
    }
    if(!is.null(SMA_out)& ncol(d)>3){
      if(ncol(SMA_out)==3){
        SMA_out<-d
      } else {
        if(SMA_out[nrow(SMA_out),7]>d[nrow(d),7])  SMA_out<-d
      }
    }
    cat(ini,"*")
  }
  list(
    client = m,
    client_name = if (!is.null(names(X))) names(X)[m] else as.character(m),
    out = SMA_out
  )
})
stopCluster(cl)
out_matrix <- do.call(
  rbind,
  lapply(results, function(res) {
    c(
      client = res$client,
      res$out
    )
  })
)
out_table <- as.data.frame(out_matrix)
names(out_table) <- c("client","ini","K_hat","iters","iters","MSE","Acc","l","ini")
out_table$client_name <- vapply(
  results,
  function(res) res$client_name,
  character(1)
)
out_table[which.min(out_table$l),]

x <- as.matrix(out_table[, 1:9])

summary_table <- rbind(
  Mean   = colMeans(x, na.rm = TRUE),
  Max    = matrixStats::colMaxs(x, na.rm = TRUE),
  Min    = matrixStats::colMins(x, na.rm = TRUE),
  Median = matrixStats::colMedians(x, na.rm = TRUE)
)

cat("SC summary statistics (five significant digits):\n")
cat("Statistic & ", paste(colnames(summary_table), collapse = " & "), "\n", sep = "")
for (i in seq_len(nrow(summary_table))) {
  formatted_values <- vapply(seq_len(ncol(summary_table)), function(j) {
    if (colnames(summary_table)[j] %in% c("MSE", "Acc")) {
      formatC(summary_table[i, j], format = "f", digits = 5L)
    } else {
      format(summary_table[i, j], digits = 5, trim = TRUE)
    }
  }, character(1))
  cat(
    rownames(summary_table)[i], " & ",
    paste(formatted_values, collapse = " & "),
    "\n", sep = ""
  )
}

section("LATEX TABLE")

E_FICR_mse <- as.numeric(E_FICR_values["MSE"])
E_FICR_acc <- as.numeric(E_FICR_values["Accuracy"])
aa_mse <- as.numeric(aa_score["mse"])
aa_acc <- as.numeric(aa_score["acc"])

# Optimize SC MSE and Acc independently.
sc_best_mse <- min(out_table$MSE, na.rm = TRUE)
sc_best_acc <- max(out_table$Acc, na.rm = TRUE)
sc_worst_mse <- max(out_table$MSE, na.rm = TRUE)
sc_worst_acc <- min(out_table$Acc, na.rm = TRUE)
sc_median_mse <- stats::median(out_table$MSE, na.rm = TRUE)
sc_median_acc <- stats::median(out_table$Acc, na.rm = TRUE)
sc_mean_mse <- mean(out_table$MSE, na.rm = TRUE)
sc_mean_acc <- mean(out_table$Acc, na.rm = TRUE)

latex_number <- function(x) {
  if (!is.finite(x)) return(as.character(x))

  integer_digits <- if (abs(x) < 1) 1L else floor(log10(abs(x))) + 1L
  decimal_digits <- max(0L, 5L - integer_digits)
  rounded_x <- round(x, decimal_digits)

  rounded_integer_digits <- if (abs(rounded_x) < 1) {
    1L
  } else {
    floor(log10(abs(rounded_x))) + 1L
  }
  decimal_digits <- max(0L, 5L - rounded_integer_digits)

  formatC(rounded_x, format = "f", digits = decimal_digits)
}

cat("The following block can be copied directly into the LaTeX manuscript.\n\n")
cat("\\begin{table}\n")
cat("\\TABLE\n")
cat("{Estimation results obtained by different algorithms on the PISA 2022 student questionnaire data.\\label{rd\\_table2}}\n")
cat("{\n")
cat("\\begin{tabular}{lcc|lcc}\n")
cat("\\hline\\hline\n")
cat("Alg & MSE & Acc (\\%) & Alg & MSE & Acc (\\%) \\\\\n")
cat("\\hline\n")
cat(
  "E-FICR & ", latex_number(E_FICR_mse), " & ", latex_number(E_FICR_acc),
  " & SC: best & ", latex_number(sc_best_mse), " & ", latex_number(sc_best_acc),
  " \\\\\n", sep = ""
)
cat(
  "AA & ", latex_number(aa_mse), " & ", latex_number(aa_acc),
  " & SC: worst & ", latex_number(sc_worst_mse), " & ", latex_number(sc_worst_acc),
  " \\\\\n", sep = ""
)
cat(
  "SC: median & ", latex_number(sc_median_mse), " & ", latex_number(sc_median_acc),
  " & SC: mean & ", latex_number(sc_mean_mse), " & ", latex_number(sc_mean_acc),
  " \\\\\n", sep = ""
)
cat("\\hline\\hline\n")
cat("\\end{tabular}\n")
cat("}\n")
cat("{}\n")
cat("\\end{table}\n")

save(cicr, z0, beta_init, E_FICR_out, E_FICR_values, beta_hat, aa_score, results, out_table, summary_table, iters, shred, K_hat, num_init,
     file = sub("[.]txt$", ".Rdata", out_path))
cat("\nFinished:", format(Sys.time()), "\n")
cat("Results file:", normalizePath(out_path, mustWork = FALSE), "\n")

sink()
