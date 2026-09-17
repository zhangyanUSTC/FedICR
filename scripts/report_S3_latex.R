print_s3_latex <- function(output) {
  files <- list.files(output, "^S3_n[0-9].*[.]Rdata$", full.names = TRUE)
  if (!length(files)) stop("No S3 result files found in ", output)
  methods <- c("Oracle", "CICR", "E_FICR", "AA", "SC")
  results <- lapply(files, function(file) {
    e <- new.env(parent = emptyenv())
    load(file, envir = e)
    summary <- e$selected_summary
    required <- c("Method", "MSE2_mean", "MSE2_SD", "ARI_mean", "ARI_SD", "Successful_N")
    if (!all(required %in% names(summary))) stop("Missing selected-K summary in ", basename(file))
    summary$Method <- sub("^(FedICR|FICR|DKR|E-FICR)$", "E_FICR", summary$Method)
    if (anyDuplicated(summary$Method) || !all(methods %in% summary$Method)) {
      stop("Missing or duplicate methods in ", basename(file))
    }
    list(M = e$M, requested = e$seed_all, summary = summary[match(methods, summary$Method), ])
  })
  clients <- vapply(results, function(e) e$M, numeric(1))
  if (anyDuplicated(clients)) stop("Expected one S3 result file per M")
  results <- results[order(clients)]

  cell <- function(summary, method, metric, M) {
    row <- summary[summary$Method == method, ]
    digits <- if (metric == "MSE2") 2L else 3L
    value <- row[[paste0(metric, "_mean")]]
    if (!is.finite(value)) return("--")
    candidates <- summary[[paste0(metric, "_mean")]]
    candidates <- as.numeric(sprintf(paste0("%.", digits, "f"), candidates[is.finite(candidates)]))
    best <- if (metric == "MSE2") min(candidates) else max(candidates)
    mean_text <- sprintf(paste0("%.", digits, "f"), value)
    if (as.numeric(mean_text) == best) mean_text <- paste0("\\textbf{", mean_text, "}")
    sd_digits <- if (metric == "MSE2" && method == "E_FICR" && M %in% c(40, 50)) 3L else digits
    sd <- row[[paste0(metric, "_SD")]]
    sd_text <- if (is.finite(sd)) sprintf(paste0("%.", sd_digits, "f"), sd) else "\\mathrm{NA}"
    paste0("$", mean_text, "_{(", sd_text, ")}$")
  }
  rows <- vapply(methods, function(method) {
    cells <- unlist(lapply(seq_along(results), function(j) {
      if (method %in% c("Oracle", "CICR") && j > 1L) return(c("-", "-"))
      e <- results[[j]]
      vapply(c("MSE2", "ARI"), function(metric) cell(e$summary, method, metric, e$M), character(1))
    }), use.names = FALSE)
    label <- if (method == "E_FICR") "E-FICR" else method
    paste0(paste(c(label, cells), collapse = " & "), " \\\\")
  }, character(1), USE.NAMES = FALSE)
  headers <- vapply(seq_along(results), function(j) {
    paste0("& \\multicolumn{2}{", if (j < length(results)) "c|" else "c",
           "}{$M=", results[[j]]$M, "$}", if (j == length(results)) " \\\\" else "")
  }, character(1))
  metrics <- rep("& $\\mathrm{MSE}_2$ & ARI", length(results))
  metrics[length(metrics)] <- paste0(metrics[length(metrics)], " \\\\")
  latex <- c(
    "\\begin{table}", "\\TABLE",
    "{$\\mathrm{MSE}_2$ and ARI of different algorithms under Setting 3, using the selected $\\hat{k}$.\\label{simu2_table2}}",
    "{", "\\renewcommand{\\arraystretch}{1.2}", "\\setlength{\\tabcolsep}{6pt}",
    paste0("\\begin{tabular}{l", paste(rep("cc", length(results)), collapse = "|"), "}"),
    "\\hline\\hline", "\\multirow{2}{*}{Alg}", headers, metrics, "\\hline",
    rows, "\\hline\\hline", "\\end{tabular}", "}", "{}", "\\end{table}"
  )
  cat("\nLaTeX table: MSE2 and ARI at the selected K\n")
  cat("Entries are means with standard deviations; bold indicates the best mean at the reported precision, including ties.\n")
  cat("Unavailable means are --; unavailable standard deviations are NA.\n")
  cat("Oracle and CICR are displayed only in the first available M column.\n")
  for (e in results) {
    incomplete <- e$summary[e$summary$Successful_N < e$requested, ]
    for (i in seq_len(nrow(incomplete))) {
      method <- sub("^E_FICR$", "E-FICR", incomplete$Method[i])
      cat(sprintf("Valid repetitions: %s, M = %d: %d of %d.\n",
                  method, e$M, incomplete$Successful_N[i], e$requested))
    }
  }
  cat("\n", paste(latex, collapse = "\n"), "\n", sep = "")
  invisible(latex)
}
