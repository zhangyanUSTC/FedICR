source("scripts/report_S3_latex.R")
make_results <- function(output, clients = c(20L, 40L, 50L), legacy = FALSE) {
  dir.create(output)
  for (M in clients) {
    n <- 12000 / M
    seed_all <- 1000L
    selected_summary <- data.frame(
      Method = c("Oracle", "CICR", if (legacy) "FICR" else "E_FICR", "AA", "SC"),
      MSE2_mean = c(27.297046, 27.375893, 27.377302, 29.576860, 39.247136),
      MSE2_SD = c(.589179, .595323, .594794, 9.592724, 15.039363),
      ARI_mean = c(.669227, .668855, .668842, .643300, .529890),
      ARI_SD = c(.006429, .006456, .006447, .096192, .127613),
      Successful_N = c(1000L, 1000L, 1000L, 1000L, if (M == 50L) 996L else 1000L)
    )
    save(M, n, seed_all, selected_summary,
         file = file.path(output, sprintf("S3_n%d_m%d_fixture.Rdata", n, M)))
  }
}
out <- tempfile("s3_latex_")
make_results(out)
set.seed(42)
before <- .Random.seed
hashes <- tools::md5sum(list.files(out, full.names = TRUE))
log <- capture.output(latex <- print_s3_latex(out))
stopifnot(identical(before, .Random.seed), identical(hashes, tools::md5sum(names(hashes))),
          is.null(names(latex)), latex[1] == "\\begin{table}", tail(latex, 1) == "\\end{table}",
          "\\begin{tabular}{lcc|cc|cc}" %in% latex,
          any(grepl("SC, M = 50: 996 of 1000", log, fixed = TRUE)))
rows <- latex[grepl("^(Oracle|CICR|E-FICR|AA|SC) &", latex)]
stopifnot(length(rows) == 5L, all(lengths(strsplit(rows, "&", fixed = TRUE)) == 7L),
          all(endsWith(rows, "\\\\")),
          any(grepl("$27.38_{(0.595)}$", rows, fixed = TRUE)),
          any(grepl("$\\textbf{0.669}_{(0.006)}$", rows, fixed = TRUE)))
legacy <- tempfile("s3_latex_legacy_")
make_results(legacy, legacy = TRUE)
invisible(capture.output(legacy_latex <- print_s3_latex(legacy)))
stopifnot(identical(latex, legacy_latex), identical(before, .Random.seed))

partial <- tempfile("s3_latex_partial_")
make_results(partial, clients = 40L)
invisible(capture.output(partial_latex <- print_s3_latex(partial)))
stopifnot("\\begin{tabular}{lcc}" %in% partial_latex,
          any(grepl("\\multicolumn{2}{c}{$M=40$}", partial_latex, fixed = TRUE)),
          !any(grepl("$M=20$", partial_latex, fixed = TRUE)))
partial_rows <- partial_latex[grepl("^(Oracle|CICR|E-FICR|AA|SC) &", partial_latex)]
stopifnot(all(lengths(strsplit(partial_rows, "&", fixed = TRUE)) == 3L),
          !any(grepl(" & - & -", partial_rows, fixed = TRUE)))
file <- list.files(partial, full.names = TRUE)
e <- new.env()
load(file, envir = e)
e$selected_summary[5L, c("MSE2_mean", "MSE2_SD", "ARI_mean", "ARI_SD")] <- NA_real_
e$selected_summary$Successful_N[5L] <- 0L
e$selected_summary[4L, c("MSE2_SD", "ARI_SD")] <- NA_real_
e$selected_summary$Successful_N[4L] <- 1L
save(list = ls(e), file = file, envir = e)
missing_log <- capture.output(missing_latex <- print_s3_latex(partial))
stopifnot("SC & -- & -- \\\\" %in% missing_latex,
          any(grepl("\\mathrm{NA}", missing_latex, fixed = TRUE)),
          any(grepl("SC, M = 40: 0 of 1000", missing_log, fixed = TRUE)),
          identical(before, .Random.seed))
cat("PASS: full/partial LaTeX tables, precision, ties, legacy names, missing results and unchanged files/RNG.\n")

args <- commandArgs(trailingOnly = TRUE)
if (length(args)) {
  files <- list.files(args[1L], full.names = TRUE)
  hashes <- tools::md5sum(files)
  saved <- readLines(file.path(args[1L], "output.txt"))
  start <- which(saved == "\\begin{table}")
  invisible(capture.output(generated <- print_s3_latex(args[1L])))
  stopifnot(length(start) == 1L,
            identical(generated, saved[seq.int(start, length.out = length(generated))]),
            identical(hashes, tools::md5sum(files)), identical(before, .Random.seed))
  cat("PASS: generated table exactly matches the approved table in the existing S3 report.\n")
}
