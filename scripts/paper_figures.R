# Plot saved summaries without refitting or drawing random numbers.
paper_figures <- function(output) {
  read_result <- function(path) {
    e <- new.env(parent = emptyenv())
    load(path, envir = e)
    e
  }
  ss_files <- list.files(output, "^S1_SS_n[0-9].*[.]Rdata$", recursive = TRUE, full.names = TRUE)
  if (length(ss_files)) {
    dat <- do.call(rbind, lapply(ss_files, function(f) {
      e <- read_result(f)
      data.frame(Method = c("Oracle", "CICR", "E-FICR", "AA", "SC"),
                 SS = e$SS, n = e$n, M = e$M, MSE = e$all[,1], Acc = e$all[,2], ell = e$all[,3])
    }))
    groups <- unique(dat[c("n", "M")])
    for (i in seq_len(nrow(groups))) {
      n <- groups$n[i]
      m <- groups$M[i]
      d <- dat[dat$n == n & dat$M == m, ]
      group_dir <- file.path(output, sprintf("n%d_m%d", n, m))
      dir.create(group_dir, showWarnings = FALSE)
      save(list = "dat", envir = list2env(list(dat = d)),
           file = file.path(group_dir, "S1_SS_figure_data.Rdata"))
    }
  }
  if (length(ss_files)) {
    source(file.path(getOption("E_FICR.root", getwd()), "scripts", "render_S1SS_ijoc.R"), local = TRUE)
    render_S1SS_ijoc(dat, output)
  }
  s3_files <- list.files(output, "^S3_n[0-9].*[.]Rdata$", full.names = TRUE)
  if (length(s3_files)) {
    dat <- do.call(rbind, lapply(s3_files, function(f) {
      e <- read_result(f)
      e$K_out <- E_FICR_result_names(e$K_out)
      do.call(rbind, lapply(c("CICR", "E_FICR", "AA", "SC"), function(method) {
        k <- e$K_out[[method]]
        data.frame(M = e$M, n = e$n, Method = if (method == "E_FICR") "E-FICR" else method,
          K = factor(e$K_choose), Count = as.integer(table(factor(k, levels = e$K_choose))),
          Successful = sum(is.finite(k)), Requested = e$seed_all)
      }))
    }))
    dat$Percent <- ifelse(dat$Successful > 0, 100*dat$Count/dat$Successful, NA_real_)
    source(file.path(getOption("E_FICR.root", getwd()), "scripts", "render_figure3_ijoc.R"), local = TRUE)
    render_figure3_ijoc(dat, output)
    save(dat, file = file.path(output,"S3_figure_data.Rdata"))
  }
  invisible(NULL)
}
