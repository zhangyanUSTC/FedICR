source("scripts/control.R")
options(device = function(...) grDevices::pdf(file = NULL))
find_hook <- function(x) {
  if (missing(x)) return(NULL)
  if (!is.call(x)) return(NULL)
  if (identical(x[[1L]], as.name("if")) &&
      grepl("setequal(settings$M[s1_rows]", paste(deparse(x[[2L]]), collapse = " "), fixed = TRUE)) return(x)
  for (part in as.list(x)[-1L]) {
    found <- find_hook(part)
    if (!is.null(found)) return(found)
  }
  NULL
}
hook <- find_hook(body(run_project))
stopifnot(!is.null(hook))
for (task in c("S1_M", "all")) {
  for (clients in list(NULL, 10, 25, 40, c(10, 25))) {
    settings <- simulation_settings(task)
    if (!is.null(clients)) settings <- settings[settings$M %in% clients, ]
    s1_rows <- which(settings$task == "S1_M")
    triggered <- integer()
    for (j in seq_len(nrow(settings))) {
      cfg <- as.list(settings[j, ])
      if (eval(hook[[2L]])) triggered <- c(triggered, j)
    }
    expected <- if (is.null(clients)) max(s1_rows) else integer()
    stopifnot(identical(triggered, expected))
  }
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args)) {
  suppressPackageStartupMessages(source("simulation/S0_basic.R"))
  source("scripts/render_S1M_ijoc.R")
  out <- tempfile("s1_m_figure_")
  dir.create(out)
  files <- list.files(args[1L], "^S1_M_RES_.*[.]Rdata$", full.names = TRUE)
  stopifnot(length(files) == 3L, all(file.copy(files, out)))
  set.seed(42)
  before <- .Random.seed
  paths <- render_S1M_ijoc(out)
  stopifnot(identical(names(paths), "png"),
            identical(unname(basename(paths)), "S1M_combined.png"),
            length(list.files(out, "[.](pdf|rds)$")) == 0L)
  stopifnot(all(file.exists(paths)), identical(before, .Random.seed))
  image_hash <- tools::md5sum(paths)
  for (file in list.files(out, "^S1_M_RES_.*[.]Rdata$", full.names = TRUE)) {
    saved <- new.env(parent = emptyenv())
    load(file, envir = saved)
    saved$res_list <- E_FICR_result_names(saved$res_list)
    save(list = ls(saved), file = file, envir = saved)
  }
  updated_paths <- render_S1M_ijoc(out)
  stopifnot(identical(image_hash, tools::md5sum(updated_paths)),
            identical(before, .Random.seed))
}
cat("PASS: S1_M figure gating and optional old/new result-key rendering parity.\n")
