parse_options <- function(args) {
  opt <- list(task = "help", seeds = "1000", workers = "10", clients = NULL,
              output = NULL, data = NULL, ss_count = "20")
  if (!length(args) || identical(args, "--help")) return(opt)
  if (length(args) %% 2L) stop("Options require a value; use --help.")
  for (i in seq.int(1L, length(args), by = 2L)) {
    key <- gsub("-", "_", sub("^--", "", args[i]), fixed = TRUE)
    if (!startsWith(args[i], "--") || !key %in% names(opt)) stop("Unknown option: ", args[i])
    opt[[key]] <- args[i + 1L]
  }
  for (key in c("seeds", "workers", "ss_count", if (!is.null(opt$clients)) "clients")) {
    value <- suppressWarnings(as.numeric(opt[[key]]))
    if (!is.finite(value) || value < 1 || value > .Machine$integer.max || value != floor(value))
      stop("--", key, " must be a positive integer.")
    opt[[key]] <- as.integer(value)
  }
  opt
}

simulation_settings <- function(task, ss_count = 20L) {
  settings <- rbind(
    data.frame(task = "S1_M", n = c(400,160,100), M = c(10,25,40), SS = 10,
               script = "S1_simu.R"),
    data.frame(task = "S1_N", n = 200, M = 20, SS = c(5,10,20),
               script = "S12_simu_ppfl.R"),
    do.call(rbind, lapply(c(10,40,25), function(m)
      data.frame(task = "S1_SS", n = 4000/m, M = m, SS = seq(110, by = 10, length.out = ss_count),
                 script = "S12_simu.R"))),
    data.frame(task = "S2", n = c(400,100,160), M = c(10,40,25), SS = 10,
               script = "S2_simu_unba.R"),
    data.frame(task = "S3", n = c(600,300,240), M = c(20,40,50), SS = 10,
               script = "S3_simu_cva.R")
  )
  if (task != "all") settings <- settings[settings$task == task, , drop = FALSE]
  settings
}

# Track workers for cleanup on errors and interrupts.
worker_registry <- new.env(parent = emptyenv())
worker_registry$clusters <- list()
start_workers <- function(n) {
  cl <- parallel::makeCluster(n)
  worker_registry$clusters[[length(worker_registry$clusters) + 1L]] <- cl
  parallel::clusterCall(cl, function(path, kind) {
    do.call(RNGkind, as.list(kind))
    suppressPackageStartupMessages(source(path, local = .GlobalEnv))
    NULL
  }, getOption("E_FICR.core"), RNGkind())
  cl
}

run_project <- function(root, args) {
  opt <- parse_options(args)
  tasks <- c("S1_M", "S1_N", "S1_SS", "S2", "S3", "all",
             "pisa-cva", "pisa-compare", "pisa-sensitivity", "pisa-all")
  if (opt$task == "help") {
    cat("Rscript run.R --task TASK [--seeds 1000] [--workers 10] [--clients M]\n",
        "Tasks: ", paste(tasks, collapse = ", "), "\n",
        "Optional: --output NEW_DIRECTORY; --data PREPARED_DATA_DIRECTORY\n",
        "S1_SS or all: --ss-count 20 (SS starts at 110, step 10; default ends at 300)\n",
        "--seeds and --clients apply only to simulations. No arguments: help only.\n", sep = "")
    return(invisible(NULL))
  }
  if (!opt$task %in% tasks) stop("Unknown task: ", opt$task)
  if (any(c("--ss-count", "--ss_count") %in% args) && !opt$task %in% c("S1_SS", "all"))
    stop("--ss-count applies only to S1_SS or all.")
  is_pisa <- startsWith(opt$task, "pisa-")
  if (is_pisa && any(c("--seeds", "--clients") %in% args))
    stop("PISA uses the fixed paper design; --seeds/--clients are simulation-only.")
  settings <- if (!is_pisa) simulation_settings(opt$task) else NULL
  if (!is_pisa && opt$ss_count != 20L) settings <- simulation_settings(opt$task, ss_count = opt$ss_count)
  if (!is.null(opt$clients)) settings <- settings[settings$M == opt$clients, , drop = FALSE]
  if (!is_pisa && !nrow(settings)) stop("No setting matches --clients.")
  output <- if (is.null(opt$output)) {
    suffix <- if (is_pisa) {
      paste0(format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), "_", Sys.getpid())
    } else {
      format(Sys.time(), "%y-%m-%d_%H-%M")
    }
    file.path(root, "output", paste0(opt$task, "_", suffix))
  } else opt$output
  if (dir.exists(output) || file.exists(output)) stop("Output already exists; choose a new --output directory.")
  dir.create(output, recursive = TRUE)
  output <- normalizePath(output)
  data <- if (is.null(opt$data)) file.path(root, "data") else normalizePath(opt$data)
  old_options <- options(E_FICR.root = root, E_FICR.output = output,
                        E_FICR.workers = opt$workers, E_FICR.data = data,
                        E_FICR.core = file.path(root, if (is_pisa) "pisa/pisa_federated.R" else "simulation/S0_basic.R"))
  sink_level <- sink.number()
  oldwd <- getwd()
  on.exit({
    for (cl in worker_registry$clusters) try(parallel::stopCluster(cl), silent = TRUE)
    worker_registry$clusters <- list()
    if (requireNamespace("foreach", quietly = TRUE)) foreach::registerDoSEQ()
    while (sink.number() > sink_level) sink()
    options(old_options)
    setwd(oldwd)
  }, add = TRUE)
  suppressPackageStartupMessages(source(getOption("E_FICR.core"), local = .GlobalEnv))
  info <- list(options = opt, settings = settings, started = Sys.time(),
               session = sessionInfo(), RNGkind = RNGkind(),
               code_md5 = tools::md5sum(list.files(root, "[.]R$", recursive = TRUE, full.names = TRUE)))
  if (is_pisa) info$data_md5 <- tools::md5sum(file.path(data,
    c("pisa_federated_data_pv1.rds", "pisa_federated_data_pv.rds")))
  saveRDS(info, file.path(output, "run_info.rds"))
  if (is_pisa) {
    scripts <- c("pisa-cva" = "pisa_cva.R", "pisa-compare" = "pisa_all_alg.R",
                 "pisa-sensitivity" = "pisa_pv_sensitivity.R")
    if (opt$task != "pisa-all") scripts <- scripts[opt$task]
    for (s in scripts) source(file.path(root, "pisa", s), local = .GlobalEnv)
  } else {
    for (j in seq_len(nrow(settings))) {
      cfg <- as.list(settings[j, ])
      task_rows <- which(settings$task == cfg$task)
      if (j == min(task_rows)) {
        task_output <- if (opt$task == "all") file.path(output, cfg$task) else output
        dir.create(task_output, showWarnings = FALSE)
        options(E_FICR.output = task_output)
        sink(file.path(task_output, "output.txt"), split = TRUE)
      }
      cfg <- c(cfg, list(K = 4L, p = 4L, num_init = 3L, iters = 30L,
                         seed_all = opt$seeds, K_choose = 2:6))
      cfg$file_path <- task_output
      if (cfg$task == "S1_SS") {
        cfg$file_path <- file.path(task_output, sprintf("n%d_m%d", cfg$n, cfg$M))
        dir.create(cfg$file_path, showWarnings = FALSE)
      }
      cfg$beta0 <- cfg$SS * diag(cfg$p)
      cat(sprintf("\nRunning %s: n=%d, M=%d, K=%d, SS=%d, iters=%d, seeds=%d, workers=%d\n",
                  cfg$task, cfg$n, cfg$M, cfg$K, cfg$SS, cfg$iters, cfg$seed_all, min(cfg$seed_all, opt$workers)))
      list2env(cfg, envir = .GlobalEnv)
      source(file.path(root, "simulation", cfg$script), local = .GlobalEnv)
      cat(format(Sys.time(), "Finished: %Y-%m-%d %H:%M:%S\n"))
      s1_rows <- which(settings$task == "S1_M")
      if (cfg$task == "S1_M" && length(s1_rows) == 3L &&
          setequal(settings$M[s1_rows], c(10, 25, 40)) && j == max(s1_rows)) {
        source(file.path(root, "scripts", "render_S1M_ijoc.R"), local = .GlobalEnv)
        render_S1M_ijoc(task_output)
      }
      if (j == max(task_rows)) {
        source(file.path(root, "scripts", "paper_figures.R"), local = .GlobalEnv)
        paper_figures(task_output)
        cat("\nTask results directory:", task_output, "\n")
        if (cfg$task == "S3") {
          source(file.path(root, "scripts", "report_S3_latex.R"), local = .GlobalEnv)
          print_s3_latex(task_output)
        }
        sink()
      }
    }
  }
  cat("\nResults directory:", output, "\n")
  invisible(output)
}
