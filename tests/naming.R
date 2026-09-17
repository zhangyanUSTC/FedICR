root <- normalizePath(".")
options(device = function(...) grDevices::pdf(file = NULL))
files <- list.files(c("simulation", "pisa", "scripts", "tests"),
                    "[.]R$", recursive = TRUE, full.names = TRUE)
legacy <- "^(FedICR|FICR|DKR)(_|$)|^(fed|eficr)_|^fedicr[._]|^FEDICR_"
for (file in files) {
  tokens <- getParseData(parse(file, keep.source = TRUE))
  symbols <- tokens$text[grepl("SYMBOL", tokens$token)]
  stopifnot(!any(grepl(legacy, symbols)))
}
suppressPackageStartupMessages(source("simulation/S0_basic.R"))
stopifnot(is.function(E_FICR_alg))
set.seed(42)
before <- .Random.seed
for (key in c("FedICR", "FICR", "DKR", "E-FICR", "E_FICR")) {
  old <- setNames(list(matrix(1:4, 2)), key)
  updated <- E_FICR_result_names(old)
  stopifnot(identical(names(updated), "E_FICR"),
            identical(updated[[1]], old[[1]]), identical(names(old), key))
}
stopifnot(identical(before, .Random.seed))
cat("PASS: current identifiers and legacy result-key compatibility.\n")

source("scripts/paper_figures.R")
figure_hashes <- character()
for (key in c("FICR", "E_FICR")) {
  out <- tempfile("s3_naming_")
  dir.create(out)
  K_choose <- 2:6
  M <- 10L
  n <- 120L
  seed_all <- 3L
  K_out <- data.frame(seed = 1:3, CICR = c(4L, 4L, 3L),
                      E_FICR = c(4L, 4L, 4L), AA = c(4L, 3L, 5L), SC = c(2L, 4L, 3L))
  names(K_out)[3L] <- key
  file <- file.path(out, "S3_n120_m10_fixture.Rdata")
  save(K_out, K_choose, M, n, seed_all, file = file)
  file_hash <- tools::md5sum(file)
  paper_figures(out)
  stopifnot(identical(file_hash, tools::md5sum(file)))
  figure_hashes <- c(figure_hashes, unname(tools::md5sum(file.path(out, "S3_selected_K.png"))))
}
stopifnot(length(unique(figure_hashes)) == 1L, identical(before, .Random.seed))
cat("PASS: identical S3 figures from old/new keys; saved inputs and RNG unchanged.\n")

# Optional pre-rename source snapshot for exact regression checks.
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) quit(status = 0L)
baseline <- normalizePath(args[1L])
rename_text <- function(x) {
  x <- gsub("(?<![A-Za-z0-9_-])(FedICR|FICR|fedicr|FEDICR|eficr)(?=[_.]|\\b)",
            "E_FICR", x, perl = TRUE)
  x <- gsub("DKR", "E_FICR", x, fixed = TRUE)
  gsub("\\bfed_", "E_FICR_", x, perl = TRUE)
}
normalize_result <- function(x) {
  attrs <- attributes(x)
  if (is.list(x)) x <- lapply(x, normalize_result)
  if (is.character(x)) x <- rename_text(x)
  if (!is.null(attrs)) attributes(x) <- lapply(attrs, normalize_result)
  x
}
read_core <- function(path) {
  e <- new.env(parent = globalenv())
  suppressPackageStartupMessages(sys.source(path, envir = e))
  e
}
check_bodies <- function(old, current) {
  for (name in ls(old)) {
    if (!is.function(old[[name]]) || name == "print_simulation_table") next
    new_name <- rename_text(name)
    stopifnot(is.function(current[[new_name]]))
    for (part in c("formals", "body")) {
      extract <- get(part, envir = baseenv())
      a <- rename_text(paste(deparse(extract(old[[name]])), collapse = "\n"))
      b <- paste(deparse(extract(current[[new_name]])), collapse = "\n")
      if (!identical(parse(text = a, keep.source = FALSE),
                     parse(text = b, keep.source = FALSE))) {
        stop("Function differs beyond renaming: ", name, " / ", part)
      }
    }
  }
}
check_fit <- function(old, current, old_name, arguments) {
  set.seed(731)
  a <- do.call(old[[old_name]], arguments)
  rng <- .Random.seed
  set.seed(731)
  b <- do.call(current[[rename_text(old_name)]], arguments)
  stopifnot(identical(rng, .Random.seed), identical(normalize_result(a), b))
}
old <- read_core(file.path(baseline, "simulation", "S0_basic.R"))
current <- read_core(file.path(root, "simulation", "S0_basic.R"))
check_bodies(old, current)
cat("PASS: all simulation function bodies/formals differ only by names (excluding reporting).\n")
for (seed in 1:2) {
  for (task in c("S1_M", "S1_N", "S1_SS", "S2", "S3")) {
    ss <- switch(task, S1_N = 5, S1_SS = 110, 10)
    fun <- switch(task, S1_M = "S1_main", S1_N = "S12_main", S1_SS = "S12_main",
                  S2 = "S2_main", S3 = "S3_main")
    arguments <- list(beta0 = ss * diag(4), num_init = 3, seed = seed,
                      iters = 30, n = 120, M = 10, shred = .0001,
                      balance = task != "S2")
    if (task == "S3") {
      arguments <- c(arguments, list(K_choose = 2:6, evaluate_selected = TRUE))
    } else arguments$K_hat <- 4
    if (task %in% c("S1_N", "S1_SS")) arguments$run_sc2 <- FALSE
    check_fit(old, current, fun, arguments)
    cat("PASS: pre/post numerical results and RNG state:", task, "seed", seed, "\n")
  }
}
old <- read_core(file.path(baseline, "pisa", "pisa_federated.R"))
current <- read_core(file.path(root, "pisa", "pisa_federated.R"))
check_bodies(old, current)
set.seed(81)
beta <- rbind(c(1, 2, 3), c(6, -2, -1))
X <- lapply(1:6, function(i) cbind(1, matrix(rnorm(160), 80, 2)))
z <- lapply(X, function(x) diag(2)[rep(1:2, each = 40), ])
Y <- Map(function(x, w) rowSums(x * (w %*% beta)) + rnorm(80), X, z)
check_fit(old, current, "FedICR_alg",
          list(X = X, Y = Y, beta_init = beta + .2, z0 = z,
               beta0 = beta, iters = 40, shred = 1e-4))
cat("PASS: all PISA function bodies/formals, E_FICR estimates and RNG state.\n")
