suppressPackageStartupMessages(source("simulation/S0_basic.R"))
source("scripts/render_S1SS_ijoc.R")
methods <- c("Oracle", "CICR", "E-FICR", "AA", "SC")
dat <- expand.grid(Method = methods, SS = seq(110, 300, 10), M = c(10, 25, 40),
                   stringsAsFactors = FALSE)
dat$n <- 4000 / dat$M
dat$MSE <- match(dat$Method, methods) * 100 / dat$SS
dat$Acc <- 100 - dat$MSE
dat$ell <- dat$MSE * 4000
out <- tempfile("ss_figures_")
dir.create(out)
set.seed(42)
before <- .Random.seed
render_S1SS_ijoc(dat, out)
stopifnot(identical(before, .Random.seed),
          file.exists(file.path(out, "n400_m10", "S1_SS_M10.png")),
          file.exists(file.path(out, "S1_SS_M25_M40_combined.png")),
          length(list.files(out, "[.]png$", recursive = TRUE)) == 4L,
          length(list.files(out, "[.](pdf|rds)$", recursive = TRUE)) == 0L)
for (m in c(10, 25)) {
  partial <- tempfile("ss_partial_")
  dir.create(partial)
  render_S1SS_ijoc(dat[dat$M == m & dat$SS <= 130, ], partial)
  stopifnot(length(list.files(partial, "[.]png$", recursive = TRUE)) == 1L,
            !file.exists(file.path(partial, "S1_SS_M25_M40_combined.png")),
            identical(before, .Random.seed))
}
cat("PASS: full/partial SS figures, combined-figure gating, PNG-only output, unchanged RNG.\n")
