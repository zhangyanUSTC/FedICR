source("scripts/render_figure3_ijoc.R")
methods <- c("CICR", "E-FICR", "AA", "SC")
dat <- expand.grid(Method = methods, K = 2:6, M = c(20L, 40L, 50L),
                   stringsAsFactors = FALSE)
dat$n <- 12000 / dat$M
dat$Percent <- 0
dat$Percent[dat$K == 4] <- 100
idx <- dat$Method == "SC" & dat$M == 50L
dat$Percent[idx] <- 100 * c(63, 41, 100, 247, 545) / 996
set.seed(42)
before <- .Random.seed
out <- tempfile("s3_figure_")
dir.create(out)
rendered <- render_figure3_ijoc(dat, out)
stopifnot(identical(before, .Random.seed), file.exists(rendered$png),
          identical(list.files(out), "S3_selected_K.png"))
plot <- ggplot2::ggplot_build(rendered$plot)
bars <- plot$data[[1L]]
labels <- plot$data[[2L]]
stopifnot(all(is.finite(bars$ymin)), all(is.finite(bars$ymax)),
          all(abs(aggregate(ymax ~ PANEL + x, bars, max)$ymax - 100) < 1e-10),
          all(labels$colour == "#111111"), nrow(labels) == 12L,
          identical(as.character(rendered$data$Scenario[1]), "M == 20 * ',' ~~ n[m] == 600"),
          identical(unname(plot$plot$scales$get_scales("fill")$map(as.character(2:6))),
                    c("#E8EEF3", "#B9CCDC", "#7898B4", "#416C91", "#183B59")))
segments <- bars[bars$fill == "#7898B4", c("PANEL", "x", "ymin", "ymax")]
centres <- merge(labels[c("PANEL", "x", "y")], segments, by = c("PANEL", "x"))
stopifnot(all(abs(centres$y - (centres$ymin + centres$ymax) / 2) < 1e-10))
size <- readBin(rendered$png, "raw", n = 24L)
png_number <- function(bytes) sum(as.integer(bytes) * 256^(3:0))
stopifnot(png_number(size[17:20]) == 2472, png_number(size[21:24]) == 786)

partial <- dat[dat$M == 50L, ]
partial$Percent[partial$Method == "SC"] <- NA_real_
rendered_partial <- render_figure3_ijoc(partial, out)
stopifnot(length(unique(rendered_partial$data$Scenario)) == 1L,
          identical(before, .Random.seed))
cat("PASS: IJOC palette, centred black labels, complete 100% stacks, PNG dimensions, partial/missing data and unchanged RNG.\n")
