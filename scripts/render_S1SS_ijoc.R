render_S1SS_ijoc <- function(dat, output, root = getOption("E_FICR.root", getwd())) {
  methods <- c("Oracle", "CICR", "E-FICR", "AA", "SC")
  read_ss_results <- function(d) {
    d <- d[order(d$SS, match(d$Method, methods)), ]
    stopifnot(!anyDuplicated(d[c("SS", "Method")]),
              all(table(d$SS) == 5L), all(d$Method %in% methods))
    data.frame(MSE = d$MSE, Acc = d$Acc, ell = log(d$ell), SS = d$SS)
  }
  candidates <- unique(c(Sys.getenv("E_FICR_PYTHON"), Sys.which("python3"),
    file.path(root, ".venv", "bin", "python"),
    path.expand("~/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3")))
  candidates <- candidates[nzchar(candidates) & file.exists(candidates)]
  usable <- vapply(candidates, function(p) {
    identical(suppressWarnings(system2(p, c("-c", shQuote("from PIL import Image")),
                                      stdout = FALSE, stderr = FALSE)), 0L)
  }, logical(1))
  if (!any(usable)) stop("Formula overlay requires Python with Pillow. Install Pillow and set E_FICR_PYTHON to that Python executable.")
  python <- candidates[which(usable)[1L]]
  overlay_formula <- function(path, layout) {
    asset <- "ss_formula_log_mean.png"
    args <- c(file.path(root, "scripts", "overlay_ss_formula.py"),
              path, file.path(root, "assets", asset), layout)
    height_fraction <- if (layout == "single") 0.17 else 0.17 * 2.500875 / 4.4
    args <- c(args, "--height-fraction", as.character(height_fraction))
    status <- system2(python, shQuote(args))
    if (status != 0L) stop("Formula overlay failed: ", path)
    message("Created: ", path)
  }
plot_mse2_ijoc <- function(iter_my, iter_non, iter_sin, iter_ave,
                           oracle, name, iters, ss_values) {
  plot_data <- data.frame(
    SS = rep(ss_values, 5L),
    value = c(oracle, iter_my, iter_non, iter_sin, iter_ave),
    algorithm = factor(
      rep(c("Oracle", "E-FICR", "CICR", "SC", "AA"), each = iters),
      levels = c("Oracle", "CICR", "E-FICR", "AA", "SC")
    )
  )

  ggplot(plot_data, aes(
    x = SS, y = value, colour = algorithm,
    shape = algorithm, linetype = algorithm
  )) +
    geom_line(linewidth = 0.55) +
    geom_point(size = 1.3, stroke = 0.32) +
    labs(x = "SS", y = name, colour = NULL, shape = NULL, linetype = NULL) +
    scale_colour_manual(values = c(
      "Oracle" = "#D89000", "CICR" = "#4477AA",
      "E-FICR" = "#EE6677", "AA" = "#7B6BAF", "SC" = "#66A99A"
    )) +
    scale_shape_manual(values = c(
      "Oracle" = 13, "CICR" = 17, "E-FICR" = 15,
      "AA" = 16, "SC" = 18
    )) +
    scale_linetype_manual(values = c(
      "Oracle" = "dashed", "CICR" = "solid", "E-FICR" = "dashed",
      "AA" = "dotdash", "SC" = "longdash"
    )) +
    scale_x_continuous(
      breaks = seq(120, 300, by = 40),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    theme_bw(base_size = 14.1, base_family = "sans") +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "#E6E6E6", linewidth = 0.28),
      panel.border = element_rect(colour = "#555555", linewidth = 0.35),
      axis.text = element_text(colour = "#222222", size = 9.84),
      axis.text.x = element_text(margin = margin(t = 3)),
      axis.title = element_text(colour = "#111111", size = 12),
      axis.title.x = element_text(margin = margin(t = 4)),
      axis.title.y = element_text(margin = margin(r = 3)),
      axis.ticks = element_line(linewidth = 0.3),
      axis.ticks.length = unit(1.5, "mm"),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.key.width = unit(16, "mm"),
      legend.key.height = unit(3.2, "mm"),
      legend.text = element_text(size = 13.8),
      plot.margin = margin(1.5, 2, 3, 2)
    ) +
    guides(
      colour = guide_legend(nrow = 1, override.aes = list(linewidth = 0.8)),
      shape = guide_legend(nrow = 1),
      linetype = guide_legend(nrow = 1)
    )
}

render_one_m <- function(results, M_value, folder_path) {
  part_data <- read_ss_results(results)
  iters <- nrow(part_data) / 5L

  make_metric_plot <- function(i, label) {
    plot_mse2_ijoc(
      iter_my = part_data[(seq_len(iters) * 5) - 2, i],
      iter_non = part_data[(seq_len(iters) * 5) - 3, i],
      iter_sin = part_data[seq_len(iters) * 5, i],
      iter_ave = part_data[(seq_len(iters) * 5) - 1, i],
      oracle = part_data[(seq_len(iters) * 5) - 4, i],
      name = label,
      iters = iters,
      ss_values = unique(part_data$SS)
    )
  }

  p1 <- make_metric_plot(1, "MSE")
  p2 <- make_metric_plot(2, "Acc")
  p3 <- make_metric_plot(3, " ")

  legend <- simulation_legend(
    p3 + theme(legend.position = "bottom"),
    position = "bottom"
  )
  panel <- cowplot::plot_grid(
    p1 + theme(legend.position = "none"),
    p2 + theme(legend.position = "none"),
    p3 + theme(legend.position = "none"),
    ncol = 3,
    align = "hv",
    axis = "tblr"
  )
  combined_plot <- cowplot::plot_grid(
    panel, legend, ncol = 1, rel_heights = c(1, 0.10)
  )

  png_path <- file.path(folder_path, sprintf("S1_SS_M%d.png", M_value))
  ggsave(filename = png_path, plot = combined_plot,
         width = 8.24, height = 2.500875, units = "in", dpi = 300, bg = "white")
  overlay_formula(png_path, "single")
  invisible(png_path)
}
plot_snr_ijoc <- function(iter_my, iter_non, iter_sin, iter_ave,
                          oracle, y_label, iters, x_label, panel_title, ss_values) {
  plot_data <- data.frame(
    SS = rep(ss_values, 5L),
    value = c(oracle, iter_my, iter_non, iter_sin, iter_ave),
    algorithm = factor(
      rep(c("Oracle", "E-FICR", "CICR", "SC", "AA"), each = iters),
      levels = c("Oracle", "CICR", "E-FICR", "AA", "SC")
    )
  )

  ggplot(plot_data, aes(
    x = SS, y = value, colour = algorithm,
    shape = algorithm, linetype = algorithm
  )) +
    geom_line(linewidth = 0.55) +
    geom_point(size = 1.3, stroke = 0.32) +
    labs(
      x = x_label, y = y_label, title = panel_title,
      colour = NULL, shape = NULL, linetype = NULL
    ) +
    scale_colour_manual(values = c(
      "Oracle" = "#D89000", "CICR" = "#4477AA",
      "E-FICR" = "#EE6677", "AA" = "#7B6BAF", "SC" = "#66A99A"
    )) +
    scale_shape_manual(values = c(
      "Oracle" = 13, "CICR" = 17, "E-FICR" = 15,
      "AA" = 16, "SC" = 18
    )) +
    scale_linetype_manual(values = c(
      "Oracle" = "dashed", "CICR" = "solid", "E-FICR" = "dashed",
      "AA" = "dotdash", "SC" = "longdash"
    )) +
    scale_x_continuous(
      breaks = seq(120, 300, by = 40),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    theme_bw(base_size = 14.1, base_family = "sans") +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "#E6E6E6", linewidth = 0.28),
      panel.border = element_rect(colour = "#555555", linewidth = 0.35),
      axis.text = element_text(colour = "#222222", size = 9.84),
      axis.text.x = element_text(margin = margin(t = 2)),
      axis.title = element_text(colour = "#111111", size = 12),
      axis.title.x = element_text(margin = margin(t = 2)),
      axis.title.y = element_text(margin = margin(r = 3)),
      axis.ticks = element_line(linewidth = 0.3),
      axis.ticks.length = unit(1.4, "mm"),
      plot.title = element_text(
        colour = "#111111", size = 13.5, hjust = 0,
        margin = margin(b = 2)
      ),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.key.width = unit(16, "mm"),
      legend.key.height = unit(3.2, "mm"),
      legend.text = element_text(size = 13.8),
      plot.margin = margin(1, 2, 1, 2)
    ) +
    guides(
      colour = guide_legend(nrow = 1, override.aes = list(linewidth = 0.8)),
      shape = guide_legend(nrow = 1),
      linetype = guide_legend(nrow = 1)
    )
}

make_row <- function(folder_path, M_value, top_row = FALSE) {
  part_data <- read_ss_results(folder_path)
  iters <- nrow(part_data) / 5L
  seq_i <- seq_len(iters)
  x_label <- if (top_row) NULL else "SS"
  row_title <- sprintf("(%s)  M = %d", if (M_value == 25) "a" else "b", M_value)

  make_metric <- function(i, y_label, title = " ") {
    plot_snr_ijoc(
      iter_my = part_data[(seq_i * 5) - 2, i],
      iter_non = part_data[(seq_i * 5) - 3, i],
      iter_sin = part_data[seq_i * 5, i],
      iter_ave = part_data[(seq_i * 5) - 1, i],
      oracle = part_data[(seq_i * 5) - 4, i],
      y_label = y_label,
      iters = iters,
      x_label = x_label,
      panel_title = title,
      ss_values = unique(part_data$SS)
    )
  }

  list(
    make_metric(1, "MSE", row_title),
    make_metric(2, "Acc"),
    make_metric(3, " ")
  )
}

render_pair <- function(data25, data40, output_folder) {
plots_25 <- make_row(data25, 25, top_row = TRUE)
plots_40 <- make_row(data40, 40, top_row = FALSE)
all_plots <- c(plots_25, plots_40)

legend <- simulation_legend(
  plots_40[[3]] + theme(legend.position = "bottom"),
  position = "bottom"
)
panel <- cowplot::plot_grid(
  plotlist = lapply(all_plots, function(p) p + theme(legend.position = "none")),
  ncol = 3,
  nrow = 2,
  align = "hv",
  axis = "tblr"
)
combined_plot <- cowplot::plot_grid(
  panel,
  legend,
  ncol = 1,
  rel_heights = c(1, 0.10)
)

png_path <- file.path(output_folder, "S1_SS_M25_M40_combined.png")
png(filename = png_path, width = 8.24, height = 4.4, units = "in",
    res = 300, type = "cairo", bg = "white")
tryCatch(print(combined_plot), finally = dev.off())
overlay_formula(png_path, "combined")
invisible(png_path)
}

  groups <- unique(dat[c("n", "M")])
  for (j in seq_len(nrow(groups))) {
    d <- dat[dat$n == groups$n[j] & dat$M == groups$M[j], ]
    group_dir <- file.path(output, sprintf("n%d_m%d", groups$n[j], groups$M[j]))
    dir.create(group_dir, showWarnings = FALSE)
    render_one_m(d, groups$M[j], group_dir)
  }
  if (all(c(25, 40) %in% groups$M)) {
    stopifnot(sum(groups$M == 25) == 1L, sum(groups$M == 40) == 1L)
    data25 <- dat[dat$M == 25, ]
    data40 <- dat[dat$M == 40, ]
    if (setequal(data25$SS, data40$SS)) {
      render_pair(data25, data40, output)
    } else {
      message("Skipping M25/M40 combined figure: SS settings differ.")
    }
  }
  invisible(NULL)
}
