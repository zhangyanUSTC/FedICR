render_figure3_ijoc <- function(dat, target_folder, label_colour = "#111111") {
  algorithm_names <- c("CICR", "E-FICR", "AA", "SC")
  K_choose <- 2:6
  stopifnot(all(c("Method", "K", "Percent", "M", "n") %in% names(dat)),
            all(dat$Method %in% algorithm_names),
            all(as.character(dat$K) %in% as.character(K_choose)))
  selection_data <- data.frame(
    Algorithm = factor(dat$Method, levels = algorithm_names),
    K = factor(dat$K, levels = K_choose),
    Percentage = dat$Percent,
    M = as.integer(dat$M),
    n_m = as.integer(dat$n)
  )

  # Match ggplot's reversed stacking order when centring the K = 4 labels.
  selection_data$label_y <- NA_real_
  for (panel_M in unique(selection_data$M)) {
    for (algorithm in algorithm_names) {
      idx <- which(selection_data$M == panel_M & selection_data$Algorithm == algorithm)
      idx <- idx[order(as.numeric(as.character(selection_data$K[idx])), decreasing = TRUE)]
      selection_data$label_y[idx] <- cumsum(selection_data$Percentage[idx]) -
        selection_data$Percentage[idx] / 2
    }
  }

  k_colours <- c(
    `2` = "#E8EEF3",
    `3` = "#B9CCDC",
    `4` = "#7898B4",
    `5` = "#416C91",
    `6` = "#183B59"
  )
  panels <- unique(selection_data[c("M", "n_m")])
  panels <- panels[order(panels$M), , drop = FALSE]
  if (anyDuplicated(panels$M)) stop("Expected one sample size per M in Figure 3")
  scenario_labels <- sprintf("M == %d * ',' ~~ n[m] == %d", panels$M, panels$n_m)
  selection_data$Scenario <- factor(
    selection_data$M, levels = panels$M, labels = scenario_labels
  )

  combined_plot <- ggplot2::ggplot(
    selection_data,
    ggplot2::aes(x = Algorithm, y = Percentage, fill = K)
  ) +
    ggplot2::geom_col(width = 0.68, colour = "white", linewidth = 0.22, na.rm = TRUE) +
    ggplot2::geom_text(
      data = subset(selection_data, as.character(K) == "4"),
      ggplot2::aes(y = label_y, label = sprintf("%.1f%%", Percentage)),
      size = 3.15,
      family = "sans",
      fontface = "plain",
      colour = label_colour,
      show.legend = FALSE,
      na.rm = TRUE
    ) +
    ggplot2::facet_wrap(
      ggplot2::vars(Scenario), nrow = 1,
      labeller = ggplot2::label_parsed
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Selection frequency (%)",
      fill = expression(k * ":")
    ) +
    ggplot2::scale_fill_manual(
      values = k_colours,
      breaks = as.character(K_choose),
      drop = FALSE
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 100),
      breaks = seq(0, 100, by = 25),
      oob = scales::squish,
      expand = ggplot2::expansion(mult = c(0, 0.005))
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_bw(base_size = 14.1, base_family = "sans") +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(colour = "#E3E3E3", linewidth = 0.28),
      panel.border = ggplot2::element_rect(colour = "#555555", linewidth = 0.38),
      panel.spacing.x = grid::unit(3.5, "mm"),
      axis.text = ggplot2::element_text(colour = "#222222", size = 9.8),
      axis.text.x = ggplot2::element_text(size = 10.2, margin = ggplot2::margin(t = 3)),
      axis.title.y = ggplot2::element_text(
        colour = "#111111", size = 11.5, margin = ggplot2::margin(r = 4)
      ),
      axis.ticks = ggplot2::element_line(colour = "#555555", linewidth = 0.3),
      axis.ticks.length = grid::unit(1.4, "mm"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(
        colour = "#111111", size = 12.2, face = "plain",
        margin = ggplot2::margin(b = 3)
      ),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.key.width = grid::unit(12, "mm"),
      legend.key.height = grid::unit(3.8, "mm"),
      legend.title = ggplot2::element_text(size = 11.5),
      legend.text = ggplot2::element_text(size = 11.5),
      legend.spacing.x = grid::unit(1.0, "mm"),
      legend.margin = ggplot2::margin(t = -1, b = -1),
      plot.margin = ggplot2::margin(2.5, 3, 1.5, 2.5)
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(
      nrow = 1,
      title.position = "left",
      byrow = TRUE
    ))

  png_path <- file.path(target_folder, "S3_selected_K.png")
  grDevices::png(
    filename = png_path,
    width = 8.24, height = 2.62, units = "in", res = 300,
    type = "cairo", bg = "white"
  )
  tryCatch(print(combined_plot), finally = grDevices::dev.off())
  message("Created: ", png_path)
  invisible(list(png = png_path, plot = combined_plot, data = selection_data))
}
