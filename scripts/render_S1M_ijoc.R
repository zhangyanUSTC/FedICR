render_S1M_ijoc <- function(target_folder) {
  file_names <- list.files(
    path = target_folder,
    pattern = "^S1_M_RES_.*\\.Rdata$",
    full.names = TRUE,
    recursive = FALSE
  )
  n_values <- as.numeric(sub(".*_n([0-9]+)_.*", "\\1", basename(file_names)))
  file_names <- file_names[order(n_values, decreasing = TRUE)]
  
  if (length(file_names) != 3L) {
    stop(sprintf("Expected 3 RES files, found %d.", length(file_names)))
  }
  
  plot_mse_ijoc <- function(iter_my, iter_non, iter_sin, iter_ave,
                            oracle, name, iters, na.rm = FALSE) {
    plot_data <- data.frame(
      iteration = rep(seq_len(iters), 4L),
      value = c(iter_my, iter_non, iter_sin, iter_ave),
      algorithm = factor(
        rep(c("E-FICR", "CICR", "SC", "AA"), each = iters),
        levels = c("CICR", "E-FICR", "AA", "SC")
      )
    )
    
    p <- ggplot(plot_data, aes(
      x = iteration, y = value, colour = algorithm,
      shape = algorithm, linetype = algorithm
    )) +
      geom_line(linewidth = 0.55, na.rm = na.rm) +
      geom_point(size = 1.3, stroke = 0.32, na.rm = na.rm) +
      labs(x = NULL, y = name, colour = NULL, shape = NULL, linetype = NULL) +
      scale_colour_manual(values = c(
        "CICR" = "#4477AA", "E-FICR" = "#EE6677",
        "AA" = "#7B6BAF", "SC" = "#66A99A"
      )) +
      scale_shape_manual(values = c("CICR" = 17, "E-FICR" = 15, "AA" = 16, "SC" = 18)) +
      scale_linetype_manual(values = c("CICR" = "solid", "E-FICR" = "dashed",
                                       "AA" = "dotdash", "SC" = "longdash")) +
      scale_x_continuous(breaks = seq(5, iters, by = 5), expand = expansion(mult = c(0.02, 0.02))) +
      theme_bw(base_size = 14.1, base_family = "sans") +
      theme(
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "#E6E6E6", linewidth = 0.28),
        panel.border = element_rect(colour = "#555555", linewidth = 0.35),
        axis.text = element_text(colour = "#222222", size = 9.84),
        axis.title.y = element_text(colour = "#111111", size = 12, margin = margin(r = 3)),
        axis.ticks = element_line(linewidth = 0.3),
        axis.ticks.length = unit(1.5, "mm"),
        legend.position = "bottom",
        legend.key.width = unit(12, "mm"),
        legend.key.height = unit(3.2, "mm"),
        legend.text = element_text(size = 13.8),
        plot.margin = margin(2.5, 3.5, 2.5, 2.5)
      )
    
    if (!is.null(oracle)) {
      p <- p + geom_hline(
        yintercept = oracle, colour = "#D89000",
        linetype = "dashed", linewidth = 0.55
      )
    }
    p
  }
  
  all_plots <- list()
  plot_counter <- 1L
  legend <- NULL
  
  for (file_path in file_names) {
    saved <- new.env(parent = emptyenv())
    load(file_path, envir = saved)
    res_list <- E_FICR_result_names(saved$res_list)
    Oracle_res <- saved$Oracle_res
    K <- saved$K
    M <- saved$M
    iters <- saved$iters
    
    iter_list <- lapply(res_list, function(res) {
      it <- mean_index(
        res = res[res[, 2] == K, c(1:6, 15:18), drop = FALSE],
        iters
      )
      gap <- it[, 11] - tail(it[, 11], 1L)
      positive <- gap > 0
      it[, 11] <- NA_real_
      it[positive, 11] <- log(gap[positive])
      it
    })
    oracle_index <- colMeans(Oracle_res[, c(3, 4, 13:16), drop = FALSE])
    
    metric_idx <- c(1, 2, 6)
    metric_names <- c("MSE", "Acc", "Ell")
    
    for (idx in seq_along(metric_idx)) {
      i <- metric_idx[idx]
      m_name <- metric_names[idx]
      curr_iters <- if (i == 6) min(20, iters) else min(25, iters)
      curr_oracle <- if (i == 6) NULL else oracle_index[i]
      dynamic_name <- bquote(.(as.name(m_name)) ~ scriptstyle((M == .(M))))
      
      p_tmp <- plot_mse_ijoc(
        iter_my = iter_list$E_FICR[seq_len(curr_iters), i + 5],
        iter_non = iter_list$KR[seq_len(curr_iters), i + 5],
        iter_sin = iter_list$SMA[seq_len(curr_iters), i + 5],
        iter_ave = iter_list$AA[seq_len(curr_iters), i + 5],
        oracle = curr_oracle,
        name = dynamic_name,
        iters = curr_iters,
        na.rm = i == 6L
      )
      
      all_plots[[plot_counter]] <- p_tmp + theme(legend.position = "none")
      plot_counter <- plot_counter + 1L
    }
  }
  
  legend_levels <- c("Oracle", "CICR", "E-FICR", "AA", "SC")
  legend_data <- data.frame(
    x = rep(1:2, length(legend_levels)),
    y = rep(seq_along(legend_levels), each = 2),
    algorithm = factor(rep(legend_levels, each = 2), levels = legend_levels)
  )
  legend_plot <- ggplot(
    legend_data,
    aes(x = x, y = y, colour = algorithm, linetype = algorithm, group = algorithm)
  ) +
    geom_line(linewidth = 0.55) +
    geom_point(
      data = subset(legend_data, algorithm != "Oracle"),
      aes(shape = algorithm), size = 1.3, stroke = 0.32
    ) +
    scale_colour_manual(values = c(
      "Oracle" = "#D89000", "CICR" = "#4477AA", "E-FICR" = "#EE6677",
      "AA" = "#7B6BAF", "SC" = "#66A99A"
    )) +
    scale_linetype_manual(values = c(
      "Oracle" = "dashed", "CICR" = "solid", "E-FICR" = "dashed",
      "AA" = "dotdash", "SC" = "longdash"
    )) +
    scale_shape_manual(values = c(
      "Oracle" = NA, "CICR" = 17, "E-FICR" = 15, "AA" = 16, "SC" = 18
    )) +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.key.width = unit(16, "mm"),
      legend.key.height = unit(3.2, "mm"),
      legend.text = element_text(size = 13.8)
    ) +
    guides(
      colour = guide_legend(
        nrow = 1, title = NULL,
        override.aes = list(
          linewidth = 0.8,
          shape = c(NA, 17, 15, 16, 18),
          linetype = c("dashed", "solid", "dashed", "dotdash", "longdash")
        )
      ),
      linetype = "none",
      shape = "none"
    )
  legend <- simulation_legend(legend_plot, position = "bottom")
  
  panel <- cowplot::plot_grid(
    plotlist = all_plots,
    ncol = 3,
    nrow = 3,
    align = "hv",
    axis = "tblr"
  )
  
  combined_plot <- cowplot::plot_grid(
    panel,
    legend,
    ncol = 1,
    rel_heights = c(1, 0.055)
  )
  
  png_path <- file.path(target_folder, "S1M_combined.png")
  
  ggsave(
    filename = png_path,
    plot = combined_plot,
    width = 8.24,
    height = 5.04,
    units = "in",
    dpi = 300,
    bg = "white"
  )
  
  message("Created: ", png_path)
  invisible(c(png = png_path))
  
}
