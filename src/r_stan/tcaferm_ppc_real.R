library(dplyr)
library(tidyr)
library(ggplot2)
library(posterior)
library(patchwork)
library(scales)
library(grid)


# ============================================================
# Observed data
# ============================================================

y_obs <- c(
  pyruvate       = 3.66,
  oxaloacetate   = 0.000487,
  acetyl_coa     = 0.606,
  citrate        = 1.96,
  cis_aconitate  = 0.0161,
  isocitrate     = 0.0367,
  aketo          = 0.443,
  succinyl_coa   = 0.233,
  succinate      = 0.569,
  fumarate       = 0.115,
  
  # IMPORTANT:
  # If the experimental value is 1.07 µM, this should be
  # 0.00107 mM rather than 1.07 mM.
  acetyl_p       = 1.07,
  
  nad            = 2.8,
  nadh           = 0.0545
)

f_obs_ferm <- 1.0833
f_obs_tca  <- 0.3


met_labels <- c(
  "Pyruvate",
  "Oxaloacetate",
  "Acetyl-CoA",
  "Citrate",
  "cis-Aconitate",
  "Isocitrate",
  "\u03b1-Ketoglutarate",
  "Succinyl-CoA",
  "Succinate",
  "Fumarate",
  "Acetyl-phosphate",
  "NAD\u207A",
  "NADH"
)

flux_labels <- c(
  "vPK", "vPDH", "vCS", "vACO1", "vACO2",
  "vIDH", "vAKGDH", "vSCS", "vSDH", "vFUM",
  "vMDH", "vPTA", "vACK"
)


draws <- tca_ferm_fit_real$draws()


# ============================================================
# Final journal theme
# ============================================================

theme_ppc_final <- function(base_size = 8.5,
                            base_family = "sans") {
  
  theme_classic(
    base_size = base_size,
    base_family = base_family
  ) +
    theme(
      axis.title = element_text(
        size = 9,
        colour = "black"
      ),
      
      axis.text = element_text(
        size = 7.5,
        colour = "black"
      ),
      
      axis.line = element_line(
        linewidth = 0.35,
        colour = "black"
      ),
      
      axis.ticks = element_line(
        linewidth = 0.30,
        colour = "black"
      ),
      
      axis.ticks.length = unit(1.5, "mm"),
      
      strip.background = element_blank(),
      
      strip.text = element_text(
        size = 8.5,
        face = "plain",
        colour = "black",
        margin = margin(b = 2.5)
      ),
      
      panel.spacing.x = unit(0.48, "cm"),
      panel.spacing.y = unit(0.40, "cm"),
      
      plot.margin = margin(
        t = 4,
        r = 5,
        b = 4,
        l = 5
      )
    )
}


# ============================================================
# Select a visually sensible common tick interval
# ============================================================

nice_tick_step <- function(target_step) {
  
  magnitude <- 10^floor(log10(target_step))
  
  candidates <- magnitude * c(
    1, 2, 2.5, 5, 10
  )
  
  candidates[
    which.min(abs(candidates - target_step))
  ]
}


# ============================================================
# Reusable PPC histogram function
# ============================================================

# ============================================================
# Format posterior summaries on the original scale
# ============================================================

format_ppc_number <- function(x) {
  
  vapply(
    x,
    function(z) {
      
      if (!is.finite(z)) {
        return("NA")
      }
      
      az <- abs(z)
      
      if (az > 0 && az < 1e-3) {
        formatC(
          z,
          format = "e",
          digits = 1
        )
        
      } else if (az < 0.01) {
        formatC(
          z,
          format = "f",
          digits = 4
        )
        
      } else if (az < 0.1) {
        formatC(
          z,
          format = "f",
          digits = 3
        )
        
      } else if (az < 10) {
        formatC(
          z,
          format = "f",
          digits = 2
        )
        
      } else {
        formatC(
          z,
          format = "f",
          digits = 1
        )
      }
    },
    character(1)
  )
}


make_ppc_histogram <- function(
    draws,
    prefix,
    panel_labels,
    observed,
    x_title,
    ncol = 5,
    interval_prob = 0.90,
    display_prob = 0.998,
    number_of_bins = 30,
    posterior_colour = "#9ECAE1",
    interval_colour = "#4292C6",
    observed_colour = "#D55E00") {
  
  n_panels <- length(panel_labels)
  
  variable_names <- sprintf(
    "%s[%d]",
    prefix,
    seq_len(n_panels)
  )
  
  draws_df <- posterior::as_draws_df(draws)
  
  missing_variables <- setdiff(
    variable_names,
    names(draws_df)
  )
  
  if (length(missing_variables) > 0) {
    stop(
      "The following variables are missing from the posterior draws: ",
      paste(missing_variables, collapse = ", ")
    )
  }
  
  # ----------------------------------------------------------
  # Posterior draws in long format
  # ----------------------------------------------------------
  
  posterior_long <- draws_df |>
    dplyr::select(all_of(variable_names)) |>
    mutate(draw = row_number()) |>
    pivot_longer(
      cols = all_of(variable_names),
      names_to = "variable",
      values_to = "value"
    ) |>
    mutate(
      panel_id = match(variable, variable_names),
      
      panel = factor(
        panel_labels[panel_id],
        levels = panel_labels
      )
    ) |>
    filter(
      is.finite(value),
      value > 0
    ) |>
    mutate(
      log_value = log10(value)
    )
  
  # ----------------------------------------------------------
  # Observed values
  # ----------------------------------------------------------
  
  observed_df <- tibble(
    panel_id = seq_len(n_panels),
    observed = as.numeric(observed)
  ) |>
    mutate(
      panel = factor(
        panel_labels[panel_id],
        levels = panel_labels
      ),
      
      log_observed = log10(observed)
    ) |>
    filter(
      is.finite(observed),
      observed > 0
    )
  
  # ----------------------------------------------------------
  # Central 95% posterior predictive intervals
  # ----------------------------------------------------------
  
  
  interval_tail <- (
    1 - interval_prob
  ) / 2
  
  interval_df <- posterior_long |>
    group_by(panel) |>
    summarise(
      interval_low = quantile(
        log_value,
        interval_tail,
        na.rm = TRUE
      ),
      
      posterior_median = median(
        log_value,
        na.rm = TRUE
      ),
      
      interval_high = quantile(
        log_value,
        1 - interval_tail,
        na.rm = TRUE
      ),
      
      .groups = "drop"
    ) |>
    mutate(
      # Back-transform to the original scale
      interval_low_original =
        10^interval_low,
      
      posterior_median_original =
        10^posterior_median,
      
      interval_high_original =
        10^interval_high,
      
      numerical_summary = paste0(
        format_ppc_number(
          posterior_median_original
        ),
        " [",
        format_ppc_number(
          interval_low_original
        ),
        "\u2013",
        format_ppc_number(
          interval_high_original
        ),
        "]"
      ),
      
      strip_label = paste0(
        as.character(panel),
        "\n",
        numerical_summary
      )
    )
  
  facet_labels <- setNames(
    interval_df$strip_label,
    as.character(interval_df$panel)
  )
  
  # ----------------------------------------------------------
  # Display range
  #
  # Use the central 99.8% of the draws to prevent a small
  # number of extreme draws from setting the complete x-range.
  # The observed value is always included.
  # ----------------------------------------------------------
  
  display_tail <- (1 - display_prob) / 2
  
  range_df <- posterior_long |>
    group_by(panel) |>
    summarise(
      posterior_low = quantile(
        log_value,
        display_tail,
        na.rm = TRUE
      ),
      
      posterior_high = quantile(
        log_value,
        1 - display_tail,
        na.rm = TRUE
      ),
      
      .groups = "drop"
    ) |>
    left_join(
      observed_df |>
        select(panel, log_observed),
      by = "panel"
    ) |>
    mutate(
      range_low = pmin(
        posterior_low,
        log_observed
      ),
      
      range_high = pmax(
        posterior_high,
        log_observed
      ),
      
      panel_midpoint = (
        range_low + range_high
      ) / 2,
      
      required_width = (
        range_high - range_low
      )
    )
  
  # Same log10 width for every panel
  common_width <- max(
    range_df$required_width,
    na.rm = TRUE
  ) * 1.00
  
  if (!is.finite(common_width) || common_width <= 0) {
    stop("Unable to calculate a valid plotting range.")
  }
  
  # Same major-tick spacing for every panel
  tick_step <- nice_tick_step(
    common_width / 4
  )
  
  x_break_function <- function(limits) {
    
    first_break <- ceiling(
      limits[1] / tick_step
    ) * tick_step
    
    last_break <- floor(
      limits[2] / tick_step
    ) * tick_step
    
    if (first_break > last_break) {
      return(numeric(0))
    }
    
    seq(
      first_break,
      last_break,
      by = tick_step
    )
  }
  
  bounds_df <- range_df |>
    mutate(
      xmin = panel_midpoint - common_width / 2,
      xmax = panel_midpoint + common_width / 2
    ) |>
    select(
      panel,
      xmin,
      xmax
    )
  
  # Clip only the extreme plotting tails
  posterior_clipped <- posterior_long |>
    left_join(
      bounds_df,
      by = "panel"
    ) |>
    filter(
      log_value >= xmin,
      log_value <= xmax
    )
  
  # Invisible anchors force every facet to have the same width
  xlim_df <- bounds_df |>
    pivot_longer(
      cols = c(xmin, xmax),
      names_to = "bound",
      values_to = "log_value"
    ) |>
    mutate(y = 0)
  
  histogram_binwidth <- (
    common_width / number_of_bins
  )
  
  # ----------------------------------------------------------
  # Plot
  # ----------------------------------------------------------
  
  ggplot(
    posterior_clipped,
    aes(x = log_value)
  ) +
    
    geom_blank(
      data = xlim_df,
      aes(
        x = log_value,
        y = y
      ),
      inherit.aes = FALSE
    ) +
    
    geom_rect(
      data = interval_df,
      aes(
        xmin = interval_low,
        xmax = interval_high,
        ymin = 0,
        ymax = Inf
      ),
      inherit.aes = FALSE,
      fill = interval_colour,
      alpha = 0.11
    ) +
    
    geom_histogram(
      aes(
        y = after_stat(density)
      ),
      binwidth = histogram_binwidth,
      boundary = 0,
      fill = posterior_colour,
      alpha = 0.80,
      colour = NA
    ) +
    
    # Posterior predictive median
    geom_vline(
      data = interval_df,
      aes(
        xintercept = posterior_median
      ),
      inherit.aes = FALSE,
      linewidth = 0.50,
      linetype = "dashed",
      colour = interval_colour
    ) +
    
    geom_vline(
      data = observed_df,
      aes(
        xintercept = log_observed
      ),
      linewidth = 0.65,
      colour = observed_colour
    ) +
    
    facet_wrap(
      vars(panel),
      scales = "free_x",
      ncol = ncol,
      drop = FALSE
    ) +
    
    scale_x_continuous(
      breaks = x_break_function,
      labels = label_number(
        accuracy = tick_step,
        trim = TRUE
      ),
      expand = expansion(mult = 0)
    ) +
    
    scale_y_continuous(
      limits = c(0, NA),
      breaks = breaks_extended(n = 3),
      labels = label_number(
        accuracy = 0.1,
        trim = TRUE
      ),
      expand = expansion(
        mult = c(0, 0.08)
      )
    ) +
    
    labs(
      x = x_title,
      y = "Posterior predictive density"
    ) +
    
    theme_ppc_final()
}


# ============================================================
# (a) Metabolite PPC
# ============================================================

p_met_ppc <- make_ppc_histogram(
  draws = draws,
  prefix = "y_rep",
  panel_labels = met_labels,
  observed = y_obs,
  
  x_title = expression(
    log[10] * " concentration (mM)"
  ),
  
  ncol = 5
)


# ============================================================
# (b) Flux PPC
# ============================================================

f_obs <- c(
  f_obs_ferm + f_obs_tca,  # vPK
  f_obs_ferm + f_obs_tca,  # vPDH
  
  rep(
    f_obs_tca,
    9
  ),                       # vCS to vMDH
  
  f_obs_ferm,              # vPTA
  f_obs_ferm               # vACK
)

stopifnot(
  length(f_obs) == length(flux_labels)
)

p_flux_ppc <- make_ppc_histogram(
  draws = draws,
  prefix = "f_rep",
  panel_labels = flux_labels,
  observed = f_obs,
  
  x_title = expression(
    log[10] * " flux (mM " * s^{-1} * ")"
  ),
  
  ncol = 5
)


# ============================================================
# Combined figure
# ============================================================

final_fig <- (
  p_met_ppc /
    p_flux_ppc
) +
  plot_layout(
    heights = c(1, 1)
  ) +
  plot_annotation(
    tag_levels = "A",
    tag_prefix = "(",
    tag_suffix = ")"
  )

final_fig <- final_fig &
  theme(
    plot.tag = element_text(
      size = 10,
      colour = "black"
    ),
    
    plot.tag.position = c(0, 1)
  )

final_fig


save_journal_plot(p_met_ppc, "tcaferm_real_met_ppc_final1024", width = 20, height = 10)
save_journal_plot(p_flux_ppc, "tcaferm_real_flux_ppc_final1024", width = 20, height = 10)
save_journal_plot(final_fig, "tcaferm_real_final_fig_final1024", width = 20, height = 15)
