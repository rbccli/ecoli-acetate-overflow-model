library(dplyr)
library(tidyr)
library(ggplot2)

# -------------------------------
# 0. Parameters to extract / plot
# -------------------------------

raw_params <- c(
  "pdh_raw", "akd_raw", "scs_raw", "aco_raw",
  "sd_raw", "pk_raw", "cs_raw", "id_raw",
  "pta_raw", "md_raw",
  "f_raw", "ack_raw","alpha_raw","k_raw"
)

# 如果你之后 posterior 里也有 f_raw 和 ack_raw，
# 可以改成下面这个：
# raw_params <- c(
#   "pdh_raw", "akd_raw", "scs_raw", "aco_raw",
#   "sd_raw", "pk_raw", "cs_raw", "id_raw",
#   "pta_raw", "md_raw", "acetate_raw",
#   "f_raw", "ack_raw",
#   "k_raw", "alpha_raw"
# )


# -------------------------------
# 1. Posterior samples
# -------------------------------

posterior_raw <- tca_ferm_fit_real$draws(
  variables = raw_params,
  format = "df"
)

df_post_raw <- posterior_raw %>%
  pivot_longer(
    cols = all_of(raw_params),
    names_to = "param_name",
    values_to = "value"
  ) %>%
  mutate(type = "Posterior")


# -------------------------------
# 2. Prior specification
# -------------------------------

# -------------------------------
# 2. Prior specification
# -------------------------------

prior_spec <- tibble::tribble(
  ~param_name,     ~prior_mean,                 ~prior_sd,
  
  "pdh_raw",       -2,                           1,
  "aco_raw",       -2,                           1,
  "akd_raw",       -2,                           1,
  #"scs_raw",       -2,                           1,
  "scs_raw",       log10(0.0061704935),                       0.5,
  "sd_raw",        -2,                           1,
  "pk_raw",        -2,                           1,
  "cs_raw",        -2,                           1,
  "id_raw",        -2,                           1,
  "md_raw",        log10(0.023603447),           0.5,
  
  "f_raw",         log10(0.001278269),           0.5,
  #"ack_raw",       log10(0.002828854),           0.3,
  
  #"md_raw",        -2,                           1,
  
  #"f_raw",         -2,                           1,
  "ack_raw",       -2,                           1,
  
 # "pta_raw",       -2,                           1,
  "pta_raw",       log10(0.000887043),                           0.5,
  
  "k_raw",          2,                   0.5,
  "alpha_raw",      log10(0.0545/2.8),           0.5
  #"k_raw",          1,                           1,
  #"alpha_raw",      -2,                          1
)



# check whether every plotted parameter has a prior
missing_prior <- setdiff(raw_params, prior_spec$param_name)

if (length(missing_prior) > 0) {
  stop(
    "Missing prior specification for: ",
    paste(missing_prior, collapse = ", ")
  )
}

set.seed(1)

df_prior <- lapply(raw_params, function(p) {
  
  spec <- prior_spec %>%
    filter(param_name == p)
  
  data.frame(
    param_name = p,
    value = rnorm(
      5000,
      mean = spec$prior_mean,
      sd   = spec$prior_sd
    ),
    type = "Prior"
  )
  
}) %>% bind_rows()


# -------------------------------
# 3. Combine prior + posterior
# -------------------------------

df_all <- bind_rows(df_prior, df_post_raw)


# -------------------------------
# 4. True values
# -------------------------------

true_spec <- tibble::tribble(
  ~param_name,      ~true_value,
  
  "pdh_raw",        -1,
  "akd_raw",        -1,
  "scs_raw",        -1,
  "aco_raw",        -1,
  "sd_raw",         -1,
  "pk_raw",         -1,
  "cs_raw",         -1,
  "id_raw",         -1,
  "pta_raw",        -1,
  "md_raw",         -1,
  "f_raw",          -1,
  "ack_raw",        -1,
  # "acetate_raw",    log10(0.005),
  "alpha_raw",      log10(0.0545 / 2.8),
  "k_raw",              2
  
  # 如果之后加 f_raw 和 ack_raw，也可以在这里加 true value。
  # "f_raw",       ...,
  # "ack_raw",     ...
)

df_true <- true_spec %>%
  filter(param_name %in% raw_params)

df_ci <- df_post_raw %>%
  group_by(param_name) %>%
  summarise(
    ci_lower = quantile(value, 0.05),
    ci_upper = quantile(value, 0.95),
    median_value = quantile(value, 0.5),
    .groups = "drop"
  )


# -------------------------------
# 5. Facet labels
# -------------------------------

param_labels <- c(
  pdh_raw     = "E[PDH]",
  akd_raw     = "E[AKGDH]",
  scs_raw     = "E[SCS]",
  sd_raw      = "E[SDH]",
  pk_raw      = "E[PK]",
  cs_raw      = "E[CS]",
  id_raw      = "E[IDH]",
  pta_raw     = "E[PTA]",
  aco_raw     = "E[ACO]",
  md_raw      = "E[MDH]",
  
  f_raw       = "E[FUM]",
  ack_raw     = "E[ACKA]",
  
  #acetate_raw = '"[Acetate]"',
  alpha_raw   = "alpha",
  k_raw       = "k"
)


# -------------------------------
# 6. Plot
# -------------------------------

prior_fill <- "#BDBDBD"

posterior_fill <- "#6BAED6"

ci_colour <- "#2171B5"


# ============================================================
# Panel-specific vertical position for median and 95% CI
# ============================================================

density_heights <- df_all %>%
  filter(
    is.finite(value)
  ) %>%
  group_by(
    param_name,
    type
  ) %>%
  summarise(
    density_max = max(
      stats::density(
        value,
        adjust = 1
      )$y
    ),
    .groups = "drop"
  ) %>%
  group_by(
    param_name
  ) %>%
  summarise(
    density_max = max(density_max),
    .groups = "drop"
  )

df_ci_plot <- df_ci %>%
  left_join(
    density_heights,
    by = "param_name"
  ) %>%
  mutate(
    # Place the CI bar slightly above the density baseline
    ci_y = 0.06 * density_max
  )



# ============================================================
# 11. Construct the main plot
# ============================================================

prior_posterior_base <- ggplot() +
  
  # Prior density
  geom_density(
    data = df_all %>%
      filter(type == "Prior"),
    mapping = aes(
      x = value,
      y = after_stat(density)
    ),
    fill = scales::alpha(
      prior_fill,
      0.30
    ),
    colour = prior_fill,
    linewidth = 0.35,
    adjust = 1
  ) +
  
  # Posterior density
  geom_density(
    data = df_all %>%
      filter(type == "Posterior"),
    mapping = aes(
      x = value,
      y = after_stat(density)
    ),
    fill = scales::alpha(
      posterior_fill,
      0.45
    ),
    colour = posterior_fill,
    linewidth = 0.45,
    adjust = 1
  ) +
  
  # True/reference value
  #  geom_vline(
    #    data = df_true,
    #    mapping = aes(
      #      xintercept = true_value
      #    ),
    #    colour = "black",
  #  linetype = "dashed",
  #  linewidth = 0.45
  #) +
  
  # Posterior 95% credible interval
# Posterior 95% credible interval
geom_segment(
  data = df_ci_plot,
  mapping = aes(
    x = ci_lower,
    xend = ci_upper,
    y = ci_y,
    yend = ci_y
  ),
  colour = ci_colour,
  linewidth = 1.2,
  lineend = "round",
  inherit.aes = FALSE
) +
  
  # Posterior median
  geom_point(
    data = df_ci_plot,
    mapping = aes(
      x = median_value,
      y = ci_y
    ),
    colour = ci_colour,
    fill = "white",
    shape = 21,
    stroke = 0.7,
    size = 2.1,
    inherit.aes = FALSE
  ) +
  
  # Four columns produce two empty spaces in the bottom-right
  facet_wrap(
    facets = vars(param_name),
    ncol = 5,
    scales = "free",
    labeller = labeller(
      param_name = as_labeller(
        param_labels,
        label_parsed
      )
    )
  ) +
  
  labs(
    x = expression(log[10] * " parameter value"),
    y = "Density"
  ) +
  
  scale_y_continuous(
    expand = expansion(
      mult = c(0.03, 0.08)
    )
  ) +
  
  coord_cartesian(
    clip = "off"
  ) +
  
  theme_classic(
    base_size = 9,
    base_family = "sans"
  ) +
  
  theme(
    legend.position = "none",
    
    strip.background = element_blank(),
    
    strip.text = element_text(
      size = 9,
      face = "plain",
      colour = "black"
    ),
    
    axis.title = element_text(
      size = 9,
      colour = "black"
    ),
    
    axis.text = element_text(
      size = 7,
      colour = "black"
    ),
    
    axis.ticks = element_line(
      linewidth = 0.35,
      colour = "black"
    ),
    
    axis.line = element_line(
      linewidth = 0.4,
      colour = "black"
    ),
    
    panel.spacing = grid::unit(
      1,
      "lines"
    ),
    
    plot.margin = margin(
      t = 5,
      r = 5,
      b = 5,
      l = 5
    )
  )

# ============================================================
# Manually constructed legend grob
# ============================================================
# ============================================================
# Compact manual legend
# ============================================================

# ============================================================
# Compact manual legend
# ============================================================

manual_legend_grob <- grid::grobTree(
  
  # Prior
  grid::rectGrob(
    x = 0.12,
    y = 0.75,
    width = 0.18,
    height = 0.14,
    gp = grid::gpar(
      fill = scales::alpha(
        prior_fill,
        0.30
      ),
      col = NA
    )
  ),
  
  grid::textGrob(
    label = "Prior",
    x = 0.26,
    y = 0.75,
    just = "left",
    gp = grid::gpar(
      fontsize = 7,
      fontfamily = "sans",
      col = "black"
    )
  ),
  
  # Posterior
  grid::rectGrob(
    x = 0.12,
    y = 0.50,
    width = 0.18,
    height = 0.14,
    gp = grid::gpar(
      fill = scales::alpha(
        posterior_fill,
        0.45
      ),
      col = NA
    )
  ),
  
  grid::textGrob(
    label = "Posterior",
    x = 0.26,
    y = 0.50,
    just = "left",
    gp = grid::gpar(
      fontsize = 7,
      fontfamily = "sans",
      col = "black"
    )
  ),
  
  # Posterior 95% CI
  grid::segmentsGrob(
    x0 = 0.03,
    x1 = 0.21,
    y0 = 0.25,
    y1 = 0.25,
    gp = grid::gpar(
      col = ci_colour,
      lwd = 1.6,
      lineend = "round"
    )
  ),
  
  # Posterior median
  grid::pointsGrob(
    x = 0.12,
    y = 0.25,
    pch = 21,
    size = grid::unit(
      1.7,
      "mm"
    ),
    gp = grid::gpar(
      col = ci_colour,
      fill = "white",
      lwd = 0.8
    )
  ),
  
  grid::textGrob(
    label = "Median and 90% CI",
    x = 0.26,
    y = 0.25,
    just = "left",
    gp = grid::gpar(
      fontsize = 6.5,
      fontfamily = "sans",
      col = "black"
    )
  )
)
# ============================================================
# Convert the ggplot to a gtable
# ============================================================
prior_posterior_gtable <- ggplot2::ggplotGrob(
  prior_posterior_base
)



# Find facet-panel positions
panel_positions <- prior_posterior_gtable$layout[
  grepl(
    "^panel",
    prior_posterior_gtable$layout$name
  ),
]


panel_rows <- sort(
  unique(
    panel_positions$t
  )
)

panel_columns <- sort(
  unique(
    panel_positions$l
  )
)


# Bottom facet row
legend_top <- max(
  panel_rows
)

legend_bottom <- max(
  panel_positions$b[
    panel_positions$t == legend_top
  ]
)


# Span the two unused bottom-right facet columns
# ============================================================
# Put the legend in the single empty bottom-right facet
# ============================================================

legend_top <- max(
  panel_rows
)

legend_bottom <- max(
  panel_positions$b[
    panel_positions$t == legend_top
  ]
)

# Use only the final empty column
legend_left <- max(
  panel_columns
)

legend_right <- legend_left

# ============================================================
# Add the manual legend to the empty bottom-right facet area
# ============================================================

prior_posterior_fig <- gtable::gtable_add_grob(
  prior_posterior_gtable,
  grobs = manual_legend_grob,
  t = legend_top,
  b = legend_bottom,
  l = legend_left,
  r = legend_right,
  z = Inf,
  clip = "off",
  name = "manual-legend"
)


# ============================================================
# Display final figure
# ============================================================

grid::grid.newpage()
grid::grid.draw(
  prior_posterior_fig
)

prior_posterior_fig

save_journal_plot(
  prior_posterior_fig,
  "tcaferm_real_prior_posterior_fig_final1024",
  width = 20,
  height = 10
)






