# ============================================================================
# PCA OF NCp: EXTRACT PC1 AND REFIT MIXED EFFECTS MODELS
# ============================================================================
# Reads analysis_data_struct.csv directly and runs end-to-end.
# Principal Component Analysis confirms that 44 inter-correlated NCp ROIs
# are dominated by a single latent factor. PC1 is extracted and used to refit
# the group moderation model, demonstrating that ROI-level findings reflect a
# spatially distributed structural signal rather than region-specific effects.

# 1. PACKAGES =================================================================
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, lme4, lmerTest, broom.mixed, rio)

# 2. DATA IMPORT & SETUP ======================================================
df <- rio::import("analysis_data_struct.csv")
if (nrow(df) == 0) stop("Data file is empty.")

# Factor coding
df <- df %>%
  mutate(
    participant = factor(participant),
    grp         = factor(grp, levels = c("Orthopedic", "Concussion")),
    sex         = factor(sex),
    hemisphere  = factor(hemisphere),
    session     = factor(session, levels = c("B0", "M3", "M6")),
    metric      = factor(metric),
    roi         = factor(roi)
  )

# Scale predictors
df <- df %>%
  mutate(
    time_c = as.numeric(scale(time, center = TRUE, scale = FALSE)),
    time2_c = time_c^2,
    age_z = as.numeric(scale(age))
  ) %>%
  group_by(metric, roi) %>%
  mutate(
    metric_value_z = as.numeric(scale(metric_value))
  ) %>%
  ungroup()

cat("Data loaded:", nrow(df), "rows,", n_distinct(df$participant), "participants\n")

# 3. BUILD WIDE-FORMAT NCp MATRIX =============================================
# Average across hemispheres (hemisphere is a covariate in the main models,
# not a structural distinction), then pivot to participant × ROI format

ncp_wide <- df %>%
  filter(metric == "NCp") %>%
  group_by(participant, roi, session) %>%
  summarise(metric_value = mean(metric_value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = roi, values_from = metric_value)

cat("Wide NCp matrix:", nrow(ncp_wide), "rows (participant × session),",
    ncol(ncp_wide) - 2, "ROIs\n")

# 4. COMPUTE PCA ON NCp MATRIX ================================================
ncp_matrix <- ncp_wide %>%
  select(-participant, -session) %>%
  as.matrix()

pca_result <- prcomp(ncp_matrix, center = TRUE, scale. = TRUE)

# Variance explained
var_exp <- summary(pca_result)$importance[2, ] * 100
cum_var  <- summary(pca_result)$importance[3, ] * 100

cat("\n=== PCA VARIANCE EXPLAINED ===\n")
cat("PC1:", round(var_exp[1], 1), "%\n")
cat("PC2:", round(var_exp[2], 1), "%\n")
cat("Components explaining >80% variance:", min(which(cum_var >= 80)), "\n")
cat("Components explaining >90% variance:", min(which(cum_var >= 90)), "\n")

# 5. EXTRACT PC1 SCORES AND MERGE =============================================
pc1_scores <- ncp_wide %>%
  select(participant, session) %>%
  mutate(pc1 = pca_result$x[, 1])

# Join PC1 to main dataframe
# Use a single representative NCp row per participant × session to avoid
# duplicating rows across ROIs and hemispheres
df_pc1 <- df %>%
  filter(metric == "NCp") %>%
  select(participant, session, grp, sex, age_z,
         time, time_c, time2_c, optimal_functioning) %>%
  distinct() %>%
  left_join(pc1_scores, by = c("participant", "session")) %>%
  mutate(pc1_z = as.numeric(scale(pc1)))

cat("\nPC1 dataframe:", nrow(df_pc1), "rows,",
    sum(is.na(df_pc1$pc1)), "missing PC1 scores\n")

# 6. FIT GROUP MODERATION MODEL: PC1 × GROUP × TIME ==========================
cat("\n=== FITTING PC1 GROUP MODERATION MODEL ===\n")

model_pc1_moderation <- lmer(
  optimal_functioning ~ pc1_z * grp * (time_c + time2_c) +
    age_z + sex + (1 | participant),
  data = df_pc1,
  REML = FALSE
)

# 7. EXTRACT AND DISPLAY RESULTS ==============================================
cat("\n")
cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("PC1 GROUP MODERATION MODEL: pc1_z × grp × time\n")
cat("═══════════════════════════════════════════════════════════════════════════\n")

results_pc1_mod <- broom.mixed::tidy(model_pc1_moderation, effects = "fixed") %>%
  mutate(
    ci_lower = estimate - 1.96 * std.error,
    ci_upper = estimate + 1.96 * std.error
  ) %>%
  filter(grepl("pc1_z", term)) %>%
  select(term, estimate, std.error, statistic, p.value, ci_lower, ci_upper)

print(results_pc1_mod, n = Inf)

# Highlight key three-way interaction
cat("\n=== KEY TERM: pc1_z × grpConcussion × time_c ===\n")
results_pc1_mod %>%
  filter(term == "pc1_z:grpConcussion:time_c") %>%
  print()

# 8. SIMPLE SLOPES PLOT =======================================================
cat("\nGenerating simple slopes plot...\n")

q25 <- quantile(df_pc1$pc1_z, 0.25, na.rm = TRUE)
q50 <- quantile(df_pc1$pc1_z, 0.50, na.rm = TRUE)
q75 <- quantile(df_pc1$pc1_z, 0.75, na.rm = TRUE)

pc1_levels <- c(q25, q50, q75)
pc1_labels <- c("Low (Q25)", "Medium (Q50)", "High (Q75)")

time_mean    <- mean(df_pc1$time, na.rm = TRUE)
time_seq_raw <- seq(min(df_pc1$time, na.rm = TRUE),
                    max(df_pc1$time, na.rm = TRUE),
                    length.out = 50)

pred_list <- list()

for (g in levels(df_pc1$grp)) {
  for (i in seq_along(pc1_levels)) {
    
    pred_grid <- data.frame(
      pc1_z       = pc1_levels[i],
      time_c      = time_seq_raw - time_mean,
      time2_c     = (time_seq_raw - time_mean)^2,
      grp         = factor(g, levels = levels(df_pc1$grp)),
      age_z       = 0,
      sex         = factor(df_pc1$sex[1], levels = levels(df_pc1$sex)),
      participant = factor(df_pc1$participant[1], levels = levels(df_pc1$participant)),
      pc1_label   = factor(pc1_labels[i], levels = pc1_labels),
      time_raw    = time_seq_raw
    )
    
    tryCatch({
      pred_grid$pred <- predict(model_pc1_moderation, newdata = pred_grid,
                                re.form = NA, allow.new.levels = TRUE)
      pred_list[[paste0(g, "_", i)]] <- pred_grid
    }, error = function(e) {
      cat("Prediction error for", g, pc1_labels[i], ":", e$message, "\n")
    })
  }
}

pred_data <- bind_rows(pred_list)

p_pc1 <- ggplot(pred_data, aes(x = time_raw, y = pred,
                               color = pc1_label, linetype = pc1_label)) +
  geom_line(size = 1.2, alpha = 0.9) +
  geom_point(size = 2.5, alpha = 0.7) +
  facet_wrap(~ grp) +
  scale_color_manual(
    name   = "NCp PC1\n(z-scored)",
    values = c("Low (Q25)"    = "#2166AC",
               "Medium (Q50)" = "#92C5DE",
               "High (Q75)"   = "#B2182B")
  ) +
  scale_linetype_manual(
    name   = "NCp PC1\n(z-scored)",
    values = c("Low (Q25)"    = "solid",
               "Medium (Q50)" = "dashed",
               "High (Q75)"   = "solid")
  ) +
  scale_y_continuous(limits = c(0, 11)) +
  labs(
    title    = "Structural NCp PC1 × Group × Time",
    subtitle = "Three-Way Interaction: PC1 (z) × Group × Time (Simple Slopes)",
    x        = "Days Post-Injury",
    y        = "Predicted Optimal Functioning"
  ) +
  theme_minimal() +
  theme(
    plot.title       = element_text(size = 13, face = "bold"),
    plot.subtitle    = element_text(size = 10, color = "gray50"),
    legend.position  = "right",
    panel.grid.minor = element_blank(),
    strip.text       = element_text(size = 11, face = "bold")
  )

print(p_pc1)

if (!dir.exists("plots")) dir.create("plots")
ggsave("plots/PC1_NCp_threeway_interaction.png", p_pc1,
       width = 12, height = 6, dpi = 300)

# 9. EXPORT ===================================================================
write.csv(results_pc1_mod, "pc1_moderation_results.csv", row.names = FALSE)

cat("\n✓ PC1 analysis complete\n")
cat("Analysis completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")