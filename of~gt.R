# ============================================================================
# Quadratic Polynomial Mixed Effects Analysis
# Graph-Theoretical Metrics & Optimal Functioning
# Olivier Brown, 2026-02-26, R version 4.4.1
# ============================================================================

rm(list=ls())

# CONFIGURATION — change these to switch between structural and functional
data_file     <- "analysis_data_struct.csv"      # or "analysis_data.csv"
global_metrics <- c("Sigma", "Eg", "Lambda", "Gamma", "degree")
output_label  <- "structural"                    # or "functional"

# 1. PACKAGES ===============================================================
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, lme4, lmerTest, broom.mixed, rio, effectsize)

# 2. DATA IMPORT & SETUP ====================================================
df <- rio::import(data_file)
if (nrow(df) == 0) stop("Data file is empty.")

# Factor coding
df <- df %>%
  mutate(
    participant = factor(participant),
    grp = factor(grp, levels = c("Orthopedic", "Concussion")),
    sex = factor(sex),
    hemisphere = factor(hemisphere),
    session = factor(session, levels = c("B0", "M3", "M6")),
    metric = factor(metric),
    roi = factor(roi)
  )

# Scale predictors
# Mean-center time (no variance scaling) to reduce collinearity with time²
# Scale metric_value WITHIN each metric × roi to make graph-theory metrics comparable
df <- df %>%
  mutate(
    time_c  = as.numeric(scale(time, center = TRUE, scale = FALSE)),
    time2_c = time_c^2,
    age_z   = as.numeric(scale(age))
  ) %>%
  group_by(metric, roi) %>%
  mutate(metric_value_z = as.numeric(scale(metric_value))) %>%
  ungroup()

# 3. MIXED EFFECTS MODEL FITTING ============================================
fit_models <- function(data, group_filter = NULL, include_group_interaction = TRUE) {
  
  # Subset to single group if requested (Concussion-only analyses)
  if (!is.null(group_filter)) {
    data <- data %>%
      filter(grp == group_filter) %>%
      mutate(grp = droplevels(grp))
  }
  
  # Distinguish global (whole-brain) from nodal (ROI-level) metrics
  nodal_metrics <- setdiff(unique(data$metric), global_metrics)
  results <- list()
  
  # ─ GLOBAL MODELS ────────────────────────────────────────────────────────
  # Whole-brain metrics only; no ROI or hemisphere breakdown
  for (gm in global_metrics) {
    df_gm <- data %>%
      filter(metric == gm, roi == "wb", hemisphere == "both")
    if (nrow(df_gm) == 0) next
    
    # Include 3-way group interaction only when comparing both groups
    if (include_group_interaction && length(unique(data$grp)) > 1) {
      formula <- optimal_functioning ~ metric_value_z * grp * (time_c + time2_c) +
        age_z + sex + (1 | participant)
    } else {
      formula <- optimal_functioning ~ metric_value_z * (time_c + time2_c) +
        age_z + sex + (1 | participant)
    }
    
    model <- lmer(formula, data = df_gm, REML = FALSE)
    results[[paste0(gm, "_global")]] <- list(
      data = df_gm, model = model, type = "global", metric = gm
    )
  }
  
  # ─ NODAL MODELS ─────────────────────────────────────────────────────────
  # Looped over each metric × ROI combination
  for (nm in nodal_metrics) {
    df_nm <- data %>% filter(metric == nm)
    if (nrow(df_nm) == 0) next
    
    for (r in unique(df_nm$roi)) {
      df_roi <- df_nm %>%
        filter(roi == r) %>%
        mutate(
          grp        = droplevels(grp),
          hemisphere = droplevels(hemisphere)
        )
      
      if (include_group_interaction && length(unique(data$grp)) > 1) {
        formula <- optimal_functioning ~ metric_value_z * grp * (time_c + time2_c) +
          age_z + sex + hemisphere + (1 | participant)
      } else {
        formula <- optimal_functioning ~ metric_value_z * (time_c + time2_c) +
          age_z + sex + hemisphere + (1 | participant)
      }
      
      model <- lmer(formula, data = df_roi, REML = FALSE)
      results[[paste0(nm, "_", r, "_nodal")]] <- list(
        data = df_roi, model = model, type = "nodal", metric = nm, roi = r
      )
    }
  }
  
  return(results)
}

# 4. EXTRACT & APPLY BY-METRIC FDR ==========================================
process_results <- function(models_list, analysis_type = "both") {
  
  # Extract fixed-effect estimates
  summary <- models_list %>%
    imap_dfr(function(result, name) {
      fixed <- broom.mixed::tidy(result$model, effects = "fixed")
      parts <- str_split(name, "_")[[1]]
      
      if (result$type == "global") {
        fixed %>% mutate(metric = parts[1], roi = "wb", model_type = "global")
      } else {
        fixed %>% mutate(metric = parts[1], roi = parts[2], model_type = "nodal")
      }
    })
  
  # Flag relevant test terms based on analysis type
  # Metric × time interactions capture differential trajectories and are robust to covariates
  if (analysis_type == "interaction") {
    summary <- summary %>%
      mutate(is_test_term = grepl("metric_value_z:.*:time", term))
  } else if (analysis_type == "main") {
    summary <- summary %>%
      mutate(is_test_term = grepl("^metric_value_z$", term))
  } else {
    summary <- summary %>%
      mutate(is_test_term = grepl("metric_value_z:.*time|^metric_value_z$", term))
  }
  
  # Initialize FDR columns
  summary <- summary %>%
    mutate(p_fdr = NA_real_, sig = FALSE)
  
  # Apply by-metric FDR correction
  # This controls false discovery rate within each graph-theory metric separately,
  # reducing multiple-comparison burden while maintaining interpretability per metric
  for (m in unique(summary$metric)) {
    test_idx <- which(summary$metric == m & summary$is_test_term == TRUE)
    
    if (length(test_idx) > 0) {
      fdr_adj <- p.adjust(summary$p.value[test_idx], method = "BH")
      summary$p_fdr[test_idx] <- fdr_adj
      summary$sig[test_idx]   <- fdr_adj < 0.05
    }
  }
  
  return(summary)
}

# 5. EFFECT SIZE COMPUTATION ================================================
# Cohen's d with 95% CI (pooled SD), computed on model residuals
# Residuals partial out covariates, so effect size reflects group differences
# above and beyond age, sex, and longitudinal trajectories
compute_cohens_d_ci <- function(y1, y2) {
  y1 <- na.omit(y1)
  y2 <- na.omit(y2)
  
  if (length(y1) < 2 || length(y2) < 2) {
    return(data.frame(cohens_d = NA, ci_lower = NA, ci_upper = NA))
  }
  
  d_result <- cohens_d(y1, y2, ci = 0.95)
  
  return(data.frame(
    cohens_d = as.numeric(d_result$Cohens_d),
    ci_lower = as.numeric(d_result$CI_low),
    ci_upper = as.numeric(d_result$CI_high)
  ))
}

# 6. ADD EFFECT SIZES TO SIGNIFICANT FINDINGS ===============================
add_effect_sizes <- function(summary_df, models_list) {
  
  sig_df <- summary_df %>%
    filter(sig == TRUE, grepl("metric_value_z:.*:time", term))
  
  if (nrow(sig_df) == 0) return(sig_df)
  
  effect_sizes <- list()
  counter <- 1
  
  for (i in seq_len(nrow(sig_df))) {
    row          <- sig_df[i, ]
    metric_i     <- row$metric
    roi_i        <- row$roi
    model_type_i <- row$model_type
    
    key <- if (model_type_i == "global") {
      paste0(metric_i, "_global")
    } else {
      paste0(metric_i, "_", roi_i, "_nodal")
    }
    
    if (!(key %in% names(models_list))) next
    
    data_i  <- models_list[[key]]$data
    model_i <- models_list[[key]]$model
    
    # Compute effect size on residuals if model includes group term
    if ("grpConcussion" %in% names(fixef(model_i))) {
      residuals_ortho <- residuals(model_i)[data_i$grp == "Orthopedic"]
      residuals_conc  <- residuals(model_i)[data_i$grp == "Concussion"]
      
      d_result <- compute_cohens_d_ci(residuals_conc, residuals_ortho)
      
      effect_sizes[[counter]] <- data.frame(
        metric     = metric_i,
        roi        = roi_i,
        model_type = model_type_i,
        term       = row$term,
        cohens_d   = d_result$cohens_d,
        ci_lower   = d_result$ci_lower,
        ci_upper   = d_result$ci_upper,
        stringsAsFactors = FALSE
      )
      counter <- counter + 1
    }
  }
  
  if (length(effect_sizes) > 0) {
    effect_df <- bind_rows(effect_sizes)
    sig_df <- sig_df %>%
      left_join(effect_df, by = c("metric", "roi", "model_type", "term"))
  } else {
    sig_df <- sig_df %>%
      mutate(cohens_d = NA, ci_lower = NA, ci_upper = NA)
  }
  
  return(sig_df)
}

# 7. ANALYSIS 1: GROUP × METRIC × TIME INTERACTIONS =========================
cat("\n")
cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("ANALYSIS 1: GROUP × METRIC × TIME INTERACTIONS (Moderation Test)\n")
cat("═══════════════════════════════════════════════════════════════════════════\n")

models_all  <- fit_models(df, group_filter = NULL, include_group_interaction = TRUE)
results_all <- process_results(models_all, analysis_type = "interaction")

sig_findings_1 <- results_all %>%
  filter(sig == TRUE) %>%
  select(model_type, metric, roi, term, estimate, std.error, statistic, p.value, p_fdr)

sig_findings_1 <- add_effect_sizes(results_all, models_all)

cat("\nSignificant three-way interactions (BY-METRIC FDR < 0.05):", nrow(sig_findings_1), "\n")
if (nrow(sig_findings_1) > 0) {
  print(sig_findings_1, n = Inf)
} else {
  cat("None\n")
}

# 8. ANALYSIS 2: CONCUSSION-ONLY ASSOCIATIONS ===============================
cat("\n")
cat("═══════════════════════════════════════════════════════════════════════════\n")
cat("ANALYSIS 2: METRIC EFFECTS (Concussion Group Only)\n")
cat("═══════════════════════════════════════════════════════════════════════════\n")

models_conc  <- fit_models(df, group_filter = "Concussion", include_group_interaction = FALSE)
results_conc <- process_results(models_conc, analysis_type = "both")

sig_main <- results_conc %>%
  filter(sig == TRUE, grepl("^metric_value_z$", term)) %>%
  select(model_type, metric, roi, term, estimate, std.error, statistic, p.value, p_fdr)

sig_interact <- results_conc %>%
  filter(sig == TRUE, grepl("metric_value_z:time", term)) %>%
  select(model_type, metric, roi, term, estimate, std.error, statistic, p.value, p_fdr)

cat("\nMain effects (BY-METRIC FDR < 0.05):", nrow(sig_main), "\n")
if (nrow(sig_main) > 0) print(sig_main, n = Inf) else cat("None\n")

cat("\nMetric × Time interactions (BY-METRIC FDR < 0.05):", nrow(sig_interact), "\n")
if (nrow(sig_interact) > 0) print(sig_interact, n = Inf) else cat("None\n")

cat("\n")
cat("═══════════════════════════════════════════════════════════════════════════\n")

# 9. EXPORT ALL SIGNIFICANT FINDINGS =========================================
all_sig <- bind_rows(
  sig_findings_1 %>% mutate(analysis = "moderation"),
  sig_main       %>% mutate(analysis = "concussion_main"),
  sig_interact   %>% mutate(analysis = "concussion_interaction")
)

write.csv(all_sig, paste(output_label, "_all_significant.csv"), row.names = FALSE)