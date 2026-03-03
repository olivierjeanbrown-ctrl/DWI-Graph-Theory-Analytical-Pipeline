# A-CAP fMRI Graph-Theoretical Analysis Pipeline

Data processing and analysis pipeline for functional graph-theoretical metrics and optimal functioning in pediatric concussion.

## Research Questions

**Primary:** Are graph-theoretical metrics associated with optimal functioning across time in concussion?

**Secondary:** Does injury type moderate these associations?

## Data Processing (SQL)

`analysis_pipeline.sql` using `analysis_pipeline.db` creates `analysis_data.csv` through 5 steps:

1. Reshape optimal functioning to long format (B0, M3, M6)
2. Extract demographics and create time variables
3. Join demographics with outcomes
4. Join with graph-theoretical metrics (weighted, non-null)
5. Create performance indexes

**Columns:** participant, session, time, time2, age, sex, grp, metric, metric_value, roi, hemisphere, optimal_functioning

## Statistical Analysis (R)

### Script 1: `of~funcgt_concussion_v02.R`
Concussion-only analysis of metric × time associations

**Model:** `optimal_functioning ~ metric_z × (time_c + time2_c) + age_z + sex + (1 | participant)`

**Steps:**
1. Standardize predictors
2. Fit LMM for each global metric (Sigma, Eg, Lambda, Gamma)
3. Fit LMM for each nodal metric × ROI
4. Extract significant interactions
5. Generate interpretations and plots

**Output:** significant_interactions_concussion.csv, interpretations CSV, plots

### Script 2: `of~funcgt_v02.R`
Full-sample moderation analysis

**Model:** `optimal_functioning ~ metric_z × grp × (time_c + time2_c) + age_z + sex + (1 | participant)`

**Output:** significant_threeway_interactions.csv, moderation interpretations, faceted plots

## Methods

- Standardization: z-scored (time: centered only)
- FDR Correction: Metric-level families, α = .05
- Random effects: Participant intercepts
- Visualization: Z-scored metrics, raw time scale

## Pipeline

SQL (5 steps) → analysis_data.csv → Script 1 (2-way) + Script 2 (3-way)

