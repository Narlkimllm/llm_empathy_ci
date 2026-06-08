library(lme4)       # lmer
library(lmerTest)   # Satterthwaite df + p-values for lmer
library(dplyr)      # data wrangling
library(car)        # vif()
library(ggplot2)    # coefficient plot
library(tidyverse)
library(readxl)
setwd("C:/Users/user/Desktop/초거대AI/manuscript/publication/SCI_REP/upload")
df <- read_excel("df_pathology_magnitude.xlsx")
df <- read_excel("df_pathology_redundancy.xlsx")
cat("Data shape:", nrow(df), "rows,", ncol(df), "cols\n")
cat("Participants:", n_distinct(df$participant), "\n")
cat("Emotions    :", n_distinct(df$emotion), "\n")

# =============================================================================
# Start with magnitude outcome, then change df <- read_excel("df_pathology_redundancy.xlsx")
# 1. Standardize predictors
# =============================================================================

predictors <- c("k_maia", "phq9", "gad7", "stai_x1", "erq", "pfq2")

df <- df %>%
  mutate(across(all_of(predictors), ~ as.numeric(scale(.)), .names = "{.col}_z"))
#change magnitude_z for redundancy_z
df <- df %>%
  mutate(
    across(all_of(predictors),
           ~ as.numeric(scale(.)),
           .names = "{.col}_z"),
    redundancy_z  = as.numeric(scale(redundancy))
    )#,
#    redundancy_z = as.numeric(scale(redundancy))
#  )

# =============================================================================
# 2. VIF check  (on participant-level data, one row per participant)
# =============================================================================

df_participant <- df %>%
  select(participant, all_of(paste0(predictors, "_z"))) %>%
  distinct()

lm_vif <- lm(
  k_maia_z ~ phq9_z + gad7_z + stai_x1_z + erq_z + pfq2_z,
  data = df_participant
)

predictors_z <- c("k_maia_z", "phq9_z", "gad7_z", "stai_x1_z", "erq_z", "pfq2_z")

vif_results <- data.frame(predictor = predictors, VIF = NA)
for (i in seq_along(predictors_z)) {
  pred <- predictors_z[i]
  
  formula_reduced <- as.formula(paste(
    pred, "~",
    paste(setdiff(predictors_z, pred), collapse = " + ")
  ))
  
  lm_vif <- lm(formula_reduced, data = df_participant)
  
  r2   <- summary(lm_vif)$r.squared
  vif  <- 1 / (1 - r2)
  
  vif_results$VIF[i] <- round(vif, 3)
}

cat("\n[VIF — Multicollinearity Check]\n")
print(vif_results)
cat("\nVIF > 10: severe | VIF 5-10: moderate | VIF < 5: acceptable\n")


# =============================================================================
# # Start with magnitude outcome, then change magnitude_z to redundancy_z
# 3. Fit LMM
# =============================================================================

model <- lmer(
  redundancy_z ~ k_maia_z + phq9_z + gad7_z + stai_x1_z + erq_z + pfq2_z
  + (1 | participant)
  + (1 | emotion),
  data = df,
  REML = FALSE   # REML for variance component estimation
)

model_min1 <- lmer(
  redundancy_z ~ k_maia_z 
  + (1 | participant)
  + (1 | emotion),
  data = df,
  REML = FALSE   # REML for variance component estimation
)

model_min2 <- lmer(
  redundancy_z ~ phq9_z
  + (1 | participant)
  + (1 | emotion),
  data = df,
  REML = FALSE
)

model_min3 <- lmer(
  redundancy_z ~ gad7_z 
  + (1 | participant)
  + (1 | emotion),
  data = df,
  REML = FALSE
)

model_min4 <- lmer(
  redundancy_z ~ stai_x1_z
  + (1 | participant)
  + (1 | emotion),
  data = df,
  REML = FALSE
)
model_min5 <- lmer(
  redundancy_z ~ erq_z
  + (1 | participant)
  + (1 | emotion),
  data = df,
  REML = FALSE
)
model_min6 <- lmer(
  redundancy_z ~ pfq2_z
  + (1 | participant)
  + (1 | emotion),
  data = df,
  REML = FALSE
)

predictors_z <- c("phq9_z","gad7_z")
for (pred in predictors_z) {
  formula_reduced <- as.formula(paste(
    "redundancy_z ~",
    paste(setdiff(predictors_z, pred), collapse = " + "),
    "+ (1 | participant) + (1 | emotion)"
  ))
  
  model_reduced <- lmer(formula_reduced, data = df, REML = FALSE)
  AIC(model_reduced)
  
  cat("\n[Dropped:", pred, "]\n")
  print(AIC(model_reduced))
  #print(anova(model_reduced, model_min1))
}
model_comp <- anova(model_min6 , model_min5, model_min4, model_min3, model_min2, model_min1, model)
cat("\n[Model Comparisons]\n")
print(model_comp)
cat("\n[Model Summary]\n")
print(summary(model))
# =============================================================================
# 4. Random effects variance components
# =============================================================================
cat("\n[Random Effects Variance Components]\n")
print(as.data.frame(VarCorr(model)))

# =============================================================================
# 5. Manual marginal effect computation
#    - Vary one predictor across its range
#    - Hold all others at 0 (= mean, since z-scored)
#    - Compute predicted magnitude + 95% CI from fixed effects only
# =============================================================================
fe       <- fixef(model)
vcov_mat <- as.matrix(vcov(model))

predict_marginal <- function(predictor_name, n_points = 100) {
  
  # Range: -2.5 SD to +2.5 SD (covers ~99% of z-scored data)
  z_seq <- seq(-2.5, 2.5, length.out = n_points)
  
  all_predictors_z <- c("k_maia_z", "phq9_z", "gad7_z", "stai_x1_z", "erq_z", "pfq2_z")
  
  preds <- lapply(z_seq, function(z_val) {
    
    # Design vector: intercept + all predictors at 0 except focal
    x <- setNames(rep(0, length(all_predictors_z)), all_predictors_z)
    x[predictor_name] <- z_val
    x_vec <- c("(Intercept)" = 1, x)
    
    # Predicted value
    y_hat <- sum(fe * x_vec)
    
    # SE from delta method: sqrt(x' * vcov * x)
    se <- sqrt(as.numeric(t(x_vec) %*% vcov_mat %*% x_vec))
    
    data.frame(
      z_val  = z_val,
      y_hat  = y_hat,
      lower  = y_hat - 1.96 * se,
      upper  = y_hat + 1.96 * se
    )
  })
  
  bind_rows(preds) %>% mutate(predictor = predictor_name)
}

# =============================================================================
# 6. Compute marginal effects for all predictors
# =============================================================================

marginal_df <- bind_rows(lapply(
  c("k_maia_z", "phq9_z", "gad7_z", "stai_x1_z", "erq_z", "pfq2_z"),
  predict_marginal
))

# Predictor labels for plot
predictor_labels <- c(
  k_maia_z   = "K-MAIA (z)",
  phq9_z     = "PHQ-9 (z)",
  gad7_z     = "GAD-7 (z)",
  stai_x1_z  = "STAI-X1 (z)",
  erq_z      = "ERQ (z)",
  pfq2_z     = "PFQ-2 (z)"
)

# Significance labels from model summary
fe_summary <- as.data.frame(coef(summary(model))) %>%
  filter(row.names(.) != "(Intercept)") %>%
  mutate(
    predictor = row.names(.),
    sig_label = case_when(
      `Pr(>|t|)` < .001 ~ "p < .001",
      `Pr(>|t|)` < .01  ~ "p < .01",
      `Pr(>|t|)` < .05  ~ "p < .05",
      `Pr(>|t|)` < .10  ~ "p < .10 (trend)",
      TRUE               ~ "ns"
    ),
    color = case_when(
      `Pr(>|t|)` < .05  ~ "significant",
      `Pr(>|t|)` < .10  ~ "trend",
      TRUE               ~ "ns"
    )
  )

marginal_df <- marginal_df %>%
  left_join(fe_summary %>% select(predictor, sig_label, color),
            by = "predictor") %>%
  mutate(
    predictor_label = predictor_labels[predictor],
    predictor_label = factor(predictor_label,
                             levels = predictor_labels)
  )

# =============================================================================
# 7. Plot A: Faceted marginal effects (all 6 predictors)
# =============================================================================

color_values <- c(
  "significant" = "#2ca02c",
  "trend"       = "#ff7f0e",
  "ns"          = "#aaaaaa"
)

p_all <- ggplot(marginal_df,
                aes(x = z_val, y = y_hat,
                    color = color, fill = color)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15,
              color = NA) +
  geom_line(linewidth = 1.0) +
  facet_wrap(~ predictor_label, nrow = 2, scales = "free_x") +
  geom_text(
    data = marginal_df %>% group_by(predictor_label) %>% slice(1),
    aes(x = -2.3, y = Inf, label = sig_label),
    vjust = 1.5, hjust = 0, size = 3.2, color = "black"
  ) +
  scale_color_manual(values = color_values) +
  scale_fill_manual(values  = color_values) +
  labs(
    x        = "Predictor (z-score)",
    y        = "Predicted Core Affect Response Redundancy",
    title    = "Marginal Effects of Psychopathology Ratings on Core Affect Response Redundancy",
    subtitle = "Other predictors held at mean (z = 0); shading = 95% CI",
    color    = NULL, fill = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position  = "none",
    strip.background = element_rect(fill = "gray95"),
    strip.text       = element_text(face = "bold")
  )

ggsave("FIG6.pdf", p_all,
       width = 12, height = 7, dpi = 300)
cat("Saved: marginal_effects_all.png\n")

# =============================================================================
# 8. Coefficient summary table (no FDR correction)
# =============================================================================

fe_table <- as.data.frame(coef(summary(model))) %>%
  filter(row.names(.) != "(Intercept)") %>%
  rename(
    estimate = Estimate,
    se       = `Std. Error`,
    df       = df,
    t        = `t value`,
    p        = `Pr(>|t|)`
  ) %>%
  mutate(
    predictor = predictor_labels[row.names(.)],
    sig       = case_when(
      p < .001 ~ "***",
      p < .01  ~ "**",
      p < .05  ~ "*",
      p < .10  ~ ".",
      TRUE     ~ "ns"
    )
  ) %>%
  select(predictor, estimate, se, df, t, p, sig) %>%
  mutate(across(where(is.numeric), ~ round(., 3)))

cat("\n[Fixed Effects Summary (no FDR correction)]\n")
print(fe_table)

# ========
# Individual associations
# ========
predictors_z <- c("k_maia_z", "phq9_z", "gad7_z", "stai_x1_z", "erq_z", "pfq2_z")
for (pred in predictors_z) {
  formula_reduced <- as.formula(paste(
    "redundancy_z ~",
    pred,
    "+ (1 | participant) + (1 | emotion)"
  ))
  
  model_reduced <- lmer(formula_reduced, data = df, REML = FALSE)
  
  cat("\n[Random Effects Variance Components]\n")
  print(as.data.frame(VarCorr(model_reduced)))
  
  fe       <- fixef(model_reduced)
  vcov_mat <- as.matrix(vcov(model_reduced))
  
  predict_marginal_reduced <- function(predictor_name, n_points = 100) {
    
    # Range: -2.5 SD to +2.5 SD (covers ~99% of z-scored data)
    z_seq <- seq(-2.5, 2.5, length.out = n_points)
    
    all_predictors_z <- pred
    
    preds <- lapply(z_seq, function(z_val) {
      
      # Design vector: intercept + all predictors at 0 except focal
      x <- setNames(rep(0, length(all_predictors_z)), all_predictors_z)
      x[predictor_name] <- z_val
      x_vec <- c("(Intercept)" = 1, x)
      
      # Predicted value
      y_hat <- sum(fe * x_vec)
      
      # SE from delta method: sqrt(x' * vcov * x)
      se <- sqrt(as.numeric(t(x_vec) %*% vcov_mat %*% x_vec))
      
      data.frame(
        z_val  = z_val,
        y_hat  = y_hat,
        lower  = y_hat - 1.96 * se,
        upper  = y_hat + 1.96 * se
      )
    })
    
    bind_rows(preds) %>% mutate(predictor = predictor_name)
  }
  
  marginal_df_reduced <- bind_rows(lapply(
    pred,
    predict_marginal_reduced
  ))
  
  predictor_labels_reduced <- pred
  
  # Significance labels from model summary
  fe_summary_reduced <- as.data.frame(coef(summary(model_reduced))) %>%
    filter(row.names(.) != "(Intercept)") %>%
    mutate(
      predictor = row.names(.),
      sig_label = case_when(
        `Pr(>|t|)` < .001 ~ "p < .001",
        `Pr(>|t|)` < .01  ~ "p < .01",
        `Pr(>|t|)` < .05  ~ "p < .05",
        `Pr(>|t|)` < .10  ~ "p < .10 (trend)",
        TRUE               ~ "ns"
      ),
      color = case_when(
        `Pr(>|t|)` < .05  ~ "significant",
        `Pr(>|t|)` < .10  ~ "trend",
        TRUE               ~ "ns"
      )
    )
  
  marginal_df_reduced <- marginal_df_reduced %>%
    left_join(fe_summary_reduced %>% select(predictor, sig_label, color),
              by = "predictor") %>%
    mutate(
      predictor_label = predictor_labels_reduced[predictor],
      predictor_label = factor(predictor_label,
                               levels = predictor_labels_reduced)
    )
  
  color_values <- c(
    "significant" = "#2ca02c",
    "trend"       = "#ff7f0e",
    "ns"          = "#aaaaaa"
  )
  
  p_all_reduced <- ggplot(marginal_df_reduced,
                  aes(x = z_val, y = y_hat,
                      color = color, fill = color)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15,
                color = NA) +
    geom_line(linewidth = 1.0) +
    facet_wrap(~ predictor_label, nrow = 2, scales = "free_x") +
    geom_text(
      data = marginal_df_reduced %>% group_by(predictor_label) %>% slice(1),
      aes(x = -2.3, y = Inf, label = sig_label),
      vjust = 1.5, hjust = 0, size = 3.2, color = "black"
    ) +
    scale_color_manual(values = color_values) +
    scale_fill_manual(values  = color_values) +
    labs(
      x        = "Predictor (z-score)",
      y        = "Response Magnitude",
      #title    = "Marginal Effects of Psychopathology Ratings on Core Affect Response Redundancy",
      #subtitle = "Other predictors held at mean (z = 0); shading = 95% CI",
      color    = NULL, fill = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position  = "none",
      strip.background = element_rect(fill = "gray95"),
      strip.text       = element_text(face = "bold")
    )
  name_reduced <- paste("marginal_effects_magn", pred, ".png")
  ggsave(name_reduced, p_all_reduced,
         width = 4, height = 4, dpi = 150)
  cat("Saved: marginal_effects_all.png\n")
  
  fe_table <- as.data.frame(coef(summary(model_reduced))) %>%
    filter(row.names(.) != "(Intercept)") %>%
    rename(
      estimate = Estimate,
      se       = `Std. Error`,
      df       = df,
      t        = `t value`,
      p        = `Pr(>|t|)`
    ) %>%
    mutate(
      predictor = predictor_labels[row.names(.)],
      sig       = case_when(
        p < .001 ~ "***",
        p < .01  ~ "**",
        p < .05  ~ "*",
        p < .10  ~ ".",
        TRUE     ~ "ns"
      )
    ) %>%
    select(predictor, estimate, se, df, t, p, sig) %>%
    mutate(across(where(is.numeric), ~ round(., 3)))
  
  cat("\n[Fixed Effects Summary (no FDR correction)]\n")
  print(fe_table)
  name_table <- paste("lmm_fe_table_magn_", pred, "csv")
  #write_csv(fe_table, name_table)
  cat("Saved: lmm_fe_table.csv\n")
}

predictors_z <- c("k_maia_z", "phq9_z", "gad7_z", "stai_x1_z", "erq_z", "pfq2_z")
for (pred in predictors_z) {
  formula_reduced <- as.formula(paste(
    "redundancy_z ~",
    pred,
    "+ (1 | participant) + (1 | emotion)"
  ))
  cat("Dropped : ", pred)
  model_reduced <- lmer(formula_reduced, data = df, REML = FALSE)
  print(anova(model_reduced, model))
}