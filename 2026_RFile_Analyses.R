# =============================================================================
# ANALYSIS SCRIPT: It Could Be Worse: How Misaligned Social Support
#                  Relates to Test Anxiety
# Authors: Blinded for Review
# Date:    2026
# =============================================================================
#


# -----------------------------------------------------------------------------
# 0. LIBRARIES
# -----------------------------------------------------------------------------

library(tidyverse)
library(ggplot2)
library(jtools)
library(interactions)
library(sandwich)
library(dplyr)
library(parameters)
library(report)
library(broom)
library(modelsummary)
library(Hmisc)
library(psych)
library(lm.beta)
library(purrr)
library(readr)
library(flextable)
library(officer)

options(scipen = 999)


# -----------------------------------------------------------------------------
# 1. LOAD DATA
# -----------------------------------------------------------------------------

setwd("enter your pathway for data here anonymised for peer review")

data_raw <- read.csv("2025_FullDataSet.csv")


# -----------------------------------------------------------------------------
# 2. PREPARE ANALYSIS SAMPLE
# -----------------------------------------------------------------------------

# Recode Gender as factor (1 = Male, 2 = Female).
# Participants with Gender codes other than 1 or 2 (n = 4) are excluded
# from the analysis sample, consistent with the reported sample (N = 244).

data <- data_raw |>
  filter(Gender %in% c(1, 2)) |>
  mutate(Gender = factor(Gender, levels = c(1, 2), labels = c("Male", "Female")))

cat("Analysis sample N =", nrow(data), "\n")
cat("Female:", sum(data$Gender == "Female"),
    sprintf("(%.1f%%)\n", mean(data$Gender == "Female") * 100))
cat("Male:  ", sum(data$Gender == "Male"),
    sprintf("(%.1f%%)\n", mean(data$Gender == "Male")   * 100))


# -----------------------------------------------------------------------------
# 3. SCORE COMPOSITES
# -----------------------------------------------------------------------------

# ROES subscales (1-6 scale)
data <- data |>
  mutate(
    roess_reappraisal = rowMeans(pick(ROES_1, ROES_2, ROES_3, ROES_4),   na.rm = TRUE),
    roess_downward    = rowMeans(pick(ROES_5, ROES_6, ROES_7, ROES_8),   na.rm = TRUE),
    roess_humour      = rowMeans(pick(ROES_9, ROES_10, ROES_11, ROES_12), na.rm = TRUE)
  )

# SAM subscales (1-5 scale; all 7 computed, 5 used in primary analyses)
data <- data |>
  mutate(
    sams_threat         = rowMeans(pick(SAM_5,  SAM_11, SAM_20, SAM_28), na.rm = TRUE),
    sams_challenge      = rowMeans(pick(SAM_7,  SAM_8,  SAM_10, SAM_19), na.rm = TRUE),
    sams_centrality     = rowMeans(pick(SAM_6,  SAM_9,  SAM_13, SAM_27), na.rm = TRUE),  # not in primary analyses
    sams_control_self   = rowMeans(pick(SAM_12, SAM_14, SAM_22, SAM_25), na.rm = TRUE),
    sams_control_others = rowMeans(pick(SAM_4,  SAM_15, SAM_17, SAM_23), na.rm = TRUE),
    sams_uncontrollable = rowMeans(pick(SAM_1,  SAM_3,  SAM_18, SAM_21), na.rm = TRUE),
    sams_stressfulness  = rowMeans(pick(SAM_2,  SAM_16, SAM_24, SAM_26), na.rm = TRUE)   # not in primary analyses
  )

# MTAS subscales (1-5 scale)
data <- data |>
  mutate(
    mtas_worry         = rowMeans(pick(MTA_1, MTA_5, MTA_9,  MTA_13), na.rm = TRUE),
    mtas_cog_interf    = rowMeans(pick(MTA_2, MTA_6, MTA_10, MTA_14), na.rm = TRUE),
    mtas_tension       = rowMeans(pick(MTA_3, MTA_7, MTA_11, MTA_15), na.rm = TRUE),
    mtas_physiological = rowMeans(pick(MTA_4, MTA_8, MTA_12, MTA_16), na.rm = TRUE)
  )


# -----------------------------------------------------------------------------
# 4. DESCRIPTIVE STATISTICS
# -----------------------------------------------------------------------------

# --- 4a. Means and SDs ---

desc_vars <- c(
  "roess_reappraisal", "roess_downward", "roess_humour",
  "sams_threat", "sams_challenge", "sams_control_self",
  "sams_control_others", "sams_uncontrollable",
  "mtas_worry", "mtas_cog_interf", "mtas_tension", "mtas_physiological",
  "Age"
)

summary_stats <- data |>
  select(all_of(desc_vars)) |>
  summarise(across(everything(), list(
    mean = ~mean(.x, na.rm = TRUE),
    sd   = ~sd(.x,   na.rm = TRUE)
  ), .names = "{.col}__{.fn}")) |>
  pivot_longer(everything(),
               names_to  = c("Variable", "Statistic"),
               names_sep = "__",
               values_to = "Value") |>
  pivot_wider(names_from = Statistic, values_from = Value) |>
  mutate(across(c(mean, sd), ~round(.x, 2)))

print(summary_stats)


# --- 4b. Correlation matrix ---

cor_vars <- data |>
  select(
    roess_reappraisal, roess_downward, roess_humour,
    sams_threat, sams_challenge, sams_control_self,
    sams_control_others, sams_uncontrollable,
    mtas_worry, mtas_cog_interf, mtas_tension, mtas_physiological,
    Age
  ) |>
  mutate(Age = as.numeric(Age)) |>
  as.matrix()

Corrs1 <- rcorr(cor_vars, type = "pearson")

annotate_correlation_matrix <- function(corrs_object) {
  cor_matrix  <- corrs_object$r
  pval_matrix <- corrs_object$P
  annotated   <- matrix("", nrow = nrow(cor_matrix), ncol = ncol(cor_matrix))
  
  for (i in seq_len(nrow(cor_matrix))) {
    for (j in seq_len(ncol(cor_matrix))) {
      coeff <- cor_matrix[i, j]
      pval  <- pval_matrix[i, j]
      
      if (is.na(pval) | is.na(coeff)) {
        annotated[i, j] <- "—"
      } else {
        stars <- if      (pval < .001) "***"
        else if (pval < .01)  "**"
        else if (pval < .05)  "*"
        else                  ""
        annotated[i, j] <- paste0(round(coeff, 2), stars)
      }
    }
  }
  
  rownames(annotated) <- rownames(cor_matrix)
  colnames(annotated) <- colnames(cor_matrix)
  annotated
}

annotated_matrix <- annotate_correlation_matrix(Corrs1)
print(annotated_matrix)


# --- 4c. Cronbach's alpha ---

results_alpha <- data |>
  dplyr::select(
    # ROES
    Reapp_1  = ROES_1,  Reapp_2  = ROES_2,  Reapp_3  = ROES_3,  Reapp_4  = ROES_4,
    Down_1   = ROES_5,  Down_2   = ROES_6,  Down_3   = ROES_7,  Down_4   = ROES_8,
    Humour_1 = ROES_9,  Humour_2 = ROES_10, Humour_3 = ROES_11, Humour_4 = ROES_12,
    # SAM (five primary subscales)
    Threat_1       = SAM_5,  Threat_2       = SAM_11, Threat_3       = SAM_20, Threat_4       = SAM_28,
    Challenge_1    = SAM_7,  Challenge_2    = SAM_8,  Challenge_3    = SAM_10, Challenge_4    = SAM_19,
    ControlSelf_1  = SAM_12, ControlSelf_2  = SAM_14, ControlSelf_3  = SAM_22, ControlSelf_4  = SAM_25,
    ControlOther_1 = SAM_4,  ControlOther_2 = SAM_15, ControlOther_3 = SAM_17, ControlOther_4 = SAM_23,
    Uncontrol_1    = SAM_1,  Uncontrol_2    = SAM_3,  Uncontrol_3    = SAM_18, Uncontrol_4    = SAM_21,
    # MTAS
    Worry_1  = MTA_1, Worry_2  = MTA_5,  Worry_3  = MTA_9,  Worry_4  = MTA_13,
    CogInt_1 = MTA_2, CogInt_2 = MTA_6,  CogInt_3 = MTA_10, CogInt_4 = MTA_14,
    Tension_1= MTA_3, Tension_2= MTA_7,  Tension_3= MTA_11, Tension_4= MTA_15,
    Phys_1   = MTA_4, Phys_2   = MTA_8,  Phys_3   = MTA_12, Phys_4   = MTA_16
  ) |>
  pivot_longer(cols = everything(),
               names_to  = c("Factor", "Item"),
               names_pattern = "(.*)_(.*)",
               values_to = "Score") |>
  pivot_wider(names_from  = Item,
              names_prefix = "Item_",
              values_from = Score,
              values_fn   = list) |>
  mutate(Factor = as.factor(Factor))

alpha_results <- results_alpha |>
  group_by(Factor) |>
  group_split() |>
  map_df(~{
    df_items <- select(.x, starts_with("Item_")) |> map_df(unlist)
    result   <- tryCatch(
      psych::alpha(df_items, check.keys = TRUE)$total,
      error = function(e) tibble(raw_alpha = NA, std.alpha = NA)
    )
    result        <- as_tibble(result)
    result$Factor <- unique(.x$Factor)
    result
  }) |>
  relocate(Factor) |>
  arrange(Factor)

print(alpha_results)


# -----------------------------------------------------------------------------
# 5. MEAN-CENTRING FOR MODERATION ANALYSES
# -----------------------------------------------------------------------------

vars_to_center <- c(
  "sams_threat", "sams_challenge", "sams_control_self",
  "sams_control_others", "sams_uncontrollable",
  "mtas_worry", "mtas_cog_interf", "mtas_tension", "mtas_physiological",
  "roess_reappraisal", "roess_downward", "roess_humour"
)

data <- data |>
  mutate(across(
    all_of(vars_to_center),
    ~ .x - mean(.x, na.rm = TRUE),
    .names = "{.col}_c"
  ))


# -----------------------------------------------------------------------------
# 6. MODEL 1: TEST ANXIETY ~ IEER STRATEGY + GENDER (H1)
# -----------------------------------------------------------------------------

modelstest <- list(
  reapp_worry  = lm(mtas_worry        ~ roess_reappraisal + Gender, data = data),
  reapp_cog    = lm(mtas_cog_interf   ~ roess_reappraisal + Gender, data = data),
  reapp_tens   = lm(mtas_tension      ~ roess_reappraisal + Gender, data = data),
  reapp_phys   = lm(mtas_physiological ~ roess_reappraisal + Gender, data = data),
  
  down_worry   = lm(mtas_worry        ~ roess_downward + Gender, data = data),
  down_cog     = lm(mtas_cog_interf   ~ roess_downward + Gender, data = data),
  down_tens    = lm(mtas_tension      ~ roess_downward + Gender, data = data),
  down_phys    = lm(mtas_physiological ~ roess_downward + Gender, data = data),  # Gender reinstated
  
  humour_worry = lm(mtas_worry        ~ roess_humour + Gender, data = data),
  humour_cog   = lm(mtas_cog_interf   ~ roess_humour + Gender, data = data),
  humour_tens  = lm(mtas_tension      ~ roess_humour + Gender, data = data),
  humour_phys  = lm(mtas_physiological ~ roess_humour + Gender, data = data)
)

# Coefficient-level output with standardised betas
coeffs1 <- imap_dfr(modelstest, function(mod, name) {
  mod_beta <- lm.beta(mod)
  tidy(mod_beta, conf.int = TRUE) |>
    mutate(standardised_beta = lm.beta(mod)$standardized.coefficients,
           model = name)
})

# Model-level fit indices (fixed: was referencing undefined `models1`)
modelsum1 <- imap_dfr(modelstest, function(mod, name) {
  glance(mod) |> mutate(model = name)
})

write_csv(coeffs1,   "model1_coefficients.csv")
write_csv(modelsum1, "model1_modelsummary.csv")


# -----------------------------------------------------------------------------
# 7. MODEL 2: IEER STRATEGY ~ APPRAISAL + GENDER (H2)
#    Five primary appraisal subscales only (centrality and stressfulness
#    are computed but not entered into primary analyses)
# -----------------------------------------------------------------------------

models2 <- list(
  reapp_threat          = lm(roess_reappraisal ~ sams_threat         + Gender, data = data),
  reapp_challenge       = lm(roess_reappraisal ~ sams_challenge      + Gender, data = data),
  reapp_control_self    = lm(roess_reappraisal ~ sams_control_self   + Gender, data = data),
  reapp_control_others  = lm(roess_reappraisal ~ sams_control_others + Gender, data = data),
  reapp_uncontrollable  = lm(roess_reappraisal ~ sams_uncontrollable + Gender, data = data),
  
  down_threat           = lm(roess_downward ~ sams_threat         + Gender, data = data),
  down_challenge        = lm(roess_downward ~ sams_challenge      + Gender, data = data),
  down_control_self     = lm(roess_downward ~ sams_control_self   + Gender, data = data),
  down_control_others   = lm(roess_downward ~ sams_control_others + Gender, data = data),
  down_uncontrollable   = lm(roess_downward ~ sams_uncontrollable + Gender, data = data),
  
  humour_threat         = lm(roess_humour ~ sams_threat         + Gender, data = data),
  humour_challenge      = lm(roess_humour ~ sams_challenge      + Gender, data = data),
  humour_control_self   = lm(roess_humour ~ sams_control_self   + Gender, data = data),
  humour_control_others = lm(roess_humour ~ sams_control_others + Gender, data = data),
  humour_uncontrollable = lm(roess_humour ~ sams_uncontrollable + Gender, data = data)
)

coeffs2 <- imap_dfr(models2, function(mod, name) {
  mod_beta <- lm.beta(mod)
  tidy(mod_beta, conf.int = TRUE) |>
    mutate(standardised_beta = lm.beta(mod)$standardized.coefficients,
           model = name)
})

modelsum2 <- imap_dfr(models2, function(mod, name) {
  glance(mod) |> mutate(model = name)
})

print(coeffs2)

write_csv(coeffs2,   "model2_coefficients.csv")
write_csv(modelsum2, "model2_modelsummary.csv")


# -----------------------------------------------------------------------------
# 8. MODEL 3: TEST ANXIETY ~ APPRAISAL + GENDER (H3)
#    Five primary appraisal subscales only
# -----------------------------------------------------------------------------

models3 <- list(
  # Worry
  worry_threat          = lm(mtas_worry ~ sams_threat         + Gender, data = data),
  worry_challenge       = lm(mtas_worry ~ sams_challenge      + Gender, data = data),
  worry_control_self    = lm(mtas_worry ~ sams_control_self   + Gender, data = data),
  worry_control_others  = lm(mtas_worry ~ sams_control_others + Gender, data = data),
  worry_uncontrollable  = lm(mtas_worry ~ sams_uncontrollable + Gender, data = data),
  
  # Cognitive Interference
  cog_threat            = lm(mtas_cog_interf ~ sams_threat         + Gender, data = data),
  cog_challenge         = lm(mtas_cog_interf ~ sams_challenge      + Gender, data = data),
  cog_control_self      = lm(mtas_cog_interf ~ sams_control_self   + Gender, data = data),
  cog_control_others    = lm(mtas_cog_interf ~ sams_control_others + Gender, data = data),
  cog_uncontrollable    = lm(mtas_cog_interf ~ sams_uncontrollable + Gender, data = data),
  
  # Tension
  tens_threat           = lm(mtas_tension ~ sams_threat         + Gender, data = data),
  tens_challenge        = lm(mtas_tension ~ sams_challenge      + Gender, data = data),
  tens_control_self     = lm(mtas_tension ~ sams_control_self   + Gender, data = data),
  tens_control_others   = lm(mtas_tension ~ sams_control_others + Gender, data = data),
  tens_uncontrollable   = lm(mtas_tension ~ sams_uncontrollable + Gender, data = data),
  
  # Physiological Arousal
  phys_threat           = lm(mtas_physiological ~ sams_threat         + Gender, data = data),
  phys_challenge        = lm(mtas_physiological ~ sams_challenge      + Gender, data = data),
  phys_control_self     = lm(mtas_physiological ~ sams_control_self   + Gender, data = data),
  phys_control_others   = lm(mtas_physiological ~ sams_control_others + Gender, data = data),
  phys_uncontrollable   = lm(mtas_physiological ~ sams_uncontrollable + Gender, data = data)
)

coeffs3 <- imap_dfr(models3, function(mod, name) {
  mod_beta <- lm.beta(mod)
  tidy(mod_beta, conf.int = TRUE) |>
    mutate(standardised_beta = lm.beta(mod)$standardized.coefficients,
           model = name)
})

modelsum3 <- imap_dfr(models3, function(mod, name) {
  glance(mod) |> mutate(model = name)
})

print(coeffs3)

write_csv(coeffs3,   "model3_coefficients.csv")
write_csv(modelsum3, "model3_modelsummary.csv")


# -----------------------------------------------------------------------------
# 9. MODERATION ANALYSES (H4)
#    Each loop: MTA outcome ~ ROES predictor * SAM moderator + Gender
# -----------------------------------------------------------------------------

roes_vars    <- c("roess_reappraisal_c", "roess_downward_c", "roess_humour_c")
mta_outcomes <- c("mtas_worry_c", "mtas_cog_interf_c", "mtas_tension_c", "mtas_physiological_c")

run_moderation_block <- function(sam_mod_var, label) {
  models_block <- list()
  for (mta in mta_outcomes) {
    for (roes in roes_vars) {
      mod_name <- paste(mta, roes, sam_mod_var, sep = "_mod_")
      formula  <- as.formula(paste(
        mta, "~", roes, "+", sam_mod_var, "+",
        paste0(roes, ":", sam_mod_var), "+ Gender"
      ))
      models_block[[mod_name]] <- lm(formula, data = data)
    }
  }
  
  purrr::iwalk(models_block, ~{
    cat("\n\n### Model:", .y, "\n")
    print(summary(.x))
  })
  
  results   <- purrr::imap_dfr(models_block, ~ tidy(.x)   |> mutate(model = .y))
  model_fit <- purrr::imap_dfr(models_block, ~ glance(.x) |> mutate(model = .y))
  
  write_csv(results,    paste0("moderation_results_", label, ".csv"))
  write_csv(model_fit,  paste0("model_fit_",          label, ".csv"))
  
  invisible(models_block)
}

models_threat        <- run_moderation_block("sams_threat_c",         "threat")
models_challenge     <- run_moderation_block("sams_challenge_c",      "challenge")
models_control_self  <- run_moderation_block("sams_control_self_c",   "control_self")
models_control_others<- run_moderation_block("sams_control_others_c", "control_others")
models_uncontrollable<- run_moderation_block("sams_uncontrollable_c", "uncontrollable")


# -----------------------------------------------------------------------------
# 10. SIMPLE SLOPES ANALYSES
#     Probing significant interactions identified in H4
# -----------------------------------------------------------------------------

# --- Downward Comparison x Challenge ---

# Cognitive Interference
mod1 <- lm(mtas_cog_interf_c ~ roess_downward_c * sams_challenge_c + Gender, data = data)
cat("\n\n### Simple Slopes: Downward x Challenge -> Cognitive Interference\n")
sim_slopes(mod1, pred = "roess_downward_c", modx = "sams_challenge_c", modx.values = "plus-minus")
interact_plot(mod1, pred = "roess_downward_c", modx = "sams_challenge_c",
              modx.values = "plus-minus", plot.points = FALSE, interval = TRUE) +
  theme_minimal() +
  labs(x = "Downward Social Comparison (centred)",
       y = "Cognitive Interference (centred)",
       colour = "Challenge Appraisal")
ggsave("simple_slope_downward_challenge_coginterf.png", dpi = 300, width = 6, height = 4)

# Physiological Arousal
mod2 <- lm(mtas_physiological_c ~ roess_downward_c * sams_challenge_c + Gender, data = data)
cat("\n\n### Simple Slopes: Downward x Challenge -> Physiological Arousal\n")
sim_slopes(mod2, pred = "roess_downward_c", modx = "sams_challenge_c", modx.values = "plus-minus")
interact_plot(mod2, pred = "roess_downward_c", modx = "sams_challenge_c",
              modx.values = "plus-minus", plot.points = FALSE, interval = TRUE) +
  theme_minimal() +
  labs(x = "Downward Social Comparison (centred)",
       y = "Physiological Arousal (centred)",
       colour = "Challenge Appraisal")
ggsave("simple_slope_downward_challenge_physio.png", dpi = 300, width = 6, height = 4)

# Worry (non-significant; retained for completeness)
mod3 <- lm(mtas_worry_c ~ roess_downward_c * sams_challenge_c + Gender, data = data)
cat("\n\n### Simple Slopes: Downward x Challenge -> Worry\n")
sim_slopes(mod3, pred = "roess_downward_c", modx = "sams_challenge_c", modx.values = "plus-minus")


# --- Reappraisal x Control-Others ---

# Cognitive Interference
mod4 <- lm(mtas_cog_interf_c ~ roess_reappraisal_c * sams_control_others_c + Gender, data = data)
cat("\n\n### Simple Slopes: Reappraisal x Control-Others -> Cognitive Interference\n")
sim_slopes(mod4, pred = "roess_reappraisal_c", modx = "sams_control_others_c", modx.values = "plus-minus")
interact_plot(mod4, pred = "roess_reappraisal_c", modx = "sams_control_others_c",
              modx.values = "plus-minus", plot.points = FALSE, interval = TRUE) +
  theme_minimal() +
  labs(x = "Reappraisal (centred)",
       y = "Cognitive Interference (centred)",
       colour = "Control-Others Appraisal")
ggsave("simple_slope_reappraisal_controlothers_coginterf.png", dpi = 300, width = 6, height = 4)

# Physiological Arousal
mod5 <- lm(mtas_physiological_c ~ roess_reappraisal_c * sams_control_others_c + Gender, data = data)
cat("\n\n### Simple Slopes: Reappraisal x Control-Others -> Physiological Arousal\n")
sim_slopes(mod5, pred = "roess_reappraisal_c", modx = "sams_control_others_c", modx.values = "plus-minus")
interact_plot(mod5, pred = "roess_reappraisal_c", modx = "sams_control_others_c",
              modx.values = "plus-minus", plot.points = FALSE, interval = TRUE) +
  theme_minimal() +
  labs(x = "Reappraisal (centred)",
       y = "Physiological Arousal (centred)",
       colour = "Control-Others Appraisal")
ggsave("simple_slope_reappraisal_controlothers_physio.png", dpi = 300, width = 6, height = 4)


# -----------------------------------------------------------------------------
# 11. TABLE GENERATION
# -----------------------------------------------------------------------------

# Helper: format p value for table display (APA style, no leading zero)
fmt_p <- function(p) {
  case_when(
    p < .001 ~ "< .001",
    TRUE     ~ sub("^0", "", sprintf("%.3f", p))
  )
}

# Helper: extract predictor row (row 2 = strategy or appraisal, after intercept)
extract_predictor_row <- function(mod, strategy_label, outcome_label) {
  mb  <- lm.beta(mod)
  td  <- broom::tidy(mb)
  gl  <- broom::glance(mod)
  pr  <- td[2, ]
  tibble(
    Outcome  = outcome_label,
    Predictor = strategy_label,
    B        = round(pr$estimate, 2),
    SE       = round(pr$std.error, 2),
    beta     = round(mb$standardized.coefficients[[2]], 2),
    t        = round(pr$statistic, 2),
    p        = pr$p.value,
    R2       = round(gl$r.squared, 3)
  )
}


# ---- TABLE 1: Descriptive Statistics and Intercorrelation Matrix ----

# Gender as numeric for correlation (1 = Male, 2 = Female)
data <- data |> mutate(Gender_num = as.numeric(Gender))

table_vars <- c(
  "roess_reappraisal", "roess_downward", "roess_humour",
  "sams_threat", "sams_challenge", "sams_control_self",
  "sams_control_others", "sams_uncontrollable",
  "mtas_worry", "mtas_cog_interf", "mtas_tension", "mtas_physiological",
  "Gender_num"
)

var_labels <- c(
  "1. Reappraisal", "2. Downward Comparison", "3. Humour",
  "4. Threat", "5. Challenge", "6. Control-Self",
  "7. Control-Others", "8. Uncontrollability",
  "9. Worry", "10. Cognitive Interference", "11. Tension", "12. Physiological Arousal",
  "13. Gender"
)

# Means and SDs for continuous variables only (Gender reported as n/% in text)
continuous_vars <- table_vars[table_vars != "Gender_num"]

desc_tbl <- data |>
  select(all_of(continuous_vars)) |>
  summarise(across(everything(), list(
    M  = ~round(mean(.x, na.rm = TRUE), 2),
    SD = ~round(sd(.x,   na.rm = TRUE), 2)
  ), .names = "{.col}__{.fn}")) |>
  pivot_longer(everything(), names_to = c("var", "stat"), names_sep = "__", values_to = "val") |>
  pivot_wider(names_from = stat, values_from = val) |>
  # Add a blank row for Gender (M and SD not applicable)
  add_row(var = "Gender_num", M = NA, SD = NA)

# Pearson correlations with significance (all 13 variables including Gender)
cor_mat  <- data |> select(all_of(table_vars)) |> as.matrix()
cor_res  <- rcorr(cor_mat, type = "pearson")
r_vals   <- cor_res$r
p_vals_c <- cor_res$P
n_vars   <- length(table_vars)

sig_stars <- function(p) {
  case_when(is.na(p) ~ "", p < .001 ~ "***", p < .01 ~ "**", p < .05 ~ "*", TRUE ~ "")
}

# Lower-triangle character matrix (n_vars rows x n_vars-1 cols)
cor_display <- matrix("", nrow = n_vars, ncol = n_vars - 1)
for (i in seq_len(n_vars)) {
  for (j in seq_len(n_vars - 1)) {
    if (j >= i) {
      cor_display[i, j] <- ""
    } else {
      r <- round(r_vals[i, j], 2)
      s <- sig_stars(p_vals_c[i, j])
      cor_display[i, j] <- paste0(sprintf("%.2f", r), s)
    }
  }
}

col_headers <- as.character(seq_len(n_vars - 1))

# Format M and SD: display "—" for Gender (categorical covariate)
M_display  <- ifelse(is.na(desc_tbl$M),  "—", sprintf("%.2f", desc_tbl$M))
SD_display <- ifelse(is.na(desc_tbl$SD), "—", sprintf("%.2f", desc_tbl$SD))

table1_df <- tibble(
  Variable = var_labels,
  M        = M_display,
  SD       = SD_display
) |>
  bind_cols(
    as_tibble(cor_display, .name_repair = "minimal") |> setNames(col_headers)
  )

ft1 <- flextable(table1_df) |>
  set_header_labels(
    Variable = "Variable", M = "M", SD = "SD",
    .list = setNames(as.list(col_headers), col_headers)
  ) |>
  bold(part = "header") |>
  align(j = seq(4, ncol(table1_df)), align = "center", part = "all") |>
  fontsize(size = 10, part = "all") |>
  font(fontname = "Times New Roman", part = "all") |>
  add_footer_lines(paste0(
    "Note. N = ", nrow(data), ". Correlations are Pearson's r. ",
    "Gender coded 1 = Male, 2 = Female. ",
    "* p < .05. ** p < .01. *** p < .001."
  )) |>
  autofit()

save_as_docx(ft1, path = "Table1_Descriptives_Correlations.docx")
cat("Table 1 saved.\n")

# ---- TABLE 2: IEER Strategies Predicting Situational Appraisals (H2) ----

appraisal_map <- list(
  "Threat"            = "threat",
  "Challenge"         = "challenge",
  "Control-Self"      = "control_self",
  "Control-Others"    = "control_others",
  "Uncontrollability" = "uncontrollable"
)

strategy_map <- list(
  "Reappraisal"         = "reapp",
  "Downward Comparison" = "down",
  "Humour"              = "humour"
)

# Extract all model results into a long tibble
table2_long <- imap_dfr(appraisal_map, function(appr_suffix, appr_label) {
  imap_dfr(strategy_map, function(strat_prefix, strat_label) {
    key     <- paste(strat_prefix, appr_suffix, sep = "_")
    matched <- models2[names(models2) == key]
    if (length(matched) == 0) return(NULL)
    extract_predictor_row(matched[[1]], strat_label, appr_label) |>
      mutate(p_fmt = fmt_p(as.numeric(p))) |>
      select(-p)
  })
})

# Build wide table manually to guarantee column order.
# pivot_wider sorts alphabetically by both names_from and values_from,
# which swaps t and p and reorders strategy groups — so we do this by hand.
strat_names   <- names(strategy_map)   # Reappraisal, Downward Comparison, Humour
appraisal_names <- names(appraisal_map)

table2_wide <- map_dfc(strat_names, function(s) {
  sub_df <- table2_long |>
    filter(Predictor == s) |>
    arrange(match(Outcome, appraisal_names)) |>
    select(B, SE, beta, t, p_fmt, R2)
  names(sub_df) <- paste(s, c("B", "SE", "beta", "t", "p", "R2"), sep = "||")
  sub_df
}) |>
  mutate(Appraisal = appraisal_names, .before = 1)

# For flextable, column names must be unique — use the Predictor||stat names as-is,
# and map them to display labels via set_header_labels
col_keys     <- names(table2_wide)
stat_display <- c("B", "SE", "\u03B2", "t", "p", "R\u00B2")
display_map  <- setNames(
  as.list(c("Appraisal", rep(stat_display, times = length(strat_names)))),
  col_keys
)

ft2 <- flextable(table2_wide, col_keys = col_keys) |>
  set_header_labels(values = display_map) |>
  add_header_row(
    values    = c("", strat_names),
    colwidths = c(1, rep(6, length(strat_names)))
  ) |>
  merge_h(part = "header") |>
  bold(part = "header") |>
  align(part = "header", align = "center") |>
  align(j = -1, align = "center", part = "body") |>
  fontsize(size = 10, part = "all") |>
  font(fontname = "Times New Roman", part = "all") |>
  add_footer_lines(paste0(
    "Note. N = ", nrow(data), ". B = unstandardised coefficient. ",
    "\u03B2 = standardised coefficient. ",
    "Gender included as covariate in all models (results not shown). ",
    "* p < .05. ** p < .01. *** p < .001."
  )) |>
  autofit()

save_as_docx(ft2, path = "Table2_H2_StrategiesOnAppraisals.docx")
cat("Table 2 saved.\n")


# ---- TABLE 3: Appraisals Predicting Test Anxiety (H3) ----

anxiety_map <- list(
  "Worry"                  = "worry",
  "Cognitive Interference" = "cog",
  "Tension"                = "tens",
  "Physiological Arousal"  = "phys"
)

table3_rows <- imap_dfr(anxiety_map, function(anx_prefix, anx_label) {
  imap_dfr(appraisal_map, function(appr_suffix, appr_label) {
    key     <- paste(anx_prefix, appr_suffix, sep = "_")
    matched <- models3[names(models3) == key]
    if (length(matched) == 0) return(NULL)
    extract_predictor_row(matched[[1]], appr_label, anx_label)
  })
}) |>
  mutate(p_fmt = fmt_p(p)) |>
  select(-p) |>
  rename(
    `Test Anxiety Component` = Outcome,
    Appraisal                = Predictor,
    p                        = p_fmt,
    beta_std                 = beta
  )

# Rename beta column to the Greek symbol after the tibble is built
# (unicode in backtick column names is not supported in R; set via flextable header instead)

ft3 <- flextable(table3_rows) |>
  set_header_labels(
    `Test Anxiety Component` = "Test Anxiety Component",
    Appraisal = "Appraisal",
    B         = "B",
    SE        = "SE",
    beta_std  = "\u03B2",
    t         = "t",
    p         = "p",
    R2        = "R\u00B2"
  ) |>
  merge_v(j = "Test Anxiety Component") |>
  bold(part = "header") |>
  align(j = -1, align = "center", part = "body") |>
  fontsize(size = 10, part = "all") |>
  font(fontname = "Times New Roman", part = "all") |>
  add_footer_lines(paste0(
    "Note. N = ", nrow(data), ". B = unstandardised coefficient. ",
    "\u03B2 = standardised coefficient. ",
    "Gender included as covariate in all models (results not shown). ",
    "* p < .05. ** p < .01. *** p < .001."
  )) |>
  autofit()

save_as_docx(ft3, path = "Table3_H3_AppraisalsOnAnxiety.docx")
cat("Table 3 saved.\n")

cat("\nAll analyses complete.\n")