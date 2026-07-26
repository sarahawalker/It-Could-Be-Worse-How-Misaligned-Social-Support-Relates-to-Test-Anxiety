library(tidyverse)
library(interactions)

# ---- Load and prepare data ----
setwd("anonymised for review")
df <- read_csv("2025_FullDataSet.csv") %>%
  filter(Gender %in% c(1, 2)) %>%
  mutate(Gender = factor(Gender, levels = c(1, 2), labels = c("Male", "Female")))

# ---- Build composite scores ----
df <- df %>%
  mutate(
    roess_downward      = rowMeans(select(., ROES_5, ROES_6, ROES_7, ROES_8)),
    sams_challenge       = rowMeans(select(., SAM_7, SAM_8, SAM_10, SAM_19)),
    mtas_worry           = rowMeans(select(., MTA_1, MTA_5, MTA_9, MTA_13)),
    mtas_cog_interf      = rowMeans(select(., MTA_2, MTA_6, MTA_10, MTA_14)),
    mtas_tension         = rowMeans(select(., MTA_3, MTA_7, MTA_11, MTA_15)),
    mtas_physiological   = rowMeans(select(., MTA_4, MTA_8, MTA_12, MTA_16))
  )

# ---- Mean-centre predictors ----
df <- df %>%
  mutate(
    roess_downward_c  = roess_downward - mean(roess_downward),
    sams_challenge_c  = sams_challenge - mean(sams_challenge)
  )

# ---- Correlation: downward comparison x challenge appraisal ----
cor.test(df$roess_downward, df$sams_challenge)

# ---- Main effect models (no moderator) ----
m_worry_main   <- lm(mtas_worry ~ roess_downward + Gender, data = df)
m_cog_main     <- lm(mtas_cog_interf ~ roess_downward + Gender, data = df)
m_tension_main <- lm(mtas_tension ~ roess_downward + Gender, data = df)
m_physio_main  <- lm(mtas_physiological ~ roess_downward + Gender, data = df)

summary(m_worry_main)
summary(m_cog_main)
summary(m_tension_main)
summary(m_physio_main)

# ---- Moderation models (downward comparison x challenge appraisal) ----
m_worry   <- lm(mtas_worry ~ roess_downward_c * sams_challenge_c + Gender, data = df)
m_cog     <- lm(mtas_cog_interf ~ roess_downward_c * sams_challenge_c + Gender, data = df)
m_tension <- lm(mtas_tension ~ roess_downward_c * sams_challenge_c + Gender, data = df)
m_physio  <- lm(mtas_physiological ~ roess_downward_c * sams_challenge_c + Gender, data = df)

summary(m_worry)
summary(m_cog)
summary(m_tension)
summary(m_physio)

# ---- Simple slopes at +/-1 SD of challenge appraisal ----
sim_slopes(m_worry, pred = roess_downward_c, modx = sams_challenge_c,
           johnson_neyman = FALSE)
sim_slopes(m_cog, pred = roess_downward_c, modx = sams_challenge_c,
           johnson_neyman = FALSE)
sim_slopes(m_physio, pred = roess_downward_c, modx = sams_challenge_c,
           johnson_neyman = FALSE)

# ---- Bonferroni threshold for this 4-test family ----
0.05 / 4
