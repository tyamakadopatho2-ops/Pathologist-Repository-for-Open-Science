# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Tetsuhiro Yamakado
#
# DMD Research Script Commentary
# Annotated companion R script for the doctoral dissertation:
# "A Systematic Histopathological Study of Duchenne Muscular Dystrophy
# Using Semi-Quantitative Image Analysis, Digital Restoration Techniques,
# and Exploratory Statistical Approaches"
# (Hokkaido University, 2026)
#
# Provenance:
# This file is an annotated and expanded dissertation version of the DMD
# analysis scripts developed for the related study and dissertation project.
#
# Scope:
# This script contains:
#   (1) confirmatory analyses reported in the doctoral dissertation and/or
#       associated publication(s), and
#   (2) clearly marked exploratory or preliminary analyses retained for
#       transparency, reproducibility, and research traceability.
#
# Interpretation note:
# Exploratory/preliminary sections are provided for transparency and should
# not be interpreted as the primary inferential basis of the peer-reviewed
# article unless explicitly stated in the dissertation text.
#
# Changes from the published-paper version:
# - added explanatory comments and annotations
# - reorganized sections for dissertation-level reproducibility
# - added exploratory/preliminary analytical sections
#
# Research-use notice:
# This script is shared for transparency and reproducibility of academic
# research. It is not validated for clinical decision-making.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


### NOTE FOR RSTUDIO USERS ###

# FIRST, PRESS "Alt + O" (Windows) or "Option + O" (macOS) 
# TO FOLD ALL SECTIONS AT ONCE.

# Expand and run sections sequentially as needed, starting from the top.
# Some sections depend on the execution of earlier ones.
# (e.g., data preparation sections must be executed first)
# If you want to unfold all sections,
# press "Shift + Alt + O" (Windows) or "Shift + Option + O"(macOS).

# NOTE 1:
# We acknowledge that the statistical methods used in this study
# are extensive and varied.
#
# The rationale for this comprehensive approach is threefold:
# (1) To address concerns arising from the small sample size 
#     and ensure robust, reliable findings.
# (2) To provide careful validation of results by utilizing
#     multiple statistical analyses rather than relying on a simple approach.
# (3) To accommodate the exploratory nature of our study,
#     given that the novel biomarker introduced here (MFD) currently lacks 
#     established validation frameworks or standard evaluation criteria.

# We greatly appreciate the reader's understanding regarding the necessity 
# and appropriateness of this thorough analytical strategy.


## Parameters ---
#  (1) Mean: The mean myofiber size (Mean, μm²);
#  (2) Sd: The standard deviation of myofiber size (Sd, μm²);
#  (3) Cov: The coefficient of variation of myofiber size (Cov); 
#  (4) MFD: Myofiber density (MFD, fibers/mm²);
#  (5) MFA: Myofiber area (MFA, %);
#  (6) CFA: Connective/Fibrotic tissue area (CFA, %);
#  (7) NFA: Necrotic fiber area (NFA, %);
#  (8) RFA: Regenerative fiber area (RFA, %);
#  (9) Fat: Fatty degeneration area (Fat, %);
# (10) Opaque: Percentage of opaque fibers relative to MFD (Opaque, count/MFD, %); and
# (11) IntN: Percentage of internally nucleated fibers relative to MFD (IntN, count/MFD, %)

################################################################################
# Preparation 1--------------------------------------------------------------
{
  # NOTE: The initial run may take some time.
  packages <- c("dplyr", "ggplot2", "ggfortify", "scales","car", "PMCMRplus", 
                "pwr", "effsize", "multcomp", "MBESS", "compute.es",
                "cvTools", "boot", "segmented","purrr", "brms",
                "stan", "posterior", "hdi","HDInterval")
  
  # Install missing packages
  install_if_missing <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
    }
  }
  
  # Check each package
  sapply(packages, install_if_missing)
}


# Essential libraries for DMD study analysis
{
  # Statistical analysis
  library(car)          # Levene's test
  library(PMCMRplus)    # Post-hoc tests: Games-Howell, Conover
  library(pwr)          # Power calculations
  library(effsize)      # Effect size
  library(multcomp)     # Multiple comparisons
  library(MBESS)        # Confidence intervals for effect sizes
  library(compute.es)   # Effect sizes for various test statistics
  library(segmented)    # Segmented regression analysis
  library(purrr)        # Functional programming helpers (map/reduce) for tidy iteration
  
  # Data manipulation and visualization
  library(dplyr)
  library(ggplot2)
  library(ggfortify)
  library(scales)
  
  # Cross-validation and modeling
  library(cvTools)      # Cross-validation tools
  
  # Bootstrap method
  library(boot)         # Bootstrapping tools
  
  # Bayesian analysis
  library(brms)
  library(rstan)
  library(posterior)
  library(hdi)
  library(HDInterval)
}

# rm(list=ls()) # Clear all variables from workspace

L <- read.csv("File_Path_of_.csv")
# or
L <- read.csv(file.choose()) # Choose source.csv 
# (covert the provided Excel file (eData 1) to CSV format before this execution)

{
  ABx <- L$ABx
  Mean <- L$Mean
  Cov <- L$Cov
  Sd <- L$Sd
  MFD <- L$MFD
  MFA <- L$MFA*100
  CFA <- L$CFA*100
  Fat <- L$Fat*100
  IntN <- L$IntN*100
  Opaque <- L$Opaque*100
  NFA <- L$NFA*100
  RFA <- L$RFA*100
}

{
  dataAll <- data.frame(
    ABx = L$ABx,
    Mean = L$Mean,
    Cov = L$Cov,
    Sd = L$Sd,
    MFD = L$MFD,
    MFA = L$MFA*100, 
    CFA = L$CFA*100,
    Fat = L$Fat*100,
    IntN = L$IntN*100,
    Opaque = L$Opaque*100, 
    NFA = L$NFA*100, 
    RFA = L$RFA*100
  )
}

# Preparation 2--------------------------------------------------------------
{
  # log transformation
  {
    LABx <- log(L$ABx); ABx <- L$ABx
    LMean <- log(L$Mean); Mean <- L$Mean
    LCov <- log(L$Cov); Cov <- L$Cov
    LSd <- log(L$Sd); Sd <- L$Sd
    LMFD <- log(L$MFD); MFD <- L$MFD
    LIntN <- log(L$IntN); IntN <- L$IntN
    LMFA <- log(L$MFA); MFA <- L$MFA
    LCFA <- log(L$CFA); CFA <- L$CFA
    LNFA <- log(L$NFA); NFA <- L$NFA
    LRFA <- log(L$RFA); RFA <- L$RFA
    LFat <- log(L$Fat); Fat <- L$Fat
    LOpaque <- log(L$Opaque); Opaque <- L$Opaque
  }
  
  # Data Frame
  {
    L0 <- data.frame(LABx, LCov, LSd, LMFD, LIntN, LMFA, LCFA, LMean,
                     LOpaque, LNFA, LRFA, LFat)
    L1 <- as.data.frame(scale(L0))  
  }
  
  # z score
  {
    L1$ABx.C <- L1$LABx; ABx.C <- L1$ABx.C
    L1$Mean.C <- L1$LMean; Mean.C <- L1$Mean.C
    L1$Cov.C <- L1$LCov; Cov.C <- L1$Cov.C
    L1$CFA.C <- L1$LCFA; CFA.C <- L1$CFA.C
    L1$IntN.C <- L1$LIntN; IntN.C <- L1$IntN.C
    L1$MFD.C <- L1$LMFD; MFD.C <- L1$MFD.C
    L1$Sd.C <- L1$LSd; Sd.C <- L1$Sd.C
    L1$MFA.C <- L1$LMFA; MFA.C <- L1$MFA.C
    L1$Opaque.C <- L1$LOpaque; Opaque.C <- L1$Opaque.C
    L1$NFA.C <- L1$LNFA; NFA.C <- L1$NFA.C
    L1$RFA.C <- L1$LRFA; RFA.C <- L1$RFA.C
    L1$Fat.C <- L1$LFat; Fat.C <- L1$Fat.C
  }
  
  # centering (log-transform only)
  {
    ABx.c = LABx - mean(log(ABx))
    Mean.c = LMean - mean(log(Mean))
    IntN.c = LIntN - mean(log(IntN)) 
    CFA.c = LCFA - mean(log(CFA))
    Cov.c = LCov - mean (log(Cov))
    MFA.c = LMFA - mean(log(MFA))
    MFD.c = LMFD - mean(log(MFD))
    Sd.c = LSd - mean(log(Sd))
    Fat.c = LFat - mean(log(Fat))
    NFA.c = LNFA - mean(log(NFA))
    Opaque.c = LOpaque - mean(log(Opaque))
    RFA.c = LRFA - mean(log(RFA))
  }
}

# Correlation table-------------------------------------------------------------
{ df_cor <- L %>% 
  dplyr::select(ABx, Mean, Cov, Sd, MFD, MFA, CFA, IntN, Opaque, NFA, RFA, Fat) 

  
  df_log <- df_cor %>% 
    mutate(
      log_ABx = log(ABx),
      log_Mean = log(Mean),
      log_Sd = log(Sd),
      log_Cov = log(Cov),
      log_MFD = log(MFD),
      log_MFA = log(MFA),
      log_CFA = log(CFA),
      log_IntN = log(IntN),
      log_Opaque = log(Opaque),
      log_NFA = log(NFA),
      log_RFA = log(RFA),
      log_Fat = log(Fat)
    )
  
  df_log <- dplyr::select(df_log, log_ABx, log_Mean, log_Sd, log_Cov, log_MFD,
                          log_MFA, log_CFA, log_IntN, log_Opaque, log_NFA, log_RFA, log_Fat)
  cor(df_log, method = "pearson")
}

# examples
cor.test(log(ABx), log(MFD), method = "pearson") # -0.8597013
shapiro.test(MFD) # p-value = 0.002265
shapiro.test(log(MFD)) # p-value = 0.8188
qplot(log(ABx), log(MFD)) + 
  geom_point(size = 2, colour = "black") +
  geom_smooth(method = "lm", colour = "black", se = F)+
  scale_x_continuous(limits = c(0, 3), breaks = 0:3) +
  theme_bw()+
  theme(
    axis.title.x = element_text(size = 20, face = "bold"), 
    axis.title.y = element_text(size = 20, face = "bold"), 
    axis.text.x = element_text(size = 16), 
    axis.text.y = element_text(size = 16)  
  )

cor.test(log(ABx), log(Sd), method = "pearson") # 0.7444486
shapiro.test(Sd) # p-value = 9.15e-06
shapiro.test(log(Sd)) # p-value = 0.5848
qplot(log(ABx), log(Sd)) + 
  geom_point(size = 2, colour = "black") +
  geom_smooth(method = "lm", colour = "black", se = F)+
  scale_x_continuous(limits = c(0, 3), breaks = 0:3) +
  theme_bw()+
  theme(
    axis.title.x = element_text(size = 20, face = "bold"), 
    axis.title.y = element_text(size = 20, face = "bold"), 
    axis.text.x = element_text(size = 16), 
    axis.text.y = element_text(size = 16)  
  )

###------------------------------### 0. Batch Effect Validation-----------------------------

# PART 1: Preliminary comparative analysis of MFA values using H&E and Gomori staining--------------------------

# MFA data obtained by HE and Gomori staining
HE <- c(55.74, 62.03, 67.74, 63.06, 46.33, 58.81, 68.59, 38.23, 60.54, 73.03,
        62.090, 39.35, 70.43, 50.20, 58.25, 51.86, 30.70, 59.25, 67.42, 27.66,
        50.96, 60.48, 43.13, 42.25, 38.62, 56.94, 31.04, 70.90)
Gomori <- c(51.42, 56.25, 63.66, 59.93, 48.14, 51.39, 69.42, 36.54, 61.07,
            66.46, 62.87, 36.77, 73.02, 49.04, 61.14, 53.31, 34.83, 60.12,
            69.29, 23.83, 50.72, 62.091, 42.74, 46.90, 41.25, 58.59, 
            32.14, 72.06)

# Wilcoxon Signed-Rank Test (Non-parametric Test)
wilcox.test(HE, Gomori) # p-value = 0.9288 (No statistically significant difference)

# Nonparametric effect size：Cliff's Delta
cliff.delta(HE, Gomori) 
# effect size = 0.01530612 (negligible) 
# 95％ CI (-0.2848065 to 0.3126860)

# Normality Test (Shapiro-Wilk test)
shapiro.test(HE)    # p-value = 0.13（It is thought to follow a normal distribution）
shapiro.test(Gomori) # p-value = 0.401（It is thought to follow a normal distribution）

# Test of equal variances (F-test)
var.test(HE, Gomori) # p-value = 0.9287（It can be considered to have equal variances）

# t-test assuming equal variances
t.test(HE, Gomori, var.equal = T) # p-value = 0.9137（not significant）

# Calculation of effect size (Hedges' d) and estimation of 95% CI by bootstrap
{
  # A function to calculate the confidence interval of Hedges' d using bootstrap
  bootstrap_ci_hedges <- function(data1, data2,
                                  n_iter = 5000,        # Bootstrap iterations
                                  conf_level = 0.95,    # Confidence interval level (95%)
                                  seed = 1234) {        # Random seed for reproducibility
    # the seed is fixed to ensure reproducibility
    if (!is.null(seed)) set.seed(seed)
    
    # Sample size for each group
    n1 <- length(data1)
    n2 <- length(data2)
    
    # Vector for storing the results
    boot_d <- numeric(n_iter)
    
    # Bootstrap processing
    for (i in seq_len(n_iter)) {
      s1 <- sample(data1, n1, replace = TRUE) # Recovery and extraction from data1
      s2 <- sample(data2, n2, replace = TRUE) # Recovery and extraction from data1
      # Calculation of Hedges' d (with correction)
      boot_d[i] <- cohen.d(s1, s2, hedges.correction = TRUE)$estimate
    }
    
    # Hedges' d (without bootstrap) based on original data
    original_d <- cohen.d(data1, data2, hedges.correction = TRUE)$estimate
    
    # Calculation of confidence intervals (CI)
    alpha <- 1 - conf_level
    ci_lower <- quantile(boot_d, alpha / 2)         # lower limit
    ci_upper <- quantile(boot_d, 1 - alpha / 2)     # upper limit
    
    # Returns the results as a list
    list(
      original_d = original_d,  # Original Hedges' d value
      boot_mean  = mean(boot_d), # Mean value of Hedges' d via bootstrap
      ci_lower   = ci_lower,     # lower limit of CI
      ci_upper   = ci_upper,     # upper limit of CI
      boot_dist  = boot_d        # Bootstrap distribution (use when necessary)
    )
  }
  
  # Perform calculation of Hedges' d and estimation of CI
  res <- bootstrap_ci_hedges(HE, Gomori)
  
  # Calculate the measured detection power (Observed Power)
  power_obs <- pwr.t2n.test(
    n1        = length(HE),
    n2        = length(Gomori),
    d         = res$original_d,
    sig.level = 0.05
  )$power * 100 # % Notation
  
  # Calculate the sample size for each group required to obtain 80% detection power
  required_N <- ceiling(
    pwr.t.test(
      d         = abs(res$original_d),
      power     = 0.80,
      sig.level = 0.05,
      type      = "two.sample"
    )$n
  )
  
  # Displaying Results
  cat("HE vs Gomori  (n =", length(HE), "vs", length(Gomori), ")\n")
  cat("Hedges' d (Original)        :", sprintf("%.4f", res$original_d), "\n")
  cat("Hedges' d (Bootstrap mean)  :", sprintf("%.4f", res$boot_mean),  "\n")
  # Display of confidence interval (explicitly stated as 95%)
  cat(sprintf("%d%% Bootstrap CI            : %.4f  –  %.4f\n",
              100 * 0.95, res$ci_lower, res$ci_upper))
  cat("Observed power (%)          :", sprintf("%.2f", power_obs), "\n")
  cat("Required N per group (80%%)  :", required_N, "\n")
}



# Correlation analysis
cor.test(HE, Gomori) # 0.9699253(0.97)

# plot
dataHEGomori <- data.frame(HE, Gomori)

ggplot(dataHEGomori, aes(x = HE, y = Gomori)) +
  geom_smooth(method = "lm", colour = "red", se = F) +
  geom_point(size = 3, colour = 'black') +
  theme_bw()+
  theme(
    axis.title.x = element_text(size = 30, face = "bold"), 
    axis.title.y = element_text(size = 30, face = "bold"), 
    axis.text.x = element_text(size = 16), 
    axis.text.y = element_text(size = 16)  
  ) +
  ylab("MFA (G-T)") +
  xlab("MFA (H&E)") + 
  coord_cartesian(clip = "off")+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black")


# Comparison between MFA and CFA
cor.test(MFA, CFA, method = "pearson") # -0.9625728, data: all dataset (n = 38)
shapiro.test(MFA) # p-value = 0.06048
shapiro.test(CFA) # p-value = 0.07888

wilcox.test(MFA, CFA) # p-value = 0.3081
cliff.delta(MFA, CFA) # 0.1371191 (negligible) (-0.1238352 to 0.3803345)

qplot(MFA, CFA) +
  geom_point(size = 3, colour = 'black') +
  ylab("CFA") +
  xlab("MFA") + 
  theme_bw() +
  coord_cartesian(clip = "off") +
  geom_smooth(method = "lm", se = F, color = "red")


# PART 2: Bayes Approach--------------------------------------------------------
{
  ## 1. JZS Bayes factor (paired t) ---
  # Function definition for performing the JZS Bayes t-test (Rouder et al., 2009)
  bf_ttest <- function(t, n, r = sqrt(2)/2){
    v <- n - 1
    j <- 1 + n * r^2
    logBF10 <- lgamma((v + 1)/2) - lgamma(v/2) -
      0.5 * log(pi * v) - 0.5 * log(j) -
      (v + 1)/2 * log(1 + t^2 / (v * j))
    exp(logBF10)
  }
  
  d   <- HE - Gomori                               # Difference between HE and Gomori
  n   <- length(d)                                 # sample size
  tstat <- t.test(HE, Gomori, paired = TRUE)$statistic # Paired t-test statistic
  BF10  <- bf_ttest(tstat, n = length(HE), r = 0.707) # Bayes factor BF10 (standard r=0.707)
  BF01  <- 1 / BF10                                # Bayes factor BF01 (null hypothesis)
  
  ## 2. Posterior draws (normal approx) ---
  # Sampling of the posterior distribution (effect size δ) 
  # from an approximate standard normal distribution
  post_mean <- as.numeric(tstat) / sqrt(n)         # Posterior mean
  post_sd   <- sqrt(1/n)                           # Posterior standard deviation
  
  set.seed(123)
  delta_draws <- rnorm(10000, post_mean, post_sd)  # Posterior distribution sample of δ (10000 times)
  CrI95 <- quantile(delta_draws, c(.025, .975))    # 95% credible interval
  
  ## 3. ROPE decision ---
  rope_fix <- 0.30                                 # ROPE range setting (Lakens, 2022; recommended)
  rope_pct <- mean(abs(delta_draws) < rope_fix) * 100 # robability of falling within the ROPE range (%)
  
  # For reference: δ_error based on TEM
  sd_d      <- sd(d)
  sd_pool   <- sqrt((var(HE) + var(Gomori)) / 2)
  delta_err <- sd_d / sd_pool                      # δ_error (apporoximately 0.246)
  
  # Width setting for ROPE sensitivity analysis
  rope_grid <- c(0.25, 0.30, 0.40)
  grid_pct  <- sapply(rope_grid,
                      function(r) mean(abs(delta_draws) < r) * 100)
  
  ## 4. Batch (specimen) random effect τ² (Test using a linear mixed-effects model) ---
  library(Matrix)
  library(lme4)
  
  # Data preparation for the mixed-effects model
  df <- data.frame(
    value = c(HE, Gomori),
    stain = factor(rep(c("HE", "Gomori"), each = length(HE))),
    specimen = factor(rep(seq_along(HE), times = 2))
  )
  
  # Fit the mixed-effects model (unadjusted and stain-adjusted)
  fit0 <- lmer(value ~ 1 + (1|specimen), data = df, REML = TRUE)  # unadjusted
  fit1 <- lmer(value ~ stain + (1|specimen), data = df, REML = TRUE) # Stain-adjusted
  
  # Sample-level variance (τ²) extraction
  tau2_0 <- as.data.frame(VarCorr(fit0))$vcov[1]   # unadjusted τ²
  tau2_1 <- as.data.frame(VarCorr(fit1))$vcov[1]   # stain-adjusted τ²
  delta_tau <- 100 * (tau2_0 - tau2_1) / tau2_0    # τ² reduction percentage (%)
  
  ## 5. Output of Results ---
  cat("\n===========  Bayesian comparison: HE  vs  Gomori  ===========\n")
  cat(sprintf("t statistic (paired)      :  %.3f\n", tstat))
  cat(sprintf("Bayes Factor  BF10        :  %.3f\n", BF10))
  cat(sprintf("Bayes Factor  BF01        :  %.3f   (%s)\n",
              BF01, ifelse(BF10 < 1,
                           "supports null model (no difference)",
                           "supports alternative model (difference)")))
  cat(sprintf("Posterior mean (delta)    :  %.3f\n",  mean(delta_draws)))
  cat(sprintf("Posterior 95%% HDI         : [%.3f , %.3f]\n",
              CrI95[1], CrI95[2]))
  
  cat("\n--  ROPE decision (Lakens 2022 : |d| < 0.30) --------------\n")
  cat(sprintf("ROPE ±%.2f → Posterior mass  %.1f %%\n", rope_fix, rope_pct))
  
  # ROPE sensitivity analysis results
  cat("\n--  ROPE sensitivity analysis -----------------------------\n")
  for(i in seq_along(rope_grid)){
    cat(sprintf("ROPE ±%.3f →  %.1f %% inside\n",
                rope_grid[i], grid_pct[i]))
  }
  
  # Results of reduction at sample level τ²
  cat("\n--  Batch effect (specimen‑level τ², mixed-effects model) -----------------------\n")
  cat(sprintf("τ²_0 (unadjusted)         :  %.4f\n", tau2_0))
  cat(sprintf("τ²_1 (stain‑adjusted)     :  %.4f\n", tau2_1))
  cat(sprintf("Δτ² reduction             :  %.2f %%  →  %s\n",
              delta_tau,
              ifelse(delta_tau < 5, "stain effect negligible", "caution: non‑negligible")))
  cat("============================================================\n\n")
}


## 6. Visualizations-------------------------------

## 6-A: Density Plot of Posterior Distribution (delta) and Display of ROPE Range
ggplot(data.frame(delta = delta_draws), aes(delta)) +
  geom_density(fill = "#3182BD", alpha = .4) +  # Draw density curve (blue with transparency 0.4)
  annotate("rect", xmin = -rope_fix, xmax = rope_fix,  # Display ROPE range (±0.30) in gray
           ymin = 0, ymax = Inf,
           fill = "grey70", alpha = .15) +
  geom_vline(xintercept = 0, linetype = "dashed") +    # Display center line (delta=0) as a dashed line
  labs(
    x = "Standardised effect size (delta)",            # X-axis label: Standardized effect size (delta)
    y = "Posterior density",                           # Y-axis label: Posterior density
    title = "Posterior of delta (HE – Gomori)",        # Title: Posterior distribution of 
    subtitle = sprintf(                                # delta (comparison of HE and Gomori)
      "ROPE = ±%.2f  (%.1f%% inside)",                 # Subtitle: Displays the ROPE range 
      rope_fix, rope_pct                               # and the proportion within the posterior distribution
    )
  ) +
  theme_bw(base_size = 13)                             # Simple theme (font size 13)

## 6-B: Bar graph showing sensitivity analysis results with varying ROPE widths
ggplot(
  data.frame(
    width = factor(rope_grid),                         # Set ROPE width as category
    pct   = grid_pct                                   # Percentage of posterior distributions falling within each ROPE width (%)
  ),
  aes(width, pct)
) +
  geom_col(fill = "#3182BD", alpha = .6, width = .6) +  # Display bar graph (blue with transparency 0.6)
  geom_text(
    aes(label = sprintf("%.1f%%", pct)),               # Display percentage above bar graph
    vjust = -0.5, size = 4
  ) +
  labs(
    x = "ROPE width (Cohen's d)",                      # X-axis label: ROPE width (Cohen's d)
    y = "Posterior mass (%)",                          # Y-axis label: Percentage of posterior distribution (%)
    title = "Sensitivity analysis across ROPE widths"  # Title: Sensitivity analysis across ROPE widths
  ) +
  ylim(0, 100) +                                       # Set Y-axis range to 0-100%
  theme_bw(base_size = 12)                             # Simple theme (font size 12)




###  [Result Summary: H&E vs Gomori (Staining Method Comparison)] ###
#
# Study phase & dataset
# - Pilot analysis conducted at the outset of the project.
# - Limited subset of n = 28 paired sections from the same specimens.

# Why Myofiber Area (MFA) was chosen
# - Parameters most likely to reveal staining artefacts are MFA 
#   and CFA because trichrome stains collagen differently from H&E.
# - MFA is objective and readily quantifiable; at this stage,
#   LUT-based digital restoration and subjective indices 
#   (e.g., internal nuclei) were not yet feasible.

# Statistical comparison (MFA, paired H&E vs Gomori)
#   ---Test---	                      ---Result---	       ---Interpretation---
#   Wilcoxon signed-rank	            p = 0.9288	         No location shift
#   Paired t-test	                    p = 0.9137	         Means equivalent
#   Hedges’ d (bootstrapped 95 % CI)	0.03 [–0.48, 0.57]	 Trivial effect

# Implications for CFA
# - MFA and CFA show a very strong inverse correlation (r = –0.96).
# - Therefore CFA is likewise expected to exhibit 
#   minimal stain-dependent variability; dedicated CFA tests deemed unnecessary.

# Conclusion & impact on main study
# - No meaningful difference exists between H&E and G-T for MFA.
# - Data from both staining techniques can be pooled on a common scale
#   throughout the study.
# - The finding validates MFA as a consistent, reliable histological metric
#   across stains while other parameters awaited further methodological development.
#  
#
# ---[Bayesian approach summary]---
#  ● JZS Bayesian paired t‑test
#     - t = 0.626
#     - BF10 = 0.101 ← The no difference model is favoured over the difference model
#       by 1/0.101 ≈ 9.9 times (almost very strong)
#       (Jeffreys scale: BF01 between 3 and 10 = “moderate evidence for the null”)
#
#  ● Posterior distribution (effect size d)
#     - Mean d = 0.12 (essentially zero)
#     - 95% HDI = –0.26 to 0.49
#
#  ● ROPE decision (Lakens, 2022 recommendation: |d| < 0.30)
#     - 82.1 % of the posterior lies within ROPE ±0.30
#       → Practical equivalence (no meaningful difference)
#     - Sensitivity: varying ROPE width yields the same conclusion
#       (±0.25 → 73 %, ±0.40 → 93 %)
#
#  ● Batch effect (specimen‑level variance τ²)
#     - τ²_0 (no stain adjustment)  = 165.3584   (unit: (%pt)²)
#       sqrt(τ²) ≈ 12.9 %pt → natural inter‑specimen variability
#     - τ²_1 (stain as fixed effect) = 165.3025
#       → Δτ² = 0.03 %: staining does not alter inter‑specimen variability
#
#  Conclusion:
#  - Bayes factor and ROPE both indicate no substantive difference
#    between H&E and Gomori/Masson trichrome for MFA measurements.
#  - Stain‑related batch effects are negligible, so both staining methods
#    can be treated on a common scale in this study.
###------------------------------### 1. Correlation Analysis--------------------------------

# Test the relationship between each parameter and ABx using cor.test(),
# and draw scatter plots with ggplot2 + geom_smooth() where needed.
#
# ★Comment style★
#   ・ For cor.test() lines, attach the result as a comment (e.g., rho = -0.845..).
#   ・ In the ggplot section, explain the main geoms and theme settings line by line.

# --- 7.1 Representative example: MFD vs ABx ---
cor.test(ABx, MFD, method = "spearman")  # rho = -0.8454973
cor.test(ABx, MFD)                        # Pearson r = -0.800855
cor.test(ABx, log(MFD))                   # after log-transformation (MFD)
cor.test(log(ABx), log(MFD))              # transform both parameter

# qplot() has been deprecated since 3.4.0, but the legacy code is retained 
# -> show the correlation coefficient with annotate()

qplot(ABx, MFD) +                               # Scatter plot (simple API)
  geom_point(size = 2, colour = "black") +      # Point size and color settings
  annotate("text", x = 11, y = 1500,            # In-plot text (displaying the correlation coefficient)
           size = 10, fontface = "bold", 
           label = "rho = -0.85") +
  geom_smooth(span = 1, colour = "red", se = FALSE) +  # LOESS curve (span = 1)
  scale_x_continuous(limits = c(0, 16.5), breaks = 0:16) +  # X-axis range and tick marks
  theme_bw() +                                             # White-background theme
  theme(
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.text.x  = element_text(size = 16),
    axis.text.y  = element_text(size = 16)
  )

# --- 7.2 MFD (subset with age < 11 years) ---
L_sub  <- subset(L, ABx < 11)                 # Filter to age < 11 years
ABx_sub <- L_sub$ABx
MFD_sub <- L_sub$MFD

cor.test(ABx_sub, MFD_sub, method = "spearman")  # rho = -0.810084

data_sub <- data.frame(ABx_sub, MFD_sub)

ggplot(data_sub, aes(x = ABx_sub, y = MFD_sub)) +
  geom_point(size = 2, colour = "black") +
  scale_x_continuous(limits = c(0, 11.5), breaks = 0:11) +
  # geom_smooth(span = 0.9, colour = "red", se = FALSE) + optional smooth line
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 40, face = "bold"),
    axis.title.y = element_text(size = 40, face = "bold"),
    axis.text.x  = element_text(size = 16),
    axis.text.y  = element_text(size = 16)
  )

# PART 2: Mean-------
# Detailed comments are provided here for Mean. 
# The remaining parameters are evaluated with the same structure.

# (a) Spearman correlation between Mean and ABx
cor.test(ABx, Mean, method = "spearman")  # rho = 0.6124302

# (b) Scatter plot
#    - geom_point(): data points
#    - scale_x_continuous(): fix the x-axis to 0-16.5 years, with ticks every 1 year
#    - theme_bw(): white background
#    - theme(): adjust font size and boldface

ggplot(dataAll, aes(x = ABx, y = Mean)) +
  geom_point(size = 2, colour = "black") +
  scale_x_continuous(limits = c(0, 16.5), breaks = 0:16) +
  # geom_smooth(span = 1, colour = "red", se = FALSE)  # Smooth line (optional)
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 40, face = "bold"),
    axis.title.y = element_text(size = 40, face = "bold"),
    axis.text.x  = element_text(size = 16),
    axis.text.y  = element_text(size = 16)
  )

# >>> The same pattern is used below for Sd, Cov, MFA, CFA, Fat, IntN, Opaque, NFA, and RFA <<<
#   ・ Check the statistics with cor.test()
#   ・ Draw scatter plots with ggplot() (with smooth lines where needed)
#   ・ Use comments to explain the main settings (size, color, scale, etc.)

# By copying this template and replacing only the variable names,
# any parameter can be visualized quickly.。

# PART 3: Sd--------------------------------------------------------------------
cor.test(ABx, Sd, method = "spearman") # 0.7369515 
ggplot(dataAll, aes(x = ABx, y = Sd)) +
  geom_point(size = 2, colour = "black") +
  scale_x_continuous(limits = c(0, 16.5), breaks = 0:16.5) +
  #geom_smooth(span = 1, colour = "red", se = FALSE) +
  theme_bw()+
  theme(
    axis.title.x = element_text(size = 40, face = "bold"), 
    axis.title.y = element_text(size = 40, face = "bold"), 
    axis.text.x = element_text(size = 16), 
    axis.text.y = element_text(size = 16)  
  )

# PART 4: Cov-------------------------------------------------------------------
cor.test(ABx, Cov, method = "spearman") # 0.6575118 
ggplot(dataAll, aes(x = ABx, y = Cov)) +
  geom_point(size = 2, colour = "black") +
  scale_x_continuous(limits = c(0, 16.5), breaks = 0:16.5) +
  #geom_smooth(span = 1, colour = "red", se = FALSE) +
  theme_bw()+
  theme(
    axis.title.x = element_text(size = 40, face = "bold"), 
    axis.title.y = element_text(size = 40, face = "bold"), 
    axis.text.x = element_text(size = 16), 
    axis.text.y = element_text(size = 16)  
  )

# PART 5: MFD-------------------------------------------------------------------
cor.test(ABx, MFD, method = "spearman") # -0.8454973

ggplot(dataAll, aes(x = ABx, y = MFD)) +
  geom_point(size = 2, colour = "black") +
  scale_x_continuous(limits = c(0, 16.5), breaks = 0:16.5) +
  geom_smooth(span = 1, colour = "red", se = FALSE) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 20, face = "bold.italic"), 
    axis.title.y = element_text(size = 20, face = "bold.italic"), 
    axis.text.x = element_text(size = 16), 
    axis.text.y = element_text(size = 16)  
  )

ggplot(dataAll, aes(x = ABx, y = MFD)) +
  scale_x_continuous(limits = c(0, 16.5), breaks = 0:16.5) +
  geom_smooth(span = 1, colour = "black", se = FALSE) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 20, face = "bold.italic"), 
    axis.title.y = element_text(size = 20, face = "bold.italic"), 
    axis.text.x = element_text(size = 16), 
    axis.text.y = element_text(size = 16)  
  )

L_sub <- subset(L, ABx < 11)

ABx_sub <- L_sub$ABx
MFD_sub <- L_sub$MFD

cor.test(ABx_sub, MFD_sub, method = "spearman") # rho = -0.810084

data_sub <- data.frame(ABx_sub, MFD_sub)
ggplot(data_sub, aes(x = ABx_sub, y = MFD_sub)) +
  geom_point(size = 2, colour = "black") +
  scale_x_continuous(limits = c(0, 11.5), breaks = 0:11.5) +
  #geom_smooth(span = 0.9, colour = "red", se = FALSE) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 40, face = "bold"), 
    axis.title.y = element_text(size = 40, face = "bold"), 
    axis.text.x = element_text(size = 16), 
    axis.text.y = element_text(size = 16)  
  )

ggplot(data_sub, aes(x = ABx_sub, y = MFD_sub)) +
  scale_x_continuous(limits = c(0, 11.5), breaks = 0:11.5) +
  geom_smooth(span = 0.9, colour = "black", se = FALSE) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 40, face = "bold"), 
    axis.title.y = element_text(size = 40, face = "bold"), 
    axis.text.x = element_text(size = 16), 
    axis.text.y = element_text(size = 16)  
  )
# PART 6: MFA-------------------------------------------------------------------
cor.test(ABx, MFA, method = "spearman") # -0.6032389 
ggplot(dataAll, aes(x = ABx, y = MFA)) +
  geom_point(size = 2, colour = "black") +
  scale_x_continuous(limits = c(0, 16.5), breaks = 0:16.5) +
  #geom_smooth(span = 1, colour = "red", se = FALSE) +
  theme_bw()+
  theme(
    axis.title.x = element_text(size = 40, face = "bold"), 
    axis.title.y = element_text(size = 40, face = "bold"), 
    axis.text.x = element_text(size = 16), 
    axis.text.y = element_text(size = 16)  
  )

# PART 7: CFA-------------------------------------------------------------------
cor.test(ABx, CFA, method = "spearman") # -0.483751  
ggplot(dataAll, aes(x = ABx, y = CFA)) +
  geom_point(size = 2, colour = "black") +
  scale_x_continuous(limits = c(0, 16.5), breaks = 0:16.5) +
  #geom_smooth(span = 1, colour = "red", se = FALSE) +
  theme_bw()+
  theme(
    axis.title.x = element_text(size = 40, face = "bold"), 
    axis.title.y = element_text(size = 40, face = "bold"), 
    axis.text.x = element_text(size = 16), 
    axis.text.y = element_text(size = 16)  
  )

# PART 8: Fat-------------------------------------------------------------------
cor.test(ABx, Fat, method = "spearman") # 0.690557   
ggplot(dataAll, aes(x = ABx, y = Fat)) +
  geom_point(size = 2, colour = "black") +
  scale_x_continuous(limits = c(0, 16.5), breaks = 0:16.5) +
  #geom_smooth(span = 1, colour = "red", se = FALSE) +
  theme_bw()+
  theme(
    axis.title.x = element_text(size = 40, face = "bold"), 
    axis.title.y = element_text(size = 40, face = "bold"), 
    axis.text.x = element_text(size = 16), 
    axis.text.y = element_text(size = 16)  
  )

# PART 9: IntN------------------------------------------------------------------
cor.test(ABx, IntN, method = "spearman") # -0.3034249 NS 
ggplot(dataAll, aes(x = ABx, y = IntN)) +
  geom_point(size = 2, colour = "black") +
  scale_x_continuous(limits = c(0, 16.5), breaks = 0:16.5) +
  #geom_smooth(span = 1, colour = "red", se = FALSE) +
  theme_bw()+
  theme(
    axis.title.x = element_text(size = 40, face = "bold"), 
    axis.title.y = element_text(size = 40, face = "bold"), 
    axis.text.x = element_text(size = 16), 
    axis.text.y = element_text(size = 16)  
  )

# PART 10: Opaque---------------------------------------------------------------
cor.test(ABx, Opaque, method = "spearman") # -0.3439107  
ggplot(dataAll, aes(x = ABx, y = Opaque)) +
  geom_point(size = 2, colour = "black") +
  scale_x_continuous(limits = c(0, 16.5), breaks = 0:16.5) +
  #geom_smooth(span = 1, colour = "red", se = FALSE) +
  theme_bw()+
  theme(
    axis.title.x = element_text(size = 40, face = "bold"), 
    axis.title.y = element_text(size = 40, face = "bold"), 
    axis.text.x = element_text(size = 16), 
    axis.text.y = element_text(size = 16)  
  )

# PART 11: NFA------------------------------------------------------------------
cor.test(ABx, NFA, method = "spearman") # 0.2719116 NS  
ggplot(dataAll, aes(x = ABx, y = NFA)) +
  geom_point(size = 2, colour = "black") +
  scale_x_continuous(limits = c(0, 16.5), breaks = 0:16.5) +
  #geom_smooth(span = 1, colour = "red", se = FALSE) +
  theme_bw()+
  theme(
    axis.title.x = element_text(size = 40, face = "bold"), 
    axis.title.y = element_text(size = 40, face = "bold"), 
    axis.text.x = element_text(size = 16), 
    axis.text.y = element_text(size = 16)  
  )

# PART 12: RFA------------------------------------------------------------------
cor.test(ABx, RFA, method = "spearman") # -0.2119488 NS
ggplot(dataAll, aes(x = ABx, y = RFA)) +
  geom_point(size = 2, colour = "black") +
  scale_x_continuous(limits = c(0, 16.5), breaks = 0:16.5) +
  #geom_smooth(span = 1, colour = "red", se = FALSE) +
  theme_bw()+
  theme(
    axis.title.x = element_text(size = 40, face = "bold"), 
    axis.title.y = element_text(size = 40, face = "bold"), 
    axis.text.x = element_text(size = 16), 
    axis.text.y = element_text(size = 16)  
  )


### -----------------------------### 2. Multiple Regression------------------------
# PART 1: Main multiple regression----------------------------------------------

# We initially hypothesized that if a critical age window (Peverelli et al., 2015,
# Neurology) did not exist, ABx could be predicted smoothly and continuously from 
# morphological variables across all samples.

# Given the lack of established quantitative markers, we aimed to identify 
# potential predictive parameters through an unbiased, hypothesis-generating approach.

Simple <- lm(ABx.C ~ MFD.C)
summary(Simple) # Adjusted R-squared: 0.7318

Attempt <- lm(ABx.C ~ (MFD.C + CFA.C + Sd.C + MFA.C + Fat.C + IntN.C + Cov.C)^2)
step(Attempt)
Attempt_Reg <- lm(ABx.C ~ MFD.C + CFA.C + Sd.C + MFA.C + Fat.C + IntN.C + 
                    Cov.C + MFD.C:CFA.C + MFD.C:Sd.C + MFD.C:MFA.C + MFD.C:Fat.C + 
                    CFA.C:Sd.C + CFA.C:Fat.C + CFA.C:IntN.C + CFA.C:Cov.C + Sd.C:MFA.C + 
                    Sd.C:Fat.C + Sd.C:IntN.C + Sd.C:Cov.C + MFA.C:Fat.C + MFA.C:IntN.C + 
                    MFA.C:Cov.C + Fat.C:IntN.C + Fat.C:Cov.C + IntN.C:Cov.C)
summary(Attempt_Reg)

Attempt2 <- lm(ABx.C ~ (MFD.C + Sd.C + Cov.C)^2)
step(Attempt2) 
Attempt_Reg2 <- lm(ABx.C ~ MFD.C + Sd.C + MFD.C:Sd.C)
summary(Attempt_Reg2)


# Sd log-transformed linear regression
SD_Lreg <- lm(ABx.C ~ MFD.C + MFD.C:Sd.C)
summary(SD_Lreg) # RSE: 0.4609, Adjusted R-squared:  0.7876, p-value: 6.343e-13
confint(SD_Lreg)
AIC(SD_Lreg) # 53.83862
autoplot(SD_Lreg, smooth.colour = NA)
vif(SD_Lreg, type="predictor") 
test_anovaSd <- lm(ABx.C ~ MFD.C) 
anova(test_anovaSd, SD_Lreg) #Pr(>F): 0.002669**


# ABx prediction
{
  AnswerSD_reg <- lm(ABx.c ~ MFD.c + Sd.c:MFD.c)
  summary(AnswerSD_reg) 
  AnswerSD <- predict(AnswerSD_reg); AnswerSD 
  Log_ABxSD = AnswerSD + mean(log(ABx)); Log_ABxSD
  ABx_PredictSD = exp(Log_ABxSD); ABx_PredictSD
}


# Scatter plot comparing ABx (observed age) and ABx_PredictSD (predicted age) 
# with standard error

plot_ABx_PredictSD <- qplot(ABx, ABx_PredictSD) +  # Create a scatter plot with ABx on the x-axis and ABx_PredictSD on the y-axis
  geom_point(size = 6, colour = "#0072B2", alpha = 0.8, shape = 16) + # Display data points in blue with size 6 (alpha = 0.8)
  scale_x_continuous(limits = c(0, 17), breaks = 0:17) + # Set the x-axis range to 0-17 with tick marks every 1 unit
  scale_y_continuous(limits = c(0, 17), breaks = 0:17) + # Apply the same settings to the y-axis (0-17)
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", 
              color = "black", size = 0.7) + # Add the identity line (y = x) as a black dashed line
  theme_bw() + # Apply a simple white-background theme
  theme(legend.position = "top") + # Place the legend at the top (although this plot has no legend)
  coord_fixed(); print(plot_ABx_PredictSD) + # Fix the x/y aspect ratio and print the plot
  
  geom_smooth(method = "lm", colour = "red", alpha = 0.5, se = T, 
              linetype = "solid") + # Add a red linear regression line (with a semi-transparent 95% CI)）
  
  theme(
    axis.title.x = element_text(size = 30, face = "bold"), # Display the x-axis title in bold with size 30
    axis.title.y = element_text(size = 30, face = "bold"), # Apply the same settings to the y-axis title
    axis.text.x = element_text(size = 15, face = "bold"),  # Display x-axis tick labels in bold with size 15
    axis.text.y = element_text(size = 15, face = "bold")   # Apply the same settings to the y-axis tick labels
  )



# Plot showing the 95% prediction interval of ABx_PredictSD (predicted age) against ABx (observed age)
{
  # Calculate the predicted values and 95% prediction interval of ABx_PredictSD using a linear model
  predict_PredictABx <- predict(lm(ABx_PredictSD ~ ABx, data = L),
                                interval = "predict", level = 0.95)
  
  # Extract the lower and upper bounds of the 95% prediction interval
  predict_PredictABx_lower <- predict_PredictABx[, "lwr"]
  predict_PredictABx_upper <- predict_PredictABx[, "upr"]
  
  # Calculate half of the prediction interval width for each data point
  predict_intervals_half_width <- (predict_PredictABx_upper - predict_PredictABx_lower) / 2
  
  # Compute the mean half-width of the prediction interval
  mean_predict_interval_half_width <- mean(predict_intervals_half_width)
  
  # Draw the scatter plot, regression line, and prediction interval with ggplot
  plotABx_predictSD <- ggplot(L, aes(x = ABx)) +
    geom_point(aes(y = ABx_PredictSD, color = "Predicted ABx"),
               size = 6, alpha = 0.8, shape = 16) +  # Display predicted age as blue points
    geom_smooth(aes(y = ABx_PredictSD), method = "lm",
                color = "red", size = 2, se = FALSE) +  # Display the regression line in red (without CI)
    geom_ribbon(aes(ymin = predict_PredictABx_lower, ymax = predict_PredictABx_upper),
                fill = "gray", alpha = 0.3) +  # Display the 95% prediction interval as a gray ribbon
    scale_color_manual(values = c("#0072B2")) +  # Set the point color
    scale_shape_manual(values = c("ABx_predictSD" = 15)) +  # Set the point shape (has no effect in this plot)
    theme_bw() +  # Apply a white-background theme
    theme(legend.position = "top") +  # Place the legend at the top (only one legend in this case)
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +  # Add the reference line y = x
    scale_x_continuous(limits = c(0, 16.5), breaks = seq(0, 16, 1)) +  # Set the x-axis range to 0-16.5 with 1-unit ticks
    scale_y_continuous(limits = c(-1, 18), breaks = seq(0, 18, 1)) +  # Set the y-axis range to -1-18 with 1-unit ticks
    coord_cartesian(clip = "off") +  # Allow drawing outside the plot region
    labs(x = "ABx", y = "Value")  # Set axis labels
  
  # Print the plot
  print(plotABx_predictSD)
  
  # Display the mean half-width of the prediction interval
  cat("Mean prediction interval half-width (ABx scale):",
      mean_predict_interval_half_width, "\n") # Example output: 3.364948
  }


# exploratory multiple regression (replace: MFD only)
LMFD_ATP <- log(MFD_ATP)
MFD_ATP.c <- LMFD_ATP - mean(LMFD_ATP)

ATP_Reg <- lm(ABx.c ~ MFD_ATP.c + Sd.c:MFD_ATP.c)
summary(ATP_Reg) # Adjusted R-squared: 0.8113
autoplot(ATP_Reg, smooth.colour = NA)
confint(ATP_Reg)
test_anova_ATP <- lm(ABx.c ~ MFD_ATP.c) 
anova(test_anova_ATP, ATP_Reg) #Pr(>F): 0.0008454***


# PART 2: MSE calculation by LOOCV (Leave-One-Out Cross-Validation) for 38 cases ------
{
  # Data preparation
  dataL17 <- L1 # Use all 38 cases here (data for age < 7 years will be calculated later)
  
  # Define a function for LOOCV (calculates MSE)
  loocv_mse <- function(dataL17) {
    n <- nrow(dataL17)                  # Number of samples (38 here)
    mse_values <- numeric(n)            # Prepare a vector to store the MSE for each case
    
    # LOOCV loop (build the model while leaving out one case at a time)
    for (i in 1:n) {
      training_data <- dataL17[-i, ]    # Use all remaining data except the excluded case as training data
      testing_data  <- dataL17[i, , drop = FALSE] # Use the excluded case as test data
      
      # Build a linear regression model 
      # (predict ABx.C from MFD.C and the interaction between MFD.C and Sd.C)
      model <- lm(ABx.C ~ MFD.C + MFD.C:Sd.C, data = training_data)
      
      # Calculate the predicted value for the excluded case
      predicted_value <- predict(model, newdata = testing_data)
      
      # Compute and store the squared prediction error (MSE)
      mse_values[i] <- (testing_data$ABx.C - predicted_value)^2
    }
    
    # Return the MSE for each case
    return(mse_values)
  }
  
  # Run LOOCV and obtain the MSE for each case
  loocv_mse_values <- loocv_mse(dataL17)
  
  # Calculate the mean and standard deviation of the MSE values obtained by LOOCV
  mean_loocv_mse <- mean(loocv_mse_values)
  sd_loocv_mse   <- sd(loocv_mse_values)
  
  # Display the results
  cat("LOOCV MSE's mean:", mean_loocv_mse, "\n")          # Mean MSE (example: 0.2379219)
  cat("LOOCV RMSE's mean:", sqrt(mean_loocv_mse), "\n")   # Mean RMSE (square root of MSE; example: 0.4877724)
  cat("LOOCV MSE's SD:", sd_loocv_mse, "\n")              # Standard deviation of MSE (example: 0.3199825)
}


# PART 3: Additional model(Mean)------------------------------------------------
RegMean <- lm(ABx.C ~ MFD.C + MFD.C:Mean.C)
summary(RegMean) # RSE: 0.4544, Adjusted R-squared:  0.7935
autoplot(RegMean, smooth.colour = NA) # Cook's distance > 0.5
AIC(RegMean) # 52.76904
confint(RegMean)
vif(RegMean, type="predictor") 
test_anovaMean <- lm(ABx.C ~ MFD.C) 
anova(test_anovaMean, RegMean) #Pr(>F): 0.001572**

###------------------------------### 3. ABx-limited Multiple Regression Models-----
# PART 1: Data frame------------------------------------------------------------
{
  # data frame
  data_c <- data.frame(ABx.c = ABx.c, MFD.c = MFD.c, Sd.c = Sd.c, Cov.c = Cov.c, 
                       CFA.c = CFA.c, MFA.c = MFA.c, Opaque.c = Opaque.c,
                       NFA.c = NFA.c, RFA.c = RFA.c, Fat.c = Fat.c, 
                       Mean.c = Mean.c, IntN.c = IntN.c)
  
  Data_C <- data.frame(ABx.C = ABx.C, MFD.C = MFD.C, Sd.C = Sd.C, Cov.C = Cov.C, 
                       CFA.C = CFA.C, MFA.C = MFA.C, Opaque.C = Opaque.C,
                       NFA.C = NFA.C, RFA.C = RFA.C, Fat.C = Fat.C, 
                       Mean.C = Mean.C, IntN.C = IntN.C)
}

# PART 2: < 11yo----------------------------------------------------------------
{
  data_under_11 <- data_c[data_c$ABx.c < 0.6, ] # see "ABx" and "ABx.c"; ABx < 11 vs ABx.c < 0.6 
  Data_under_11 <- Data_C[Data_C$ABx.C < 1.2, ] # see "ABx" and "ABx.C"; ABx < 11 vs ABx.c < 1.2
  data_under11ABx <- L[L$ABx < 11, ]
  ABx11 = data_under11ABx$ABx
}

RegSD11 <- lm(ABx.C ~ MFD.C + MFD.C:Sd.C, data = Data_under_11)
summary(RegSD11) # Adjusted R-squared:  0.7644
confint(RegSD11)

# ABx prediction < 11yo
{
  AnswerSD_reg11 <- lm(ABx.c ~ MFD.c + Sd.c:MFD.c, data = data_under_11)
  summary(AnswerSD_reg11) 
  AnswerSD11 <- predict(AnswerSD_reg11); AnswerSD11 
  Log_ABxSD11 = AnswerSD11 + mean(log(ABx)); Log_ABxSD11
  ABx_PredictSD11 = exp(Log_ABxSD11); ABx_PredictSD11
}

# Validate ABx prediction < 9yo (1-6,6-7,7-9yo) in the same manner.
# PART 3: < 9yo-----------------------------------------------------------------
{
  data_under_9 <- data_c[data_c$ABx.c < 0.39, ]   
  Data_under_9 <- Data_C[Data_C$ABx.C < 0.69, ]  
  data_under9ABx <- L[L$ABx < 9, ]  
  ABx9 = data_under9ABx$ABx  
}

RegSD9 <- lm(ABx.C ~ MFD.C + MFD.C:Sd.C, data = Data_under_9)
summary(RegSD9) # Adjusted R-squared: 0.7684 

{
  AnswerSD_reg9 <- lm(ABx.c ~ MFD.c + Sd.c:MFD.c, data = data_under_9)
  summary(AnswerSD_reg9)  # Adjusted R-squared:  0.7684
  AnswerSD9 <- predict(AnswerSD_reg9); AnswerSD9 
  Log_ABxSD9 = AnswerSD9 + mean(log(ABx)); Log_ABxSD9
  ABx_PredictSD9 = exp(Log_ABxSD9); ABx_PredictSD9
}


# PART 4: < 7yo-----------------------------------------------------------------
{
  data_under_7 <- data_c[data_c$ABx.c < 0.14, ]
  Data_under_7 <- Data_C[Data_C$ABx.C < 0.24, ]
  data_under7ABx <- L[L$ABx < 7, ]
  ABx7 = data_under7ABx$ABx
}

RegSD7 <- lm(ABx.C ~ MFD.C + MFD.C:Sd.C, data = Data_under_7)
summary(RegSD7) # RSE: 0.4321, Adjusted R-squared:  0.7449
AIC(RegSD7) # 29.11788
confint(RegSD7)

# ABx prediction < 7yo
{
  AnswerSD_reg7 <- lm(ABx.c ~ MFD.c + Sd.c:MFD.c, data = data_under_7)
  summary(AnswerSD_reg7) 
  AnswerSD7 <- predict(AnswerSD_reg7); AnswerSD7 
  Log_ABxSD7 = AnswerSD7 + mean(log(ABx)); Log_ABxSD7
  ABx_PredictSD7 = exp(Log_ABxSD7); ABx_PredictSD7
}

### -----------------------------### 4. Residual (Predict Error) Analysis----------
# PART 1: Error analysis (Main)--------------

# Purpose:
# To test whether biopsy-age prediction accuracy
# (prediction error = predicted ABx - observed ABx)
# varies across age intervals.
# We specifically examined whether a critical histological window
# exists around 6-7 years, as proposed by Peverelli et al. (2015, Neurology).

# Statistical approach:
# First, prediction errors were compared across broader age intervals
# (1-6, 6-7, 7-11, >=11 years, etc.).
# Because the 11-17 years group was small (n = 3),
# Kruskal-Wallis tests were used for multi-group comparisons.
# These tests consistently showed significant between-group differences (p < 0.05),
# indicating age-dependent heterogeneity in predictive accuracy.

# Next, to focus on the putative critical interval,
# pairwise comparisons were performed for 1-6 years vs. 6-7 years only.
# Welch's t-test and the Wilcoxon rank-sum test both showed no significant difference
# in prediction errors between these two groups. 


# Residual definitions
errors17 <- ABx_PredictSD - ABx # Full-range prediction errors
errors11 <- ABx_PredictSD11 - ABx11 # Prediction errors: age < 11 years
errors9 <- ABx_PredictSD9 - ABx9 # Prediction errors: age < 9 years
errors7 <- ABx_PredictSD7 - ABx7 # Prediction errors: age < 7 years

Errors <- data.frame(ABx, errors17) # Full-range residual data
Errors11 <- data.frame(ABx11, errors11) # Residual data: age < 11 years
Errors9 <- data.frame(ABx9, errors9) # Residual data: age < 9 years

# Boxplots: finer age split
{
  age_groups <- cut(ABx, breaks = c(1,4,6,7,9,11,17), right = FALSE, 
                    labels = c("1–4","4–6","6–7","7–9","9–11","11–17"))
  
  boxplot(errors17 ~ age_groups,
          xlab = "Age Group (years)",
          ylab = "Prediction Errors",
          main = "Boxplot of Prediction Errors by Age Group",
          col = "lightblue")
}

# Boxplots: broader age split
{
  age_groups_2_1 <- cut(ABx, breaks = c(1,6,7,11,17), right = FALSE, 
                        labels = c("1–6","6–7","7–11","11–17"))
  
  boxplot(errors17 ~ age_groups_2_1,
          xlab = "Age Group (years)",
          ylab = "Prediction Errors",
          main = "Prediction Errors by Age Groups",
          col = "lightblue")
}


# Three boxplots shown side by side
{
  par(mfrow = c(1,3))  # 1 row x 3 columns
  {
    age_groups_2_2 <- cut(ABx11, breaks = c(1,6,7,11), right = FALSE, 
                          labels = c("1–6","6–7","7–11"))
    
    boxplot(errors11 ~ age_groups_2_2,
            xlab = "Age Group (years)",
            ylab = "Prediction Errors",
            col = "lightblue")
  }
  
  {
    age_groups_3 <- cut(ABx9, breaks = c(1,6,7,9), right = FALSE, 
                        labels = c("1–6","6–7","7–9"))
    
    boxplot(errors9 ~ age_groups_3,
            xlab = "Age Group (years)",
            ylab = "Prediction Errors",
            main = "Boxplot of Prediction Errors by Age Group",
            col = "lightblue")
  }
  
  {
    age_groups_4 <- cut(ABx7, breaks = c(1,6,7), right = FALSE, 
                        labels = c("1–6","6–7"))
    
    boxplot(errors7 ~ age_groups_4,
            xlab = "Age Group (years)",
            ylab = "Prediction Errors",
            col = "lightblue")
  }
  
  
  par(mfrow = c(1,1))  # Reset layout
}


# Global nonparametric tests by age grouping
ABxGroups17 <- cut(ABx, breaks = c(1,6,7,11,17), right = F)  # 4-group split
kruskal.test(errors17 ~ ABxGroups17) # p-value = 0.002931 (1,6,7,11,17)

ABxGroups17_2 <- cut(ABx, breaks = c(1,6,7,9,11,17), right = F)  # 5-group split
kruskal.test(errors17 ~ ABxGroups17_2) # p-value = 0.007006 (1,6,7,9,11,17)

ABxGroups11 <- cut(ABx11, breaks = c(1,6,7,11), right = F)  # 3-group split: age < 11
kruskal.test(errors11 ~ ABxGroups11) # p-value = 0.001418 (1,6,7,11)
conover_result_Errors11 <- kwAllPairsConoverTest(errors11 ~ ABxGroups11, data = Errors11, p.adjust.method = "bonferroni")
print(conover_result_Errors11)  # Post-hoc pairwise comparison

ABxGroups11_2 <- cut(ABx11, breaks = c(1,6,7,9,11), right = F)  # 4-group split: age < 11
kruskal.test(errors11 ~ ABxGroups11_2) # p-value = 0.003798 (1,6,7,9,11)
conover_result_Errors11_2 <- kwAllPairsConoverTest(errors11 ~ ABxGroups11_2, data = Errors11, p.adjust.method = "bonferroni")
print(conover_result_Errors11_2)  # Post-hoc pairwise comparison

ABxGroups9 <- cut(ABx9, breaks = c(1,6,7,9), right = F)  # 3-group split: age < 9
kruskal.test(errors9 ~ ABxGroups9) # p-value = 0.002122
conover_result_Errors9 <- kwAllPairsConoverTest(errors9 ~ ABxGroups9, data = Errors9, p.adjust.method = "bonferroni")
print(conover_result_Errors9)  # Post-hoc pairwise comparison

# PART 2: Error analysis (<7yo)----------------------------------------------
ABxGroups7 <- cut(ABx7, breaks = c(1,6,7), right =F)  # Two-group split: 1-6 vs. 6-7

group1_error7 <- errors7[ABxGroups7 == "[1,6)"]  # Prediction errors: 1-6 years
group2_error7 <- errors7[ABxGroups7 == "[6,7)"]  # Prediction errors: 6-7 years

shapiro.test(group1_error7) # p-value = 0.3177
shapiro.test(group2_error7) # p-value = 0.8853

var.test(group1_error7, group2_error7) # p-value = 0.02343*

t.test(group1_error7, group2_error7, var.equal = FALSE) # p-value = 0.5515, CI: -0.4892820 to 0.8853815
wilcox.test(group1_error7, group2_error7) # p-value = 0.3011

cliff.delta(group1_error7, group2_error7) # 0.2884615: small (-0.2404896 to 0.6853153)


# Interpretation:
# Multi-group comparisons including older patients showed significant heterogeneity,
# whereas the direct two-group comparison (1-6 vs. 6-7 years) did not.
# Notably, restricting analysis to patients under 7 years reduced
# the variability in prediction errors, suggesting more homogeneous
# histopathological dynamics before this age range.

# Effect size:
# Cliff's delta indicated only a small difference between 1-6 and 6-7 years.
# Taken together, these findings support cautious interpretation:
# the broader residual pattern is compatible with a transition around 6-7 years,
# but the isolated two-group contrast was not statistically significant.
# Therefore, the age intervals (1-6, 6-7, 7-11 years) were adopted for subsequent
# ANOVA/Kruskal-Wallis analyses and effect size comparisons of histopathological parameters.


boxplot(errors7 ~ cut(ABx7, breaks = c(1, 6, 7), 
                      labels = c("[1-6)","[6-7)"), 
                      right = FALSE), 
        xlab = "ABx Groups7",
        cex.axis = 1.5,
        ylab = "errors7",
        cex.lab = 1.5,
        main = "",
        outline = FALSE)

# Plot + prediction error: prediction interval plot using only data under 7 years
{
  # Compute the 95% prediction interval from a linear model fitted to the under-7 dataset
  predict_PredictABx7 <- predict(lm(ABx_PredictSD7 ~ ABx7), interval = "predict", level = 0.95)
  predict_PredictABx7_lower <- predict_PredictABx7[, "lwr"]  # Lower bound
  predict_PredictABx7_upper <- predict_PredictABx7[, "upr"]  # Upper bound
  
  # Compute the half-width of the prediction interval
  predict_intervals_half_width7 <- (predict_PredictABx7_upper - predict_PredictABx7_lower)/2
  mean_predict_interval_half_width7 <- mean(predict_intervals_half_width7)  # 平均半幅
  
  # Build the plot with ggplot
  plotABx_predictSD7 <- ggplot(data_under7ABx, aes(x = ABx)) +
    
    # Plot predicted values from the under-7 dataset
    geom_point(aes(y = ABx_PredictSD7, color = "ABx_PredictSD7"), size = 6, alpha = 0.8, shape = 16) +
    
    # Regression line for the under-7 dataset
    geom_smooth(aes(y = ABx_PredictSD7), method = "lm", color = "red", size = 2, se = FALSE) +
    
    # 95% prediction interval shown as a gray ribbon
    geom_ribbon(aes(ymin = predict_PredictABx7_lower, ymax = predict_PredictABx7_upper), fill = "gray", alpha = 0.3) +
    
    # Color and shape settings
    scale_color_manual(values = c("#0072B2")) +
    scale_shape_manual(values = c("ABx_predictSD7" = 15)) +
    
    # White-background theme
    theme_bw() +
    theme(legend.position = "top") +  # Place the legend at the top
    
    # Identity line (ideal agreement)
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
    
    # Axis scales
    scale_x_continuous(limits = c(0, 7), breaks = seq(0, 7, 1)) +      # X-axis: age
    scale_y_continuous(limits = c(-1, 8.5), breaks = seq(0, 8.5, 1)) + # Y-axis: predicted age
    coord_cartesian(clip = "off") + # Allow drawing beyond the panel
    
    # Axis labels
    labs(x = "ABx", y = "Value")
  
  # Print the plot
  print(plotABx_predictSD7)
  
  # Print the mean half-width of the prediction interval
  cat("Mean prediction interval half_width7 (ABx scale):", mean_predict_interval_half_width7, "\n") # 1.823573
}


# Combined plot: full-range and under-7 prediction interval plots shown together
{
  # Compute the 95% prediction interval for the full-range dataset (L)
  predict_intervals_FullRange <- predict(lm(ABx_PredictSD ~ ABx, data = L), interval = "predict", level = 0.95)
  predict_intervals_Full_lower <- predict_intervals_FullRange[, "lwr"]  # Lower bound
  predict_intervals_Full_upper <- predict_intervals_FullRange[, "upr"]  # Upper bound
  
  # Compute the 95% prediction interval for the under-7 dataset
  predict_intervals_Under7Range <- predict(lm(ABx_PredictSD7 ~ ABx7, data = data_under7ABx), interval = "predict", level = 0.95)
  predict_intervals_Under7_lower <- predict_intervals_Under7Range[, "lwr"]  # Lower bound
  predict_intervals_Under7_upper <- predict_intervals_Under7Range[, "upr"]  # Upper bound
  
  # Build the combined plot with ggplot
  Combined_plot <- ggplot(mapping = aes(x = ABx)) + 
    
    # Full-range data (L)
    geom_point(data = L, aes(y = ABx_PredictSD, color = "ABx_PredictSD", shape = "ABx_PredictSD"), size = 2, alpha = 0.8) +  # Full-range points
    geom_smooth(data = L, aes(y = ABx_PredictSD), method = "lm", color = "#0072B2", size = 0.5, se = FALSE) +  # Full-range regression line
    geom_ribbon(data = L, aes(ymin = predict_intervals_Full_lower, ymax = predict_intervals_Full_upper), fill = "#0072B2", alpha = 0.3) +  # Full-range prediction ribbon
    
    # Under-7 data
    geom_point(data = data_under7ABx, aes(x = ABx7, y = ABx_PredictSD7, color = "ABx_PredictSD7", shape = "ABx_PredictSD7"), size = 2, alpha = 0.8) +  # Under-7 points
    geom_smooth(data = data_under7ABx, aes(x = ABx7, y = ABx_PredictSD7), method = "lm", color = "#E69F00", size = 0.5, se = FALSE) +  # Under-7 regression line
    geom_ribbon(data = data_under7ABx, aes(x = ABx7, ymin = predict_intervals_Under7_lower, ymax = predict_intervals_Under7_upper), fill = "#E69F00", alpha = 0.3) +  # Under-7 prediction ribbon
    
    # Legend settings for color and shape
    scale_color_manual(
      name = "Series      ",
      values = c("ABx_PredictSD" = "blue", "ABx_PredictSD7" = "#D55E00"),
      labels = c("ABx_PredictSD" = "Full Range: Blue", "ABx_PredictSD7" = "Under 7: Orange")
    ) +
    scale_shape_manual(
      name = "Series      ",
      values = c("ABx_PredictSD" = 15, "ABx_PredictSD7" = 16),
      labels = c("ABx_PredictSD" = "Full Range: Blue", "ABx_PredictSD7" = "Under 7: Orange")
    ) +
    
    # Plot theme
    theme_bw() +
    theme(
      legend.position = "top",  # Legend at the top
      legend.title = element_text(size = 15),  # Legend title size
      legend.text = element_text(size = 15),   # Legend text size
      axis.title.x = element_text(size = 25, face = "bold"),  # X-axis title style
      axis.title.y = element_text(size = 25, face = "bold"),  # Y-axis title style
      axis.text.x = element_text(size = 15),  # X-axis tick label size
      axis.text.y = element_text(size = 15)   # Y-axis tick label size
    ) +
    
    # Identity line
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +  # Ideal agreement line
    
    # Axis settings
    scale_x_continuous(limits = c(0, 17), breaks = seq(0, 16, 1)) +  # X-axis range and ticks
    scale_y_continuous(limits = c(-1, 17.5), breaks = seq(0, 18, 1)) +  # Y-axis range and ticks
    coord_cartesian(clip = "off") +  # Do not clip elements outside the panel
    
    # Axis labels
    labs(x = "ABx", y = "Predicted ABx")  
  
  # Print the plot
  print(Combined_plot)
}



# PART 3: LOOCV-based MSE calculation for the under-7 dataset---------------------------------------------
{
  # Define the under-7 dataset separately from dataL17
  dataL7 <- Data_under_7
  
  # Function to compute squared prediction errors by LOOCV
  loocv_mse7 <- function(dataL7) {
    n <- nrow(dataL7)                  # Number of observations
    mse_values7 <- numeric(n)          # Vector to store the MSE values
    
    for (i in 1:n) {
      # Split into training and test data by leaving out the i-th case
      training_data7 <- dataL7[-i, ]
      testing_data7 <- dataL7[i, , drop = FALSE]
      
      # Fit the linear regression model
      # Predict ABx.C from MFD.C and its interaction with Sd.C
      model7 <- lm(ABx.C ~ MFD.C + MFD.C:Sd.C, data = training_data7)
      
      # Predict the left-out case
      predicted_value7 <- predict(model7, newdata = testing_data7)
      
      # Store the squared prediction error
      mse_values7[i] <- (testing_data7$ABx.C - predicted_value7)^2
    }
    
    return(mse_values7)  # Return case-wise MSE values
  }
  
  # Run LOOCV and collect the MSE values
  loocv_mse_values7 <- loocv_mse7(dataL7)
  
  # Summary metrics of predictive performance
  mean_loocv_mse7 <- mean(loocv_mse_values7)       # Mean MSE
  sd_loocv_mse7 <- sd(loocv_mse_values7)           # SD of MSE
  
  # Print the results
  cat("LOOCV MSE's mean for dataL7:", mean_loocv_mse7, "\n") # Mean MSE: 0.2252986
  cat("LOOCV RMSE's mean for dataL7:", sqrt(mean_loocv_mse7), "\n") # RMSE: 0.4746563
  cat("LOOCV MSE's SD for dataL7:", sd_loocv_mse7, "\n")     # SD of MSE: 0.3848039
}


### -----------------------------### 5. Overall Interpretation on Preliminary Analyses----

# Our preliminary analyses employing multiple regression modeling aimed at predicting age
# from morphological parameters (particularly focusing on MFD and related parameters)
# served a critical exploratory purpose: to objectively and quantitatively assess whether 
# the relationship between muscle histopathology and age progression is consistently linear 
# across all early age groups, or if distinct intervals exist that warrant separate analytical attention.
#
# Specifically:
# 1. Multiple regression analyses initially demonstrated substantial predictive accuracy 
#    for patient age across the entire sample range. However, subsequent residual (prediction error) 
#    analyses revealed considerable heterogeneity in prediction accuracy across age intervals.
#
# 2. Notably, when we restricted analyses to patients younger than 7 years, prediction errors
#    appeared reduced, suggesting a more homogeneous and predictable histological progression
#    pattern before this age. In contrast, including patients older than 7 years introduced greater
#    variability, reflecting more heterogeneous histological trajectories beyond this age.
#
# 3. Residual analyses across multiple age-group configurations suggested age-dependent
#    differences in predictive accuracy, including around the 6-7-year interval.
#    However, the direct two-group comparison between 1-6 and 6-7 years was not statistically significant.
#    These findings therefore support cautious interpretation, while remaining compatible with the possibility
#    that these intervals represent distinct pathological phases with different dynamics of muscle deterioration.
#
# Thus, these preliminary statistical explorations provided quantitative support
# for our final analytical decision to stratify patient cohorts into distinct 
# intervals (1-6, 6-7, and 7-11 years). The age-specific residual patterns uncovered 
# by these analyses are consistent with the biological hypothesis that muscle histopathological
# progression in DMD is not a continuous, linear process across early childhood, but may instead
# involve transitions or inflection points, particularly around the 6-7-year period.
#
# Consequently, these preliminary modeling and error analyses not only supported
# our chosen age-group stratifications but also reinforced the biological plausibility
# and statistical rationale for examining early-stage DMD progression in distinct
# developmental intervals rather than as a single homogeneous continuum.


###------------------------------### 6. Analysis of Variance-----------------------
# PART 1: Age grouping----------------------------------------------------------

# Divide ABx into 3 categories and create the AgeGroups factor
AgeGroups <- cut(ABx,                         # Allocation target
                 breaks = c(1, 6, 7, 11),    # Delimiter
                 labels = c("[1,6)", "[6,7)", "[7,11)"), 
                 right  = FALSE)             # Exclude rightmost point

# change age group (details are described later)
AgeGroups_rev <- cut(ABx,
                     breaks = c(1, 5, 7.5, 11),
                     labels = c("[1,5)", "[5,7.5)", "[7.5,11)"),
                     right = FALSE)


# PART 2: Mean------------------------------------------------------------------
# Data frame for analysis
dataMean <- data.frame(ABx = ABx, Mean = Mean)

{
  # ── Boxplot ───────────────────────────────
  boxplot(Mean ~ cut(ABx, breaks = c(1, 6, 7, 11), 
                     labels = c("[1-6)", "[6-7)", "[7-11)"), 
                     right = FALSE), 
          xlab = "",            # Leave axis label blank (described in the figure caption)
          cex.axis = 2.5,       # Tick-label font scaling
          ylab = "", 
          cex.lab = 2.5,
          main = "",            # No title
          data = dataMean, 
          outline = FALSE)      # Suppress plotting of outliers
  
  # ── Overlay data points with jitter ──────────────
  stripchart(Mean ~ cut(ABx, breaks = c(1, 6, 7, 11), 
                        labels = c("[1-6)", "[6-7)", "[7-11)"), 
                        right = FALSE), 
             data      = dataMean, 
             method    = "jitter",   # Scatter points with jitter to avoid overlap
             pch       = 16,         # Plot symbol: filled circle
             col       = "black",    # Point color
             vertical  = TRUE,       # Arrange along the y-axis
             add       = TRUE)       # Add to the existing boxplot
}

# ── Extract vectors by age group ─────────────
Mean_Group1 <- dataMean %>% filter(ABx < 6)             %>% pull(Mean)
Mean_Group2 <- dataMean %>% filter(ABx >= 6 & ABx < 7)   %>% pull(Mean)
Mean_Group3 <- dataMean %>% filter(ABx >= 7 & ABx < 11)  %>% pull(Mean)

# ── Summary statistics (mean and standard deviation) ───────────────────────
mean(Mean_Group1)  # 714.8444
sd(Mean_Group1)    # 186.4232

mean(Mean_Group2)  # 1271.2773
sd(Mean_Group2)    # 422.0449

mean(Mean_Group3)  # 1213.5761
sd(Mean_Group3)    # 749.4944

# ── Normality test (Shapiro-Wilk) ────────────────
shapiro.test(Mean_Group1)  # p = 0.2669
shapiro.test(Mean_Group2)  # p = 0.1191
shapiro.test(Mean_Group3)  # p = 9.87e-05

# ── Homogeneity of variance test (Levene) ───────────────────
leveneTest(Mean ~ AgeGroups, data = dataMean)  # p = 0.3581

# ── Group comparison: Kruskal-Wallis ────────────────
kruskal.test(Mean ~ AgeGroups)  # p = 0.002888

# ── Post hoc test: Conover + Bonferroni ─────────
conover_result_Mean <- kwAllPairsConoverTest(Mean ~ AgeGroups,
                                             data = dataMean,
                                             p.adjust.method = "bonferroni")
print(conover_result_Mean)

# ── Review boxplot summary statistics by group ──────────
stats_Mean <- tapply(Mean, AgeGroups, summary)
print(stats_Mean)


# PART 3: Cov-------------------------------------------------------------------
# Extract ABx and Cov only to create a data frame
dataCov <- data.frame(ABx = ABx, Cov = Cov)

{
  # Draw boxplots of Cov across 3 age intervals
  boxplot(Cov ~ cut(ABx, breaks = c(1, 6, 7, 11), 
                    labels = c("[1-6)", "[6-7)", "[7-11)"), 
                    right = FALSE), 
          xlab = "",            # Leave the axis label blank
          cex.axis = 2.5,       # Tick-label font scaling
          ylab = "",
          cex.lab = 2.5,
          main = "",
          data = dataCov, 
          outline = FALSE)      # Hide outliers
  
  # Overlay data points with jitter
  stripchart(Cov ~ cut(ABx, breaks = c(1, 6, 7, 11), 
                       labels = c("[1-6)", "[6-7)", "[7-11)"), 
                       right = FALSE), 
             data      = dataCov, 
             method    = "jitter", 
             pch       = 16,     # Filled circle
             col       = "black",
             vertical  = TRUE, 
             add       = TRUE)         
}

# Extract Cov vectors by age group
Cov_Group1 <- dataCov %>% filter(ABx < 6)            %>% pull(Cov)
Cov_Group2 <- dataCov %>% filter(ABx >= 6 & ABx < 7)  %>% pull(Cov)
Cov_Group3 <- dataCov %>% filter(ABx >= 7 & ABx < 11) %>% pull(Cov)

# Check the mean and standard deviation of each group
mean(Cov_Group1)  # 0.5804623
sd(Cov_Group1)    # 0.1934304
mean(Cov_Group2)  # 0.7849583
sd(Cov_Group2)    # 0.2088378
mean(Cov_Group3)  # 0.9155027
sd(Cov_Group3)    # 0.2822029

# Normality test (Shapiro-Wilk)
shapiro.test(Cov_Group1)  # p = 2.876e-05* 
shapiro.test(Cov_Group2)  # p = 0.9998
shapiro.test(Cov_Group3)  # p = 0.004765*

# Homogeneity of variance test (Levene)
leveneTest(Cov ~ AgeGroups, data = dataCov)  # p = 0.4404

# Group comparison: Kruskal-Wallis
kruskal.test(Cov ~ AgeGroups)  # p = 0.0006608

# Post hoc test: Conover + Bonferroni
conover_result_Cov <- kwAllPairsConoverTest(Cov ~ AgeGroups, data = dataCov,
                                            p.adjust.method = "bonferroni")
print(conover_result_Cov)

# Output boxplot summary statistics
stats_Cov <- tapply(Cov, AgeGroups, summary)
print(stats_Cov)  # 1st Qu. 0.50, 3rd Qu. 0.56 (1-6yo group)
# 1st Qu. 0.50, 3rd Qu. 0.56 (1-6yo group)




# PART 4: Sd--------------------------------------------------------------------
# Data frame for analysis
dataSd <- data.frame(ABx = ABx, Sd = Sd)   # Combine age and standard deviation

{
  # ── Boxplot ─────────────────────────────
  boxplot(Sd ~ cut(ABx, breaks = c(1, 6, 7, 11),           # across 3 age intervals
                   labels = c("[1-6)", "[6-7)", "[7-11)"),
                   right = FALSE),
          xlab = "",            # Leave the axis label blank (described in the figure caption)
          cex.axis = 2.5,       # Tick-label font scaling
          ylab = "",            # Leave the y-axis label blank
          cex.lab = 2.5,        # Label font scaling
          main = "",            # No title
          data = dataSd,
          outline = FALSE)      # Hide outliers
  
  # ── Overlay data points with jitter ────────────
  stripchart(Sd ~ cut(ABx, breaks = c(1, 6, 7, 11),
                      labels = c("[1-6)", "[6-7)", "[7-11)"),
                      right = FALSE),
             data      = dataSd,
             method    = "jitter",   # Scatter points with jitter to avoid overlap
             pch       = 16,         # Plot symbol: filled circle
             col       = "black",    # Point color
             vertical  = TRUE,       # Arrange along the y-axis
             add       = TRUE)       # Add to the existing boxplot
}

# ── Extract values by age group ───────────────────
Sd_Group1 <- dataSd %>% filter(ABx < 6)             %>% pull(Sd)  # <6yo
Sd_Group2 <- dataSd %>% filter(ABx >= 6 & ABx < 7)   %>% pull(Sd)  # 6–7yo
Sd_Group3 <- dataSd %>% filter(ABx >= 7 & ABx < 11)  %>% pull(Sd)  # 7–11yo

# ── Summary statistics (mean and standard deviation) ──────────────────────
mean(Sd_Group1); sd(Sd_Group1)   # 409.13, 137.99
mean(Sd_Group2); sd(Sd_Group2)   # 974.10, 373.86
mean(Sd_Group3); sd(Sd_Group3)   # 1144.48, 832.47

# ── Normality test (Shapiro-Wilk) ────────────────
shapiro.test(Sd_Group1)  # p = 0.4728  -> normality not rejected
shapiro.test(Sd_Group2)  # p = 0.2271  -> normality not rejected
shapiro.test(Sd_Group3)  # p = 0.00097 -> non-normal

# ── Homogeneity of variance test (Levene) ───────────────────
leveneTest(Sd ~ AgeGroups, data = dataSd)  # p = 0.09449 → Almost equal variance

# ── Group comparison: Kruskal-Wallis ───────────────
kruskal.test(Sd ~ AgeGroups)  # p = 0.0001266 -> significant group difference

# ── Post hoc test: Conover (Bonferroni correction) ───
conover_result_Sd <- kwAllPairsConoverTest(Sd ~ AgeGroups,
                                           data = dataSd,
                                           p.adjust.method = "bonferroni")
print(conover_result_Sd)

# ── Boxplot summary statistics (median, IQR, etc.) ─────────
stats_Sd <- tapply(Sd, AgeGroups, summary)
print(stats_Sd)


# PART 5: MFD-------------------------------------------------------------------
# ── Create data frame ───────────────────────
dataMFD <- data.frame(ABx = ABx, MFD = MFD)   # Combine age and myofibre density

{
  # ── Boxplot ─────────────────────────────
  boxplot(MFD ~ cut(ABx, breaks = c(1, 6, 7, 11),           # Plot across 3 age intervals
                    labels = c("[1-6)", "[6-7)", "[7-11)"),
                    right = FALSE),
          xlab = "",            # Leave the axis label blank
          cex.axis = 2.5,       # Tick-label font scaling
          ylab = "",            # Leave the y-axis label blank
          cex.lab = 2.5,        # Label font scaling
          main = "",            # No title
          data = dataMFD,
          outline = FALSE)      # Hide outliers
  
  # ── Overlay data points with jitter ────────────
  stripchart(MFD ~ cut(ABx, breaks = c(1, 6, 7, 11),
                       labels = c("[1-6)", "[6-7)", "[7-11)"),
                       right = FALSE),
             data      = dataMFD,
             method    = "jitter",   # Scatter points with jitter to avoid overlap
             pch       = 16,         # Filled circle
             col       = "black",    # Point color
             vertical  = TRUE,       # Arrange along the y-axis
             add       = TRUE)       # Add to the existing boxplot
}

# ── Extract values by age group ───────────────────
MFD_Group1 <- dataMFD %>% filter(ABx < 6)            %>% pull(MFD)  # 1–6yo
MFD_Group2 <- dataMFD %>% filter(ABx >= 6 & ABx < 7)  %>% pull(MFD) # 6–7yo
MFD_Group3 <- dataMFD %>% filter(ABx >= 7 & ABx < 11) %>% pull(MFD) # 7–11yo

# ── Summary statistics (mean and standard deviation) ──────────────────────
mean(MFD_Group1); sd(MFD_Group1)   # 894.85, 283.75
mean(MFD_Group2); sd(MFD_Group2)   # 399.26, 119.21
mean(MFD_Group3); sd(MFD_Group3)   # 382.04, 125.70

# ── Normality test (Shapiro-Wilk) ────────────────
shapiro.test(MFD_Group1)  # p = 0.5893
shapiro.test(MFD_Group2)  # p = 0.7833
shapiro.test(MFD_Group3)  # p = 0.7894

# ── Homogeneity of variance test (Levene) ───────────────────
leveneTest(MFD ~ AgeGroups, data = dataMFD)  # p = 0.01167

# ── Group comparison: Welch ANOVA ───────────────────
welch_anova_result_MFD <- oneway.test(MFD ~ AgeGroups, data = dataMFD, var.equal = FALSE)
print(welch_anova_result_MFD)  # p = 4.466e-05

# ── Post hoc test: Games-Howell ──────────────────
gh_result_MFD <- gamesHowellTest(MFD ~ AgeGroups, data = dataMFD)
print(gh_result_MFD)

# ── Boxplot summary statistics (median, IQR, etc.) ─────────
stats_MFD <- tapply(MFD, AgeGroups, summary); print(stats_MFD)

# Moving-Average plot
{
  # Load required packages
  library(ggplot2)
  library(dplyr)
  library(zoo)        # rollmean()
  library(purrr)      # map()
  library(cowplot)
  
  # Subset to ABx < 11 years and sort by age
  dat_sub <- dataMFD %>%
    filter(ABx < 11) %>%
    arrange(ABx)
  
  # Moving-average window sizes to compare
  ks <- 2:6
  
  # Compute moving averages for each window size
  ma_long <- map_dfr(ks, function(k) {
    dat_sub %>%
      mutate(MA = rollmean(MFD, k, align = "left", fill = NA),
             window = factor(paste0("k = ", k), levels = paste0("k = ", ks)))
  })
  
  # Reference moving-average window for LOESS matching
  k_ma <- 2    
  
  # Helper function for moving average
  ma_fun <- function(k) rollmean(dat_sub$MFD, k, align = "left", fill = NA)
  
  # Search LOESS span giving the closest fit to the moving average
  span_grid <- seq(0.20, 0.80, 0.02)
  rmse_vec  <- vapply(span_grid, function(sp) {
    lo <- loess(MFD ~ ABx, dat_sub, span = sp,
                control = loess.control(surface = "direct"))
    pred <- predict(lo, dat_sub$ABx)
    sqrt(mean((pred - ma_fun(k_ma))^2, na.rm = TRUE))
  }, numeric(1))
  
  # Select the span with minimum RMSE
  span_eq <- span_grid[ which.min(rmse_vec) ]  
  cat(sprintf("span_eq (k = %d) = %.2f\n", k_ma, span_eq))
  
  # Fit LOESS using the selected span
  lo <- loess(MFD ~ ABx, dat_sub, span = span_eq,
              control = loess.control(surface = "direct"))
  
  # Prediction grid for smooth curve and local slope
  grid <- data.frame(
    ABx = seq(min(dat_sub$ABx), max(dat_sub$ABx), length.out = 200)   
  )
  grid$fit   <- predict(lo, newdata = grid)
  
  # First derivative approximation of the fitted curve
  grid$slope <- c(NA, diff(grid$fit) / diff(grid$ABx))
  
  # Smooth the local slope curve
  grid$slope_smooth <- predict(
    loess(slope ~ ABx, data = grid, span = 0.15, na.action = na.exclude)
  )
  
  # Identify the age with the steepest decline
  slope_vec <- grid$slope_smooth
  slope_vec[is.na(slope_vec)] <- Inf
  bp_loess  <- grid$ABx[ which.min(slope_vec) ]
  bp_lab    <- sprintf("Candidate ≈ %.2f", bp_loess)
  
  # Color palette for moving-average windows
  pal_col <- c("k = 2" = "#d73027",
               "k = 3" = "#1a9850",
               "k = 4" = "#4575b4",
               "k = 5" = "#e6ac00",
               "k = 6" = "#56B4E9")
  
  # Line-type palette for moving-average windows
  pal_lty <- c("k = 2" = "longdash",
               "k = 3" = "dotdash",
               "k = 4" = "twodash",
               "k = 5" = "dotted",
               "k = 6" = "solid")
  
  # Top panel: raw data, moving averages, and candidate breakpoint
  p_top <- ggplot(dat_sub, aes(ABx, MFD)) +
    geom_point(size = 2) +
    geom_line(data = ma_long,
              aes(y = MA, colour = window, linetype = window),
              linewidth = 0.9) +
    geom_vline(xintercept = bp_loess, linetype = "dashed") +
    annotate("text", x = bp_loess, y = max(dat_sub$MFD),
             label = bp_lab, hjust = -0.1, vjust = 0.5, size = 4) +
    scale_colour_manual(values = pal_col, name = "Window\n(size)") +
    scale_linetype_manual(values = pal_lty,  name = "Window\n(size)") +
    guides(colour = guide_legend(
      override.aes = list(linetype = pal_lty, linewidth = 1))) +
    scale_x_continuous(breaks = seq(floor(min(dat_sub$ABx)),
                                    ceiling(max(dat_sub$ABx)), by = 1)) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
    labs(title = "Moving-average comparison (windows 2–6)",
         subtitle = sprintf("Dashed line = steepest LOESS decline (span = %.2f)", span_eq),
         x = "ABx (years)", y = "MFD") +
    theme_minimal(base_size = 14) +
    theme(legend.position = "right")
  
  # Bottom panel: smoothed local slope of the LOESS curve
  p_bottom <- ggplot(grid, aes(ABx, slope_smooth)) +
    geom_area(fill = "steelblue", alpha = 0.8) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = bp_loess, linetype = "dashed") +
    coord_cartesian(ylim = c(-400, 400)) +        
    labs(x = "ABx (years)", y = "Local slope (ΔMFD / Δyear)") +
    theme_minimal(base_size = 14)
  
  # Replot bottom panel with a wider y-axis range
  p_bottom_fix <- p_bottom + coord_cartesian(ylim = c(-1000, 1000))
  
  # Combine the two panels vertically
  cowplot::plot_grid(p_top, p_bottom_fix, ncol = 1, align = "v",
                     rel_heights = c(3, 1.15))
}


# ------- Detailed age group---------------------------------  
# # Reclassify age into 5 intervals
AgeGroups_detailed <- cut(ABx, breaks = c(1, 4, 6, 7, 9, 11), 
                          labels = c("[1,4)", "[4,6)", "[6,7)", "[7,9)", "[9,11)"), right = FALSE)  # Create factor

{  # Visualization block
  boxplot(MFD ~ cut(ABx, breaks = c(1, 4, 6, 7, 9, 11),
                    labels = c("[1-4)", "[4-6)", "[6-7)", "[7-9)", "[9-11)"), right = FALSE), 
          xlab = "", cex.axis = 2, ylab = "", cex.lab = 2.5, main = "", data = dataMFD, outline = FALSE)  # Boxplot
  
  stripchart(MFD ~ cut(ABx, breaks = c(1, 4, 6, 7, 9, 11),
                       labels = c("[1-4)", "[4-6)", "[6-7)", "[7-9)", "[9-11)"), right = FALSE), 
             data = dataMFD, method = "jitter", pch = 16, col = "black", vertical = TRUE, add = TRUE)  # Jittered points
}

{  # Extract groups
  MFD_Group4 <- dataMFD %>% filter(ABx < 4)            %>% pull(MFD)  # 1–4yo
  MFD_Group5 <- dataMFD %>% filter(ABx >= 4 & ABx < 6)  %>% pull(MFD)  # 4–6yo
  MFD_Group6 <- dataMFD %>% filter(ABx >= 6 & ABx < 7)  %>% pull(MFD)  # 6–7yo
  MFD_Group7 <- dataMFD %>% filter(ABx >= 7 & ABx < 9)  %>% pull(MFD)  # 7–9yo
  MFD_Group8 <- dataMFD %>% filter(ABx >= 9 & ABx < 11) %>% pull(MFD)  # 9–11yo
}

# Normality checks
shapiro.test(MFD_Group4)  # p = 0.977
shapiro.test(MFD_Group5)  # p = 0.2791
shapiro.test(MFD_Group6)  # p = 0.7833
shapiro.test(MFD_Group7)  # p = 0.5523
shapiro.test(MFD_Group8)  # p = 0.3354

# Homogeneity of variance
leveneTest(MFD ~ AgeGroups_detailed, data = dataMFD)  # p = 0.2488

anova_result_MFD_detailed <- aov(MFD ~ AgeGroups_detailed, data = dataMFD)  # One-way ANOVA
summary(anova_result_MFD_detailed)  # p = 6.59e-09

tukey_result_MFD_detailed <- TukeyHSD(anova_result_MFD_detailed)  # Tukey HSD
print(tukey_result_MFD_detailed)

stats_MFD_detail <- tapply(MFD, AgeGroups_detailed, summary); stats_MFD_detail  # Summary statistics
stats_MFD_detail_rounded <- lapply(stats_MFD_detail, function(x) format(round(x, 3), nsmall = 3))
print(stats_MFD_detail_rounded)


# ------- two-group comparison-----------
# Grouping
# Purpose: compare MFD (myofibre density) between 2 groups:
# ABx < 6 vs. 6 <= ABx < 11.

# Extract MFD values from the group with ABx < 6
MFD_Group1 <- dataMFD %>% filter(ABx < 6) %>% pull(MFD)

# Extract MFD values from the group with 6 <= ABx < 11
MFD_Group9 <- dataMFD %>% filter(ABx >= 6 & ABx < 11) %>% pull(MFD)

# Boxplot
# Purpose: visually compare the MFD distributions between the 2 groups.
boxplot(MFD ~ cut(ABx, breaks = c(1, 6, 11), 
                  labels = c("[1-6)", "[6-11)"), 
                  right = FALSE))

# Descriptive statistics
# Compute and display the mean and standard deviation for each group.
mean(MFD_Group1) # 894.8523
sd(MFD_Group1) # 283.7481

mean(MFD_Group9) # 388.3022
sd(MFD_Group9) # 120.7904

# Normality test (Shapiro-Wilk)
shapiro.test(MFD_Group1) # p-value = 0.5893
shapiro.test(MFD_Group9) # p-value = 0.6224

# t-test
t.test(MFD_Group1, MFD_Group9) # p-value = 2.209e-05

# Glass's Delta, confidence interval, power, and required sample size
# This block defines and applies a function for effect-size analysis.
{
  # Recreate grouped data frames while keeping all variables
  group1 <- dataMFD %>% filter(ABx < 6)
  group9 <- dataMFD %>% filter(ABx >= 6 & ABx < 11)
  
  # Variables to analyse (MFD only here)
  variables_glass <- c("MFD")
  
  # Function to compute Glass's Delta, 95% CI, power, and required sample size
  calculate_glass_delta <- function(var, ref_group, comp_group) {
    # Remove missing values and extract numeric vectors
    # na.omit() excludes missing values, and as.numeric() ensures numeric input.
    ref_var <- as.numeric(na.omit(ref_group[[var]]))  # Reference group
    comp_var <- as.numeric(na.omit(comp_group[[var]])) # Comparison group
    
    # Compute Glass's Delta (effect size)
    effect_size <- (mean(comp_var) - mean(ref_var)) / sd(ref_var)
    
    # Compute the 95% confidence interval using ci.smd() from MBESS
    ci_lower_upper <- ci.smd(
      smd = effect_size,
      n.1 = length(comp_var),
      n.2 = length(ref_var),
      conf.level = 0.95
    )
    
    # Estimate the required sample size per group for 80% power
    required_n <- pwr.t.test(
      d = effect_size,
      power = 0.8,
      sig.level = 0.05,
      type = "two.sample"
    )$n
    
    # Compute observed power only when both groups have at least 2 observations
    if (length(comp_var) > 1 && length(ref_var) > 1) {
      observed_power <- pwr.t2n.test(
        n1 = length(ref_var),
        n2 = length(comp_var),
        d = effect_size,
        sig.level = 0.05,
        alternative = "two.sided"
      )$power
    } else {
      observed_power <- NA
    }
    
    # Return results as a list
    list(
      effect_size = effect_size,
      ci = ci_lower_upper,
      required_n = required_n,
      observed_power = observed_power
    )
  }
  
  # Apply the function to each variable listed in variables_glass
  for (var in variables_glass) {
    # Display the variable name being analysed
    cat("\nAnalyzing", var, "(Glass's Delta)\n")
    
    # Run the function and store the result
    result <- calculate_glass_delta(var, group1, group9)
    
    # Print the results in a readable format
    cat("1-6yo vs 6-11yo Glass's Delta:", result$effect_size, "\n")
    cat("95% CI:", result$ci$Lower.Conf.Limit.smd, "-", result$ci$Upper.Conf.Limit.smd, "\n")
    cat("Observed power:", result$observed_power, "\n")
    cat("Required sample size (per group, for 80% power):", result$required_n, "\n\n")
  }
}

# Analyzing MFD (Glass's Delta)
# 1-6yo vs 6-11yo Glass's Delta: -1.78521 
# 95% CI: -2.584194 to -0.9666097 
# Observed power: 0.9986048 
# Required sample size (per group, for 80% power): 6.054806


# PART 6: MFA-------------------------------------------------------------------
# 【データフレームの作成】
# ABx（年齢）とMFA（筋線維面積）のデータを含むデータフレームを作成。
dataMFA <- data.frame(ABx = ABx, MFA = MFA)

# 【ボックスプロット描画】
# 年齢別のMFAの分布を視覚的に比較。
{
  boxplot(MFA ~ cut(ABx, breaks = c(1, 6, 7, 11), 
                    labels = c("[1-6)", "[6-7)", "[7-11)"), 
                    right = FALSE), 
          xlab = "",              # x軸ラベル（空白）
          cex.axis = 2.5,         # 軸の文字サイズを2.5倍に拡大
          ylab = "",              # y軸ラベル（空白）
          cex.lab = 2.5,          # ラベルの文字サイズを2.5倍に拡大
          main = "",              # グラフタイトル（空白）
          data = dataMFA,
          outline = FALSE)        # 外れ値を非表示にする
  
  # 【データポイントの追加】
  # 各データポイントをボックスプロット上に散布図として表示。
  stripchart(MFA ~ cut(ABx, breaks = c(1, 6, 7, 11), 
                       labels = c("[1-6)", "[6-7)", "[7-11)"), 
                       right = FALSE), 
             data = dataMFA, 
             method = "jitter",  # データ点が重ならないようにランダムにずらす
             pch = 16,           # 点の形状（塗りつぶした丸）
             col = "black",      # 点の色を黒に設定
             vertical = TRUE,    # 垂直方向に散布
             add = TRUE)         # 現在のグラフに追加表示
}

# 【グループ分け】
# 年齢（ABx）に基づき3つのグループに分け、それぞれのMFA値を抽出。
{
  MFA_Group1 <- dataMFA %>%
    filter(ABx < 6) %>%
    pull(MFA)
  
  MFA_Group2 <- dataMFA %>%
    filter(ABx >= 6 & ABx < 7) %>%
    pull(MFA)
  
  MFA_Group3 <- dataMFA %>%
    filter(ABx >= 7 & ABx < 11) %>%
    pull(MFA)
}

# 【各グループの平均と標準偏差】
mean(MFA_Group1) # 61.33223
sd(MFA_Group1)   # 8.209139

mean(MFA_Group2) # 48.13858
sd(MFA_Group2)   # 13.86261

mean(MFA_Group3) # 42.81239
sd(MFA_Group3)   # 14.54458

# 【正規性の検定（Shapiro-Wilk検定）】
# 各グループのデータが正規分布に従うか検定。
shapiro.test(MFA_Group1) # p-value = 0.5431
shapiro.test(MFA_Group2) # p-value = 0.5803
shapiro.test(MFA_Group3) # p-value = 0.7888

# 【等分散性の検定（Levene検定）】
# グループ間で分散が等しいか検定。
# ここで使用する「AgeGroups」は、ABxを前述の区間で分割した新たな因子変数。
leveneTest(MFA ~ AgeGroups, data=dataMFA) # Pr(>F) = 0.07045

# 【分散分析（ANOVA）】
# グループ間のMFAの平均値の差が統計的に有意か検定。
anova_result_MFA <- aov(MFA ~ AgeGroups, data=dataMFA)
summary(anova_result_MFA) # Pr(>F) = 0.00174  

# 【多重比較（TukeyのHSD検定）】
# ANOVAで有意差があった場合、どのグループ間に差があるかを特定。
tukey_result_MFA <- TukeyHSD(anova_result_MFA)
print(tukey_result_MFA)


# ------- two-group comparison-----

# 【グループ分け】
# ABx（年齢）が6未満と6以上11未満の2グループに分け、それぞれMFA値を抽出。
{
  MFA_Group4 <- dataMFA %>%
    filter(ABx < 6) %>%
    pull(MFA)
  
  MFA_Group5 <- dataMFA %>%
    filter(ABx >= 6 & ABx < 11) %>%
    pull(MFA)
}

# 【ボックスプロットの描画】
# 2グループ間のMFA値の分布を視覚的に比較。
boxplot(MFA ~ cut(ABx, breaks = c(1, 6, 11), 
                  labels = c("[1-6)", "[6-11)"), 
                  right = FALSE))

# 【各グループの平均と標準偏差（記述統計）】
mean(MFA_Group4) # 61.33223
sd(MFA_Group4)   # 8.209139

mean(MFA_Group5) # 44.74919
sd(MFA_Group5)   # 14.20883

# 【正規性の検定（Shapiro-Wilk検定）】
# 各グループが正規分布に従うかを検定。
shapiro.test(MFA_Group4) # p-value = 0.5431
shapiro.test(MFA_Group5) # p-value = 0.528

# 【t検定】
# 2グループ間の平均値の差が統計的に有意かを検定。
t.test(MFA_Group4, MFA_Group5) # p-value = 0.0001144

# 【Glass's Deltaによる効果量、信頼区間、検出力の計算】
{
  # 改めてグループ分けをデータフレームの形で行う
  group4 <- dataMFA %>% filter(ABx < 6)
  group5 <- dataMFA %>% filter(ABx >= 6 & ABx < 11)
  
  # 分析対象となる変数をリストとして指定（ここでは「MFA」のみ）
  variables_glass <- c("MFA")
  
  # Glass's Delta、95%信頼区間、検出力、必要サンプルサイズを計算するための関数
  calculate_glass_delta <- function(var, ref_group, comp_group) {
    # 【欠損値の処理】
    # na.omit()で欠損値(NA)を除外し、as.numeric()で数値型ベクトルに変換。
    ref_var <- as.numeric(na.omit(ref_group[[var]]))  # 基準グループ（年齢が若いグループ）
    comp_var <- as.numeric(na.omit(comp_group[[var]])) # 比較グループ（年齢が高いグループ）
    
    # 【効果量（Glass's Delta）の計算】
    # Glass's Deltaは比較群の平均から基準群の平均を引き、基準群の標準偏差で割った値。
    effect_size <- (mean(comp_var) - mean(ref_var)) / sd(ref_var)
    
    # 【95%信頼区間の計算】
    # MBESSパッケージのci.smd()関数を用いて効果量の95%信頼区間を計算。
    ci_lower_upper <- ci.smd(
      smd = effect_size,
      n.1 = length(comp_var),
      n.2 = length(ref_var),
      conf.level = 0.95
    )
    
    # 【必要サンプルサイズの計算】
    # pwrパッケージのpwr.t.test()関数を使用し、検出力80%を達成するために必要なサンプルサイズを算出。
    required_n <- pwr.t.test(
      d = effect_size,
      power = 0.8,
      sig.level = 0.05,
      type = "two.sample"
    )$n
    
    # 【条件分岐による検出力の計算】
    # 各グループのサンプル数が2以上の場合のみ検出力を計算可能。
    # それ以外の場合は計算できないため、NAを返す。
    if (length(comp_var) > 1 && length(ref_var) > 1) {
      observed_power <- pwr.t2n.test(
        n1 = length(ref_var),
        n2 = length(comp_var),
        d = effect_size,
        sig.level = 0.05,
        alternative = "two.sided"
      )$power
    } else {
      observed_power <- NA
    }
    
    # 結果をリスト形式で返す。
    list(
      effect_size = effect_size,
      ci = ci_lower_upper,
      required_n = required_n,
      observed_power = observed_power
    )
  }
  
  # 【forループによる計算の繰り返し】
  # 変数リスト（variables_glass）の各変数に対してGlass's Delta等を計算。
  for (var in variables_glass) {
    cat("\nAnalyzing", var, "(Glass's Delta)\n")
    
    result <- calculate_glass_delta(var, group4, group5)
    
    # 結果の表示
    cat("1-6yo vs 6-11yo Glass's Delta:", result$effect_size, "\n")
    cat("95% CI:", result$ci$Lower.Conf.Limit.smd, "-", result$ci$Upper.Conf.Limit.smd, "\n")
    cat("Observed power:", result$observed_power, "\n")
    cat("Required sample size (per group, for 80% power):", result$required_n, "\n\n")
  }
}


# PART 7: CFA-------------------------------------------------------------------
# 【データフレームの作成】
# ABx（年齢）とCFA（結合組織面積）のデータを含むデータフレームを作成。
dataCFA <- data.frame(ABx = ABx, CFA = CFA)

# 【ボックスプロットの描画】
# 年齢別のCFA値の分布を視覚的に比較。
{
  boxplot(CFA ~ cut(ABx, breaks = c(1, 6, 7, 11), 
                    labels = c("[1-6)", "[6-7)", "[7-11)"), 
                    right = FALSE), 
          xlab = "",                # x軸ラベル（空白）
          cex.axis = 2.5,           # 軸目盛りの文字サイズを2.5倍に拡大
          ylab = "",                # y軸ラベル（空白）
          cex.lab = 2.5,            # ラベルの文字サイズを2.5倍に拡大
          main = "",                # グラフタイトル（空白）
          data = dataCFA,
          outline = FALSE)          # 外れ値を非表示
  
  # 【データポイントの追加】
  # 各データポイントをボックスプロット上に散布図として表示。
  stripchart(CFA ~ cut(ABx, breaks = c(1, 6, 7, 11), 
                       labels = c("[1-6)", "[6-7)", "[7-11)"), 
                       right = FALSE), 
             data = dataCFA, 
             method = "jitter",    # 点が重ならないようにランダムにずらす
             pch = 16,             # 点の形状を塗りつぶした丸に指定
             col = "black",        # 点の色を黒に指定
             vertical = TRUE,      # 縦方向に表示
             add = TRUE)           # 既存のグラフに追加
}


# ------- two-group comparison--------------------------------------------------

# 【グループ分け】
# ABx（年齢）を6未満と6以上11未満の2グループに分け、各グループのCFA値を抽出。
{
  CFA_Group4 <- dataCFA %>%
    filter(ABx < 6) %>%
    pull(CFA)
  
  CFA_Group5 <- dataCFA %>%
    filter(ABx >= 6 & ABx < 11) %>%
    pull(CFA)
}

# 【ボックスプロットの描画】
# 2グループ間のCFA値の分布を視覚的に比較。
boxplot(CFA ~ cut(ABx, breaks = c(1, 6, 11), 
                  labels = c("[1-6)", "[6-11)"), 
                  right = FALSE))

# 【各グループの平均値と標準偏差】
mean(CFA_Group4) # 36.03298
sd(CFA_Group4)   # 7.088319

mean(CFA_Group5) # 48.77115
sd(CFA_Group5)   # 13.2095

# 【正規性の検定（Shapiro-Wilk検定）】
# 各グループが正規分布に従うかを検定。
shapiro.test(CFA_Group4) # p-value = 0.9153
shapiro.test(CFA_Group5) # p-value = 0.3401

# 【t検定】
# 2グループ間の平均値の差が統計的に有意かを検定。
t.test(CFA_Group4, CFA_Group5) # p-value = 0.0007671

# 【Glass's Deltaによる効果量、信頼区間、検出力の計算】
{
  # 再度グループ分けをデータフレーム形式で行う
  group4 <- dataCFA %>% filter(ABx < 6)
  group5 <- dataCFA %>% filter(ABx >= 6 & ABx < 11)
  
  # 分析する変数のリストを指定（ここではCFAのみ）
  variables_glass <- c("CFA")
  
  # Glass's Delta、95%信頼区間、検出力、必要サンプルサイズを計算する関数を作成
  calculate_glass_delta <- function(var, ref_group, comp_group) {
    # 【欠損値の処理】
    # na.omit()で欠損値(NA)を除外し、as.numeric()で数値型ベクトルに変換。
    ref_var <- as.numeric(na.omit(ref_group[[var]]))  # 基準グループ（若年グループ）
    comp_var <- as.numeric(na.omit(comp_group[[var]])) # 比較グループ（年齢が高いグループ）
    
    # 【効果量（Glass's Delta）の計算】
    # Glass's Deltaは、比較群の平均から基準群の平均を引き、基準群の標準偏差で割った値。
    effect_size <- (mean(comp_var) - mean(ref_var)) / sd(ref_var)
    
    # 【95%信頼区間の計算】
    # MBESSパッケージのci.smd関数で効果量の95%信頼区間を計算。
    ci_lower_upper <- ci.smd(
      smd = effect_size,
      n.1 = length(comp_var),
      n.2 = length(ref_var),
      conf.level = 0.95
    )
    
    # 【必要サンプルサイズの計算】
    # pwrパッケージのpwr.t.test関数で、検出力80%を達成するために必要なサンプルサイズを計算。
    required_n <- pwr.t.test(
      d = effect_size,
      power = 0.8,
      sig.level = 0.05,
      type = "two.sample"
    )$n
    
    # 【条件分岐による検出力の計算】
    # 各グループのサンプル数が2以上の場合のみ検出力を計算し、それ以外はNAを返す。
    if (length(comp_var) > 1 && length(ref_var) > 1) {
      observed_power <- pwr.t2n.test(
        n1 = length(ref_var),
        n2 = length(comp_var),
        d = effect_size,
        sig.level = 0.05,
        alternative = "two.sided"
      )$power
    } else {
      observed_power <- NA
    }
    
    # 結果をリスト形式で返す。
    list(
      effect_size = effect_size,
      ci = ci_lower_upper,
      required_n = required_n,
      observed_power = observed_power
    )
  }
  
  # 【forループによる計算の繰り返し】
  # 指定された各変数に対してGlass's Deltaなどの統計量を計算。
  for (var in variables_glass) {
    cat("\nAnalyzing", var, "(Glass's Delta)\n")
    
    result <- calculate_glass_delta(var, group4, group5)
    
    # 結果の詳細な表示
    cat("1-6yo vs 6-11yo Glass's Delta:", result$effect_size, "\n")
    cat("95% CI:", result$ci$Lower.Conf.Limit.smd, "-", result$ci$Upper.Conf.Limit.smd, "\n")
    cat("Observed power:", result$observed_power, "\n")
    cat("Required sample size (per group, for 80% power):", result$required_n, "\n\n")
  }
}


# PART 8: Fat-------------------------------------------------------------------
# 【データフレーム作成】
# ABx（年齢）とFat（脂肪組織面積）のデータをまとめたデータフレームを作成。
dataFat <- data.frame(ABx = ABx, Fat = Fat)

# 【グループ分け】
# 年齢に基づき3つのグループに分け、それぞれFat値を抽出。
{
  Fat_Group1 <- dataFat %>%
    filter(ABx < 6) %>%
    pull(Fat)
  
  Fat_Group2 <- dataFat %>%
    filter(ABx >= 6 & ABx < 7) %>%
    pull(Fat)
  
  Fat_Group3 <- dataFat %>%
    filter(ABx >= 7 & ABx < 11) %>%
    pull(Fat)
}

# 【各グループの平均値と標準偏差】
mean(Fat_Group1) # 1.112968
sd(Fat_Group1)   # 1.416389

mean(Fat_Group2) # 1.507363
sd(Fat_Group2)   # 1.209222

mean(Fat_Group3) # 6.53602
sd(Fat_Group3)   # 4.829781

# 【正規性の検定（Shapiro-Wilk検定）】
# 各グループが正規分布に従うかを検定。
shapiro.test(Fat_Group1) # p-value = 0.0006062（正規性なし）
shapiro.test(Fat_Group2) # p-value = 0.2932
shapiro.test(Fat_Group3) # p-value = 0.06398

# 【等分散性の検定（Levene検定）】
# グループ間で分散が等しいか検定。
leveneTest(Fat ~ AgeGroups, data=dataFat) # Pr(>F) = 0.01718（等分散性なし）

# 【Kruskal-Wallis検定】
# 正規性・等分散性のいずれかが満たされないため、ノンパラメトリック検定を使用。
kruskal.test(Fat ~ AgeGroups) # p-value = 0.0003663（有意差あり）

# 【Conover多重比較検定（Bonferroni補正）】
# Kruskal-Wallis検定で有意差が認められた場合、具体的なグループ間差を特定。
conover_result_Fat <- kwAllPairsConoverTest(Fat ~ AgeGroups, data = dataFat, p.adjust.method = "bonferroni")
print(conover_result_Fat)

# 【各グループの統計量】
# 各グループごとの記述統計量（最小値、四分位点、中央値、平均、最大値など）を表示。
stats_Fat <- tapply(Fat, AgeGroups, summary)
print(stats_Fat)


# PART 9: Conover (Mean,Sd,Cov,Fat)---------------------------------------------
# 【Conover検定とBonferroni補正による多重比較（95％信頼区間付き）】
{
  # 【データの前処理】
  # ABxが11以下のデータのみを抽出し、欠損値を除去
  filtered_dataAll <- dataAll %>% filter(ABx <= 11) %>% na.omit()
  
  # 【年齢グループの定義】
  # ABxを以下の3つの年齢区間に分類（1-6歳未満、6-7歳未満、7-11歳未満）
  filtered_dataAll$AgeGroups <- cut(filtered_dataAll$ABx, breaks = c(1, 6, 7, 11), 
                                    labels = c("1-6yo", "6-7yo", "7-11yo"), right = F)
  
  # 【分析対象変数リスト】
  variables <- c("Mean", "Sd", "Cov", "Fat")
  
  # 【信頼区間を計算する関数の作成】
  calculate_confidence_interval <- function(diff, sd1, sd2, n1, n2) {
    # 2つのグループ間の平均差(diff)の標準誤差(se_diff)を計算
    se_diff <- sqrt((sd1^2 / n1) + (sd2^2 / n2))
    alpha <- 0.05  # 有意水準（95%信頼区間）
    t_crit <- qt(1 - alpha / 2, df = n1 + n2 - 2)  # 自由度は(n1+n2-2)
    ci_lower <- diff - t_crit * se_diff  # 信頼区間の下限
    ci_upper <- diff + t_crit * se_diff  # 信頼区間の上限
    return(c(ci_lower, ci_upper))
  }
  
  # 【Conover検定の実施】
  # 指定された各変数について繰り返し分析
  for (var in variables) {
    cat("\nAnalyzing", var, "\n")
    
    # Conover検定を実行し、Bonferroni補正を用いた多重比較を実施
    conover_result <- kwAllPairsConoverTest(as.formula(paste(var, "~ AgeGroups")), 
                                            data = filtered_dataAll, p.adjust.method = "bonferroni")
    print(conover_result)
    
    # 各グループ間の比較結果とp値を表示
    comparison_pairs <- rownames(conover_result$p.value)
    for (i in seq_along(comparison_pairs)) {
      comparison <- comparison_pairs[i]
      p_value <- conover_result$p.value[i]
      cat("Comparison:", comparison, "\n")
      cat("p-value:", p_value, "\n\n")
    }
    
    # 【各グループごとのデータ抽出】
    group1 <- filtered_dataAll[[var]][filtered_dataAll$AgeGroups == "1-6yo"]
    group2 <- filtered_dataAll[[var]][filtered_dataAll$AgeGroups == "6-7yo"]
    group3 <- filtered_dataAll[[var]][filtered_dataAll$AgeGroups == "7-11yo"]
    
    # 【平均値の差および信頼区間の計算と表示】
    # グループ1 (1-6yo) とグループ2 (6-7yo) 間の比較
    diff_1_2 <- mean(group2, na.rm = TRUE) - mean(group1, na.rm = TRUE)
    sd_1 <- sd(group1, na.rm = TRUE)
    sd_2 <- sd(group2, na.rm = TRUE)
    n_1 <- length(group1)
    n_2 <- length(group2)
    ci_1_2 <- calculate_confidence_interval(diff_1_2, sd_2, sd_1, n_2, n_1)
    cat("1-6yo vs 6-7yo: Mean Difference =", diff_1_2, "Confidence Interval =", ci_1_2[1], "to", ci_1_2[2], "\n")
    
    # グループ1 (1-6yo) とグループ3 (7-11yo) 間の比較
    diff_1_3 <- mean(group3, na.rm = TRUE) - mean(group1, na.rm = TRUE)
    sd_3 <- sd(group3, na.rm = TRUE)
    n_3 <- length(group3)
    ci_1_3 <- calculate_confidence_interval(diff_1_3, sd_3, sd_1, n_3, n_1)
    cat("1-6yo vs 7-11yo: Mean Difference =", diff_1_3, "Confidence Interval =", ci_1_3[1], "to", ci_1_3[2], "\n")
    
    # グループ2 (6-7yo) とグループ3 (7-11yo) 間の比較
    diff_2_3 <- mean(group3, na.rm = TRUE) - mean(group2, na.rm = TRUE)
    ci_2_3 <- calculate_confidence_interval(diff_2_3, sd_3, sd_2, n_3, n_2)
    cat("6-7yo vs 7-11yo: Mean Difference =", diff_2_3, "Confidence Interval =", ci_2_3[1], "to", ci_2_3[2], "\n")
  }
}


# PART 10: Games-Howell (MFD,CFA)-----------------------------------------------
# 【Games-Howell検定（95％信頼区間付き）による多重比較】
{
  # 【データのフィルタリング】
  # ABxが11以下のデータのみを抽出し、欠損値を除去
  filtered_dataAll <- dataAll %>% filter(ABx <= 11) %>% na.omit()
  
  # 【年齢グループの定義】
  # 年齢を1-6歳、6-7歳、7-11歳の3つのグループに分類
  filtered_dataAll$AgeGroups <- cut(filtered_dataAll$ABx, breaks = c(1, 6, 7, 11), 
                                    labels = c("1-6yo", "6-7yo", "7-11yo"), right = FALSE)
  
  # 【分析対象変数リスト】
  variables <- c("MFD", "CFA")
  
  # 【95%信頼区間を計算する関数の作成】
  calculate_confidence_interval <- function(diff, sd1, sd2, n1, n2) {
    # 2群間の平均差(diff)の標準誤差(se_diff)を計算
    se_diff <- sqrt((sd1^2 / n1) + (sd2^2 / n2))
    alpha <- 0.05  # 有意水準（95%信頼区間）
    t_crit <- qt(1 - alpha / 2, df = n1 + n2 - 2)  # 自由度(n1+n2-2)でのt値を取得
    ci_lower <- diff - t_crit * se_diff  # 信頼区間の下限
    ci_upper <- diff + t_crit * se_diff  # 信頼区間の上限
    return(c(ci_lower, ci_upper))
  }
  
  # 【Games-Howell検定の実施】
  for (var in variables) {
    cat("\nAnalyzing", var, "\n")
    # Games-Howell検定を実行
    gh_result <- gamesHowellTest(as.formula(paste(var, "~ AgeGroups")), data = filtered_dataAll)
    print(gh_result)
    
    # グループ別のデータ抽出
    group1 <- filtered_dataAll[[var]][filtered_dataAll$AgeGroups == "1-6yo"]
    group2 <- filtered_dataAll[[var]][filtered_dataAll$AgeGroups == "6-7yo"]
    group3 <- filtered_dataAll[[var]][filtered_dataAll$AgeGroups == "7-11yo"]
    
    # 比較対象の組み合わせリストを作成
    comparisons <- list(
      "1-6yo vs 6-7yo" = list(group1, group2),
      "1-6yo vs 7-11yo" = list(group1, group3),
      "6-7yo vs 7-11yo" = list(group2, group3)
    )
    
    # 各比較ごとに平均差と95％信頼区間を計算
    for (comp in names(comparisons)) {
      g1 <- comparisons[[comp]][[2]]  # 比較群のデータ
      g2 <- comparisons[[comp]][[1]]  # 基準群のデータ
      diff <- mean(g1, na.rm = TRUE) - mean(g2, na.rm = TRUE)  # 平均差の計算
      sd1 <- sd(g1, na.rm = TRUE)
      sd2 <- sd(g2, na.rm = TRUE)
      n1 <- length(g1)
      n2 <- length(g2)
      ci <- calculate_confidence_interval(diff, sd1, sd2, n1, n2)  # 信頼区間を計算
      cat("Comparison:", comp, "\n")
      cat("Mean Difference:", diff, "\n")
      cat("Confidence Interval:", ci[1], "to", ci[2], "\n\n")
    }
  }
}


# PART 11: Tukey HSD (MFA)---------------------------------------------------
# 【TukeyのHSD検定（95％信頼区間付き）による多重比較】
{
  # 【データのフィルタリング】
  # ABx（年齢）が11以下のデータを抽出し、欠損値を除去。
  filtered_dataAll <- dataAll %>% filter(ABx <= 11) %>% na.omit()
  
  # 【年齢グループの定義】
  # 年齢を以下の3つの区間（1-6歳未満、6-7歳未満、7-11歳未満）に分類。
  filtered_dataAll$AgeGroups <- cut(filtered_dataAll$ABx, breaks = c(1, 6, 7, 11), 
                                    labels = c("1-6yo", "6-7yo", "7-11yo"), right = FALSE)
  
  # 【分析対象の変数リスト】
  variables <- c("MFA")
  
  # 【TukeyのHSD検定の実施】
  # TukeyのHSD検定は、ANOVAで有意差が見られた後に具体的なグループ間差を確認する。
  for (var in variables) {
    cat("\nAnalyzing", var, "\n")
    
    # ANOVAモデルを作成（変数を年齢グループで比較）
    model <- aov(as.formula(paste(var, "~ AgeGroups")), data = filtered_dataAll)
    
    # Tukey HSD検定を実行
    tukey_result <- TukeyHSD(model)
    
    # 【結果の整形と表示】
    cat("\nTukey HSD results for", var, ":\n")
    result <- as.data.frame(tukey_result$AgeGroups)
    colnames(result) <- c("Mean Difference", "Lower CI", "Upper CI", "Adjusted p-value")
    
    # 各比較の平均差、信頼区間、および調整済みp値を表示
    print(result)
  }
}


### -----------------------------### 7. Effect Size Analysis with 95% CI-----------
# PART 1: Cohen's d-------------------------------------------------------------
{
  # 【グループ分け】
  # 年齢（ABx）を基準に、以下の3つのグループに分割。
  group1 <- dataAll %>% filter(ABx < 6)               # 1-6歳未満のグループ
  group2 <- dataAll %>% filter(ABx >= 6 & ABx < 7)    # 6-7歳未満のグループ
  group3 <- dataAll %>% filter(ABx >= 7 & ABx < 11)   # 7-11歳未満のグループ
  
  # 【分析対象変数のリスト】
  variables <- c("Mean", "Sd", "Cov", "MFD", "MFA", "CFA", "Fat")
  
  # 【効果量（Cohen's d）および必要サンプルサイズを計算する関数を作成】
  calculate_effect_size <- function(var, group1, group2) {
    # 【欠損値処理】
    # na.omit()を使用して欠損値を除去し、as.numeric()で数値型ベクトルに変換。
    var1 <- as.numeric(na.omit(group1[[var]]))  # 基準グループ（比較元）
    var2 <- as.numeric(na.omit(group2[[var]]))  # 比較対象グループ
    
    # 【Cohen's dの計算】
    # cohen.d()関数で効果量（2群間の平均差を標準偏差で標準化した指標）を算出。
    effect_size <- cohen.d(var2, var1)$estimate
    
    # 【95%信頼区間の計算】
    # ci.smd()関数で効果量の95%信頼区間を計算。
    ci <- ci.smd(smd = effect_size, n.1 = length(var1), n.2 = length(var2), conf.level = 0.95)
    
    # 【必要なサンプルサイズの計算】
    # pwr.t.test()関数で検出力80%を達成するための必要サンプルサイズを算出。
    sample_size <- pwr.t.test(d = effect_size, power = 0.8, sig.level = 0.05, type = "two.sample")$n
    
    # 結果をリスト形式で返す。
    list(effect_size = effect_size, ci = ci, sample_size = sample_size)
  }
  
  # 【各変数に対する計算の繰り返し（forループ）】
  for (var in variables) {
    cat("\nAnalyzing", var, "\n")
    
    # 各グループ間の比較で効果量・信頼区間・必要サンプルサイズを計算。
    result_1_2 <- calculate_effect_size(var, group1, group2)
    result_1_3 <- calculate_effect_size(var, group1, group3)
    result_2_3 <- calculate_effect_size(var, group2, group3)
    
    # 【結果の表示】
    # 1-6歳未満 vs 6-7歳未満の比較
    cat("1-6yo vs 6-7yo Effect Size (Cohen's d):", result_1_2$effect_size, "\n")
    cat("1-6yo vs 6-7yo 95% CI:", result_1_2$ci$Lower.Conf.Limit.smd, "-", result_1_2$ci$Upper.Conf.Limit.smd, "\n")
    cat("Required sample size:", result_1_2$sample_size, "\n\n")
    
    # 1-6歳未満 vs 7-11歳未満の比較
    cat("1-6yo vs 7-11yo Effect Size (Cohen's d):", result_1_3$effect_size, "\n")
    cat("1-6yo vs 7-11yo 95% CI:", result_1_3$ci$Lower.Conf.Limit.smd, "-", result_1_3$ci$Upper.Conf.Limit.smd, "\n")
    cat("Required sample size:", result_1_3$sample_size, "\n\n")
    
    # 6-7歳未満 vs 7-11歳未満の比較
    cat("6-7yo vs 7-11yo Effect Size (Cohen's d):", result_2_3$effect_size, "\n")
    cat("6-7yo vs 7-11yo 95% CI:", result_2_3$ci$Lower.Conf.Limit.smd, "-", result_2_3$ci$Upper.Conf.Limit.smd, "\n")
    cat("Required sample size:", result_2_3$sample_size, "\n")
  }
}


# PART 2: Hedges' d: unbiased d-------------------------------------------------
# 注意事項:
# - Bootstrapped CI（ブートストラップ法を使った信頼区間）についての詳細は以下の論文を参照。
# - Nakagawa S et al. (Biol Rev Camb Philos Soc. 2007)
{
  # 【年齢グループの定義】
  group1 <- dataAll %>% filter(ABx < 6)               # 1-6歳未満のグループ
  group2 <- dataAll %>% filter(ABx >= 6 & ABx < 7)    # 6-7歳未満のグループ
  group3 <- dataAll %>% filter(ABx >= 7 & ABx < 11)   # 7-11歳未満のグループ
  
  # 【効果量(Hedges' d)と信頼区間をブートストラップ法で算出する関数】
  bootstrap_ci_hedges <- function(data1, data2, n_iter = 5000, conf_level = 0.95, seed = NULL) {
    # 乱数生成の再現性を保つためのシード設定（任意）
    if (!is.null(seed)) set.seed(seed)
    
    # Hedges' dの値をブートストラップで繰り返し計算し、ベクトルに格納
    boot_d <- numeric(n_iter)
    for (i in 1:n_iter) {
      sample1 <- sample(data1, length(data1), replace = TRUE)
      sample2 <- sample(data2, length(data2), replace = TRUE)
      boot_d[i] <- cohen.d(sample1, sample2, hedges.correction = TRUE)$estimate
    }
    
    # 元データによるHedges' d（ブートストラップなし）
    original_d <- cohen.d(data1, data2, hedges.correction = TRUE)$estimate
    
    # ブートストラップ分布から95%信頼区間を計算
    lower_bound <- quantile(boot_d, (1 - conf_level) / 2)
    upper_bound <- quantile(boot_d, 1 - (1 - conf_level) / 2)
    
    # 結果をリスト形式で返す
    return(list(
      original_d = original_d,
      boot_mean = mean(boot_d),
      ci_lower = lower_bound,
      ci_upper = upper_bound
    ))
  }
  
  # 【分析対象の変数リスト】
  variables <- c("Mean", "Sd", "Cov", "MFD", "MFA", "CFA", "Fat")
  
  # 【各変数に対する効果量と信頼区間の計算（繰り返し処理）】
  for (var in variables) {
    cat("\nAnalyzing", var, "\n")
    
    # 各グループの欠損値を除去し、数値型ベクトルに変換
    group1_var <- as.numeric(na.omit(group1[[var]]))
    group2_var <- as.numeric(na.omit(group2[[var]]))
    group3_var <- as.numeric(na.omit(group3[[var]]))
    
    # 1-6歳未満 vs 6-7歳未満の比較
    result_1_2 <- bootstrap_ci_hedges(group1_var, group2_var, seed = NULL)
    power_1_2 <- pwr.t2n.test(n1=length(group1_var), n2=length(group2_var),
                              d=result_1_2$original_d, sig.level=0.05)$power * 100
    requiredN_1_2 <- pwr.t.test(d=abs(result_1_2$original_d), power=0.80,
                                sig.level=0.05, type="two.sample")$n
    
    cat("1-6yo vs 6-7yo Hedges' d (Original):", result_1_2$original_d, "\n")
    cat("   Hedges' d (Bootstrap mean):", result_1_2$boot_mean, "\n")
    cat("   95% Bootstrap CI:", result_1_2$ci_lower, "-", result_1_2$ci_upper, "\n")
    cat("   Observed power (%):", round(power_1_2, 2), "\n")
    cat("   Required sample size (for 80% power, each group):", requiredN_1_2, "\n\n")
    
    # 1-6歳未満 vs 7-11歳未満の比較
    result_1_3 <- bootstrap_ci_hedges(group1_var, group3_var, seed = NULL)
    power_1_3 <- pwr.t2n.test(n1=length(group1_var), n2=length(group3_var),
                              d=result_1_3$original_d, sig.level=0.05)$power * 100
    requiredN_1_3 <- pwr.t.test(d=abs(result_1_3$original_d), power=0.80,
                                sig.level=0.05, type="two.sample")$n
    
    cat("1-6yo vs 7-11yo Hedges' d (Original):", result_1_3$original_d, "\n")
    cat("   Hedges' d (Bootstrap mean):", result_1_3$boot_mean, "\n")
    cat("   95% Bootstrap CI:", result_1_3$ci_lower, "-", result_1_3$ci_upper, "\n")
    cat("   Observed power (%):", round(power_1_3, 2), "\n")
    cat("   Required sample size (for 80% power, each group):", requiredN_1_3, "\n\n")
    
    # 6-7歳未満 vs 7-11歳未満の比較
    result_2_3 <- bootstrap_ci_hedges(group2_var, group3_var, seed = NULL)
    power_2_3 <- pwr.t2n.test(n1=length(group2_var), n2=length(group3_var),
                              d=result_2_3$original_d, sig.level=0.05)$power * 100
    requiredN_2_3 <- pwr.t.test(d=abs(result_2_3$original_d), power=0.80,
                                sig.level=0.05, type="two.sample")$n
    
    cat("6-7yo vs 7-11yo Hedges' d (Original):", result_2_3$original_d, "\n")
    cat("   Hedges' d (Bootstrap mean):", result_2_3$boot_mean, "\n")
    cat("   95% Bootstrap CI:", result_2_3$ci_lower, "-", result_2_3$ci_upper, "\n")
    cat("   Observed power (%):", round(power_2_3, 2), "\n")
    cat("   Required sample size (for 80% power, each group):", requiredN_2_3, "\n\n")
  }
}

# Note:
# Hedges' d is strictly appropriate only for MFA (myofiber area), 
# as it meets the assumption of equal variances.
# However, since the alternative method (Glass's delta), described later, 
# tends to underestimate the effect size for MFD due to its larger variance within the 1–6-year group,
# Hedges' d was computed for all parameters to allow consistent comparison across parameters,
# thereby complementing the primary analyses.
# 0.2: Small, 0.5: Medium, 0.8: Large

   
# PART 3: Cliff's Delta---------------------------------------------------------
{
  # 【グループ分け】
  group1 <- dataAll %>% filter(ABx < 6)               # 1-6歳未満
  group2 <- dataAll %>% filter(ABx >= 6 & ABx < 7)    # 6-7歳未満
  group3 <- dataAll %>% filter(ABx >= 7 & ABx < 11)   # 7-11歳未満
  
  # 正規性を満たさない可能性がある変数リスト
  non_normal_vars <- c("Mean", "Sd", "Cov", "Fat")
  
  # 【Cliff's Deltaと95％信頼区間を計算する関数】
  calculate_cliffs_delta <- function(var, groupA, groupB) {
    varA <- as.numeric(na.omit(groupA[[var]]))
    varB <- as.numeric(na.omit(groupB[[var]]))
    cliff_result <- cliff.delta(varA, varB, conf.level = 0.95)
    
    list(
      delta     = cliff_result$estimate,
      ci_lower  = cliff_result$conf.int[1],
      ci_upper  = cliff_result$conf.int[2],
      magnitude = cliff_result$magnitude
    )
  }
  
  # 【Cliff's DeltaをCohen's dに近似変換する関数】
  cliffs_to_d <- function(delta) {
    p <- (delta + 1) / 2
    p <- ifelse(p <= 0, 1e-10, ifelse(p >= 1, 1 - 1e-10, p))
    d_approx <- sqrt(2) * qnorm(p)
    return(d_approx)
  }
  
  # 【各変数に対する効果量と検出力の計算】
  for (var in non_normal_vars) {
    cat("\nAnalyzing", var, "(Non-parametric: Cliff's Delta)\n")
    
    group1_var <- as.numeric(na.omit(group1[[var]]))
    group2_var <- as.numeric(na.omit(group2[[var]]))
    group3_var <- as.numeric(na.omit(group3[[var]]))
    
    # Cliff's Delta計算
    result_1_2 <- calculate_cliffs_delta(var, group2, group1)
    result_1_3 <- calculate_cliffs_delta(var, group3, group1)
    result_2_3 <- calculate_cliffs_delta(var, group3, group2)
    
    # 結果表示
    cat("1-6yo vs 6-7yo Cliff's Delta:", result_1_2$delta, "\n")
    cat("  95% CI:", result_1_2$ci_lower, "-", result_1_2$ci_upper, "\n")
    cat("  Magnitude:", result_1_2$magnitude, "\n\n")
    
    cat("1-6yo vs 7-11yo Cliff's Delta:", result_1_3$delta, "\n")
    cat("  95% CI:", result_1_3$ci_lower, "-", result_1_3$ci_upper, "\n")
    cat("  Magnitude:", result_1_3$magnitude, "\n\n")
    
    cat("6-7yo vs 7-11yo Cliff's Delta:", result_2_3$delta, "\n")
    cat("  95% CI:", result_2_3$ci_lower, "-", result_2_3$ci_upper, "\n")
    cat("  Magnitude:", result_2_3$magnitude, "\n\n")
    
    # Cohen's d近似値計算
    d_1_2_approx <- cliffs_to_d(result_1_2$delta)
    d_1_3_approx <- cliffs_to_d(result_1_3$delta)
    d_2_3_approx <- cliffs_to_d(result_2_3$delta)
    
    # 観測された検出力と必要サンプルサイズの計算
    n1_1_2 <- length(group1_var); n2_1_2 <- length(group2_var)
    obs_power_1_2 <- ifelse(n1_1_2 > 1 && n2_1_2 > 1,
                            pwr.t2n.test(n1 = n1_1_2, n2 = n2_1_2, d = d_1_2_approx, sig.level = 0.05)$power * 100, NA)
    reqN_1_2 <- pwr.t.test(d = abs(d_1_2_approx), power = 0.80, sig.level = 0.05, type = "two.sample")$n
    
    cat("   Observed power (%):", round(obs_power_1_2, 2), "\n")
    cat("   Required sample size (80% power):", round(reqN_1_2, 2), "\n\n")
    
    n1_1_3 <- length(group1_var); n2_1_3 <- length(group3_var)
    obs_power_1_3 <- ifelse(n1_1_3 > 1 && n2_1_3 > 1,
                            pwr.t2n.test(n1 = n1_1_3, n2 = n2_1_3, d = d_1_3_approx, sig.level = 0.05)$power * 100, NA)
    reqN_1_3 <- pwr.t.test(d = abs(d_1_3_approx), power = 0.80, sig.level = 0.05, type = "two.sample")$n
    
    cat("   Observed power (%):", round(obs_power_1_3, 2), "\n")
    cat("   Required sample size (80% power):", round(reqN_1_3, 2), "\n\n")
    
    n1_2_3 <- length(group2_var); n2_2_3 <- length(group3_var)
    obs_power_2_3 <- ifelse(n1_2_3 > 1 && n2_2_3 > 1,
                            pwr.t2n.test(n1 = n1_2_3, n2 = n2_2_3, d = d_2_3_approx, sig.level = 0.05)$power * 100, NA)
    reqN_2_3 <- pwr.t.test(d = abs(d_2_3_approx), power = 0.80, sig.level = 0.05, type = "two.sample")$n
    
    cat("   Observed power (%):", round(obs_power_2_3, 2), "\n")
    cat("   Required sample size (80% power):", round(reqN_2_3, 2), "\n\n")
  }
}

# 【Bootstrap法によるCliff's Deltaと信頼区間の算出】
# Bootstrap function for Cliff's Delta with reproducible seed setting
{
  # 【ブートストラップ法を用いたCliff's Deltaと95%信頼区間を計算する関数】
  bootstrap_cliffs_delta <- function(var, groupA, groupB, n_bootstrap = 5000, conf_level = 0.95, seed = NULL) {
    
    # データを数値型ベクトルとして取得（欠損値を除去）
    data_A <- as.numeric(na.omit(groupA[[var]]))  # 比較グループ
    data_B <- as.numeric(na.omit(groupB[[var]]))  # 基準グループ
    
    # オリジナルのCliff's Deltaを計算
    original_delta <- cliff.delta(data_A, data_B)$estimate
    
    # 再現性を確保するための乱数シード設定（オプション）
    if (!is.null(seed)) set.seed(seed)
    
    # ブートストラップ法を用いたCliff's Deltaの計算（繰り返し）
    bootstrap_results <- replicate(n_bootstrap, {
      sample_A <- sample(data_A, length(data_A), replace = TRUE)
      sample_B <- sample(data_B, length(data_B), replace = TRUE)
      cliff.delta(sample_A, sample_B)$estimate
    })
    
    # ブートストラップ分布から95％信頼区間を算出
    alpha <- 1 - conf_level
    ci_lower <- quantile(bootstrap_results, alpha / 2)
    ci_upper <- quantile(bootstrap_results, 1 - alpha / 2)
    
    # 結果をリストとして返す
    list(
      original_delta = original_delta,
      delta_bootstrap_mean = mean(bootstrap_results),
      ci_lower_bootstrap = ci_lower,
      ci_upper_bootstrap = ci_upper
    )
  }
  
  # 各変数と年齢グループ間の比較に対してブートストラップ分析を実施
  for (var in non_normal_vars) {
    cat("\nBootstrap analysis for", var, "(Cliff's Delta)\n")
    
    # 各グループのデータ取得（欠損値を除去）
    group1_var <- as.numeric(na.omit(group1[[var]]))
    group2_var <- as.numeric(na.omit(group2[[var]]))
    group3_var <- as.numeric(na.omit(group3[[var]]))
    
    # 1-6歳未満 vs 6-7歳未満の比較
    boot_result_1_2 <- bootstrap_cliffs_delta(var, group2, group1, n_bootstrap = 5000, seed = NULL)
    cat("1-6yo vs 6-7yo Cliff's Delta (Original):", round(boot_result_1_2$original_delta, 4), "\n")
    cat("   Cliff's Delta (Bootstrap mean):", round(boot_result_1_2$delta_bootstrap_mean, 4), "\n")
    cat("   95% CI (bootstrap):", round(boot_result_1_2$ci_lower_bootstrap, 4), "-", round(boot_result_1_2$ci_upper_bootstrap, 4), "\n\n")
    
    # 1-6歳未満 vs 7-11歳未満の比較
    boot_result_1_3 <- bootstrap_cliffs_delta(var, group3, group1, n_bootstrap = 5000, seed = NULL)
    cat("1-6yo vs 7-11yo Cliff's Delta (Original):", round(boot_result_1_3$original_delta, 4), "\n")
    cat("   Cliff's Delta (Bootstrap mean):", round(boot_result_1_3$delta_bootstrap_mean, 4), "\n")
    cat("   95% CI (bootstrap):", round(boot_result_1_3$ci_lower_bootstrap, 4), "-", round(boot_result_1_3$ci_upper_bootstrap, 4), "\n\n")
    
    # 6-7歳未満 vs 7-11歳未満の比較
    boot_result_2_3 <- bootstrap_cliffs_delta(var, group3, group2, n_bootstrap = 5000, seed = NULL)
    cat("6-7yo vs 7-11yo Cliff's Delta (Original):", round(boot_result_2_3$original_delta, 4), "\n")
    cat("   Cliff's Delta (Bootstrap mean):", round(boot_result_2_3$delta_bootstrap_mean, 4), "\n")
    cat("   95% CI (bootstrap):", round(boot_result_2_3$ci_lower_bootstrap, 4), "-", round(boot_result_2_3$ci_upper_bootstrap, 4), "\n\n")
  }
}

# Cliff's delta is a non-parametric method used when data violate assumptions of normality.
# For parameters "Sd" and "Cov", the assumption of normality was initially violated due to a few extreme outliers.
# Therefore, these parameters were subsequently re-analyzed after excluding such outliers (sensitivity analysis).
# Effect size benchmarks for Cliff's delta are: Small (≥ 0.147), Medium (≥ 0.33), Large (≥ 0.474).


# PART 4: Glass's Delta-------------------------------------------------------------
{
  # 【グループ分け】
  # ABx（年齢）を基準として3つのグループ（1-6歳未満、6-7歳未満、7-11歳未満）に分割
  group1 <- dataAll %>% filter(ABx < 6)
  group2 <- dataAll %>% filter(ABx >= 6 & ABx < 7)
  group3 <- dataAll %>% filter(ABx >= 7 & ABx < 11)
  
  # 【分析対象変数リスト（Glass's Deltaを計算する変数）】
  variables_glass <- c("MFD", "CFA")
  
  # 【Glass's Delta、95％信頼区間、観測された検出力、必要なサンプルサイズを計算する関数】
  calculate_glass_delta <- function(var, ref_group, comp_group) {
    # 欠損値を除去し、数値ベクトルに変換
    ref_var <- as.numeric(na.omit(ref_group[[var]]))
    comp_var <- as.numeric(na.omit(comp_group[[var]]))
    
    # Glass's Deltaを計算（比較群と基準群の平均の差を基準群の標準偏差で割ったもの）
    effect_size <- (mean(comp_var) - mean(ref_var)) / sd(ref_var)
    
    # MBESSパッケージのci.smd関数を用いて95%信頼区間を算出
    ci_lower_upper <- ci.smd(
      smd = effect_size,
      n.1 = length(comp_var),
      n.2 = length(ref_var),
      conf.level = 0.95
    )
    
    # pwrパッケージを使用し、80%の検出力を達成するために必要なサンプルサイズを算出
    required_n <- pwr.t.test(
      d = effect_size,
      power = 0.8,
      sig.level = 0.05,
      type = "two.sample"
    )$n
    
    # 観測された検出力を計算（実際のサンプルサイズを使用）
    if (length(comp_var) > 1 && length(ref_var) > 1) {
      observed_power <- pwr.t2n.test(
        n1 = length(ref_var),
        n2 = length(comp_var),
        d = effect_size,
        sig.level = 0.05,
        alternative = "two.sided"
      )$power
    } else {
      observed_power <- NA
    }
    
    list(
      effect_size = effect_size,
      ci = ci_lower_upper,
      required_n = required_n,
      observed_power = observed_power
    )
  }
  
  # 【各変数についてGlass's Deltaの計算を繰り返し実施】
  for (var in variables_glass) {
    cat("\nAnalyzing", var, "(Glass's Delta)\n")
    
    # 各比較（1-2、1-3、2-3）で計算
    result_1_2 <- calculate_glass_delta(var, group1, group2)
    result_1_3 <- calculate_glass_delta(var, group1, group3)
    result_2_3 <- calculate_glass_delta(var, group2, group3)
    
    # 結果の表示
    # 1) 1-6歳未満 vs 6-7歳未満
    cat("1-6yo vs 6-7yo Glass's Delta:", result_1_2$effect_size, "\n")
    cat("95% CI:", result_1_2$ci$Lower.Conf.Limit.smd, "-", result_1_2$ci$Upper.Conf.Limit.smd, "\n")
    cat("Observed power:", result_1_2$observed_power, "\n")
    cat("Required sample size (per group, for 80% power):", result_1_2$required_n, "\n\n")
    
    # 2) 1-6歳未満 vs 7-11歳未満
    cat("1-6yo vs 7-11yo Glass's Delta:", result_1_3$effect_size, "\n")
    cat("95% CI:", result_1_3$ci$Lower.Conf.Limit.smd, "-", result_1_3$ci$Upper.Conf.Limit.smd, "\n")
    cat("Observed power:", result_1_3$observed_power, "\n")
    cat("Required sample size (per group, for 80% power):", result_1_3$required_n, "\n\n")
    
    # 3) 6-7歳未満 vs 7-11歳未満
    cat("6-7yo vs 7-11yo Glass's Delta:", result_2_3$effect_size, "\n")
    cat("95% CI:", result_2_3$ci$Lower.Conf.Limit.smd, "-", result_2_3$ci$Upper.Conf.Limit.smd, "\n")
    cat("Observed power:", result_2_3$observed_power, "\n")
    cat("Required sample size (per group, for 80% power):", result_2_3$required_n, "\n")
  }
}


# Function to compute a Bootstrap CI for Glass's Delta---
{
  # 【Glass's Deltaのブートストラップ信頼区間を計算する関数】
  bootstrap_ci_glass_delta <- function(ref_group, comp_group, n_iter = 5000, conf_level = 0.95, seed = NULL) {
    
    # 再現可能性のための乱数シード設定（オプション）
    if (!is.null(seed)) set.seed(seed)
    
    n_ref  <- length(ref_group)   # 基準グループのサンプルサイズ
    n_comp <- length(comp_group)  # 比較グループのサンプルサイズ
    boot_deltas <- numeric(n_iter)  # ブートストラップ結果格納用のベクトル
    
    # オリジナルのGlass's Deltaを計算（ブートストラップなし）
    original_delta <- (mean(comp_group) - mean(ref_group)) / sd(ref_group)
    
    # ブートストラップを用いたGlass's Deltaの計算（繰り返し）
    for (i in 1:n_iter) {
      ref_sample  <- sample(ref_group, n_ref, replace = TRUE)
      comp_sample <- sample(comp_group, n_comp, replace = TRUE)
      boot_deltas[i] <- (mean(comp_sample) - mean(ref_sample)) / sd(ref_sample)
    }
    
    # ブートストラップ分布から95％信頼区間を算出
    lower_bound <- quantile(boot_deltas, (1 - conf_level) / 2)
    upper_bound <- quantile(boot_deltas, 1 - (1 - conf_level) / 2)
    
    # 計算結果をリスト形式で返す
    return(list(
      original_glass_delta = original_delta,
      boot_mean            = mean(boot_deltas),
      ci_lower             = lower_bound,
      ci_upper             = upper_bound
    ))
  }
  
  # 【実際のデータを用いた例（group1, group2, group3を使用）】
  variables_glass <- c("MFD", "CFA") # 分析対象変数を指定
  
  for (var in variables_glass) {
    cat("\nBootstrap Analysis for", var, "(Glass's Delta)\n")
    
    # 各グループの欠損値を除去し数値ベクトルとして抽出
    group1_var <- as.numeric(na.omit(group1[[var]]))
    group2_var <- as.numeric(na.omit(group2[[var]]))
    group3_var <- as.numeric(na.omit(group3[[var]]))
    
    # 1-6歳未満 vs 6-7歳未満
    result_1_2 <- bootstrap_ci_glass_delta(group1_var, group2_var, n_iter = 5000, seed = NULL)
    cat("1-6yo vs 6-7yo Glass's Delta (Original):", result_1_2$original_glass_delta, "\n")
    cat("   Bootstrap mean:", result_1_2$boot_mean, "\n")
    cat("   95% Bootstrap CI:", result_1_2$ci_lower, "-", result_1_2$ci_upper, "\n\n")
    
    # 1-6歳未満 vs 7-11歳未満
    result_1_3 <- bootstrap_ci_glass_delta(group1_var, group3_var, n_iter = 5000, seed = NULL)
    cat("1-6yo vs 7-11yo Glass's Delta (Original):", result_1_3$original_glass_delta, "\n")
    cat("   Bootstrap mean:", result_1_3$boot_mean, "\n")
    cat("   95% Bootstrap CI:", result_1_3$ci_lower, "-", result_1_3$ci_upper, "\n\n")
    
    # 6-7歳未満 vs 7-11歳未満
    result_2_3 <- bootstrap_ci_glass_delta(group2_var, group3_var, n_iter = 5000, seed = NULL)
    cat("6-7yo vs 7-11yo Glass's Delta (Original):", result_2_3$original_glass_delta, "\n")
    cat("   Bootstrap mean:", result_2_3$boot_mean, "\n")
    cat("   95% Bootstrap CI:", result_2_3$ci_lower, "-", result_2_3$ci_upper, "\n")
  }
}


# Glass's delta is an effect size measure used when variances between groups are unequal.
# In this study, Glass's delta was specifically calculated for "MFD" and "CFA".
# Effect size benchmarks are: Small (0.2), Medium (0.5), Large (≥ 0.8).
#
# Interpretation of bootstrap-derived 95% confidence intervals (CI):
# For "MFD", bootstrap analysis notably narrowed the 95% CI,
# reflecting consistent group differences and stability of effect size despite larger variance within the 1–6-year group.
# Conversely, for "CFA", the bootstrap approach widened the 95% CI,
# indicating greater uncertainty and variability in the magnitude of effect size estimation.
# These contrasting results highlight important differences in data characteristics:
# "MFD" exhibits stable group-level effects despite higher variance at younger ages,
# while "CFA" shows less consistent group differences across individuals.


### -----------------------------### 8. Sensitivity Analysis (Effect size)---------
# PART 1: Cov with Hedges' d----------------------

# 【外れ値を除去したデータセットの作成】
dataCov_noOutlier <- dataCov[-c(8, 26), ]  # 8番目と26番目のデータを外れ値として除去

# 【外れ値を除去したグループごとのデータ抽出（年齢区分）】
Cov_Group1_noOutlier <- dataCov_noOutlier %>%
  filter(ABx < 6) %>%
  pull(Cov)

Cov_Group2_noOutlier <- dataCov_noOutlier %>%
  filter(ABx >= 6 & ABx < 7) %>%
  pull(Cov)

Cov_Group3_noOutlier <- dataCov_noOutlier %>%
  filter(ABx >= 7 & ABx < 11) %>%
  pull(Cov)

# 【正規性の検定（Shapiro-Wilk検定）】
shapiro.test(Cov_Group1_noOutlier) # p-value = 0.5169（正規性あり）
shapiro.test(Cov_Group3_noOutlier) # p-value = 0.3566（正規性あり）

{
  # 【1. データのコピーと外れ値の除去】
  dataCov_HedgesNoOutlier <- dataCov
  dataCov_HedgesNoOutlier <- dataCov_HedgesNoOutlier[-c(8, 26), ]
  
  # 【2. グループ分け】
  group1_Cov_HedgesNoOutlier <- dataCov_HedgesNoOutlier %>% filter(ABx < 6)
  group2_Cov_HedgesNoOutlier <- dataCov_HedgesNoOutlier %>% filter(ABx >= 6 & ABx < 7)
  group3_Cov_HedgesNoOutlier <- dataCov_HedgesNoOutlier %>% filter(ABx >= 7 & ABx < 11)
  
  # 【3. 対象変数とHedges' d計算用の関数定義】
  variables_Cov_HedgesNoOutlier <- c("Cov")
  
  # ブートストラップ法を用いてHedges' dの95％信頼区間を算出する関数
  bootstrap_ci_Cov_Hedges <- function(data1, data2, n_iter = 5000, conf_level = 0.95) {
    n1 <- length(data1)
    n2 <- length(data2)
    boot_vals <- numeric(n_iter)
    for (i in seq_len(n_iter)) {
      s1 <- sample(data1, n1, replace = TRUE)
      s2 <- sample(data2, n2, replace = TRUE)
      boot_vals[i] <- cohen.d(s1, s2, hedges.correction = TRUE)$estimate
    }
    c(
      quantile(boot_vals, (1 - conf_level) / 2),
      quantile(boot_vals, 1 - (1 - conf_level) / 2)
    )
  }
  
  # 各変数に対してHedges' dおよび95％ブートストラップ信頼区間の計算
  for (v in variables_Cov_HedgesNoOutlier) {
    cat("\n[Cov] Hedges' d (NoOutlier)\n")
    
    g1 <- as.numeric(na.omit(group1_Cov_HedgesNoOutlier[[v]]))
    g2 <- as.numeric(na.omit(group2_Cov_HedgesNoOutlier[[v]]))
    g3 <- as.numeric(na.omit(group3_Cov_HedgesNoOutlier[[v]]))
    
    # Hedges' d計算（外れ値なし）
    d_1_2 <- cohen.d(g2, g1, hedges.correction = TRUE)$estimate
    d_1_3 <- cohen.d(g3, g1, hedges.correction = TRUE)$estimate
    d_2_3 <- cohen.d(g3, g2, hedges.correction = TRUE)$estimate
    
    # ブートストラップ法による95％信頼区間の計算
    ci_1_2_boot <- bootstrap_ci_Cov_Hedges(g2, g1)
    ci_1_3_boot <- bootstrap_ci_Cov_Hedges(g3, g1)
    ci_2_3_boot <- bootstrap_ci_Cov_Hedges(g3, g2)
    
    # 結果の表示
    cat("1-6yo vs 6-7yo:", d_1_2, "\n",
        " 95% Boot CI:", ci_1_2_boot[1], "-", ci_1_2_boot[2], "\n")
    
    cat("1-6yo vs 7-11yo:", d_1_3, "\n",
        " 95% Boot CI:", ci_1_3_boot[1], "-", ci_1_3_boot[2], "\n")
    
    cat("6-7yo vs 7-11yo:", d_2_3, "\n",
        " 95% Boot CI:", ci_2_3_boot[1], "-", ci_2_3_boot[2], "\n")
  }
}


# Sensitivity analysis for "Cov":
# After excluding clear outliers, data stability significantly improved, 
# particularly in the youngest group (1–6 years).
# This resulted in notably improved (narrower) confidence intervals for the effect size,
# especially when comparing age groups 1–6 and 6–7.


# PART 2: Sd with Hedges' d------------------------------
{
  # 【外れ値を除去したデータセットの作成】
  dataAll_SdNoOutlier <- dataSd[-33, ]  # 33番目のデータを外れ値として除去
  
  # 【外れ値除去後のグループごとデータ抽出（年齢区分）】
  group1_SdNoOutlier <- dataAll_SdNoOutlier %>% filter(ABx < 6)
  group2_SdNoOutlier <- dataAll_SdNoOutlier %>% filter(ABx >= 6 & ABx < 7)
  group3_SdNoOutlier <- dataAll_SdNoOutlier %>% filter(ABx >= 7 & ABx < 11)
  
  # グループごとの標準偏差（Sd）の抽出
  Sd_Group1NoOutlier <- group1_SdNoOutlier %>% pull(Sd)
  Sd_Group2NoOutlier <- group2_SdNoOutlier %>% pull(Sd)
  Sd_Group3NoOutlier <- group3_SdNoOutlier %>% pull(Sd)
  
  # 【正規性の検定（Shapiro-Wilk検定）】
  shapiro.test(Sd_Group3NoOutlier) # p-value = 0.07849（正規性あり）
  
  {
    # 【1. データのコピーと外れ値の除去】
    dataSd_HedgesNoOutlier <- dataAll
    dataSd_HedgesNoOutlier <- dataSd_HedgesNoOutlier[-33, ]
    
    # 【2. グループ分け】
    group1_Sd_HedgesNoOutlier <- dataSd_HedgesNoOutlier %>% filter(ABx < 6)
    group2_Sd_HedgesNoOutlier <- dataSd_HedgesNoOutlier %>% filter(ABx >= 6 & ABx < 7)
    group3_Sd_HedgesNoOutlier <- dataSd_HedgesNoOutlier %>% filter(ABx >= 7 & ABx < 11)
    
    # 【3. 対象変数とHedges' d計算用の関数定義】
    variables_Sd_HedgesNoOutlier <- c("Sd")
    
    # 【4. ブートストラップ法を用いたHedges' dの95％信頼区間を算出する関数】
    bootstrap_ci_Sd_Hedges <- function(data1, data2, n_iter = 5000, conf_level = 0.95) {
      n1 <- length(data1)
      n2 <- length(data2)
      boot_vals <- numeric(n_iter)
      for (i in seq_len(n_iter)) {
        s1 <- sample(data1, n1, replace = TRUE)
        s2 <- sample(data2, n2, replace = TRUE)
        boot_vals[i] <- cohen.d(s1, s2, hedges.correction = TRUE)$estimate
      }
      c(
        quantile(boot_vals, (1 - conf_level) / 2),
        quantile(boot_vals, 1 - (1 - conf_level) / 2)
      )
    }
    
    # 各変数に対するHedges' dと95％ブートストラップ信頼区間の計算
    for (v in variables_Sd_HedgesNoOutlier) {
      cat("\n[Sd] Hedges' d (NoOutlier: row 33)\n")
      
      g1 <- as.numeric(na.omit(group1_Sd_HedgesNoOutlier[[v]]))
      g2 <- as.numeric(na.omit(group2_Sd_HedgesNoOutlier[[v]]))
      g3 <- as.numeric(na.omit(group3_Sd_HedgesNoOutlier[[v]]))
      
      # Hedges' d計算（外れ値なし）
      d_1_2 <- cohen.d(g2, g1, hedges.correction = TRUE)$estimate
      d_1_3 <- cohen.d(g3, g1, hedges.correction = TRUE)$estimate
      d_2_3 <- cohen.d(g3, g2, hedges.correction = TRUE)$estimate
      
      # ブートストラップ法による95％信頼区間の計算
      ci_1_2_boot <- bootstrap_ci_Sd_Hedges(g2, g1)
      ci_1_3_boot <- bootstrap_ci_Sd_Hedges(g3, g1)
      ci_2_3_boot <- bootstrap_ci_Sd_Hedges(g3, g2)
      
      # 結果の表示
      cat("1-6yo vs 6-7yo:", d_1_2, "\n",
          " 95% Boot CI:", ci_1_2_boot[1], "-", ci_1_2_boot[2], "\n")
      
      cat("1-6yo vs 7-11yo:", d_1_3, "\n",
          " 95% Boot CI:", ci_1_3_boot[1], "-", ci_1_3_boot[2], "\n")
      
      cat("6-7yo vs 7-11yo:", d_2_3, "\n",
          " 95% Boot CI:", ci_2_3_boot[1], "-", ci_2_3_boot[2], "\n")
    }
  }
}


# Sensitivity analysis for "Sd":
# Similar to "Cov", removing clear outliers resulted in stable effect sizes 
# when comparing age groups 1–6 and 7–11.
# However, for the comparison between age groups 1–6 and 6–7, 
# the confidence interval remained wide, indicating persistent variability.


# PART 3: Mean with Hedges' d--------------------------------
{
  # 【外れ値を除去したデータセットの作成】
  dataAll_MeanNoOutlier <- dataMean[-33, ]  # 33番目のデータを外れ値として除去
  
  # 【外れ値除去後のグループごとデータ抽出（年齢区分）】
  group1_MeanNoOutlier <- dataAll_MeanNoOutlier %>% filter(ABx < 6)
  group2_MeanNoOutlier <- dataAll_MeanNoOutlier %>% filter(ABx >= 6 & ABx < 7)
  group3_MeanNoOutlier <- dataAll_MeanNoOutlier %>% filter(ABx >= 7 & ABx < 11)
  
  # グループごとの平均値（Mean）の抽出
  Mean_Group1NoOutlier <- group1_MeanNoOutlier %>% pull(Mean)
  Mean_Group2NoOutlier <- group2_MeanNoOutlier %>% pull(Mean)
  Mean_Group3NoOutlier <- group3_MeanNoOutlier %>% pull(Mean)
  
  # 【正規性の検定（Shapiro-Wilk検定）】
  shapiro.test(Mean_Group3NoOutlier) # p-value = 0.8786（正規性あり）
  
  {
    # 【1. データのコピーと外れ値の除去】
    dataMean_HedgesNoOutlier <- dataAll
    dataMean_HedgesNoOutlier <- dataMean_HedgesNoOutlier[-33, ]
    
    # 【2. グループ分け】
    group1_Mean_HedgesNoOutlier <- dataMean_HedgesNoOutlier %>% filter(ABx < 6)
    group2_Mean_HedgesNoOutlier <- dataMean_HedgesNoOutlier %>% filter(ABx >= 6 & ABx < 7)
    group3_Mean_HedgesNoOutlier <- dataMean_HedgesNoOutlier %>% filter(ABx >= 7 & ABx < 11)
    
    # 【3. 対象変数】
    variables_Mean_HedgesNoOutlier <- c("Mean")
    
    # 【4. ブートストラップ法を用いたHedges' dの95％信頼区間を算出する関数】
    bootstrap_ci_Mean_Hedges <- function(data1, data2, n_iter = 5000, conf_level = 0.95) {
      n1 <- length(data1)
      n2 <- length(data2)
      boot_vals <- numeric(n_iter)
      for (i in seq_len(n_iter)) {
        s1 <- sample(data1, n1, replace = TRUE)
        s2 <- sample(data2, n2, replace = TRUE)
        boot_vals[i] <- cohen.d(s1, s2, hedges.correction = TRUE)$estimate
      }
      c(
        quantile(boot_vals, (1 - conf_level) / 2),
        quantile(boot_vals, 1 - (1 - conf_level) / 2)
      )
    }
    
  }
  
  # 【5. Hedges' dの計算と95％ブートストラップ信頼区間の表示】
  for (v in variables_Mean_HedgesNoOutlier) {
    cat("\n[Mean] Hedges' d (NoOutlier: row 33)\n")
    
    g1 <- as.numeric(na.omit(group1_Mean_HedgesNoOutlier[[v]]))
    g2 <- as.numeric(na.omit(group2_Mean_HedgesNoOutlier[[v]]))
    g3 <- as.numeric(na.omit(group3_Mean_HedgesNoOutlier[[v]]))
    
    # Hedges' d計算（外れ値なし）
    d_1_2 <- cohen.d(g2, g1, hedges.correction = TRUE)$estimate
    d_1_3 <- cohen.d(g3, g1, hedges.correction = TRUE)$estimate
    d_2_3 <- cohen.d(g3, g2, hedges.correction = TRUE)$estimate
    
    # ブートストラップ法による95％信頼区間の計算
    ci_1_2_boot <- bootstrap_ci_Mean_Hedges(g2, g1)
    ci_1_3_boot <- bootstrap_ci_Mean_Hedges(g3, g1)
    ci_2_3_boot <- bootstrap_ci_Mean_Hedges(g3, g2)
    
    # 結果の表示
    cat("1-6yo vs 6-7yo:", d_1_2, "\n",
        " 95% Boot CI:", ci_1_2_boot[1], "-", ci_1_2_boot[2], "\n")
    cat("1-6yo vs 7-11yo:", d_1_3, "\n",
        " 95% Boot CI:", ci_1_3_boot[1], "-", ci_1_3_boot[2], "\n")
    cat("6-7yo vs 7-11yo:", d_2_3, "\n",
        " 95% Boot CI:", ci_2_3_boot[1], "-", ci_2_3_boot[2], "\n")
  }
}


# Sensitivity analysis for "Mean":
# Outliers were found in the same samples identified for "Sd" and were similarly excluded.
# Post-exclusion, "Mean" demonstrated a comparable pattern of CI width to that observed with "Sd",
# suggesting consistent data characteristics across these parameters.


# PART 4: Fat--------------------------------------------------
dataAll_FatNoOutlier <- dataFat[-8, ]

group1_FatNoOutlier <- dataAll_FatNoOutlier %>% filter(ABx < 6)
group2_FatNoOutlier <- dataAll_FatNoOutlier %>% filter(ABx >= 6 & ABx < 7)
group3_FatNoOutlier <- dataAll_FatNoOutlier %>% filter(ABx >= 7 & ABx < 11)

Fat_Group1NoOutlier <- group1_FatNoOutlier %>% pull(Fat)
Fat_Group2NoOutlier <- group2_FatNoOutlier %>% pull(Fat)
Fat_Group3NoOutlier <- group3_FatNoOutlier %>% pull(Fat)

shapiro.test(Fat_Group1NoOutlier) # p-value = 0.02856 (!)
shapiro.test(Fat_Group2NoOutlier) # p-value = 0.2932
shapiro.test(Fat_Group3NoOutlier) # p-value = 0.06398

# Sensitivity analysis for "Fat":
# After removal of extreme outliers, the assumption of normality was still violated.
# Consequently, no further effect size calculations under parametric assumptions
# were performed for this parameter.


### -----------------------------### 9. Replication of Peverelli et al. (2015)---------
# PART 1: 1st simulations--------------------------------------------------------
# 【仮定条件】
#  - Group1 (1–6歳): 平均値 = 16.98, 標準偏差 = 4.90, サンプルサイズ = 24
#  - Group2 (7–10歳): 平均値 = 30.00, 標準偏差 = 8.65, サンプルサイズ = 16
# Cohen's dを計算（等分散を仮定）し、その95%信頼区間を近似的に推定する。
{
  cohen_d_ci <- function(mean1, mean2, sd1, sd2, n1, n2, alpha = 0.05) {
    # 1) 等分散仮定のもと、プールされた標準偏差を計算
    pooled_sd <- sqrt( ((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2) )
    
    # 2) Cohen's dの計算
    d_value <- (mean2 - mean1) / pooled_sd
    
    # 3) Cohen's dの近似分散を計算（Cumming 2012の式を使用）
    var_d <- ((n1 + n2) / (n1 * n2)) + (d_value^2 / (2 * (n1 + n2 - 2)))
    
    # 4) 標準誤差（SE）を算出
    se_d <- sqrt(var_d)
    
    # 5) 信頼区間を求めるための臨界z値を算出
    z_crit <- qnorm(1 - alpha/2)
    
    # 6) 95%信頼区間を算出
    lower_ci <- d_value - z_crit * se_d
    upper_ci <- d_value + z_crit * se_d
    
    return(list(
      d = d_value,
      SE = se_d,
      lower_95 = lower_ci,
      upper_95 = upper_ci
    ))
  }
  
  # 各グループのパラメータを設定
  mean1 <- 16.98; sd1 <- 4.90; n1 <- 24  # 1–6歳グループ
  mean2 <- 30.00; sd2 <- 8.65; n2 <- 16  # 7–10歳グループ
  
  # Cohen's d とその95％信頼区間を計算
  result <- cohen_d_ci(mean1, mean2, sd1, sd2, n1, n2)
  
  # 結果を小数点以下3桁で表示
  cat("Cohen's d (1st simulation):", round(result$d, 3), "\n")
  cat("95% CI: [", round(result$lower_95, 3), ",", round(result$upper_95, 3), "]\n")
}

# 【Glass's Deltaの仮定条件と計算】
#   - "対照群" (1–6歳): 平均値 = 16.98, 標準偏差 = 4.90, サンプルサイズ = 24
#   - "処置群" (7–10歳): 平均値 = 30.00, 標準偏差 = 8.65, サンプルサイズ = 16
#   - Glass's Deltaは対照群の標準偏差を用いて算出
#   - Cohen's dと類似の方法で近似的に95%信頼区間を算出（厳密ではないがよく使用される）

{
  glass_delta_ci <- function(mean_control, mean_treatment, sd_control, n_control, n_treatment, alpha = 0.05) {
    # 1) Glass's Deltaの計算
    delta_value <- (mean_treatment - mean_control) / sd_control
    
    # 2) 近似分散をCohen's dの式で代用
    var_delta <- (n_control + n_treatment) / (n_control * n_treatment) + (delta_value^2) / (2 * (n_control + n_treatment - 2))
    
    # 3) 標準誤差の算出
    se_delta <- sqrt(var_delta)
    
    # 4) 信頼区間の臨界z値を算出
    z_crit <- qnorm(1 - alpha/2)
    
    # 5) 95%信頼区間を算出
    lower_ci <- delta_value - z_crit * se_delta
    upper_ci <- delta_value + z_crit * se_delta
    
    return(list(
      delta = delta_value,
      SE = se_delta,
      lower_95 = lower_ci,
      upper_95 = upper_ci
    ))
  }
  
  # 各グループのパラメータ設定
  mean1 <- 16.98; sd1 <- 4.90; n1 <- 24   # "対照群"
  mean2 <- 30.00; sd2 <- 8.65; n2 <- 16   # "処置群"
  
  # Glass's Deltaの計算
  result_glass <- glass_delta_ci(mean_control = mean1, mean_treatment = mean2, sd_control = sd1, n_control = n1, n_treatment = n2)
  
  # 結果の表示（小数点以下3桁）
  cat("Glass's Delta (1st simulation):", round(result_glass$delta, 3), "\n")
  cat("Approx. 95% CI: [", round(result_glass$lower_95, 3), ",", round(result_glass$upper_95, 3), "]\n")
}


# PART 2: 2nd simulations-------------------------------------------------------
# 仮定条件:
# - 1歳群: SD = 7.84（1-6歳患者群のSDとして仮定）
# - 7歳群: SD = 7.75（7-10歳患者群のSDとして仮定）
# - 1-6歳の平均 = 16.98% (n = 24), 7-10歳の平均 = 30% (n = 16)

# Cohen's dの計算
{
  # Cohen's dを計算する関数
  cohen_d <- function(mean1, mean2, sd1, sd2, n1, n2) {
    # プールされた標準偏差を計算（等分散仮定）
    pooled_sd <- sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
    # Cohen's dを算出
    d <- (mean2 - mean1) / pooled_sd
    return(d)
  }
  
  # 各グループのパラメータを設定
  mean1 <- 16.98  # 1-6歳群の平均
  mean2 <- 30.00  # 7-10歳群の平均
  sd1 <- 7.84     # 1-6歳群の標準偏差
  sd2 <- 7.75     # 7-10歳群の標準偏差
  n1 <- 24        # 1-6歳群のサンプルサイズ
  n2 <- 16        # 7-10歳群のサンプルサイズ
  
  # Cohen's dを計算し結果を表示
  d_value <- cohen_d(mean1, mean2, sd1, sd2, n1, n2)
  print(paste("Cohen's d (2nd simulation):", round(d_value, 3)))
}

# Glass's Deltaの計算
{
  # Glass's Deltaを計算する関数（対照群の標準偏差を使用）
  glass_delta <- function(mean1, mean2, sd_control) {
    delta <- (mean2 - mean1) / sd_control
    return(delta)
  }
  
  # パラメータ設定
  mean1 <- 16.98     # 対照群（1-6歳群）の平均
  mean2 <- 30.00     # 処置群（7-10歳群）の平均
  sd_control <- 7.84 # 対照群の標準偏差
  
  # Glass's Deltaを計算し結果を表示
  delta_value <- glass_delta(mean1, mean2, sd_control)
  print(paste("Glass's Delta (2nd simulation):", round(delta_value, 3)))
}


### -----------------------------### 10. Logistic Regression Analysis---------------
# PART 1: MFD Threshold---------------------------------------------------------
{
  # 【データ準備】ABxが11歳未満のデータを抽出し、データフレームを作成
  data_logis <- data.frame(ABx, MFD)
  data_logis <- subset(data_logis, ABx < 11)
  options(warn = -1) # glmの完全分離に関する警告を無視（安全）
  
  # 【年齢の二値化】年齢の境界を6歳として二値化（4歳、5歳、7歳など他の年齢での調整可能）
  data_logis$AgeBinary <- ifelse(data_logis$ABx < 6, 0, 1)
  results <- data.frame()
  
  # 【MFD閾値の設定】閾値を500～600まで1ずつ増加させて試行
  thresholds <- seq(500, 600, by = 1)
  
  # 各閾値に対してロジスティック回帰を実施
  for (threshold in thresholds) {
    data_logis$MFDBinary <- ifelse(data_logis$MFD > threshold, 1, 0)
    
    # 【ロジスティック回帰モデルの適合】AgeBinaryをMFDBinaryで予測
    model <- glm(AgeBinary ~ MFDBinary, data = data_logis, family = binomial)
    
    # 【統計指標の計算】AIC、p値、オッズ比、95％信頼区間を計算
    aic_value <- AIC(model)
    p_value <- summary(model)$coefficients[2, 4]
    odds_ratio <- exp(coef(model)[2])
    conf_int <- exp(confint(model))[2, ] 
    
    # 【分類精度の計算】
    predicted <- ifelse(predict(model, type = "response") > 0.5, 1, 0)
    accuracy <- mean(predicted == data_logis$AgeBinary) * 100
    
    # 結果を保存（有意かどうかで表記を分ける）
    if (p_value < 0.05) {
      results <- rbind(results, data.frame(
        Threshold = threshold, AIC = aic_value, P_Value = format(p_value, scientific = FALSE), 
        Accuracy = accuracy, Odds_Ratio = odds_ratio, CI_Lower = conf_int[1], CI_Upper = conf_int[2]
      ))
    } else {
      results <- rbind(results, data.frame(
        Threshold = threshold, AIC = aic_value, P_Value = format(p_value, scientific = FALSE), 
        Accuracy = accuracy, Odds_Ratio = "NS", CI_Lower = "NS", CI_Upper = "NS"
      ))
    }
  }
  
  # 結果を表示
  print(results)
  options(warn = 0)
}

# Bonferroni調整済みp値: 0.0004950495

# 546-551: p = 0.0007121682 (>0.000495)
# 552,553: p = 0.0003530315
# 554-559: p = 0.000179528
# 560-563: p = 0.0001518293
# 564-571: p = 0.0002070767
#     572: P = 0.0004800417

{
  # 年齢の境界値を6歳と設定し、二値化したデータフレームを作成
  age_break_log <- 6
  dat_log <- dataAll %>%             
    filter(ABx < 11) %>%               
    mutate(AgeBin = ifelse(ABx < age_break_log, 1, 0))
  
  # ロジスティック回帰モデルを適合（MFDによる年齢区分の予測）
  fit_log  <- glm(AgeBin ~ MFD, data = dat_log, family = binomial)
  
  # ロジット関数を定義
  logit <- function(p) log(p/(1-p))
  # MFDの閾値を計算（50％の分類境界）
  cut_log    <- (logit(0.50) - coef(fit_log)[1]) / coef(fit_log)[2]
  
  # 計算された閾値を表示
  print(cut_log) # 566.7945 
}



# Empirical Threshold Determination:
#  We conducted a detailed logistic regression analysis to determine the 
#  optimal threshold of MFD for distinguishing DMD patients younger than 6 years from
#  those aged 6 years or older. Thresholds ranging from 546 to 572 fibers/mm² 
#  were assessed using iterative logistic regression modeling, 
#  evaluating classification accuracy, OR, and Bonferroni-corrected
#  statistical significance (corrected α-level: p=0.000495).

# Optimal Threshold Range:
#  The most robust and statistically significant threshold range 
#  identified was 560–563 fibers/mm², which yielded the highest classification accuracy
#  of 94.3%, along with the lowest ORs (minimum OR = 0.0040; 95% CI: 0.00010–0.043)
#  and exceptional statistical significance (p = 0.000152).
#  This clearly delineates a highly precise and clinically relevant MFD threshold range.

# Alignment with Analytical Logistic Model:
#  Separately, we derived a theoretical logistic threshold by fitting
#  a continuous logistic regression model (without binarization), 
#  which provided an analytical solution for the MFD value corresponding to 
#  a 50% probability of being classified as younger than 6 years. 
#  This analytical approach yielded a threshold of approximately 566.8 fibers/mm²,
#  very close to the empirically identified optimal range (560–563 fibers/mm²).

# Clinical and Statistical Implications:
#  The close alignment between the empirically determined optimal threshold
#  (560–563 fibers/mm²) and the analytically derived threshold (566.8 fibers/mm²)
#  strongly supports the robustness and validity of these values. 
#  The empirical approach highlights the practical threshold that maximizes
#  clinical diagnostic accuracy, whereas the analytical approach identifies
#  the theoretical boundary point at which classification probability is exactly balanced (50%).

# Sensitivity Analyses:
#  Additional sensitivity analyses using alternative critical ages (4, 5, and 7 years)
#  consistently produced broader and less precise threshold ranges, 
#  confirming that the chosen critical age of 6 years provides the sharpest
#  and most clinically meaningful distinction.


# PART 2: 1-6/6-11 vs 1-7/7-11--------------------------------------------------
{
  # 【プロット領域の設定】左右に2つのプロットを並べて表示
  par(mfrow = c(1, 2), mar = c(5, 5, 2, 2))
  
  # 【1つ目の箱ひげ図】年齢区分 [1-6) と [6-11) によるMFD分布を表示
  boxplot(MFD ~ cut(ABx, breaks = c(1, 6, 11), labels = c("[1-6)", "[6-11)"), right = FALSE), 
          xlab = "Age Groups", ylab = "MFD", size = 4, outline = FALSE)
  # 各データポイントをジッターで追加（散布点）
  stripchart(MFD ~ cut(ABx, breaks = c(1, 6, 11), labels = c("[1-6)", "[6-11)"), right = FALSE), 
             method = "jitter", pch = 16, col = "black", vertical = TRUE, add = TRUE)
  
  # 【2つ目の箱ひげ図】年齢区分 [1-7) と [7-11) によるMFD分布を表示
  boxplot(MFD ~ cut(ABx, breaks = c(1, 7, 11), labels = c("[1-7)", "[7-11)"), right = FALSE), 
          xlab = "Age Groups", ylab = "MFD", size = 4, outline = FALSE)
  # 各データポイントをジッターで追加（散布点）
  stripchart(MFD ~ cut(ABx, breaks = c(1, 7, 11), labels = c("[1-7)", "[7-11)"), right = FALSE), 
             method = "jitter", pch = 16, col = "black", vertical = TRUE, add = TRUE)
  
  # プロット設定を元に戻す
  par(mfrow = c(1, 1))
}

# 【ggplotを使ったヒストグラム】MFDの分布を年齢二値（AgeBinary）ごとに表示
# 赤色の破線はMFDの閾値（ここでは560を例として示している）
ggplot(data_logis, aes(x = MFD, fill = as.factor(AgeBinary))) + 
  geom_histogram(binwidth = 10, alpha = 0.6, position = "identity") +  # ヒストグラム表示（透明度を調整し重ねて表示）
  geom_vline(xintercept = 560, color = "red", linetype = "dashed") +   # 閾値の位置を示す垂直線
  labs(title = "MFD Distribution by Age Groups", fill = "AgeBinary") +  # タイトルと凡例を追加
  scale_y_continuous(breaks = seq(0, max(2), by = 1))  # Y軸目盛りの調整


### -----------------------------### 11. Segmented Regression Analysis--------------
# PART 1: Preparation-------------------------------------------------------------------
{
  dataMean <- data.frame(ABx = ABx, Mean = Mean)
  dataSd <- data.frame(ABx = ABx, Sd = Sd)
  dataCov <- data.frame(ABx = ABx, Cov = Cov)
  dataMFD <- data.frame(ABx = ABx, MFD = MFD)
  dataMFA <- data.frame(ABx = ABx, MFA = MFA)
  dataFat <- data.frame(ABx = ABx, Fat = Fat)
  dataCFA <- data.frame(ABx = ABx, CFA = CFA)
}

{
  # Filter data for 1-11 years age group
  dataSegment_Mean <- dataMean %>%
    filter(ABx >= 1 & ABx <= 11)
  
  dataSegment_Sd <- dataSd %>%
    filter(ABx >= 1 & ABx <= 11)
  
  dataSegment_Cov <- dataCov %>%
    filter(ABx >= 1 & ABx <= 11)
  
  dataSegment_MFD <- dataMFD %>%
    filter(ABx >= 1 & ABx <= 11)
  
  dataSegment_MFA <- dataMFA %>%
    filter(ABx >= 1 & ABx <= 11)
  
  dataSegment_CFA <- dataCFA %>%
    filter(ABx >= 1 & ABx <= 11)
  
  dataSegment_Fat <- dataFat %>%
    filter(ABx >= 1 & ABx <= 11)
  
  # Fit initial linear model (simple linear regression)
  lm_init_Mean <- lm(Mean ~ ABx, data = dataSegment_Mean)
  lm_init_Sd <- lm(Sd ~ ABx, data = dataSegment_Sd)
  lm_init_Cov <- lm(Cov ~ ABx, data = dataSegment_Cov)
  lm_init_MFD <- lm(MFD ~ ABx, data = dataSegment_MFD)
  lm_init_MFA <- lm(MFA ~ ABx, data = dataSegment_MFA)
  lm_init_CFA <- lm(CFA ~ ABx, data = dataSegment_CFA)
  lm_init_Fat <- lm(Fat ~ ABx, data = dataSegment_Fat)
  
  # Fit segmented regression model with initial breakpoint at 6 years
  segmented_model_Mean <- segmented(lm_init_Mean, seg.Z = ~ABx, psi = 6)
  segmented_model_Sd <- segmented(lm_init_Sd, seg.Z = ~ABx, psi = 6)
  segmented_model_Cov <- segmented(lm_init_Cov, seg.Z = ~ABx, psi = 6)
  segmented_model_MFD <- segmented(lm_init_MFD, seg.Z = ~ABx, psi = 6)
  segmented_model_MFA <- segmented(lm_init_MFA, seg.Z = ~ABx, psi = 6)
  segmented_model_CFA <- segmented(lm_init_CFA, seg.Z = ~ABx, psi = 6)
  segmented_model_Fat <- segmented(lm_init_Fat, seg.Z = ~ABx, psi = 6)
}

# PART 2: Mean--------------------------------------------------------------------------
# Summarize segmented regression results
summary(segmented_model_Mean)

# Extract and visualize breakpoint
plot(segmented_model_Mean, main = "Segmented Regression for Mean vs Age")
abline(v = segmented_model_Mean$psi[2], col = "red", lty = 2)

# Extract slopes (before and after the breakpoint)
slope(segmented_model_Mean)

# Confidence intervals for breakpoint and slopes
confint(segmented_model_Mean)
anova(lm_init_Mean, segmented_model_Mean)

# PART 3: Sd----------------------------------------------------------------------------
# Summarize segmented regression results
summary(segmented_model_Sd)

# Extract and visualize breakpoint
plot(segmented_model_Sd, main = "Segmented Regression for Sd vs Age")
abline(v = segmented_model_Sd$psi[2], col = "red", lty = 2)

# Extract slopes (before and after the breakpoint)
slope(segmented_model_Sd)

# Confidence intervals for breakpoint and slopes
confint(segmented_model_Sd)
anova(lm_init_Sd, segmented_model_Sd)

# PART 4: Cov---------------------------------------------------------------------------
# Summarize segmented regression results
summary(segmented_model_Cov)

# Extract and visualize breakpoint
plot(segmented_model_Cov, main = "Segmented Regression for Cov vs Age")
abline(v = segmented_model_Cov$psi[2], col = "red", lty = 2)

# Extract slopes (before and after the breakpoint)
slope(segmented_model_Cov)

# Confidence intervals for breakpoint and slopes
confint(segmented_model_Cov)
anova(lm_init_Cov, segmented_model_Cov)

# PART 5: MFD---------------------------------------------------------------------------
# Summarize segmented regression results
summary(segmented_model_MFD) # breakpoint = 6.25 (St.Err 0.575), Adjusted R-squared: 0.7988

# Normality of residuals
autoplot(lm_init_MFD, smooth.colour = NA)
shapiro.test(residuals(lm_init_MFD)) # p-value = 0.7144

# Extract and visualize breakpoint
plot(segmented_model_MFD, main = "Segmented Regression for MFD vs Age")
abline(v = segmented_model_MFD$psi[2], col = "red", lty = 2)

# Extract slopes (before and after the breakpoint)
slope(segmented_model_MFD)
#           Est. St.Err. t value CI(95%).l CI(95%).u
# slope1 -170.68  20.951 -8.1463  -213.410  -127.950
# slope2  -28.31  23.612 -1.1990   -76.467    19.846

# Confidence intervals for breakpoint and slopes
confint(segmented_model_MFD) # 6.25 (5.07629  to  7.42371)
anova(lm_init_MFD, segmented_model_MFD) # Pr(>F) = 0.0003556 ***

dataSegment_MFD$Predicted_MFD <- predict(segmented_model_MFD)
RMSE_original <- sqrt(mean((dataSegment_MFD$MFD - dataSegment_MFD$Predicted_MFD)^2))
print(RMSE_original) # 132.9047

AIC_init <- AIC(lm_init_MFD)
AIC_seg <- AIC(segmented_model_MFD)
AIC_seg - AIC_init

# Sensitivity Analysis
{
  psi_values <- c(4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8)
  sensitivity_results <- sapply(psi_values, function(psi){
    seg_model <- segmented(lm_init_MFD, seg.Z = ~ABx, psi = psi)
    breakpoint <- seg_model$psi[2]
    return(breakpoint)
  })
  names(sensitivity_results) <- psi_values
  print(sensitivity_results)
} # 6.250000 to 6.250002


# ------- MFD Segmented model plot--------------------------------------------------
{
  bp   <- segmented_model_MFD$psi[2]            # 6.25 y
  xseq <- seq(0, 11, 0.1)
  pred <- predict(segmented_model_MFD,
                  newdata = data.frame(ABx = xseq),
                  se.fit  = FALSE)
  pred_df <- data.frame(ABx = xseq, fit = pred)
  
  # color setting
  pal_age <- gradient_n_pal(c("#3182BD", "#FEE08B", "#D73027"))
  
  col_line_pre  <- "#2AAA5E"   # regression line <6
  col_line_post <- "#666666"   # regression line >6
  col_bp        <- "#D72631"   # breakpoint
  fill_pre      <- "#E5F4EA"   # background
  fill_post     <- "#F2F2F2"   # background
  
  # plot and regression line
  ggplot(dataAll %>% filter(ABx <= 11), aes(x = ABx, y = MFD)) +
    # background
    annotate("rect", xmin = 0,  xmax = bp,
             ymin = -Inf, ymax = Inf, fill = fill_pre, alpha = 1) +
    annotate("rect", xmin = bp, xmax = 11,
             ymin = -Inf, ymax = Inf, fill = fill_post, alpha = 1) +
    
    # plot: age grad
    geom_point(aes(colour = ABx), size = 3, stroke = 0.1) +
    scale_colour_gradientn(
      colours = c("#3182BD", "#FEE08B", "#D73027"),
      limits = c(0, 11), name = "Age"
    ) +
    
    # segmented regression line
    geom_line(
      data = pred_df %>% filter(ABx <= bp),
      aes(x = ABx, y = fit),
      colour = col_line_pre, size = 1.4
    ) +
    geom_line(
      data = pred_df %>% filter(ABx >= bp),
      aes(x = ABx, y = fit),
      colour = col_line_post, size = 1.4
    ) +
    
    # breakpoint and annotation
    geom_vline(xintercept = bp,
               colour = col_bp, linetype = "dashed", size = 1) +
    annotate("text", x = bp + 0.15,
             y = max(dataAll$MFD[dataAll$ABx <= 11]) * 0.92,
             label = sprintf("Breakpoint  %.2f y", bp),
             hjust = 0, vjust = 1, angle = 0,
             colour = col_bp, size = 4.2, fontface = "bold") +
    
    # axis
    scale_x_continuous(limits = c(0, 11), breaks = 0:11) +
    scale_y_continuous(
      limits = c(0, max(dataAll$MFD[dataAll$ABx <= 11]) * 1.05),
      expand = expansion(mult = c(0, 0.02))
    ) +
    
    # label and theme
    labs(
      x = "ABx",
      y = "MFD"
    ) +
    theme_bw() +
    theme(
      axis.title = element_text(size = 20, face = "bold.italic"),
      axis.text  = element_text(size = 12),
      legend.position = "none",
      plot.margin = margin(8, 12, 8, 8)
    )
}

# ------- MFD LOOCV analysis----------------------------------------------------------
{
  # Prepare vectors to store LOOCV results
  n <- nrow(dataSegment_MFD)
  breakpoints <- numeric(n)
  CI_lower <- numeric(n)
  CI_upper <- numeric(n)
  
  # Perform LOOCV for segmented regression
  for (i in 1:n) {
    
    # Create LOOCV dataset by excluding one observation
    data_LOOCV <- dataSegment_MFD[-i, ]
    
    # Segmented regression model for LOOCV
    lm_LOOCV <- lm(MFD ~ ABx, data = data_LOOCV)
    
    # Perform segmented regression with initial breakpoint at age 6
    segmented_LOOCV <- tryCatch({
      segmented(lm_LOOCV, seg.Z = ~ABx, psi = 6)
    }, error = function(e) {
      return(NULL)
    })
    
    # If successful, save breakpoint and confidence intervals
    if (!is.null(segmented_LOOCV)) {
      breakpoints[i] <- segmented_LOOCV$psi[2]
      
      # Extract correct CI (lower and upper)
      CI_vals <- confint(segmented_LOOCV)[1, ]
      CI_lower[i] <- CI_vals[2]
      CI_upper[i] <- CI_vals[3]
    } else {
      # If model fails, assign NA
      breakpoints[i] <- NA
      CI_lower[i] <- NA
      CI_upper[i] <- NA
    }
  }
  
  # Remove NA results from LOOCV
  valid_breakpoints <- breakpoints[!is.na(breakpoints)]
  valid_CI_lower <- CI_lower[!is.na(CI_lower)]
  valid_CI_upper <- CI_upper[!is.na(CI_upper)]
  
  # Calculate LOOCV mean and median for breakpoint and confidence intervals
  mean_breakpoint <- mean(valid_breakpoints)
  median_breakpoint <- median(valid_breakpoints)
  mean_CI_lower <- mean(valid_CI_lower)
  mean_CI_upper <- mean(valid_CI_upper)
  
  # Output the summarized LOOCV results
  cat("LOOCV Mean Breakpoint:", mean_breakpoint, "\n")
  cat("LOOCV Median Breakpoint:", median_breakpoint, "\n")
  cat("LOOCV Mean CI:", mean_CI_lower, "to", mean_CI_upper, "\n")
  
  # Calculate range of breakpoints
  breakpoint_min <- min(breakpoints, na.rm = TRUE)
  breakpoint_max <- max(breakpoints, na.rm = TRUE)
  
  # Display range
  cat("Breakpoint Range:", breakpoint_min, "to", breakpoint_max, "\n")
  
  # Prediction by original data
  dataSegment_MFD$Predicted_MFD <- predict(segmented_model_MFD)
  
  # RMSE（original）
  RMSE_original <- sqrt(mean((dataSegment_MFD$MFD - dataSegment_MFD$Predicted_MFD)^2))
  cat("Original Data RMSE:", RMSE_original, "\n")
  
  # LOOCV analysis (RMSE calculation)
  {
    # Prepare vector to store predictions
    LOOCV_predictions <- numeric(n)
    
    # Perform LOOCV to get predicted values
    for (i in 1:n) {
      # Exclude one observation
      data_LOOCV <- dataSegment_MFD[-i, ]
      
      # Fit model
      lm_LOOCV <- lm(MFD ~ ABx, data = data_LOOCV)
      segmented_LOOCV <- tryCatch({
        segmented(lm_LOOCV, seg.Z = ~ABx, psi = 6)
      }, error = function(e) {
        return(NULL)
      })
      
      # If successful, predict the excluded observation
      if (!is.null(segmented_LOOCV)) {
        LOOCV_predictions[i] <- predict(segmented_LOOCV, 
                                        newdata = dataSegment_MFD[i,])
      } else {
        LOOCV_predictions[i] <- NA
      }
    }
    
    # Remove any NA predictions (just in case)
    valid_indices <- !is.na(LOOCV_predictions)
    LOOCV_actuals <- dataSegment_MFD$MFD[valid_indices]
    LOOCV_preds   <- LOOCV_predictions[valid_indices]
    
    # RMSE (LOOCV data)
    RMSE_LOOCV <- sqrt(mean((LOOCV_actuals - LOOCV_preds)^2))
    cat("LOOCV RMSE:", RMSE_LOOCV, "\n")
    
    # RMSE difference and ratio
    cat("Difference in RMSE (LOOCV - Original):", RMSE_LOOCV - RMSE_original, "\n")
    cat("RMSE Ratio (LOOCV / Original):", RMSE_LOOCV / RMSE_original, "\n")
  }
}

# LOOCV Mean Breakpoint: 6.256918 
# LOOCV Median Breakpoint: 6.250007 
# LOOCV Mean CI: 5.064994 to 7.448842
# Breakpoint Range: 5.504903 to 6.404331
# LOOCV RMSE: 151.4418 
# Difference in RMSE (LOOCV - Original): 18.53717 
# RMSE Ratio (LOOCV / Original): 1.139477 

# PART 6: MFA---------------------------------------------------------------------------
# Summarize segmented regression results
summary(segmented_model_MFA)

# Extract and visualize breakpoint
plot(segmented_model_MFA, main = "Segmented Regression for MFA vs Age")
abline(v = segmented_model_MFA$psi[2], col = "red", lty = 2)

# Extract slopes (before and after the breakpoint)
slope(segmented_model_MFA)

# Confidence intervals for breakpoint and slopes
confint(segmented_model_MFA)
anova(lm_init_MFA, segmented_model_MFA)

# PART 7: CFA---------------------------------------------------------------------------
# Summarize segmented regression results
summary(segmented_model_CFA)

# Extract and visualize breakpoint
plot(segmented_model_CFA, main = "Segmented Regression for CFA vs Age")
abline(v = segmented_model_CFA$psi[2], col = "red", lty = 2)

# Extract slopes (before and after the breakpoint)
slope(segmented_model_CFA)

# Confidence intervals for breakpoint and slopes
confint(segmented_model_CFA)
anova(lm_init_CFA, segmented_model_CFA)

# PART 8: Fat---------------------------------------------------------------------------
# Summarize segmented regression results
summary(segmented_model_Fat)

# Extract and visualize breakpoint
plot(segmented_model_Fat, main = "Segmented Regression for Fat vs Age")
abline(v = segmented_model_Fat$psi[2], col = "red", lty = 2)

# Extract slopes (before and after the breakpoint)
slope(segmented_model_Fat)

# Confidence intervals for breakpoint and slopes
confint(segmented_model_Fat)
anova(lm_init_Fat, segmented_model_Fat)

### -----------------------------### 12. Replacement by ATPase----------------------
# PART 1: Preparation and Correlation Analysis----------------------------------
{
  # Copy original MFD values into a new column
  L$MFD_ATP <- L$MFD
  
  # Indices of rows (subject numbers) to replace with new values
  row_indices <- c(1, 4, 7, 12, 13, 19, 20, 21, 29, 30, 32, 33, 34, 38)
  
  # New ATP-derived MFD values for specified subjects
  new_values <- c(1357.017544, 1116.733379, 
                  644.9190409, 638.5964912, 961.1027678,
                  249.2767811, 350.5436414, 324.6866589, 395.9397818, 308.743266, 
                  362.6611691, 184.5928212, 235.3522303, 90.07487949)
  
  # Replace original MFD values with ATP-derived values
  L$MFD_ATP[row_indices] <- new_values
  MFD_ATP <- L$MFD_ATP
}

# 【Correlation Analysis】再計算したMFD（ATP由来）との相関を再確認
cor.test(ABx, MFD_ATP, method = "spearman") # Spearman相関係数: -0.8665062
cor.test(ABx, MFD_ATP, method = "pearson")  # Pearson相関係数: -0.8280345
cor.test(log(ABx), log(MFD_ATP), method = "pearson") # 対数変換後のPearson相関係数: -0.8680387

# 散布図を描画
qplot(ABx, MFD_ATP) +
  geom_point(size = 2, colour = "black") +
  geom_smooth(span =1, colour = "red", se = FALSE) +
  scale_x_continuous(limits = c(0, 16.5), breaks = 0:16.5) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 40, face = "bold"), 
    axis.title.y = element_text(size = 40, face = "bold"), 
    axis.text.x = element_text(size = 16), 
    axis.text.y = element_text(size = 16)
  )

# PART 2: Box plot analysis-----------------------------------------------------
# MFD（ATP由来）のボックスプロットを作成
{
  dataMFD_ATP <- data.frame(ABx = ABx, MFD_ATP = MFD_ATP)
  
  # Age groups: [1-6), [6-7), [7-11)
  boxplot(MFD_ATP ~ cut(ABx, breaks = c(1, 6, 7, 11), labels = c("[1-6)", "[6-7)", "[7-11)"), right = FALSE), 
          xlab = "ABx", cex.axis = 2.5, ylab = "ATPase (replaced)", cex.lab = 2.5, data = dataMFD_ATP, outline = FALSE)
  stripchart(MFD_ATP ~ cut(ABx, breaks = c(1, 6, 7, 11), labels = c("[1-6)", "[6-7)", "[7-11)"), right = FALSE), 
             data = dataMFD_ATP, method = "jitter", pch = 16, col = "black", vertical = TRUE, add = TRUE)
  
  # Detailed age groups: [1-4), [4-6), [6-7), [7-9), [9-11)
  AgeGroups_detailed_ATP <- cut(ABx, breaks = c(1, 4, 6, 7, 9, 11), labels = c("[1,4)", "[4,6)", "[6,7)", "[7,9)", "[9,11)"), right = FALSE)
  
  boxplot(MFD_ATP ~ AgeGroups_detailed_ATP, xlab = "ABx", cex.axis = 2,
          ylab = "ATPase (replaced)", cex.lab = 2.5, data = dataMFD_ATP, outline = FALSE)
  stripchart(MFD_ATP ~ AgeGroups_detailed_ATP, data = dataMFD_ATP,
             method = "jitter", pch = 16, col = "black", vertical = TRUE, add = TRUE)
  
  # Side-by-side boxplots for simplified age groups
  par(mfrow = c(1, 2), mar = c(5, 5, 2, 2))
  
  # [1-6) vs [6-11)
  boxplot(MFD_ATP ~ cut(ABx, breaks = c(1, 6, 11), labels = c("[1-6)", "[6-11)"), right = FALSE), 
          xlab = "Age Groups", ylab = "MFD (ATPase replaced)", outline = FALSE)
  stripchart(MFD_ATP ~ cut(ABx, breaks = c(1, 6, 11), labels = c("[1-6)", "[6-11)"), right = FALSE), 
             method = "jitter", pch = 16, col = "black", vertical = TRUE, add = TRUE)
  
  # [1-7) vs [7-11)
  boxplot(MFD_ATP ~ cut(ABx, breaks = c(1, 7, 11), labels = c("[1-7)", "[7-11)"), right = FALSE), 
          xlab = "Age Groups", ylab = "MFD (ATPase replaced)", outline = FALSE)
  stripchart(MFD_ATP ~ cut(ABx, breaks = c(1, 7, 11), labels = c("[1-7)", "[7-11)"), right = FALSE), 
             method = "jitter", pch = 16, col = "black", vertical = TRUE, add = TRUE)
  
  par(mfrow = c(1, 1)) # Reset plot area
}

# PART 3: Logistic Regression Analysis------------------------------
# MFD（ATP由来）を使ったロジスティック回帰分析
{
  data_logis_ATP <- subset(data.frame(ABx, MFD_ATP), ABx < 11)
  options(warn = -1) # 警告抑制
  
  data_logis_ATP$AgeBinary_ATP <- ifelse(data_logis_ATP$ABx < 6, 0, 1)
  results_ATP <- data.frame()
  
  thresholds_ATP <- seq(500, 600, by = 1)
  
  for (threshold_ATP in thresholds_ATP) {
    data_logis_ATP$MFDBinary_ATP <- ifelse(data_logis_ATP$MFD_ATP > threshold_ATP, 1, 0)
    model_ATP <- glm(AgeBinary_ATP ~ MFDBinary_ATP, data = data_logis_ATP, family = binomial)
    
    aic_value_ATP <- AIC(model_ATP)
    p_value_ATP <- summary(model_ATP)$coefficients[2, 4]
    odds_ratio_ATP <- exp(coef(model_ATP)[2])
    conf_int_ATP <- exp(confint(model_ATP))[2, ]
    
    predicted_ATP <- ifelse(predict(model_ATP, type = "response") > 0.5, 1, 0)
    accuracy_ATP <- mean(predicted_ATP == data_logis_ATP$AgeBinary_ATP) * 100
    
    results_ATP <- rbind(results_ATP, data.frame(
      Threshold_ATP = threshold_ATP, AIC_ATP = aic_value_ATP, P_Value_ATP = format(p_value_ATP, scientific = FALSE),
      Accuracy_ATP = accuracy_ATP, Odds_Ratio_ATP = ifelse(p_value_ATP < 0.05, odds_ratio_ATP, "NS"),
      CI_Lower_ATP = ifelse(p_value_ATP < 0.05, conf_int_ATP[1], "NS"), CI_Upper_ATP = ifelse(p_value_ATP < 0.05, conf_int_ATP[2], "NS")
    ))
  }
  
  print(results_ATP)
  options(warn = 0)
}


# Logistic regression analysis using ATPase‑based MFD
#
# • Significant thresholds (≈ 546–572 fibers/mm²) were virtually identical
#   to those obtained with the original H&E/G‑T/M‑T dataset, even though
#   a few cut‑offs narrowly missed the strict Bonferroni criterion
#   (p < 0.000495).  This concordance underscores the robustness and
#   reproducibility of the MFD threshold across different staining
#   methods and analytic conditions.
#
# • Thresholds ≥ 573 fibers/mm² delivered near‑100% classification accuracy
#   but failed to reach statistical significance. In this zone, the two
#   age groups are (quasi‑)completely separated, producing spuriously low
#   AIC values and probably rendering logistic coefficients and p‑values
#   non‑informative. Consequently, formal inference should be confined
#   to thresholds ≤ 572 fibers/mm².


# PART 4: Segmented regression analysis-----------------------------
{
  # 【データ準備】ATP由来のMFDデータを使用した新たなデータフレームを作成
  dataMFD_ATP <- data.frame(ABx = ABx, MFD_ATP = MFD_ATP)
  dataSegment_MFD_ATP <- dataMFD_ATP %>%
    filter(ABx >= 1 & ABx <= 11)
  
  # 【線形モデルの初期フィット】単純線形回帰を実施
  lm_init_MFD_ATP <- lm(MFD_ATP ~ ABx, data = dataSegment_MFD_ATP)
  
  # 【分節回帰モデルの構築】初期ブレークポイントを6歳に設定
  segmented_model_MFD_ATP <- segmented(lm_init_MFD_ATP, seg.Z = ~ABx, psi = 6)
}

# 【分節回帰結果の要約】
summary(segmented_model_MFD_ATP) # 調整済み決定係数（Adjusted R-squared）: 0.8223

# 【ブレークポイントの抽出と可視化】
plot(segmented_model_MFD_ATP, main = "Segmented Regression for MFD_ATP vs Age")
abline(v = segmented_model_MFD_ATP$psi[2], col = "red", lty = 2)

# 【ブレークポイント前後の傾きを抽出】
slope(segmented_model_MFD_ATP)

# 【ブレークポイントおよび傾きの信頼区間を計算】
confint(segmented_model_MFD_ATP) # ブレークポイントの信頼区間: 6.42 [5.28036 to 7.55964]
anova(lm_init_MFD_ATP, segmented_model_MFD_ATP) # 有意差検定の結果: Pr(>F) = 0.0002935 ***

# 【予測値とRMSE（残差平方和平方根）の計算】
dataSegment_MFD_ATP$Predicted_MFD_ATP <- predict(segmented_model_MFD_ATP)
RMSE_original_ATP <- sqrt(mean((dataSegment_MFD_ATP$MFD_ATP - dataSegment_MFD_ATP$Predicted_MFD_ATP)^2))
print(RMSE_original_ATP) # RMSE（ATP由来データ）: 121.2629

# ------- Segmented model plot--------------------------------------------------
{
  # 【プロット用のデータ準備】
  bp   <- segmented_model_MFD_ATP$psi[2]            # ブレークポイント（6.42歳）
  xseq <- seq(0, 11, 0.1)
  pred <- predict(segmented_model_MFD_ATP,
                  newdata = data.frame(ABx = xseq),
                  se.fit  = FALSE)
  pred_df <- data.frame(ABx = xseq, fit = pred)
  
  # 【カラー設定】
  pal_age <- gradient_n_pal(c("#3182BD", "#FEE08B", "#D73027"))
  
  col_line_pre  <- "#2AAA5E"   # ブレークポイント前の回帰直線
  col_line_post <- "#666666"   # ブレークポイント後の回帰直線
  col_bp        <- "#D72631"   # ブレークポイントライン
  fill_pre      <- "#E5F4EA"   # 背景色（前半）
  fill_post     <- "#F2F2F2"   # 背景色（後半）
  
  # 【プロット描画】
  ggplot(dataMFD_ATP %>% filter(ABx <= 11), aes(x = ABx, y = MFD_ATP)) +
    # 背景設定
    annotate("rect", xmin = 0, xmax = bp,
             ymin = -Inf, ymax = Inf, fill = fill_pre, alpha = 1) +
    annotate("rect", xmin = bp, xmax = 11,
             ymin = -Inf, ymax = Inf, fill = fill_post, alpha = 1) +
    
    # 散布図（年齢に基づくカラーグラデーション）
    geom_point(aes(colour = ABx), size = 4, stroke = 0.1) +
    scale_colour_gradientn(
      colours = c("#3182BD", "#FEE08B", "#D73027"),
      limits = c(0, 11), name = "Age"
    ) +
    
    # 分節回帰の直線を追加
    geom_line(
      data = pred_df %>% filter(ABx <= bp),
      aes(x = ABx, y = fit),
      colour = col_line_pre, size = 1.4
    ) +
    geom_line(
      data = pred_df %>% filter(ABx >= bp),
      aes(x = ABx, y = fit),
      colour = col_line_post, size = 1.4
    ) +
    
    # ブレークポイントのラインと注釈
    geom_vline(xintercept = bp,
               colour = col_bp, linetype = "dashed", size = 1) +
    annotate("text", x = bp + 0.15,
             y = max(dataMFD_ATP$MFD_ATP[dataMFD_ATP$ABx <= 11]) * 0.92,
             label = sprintf("Breakpoint %.2f y", bp),
             hjust = 0, vjust = 1, angle = 0,
             colour = col_bp, size = 5.2, fontface = "bold") +
    
    # 軸設定
    scale_x_continuous(limits = c(0, 11), breaks = 0:11) +
    scale_y_continuous(
      limits = c(0, max(dataMFD_ATP$MFD_ATP[dataMFD_ATP$ABx <= 11]) * 1.05),
      expand = expansion(mult = c(0, 0.02))
    ) +
    
    # 軸ラベルおよびテーマ設定
    labs(x = "ABx", y = "MFD (ATPase replaced)") +
    theme_bw() +
    theme(
      axis.title = element_text(size = 28, face = "bold.italic"),
      axis.text  = element_text(size = 12),
      legend.position = "none",
      plot.margin = margin(8, 12, 8, 8)
    )
}


# ------- LOOCV analysis----------------------------------------------------------
{
  # 【LOOCV準備】LOOCV結果を格納するベクトルの初期化
  n <- nrow(dataSegment_MFD_ATP)
  breakpoints_ATP <- numeric(n)
  CI_lower_ATP <- numeric(n)
  CI_upper_ATP <- numeric(n)
  
  # 【LOOCV実施】データを1つずつ除外して分節回帰を繰り返す
  for (i in 1:n) {
    
    # 1観測値を除外して新たなデータセットを作成
    data_LOOCV_ATP <- dataSegment_MFD_ATP[-i, ]
    
    # 除外データで線形回帰モデルをフィット
    lm_LOOCV_ATP <- lm(MFD_ATP ~ ABx, data = data_LOOCV_ATP)
    
    # 初期ブレークポイントを6歳に設定して分節回帰モデルをフィット
    segmented_LOOCV_ATP <- tryCatch({
      segmented(lm_LOOCV_ATP, seg.Z = ~ABx, psi = 6)
    }, error = function(e) {
      return(NULL)
    })
    
    # フィット成功時、ブレークポイントと信頼区間を保存
    if (!is.null(segmented_LOOCV_ATP)) {
      breakpoints_ATP[i] <- segmented_LOOCV_ATP$psi[2]
      
      # 信頼区間（下限・上限）を抽出
      CI_vals_ATP <- confint(segmented_LOOCV_ATP)[1, ]
      CI_lower_ATP[i] <- CI_vals_ATP[2]
      CI_upper_ATP[i] <- CI_vals_ATP[3]
    } else {
      # モデル失敗時はNAを設定
      breakpoints_ATP[i] <- NA
      CI_lower_ATP[i] <- NA
      CI_upper_ATP[i] <- NA
    }
  }
  
  # NA値を除去したLOOCV結果
  valid_breakpoints_ATP <- breakpoints_ATP[!is.na(breakpoints_ATP)]
  valid_CI_lower_ATP <- CI_lower_ATP[!is.na(CI_lower_ATP)]
  valid_CI_upper_ATP <- CI_upper_ATP[!is.na(CI_upper_ATP)]
  
  # ブレークポイントおよび信頼区間の平均値・中央値を計算
  mean_breakpoint_ATP <- mean(valid_breakpoints_ATP)
  median_breakpoint_ATP <- median(valid_breakpoints_ATP)
  mean_CI_lower_ATP <- mean(valid_CI_lower_ATP)
  mean_CI_upper_ATP <- mean(valid_CI_upper_ATP)
  
  # 【結果表示】LOOCVで得られたブレークポイントと信頼区間
  cat("LOOCV Mean Breakpoint_ATP:", mean_breakpoint_ATP, "\n")
  cat("LOOCV Median Breakpoint_ATP:", median_breakpoint_ATP, "\n")
  cat("LOOCV Mean CI_ATP:", mean_CI_lower_ATP, "to", mean_CI_upper_ATP, "\n")
  
  # ブレークポイントの範囲を計算して表示
  breakpoint_min_ATP <- min(breakpoints_ATP, na.rm = TRUE)
  breakpoint_max_ATP <- max(breakpoints_ATP, na.rm = TRUE)
  cat("Breakpoint Range_ATP:", breakpoint_min_ATP, "to", breakpoint_max_ATP, "\n")
  
  # オリジナルデータの予測とRMSEを計算
  dataSegment_MFD_ATP$Predicted_MFD_ATP <- predict(segmented_model_MFD_ATP)
  RMSE_original_ATP <- sqrt(mean((dataSegment_MFD_ATP$MFD_ATP - dataSegment_MFD_ATP$Predicted_MFD_ATP)^2))
  cat("Original Data RMSE_ATP:", RMSE_original_ATP, "\n")
  
  # LOOCVでの予測値を取得してRMSEを計算
  LOOCV_predictions_ATP <- numeric(n)
  for (i in 1:n) {
    data_LOOCV_ATP <- dataSegment_MFD_ATP[-i, ]
    lm_LOOCV_ATP <- lm(MFD_ATP ~ ABx, data = data_LOOCV_ATP)
    segmented_LOOCV_ATP <- tryCatch({
      segmented(lm_LOOCV_ATP, seg.Z = ~ABx, psi = 6)
    }, error = function(e) {
      return(NULL)
    })
    
    if (!is.null(segmented_LOOCV_ATP)) {
      LOOCV_predictions_ATP[i] <- predict(segmented_LOOCV_ATP, newdata = dataSegment_MFD_ATP[i, ])
    } else {
      LOOCV_predictions_ATP[i] <- NA
    }
  }
  
  # NA値を除去しRMSEを計算
  valid_indices_ATP <- !is.na(LOOCV_predictions_ATP)
  LOOCV_actuals_ATP <- dataSegment_MFD_ATP$MFD_ATP[valid_indices_ATP]
  LOOCV_preds_ATP <- LOOCV_predictions_ATP[valid_indices_ATP]
  
  RMSE_LOOCV_ATP <- sqrt(mean((LOOCV_actuals_ATP - LOOCV_preds_ATP)^2))
  cat("LOOCV RMSE_ATP:", RMSE_LOOCV_ATP, "\n")
  
  # RMSEの差分および比率を計算して表示
  cat("Difference in RMSE (LOOCV_ATP - Original_ATP):", RMSE_LOOCV_ATP - RMSE_original_ATP, "\n")
  cat("RMSE Ratio_ATP (LOOCV_ATP / Original_ATP):", RMSE_LOOCV_ATP / RMSE_original_ATP, "\n")
  
  # ブレークポイントの分布をヒストグラムで可視化
  hist(valid_breakpoints_ATP,
       main = "LOOCV Breakpoints Distribution_ATP",
       xlab = "Breakpoint Age_ATP",
       col = "lightblue",
       border = "black")
  abline(v = mean_breakpoint_ATP, col = "red", lty = 2, lwd = 2)
  abline(v = median_breakpoint_ATP, col = "blue", lty = 2, lwd = 2)
  legend("topright", legend = c("Mean_ATP", "Median_ATP"),
         col = c("red", "blue"), lty = 2, lwd = 2, bty = "n")
}


# LOOCV Mean Breakpoint_ATP: 6.397694 
# LOOCV Median Breakpoint_ATP: 6.419996 
# LOOCV Mean CI_ATP: 5.231917 to 7.56347 
# Breakpoint Range_ATP: 6.250004 to 6.420017 
# Original Data RMSE_ATP: 121.2629 
# LOOCV RMSE_ATP: 132.3526 
# Difference in RMSE (LOOCV_ATP - Original_ATP): 11.08972 
# RMSE Ratio_ATP (LOOCV_ATP / Original_ATP): 1.091452 

# Our exploratory analysis using ATPase-stained sections,
# which allowed selective analysis of better-preserved tissue areas,
# yielded more stable results with narrower confidence intervals,
# suggesting that methodological adjustments may mitigate the impact of influential outliers.
#
# Notably, in contrast to analyses using original MFD data,
# no influential outliers (such as patient #13 in the original dataset) were detected 
# in ATPase-derived MFD values, highlighting the context-dependent nature of outlier detection.
# This implies that previous outlier-driven variability in breakpoint estimates 
# was not due to inherent unreliability of the MFD metric itself,
# but rather due to methodological conditions (e.g., tissue preservation, staining methods, analytical approach).
#
# Thus, methodological optimization (e.g., selective use of better-preserved tissue regions via ATPase staining)
# can effectively enhance the stability and robustness of breakpoint estimation, 
# reinforcing the overall reliability and sensitivity of the MFD metric.


# PART 5: Overall interpretation on MFD_ATP results-----------------------------

# The additional analyses using MFD derived from ATPase-stained sections (MFD_ATP)
# provided highly consistent and robust results across all statistical methods employed:
#
# 1. Correlation analysis:
#    MFD_ATP maintained a similarly strong negative Spearman correlation with age at biopsy (rho ≈ -0.87),
#    corroborating the original MFD results (rho = -0.85). This reinforces the validity and
#    robustness of MFD as an age-associated pathological marker, irrespective of the staining method.
#
# 2. Box plot visualization:
#    Box plot distributions of MFD_ATP across defined age groups demonstrated nearly complete 
#    separation between patients younger than 6 years and those aged 6 years or older,
#    consistent with the original MFD data. This visually confirms that the observed critical
#    threshold around age 6 is robust across staining conditions and analytical strategies.
#
# 3. Logistic regression analysis:
#    Logistic regression models using MFD_ATP identified threshold ranges of 546–572 fibers/mm² 
#    as statistically significant discriminatory points, precisely replicating the findings from 
#    original MFD analyses. Despite slight variations in Bonferroni-adjusted significance levels,
#    these results strongly reinforce the consistency and biological relevance of this critical threshold.
#
# 4. Segmented regression analysis:
#    Segmented regression analysis of MFD_ATP confirmed a stable breakpoint around 6.40 years
#    (95% CI approximately 5.23 to 7.56), closely aligning with the original MFD breakpoint 
#    estimation (6.25 years). Furthermore, the use of ATPase staining allowed better selection 
#    of well-preserved tissue regions, thereby improving model stability and mitigating the impact
#    of influential outliers seen in original analyses.
#
# These convergent findings across multiple analytical methods clearly demonstrate that MFD, 
# when derived from optimally selected and digitally restored ATPase-stained sections, 
# consistently and robustly identifies the early-stage critical transition around age 6 in DMD patients.
#
# Therefore, concerns regarding potential instability or unreliability of the MFD metric due to 
# staining or sample quality variability are effectively addressed. The comprehensive evidence 
# presented strongly supports the scientific validity, reproducibility, and biological significance 
# of MFD as an early and sensitive biomarker for assessing Duchenne muscular dystrophy progression.


### -----------------------------### 13. Robustness Check----
# PART 1: Bayesian segmented model  (fixed BP = 6.25 y)---------------
# 1: Build the design matrix for Bayesian segmented regression (BP = 6.25 y)
{
  idx <- ABx < 11        
  ABx_sub <- ABx[idx]
  MFD_sub <- MFD[idx]  
  
  bp <- 6.25
  below <- pmin(ABx_sub, bp)           # slope before BP  (β1)
  above <- pmax(ABx_sub - bp, 0)       # slope after BP   (β2)
  
  X <- cbind(1, below, above)      # (Intercept, β1, β2)
  y <-  MFD_sub
  
  n <- nrow(X)
  p <- ncol(X)
  
  # 2: Compute ordinary-least-squares (OLS) estimates for initialization
  XtX_inv  <- solve(t(X) %*% X)
  beta_hat <- XtX_inv %*% t(X) %*% y
  resid    <- y - X %*% beta_hat
  s2_hat   <- sum(resid^2) / (n - p)       # unbiased residual variance
  
  # 3: Sample from conjugate posterior distribution
  set.seed(123)
  draws <- 10000
  
  # sample sigma² from scaled inverse chi-squared distribution
  sigma2 <- (n - p) * s2_hat / rchisq(draws, df = n - p)  # length = draws
  
  # Cholesky of (X'X)^{-1}
  chol_XtX <- chol(XtX_inv)                  # 3 x 3 upper‑triangular
  
  # generate standard normals, then scale
  Z  <- matrix(rnorm(draws * p), nrow = draws, ncol = p)   # draws × 3
  beta_draws <- sweep(Z %*% chol_XtX, 1, sqrt(sigma2), "*")       # scale rows
  beta_draws <- sweep(beta_draws, 2, beta_hat, "+")        # add mean
  
  colnames(beta_draws) <- c("Intercept", "beta1", "beta2")
  slope_pre  <- beta_draws[, "beta1"]   # slope before BP
  slope_post <- beta_draws[, "beta2"]   # slope after BP
  
  # 4:  Define ROPE using measurement error (differences between actual MFD and ATPase-replaced MFD)
  {
    # MFD value copy
    L$MFD_ATP <- L$MFD
    # rows
    row_indices <- c(1, 4, 7, 12, 13, 19, 20, 21, 29, 30, 32, 33, 34, 38) # subject number
    
    # replace values
    new_values <- c(1357.017544, 1116.733379, 
                    644.9190409, 638.5964912, 961.1027678,
                    249.2767811, 350.5436414, 324.6866589, 395.9397818, 308.743266, 
                    362.6611691, 184.5928212, 235.3522303, 90.07487949)
    
    # new parameter
    L$MFD_ATP[row_indices] <- new_values; MFD_ATP <- L$MFD_ATP; MFD_ATP
  }
  
  # subject (row) indices already defined
  row_indices <- c(1, 4, 7, 12, 13, 19, 20, 21, 29, 30, 32, 33, 34, 38)
  
  # vector of differences
  abs_diff <- abs(L$MFD[row_indices] - MFD_ATP[row_indices])
  
  # results
  abs_diff                       # individual differences
  mean_diff <- mean(abs_diff)
  rope <- mean_diff; rope
  
  # 5: Summarize posterior slopes using ROPE criteria 
  cat("ROPE  (±", round(rope,2), " MFD/yr)\n")
  cat("  fraction inside  <BP :",
      round(mean(abs(slope_pre)  < rope) * 100, 2), "%\n")
  cat("  fraction inside  ≥BP :", 
      round(mean(abs(slope_post) < rope) * 100, 2), "%\n\n")
  
  cat("Posterior slope < 6.25 y  :",
      "mean",  round(mean(slope_pre),  2),
      " 95% CrI",
      round(quantile(slope_pre,  c(.025,.975)), 2), "\n")
  
  cat("Posterior slope ≥ 6.25 y  :",
      "mean",  round(mean(slope_post), 2),
      " 95% CrI",
      round(quantile(slope_post, c(.025,.975)), 2), "\n\n")
  
  # 6  Model comparison (BIC → BF approximation)
  fit_seg  <- lm(y ~ below + above)   # segmented
  fit_null <- lm(y ~ ABx_sub)         # single-line (null) model
  
  bic_seg  <- BIC(fit_seg)
  bic_null <- BIC(fit_null)
  bf10     <- exp((bic_null - bic_seg) / 2)
  
  cat("BIC  single‑line :", round(bic_null,2),
      "|  segmented :", round(bic_seg,2),
      "  →  BF10 ≈", format(bf10, digits = 3), "\n")
}
# ROPE  (± 62  MFD/yr)
# fraction inside  <BP : 0 %
# fraction inside  ≥BP : 96.26 %

# Posterior slope < 6.25 y  : mean -170.44  95% CrI -207.47 -133.58 
# Posterior slope ≥ 6.25 y  : mean -28.42  95% CrI -65.22 8.77 
# BIC  single‑line : 470.2 |  segmented : 455.82   →  BF10 ≈ 1324

# 7: Sensitivity analysis of ROPE criterion (varying by SD of measurement error)
{
  abs_diff <- abs(L$MFD[row_indices] - MFD_ATP[row_indices])
  sd_diff   <- sd(abs_diff); sd_diff # 38.48842, SD of absolute differences
  
  rope_grid <- c(rope-sd_diff, rope, rope+sd_diff, rope+1.5*sd_diff, rope+2*sd_diff)   # unit: MFD/year
  
  for (r in rope_grid) {
    cat("±", r, "→  anterior", round(mean(abs(slope_pre)  < r)*100,2), "%,",
        "posterior", round(mean(abs(slope_post) < r)*100,2), "%\n")
  }
}
# ± 23.51481 →  anterior 0 %, posterior 39.23 %
# ± 62.00323 →  anterior 0 %, posterior 96.26 %
# ± 100.4917 →  anterior 0.01 %, posterior 99.98 %
# ± 119.7359 →  anterior 0.34 %, posterior 100 %
# ± 138.9801 →  anterior 4.49 %, posterior 100 %


# ------- Bayesian segmented Plot------------------------------------------------------------------
# Data preparation
{
  sd_pre  <- sd(slope_pre)
  sd_post <- sd(slope_post) 
  
  set.seed(123)
  slope_pre  <- rnorm(10000, mean = -170.4, sd = sd_pre)
  slope_post <- rnorm(10000, mean = -28.4, sd = sd_post)
  rope <- 62  # Region of Practical Equivalence (ROPE)
  
  # Define scale range
  slope_min <- min(c(slope_pre, slope_post, -rope)) - 10
  slope_max <- max(c(slope_pre, slope_post, rope)) + 10
  
  # Create a single data frame for plotting
  posterior_data <- data.frame(
    Slope = c(slope_pre, slope_post),
    Group = factor(rep(c("Before 6.25 yr", "After 6.25 yr"), each = 10000),
                   levels = c("Before 6.25 yr", "After 6.25 yr"))
  )
  
  # Plot
  ggplot(posterior_data, aes(x = Slope, fill = Group)) +
    geom_density(alpha = 0.6) +
    facet_wrap(~ Group, scales = "fixed") +
    scale_x_continuous(limits = c(slope_min, slope_max)) +
    geom_vline(xintercept = c(-rope, rope), linetype = "dashed", color = "red", linewidth = 1) +
    labs(
      x = "Regression Slope (MFD/year)",
      y = "Density"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 14),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 13, face = "bold"),
      strip.text = element_text(size = 12, face = "bold")
    ) +
    scale_fill_manual(values = c("#F8766D", "#00BFC4"))
}

# ROPE ±62 MFD/yr was derived from the mean absolute duplicate error.
# Under this margin, the pre-break slope (≈ −171 MFD/yr) lies entirely outside the ROPE,
# indicating a statistically meaningful decline,
# whereas the post-break slope (≈ −29 MFD/yr) lies almost entirely inside,
# indicating practical flatness.
# Importantly, slopes obtained with the ATP method (−170 / −29 MFD/yr)
# were virtually identical to those from the original staining.
# This consistency confirms that even when considering typical measurement variability,
# the age-dependent trajectory estimation remains robust.

# Sensitivity checks with broader or narrower margins (±24 to 138 MFD/yr)
# produced the same qualitative conclusion, further reinforcing the breakpoint effect.

# PART 2: Bayesian segmented model (brms version)-------------------------------

# Note 1: The following Bayesian segmented regression using brms 
# involves MCMC sampling, which may require significant computational time.

# NOTE 2: This Bayesian analysis using 'brms' and 'Stan' requires Rtools
# to be installed and correctly configured.
# Ensure Rtools (e.g., Rtools43 for R version 4.3.x) is installed and
# its paths are properly set before running this section.
# Refer to: https://cran.r-project.org/bin/windows/Rtools/

{
  # Define ROPE using measurement error (differences between actual MFD and ATPase-replaced MFD)
  # subject (row) indices already defined
  row_indices <- c(1, 4, 7, 12, 13, 19, 20, 21, 29, 30, 32, 33, 34, 38)
  
  # vector of differences
  abs_diff <- abs(L$MFD[row_indices] - MFD_ATP[row_indices])
  
  # results
  abs_diff                       # individual differences
  mean_diff <- mean(abs_diff)
  rope <- mean_diff # 62.00323
  
  {
    mu_x <- mean(ABx_sub)
    df   <- data.frame(x = ABx_sub - mu_x,   
                       y = MFD_sub)
    
    form <- brms::bf(
      y ~ alpha + beta1*x + betaDelta * (x - bp) * step(x - bp),
      nl = TRUE,
      alpha     ~ 1,
      beta1     ~ 1,
      betaDelta ~ 1,
      bp        ~ 1
    )
    
    # Intercept SD = 300 was chosen to match the empirical SD of MFD (~320).
    # Slope SD = 40 was set to be ~2 × the empirical SEs of slope1/2 (≈21–24), giving weakly-informative priors
    # breakpoint: SD = 2 y  →  95 % of the prior mass lies within ±4 y
    # (covers almost the full 1–11 y analysis range (~3.5 × SE 0.58) → weakly-informative)
    
    pri <- c(
      prior(normal(1500,300), nlpar="alpha"),
      prior(normal(-171, 40), nlpar="beta1"),
      prior(normal( 142, 40), nlpar="betaDelta"),
      prior(normal(   0,  2), nlpar="bp")        
    )
    
    
    fit_brm_pw <- brm(
      form, data=df, prior=pri, family=gaussian(),
      chains=4, cores=4, iter=8000, warmup=4000,
      control=list(adapt_delta=0.99, max_treedepth=15)
    )
    
    check_hmc_diagnostics(fit_brm_pw$fit)   # divergence 0
    pp_check(fit_brm_pw, type = "dens_overlay")
    
    draws <- as_draws_df(fit_brm_pw)           
    
    mu_x <- mean(ABx_sub)
    
    draws$bp_raw        <- draws$b_bp_Intercept         + mu_x
    draws$slope_pre     <- draws$b_beta1_Intercept
    draws$slope_post    <- draws$b_beta1_Intercept      + draws$b_betaDelta_Intercept
    draws$diff_slope    <- draws$b_betaDelta_Intercept
    draws$intercept_pre_raw <-
      draws$b_alpha_Intercept -
      draws$b_beta1_Intercept      * mu_x -
      draws$b_betaDelta_Intercept  * pmax(0, -draws$b_bp_Intercept)
    
    pars <- c("intercept_pre_raw","slope_pre","slope_post",
              "diff_slope","bp_raw","sigma")
    
    summ <- t(apply(draws[, pars], 2, function(v) {
      n <- length(v)
      c(mean    = mean(v),
        se_mean = sd(v) / sqrt(n),     # Monte-Carlo
        sd      = sd(v),
        quantile(v, probs = c(.025, .50, .975)))
    }))
    
    colnames(summ) <- c("mean","se_mean","sd","2.5%","50%","97.5%")
    
    print(round(summ, 3))
    prob_slope_pre_in_rope <- mean(abs(draws$slope_pre) < rope)
    prob_slope_post_in_rope <- mean(abs(draws$slope_post) < rope)
    
    cat(" fraction inside (< BP)  :", round(prob_slope_pre_in_rope * 100, 2), "%\n")
    cat(" fraction inside (≥ BP)  :", round(prob_slope_post_in_rope * 100, 2), "%\n")
    
  }
  
  
  # The estimates vary slightly from run to run because the sampler is stochastic.
  #
  #                       mean      sd     2.5%      50%    97.5%  Rhat 
  # intercept_pre_raw 1459.611  81.669 1292.298 1460.278 1617.477  1.00
  # slope_pre         -165.137  17.339 -200.319 -164.869 -131.151  1.00
  # slope_post         -28.723  21.087  -67.841  -29.730   15.941  1.00
  # diff_slope         136.414  24.485   87.635  136.422  184.847   NA   
  # bp_raw               6.373   0.577    5.236    6.367    7.657  1.00
  # sigma              145.528  19.262  113.883  143.601  188.628  1.00
  
  # diff_slope is a deterministic transform of sampled parameters; therefore Rhat is not reported.
  
  # ROPE
  # fraction inside (< BP)  : 0 %
  # fraction inside (≥ BP)  : 95.01 %
  
  # Sensitivity analysis
  abs_diff <- abs(L$MFD[row_indices] - MFD_ATP[row_indices])
  sd_diff   <- sd(abs_diff); sd_diff # 38.48842, SD of absolute differences
  
  rope_grid <- c(rope - sd_diff, rope, rope + sd_diff,
                 rope + 1.5*sd_diff, rope + 2*sd_diff)
  
  for (r in rope_grid) {
    cat("±", round(r,1), "fib/mm² →  anterior",
        round(mean(abs(slope_pre)  < r)*100,2), "% ,",
        "posterior", round(mean(abs(slope_post) < r)*100,2), "%\n")
  }
}
#  ± 23.5 fib/mm² →  anterior 0 % , posterior 38.87 %
#  ± 62 fib/mm² →  anterior 0 % , posterior 95.01 %
#  ± 100.5 fib/mm² →  anterior 0.01 % , posterior 100 %
#  ± 119.7 fib/mm² →  anterior 0.37 % , posterior 100 %
#  ± 139 fib/mm² →  anterior 4.48 % , posterior 100 %

# Divergences:
# 0 of 16000 iterations ended with a divergence.

# Tree depth:
# 0 of 16000 iterations saturated the maximum tree depth of 15.

# Energy:
# E-BFMI indicated no pathological behavior.

# Rhat, Bulk_ESS, Tail_ESS etc.
fit_brm_pw

# plot
pp_check(fit_brm_pw, type = "dens_overlay")


# HDI
{
  hdi_beta1 <- hdi(as.vector(slope_pre), credMass = 0.95)
  print(hdi_beta1)
  
  hdi_beta2 <- hdi(as_vector(slope_post), credMass = 0.95)
  print(hdi_beta2)
  
  hdi_bp <- hdi(draws$bp_raw, credMass = 0.95)
  print(hdi_bp)
}

# slope_pre
# lower     upper 
# -197.4733 -133.1994 

# slope_post
# lower     upper 
# -67.03500  10.09105 

# breakpoint
# lower    upper 
# 5.252417 7.776177 


# PART 3: Segmented MFD cut-off & Monte-Carlo misclassification analysis---------------- 

## 1. Data preparation

{
  dat <- dataAll %>% 
    mutate(
      ABx = as.numeric(ABx),          
      MFD = as.numeric(MFD)           
    ) %>% 
    arrange(ABx)
  options(warn = -1)
  
  ## 2. Restrict analysis to ABx < 11 y
  
  dat_sub <- dat %>% filter(ABx < 11)
  
  
  ## 3. Segmented regression (MFD ~ Age)
  
  lm0 <- lm(MFD ~ ABx, data = dat_sub)          # initial linear model
  seg <- segmented(lm0, seg.Z = ~ABx, psi = 6)  # starting breakpoint = 6 y
  
  # Breakpoint & slopes
  bp     <- seg$psi[2]                 # estimated breakpoint ≈ 6.25 y
  slope1 <- slope(seg)$ABx[1, "Est."]  # slope before breakpoint
  slope2 <- slope(seg)$ABx[2, "Est."]  # slope after  breakpoint
  
  
  ## 4. Error metrics
  
  # 4-A) RMSE on fitted data
  rmse_fit <- sqrt(mean(residuals(seg)^2))
  
  # 4-B) Leave-One-Out CV RMSE
  get_pred_err <- function(i) {
    train <- dat_sub[-i, ]
    test  <- dat_sub[i , , drop = FALSE]
    m0 <- lm(MFD ~ ABx, data = train)
    s  <- segmented(m0, seg.Z = ~ABx, psi = bp)   # use bp as starting value
    pred <- predict(s, newdata = test)
    (test$MFD - pred)^2                          # squared error
  }
  err_vec <- map_dbl(seq_len(nrow(dat_sub)), get_pred_err)
  rmse_cv <- sqrt(mean(err_vec))
  
  # 4-C) Convert RMSE (vertical error) to “age-equivalent” error (horizontal)
  age_err_pre  <- rmse_cv / abs(slope1)
  age_err_post <- rmse_cv / abs(slope2)
  
  
  ## 5. Monte-Carlo simulation of mis-classification error
  age_break <- 6.25 # reference point aligned to the segmented regression breakpoint
  dat_log2 <- dataAll %>%             
    filter(ABx < 11) %>%               
    mutate(AgeBin = ifelse(ABx < age_break, 1, 0))
  
  # θ(MFD) = P(Age < 6.25 y | MFD)
  fit2  <- glm(AgeBin ~ MFD, data = dat_log2, family = binomial)
  
  ## P = 0.50（odds = 1）
  logit <- function(p) log(p/(1-p))
  cut    <- (logit(0.50) - coef(fit2)[1]) / coef(fit2)[2]     # New cut corresponding to 6.25 y
  cut # 510.9353 (consistent with BP = 6.25 y)
  set.seed(123)
  S        <- 2e6                          # number of simulations
  sig      <- rmse_cv                     # vertical SD ≈ RMSE
  DeltaAge <- seq(0.25, 2, by = 0.25)      # ±0.25 – ±2 y
  
  get_mis <- function(dA) {
    dM <- abs(slope1) * dA                 # ΔMFD corresponding to ΔAge
    z  <- rnorm(S, 0, sig)                 # add noise
    mis_low  <- mean((cut + dM + z) <  cut)   # false “<6.25 y”
    mis_high <- mean((cut - dM + z) >= cut)   # false “≥6.25 y”
    c(mis_low, mis_high, (mis_low + mis_high) / 2)
  }
  mc <- t(sapply(DeltaAge, get_mis))
  colnames(mc) <- c("Mis_under_BP", "Mis_over_BP", "MeanErr")
  mcTab <- data.frame(
    Δage = DeltaAge,
    ΔMFD = abs(slope1) * DeltaAge,
    round(mc, 4)
  )
  
  
  ## 6. Logistic regression to derive 80/50/20 % thresholds
  
  dat_log <- dataAll %>% 
    filter(ABx < 11) %>% 
    transmute(
      AgeBin = ifelse(ABx < 6.25, 1, 0),  # 1 = ABx < 6.25 y, 0 = ABx ≥ 6.25 y
      MFD
    )
  
  fit <- glm(AgeBin ~ MFD, data = dat_log, family = binomial)
  logit <- function(p) log(p / (1 - p))
  
  th80 <- (logit(0.80) - coef(fit)[1]) / coef(fit)[2]
  th50 <- (logit(0.50) - coef(fit)[1]) / coef(fit)[2]
  th20 <- (logit(0.20) - coef(fit)[1]) / coef(fit)[2]
  
  thr_vec <- c(hi80 = unname(th80),
               mid50 = unname(th50),
               lo20 = unname(th20))
  
  
  ## 7. Categorise MFD into High / Grey / Low zones
  
  dat_zone <- dat_log %>% 
    mutate(
      MFD_zone = case_when(
        MFD >= thr_vec["hi80"] ~ "High (≥80% <BP)",
        MFD <  thr_vec["lo20"] ~ "Low  (≥80% ≥BP)",
        TRUE                   ~ "Grey (≈50%)"
      )
    )
  
  
  ## 8. Zone-wise summary table
  
  zone_tbl <- dat_zone %>% 
    count(MFD_zone, name = "Cases") %>% 
    mutate(
      Prob_ltBP = case_when(
        MFD_zone == "High (≥80% <BP)" ~ "≥80%",
        MFD_zone == "Grey (≈50%)"     ~ "20–80%",
        MFD_zone == "Low  (≥80% ≥BP)" ~ "≤20%"
      ),
      Window = case_when(
        MFD_zone == "High (≥80% <BP)" ~ "Large therapeutic window",
        MFD_zone == "Grey (≈50%)"     ~ "Requires individual judgement",
        MFD_zone == "Low  (≥80% ≥BP)" ~ "Small therapeutic window"
      ),
      MFD_zone = factor(
        MFD_zone,
        levels = c("High (≥80% <BP)", "Grey (≈50%)", "Low  (≥80% ≥BP)")
      )
    ) %>% 
    arrange(MFD_zone) %>% 
    dplyr::select(
      `MFD zone`                = MFD_zone,
      Cases,
      `Prob(<BP)`               = Prob_ltBP,
      `Presumed therapeutic window`      = Window
    )
  
  
  ## 9. Output
  
  {
    cat("\n--- Recheck: Segmented regression (ABx < 11 y) ---\n")
    cat(sprintf("  Breakpoint (BP)     : %.3f y\n", bp))
    cat(sprintf("  Slope < BP          : %.2f | Slope ≥ BP : %.2f\n", slope1, slope2))
    cat(sprintf("  RMSE (fit)          : %.2f | RMSE (LOOCV) : %.2f\n", rmse_fit, rmse_cv))
    cat(sprintf("  Age-error estimate  : %.2f y (pre) | %.2f y (post)\n",
                age_err_pre, age_err_post))
    cat("\n--- Monte-Carlo mis-classification table ---\n")
    print(mcTab, row.names = FALSE)
    
    cat("\n--- Thresholds derived from logistic model (reference point = 6.25 yr) ---\n")
    cat(sprintf("  High  zone (≈80%% <BP) :  MFD ≥ %.2f fibers/mm²\n", thr_vec["hi80"]))
    cat(sprintf("  Grey  zone (≈50%%)      :  MFD  = %.2f – %.2f fibers/mm²\n",
                thr_vec["lo20"], thr_vec["hi80"]))
    cat(sprintf("  Low   zone (≈80%% ≥BP) :  MFD ≤ %.2f fibers/mm²\n\n", thr_vec["lo20"]))
    
    cat("\n--- Zone summary table ---\n")
    print(zone_tbl, row.names = FALSE)
    options(warn = 0)
    }
}

# --- Recheck: Segmented regression (ABx < 11 y) ---
# Breakpoint (BP): 6.250 y
# Slope < BP: -170.68 | Slope ≥ BP : -28.31
# RMSE (fit): 132.90 | RMSE (LOOCV) : 151.44
# Age-error estimate: 0.89 y (pre) | 5.35 y (post)

# --- Monte-Carlo mis-classification table ---
#   Δage  ΔMFD   Mis_under_BP  Mis_over_BP  MeanErr
#   0.25  42.67       0.3896      0.3887    0.3892
#   0.50  85.34       0.2867      0.2865    0.2866
#   0.75 128.01       0.1993      0.1985    0.1989
#   1.00 170.68       0.1297      0.1301    0.1299
#   1.25 213.35       0.0793      0.0797    0.0795
#   1.50 256.02       0.0456      0.0456    0.0456
#   1.75 298.69       0.0243      0.0242    0.0243
#   2.00 341.36       0.0121      0.0122    0.0121

# --- Thresholds derived from logistic model (reference point = 6.25 yr) ---
# High  zone (≈80% <BP) :  MFD ≥ 595.51 fibers/mm²
# Grey  zone (≈50%)     :  MFD  = 426.36 – 595.51 fibers/mm²
# Low   zone (≈80% ≥BP) :  MFD ≤ 426.36 fibers/mm²


# --- Zone summary table ---
#      MFD zone      Cases  Prob(<BP)   Presumed therapeutic window
#   High (≥80% <BP)    10      ≥80%      Large therapeutic window
#        Grey (≈50%)   11    20–80%    Requires individual judgement
#   Low  (≥80% ≥BP)    14      ≤20%      Small therapeutic window


## Interpretation
# 1) Segmented regression
#    • Model LOOCV RMSE = 151 is adopted as the empirical SD of
#      vertical measurement error.
#    • ΔMFD = |slope_pre| × ΔAge converts a horizontal age shift into the
#      corresponding vertical shift in MFD on the steep limb.

# 2) Monte-Carlo simulation
#    • Vertical SD is projected to the horizontal (age) axis via the steep
#      pre-breakpoint slope.
#    • For each ΔAge (0.25–2 y) we simulate 2 M trials and count how often a
#      single MFD reading crosses the optimal cut-off.
#    • Outcome: at ±1.25 y the mean mis-classification rate drops to ~8 %,

# 3) Logistic cut-offs
#    • Posterior probability thresholds P(Age ≥ 6.25 y | MFD) of 0.20 / 0.50 / 0.80 yield
#      estimated MFD cut-offs of 426.4 / 595.5 fibers/mm².
#    • “High” (≥80% probability Age < 6.25 y), “Grey” (≈50%), and “Low” (≥80% probability Age ≥ 6.25 y)
#      zones are presumed categories translating these probabilities into
#      a qualitative presumed therapeutic window.

# 4) Practical takeaway
#    • Outside ±1.25 y of the 6.25-year mark, the single-cut strategy is highly
#      reliable.
#    • Inside that band, rely on the probabilistic zones plus clinical context
#      rather than a strict binary rule.


# PART 4: Sample size estimation------------------------------------------------
# NOTE: requires significant computational time.

{
  library(parallel)
  library(pbapply)
  library(data.table)
  library(dplyr)
  library(segmented)
  
  # 1) Posterior draws (10000 rows) 
  set.seed(123)
  post_draws <- draws %>%
    dplyr::select(
      beta1   = tidyselect::matches("b_beta1"),
      betaDel = tidyselect::matches("b_betaDelta"),
      sigma   = tidyselect::matches("^sigma$"),
      bp      = tidyselect::matches("b_bp")
    ) %>%
    slice_sample(n = 10000)                   # posterior subsample
  
  # 2) Data subset and centering constant 
  dat_sub <- dataAll %>% filter(ABx < 11)
  mu_x    <- mean(dat_sub$ABx)
  
  # 3) Single simulation 
  one_run <- function(N, tol_age) {
    p      <- post_draws[sample.int(nrow(post_draws), 1), ]
    beta1  <- p$beta1
    beta2  <- p$beta1 + p$betaDel
    bp     <- p$bp
    sig    <- p$sigma
    
    x_raw  <- sort(sample(dat_sub$ABx, N, TRUE))
    x      <- x_raw - mu_x                 # centred age
    
    y <- ifelse(x < bp,
                beta1 * (x - bp) + rnorm(N, 0, sig),
                beta2 * (x - bp) + rnorm(N, 0, sig))
    
    seg_fit <- try(
      segmented(lm(y ~ x), seg.Z = ~x, psi = bp,
                control = seg.control(display = FALSE)),
      silent = TRUE)
    
    if (inherits(seg_fit, "segmented") &&
        abs(seg_fit$psi[2] - bp) <= tol_age) 1 else 0
  }
  
  # 4) Wrapper with parallel execution 
  design_power_mc <- function(N, tol_age,
                              iter    = 10000,
                              ncore   = max(1L, detectCores() - 1L)) {
    cl <- makeCluster(ncore, outfile = "")
    on.exit(try(stopCluster(cl), silent = TRUE), add = TRUE)
    
    clusterCall(cl, function() { suppressMessages(library(segmented)); NULL })
    clusterExport(cl,
                  c("dat_sub", "post_draws", "one_run",
                    "tol_age", "N", "mu_x"),
                  envir = environment())
    
    mean(parSapply(cl, seq_len(iter),
                   function(i) one_run(N, tol_age)))
  }
  
  # 5) Grid and progress-bar execution 
  sizes        <- seq(20, 70, 5)
  tol_vec      <- c(1, 1.25, 1.5)
  grid         <- expand.grid(N = sizes, tol = tol_vec, KEEP.OUT.ATTRS = FALSE)
  
  pboptions(type = "timer")
  res_vec <- pblapply(seq_len(nrow(grid)), function(i) {
    with(grid[i, ], design_power_mc(N, tol_age = tol))
  })
  
  # 6) Assemble results 
  out_dt <- data.table(grid, power = unlist(res_vec))
  result_wide <- dcast(out_dt, N ~ tol, value.var = "power")
  setnames(result_wide,
           c("N", "1", "1.25", "1.5"),
           c("N", "Prec_1yo", "Prec_1.25yo", "Prec_1.5yo"))
  
  print(result_wide, digits = 5)
}


#  N Prec_1yo Prec_1.25yo Prec_1.5yo
# 20   0.5818      0.6543     0.7266
# 25   0.6565      0.7101     0.7677
# 30   0.6980      0.7643     0.8208
# 35   0.7243      0.7947     0.8454 : Our study
# 40   0.7511      0.8164     0.8644
# 45   0.7700      0.8408     0.8816
# 50   0.7971      0.8614     0.8996
# 55   0.8159      0.8763     0.9158
# 60   0.8342      0.8901     0.9240
# 65   0.8451      0.8981     0.9387
# 70   0.8561      0.9065     0.9401


#   N Prec_0.25yo Prec_0.5yo Prec_0.75yo
# 100      0.4265     0.6965      0.8498
# 110      0.4548     0.7204      0.8639
# 120      0.4794     0.7350      0.8658
# 130      0.4882     0.7563      0.8807
# 140      0.4980     0.7665      0.8905
# 150      0.5230     0.7831      0.9001
# 160      0.5311     0.8037      0.9137
# 170      0.5432     0.8110      0.9216
# 180      0.5585     0.8139      0.9311
# 190      0.5684     0.8287      0.9307
# 200      0.5768     0.8410      0.9408
# 210      0.5865     0.8535      0.9410
# 220      0.6071     0.8581      0.9444
# 230      0.5992     0.8644      0.9510
# 240      0.6169     0.8707      0.9535
# 250      0.6272     0.8758      0.9575
# 260      0.6295     0.8809      0.9596
# 270      0.6448     0.8875      0.9627
# 280      0.6451     0.8971      0.9683
# 290      0.6475     0.8995      0.9661
# 300      0.6613     0.8949      0.9712


###------------------------------### 14. Revised Analysis of Variance-----------------------
# PART 1: Age grouping----------------------------------------------------------
AgeGroups_rev <- cut(ABx,
                     breaks = c(1, 5, 7.5, 11),
                     labels = c("[1,5)", "[5,7.5)", "[7.5,11)"),
                     right = FALSE)

# PART 2: Mean analysis  -----------------------------------------------------
dataMean_rev <- data.frame(ABx         = ABx,
                           Mean        = Mean,
                           AgeGroups_rev = AgeGroups_rev)

boxplot(Mean ~ AgeGroups_rev,
        data      = dataMean_rev,
        xlab      = "",          
        ylab      = "",
        cex.axis  = 2.5,
        cex.lab   = 2.5,
        main      = "",
        outline   = FALSE)

stripchart(Mean ~ AgeGroups_rev,
           data      = dataMean_rev,
           method    = "jitter",
           pch       = 16,
           col       = "black",
           vertical  = TRUE,
           add       = TRUE)

Mean_Group1_rev <- dataMean_rev %>% filter(ABx < 5)              %>% pull(Mean)
Mean_Group2_rev <- dataMean_rev %>% filter(ABx >= 5  & ABx < 7.5) %>% pull(Mean)
Mean_Group3_rev <- dataMean_rev %>% filter(ABx >= 7.5 & ABx < 11) %>% pull(Mean)

mean(Mean_Group1_rev) # 696.4769
sd(Mean_Group1_rev) # 178.7592

mean(Mean_Group2_rev) # 1100.194
sd(Mean_Group2_rev) # 405.3749

mean(Mean_Group3_rev) # 1317.226
sd(Mean_Group3_rev) # 819.061

shapiro.test(Mean_Group1_rev) # 0.4544
shapiro.test(Mean_Group2_rev) # 0.2054
shapiro.test(Mean_Group3_rev) # 0.0002137

leveneTest(Mean ~ AgeGroups_rev, data = dataMean_rev) # 0.2433

kruskal.test(Mean ~ AgeGroups_rev, data = dataMean_rev)
# p-value = 0.00419

conover_result_Mean_rev <- kwAllPairsConoverTest(
  Mean ~ AgeGroups_rev,
  data = dataMean_rev,
  p.adjust.method = "bonferroni")
print(conover_result_Mean_rev)

stats_Mean_rev <- tapply(dataMean_rev$Mean,
                         dataMean_rev$AgeGroups_rev,
                         summary)
print(stats_Mean_rev)

## ─────────────────────────────────────────────
## Figure: Mean vs AgeGroups_rev
## ─────────────────────────────────────────────
{
  alpha  <- 0.05
  grp_lv <- levels(dataMean_rev$AgeGroups_rev)      # c("[1,5)", "[5,7.5)", "[7.5,11)")
  p_mat  <- conover_result_Mean_rev$p.value         
  
  sig_pairs <- list()      
  ctr <- 1
  for (k in 1:(length(grp_lv) - 1)) {
    for (l in (k + 1):length(grp_lv)) {
      g1 <- grp_lv[k]; g2 <- grp_lv[l]
      if (g1 %in% rownames(p_mat) && g2 %in% colnames(p_mat)) {
        p_val <- p_mat[g1, g2]
      } else if (g2 %in% rownames(p_mat) && g1 %in% colnames(p_mat)) {
        p_val <- p_mat[g2, g1]
      } else {
        next                              
      }
      if (!is.na(p_val) && p_val < alpha) {
        sig_pairs[[ctr]] <- c(k, l)       
        ctr <- ctr + 1
      }
    }
  }
  
  bp0   <- boxplot(Mean ~ AgeGroups_rev, data = dataMean_rev, plot = FALSE)
  y_max <- max(bp0$stats[5, ])
  y_lim <- c(min(dataMean_rev$Mean), y_max * 1.15)
  
  par(mar = c(5, 6, 4, 2) + 0.1)
  bp <- boxplot(Mean ~ AgeGroups_rev, data = dataMean_rev,
                ylim = y_lim,
                xlab = "", ylab = "",
                cex.axis = 2.2, cex.lab = 2.2,
                outline = FALSE)
  stripchart(Mean ~ AgeGroups_rev, data = dataMean_rev,
             method = "jitter", pch = 16, col = "black",
             vertical = TRUE, add = TRUE)
  
  mtext(side = 2,
        text = expression(bold(italic("Mean"))),
        line = 4,
        cex  = 2.2)
  
  step   <- 0.05 * y_max
  height <- y_max + step
  for (pair in sig_pairs) {
    i <- pair[1]; j <- pair[2]
    segments(i, height, j, height, lwd = 1.5)
    text((i + j) / 2, height * 1.02, "*", cex = 2.8)
    height <- height + step
  }
  
  mtext(side = 1,
        text = "Age Groups (years)",
        line = 3,
        cex  = 2.0)
}

# PART 3: Cov----------------------------------------------------------
dataCov_rev <- data.frame(ABx          = ABx,
                          Cov          = Cov,
                          AgeGroups_rev = AgeGroups_rev)

Cov_Group1_rev <- dataCov_rev %>% filter(ABx < 5)              %>% pull(Cov)
Cov_Group2_rev <- dataCov_rev %>% filter(ABx >= 5  & ABx < 7.5) %>% pull(Cov)
Cov_Group3_rev <- dataCov_rev %>% filter(ABx >= 7.5 & ABx < 11) %>% pull(Cov)

lapply(list(Cov_Group1_rev, Cov_Group2_rev, Cov_Group3_rev), shapiro.test)
# p-value = 2.367e-05*
# p-value = 0.9867
# p-value = 0.014*

leveneTest(Cov ~ AgeGroups_rev, data = dataCov_rev) # 0.5162

# Kruskal-Wallis
kruskal.test(Cov ~ AgeGroups_rev, data = dataCov_rev)
# p-value = 0.002901

conover_result_Cov_rev <- kwAllPairsConoverTest(
  Cov ~ AgeGroups_rev,
  data = dataCov_rev,
  p.adjust.method = "bonferroni")
print(conover_result_Cov_rev)

## ─────────────────────────────────────────────
## Figure: Cov vs AgeGroups_rev
## ─────────────────────────────────────────────）
{
  alpha  <- 0.05
  grp_lv <- levels(dataCov_rev$AgeGroups_rev)          
  p_mat  <- conover_result_Cov_rev$p.value             
  sig_pairs <- list()                                  
  ctr <- 1
  for (i in 1:(length(grp_lv)-1)) {
    for (j in (i+1):length(grp_lv)) {
      g1 <- grp_lv[i]; g2 <- grp_lv[j]
      if (g1 %in% rownames(p_mat) && g2 %in% colnames(p_mat)) {
        p_val <- p_mat[g1, g2]
      } else if (g2 %in% rownames(p_mat) && g1 %in% colnames(p_mat)) {
        p_val <- p_mat[g2, g1]
      } else next
      if (!is.na(p_val) && p_val < alpha) {
        sig_pairs[[ctr]] <- c(i, j)
        ctr <- ctr + 1
      }
    }
  }
  
  bp0   <- boxplot(Cov ~ AgeGroups_rev, data = dataCov_rev, plot = FALSE)
  y_max <- max(bp0$stats[5, ])          
  y_lim <- c(min(dataCov_rev$Cov), y_max * 1.15)
  
  par(mar = c(5, 6, 4, 2) + 0.1)
  
  bp_cov <- boxplot(Cov ~ AgeGroups_rev,
                    data   = dataCov_rev,
                    plot   = FALSE)        
  
  out_vals_cov <- bp_cov$out               
  
  boxplot(Cov ~ AgeGroups_rev,
          data     = dataCov_rev,
          ylim     = y_lim,                
          xlab     = "", ylab = "",
          cex.axis = 2.2, cex.lab = 2.2,
          outline  = FALSE,                
          col      = "lightgrey")
  
  keep_idx_cov <- !(dataCov_rev$Cov %in% out_vals_cov)
  
  stripchart(Cov ~ AgeGroups_rev,
             data      = dataCov_rev[keep_idx_cov, ],
             method    = "jitter",
             pch       = 16,
             col       = "black",
             vertical  = TRUE,
             add       = TRUE)
  
  mtext(side = 2,
        text = expression(bold(italic("Cov"))),
        line = 4,
        cex  = 2.2)
  
  if (length(sig_pairs) > 0) {
    step   <- 0.05 * y_max
    height <- y_max + step
    for (pair in sig_pairs) {
      i <- pair[1]; j <- pair[2]
      segments(i, height, j, height, lwd = 1.5)
      text((i + j) / 2, height * 1.02, "*", cex = 2.8)
      height <- height + step
    }
  }
  
  mtext(side = 1,
        text = "Age group (years)",
        line = 3,
        cex  = 2.0)
  
  stats_Cov_rev <- tapply(dataCov_rev$Cov,
                          dataCov_rev$AgeGroups_rev,
                          summary)
  print(stats_Cov_rev)
}

# $`[1,5)`
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.4346  0.5053  0.5311  0.5925  0.5712  1.2006

# PART 4: Sd----------------------------------------------------------
dataSd_rev <- data.frame(ABx          = ABx,
                         Sd           = Sd,
                         AgeGroups_rev = AgeGroups_rev)

Sd_Group1_rev <- dataSd_rev %>% filter(ABx < 5)              %>% pull(Sd)
Sd_Group2_rev <- dataSd_rev %>% filter(ABx >= 5  & ABx < 7.5) %>% pull(Sd)
Sd_Group3_rev <- dataSd_rev %>% filter(ABx >= 7.5 & ABx < 11) %>% pull(Sd)

lapply(list(Sd_Group1_rev, Sd_Group2_rev, Sd_Group3_rev), shapiro.test)
# p-value = 0.4462
# p-value = 0.1277
# p-value = 0.004873*

leveneTest(Sd ~ AgeGroups_rev, data = dataSd_rev) # 0.08133

kruskal.test(Sd ~ AgeGroups_rev, data = dataSd_rev) # p-value = 0.0004755

conover_result_Sd_rev <- kwAllPairsConoverTest(
  Sd ~ AgeGroups_rev,
  data = dataSd_rev,
  p.adjust.method = "bonferroni")
print(conover_result_Sd_rev)

{
  alpha  <- 0.05
  grp_lv <- levels(dataSd_rev$AgeGroups_rev)      
  p_mat  <- conover_result_Sd_rev$p.value         
  sig_pairs <- list(); ctr <- 1
  for (i in 1:(length(grp_lv) - 1)) {
    for (j in (i + 1):length(grp_lv)) {
      g1 <- grp_lv[i]; g2 <- grp_lv[j]
      if (g1 %in% rownames(p_mat) && g2 %in% colnames(p_mat)) {
        p_val <- p_mat[g1, g2]
      } else if (g2 %in% rownames(p_mat) && g1 %in% colnames(p_mat)) {
        p_val <- p_mat[g2, g1]
      } else next
      if (!is.na(p_val) && p_val < alpha) {
        sig_pairs[[ctr]] <- c(i, j); ctr <- ctr + 1
      }
    }
  }
  
  bp0   <- boxplot(Sd ~ AgeGroups_rev, data = dataSd_rev, plot = FALSE)
  y_max <- max(bp0$stats[5, ])          
  y_lim <- c(min(dataSd_rev$Sd), y_max * 1.15)
  
  par(mar = c(5, 6, 4, 2) + 0.1)
  
  boxplot(Sd ~ AgeGroups_rev,
          data     = dataSd_rev,
          ylim     = y_lim,
          xlab     = "", ylab = "",
          cex.axis = 2.2, cex.lab = 2.2,
          outline  = FALSE,
          col      = "lightgrey")
  
  out_vals_sd <- bp0$out
  keep_idx_sd <- !(dataSd_rev$Sd %in% out_vals_sd)
  stripchart(Sd ~ AgeGroups_rev,
             data      = dataSd_rev[keep_idx_sd, ],
             method    = "jitter",
             pch       = 16, col = "black",
             vertical  = TRUE, add = TRUE)
  
  mtext(side = 2,
        text = expression(bold(italic("Sd"))),
        line = 4,
        cex  = 2.2)
  
  if (length(sig_pairs) > 0) {
    step   <- 0.05 * y_max
    height <- y_max + step
    for (pair in sig_pairs) {
      i <- pair[1]; j <- pair[2]
      segments(i, height, j, height, lwd = 1.5)
      text((i + j) / 2, height * 1.02, "*", cex = 2.8)
      height <- height + step
    }
  }
  
  mtext(side = 1,
        text = "Age group (years)",
        line = 3,
        cex  = 2.0)
  
  stats_Sd_rev <- tapply(dataSd_rev$Sd,
                         dataSd_rev$AgeGroups_rev,
                         summary)
  print(stats_Sd_rev)
}

# PART 5: MFD------------------------------
dataMFD_rev <- data.frame(ABx          = ABx,
                          MFD          = MFD,
                          AgeGroups_rev = AgeGroups_rev)

MFD_1 <- dataMFD_rev %>% filter(ABx < 5)              %>% pull(MFD)
MFD_2 <- dataMFD_rev %>% filter(ABx >= 5  & ABx < 7.5) %>% pull(MFD)
MFD_3 <- dataMFD_rev %>% filter(ABx >= 7.5 & ABx < 11) %>% pull(MFD)

lapply(list(MFD_1, MFD_2, MFD_3), shapiro.test)
# p-value = 0.7836
# p-value = 0.06719
# p-value = 0.6483

leveneTest(MFD ~ AgeGroups_rev, data = dataMFD_rev) # 0.0503: it requires caution

# two approaches
# 1. Welch's ANOVA
welch_anova_MFD <- oneway.test(MFD ~ AgeGroups_rev,
                               data = dataMFD_rev,
                               var.equal = FALSE)
print(welch_anova_MFD) # p-value = 3.235e-05

gh_MFD_rev <- gamesHowellTest(MFD ~ AgeGroups_rev, data = dataMFD_rev)
print(gh_MFD_rev)

# 2. ANOVA
aov_MFD_rev <- aov(MFD ~ AgeGroups_rev, data = dataMFD_rev)
summary(aov_MFD_rev) # Pr(>F) = 7.01e-07 

tukey_MFD_rev <- TukeyHSD(aov_MFD_rev)
print(tukey_MFD_rev)


# Myofiber Density (MFD) box-plot with Games–Howell asterisks

{
  dataMFD <- data.frame(
    ABx           = ABx,
    MFD           = MFD,
    AgeGroups_rev = AgeGroups_rev          # "[1,5)","[5,7.5)","[7.5,11)"
  )
  
  aov_MFD <- aov(MFD ~ AgeGroups_rev, data = dataMFD)
  print(summary(aov_MFD))                 
  
  tukey <- TukeyHSD(aov_MFD)
  p_mat <- tukey$AgeGroups_rev[ , "p adj"]  
  grp_lv <- levels(dataMFD$AgeGroups_rev)
  
  alpha <- 0.05
  sig_pairs <- list()
  for (rn in names(p_mat)) {
    if (p_mat[rn] < alpha) {
      g <- strsplit(rn, "-", fixed = TRUE)[[1]]       # c("B", "A")
      idx <- sort(c(match(g[1], grp_lv), match(g[2], grp_lv)))
      if (!any(vapply(sig_pairs, identical, logical(1), idx)))
        sig_pairs[[length(sig_pairs) + 1]] <- idx    
    }
  }
  
  bp0   <- boxplot(MFD ~ AgeGroups_rev, data = dataMFD, plot = FALSE)
  y_max <- max(bp0$stats[5, ])
  y_lim <- c(min(dataMFD$MFD), y_max * 1.25)
  
  par(mar = c(5, 6, 4, 2) + 0.1)
  boxplot(MFD ~ AgeGroups_rev, data = dataMFD,
          outline = FALSE, col = "lightgrey",
          ylim = y_lim, xlab = "", ylab = "",
          cex.axis = 2.2, cex.lab = 2.2)
  
  stripchart(MFD ~ AgeGroups_rev,
             data = subset(dataMFD, !(MFD %in% bp0$out)),
             method = "jitter", pch = 16, col = "black",
             vertical = TRUE, add = TRUE)
  
  mtext(side = 2, text = expression(bold(italic("MFD"))),
        line = 4, cex = 2.2)
  
  if (length(sig_pairs) > 0) {
    step   <- 0.05 * y_max
    height <- y_max + step
    for (pair in sig_pairs) {
      segments(pair[1], height, pair[2], height, lwd = 1.5)
      text(mean(pair), height * 1.02, "*", cex = 2.8)
      height <- height + step
    }
  }
  
  mtext(side = 1, text = "Age group (years)",
        line = 3, cex = 2.0)
  
  print(round(tukey$AgeGroups_rev, 4))
}

# PART 6: MFA----------------------------------------------------------

dataMFA_rev <- data.frame(
  ABx           = ABx,
  MFA           = MFA,
  AgeGroups_rev = AgeGroups_rev   
)

MFA_g1 <- dataMFA_rev %>% filter(ABx < 5)              %>% pull(MFA)
MFA_g2 <- dataMFA_rev %>% filter(ABx >= 5 & ABx < 7.5)  %>% pull(MFA)
MFA_g3 <- dataMFA_rev %>% filter(ABx >= 7.5 & ABx < 11) %>% pull(MFA)

mean(MFA_g1) # 0.6117417
sd(MFA_g1)# 0.0890273

mean(MFA_g2)# 0.4939644
sd(MFA_g2)# 0.1228538

mean(MFA_g3)# 0.4243013
sd(MFA_g3)# 0.1646846

shapiro.test(MFA_g1) # p-value = 0.3223
shapiro.test(MFA_g2) # p-value = 0.5591
shapiro.test(MFA_g3) # p-value = 0.3185

leveneTest(MFA ~ AgeGroups_rev, data = dataMFA_rev) # 0.08111

aov_MFA_rev <- aov(MFA ~ AgeGroups_rev, data = dataMFA_rev)
summary(aov_MFA_rev) # Pr(>F) = 0.00631 

tukey_MFA_rev <- TukeyHSD(aov_MFA_rev)
print(tukey_MFA_rev)


# Plot
{
  alpha <- 0.05
  tk_mat <- tukey_MFA_rev$AgeGroups_rev           
  sig_pairs <- list()
  
  for (k in seq_len(nrow(tk_mat))) {
    if (tk_mat[k, "p adj"] < alpha) {
      lhs_rhs <- strsplit(rownames(tk_mat)[k], "-", fixed = TRUE)[[1]]
      g1 <- match(lhs_rhs[2], levels(dataMFA_rev$AgeGroups_rev))  # A
      g2 <- match(lhs_rhs[1], levels(dataMFA_rev$AgeGroups_rev))  # B
      sig_pairs[[length(sig_pairs)+1]] <- sort(c(g1, g2))        
    }
  }
  
  bp0   <- boxplot(MFA ~ AgeGroups_rev, data = dataMFA_rev, plot = FALSE)
  y_max <- max(bp0$stats[5, ])
  y_lim <- c(min(dataMFA_rev$MFA), y_max * 1.25)
  
  par(mar = c(5, 6, 4, 2) + 0.1)
  
  boxplot(MFA ~ AgeGroups_rev,
          data     = dataMFA_rev,
          outline  = FALSE,
          col      = "lightgrey",
          ylim     = y_lim,
          xlab     = "", ylab = "",
          cex.axis = 2.2, cex.lab = 2.2)
  
  stripchart(MFA ~ AgeGroups_rev,
             data      = subset(dataMFA_rev, !(MFA %in% bp0$out)),
             method    = "jitter",
             pch       = 16, col = "black",
             vertical  = TRUE, add = TRUE)
  
  mtext(side = 2,
        text = expression(bold(italic("MFA"))),
        line = 4, cex = 2.2)
  
  if (length(sig_pairs) > 0) {
    height <- y_max + 0.05 * y_max         
    for (pair in sig_pairs) {
      i <- pair[1]; j <- pair[2]           # i < j
      segments(i, height, j, height, lwd = 1.5)
      text((i + j)/2, height * 1.02, "*", cex = 2.8)
    }
  }
  
  mtext(side = 1,
        text = "Age group (years)",
        line = 3, cex = 2.0)
}

# PART 7: CFA------------------------------------------------------------
dataCFA_rev <- data.frame(
  ABx           = ABx,
  CFA           = CFA,
  AgeGroups_rev = AgeGroups_rev     # "[1,5)" "[5,7.5)" "[7.5,11)"
)

{
  CFA_Group1_rev <- dataCFA_rev %>%
    filter(ABx < 5) %>%
    pull(CFA)
  
  CFA_Group2_rev <- dataCFA_rev %>%
    filter(ABx >= 5 & ABx < 7.5) %>%
    pull(CFA)
  
  CFA_Group3_rev <- dataCFA_rev %>%
    filter(ABx >= 7.5 & ABx < 11) %>%
    pull(CFA)
}

mean(CFA_Group1_rev)  # 0.3611023
sd(CFA_Group1_rev) # 0.07727264

mean(CFA_Group2_rev) # 0.4630117
sd(CFA_Group2_rev) # 0.1139217

mean(CFA_Group3_rev) # 0.4929692
sd(CFA_Group3_rev) # 0.1541594

shapiro.test(CFA_Group1_rev) # p-value = 0.6639
shapiro.test(CFA_Group2_rev) # p-value = 0.6488
shapiro.test(CFA_Group3_rev) # p-value = 0.3168

leveneTest(CFA ~ AgeGroups_rev, data = dataCFA_rev) # Pr(>F) = 0.09644

aov_CFA_rev <- aov(CFA ~ AgeGroups_rev, data = dataCFA_rev)
summary(aov_CFA_rev) # Pr(>F) = 0.00631 

tukey_CFA_rev <- TukeyHSD(aov_CFA_rev)
print(tukey_CFA_rev)


# plot
{
  bp0   <- boxplot(CFA ~ AgeGroups_rev, data = dataCFA_rev, plot = FALSE)
  y_max <- max(bp0$stats[5, ])             
  y_lim <- c(min(dataCFA_rev$CFA), y_max * 1.20)  
  
  par(mar = c(5, 6, 4, 2) + 0.1)           
  
  boxplot(CFA ~ AgeGroups_rev,
          data     = dataCFA_rev,
          outline  = FALSE,
          col      = "lightgrey",
          ylim     = y_lim,
          xlab     = "",
          ylab     = "",
          cex.axis = 2.5,
          cex.lab  = 2.5)
  
  stripchart(CFA ~ AgeGroups_rev,
             data      = subset(dataCFA_rev, !(CFA %in% bp0$out)),
             method    = "jitter",
             pch       = 16,
             col       = "black",
             vertical  = TRUE,
             add       = TRUE)
  
  mtext(side = 2,
        text = expression(bold(italic("CFA"))),
        line = 4,
        cex  = 2.5)
  
  height <- y_max + 0.05 * y_max           
  segments(1, height, 3, height, lwd = 1.5)
  text(2, height * 1.02, "*", cex = 3)
  
  mtext(side = 1,
        text = "Age group (years)",
        line = 3,
        cex  = 2.5)
}


# PART 8: Fat --------------------------------------------------------------------------
dataFat_rev <- data.frame(
  ABx           = ABx,
  Fat           = Fat,
  AgeGroups_rev = AgeGroups_rev     # "[1,5)" "[5,7.5)" "[7.5,11)"
)

{
  Fat_Group1_rev <- dataFat_rev %>%
    filter(ABx < 5) %>%
    pull(Fat)
  
  Fat_Group2_rev <- dataFat_rev %>%
    filter(ABx >= 5 & ABx < 7.5) %>%
    pull(Fat)
  
  Fat_Group3_rev <- dataFat_rev %>%
    filter(ABx >= 7.5 & ABx < 11) %>%
    pull(Fat)
}

mean(Fat_Group1_rev) # 0.01153532
sd(Fat_Group1_rev) # 0.01504892

mean(Fat_Group2_rev) # 0.02657025
sd(Fat_Group2_rev) # 0.03444286

mean(Fat_Group3_rev) # 0.06436508
sd(Fat_Group3_rev) # 0.04927252

shapiro.test(Fat_Group1_rev) # p-value = 0.0003974
shapiro.test(Fat_Group2_rev) # p-value = 0.0002673
shapiro.test(Fat_Group3_rev) # p-value = 0.08485

car::leveneTest(Fat ~ AgeGroups_rev, data = dataFat_rev) # Pr(>F) = 0.1715

kruskal.test(Fat ~ AgeGroups_rev, data = dataFat_rev) # p-value = 0.004226

conover_Fat_rev <- PMCMRplus::kwAllPairsConoverTest(
  Fat ~ AgeGroups_rev,
  data = dataFat_rev,
  p.adjust.method = "bonferroni")
print(conover_Fat_rev)

stats_Fat_rev <- tapply(dataFat_rev$Fat,
                        dataFat_rev$AgeGroups_rev,
                        summary)
print(stats_Fat_rev)


{
  bp0   <- boxplot(Fat ~ AgeGroups_rev, data = dataFat_rev, plot = FALSE)
  y_max <- max(bp0$stats[5, ])                       
  y_lim <- c(min(dataFat_rev$Fat), y_max * 1.20)     
  
  boxplot(Fat ~ AgeGroups_rev,
          data     = dataFat_rev,
          outline  = FALSE,
          col      = "lightgrey",
          ylim     = y_lim,
          xlab     = "",
          ylab     = "",
          cex.axis = 2.5,
          cex.lab  = 2.5)
  
  stripchart(Fat ~ AgeGroups_rev,
             data      = subset(dataFat_rev, !(Fat %in% bp0$out)),
             method    = "jitter",
             pch       = 16,
             col       = "black",
             vertical  = TRUE,
             add       = TRUE)
  
  mtext(side = 2,
        text = expression(bold(italic("Fat"))),
        line = 4,
        cex  = 2.5)
  
  p_mat <- conover_Fat_rev$p.value
  lv    <- levels(dataFat_rev$AgeGroups_rev)
  
  sig_pairs <- list()
  for (r in rownames(p_mat)) {
    for (c in colnames(p_mat)) {
      p <- p_mat[r, c]
      if (!is.na(p) && p < 0.05) {
        idx <- sort(c(match(r, lv), match(c, lv)))     
        if (!any(vapply(sig_pairs, identical, logical(1), idx)))
          sig_pairs[[length(sig_pairs) + 1]] <- idx
      }
    }
  }
  
  if (length(sig_pairs) > 0) {
    step   <- 0.05 * y_max
    height <- y_max + step
    for (pair in sig_pairs) {
      segments(pair[1], height, pair[2], height, lwd = 1.5)
      text(mean(pair), height * 1.02, "*", cex = 3)
      height <- height + step
    }
  }
  
  mtext(side = 1,
        text = "Age group (years)",
        line = 3,
        cex  = 2.5)
}

# PART 9: Conover (Mean,Sd,Cov,Fat)------------------------------------------------------------------
{
  filtered_dataAll_rev <- dataAll %>%
    filter(ABx <= 11) %>%
    na.omit()
  
  filtered_dataAll_rev$AgeGroups_rev <- cut(
    filtered_dataAll_rev$ABx,
    breaks  = c(1, 5, 7.5, 11),
    labels  = c("1-5yo", "5-7.5yo", "7.5-11yo"),
    right   = FALSE
  )
  
  variables <- c("Mean", "Sd", "Cov", "Fat")
  pairs     <- list(c("1-5yo",  "5-7.5yo"),
                    c("1-5yo",  "7.5-11yo"),
                    c("5-7.5yo","7.5-11yo"))
  
  conf_int <- function(diff, sd1, sd2, n1, n2, alpha = 0.05) {
    se  <- sqrt((sd1^2 / n1) + (sd2^2 / n2))
    t_c <- qt(1 - alpha / 2, df = n1 + n2 - 2)
    c(diff - t_c * se, diff + t_c * se)
  }
  
  for (var in variables) {
    cat("\n=========== ", var, " ===========\n", sep = "")
    
    con <- kwAllPairsConoverTest(
      as.formula(paste(var, "~ AgeGroups_rev")),
      data            = filtered_dataAll_rev,
      p.adjust.method = "none"
    )
    
    raw_p_mat <- con$p.value                             # matrix
    bh_p_mat  <- matrix(
      p.adjust(as.vector(raw_p_mat), method = "BH"),     # vector → matrix
      nrow = nrow(raw_p_mat),
      dimnames = dimnames(raw_p_mat)
    )
    
    g <- split(filtered_dataAll_rev[[var]],
               filtered_dataAll_rev$AgeGroups_rev)
    
    tidy <- data.frame(Pair = character(0),
                       Diff = numeric(0),
                       CI_low = numeric(0),
                       CI_high = numeric(0),
                       Raw_p = numeric(0),
                       BH_p = numeric(0),
                       stringsAsFactors = FALSE)
    
    for (pr in pairs) {
      lv <- pr[1]; hv <- pr[2]            # low vs high
      
      # fetch p-value regardless of matrix orientation
      if (lv %in% rownames(raw_p_mat) && hv %in% colnames(raw_p_mat)) {
        p_raw <- raw_p_mat[lv, hv]; p_bh <- bh_p_mat[lv, hv]
      } else {
        p_raw <- raw_p_mat[hv, lv]; p_bh <- bh_p_mat[hv, lv]
      }
      
      diff <- mean(g[[hv]]) - mean(g[[lv]])
      ci   <- conf_int(diff,
                       sd(g[[hv]]), sd(g[[lv]]),
                       length(g[[hv]]), length(g[[lv]]))
      
      tidy <- rbind(tidy,
                    data.frame(
                      Pair   = paste(lv, hv, sep = " vs "),
                      Diff   = round(diff, 3),
                      CI_low = round(ci[1], 3),
                      CI_high= round(ci[2], 3),
                      Raw_p  = signif(p_raw, 3),
                      BH_p   = signif(p_bh, 3),
                      stringsAsFactors = FALSE
                    )
      )
    }
    
    print(tidy, row.names = FALSE)
  }
}

# =========== Mean
#  Pair                 Diff   CI_low  CI_high   Raw_p  BH_p(Benjamini–Hochberg)
#1-5yo vs 5-7.5yo    403.717  145.142  662.292 0.00356   0.00534**
#1-5yo vs 7.5-11yo   620.750   93.482 1148.017 0.00105   0.00315**
#5-7.5yo vs 7.5-11yo 217.032 -345.703  779.767 0.55000   0.55000

#=========== Sd 
#  Pair                 Diff   CI_low  CI_high    Raw_p      BH_p
#1-5yo vs 5-7.5yo    414.100  180.253  647.947 7.10e-04  1.07e-03**
#1-5yo vs 7.5-11yo   861.902  288.900 1434.904 2.46e-05  7.38e-05***
#5-7.5yo vs 7.5-11yo 447.802 -154.308 1049.912 1.77e-01  1.77e-01

#=========== Cov 
#  Pair                Diff   CI_low   CI_high    Raw_p     BH_p
#1-5yo vs 5-7.5yo     0.156   -0.019    0.331  0.025900 0.038900 (Mean-based CI includes zero)
#1-5yo vs 7.5-11yo    0.353    0.119    0.586  0.000275 0.000824***
#5-7.5yo vs 7.5-11yo  0.197   -0.028    0.421  0.064200 0.064200

#=========== Fat 
#  Pair               Diff  CI_low  CI_high    Raw_p    BH_p
#1-5yo vs 5-7.5yo    1.503  -0.690   3.697  0.147000  0.14700
#1-5yo vs 7.5-11yo   5.283   2.043   8.523  0.000537  0.00161**
#5-7.5yo vs 7.5-11yo 3.779   0.117   7.442  0.017000  0.02550*

# PART 10: Tukey's HSD (MFD, MFA, CFA)------------------------------------------
{
  filtered_dataAll_rev <- dataAll %>% 
    filter(ABx <= 11) %>% 
    na.omit()
  
  filtered_dataAll_rev$AgeGroups_rev <- cut(
    filtered_dataAll_rev$ABx,
    breaks  = c(1, 5, 7.5, 11),
    labels  = c("1-5yo", "5-7.5yo", "7.5-11yo"),
    right   = FALSE
  )
  
  fmt_num <- function(x, digits = 4, sci_cut = 1e-20) {
    ifelse(abs(x) < sci_cut & x != 0,
           formatC(x, format = "e", digits = 2),     
           sprintf(paste0("%.", digits, "f"), x))     
  }
  
  variables <- c("MFD", "MFA", "CFA")
  
  for (var in variables) {
    cat("\nAnalyzing", var, "\n")
    
    model <- aov(as.formula(paste(var, "~ AgeGroups_rev")),
                 data = filtered_dataAll_rev)
    tk <- TukeyHSD(model)$AgeGroups_rev
    tk_df <- as.data.frame(tk)
    tk_df$BH_p <- p.adjust(tk_df$`p adj`, method = "BH")
    
    out <- cbind(Comparison = rownames(tk_df),
                 Diff       = fmt_num(tk_df$diff, 3),
                 Lower_CI   = fmt_num(tk_df$lwr, 3),
                 Upper_CI   = fmt_num(tk_df$upr, 3),
                 Tukey_p    = fmt_num(tk_df$`p adj`, 8),
                 BH_p       = fmt_num(tk_df$BH_p, 8))
    print(out)
  }
}

# Analyzing MFD 
# Comparison            Diff     Lower_CI   Upper_CI    Tukey_p   BH_p    
# "5-7.5yo-1-5yo"    "-431.248"  "-641.036" "-221.459"     "0"    "1e-04" ***
# "7.5-11yo-1-5yo"   "-573.517"  "-791.871" "-355.162"     "0"        "0" ***    
# "7.5-11yo-5-7.5yo" "-142.269"  "-352.058"    "67.52"  "0.2335"  "0.2335"

# Analyzing MFA 
# Comparison            Diff     Lower_CI   Upper_CI    Tukey_p   BH_p    
# "5-7.5yo-1-5yo"    "-11.778"   "-24.753"   "1.198"    "0.0812"  "0.1218"
# "7.5-11yo-1-5yo"   "-18.744"   "-32.249"  "-5.239"    "0.0049"  "0.0147" *
# "7.5-11yo-5-7.5yo"  "-6.966"   "-19.942"   "6.009"    "0.3951"  "0.3951"

# Analyzing CFA 
# Comparison            Diff     Lower_CI   Upper_CI    Tukey_p   BH_p    
# "5-7.5yo-1-5yo"    "10.191"    "-1.788"   "22.17"     "0.1078"  "0.1617"
# "7.5-11yo-1-5yo"   "13.187"    "0.718"    "25.655"    "0.0364"  "0.1092"
# "7.5-11yo-5-7.5yo"  "2.996"    "-8.984"   "14.975"    "0.8133"  "0.8133"

### -----------------------------### 15. Effect Size Analysis with 95% CI-----------
# PART 1: Hedges' d: unbiased d-----------------------------------------------------------
{
  group1_rev <- dataAll %>% filter(ABx < 5)                # 1–<5
  group2_rev <- dataAll %>% filter(ABx >= 5 & ABx < 7.5)   # 5–<7.5
  group3_rev <- dataAll %>% filter(ABx >= 7.5 & ABx < 11)  # 7.5–<11
  
  signed_d <- function(low_vec, high_vec) {
    effsize::cohen.d(high_vec, low_vec, hedges.correction = TRUE)$estimate
  }
  
  # ──────────────────────────────────────────────
  bootstrap_ci_hedges <- function(data_low, data_high,
                                  n_iter = 5000,
                                  conf_level = 0.95,
                                  seed = NULL) {
    
    if (!is.null(seed)) set.seed(seed)
    
    boot_d <- replicate(n_iter, {
      s_low  <- sample(data_low,  length(data_low),  replace = TRUE)
      s_high <- sample(data_high, length(data_high), replace = TRUE)
      signed_d(s_low, s_high)
    })
    
    list(
      original_d = signed_d(data_low, data_high),
      boot_mean  = mean(boot_d),
      ci_lower   = quantile(boot_d, (1 - conf_level) / 2),
      ci_upper   = quantile(boot_d, 1 - (1 - conf_level) / 2)
    )
  }
  
  variables <- c("MFD", "MFA", "CFA")
  
  for (var in variables) {
    cat("\nAnalyzing", var, "\n")
    
    g1 <- as.numeric(na.omit(group1_rev[[var]]))  
    g2 <- as.numeric(na.omit(group2_rev[[var]]))  
    g3 <- as.numeric(na.omit(group3_rev[[var]]))  
    
    ## 1–5 vs 5–7.5
    r12 <- bootstrap_ci_hedges(g1, g2)
    p12 <- pwr.t2n.test(n1 = length(g1), n2 = length(g2),
                        d = abs(r12$original_d), sig.level = 0.05)$power * 100
    n12 <- pwr.t.test(d = abs(r12$original_d), power = 0.80,
                      sig.level = 0.05, type = "two.sample")$n
    
    cat("1–5yo vs 5–7.5yo  Hedges' d:",  r12$original_d, "\n",
        "  Bootstrap mean:",            r12$boot_mean,  "\n",
        "  95% CI:",                    r12$ci_lower, "-", r12$ci_upper, "\n",
        "  Observed power (%):",        round(p12, 2), "\n",
        "  Required N (80% power):",    n12, "\n\n")
    
    ## 1–5 vs 7.5–11
    r13 <- bootstrap_ci_hedges(g1, g3)
    p13 <- pwr.t2n.test(n1 = length(g1), n2 = length(g3),
                        d = abs(r13$original_d), sig.level = 0.05)$power * 100
    n13 <- pwr.t.test(d = abs(r13$original_d), power = 0.80,
                      sig.level = 0.05, type = "two.sample")$n
    
    cat("1–5yo vs 7.5–11yo Hedges' d:",  r13$original_d, "\n",
        "  Bootstrap mean:",             r13$boot_mean,  "\n",
        "  95% CI:",                     r13$ci_lower, "-", r13$ci_upper, "\n",
        "  Observed power (%):",         round(p13, 2), "\n",
        "  Required N (80% power):",     n13, "\n\n")
    
    ## 5–7.5 vs 7.5–11
    r23 <- bootstrap_ci_hedges(g2, g3)
    p23 <- pwr.t2n.test(n1 = length(g2), n2 = length(g3),
                        d = abs(r23$original_d), sig.level = 0.05)$power * 100
    n23 <- pwr.t.test(d = abs(r23$original_d), power = 0.80,
                      sig.level = 0.05, type = "two.sample")$n
    
    cat("5–7.5yo vs 7.5–11yo Hedges' d:", r23$original_d, "\n",
        "  Bootstrap mean:",              r23$boot_mean,  "\n",
        "  95% CI:",                      r23$ci_lower, "-", r23$ci_upper, "\n",
        "  Observed power (%):",          round(p23, 2), "\n",
        "  Required N (80% power):",      n23, "\n\n")
  }
}

# Analyzing MFD 
# 1–5yo vs 5–7.5yo  Hedges' d: -1.732657 
#    Bootstrap mean: -1.887192 
#    95% CI: -3.056972 to -1.010826 
#    Observed power (%): 98.12 
#    Required N (80% power, per group): 6.349431 

# 1–5yo vs 7.5–11yo Hedges' d: -2.507261 
#    Bootstrap mean: -2.70815 
#    95% CI: -3.974917 to -1.886407 
#    Observed power (%): 99.99 
#    Required N (80% power, per group): 3.748474 

# 5–7.5yo vs 7.5–11yo Hedges' d: -0.8765553 
#    Bootstrap mean: -0.9397959 
#    95% CI: -1.68972 to -0.255354 
#    Observed power (%): 53.42 
#    Required N (80% power, per group): 21.43455 


# Analyzing MFA 
# 1–5yo vs 5–7.5yo  Hedges' d: -1.045281 
#   Bootstrap mean: -1.119667 
#   95% CI: -2.075234 to -0.3285498 
#   Observed power (%): 68.4 
#   Required N (80% power, per group): 15.38927 

# 1–5yo vs 7.5–11yo Hedges' d: -1.362198 
#    Bootstrap mean: -1.462821 
#    95% CI: -2.543632 to -0.6491953 
#    Observed power (%): 85.95 
#    Required N (80% power, per group): 9.522869 

# 5–7.5yo vs 7.5–11yo Hedges' d: -0.4690813 
#    Bootstrap mean: 0.4991922 
#    95% CI: -1.424109 to 0.3315161 
#    Observed power (%): 19.48 
#    Required N (80% power, per group): 72.3142 


# Analyzing CFA 
# 1–5yo vs 5–7.5yo  Hedges' d: 0.9942917 
#    Bootstrap mean: 1.056281 
#    95% CI: 0.313086 to 1.915683 
#    Observed power (%): 64.07 
#    Required N (80% power, per group): 16.89485 

# 1–5yo vs 7.5–11yo Hedges' d: 1.040386 
#    Bootstrap mean: 1.118257  
#    95% CI: 0.3116748 to 2.161084 
#    Observed power (%): 64.12 
#    Required N (80% power, per group): 15.52422 

# 5–7.5yo vs 7.5–11yo Hedges' d: 0.2163049 
#    Bootstrap mean: -0.2364516 
#    95% CI: -0.6194971 to 1.130681
#    Observed power (%): 7.97 
#    Required N (80% power, per group): 336.4719

# PART 2: Cliff's Delta-----------------------------------------------------------
{
  group1_rev <- dataAll %>% filter(ABx < 5)
  group2_rev <- dataAll %>% filter(ABx >= 5 & ABx < 7.5)
  group3_rev <- dataAll %>% filter(ABx >= 7.5 & ABx < 11)
  
  signed_delta <- function(low_vec, high_vec) {
    effsize::cliff.delta(high_vec, low_vec)$estimate   
  }
  
  boot_cliff <- function(x_low, x_high,
                         n_boot = 5000,
                         conf   = 0.95,
                         seed   = NULL) {
    if (!is.null(seed)) set.seed(seed)
    
    boot_d <- replicate(n_boot, {
      signed_delta(sample(x_low,  length(x_low),  TRUE),
                   sample(x_high, length(x_high), TRUE))
    })
    
    delta_orig <- signed_delta(x_low, x_high)
    ci <- quantile(boot_d, probs = c((1-conf)/2, 1-(1-conf)/2))
    
    list(delta = delta_orig,
         mean  = mean(boot_d),
         ci_lo = ci[1], ci_hi = ci[2])
  }
  
  
  cliff_to_d <- function(delta) {
    p <- (delta + 1) / 2
    p <- ifelse(p <= 0, 1e-10, ifelse(p >= 1, 1-1e-10, p))
    sqrt(2) * qnorm(p)
  }
  
  vars <- c("Mean", "Sd", "Cov", "Fat")
  
  for (v in vars) {
    cat("\n=====  ", v, "  (Bootstrap Cliff’s Δ + Power)  =====\n")
    
    g1 <- as.numeric(na.omit(group1_rev[[v]]))
    g2 <- as.numeric(na.omit(group2_rev[[v]]))
    g3 <- as.numeric(na.omit(group3_rev[[v]]))
    
    pairs <- list(
      `1–5 vs 5–7.5`   = list(a = g1, b = g2),
      `1–5 vs 7.5–11`  = list(a = g1, b = g3),
      `5–7.5 vs 7.5–11`= list(a = g2, b = g3)
    )
    
    for (lbl in names(pairs)) {
      a <- pairs[[lbl]]$a
      b <- pairs[[lbl]]$b
      
      res <- boot_cliff(a, b)
      d   <- cliff_to_d(res$delta)
      
      pow <- if (length(a) > 1 && length(b) > 1) {
        pwr.t2n.test(n1 = length(a), n2 = length(b),
                     d = d, sig.level = 0.05)$power * 100
      } else NA
      
      reqN <- pwr.t.test(d = abs(d), power = 0.80,
                         sig.level = 0.05, type = "two.sample")$n
      
      cat(lbl, "\n",
          "  Cliff’s Δ :", round(res$delta, 4),
          "  95% CI [", round(res$ci_lo, 4), ",", round(res$ci_hi, 4), "]\n",
          "  Cohen’s d (approx):", round(d, 4), "\n",
          "  Observed power (%) :", ifelse(is.na(pow), "NA", round(pow, 2)), "\n",
          "  Required N (80% power, per group):", round(reqN, 2), "\n\n")
    }
  }
}

# =====   Mean   (Bootstrap Cliff’s Δ + Power)
#   1–5 vs 5–7.5 
# Cliff’s Δ : 0.6084   95% CI [ 0.2168 , 0.9161 ]
# Cohen’s d (approx): 1.2116 
# Observed power (%) : 80.68 
# Required N (80% power, per group): 11.74 

# 1–5 vs 7.5–11 
# Cliff’s Δ : 0.8017   95% CI [ 0.4876 , 1 ]
# Cohen’s d (approx): 1.8191 
# Observed power (%) : 98.19 
# Required N (80% power, per group): 5.88 

# 5–7.5 vs 7.5–11 
# Cliff’s Δ : 0.0909   95% CI [ -0.3846 , 0.5664 ]
# Cohen’s d (approx): 0.1615 
# Observed power (%) : 6.65 
# Required N (80% power, per group): 602.95 


# =====   Sd   (Bootstrap Cliff’s Δ + Power)
#   1–5 vs 5–7.5 
# Cliff’s Δ : 0.7343   95% CI [ 0.3846 , 0.972 ]
# Cohen’s d (approx): 1.5739 
# Observed power (%) : 95.63 
# Required N (80% power, per group): 7.43 

# 1–5 vs 7.5–11 
# Cliff’s Δ : 0.8843   95% CI [ 0.6529 , 1 ]
# Cohen’s d (approx): 2.2247 
# Observed power (%) : 99.86 
# Required N (80% power, per group): 4.38 

# 5–7.5 vs 7.5–11 
# Cliff’s Δ : 0.3007   95% CI [ -0.1612 , 0.7343 ]
# Cohen’s d (approx): 0.5463 
# Observed power (%) : 24.75 
# Required N (80% power, per group): 53.58 


# =====   Cov   (Bootstrap Cliff’s Δ + Power)
#   1–5 vs 5–7.5 
# Cliff’s Δ : 0.4825   95% CI [ 0.007 , 0.8881 ]
# Cohen’s d (approx): 0.9153 
# Observed power (%) : 57 
# Required N (80% power, per group): 19.74 

# 1–5 vs 7.5–11 
# Cliff’s Δ : 0.8347   95% CI [ 0.4876 , 1 ]
# Cohen’s d (approx): 1.9622 
# Observed power (%) : 99.21 
# Required N (80% power, per group): 5.24 

# 5–7.5 vs 7.5–11 
# Cliff’s Δ : 0.3986   95% CI [ -0.0629 , 0.7902 ]
# Cohen’s d (approx): 0.7388 
# Observed power (%) : 40.71 
# Required N (80% power, per group): 29.75 


# =====   Fat   (Bootstrap Cliff’s Δ + Power)
#   1–5 vs 5–7.5 
# Cliff’s Δ : 0.3566   95% CI [ -0.1049 , 0.7762 ]
# Cohen’s d (approx): 0.6548 
# Observed power (%) : 33.34 
# Required N (80% power, per group): 37.6 

# 1–5 vs 7.5–11 
# Cliff’s Δ : 0.7521   95% CI [ 0.3554 , 1 ]
# Cohen’s d (approx): 1.634 
# Observed power (%) : 95.37 
# Required N (80% power, per group): 6.98 

# 5–7.5 vs 7.5–11 
# Cliff’s Δ : 0.5664   95% CI [ 0.1189 , 0.9441 ]
# Cohen’s d (approx): -1.1075 
# Observed power (%) : 73.36 
# Required N (80% power, per group): 13.83


# PART 3: Glass's Delta--------------------------------------------------------
# if var = FALSE
{
  group1_rev <- dataAll %>% filter(ABx < 5)                 # 1–5 
  group2_rev <- dataAll %>% filter(ABx >= 5 & ABx < 7.5)    # 5–7.5 
  group3_rev <- dataAll %>% filter(ABx >= 7.5 & ABx < 11)   # 7.5–11 
  
  calc_glass_delta_rev <- function(x_ref, x_cmp) {
    ## Glass’s Δ
    d_glass <- (mean(x_cmp) - mean(x_ref)) / sd(x_ref)
    
    ## 95 % CI（MBESS::ci.smd）
    ci <- MBESS::ci.smd(
      smd      = d_glass,
      n.1      = length(x_cmp),
      n.2      = length(x_ref),
      conf.lvl = 0.95)
    
    reqN <- pwr::pwr.t.test(d = abs(d_glass),
                            power = 0.80,
                            sig.level = 0.05,
                            type = "two.sample")$n
    
    obsP <- pwr::pwr.t2n.test(n1 = length(x_ref),
                              n2 = length(x_cmp),
                              d  = d_glass,
                              sig.level = 0.05)$power * 100
    
    list(d = d_glass,
         ci_lo = ci$Lower.Conf.Limit.smd,
         ci_hi = ci$Upper.Conf.Limit.smd,
         power = obsP,
         reqN  = reqN)
  }
  
  var <- "MFD"
  cat("\n=====  Glass’s Δ for", var, "=====\n")
  
  g1 <- as.numeric(na.omit(group1_rev[[var]]))
  g2 <- as.numeric(na.omit(group2_rev[[var]]))
  g3 <- as.numeric(na.omit(group3_rev[[var]]))
  
  ### 1–5  vs 5–7.5
  r12 <- calc_glass_delta_rev(g1, g2)
  cat("1–5yo vs 5–7.5yo  Glass’s Δ :", round(r12$d, 4),
      "  95% CI [", round(r12$ci_lo, 4), ",", round(r12$ci_hi, 4), "]\n",
      "  Observed power (%)         :", round(r12$power, 2), "\n",
      "  Required N (each, 80% pow) :", round(r12$reqN, 2), "\n\n")
  
  ### 1–5  vs 7.5–11
  r13 <- calc_glass_delta_rev(g1, g3)
  cat("1–5yo vs 7.5–11yo Glass’s Δ :", round(r13$d, 4),
      "  95% CI [", round(r13$ci_lo, 4), ",", round(r13$ci_hi, 4), "]\n",
      "  Observed power (%)         :", round(r13$power, 2), "\n",
      "  Required N (each, 80% pow) :", round(r13$reqN, 2), "\n\n")
  
  ### 5–7.5 vs 7.5–11
  r23 <- calc_glass_delta_rev(g2, g3)
  cat("5–7.5yo vs 7.5–11yo Glass’s Δ :", round(r23$d, 4),
      "  95% CI [", round(r23$ci_lo, 4), ",", round(r23$ci_hi, 4), "]\n",
      "  Observed power (%)          :", round(r23$power, 2), "\n",
      "  Required N (each, 80% pow)  :", round(r23$reqN, 2), "\n")
}

# =====  Glass’s Δ for MFD
#   1–5yo vs 5–7.5yo  
# Glass’s Δ : -1.4797   
# 95% CI [ -2.3797 , -0.5537 ]
# Observed power (%)         : 93.19 
# Required N (each, 80% pow) : 8.25 

# 1–5yo vs 7.5–11yo 
# Glass’s Δ : -1.9679   
# 95% CI [ -2.9839 , -0.9195 ]
# Observed power (%)         : 99.23 
# Required N (each, 80% pow) : 5.21 

# 5–7.5yo vs 7.5–11yo 
# Glass’s Δ : -0.7594   
# 95% CI [ -1.5848 , 0.0819 ]
# Observed power (%)          : 42.59 
# Required N (each, 80% pow)  : 28.21


### -----------------------------### 16. PROVENANCE AND REPRODUCIBILITY NOTE -----------------------------

# This annotated R script accompanies the doctoral dissertation:
# "A Systematic Histopathological Study of Duchenne Muscular Dystrophy
# Using Semi-Quantitative Image Analysis, Digital Restoration Techniques,
# and Exploratory Statistical Approaches" (Hokkaido University, 2026).
#
# The script was prepared to enhance transparency, reproducibility, and
# research traceability of the statistical workflow used in the dissertation.
#
# Scope of this script:
# - exploratory analyses used to refine the statistical strategy,
# - confirmatory analyses used for the main inferential framework, and
# - supplementary visualization and robustness checks retained for documentation.
#
# Relation to the published article:
# Part of this workflow corresponds to analyses reported in the related
# peer-reviewed article:
# "Myofibre Density Reveals a Critical Threshold Around Age 6 in
# Steroid-Naïve Duchenne Muscular Dystrophy: A Retrospective
# Observational Study" (Neuropathology and Applied Neurobiology, 2026).
#
# Interpretation note:
# Exploratory or preliminary sections are retained for transparency and to
# document how the final analytical framework was established. They should
# not be interpreted as fully independent confirmatory evidence unless
# explicitly stated in the dissertation text.
#
# Research-use notice:
# This script is provided for academic transparency and methodological
# reproducibility. It is not intended for clinical decision-making.


### Statistics
# Yamakado T (MD): drafted, wrote, and executed the R-based statistical analyses.
#  1) Department of Cancer Pathology, Faculty of Medicine, Hokkaido University, Japan
#  2) Department of Pathology, National Hospital Organization Hokkaido Medical
#     Center, Japan
#  3) Center for Neuromuscular Disease, Child Health and Development, 
#     National Hospital Organization Hokkaido Medical Center, Japan


# Record the R session information used for this analysis
sessionInfo() 


### The R version used during development of this script is as follows:
# R version 4.3.2 (2022-10-31 ucrt)