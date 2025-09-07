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
###------------------------ Main Source Code------------------------------------
# "Myofibre Density Reveals a Critical Threshold around Age 6 in Steroid-Naïve
#  Duchenne Muscular Dystrophy: A Retrospective Observational Study"

# Description:
# This R script is provided to ensure reproducibility of the primary analyses and results 
# reported in the above manuscript. It includes essential data preparation, modeling, 
# and visualization steps. For additional exploratory or supplementary code, please refer 
# to the supplementary section.

# Note:
# - The script relies on R (v4.3.2) and standard packages for data manipulation, 
#   statistical analysis, and graphics.
# - All data were de-identified and handled in accordance with institutional review board
#   and ethical guidelines.
# - For questions or clarifications, please contact the corresponding author.

################################################################################
# Preparation--------------------------------------------------------------
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

### -----------------------------### 1. Correlation Analysis#----------------------
# PART 1: Overview--------------------------------------------------------------
df_cor <- L %>% 
  dplyr::select(ABx, Mean, Cov, Sd, MFD, MFA, CFA, IntN, Opaque, NFA, RFA, Fat)  

cor(df_cor, method = "pearson")
cor(df_cor, method = "spearman")

# examples (showed the strongest correlation with ABx)
cor.test(ABx, MFD, method = "spearman") # -0.8454973
cor.test(ABx, MFD) # -0.800855
cor.test(ABx, log(MFD)) # -0.8495001
cor.test(log(ABx), log(MFD)) # -0.8597013

qplot(ABx, MFD) + # `qplot()` was deprecated in ggplot2 3.4.0.
  geom_point(size = 2, colour = "black") +
  annotate("text", x = c(11), y = c(1500), size = 10, fontface = "bold", label = c("rho = -0.85"))+
  geom_smooth(span = 1, colour = "red", se = F)+
  scale_x_continuous(limits = c(0, 16.5), breaks = 0:16.5) +
  theme_bw()+
  theme(
    axis.title.x = element_text(size = 20, face = "bold"), 
    axis.title.y = element_text(size = 20, face = "bold"), 
    axis.text.x = element_text(size = 16), 
    axis.text.y = element_text(size = 16)  
  )

# PART 2: Mean------------------------------------------------------------------
cor.test(ABx, Mean, method = "spearman") # 0.6124302 
SpearmanRho(ABx, Mean, conf.level = 0.95) # 0.3640151 to 0.7795001
ggplot(dataAll, aes(x = ABx, y = Mean)) +
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

# PART 3: Sd--------------------------------------------------------------------
cor.test(ABx, Sd, method = "spearman") # 0.7369515 
SpearmanRho(ABx, Sd, conf.level = 0.95) # 0.5458704 to 0.8551655
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
SpearmanRho(ABx, Cov, conf.level = 0.95) # 0.4277368 to 0.8074688
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
SpearmanRho(ABx, MFD, conf.level = 0.95) # -0.9172546 to -0.7205808
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

L_sub <- subset(L, ABx < 11)

ABx_sub <- L_sub$ABx
MFD_sub <- L_sub$MFD

cor.test(ABx_sub, MFD_sub, method = "spearman") # rho = -0.810084
SpearmanRho(ABx_sub, MFD_sub, conf.level = 0.95) # -0.9002903 to -0.6531642

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

# PART 6: MFA-------------------------------------------------------------------
cor.test(ABx, MFA, method = "spearman") # -0.6032389 
SpearmanRho(ABx, MFA, conf.level = 0.95) # -0.7737149 to -0.3513026
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
cor.test(ABx, CFA, method = "spearman") # 0.483751 
SpearmanRho(ABx, CFA, conf.level = 0.95) # 0.1940818 to 0.6958267
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
SpearmanRho(ABx, Fat, conf.level = 0.95) # 0.4759424 to 0.8275506
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
SpearmanRho(ABx, Opaque, conf.level = 0.95) # 0.02721999 to 0.59786345
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


# PART 13: Subgroup (1-11 yrs)--------------------------------------------------
{
  L_sub <- subset(L, ABx < 11)
  
  ABx_sub <- L_sub$ABx
  Mean_sub <- L_sub$Mean
  MFD_sub <- L_sub$MFD
  Sd__sub <- L_sub$Sd
  Cov_sub <- L_sub$Cov
  MFA_sub <- L_sub$MFA
  Fat_sub <- L_sub$Fat
  CFA_sub <- L_sub$CFA
  Opaque_sub <- L_sub$Opaque
  }

cor.test(ABx_sub, MFD_sub, method = "spearman") # rho = -0.810084, p-value = 2.567e-07***
cor.test(ABx_sub, MFD_sub, method = "pearson") # rho = -0.832894, p-value = 25.428e-10***

cor.test(ABx_sub, Mean_sub, method = "spearman") # rho = 0.5778711, p-value = 0.0003552*** 
cor.test(ABx_sub, Mean_sub, method = "pearson") # rho = 0.5062608, p-value = 0.001916** 

cor.test(ABx_sub, Sd__sub, method = "spearman") # rho = 0.6966387, p-value = 6.459e-06*** 
cor.test(ABx_sub, Sd__sub, method = "pearson") # rho = 0.5713054, p-value = 0.0003377*** 

cor.test(ABx_sub, Cov_sub, method = "spearman") # rho = 0.6156863, p-value = 0.0001151*** 
cor.test(ABx_sub, Cov_sub, method = "pearson") # rho = 0.4968486, p-value = 0.002396**

cor.test(ABx_sub, MFA_sub, method = "spearman") # rho = -0.512605, p-value = 0.001884**
cor.test(ABx_sub, MFA_sub, method = "pearson") # rho = -0.5353767, p-value = 0.0009195***

cor.test(ABx_sub, Fat_sub, method = "spearman") # rho = 0.6518207, p-value = 3.436e-05***
cor.test(ABx_sub, Fat_sub, method = "pearson") # rho = 0.542986, p-value = 0.0007507***

cor.test(ABx_sub, CFA_sub, method = "spearman") # rho = 0.3977591, p-value = 0.01863* 
cor.test(ABx_sub, CFA_sub, method = "pearson") # rho = 0.427533, p-value = 0.01041* 

cor.test(ABx_sub, Opaque_sub, method = "spearman") # rho = 0.2980392, p-value = 0.08233 NS
cor.test(ABx_sub, Opaque_sub, method = "pearson") # rho = 0.2642391, p-value = 0.1251 NS


### -----------------------------### 2. Segmented Regression Analysis--------------
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

### -----------------------------### 3. Robustness Check (Supplementary)----
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


###------------------------------### 4. Analysis of Variance-----------------------
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


# Muscle Fiber Density (MFD) box-plot with Games–Howell asterisks

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

### -----------------------------### 5. Effect Size Analysis with 95% CI-----------
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

### -----------------------------### 6. Replacement by ATPase----------------------
# PART 1: Preparation and Correlation Analysis----------------------------------
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

# re-check correlation coefficient and plot
cor.test(ABx, MFD_ATP, method = "spearman") # -0.8665062
cor.test(ABx, MFD_ATP, method = "pearson") # -0.8280345
cor.test(log(ABx), log(MFD_ATP), method = "pearson") # -0.8680387 

qplot(ABx, MFD_ATP) +
  geom_point(size = 2, colour = "black") +
  geom_smooth(span =1, colour = "red", se = F) +
  scale_x_continuous(limits = c(0, 16.5), breaks = 0:16.5) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 40, face = "bold"), 
    axis.title.y = element_text(size = 40, face = "bold"), 
    axis.text.x = element_text(size = 16), 
    axis.text.y = element_text(size = 16)  
  )

# PART 2: (MFD_ATP)Segmented regression analysis -----------------------------
{
  dataMFD_ATP <- data.frame(ABx = ABx, MFD_ATP = MFD_ATP)
  dataSegment_MFD_ATP <- dataMFD_ATP %>%
    filter(ABx >= 1 & ABx <= 11)
  
  lm_init_MFD_ATP <- lm(MFD_ATP ~ ABx, data = dataSegment_MFD_ATP)
  segmented_model_MFD_ATP <- segmented(lm_init_MFD_ATP, seg.Z = ~ABx, psi = 6)
}

# Summarize segmented regression results
summary(segmented_model_MFD_ATP) # Adjusted R-squared: 0.8223

# Extract and visualize breakpoint
plot(segmented_model_MFD_ATP, main = "Segmented Regression for MFD_ATP vs Age")
abline(v = segmented_model_MFD_ATP$psi[2], col = "red", lty = 2)

# Extract slopes (before and after the breakpoint)
slope(segmented_model_MFD_ATP)
#            Est. St.Err. t value CI(95%).l CI(95%).u
# slope1 -161.750  19.116 -8.4615  -200.740  -122.760
# slope2  -28.916  21.544 -1.3422   -72.854    15.023

# Confidence intervals for breakpoint and slopes
confint(segmented_model_MFD_ATP) # 6.42 [5.28036  to  7.55964]
anova(lm_init_MFD_ATP, segmented_model_MFD_ATP) # Pr(>F) = 0.0002935 ***

dataSegment_MFD_ATP$Predicted_MFD_ATP <- predict(segmented_model_MFD_ATP)
RMSE_original_ATP <- sqrt(mean((dataSegment_MFD_ATP$MFD_ATP - dataSegment_MFD_ATP$Predicted_MFD_ATP)^2))
print(RMSE_original_ATP) # 121.2629


# ------- (MFD_ATP)Segmented model plot--------------------------------------------------
{
  bp   <- segmented_model_MFD_ATP$psi[2]            # 6.25 y
  xseq <- seq(0, 11, 0.1)
  pred <- predict(segmented_model_MFD_ATP,
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
  ggplot(dataMFD_ATP %>% filter(ABx <= 11), aes(x = ABx, y = MFD_ATP)) +
    # background
    annotate("rect", xmin = 0,  xmax = bp,
             ymin = -Inf, ymax = Inf, fill = fill_pre, alpha = 1) +
    annotate("rect", xmin = bp, xmax = 11,
             ymin = -Inf, ymax = Inf, fill = fill_post, alpha = 1) +
    
    # plot: age grad
    geom_point(aes(colour = ABx), size = 4, stroke = 0.1) +
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
             y = max(dataMFD_ATP$MFD_ATP[dataMFD_ATP$ABx <= 11]) * 0.92,
             label = sprintf("Breakpoint  %.2f y", bp),
             hjust = 0, vjust = 1, angle = 0,
             colour = col_bp, size = 5.2, fontface = "bold") +
    
    # axis
    scale_x_continuous(limits = c(0, 11), breaks = 0:11) +
    scale_y_continuous(
      limits = c(0, max(dataMFD_ATP$MFD_ATP[dataMFD_ATP$ABx <= 11]) * 1.05),
      expand = expansion(mult = c(0, 0.02))
    ) +
    
    # label and theme
    labs(
      x = "ABx",
      y = "MFD (ATPase replaced)"
    ) +
    theme_bw() +
    theme(
      axis.title = element_text(size = 28, face = "bold.italic"),
      axis.text  = element_text(size = 12),
      legend.position = "none",
      plot.margin = margin(8, 12, 8, 8)
    )
}


# ------- (MFD_ATP)LOOCV analysis-----------------------------------------------
{
  # Prepare vectors to store LOOCV results
  n <- nrow(dataSegment_MFD_ATP)
  breakpoints_ATP <- numeric(n)
  CI_lower_ATP <- numeric(n)
  CI_upper_ATP <- numeric(n)
  
  # Perform LOOCV for segmented regression
  for (i in 1:n) {
    
    # Create LOOCV dataset by excluding one observation
    data_LOOCV_ATP <- dataSegment_MFD_ATP[-i, ]
    
    # Segmented regression model for LOOCV
    lm_LOOCV_ATP <- lm(MFD_ATP ~ ABx, data = data_LOOCV_ATP)
    
    # Perform segmented regression with initial breakpoint at age 6
    segmented_LOOCV_ATP <- tryCatch({
      segmented(lm_LOOCV_ATP, seg.Z = ~ABx, psi = 6)
    }, error = function(e) {
      return(NULL)
    })
    
    # If successful, save breakpoint and confidence intervals
    if (!is.null(segmented_LOOCV_ATP)) {
      breakpoints_ATP[i] <- segmented_LOOCV_ATP$psi[2]
      
      # Extract correct CI (lower and upper)
      CI_vals_ATP <- confint(segmented_LOOCV_ATP)[1, ]
      CI_lower_ATP[i] <- CI_vals_ATP[2]
      CI_upper_ATP[i] <- CI_vals_ATP[3]
    } else {
      # If model fails, assign NA
      breakpoints_ATP[i] <- NA
      CI_lower_ATP[i] <- NA
      CI_upper_ATP[i] <- NA
    }
  }
  
  # Remove NA results from LOOCV
  valid_breakpoints_ATP <- breakpoints_ATP[!is.na(breakpoints_ATP)]
  valid_CI_lower_ATP <- CI_lower_ATP[!is.na(CI_lower_ATP)]
  valid_CI_upper_ATP <- CI_upper_ATP[!is.na(CI_upper_ATP)]
  
  # Calculate LOOCV mean and median for breakpoint and confidence intervals
  mean_breakpoint_ATP <- mean(valid_breakpoints_ATP)
  median_breakpoint_ATP <- median(valid_breakpoints_ATP)
  mean_CI_lower_ATP <- mean(valid_CI_lower_ATP)
  mean_CI_upper_ATP <- mean(valid_CI_upper_ATP)
  
  # Output the summarized LOOCV results
  cat("LOOCV Mean Breakpoint_ATP:", mean_breakpoint_ATP, "\n")
  cat("LOOCV Median Breakpoint_ATP:", median_breakpoint_ATP, "\n")
  cat("LOOCV Mean CI_ATP:", mean_CI_lower_ATP, "to", mean_CI_upper_ATP, "\n")
  
  # Calculate range of breakpoints
  breakpoint_min_ATP <- min(breakpoints_ATP, na.rm = TRUE)
  breakpoint_max_ATP <- max(breakpoints_ATP, na.rm = TRUE)
  
  # Display range
  cat("Breakpoint Range_ATP:", breakpoint_min_ATP, "to", breakpoint_max_ATP, "\n")
  
  # Prediction by original data
  dataSegment_MFD_ATP$Predicted_MFD_ATP <- predict(segmented_model_MFD_ATP)
  
  # RMSE（original）
  RMSE_original_ATP <- sqrt(mean((dataSegment_MFD_ATP$MFD_ATP - dataSegment_MFD_ATP$Predicted_MFD_ATP)^2))
  cat("Original Data RMSE_ATP:", RMSE_original_ATP, "\n")
  
  {
    # Prepare vector to store predictions
    LOOCV_predictions_ATP <- numeric(n)
    
    # Perform LOOCV to get predicted values
    for (i in 1:n) {
      # Exclude one observation
      data_LOOCV_ATP <- dataSegment_MFD_ATP[-i, ]
      
      # Fit model
      lm_LOOCV_ATP <- lm(MFD_ATP ~ ABx, data = data_LOOCV_ATP)
      segmented_LOOCV_ATP <- tryCatch({
        segmented(lm_LOOCV_ATP, seg.Z = ~ABx, psi = 6)
      }, error = function(e) {
        return(NULL)
      })
      
      # If successful, predict the excluded observation
      if (!is.null(segmented_LOOCV_ATP)) {
        LOOCV_predictions_ATP[i] <- predict(segmented_LOOCV_ATP, 
                                            newdata = dataSegment_MFD_ATP[i,])
      } else {
        LOOCV_predictions_ATP[i] <- NA
      }
    }
    
    # Remove any NA predictions (just in case)
    valid_indices_ATP <- !is.na(LOOCV_predictions_ATP)
    LOOCV_actuals_ATP <- dataSegment_MFD_ATP$MFD_ATP[valid_indices_ATP]
    LOOCV_preds_ATP   <- LOOCV_predictions_ATP[valid_indices_ATP]
    
    # RMSE (LOOCV data)
    RMSE_LOOCV_ATP <- sqrt(mean((LOOCV_actuals_ATP - LOOCV_preds_ATP)^2))
    cat("LOOCV RMSE_ATP:", RMSE_LOOCV_ATP, "\n")
    
    # RMSE difference and ratio
    cat("Difference in RMSE (LOOCV_ATP - Original_ATP):", RMSE_LOOCV_ATP - RMSE_original_ATP, "\n")
    cat("RMSE Ratio_ATP (LOOCV_ATP / Original_ATP):", RMSE_LOOCV_ATP / RMSE_original_ATP, "\n")
  }
  
  # Visualize distribution of LOOCV breakpoints
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
# LOOCV Mean CI_ATP: 5.230753 to 7.564639 
# Breakpoint Range_ATP: 6.250004 to 6.420001 
# Original Data RMSE_ATP: 121.2629 
# LOOCV RMSE_ATP: 132.3527 
# Difference in RMSE (LOOCV_ATP - Original_ATP): 11.08983  
# RMSE Ratio_ATP (LOOCV_ATP / Original_ATP): 1.091453  

###------------------------------### 7. Batch Effect Validation--------------------
# PART 1: Basic comparison of MFA between H&E and Gomori---------------------------------------

HE <- c(55.74, 62.03, 67.74, 63.06, 46.33, 58.81, 68.59, 38.23, 60.54, 73.03,
        62.090, 39.35, 70.43, 50.20, 58.25, 51.86, 30.70, 59.25, 67.42, 27.66,
        50.96, 60.48, 43.13, 42.25, 38.62, 56.94, 31.04, 70.90)
Gomori <- c(51.42, 56.25, 63.66, 59.93, 48.14, 51.39, 69.42, 36.54, 61.07,
            66.46, 62.87, 36.77, 73.02, 49.04, 61.14, 53.31, 34.83, 60.12,
            69.29, 23.83, 50.72, 62.091, 42.74, 46.90, 41.25, 58.59, 
            32.14, 72.06)

# Wilcoxon signed-rank sum test and effect size
wilcox.test(HE, Gomori) # p-value = 0.9288
cliff.delta(HE, Gomori) # 0.01530612: negligible (-0.2848065 to 0.3126860)

# t-test
shapiro.test(HE) # p-value = 0.13
shapiro.test(Gomori) # p-value = 0.401
var.test(HE, Gomori) # p-value = 0.9287
t.test(HE, Gomori, var.equal = T) # p-value = 0.9137

# Effect size (Hedges' d with bootstrapped CI)
{
  # bootstrap function
  bootstrap_ci_hedges <- function(data1, data2,
                                  n_iter = 5000,
                                  conf_level = 0.95,
                                  seed = 1234) {
    if (!is.null(seed)) set.seed(seed)
    
    n1 <- length(data1)
    n2 <- length(data2)
    boot_d <- numeric(n_iter)
    
    for (i in seq_len(n_iter)) {
      s1 <- sample(data1, n1, replace = TRUE)
      s2 <- sample(data2, n2, replace = TRUE)
      boot_d[i] <- cohen.d(s1, s2, hedges.correction = TRUE)$estimate
    }
    
    # Hedges d (original data)
    original_d <- cohen.d(data1, data2, hedges.correction = TRUE)$estimate
    
    # CI
    alpha <- 1 - conf_level
    ci_lower <- quantile(boot_d, alpha / 2)
    ci_upper <- quantile(boot_d, 1 - alpha / 2)
    
    list(
      original_d = original_d,
      boot_mean  = mean(boot_d),
      ci_lower   = ci_lower,
      ci_upper   = ci_upper,
      boot_dist  = boot_d           # if needed
    )
  }
  res <- bootstrap_ci_hedges(HE, Gomori)
  
  # Observed power (two‑sample）
  power_obs <- pwr.t2n.test(
    n1        = length(HE),
    n2        = length(Gomori),
    d         = res$original_d,
    sig.level = 0.05
  )$power * 100
  
  # required sample size (80%)
  required_N <- ceiling(
    pwr.t.test(
      d         = abs(res$original_d),
      power     = 0.80,
      sig.level = 0.05,
      type      = "two.sample"
    )$n
  )
  cat("HE vs Gomori  (n =", length(HE), "vs", length(Gomori), ")\n")
  cat("Hedges' d (Original)        :", sprintf("%.4f", res$original_d), "\n")
  cat("Hedges' d (Bootstrap mean)  :", sprintf("%.4f", res$boot_mean),  "\n")
  cat(sprintf("%d%% Bootstrap CI            : %.4f  –  %.4f\n",
              100 * 0.95, res$ci_lower, res$ci_upper))
  cat("Observed power (%)          :", sprintf("%.2f", power_obs), "\n")
  cat("Required N per group (80%%)  :", required_N, "\n")
}
# HE vs Gomori  (n = 28 vs 28 )
# Hedges' d (Original): 0.0287 
# Hedges' d (Bootstrap mean): 0.0396 
# 95% Bootstrap CI: -0.4799  to  0.5708
# Observed power (%): 5.13 
# Required N per group (80%): 19063


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
  #  HE vs Gomori : Bayesian paired‑samples analysis
  #  ────────────────────────────────────────────────────────
  #  * JZS Bayes t‑test  (Rouder et al., 2009)
  #  * Posterior ~ N(t/√n,1/n)  approximation
  #  * ROPE = ±0.30 (Lakens,2022)  + optional sensitivity grid
  #  * τ² reduction to check batch (specimen) effect
  
  
  ## 1. JZS Bayes factor (paired t) ---
  bf_ttest <- function(t, n, r = sqrt(2)/2){
    v <- n - 1
    j <- 1 + n * r^2
    logBF10 <- lgamma((v + 1)/2) - lgamma(v/2) -
      0.5 * log(pi * v) - 0.5 * log(j) -
      (v + 1)/2 * log(1 + t^2 / (v * j))
    exp(logBF10)
  }
  
  d   <- HE - Gomori
  n   <- length(d)
  tstat <- t.test(HE, Gomori, paired = TRUE)$statistic
  BF10  <- bf_ttest(tstat, n = length(HE), r = 0.707)      # default r = 0.707
  BF01  <- 1 / BF10
  
  
  ## 2. Posterior draws (normal approx) ---
  post_mean <- as.numeric(tstat) / sqrt(n)
  post_sd   <- sqrt(1/n)
  
  set.seed(123)
  delta_draws <- rnorm(10000, post_mean, post_sd)
  CrI95 <- quantile(delta_draws, c(.025, .975))
  
  
  ## 3. ROPE decision ---
  rope_fix <- 0.30                               # Lakens
  rope_pct <- mean(abs(delta_draws) < rope_fix) * 100
  
  # optional:  TEM‑based  δ_error  (kept for reference)
  sd_d      <- sd(d)
  sd_pool   <- sqrt((var(HE) + var(Gomori)) / 2)
  delta_err <- sd_d / sd_pool                    # ≈ 0.246
  # rope_tem  <- delta_err
  # pct_tem   <- mean(abs(delta_draws) < rope_tem) * 100
  
  #  Sensitivity grid (edit as needed)
  rope_grid <- c(0.25, 0.30, 0.40)
  grid_pct  <- sapply(rope_grid,
                      function(r) mean(abs(delta_draws) < r) * 100)
  
  
  ## 4. Batch (specimen) random effect τ² (corrected by linear mixed-effects model) ---
  library(Matrix)
  library(lme4)
  
  # Data preparation for mixed-effects model
  df <- data.frame(
    value = c(HE, Gomori),
    stain = factor(rep(c("HE", "Gomori"), each = length(HE))),
    specimen = factor(rep(seq_along(HE), times = 2))
  )
  
  # Fit models
  fit0 <- lmer(value ~ 1 + (1|specimen), data = df, REML = TRUE)
  fit1 <- lmer(value ~ stain + (1|specimen), data = df, REML = TRUE)
  
  # Extract specimen-level variance τ² from the models
  tau2_0 <- as.data.frame(VarCorr(fit0))$vcov[1]
  tau2_1 <- as.data.frame(VarCorr(fit1))$vcov[1]
  delta_tau <- 100 * (tau2_0 - tau2_1) / tau2_0      # % reduction
  
  ## 5. Output (unchanged) ---
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
  
  # List the sensitivity analysis
  cat("\n--  ROPE sensitivity analysis -----------------------------\n")
  for(i in seq_along(rope_grid)){
    cat(sprintf("ROPE ±%.3f →  %.1f %% inside\n",
                rope_grid[i], grid_pct[i]))
  }
  
  cat("\n--  Batch effect (specimen‑level τ², mixed-effects model) -----------------------\n")
  cat(sprintf("τ²_0 (unadjusted)         :  %.4f\n", tau2_0))
  cat(sprintf("τ²_1 (stain‑adjusted)     :  %.4f\n", tau2_1))
  cat(sprintf("Δτ² reduction             :  %.2f %%  →  %s\n",
              delta_tau,
              ifelse(delta_tau < 5, "stain effect negligible", "caution: non‑negligible")))
  cat("============================================================\n\n")
}

# ===========  Bayesian comparison: HE  vs  Gomori
# t statistic (paired)      :  0.626
# Bayes Factor  BF10        :  0.101
# Bayes Factor  BF01        :  9.931   (supports null model (no difference))
# Posterior mean (delta)    :  0.118
# Posterior 95% HDI         : [-0.255 , 0.487]

# --  ROPE decision (Lakens 2022 : |d| < 0.30)
#   ROPE ±0.30 → Posterior mass  82.1 %

# --  ROPE sensitivity analysis 
#   ROPE ±0.250 →  73.0 % inside
# ROPE ±0.300 →  82.1 % inside
# ROPE ±0.400 →  93.1 % inside

# --  Batch effect (specimen‑level τ², mixed-effects model)
#   τ²_0 (unadjusted)         :  165.3584
# τ²_1 (stain‑adjusted)     :  165.3025
# Δτ² reduction             :  0.03 %  →  stain effect negligible


## 6. Visualisations

## 6‑A  posterior density + ROPE
ggplot(data.frame(delta = delta_draws), aes(delta)) +
  geom_density(fill = "#3182BD", alpha = .4) +
  annotate("rect", xmin = -rope_fix, xmax = rope_fix,
           ymin = 0, ymax = Inf,
           fill = "grey70", alpha = .15) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Standardised effect size (delta)",
       y = "Posterior density",
       title = "Posterior of delta (HE – Gomori)",
       subtitle = sprintf("ROPE = ±%.2f  (%.1f%% inside)",
                          rope_fix, rope_pct)) +
  theme_bw(base_size = 13)

## 6‑B  ROPE sensitivity bar chart
ggplot(data.frame(width = factor(rope_grid),
                  pct   = grid_pct),
       aes(width, pct)) +
  geom_col(fill = "#3182BD", alpha = .6, width = .6) +
  geom_text(aes(label = sprintf("%.1f%%", pct)),
            vjust = -0.5, size = 4) +
  labs(x = "ROPE width (Cohen's d)", y = "Posterior mass (%)",
       title = "Sensitivity analysis across ROPE widths") +
  ylim(0, 100) +
  theme_bw(base_size = 12)


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
### -----------------------------### 8. AFTERWORD----------------------------------

# This R script was developed to enhance reproducibility and transparency of
# our statistical analyses for the research manuscript entitled 
# "Myofibre Density Reveals a Critical Threshold around Age 6 in Steroid-Naïve
# Duchenne Muscular Dystrophy: A Retrospective Observational Study" (2025).



### Statistics — Author Contributions and Affiliations(Sapporo,Japan: unless noted)

# Yamakado T., MD., drafted and executed the R-based statistical analyses.
#  1) Department of Cancer Pathology, Faculty of Medicine, Hokkaido University
#  2) Department of Pathology, National Hospital Organization Hokkaido Medical
#     Center
#  3) Center for Neuromuscular Disease, Child Health and Development, 
#     National Hospital Organization Hokkaido Medical Center


# Tanaka S., MD., PhD.: Overall supervision of statistical analysis
#  1) Department of Cancer Pathology, Faculty of Medicine, Hokkaido University
#  2) Department of Surgical Pathology, Hokkaido University Hospital
#  3) Institute for Chemical Reaction Design and Discovery (WPI-ICReDD), 
#     Hokkaido University


# Tanei ZI.(corresponding author), MD., PhD., supervised statistical methodology
#  1) Department of Cancer Pathology, Faculty of Medicine, Hokkaido University, 
#     Sapporo, Japan
#  2) Department of Surgical Pathology, Hokkaido University Hospital
#     Sapporo, Japan
#  Email: tanei@med.hokudai.ac.jp


### The R version and package information used in this analysis are available below:
sessionInfo() 
# R version 4.3.2 (2023-10-31 ucrt)
