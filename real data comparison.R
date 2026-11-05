# ==============================================================================
# Hard-Truncated EWMA Group Sampling Plan: Real Data Empirical Benchmark
# Comparative Analysis: Proposed Plan (K=3) vs. SASP (K=1) vs. Unbounded RASP (K=Inf)
# ==============================================================================

# Install required packages if not available
if (!require("fitdistrplus")) install.packages("fitdistrplus")
library(fitdistrplus)

# ------------------------------------------------------------------------------
# STEP 1: Real-World Industrial Dataset Initialization
# Real dataset: Ball bearing fatigue failure times (in millions of revolutions)
# Source: Lieblein & Zelen (1956) / Standard Reliability Engineering Benchmark
# ------------------------------------------------------------------------------
real_failure_times <- c(
  17.88, 28.92, 33.00, 41.52, 42.12, 45.60, 48.48, 51.84, 
  51.96, 54.12, 55.56, 67.80, 68.64, 68.64, 68.88, 84.12, 
  93.12, 98.64, 105.84, 127.92, 128.04, 173.40, 173.40
)

cat("========================================================\n")
cat("STEP 1: Real Data Lifetime Analysis & Model Fitting\n")
cat("========================================================\n")
cat("Sample Size:", length(real_failure_times), "items\n")

# Fit 2-parameter Weibull distribution using Maximum Likelihood Estimation
fit_weibull <- fitdist(real_failure_times, "weibull", method = "mle")
shape_hat <- fit_weibull$estimate["shape"]  # omega (shape/dispersion parameter)
scale_hat <- fit_weibull$estimate["scale"]  # xi (scale/characteristic lifetime)

# Derive baseline median lifetime M0
median_hat <- scale_hat * (log(2))^(1 / shape_hat)

cat(sprintf("Fitted Weibull Shape (omega): %.4f\n", shape_hat))
cat(sprintf("Fitted Weibull Scale (xi)   : %.4f\n", scale_hat))
cat(sprintf("Estimated Median Life (M0)  : %.4f million revolutions\n\n", median_hat))

# ------------------------------------------------------------------------------
# STEP 2: Acceptance Sampling Probability Metrics Mapping
# Parameter specification under lifetime ratio constraints:
# Producer risk alpha = 0.05, Consumer risk beta = 0.05
# Acceptable quality ratio r2 = M2/M0 = 2.0, Unacceptable r1 = M1/M0 = 1.0
# Test termination ratio d0 = t0 / M0 = 1.0
# ------------------------------------------------------------------------------
alpha <- 0.05
beta  <- 0.05
r1    <- 1.0
r2    <- 2.0
d0    <- 1.0

# Compute unreliability failure fractions under empirical Weibull shape (omega_hat)
p1 <- 1 - exp(-log(2) * (d0 / r1)^shape_hat) # Unacceptable defect level (at r1=1.0)
p2 <- 1 - exp(-log(2) * (d0 / r2)^shape_hat) # Acceptable defect level (at r2=2.0)

cat("Mapped Defect Probabilities for Life Test:\n")
cat(sprintf("Unacceptable Quality Level (p1): %.4f\n", p1))
cat(sprintf("Acceptable Quality Level   (p2): %.4f\n\n", p2))

# ------------------------------------------------------------------------------
# STEP 3: Closed-Form Evaluation Functions
# ------------------------------------------------------------------------------
# Tester group failure probability
p_tester <- function(p, r, c1) {
  1 - pbinom(c1, size = r, prob = p)
}

# Single-cycle Accept (Pa), Reject (Pr), and Indecisive (Pi) Probabilities
calc_single_cycle <- function(g, r, c1, c2, p) {
  pg <- p_tester(p, r, c1)
  Pa <- pbinom(c1, size = g, prob = pg)
  Pr <- 1 - pbinom(c2 - 1, size = g, prob = pg)
  Pi <- 1 - Pa - Pr
  return(list(Pa = Pa, Pr = Pr, Pi = Pi))
}

# Finite K-Truncated Operating Characteristic (OC) and ASN
calc_truncated_metrics <- function(g, r, c1, c2, p, K) {
  probs <- calc_single_cycle(g, r, c1, c2, p)
  Pa <- probs$Pa
  Pr <- probs$Pr
  Pi <- probs$Pi
  
  if (is.infinite(K)) {
    OC  <- Pa / (Pa + Pr)
    ASN <- (g * r) / (Pa + Pr)
  } else if (K == 1) {
    OC  <- Pa
    ASN <- g * r
  } else {
    if (abs(Pi - 1) < 1e-10) {
      OC  <- Pa * K
      ASN <- (g * r) * K
    } else {
      factor <- (1 - Pi^K) / (1 - Pi)
      OC     <- Pa * factor
      ASN    <- (g * r) * factor
    }
  }
  return(list(OC = OC, ASN = ASN, n_single = g * r))
}

# ------------------------------------------------------------------------------
# STEP 4: Grid Optimization Engine for Sampling Plans
# ------------------------------------------------------------------------------
find_optimal_plan <- function(K_val, p1, p2, alpha_val, beta_val) {
  best_asn <- Inf
  best_params <- NULL
  
  for (g in 2:20) {
    for (r in 2:10) {
      for (c1 in 0:(g - 1)) {
        for (c2 in (c1 + 1):g) {
          metrics_p1 <- calc_truncated_metrics(g, r, c1, c2, p1, K_val)
          metrics_p2 <- calc_truncated_metrics(g, r, c1, c2, p2, K_val)
          
          # Check dual risk constraints: OC(p1) <= beta and OC(p2) >= 1 - alpha
          if (metrics_p1$OC <= beta_val && metrics_p2$OC >= (1 - alpha_val)) {
            if (metrics_p1$ASN < best_asn) {
              best_asn <- metrics_p1$ASN
              best_params <- list(
                g = g, r = r, c1 = c1, c2 = c2,
                ASN_p1 = metrics_p1$ASN,
                P_accept_p2 = metrics_p2$OC,
                P_accept_p1 = metrics_p1$OC,
                n_single = g * r,
                max_risk = ifelse(is.infinite(K_val), "Unbounded (Inf)", as.character(g * r * K_val))
              )
            }
          }
        }
      }
    }
  }
  return(best_params)
}

# ------------------------------------------------------------------------------
# STEP 5: Execute Comparative Analysis Across Plan Architectures
# ------------------------------------------------------------------------------
cat("========================================================\n")
cat("STEP 2: Optimization Results on Real Ball Bearing Data\n")
cat("========================================================\n")

plan_sasp     <- find_optimal_plan(K_val = 1,        p1, p2, alpha, beta)
plan_unbounded<- find_optimal_plan(K_val = Inf,      p1, p2, alpha, beta)
plan_proposed <- find_optimal_plan(K_val = 3,        p1, p2, alpha, beta)

# Compute percent ASN savings relative to SASP
savings_unbounded <- ((plan_sasp$ASN_p1 - plan_unbounded$ASN_p1) / plan_sasp$ASN_p1) * 100
savings_proposed  <- ((plan_sasp$ASN_p1 - plan_proposed$ASN_p1)  / plan_sasp$ASN_p1) * 100

comparison_df <- data.frame(
  Architecture = c("Single Group (SASP)", "Unbounded RASP", "Proposed Hard-Truncated Plan"),
  Truncation_K = c("1", "Inf", "3"),
  Optimal_g_r  = c(sprintf("(%d, %d)", plan_sasp$g, plan_sasp$r),
                   sprintf("(%d, %d)", plan_unbounded$g, plan_unbounded$r),
                   sprintf("(%d, %d)", plan_proposed$g, plan_proposed$r)),
  Optimal_c1_c2= c(sprintf("(%d, %d)", plan_sasp$c1, plan_sasp$c2),
                   sprintf("(%d, %d)", plan_unbounded$c1, plan_unbounded$c2),
                   sprintf("(%d, %d)", plan_proposed$c1, plan_proposed$c2)),
  ASN_p1       = c(round(plan_sasp$ASN_p1, 2), round(plan_unbounded$ASN_p1, 2), round(plan_proposed$ASN_p1, 2)),
  Max_Risk_nK  = c(as.character(plan_sasp$n_single), "Unbounded (Inf)", as.character(plan_proposed$n_single * 3)),
  ASN_Savings  = c("0.00%", sprintf("%.2f%%", savings_unbounded), sprintf("%.2f%%", savings_proposed))
)

print(comparison_df)

# ------------------------------------------------------------------------------
# STEP 6: Graphical Comparison across Defect Spectrum (using windows())
# ------------------------------------------------------------------------------
p_seq <- seq(0.001, 0.60, length.out = 300)

metrics_sasp_curve <- t(sapply(p_seq, function(p) {
  res <- calc_truncated_metrics(plan_sasp$g, plan_sasp$r, plan_sasp$c1, plan_sasp$c2, p, 1)
  c(OC = res$OC, ASN = res$ASN)
}))

metrics_rasp_curve <- t(sapply(p_seq, function(p) {
  res <- calc_truncated_metrics(plan_unbounded$g, plan_unbounded$r, plan_unbounded$c1, plan_unbounded$c2, p, Inf)
  c(OC = res$OC, ASN = res$ASN)
}))

metrics_prop_curve <- t(sapply(p_seq, function(p) {
  res <- calc_truncated_metrics(plan_proposed$g, plan_proposed$r, plan_proposed$c1, plan_proposed$c2, p, 3)
  c(OC = res$OC, ASN = res$ASN)
}))

# Open Windows graphics device
windows(width = 11, height = 5.5, title = "Comparative Efficiency Analysis - Real Data")

# Configure 1x2 panel layout
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1), cex.main = 1.05, font.main = 2)

col_sasp <- "#E41A1C" # Red
col_rasp <- "#377EB8" # Blue
col_prop <- "#4DAF4A" # Green

# Left Panel: ASN(p) Curves
plot(p_seq, metrics_sasp_curve[, "ASN"], type = "l", lty = 2, lwd = 2.2, col = col_sasp,
     ylim = c(0, max(metrics_sasp_curve[, "ASN"]) * 1.15),
     xlab = expression(paste("Lot Failure Probability (", p, ")")),
     ylab = expression(paste("Average Sample Number (", ASN(p), ")")),
     main = "")
grid(col = "lightgray", lty = "dotted")

lines(p_seq, metrics_rasp_curve[, "ASN"], lty = 3, lwd = 2.2, col = col_rasp)
lines(p_seq, metrics_prop_curve[, "ASN"], lty = 1, lwd = 2.5, col = col_prop)

legend("topright", 
       legend = c(sprintf("Single Group (SASP, K=1)"), 
                  sprintf("Unbounded RASP (K=\u221E)"), 
                  sprintf("Proposed Plan (K=3)")),
       col = c(col_sasp, col_rasp, col_prop),
       lty = c(2, 3, 1), lwd = c(2.2, 2.2, 2.5),
       bty = "o", bg = "white", cex = 0.85)

# Right Panel: OC(p) Curves
plot(p_seq, metrics_sasp_curve[, "OC"], type = "l", lty = 2, lwd = 2.2, col = col_sasp,
     ylim = c(0, 1.05),
     xlab = expression(paste("Lot Failure Probability (", p, ")")),
     ylab = expression(paste("Acceptance Probability (", OC(p), ")")),
     main = "")
grid(col = "lightgray", lty = "dotted")

lines(p_seq, metrics_rasp_curve[, "OC"], lty = 3, lwd = 2.2, col = col_rasp)
lines(p_seq, metrics_prop_curve[, "OC"], lty = 1, lwd = 2.5, col = col_prop)

legend("bottomleft", 
       legend = c("Single Group (SASP, K=1)", 
                  "Unbounded RASP (K=\u221E)", 
                  "Proposed Plan (K=3)"),
       col = c(col_sasp, col_rasp, col_prop),
       lty = c(2, 3, 1), lwd = c(2.2, 2.2, 2.5),
       bty = "o", bg = "white", cex = 0.85)

# Save high-resolution PNG file for manuscript insertion
dev.copy(png, filename = "Real_Data_Comparative_Curves.png", width = 11, height = 5.5, units = "in", res = 300)
dev.off()