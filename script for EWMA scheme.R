# ==============================================================================
# Hard-Truncated EWMA Group Acceptance Sampling Scheme for Weibull Models
# Fully Corrected & Comprehensive R Execution Engine
# Theoretical Specifications based on Obulezi et al.
# ==============================================================================

suppressPackageStartupMessages({
  if (!require("ggplot2")) install.packages("ggplot2", quietly = TRUE)
  library(ggplot2)
})

output_dir <- getwd()

# ------------------------------------------------------------------------------
# 1. CORE CLOSED-FORM MATHEMATICAL METRICS (THEOREMS 2.1 - 2.6)
# ------------------------------------------------------------------------------

# Theorem 2.1: Hard-Truncated Acceptance Probability P_accept^(K)(p)
oc_k <- function(Pa, Pr, Pi, K) {
  Pa <- as.numeric(Pa); Pr <- as.numeric(Pr); Pi <- as.numeric(Pi)
  if (is.infinite(K)) {
    if ((Pa + Pr) < 1e-12) return(0)
    return(Pa / (Pa + Pr))
  }
  if (abs(Pi) < 1e-12) return(Pa)
  if (abs(Pi - 1) < 1e-12) return(0)
  Pa * (1 - Pi^K) / (1 - Pi)
}

# Theorem 2.2: Hard-Truncated Rejection Probability P_reject^(K)(p)
reject_k <- function(Pa, Pr, Pi, K) {
  Pa <- as.numeric(Pa); Pr <- as.numeric(Pr); Pi <- as.numeric(Pi)
  if (is.infinite(K)) {
    if ((Pa + Pr) < 1e-12) return(1)
    return(Pr / (Pa + Pr))
  }
  if (abs(Pi) < 1e-12) return(Pr)
  if (abs(Pi - 1) < 1e-12) return(1)
  Pr * (1 - Pi^K) / (1 - Pi) + Pi^K
}

# Theorem 2.3: Average Sample Number ASN^(K)(p)
asn_k <- function(n_stage, Pi, K) {
  n_stage <- as.numeric(n_stage); Pi <- as.numeric(Pi)
  if (is.infinite(K)) {
    if (abs(1 - Pi) < 1e-12) return(Inf)
    return(n_stage / (1 - Pi))
  }
  if (abs(Pi) < 1e-12) return(n_stage)
  n_stage * (1 - Pi^K) / (1 - Pi)
}

# Theorem 2.4: Variance of Inspection Burden Var^(K)(ASN(p))
var_asn_k <- function(n_stage, Pi, K) {
  n_stage <- as.numeric(n_stage); Pi <- as.numeric(Pi)
  if (is.infinite(K)) {
    if (abs(1 - Pi) < 1e-12) return(Inf)
    return((n_stage^2) * Pi / ((1 - Pi)^2))
  }
  if (abs(Pi) < 1e-12) return(0)
  num <- Pi * (1 + Pi^K - 2 * K * Pi^(K - 1) * (1 - Pi)) - Pi^(2 * K)
  den <- (1 - Pi)^2
  (n_stage^2) * (num / den)
}

# Theorem 2.5: Average Total Inspection ATI^(K)(p)
ati_k <- function(N_lot, n_stage, Pa, Pr, Pi, K) {
  ASN_val <- asn_k(n_stage, Pi, K)
  P_accept_val <- oc_k(Pa, Pr, Pi, K)
  ASN_val + (N_lot - ASN_val) * (1 - P_accept_val)
}

# Theorem 2.6: Average Outgoing Quality AOQ^(K)(p)
aoq_k <- function(p_defect, N_lot, n_stage, Pa, Pr, Pi, K) {
  p_defect <- as.numeric(p_defect)
  ASN_val <- asn_k(n_stage, Pi, K)
  P_accept_val <- oc_k(Pa, Pr, Pi, K)
  p_defect * P_accept_val * (N_lot - ASN_val) / N_lot
}

# ------------------------------------------------------------------------------
# 2. BIVARIATE WEIBULL GROUP PROBABILITIES (DEFINITION 2.2 & SECTION 4.1)
# ------------------------------------------------------------------------------

# Compute Weibull unreliability p(r) given median ratio r and termination ratio d0
weibull_p <- function(r_ratio, d0, omega = 1.5) {
  r_ratio <- as.numeric(r_ratio); d0 <- as.numeric(d0); omega <- as.numeric(omega)
  as.numeric(1 - exp(-log(2) * (d0 / r_ratio)^omega))
}

# Trinomial group failure probabilities under hierarchical location/dispersion rules
compute_weibull_group_probs <- function(p_defect, g, r, c1, c2, d1, d2, omega_val = 1.5) {
  p_defect <- as.numeric(p_defect)
  g <- as.integer(g); r <- as.integer(r)
  c1 <- as.integer(c1); c2 <- as.integer(c2)
  d1 <- as.integer(d1); d2 <- as.integer(d2)
  
  # 1. Location failure probability (pM): probability group exceeds item failure limit c_g
  pM <- as.numeric(1 - pbinom(c1, size = r, prob = p_defect))
  
  # 2. Conditional dispersion failure probability (p_omega): evaluated on groups meeting location standard
  # Calculated via shape parameter variance shift on truncated observations
  disp_factor <- pmin(0.99, pmax(0.01, p_defect * (1.5 / as.numeric(omega_val))))
  p_omega <- as.numeric((1 - pM) * (1 - pbinom(d1, size = r, prob = disp_factor)))
  
  p_conf <- max(0, 1 - pM - p_omega)
  prob_vec <- c(pM, p_omega, p_conf)
  prob_vec <- prob_vec / sum(prob_vec) # Conserve trinomial probability measure
  
  # Single-cycle acceptance probability Pa with upper index bounds min(d1, g - yM)
  Pa <- 0
  for (yM in 0:c1) {
    max_yw <- min(d1, g - yM)
    if (max_yw >= 0) {
      for (yw in 0:max_yw) {
        Pa <- Pa + dmultinom(c(yM, yw, g - yM - yw), prob = prob_vec)
      }
    }
  }
  
  # Single-cycle non-rejection probability (1 - Pr) with upper index bounds min(d2, g - yM)
  P_not_reject <- 0
  for (yM in 0:c2) {
    max_yw <- min(d2, g - yM)
    if (max_yw >= 0) {
      for (yw in 0:max_yw) {
        P_not_reject <- P_not_reject + dmultinom(c(yM, yw, g - yM - yw), prob = prob_vec)
      }
    }
  }
  
  Pr <- max(0, 1 - P_not_reject)
  Pi <- max(0, 1 - Pa - Pr)
  
  list(Pa = Pa, Pr = Pr, Pi = Pi)
}

# ------------------------------------------------------------------------------
# 3. OPTIMIZATION SOLVERS (CLASSICAL TWO-POINT & WEIGHTED RISK)
# ------------------------------------------------------------------------------

# Classical (alpha, beta) Optimization Model (Section 3.1)
solve_weibull_two_point <- function(d0 = 0.5, omega = 1.5, r1 = 1.0, r2 = 2.0, 
                                    alpha = 0.05, beta = 0.05, K_max = 3) {
  p1 <- weibull_p(r1, d0, omega) # Unacceptable quality level
  p2 <- weibull_p(r2, d0, omega) # Acceptable quality level
  
  best_plan <- NULL
  min_asn <- Inf
  
  for (g in 2:30) {
    for (r in 2:10) {
      n_stage <- g * r
      for (c1 in 0:min(4, g - 2)) {
        for (c2 in (c1 + 1):min(c1 + 4, g - 1)) {
          for (d1 in 0:min(4, g - 2)) {
            for (d2 in (d1 + 1):min(d1 + 4, g - 1)) {
              
              probs1 <- compute_weibull_group_probs(p1, g, r, c1, c2, d1, d2, omega)
              probs2 <- compute_weibull_group_probs(p2, g, r, c1, c2, d1, d2, omega)
              
              P_accept_p1 <- oc_k(probs1$Pa, probs1$Pr, probs1$Pi, K_max)
              P_accept_p2 <- oc_k(probs2$Pa, probs2$Pr, probs2$Pi, K_max)
              
              # Risk limits: Consumer risk P_accept(p1) <= beta and Producer risk 1 - P_accept(p2) <= alpha
              if (P_accept_p1 <= beta && (1 - P_accept_p2) <= alpha) {
                current_asn <- asn_k(n_stage, probs1$Pi, K_max)
                if (current_asn < min_asn) {
                  min_asn <- current_asn
                  best_plan <- list(
                    g = g, r = r, n_stage = n_stage,
                    c1 = c1, c2 = c2, d1 = d1, d2 = d2, K = K_max,
                    p1 = p1, p2 = p2, P_accept_p2 = P_accept_p2, ASN_p1 = current_asn
                  )
                }
              }
            }
          }
        }
      }
    }
  }
  return(best_plan)
}

# Adaptive Weighted-Risk (gamma-WR) Optimization Model (Section 3.2)
solve_weibull_weighted_risk <- function(d0 = 1.0, omega = 1.5, r1 = 1.0, r2 = 2.0,
                                        w0 = 0.5, w1 = 0.5, gamma_limit = 0.05, 
                                        Vmax = 5000, K_max = 3) {
  p1 <- weibull_p(r1, d0, omega)
  p2 <- weibull_p(r2, d0, omega)
  
  best_plan <- NULL
  min_asn <- Inf
  
  for (g in 2:25) {
    for (r in 2:6) {
      n_stage <- g * r
      for (c1 in 0:min(3, g - 2)) {
        for (c2 in (c1 + 1):min(c1 + 3, g - 1)) {
          for (d1 in 0:min(3, g - 2)) {
            for (d2 in (d1 + 1):min(d1 + 3, g - 1)) {
              
              probs1 <- compute_weibull_group_probs(p1, g, r, c1, c2, d1, d2, omega)
              probs2 <- compute_weibull_group_probs(p2, g, r, c1, c2, d1, d2, omega)
              
              P_accept_p1 <- oc_k(probs1$Pa, probs1$Pr, probs1$Pi, K_max)
              P_accept_p2 <- oc_k(probs2$Pa, probs2$Pr, probs2$Pi, K_max)
              
              WR <- w0 * (1 - P_accept_p2) + w1 * P_accept_p1
              var_asn_val <- var_asn_k(n_stage, probs1$Pi, K_max)
              
              if (WR <= gamma_limit && var_asn_val <= Vmax) {
                current_asn <- asn_k(n_stage, probs1$Pi, K_max)
                if (current_asn < min_asn) {
                  min_asn <- current_asn
                  best_plan <- list(
                    g = g, r = r, n_stage = n_stage,
                    c1 = c1, c2 = c2, d1 = d1, d2 = d2, K = K_max,
                    WR = WR, ASN_p1 = current_asn, P_accept_p2 = P_accept_p2
                  )
                }
              }
            }
          }
        }
      }
    }
  }
  return(best_plan)
}

# ------------------------------------------------------------------------------
# 4. EWMA PROCESS MEMORY CONTROLLER SIMULATION (SECTION 2.4)
# ------------------------------------------------------------------------------

simulate_weibull_ewma <- function(num_lots = 50, p0 = 0.05, lambda_ewma = 0.2, L_width = 3.0, n_sample = 20) {
  set.seed(123)
  
  sigma_z <- sqrt((p0 * (1 - p0) / n_sample) * (lambda_ewma / (2 - lambda_ewma)))
  UCL <- p0 + L_width * sigma_z
  LCL <- max(0, p0 - L_width * sigma_z)
  
  Z <- numeric(num_lots)
  p_observed <- numeric(num_lots)
  K_adaptive <- numeric(num_lots)
  
  # Simulated process shift starting at Lot 25
  true_p <- c(rep(0.04, 24), rep(0.12, 26))
  
  Z_prev <- p0
  for (i in 1:num_lots) {
    p_obs <- rbinom(1, size = n_sample, prob = true_p[i]) / n_sample
    p_observed[i] <- p_obs
    
    Z_curr <- lambda_ewma * p_obs + (1 - lambda_ewma) * Z_prev
    Z[i] <- Z_curr
    
    if (Z_curr < LCL) {
      K_adaptive[i] <- 5 # State 1: In-Control Good Quality (High Resampling Allowed)
    } else if (Z_curr <= UCL) {
      K_adaptive[i] <- 3 # State 2: Normal Quality (Standard Truncation)
    } else {
      K_adaptive[i] <- 1 # State 3: Process Shift / Poor Quality (Immediate Hard Stop)
      Z_curr <- p0       # Reset controller post corrective action
    }
    Z_prev <- Z_curr
  }
  
  data.frame(
    Lot = 1:num_lots, True_p = true_p, Observed_p = p_observed,
    EWMA_Z = Z, LCL = LCL, UCL = UCL, Truncation_K = K_adaptive
  )
}

# ------------------------------------------------------------------------------
# 5. AUTOMATED LATEX TABLE GENERATORS
# ------------------------------------------------------------------------------

export_weibull_tables <- function() {
  cat("\n=================================================================\n")
  cat("TABLE 1: CLASSICAL (ALPHA, BETA) WEIBULL GROUP SAMPLING PLANS\n")
  cat("=================================================================\n\n")
  
  beta_list <- c(0.25, 0.10, 0.05, 0.01)
  k_list <- c(2, 3, Inf)
  
  for (beta_val in beta_list) {
    for (k_val in k_list) {
      res_05 <- solve_weibull_two_point(d0 = 0.5, omega = 1.5, alpha = 0.05, beta = beta_val, K_max = k_val)
      res_10 <- solve_weibull_two_point(d0 = 1.0, omega = 1.5, alpha = 0.05, beta = beta_val, K_max = k_val)
      
      k_str <- if (is.infinite(k_val)) "\\infty" else as.character(k_val)
      
      str_05_params  <- if (!is.null(res_05)) sprintf("(%d, %d)", res_05$g, res_05$r) else "--"
      str_05_cutoffs <- if (!is.null(res_05)) sprintf("(%d, %d)", res_05$c1, res_05$c2) else "--"
      str_05_pacc    <- if (!is.null(res_05)) sprintf("%.4f", res_05$P_accept_p2) else "--"
      str_05_asn     <- if (!is.null(res_05)) sprintf("%.2f", res_05$ASN_p1) else "--"
      
      str_10_params  <- if (!is.null(res_10)) sprintf("(%d, %d)", res_10$g, res_10$r) else "--"
      str_10_cutoffs <- if (!is.null(res_10)) sprintf("(%d, %d)", res_10$c1, res_10$c2) else "--"
      str_10_pacc    <- if (!is.null(res_10)) sprintf("%.4f", res_10$P_accept_p2) else "--"
      str_10_asn     <- if (!is.null(res_10)) sprintf("%.2f", res_10$ASN_p1) else "--"
      
      cat(sprintf("%.2f & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\\n",
                  beta_val, k_str,
                  str_05_params, str_05_cutoffs, str_05_pacc, str_05_asn,
                  str_10_params, str_10_cutoffs, str_10_pacc, str_10_asn))
    }
    cat("[0.5ex]\n")
  }
}

# Execute Table Exporter and Plots
export_weibull_tables()

# ------------------------------------------------------------------------------
# 6. DIAGNOSTIC PLOTS GENERATION
# ------------------------------------------------------------------------------

opt_demo <- solve_weibull_two_point(d0 = 0.5, omega = 1.5, alpha = 0.05, beta = 0.05, K_max = 3)
p_seq <- seq(0.001, 0.25, length.out = 150)

oc_vals <- sapply(p_seq, function(p) {
  probs <- compute_weibull_group_probs(p, g = opt_demo$g, r = opt_demo$r, 
                                       c1 = opt_demo$c1, c2 = opt_demo$c2, 
                                       d1 = opt_demo$d1, d2 = opt_demo$d2)
  oc_k(probs$Pa, probs$Pr, probs$Pi, K = opt_demo$K)
})

asn_vals <- sapply(p_seq, function(p) {
  probs <- compute_weibull_group_probs(p, g = opt_demo$g, r = opt_demo$r, 
                                       c1 = opt_demo$c1, c2 = opt_demo$c2, 
                                       d1 = opt_demo$d1, d2 = opt_demo$d2)
  asn_k(n_stage = opt_demo$n_stage, probs$Pi, K = opt_demo$K)
})

plot_df <- data.frame(
  p = rep(p_seq, 2),
  Value = c(oc_vals, asn_vals),
  Metric = rep(c("Operating Characteristic OC(p)", "Average Sample Number ASN(p)"), each = 150)
)

p1 <- ggplot(plot_df, aes(x = p, y = Value, color = Metric)) +
  geom_line(size = 1.2) +
  facet_wrap(~Metric, scales = "free_y") +
  theme_bw() +
  labs(title = "Weibull Group Sampling Performance Diagnostic Curves",
       subtitle = paste0("Optimal Plan: g = ", opt_demo$g, ", r = ", opt_demo$r, ", K = ", opt_demo$K),
       x = "Lot Failure Probability (p)", y = "Metric Value") +
  theme(plot.title = element_text(face = "bold"), legend.position = "none")

ggsave(file.path(output_dir, "Weibull_OC_ASN_Curves.png"), p1, width = 8, height = 4, dpi = 300)

ewma_df <- simulate_weibull_ewma()
p2 <- ggplot(ewma_df, aes(x = Lot)) +
  geom_line(aes(y = EWMA_Z, color = "EWMA Statistic Z_i"), size = 1) +
  geom_point(aes(y = EWMA_Z, color = "EWMA Statistic Z_i"), size = 2) +
  geom_hline(aes(yintercept = UCL, linetype = "UCL"), color = "red", size = 0.8) +
  geom_hline(aes(yintercept = LCL, linetype = "LCL"), color = "blue", size = 0.8) +
  geom_col(aes(y = Truncation_K / 50, fill = "Adaptive K Depth"), alpha = 0.2) +
  scale_y_continuous(sec.axis = sec_axis(~ . * 50, name = "Truncation Depth (K)")) +
  theme_light() +
  labs(title = "EWMA Process-Memory Control Chart under Weibull Lot Process",
       x = "Submitted Lot Sequence", y = "EWMA Statistic (Z_i)",
       color = "Metric", linetype = "Control Limits", fill = "Process State") +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

ggsave(file.path(output_dir, "Weibull_EWMA_Tracking.png"), p2, width = 8, height = 5, dpi = 300)

cat("\nExecution completed. All functions updated and diagnostic figures saved to working directory.\n")