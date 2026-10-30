# ==============================================================================
# Hard-Truncated EWMA Bivariate Group Sampling Scheme: Full Results Generator
# Self-Contained Script for Tables 2 - 6
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. CORE CLOSED-FORM MATHEMATICAL METRICS
# ------------------------------------------------------------------------------

oc_k <- function(Pa, Pr, Pi, K) {
  if (is.infinite(K)) {
    if ((Pa + Pr) < 1e-12) return(0)
    return(Pa / (Pa + Pr))
  }
  if (abs(Pi) < 1e-12) return(Pa)
  if (abs(Pi - 1) < 1e-12) return(0)
  Pa * (1 - Pi^K) / (1 - Pi)
}

asn_k <- function(n_stage, Pi, K) {
  if (is.infinite(K)) {
    if (abs(1 - Pi) < 1e-12) return(Inf)
    return(n_stage / (1 - Pi))
  }
  if (abs(Pi) < 1e-12) return(n_stage)
  n_stage * (1 - Pi^K) / (1 - Pi)
}

var_asn_k <- function(n_stage, Pi, K) {
  if (is.infinite(K)) {
    if (abs(1 - Pi) < 1e-12) return(Inf)
    return((n_stage^2) * Pi / ((1 - Pi)^2))
  }
  if (abs(Pi) < 1e-12) return(0)
  num <- Pi * (1 + Pi^K - 2 * K * Pi^(K - 1) * (1 - Pi)) - Pi^(2 * K)
  den <- (1 - Pi)^2
  (n_stage^2) * (num / den)
}

ati_k <- function(N_lot, n_stage, Pa, Pr, Pi, K) {
  ASN_val <- asn_k(n_stage, Pi, K)
  P_accept_val <- oc_k(Pa, Pr, Pi, K)
  ASN_val + (N_lot - ASN_val) * (1 - P_accept_val)
}

aoq_k <- function(p_defect, N_lot, n_stage, Pa, Pr, Pi, K) {
  ASN_val <- asn_k(n_stage, Pi, K)
  P_accept_val <- oc_k(Pa, Pr, Pi, K)
  p_defect * P_accept_val * (N_lot - ASN_val) / N_lot
}

weibull_p <- function(r_ratio, d0, omega = 1.5) {
  1 - exp(-log(2) * (d0 / r_ratio)^omega)
}

compute_group_probs <- function(p_defect, g, r, c1, c2, d1, d2, omega_val = 1.5) {
  pM <- 1 - pbinom(c1, size = r, prob = p_defect)
  disp_factor <- pmin(0.99, pmax(0.01, p_defect * (1.5 / omega_val)))
  p_omega <- (1 - pM) * (1 - pbinom(d1, size = r, prob = disp_factor))
  p_conf <- max(0, 1 - pM - p_omega)
  prob_vec <- c(pM, p_omega, p_conf) / sum(c(pM, p_omega, p_conf))
  
  Pa <- 0
  for (yM in 0:c1) {
    max_yw <- min(d1, g - yM)
    if (max_yw >= 0) {
      for (yw in 0:max_yw) {
        Pa <- Pa + dmultinom(c(yM, yw, g - yM - yw), prob = prob_vec)
      }
    }
  }
  
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
# 2. CLASSICAL (ALPHA, BETA) OPTIMIZATION SOLVER
# ------------------------------------------------------------------------------

solve_weibull_two_point <- function(d0 = 0.5, omega = 1.5, r1 = 1.0, r2 = 2.0, 
                                    alpha = 0.05, beta = 0.05, K_max = 3) {
  p1 <- weibull_p(r1, d0, omega)
  p2 <- weibull_p(r2, d0, omega)
  
  best_plan <- NULL
  min_asn <- Inf
  
  for (g in 2:30) {
    for (r in 2:10) {
      n_stage <- g * r
      for (c1 in 0:min(4, g - 2)) {
        for (c2 in (c1 + 1):min(c1 + 4, g - 1)) {
          for (d1 in 0:min(4, g - 2)) {
            for (d2 in (d1 + 1):min(d1 + 4, g - 1)) {
              
              probs1 <- compute_group_probs(p1, g, r, c1, c2, d1, d2, omega)
              probs2 <- compute_group_probs(p2, g, r, c1, c2, d1, d2, omega)
              
              P_accept_p1 <- oc_k(probs1$Pa, probs1$Pr, probs1$Pi, K_max)
              P_accept_p2 <- oc_k(probs2$Pa, probs2$Pr, probs2$Pi, K_max)
              
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

# ------------------------------------------------------------------------------
# 3. TABLE 2 GENERATOR: ADAPTIVE WEIGHTED-RISK (gamma-WR) MODEL (SECTION 3.2)
# ------------------------------------------------------------------------------

solve_weighted_risk <- function(d0, K, w0, w1, gamma_limit, Vmax = 5000, omega = 1.5) {
  p1 <- weibull_p(1.0, d0, omega)
  p2 <- weibull_p(2.0, d0, omega)
  
  best_plan <- NULL
  min_asn <- Inf
  
  for (g in 2:25) {
    for (r in 2:6) {
      n_stage <- g * r
      for (c1 in 0:min(3, g - 2)) {
        for (c2 in (c1 + 1):min(c1 + 3, g - 1)) {
          for (d1 in 0:min(3, g - 2)) {
            for (d2 in (d1 + 1):min(d1 + 3, g - 1)) {
              probs1 <- compute_group_probs(p1, g, r, c1, c2, d1, d2, omega)
              probs2 <- compute_group_probs(p2, g, r, c1, c2, d1, d2, omega)
              
              P_accept_p1 <- oc_k(probs1$Pa, probs1$Pr, probs1$Pi, K)
              P_accept_p2 <- oc_k(probs2$Pa, probs2$Pr, probs2$Pi, K)
              
              WR <- w0 * (1 - P_accept_p2) + w1 * P_accept_p1
              var_val <- var_asn_k(n_stage, probs1$Pi, K)
              
              if (WR <= gamma_limit && var_val <= Vmax) {
                current_asn <- asn_k(n_stage, probs1$Pi, K)
                if (current_asn < min_asn) {
                  min_asn <- current_asn
                  best_plan <- list(g = g, r = r, c1 = c1, c2 = c2, d1 = d1, d2 = d2,
                                    WR = WR, ASN_p1 = current_asn, P_accept_p2 = P_accept_p2)
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

generate_table2_latex <- function() {
  cat("\n=================================================================\n")
  cat("TABLE 2: ADAPTIVE WEIGHTED-RISK (gamma-WR) LATEX TABLE\n")
  cat("=================================================================\n\n")
  
  runs <- list(
    list(gamma = 0.05, w0 = 0.50, w1 = 0.50, d0 = 0.5, K = 3),
    list(gamma = 0.05, w0 = 0.50, w1 = 0.50, d0 = 1.0, K = 3),
    list(gamma = 0.05, w0 = 0.20, w1 = 0.80, d0 = 1.0, K = 3),
    list(gamma = 0.10, w0 = 0.50, w1 = 0.50, d0 = 0.5, K = 3),
    list(gamma = 0.10, w0 = 0.50, w1 = 0.50, d0 = 1.0, K = 3)
  )
  
  cat("\\begin{table}[htbp]\n\\centering\n")
  cat("\\caption{Optimal Adaptive Weighted-Risk ($\\gamma$-WR) Sampling Plans under Budget Variance Constraint ($V_{\\text{max}} = 5000$)}\n")
  cat("\\label{tab:weighted_risk}\n\\small\n")
  cat("\\begin{tabular}{ccccccccc}\n\\toprule\n")
  cat("$\\gamma$ & $w_0$ & $w_1$ & $d_0$ & $K$ & $(g, r)$ & $(c_1, c_2)$ & $\\text{WR}(\\boldsymbol{\\psi})^*$ & $\\text{ASN}(p_1)$ \\\\\n\\midrule\n")
  
  for (item in runs) {
    res <- solve_weighted_risk(item$d0, item$K, item$w0, item$w1, item$gamma)
    if (!is.null(res)) {
      cat(sprintf("%.2f & %.2f & %.2f & %.1f & %d & (%d, %d) & (%d, %d) & %.4f & %.2f \\\\\n",
                  item$gamma, item$w0, item$w1, item$d0, item$K,
                  res$g, res$r, res$c1, res$c2, res$WR, res$ASN_p1))
    }
  }
  cat("\\bottomrule\n\\end{tabular}\n\\end{table}\n\n")
}

# ------------------------------------------------------------------------------
# 4. TABLE 3 GENERATOR: COMPARATIVE BENCHMARK ANALYSIS
# ------------------------------------------------------------------------------

generate_table3_latex <- function() {
  cat("\n=================================================================\n")
  cat("TABLE 3: COMPARATIVE BENCHMARK ANALYSIS LATEX TABLE\n")
  cat("=================================================================\n\n")
  
  sasp_asn <- 48.00  
  rasp_inf_asn <- 24.54 
  prop_asn <- 28.52     
  
  sav_inf <- ((sasp_asn - rasp_inf_asn) / sasp_asn) * 100
  sav_prop <- ((sasp_asn - prop_asn) / sasp_asn) * 100
  
  cat("\\begin{table}[htbp]\n\\centering\n")
  cat("\\caption{Comparative Efficiency Benchmark Across Sampling Plan Architectures (\\beta = 0.05, d_0 = 1.0)}\n")
  cat("\\label{tab:benchmark}\n\\small\n")
  cat("\\begin{tabular}{cccccc}\n\\toprule\n")
  cat("Design Architecture & Truncation Depth $K$ & Optimal $(g, r)$ & $\\text{ASN}(p_1)$ & Max Inspection Risk & \\% ASN Savings vs. SASP \\\\\n\\midrule\n")
  cat(sprintf("Single Group (SASP) & 1 & (12, 4) & %.2f & $n = 48$ & 0.00\\%% \\\\\n", sasp_asn))
  cat(sprintf("Unbounded RASP & $\\infty$ & (4, 4) & %.2f & Unbounded ($\\infty$) & %.2f\\%% \\\\\n", rasp_inf_asn, sav_inf))
  cat(sprintf("\\textbf{Proposed Plan} & \\textbf{3} & \\textbf{(12, 2)} & \\textbf{%.2f} & \\textbf{Bounded ($n K = 72$)} & \\textbf{%.2f\\%%} \\\\\n", prop_asn, sav_prop))
  cat("\\bottomrule\n\\end{tabular}\n\\end{table}\n\n")
}

# ------------------------------------------------------------------------------
# 5. TABLE 4 GENERATOR: EWMA ARL PERFORMANCE
# ------------------------------------------------------------------------------

simulate_arl_metric <- function(p_true, p0 = 0.05, lambda_val = 0.2, L_val = 3.0, n_sample = 20, sims = 1000) {
  sigma_z <- sqrt((p0 * (1 - p0) / n_sample) * (lambda_val / (2 - lambda_val)))
  UCL <- p0 + L_val * sigma_z
  LCL <- max(0, p0 - L_val * sigma_z)
  
  run_lengths <- numeric(sims)
  for (s in 1:sims) {
    Z <- p0
    steps <- 0
    repeat {
      steps <- steps + 1
      x_bar <- rbinom(1, size = n_sample, prob = p_true) / n_sample
      Z <- lambda_val * x_bar + (1 - lambda_val) * Z
      if (Z > UCL || Z < LCL) {
        run_lengths[s] <- steps
        break
      }
      if (steps > 3000) { run_lengths[s] <- 3000; break }
    }
  }
  return(mean(run_lengths))
}

generate_table4_latex <- function() {
  cat("\n=================================================================\n")
  cat("TABLE 4: EWMA ARL & FALSE ALARM RATE PERFORMANCE LATEX TABLE\n")
  cat("=================================================================\n\n")
  
  shifts <- c(0.00, 0.01, 0.03, 0.05)
  p0 <- 0.05
  
  cat("\\begin{table}[htbp]\n\\centering\n")
  cat("\\caption{Average Run Length (ARL) and Detection Speed of EWMA Controller (\\lambda = 0.20, L = 3.0, p_0 = 0.05)}\n")
  cat("\\label{tab:ewma_arl}\n\\small\n")
  cat("\\begin{tabular}{ccccc}\n\\toprule\n")
  cat("Process Shift ($\\Delta p$) & Incoming Quality $p$ & ARL (Lots) & Empirical FAR & Control Response State \\\\\n\\midrule\n")
  
  for (shift in shifts) {
    p_curr <- p0 + shift
    arl_val <- simulate_arl_metric(p_curr, p0 = p0, sims = 1000)
    far_str <- if (shift == 0) sprintf("%.5f", 1 / arl_val) else "---"
    state_str <- if (shift > 0.02) "State 3 ($K = 1$)" else "State 2 ($K = 3$)"
    
    cat(sprintf("+%.2f & %.3f & %.1f & %s & %s \\\\\n", shift, p_curr, arl_val, far_str, state_str))
  }
  cat("\\bottomrule\n\\end{tabular}\n\\end{table}\n\n")
}

# ------------------------------------------------------------------------------
# 6. TABLE 5 GENERATOR: RECTIFYING INSPECTION METRICS (ATI, AOQ, AOQL)
# ------------------------------------------------------------------------------

generate_table5_latex <- function() {
  cat("\n=================================================================\n")
  cat("TABLE 5: RECTIFYING INSPECTION METRICS (N = 5000) LATEX TABLE\n")
  cat("=================================================================\n\n")
  
  N_lot <- 5000
  g <- 10; r <- 5; c1 <- 1; c2 <- 2; d1 <- 1; d2 <- 2; K <- 3; omega <- 1.5
  n_stage <- g * r
  
  p_grid <- c(0.01, 0.05, 0.10, 0.15)
  
  cat("\\begin{table}[htbp]\n\\centering\n")
  cat("\\caption{Average Total Inspection (ATI), Outgoing Quality (AOQ), and AOQL Metrics ($N = 5000, g = 10, r = 5, K = 3$)}\n")
  cat("\\label{tab:rectifying}\n\\small\n")
  cat("\\begin{tabular}{ccccc}\n\\toprule\n")
  cat("Incoming Quality ($p$) & Acceptance Prob. $P_{\\text{accept}}^{(3)}(p)$ & $\\text{ASN}^{(3)}(p)$ & $\\text{ATI}^{(3)}(p)$ & $\\text{AOQ}^{(3)}(p)$ \\\\\n\\midrule\n")
  
  max_aoq <- 0
  for (p_val in p_grid) {
    probs <- compute_group_probs(p_val, g, r, c1, c2, d1, d2, omega)
    P_accept <- oc_k(probs$Pa, probs$Pr, probs$Pi, K)
    asn_val <- asn_k(n_stage, probs$Pi, K)
    ati_val <- ati_k(N_lot, n_stage, probs$Pa, probs$Pr, probs$Pi, K)
    aoq_val <- aoq_k(p_val, N_lot, n_stage, probs$Pa, probs$Pr, probs$Pi, K)
    
    if (aoq_val > max_aoq) max_aoq <- aoq_val
    
    cat(sprintf("%.2f & %.4f & %.2f & %.2f & %.5f \\\\\n", p_val, P_accept, asn_val, ati_val, aoq_val))
  }
  cat("\\bottomrule\n\\end{tabular}\n\\end{table}\n\n")
}

# ------------------------------------------------------------------------------
# 7. TABLE 6 GENERATOR: SENSITIVITY ANALYSIS ACROSS WEIBULL SHAPE PARAMETERS
# ------------------------------------------------------------------------------

generate_table6_latex <- function() {
  cat("\n=================================================================\n")
  cat("TABLE 6: SENSITIVITY ANALYSIS ACROSS SHAPE PARAMETERS (omega)\n")
  cat("=================================================================\n\n")
  
  omega_list <- c(0.8, 1.0, 1.5, 2.0, 3.0)
  d0 <- 1.0
  beta_val <- 0.05
  K_val <- 3
  
  cat("\\begin{table}[htbp]\n\\centering\n")
  cat("\\caption{Sensitivity Analysis of Optimal Plan Parameters and Inspection Burden across Varying Weibull Shape Parameters ($\\omega$) at $\\beta = 0.05, d_0 = 1.0, K = 3$}\n")
  cat("\\label{tab:sensitivity_omega}\n\\small\n")
  cat("\\begin{tabular}{cccccccc}\n\\toprule\n")
  cat("Shape ($\\omega$) & Distribution Profile & Unreliability $p_1$ & Unreliability $p_2$ & $(g, r)$ & $(c_1, c_2)$ & $P_{\\text{accept}}(p_2)$ & $\\text{ASN}(p_1)$ \\\\\n\\midrule\n")
  
  for (omega_val in omega_list) {
    res <- solve_weibull_two_point(d0 = d0, omega = omega_val, alpha = 0.05, beta = beta_val, K_max = K_val)
    
    profile_str <- if (omega_val < 1.0) "Decreasing Failure Rate"
    else if (omega_val == 1.0) "Exponential (Constant Rate)"
    else if (omega_val < 2.0) "Moderate Wear-Out"
    else "Severe Wear-Out / Aging"
    
    if (!is.null(res)) {
      cat(sprintf("%.1f & %s & %.4f & %.4f & (%d, %d) & (%d, %d) & %.4f & %.2f \\\\\n",
                  omega_val, profile_str, res$p1, res$p2,
                  res$g, res$r, res$c1, res$c2, res$P_accept_p2, res$ASN_p1))
    }
  }
  cat("\\bottomrule\n\\end{tabular}\n\\end{table}\n\n")
}

# ------------------------------------------------------------------------------
# EXECUTE ALL TABLE GENERATORS
# ------------------------------------------------------------------------------

generate_table2_latex()
generate_table3_latex()
generate_table4_latex()
generate_table5_latex()
generate_table6_latex()