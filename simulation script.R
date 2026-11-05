
setwd("C:\\Users\\Dr. O. J. Obulezi\\Documents\\R projects\\hard-truncated EMWA")
# ==============================================================================
# Step 1: Theoretical Grid Optimization for Hard-Truncated Group Sampling
# ==============================================================================

# --- 1. Input Parameters ---
alpha <- 0.05    # Maximum Producer Risk (e.g., 5%)
beta  <- 0.10    # Maximum Consumer Risk (e.g., 10%)
r1    <- 1.0     # Unacceptable Quality Ratio (M1/M0)
r2    <- 2.0     # Acceptable Quality Ratio (M2/M0)
d0    <- 1.0     # Test termination ratio (t0/M0)
w     <- 1.5     # Weibull shape parameter
K     <- 3       # Truncation depth

# --- 2. Calculate Item Unreliability (p1 and p2) ---
# Unreliability (failure probability) under Weibull distribution (Eq. 19 & 20)
p1 <- 1 - exp(-log(2) * (d0)^w)       # Defect rate at unacceptable quality
p2 <- 1 - exp(-log(2) * (d0 / r2)^w)  # Defect rate at acceptable quality

cat(sprintf("Target Unreliability p1 (Unacceptable): %.4f\n", p1))
cat(sprintf("Target Unreliability p2 (Acceptable):   %.4f\n\n", p2))

# --- 3. Single-Cycle Binomial/Multinomial Helper Functions ---
# Item failure probability within a tester of size r
calc_item_p <- function(p_lot, r) {
  # Probability that a single item fails before time t0
  return(p_lot)
}

# Probability of a group tester failing location limit (c)
calc_pg <- function(p_lot, r, c_limit) {
  # P(Y > c_limit) where Y ~ Binomial(r, p_lot)
  return(1 - pbinom(c_limit, size = r, prob = p_lot))
}

# Compute single-cycle Pa, Pr, and Pi (Eq. 8, 9, 10)
calc_single_cycle <- function(g, r, c1, c2, d1, d2, p_lot) {
  # Note: Under standard single-variable failure mapping:
  # Acceptance occurs when failure count <= c1
  # Rejection occurs when failure count > c2
  # Indecisive region occurs when c1 < failure count <= c2
  
  Pa <- pbinom(c1, size = g, prob = calc_pg(p_lot, r, c1))
  Pr <- 1 - pbinom(c2 - 1, size = g, prob = calc_pg(p_lot, r, c1))
  
  # Ensure valid probability bounds
  Pa <- max(0, min(1, Pa))
  Pr <- max(0, min(1, Pr))
  
  # If Pa + Pr > 1 due to boundary definitions, constrain them
  if (Pa + Pr > 1) {
    Pr <- 1 - Pa
  }
  
  Pi <- 1 - Pa - Pr
  return(list(Pa = Pa, Pr = Pr, Pi = Pi))
}

# --- 4. Theorems 2.1 & 2.3: Closed-Form Truncated OC and ASN ---
calc_metrics <- function(Pa, Pr, K, g, r) {
  Pi <- 1 - Pa - Pr
  
  if (abs(Pi - 1) < 1e-8) {
    P_accept <- Pa * K
    ASN      <- (g * r) * K
  } else if (abs(Pi) < 1e-8) {
    P_accept <- Pa
    ASN      <- (g * r)
  } else {
    # Theorem 2.1: Exact Finite K-Truncated OC
    P_accept <- Pa * (1 - Pi^K) / (1 - Pi)
    # Theorem 2.3: Exact Finite K-Truncated ASN
    ASN      <- (g * r) * (1 - Pi^K) / (1 - Pi)
  }
  
  return(c(P_accept = P_accept, ASN = ASN))
}

# --- 5. Grid Search Optimization Routine ---
best_asn <- Inf
best_plan <- NULL

cat("Running Grid Search Optimization...\n")

# Bounded parameter search space
for (g in 2:15) {
  for (r in 2:10) {
    for (c1 in 0:(g - 1)) {
      for (c2 in (c1 + 1):g) {
        
        # Simplify dispersion boundaries for grid search mapping
        d1 <- c1
        d2 <- c2
        
        # --- Evaluate Consumer Risk at p1 ---
        cycle_p1 <- calc_single_cycle(g, r, c1, c2, d1, d2, p1)
        res_p1   <- calc_metrics(cycle_p1$Pa, cycle_p1$Pr, K, g, r)
        
        # Check Consumer Risk Constraint: P_accept(p1) <= beta
        if (res_p1["P_accept"] > beta) next
        
        # --- Evaluate Producer Risk at p2 ---
        cycle_p2 <- calc_single_cycle(g, r, c1, c2, d1, d2, p2)
        res_p2   <- calc_metrics(cycle_p2$Pa, cycle_p2$Pr, K, g, r)
        
        # Check Producer Risk Constraint: P_accept(p2) >= 1 - alpha
        if (res_p2["P_accept"] < (1 - alpha)) next
        
        # --- Select Plan with Minimal ASN(p1) ---
        current_asn <- res_p1["ASN"]
        
        if (current_asn < best_asn) {
          best_asn  <- current_asn
          best_plan <- data.frame(
            g         = g,
            r         = r,
            c1        = c1,
            c2        = c2,
            n_initial = g * r,
            ASN_p1    = round(res_p1["ASN"], 2),
            P_acc_p2  = round(res_p2["P_accept"], 4),
            P_acc_p1  = round(res_p1["P_accept"], 4)
          )
        } else if (abs(current_asn - best_asn) < 1e-5) {
          # Tie-breaking rule: pick smaller initial group size g
          if (g < best_plan$g) {
            best_plan <- data.frame(
              g         = g,
              r         = r,
              c1        = c1,
              c2        = c2,
              n_initial = g * r,
              ASN_p1    = round(res_p1["ASN"], 2),
              P_acc_p2  = round(res_p2["P_accept"], 4),
              P_acc_p1  = round(res_p1["P_accept"], 4)
            )
          }
        }
        
      }
    }
  }
}

# --- 6. Display Optimal Plan Results ---
cat("\n======================================================\n")
cat("            OPTIMAL SAMPLING PLAN RESULTS             \n")
cat("======================================================\n")
if (!is.null(best_plan)) {
  print(best_plan, row.names = FALSE)
} else {
  cat("No valid plan found within the given parameter ranges.\n")
}

# ==============================================================================
# Step 2: Monte Carlo Simulation & EWMA Process Monitoring
# ==============================================================================

set.seed(42)

# --- 1. Simulation Parameters ---
N_lots  <- 50       # Total submitted lots in sequence
shift_lot <- 25     # Lot number where process drift occurs
p_base  <- 0.05     # Baseline in-control defect rate (p0)
p_shift <- 0.10     # Out-of-control defect rate after drift (Delta p = +0.05)

# Optimal plan parameters from Step 1
g        <- 5
r        <- 4
n_sample <- g * r   # Nominal stage sample size (n = 20)
K_base   <- 3       # Normal truncation depth
K_max    <- 4       # Expanded truncation depth

# EWMA Control Parameters (Table 4 settings)
lambda <- 0.20
L      <- 3.0
Z0     <- p_base

# --- 2. Calculate Asymptotic EWMA Control Limits (Eq. 12) ---
sigma_Z <- sqrt((p_base * (1 - p_base) / n_sample) * (lambda / (2 - lambda)))
UCL     <- p_base + L * sigma_Z
LCL     <- max(0, p_base - L * sigma_Z)

cat(sprintf("EWMA Control Limits -> LCL: %.4f | UCL: %.4f\n\n", LCL, UCL))

# --- 3. Run Sequential Lot Simulation ---
Z          <- numeric(N_lots)
p_true     <- numeric(N_lots)
K_applied  <- numeric(N_lots)
p_hat      <- numeric(N_lots)
state      <- character(N_lots)

Z_prev <- Z0

for (i in 1:N_lots) {
  # Define true defect rate for current lot
  p_current <- if (i < shift_lot) p_base else p_shift
  p_true[i] <- p_current
  
  # Determine Adaptive Truncation Depth (K_i) based on previous Z_{i-1}
  if (i == 1) {
    K_curr <- K_base
  } else if (Z_prev > UCL) {
    K_curr <- 1         # State 3: Hard stop due to out-of-control signal
  } else if (Z_prev < LCL) {
    K_curr <- K_max       # State 1: In-control high quality
  } else {
    K_curr <- K_base      # State 2: Normal operation
  }
  K_applied[i] <- K_curr
  
  # Simulate group testing inspection to obtain sample failure fraction (p_hat)
  # Total items inspected in single stage = g * r
  failures  <- rbinom(1, size = n_sample, prob = p_current)
  p_hat[i]  <- failures / n_sample
  
  # Update EWMA Tracking Statistic (Eq. 11)
  Z[i]   <- lambda * p_hat[i] + (1 - lambda) * Z_prev
  
  # Record Process State
  if (Z[i] > UCL) {
    state[i] <- "State 3 (Shift Detected -> K=1)"
    # Optional reset after action: Z_prev <- p_base
    Z_prev <- Z[i]
  } else if (Z[i] < LCL) {
    state[i] <- "State 1 (High Quality -> K_max)"
    Z_prev <- Z[i]
  } else {
    state[i] <- "State 2 (Normal -> K_base)"
    Z_prev <- Z[i]
  }
}

# --- 4. Summarize Simulation Results ---
sim_results <- data.frame(
  Lot       = 1:N_lots,
  True_p    = p_true,
  p_hat     = round(p_hat, 4),
  EWMA_Z    = round(Z, 4),
  Applied_K = K_applied,
  State     = state
)

cat("First 10 Lots (In-Control Phase):\n")
print(head(sim_results, 10))

cat("\nLots Around Process Shift (Lots 22-32):\n")
print(sim_results[22:32, ])
windows()
# --- 5. Plot EWMA Tracking Chart ---
plot(1:N_lots, Z, type = "o", col = "blue", pch = 16,
     ylim = c(0, max(Z, UCL) * 1.2),
     xlab = "Submitted Lot Sequence", ylab = "EWMA Statistic (Z_i)",
     main = "")

abline(h = UCL, col = "red", lty = 2, lwd = 2)
abline(h = LCL, col = "darkgreen", lty = 2, lwd = 2)
abline(h = p_base, col = "gray", lty = 3)
abline(v = shift_lot, col = "purple", lty = 4)
text(shift_lot + 1, max(Z)*1.1, "Shift Introduced", pos = 4, col = "purple")
legend("topleft", legend = c("EWMA Z_i", "UCL", "LCL", "Baseline p0"),
       col = c("blue", "red", "darkgreen", "gray"), lty = c(1, 2, 2, 3), pch = c(16, NA, NA, NA))

# ==============================================================================
# Step 3: Empirical ARL Performance Benchmarking (10,000 Monte Carlo Runs)
# ==============================================================================

set.seed(123)

# System Setup
g <- 5; r <- 4; n_sample <- g * r
lambda <- 0.20; L <- 3.0; p0 <- 0.05
sigma_Z <- sqrt((p0 * (1 - p0) / n_sample) * (lambda / (2 - lambda)))
UCL <- p0 + L * sigma_Z
LCL <- max(0, p0 - L * sigma_Z)

# Vector of process shifts (Delta p) matching Table 4 in the paper
shifts <- c(0.00, 0.01, 0.03, 0.05)
N_sims <- 10000
max_lots <- 2000

arl_results <- data.frame(Shift = character(), Incoming_p = numeric(), 
                          ARL = numeric(), State_Trigger = character())

cat("Running 10,000 Monte Carlo Simulations per shift level...\n\n")

for (delta in shifts) {
  p_test <- p0 + delta
  rl_vec <- numeric(N_sims)
  
  for (s in 1:N_sims) {
    Z_prev <- p0
    run_length <- 0
    
    repeat {
      run_length <- run_length + 1
      failures <- rbinom(1, size = n_sample, prob = p_test)
      p_hat <- failures / n_sample
      Z_curr <- lambda * p_hat + (1 - lambda) * Z_prev
      
      # Stop run when process shift is detected
      if (Z_curr > UCL || run_length >= max_lots) {
        rl_vec[s] <- run_length
        break
      }
      Z_prev <- Z_curr
    }
  }
  
  mean_arl <- mean(rl_vec)
  state_action <- if (delta == 0) "State 2 (In-Control)" else "State 3 (Hard Stop K=1)"
  
  arl_results <- rbind(arl_results, data.frame(
    Shift         = sprintf("+%.2f", delta),
    Incoming_p    = p_test,
    ARL           = round(mean_arl, 1),
    State_Trigger = state_action
  ))
}

cat("======================================================\n")
cat("      EWMA CONTROLLER ARL PERFORMANCE BENCHMARK        \n")
cat("======================================================\n")
print(arl_results, row.names = FALSE)