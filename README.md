# Hard-Truncated EWMA Bivariate Group Acceptance Sampling Scheme

[![R](https://img.shields.io/badge/Language-R-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

This repository contains the official R implementation for the **Hard-Truncated Adaptive Group Acceptance Sampling Plan** under Weibull lifetime processes. The proposed methodology introduces a finite truncation depth ($K$) to eliminate unbounded re-testing risks inherent in traditional Repetitive Acceptance Sampling Plans (RASP), evaluates bivariate quality standards (median lifetime $M$ and tail-risk dispersion $\omega$), and integrates an Exponentially Weighted Moving Average (EWMA) process-memory controller for real-time drift detection.

---

## 📌 Features & Key Capabilities

- **Exact Closed-Form Metrics**: Derivations for Operating Characteristic ($OC^{(K)}$), Average Sample Number ($ASN^{(K)}$), Inspection Burden Variance ($\text{Var}(ASN)$), Average Total Inspection ($ATI^{(K)}$), Average Outgoing Quality ($AOQ^{(K)}$), and $AOQL^{(K)}$.
- **Bivariate Quality Constraints**: Evaluates location ($Y_M$) and failure dispersion ($Y_\omega$) simultaneously using trinomial group multinomial probabilities.
- **Dual Optimization Solvers**:
  - **Classical $(\alpha, \beta)$ Solver**: Minimizes $ASN(p_1)$ subject to producer risk $\alpha$ and consumer risk $\beta$.
  - **Adaptive Weighted-Risk ($\gamma$-WR) Solver**: Minimizes $ASN(p_1)$ under combined risk limits $\gamma$ and budget variance constraint $V_{\text{max}}$.
- **EWMA Process-Memory Controller**: Simulates lot-by-lot process monitoring with dynamic truncation depth adaptation ($K_i \in \{1, K_{\text{base}}, K_{\text{max}}\}$).
- **Automated Exporters**: Formats numerical optimization results directly into APA-compliant LaTeX tables and generates diagnostic visualizations.

---

## 📁 Repository Structure

```text
.
├── R/
│   ├── comprehensive_script.R     # Self-contained solver for Tables 2–6 (WR, Benchmark, ARL, ATI/AOQ, Sensitivity)
│   └── script_for_EWMA_scheme.R   # Core engine with diagnostic plot exporters and Table 1 solver
├── output/                        # Saved diagnostic performance charts (PNG)
├── README.md                      # Project documentation
└── LICENSE                        # MIT License
