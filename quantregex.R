###############################################################################
# Expectile / Quantile-Type Nonlinear Regression
# LS initialization + Sandwich Variance Confidence Intervals
###############################################################################

###############################################################################
# 1. Read Data
###############################################################################

library(readxl)

a <- read_excel(
  "quantregexdata.xlsx"
)

names(a) <- c("cd19pos_bcells_pct_tot_c",
              "AGE",
              "Charletson_Comorbidity_Index__CC")

###############################################################################
# 2. Linear Regression (LS Starting Values)
###############################################################################

fit_lm <- lm(cd19pos_bcells_pct_tot_c ~ AGE +
               Charletson_Comorbidity_Index__CC,
             data = a)

rmse <- sqrt(mean(residuals(fit_lm)^2))
cat("Root MSE =", rmse, "\n")

###############################################################################
# 3. Nonlinear Optimization (PROC NLP Equivalent)
###############################################################################

n   <- nrow(a)
tau <- rmse / sqrt(n)

# Quantile (expectile) of interest
quantile_of_interest <- 0.5

# Confidence level
alpha <- 0.05
zcrit <- qnorm(1 - alpha/2)

# Extract variables
y  <- a$cd19pos_bcells_pct_tot_c
x1 <- a$AGE
x2 <- a$Charletson_Comorbidity_Index__CC

objective <- function(par, y, x1, x2, tau, q) {
  
  beta0 <- par[1]
  beta1 <- par[2]
  beta2 <- par[3]
  
  theta <- beta0 + beta1*x1 + beta2*x2
  
  loss <- tau * (
            log1p(exp(-(y - theta)/tau)) +
            log1p(exp((y - theta)/tau))
          ) +
          (2*q - 1)*(y - theta)
  
  sum(loss)
}

###############################################################################
# 4. Starting Values from LS
###############################################################################

init <- coef(fit_lm)

###############################################################################
# 5. Run Optimization
###############################################################################

fit_nlp <- optim(
  par = init,
  fn  = objective,
  y   = y,
  x1  = x1,
  x2  = x2,
  tau = tau,
  q   = quantile_of_interest,
  method  = "BFGS",
  hessian = TRUE
)

###############################################################################
# 6. Proper Sandwich Variance Estimator
###############################################################################

estimates <- fit_nlp$par

beta0 <- estimates[1]
beta1 <- estimates[2]
beta2 <- estimates[3]

theta_hat <- beta0 + beta1*x1 + beta2*x2

# Score ψ_i
psi_i <- 2 * ( 1/(1 + exp((y - theta_hat)/tau)) - quantile_of_interest )

# Derivative ψ'_i
dpsi_i <- (2/tau) *
          exp((y - theta_hat)/tau) /
          (1 + exp((y - theta_hat)/tau))^2

# Design matrix
X <- cbind(1, x1, x2)

###############################################################################
# A matrix (sensitivity)
###############################################################################

A_hat <- matrix(0, 3, 3)
for(i in 1:n){
  A_hat <- A_hat + dpsi_i[i] * (X[i,] %*% t(X[i,]))
}
A_hat <- A_hat / n

###############################################################################
# B matrix (variability)
###############################################################################

B_hat <- matrix(0, 3, 3)
for(i in 1:n){
  B_hat <- B_hat + psi_i[i]^2 * (X[i,] %*% t(X[i,]))
}
B_hat <- B_hat / n

###############################################################################
# Sandwich Covariance
###############################################################################

cov_matrix <- solve(A_hat) %*% B_hat %*% solve(A_hat) / n

std_errors <- sqrt(diag(cov_matrix))

###############################################################################
# Wald Confidence Intervals
###############################################################################

lower_ci <- estimates - zcrit * std_errors
upper_ci <- estimates + zcrit * std_errors

results <- data.frame(
  Parameter = c("Intercept","AGE","Charletson_Comorbidity_Index__CC"),
  Estimate  = as.numeric(estimates),
  StdError  = std_errors,
  Lower_CI  = lower_ci,
  Upper_CI  = upper_ci,
  row.names = NULL
)

###############################################################################
# 7. Final Summary
###############################################################################

cat("\nQuantile of Interest =", quantile_of_interest, "\n")
cat("Confidence Level =", 100*(1-alpha), "%\n")
cat("Objective Function Value =", fit_nlp$value, "\n")
cat("Convergence Code =", fit_nlp$convergence, "\n\n")
print(results)

