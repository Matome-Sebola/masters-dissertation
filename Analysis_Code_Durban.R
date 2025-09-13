# =======================================================================
# MAIN SCRIPT: Time Varying Threshold Generalised Pareto Regression Tree
# =======================================================================

# Load packages
library(rpart)
library(evir)
library(testit)
library(ggplot2)
library(tseries)
library(fBasics)
library(treeClust)
library(readxl)
library(mgcv)
library(testit)
library(rpart.plot)
library(dplyr)
library(corrplot)
library(glmnet)
library(testit)
library(evmix)

# Source GP regression tree method functions
source("./static_gprt_func.R")
source("./GPRT_func.R")
#--------------------------
# STEP 1: Load the dataset
# --------------------------
rainfall <- read_excel("C:/Users/matom/OneDrive/Data Analysis/durban_station_rainfall.xlsx")
options(scipen = 999)
R <-ts(rainfall$PRECTOTCORR)

#-----------------------------------
# STEP 2: Exploratory Data Analysis
# ----------------------------------
win.graph()
par(mfrow =c(2 ,2))

ggplot(rainfall,aes(x=YEAR, y=PRECTOTCORR)) + geom_line() + xlab("")
plot (R , xlab =" Observation number " , ylab =" Rainfall (mm) ", main = "(a) Plot of Rainfall",col = "red")
qqnorm (R, col = "red",main ="(b) Normal QQ plot ")
qqline (R)
boxplot (rainfall$PRECTOTCORR, main ="(c) Box plot ",varwidth =TRUE ,xlab = " Rainfall (mm) ", col = "red" ,horizontal = TRUE )
plot ( density (R),xlab =" Rainfall (mm) " ,main ="(d) Density plot ",col = "red")

#-------------------------------
# STEP 3: Descriptive Statistics
# ------------------------------
summary(R)
basicStats(R)
#------------------------------

#--------------------------
# STEP 4: Data Preparations
# -------------------------
# Separate features (X) and target (y)
rainfall <- rainfall[, !(names(rainfall) %in% c("YEAR","DOY"))]

# Set split index (80% for training)
split_index <- floor(nrow(rainfall) * 0.8)

# Train and test datasets
train_data <- rainfall[1:split_index, ]
test_data  <- rainfall[(split_index + 1):nrow(rainfall), ]

# ---------------------------
# STEP 5: Threshold Selection
# ---------------------------
# Selecting the best variables to estimate a time-varying threshold is a complex task that requires careful consideration of the underlying data characteristics, the purpose of the threshold, and the chosen estimation methodology.
# We use LASSO regression to select to select variables for estimating a time-varying threshold

# --------------------------
# STEP 5A: LASSO Regression
# --------------------------
y <- train_data$PRECTOTCORR
x <- model.matrix(PRECTOTCORR ~ . , data = train_data)[, -1]


set.seed(123)  # for reproducibility

cv_fit <- cv.glmnet(x, y, alpha = 1)  # alpha=1 for LASSO

# Plot cross-validation curve
plot(cv_fit)

best_lambda_1se <- cv_fit$lambda.1se
final_model_1se <- glmnet(x, y, alpha = 1, lambda = best_lambda_1se)
coef(final_model_1se)

# -----------------------
# STEP 6B: Fit Threshold
# -----------------------

fit <- gam(PRECTOTCORR ~ s(T2M_MAX) + s(WS2M) + s(WD2M) + s(RH2M) + s(WS2M_MAX) + s(WS2M_MIN), 
           family = Gamma(link = "log"), data = train_data)

threshold <- predict(fit, type = "response")  # Returns values on original scale
# -----------------------
# STEP 6C: Plot Threshold
# -----------------------

start_year <- 1981
end_year <-2025
years <- seq(start_year, end_year, by = 5)  

break_positions <- (years - start_year) * 295 +1

# Plot the threshold with legend
ggplot(train_data, aes(x = 1:nrow(train_data))) +
  geom_point(aes(y = PRECTOTCORR, colour = "Observed Rainfall")) +
  geom_line(aes(y = threshold, colour = "Time-Varying Threshold")) +
  scale_color_manual(
    values = c("Observed Rainfall" = "black", "Time-Varying Threshold" = "red"),
    name = ""  # legend title
  ) +
  scale_x_continuous(breaks = break_positions, labels = years) +
  labs(x = "Years", y = "Rainfall (mm)")


# Observed response
observed <- train_data$PRECTOTCORR

# Count exceedances
exceedances <- sum(observed > threshold)


cat("Number of exceedances above threshold:", exceedances, "\n")

parms = list(threshold = threshold)

# Ensure the threshold vector has the same length as the training data
stopifnot(length(threshold) == nrow(train_data))

# Add the threshold as a column to the training data frame
train_data$threshold_col <- threshold

# Create the matrix response variable
response_matrix <- cbind(train_data$PRECTOTCORR, train_data$threshold_col)


# -------------------------------------------
# STEP 6C: Plot Threshold for static GPRT
# ----------------------------------------
mrlplot(R)
tcplot(R)


# ----------------------------------------
# STEP 7: Fit Time-Varying GPRT (Unpruned)
# ----------------------------------------
# Grow the tree with the updated data
gprt_model <- rpart(
  formula = response_matrix ~ T2M + T2M_MAX + T2M_MIN + PS + 
    WS10M + WS2M + WD10M + WD2M + RH2M + WS10M_MAX + WS10M_MIN + 
    WS2M_MAX + WS2M_MIN,
  data = train_data,
  method = GPRT_method_xi_in_R_plus,
  parms = list(),
  control = rpart.control(minsplit = 100, minbucket = 70, cp = 0.0015)
)
printcp(gprt_model)
print(gprt_model)


# -------------------------------------------
# STEP 7A: Plot Time-Varying GPRT (Unpruned)
# -------------------------------------------

# Base tree plot
par(bg = "white", mar = c(1, 1, 1, 1))
plot(gprt_model, uniform = TRUE, margin = 0.1)

# Extract node positions and frame info
coords <- rpart:::rpartco(gprt_model)
frame <- gprt_model$frame
yvals <- frame$yval2

# Draw shapes first
for (i in seq_len(nrow(frame))) {
  x <- coords$x[i]
  y <- coords$y[i]
  
  if (frame$var[i] != "<leaf>") {
    # Rectangle for split nodes
    rect(x - 0.7, y - 0.03, x + 0.7, y + 0.03, col = "lightblue", border = "blue")
  } else {
    # Circle for terminal nodes
    symbols(x, y, circles = 0.4, inches = FALSE, add = TRUE, bg = "lightgreen")
  }
}

# Add default node labels (centered inside shapes)
text(gprt_model, use.n = TRUE, all = TRUE, cex = 0.8, pos=1)

# Annotate GPD parameters and sample sizes
for (i in seq_len(nrow(frame))) {
  x <- coords$x[i]
  y <- coords$y[i]
  
  shape <- round(yvals[i, 1], 3)
  scale <- round(yvals[i, 2], 3)
  gpd_label <- paste0("ξ = ", shape, "\nσ = ", scale)
  
  if (frame$var[i] != "<leaf>") {
    # GPD parameters below rectangle
    text(x, y - 0.08, gpd_label, cex = 0.75, col = "blue", adj = 0.5)
  } else {
    # GPD parameters inside circle
    text(x, y, gpd_label, cex = 0.75, col = "blue", adj = 0.5)
    
    # Sample size below circle
    n_obs <- frame$n[i]
    text(x, y - 0.10, paste0("n = ", n_obs), cex = 0.75, adj = 0.5)
  }
}

# ----------------------------------------------------
# STEP 7B: Penalized Pruning of time-varying GPTR model
# ----------------------------------------------------
A_candidates <- seq(20, 200, by = 20)
pruning_result <- cross_validate_penalty(gprt_model, A_candidates)
pruned_tree <- pruning_result$tree
optimal_A <- pruning_result$A
cat("Optimal penalty A:", optimal_A, "\n")

# -------------------------------------------
# STEP 7C: Plot pruned time-varying GPRT
# -------------------------------------------
# Base tree plot
par(bg = "white", mar = c(2, 2, 2, 2))
plot(pruned_tree, uniform = TRUE, margin = 0.1)


# Extract node positions and frame info
coords <- rpart:::rpartco(pruned_tree)
frame <- pruned_tree$frame
yvals <- frame$yval2

# Draw shapes first
for (i in seq_len(nrow(frame))) {
  x <- coords$x[i]
  y <- coords$y[i]
  
  if (frame$var[i] != "<leaf>") {
    # Rectangle for split nodes
    rect(x - 0.5, y - 0.03, x + 0.5, y + 0.03, col = "lightblue", border = "blue")
  } else {
    # Circle for terminal nodes
    symbols(x, y, circles = 0.3, inches = FALSE, add = TRUE, bg = "lightgreen")
  }
}

# Add default node labels (centered inside shapes)
text(pruned_tree, use.n = TRUE, all = TRUE, cex = 0.8, pos=1)

# Annotate GPD parameters and sample sizes
for (i in seq_len(nrow(frame))) {
  x <- coords$x[i]
  y <- coords$y[i]
  
  shape <- round(yvals[i, 1], 3)
  scale <- round(yvals[i, 2], 3)
  gpd_label <- paste0("ξ = ", shape, "\nσ = ", scale)
  
  if (frame$var[i] != "<leaf>") {
    # GPD parameters below rectangle
    text(x, y - 0.08, gpd_label, cex = 0.75, col = "blue", adj = 0.5)
  } else {
    # GPD parameters inside circle
    text(x, y, gpd_label, cex = 0.75, col = "blue", adj = 0.5)
    
    # Sample size below circle
    n_obs <- frame$n[i]
    text(x, y - 0.10, paste0("n = ", n_obs), cex = 0.75, adj = 0.5)
  }
}

# ----------------------------
# STEP 8: Benchmarking models
# ----------------------------

# ----------------------------
# STEP 8A: Static GPRT
# ----------------------------
# Static Generalised Pareto regression trees
static_gprt <- rpart(
  formula = PRECTOTCORR ~ T2M + T2M_MAX + T2M_MIN + PS + 
    WS10M + WS2M + WD10M + WD2M + RH2M + WS10M_MAX + WS10M_MIN + 
    WS2M_MAX + WS2M_MIN,
  data = train_data,
  method = GPRT_method,
  control = rpart.control(
    minsplit = 100, # Minimum number of observations in a node to be considered for a split
    minbucket = 70,# Minimum number of observations in any terminal node (leaf)
    xval = 10,
    cp = 0.0015# Complexity parameter; setting to 0 ensures a full tree is grown
  )
)

cv_results <- cross_val_gprt(
  data = train_data,
  formula = PRECTOTCORR ~ T2M + T2M_MAX + T2M_MIN + PS + 
    WS10M + WS2M + WD10M + WD2M + RH2M + WS10M_MAX + WS10M_MIN + 
    WS2M_MAX + WS2M_MIN,
  arbre = static_gprt,
  n_fold = 10, # Number of folds for cross-validation
  seed = 42 # For reproducibility
)

optimal_cp <- cv_results[[1]]
cat("Optimal Cost parameter CP:", optimal_cp, "\n")

pruned_gprt <- prune(static_gprt, cp = 0.0057574 )

printcp(static_gprt)
# -------------------------------
# STEP 8B: Fit time varying GPD
# -------------------------------

# Observed rainfall values
observed <- train_data$PRECTOTCORR

# Your vector of thresholds (same length as observed)
threshold_vec <- threshold  

# Step 1: Extract exceedances (elementwise)
exceedances <- observed[observed > threshold_vec]

# Step 2: Fit GPD with threshold = 0
gpd_model <- gpd(exceedances, threshold = 0)

# Step 3: Inspect results
print(gpd_model)


# ------------------------------------------
# STEP 9: Model Evaluation Accuracy Metrics
# ------------------------------------------
#Static threshold gprt
gprt_result <- compute_bic_gprt(static_gprt, train_data)
print(gprt_result)

#Time varying threshold gprt
time_gprt_result <- compute_bic_gprt(pruned_tree, train_data)
print(time_gprt_result)

#Time varying gpd
gpd_result <- compute_bic_gpd(gpd_model, length(exceedances))
print(gpd_result)



# -----------------------------------------------------
# STEP 9: Model Evaluation Accuracy Metrics on Test Set
# -----------------------------------------------------
## -------------------------------
# BIC and Log-Likelihood on Test Set
# -------------------------------

# Helper: Compute log-likelihood and BIC for GPRT models
compute_bic_gprt <- function(tree, data) {
  loglik <- -sum(tree$frame[tree$frame$var == "<leaf>", "dev"], na.rm = TRUE)
  n_params <- sum(tree$frame$var == "<leaf>") * 2  # shape and scale per leaf
  n_obs <- nrow(data)
  bic <- -2 * loglik + log(n_obs) * n_params
  return(list(logLik = loglik, BIC = bic))
}

# Helper: Compute log-likelihood and BIC for GPD model
compute_bic_gpd <- function(gpd_model, n_obs) {
  loglik <- -gpd_model$nllh.final
  n_params <- 2  # shape and scale
  bic <- -2 * loglik + log(n_obs) * n_params
  return(list(logLik = loglik, BIC = bic))
}

# -------------------------------
# 1. Static GPRT Evaluation
# -------------------------------
static_gprt_test_result <- compute_bic_gprt(static_gprt, test_data)
cat("Static GPRT - LogLik:", static_gprt_test_result$logLik, "\n")
cat("Static GPRT - BIC:", static_gprt_test_result$BIC, "\n")

# -------------------------------
# 2. Pruned Time-Varying GPRT Evaluation
# -------------------------------
pruned_tree_test_result <- compute_bic_gprt(gprt_model, test_data)
cat("Pruned Time-Varying GPRT - LogLik:", pruned_tree_test_result$logLik, "\n")
cat("Pruned Time-Varying GPRT - BIC:", pruned_tree_test_result$BIC, "\n")

# -------------------------------
# 3. Time-Varying GPD Evaluation
# -------------------------------
# Use the GAM model fitted on training data to predict threshold for test data
threshold_vec_test <- predict(fit, newdata = test_data, type = "response")

# Ensure alignment
stopifnot(length(threshold_vec_test) == nrow(test_data))

observed_test <- test_data$PRECTOTCORR
exceedances_test <- observed_test[observed_test > threshold_vec_test]

loglik_test <- -sum(log(evir::dgpd(exceedances_test,
                                   xi = gpd_model$par.ests[1],
                                   beta = gpd_model$par.ests[2],
                                   mu = 0)), na.rm = TRUE)

bic_test <- -2 * loglik_test + log(length(exceedances_test)) * 2

cat("Time-Varying GPD (Original Fit) - LogLik:", loglik_test, "\n")
cat("Time-Varying GPD (Original Fit) - BIC:", bic_test, "\n")


# Compute log-likelihood using original parameters
loglik_test <- -sum(log(evir::dgpd(exceedances_test,
                                   xi = gpd_model$par.ests[1],
                                   beta = gpd_model$par.ests[2],
                                   mu = 0)), na.rm = TRUE)

bic_test <- -2 * loglik_test + log(length(exceedances_test)) * 2

cat("Time-Varying GPD (Original Fit) - LogLik:", loglik_test, "\n")
cat("Time-Varying GPD (Original Fit) - BIC:", bic_test, "\n")
