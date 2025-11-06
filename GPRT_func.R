# Generalized Pareto Regression Trees with Time-Varying Thresholds

library(rpart)
library(testit)
library(evir)
library(treeClust)

# Helper function to handle negative log-likelihood values
positive_nll <- function(x){
  if(x > 0){
    return(x + 1)
  } else {
    return(exp(x))
  }
}

# -------------------------------------
# Functions for growing the tree
# -------------------------------------

# Initialisation function (updated for a matrix response)
itempgprt <- function(y, offset, parms, wt) {
  if (!is.matrix(y) || ncol(y) != 2)
    stop("Matrix response with 2 columns required: (response, threshold)")
  
  sfun <- function(yval, dev, wt, ylevel, digits ) {
    paste(" Shape=", format(signif(yval[1], digits)),
          ", Scale=" , format(signif(yval[2], digits)),
          ", Deviance=" , format(signif(dev/wt, digits)),
          sep = '')
  }
  environment(sfun) <- .GlobalEnv
  list(y = y, parms = parms, numresp = 2, numy = 2, summary = sfun)
}

# Evaluation function with time-varying threshold (for a matrix response)
etempgprt <- function(y, wt, parms) {
  rainfall <- y[, 1]
  threshold <- y[, 2]
  node_size <- length(rainfall)
  
  cat("\n--- Evaluating Node ---\n")
  cat("Node size:", node_size, "\n")
  cat("Rainfall summary:", summary(rainfall), "\n")
  cat("Threshold summary:", summary(threshold), "\n")
  
  exceedances <- rainfall - threshold
  positive_exceedances <- exceedances[exceedances > 0]
  cat("Positive exceedances:", length(positive_exceedances), "\n")
  
  leval <- FALSE
  fit_try <- try(evir::gpd(positive_exceedances, threshold = 0, method = "ml"), silent = TRUE)
  
  sh <- NA
  sc <- NA
  rss <- Inf
  
  if (!has_warning(fit_try) && !has_error(fit_try) && !is.null(fit_try$par.ests) &&
      !is.na(fit_try$par.ests[1]) && !is.na(fit_try$par.ests[2]) &&
      is.finite(fit_try$nllh.final)) {
    
    leval <- TRUE
    fit <- fit_try
    sh <- fit$par.ests[1]
    sc <- fit$par.ests[2]
    rss <- positive_nll(fit$nllh.final)
    
    cat("GPD fit successful\n")
    cat("Shape (ξ):", round(sh, 4), "\n")
    cat("Scale (σ):", round(sc, 4), "\n")
    cat("Deviance:", round(rss, 4), "\n")
    
  } else {
    cat("GPD fit failed or returned invalid values\n")
  }
  
  list(label = c(sh, sc), deviance = rss)
}

# Split function with time-varying threshold (revised for matrix response)
stempgprt <- function(y, wt, x, parms, continuous) {
  rainfall <- y[, 1]
  threshold <- y[, 2]
  
  # Define a minimum number of exceedances for a stable fit
  min_exceedances <- 70 # 
  
  stopifnot(all(!is.na(rainfall)), all(!is.na(threshold)))
  
  parent_exceedances <- rainfall - threshold
  parent_exceedances <- parent_exceedances[parent_exceedances > 0]
  
  # If the parent node itself is too small to fit, don't split
  if (length(parent_exceedances) < min_exceedances) {
    return(list(goodness = 0, direction = 0))
  }
  
  cat("Number of exceedances in parent node:", length(parent_exceedances), "\n")
  
  failed_deviance <- Inf
  
  lparent <- FALSE
  fitparent_try <- try(evir::gpd(parent_exceedances, threshold = 0, method = "ml"), silent = TRUE)
  
  parent_nll <- failed_deviance
  if (!has_warning(fitparent_try) && !has_error(fitparent_try) && !is.null(fitparent_try$par.ests)) {
    fitparent <- fitparent_try
    if (fitparent$par.ests[1] > 0 && is.finite(fitparent$nllh.final)) {
      lparent <- TRUE
      parent_nll <- positive_nll(fitparent$nllh.final)
    }
  }
  
  # If parent fit is invalid, no split is possible
  if (!lparent) {
    return(list(goodness = 0, direction = 0))
  }
  
  if (continuous) {
    n <- length(rainfall)
    goodness <- double(n - 1)
    direction <- goodness
    
    for (i in 1:(n - 1)) {
      child1_indices <- 1:i
      child2_indices <- (i + 1):n
      
      child1_exceedances <- rainfall[child1_indices] - threshold[child1_indices]
      child1_exceedances <- child1_exceedances[child1_exceedances > 0]
      
      child2_exceedances <- rainfall[child2_indices] - threshold[child2_indices]
      child2_exceedances <- child2_exceedances[child2_exceedances > 0]
      
      if (length(child1_exceedances) >= min_exceedances &&
          length(child2_exceedances) >= min_exceedances) {
        
        fitmoins_try <- try(evir::gpd(child1_exceedances, threshold = 0, method = "ml"), silent = TRUE)
        fitplus_try <- try(evir::gpd(child2_exceedances, threshold = 0, method = "ml"), silent = TRUE)
        
        if (!has_warning(fitmoins_try) && !has_error(fitmoins_try) && !is.null(fitmoins_try$par.ests) &&
            !has_warning(fitplus_try) && !has_error(fitplus_try) && !is.null(fitplus_try$par.ests)) {
          
          fitmoins <- fitmoins_try
          fitplus <- fitplus_try
          
          if (fitmoins$par.ests[1] > 0 && fitplus$par.ests[1] > 0 &&
              is.finite(fitmoins$nllh.final) && is.finite(fitplus$nllh.final)) {
            
            goodness[i] <- max(0, parent_nll - (positive_nll(fitmoins$nllh.final) + positive_nll(fitplus$nllh.final)))
            direction[i] <- sign(positive_nll(fitmoins$nllh.final) - positive_nll(fitplus$nllh.final))
          } else {
            goodness[i] <- 0 # If fit is bad, goodness is zero
          }
        } else {
          goodness[i] <- 0 # If fit failed, goodness is zero
        }
      } else {
        goodness[i] <- 0 # If not enough exceedances, goodness is zero
      }
      
      if (goodness[i] < 0) print("Warning: negative goodness")
      if (is.infinite(goodness[i])) print("Warning: infinite goodness")
    }
    
    list(goodness = goodness, direction = direction)
    
  } else {
    ux <- sort(unique(x))
    ymean <- tapply(rainfall, x, mean)
    ord <- order(ymean)
    n <- length(ord)
    
    goodness <- double(n - 1)
    direction <- goodness
    
    for (i in 1:(n - 1)) {
      ymoins_indices <- x %in% ux[ord[1:i]]
      yplus_indices <- x %in% ux[ord[(i + 1):n]]
      
      ymoins_exceedances <- rainfall[ymoins_indices] - threshold[ymoins_indices]
      ymoins_exceedances <- ymoins_exceedances[ymoins_exceedances > 0]
      
      yplus_exceedances <- rainfall[yplus_indices] - threshold[yplus_indices]
      yplus_exceedances <- yplus_exceedances[yplus_exceedances > 0]
      
      if (length(ymoins_exceedances) >= min_exceedances &&
          length(yplus_exceedances) >= min_exceedances) {
        
        fitmoins_try <- try(evir::gpd(ymoins_exceedances, threshold = 0, method = "ml"), silent = TRUE)
        fitplus_try <- try(evir::gpd(yplus_exceedances, threshold = 0, method = "ml"), silent = TRUE)
        
        if (!has_warning(fitmoins_try) && !has_error(fitplus_try) && !is.null(fitmoins_try$par.ests) &&
            !has_warning(fitplus_try) && !has_error(fitplus_try) && !is.null(fitplus_try$par.ests)) {
          
          fitmoins <- fitmoins_try
          fitplus <- fitplus_try
          
          if (fitmoins$par.ests[1] > 0 && fitplus$par.ests[1] > 0 &&
              is.finite(fitmoins$nllh.final) && is.finite(fitplus$nllh.final)) {
            
            goodness[i] <- max(0, parent_nll - (positive_nll(fitmoins$nllh.final) + positive_nll(fitplus$nllh.final)))
            direction[i] <- sign(positive_nll(fitmoins$nllh.final) - positive_nll(fitplus$nllh.final))
          } else {
            goodness[i] <- 0
          }
        } else {
          goodness[i] <- 0
        }
      } else {
        goodness[i] <- 0
      }
      
      if (goodness[i] < 0) print("Warning: negative goodness")
      if (is.infinite(goodness[i])) print("Warning: infinite goodness")
    }
    
    list(goodness = goodness, direction = ux[ord])
  }
}

GPRT_method_xi_in_R_plus <- list(eval = etempgprt, split = stempgprt, init = itempgprt)


# ================================================================
# Custom Penalized Pruning for GPRT
# ================================================================

# Compute penalized deviance criterion: sum of leaf deviances + A * (#leaves)
get_penalized_criterion <- function(tree, A) {
  frame <- tree$frame
  leaves <- frame[frame$var == "<leaf>", ]
  total_deviance <- sum(leaves$dev, na.rm = TRUE)
  n_leaves <- nrow(leaves)
  return(total_deviance + A * n_leaves)
}

# Generate all possible subtrees by pruning one internal node at a time
# until only the root remains
# Generate all possible pruned subtrees (using cp values)
generate_subtrees <- function(tree) {
  # All possible cp values give valid subtrees
  cp_values <- sort(unique(tree$cptable[, "CP"]))
  subtrees <- lapply(cp_values, function(cp) {
    prune(tree, cp = cp)
  })
  return(subtrees)
}


# Select the best subtree given penalty A
select_best_subtree <- function(tree, A) {
  subtrees <- generate_subtrees(tree)
  scores <- sapply(subtrees, function(st) get_penalized_criterion(st, A))
  best_index <- which.min(scores)
  return(list(tree = subtrees[[best_index]], A = A, score = scores[best_index]))
}

# Cross-validate penalty A over candidate grid
cross_validate_penalty <- function(tree, A_values) {
  results <- lapply(A_values, function(A) select_best_subtree(tree, A))
  scores <- sapply(results, function(res) res$score)
  best_result <- results[[which.min(scores)]]
  return(best_result)
}


# -------------------------------
# Function to Compute BIC
# -------------------------------
compute_bic <- function(tree, n_obs, params_per_leaf = 2) {
  frame <- tree$frame
  leaves <- frame[frame$var == "<leaf>", ]
  
  # Total deviance from terminal nodes
  total_deviance <- sum(leaves$dev, na.rm = TRUE)
  
  # Number of leaves (terminal nodes)
  n_leaves <- nrow(leaves)
  
  # Total number of parameters (e.g., 2 per leaf: shape and scale)
  k <- params_per_leaf * n_leaves
  
  # BIC formula
  bic <- total_deviance + log(n_obs) * k
  return(bic)
}


#---------------------------------------------------------------------
# STEP: Fit GPD per terminal node using node-specific threshold
#---------------------------------------------------------------------

gpd_return_levels_by_node <- function(pruned_tree, data, response_col = "PRECTOTCORR",
                                      q_thresh = 0.98, min_exc = 10,
                                      return_periods = c(10, 20, 50, 100, 150, 200),
                                      obs_per_year = 365) {
  
  terminal_nodes <- unique(pruned_tree$where)
  results <- list()
  
  for (node in terminal_nodes) {
    node_data <- data[pruned_tree$where == node, ]
    if (nrow(node_data) <= min_exc) next
    
    thr <- quantile(node_data[[response_col]], probs = q_thresh, na.rm = TRUE)
    exc_data <- node_data[[response_col]][node_data[[response_col]] > thr]
    n_exc <- length(exc_data)
    if (n_exc < min_exc) next
    
    # Fit GPD
    gpd_fit <- evir::gpd(node_data[[response_col]], threshold = thr)
    
    shape <- gpd_fit$par.ests["xi"]
    scale <- gpd_fit$par.ests["beta"]
    se_shape <- gpd_fit$par.ses["xi"]
    se_scale <- gpd_fit$par.ses["beta"]
    
    # 95% CI for shape parameter
    shape_lower95 <- shape - 1.96 * se_shape
    shape_upper95 <- shape + 1.96 * se_shape
    
    n_obs <- nrow(node_data)
    zeta_u <- n_exc / n_obs
    
    # Compute return levels + CI
    rl_df <- data.frame()
    for (R in return_periods) {
      T <- R * obs_per_year
      
      if (abs(shape) > 1e-6) {
        rl <- thr + (scale / shape) * ((T * zeta_u)^shape - 1)
        
        # Delta method: partial derivatives
        d_scale <- ((T * zeta_u)^shape - 1) / shape
        d_shape <- -(scale / shape^2) * ((T * zeta_u)^shape - 1) +
          (scale / shape) * ((T * zeta_u)^shape * log(T * zeta_u))
        
      } else {
        rl <- thr + scale * log(T * zeta_u)
        
        d_scale <- log(T * zeta_u)
        d_shape <- 0  # approx
      }
      
      # Variance approximation (ignoring covariance term)
      var_rl <- (d_scale^2) * (se_scale^2) + (d_shape^2) * (se_shape^2)
      se_rl <- sqrt(var_rl)
      
      rl_df <- rbind(rl_df, data.frame(
        node = node,
        n_obs = n_obs,
        threshold = thr,
        n_exceedances = n_exc,
        shape = shape,
        shape_lower95 = shape_lower95,
        shape_upper95 = shape_upper95,
        scale = scale,
        return_period = R,
        return_level = rl,
        rl_lower95 = rl - 1.96 * se_rl,
        rl_upper95 = rl + 1.96 * se_rl
      ))
    }
    
    results[[as.character(node)]] <- rl_df
  }
  
  return(do.call(rbind, results))
}