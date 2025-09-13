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







































































#' 
#' 
#' 
#' prune_gprt_custom <- function(tree, train_data, test_data, A_values = seq(0, 100, by = 5), verbose = TRUE) {
#'   compute_CA <- function(tree, data, A) {
#'     frame <- tree$frame
#'     leaf_nodes <- which(frame$var == "<leaf>")
#'     preds <- rep(NA, nrow(data))
#'     
#'     for (i in leaf_nodes) {
#'       node_id <- as.numeric(rownames(frame)[i])
#'       node_obs <- which(tree$where == node_id)
#'       
#'       shape <- frame$yval2[i, 1]
#'       scale <- frame$yval2[i, 2]
#'       threshold <- data$threshold_col[node_obs]
#'       rainfall <- data$PRECTOTCORR[node_obs]
#'       
#'       exceedances <- rainfall - threshold
#'       exceedances <- exceedances[exceedances > 0]
#'       
#'       if (length(exceedances) > 0 && shape > 0 && scale > 0) {
#'         nll <- -sum(log(1 + shape * exceedances / scale)) / shape + length(exceedances) * log(scale)
#'       } else {
#'         nll <- Inf
#'       }
#'       
#'       preds[node_obs] <- nll
#'     }
#'     
#'     total_loss <- sum(preds, na.rm = TRUE)
#'     n_leaves <- length(leaf_nodes)
#'     CA <- total_loss + A * n_leaves
#'     return(CA)
#'   }
#'   
#'   best_tree <- NULL
#'   best_CA <- Inf
#'   best_A <- NA
#'   
#'   for (A in A_values) {
#'     current_tree <- tree
#'     current_CA <- compute_CA(current_tree, test_data, A)
#'     
#'     repeat {
#'       frame <- current_tree$frame
#'       leaf_nodes <- which(frame$var == "<leaf>")
#'       if (length(leaf_nodes) <= 1) break
#'       
#'       leaf_deviances <- frame$dev[leaf_nodes]
#'       worst_leaf <- leaf_nodes[which.min(leaf_deviances)]
#'       
#'       new_cp <- frame$complexity[worst_leaf] + 1e-4
#'       pruned <- try(prune(current_tree, cp = new_cp), silent = TRUE)
#'       if (inherits(pruned, "try-error")) break
#'       
#'       new_CA <- compute_CA(pruned, test_data, A)
#'       
#'       if (new_CA < current_CA) {
#'         current_tree <- pruned
#'         current_CA <- new_CA
#'       } else {
#'         break
#'       }
#'     }
#'     
#'     if (verbose) cat("A =", A, "→ C_A =", round(current_CA, 2), "\n")
#'     
#'     if (current_CA < best_CA) {
#'       best_CA <- current_CA
#'       best_tree <- current_tree
#'       best_A <- A
#'     }
#'   }
#'   
#'   if (verbose) cat("Selected A:", best_A, "\n")
#'   return(best_tree)
#' }
#' 
#' 
#' 
#' #' 
#' #' 
#' #' 
#' #' #Functions for pruning
#' #' 
#' #' #'@param y a numeric vector
#' #' #'@param sh a shape parameter
#' #' #'@param sc a scale parameter
#' #' #'@return goodness the negative log likelihood computed on y according to a generalized pareto distribution with paramters shape, scale and location=0
#' #' LL=function(y,sh,sc){
#' #'   temp=evir::dgpd(y,xi=sh,beta=sc,mu=0)
#' #'   goodness=-(sum(log(temp)))
#' #'   if(is.na(goodness)){
#' #'     #print("Warning: unable to compute the LL for a subset of the test sample in a leaf")
#' #'     goodness=+Inf
#' #'   }
#' #'   return(goodness)
#' #' }
#' #' 
#' #' #'@param arbre a tree of the class rpart
#' #' #'@param train_data a test data sample
#' #' #'@return a vector of the x relative errors of the tree, ie the ratio between the sum of the negative log likelihood at each leaf and the root negative log likelihood, for each tree pruned thanks to the complexity parameter and with the test sample
#' #' construct_x_rel_error_alpha_gprt=function(arbre,test_data,name_response,alpha){
#' #'   #tree=prune(arbre,cp=alpha/arbre$frame[1,"dev"])
#' #'   tree=prune_tree(rpart_object=arbre,CP=alpha)
#' #' 
#' #'   y_val_tree=matrix(c(1:dim(tree$frame)[1],
#' #'                       as.numeric(row.names(tree$frame)),
#' #'                       as.numeric(as.matrix(tree$frame)[,c("yval2.1","yval2.2")])),
#' #'                     nrow=dim(tree$frame)[1],ncol=4,byrow=FALSE)
#' #'   colnames(y_val_tree)=c("where","node","sh","sc")
#' #' 
#' #'   number_of_terminal_nodes=sum(tree$frame[,"var"]=="<leaf>")
#' #' 
#' #'   where_test_data=as.vector(rpart.predict.leaves(tree, test_data, type = "where"))
#' #'   unique_where_test_data=sort(unique(where_test_data))
#' #'   unique_where_test_data=matrix(c(unique_where_test_data,rep(0,times=length(unique_where_test_data))),ncol=2,byrow=FALSE)
#' #'   colnames(unique_where_test_data)=c("where","NLL")
#' #' 
#' #'   for(i in 1:dim(unique_where_test_data)[1]){
#' #'     where=unique_where_test_data[i]
#' #'     y=test_data[where_test_data==where,name_response]
#' #'     # unique_where_test_data[i,2]=positive_nll(LL(y,
#' #'     #                                sh=y_val_tree[y_val_tree[,"where"]==where,"sh"],
#' #'     #                                sc=y_val_tree[y_val_tree[,"where"]==where,"sc"]))
#' #'     unique_where_test_data[i,2]=LL(y,
#' #'                                    sh=y_val_tree[y_val_tree[,"where"]==where,"sh"],
#' #'                                    sc=y_val_tree[y_val_tree[,"where"]==where,"sc"])
#' #'   }
#' #' 
#' #'   return(sum(unique_where_test_data[,"NLL"]))
#' #' 
#' #' }
#' #' 
#' #' 
#' #' #'@param arbre a tree of the class rpart
#' #' #'@param train_data a test data sample
#' #' #'@return a vector of the x relative errors of the tree, ie the ratio between the sum of the negative log likelihood at each leaf and the root negative log likelihood, for each tree pruned thanks to the complexity parameter and with the test sample
#' #' cross_val_gprt=function(data,formula,arbre,n_fold,seed,choice="first increase"){
#' #'   
#' #'   #On the tree
#' #'   
#' #'   alpha=as.numeric(get_table(arbre)[,"CP"])
#' #'   
#' #'   if(alpha[length(alpha)]==0){
#' #'     alpha=alpha[1:(length(alpha)-1)]
#' #'   }
#' #'   
#' #'   beta=sort(c(1,seq(alpha[length(alpha)],alpha[1],length.out=1000),0),decreasing=TRUE)
#' #'   
#' #'   #On the n_fold trees
#' #'   
#' #'   set.seed(seed)
#' #'   xgroup=sample(1:n_fold,size=dim(data)[1],replace=TRUE)
#' #'   
#' #'   
#' #'   list_train_tree=list()
#' #'   Risk_matrix=matrix(nrow = length(beta),ncol = n_fold)
#' #'   
#' #'   for(k in 1:n_fold){
#' #'     train_data=data[xgroup!=k,]
#' #'     test_data=data[xgroup==k,]
#' #'     
#' #'     tree_train=rpart(data = train_data,
#' #'                      formula = formula,
#' #'                      method = GPRT_method_xi_in_R_plus,
#' #'                      control = rpart.control(cp=0,minsplit=50, minbucket=30, maxcompete=0, maxsurrogate=0),
#' #'                      xval=0)
#' #'     
#' #'     list_train_tree[[k]]=tree_train
#' #'     
#' #'     # for(l in 1:length(beta)){
#' #'     #   Risk_matrix[l,k]=construct_x_rel_error_alpha_gprt(arbre = tree_train,test_data = test_data,name_response = as.character(formula[[2]]),alpha = beta[l])
#' #'     # }
#' #'     Risk_matrix[,k]=sapply(X = beta,FUN="construct_x_rel_error_alpha_gprt",arbre = tree_train,test_data = test_data,name_response = as.character(formula[[2]]))
#' #'     
#' #'     
#' #'   }
#' #'   
#' #'   Risk_beta=rowSums(Risk_matrix)
#' #'   #res=beta[Risk_beta==min(Risk_beta)]
#' #'   #We choose the deepth just before which the test error begin to increase
#' #'   pos=match(TRUE,Risk_beta[-1]>Risk_beta[-length(Risk_beta)])
#' #'   if(is.na(pos) | choice=="min"){
#' #'     #If the test error does not increase, we choose the minimal deepth that minimize the test error
#' #'     depth=beta[match(TRUE,Risk_beta==min(Risk_beta))][1]
#' #'   }else{
#' #'     depth=beta[pos][1]
#' #'   }
#' #'   
#' #'   res=list(depth,beta,Risk_beta)
#' #'   
#' #'   return(res)
#' #' }
#' #' 
#' #' 
#' #' 
#' #' #'@param arbre a tree of the class rpart
#' #' #'@param train_data a test data sample
#' #' #'@return a vector of the x relative errors of the tree, ie the ratio between the sum of the negative log likelihood at each leaf and the root negative log likelihood, for each tree pruned thanks to the complexity parameter and with the test sample
#' #' prune_train_test_gprt=function(data_train,data_test,formula,arbre_apprentissage){
#' #'   
#' #'   #On the tree
#' #'   
#' #'   alpha=c(as.numeric(arbre_apprentissage$cptable[,"CP"])*arbre_apprentissage$frame[1,"dev"])
#' #'   
#' #'   Risk_alpha=double(length(alpha))
#' #'   
#' #'   for(l in 1:length(alpha)){
#' #'     Risk_alpha[l]=construct_x_rel_error_alpha_gprt(arbre = arbre_apprentissage,test_data = data_test,name_response = as.character(formula[[2]]),alpha = alpha[l])
#' #'   }
#' #'   
#' #'   #Compute the risk for each beta and each fold
#' #'   
#' #'   
#' #'   h=Risk_alpha[1]-min(Risk_alpha,na.rm=TRUE)
#' #'   #res = max(match(TRUE, Risk_alpha <= min(Risk_alpha,na.rm=TRUE)+1/2*h)-1,1)
#' #'   res = match(TRUE, Risk_alpha == min(Risk_alpha,na.rm=TRUE))
#' #'   
#' #'   #res=alpha[match(TRUE,Risk_alpha==min(Risk_alpha))]
#' #'   
#' #'   return(res)
#' #' }
#' #' 
#' #' 
#' #' positive_nll=function(x){
#' #'   if(x>0){
#' #'     return(x+1)
#' #'   }else{
#' #'     return(exp(x))
#' #'   }
#' #' }
#' #' 
#' #' ll_ratio_test_prune_gprt=function(arbre){
#' #'   ll=c()
#' #'   n_param=c()
#' #'   depth=dim(arbre$cptable)[1]
#' #'   for(i in 1:depth){
#' #'     sub_tree=prune(arbre,cp=arbre$cptable[i,"CP"])
#' #'     ll=c(ll,-sum(sub_tree$frame[sub_tree$frame$var=="<leaf>","dev"]))
#' #'   }
#' #'   n_param=(arbre$cptable[,"nsplit"]+1)*2
#' #'   valeurs=ll[2:length(ll)]-ll[1:(length(ll)-1)]
#' #'   seuils=qchisq(p = 0.95,df = n_param[2:length(n_param)]-n_param[1:(length(n_param)-1)])
#' #'   decision=valeurs>seuils
#' #'   level=sum(cumsum(decision)==1:length(decision))+1
#' #'   alpha=arbre$cptable[level,"CP"]*arbre$frame[1,"dev"]
#' #'   return(alpha)
#' #' }
#' #' 
#' #' 
#' #' get_table=function(rpart_object){
#' #'   n_split=1:nrow(rpart_object$splits)
#' #'   var=rownames(rpart_object$splits)
#' #'   improve=rpart_object$splits[,"improve"]
#' #'   n=rpart_object$splits[,"count"]
#' #'   CP=(rpart_object$frame[1,"dev"]-cumsum(rpart_object$splits[,"improve"]))/rpart_object$frame[1,"dev"]
#' #'   return(cbind(n_split,CP,improve,var,n))
#' #' }
#' #' 
#' #' match_id_id_bis=function(id_bis,table,frame){
#' #'   frame$improve=NA
#' #'   for(i in 1:nrow(frame)){
#' #'     if(frame$var[i]!="<leaf>"){
#' #'       id=row.names(frame)[i]
#' #'       frame$improve[i]=frame$dev[row.names(frame) == as.character(as.numeric(id))]-(frame$dev[row.names(frame) == as.character(2*as.numeric(id))] + frame$dev[row.names(frame) == as.character(2*as.numeric(id) + 1)])
#' #'     }
#' #'   }
#' #'   
#' #'   id=rownames(frame)[frame$var==table[id_bis,"var"] & frame$n==table[id_bis,"n"]  &  round(as.numeric(frame$improve),digits=6)==round(as.numeric(table[id_bis,"improve"]),digits=6)]
#' #'   return(as.numeric(id))
#' #' }
#' #' 
#' #' prune_tree=function(rpart_object,CP){
#' #'   frame <- rpart_object$frame
#' #'   
#' #'   COMPUTE=FALSE
#' #'   if(length(rpart_object$splits)>0){
#' #'     if(nrow(rpart_object$splits)>0){
#' #'       COMPUTE=TRUE
#' #'     }
#' #'   }
#' #'   if(COMPUTE){
#' #'     
#' #'     table=get_table(rpart_object)
#' #'     id_bis_splits_to_remove=table[,"n_split"][table[,"CP"]<CP]
#' #'     
#' #'     if(length(id_bis_splits_to_remove)==nrow(table)){
#' #'       return(prune(rpart_object,cp=Inf))
#' #'     }else{
#' #'       if(length(id_bis_splits_to_remove)==0){
#' #'         return(rpart_object)
#' #'       }else{
#' #'         id_splits_to_remove=sapply(as.numeric(id_bis_splits_to_remove),FUN="match_id_id_bis",table=table,frame=frame)
#' #'         return(snip.rpart(rpart_object, id_splits_to_remove))
#' #'       }
#' #'     }
#' #'     
#' #'   }else{
#' #'     return(rpart_object)
#' #'   }
#' #'   
#' #' }
