#########################################################################################################################################
###                                                        SET OF FUNCTIONS FOR THE SIMULATION                                        ###
#########################################################################################################################################

## Packages -----------------------------------------------------------------------------------------------------------------------------
if(!require(FSelector)){
  # Package used for feature selection
  devtools::install_version('FSelector', version = 'FSelctor', version = 0.34)
  library(FSelector)
}

if(!require(caret)){
  # Package used to calibrate models (and compute some indicators)
  devtools::install_version('caret', version = '6.0.94')
  library(caret)
}

if(!require(reshape2)){
  # Package use to format data table
  devtools::install_version('reshape2', version = '1.4.4')
}

if(!require(mlr3)){
  # Dependency package for mlr3verse used for the SVM-RFE method
  devtools::install_version('mlr3', version = '0.18.0')
  library(mlr3)
}

if(!require(mlr3verse)){
  # Package used for the SVM-RFE method
  devtools::install_version('mlr3verse')
  library(mlr3verse)
}

if(!require(glmnet)){
  # Package used for the LASSO feature selection method
  devtools::install_version('glmnet', version = '4.1-8')
  library(glmnet)
}

## Sources ------------------------------------------------------------------------------------------------------------------------------
source('functionGammaMetric.R')
source('functionGeneration.R')

## Feature selection methods ------------------------------------------------------------------------------------------------------------
# To facilitate parallel work, we put all the feature selection methods function into a list and R we go through that list to execute a 
# method and call the corresponding function. All the functions in that list have the same inputs and outputs. 
# Inputs:
#      f: The formula object we want to do selection for.
#   data: The dataset with the first column being the target feature y
# Outputs:

approachFunctions <- list(
  'FULL' <- function(f, data) colnames(data[, -1]),
  'BEST' <- function(f, data) colnames(data[, which(beta != 0)+1]),
  'CFS' <- function(f, data) FSelector::cfs(f, data),
  'CHI2' <- function(f, data) FSelector::cutoff.biggest.diff(FSelector::chi.squared(f, data)),
  'CONS' <- function(f, data) FSelector::consistency(f, data),
  'IG' <- function(f, data) FSelector::cutoff.biggest.diff(FSelector::information.gain(f, data)),
  'IGR' <- function(f, data) FSelector::cutoff.biggest.diff(FSelector::gain.ratio(f, data)),
  'GAMMA_BACK' <- function(f, data){
    # Parallel
    inner_cl <- parallel::makeCluster(2)
    doParallel::registerDoParallel(inner_cl)
    
    estim <- ifelse(min(table(data[, 1])) < ncol(data[, -1]), 'shrink', 'empiric')
    parallel::clusterEvalQ(inner_cl, source('functionGammaMetric.R'))
    
    evaluator <- function(attributes, data, dependent = 'y'){
      gammaMetric(X = data[,attributes], class = data[, dependent], covEstim = estim, plot = FALSE)
    }
    
    result <- FSelectorRcpp::feature_search(
      attributes = names(data[, -1]), 
      fun = evaluator, 
      data = data, 
      mode = 'greedy',
      type = 'backward'
    )
    
    var_select <- names(data[,-1])[which(result$best == 1)]
    
    parallel::stopCluster(inner_cl)
    return(var_select)
  },
  'GAMMA_BF' <- function(f, data){
    # Which estimator of the variance-covariance matrix to use
    estim <- ifelse(min(table(data[, 1])) < ncol(data[, -1]), 'shrink', 'empiric')
    
    # Feature selection with the best-first search direction
    var_select <- FSelector::best.first.search(colnames(data[, -1]), eval.fun = function(subset){
      # Computation of the gamma-metric value for the candidate subset of features
      g <- gammaMetric(data[, subset], class = data[, 1], covEstim = estim, plot = FALSE)
      return(g)
    })
    
    return(var_select)
  },
  'GAMMA_FORW' <- function(f, data){
    # Which estimator of the variance-covariance matrix to use
    estim <- ifelse(min(table(data[, 1])) < ncol(data[, -1]), 'shrink', 'empiric')
    
    # Feature selection with the forward search direction
    var_select <- FSelector::forward.search(colnames(data[, -1]), eval.fun = function(subset){
      # Computation of the gamma-metric value for the candidate subset of features
      g <- gammaMetric(data[, subset], class = data[, 1], covEstim = estim, plot = FALSE)
      return(g)
    })
    
    return(var_select)
  },
  'GAMMA_HC' <- function(f, data){
    # Which estimator of the variance-covariance matrix to use
    estim <- ifelse(min(table(data[, 1])) < ncol(data[, -1]), 'shrink', 'empiric')
    
    # Feature selection with the hill-climbing search direction
    var_select <- FSelector::hill.climbing.search(colnames(data[, -1]), eval.fun = function(subset){
      # Computation of the gamma-metric value for the candidate subset of features
      g <- gammaMetric(data[, subset], class = data[, 1], covEstim = estim, plot = FALSE)
      return(g)
    })
    
    return(var_select)
  },
  'LASSO' <- function(f, data){
    # Calibration of the lambda values
    cv.lasso <- glmnet::cv.glmnet(x = as.matrix(data[,-1]), y = data[,'y'], family = 'binomial', alpha = 1, parallel = FALSE, standardize = TRUE, type.measure = 'auc')
    
    # Get the beta coefficient values
    df_coef <- round(as.matrix(coef(cv.lasso, s = cv.lasso$lambda.min)), 2)
    
    # Get all the contributing variables
    var_select <- names(df_coef[df_coef[,1] != 0,])[-1]
    
    return(var_select)
  },
  'ONER' <- function(f, data) FSelector::cutoff.biggest.diff(FSelector::oneR(f, data)),
  'RFI' <- function(f, data) FSelector::cutoff.biggest.diff(FSelector::random.forest.importance(f, data)),
  'RELIEF' <- function(f, data) FSelector::cutoff.biggest.diff(FSelector::relief(f, data)),
  'STEP' <- function(f, data){
    df <- data

    # Model
    all <- glm(f, family = 'binomial', data = df)
    intercept_only <- glm(y~1, family = 'binomial', data = df)
    
    
    # Stepwise selection with AIC
    both <- step(intercept_only, direction = 'both', scope = formula(all), trace = 0)
    
    
    # Selected features
    var_select <- names(both$coefficients)[-1]
    
    return(var_select)
  },
  'SU' <- function(f, data) FSelector::cutoff.biggest.diff(FSelector::symmetrical.uncertainty(f, data)),
  'SVM-RFE' <- function(f, data){
    # To ensure that y is a factor
    data$y <- factor(data$y)
    
    # Definition of the feature selection
    optimizer <- mlr3verse::fs('rfe', n_features = 1, feature_fraction = 0.9)
    
    # Definition of the learner
    learner <- mlr3::lrn("classif.svm", type = 'C-classification', kernel = 'linear', predict_type = 'prob')
    
    # Definition of the task
    task <- mlr3::TaskClassif$new(id = 'simulated_data', 
                                  backend = data, 
                                  target = 'y', 
                                  positive = '1')
    
    # Add the task to the dictionary
    mlr3::mlr_tasks$add('simulated_data', task)
    
    task_gc <- mlr3::tsk('simulated_data')
    task_gc$col_roles$stratum <- 'y'
    cv5 <- mlr3::rsmp('cv', folds = 4)
    cv5$instantiate(task_gc)
    
    # Instance
    instance <- mlr3verse::fsi(task = task_gc,
                               learner = learner, 
                               resampling = cv5,
                               measures = mlr3::msr('classif.auc'),
                               terminator = mlr3verse::trm('none'),
                               callback = mlr3verse::clbk('mlr3fselect.svm_rfe')
                               )
    
    optimizer$optimize(inst = instance)
    var_select <- instance$result_feature_set
    
    return(var_select)
  }

)


## Simulation iteration -----------------------------------------------------------------------------------------------------------------
# A single function to produce the results of the simulation for a given feature selection method and given train/test data set.
performSimulationIteration <- function(X_train, Y_train, X_test, Y_test, beta, fs){
  ## Inputs: 
  ##    X_train: Matrix of observations X for the train data set (without intercept)
  ##    Y_train: Vectors of labels of the target feature for the train data set
  ##     X_test: Matrix of observations X for the test data set (without intercept)
  ##     Y_test: Vectors of labels of the target feature for the test data set
  ##       beta: Beta vectors without beta0. It was only used to identify quickly informative features
  ##         fs: Name of the feature selection method to apply
  ## Outputs: 
  ##    A dataframe resIteration with multiples column corresponding to the results of one iteration of the simulation on the train and test data set.
  ##      * Approach: Name of the feature selection methods 
  ##      * NbVarSelected: Total number of features selected by the method
  ##      * NbVarInf: Total number of informative features selected 
  ##      * Accuracy_selection: The accuracy of the selection while we defined a True Positive as the feature is informative and is selected by the method
  ##      * Specificity_selection: The specificity of the selection 
  ##      * Sensitivity_selection: The sensitivity of the selection
  ##      * MCC_selection: The MCC of the selection
  ##      * Accuracy_training: The accuracy computed on the adjusted values of the training sample
  ##      * Specificity_training: The specificity computed on the adjusted values of the training sample
  ##      * Sensitivity_training: The sensitivity computed on the adjusted values of the training sample
  ##      * MCC_training: The MCC computed on the adjusted values of the training sample
  ##      * Accuracy_test: The accuracy computed on the predicted values of the test sample
  ##      * Specificity_test: The specificity computed on the predicted values of the test sample
  ##      * Sensitivity_test: The sensitivity computed on the predicted values of the test sample
  ##      * MCC_test: The MCC computed on the predicted values of the test sample
  ##      * VarSelected: The complete list of features selected by the method. To fit everything in a dataframe we pasted all the names of the features
  ##                        in a string separated by ';'
  ##      * coef_estimate: This is the beta coefficients estimated by the logistic regression model. As for VarSelected, we pasted all the values 
  ##                        in a string separated by ';'
  ##      * std_estimate: This is the standard deviation of the estimation of the beta coefficient computed with the logistic regression model. As 
  ##                        for VarSelected, we pasted all the values in a string separated by ';'
  
  # Feature selection step
  runtime <- as.numeric(system.time({var_select <- do.call(fs, list(f = FSelector::as.simple.formula(colnames(X_train), 'y'), data = data.frame(y = factor(Y_train), X_train)))})['elapsed'])
  
  # Binary vector with 1 being the feature is informative and 0 the feature is non-informative
  non_zero_coefficients <- replace(rep(0, length(beta)), beta != 0, 1)
  
  # Type of features selected by the method
  selected <- (colnames(X_train) %in% var_select)*1
  
  # Computation of the 2x2 confusion matrix for selection
  cm_sel <- caret::confusionMatrix(data = factor(selected, levels = c('0', '1')), reference = factor(non_zero_coefficients), positive = '1')
  
  # Computation of the performance indicators
  spe_sel <- cm_sel$byClass['Specificity']
  sen_sel <- cm_sel$byClass['Sensitivity']
  acc_sel <- cm_sel$overall['Accuracy']
  mcc_sel <- mltools::mcc(preds = factor(selected, levels = c('0', '1')),
                          actuals = factor(non_zero_coefficients))
  
  # Calibration of the Logistic model with the selected features
  training_time <- as.numeric(system.time({
    ctr <- caret::trainControl(method = 'none', allowParallel = FALSE)
    mod.fit <- caret::train(FSelector::as.simple.formula(var_select, 'y'), 
                            data = data.frame(y = factor(Y_train), X_train), 
                            trControl = ctr, 
                            method = 'glm', 
                            family = 'binomial')
  })['elapsed'])
  
  # Prediction of the test sample
  test_time <- as.numeric(system.time({
    preds <- predict(mod.fit, data.frame(y = Y_test, X_test))
  })['elapsed'])
  
  # Classification performances (test sample)
  cm <- caret::confusionMatrix(preds, factor(Y_test), positive = '1')
  acc_test <- cm$overall['Accuracy']
  spe_test <- cm$byClass['Specificity']
  sen_test <- cm$byClass['Sensitivity']
  mcc_test <- mltools::mcc(preds, factor(Y_test))
  
  # AUC 
  results_roc <- pROC::roc(as.numeric(Y_test), as.numeric(preds), quiet = TRUE)
  auc <- as.numeric(results_roc$auc)
  
  # True Positive Rate (TPR)
  TP <- cm$table[2, 2]
  FN <- cm$table[1, 2]
  TPR <- TP/(TP+FN)
  
  # True Negative Rate (TNR)
  TN <- cm$table[1, 1]
  FP <- cm$table[2, 1]
  TNR <- TN/(TN + FP)
  
  # Mean Absolute Error (MAE)
  proba_preds <- predict(mod.fit, X_test, type = 'prob')[,2]
  MAE <- mean(abs(as.numeric(Y_test) - proba_preds))
  
  # Maximum Youden index J and corresponding sensitivity/specificity
  roc_preds <- pROC::roc(as.numeric(Y_test), proba_preds, quiet = TRUE)
  youden_se_spe <- pROC::coords(roc_preds, "best", best.method = 'youden')
  
  # In case of equal values we return the performance for the youden index which maximize the sensitivity
  youden_se_spe <- youden_se_spe[which.max(youden_se_spe[, 'sensitivity']),]
  
  # Aggregation of the results
  localRes <- data.frame(Approach = fs, NbVarSelected = length(var_select), NbVarInf = sum(non_zero_coefficients & selected), 
                         Accuracy_selection = acc_sel, Specificity_selection = spe_sel, Sensitivity_selection = sen_sel, MCC_selection = mcc_sel, 
                         Accuracy_test = acc_test, Specificity_test = spe_test, Sensitivity_test = sen_test, MCC_test = mcc_test, AUC = auc,
                         Youden = as.numeric(youden_se_spe['threshold']), YoudenSensitivity = as.numeric(youden_se_spe['sensitivity']), YoudenSpecificity = as.numeric(youden_se_spe['specificity']),
                         TPR = TPR, TNR = TNR, MAE = MAE, Runtime = runtime, TrainingTime = training_time, TestTime = test_time,
                         VarSelected = paste(var_select, collapse = ';'))
  return(localRes)
}

## Indicator function -------------------------------------------------------------------------------------------------------------------
## Function to compute the jaccard index between two sets
jaccardFunction <- function(set1, set2){
  # Inputs: 
  # - set1: the first vector of variable names 
  # - set2: the second vector of variable names
  # Outputs: 
  # A numeric value corresponding to the jaccard index between set1 and set2
  
  intersection <- length(intersect(set1, set2))
  union <- length(set1) + length(set2) - intersection
  return(intersection/union)
}

## Function to compute pairwise jaccard index for all sets of selected features from a single method. 
stabilityFunction <- function(list_of_selected_var){
  ## Inputs: 
  # list_of_selected_var: A list with the features selected by a single method for each repetition of the FSM process 
  # Outputs: 
  # The jaccard index score for pairwise sets are computed then averaged other all combinations
  
  # Number of sets (should be equal to the number of repetitions)
  m <- length(list_of_selected_var)
  
  # Jaccard index for all combinations of 2 sets of selected features
  jaccard_index <- combn(x = 1:m, 2, FUN = function(combinaison){
    # Create all combination of 2 from 1 to m. 
    jaccardFunction(set1 = list_of_selected_var[[combinaison[1]]], set2 = list_of_selected_var[[combinaison[2]]])
  })
  
  # Stability index
  stability <- as.numeric((2/(m*(m-1)))*sum(jaccard_index))
  return(stability)
}

## TOPSIS function ----------------------------------------------------------------------------------------------------------------------
# Implementation of the Technique for Order preference by similarity to ideal solution (TOPSIS)

topsisScore <- function(xij, benefit_criterion, cost_criterion, wi){
  ## Inputs: 
  # xij: the matrix of criterion values for each feature selection method. Feature selection methods as rows and columns as the criterion.
  #       Intersection is the value of the ith FSM at the jth criteria
  # benefit_criterion: a vector of criterion that are considered as benefit (the greater the better).
  # cost_criterion: a vector of criterion that are considered as cost (the lower the better).
  # Output: 
  # A data frame with the score of each FSM and their ranking
  
  # Step 1: Calculate the normalized criteria values
  rij <- t(apply(xij, 1, FUN = function(values){
    values/sqrt(sum(values^2))
  }))
  
  # Step 2: Calculate the weighted criteria values:
  vij <- wi*rij
  
  # Step 3: Find the ideal solution S+ and the negative ideal solution S-
  S_plus <- sapply(rownames(vij),FUN = function(criteria){
    if(criteria %in% benefit_criterion){
      max(vij[criteria,])
    }else{
      min(vij[criteria,])
    }
  })
  
  S_moins <- sapply(rownames(vij), FUN = function(criteria){
    if(criteria %in% benefit_criterion){
      min(vij[criteria,])
    }else{
      max(vij[criteria,])
    }
  })
  
  # Step 4: Calculate the Euclidean distance between the real and ideal solutions, and between the real solution and negative ideal solution
  D_plus <- apply(vij, 2, FUN = function(values){
    sqrt(sum((values-S_plus)^2))
  })
  
  D_moins <- apply(vij, 2, FUN = function(values){
    sqrt(sum((values - S_moins)^2))
  })
  
  # Step 5: Calculate R_j
  R_j <- D_moins/(D_plus + D_moins)
  
  # Step 6: Rank feature selection methods by maximizing Rj
  topsis <- data.frame(t(xij),
             Rank = rank(-R_j),
             Score = R_j)
  
  return(topsis)
}
topsisRanking <- function(RES, separation, repartition, criterion, benefit = benefit_criterion, cost = cost_criterion, rankings = Rankings){
  
  # Creation of the xij matrix
  xij <- t(RES[which((RES$Separation == separation) & (RES$Repartition == repartition)), criterion])
  colnames(xij) <- RES[which((RES$Separation == separation) & (RES$Repartition == repartition)), "Approach"]
  
  # Weights
  wi <- sapply(rankings, FUN = function(rank){
    n <- length(rankings)
    wi <- (n-rank+1)/sum(n - rankings + 1)
    return(wi)
  })
  
  # Topsis score and ranks
  topsis <- topsisScore(xij, benefit, cost, wi)
  return(topsis)
}

## Formating functions ------------------------------------------------------------------------------------------------------------------
# dataframeSelectedFeaturesScenario2 <- function(res_table){
#   ## Inputs: 
#   ##    res_table: The results table of the simulations.
#   ## Outputs: 
#   ##    VARSELECTED_LONG: The dataframe contains for each method the number of times each feature is selected 
#   
#   # Number of non-informative features
#   p_prime <- res_table$NbVar[1] - n_var_inf
#   
#   # data.frame with how much time each feature is selected
#   VARSELECTED <- data.frame(sapply(Approach, FUN = function(x){
#     approach_fs <- which(res_table$Approach == x)
#     list_var_selected <- unlist(strsplit(res_table$VarSelected[approach_fs], ';'))
#     indices <- match(c(paste0('x', 1:n_var_inf), paste0('Noise', 1:p_prime)), names(table(list_var_selected)))
#     table(list_var_selected)[indices]
#   }))
#   
#   # The name for SVM-RFE is modified automatically to SVM.RFE in the results, we need to fix this
#   colnames(VARSELECTED)[colnames(VARSELECTED) == 'SVM.RFE'] <- 'SVM-RFE'
#   
#   # Name of the features
#   VARSELECTED$Feature <- c(paste0('x', 1:n_var_inf), paste0('N', 1:p_prime))
#   
#   # Long format
#   VARSELECTED_LONG <- reshape2::melt(VARSELECTED, id.vars = 'Feature')
#   
#   return(VARSELECTED_LONG)
# }
# 
# dataframeSelectedFeaturesScenario3 <- function(res_table){
#   # Function to return how much each feature of Scenario 3 is selected by the feature selection methods
#   VARSELECTED <- data.frame(sapply(Approach, FUN = function(x){
#     approach_fs <- which(res_table$Approach == x)
#     list_var_selected <- unlist(strsplit(res_table$VarSelected[approach_fs], ';'))
#     indices <- match(paste0('X', 1:(s_g*n_g)), names(table(list_var_selected)))
#     table(list_var_selected)[indices]
#   }))
#   
#   # The name for SVM-RFE is modified automatically to SVM.RFE in the results, we need to change that
#   colnames(VARSELECTED)[colnames(VARSELECTED) == 'SVM.RFE'] <- 'SVM-RFE'
#   
#   # Name of the features
#   VARSELECTED$Feature <- paste0('X', 1:(s_g*n_g))
#   
#   # Long format 
#   VARSELECTED_LONG <- reshape2::melt(VARSELECTED, id.vars = 'Feature')
#   
#   return(VARSELECTED_LONG)
# }


