###########################################################################################################
#                                                 SCENARIO 1                                              #
###########################################################################################################

## Environment --------------------------------------------------------------------------------------------

if(!require(ggplot2)){
  # Package used to create figures
  devtools::install_version('ggplot2', version = '3.4.3')
  library(ggplot2)
}

if(!require(doRNG)){
  # Package used for reproductibility of parallel executions (seed)
  devtools::install_version('doRNG', version = '1.8.6')
  library(doRNG)
}
# 
if(!require(ggpubr)){
  # Package used for grids in plots
  devtools::install_version('ggpubr', version = '0.6.0')
  library(ggpubr)
}
#### Sources --------------------------------------------------------------------------------------------
before <- ls()
source('functionGammaMetric.R')
source('functionGeneration.R')
source('functions.R')
after <- ls()

new_functions <- setdiff(after, before)
new_functions <- new_functions[-which(new_functions == 'before')]

#### Parameters of the generation -----------------------------------------------------------------------
N <- 2000                                 # Number of observations
n_var_inf <- 3                            # Number of informative features
n_noise <- 22                             # Number of non-informative features
beta <- c(3, -2, 0.5, rep(0, n_noise))    # Beta coefficients for the generation process
beta0 <- 0                                # Intercept value
R <- 50                                   # Number of repetitions
set.seed(123)                             # Seed of the simulation

# Feature selection methods 
Approach <- c('FULL', 'CFS', 'CHI2', 'GAMMA_BACK', 'GAMMA_BF', 'GAMMA_FORW', 'LASSO', 'RFI', 'STEP', 'SU', 'SVM-RFE')



# EXECUTE TO RUN THE SIMULATIONS OR GO TO UPLOADING THE RESULTS
###########################################################################################################################################################################
###########################################################################################################################################################################

#### Parameters of the parallel computations ------------------------------------------------------------

# Number of clusters
n_clusters <- parallel::detectCores() - 1
cl <- parallel::makeCluster(n_clusters)
doSNOW::registerDoSNOW(cl)
time <- vector('numeric', R)
res <- NULL

# Export functions
parallel::clusterExport(cl, varlist = c("N", "n_var_inf", "n_noise", "beta", "beta0", new_functions))

# Progress bar
pb <- txtProgressBar(min = 0, max = R*length(Approach), style = 3)
progress <- function(n) setTxtProgressBar(pb, n + iteration_count)
opts <- list(progress = progress)
iteration_count <<- 0

#### BEGINING OF THE SIMULATION -------------------------------------------------------------------------

for(i in 1:R){
  # Time of computation
  t0 <- Sys.time()
  
  # Generation of the training data
  dat_train <- generation(n = N, n_var_inf = n_var_inf, n_noise = n_noise, beta = c(beta0, beta))
  X_train <- as.matrix(dat_train[, -1])
  Y_train <- as.matrix(dat_train[, 1])
  
  # Generation of the test data
  dat_test <- generation(n = N, n_var_inf = n_var_inf, n_noise = n_noise, beta = c(beta0, beta))
  X_test <- as.matrix(dat_test[, -1])
  Y_test <- as.matrix(dat_test[, 1])
  
  # Perform feature selection, modelisation and prediction in parallel for each feature selection methods
  resIteration <- foreach(fs = Approach, .combine = 'rbind', .export = ls()[!ls() %in% c('performSimulationIteration', 'beta', 'X_train', 'Y_train', 'X_test', 'Y_test')], .options.snow = opts, .verbose = FALSE) %dorng% {
    localRes <- performSimulationIteration(X_train, Y_train, X_test, Y_test, beta = beta, fs)
    row.names(localRes) <- NULL
    return(localRes)
  }
  
  # Updates
  t1 <- Sys.time()
  iteration_count <<- i*length(Approach)
  res <- rbind(res, resIteration)
  time[i] <- as.numeric(difftime(t1, t0, units = 'secs'))
}
close(pb)
parallel::stopCluster(cl)

#### Saving the results ----------------------------------------------------------------------------------
file_path <- paste0('Scenario1/res_', R, '_iterations.txt')
write.table(res, file_path)

############################################################################################################################################################################
############################################################################################################################################################################

#### Uploading the results for graphics and tables -------------------------------------------------------
res <- read.table('Scenario1/res_50_iterations.txt')

# Formating
res <- data.frame(res)

# Aggregated results
RES <- data.frame(Approach = Approach)

for(i in 1:nrow(RES)){
  # Line
  lignes <- res[which(res$Approach == RES[i, 'Approach']), ]
  
  # Number of features
  RES[i, 'NbVarSelected'] <- mean(lignes[, 'NbVarSelected'])
  RES[i, 'NbVarInf'] <- mean(lignes[, 'NbVarInf'])
  RES[i, 'NbVarNonInf'] <- RES[i, 'NbVarSelected'] - RES[i, 'NbVarInf']

  # Time indicators
  RES[i, 'Runtime'] <- mean(lignes[, 'Runtime'])
  RES[i, 'TrainingTime'] <- mean(lignes[, 'TrainingTime'])
  RES[i, 'TestTime'] <- mean(lignes[, 'TestTime'])
  
  # Global indicator
  RES[i, 'AUC'] <- mean(lignes[, 'AUC'])*100
  RES[i, 'Youden'] <- mean(lignes[,'Youden'])*100
  RES[i, 'YoudenSensitivity'] <- mean(lignes[, 'YoudenSensitivity'])*100
  RES[i, 'YoudenSpecificity'] <- mean(lignes[, 'YoudenSpecificity'])*100
  
  # Selection indicators
  RES[i, 'Accuracy_selection'] <- mean(lignes[, 'Accuracy_selection'])*100
  RES[i, 'Specificity_selection'] <- mean(lignes[, 'Specificity_selection'])*100
  RES[i, 'Sensitivity_selection'] <- mean(lignes[, 'Sensitivity_selection'])*100
  RES[i, 'MCC_selection'] <- mean(lignes[, 'MCC_selection'])
  
  # Classification indicators (test data)
  RES[i, 'Accuracy_test'] <- mean(lignes[, 'Accuracy_test'])*100
  RES[i, 'Specificity_test'] <- mean(lignes[, 'Specificity_test'])*100
  RES[i, 'Sensitivity_test'] <- mean(lignes[, 'Sensitivity_test'])*100
  RES[i, 'MCC_test'] <- mean(lignes[, 'MCC_test'])
  RES[i, 'TPR'] <- mean(lignes[, 'TPR'])
  RES[i, 'TNR'] <- mean(lignes[, 'TNR'])
  RES[i, 'MAE'] <- mean(lignes[, 'MAE'])
  
  # Stability
  list_of_selected_var <- strsplit(lignes$VarSelected, ";")
  RES[i, 'Stability'] <- stabilityFunction(list_of_selected_var)
}

## Computation of the PROMETHEE score --------------------------------------------------------------------

# Define the benefit and cost criterion
benefit_criterion <- c('AUC', 'YoudenSensitivity', 'YoudenSpecificity', 'NbVarInf','Stability')
cost_criterion <- c('MAE', 'Runtime', 'TrainingTime', 'TestTime', 'NbVarSelected', 'NbVarNonInf')

# Define the criterions to use for the TOPSIS computation
criterion <- c("AUC", "YoudenSensitivity", "YoudenSpecificity", 'NbVarInf', 'NbVarNonInf', 'Stability', 'Runtime')
Rankings <- c(AUC = 1, YoudenSensitivity = 2, YoudenSpecificity = 3,  NbVarInf = 4, NbVarNonInf = 5, Stability = 6, Runtime = 7)


# Compute the matrix xij the value for the ith criterion for jth method
xij <- t(RES[, criterion])
colnames(xij) <- RES$Approach

# Weight vector
wi <- sapply(Rankings, FUN = function(rank){
  n <- length(Rankings)
  wi <- (n-rank+1)/sum(n - Rankings + 1)
  return(wi)
})

# Topsis ranks
ranks <- topsisScore(xij, benefit_criterion = benefit_criterion, cost_criterion = cost_criterion, wi = wi)

# Table 3
print(ranks)


# Code for Figure 1
plot_01 <- ggplot(data = ranks)+geom_col(aes(x = row.names(ranks), y = Score), fill = 'darkgrey')+
  scale_x_discrete(limits = rev(row.names(ranks)))+
  scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.2))+
  geom_text(aes(x = row.names(ranks), y = Score + 0.05, label = round(Score, 2)), size = 13)+
  labs(x = "", y = '')+
  grids(axis = 'xy', color = "lightgrey", size = 1.2, linetype = "dashed")+
  theme(axis.text = element_text(size = 37), 
        axis.title = element_text(size = 35), 
        panel.background = element_rect(fill = 'white', colour = 'black'))+
  coord_flip()

# Save figure 1
setEPS()
postscript('Figures/Figure_01.eps', horizontal = FALSE, width = 18, height = 10)
plot_01
dev.off()

