####################################################################################################################################
#                                                           SCENARIO 3                                                             #  
####################################################################################################################################

## Environment ---------------------------------------------------------------------------------------------------------------------
if(!require(doRNG)){
  # Package for reproductibility of parallel execution
  devtools::install_version('doRNG', version = '1.8.6')
  library(doRNG)
}

if(!require(ggplot2)){
  # Package use for plots
  devtools::install_version('ggplot2', version = '3.4.3')
  library(ggplot2)
}

if(!require(ggpubr)){
  # Package used for grids in ggplots
  devtools::install_version('ggpubr', version = '0.6.0')
  library(ggpubr)
}

## Sources ------------------------------------------------------------------------------------------------------------------------
source('functionGammaMetric.R')
source('functionGeneration.R')
source('functions.R')

## Parameters of the generation ---------------------------------------------------------------------------------------------------
N <- 2000                                                                   # Number of observations
s_g <- 10                                                                   # Size of each groups
n_g <- 10                                                                   # Number of groups
non_zero_coeff <- 1                                                         # Informative features are the first feature of each group 
non_zero_group <- c(1, 2, 3, 4, 5, 10)                                      # Groups with non-zero coefficients
value <- 1.5                                                                # Beta coefficients value (when different from 0)
n_ind <- 1                                                                  # Number of independent group
simulation <- 1:6                                                           # Cases for this scenario 
R <- 50                                                                    # Number of repetitions
set.seed(123)                                                               # Seed
beta_vec <- produceBeta(s_g, n_g, non_zero_coeff, non_zero_group, value)    # Beta vector
beta <- beta_vec[-1]
beta0 <- beta_vec[1]

# Feature selection methods 
Approach <- c('FULL', 'CFS', 'CHI2', 'GAMMA_BACK', 'GAMMA_BF', 'GAMMA_FORW', 'LASSO', 'RFI', 'STEP', 'SU', 'SVM-RFE')


# EXECUTE TO RUN THE SIMULATIONS OR GO TO UPLOADING THE RESULTS
###########################################################################################################################################################################
###########################################################################################################################################################################

## Parameters of the parallel computations ---------------------------------------------------------------------------------------

# Number of clusters
n_clusters <- parallel::detectCores() - 1                             # Number of clusters available
cl <- parallel::makeCluster(n_clusters)                               # Initializing n_clusters
doSNOW::registerDoSNOW(cl)                                            # Parallel work

# Time vector
time <- expand.grid(Repetition = 1:R, 
                    Simulation = simulation, 
                    Time = 0)

# Initialisations
iteration_count <<- 0                                                 # Counter
i <- 1                                                                # Counter for progress bar
res <- NULL                                                           # Matrix of results

# Progress bar
total_iteration <- R*length(simulation)*length(Approach)              # Total number of iterations
pb <- txtProgressBar(min = 0, max = total_iteration, style = 3)       # A progress bar
progress <- function(n) setTxtProgressBar(pb, n + iteration_count)    # Function to update the progress bar
opts <- list(progress = progress)                                     # List of options for the progress bar

## BEGINING OF THE SIMULATION -----------------------------------------------------------------------------------------------------

for(simu in simulation){
  # We defined a set of parameters for generation of data with regard to the level and type of dependence we considered
  
  if(simu == 1){
    # First case with a constant correlation level alpha_max = 0.9
    alpha_max <- 0.9
    constant_cor <- TRUE
    
  }else if(simu == 2){
    # Second case with a constant correlation level alpha_max = 0.6
    alpha_max <- 0.6
    constant_cor <- TRUE
    
  }else if(simu == 3){
    # Third case with a constant correlation level alpha_max = 0.3
    alpha_max <- 0.3
    constant_cor <- TRUE
    
  }else if(simu == 4){
    # Fourth case with a non-constant correlation level:
    # Maximum level alpha_max = 0.9 and minimum level c = 0.35
    alpha_max <- 0.9
    c <- 0.35
    constant_cor <- FALSE
    
  }else if(simu == 5){
    # Fith case with a non-constant correlation level:
    alpha_max <- 0.6
    c <- 0.25
    constant_cor <- FALSE
  }else if(simu == 6){
    # Sixth case with a non-constant correlation level:
    alpha_max <- 0.3
    c <- 0.1
    constant_cor <- FALSE
  }
  
  ## Beginning of the repetitions
  for(r in 1:R){
    # Time counter
    t0 <- Sys.time()
    
    # Generation of the data (train data set)
    dat_train <- genData(n = N, s_g = s_g, n_g = n_g, beta0 = beta0, beta = beta, alpha_max = alpha_max, c = c, n_ind = n_ind, constant_cor = constant_cor)
    X_train <- dat_train[, -1]
    Y_train <- dat_train[, 1]
    
    # Generation of the data (test data set)
    dat_test <- genData(n = N, s_g = s_g, n_g = n_g, beta0 = beta0, beta = beta, alpha_max = alpha_max, c = c, n_ind = n_ind, constant_cor = constant_cor)
    X_test <- dat_test[, -1]
    Y_test <- dat_test[, 1]
    
    # Perform feature selection, modelisation and prediction in parallel for each feature selection method
    resIteration <- foreach(fs = Approach, .combine = 'rbind', .export = ls()[!ls() %in% c('performSimulationIteration', 'X_train', 'X_test', 'Y_train', 'Y_test', 'beta', 'c')], .options.snow = opts, .verbose = FALSE) %dorng%{
      localRes <- performSimulationIteration(X_train, Y_train, X_test, Y_test, beta, fs)
      row.names(localRes) <- NULL
      return(localRes)
    }
    
    # Filling the common columns
    resIteration$Simulation <- simu
    
    # Updates
    t1 <- Sys.time()
    iteration_count <- i*length(Approach)
    res <- rbind(res, resIteration)
    time[i, 'Time'] <- as.numeric(difftime(t1, t0, units = 'secs'))
    i <- i+1
  }
  
  # Intermediate save
  file_path <- paste0('Scenario3/res_simulation_', simu, '_', R, '_repetitions.txt')
  write.table(res, file_path)
  res <- NULL
}
close(pb)
parallel::stopCluster(cl)

###########################################################################################################################################################################
###########################################################################################################################################################################

## Uploading the results for graphics and tables --------------------------------------------------------------------------------------------------------------------------

# All the results files
files <- dir('Scenario3/')

# Aggregations of the results files
patterns <- paste0('res_simulation_', 1:6)
res <- data.frame(do.call('rbind', lapply(
  paste0('Scenario3/', files[grep(paste(patterns, collapse = '|'), files)]), read.table)
))

## Aggregation ------------------------------------------------------------------------------------------------------------------------------------------------------------

RES <- data.frame(expand.grid(
  Approach = Approach, 
  Simulation = simulation
))

for(i in 1:nrow(RES)){
  # Line
  lignes <- res[which(res$Approach == RES[i, 'Approach'] & res$Simulation == RES[i, 'Simulation']),]
  
  # Number of features
  RES[i, 'NbVarSelected'] <- mean(lignes[,'NbVarSelected'])
  RES[i, 'NbVarInf'] <- mean(lignes[, 'NbVarInf'])
  RES[i, 'NbVarNonInf'] <- RES[i, 'NbVarSelected'] - RES[i, 'NbVarInf']
  
  # Time indicators
  RES[i, 'Runtime'] <- mean(lignes[, 'Runtime'])
  RES[i, 'TrainingTime'] <- mean(lignes[, 'TrainingTime'])
  RES[i, 'TestTime'] <- mean(lignes[, 'TestTime'])
  
  # Global indicators
  RES[i, 'AUC'] <- mean(lignes[, 'AUC'])*100
  RES[i, 'Youden'] <- mean(lignes[, 'Youden'])*100
  RES[i, 'YoudenSensitivity'] <- mean(lignes[, 'YoudenSensitivity'])*100
  RES[i, 'YoudenSpecificity'] <- mean(lignes[, 'YoudenSpecificity'])*100
  
  # Selection indicators
  RES[i, 'Accuracy_selection'] <- mean(lignes[,'Accuracy_selection'])*100
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
  list_of_selected_var <- strsplit(lignes$VarSelected, ';')
  RES[i, 'Stability'] <- stabilityFunction(list_of_selected_var)
}

# Define the benefit and cost criterion
benefit_criterion <- c('AUC', 'YoudenSensitivity', 'YoudenSpecificity', 'MCC_test', 'NbVarInf', 'TPR', 'TNR', 'Stability')
cost_criterion <- c('MAE', 'Runtime', 'TrainingTime', 'TestTime', 'NbVarNonInf', 'NbVarSelected')

# TOPSIS score
criterion <- c('AUC', 'YoudenSensitivity', 'YoudenSpecificity', 'NbVarInf', 'NbVarNonInf', 'Stability', 'Runtime')
Rankings <- c(AUC = 1, YoudenSensitivity = 2, YoudenSpecificity = 3, NbVarInf = 4, NbVarNonInf = 5, Stability = 6, Runtime = 7)
wi <- (length(Rankings) - Rankings + 1)/sum(length(Rankings) - Rankings + 1)

for(simu in simulation){
  # Creation of xij
  xij <- t(RES[which(RES$Simulation == simu), criterion])
  colnames(xij) <- RES[which(RES$Simulation == simu), 'Approach']
  
  topsis <- topsisScore(xij, benefit_criterion = benefit_criterion, cost_criterion = cost_criterion, wi = wi)
  RES[which(RES$Simulation == simu), c('Score', 'Rank')] <- topsis[,c('Score', 'Rank')]
}


#### Table ----------------------------------------------------------------------------------------------------------------------------------------------------------------

# Code for Table 
cbind(RES[which(RES$Simulation == 1),c('Approach', criterion, 'Rank')], RES[which(RES$Simulation == 4),c('Approach', criterion, 'Rank')])
cbind(RES[which(RES$Simulation == 2),c('Approach', criterion, 'Rank')], RES[which(RES$Simulation == 5),c('Approach', criterion, 'Rank')])
cbind(RES[which(RES$Simulation == 3),c('Approach', criterion, 'Rank')], RES[which(RES$Simulation == 6),c('Approach', criterion, 'Rank')])

# RES 
RES[which(RES$Simulation == 1), c('Correlation', 'Alpha')] <- rep(c('Constant', '0.9'), each = length(which(RES$Simulation == 1)))
RES[which(RES$Simulation == 2), c('Correlation', 'Alpha')] <- rep(c('Constant', '0.6'), each = length(which(RES$Simulation == 2)))
RES[which(RES$Simulation == 3), c('Correlation', 'Alpha')] <- rep(c('Constant', '0.3'), each = length(which(RES$Simulation == 3)))
RES[which(RES$Simulation == 4), c('Correlation', 'Alpha')] <- rep(c('Non-constant', '0.9'), each = length(which(RES$Simulation == 4)))
RES[which(RES$Simulation == 5), c('Correlation', 'Alpha')] <- rep(c('Non-constant', '0.6'), each = length(which(RES$Simulation == 5)))
RES[which(RES$Simulation == 6), c('Correlation', 'Alpha')] <- rep(c('Non-constant', '0.3'), each = length(which(RES$Simulation == 6)))

RES$Correlation <- factor(RES$Correlation)

RES$Alpha <- factor(RES$Alpha, 
                    levels = c('0.9', '0.6', '0.3'), 
                    labels = c(
                      expression(paste(alpha[max], ' = 0.9')),
                      expression(paste(alpha[max], ' = 0.6')),
                      expression(paste(alpha[max], ' = 0.3'))
                    ))

plot_03 <- ggplot(data = RES)+
  geom_col(aes(x = Approach, y = Score), fill = 'darkgrey')+
  geom_text(aes(x = Approach, y = Score+0.05, label = round(Score, 2)), size = 8)+
  scale_x_discrete(limits = rev(Approach))+
  scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.2))+
  labs(x = "", y = "")+
  grids(axis = 'xy', color = 'lightgrey', size = 1.1, linetype = 'dashed')+
  theme(axis.text = element_text(size = 20), 
        strip.text = element_text(size = 25),
        panel.background = element_rect(fill = 'white', colour = 'black'))+
  facet_grid(Alpha~Correlation, labeller = label_parsed)+
  coord_flip()

setEPS()
postscript('Figures/Figure_03.eps', width = 18, height = 10, horizontal = FALSE)
plot_03
dev.off()

