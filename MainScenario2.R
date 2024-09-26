####################################################################################################################################
#                                                           SCENARIO 2                                                             #  
####################################################################################################################################

## Environment ---------------------------------------------------------------------------------------------------------------------
if(!require(doRNG)){
  # Package used for reproductibility of parallel looping (define a seed for parallel loop)
  devtools::install_version('doRNG', version = '1.8.6')
  library(doRNG)
}

if(!require(ggplot2)){
  # Package used for graphics
  devtools::install_version('ggplot2', version = '3.4.3')
  library(ggplot2)
}

if(!require(ggpubr)){
  # Package used for grids in ggplot graphs
  devtools::install_version('ggpubr', version = '0.6.0')
  library(ggpubr)
}


## Sources -------------------------------------------------------------------------------------------------------------------------
source('functionGammaMetric.R')
source('functionGeneration.R')
source('functions.R')

## Parameters of the generation ----------------------------------------------------------------------------------------------------
N <- 100                                           # Number of observations
n_var_inf <- 3                                     # Number of informative features
n_noise <- 197                    # Number of non-informative features
separation <- c('faible', 'forte')                 # Separation between classes
repartition <- c('equilibre', 'desequilibre')      # Repartition between classes
R <- 50                                            # Number of repetition
set.seed(123)                                      # Seed of the simulation

# Feature selection methods
Approach <- c('FULL', 'CFS', 'CHI2', 'GAMMA_BACK', 'GAMMA_BF', 'GAMMA_FORW', 'LASSO', 'RFI', 'STEP', 'SU', 'SVM-RFE')


# EXECUTE TO RUN THE SIMULATIONS OR GO TO UPLOADING THE RESULTS
###########################################################################################################################################################################
###########################################################################################################################################################################
## Parameters of the parallel computations -----------------------------------------------------------------------------------------

# Number of clusters
n_clusters <- parallel::detectCores()-1
cl <- parallel::makeCluster(n_clusters)
doSNOW::registerDoSNOW(cl)

# Time vector
time <- expand.grid(1:R, repartition, separation, 0)
colnames(time) <- c('Repetition', 'Repartition', 'Separation', 'Time')

# Initialisations
iteration_count <<- 0
i <- 1
res <- NULL

# Progress bar
total_iteration <- R*length(Approach)*length(repartition)*length(separation)
pb <- txtProgressBar(min = 0, max = total_iteration, style = 3)
progress <- function(n) setTxtProgressBar(pb, n + iteration_count)
opts <- list(progress = progress)

## BEGINING OF THE SIMULATION -----------------------------------------------------------------------------------------------------
for(sep in separation){
  for(repa in repartition){
    # Time starter
    t0 <- Sys.time()
    
    # Definition of the beta vector
    if(sep == 'faible'){
      # In the case of weak separation, balanced and unbalanced classes
      beta0 <- ifelse(repa == 'equilibre', 0.5, -2.75)
      beta <- c(0.6, -2.5, -1, rep(0, n_noise))
    }else{
      if(repa == 'equilibre'){
        # In the case of strong separation and balanced classes
        beta0 <- 0
        beta <- c(3.6, -4, -1, rep(0, n_noise))
      }else{
        # In the case of strong separation and unbalanced classes
        beta0 <- -2.65
        beta <- c(3.6, -2.2, -1, rep(0, n_noise))
      }
    }
    
    for(r in 1:R){
      # Generation of the training data
      dat_train <- generation(n = N, n_var_inf = n_var_inf, n_noise = n_noise, beta = c(beta0, beta))
      X_train <- as.matrix(dat_train[, -1])
      Y_train <- as.matrix(dat_train[, 1])
      
      # Generation of the validation data
      dat_test <- generation(n = N, n_var_inf = n_var_inf, n_noise = n_noise, beta = c(beta0, beta))
      X_test <- as.matrix(dat_test[, -1])
      Y_test <- as.matrix(dat_test[, 1])
      
      # Perform feature selection, modelisation and prediction in parallel for each feature selection methods
      resIteration <- foreach(fs = Approach, .combine = 'rbind', .export = ls()[!ls() %in% c('performSimulationIteration', 'X_train', 'Y_train', 'X_test', 'Y_test', 'beta')], .options.snow = opts, .verbose = FALSE) %dorng% {
        localRes <- performSimulationIteration(X_train = X_train, Y_train = Y_train, X_test = X_test, Y_test = Y_test, beta = beta, fs = fs)
        row.names(localRes) <- NULL
        return(localRes)
      }
      
      # Filling the common column 
      resIteration$Separation <- sep
      resIteration$Repartition <- repa
      resIteration$NbVar <- n_var_inf + n_noise
      
      ## Updates 
      
      # Counter
      iteration_count <<- i*length(Approach)
      
      # Append the results data.frames
      res <- rbind(res, resIteration)
      t1 <- Sys.time()
      
      # Time
      time[i, 'Time'] <- as.numeric(difftime(t1, t0, units = 'secs'))
      i <- i+1
    }
    
    # Intermediate saves
    file_path <- paste0('Scenario2/res_', R, '_repetitions_', sep, '_', repa, '.txt')
    write.table(res, file_path)
    res <- NULL
  }
}
close(pb)
parallel::stopCluster(cl)

###########################################################################################################################################################################
###########################################################################################################################################################################

## Uploading the results for graphics and tables --------------------------------------------------------------------------------------------------------------------------
# To upload the results of scenario 2 we have to gather all the files corresponding to the same situation. Hence we make a distinction with the number of features first.

# All the files in scenario 2
files <- dir('Scenario2/')

# Get files related to the same number of features 
files_200 <- paste('Scenario2/', files[grep(paste0(c('res_',R,'_repetitions_'), collapse = ''), files)], sep = '')

# Aggregating the files in one dataframe
res_200 <- data.frame(do.call('rbind', lapply(files_200, read.table)))

## Aggregation ------------------------------------------------------------------------------------------------------------------------------------------------------------

# Define the benefit and cost criterion
benefit_criterion <- c('AUC', 'MCC_test', 'NbVarInf', 'YoudenSensitivity', 'YoudenSpecificity', 'Stability')
cost_criterion <- c('MAE', 'Runtime', 'TrainingTime', 'TestTime', 'NbVarNonInf', 'NbVarSelected')

RES_200 <- data.frame(Approach = rep(Approach, each = length(separation)*length(repartition)), 
                      Separation = rep(separation, each = length(repartition)), 
                      Repartition = rep(repartition, each = 1))

## Filling the rows with mean results
for(i in 1:nrow(RES_200)){
  # For 200 features -------------------------------------
  lignes <- res_200[which(res_200$Approach == RES_200[i, 'Approach'] &
                            res_200$Separation == RES_200[i, 'Separation'] &
                            res_200$Repartition == RES_200[i, 'Repartition']),]

  # Number of features
  RES_200[i, 'NbVarSelected'] <- mean(lignes[, 'NbVarSelected'])
  RES_200[i, 'NbVarInf'] <- mean(lignes[, 'NbVarInf'])
  RES_200[i, 'NbVarNonInf'] <- RES_200[i, 'NbVarSelected'] - RES_200[i, 'NbVarInf']

  # Time indicators
  RES_200[i, 'Runtime'] <- mean(lignes[, 'Runtime'])
  RES_200[i, 'TrainingTime'] <- mean(lignes[, 'TrainingTime'])
  RES_200[i, 'TestTime'] <- mean(lignes[, 'TestTime'])

  # Global indicator
  RES_200[i, 'AUC'] <- mean(lignes[, 'AUC'])*100
  RES_200[i, 'Youden'] <- mean(lignes[, 'Youden'])*100
  RES_200[i, 'YoudenSensitivity'] <- mean(lignes[, 'YoudenSe'])*100
  RES_200[i, 'YoudenSpecificity'] <- mean(lignes[, 'YoudenSpe'])*100

  # Selection indicators
  RES_200[i, 'Accuracy_selection'] <- mean(lignes[, 'Accuracy_selection'])*100
  RES_200[i, 'Specificity_selection'] <- mean(lignes[, 'Specificity_selection'])*100
  RES_200[i, 'Sensitivity_selection'] <- mean(lignes[, 'Sensitivity_selection'])*100
  RES_200[i, 'MCC_selection'] <- mean(lignes[, 'MCC_selection'])

  # Classification indicators (test data)
  RES_200[i, 'Accuracy_test'] <- mean(lignes[, 'Accuracy_test'])*100
  RES_200[i, 'Specificity_test'] <- mean(lignes[, 'Specificity_test'])*100
  RES_200[i, 'Sensitivity_test'] <- mean(lignes[, 'Sensitivity_test'])*100
  RES_200[i, 'MCC_test'] <- mean(lignes[, 'MCC_test'])

  RES_200[i, 'TPR'] <- mean(lignes[, 'TPR'])
  RES_200[i, 'TNR'] <- mean(lignes[, 'TNR'])
  RES_200[i, 'MAE'] <- mean(lignes[, 'MAE'])

  # Stability
  list_of_selected_var <- strsplit(lignes$VarSelected, ';')
  RES_200[i, 'Stability'] <- stability_function(list_of_selected_var)
}

# 
criterion <- c('AUC', 'YoudenSensitivity', 'YoudenSpecificity', 'NbVarInf', 'NbVarNonInf', 'Stability', 'Runtime')
Rankings <- c(AUC = 1, YoudenSensitivity = 2, YoudenSpecificity = 3, NbVarInf = 4, NbVarNonInf = 5, Stability = 6, Runtime = 7)

# Topsis score
RES_200[which((RES_200$Separation == 'faible') & (RES_200$Repartition == 'desequilibre')), c('TopsisScore', 'TopsisRank')] <- topsisRanking(RES_200, separation = "faible", repartition = 'desequilibre', criterion = criterion, rankings = Rankings)[c('Score', 'Rank')]
RES_200[which((RES_200$Separation == 'faible') & (RES_200$Repartition == 'equilibre')), c('TopsisScore', 'TopsisRank')] <- topsisRanking(RES_200, separation = 'faible', repartition = 'equilibre', criterion = criterion, rankings = Rankings)[c('Score', 'Rank')]
RES_200[which((RES_200$Separation == 'forte') & (RES_200$Repartition == 'desequilibre')), c('TopsisScore', 'TopsisRank')] <- topsisRanking(RES_200, separation = 'forte', repartition = 'desequilibre', criterion = criterion, rankings = Rankings)[c('Score', 'Rank')]
RES_200[which((RES_200$Separation == 'forte') & (RES_200$Repartition == 'equilibre')), c('TopsisScore', 'TopsisRank')] <- topsisRanking(RES_200, separation = 'forte', repartition = 'equilibre', criterion = criterion, rankings = Rankings)[c('Score', 'Rank')]

# Topsis rank
RES_200[which((RES_200$Separation == 'faible') & RES_200$Repartition == 'desequilibre'), c('Approach', criterion, 'TopsisScore', 'TopsisRank')]
RES_200[which((RES_200$Separation == 'faible') & RES_200$Repartition == 'equilibre'), c('Approach', criterion, 'TopsisScore', 'TopsisRank')]
RES_200[which((RES_200$Separation == 'forte') & RES_200$Repartition == 'desequilibre'), c('Approach', criterion, 'TopsisScore', 'TopsisRank')]
RES_200[which((RES_200$Separation == 'forte') & RES_200$Repartition == 'equilibre'), c('Approach', criterion, 'TopsisScore', 'TopsisRank')]

# Print for upper part of Table 4
print(cbind(RES_200[which((RES_200$Separation == 'faible') & RES_200$Repartition == 'desequilibre'), c('Approach', criterion, 'TopsisRank', "TopsisScore")],
            RES_200[which((RES_200$Separation == 'faible') & RES_200$Repartition == 'equilibre'), c('Approach', criterion, 'TopsisRank', 'TopsisScore')]))

# Print for lower part of Table 4
print(cbind(RES_200[which((RES_200$Separation == 'forte') & RES_200$Repartition == 'desequilibre'), c('Approach', criterion, 'TopsisRank')],
            RES_200[which((RES_200$Separation == 'forte') & RES_200$Repartition == 'equilibre'), c('Approach', criterion, 'TopsisRank')]))

# Code Figure 2
separation.labs <- c('Weak', 'Strong')
names(separation.labs) <- c('faible', 'forte')

repartition.labs <-  c('Unbalanced', 'Balanced')
names(repartition.labs) <- c('desequilibre', 'equilibre')

plot_02 <- ggplot(data = RES_200)+geom_col(aes(x = Approach, y = TopsisScore), fill = 'darkgrey')+
  scale_x_discrete(limits = rev(Approach))+
  scale_y_continuous(limits = c(0, 1))+
  geom_text(aes(x = Approach, y = TopsisScore + 0.07, label = round(TopsisScore, 2)), size = 11)+
  labs(x = '', y = '')+
  grids(axis = 'xy', color = 'lightgrey', size = 1.2, linetype = 'dashed')+
  theme(axis.text = element_text(size = 30),
        axis.title = element_text(size = 25),
        strip.text = element_text(size = 25),
        panel.background = element_rect(fill = 'white', colour = 'black'))+
  facet_grid(Separation~Repartition, labeller = labeller(Separation = separation.labs, Repartition = repartition.labs))+
  coord_flip()

setEPS()
postscript('Figures/Figure_02.eps', width = 18, height = 10, horizontal = FALSE)
plot_02
dev.off()

