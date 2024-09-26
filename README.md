# Multivariate-filter-methods-for-feature-selection-with-the-gamma-metric

Repository with the R scripts used for the simulation study of the feature selection with the γ-metric. Each scenario can be run independently of any other scenario. If one wish to execute of the Scenario, first make sure that the scripts 'function.R', 'functionGammaMetric.R', and 'functionGeneration.R' 
are installed in the same directory and set it as the working directory of your R environment. Second for Scenario 2 and Scenario 3, intermediate saves are executed inside the for loops. Make sure you that running the script won't erase previous results. 
Third, running all the iterations of the simulation might take some times (depending on the configuration of your device), you might consider lowering the number of repetition (named R in all scripts) if you won't to execute the code. 

# Folder Scenario 1
Folder with all the results of Scenario 1 of the simulation study (1 file). 

# Folder Scenario 2
Folder with all the results of Scenario 2 of the simulation study (4 files). One file for each combination separation of the classes and balance of the classes. 

# Folder Scenario 3
Folder with all the results of Scenario 3 of the simulation study (6 files). One for each combination levels of correlation and type of correlation (constant or non-constant).

# File functionGammaMetric.R
R scripts with all the functions needed to compute the γ-metric value of a set of features. The main function to call to compute the γ-metric value of set of features is gammaMetric().

# File functionGeneration.R
This file contains the functions to generate data for all the scenarios. With the following description:
- generation(): to generate the data of scenario 1 and 2.
- genData(): to generate the data of scenario 3.
Two other functions, produceRho() and produceBeta() are used to generate the variance-covariance matrix (Σ) and β vector of scenario 3.

# File functions.R
This file contains functions that are used in the different scenarios:
- approachFunctions: it is actually a list of functions with each one being a feature selection method.
- performSimulationIteration(): This function is called inside the parallel loop of the simulation in each scneario. It is used to execute the feature selection (on train), apply the logistic regression model (on train), compute the prediction of the model (on validation) and the performance indicators.
- jaccardFunction(): This function is used to compute the jaccard index between two sets of features.
- stabilityFunction(): This function is used to compute the stability index, which is the pairwise jaccard index averaged over all sets of selected features.
- topsisScore(): Compute the TOPSIS score.
- topsisRanking(): Give the TOPSIS rank.

# File appendixFigures.R
This script contains the R code necessary to generate the figures and tables of the appendix and the supporting information.

# File sessionInfo.R
This file is the output of the R command sessionInfo(). This file list information about the R version, the platform used and the operating system. Also some details about the packages (and their versions) loaded and attached in the environment at the time of the computation of the simulation study.
