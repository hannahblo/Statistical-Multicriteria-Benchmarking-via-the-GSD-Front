# TODOS:
- @Julian bitte nochmal pmlb durchlesen und erklären was die einzelnen Files machen--> siehe TODOOOOS und die Laufzeit eintragen
- @alle: Achtung, wenn eine der Permutationsläufe im PMLB nicht verwendet wird, dann müssen wir das in der Erklärung (direkt nach Introduction) noch rausschmeißen.

# Statistical Multicriteria Benchmarking via the GSD-Front

## Introduction
This anonymous repository contains R-code and data sets corresponding to the "Statistical Multicriteria Benchmarking via the GSD-Front" article. We apply the introduced tests on two data sets: This are given by the OpenML and the PLMB suite.

The structure of the repository is as follows:
- File _setup_session.R installs all (except for gurobi) needed R-packages.
- File openml_permutation_tests.R main file to compute the permutation test and the main plots based on the OpenML suite
- File pmlb_permutation_tests.R main file to compute the permutation test and the main plots based on the PMLB suite
- Folder R/ contains the needed functions (permutation test, sampling, definition of the constraint matrix, plotting etc) for the analysis.
- Folder pmbl_results/ contains the results as well as the code of the evaluations (robustified) accuracy on the PMLB suite
- Folder openml_permutation_results contains all the results provided by the file openml_permutation_tests.R (Note that they are not automatically saved in this folder, but copied to it.)
- Folder pmlb_permutation_results contains all the results provided by the file pmlb_permutation_tests.R (Note that they are not automatically saved in this folder, but copied to it.)

The code was tested with
- R version 4.2.1
- R version 4.2.2
- R version 4.3.2

on

- Linux Ubuntu 20.04.5
- Windows 10/11

## Notation
To clarify the different notations between the article in the Code.

For OpenML:
- "classif.svm" is renamed to "SVM" in the article
- "classif.multinom" is renamed to "LR" in the article
- "classif.ranger" is renamed to "RF" in the article
- "classif.xgboost" is renamed to "xGBoost" in the article
- "classif.glmnet" is renamed to "GLMNet" in the article
- "classif.kknn" is renamed to "kNN" in the article
- "classif.rpart" is renamend to "CART" in the article

For PMLB:
- "cre" is renamed to "CRE" in the article
- "svmRadial" is renamed to "SVM" in the article
- "J48" is renamed to "CART" in the article
- "ranger" is renamed to "RF" in the article
- "knn" is fenamed to "kNN" in the article
- "glmnet is renamed to "GLMNet" in the article

## Setup
First, please install all necessary R-packages (can be found in _setup_session.R):
- For the computation of the linear programs, we used the R interface of gurobi optimizer, see [here](https://www.gurobi.com/) (accessed: 08.02.2023). This is a commercial
solver that offers a free academic licenses which can be found [here](https://www.gurobi.com/features/academic-named-user-license/) (accessed: 08.02.2023). To install this package, please follow the instructions there. A documentation can be found [here](https://www.gurobi.com/wp-content/plugins/hd_documentations/documentation/9.0/refman.pdf) (page 643ff) (accessed: 08.02.2023).
- Afterwards, please install all dependencies by sourcing the file _setup_session.R and source the files 

Then download the following files and save them in a folder named 'R' (these files can already be found in the folder R/):
- test_two_items.R
- constraints_r1_r2.R
- sample_permutation_test.R
- plotting_permutationtest.R

Afterwards, please decided which data set you are interested OpenML vs PMLB\
\
**OpemML** If you are interested in OpenML, download openml_permutation_tests.R and run it. Note that the 'R/' folder and the openml_permutation_tests.R file must be in the same folder structure. All results produced by this file will be stored in the current that contains openml_permutation_tests.R. For tidiness, we decided to copy the results to the openml_permutation_results/ folder after the entire computation. Therefore, if you are interested in running only parts of the code, please change your working directory (after you have obtained all the necessary files) to the openml_permutation_results/ folder.\
\
**PMLB** If you are interested in PMLB, download the following files and store them in a folder named 'pmlb_results' (these files can already be found in the folder pmlb_evaluations/):
   - main_pmlb_experiments.R TODOOOOOOOOO
   - helper_pmlb_experiments.R TODOOOOOOO
   - make_method_cre.R TODOOOOOOOOOOOOOOO
   - final_results.RDS TODOOOOOOOOOOOOOOOO
     
Either you can now reproduce the evaluations (computing accuracy and a robustified accuracy) by running main_pmlb_experiments.R, or you can immediately go to the next step and use the already computed results stored in final_results.RDS.\
Now that the evaluation is computed and saved, go back to the main folder and save the file pmlb_permutation_tests.R to run the permutation tests. Note that the folder 'R/', 'pmlb_results/' and the file pmlb_permutation_tests.R must be in the same folder structure. All results produced by this file will be stored in the current that contains pmlb_permutation_tests.R. For tidiness, we decided to copy the results after the entire computation to the pmlb_permutation_results/ folder. Thus, if you are interested in running only parts of the code, please change your working directory (after you have obtained all necessary files) to the folder pmlb_permutation_results/.

## Total Computation time and produced objects

The computation times of the main files is as follows
- openml_permutation_tests.R: approximately 2 weeks
- pmlb_permutation_tests.R: approximately 5 days 
- main_pmlb_experiments.R: approximately TODOOOOOOOOOOOOOOOO

The files openml_permutation_tests.R and pmlb_permutation_tests.R produce the following objects. \
\
The **OpenML** test compares the algorithm named "classif.svm" (classifier of interest in the following) against "classif.multinom", "classif.ranger", "classif.xgboost", "classif.glmnet", "classif.kknn", and "classif.rpart". In the following # is a placeholder for "classif.multinom", "classif.ranger", "classif.xgboost", "classif.glmnet", "classif.kknn", and "classif.rpart". * and ° are two placeholders for "classif.svm", "classif.multinom", "classif.ranger", "classif.xgboost", "classif.glmnet", "classif.kknn", and "classif.rpart".\
\
The **PMLB** test compares algorithm named "cre" (classifier of interest in the following against "svmRadial", "J48", "ranger", "knn", and "glmnet". In the following # is a placeholder for "svmRadial", "J48", "ranger", "knn", or "glmnet". * and ° are two placeholders for "cre" "svmRadial", "J48", "ranger", "knn", and "glmnet".
- #_result.rds stores the result of the empirical as well as the permutated test statistic where we compute the expectation of # minus the classifier of interest.
- #_computation_time.rds stores the computation time to obtain #_result.rds
- #_dat_set.rds stores the used sorted data set (here group_a corresponds to the classifier of interest and # to group_b)
- *_°_result_all.rds stores the result of the empirical as well as the permutated test statistic where we compute the expectation of * minus °
- *_°_computation_time.rds stores the time to compute *_°_result_all.rds
- proportion_below_df.rds stores the number of permutation results that are below the observed one
- dat_openml_filter.rds stores the downloaded OpenMl data (download was done 26.01.2024)
- dat_final.rds stores the rescaled OpenML data (saved 26.01.2024)
-  fig_2.pdf, fig_3.pdf, fig_4.pdf are the Figures 2 to 4 in the paper visualizing the test results for OpenML
-  fig_5.pdf, fig_6.pdf, fig_7.pdf are the Figures 5 to 7 in the paper visualizing the test results for PMLB 

## References:
- Blocher, H., Schollmeyer, G., Jansen, C., and Nalenz, M. Depth functions for partial orders with a descriptive analysis of machine learning algorithms. In Miranda, E., Montes, I., Quaeghebeur, E., and Vantaggi, B. (eds.), Proceedings of the Thirteenth International Symposium on Imprecise Probability: Theories and Applications, volume 215 of Proceedings of Machine Learning Research, pp. 59–71. PMLR, 11–14 Jul 2023
- Jansen, C., Schollmeyer, G., Blocher, H., Rodemann, J., and Augustin, T. Robust statistical comparison of random variables with locally varying scale of measurement. In Evans, R. J. and Shpitser, I. (eds.), Proceedings of the Thirty-Ninth Conference on Uncertainty in Artificial Intelligence, volume 216 of Proceedings of Machine Learning Research, pp. 941–952. PMLR, 31 Jul–04 Aug 2023
- Olson, R. S., La Cava, W., Orzechowski, P., Urbanowicz, R. J., and Moore, J. H. Pmlb: a large benchmark suite for machine learning evaluation and comparison. BioData Mining, 10:36, 2017
- Van Rijn, J. N., Bischl, B., Torgo, L., Gao, B., Umaashankar, V., Fischer, S., Winter, P., Wiswedel, B., Berthold, M. R., and Vanschoren, J. Openml: A collaborative science platform. In Machine Learning and Knowledge Discovery in Databases: European Conference, ECML PKDD 2013, Prague, Czech Republic, September 23-27, 2013, Proceedings, Part III 13, pp. 645–649. Springer, 2013


