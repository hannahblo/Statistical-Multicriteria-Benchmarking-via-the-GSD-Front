#############
## PMLB datasets
# This file (and helpers) corresponds to the pmlb analysis (see section 5.1 in the paper)
# We select all datasets from PMLB benchmarking suite (https://epistasislab.github.io/pmlb/)
# for binary classification tasks with 40 to 1000 observations and less than
# 100 features (see lines 36-49 in file "main_pmlb_experiments.R") (62 in total).
# The latter file runs 10-fold cross-validation on each considered classifier on 
# each of the 62 datasets from the PMLB benchmark suite. It returns i) classical 
# accuracies, ii) accuracies with perturbed x and iii) accuracies with perturbed y.
# the resulting data frames has i) to iii) as columns and is loaded in line 56
###########

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# TODO
# TODO
# TODO
# TODO @Hannah in den files noch erklären wie umgang mit na in eps_0 also bei
# beschreibung

# This code is with adaptation copied from
# 1.
# Hannah Blocher, Georg Schollmeyer, Christoph Jansen, and Malte Nalenz. Depth functions for
# Christoph Jansen, Georg Schollmeyer, Hannah Blocher, Julian Rodemann and Thomas Augustin (2023):
# Robust statistical comparison of random variables with locally varying scale of measurement.
# In: Proceedings of the Thirty-Ninth Conference on Uncertainty in Artificial Intelligence (UAI 2023).
# Proceedings of Machine Learning Research, vol. 216. PMLR.
# and the code provided by
# https://github.com/hannahblo/Robust_GSD_Tests (accessed: 18.01.2024)


################################################################################
# Session Settings
################################################################################
library(purrr)# for detect_index
library(dplyr) # for group_by
library(slam) # simple triplet matrix
library(gurobi) # solving LP
library(readr)  # data preparation
library(forcats) # data preparation
library(ggplot2) # visualization
library(reshape2) # visualization
library(tidyverse) # data wrangling
library(ggridges) # visualization
library(latex2exp) # for gamma (and epsilon) symbols
library(RColorBrewer) # color palettes
library(rcartocolor) # color gradients
library(dplyr)
library(farff)
library(reshape2)
library(parallel)

source("R/constraints_r1_r2.R") # contains the functions compute_constraints...
source("R/sample_permutation_test.R") # permutation test, sample etc
source("R/plotting_permutationtest.R") # plot function
source("R/test_two_items.R") # main function summarizing the computation
################################################################################
# Prepare Data Set: PMLB
################################################################################

data_pmlb <- as.data.frame(readRDS("pmlb_results/final_results.RDS"))
# head(data_pmlb_filter)
data_pmlb$classifier <- rep(NA, dim(data_pmlb)[1])
classifier_all <- c("cre", "svmRadial", "J48", "ranger", "knn", "glmnet")
for (classifier in classifier_all) {
  row_classifier <- which(grepl(classifier, rownames(data_pmlb)))
  data_pmlb$classifier[row_classifier] <- classifier
}
# data_pmlb$classifier

cutting_noisy_x <- quantile(data_pmlb$results_noisy_x_stacked)
data_pmlb$results_noisy_x_stacked <- findInterval(data_pmlb$results_noisy_x_stacked,
                                                         cutting_noisy_x)

cutting_noisy_y <- quantile(data_pmlb$results_noisy_y_stacked)
data_pmlb$results_noisy_y_stacked <- findInterval(data_pmlb$results_noisy_y_stacked,
                                                         cutting_noisy_y)


################################################################################
# Conduct the permutation test and plotting the single comparisons
################################################################################

classifier_of_interest <- "cre"
classifiers_comparison <- list("svmRadial", "J48", "ranger", "knn", "glmnet")

for (classifier in classifiers_comparison) {

  data_pmlb_selected <- data_pmlb

  # Now we select two classifiers and compare them
  # here we choose classif.rpart and classif.svm
  data_pmlb_selected <- data_pmlb_selected[data_pmlb_selected$classifier %in%
                                                 c(classifier_of_interest,
                                                   classifier), ]


  # now, we need to convert the data_final_selected into a dataframe with columns:
  # ordinal_1, ordinal_2, numeric, count_group_a, count_group_b, count_all, ID

  # Step 1: Converting the variables of interest into numeric and order modes
  data_pmlb_selected[["results_clean_stacked"]] <- as.numeric(as.character(
    data_pmlb_selected[["results_clean_stacked"]]))

  data_pmlb_selected[["results_noisy_x_stacked"]] <- as.ordered(as.character(
    data_pmlb_selected[["results_noisy_x_stacked"]]))

  data_pmlb_selected[["results_noisy_y_stacked"]] <- as.ordered(as.character(
    data_pmlb_selected[["results_noisy_y_stacked"]]))

  # data_pmlb_selected <- data_pmlb_selected[, !names(data_pmlb_selected) %in% c("classifier")]


  # Step 2: duplication handling
  data_count <- data_pmlb_selected %>% group_by_all() %>% count()


  data_algorithm_1 <- data_count[which(data_count$classifier == classifier_of_interest),
                                 c(1, 2, 3, 5)]
  data_algorithm_1 <- matrix(as.numeric(as.matrix(data_algorithm_1)), ncol = 4)
  colnames(data_algorithm_1) <-  c("numeric", "ordinal_1", "ordinal_2", "count_group_a")


  data_algorithm_2 <- data_count[which(data_count$classifier == classifier),
                                 c(1, 2, 3, 5)]
  data_algorithm_2 <- matrix(as.numeric(as.matrix(data_algorithm_2)), ncol = 4)
  colnames(data_algorithm_2) <-  c("numeric", "ordinal_1", "ordinal_2", "count_group_b")

  dat_final <- merge(x = data_algorithm_1, y = data_algorithm_2,
                     by = c("ordinal_1", "ordinal_2", "numeric"),
                     all.x = TRUE, all.y = TRUE)
  dat_final[is.na(dat_final)] <- 0
  dat_final$count_all <- dat_final$count_group_a + dat_final$count_group_b
  dat_final$ID <- seq(1:dim(dat_final)[1])


  # View(dat_final)
  # dim(dat_final)
  # min(dat_final$numeric)
  # max(dat_final$numeric)

  index_max <- which(dat_final$numeric == max(dat_final$numeric))[1]
  # dat_final[index_max, ]
  index_min <- which(dat_final$numeric == min(dat_final$numeric))[1]
  # dat_final[index_min, ]


  # Add minimal and maximal at the bottom of the matrix

  # ATTENTION: It is very important for the following analysis that the
  # the input at the second largest row is the minimal value and the largest row
  # represents the maximal value
  dat_final[dim(dat_final)[1] + 1, ] <- c(min(dat_final$ordinal_1),
                                          min(dat_final$ordinal_2),
                                          dat_final[index_min, 3],
                                          0, 0, 0,
                                          max(dat_final$ID) + 1)
  dat_final[dim(dat_final)[1] + 1, ] <- c(max(dat_final$ordinal_1),
                                          max(dat_final$ordinal_2),
                                          dat_final[index_max, 3],
                                          0, 0, 0,
                                          max(dat_final$ID) + 1)


  # Note that we can have now the problem that these two added elements are not
  # allowed to occur already in the data, thus we check this and eventually delete
  # this row

  if (all((dat_final[dim(dat_final)[1] - 1, ] == dat_final[index_min, ])[c(1,2,3)])) {
    dat_final[dim(dat_final)[1] - 1, ] <- dat_final[index_min, ]
    dat_final <- dat_final[-c(index_min), ]
    dat_final$ID <- seq(1:dim(dat_final)[1])

    # We have to update index_max as now the data frame changed
    # note that the added maximum is by default behind [1]
    index_max <- which(dat_final$numeric == max(dat_final$numeric))[1]
  }
  if (all((dat_final[dim(dat_final)[1], ] == dat_final[index_max, ])[c(1,2,3)])) {
    dat_final[dim(dat_final)[1], ] <- dat_final[c(index_max), ]
    dat_final <- dat_final[-c(index_max), ]
    dat_final$ID <- seq(1:dim(dat_final)[1])
  }


  dat_set <- dat_final
  saveRDS(dat_set, paste0(classifier, "_dat_set.rds"))
  # dat_final <- readRDS("dat_final.rds")
  start_time <- Sys.time()
  result_inner <- test_two_items(dat_set, iteration_number = 200,
                                 seed_number = 2983754,
                                 eps_0 = 0,
                                 eps_1 = 0.5,
                                 eps_2 = NA,
                                 eps_3 = NA,
                                 eps_4 = NA)

  # result[[classifier]] <- result_inner
  # plotting_permutationtest(result_inner$permutation_test, result_inner$d_observed,
  # add_name_file = classifier)
  total_time <- Sys.time() - start_time

  saveRDS(result_inner, paste0(classifier, "_result.rds"))
  saveRDS(total_time, paste0(classifier, "_computation_time.rds"))
}


################################################################################
# Result and Plotting
################################################################################
# # plotting the result of the pairwise comparisons
# classifier_of_interest <- "cre"
# classifiers_comparison <- list("svmRadial", "J48", "ranger", "knn", "glmnet")
#
# for (classifier in classifiers_comparison) {
#   result_plot <- readRDS(paste0(classifier, "_result.rds"))
#   plotting_permutationtest(result_plot$permutation_test, result_plot$d_observed,
#                          add_name_file = classifier)
# }
#
#
# # computing the test statistic values
# classifier_of_interest <- "cre"
# classifiers_comparison <- list("svmRadial", "J48", "ranger", "knn", "glmnet")
# all_eps_values <- c("result_eps_0", "result_eps_1")
#
# proportion_above_df <- as.data.frame(matrix(rep(0, 6 * 5), nrow = 5, ncol = 6), row.names = all_eps_values)
# colnames(proportion_above_df) <- unlist(classifiers_comparison)
# for (classifier in classifiers_comparison) {
#   result_classifier <- readRDS(paste0(classifier, "_result.rds"))
#
#   for (eps_value in all_eps_values) {
#     base_value <- result_classifier$d_observed[[eps_value]]
#     proportion_above_df[eps_value, classifier] <-  sum(result_classifier$permutation_test[eps_value, ] < base_value)
#   }
# }
#
# saveRDS(proportion_above_df, "final_result.rds")



# Constructing the graph representing all pairwise comparisons at once
# therfore we have to compute the reverse (switching classifier and classifier_of_interest
# position) observed minimal difference
