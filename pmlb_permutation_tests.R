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
library(farff) # for data preparation

source("R/constraints_r1_r2.R") # contains the functions compute_constraints...
source("R/sample_permutation_test.R") # permutation test, sample etc
source("R/plotting_permutationtest_pmlb.R") # plot function
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
# Conducting the permutation test and plotting the results
################################################################################

classifier_of_interest <- "cre"
classifiers_comparison <- list("svmRadial", "J48", "ranger", "knn", "glmnet")  #

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

  index_max <- which(dat_final$numeric == max(dat_final$numeric))
  # dat_final[index_max, ]
  index_min <- which(dat_final$numeric == min(dat_final$numeric))
  # dat_final[index_min, ]


  # Add minimal and maximal at the bottom of the matrix

  # ATTENTION: It is very important for the following analysis that the
  # the input at the second largest row is the minimal value and the largest row
  # represents the maximal value
  dat_final[dim(dat_final)[1] + 1, ] <- c(min(dat_final$ordinal_1),
                                          min(dat_final$ordinal_2),
                                          dat_final[index_min[1], 3],
                                          0, 0, 0,
                                          max(dat_final$ID) + 1)
  dat_final[dim(dat_final)[1] + 1, ] <- c(max(dat_final$ordinal_1),
                                          max(dat_final$ordinal_2),
                                          dat_final[index_max[1], 3],
                                          0, 0, 0,
                                          max(dat_final$ID) + 1)


  # Note that we can have now the problem that these two added elements are not
  # allowed to occur already in the data, thus we check this and eventually delete
  # this row
  for (i in seq(1, length(index_min))) {
    if (all((dat_final[dim(dat_final)[1] - 1, ] == dat_final[index_min[i], ])[c(1,2,3)])) {
      dat_final[dim(dat_final)[1] - 1, ] <- dat_final[index_min[i], ]
      dat_final <- dat_final[-c(index_min[i]), ]
      dat_final$ID <- seq(1:dim(dat_final)[1])

      # We have to update index_max as now the data frame changed
      # note that we want to compare to the last row and therefore we have
      # to delete this one in index_max
      index_max <- which(dat_final$numeric == max(dat_final$numeric))
      index_max <- index_max[-length(index_max)]
    }
  }
  for (i in seq(1, length(index_max))) {
    if (all((dat_final[dim(dat_final)[1], ] == dat_final[index_max[i], ])[c(1,2,3)])) {
      dat_final[dim(dat_final)[1], ] <- dat_final[c(index_max[i]), ]
      dat_final <- dat_final[-c(index_max[i]), ]
      dat_final$ID <- seq(1:dim(dat_final)[1])
    }
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


# computing the test statistic values
classifier_of_interest <- "cre"
classifiers_comparison <- list("svmRadial", "J48", "ranger", "knn", "glmnet")
all_eps_values <- c("result_eps_0", "result_eps_1")

proportion_below_df <- as.data.frame(matrix(rep(0, 5 * 2), nrow = 2, ncol = 5), row.names = all_eps_values)
colnames(proportion_below_df) <- unlist(classifiers_comparison)
for (classifier in classifiers_comparison) {
  result_classifier <- readRDS(paste0(classifier, "_result.rds"))

  for (eps_value in all_eps_values) {
    base_value <- result_classifier$d_observed[[eps_value]]
    proportion_below_df[eps_value, classifier] <-  sum(result_classifier$permutation_test[eps_value, ] < base_value)
  }
}

saveRDS(proportion_below_df, "proportion_below_df.rds")


################################################################################
# Result and Plotting
################################################################################
# plotting the test results (of the pairwise comparisons) as in figure 2, 3, and 4 (appendix)
# in the paper
classifiers_comparison <- list("svmRadial", "J48", "ranger", "knn", "glmnet")

results_plots = list()
for (classifier in classifiers_comparison) {
  # if(classifier == "classif.xgboost")
  #   debugonce(plotting_permutationtest)
  result_plot <- readRDS(paste0(classifier, "_result.rds"))
  results_plots[[classifier]] = result_plot
}

plotting_permutationtest_pmlb(results_plots)




# Constructing the graph representing all pairwise comparisons at once
# therfore we have to compute the reverse (switching classifier and classifier_of_interest
# position) observed minimal difference
classifiers_all <- list("cre", "svmRadial", "J48", "ranger", "knn", "glmnet")

for (k in seq(1, length(classifiers_all))) {
  for (m in seq(1, length(classifiers_all))[-k]) {

    classifier <- classifiers_all[m]
    classifier_of_interest <- classifiers_all[k]


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

    index_max <- which(dat_final$numeric == max(dat_final$numeric))
    # dat_final[index_max, ]
    index_min <- which(dat_final$numeric == min(dat_final$numeric))
    # dat_final[index_min, ]


    # Add minimal and maximal at the bottom of the matrix

    # ATTENTION: It is very important for the following analysis that the
    # the input at the second largest row is the minimal value and the largest row
    # represents the maximal value
    dat_final[dim(dat_final)[1] + 1, ] <- c(min(dat_final$ordinal_1),
                                            min(dat_final$ordinal_2),
                                            dat_final[index_min[1], 3],
                                            0, 0, 0,
                                            max(dat_final$ID) + 1)
    dat_final[dim(dat_final)[1] + 1, ] <- c(max(dat_final$ordinal_1),
                                            max(dat_final$ordinal_2),
                                            dat_final[index_max[1], 3],
                                            0, 0, 0,
                                            max(dat_final$ID) + 1)


    # Note that we can have now the problem that these two added elements are not
    # allowed to occur already in the data, thus we check this and eventually delete
    # this row
    for (i in seq(1, length(index_min))) {
      if (all((dat_final[dim(dat_final)[1] - 1, ] == dat_final[index_min[i], ])[c(1,2,3)])) {
        dat_final[dim(dat_final)[1] - 1, ] <- dat_final[index_min[i], ]
        dat_final <- dat_final[-c(index_min[i]), ]
        dat_final$ID <- seq(1:dim(dat_final)[1])

        # We have to update index_max as now the data frame changed
        # note that we want to compare to the last row and therefore we have
        # to delete this one in index_max
        index_max <- which(dat_final$numeric == max(dat_final$numeric))
        index_max <- index_max[-length(index_max)]
      }
    }
    for (i in seq(1, length(index_max))) {
      if (all((dat_final[dim(dat_final)[1], ] == dat_final[index_max[i], ])[c(1,2,3)])) {
        dat_final[dim(dat_final)[1], ] <- dat_final[c(index_max[i]), ]
        dat_final <- dat_final[-c(index_max[i]), ]
        dat_final$ID <- seq(1:dim(dat_final)[1])
      }
    }


    dat_set <- dat_final
    start_time <- Sys.time()
    result_inner <- test_two_items(dat_set, iteration_number = 1,
                                   eps_0 = 0,
                                   eps_1 = NA,
                                   eps_2 = NA,
                                   eps_3 = NA,
                                   eps_4 = NA)

    total_time <- Sys.time() - start_time

    saveRDS(result_inner, paste0(classifier, "-", classifier_of_interest, "_result_all.rds"))
    saveRDS(total_time, paste0(classifier, "-", classifier_of_interest, "_computation_time_all.rds"))
  }
}




df_eps_0 <- as.data.frame(matrix(rep(NA, length(classifiers_all) * length(classifiers_all)),
                                 nrow = length(classifiers_all)))
colnames(df_eps_0) <- rownames(df_eps_0) <- classifiers_all

for (k in seq(1, length(classifiers_all))) {
  for (m in seq(1, length(classifiers_all))[-k]) {
    classifier <- classifiers_all[k]
    classifier_of_interest <- classifiers_all[m]

    result_inner <- readRDS(paste0(classifier, "-", classifier_of_interest, "_result_all.rds"))
    df_eps_0[unlist(classifier), unlist(classifier_of_interest)] <- result_inner$d_observed$result_eps_0

  }
}


# df_eps_p explanation:
# an entry at with row x (correpsonding to algorithm x) and column y (corresponding to algorithm y)
# gives us an objectiv where algorithm x is included with minus in objective
# and algorihtm y is included with plus in objectiv
# all in all row gives - ; and column gives +
# With this we get that here the first row corresponds to the already computed
# comparisons needed for evaluating the test


# The Hasse graph of the empirical GSD relation for the PMLB datasets was created
# by first computing the value d_{62}(C,C') for all distinct pairs (C,C')  from the set
# {SVM, RF, CART, CRE, GLMNet, kNN } and then drawing a top-down edge from C to C',
# whenever d_{62}(C,C') >= 0. The resulting figure was created manually.

