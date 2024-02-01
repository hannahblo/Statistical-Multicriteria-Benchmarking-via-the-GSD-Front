# This file corresponds to the openml analysis

# Analyzing classifiers based on data sets, evaluations, and classifiers
# provided by the openly available OpenML repository
# see: https://www.openml.org/ (Accessed: 18.01.2024)




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
library(OpenML) # downloading data from openml
library(farff) # for openml data preparation

source("R/constraints_r1_r2.R") # contains the functions compute_constraints...
source("R/sample_permutation_test.R") # permutation test, sample etc
source("R/plotting_permutationtest_openml.R") # plot function
source("R/test_two_items.R") # main function summarizing the computation
################################################################################
# Prepare Data Set: OpenML
################################################################################

# Load the data set from OpenML (see https://www.openml.org/)
data_all <- OpenML::listOMLDataSets()
rel_datasets <- data_all %>% filter(status == "active" &
                                      number.of.classes == 2 &
                                      number.of.features < 100 &
                                      number.of.instances < 1000 &
                                      number.of.instances > 100 &
                                      number.of.instances.with.missing.values == 0 &
                                      max.nominal.att.distinct.values < 5
)

#### flows 2333 (rpart), 2330 (ranger), 2409 (knn), 2337 (xgboost), 2408 (glmnet), 2336(svm), 2317(logit), 2313 lda
# test = getOMLFlow(flow.id = 2333)
#test = getOMLTask(task.id = 3729)
#### 4689


flows_for_paper <- c(2333, 2330, 2317, 2337, 2408, 2336,  2409)
outls <- list()
for (i in 1:length(flows_for_paper)) {
  temp <- listOMLRunEvaluations(flow.id = flows_for_paper[i], limit = 10000)
  if (i == 1) {
    datasets = temp$data.name
  } else {
    datasets = intersect(datasets,temp$data.name)
  }
  outls[[i]] = temp
}
data_openml <- do.call('rbind', outls)
data_openml <-  data_openml[data_openml$data.name %in% datasets,]
data_openml <- data_openml %>% group_by(flow.id, data.name) %>% slice(n())



extractDataId <- function(taskid){
  print(taskid)
  if (length(tryCatch({res = getOMLTask(task.id = taskid)}, error = function(e){return(NA)})) <= 1) {
    return(NA)
  }
  return(res$input$data.set$desc$id)
}
data_openml$data.id = sapply(data_openml$task.id,function(x)extractDataId(taskid = x))


# data_final[which(is.na(data_final$data.id)), ]
# taskid = 3019: Data set has been deactivated. Thus delete this one
if (length(which(is.na(data_openml$data.id))) > 0) {
  data_openml <- data_openml[-which(is.na(data_openml$data.id)), ]
}


# We are interested in binary classification with sample size between 450 and
# 10000
data_openml_filter = data_openml %>%
  group_by(data.id) %>%
  dplyr::mutate(count = n()) %>%
  left_join(data_all, by = "data.id") %>%
  filter((count == length(flows_for_paper)) &
           (number.of.instances.x > 450) &
           (number.of.instances.x < 10000) & # 10000 500
           (number.of.classes == 2)
  )
data_openml_filter <- data_openml_filter[order(data_openml_filter$data.name), ]

# saveRDS(data_openml_filter, "data_openml_filter.rds")
# data_openml_filter <- readRDS("data_openml_filter.rds")


# We select those performance measures using for the evaluation
# numeric information: predictive.accuracy
# ordinal information: usercpu.time.accuracy and root.mean.square.error
data_openml_filter <- data_openml_filter[, c("learner.name",
                                               "predictive.accuracy",
                                               "usercpu.time.millis.training",
                                               "usercpu.time.millis.testing"
)]

# We discretize both cpu times
cuttings_test <- quantile(data_openml_filter$usercpu.time.millis.testing, probs = seq(0, 1, 0.1))
data_openml_filter$usercpu.time.millis.testing <- findInterval(data_openml_filter$usercpu.time.millis.testing, cuttings_test)
data_openml_filter$usercpu.time.millis.testing <- max(data_openml_filter$usercpu.time.millis.testing) + 1 - data_openml_filter$usercpu.time.millis.testing

cutting_train <- quantile(data_openml_filter$usercpu.time.millis.training, probs = seq(0, 1, 0.1))
data_openml_filter$usercpu.time.millis.training <- findInterval(data_openml_filter$usercpu.time.millis.training, cutting_train)
data_openml_filter$usercpu.time.millis.training <- max(data_openml_filter$usercpu.time.millis.training) + 1 - data_openml_filter$usercpu.time.millis.training

################################################################################
# Conducting the permutation test and plotting the results
################################################################################

classifier_of_interest <- "classif.svm"
classifiers_comparison <- list( "classif.multinom", "classif.ranger", "classif.xgboost",
                                "classif.glmnet", "classif.kknn", "classif.rpart")


for (classifier in classifiers_comparison) {

  data_openml_selected <- data_openml_filter

  # Now we select two classifiers and compare them
  # here we choose classif.rpart and classif.svm
  data_openml_selected <- data_openml_selected[data_openml_selected$learner.name %in%
                                                 c(classifier_of_interest,
                                                   classifier), ]


  # now, we need to convert the data_final_selected into a dataframe with columns:
  # ordinal_1, ordinal_2, numeric, count_group_a, count_group_b, count_all, ID

  # Step 1: Converting the variables of interest into numeric and order modes
  data_openml_selected[["predictive.accuracy"]] <- as.numeric(as.character(
    data_openml_selected[["predictive.accuracy"]]))

  data_openml_selected[["usercpu.time.millis.training"]] <- as.ordered(as.character(
    data_openml_selected[["usercpu.time.millis.training"]]))

  data_openml_selected[["usercpu.time.millis.testing"]] <- as.ordered(as.character(
    data_openml_selected[["usercpu.time.millis.testing"]]))


  # Step 2: duplication handling
  data_count <- data_openml_selected %>% group_by_all() %>% count()


  data_algorithm_1 <- data_count[which(data_count$learner.name == classifier_of_interest),
                                 c(2, 3, 4, 5)]
  data_algorithm_1 <- matrix(as.numeric(as.matrix(data_algorithm_1)), ncol = 4)
  colnames(data_algorithm_1) <-  c("numeric", "ordinal_1", "ordinal_2", "count_group_a")


  data_algorithm_2 <- data_count[which(data_count$learner.name == classifier),
                                 c(2, 3, 4, 5)]
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
  # input at the second largest row is the minimal value and the largest row
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
  start_time <- Sys.time()
  result_inner <- test_two_items(dat_set, iteration_number = 1000)
  total_time <- Sys.time() - start_time

  saveRDS(result_inner, paste0(classifier, "_result.rds"))
  saveRDS(total_time, paste0(classifier, "_computation_time.rds"))
}

# computing the test statistic values
classifier_of_interest <- "classif.svm"
classifiers_comparison <- list( "classif.multinom", "classif.ranger", "classif.xgboost", "classif.glmnet", "classif.kknn", "classif.rpart")
all_eps_values <- c("result_eps_0", "result_eps_1", "result_eps_2", "result_eps_3", "result_eps_4")

proportion_below_df <- as.data.frame(matrix(rep(0, 6 * 5), nrow = 5, ncol = 6), row.names = all_eps_values)
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
classifiers_comparison <- list( "classif.multinom", "classif.ranger", "classif.xgboost", "classif.glmnet", "classif.kknn", "classif.rpart")

results_plots = list()
for (classifier in classifiers_comparison) {
  # if(classifier == "classif.xgboost")
  #   debugonce(plotting_permutationtest)
  result_plot <- readRDS(paste0(classifier, "_result.rds"))
  results_plots[[classifier]] = result_plot
}
plotting_permutationtest_openml(results_plots)

#
# plotting_permutationtest_openml(result_plot)
#
#
# for (classifier in classifiers_comparison) {
#   result_plot <- readRDS(paste0(classifier, "_result.rds"))
#   plotting_permutationtest_openml(result_plot)
# }









# Constructing the graph representing all pairwise comparisons at once
# therfore we have to compute the reverse (switching classifier and classifier_of_interest
# position) observed minimal difference
classifiers_all <- list( "classif.svm",  "classif.multinom", "classif.ranger", "classif.xgboost", "classif.glmnet", "classif.kknn", "classif.rpart")

for (k in seq(1, length(classifiers_all))) {
  for (m in seq(1, length(classifiers_all))[-k]) {

    classifier <- classifiers_all[m]
    classifier_of_interest <- classifiers_all[k]


    data_openml_selected <- data_openml_filter

    # Now we select two classifiers and compare them
    # here we choose classif.rpart and classif.svm
    data_openml_selected <- data_openml_selected[data_openml_selected$learner.name %in%
                                                   c(classifier_of_interest,
                                                     classifier), ]


    # now, we need to convert the data_final_selected into a dataframe with columns:
    # ordinal_1, ordinal_2, numeric, count_group_a, count_group_b, count_all, ID

    # Step 1: Converting the variables of interest into numeric and order modes
    data_openml_selected[["predictive.accuracy"]] <- as.numeric(as.character(
      data_openml_selected[["predictive.accuracy"]]))

    data_openml_selected[["usercpu.time.millis.training"]] <- as.ordered(as.character(
      data_openml_selected[["usercpu.time.millis.training"]]))

    data_openml_selected[["usercpu.time.millis.testing"]] <- as.ordered(as.character(
      data_openml_selected[["usercpu.time.millis.testing"]]))


    # Step 2: duplication handling
    data_count <- data_openml_selected %>% group_by_all() %>% count()


    data_algorithm_1 <- data_count[which(data_count$learner.name == classifier),
                                   c(2, 3, 4, 5)]
    data_algorithm_1 <- matrix(as.numeric(as.matrix(data_algorithm_1)), ncol = 4)
    colnames(data_algorithm_1) <-  c("numeric", "ordinal_1", "ordinal_2", "count_group_a")


    data_algorithm_2 <- data_count[which(data_count$learner.name == classifier_of_interest),
                                   c(2, 3, 4, 5)]
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
    result_inner <- test_two_items(dat_set, iteration_number = 1)

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



# The Hasse graph of the empirical GSD relation for the OpenML datasets was created
# by first computing the value d_{80}(C,C') for all distinct pairs (C,C')  from the
# set {SVM, RF, CART, LR, GLMNet, xGBoost, kNN } and then drawing a top-down edge
# from C to C', whenever d_{80}(C,C') >= 0. The resulting figure was created manually.


