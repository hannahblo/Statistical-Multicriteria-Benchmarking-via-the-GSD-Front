#' This function randomly samples a subset of a data set
#'
#' @param data_final: A data frame with column name count_group_a, count_group_b
#' and count_all with maximal and minimal value at the end
#' @param size_each_group: numeric value, not NULL, when both group sample have
#' the same size
#' @param size_group_a: numeric value, sample size of group a
#' @param size_group_b: numeric value, sample size of group b
#'
#' @return data_set, same structure as data_final but only a subset
sample_random_subset <- function(dat_final, size_each_group = 100,
                                 size_group_a = NULL, size_group_b = NULL) {

  if (!is.null(size_each_group)) {
    size_group_a <- size_each_group
    size_group_b <- size_each_group
  } else {
    if (is.null(size_group_a) || is.null(size_group_b)) {
      stop("Input error")
    }
  }

  dat_set <- dat_final
  cumsum_group_a <- cumsum(dat_set$count_group_a)
  cumsum_group_b <- cumsum(dat_set$count_group_b)
  sum_group_a <- sum(dat_set$count_group_a)
  sum_group_b <- sum(dat_set$count_group_b)
  dat_set$count_group_a <- 0
  dat_set$count_group_b <- 0

  # Compute a random sample
  group_a_sample <- sample(seq(0, sum_group_a), size = size_group_a, replace = FALSE)
  for (i in group_a_sample) {
    index_sample_value <- which.max(cumsum_group_a >= i)
    dat_set$count_group_a[index_sample_value] <-
      dat_set$count_group_a[index_sample_value] + 1
  }

  group_b_sample <- sample(seq(0, sum_group_b), size = size_group_b, replace = FALSE)
  for (i in group_b_sample) {
    index_sample_value <- which.max(cumsum_group_b >= i)
    dat_set$count_group_b[index_sample_value] <-
      dat_set$count_group_b[index_sample_value] + 1
  }


  dat_set$count_all <- dat_set$count_group_a + dat_set$count_group_b
  dat_set <- dat_set[-which(dat_set$count_all == 0), ]
  # add minimal and maximal element row (the last two rows of dat_final)
  # Note that the minimal and maximal value has neither a value for male nor female
  # -> does not increase the sample size
  dat_set[c(dim(dat_set)[1] + 1, dim(dat_set)[1] + 2), ] <-
    dat_final[c(dim(dat_final)[1] - 1, dim(dat_final)[1]), ]
  dat_set$ID <- seq(1, dim(dat_set)[1])

  return(dat_set)
}



#' This function computes xi
#'
#' @param constraint_r1_values: three vectors representing the constraints given
#' by R_1, this is a output of function compute_constraints_r1 from file
#' constraints_r1_r2.R
#' @param constraint_r2_values: three vectors representing the constraints given
#' by R_2, this is a output of function compute_constraints_r1 from file
#' constraints_r1_r2.R
#' @param number_observations: integer
#'
#' @return numeric value, the value of xi
compute_xi <- function(constraint_r1_values,
                       constraint_r2_values,
                       number_observations) {

  dim_r2 <- c(max(constraint_r2_values$r_2_i), number_observations)
  dim_r1 <- c(max(constraint_r1_values$r_1_i), number_observations)

  xi_r2_constraint <- rep(-1, dim_r2[1])
  gurobi_r2_sense <-  rep(">",  dim_r2[1])
  if (!any(is.na(constraint_r2_values$xi))) {
    xi_r2_constraint[constraint_r2_values$xi] <- 0
    gurobi_r2_sense[constraint_r2_values$xi] <- "="
  }
  gurobi_sense <- c(rep(">", dim_r1[1]), gurobi_r2_sense)

  # simple_triplet matrix i=row, j=column, v=value
  xi_A_r_i <- c(constraint_r1_values$r_1_i, seq(1, dim_r1[1]), # r_1 part
                constraint_r2_values$r_2_i + dim_r1[1], seq(dim_r1[1] + 1, dim_r1[1] + dim_r2[1])) # r_2 part



  xi_A_r_j <- c(constraint_r1_values$r_1_j, rep(dim_r1[2] + 1, dim_r1[1]),
                constraint_r2_values$r_2_j, rep(dim_r2[2] + 1, dim_r2[1]))

  xi_A_r_v <- c(constraint_r1_values$r_1_v, rep(-1, dim_r1[1]),
                constraint_r2_values$r_2_v, xi_r2_constraint)


  # Note that construcing the simple_triplet_matrix needs some minutes
  # set min and max to default zero and one -> by lower and upper bound constraints
  gurobi_model_xi <- list()
  gurobi_model_xi$A <- slam::simple_triplet_matrix(xi_A_r_i, xi_A_r_j, xi_A_r_v)
  gurobi_model_xi$rhs <- rep(0, max(xi_A_r_i))
  gurobi_model_xi$lb <- c(rep(0, dim_r1[2] - 1), 1, 0)
  gurobi_model_xi$ub <- c(rep(1, dim_r1[2] - 2), 0, 1, 1)
  gurobi_model_xi$vtypes <- c(rep("C", dim_r1[2] + 1))

  gurobi_model_xi$obj <- c(rep(0, dim_r1[2]), 1)
  gurobi_model_xi$modelsense <- "max"
  gurobi_model_xi$sense <- gurobi_sense

  if (dim(gurobi_model_xi$A)[2] != number_observations + 1) {
    stop("constraint matrix dimension and number_observations do not match. Hint: check if the ID corrsponds to the sequence seq(1, number different observations).")
  }
  xi_gurobi <- gurobi::gurobi(gurobi_model_xi)

  if (xi_gurobi$status == "OPTIMAL") {
    return(xi_gurobi$objval)
  } else {
    return(xi_gurobi)
  }
}


#' This function computes the gurobi model used in the permutation test
#'
#' @param eps_0 numeric between 0 and 1, regularization value
#' @param eps_1 numeric between 0 and 1, regularization value
#' @param eps_2 numeric between 0 and 1, regularization value
#' @param eps_3 numeric between 0 and 1, regularization value
#' @param eps_4 numeric between 0 and 1, regularization value
#' @param constraint_r1_values: three vectors representing the constraints given
#' by R_1, this is a output of function compute_constraints_r1 from file
#' constraints_r1_r2.R
#' @param constraint_r2_values: three vectors representing the constraints given
#' by R_2, this is a output of function compute_constraints_r1 from file
#' constraints_r1_r2.R
#' @param number_observations: integer
#'
#' @return numeric value, the value of xi
compute_gurobi_permu <- function(eps_0, eps_1, eps_2, eps_3, eps_4,
                                 xi,
                                 constraint_r1_values,
                                 constraint_r2_values,
                                 number_observations) {

  dim_r2 <- c(max(constraint_r2_values$r_2_i), number_observations)
  dim_r1 <- c(max(constraint_r1_values$r_1_i), number_observations)

  permu_A_r_i <- c(constraint_r1_values$r_1_i, # r_1 part
                   constraint_r2_values$r_2_i + dim_r1[1]) # r_2 part

  permu_A_r_j <- c(constraint_r1_values$r_1_j,
                   constraint_r2_values$r_2_j)

  permu_A_r_v <- c(constraint_r1_values$r_1_v,
                   constraint_r2_values$r_2_v)

  # Note that construction the simple_triplet_matrix needs some minutes
  # permu_A <- slam::simple_triplet_matrix(permu_A_r_i, permu_A_r_j, permu_A_r_v)

  # The right-hand side vector for the linear constraints
  eps_0_xi_r2_constraint <- rep(xi * eps_0, dim_r2[1])
  eps_1_xi_r2_constraint <- rep(xi * eps_1, dim_r2[1])
  eps_2_xi_r2_constraint <- rep(xi * eps_2, dim_r2[1])
  eps_3_xi_r2_constraint <- rep(xi * eps_3, dim_r2[1])
  eps_4_xi_r2_constraint <- rep(xi * eps_4, dim_r2[1])
  if (!any(is.na(constraint_r2_values$xi))) {
    eps_0_xi_r2_constraint[constraint_r2_values$xi] <- 0
    eps_1_xi_r2_constraint[constraint_r2_values$xi] <- 0
    eps_2_xi_r2_constraint[constraint_r2_values$xi] <- 0
    eps_3_xi_r2_constraint[constraint_r2_values$xi] <- 0
    eps_4_xi_r2_constraint[constraint_r2_values$xi] <- 0
  }
  permu_rhs_0 <- c(rep(xi * eps_0, dim_r1[1]), eps_0_xi_r2_constraint)
  permu_rhs_1 <- c(rep(xi * eps_1, dim_r1[1]), eps_1_xi_r2_constraint)
  permu_rhs_2 <- c(rep(xi * eps_2, dim_r1[1]), eps_2_xi_r2_constraint)
  permu_rhs_3 <- c(rep(xi * eps_3, dim_r1[1]), eps_3_xi_r2_constraint)
  permu_rhs_4 <- c(rep(xi * eps_4, dim_r1[1]), eps_4_xi_r2_constraint)

  # Gurobi Model (no regularization)
  # set minimales and maximales Element to 0 and 1 resp (using lower and upper bounds)
  gurobi_model_permu <- list()
  gurobi_model_permu$A <- slam::simple_triplet_matrix(permu_A_r_i, permu_A_r_j, permu_A_r_v)
  gurobi_model_permu$lb <- c(rep(0, dim_r1[2] - 1), 1)
  gurobi_model_permu$ub <- c(rep(1, dim_r1[2] - 2), 0, 1)
  gurobi_model_permu$vtypes <- c(rep("C", dim_r1[2]))


  gurobi_r2_sense <-  rep(">",  dim_r2[1])
  if (!any(is.na(constraint_r2_values$xi))) {
    gurobi_r2_sense[constraint_r2_values$xi] <- "="
  }
  gurobi_sense <- c(rep(">", dim_r1[1]), gurobi_r2_sense)
  gurobi_model_permu$modelsense <- "min"
  gurobi_model_permu$sense <- gurobi_sense

  if (dim(gurobi_model_permu$A)[2] != number_observations) {
    stop("constraint matrix dimension and number_observations do not match. Hint: check if the ID corrsponds to the sequence seq(1, number different observations).")
  }


  return(list(permu_rhs_0 = permu_rhs_0,
              permu_rhs_1 = permu_rhs_1,
              permu_rhs_2 = permu_rhs_2,
              permu_rhs_3 = permu_rhs_3,
              permu_rhs_4 = permu_rhs_4,
         gurobi_model_permu = gurobi_model_permu))
}





#' This function computes the test statistics (based on the sample and on a
#' permutation of the sample)
#'
#' @param index integer value (not needed in the function, only added for lapply)
#' @param dat_set_permu: A data frame with column name count_group_a, count_group_b
#' and count_all with maximal and minimal value at the end
#' @param gurobi_permu: list containing the gurobi model and five different right
#' hand sides of the constraints. Note this must be an object produced by
#' function compute_gurobi_permu()
#' @param permutate_obs: logical. If TRUE, then we permute the data set. If false,
#' then we use the sample
#' @param reverse_objective: logical. if the objective is (-1)*objective
#'
#' @return five numeric values, the test statistics based on the five different
#' right-hand side constraints
compute_d <- function(index,
                      dat_set_permu,
                      gurobi_permu,
                      permutate_obs = TRUE,
                      reverse_objective = FALSE) {



  all_obs <- sum(dat_set_permu$count_all)

  # Compute a random sample
  if (permutate_obs) {
    print(paste0("Now at index (permutation count) ", index))
    group_a_index <- sample(seq(1, all_obs),
                                 size = sum(dat_set_permu$count_group_a),
                                 replace = FALSE)
    group_a_sample <- rep(0, all_obs)
    group_a_sample[group_a_index] <- 1
    dat_set_permu <- dat_set_permu
    dat_set_permu$count_group_a <- 0

    index_sample <- 1
    for (i in 1:(dim(dat_set_permu)[1] - 2)) {
      # Note that for each i dat_set_permu$dup_all[i] must >= 1, since there is
      # at least one observation
      group_a_sample_i <- group_a_sample[seq(index_sample,
                                         index_sample + dat_set_permu$count_all[i] - 1)]
      dat_set_permu$count_group_a[i] <- sum(group_a_sample_i)
      index_sample <- index_sample + dat_set_permu$count_all[i]
    }
    # min und max never observed
    dat_set_permu$count_group_b <- dat_set_permu$count_all - dat_set_permu$count_group_a
    dat_set_permu[c(dim(dat_set_permu)[1] - 1, dim(dat_set_permu)[1]), ] <- 0
  }


  # compute objective function
  if (reverse_objective) {
    gurobi_permu$gurobi_model_permu$obj <-
       -(dat_set_permu$count_group_b / all_obs) + (dat_set_permu$count_group_a / all_obs)
  } else {
    gurobi_permu$gurobi_model_permu$obj <-
      (dat_set_permu$count_group_b / all_obs) - (dat_set_permu$count_group_a / all_obs)
  }



  # storing return list
  d_return <- list()

  # Computing d(x,y) for eps_0
  if (!any(is.na(gurobi_permu$permu_rhs_0))) {
    gurobi_permu$gurobi_model_permu$rhs <- gurobi_permu$permu_rhs_0
    result <- gurobi::gurobi(gurobi_permu$gurobi_model_permu)
    # params = list(Method = 0))
    if (result$status == "OPTIMAL") {
      d_return$result_eps_0 <- result$objva
    } else {
      d_return$result_eps_0 <- result
    }
  }


  # Computing d(x,y) for eps_1
  if (!any(is.na(gurobi_permu$permu_rhs_1))) {
    gurobi_permu$gurobi_model_permu$rhs <- gurobi_permu$permu_rhs_1
    result <- gurobi::gurobi(gurobi_permu$gurobi_model_permu)
    if (result$status == "OPTIMAL") {
      d_return$result_eps_1 <- result$objva
    } else {
      d_return$result_eps_1 <- result
    }
  }


  # Computing d(x,y) for eps_2
  if (!any(is.na(gurobi_permu$permu_rhs_2))) {
    gurobi_permu$gurobi_model_permu$rhs <- gurobi_permu$permu_rhs_2
    result <- gurobi::gurobi(gurobi_permu$gurobi_model_permu)
    if (result$status == "OPTIMAL") {
      d_return$result_eps_2 <- result$objva
    } else {
      d_return$result_eps_2 <- result
    }
  }


  # Computing d(x,y) for eps_3
  if (!any(is.na(gurobi_permu$permu_rhs_3))) {
    gurobi_permu$gurobi_model_permu$rhs <- gurobi_permu$permu_rhs_3
    result <- gurobi::gurobi(gurobi_permu$gurobi_model_permu)
    if (result$status == "OPTIMAL") {
      d_return$result_eps_3 <- result$objva
    } else {
      d_return$result_eps_3 <- result
    }
  }

  # Computing d(x,y) for eps_4
  if (!any(is.na(gurobi_permu$permu_rhs_4))) {
    gurobi_permu$gurobi_model_permu$rhs <- gurobi_permu$permu_rhs_4
    result <- gurobi::gurobi(gurobi_permu$gurobi_model_permu)
    if (result$status == "OPTIMAL") {
      d_return$result_eps_4 <- result$objva
    } else {
      d_return$result_eps_4 <- result
    }
  }


  return(d_return)
}

