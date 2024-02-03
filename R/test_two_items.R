test_two_items <- function(dat_set,
                           iteration_number = 1000,
                           seed_number = 2893,
                           eps_0 = 0,
                           eps_1 = 0.25,
                           eps_2 = 0.5,
                           eps_3 = 0.75,
                           eps_4 = 1,
                           mc.cores = 1,
                           reverse = FALSE) {

  ### Step 1: Compute the constraints given by R1 and R2
  # start_time_r1 <- Sys.time()
  constraint_r1_values <- compute_constraints_r1(dat_set)
  # total_time_r1 <- Sys.time() - start_time_r1


  # start_time_r2 <- Sys.time()
  constraint_r2_values <- compute_constraints_r2(dat_set,
                                                 constraint_r1_values$df_r1_values)
   # total_time_r2 <- Sys.time() - start_time_r2
  # saveRDS(constraint_r1_values, file = "constraint_r1_values_derma.rds")
  # saveRDS(constraint_r2_values, file = "constraint_r2_values_derma.rds")




  ### Step 2: Compute xi
  # Informations about gurobi (page 643ff):
  # https://www.gurobi.com/wp-content/plugins/hd_documentations/documentation/9.0/refman.pdf
  # (accessed: 07.02.2023)
  xi <- compute_xi(constraint_r1_values, constraint_r2_values, dim(dat_set)[1])
  # saveRDS(xi, "xi_derma.rds")


  ### Step 3: Compute the permutation test based on four different eps values


  # the general gurobi model
  gurobi_permu <- compute_gurobi_permu(eps_0, eps_1, eps_2, eps_3, eps_4,
                                       xi,
                                       constraint_r1_values,
                                       constraint_r2_values,
                                       dim(dat_set)[1])
  # Reduce dat_set only to the part needed in the permutation test
  dat_set_permu <- dat_set[, c("count_group_a", "count_group_b", "count_all")]

  # Computation of the test statistic values based on the input
  # start_time_d_obs  <- Sys.time()
  if (!reverse) {
    d_observed <- compute_d(1,
                            dat_set_permu,
                            gurobi_permu,
                            permutate_obs = FALSE,
                            reverse_objective = FALSE)
  }
  if (reverse) {
    d_observed <- compute_d(1,
                            dat_set_permu,
                            gurobi_permu,
                            permutate_obs = FALSE,
                            reverse_objective = TRUE)
    # total_time_d_obs <- Sys.time() - start_time_d_obs
    # saveRDS(d_observed, file = "d_observed_derma.rds")
  }


  ### Test statistic computation based on iteration_number permuted observations
  # Note that gurobi already parallels, thus parallelism does not necessarily
  # help to reduce the computation time
  # iteration_seq <- seq(1, iteration_number)
  # set.seed(2893)
  # RNGkind("L'Ecuyer-CMRG") # reproducibilty (note that one needs to use the same
  # # number of course as we did)
  #
  # start_time <- Sys.time()
  # permutation_test <- mclapply(X = iteration_seq, FUN = compute_d,
  #                            dat_set_permu = dat_set_permu,
  #                            gurobi_permu = gurobi_permu,
  #                            mc.cores = mc.cores)
  # total_time_permu <- Sys.time() - start_time



  iteration_seq <- seq(1, iteration_number)
  set.seed(seed_number)
  start_time <- Sys.time()
  permutation_test <- sapply(iteration_seq, FUN = compute_d,
                             dat_set_permu = dat_set_permu,
                             gurobi_permu = gurobi_permu)
  total_time_permu <- Sys.time() - start_time

  # saveRDS(permutation_test, file = "permutation_test_derma.rds")
  # saveRDS(total_time_permu, file = "total_time_permu_derma.rds")

  return(c(permutation_test = list(permutation_test),
           d_observed = list(d_observed)))

}
