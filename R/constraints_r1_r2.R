#' This functions computes the constraints given by R_1
#' Note that, in contrast to the article, the implementation was done with having
#' the 'less or equal' and not the 'greater or equal' relation in mind.
#'
#' @param data_set: This is a data_frame which has as column names ordinal_1,
#' ordinal_2, numeric
#'
#' @returns the row, column and values as vectors to then use a simple_triplet_matrix
#' and a data frame representing all pairs of R_1
compute_constraints_r1 <- function(data_set) {


  # Unique combination of ordinal observations
  numb_unique_ordinal <- dim(unique(data_set[,c("ordinal_1", "ordinal_2")])) [1]

  # Sorted first by ordinal_1, then within a ordinal_1 categories is sorted by
  # ordinal_2 and analogously then by numeric.
  sort_dat <- data_set[
    order(data_set[["ordinal_1"]], data_set[["ordinal_2"]], data_set[["numeric"]]),
  ]


  # Computation when in the sorted data a switch in the ordinal data is
  value_bevor <- sort_dat[1, c("ordinal_1", "ordinal_2")]
  index_switch_ordinal <- 1
  i_start_ordinal_groups <- data.frame(ordinal_1 = rep(NA, numb_unique_ordinal),
                                       ordinal_2 = rep(NA, numb_unique_ordinal),
                                       from = rep(NA, numb_unique_ordinal),
                                       to = rep(NA, numb_unique_ordinal))
  i_start_ordinal_groups[1, ] <- c(sort_dat[1,c("ordinal_1", "ordinal_2")],
                                   1, NA)

  # We go through sort_dat and every time there is a change in the ordinal part,
  # we store the value as well as from where to where that combination goes
  # Note that this is based on sort_dat and not on dat
  for (i in 2:dim(sort_dat)[1]) {
    value_current <- sort_dat[i, c("ordinal_1", "ordinal_2")]
    if (any(value_current != value_bevor)) {
      i_start_ordinal_groups[index_switch_ordinal + 1, ] <- c(value_current,
                                                              i, NA)
      i_start_ordinal_groups[index_switch_ordinal, 4] <- i - 1
      index_switch_ordinal <- index_switch_ordinal + 1
      value_bevor <- value_current
    }
  }
  # We have to add for the last group the "to" part
  i_start_ordinal_groups[numb_unique_ordinal, 4] <- dim(sort_dat)[1]


  # Now, we compute r_1
  n <- dim(sort_dat)[1]
  df_index <- 1
  r_1_i <- c()
  r_1_j <- c()
  r_1_v <- c()

  # parallel to the calculation of r_1 we calculate df_r1_values, which stores the
  # relation r_1, but which can then be sorted again and used for the calculation
  # of r_2
  df_r1_values <- data.frame(matrix(ncol = 7, nrow = 0))
  colnames(df_r1_values) <- c("ID_lower", "ID_upper",
                              "ordinal_1_lower", "ordinal_2_lower",
                              "ordinal_1_upper", "ordinal_2_upper",
                              "difference_numeric")



  # We go through every observation
  for (i in 1:n) {
    # Value comparing which value is above basis_value
    basis_value <- sort_dat[i, c("ordinal_1", "ordinal_2", "numeric")]
    # First, we consider only the ordinal part
    larger_ordinal <- intersect(
      which(i_start_ordinal_groups[, 1] >= basis_value[1, 1]),
      which(i_start_ordinal_groups[, 2] >= basis_value[1, 2]))
    # base_entry corresponds to the
    base_entry <- 0
    index_below <- sort_dat[i, "ID"]

    # We go through every group of ordinal combinations which fulfill being larger
    # then basis value.
    for (j in larger_ordinal) {
      # First entry where also the numeric part is larger
      first_income_above <-
        purrr::detect_index(sort_dat[seq(from = i_start_ordinal_groups[j, 3],
                                         to = i_start_ordinal_groups[j,  4]),
                                     "numeric"],
                            function(x) {x >= as.numeric(basis_value[3])})

      # if there does not exist a first entry, because every numeric value is
      # below, then first_income_above = 0
      if (!(first_income_above == 0)) {
        # Based on sort_dat every observation which lies in the same ordinal group
        # but has a higher index than first_income_above is larger then base_value
        # (in all components)
        index_above <- seq(i_start_ordinal_groups[j,  3] + first_income_above - 1,
                           i_start_ordinal_groups[j,  4], 1)
        # deleting the reflexive comparison
        index_above <- index_above[!(sort_dat[index_above, "ID"] == index_below)]

        if (!(length(index_above) == 0)) {
          # R_1 is sorted by the order of the original saving of the data, thus we
          # have to use the "ID" of sort_dat and not the index of sort_dat

          # Further, we save df_r1_values. Here we serve the combinations with
          # value, ID and difference in the numeric part
          for (index_above_df in index_above) {
            r_1_i <- c(r_1_i, c(df_index, df_index))
            r_1_j <- c(r_1_j, c(index_below, sort_dat[index_above_df, "ID"]))
            r_1_v <- c(r_1_v, c(-1, 1))
            df_r1_values[df_index, ] <-
              c(index_below, sort_dat[index_above_df, "ID"],
                basis_value[["ordinal_1"]], basis_value[["ordinal_2"]],
                sort_dat[index_above_df, "ordinal_1"],
                sort_dat[index_above_df, "ordinal_2"],
                as.numeric(sort_dat[index_above_df, "numeric"]) -
                  as.numeric(basis_value[3]))
            df_index <- df_index + 1
          }
        }
      }
    }
  }
  return(list(r_1_i = r_1_i, r_1_j = r_1_j, r_1_v = r_1_v,
              df_r1_values = df_r1_values))
}



#' This functions computes the constraints given by R_2
#' Note that, in contrast to the article, the implementation was done with having
#' the 'less or equal' and not the 'greater or equal' relation in mind.
#'
#' @param data_set: This is a data_frame which has as column names ordinal_1,
#' ordinal_2, numeric
#' @param df_r1_values data frame representing all pairs of R_1 and computed by
#' the function compute_constraints_r1
#'
#' @returns the row, column and values as vectors to be later used for a simple_triplet_matrix
#' further, we return a vector which gives all equality constraints
compute_constraints_r2 <- function(data_set, df_r1_values) {
  # Sorted first by health_lower, then within a health_lower category is sorted by
  # education_lower, then within this group by health_upper, then
  # by education_lower and finally by difference_numeric
  #TODO
  # wie order?
  sort_df_r1 <- df_r1_values[
    order(df_r1_values[["ordinal_1_lower"]],
          df_r1_values[["ordinal_2_lower"]],
          as.numeric(as.factor(df_r1_values[["ordinal_1_upper"]])),
          as.numeric(as.factor(df_r1_values[["ordinal_2_upper"]])),
          -df_r1_values[["difference_numeric"]]
    ), ]

  # Unique combination of ordinal observations
  numb_unique_ordinal <- dim(unique(sort_df_r1[, c("ordinal_1_lower",
                                                   "ordinal_2_lower",
                                                   "ordinal_1_upper",
                                                   "ordinal_2_upper")])) [1]

  value_bevor <- sort_df_r1[1, c("ordinal_1_lower", "ordinal_2_lower",
                                 "ordinal_1_upper", "ordinal_2_upper")]


  # Computation when in the sorted data a switch in the ordinal data is
  index_switch_ordinal <- 1
  i_df_ordinal_groups <- data.frame("ordinal_1_lower" = rep(NA, numb_unique_ordinal),
                                    "ordinal_2_lower" = rep(NA, numb_unique_ordinal),
                                    "ordinal_1_upper" = rep(NA, numb_unique_ordinal),
                                    "ordinal_2_upper" = rep(NA, numb_unique_ordinal),
                                    from = rep(NA, numb_unique_ordinal),
                                    to = rep(NA, numb_unique_ordinal))
  i_df_ordinal_groups[1, ] <- c(sort_df_r1[1, c("ordinal_1_lower", "ordinal_2_lower",
                                                "ordinal_1_upper", "ordinal_2_upper")], 1, NA)

  # We go through sort_df_r1 and every time there is a change in the ordinal part,
  # we store the value as well as from where to where that combination goes
  for (i in 2:dim(sort_df_r1)[1]) {
    # This code is analogously to the one in the r_1 part
    value_current <- sort_df_r1[i, c("ordinal_1_lower", "ordinal_2_lower",
                                     "ordinal_1_upper", "ordinal_2_upper")]
    if (any(value_current != value_bevor)) {
      i_df_ordinal_groups[index_switch_ordinal + 1, ] <- c(value_current, i, NA)
      i_df_ordinal_groups[index_switch_ordinal, 6] <- i - 1
      index_switch_ordinal <- index_switch_ordinal + 1
      value_bevor <- value_current
    }
  }
  # We have to add for the last group the "to" part
  i_df_ordinal_groups[numb_unique_ordinal, 6] <- dim(sort_df_r1)[1]


  # Now, we compute r_2, each row represents one pair and every other pair which
  # lies above
  r_2_i <- rep(NA,  dim(sort_df_r1)[1] * dim(sort_df_r1)[1] * 4)
  r_2_j <- rep(NA,  dim(sort_df_r1)[1] * dim(sort_df_r1)[1] * 4)
  r_2_v <- rep(NA,  dim(sort_df_r1)[1] * dim(sort_df_r1)[1] * 4)
  later_delta <-  rep(NA,  dim(sort_df_r1)[1] * dim(sort_df_r1)[1] * 4)
  i_row <- 1
  i_saving <- 1
  delta_saving <- 1

  start_time <- Sys.time()
  # We go through every existing relation in r_1 which is given by sort_df_r1
  # This code is analogously to the one in the r_1 part
  for (i in 1:dim(sort_df_r1)[1]) {
    basis_value <- sort_df_r1[i, c("ordinal_1_lower", "ordinal_2_lower",
                                   "ordinal_1_upper", "ordinal_2_upper",
                                   "difference_numeric")]
    print(paste0("Now at row ", i, " of sort_df_r1 with total number of ", dim(sort_df_r1)[1], " rows."))
    basis_value_id <- sort_df_r1[i, c("ID_lower", "ID_upper")]
    larger_ordinal <- Reduce(intersect,
                             list(which(i_df_ordinal_groups[, 1] >= basis_value[1, 1]),
                                  which(i_df_ordinal_groups[, 2] >= basis_value[1, 2]),
                                  which(i_df_ordinal_groups[, 3] <= basis_value[1, 3]),
                                  which(i_df_ordinal_groups[, 4] <= basis_value[1, 4])))
    # due to reflexivity,  larger_ordinal is never empty
    for (j in larger_ordinal) {
      difference_above_1 <-
        purrr::detect_index(sort_df_r1[seq(from = i_df_ordinal_groups[j, "from"],
                                           to = i_df_ordinal_groups[j,  "to"]),
                                       "difference_numeric"],
                            function(x) {x <= as.numeric(basis_value[5])})

      if (!(difference_above_1 == 0)) {

        index_above <- seq(i_df_ordinal_groups[j,  "from"] + difference_above_1 - 1,
                           i_df_ordinal_groups[j,  "to"],
                           1)


        # delete the relation with itself
        delete_reflexiv <- intersect(which(sort_df_r1[index_above, 1] == basis_value_id[["ID_lower"]]),
                                     which(sort_df_r1[index_above, 2] == basis_value_id[["ID_upper"]]))
        if (!(length(delete_reflexiv) == 0)) {
          index_above <- index_above[-delete_reflexiv]
        }
        if (!(length(index_above) == 0)) {
          for (k in index_above) {

            if (sort_df_r1[k, 1] == basis_value_id[["ID_upper"]]) {

              r_2_i[seq(i_saving, i_saving + 2)] <- c(i_row, i_row, i_row)
              r_2_j[seq(i_saving, i_saving + 2)] <-  c(sort_df_r1[k, 1], sort_df_r1[k, 2],
                                                       basis_value_id[["ID_lower"]])
              r_2_v[seq(i_saving, i_saving + 2)] <- c(2,-1,-1)
              i_saving <- i_saving + 3
              if (all(basis_value == sort_df_r1[k, c(3,4,5,6,7)])) {
                later_delta[delta_saving] <- i_row
                delta_saving <- delta_saving + 1
              }
              i_row <- i_row + 1

            } else if (sort_df_r1[k, 2] == basis_value_id[["ID_lower"]]) {
              r_2_i[seq(i_saving, i_saving + 2)] <- c(i_row, i_row, i_row)
              r_2_j[seq(i_saving, i_saving + 2)] <-  c(sort_df_r1[k, 1], sort_df_r1[k, 2],
                                                       basis_value_id[["ID_upper"]])
              r_2_v[seq(i_saving, i_saving + 2)] <- c(1,-2,1)
              i_saving <- i_saving + 3
              if (all(basis_value == sort_df_r1[k, c(3,4,5,6,7)])) {
                later_delta[delta_saving] <- i_row
                delta_saving <- delta_saving + 1
              }
              i_row <- i_row + 1

            } else if (!(sort_df_r1[k, 1] == basis_value_id[["ID_lower"]]) & !(sort_df_r1[k, 2] == basis_value_id[["ID_upper"]])) {
              r_2_i[seq(i_saving, i_saving + 3)] <- c(i_row, i_row, i_row, i_row)
              r_2_j[seq(i_saving, i_saving + 3)] <-  c(sort_df_r1[k, 1], sort_df_r1[k, 2],
                                                       basis_value_id[["ID_lower"]], basis_value_id[["ID_upper"]])
              r_2_v[seq(i_saving, i_saving + 3)] <- c(1,-1,-1,1)
              i_saving <- i_saving + 4
              if (all(basis_value == sort_df_r1[k, c(3,4,5,6,7)])) {
                later_delta[delta_saving] <- i_row
                delta_saving <- delta_saving + 1
              }
              i_row <- i_row + 1

            }
          }
        }
      }
    }
  }


  first_na_const <- purrr::detect_index(r_2_i, function(x) {is.na(x)})
  first_na_delta <- purrr::detect_index(later_delta, function(x) {is.na(x)})
  return(list(r_2_i = r_2_i[seq(1, first_na_const - 1)], r_2_j =  r_2_j[seq(1, first_na_const - 1)],
         r_2_v = r_2_v[seq(1, first_na_const - 1)], xi = later_delta[seq(1, first_na_delta - 1)]))
}

