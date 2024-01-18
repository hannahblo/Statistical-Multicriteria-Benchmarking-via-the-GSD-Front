################################################################################
# Visualization of test statistics
################################################################################

plotting_permutationtest <- function(permutation_test, d_observed,
                                     add_name_file = "") {
  # prepare data frames for ggplot
  all_test_results = data.frame("reg0" = permutation_test[1,] %>% unlist(), "reg0.25" = permutation_test[2,] %>% unlist(),
                                "reg0.5" = permutation_test[3,] %>% unlist(),"reg0.75" = permutation_test[4,] %>% unlist(),
                                "reg1" = permutation_test[5,] %>% unlist())
  all_test_results = all_test_results %>% rev()
  d_observed = d_observed %>% unlist()
  names(d_observed) = c("eps_0","eps_0.25","eps_0.5","eps_0.75","eps_1")


  ## theme for horizontal charts
  theme_flip <-
    theme(
      axis.text.x = element_text(face = "plain", family = "Roboto Mono", size = 22),
      axis.text.y = element_text(face = "bold", family = "Roboto", size = 26),
      panel.grid.major.x = element_line(color = "grey90", size = .9),
      #panel.grid.major.y = element_blank(),
      legend.position = "top",
      legend.text = element_text(family = "Roboto Mono", size = 18),
      legend.title = element_text(face = "bold", size = 18, margin = 25)#margin(b = 25))
    )


  # nice colors
  my_pal <- carto_pal(n = 8, name = "Bold")[c(1, 3, 7, 2)]
  # colorblind friendly
  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  # for geomsegment
  d_observed = unlist(d_observed)
  d_observed_geom = data.frame(test = seq(1,5), d = d_observed )

  legend_title = latex2exp::TeX("Distribution of test statistics ${d}_{I}$, ${d}^{\\epsilon}_{I}$, $\\underline{d}_{I}$ and $\\underline{d}^{\\epsilon}_{I}$")
  legend_title = latex2exp::TeX("               Distribution of Test Statistics")# ${d}^{\\epsilon}_{I}$")

  y_title = c(latex2exp::TeX("$\\epsilon = 0$"), latex2exp::TeX("$\\epsilon = 0.25$"),latex2exp::TeX("$\\epsilon = 0.5$"),
              latex2exp::TeX("$\\epsilon = 0.75$"),latex2exp::TeX("$\\epsilon = 1$"))


  df = stack(all_test_results )

  # visualize test statistics including observed ones (see Figure 2 in paper)
  figure_2 <- ggplot(df, aes(x = values, y = fct_rev(ind), fill = stat(x))) +
    geom_density_ridges_gradient(scale = 1.3, rel_min_height = 0,
                                 quantile_lines = TRUE, quantiles = 2,
                                 color = "black", alpha = 2.8, size = 0.8,
                                 jittered_points = TRUE,
                                 position = position_points_jitter(width = 0, height = 0),
                                 point_shape = '|', point_size = 2, point_alpha = 1, alpha = 1) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_discrete( name = "Regularization Strength                 ", labels = c("0", "0.25","0.5", "0.75","1")) +
    ggtitle(y_title) +
    xlim(c(-0.125,0.025)) +
    theme(axis.title.x = element_text(margin = 20)) +
    theme(plot.title = element_text(hjust = 0.8)) +
    scale_fill_viridis_c(name = "Value of Test Statistic", option = "F") +
    coord_cartesian(clip = "off") +
    labs(title = legend_title) +
    theme_ridges(font_size = 18, grid = TRUE) +
    #theme(axis.title.y = element_blank()) +
    theme(axis.title.x = element_blank()) +
    theme(legend.position = "bottom", legend.key.size = unit(1.2, 'cm')) +
    theme(axis.text=element_text(size=13), axis.title=element_text(size=28), plot.title = element_text(size = 20, face = "bold")) +
    theme(axis.title.y = element_text(size = 20))+
    geom_segment(data = d_observed_geom, aes(x = d, xend = d, y = test,
                                             yend = test +0.99), color = "red", size = 1)


  jpeg(file= paste0(add_name_file, "_values_teststatistic.jpeg"))
  print(figure_2)
  dev.off()


  # function to compute shares of rejected resampled test statistics
  # (see supp. C)
  reject_share <- function(gamma, d_res, d_obs){
    d_obs <- rep(d_obs, length(d_res))
    sum(d_obs - d_res > 2*gamma/(1 - gamma))
  }

  # check for gamma = 0
  reject_share(0, d_res = all_test_results$reg1, d_obs = d_observed["eps_1"] )
  reject_share(0, d_res = all_test_results$reg0.75, d_obs = d_observed["eps_0.75"] )
  reject_share(0, d_res = all_test_results$reg0.5, d_obs = d_observed["eps_0.5"] )
  reject_share(0, d_res = all_test_results$reg0.25, d_obs = d_observed["eps_0.25"] )
  reject_share(0, d_res = all_test_results$reg0, d_obs = d_observed["eps_0"] )

  # define grid of gamma values...
  gamma_grid = seq(0,0.04, by = 0.00001) %>% as.list()
  #... and lapply share computation over theses grids for all regularizations
  rejection_shares = lapply(gamma_grid,
                            FUN = reject_share,
                            d_res = all_test_results$reg1, d_obs = d_observed["eps_1"])
  rej_shares_reg1 =  rejection_shares %>% unlist()
  df_rej_1 = data.frame("gamma" = gamma_grid %>% unlist, "Rejection share" = rej_shares_reg1, "method" = "eps = 1")
  #
  rejection_shares = lapply(gamma_grid,
                            FUN = reject_share,
                            d_res = all_test_results$reg0.75, d_obs = d_observed["eps_0.75"])
  rej_shares_reg075 =  rejection_shares %>% unlist()
  df_rej_075 = data.frame("gamma" = gamma_grid %>% unlist, "Rejection share" = rej_shares_reg075, "method" = "eps = 0.75")
  #
  rejection_shares = lapply(gamma_grid,
                            FUN = reject_share,
                            d_res = all_test_results$reg0.5, d_obs = d_observed["eps_0.5"] )
  rej_shares_reg05 =  rejection_shares %>% unlist()
  df_rej_05 = data.frame("gamma" = gamma_grid %>% unlist, "Rejection share" = rej_shares_reg05, "method" = "eps = 0.5")
  #
  rejection_shares = lapply(gamma_grid,
                            FUN = reject_share,
                            d_res = all_test_results$reg0.25, d_obs = d_observed["eps_0.25"])
  rej_shares_reg025 =  rejection_shares %>% unlist()
  df_rej_025 = data.frame("gamma" = gamma_grid %>% unlist, "Rejection share" = rej_shares_reg025, "method" = "eps = 0.25")
  #
  rejection_shares = lapply(gamma_grid,
                            FUN = reject_share,
                            d_res = all_test_results$reg0, d_obs = d_observed["eps_0"] )
  rej_shares_reg0 =  rejection_shares %>% unlist()
  df_rej_0 = data.frame("gamma" = gamma_grid %>% unlist, "Rejection share" = rej_shares_reg0, "method" = "eps = 0")




  ################################################################################
  # Visualization of rejection shares
  ################################################################################

  # only visualize regularized test statistics
  df_rej_all = rbind(df_rej_025, df_rej_05, df_rej_075, df_rej_1)
  method = latex2exp::TeX("$\\epsilon = 0$")

  par(mar=c(3,4,2,2))

  # color grid for regulaizations
  pal = brewer.pal(5, "YlGnBu")[c(2,3,4,5)]

  # compute p-values
  df_rej_all$Rejection_share = 1 - df_rej_all$Rejection.share/nrow(all_test_results)

  # Eventually make plot of p-avlues (1 - rejection shares) as function of contamination parameter gamma
  # (corresponds to figure 3 in paper)
  figure_3 <- ggplot(data = df_rej_all, aes(x = gamma, group = method)) +
    #geom_point(data = df_rej_all, aes(x = gamma, y = Rejection.share, colour = method))  +
    geom_line(data = df_rej_all, aes(x = gamma, y = Rejection_share, colour = method), size = 1.7) +
    labs(color = "Regularization") +
    ylab("p-values") +
    xlab(unname(TeX(c("$\\gamma$")))) +
    theme_minimal() +
    theme(axis.text=element_text(size=13), axis.title=element_text(size=28)) +
    theme(legend.key.size = unit(1, 'cm'),
          legend.key.height = unit(2, 'cm'),
          legend.key.width = unit(1, 'cm'),
          legend.title = element_text(size=20),
          legend.text = element_text(size=20),
          legend.text.align = 0,
          axis.title.x = element_text(size = 30),
          axis.title.y = element_text(size = 20)) +
    ylim(c(0, 0.2)) +
    xlim(c(0,0.012)) +
    geom_hline(yintercept=0.05, linetype="dashed", color = "red", size = 2) +
    scale_color_manual(values=pal,
                       labels = unname(TeX(c(#"$$\\epsilon = 0$$",
                         "$$\\epsilon = 0.25$$",
                         "$$\\epsilon = 0.5$$",
                         "$$\\epsilon = 0.75$$",
                         "$$\\epsilon = 1$$"))))


  jpeg(file= paste0(add_name_file, "_regularization.jpeg"))
  print(figure_3)
  dev.off()

}

