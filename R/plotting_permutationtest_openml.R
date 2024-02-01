################################################################################
# Visualization of test statistics
################################################################################
library(ggplot2)
library(dplyr)
library(ggridges)
library(forcats)
library(latex2exp)
library(RColorBrewer)
library(ggthemes)

plotting_permutationtest_openml <- function(results_plots) {
  
  
  # prepare data frames for ggplot
  all_test_results = data.frame("rpart" = results_plots$classif.rpart$permutation_test[1,] %>% unlist(),
                                "knn" = results_plots$classif.kknn$permutation_test[1,] %>% unlist(),
                                "xgboost" = results_plots$classif.xgboost$permutation_test[1,] %>% unlist(),
                                "ranger" = results_plots$classif.ranger$permutation_test[1,] %>% unlist(),
                                "GLMnet" = results_plots$classif.glmnet$permutation_test[1,] %>% unlist(),
                                "Multinom" = results_plots$classif.multinom$permutation_test[1,] %>% unlist())
  all_test_results = all_test_results %>% rev()
  d_observed = c(
  results_plots$classif.rpart$d_observed$result_eps_0,
  results_plots$classif.kknn$d_observed$result_eps_0,
  results_plots$classif.xgboost$d_observed$result_eps_0,
  results_plots$classif.ranger$d_observed$result_eps_0,
  results_plots$classif.glmnet$d_observed$result_eps_0,
  results_plots$classif.multinom$d_observed$result_eps_0
  )
  d_observed = d_observed %>% unlist()
  names(d_observed) = c("rpart", "knn", "xgboost", "ranger", "GLMnet", "Multinom")
  
  # get critical values (0.05/6) (adapt manually to quantiles
  # from stat_density_ridges {ggridges}, see documentation in appendix of paper)
  cv_rpart = quantile(all_test_results[,1], probs = 0.05/6) +0.0015
  cv_knn = quantile(all_test_results[,2], probs = 0.05/6)
  cv_xgboost = quantile(all_test_results[,3], probs = 0.05/6) -0.0043
  cv_ranger = quantile(all_test_results[,4], probs = 0.05/6) +0.0045
  cv_GLMnet = quantile(all_test_results[,5], probs = 0.05/6)
  cv_Multinom = quantile(all_test_results[,6], probs = 0.05/6) -0.002
  
  # get critical values (0.05)
  cv_rpart_0.5 = quantile(all_test_results[,1], probs = 0.05) +0.0059
  cv_knn_0.5 = quantile(all_test_results[,2], probs = 0.05) -0.005
  cv_xgboost_0.5 = quantile(all_test_results[,3], probs = 0.05)
  cv_ranger_0.5 = quantile(all_test_results[,4], probs = 0.05) +0.0005
  cv_GLMnet_0.5 = quantile(all_test_results[,5], probs = 0.05) +0.0057
  cv_Multinom_0.5 = quantile(all_test_results[,6], probs = 0.05) -0.006
  
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
 # colorblind friendly
   cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
 # for geomsegment
  d_observed = unlist(d_observed) %>% rev
  d_observed_geom = data.frame(test = seq(1,6), d = d_observed )
  legend_title = ""
  df = stack(all_test_results )
  labels = fct_rev(df$ind) %>% levels # make sure labels align with levels in factors of df
  labels_segments = fct_rev(df$ind)
  labels = c("LR", "GLMNet", "RF", "xGBoost", "kNN", "CART")
  fill_rej_reg = 0.2
  
  # visualize test statistics including observed ones (see Figure 2 in paper)
  figure_2 = ggplot(df, aes(x = values, y = fct_rev(ind), fill = factor(stat(quantile)))) +
    xlim(c(-0.45,0.05)) +
    scale_fill_manual(
      name = latex2exp::TeX("  Decision:   "), values = c("red4", "#FF0000A0", "lightsteelblue"),
      labels = c(latex2exp::TeX("  Rejection for ${\\alpha} = \\frac{0.05}{6}$           "),
                 latex2exp::TeX("  Rejection for ${\\alpha} = 0.05$            "),
                                               latex2exp::TeX("No Rejection"))
    ) + 
    stat_density_ridges(
      geom = "density_ridges_gradient",
      calc_ecdf = TRUE,
      quantiles = c(0.05/6,0.05),
      quantile_fun = quantile,
      n = 512*2,
      jittered_points = TRUE,
      quantile_lines = FALSE,
      position = position_points_jitter(width = 0, height = 0),
      point_shape = '|', point_size = 4, point_alpha = 1, alpha = 3.5
    ) +
    theme(axis.title.x = element_text(margin = 20)) +
    theme(plot.title = element_text(hjust = 0.8)) +
    scale_y_discrete( name = "SVM vs.                             ", labels = labels) +
    coord_cartesian(clip = "off") +
    labs(title = legend_title) +
    theme_ridges(font_size = 18, grid = TRUE) +
    ggthemes::theme_economist_white(gray_bg = T) +
    theme(axis.title.x = element_blank()) +
    theme(legend.position = "bottom", legend.key.size = unit(1.2, 'cm'),
          legend.box.spacing = unit(1.8, 'cm'), 
          legend.text = element_text(size = 16),
          legend.title = element_text(size=21)) +
    guides(fill = guide_legend(override.aes = list(size = 14))) +
    theme(plot.background = element_rect(fill = 'white', colour = 'white')) +
    theme(axis.text=element_text(size=13), 
          axis.title=element_text(size=28), 
          plot.title = element_text(size = 32, face = "bold")) +
    theme(axis.title.y = element_text(size = 20, hjust= 0.7, margin = margin(t = 20, r = 20, b = 0, l = 0))) +
    geom_segment(data = d_observed_geom, aes(x = d, xend = d, y = test, yend = test +0.8),
                 color = "black", inherit.aes = FALSE, lineend = "round", linejoin = "round", size = 1.4) +
    annotate("rect", xmin=-0.45, xmax=cv_rpart, ymin=1, ymax=1.8,alpha=fill_rej_reg, fill = "red4") +
    annotate("rect", xmin=-0.45, xmax=cv_knn, ymin=2, ymax=2.8,alpha=fill_rej_reg, fill = "red4") +
    annotate("rect", xmin=-0.45, xmax=cv_xgboost, ymin=3, ymax=3.8,alpha=fill_rej_reg, fill = "red4") +
    annotate("rect", xmin=-0.45, xmax=cv_ranger, ymin=4, ymax=4.8,alpha=fill_rej_reg, fill = "red4") +    
    annotate("rect", xmin=-0.45, xmax=cv_GLMnet, ymin=5, ymax=5.8,alpha=fill_rej_reg, fill = "red4") +
    annotate("rect", xmin=-0.45, xmax=cv_Multinom, ymin=6, ymax=6.8,alpha=fill_rej_reg, fill = "red4") +
    annotate("rect", xmin=cv_rpart, xmax=cv_rpart_0.5, ymin=1, ymax=1.8,alpha=fill_rej_reg, fill = "#FF0000A0") +
    annotate("rect", xmin=cv_knn, xmax=cv_knn_0.5, ymin=2, ymax=2.8,alpha=fill_rej_reg, fill = "#FF0000A0") +
    annotate("rect", xmin=cv_xgboost, xmax=cv_xgboost_0.5, ymin=3, ymax=3.8,alpha=fill_rej_reg, fill = "#FF0000A0") +
    annotate("rect", xmin=cv_ranger, xmax=cv_ranger_0.5, ymin=4, ymax=4.8,alpha=fill_rej_reg, fill = "#FF0000A0") +    
    annotate("rect", xmin=cv_GLMnet, xmax=cv_GLMnet_0.5, ymin=5, ymax=5.8,alpha=fill_rej_reg, fill = "#FF0000A0") +
    annotate("rect", xmin=cv_Multinom, xmax=cv_Multinom_0.5, ymin=6, ymax=6.8,alpha=fill_rej_reg, fill = "#FF0000A0") +
    annotate("rect", xmin=cv_rpart_0.5, xmax=0.05, ymin=1, ymax=1.8,alpha=0.2, fill = "steelblue") +
    annotate("rect", xmin=cv_knn_0.5, xmax=0.05, ymin=2, ymax=2.8,alpha=fill_rej_reg, fill = "steelblue") +
    annotate("rect", xmin=cv_xgboost_0.5, xmax=0.05, ymin=3, ymax=3.8,alpha=fill_rej_reg, fill = "steelblue") +
    annotate("rect", xmin=cv_ranger_0.5, xmax=0.05, ymin=4, ymax=4.8,alpha=fill_rej_reg, fill = "steelblue") +    
    annotate("rect", xmin=cv_GLMnet_0.5, xmax=0.05, ymin=5, ymax=5.8,alpha=fill_rej_reg, fill = "steelblue") +
    annotate("rect", xmin=cv_Multinom_0.5, xmax=0.05, ymin=6, ymax=6.8,alpha=fill_rej_reg, fill = "steelblue") 
    


pdf(file= paste0("fig_2.pdf"), width = 13, height = 8)
print(figure_2)
dev.off()



  # function to compute shares of rejected resampled test statistics
  reject_share <- function(gamma, d_res, d_obs){
    d_obs <- rep(d_obs, length(d_res))
    sum(d_res - d_obs > 2*gamma/(1 - gamma))
  }
  # define grid of gamma values...
  gamma_grid = seq(0,0.2, by = 0.001) %>% as.list()
  #... and lapply share computation over this grid
  rejection_shares = lapply(gamma_grid,
                            FUN = reject_share,
                            d_res = all_test_results$Multinom,
                            d_obs = d_observed["Multinom"] )
  rej_shares_reg0 =  rejection_shares %>% unlist()
  df_rej_multinom = data.frame("gamma" = gamma_grid %>% unlist,
                        "Rejection share" = rej_shares_reg0,
                        "test" = "svm vs. multinom")

  rejection_shares = lapply(gamma_grid,
                            FUN = reject_share,
                            d_res = all_test_results$GLMnet,
                            d_obs = d_observed["GLMnet"] )
  rej_shares_reg0 =  rejection_shares %>% unlist()
  df_rej_glmnet = data.frame("gamma" = gamma_grid %>% unlist,
                               "Rejection share" = rej_shares_reg0,
                               "test" = "svm vs. GLMNet")

  rejection_shares = lapply(gamma_grid,
                            FUN = reject_share,
                            d_res = all_test_results$ranger,
                            d_obs = d_observed["ranger"] )
  rej_shares_reg0 =  rejection_shares %>% unlist()
  df_rej_ranger = data.frame("gamma" = gamma_grid %>% unlist,
                               "Rejection share" = rej_shares_reg0,
                               "test" = "svm vs. RF")

  rejection_shares = lapply(gamma_grid,
                            FUN = reject_share,
                            d_res = all_test_results$xgboost,
                            d_obs = d_observed["xgboost"] )
  rej_shares_reg0 =  rejection_shares %>% unlist()
  df_rej_xgboost = data.frame("gamma" = gamma_grid %>% unlist,
                               "Rejection share" = rej_shares_reg0,
                               "test" = "svm vs. xGBoost")

  rejection_shares = lapply(gamma_grid,
                            FUN = reject_share,
                            d_res = all_test_results$knn,
                            d_obs = d_observed["knn"] )
  rej_shares_reg0 =  rejection_shares %>% unlist()
  df_rej_knn = data.frame("gamma" = gamma_grid %>% unlist,
                               "Rejection share" = rej_shares_reg0,
                               "test" = "svm vs. kNN")

  rejection_shares = lapply(gamma_grid,
                            FUN = reject_share,
                            d_res = all_test_results$rpart,
                            d_obs = d_observed["rpart"] )
  rej_shares_reg0 =  rejection_shares %>% unlist()
  df_rej_rpart = data.frame("gamma" = gamma_grid %>% unlist,
                               "Rejection share" = rej_shares_reg0,
                               "test" = "svm vs. CART")



  ################################################################################
  # Visualization of rejection shares
  ################################################################################

  # only visualize unregularized test statistics, see footnote 7 in paper
  df_rej_all = rbind(df_rej_glmnet, df_rej_knn, df_rej_ranger,
                     df_rej_xgboost)
  par(mar=c(3,4,2,2))
  # to do: discrete palette
  pal = c("#DB6D00","#490092", "#56B4E9", "#009E73")
  # compute p-values
  df_rej_all$Rejection_share = 1- df_rej_all$Rejection.share/nrow(all_test_results)
  # compute number of contaminations
  number_resamples = length(all_test_results$Multinom)
  number_datasets = 80
  df_rej_all$k = df_rej_all$gamma * number_datasets
  #axis limits 
  x_limits = c(0,15)
  y_limits = c(0,0.35)

  # Eventually make plot of p-values (1- rejection shares) as function of contamination parameter gamma
  # (corresponds to figure 3 in paper)
  figure_3 = ggplot(data = df_rej_all) +
    geom_line(data = df_rej_all, aes(x = k, y = Rejection_share, colour = test), 
              size = 1.7, linejoin = "round") +
    labs(color = "Test") +
    ylab("   p-values   ") +
    xlab("Contaminated Samples") +
    ggthemes::theme_economist_white(gray_bg = T) +
    theme(axis.text=element_text(size=13), axis.title=element_text(size=28)) +
    theme(legend.key.size = unit(1, 'cm'),
          legend.key.height = unit(2, 'cm'),
          legend.key.width = unit(1, 'cm'),
          legend.title = element_text(size=20),
          legend.text = element_text(size=12),
          legend.text.align = 0,
          legend.spacing.x = unit(0.2,"cm"),
          axis.title.x = element_text(size = 22),
          axis.title.y = element_text(size = 22),
          plot.title = element_text(size = 32, face = "bold")) +
    theme(legend.position = "right", legend.key.size = unit(1.2, 'cm'),
          legend.box.spacing = unit(0.6, 'cm')) +
    theme(panel.background = element_rect(fill = 'white', colour = 'black')) +
    theme(plot.background = element_rect(fill = 'white', colour = 'white')) +
    ylim(y_limits) +
    xlim(x_limits) +
    geom_hline(yintercept=0.05/6, linetype="dashed", color = "red4", size = 2) +
    geom_hline(yintercept=0.05, linetype="dashed", color = "#FF0000A0", size = 2) +
    scale_color_manual(values=pal) +
    theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1)) +
    theme(axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1)) +
    theme(axis.title.x = element_text(margin = margin(t = 14, r = 20, b = 0, l = 0))) +
    theme(axis.title.y = element_text(margin = margin(t = 14, r = 20, b = 0, l = 0))) #+
    labs(title = latex2exp::TeX("                   Effect of Contamination")) 
    
  
    pdf(file= paste0("fig_3.pdf"), width = 12, height = 9)
    print(figure_3)
    dev.off()
  
  df_cdf = df
  levels(df_cdf$ind) = c("LR","GLMnet","RF","xGBoost","kNN","CART")
  
  ## visualize CDF
  figure_4 = ggplot(df_cdf, aes(values, colour = ind)) +
    stat_ecdf(size=1.22) +
    theme(plot.title = element_text(hjust = 1.8, margin = margin(t = 0, r = 20, b = 0, l = 0))) +
    ggtitle( "CDFs of resampled test statistics (OpenML)") +
    ggtitle( "") +
    theme(axis.title.x = element_blank()) +
    xlim(c(-0.15,0))+
    theme(legend.position = "bottom", legend.key.size = unit(1, 'cm')) +
    theme(axis.text=element_text(size=13), axis.title=element_text(size=28), plot.title = element_text(size = 20, face = "bold")) +
    theme(axis.title.y = element_text(size = 20)) +
    ggthemes::theme_economist_white(gray_bg = T) +
    labs(color = "Test: SVM vs.") +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    theme(legend.position = "bottom", legend.key.size = unit(1.2, 'cm'),
          legend.box.spacing = unit(1.8, 'cm'), 
          legend.text = element_text(size = 16),
          legend.title = element_text(size=21),
          plot.title = element_text(size = 28, face = "bold",hjust = 0.1, vjust = 2, margin = margin(t = 0, r = 20, b = 0, l = 10)))
  
  

  pdf(file= paste0("fig_4.pdf"), width = 12, height = 9)
  print(figure_4)
  dev.off()

  
 }

