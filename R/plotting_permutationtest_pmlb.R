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

plotting_permutationtest_pmlb <- function(results_plots) {


# prepare data frames for ggplot
all_test_results = data.frame("rpart" = results_plots$J48$permutation_test[1,] %>% unlist(),
                              "knn" = results_plots$knn$permutation_test[1,] %>% unlist(),
                              "ranger" = results_plots$ranger$permutation_test[1,] %>% unlist(),
                              "GLMnet" = results_plots$glmnet$permutation_test[1,] %>% unlist(),
                              "SVM" = results_plots$svmRadial$permutation_test[1,] %>% unlist())
all_test_results = all_test_results %>% rev()
d_observed = c(
  results_plots$J48$d_observed$result_eps_0,
  results_plots$knn$d_observed$result_eps_0,
  results_plots$ranger$d_observed$result_eps_0,
  results_plots$glmnet$d_observed$result_eps_0,
  results_plots$svmRadial$d_observed$result_eps_0
)
d_observed = d_observed %>% unlist()
names(d_observed) = c("rpart", "knn", "ranger", "GLMnet", "SVM")

# get critical values (0.05/6) (adapt manually to quantiles from stat_density_ridges {ggridges},
# see documentation in appendix)
cv_rpart = quantile(all_test_results[,1], probs = 0.05/6) -0.005
cv_knn = quantile(all_test_results[,2], probs = 0.05/6) -0.002
cv_ranger = quantile(all_test_results[,3], probs = 0.05/6) #+0.0045
cv_GLMnet = quantile(all_test_results[,4], probs = 0.05/6) +0.0025
cv_SVM = quantile(all_test_results[,5], probs = 0.05/6) +0.006

# get critical values (0.05)
cv_rpart_0.5 = quantile(all_test_results[,1], probs = 0.05) -0.001
cv_knn_0.5 = quantile(all_test_results[,2], probs = 0.05) +0.008
cv_ranger_0.5 = quantile(all_test_results[,3], probs = 0.05) #+0.0005
cv_GLMnet_0.5 = quantile(all_test_results[,4], probs = 0.05) -0.008
cv_SVM_0.5 = quantile(all_test_results[,5], probs = 0.05) #+0.006


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

# # nice colors
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# # for geomsegment
d_observed = unlist(d_observed) 
d_observed_geom = data.frame(test = seq(1,5), d = d_observed )
legend_title = latex2exp::TeX("             ")
df = stack(all_test_results )
labels = fct_rev(df$ind) %>% levels # make sure labels align with levels in factors of df
labels_segments = fct_rev(df$ind)
labels = c("CART", "kNN", "RF","GLMNet", "SVM")
fill_rej_reg = 0.2

# visualize test statistics including observed ones (see Figure 6 in paper)
figure_6 = ggplot(df, aes(x = values, y = fct_rev(ind), fill = factor(stat(quantile)))) +
  xlim(c(-0.2,0.05)) +
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
  scale_y_discrete( name = "CRE vs.                             ", labels = labels) +
  coord_cartesian(clip = "off") +
  labs(title = legend_title) +
  theme_ridges(font_size = 18, grid = TRUE) +
  ggthemes::theme_economist_white(gray_bg = T) +
  theme(axis.title.x = element_blank()) +
  theme(legend.position = "bottom", legend.key.size = unit(1.2, 'cm'),
        legend.box.spacing = unit(1.8, 'cm'), 
        legend.text = element_text(size = 16),
        legend.title = element_text(size=21))+
  guides(fill = guide_legend(override.aes = list(size = 14))) +
  theme(plot.background = element_rect(fill = 'white', colour = 'white')) +
  theme(axis.text=element_text(size=13), 
        axis.title=element_text(size=28), 
        plot.title = element_text(size = 32, face = "bold")) +
  theme(axis.title.y = element_text(size = 20, hjust= 0.7, margin = margin(t = 20, r = 20, b = 0, l = 0))) +
  geom_segment(data = d_observed_geom, aes(x = d, xend = d, y = test, yend = test +0.8),
               color = "black", inherit.aes = FALSE, lineend = "round", linejoin = "round", size = 1.4) +
  annotate("rect", xmin=-0.2, xmax=cv_rpart, ymin=1, ymax=1.8,alpha=fill_rej_reg, fill = "red4") +
  annotate("rect", xmin=-0.2, xmax=cv_knn, ymin=2, ymax=2.8,alpha=fill_rej_reg, fill = "red4") +
  annotate("rect", xmin=-0.2, xmax=cv_ranger, ymin=3, ymax=3.8,alpha=fill_rej_reg, fill = "red4") +    
  annotate("rect", xmin=-0.2, xmax=cv_GLMnet, ymin=4, ymax=4.8,alpha=fill_rej_reg, fill = "red4") +
  annotate("rect", xmin=-0.2, xmax=cv_SVM, ymin=5, ymax=5.8,alpha=fill_rej_reg, fill = "red4") +
  annotate("rect", xmin=cv_rpart, xmax=cv_rpart_0.5, ymin=1, ymax=1.8,alpha=fill_rej_reg, fill = "#FF0000A0") +
  annotate("rect", xmin=cv_knn, xmax=cv_knn_0.5, ymin=2, ymax=2.8,alpha=fill_rej_reg, fill = "#FF0000A0") +
  annotate("rect", xmin=cv_ranger, xmax=cv_ranger_0.5, ymin=3, ymax=3.8,alpha=fill_rej_reg, fill = "#FF0000A0") +    
  annotate("rect", xmin=cv_GLMnet, xmax=cv_GLMnet_0.5, ymin=4, ymax=4.8,alpha=fill_rej_reg, fill = "#FF0000A0") +
  annotate("rect", xmin=cv_SVM, xmax=cv_SVM_0.5, ymin=5, ymax=5.8,alpha=fill_rej_reg, fill = "#FF0000A0") +
  annotate("rect", xmin=cv_rpart_0.5, xmax=0.05, ymin=1, ymax=1.8,alpha=0.2, fill = "steelblue") +
  annotate("rect", xmin=cv_knn_0.5, xmax=0.05, ymin=2, ymax=2.8,alpha=fill_rej_reg, fill = "steelblue") +
  annotate("rect", xmin=cv_ranger_0.5, xmax=0.05, ymin=3, ymax=3.8,alpha=fill_rej_reg, fill = "steelblue") +    
  annotate("rect", xmin=cv_GLMnet_0.5, xmax=0.05, ymin=4, ymax=4.8,alpha=fill_rej_reg, fill = "steelblue") +
  annotate("rect", xmin=cv_SVM_0.5, xmax=0.05, ymin=5, ymax=5.8,alpha=fill_rej_reg, fill = "steelblue") 

pdf(file= paste0("fig_6.pdf"), width = 13, height = 9)
print(figure_6)
dev.off()



## visualize CDF, see figure 7 in paper
df_cdf = df
levels(df_cdf$ind) = c("SVM", "GLMnet", "RF", "kNN","CART")

figure_7 = ggplot(df_cdf, aes(values, colour = ind)) +
  stat_ecdf(size=1.22) +
  theme(plot.title = element_text(hjust = 1.8, margin = margin(t = 0, r = 20, b = 0, l = 0))) +
  ggtitle( "CDFs of resampled test statistics (PMLB)") +
  ggtitle( "") +
  theme(axis.title.x = element_blank()) +
  xlim(c(-0.15,0))+
  theme(legend.position = "bottom", legend.key.size = unit(1, 'cm')) +
  theme(axis.text=element_text(size=13), axis.title=element_text(size=28), plot.title = element_text(size = 20, face = "bold")) +
  theme(axis.title.y = element_text(size = 20)) +
  ggthemes::theme_economist_white(gray_bg = T) +
  labs(color = "Test: CRE vs.") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.position = "bottom", legend.key.size = unit(1.2, 'cm'),
        legend.box.spacing = unit(1.8, 'cm'), 
        legend.text = element_text(size = 16),
        legend.title = element_text(size=21),
        plot.title = element_text(size = 28, face = "bold",hjust = 0.1, vjust = 2, margin = margin(t = 0, r = 20, b = 0, l = 10)))


pdf(file= paste0("fig_7.pdf"), width = 13, height = 9)
print(figure_7)
dev.off()




# function to compute shares of rejected resampled test statistics
reject_share <- function(gamma, d_res, d_obs){
  d_obs <- rep(d_obs, length(d_res))
  sum(d_res - d_obs > 2*gamma/(1 - gamma))
}
# define grid of gamma values...
gamma_grid = seq(0,0.2, by = 0.0001) %>% as.list()
#... and lapply share computation over this grid
rejection_shares = lapply(gamma_grid,
                          FUN = reject_share,
                          d_res = all_test_results$SVM,
                          d_obs = d_observed["SVM"] )
rej_shares_reg0 =  rejection_shares %>% unlist()
df_rej_SVM = data.frame("gamma" = gamma_grid %>% unlist,
                             "Rejection share" = rej_shares_reg0,
                             "test" = "CRE vs. SVM")

rejection_shares = lapply(gamma_grid,
                          FUN = reject_share,
                          d_res = all_test_results$GLMnet,
                          d_obs = d_observed["GLMnet"] )
rej_shares_reg0 =  rejection_shares %>% unlist()
df_rej_glmnet = data.frame("gamma" = gamma_grid %>% unlist,
                           "Rejection share" = rej_shares_reg0,
                           "test" = "CRE vs. GLMNet")

rejection_shares = lapply(gamma_grid,
                          FUN = reject_share,
                          d_res = all_test_results$ranger,
                          d_obs = d_observed["ranger"] )
rej_shares_reg0 =  rejection_shares %>% unlist()
df_rej_ranger = data.frame("gamma" = gamma_grid %>% unlist,
                           "Rejection share" = rej_shares_reg0,
                           "test" = "CRE vs. RF")


rejection_shares = lapply(gamma_grid,
                          FUN = reject_share,
                          d_res = all_test_results$knn,
                          d_obs = d_observed["knn"] )
rej_shares_reg0 =  rejection_shares %>% unlist()
df_rej_knn = data.frame("gamma" = gamma_grid %>% unlist,
                        "Rejection share" = rej_shares_reg0,
                        "test" = "CRE vs. kNN")

rejection_shares = lapply(gamma_grid,
                          FUN = reject_share,
                          d_res = all_test_results$rpart,
                          d_obs = d_observed["rpart"] )
rej_shares_reg0 =  rejection_shares %>% unlist()
df_rej_rpart = data.frame("gamma" = gamma_grid %>% unlist,
                          "Rejection share" = rej_shares_reg0,
                          "test" = "CRE vs. CART")



################################################################################
# Visualization of rejection shares
################################################################################

# only visualize unregularized test statistics, see footnote 7 in paper
df_rej_all = rbind(df_rej_glmnet, df_rej_knn, df_rej_ranger,
                   df_rej_SVM, df_rej_rpart)

par(mar=c(3,4,2,2))

# color grid for regulaizations
pal = brewer.pal(5, "YlGnBu")[c(2,3,4,5,6,7)]
# to do: discrete palette
pal = c("#DB6D00","#490092", "#56B4E9", "#009E73","#F0E442")
# compute p-values
df_rej_all$Rejection_share = 1- df_rej_all$Rejection.share/nrow(all_test_results)
# compute number of contaminations
number_resamples = length(all_test_results$SVM)
number_datasets = 62
df_rej_all$k = df_rej_all$gamma * number_datasets
#axis limits
x_limits = c(0,5)
y_limits = c(0,1)



# Eventually make plot of p-values (1- rejection shares) as function of contamination parameter gamma
# (corresponds to figure 8 in paper)
figure_8 = ggplot(data = df_rej_all) +
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
  scale_color_manual(values=pal)+
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1)) +
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1)) +
  theme(axis.title.x = element_text(margin = margin(t = 14, r = 20, b = 0, l = 0))) +
  theme(axis.title.y = element_text(margin = margin(t = 14, r = 20, b = 0, l = 0))) +
  labs(title = latex2exp::TeX("                   "))


pdf(file= paste0("fig_8.pdf"), width = 13, height = 9)
print(figure_8)
dev.off()


}

