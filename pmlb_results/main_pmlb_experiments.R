####
# This script runs 10-fold cross-validation on each considered classifier on 
# each of the 62 datasets from the PMLB benchmark suite. It returns i) classical 
# accuracies, ii) accuracies with perturbed x and iii) accuracies with perturbed y.
#####

#compared classifiers:
#SVM
#RPART
#Random Forest
#ELASTIC NET
#KNN

# source helper functions
source("pmlb_results/helper_pmlb_experiments.R")

# for compressed rule ensemble classifier (CRE)
source("pmlb_results/make_method_cre.R")
library(cre)

library(pmlbr)
library(caret)
library(RWeka)
library(ranger)
library(kknn)
library(RSNNS)
library(pmlbr)
library(e1071)
library(pROC)
library(randomForest)
library(kernlab)
library(ada)
library(gbm)
library(dplyr)
t <- 1
datasets <- list()
for(k in seq_len(length(classification_dataset_names))){
  index <- which(summary_stats[,1]==classification_dataset_names[k])
  metadat <- summary_stats[index,]
  if( metadat$n_instances %in% c(40:1000) &
  metadat$n_classes==2 &
  metadat$task== "classification" &
  metadat$n_features <= 100){
      dat <- fetch_data(classification_dataset_names[k])
      datasets[[t]] <- dat
      print(t)
      t <- t+1
  }
}
# optional: save downloaded data locally:
# saveRDS(datasets,"datasets.RDS")



# compared methods:
#SVM
#RPART
#Random Forest
#ELASTIC NET
#KNN
methods <- c("cre","svmRadial","J48","ranger","knn","glmnet")
results <- array(0,c(length(datasets),3*length(methods)))
colnames(results) <- rep(methods,3)

# list for accuracy matrices in each loop
clean_accuracies_list = list()
accuracies_noisy_y_list = list()
accuracies_noisy_x_list = list()

# outer loop over datasets
for(k in seq_len(length(datasets))){
  print(Sys.time())
  
  dat <- datasets[[k]]
  dim(dat)
  n_row <- nrow(dat)
  n_col <- ncol(dat)
  target_index = which(colnames(dat)== "target")
  x <- dat[,-target_index]
  y <- dat[,target_index]
  if(colnames(dat)[n_col]!="target"){print("error")}
  y_noisy <- perturbate_y(y)
  x_noisy <- perturbate_x(x)
  folds <- compute_k_folds(n_row,k=10)
  
  
  clean_accuracies <- array(0,c(10,length(methods)))
  colnames(clean_accuracies) <- methods
  accuracies_noisy_y <- clean_accuracies
  accuracies_noisy_x <- clean_accuracies
  
  
  # hyperparameter tuning
  # sets hyperparams for each dataset
  best_hyperparams <- list()
  try(
  for(i in seq_len(length(methods))){
    method <- methods[i]
    print(method)
    if(method == "cre")
      method = cre_method
    model <- caret::train(x=x,y=as.factor(y),method = method,
                          tuneLength = 10,
                          trControl = trainControl(number= 5))
    best_hyperparams[[i]] = model$bestTune
  }
)
  
  
  print("##############")
  print("k:")
  print(k)
  print("##############")
  print(max(results[k,]))
  for(l in (1:10)){
    indexs <- which(folds==l)
    x_train <- x[-indexs,]
    y_train <- y[-indexs]
    x_test <- x[indexs,]
    y_test <- y[indexs]
    
    x_train_noisy <- x_noisy[-indexs,]
    y_train_noisy <- y_noisy[-indexs]
    x_test_noisy <- x_noisy[indexs,]
    y_test_noisy <- y_noisy[indexs]
    
    
    
    ### clean data
    for(i in seq_len(length(methods))){
      method <- methods[i]
      try(
      if(method == "cre"){
        model = cre(x=x_train,y=as.factor(y_train), task = "class")
        accuracy <- compute_accuracy(predict(model,x_test),y_test)
      }
      else{
        model <- caret::train(x=x_train,y=as.factor(y_train),method = method, 
                              tuneGrid = best_hyperparams[[i]],
                              trControl = trainControl(number= 1))
        accuracy <- compute_accuracy(predict(model,x_test),y_test)
      }
      )
      clean_accuracies[l,i] <- accuracy
      
    }
    
    ### perturbation in y
    for(i in seq_len(length(methods))){
      method <- methods[i]
      try(
      if(method == "cre"){
        model = cre(x=x_train,y=as.factor(y_train_noisy), task = "class")
        #model$mat_names = names(x_test)
        accuracy <- compute_accuracy(predict(model,x_test),y_test_noisy)
      }
      else{
        model <- caret::train(x=x_train,y=as.factor(y_train_noisy),method = method, 
                              tuneGrid = best_hyperparams[[i]],
                              trControl = trainControl(number= 1))
        accuracy <- compute_accuracy(predict(model,x_test),y_test_noisy)
      }
      )
      accuracies_noisy_y[l,i] <- accuracy
      #print(c(method,accuracy))
    }

    ### perturbation in x
    for(i in seq_len(length(methods))){
      method <- methods[i]
      try(
      if(method == "cre"){
        model = cre(x=x_train_noisy,y=as.factor(y_train), task = "class")
        #model$mat_names = names(x_test_noisy)
        accuracy <- compute_accuracy(predict(model,x_test_noisy),y_test)
      }
      else{
        model <- caret::train(x=x_train_noisy,y=as.factor(y_train),method = method, 
                              tuneGrid = best_hyperparams[[i]],
                              trControl = trainControl(number= 1))
        accuracy <- compute_accuracy(predict(model,x_test_noisy),y_test)
      }
      )
      accuracies_noisy_x[l,i] <- accuracy
      #print(c(method,accuracy))
    }
      
  }
  results[k,] <- (c(colMeans(clean_accuracies),colMeans(accuracies_noisy_y),colMeans(accuracies_noisy_x)))
  clean_accuracies_list[[k]] = clean_accuracies
  accuracies_noisy_y_list[[k]] = accuracies_noisy_y
  accuracies_noisy_x_list[[k]] = accuracies_noisy_x
}



# stack results:
results_clean = results[,1:6]
results_noisy_x = results[,7:12]
results_noisy_y = results[,13:8]
results_clean_stacked = results_clean %>% t() %>% as.vector()
names(results_clean_stacked) = rep(methods,nrow(results))
results_noisy_x_stacked = results_noisy_x %>% t() %>% as.vector()
names(results_noisy_x_stacked) = rep(methods,nrow(results))
results_noisy_y_stacked = results_noisy_y %>% t() %>% as.vector()
names(results_noisy_y_stacked) = rep(methods,nrow(results))

final_res = cbind(results_clean_stacked, results_noisy_x_stacked, results_noisy_y_stacked)

# optinal: save final results
#saveRDS(final_res,file = "final_results.RDS")


