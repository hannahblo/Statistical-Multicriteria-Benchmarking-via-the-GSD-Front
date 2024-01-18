
## needed own functions:


compute_k_folds <- function(n_row,k=10,set.seed=TRUE,seed=12345678){
  result <- rep((1:k),ceiling(n_row/k))
  result <- result[(1:n_row)]
  if(set.seed){set.seed(seed)}
  return(sample(result))
}


compute_accuracy <- function(predictions,true_labels){
  if(method == "cre"){
    preds_cre = ifelse(predictions > 0.5, 1,0) %>% as.factor
    result <- mean(preds_cre==true_labels)
  }else{
    result <- mean(predictions==true_labels)
  }
  return(result)
}


perturbate_y <- function(y,p_percent=20,set.seed=TRUE,seed=12345678){
  if(set.seed){set.seed(seed)}
  indexs <- which( sample(c(1,0),size=length(y),prob=c(p_percent/100,(100-p_percent)/100),replace=TRUE)==1)
  length_indexs <- length(indexs)
  y[indexs] <- sample(y, size=length_indexs)
  return(y)
}


perturbate_x <- function(x,p_percent=20,set.seed=TRUE,seed=12345678){
  if(set.seed){set.seed(seed)}
  
  for(k in seq_len(ncol(x))){
    indexs <- which( sample(c(1,0),size=nrow(x),prob=c(p_percent/100,(100-p_percent)/100),replace=TRUE)==1)
    length_indexs <- length(indexs)
    x[indexs,k] <- sample(x[,k], size=length_indexs)
  }
  return(x)
}
########