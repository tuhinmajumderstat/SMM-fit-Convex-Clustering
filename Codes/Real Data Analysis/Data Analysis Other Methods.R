library(seqinr)
library(cvxclustr)
library(gtools)
library(dplyr)
library(Matrix)
library(parallel)

# import this functions from the location you store this
source("SMM_Fit_Other_Methods_Functions.R")

find_index <- function(input,base)
{
  m <- length(input)
  d <- base
  ind <- sum((input-1)* d^c(0:(m-1)))+1
  ind
}

log_trunc <- function(x)
{
  if(x>0)
  {
    return (log(x))
  }else{
    return (0)
  }
}

# change the working directory accordingly
setwd("/cwork/tm389/JCGS_Revision")

dat_sars <- read.fasta(file="./Data/Reference  Genes/SARS_Cov2_Ref.fasta")
dat_mers <- read.fasta(file="./Data/Reference  Genes/MERS_Ref.fasta")
dat_dengue <- read.fasta(file="./Data/Reference  Genes/Dengue_Ref.fasta")
dat_hepatitis <- read.fasta(file="./Data/Reference  Genes/Hepatitis_Ref.fasta")


####### Get Test Data #####


test_dat_sars_set_1 <- readRDS("./Data/Complete Test Genes/SARS-Cov-2_Test_Sample_Set_1.rds")
test_dat_sars_set_2 <- readRDS("./Data/Complete Test Genes/SARS-Cov-2_Test_Sample_Set_2.rds")
test_dat_sars_set_3 <- readRDS("./Data/Complete Test Genes/SARS-Cov-2_Test_Sample_Set_3.rds")
test_dat_sars_set_4 <- readRDS("./Data/Complete Test Genes/SARS-Cov-2_Test_Sample_Set_4.rds")

test_dat_mers <- readRDS("./Data/Complete Test Genes/MERS_Test_Sample_Set.rds")

test_dat_dengue <- readRDS("./Data/Complete Test Genes/Dengue_Test_Sample_Set.rds")

test_dat_hepatitis <- readRDS("./Data/Complete Test Genes/Hepatitis_Test_Sample_Set.rds")

fit_reference_model_other_methods <- function(dat,method="dpmm",m=3)
{
  dat_num <- matrix(0,ncol=length(dat[[1]]),nrow=1)
  for(i in 1:length(dat[[1]]))
  {
    if(dat[[1]][i]=="a")
    {
      dat_num[1,i] <- 1
    }else{
      if(dat[[1]][i]=="g")
      {
        dat_num[1,i] <- 2
      }else{
        if(dat[[1]][i]=="t")
        {
          dat_num[1,i] <- 3
        }else{
          dat_num[1,i] <- 4
        }
      }
    }
    
    
  }
  
  
  counts <- reformat(dat_num,m,c("1","2","3","4"))
  N <- do.call(rbind,counts$counts)
  
  dat_num_string <- apply(dat_num, 1, function(x) paste(x, collapse = ""))
  
  
  if(method=="dpmm")
  {
    fitted_model = run_dpmm(dat_num_string, 1, m, c("1","2","3","4"), 500, .1)
    fitted_model=fitted_model$output_matrix
  }
  if(method=="xiong")
  {
    fitted_model=run_xiong(dat_num_string, 1, m, c("1","2","3","4")) 
  }
  if(method=="ggl")
  {
    fitted_model <- run_ggl(dat_num_string, 1, m, c("1","2","3","4"))
  }
  if(method=="kmeans")
  {
    fitted_model <- run_kmeans(dat_num_string, 1, m, c("1","2","3","4")) 
  }
  names(fitted_model)= rownames(N)
  
  
  
  all_comb <- expand.grid(rep(list(1:d), m))
  #all_comb_mat <- t(sapply(all_comb, function(x) as.integer(strsplit(x, "")[[1]])))
  all_comb <- apply(all_comb, 1, function(x) paste0(x, collapse = ""))
  N <- N[all_comb,]
  N_data <- as.data.frame(N)
  names(N_data) <- paste0("state",c(1:d))
  
  N_clustered <- N_data
  N_clustered$cluster <- as.vector(fitted_model) 
  
  N_clustered <- N_clustered %>% group_by(cluster) %>% summarise_all(list(sum))
  N_clustered_counts <- as.data.frame(N_clustered %>% select(-cluster))
  R_hat <- N_clustered_counts / rowSums(N_clustered_counts)
  stationary_prob <- rowSums(N)/sum(N)
  
  
  fitted_model=fitted_model[all_comb]
  
  return(list("method"=method,"m"=m,"clusters"=fitted_model,"transition matrix"=R_hat,"stationary prob"=stationary_prob))
  
}

find_index <- function(input,base)
{
  m <- length(input)
  d <- base
  ind <- sum((input-1)* d^c(0:(m-1)))+1
  ind
}

log_trunc <- function(x)
{
  if(x>0)
  {
    return (log(x))
  }else{
    return (0)
  }
}



classify_virus_other_methods <- function(data,virus=c("SARS-Cov-2","MERS","Dengue","Hepatitis"),
                                         reference_model=reference_model_1,markov_order=c(4,4,3,3),
                                         sample_fraction=0.25)
{
  sample_size <- length(data)
  virus_number <- length(virus)
  log_lik_data <- matrix(0,nrow=sample_size,ncol=virus_number)
  for(sample_id in 1:sample_size)
  {
    dat <- data[[sample_id]]
    dat_num <- matrix(0,nrow=length(dat),ncol=1)
    for(i in 1:length(dat))
    {
      if(dat[i]=="a")
      {
        dat_num[i] <- 1
      }else{
        if(dat[i]=="g")
        {
          dat_num[i] <- 2
        }else{
          if(dat[i]=="t")
          {
            dat_num[i] <- 3
          }else{
            dat_num[i] <- 4
          }
        }
      }
    }
    
    seq_length <- floor(length(dat)*sample_fraction)
    set.seed(12345+sample_id)
    seq_start <- sample(c(1:(length(dat)-seq_length)),1)
    test_chain <- dat_num[seq_start:(seq_start+seq_length-1)]
    n <- length(test_chain)
    
    
    for(mdl in 1:length(virus))
    {
      model <- reference_model[[mdl]]
      d <- 4
      m <- markov_order[mdl]
      p <- d^m
      all_comb <- expand.grid(rep(list(1:d), m))
      
      N <- matrix(0,p,d,byrow=TRUE)
      for(i in (m+1):n)
      {
        prev <- test_chain[(i-m):(i-1)]
        j <- test_chain[i]
        tuple_ind <- find_index(prev,d)
        N[tuple_ind,j] <- N[tuple_ind,j]+1
      }
      
      N_data <- as.data.frame(N)
      names(N_data) <- paste0("state",c(1:d))
      cluster_id <- unique(model$cluster)
      N_clustered <- N_data %>% mutate(cluster=model$cluster)
      N_clustered <- N_clustered %>% group_by(cluster) %>% summarise_all(list(sum))
      N_clustered <- as.data.frame(N_clustered %>% select(-cluster))
      R_hat <- model$`transition matrix`
      log_R_hat <- matrix(0,nrow(R_hat),d,byrow=TRUE)
      for(i in 1:nrow(R_hat))
      {
        for(l in 1:d)
        {
          if(R_hat[i,l]==0)
          {
            log_R_hat[i,l]=0
          }else{
            log_R_hat[i,l] = log(R_hat[i,l])
          }
        }
      }
      log_lik_data[sample_id,mdl] <-  as.numeric(sum(as.matrix(N_clustered)*log_R_hat))
    }
    #print(sample_id)
  }
  
  predicted_virus <- sapply(c(1:sample_size), function(i)virus[which.max(log_lik_data[i,])]) 
  predicted_virus
  
}


#### Reference Model 1 ####

for(method in c("dpmm","xiong","ggl","kmeans"))
{
  fit_reference_sars_m_4 <- fit_reference_model_other_methods(dat_sars,method=method ,m=4)
  fit_reference_mers_m_4 <- fit_reference_model_other_methods(dat_mers,method=method,m=4)
  fit_reference_dengue_m_3 <- fit_reference_model_other_methods(dat_dengue,method=method,m=3)
  fit_reference_hepatitis_m_3 <- fit_reference_model_other_methods(dat_hepatitis,method=method,m=3)
  
  reference_model_1 <- list("SARS-Cov-2"=fit_reference_sars_m_4,"MERS"=fit_reference_mers_m_4,
                            "Dengue"=fit_reference_dengue_m_3,"Hepatitis"=fit_reference_hepatitis_m_3)
  
  
  ####################################################################
  ###################### Classification ##############################
  ####################################################################
  
  
  
  result_sars_1_percent_25 <- classify_virus_other_methods(data=test_dat_sars_set_1,reference_model=reference_model_1,markov_order=c(4,4,3,3),sample_fraction = 0.25)
  result_sars_2_percent_25 <- classify_virus_other_methods(data=test_dat_sars_set_2,reference_model=reference_model_1,markov_order=c(4,4,3,3),sample_fraction = 0.25)
  result_sars_3_percent_25 <- classify_virus_other_methods(data=test_dat_sars_set_3,reference_model=reference_model_1,markov_order=c(4,4,3,3),sample_fraction = 0.25)
  result_sars_4_percent_25 <- classify_virus_other_methods(data=test_dat_sars_set_4,reference_model=reference_model_1,markov_order=c(4,4,3,3),sample_fraction = 0.25)
  result_mers_percent_25 <- classify_virus_other_methods(data=test_dat_mers,reference_model=reference_model_1,markov_order=c(4,4,3,3),sample_fraction = 0.25)
  result_dengue_percent_25 <- classify_virus_other_methods(data=test_dat_dengue,reference_model=reference_model_1,markov_order=c(4,4,3,3),sample_fraction = 0.25)
  result_hepatitis_percent_25 <- classify_virus_other_methods(data=test_dat_hepatitis,reference_model=reference_model_1,markov_order=c(4,4,3,3),sample_fraction = 0.25)
  
  result_all_25 <- c(result_sars_1_percent_25,result_sars_2_percent_25,result_sars_3_percent_25,
                     result_sars_4_percent_25,result_mers_percent_25,result_dengue_percent_25,
                     result_hepatitis_percent_25)
  
  Confusion_MX_25 <- matrix(0,ncol=length(virus),nrow=length(virus))
  for(iter in 1:length(true_virus_label))
  {
    i <- which(virus==true_virus_label[iter])
    j <- which(virus==result_all_25[iter])
    Confusion_MX_25[i,j] <- Confusion_MX_25[i,j]+1
  }
  
  Confusion_MX_25
  
  result_sars_1_percent_10 <- classify_virus_other_methods(data=test_dat_sars_set_1,reference_model=reference_model_1,markov_order=c(4,4,3,3),sample_fraction = 0.10)
  result_sars_2_percent_10 <- classify_virus_other_methods(data=test_dat_sars_set_2,reference_model=reference_model_1,markov_order=c(4,4,3,3),sample_fraction = 0.10)
  result_sars_3_percent_10 <- classify_virus_other_methods(data=test_dat_sars_set_3,reference_model=reference_model_1,markov_order=c(4,4,3,3),sample_fraction = 0.10)
  result_sars_4_percent_10 <- classify_virus_other_methods(data=test_dat_sars_set_4,reference_model=reference_model_1,markov_order=c(4,4,3,3),sample_fraction = 0.10)
  result_mers_percent_10 <- classify_virus_other_methods(data=test_dat_mers,reference_model=reference_model_1,markov_order=c(4,4,3,3),sample_fraction = 0.10)
  result_dengue_percent_10 <- classify_virus_other_methods(data=test_dat_dengue,reference_model=reference_model_1,markov_order=c(4,4,3,3),sample_fraction = 0.10)
  result_hepatitis_percent_10 <- classify_virus_other_methods(data=test_dat_hepatitis,reference_model=reference_model_1,markov_order=c(4,4,3,3),sample_fraction = 0.10)
  
  result_all_10 <- c(result_sars_1_percent_10,result_sars_2_percent_10,result_sars_3_percent_10,
                     result_sars_4_percent_10,result_mers_percent_10,result_dengue_percent_10,
                     result_hepatitis_percent_10)
  
  Confusion_MX_10 <- matrix(0,ncol=length(virus),nrow=length(virus))
  for(iter in 1:length(true_virus_label))
  {
    i <- which(virus==true_virus_label[iter])
    j <- which(virus==result_all_10[iter])
    Confusion_MX_10[i,j] <- Confusion_MX_10[i,j]+1
  }
  
  
  result_sars_1_percent_05 <- classify_virus_other_methods(data=test_dat_sars_set_1,reference_model=reference_model_1,markov_order=c(4,4,3,3),sample_fraction = 0.05)
  result_sars_2_percent_05 <- classify_virus_other_methods(data=test_dat_sars_set_2,reference_model=reference_model_1,markov_order=c(4,4,3,3),sample_fraction = 0.05)
  result_sars_3_percent_05 <- classify_virus_other_methods(data=test_dat_sars_set_3,reference_model=reference_model_1,markov_order=c(4,4,3,3),sample_fraction = 0.05)
  result_sars_4_percent_05 <- classify_virus_other_methods(data=test_dat_sars_set_4,reference_model=reference_model_1,markov_order=c(4,4,3,3),sample_fraction = 0.05)
  result_mers_percent_05 <- classify_virus_other_methods(data=test_dat_mers,reference_model=reference_model_1,markov_order=c(4,4,3,3),sample_fraction = 0.05)
  result_dengue_percent_05 <- classify_virus_other_methods(data=test_dat_dengue,reference_model=reference_model_1,markov_order=c(4,4,3,3),sample_fraction = 0.05)
  result_hepatitis_percent_05 <- classify_virus_other_methods(data=test_dat_hepatitis,reference_model=reference_model_1,markov_order=c(4,4,3,3),sample_fraction = 0.05)
  
  result_all_05 <- c(result_sars_1_percent_05,result_sars_2_percent_05,result_sars_3_percent_05,
                     result_sars_4_percent_05,result_mers_percent_05,result_dengue_percent_05,
                     result_hepatitis_percent_05)
  
  
  Confusion_MX_05 <- matrix(0,ncol=length(virus),nrow=length(virus))
  for(iter in 1:length(true_virus_label))
  {
    i <- which(virus==true_virus_label[iter])
    j <- which(virus==result_all_05[iter])
    Confusion_MX_05[i,j] <- Confusion_MX_05[i,j]+1
  }
  
  
  Confusion_MX_05
  Confusion_MX_10
  Confusion_MX_25
  saveRDS(Confusion_MX_05,paste0("Confusion_MX_05_model_1_method_",method,".rds"))
  saveRDS(Confusion_MX_10,paste0("Confusion_MX_10_model_1_method_",method,".rds"))
  saveRDS(Confusion_MX_25,paste0("Confusion_MX_25_model_1_method_",method,".rds"))
  print(paste0("Method:",method))
  print(paste0("5% missclassification: ", 1-sum(diag(Confusion_MX_05))/sum(Confusion_MX_05)))
  print(paste0("10% missclassification: ", 1-sum(diag(Confusion_MX_10))/sum(Confusion_MX_10)))
  print(paste0("25% missclassification: ", 1-sum(diag(Confusion_MX_25))/sum(Confusion_MX_25)))
}




#### Reference Model 2 ####

for(method in c("dpmm","xiong","ggl","kmeans"))
{
  fit_reference_sars_m_3 <- fit_reference_model_other_methods(dat_sars,method=method ,m=3)
  fit_reference_mers_m_3 <- fit_reference_model_other_methods(dat_mers,method=method,m=3)
  fit_reference_dengue_m_3 <- fit_reference_model_other_methods(dat_dengue,method=method,m=3)
  fit_reference_hepatitis_m_3 <- fit_reference_model_other_methods(dat_hepatitis,method=method,m=3)
  
  reference_model_2 <- list("SARS-Cov-2"=fit_reference_sars_m_3,"MERS"=fit_reference_mers_m_3,
                            "Dengue"=fit_reference_dengue_m_3,"Hepatitis"=fit_reference_hepatitis_m_3)
  
  result_sars_1_percent_25 <- classify_virus_other_methods(data=test_dat_sars_set_1,reference_model=reference_model_2,markov_order=c(3,3,3,3),sample_fraction = 0.25)
  result_sars_2_percent_25 <- classify_virus_other_methods(data=test_dat_sars_set_2,reference_model=reference_model_2,markov_order=c(3,3,3,3),sample_fraction = 0.25)
  result_sars_3_percent_25 <- classify_virus_other_methods(data=test_dat_sars_set_3,reference_model=reference_model_2,markov_order=c(3,3,3,3),sample_fraction = 0.25)
  result_sars_4_percent_25 <- classify_virus_other_methods(data=test_dat_sars_set_4,reference_model=reference_model_2,markov_order=c(3,3,3,3),sample_fraction = 0.25)
  result_mers_percent_25 <- classify_virus_other_methods(data=test_dat_mers,reference_model=reference_model_2,markov_order=c(3,3,3,3),sample_fraction = 0.25)
  result_dengue_percent_25 <- classify_virus_other_methods(data=test_dat_dengue,reference_model=reference_model_2,markov_order=c(3,3,3,3),sample_fraction = 0.25)
  result_hepatitis_percent_25 <- classify_virus_other_methods(data=test_dat_hepatitis,reference_model=reference_model_2,markov_order=c(3,3,3,3),sample_fraction = 0.25)
  
  result_all_25 <- c(result_sars_1_percent_25,result_sars_2_percent_25,result_sars_3_percent_25,
                     result_sars_4_percent_25,result_mers_percent_25,result_dengue_percent_25,
                     result_hepatitis_percent_25)
  
  Confusion_MX_25 <- matrix(0,ncol=length(virus),nrow=length(virus))
  for(iter in 1:length(true_virus_label))
  {
    i <- which(virus==true_virus_label[iter])
    j <- which(virus==result_all_25[iter])
    Confusion_MX_25[i,j] <- Confusion_MX_25[i,j]+1
  }
  
  Confusion_MX_25
  
  result_sars_1_percent_10 <- classify_virus_other_methods(data=test_dat_sars_set_1,reference_model=reference_model_2,markov_order=c(3,3,3,3),sample_fraction = 0.10)
  result_sars_2_percent_10 <- classify_virus_other_methods(data=test_dat_sars_set_2,reference_model=reference_model_2,markov_order=c(3,3,3,3),sample_fraction = 0.10)
  result_sars_3_percent_10 <- classify_virus_other_methods(data=test_dat_sars_set_3,reference_model=reference_model_2,markov_order=c(3,3,3,3),sample_fraction = 0.10)
  result_sars_4_percent_10 <- classify_virus_other_methods(data=test_dat_sars_set_4,reference_model=reference_model_2,markov_order=c(3,3,3,3),sample_fraction = 0.10)
  result_mers_percent_10 <- classify_virus_other_methods(data=test_dat_mers,reference_model=reference_model_2,markov_order=c(3,3,3,3),sample_fraction = 0.10)
  result_dengue_percent_10 <- classify_virus_other_methods(data=test_dat_dengue,reference_model=reference_model_2,markov_order=c(3,3,3,3),sample_fraction = 0.10)
  result_hepatitis_percent_10 <- classify_virus_other_methods(data=test_dat_hepatitis,reference_model=reference_model_2,markov_order=c(3,3,3,3),sample_fraction = 0.10)
  
  result_all_10 <- c(result_sars_1_percent_10,result_sars_2_percent_10,result_sars_3_percent_10,
                     result_sars_4_percent_10,result_mers_percent_10,result_dengue_percent_10,
                     result_hepatitis_percent_10)
  
  Confusion_MX_10 <- matrix(0,ncol=length(virus),nrow=length(virus))
  for(iter in 1:length(true_virus_label))
  {
    i <- which(virus==true_virus_label[iter])
    j <- which(virus==result_all_10[iter])
    Confusion_MX_10[i,j] <- Confusion_MX_10[i,j]+1
  }
  
  
  result_sars_1_percent_05 <- classify_virus_other_methods(data=test_dat_sars_set_1,reference_model=reference_model_2,markov_order=c(3,3,3,3),sample_fraction = 0.05)
  result_sars_2_percent_05 <- classify_virus_other_methods(data=test_dat_sars_set_2,reference_model=reference_model_2,markov_order=c(3,3,3,3),sample_fraction = 0.05)
  result_sars_3_percent_05 <- classify_virus_other_methods(data=test_dat_sars_set_3,reference_model=reference_model_2,markov_order=c(3,3,3,3),sample_fraction = 0.05)
  result_sars_4_percent_05 <- classify_virus_other_methods(data=test_dat_sars_set_4,reference_model=reference_model_2,markov_order=c(3,3,3,3),sample_fraction = 0.05)
  result_mers_percent_05 <- classify_virus_other_methods(data=test_dat_mers,reference_model=reference_model_2,markov_order=c(3,3,3,3),sample_fraction = 0.05)
  result_dengue_percent_05 <- classify_virus_other_methods(data=test_dat_dengue,reference_model=reference_model_2,markov_order=c(3,3,3,3),sample_fraction = 0.05)
  result_hepatitis_percent_05 <- classify_virus_other_methods(data=test_dat_hepatitis,reference_model=reference_model_2,markov_order=c(3,3,3,3),sample_fraction = 0.05)
  
  result_all_05 <- c(result_sars_1_percent_05,result_sars_2_percent_05,result_sars_3_percent_05,
                     result_sars_4_percent_05,result_mers_percent_05,result_dengue_percent_05,
                     result_hepatitis_percent_05)
  
  
  Confusion_MX_05 <- matrix(0,ncol=length(virus),nrow=length(virus))
  for(iter in 1:length(true_virus_label))
  {
    i <- which(virus==true_virus_label[iter])
    j <- which(virus==result_all_05[iter])
    Confusion_MX_05[i,j] <- Confusion_MX_05[i,j]+1
  }
  
  
  Confusion_MX_05
  Confusion_MX_10
  Confusion_MX_25
  
  saveRDS(Confusion_MX_05,paste0("Confusion_MX_05_model_2_method_",method,".rds"))
  saveRDS(Confusion_MX_10,paste0("Confusion_MX_10_model_2_method_",method,".rds"))
  saveRDS(Confusion_MX_25,paste0("Confusion_MX_25_model_2_method_",method,".rds"))
  print(paste0("Method:",method))
  print(paste0("5% missclassification: ", 1-sum(diag(Confusion_MX_05))/sum(Confusion_MX_05)))
  print(paste0("10% missclassification: ", 1-sum(diag(Confusion_MX_10))/sum(Confusion_MX_10)))
  print(paste0("25% missclassification: ", 1-sum(diag(Confusion_MX_25))/sum(Confusion_MX_25)))
}








