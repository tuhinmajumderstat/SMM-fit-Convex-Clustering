library(seqinr)
library(cvxclustr)
library(gtools)
library(dplyr)
library(Matrix)
library(parallel)

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



## Fit the reference model

fit_reference_model <- function(dat,d=4,m=3,knn=5,phi=100,lambda_seq=seq(0,100,by=0.1))
{
  dat_num <- matrix(0,nrow=length(dat[[1]]),ncol=1)
  for(i in 1:length(dat[[1]]))
  {
    if(dat[[1]][i]=="a")
    {
      dat_num[i] <- 1
    }else{
      if(dat[[1]][i]=="g")
      {
        dat_num[i] <- 2
      }else{
        if(dat[[1]][i]=="t")
        {
          dat_num[i] <- 3
        }else{
          dat_num[i] <- 4
        }
      }
    }
    
    
  }
  
  
  chain <- dat_num
  n <- length(chain)
  p <- d^m
  all_comb <- expand.grid(rep(list(1:d), m))
  N <- matrix(0,p,d,byrow=TRUE)
  for(i in (m+1):n)
  {
    prev <- chain[(i-m):(i-1)]
    j <- chain[i]
    tuple_ind <- find_index(prev,d)
    N[tuple_ind,j] <- N[tuple_ind,j]+1
  }
  
  X <- N/rowSums(N)
  X <- t(X)
  keep_rows <- which(rowSums(N)>0)
  X <- X[,keep_rows]
  nu <- 1.8*(1/p)
  wt <- matrix(0,nrow = ncol(X),ncol=ncol(X))
  for(i in 1:ncol(X))
  {
    for(j in 1:ncol(X))
    {
      wt[i,j] <- max(abs(X[,i]-X[,j]))
    }
  }
  
  w <- array()
  start <- 1
  for(i in 1:(ncol(X)-1))
  {
    for(j in (i+1):ncol(X))
    {
      w[start] <- exp(-phi*wt[i,j])
      start <- start+1
    }
  }
  
  w <- knn_weights(w,k=knn,p)
  
  model_all <- cvxclust(X,w=w,gamma=lambda_seq,nu=nu,tol=1e-7)
  N_data <- as.data.frame(N[keep_rows,])
  names(N_data) <- paste0("state",c(1:d))
  BIC <- matrix(0,length(lambda_seq),1,byrow=TRUE)
  
  for(j in 1:length(lambda_seq))
  {
    A <- create_adjacency(model_all$V[[j]],w,ncol(X))
    A <- as(A,"dgCMatrix")
    assignment <- find_clusters(A)
    cluster_id <- unique(assignment$cluster)
    N_clustered <- N_data %>% mutate(cluster=assignment$cluster)
    N_clustered <- N_clustered %>% group_by(cluster) %>% summarise_all(list(sum))
    N_clustered <- as.data.frame(N_clustered %>% select(-cluster))
    R_hat <- N_clustered / rowSums(N_clustered)
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
    BIC[j] = as.numeric(-2*sum(as.matrix(N_clustered)*log_R_hat) + (nrow(R_hat))*(d-1)*log(n))
    print(j)
  }
  
  lambda_min <- lambda_seq[which.min(BIC)]
  #plot(lambda_seq,BIC,type="l")
  A <- create_adjacency(model_all$V[[which.min(BIC)]],w,ncol(X))
  opt_assignment <- find_clusters(A)
  unique_cluster_id <- unique(opt_assignment$cluster)
  opt_assignment
  
  N_clustered <- N_data %>% mutate(cluster=opt_assignment$cluster)
  N_clustered <- N_clustered %>% group_by(cluster) %>% summarise_all(list(sum))
  N_clustered_counts <- as.data.frame(N_clustered %>% select(-cluster))
  R_hat <- N_clustered_counts / rowSums(N_clustered_counts)
  stationary_prob <- rowSums(N)/sum(N)
  
  final_model <- list("cluster"=opt_assignment$cluster,"size"=opt_assignment$size,
                      "state space"=all_comb,"transition matrix"=R_hat,"stationary prob"=stationary_prob,
                      "lambda"=lambda_min)
  final_model
  
  
  
}


## change the working directory to which you store your datasets

setwd("/cwork/tm389/JCGS_Revision")
dat_sars <- read.fasta(file="./Data/Reference Genes/SARS_Cov2_Ref.fasta")
dat_mers <- read.fasta(file="./Data/Reference Genes/MERS_Ref.fasta")
dat_dengue <- read.fasta(file="./Data/Reference Genes/Dengue_Ref.fasta")
dat_hepatitis <- read.fasta(file="./Data/Reference Genes/Hepatitis_Ref.fasta")


#############################
######## m=3, knn=20 ########
#############################

sars_final_model_m_3_knn_20 <- fit_reference_model(dat=dat_sars,d=4,m=3,
                                                   knn=20,phi=100,lambda_seq=seq(0,20,by=0.01))

mers_final_model_m_3_knn_20 <- fit_reference_model(dat=dat_mers,d=4,m=3,
                                                   knn=20,phi=100,lambda_seq=seq(0,40,by=0.01))

dengue_final_model_m_3_knn_20 <- fit_reference_model(dat=dat_dengue,d=4,m=3,
                                                     knn=20,phi=100,lambda_seq=c(seq(0,100,by=0.1),seq(101,200,by=1)))

hepatitis_final_model_m_3_knn_20 <- fit_reference_model(dat=dat_hepatitis,d=4,m=3,
                                                        knn=20,phi=100,lambda_seq=c(seq(0,100,by=0.1),seq(110,1000,by=5)))



#############################
######## m=3, knn=5 ########
#############################

sars_final_model_m_3_knn_5 <- fit_reference_model(dat=dat_sars,d=4,m=3,
                                                  knn=5,phi=100,lambda_seq=seq(0,20,by=0.01))

mers_final_model_m_3_knn_5 <- fit_reference_model(dat=dat_mers,d=4,m=3,
                                                  knn=5,phi=100,lambda_seq=seq(0,40,by=0.01))

dengue_final_model_m_3_knn_5 <- fit_reference_model(dat=dat_dengue,d=4,m=3,
                                                    knn=5,phi=100,lambda_seq=c(seq(0,100,by=0.1),seq(101,200,by=1)))

hepatitis_final_model_m_3_knn_5 <- fit_reference_model(dat=dat_hepatitis,d=4,m=3,
                                                       knn=5,phi=100,lambda_seq=c(seq(0,100,by=0.1),seq(110,1000,by=5)))






#############################
######## m=4, knn=50 ########
#############################

sars_final_model_m_4_knn_50 <- fit_reference_model(dat=dat_sars,d=4,m=4,
                                                   knn=50,phi=100,lambda_seq=seq(0,20,by=0.01))

mers_final_model_m_4_knn_50 <- fit_reference_model(dat=dat_mers,d=4,m=4,
                                                   knn=50,phi=100,lambda_seq=seq(0,40,by=0.01))

dengue_final_model_m_4_knn_50 <- fit_reference_model(dat=dat_dengue,d=4,m=4,
                                                     knn=50,phi=100,lambda_seq=c(seq(0,100,by=0.1),seq(101,200,by=1)))




#############################
######## m=4, knn=20 ########
#############################

sars_final_model_m_4_knn_20 <- fit_reference_model(dat=dat_sars,d=4,m=4,
                                                   knn=20,phi=100,lambda_seq=seq(0,20,by=0.01))

mers_final_model_m_4_knn_20 <- fit_reference_model(dat=dat_mers,d=4,m=4,
                                                   knn=20,phi=100,lambda_seq=seq(0,40,by=0.01))

dengue_final_model_m_4_knn_20 <- fit_reference_model(dat=dat_dengue,d=4,m=4,
                                                     knn=20,phi=100,lambda_seq=c(seq(0,100,by=0.1),seq(101,200,by=1)))



#######################################################
#######################################################
#######################################################




reference_model_1 <- list("SARS-Cov-2"=sars_final_model_m_4_knn_20,"MERS"=mers_final_model_m_4_knn_20,
                          "Dengue"=dengue_final_model_m_3_knn_5,"Hepatitis"=hepatitis_final_model_m_3_knn_5)
saveRDS(reference_model_1,file="reference_smm_model_1.rds")


reference_model_2 <- list("SARS-Cov-2"=sars_final_model_m_3_knn_20,"MERS"=mers_final_model_m_3_knn_20,
                          "Dengue"=dengue_final_model_m_3_knn_20,"Hepatitis"=hepatitis_final_model_m_3_knn_20)
saveRDS(reference_model_2,file="reference_smm_model_2.rds")

reference_model_3 <- list("SARS-Cov-2"=sars_final_model_m_4_knn_50,"MERS"=mers_final_model_m_4_knn_50,
                          "Dengue"=dengue_final_model_m_3_knn_20,"Hepatitis"=hepatitis_final_model_m_3_knn_20)
saveRDS(reference_model_3,file="reference_smm_model_3.rds")

reference_model_4 <- list("SARS-Cov-2"=sars_final_model_m_3_knn_5,"MERS"=mers_final_model_m_3_knn_5,
                          "Dengue"=dengue_final_model_m_3_knn_5,"Hepatitis"=hepatitis_final_model_m_3_knn_5)
saveRDS(reference_model_4,file="reference_smm_model_4.rds")

reference_model_5 <- list("SARS-Cov-2"=sars_final_model_m_4_knn_20,"MERS"=mers_final_model_m_4_knn_20,
                          "Dengue"=dengue_final_model_m_4_knn_20,"Hepatitis"=hepatitis_final_model_m_3_knn_5)
saveRDS(reference_model_5,file="reference_smm_model_5.rds")




####################################################################
###################### Classification ##############################
####################################################################


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



classify_virus <- function(data,virus=c("SARS-Cov-2","MERS","Dengue","Hepatitis"),
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


####### Get Test Data #####

## Change the location based on your data storage
test_dat_sars_set_1 <- readRDS("./Data/Complete Test Genes/SARS-Cov-2_Test_Sample_Set_1.rds")
test_dat_sars_set_2 <- readRDS("./Data/Complete Test Genes/SARS-Cov-2_Test_Sample_Set_2.rds")
test_dat_sars_set_3 <- readRDS("./Data/Complete Test Genes/SARS-Cov-2_Test_Sample_Set_3.rds")
test_dat_sars_set_4 <- readRDS("./Data/Complete Test Genes/SARS-Cov-2_Test_Sample_Set_4.rds")

test_dat_mers <- readRDS("./Data/Complete Test Genes/MERS_Test_Sample_Set.rds")

test_dat_dengue <- readRDS("./Data/Complete Test Genes/Dengue_Test_Sample_Set.rds")

test_dat_hepatitis <- readRDS("./Data/Complete Test Genes/Hepatitis_Test_Sample_Set.rds")


###### Result: Model1 ######

virus <- c("SARS-Cov-2","MERS","Dengue","Hepatitis")
true_virus_label <- c(rep("SARS-Cov-2",200), rep("MERS",50),rep("Dengue",100),rep("Hepatitis",150))



result_sars_1_percent_25 <- classify_virus(data=test_dat_sars_set_1,reference_model=reference_model_1,markov_order=c(4,4,3,3),sample_fraction = 0.25)
result_sars_2_percent_25 <- classify_virus(data=test_dat_sars_set_2,reference_model=reference_model_1,markov_order=c(4,4,3,3),sample_fraction = 0.25)
result_sars_3_percent_25 <- classify_virus(data=test_dat_sars_set_3,reference_model=reference_model_1,markov_order=c(4,4,3,3),sample_fraction = 0.25)
result_sars_4_percent_25 <- classify_virus(data=test_dat_sars_set_4,reference_model=reference_model_1,markov_order=c(4,4,3,3),sample_fraction = 0.25)
result_mers_percent_25 <- classify_virus(data=test_dat_mers,reference_model=reference_model_1,markov_order=c(4,4,3,3),sample_fraction = 0.25)
result_dengue_percent_25 <- classify_virus(data=test_dat_dengue,reference_model=reference_model_1,markov_order=c(4,4,3,3),sample_fraction = 0.25)
result_hepatitis_percent_25 <- classify_virus(data=test_dat_hepatitis,reference_model=reference_model_1,markov_order=c(4,4,3,3),sample_fraction = 0.25)

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


result_sars_1_percent_10 <- classify_virus(data=test_dat_sars_set_1,reference_model=reference_model_1,markov_order=c(4,4,3,3),sample_fraction = 0.1)
result_sars_2_percent_10 <- classify_virus(data=test_dat_sars_set_2,reference_model=reference_model_1,markov_order=c(4,4,3,3),sample_fraction = 0.1)
result_sars_3_percent_10 <- classify_virus(data=test_dat_sars_set_3,reference_model=reference_model_1,markov_order=c(4,4,3,3),sample_fraction = 0.1)
result_sars_4_percent_10 <- classify_virus(data=test_dat_sars_set_4,reference_model=reference_model_1,markov_order=c(4,4,3,3),sample_fraction = 0.1)
result_mers_percent_10 <- classify_virus(data=test_dat_mers,reference_model=reference_model_1,markov_order=c(4,4,3,3),sample_fraction = 0.1)
result_dengue_percent_10 <- classify_virus(data=test_dat_dengue,reference_model=reference_model_1,markov_order=c(4,4,3,3),sample_fraction = 0.1)
result_hepatitis_percent_10 <- classify_virus(data=test_dat_hepatitis,reference_model=reference_model_1,markov_order=c(4,4,3,3),sample_fraction = 0.1)
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


result_sars_1_percent_05 <- classify_virus(data=test_dat_sars_set_1,reference_model=reference_model_1,markov_order=c(4,4,3,3),sample_fraction = 0.05)
result_sars_2_percent_05 <- classify_virus(data=test_dat_sars_set_2,reference_model=reference_model_1,markov_order=c(4,4,3,3),sample_fraction = 0.05)
result_sars_3_percent_05 <- classify_virus(data=test_dat_sars_set_3,reference_model=reference_model_1,markov_order=c(4,4,3,3),sample_fraction = 0.05)
result_sars_4_percent_05 <- classify_virus(data=test_dat_sars_set_4,reference_model=reference_model_1,markov_order=c(4,4,3,3),sample_fraction = 0.05)
result_mers_percent_05 <- classify_virus(data=test_dat_mers,reference_model=reference_model_1,markov_order=c(4,4,3,3),sample_fraction = 0.05)
result_dengue_percent_05 <- classify_virus(data=test_dat_dengue,reference_model=reference_model_1,markov_order=c(4,4,3,3),sample_fraction = 0.05)
result_hepatitis_percent_05 <- classify_virus(data=test_dat_hepatitis,reference_model=reference_model_1,markov_order=c(4,4,3,3),sample_fraction = 0.05)
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

saveRDS(Confusion_MX_05,"Confusion_Matrix_5_Percent_model_1.rds")
saveRDS(Confusion_MX_10,"Confusion_Matrix_10_Percent_model_1.rds")
saveRDS(Confusion_MX_25,"Confusion_Matrix_25_Percent_model_1.rds")




######### Result: Model 2 ######



virus <- c("SARS-Cov-2","MERS","Dengue","Hepatitis")
true_virus_label <- c(rep("SARS-Cov-2",200), rep("MERS",50),rep("Dengue",100),rep("Hepatitis",150))



result_sars_1_percent_25 <- classify_virus(data=test_dat_sars_set_1,reference_model=reference_model_2,markov_order=c(3,3,3,3),sample_fraction = 0.25)
result_sars_2_percent_25 <- classify_virus(data=test_dat_sars_set_2,reference_model=reference_model_2,markov_order=c(3,3,3,3),sample_fraction = 0.25)
result_sars_3_percent_25 <- classify_virus(data=test_dat_sars_set_3,reference_model=reference_model_2,markov_order=c(3,3,3,3),sample_fraction = 0.25)
result_sars_4_percent_25 <- classify_virus(data=test_dat_sars_set_4,reference_model=reference_model_2,markov_order=c(3,3,3,3),sample_fraction = 0.25)
result_mers_percent_25 <- classify_virus(data=test_dat_mers,reference_model=reference_model_2,markov_order=c(3,3,3,3),sample_fraction = 0.25)
result_dengue_percent_25 <- classify_virus(data=test_dat_dengue,reference_model=reference_model_2,markov_order=c(3,3,3,3),sample_fraction = 0.25)
result_hepatitis_percent_25 <- classify_virus(data=test_dat_hepatitis,reference_model=reference_model_2,markov_order=c(3,3,3,3),sample_fraction = 0.25)

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


result_sars_1_percent_10 <- classify_virus(data=test_dat_sars_set_1,reference_model=reference_model_2,markov_order=c(3,3,3,3),sample_fraction = 0.1)
result_sars_2_percent_10 <- classify_virus(data=test_dat_sars_set_2,reference_model=reference_model_2,markov_order=c(3,3,3,3),sample_fraction = 0.1)
result_sars_3_percent_10 <- classify_virus(data=test_dat_sars_set_3,reference_model=reference_model_2,markov_order=c(3,3,3,3),sample_fraction = 0.1)
result_sars_4_percent_10 <- classify_virus(data=test_dat_sars_set_4,reference_model=reference_model_2,markov_order=c(3,3,3,3),sample_fraction = 0.1)
result_mers_percent_10 <- classify_virus(data=test_dat_mers,reference_model=reference_model_2,markov_order=c(3,3,3,3),sample_fraction = 0.1)
result_dengue_percent_10 <- classify_virus(data=test_dat_dengue,reference_model=reference_model_2,markov_order=c(3,3,3,3),sample_fraction = 0.1)
result_hepatitis_percent_10 <- classify_virus(data=test_dat_hepatitis,reference_model=reference_model_2,markov_order=c(3,3,3,3),sample_fraction = 0.1)
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


result_sars_1_percent_05 <- classify_virus(data=test_dat_sars_set_1,reference_model=reference_model_2,markov_order=c(3,3,3,3),sample_fraction = 0.05)
result_sars_2_percent_05 <- classify_virus(data=test_dat_sars_set_2,reference_model=reference_model_2,markov_order=c(3,3,3,3),sample_fraction = 0.05)
result_sars_3_percent_05 <- classify_virus(data=test_dat_sars_set_3,reference_model=reference_model_2,markov_order=c(3,3,3,3),sample_fraction = 0.05)
result_sars_4_percent_05 <- classify_virus(data=test_dat_sars_set_4,reference_model=reference_model_2,markov_order=c(3,3,3,3),sample_fraction = 0.05)
result_mers_percent_05 <- classify_virus(data=test_dat_mers,reference_model=reference_model_2,markov_order=c(3,3,3,3),sample_fraction = 0.05)
result_dengue_percent_05 <- classify_virus(data=test_dat_dengue,reference_model=reference_model_2,markov_order=c(3,3,3,3),sample_fraction = 0.05)
result_hepatitis_percent_05 <- classify_virus(data=test_dat_hepatitis,reference_model=reference_model_2,markov_order=c(3,3,3,3),sample_fraction = 0.05)
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

saveRDS(Confusion_MX_05,"Confusion_Matrix_5_Percent_model_2.rds")
saveRDS(Confusion_MX_10,"Confusion_Matrix_10_Percent_model_2.rds")
saveRDS(Confusion_MX_25,"Confusion_Matrix_25_Percent_model_2.rds")



####### Result: Model 3 ########


virus <- c("SARS-Cov-2","MERS","Dengue","Hepatitis")
true_virus_label <- c(rep("SARS-Cov-2",200), rep("MERS",50),rep("Dengue",100),rep("Hepatitis",150))



result_sars_1_percent_25 <- classify_virus(data=test_dat_sars_set_1,reference_model=reference_model_3,markov_order=c(4,4,3,3),sample_fraction = 0.25)
result_sars_2_percent_25 <- classify_virus(data=test_dat_sars_set_2,reference_model=reference_model_3,markov_order=c(4,4,3,3),sample_fraction = 0.25)
result_sars_3_percent_25 <- classify_virus(data=test_dat_sars_set_3,reference_model=reference_model_3,markov_order=c(4,4,3,3),sample_fraction = 0.25)
result_sars_4_percent_25 <- classify_virus(data=test_dat_sars_set_4,reference_model=reference_model_3,markov_order=c(4,4,3,3),sample_fraction = 0.25)
result_mers_percent_25 <- classify_virus(data=test_dat_mers,reference_model=reference_model_3,markov_order=c(4,4,3,3),sample_fraction = 0.25)
result_dengue_percent_25 <- classify_virus(data=test_dat_dengue,reference_model=reference_model_3,markov_order=c(4,4,3,3),sample_fraction = 0.25)
result_hepatitis_percent_25 <- classify_virus(data=test_dat_hepatitis,reference_model=reference_model_3,markov_order=c(4,4,3,3),sample_fraction = 0.25)

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


result_sars_1_percent_10 <- classify_virus(data=test_dat_sars_set_1,reference_model=reference_model_3,markov_order=c(4,4,3,3),sample_fraction = 0.1)
result_sars_2_percent_10 <- classify_virus(data=test_dat_sars_set_2,reference_model=reference_model_3,markov_order=c(4,4,3,3),sample_fraction = 0.1)
result_sars_3_percent_10 <- classify_virus(data=test_dat_sars_set_3,reference_model=reference_model_3,markov_order=c(4,4,3,3),sample_fraction = 0.1)
result_sars_4_percent_10 <- classify_virus(data=test_dat_sars_set_4,reference_model=reference_model_3,markov_order=c(4,4,3,3),sample_fraction = 0.1)
result_mers_percent_10 <- classify_virus(data=test_dat_mers,reference_model=reference_model_3,markov_order=c(4,4,3,3),sample_fraction = 0.1)
result_dengue_percent_10 <- classify_virus(data=test_dat_dengue,reference_model=reference_model_3,markov_order=c(4,4,3,3),sample_fraction = 0.1)
result_hepatitis_percent_10 <- classify_virus(data=test_dat_hepatitis,reference_model=reference_model_3,markov_order=c(4,4,3,3),sample_fraction = 0.1)
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


result_sars_1_percent_05 <- classify_virus(data=test_dat_sars_set_1,reference_model=reference_model_3,markov_order=c(4,4,3,3),sample_fraction = 0.05)
result_sars_2_percent_05 <- classify_virus(data=test_dat_sars_set_2,reference_model=reference_model_3,markov_order=c(4,4,3,3),sample_fraction = 0.05)
result_sars_3_percent_05 <- classify_virus(data=test_dat_sars_set_3,reference_model=reference_model_3,markov_order=c(4,4,3,3),sample_fraction = 0.05)
result_sars_4_percent_05 <- classify_virus(data=test_dat_sars_set_4,reference_model=reference_model_3,markov_order=c(4,4,3,3),sample_fraction = 0.05)
result_mers_percent_05 <- classify_virus(data=test_dat_mers,reference_model=reference_model_3,markov_order=c(4,4,3,3),sample_fraction = 0.05)
result_dengue_percent_05 <- classify_virus(data=test_dat_dengue,reference_model=reference_model_3,markov_order=c(4,4,3,3),sample_fraction = 0.05)
result_hepatitis_percent_05 <- classify_virus(data=test_dat_hepatitis,reference_model=reference_model_3,markov_order=c(4,4,3,3),sample_fraction = 0.05)
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

saveRDS(Confusion_MX_05,"Confusion_Matrix_5_Percent_model_3.rds")
saveRDS(Confusion_MX_10,"Confusion_Matrix_10_Percent_model_3.rds")
saveRDS(Confusion_MX_25,"Confusion_Matrix_25_Percent_model_3.rds")




######## Result: Model4  (Essentially the Model 2 in our paper) ########

virus <- c("SARS-Cov-2","MERS","Dengue","Hepatitis")
true_virus_label <- c(rep("SARS-Cov-2",200), rep("MERS",50),rep("Dengue",100),rep("Hepatitis",150))



result_sars_1_percent_25 <- classify_virus(data=test_dat_sars_set_1,reference_model=reference_model_4,markov_order=c(3,3,3,3),sample_fraction = 0.25)
result_sars_2_percent_25 <- classify_virus(data=test_dat_sars_set_2,reference_model=reference_model_4,markov_order=c(3,3,3,3),sample_fraction = 0.25)
result_sars_3_percent_25 <- classify_virus(data=test_dat_sars_set_3,reference_model=reference_model_4,markov_order=c(3,3,3,3),sample_fraction = 0.25)
result_sars_4_percent_25 <- classify_virus(data=test_dat_sars_set_4,reference_model=reference_model_4,markov_order=c(3,3,3,3),sample_fraction = 0.25)
result_mers_percent_25 <- classify_virus(data=test_dat_mers,reference_model=reference_model_4,markov_order=c(3,3,3,3),sample_fraction = 0.25)
result_dengue_percent_25 <- classify_virus(data=test_dat_dengue,reference_model=reference_model_4,markov_order=c(3,3,3,3),sample_fraction = 0.25)
result_hepatitis_percent_25 <- classify_virus(data=test_dat_hepatitis,reference_model=reference_model_4,markov_order=c(3,3,3,3),sample_fraction = 0.25)

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


result_sars_1_percent_10 <- classify_virus(data=test_dat_sars_set_1,reference_model=reference_model_4,markov_order=c(3,3,3,3),sample_fraction = 0.1)
result_sars_2_percent_10 <- classify_virus(data=test_dat_sars_set_2,reference_model=reference_model_4,markov_order=c(3,3,3,3),sample_fraction = 0.1)
result_sars_3_percent_10 <- classify_virus(data=test_dat_sars_set_3,reference_model=reference_model_4,markov_order=c(3,3,3,3),sample_fraction = 0.1)
result_sars_4_percent_10 <- classify_virus(data=test_dat_sars_set_4,reference_model=reference_model_4,markov_order=c(3,3,3,3),sample_fraction = 0.1)
result_mers_percent_10 <- classify_virus(data=test_dat_mers,reference_model=reference_model_4,markov_order=c(3,3,3,3),sample_fraction = 0.1)
result_dengue_percent_10 <- classify_virus(data=test_dat_dengue,reference_model=reference_model_4,markov_order=c(3,3,3,3),sample_fraction = 0.1)
result_hepatitis_percent_10 <- classify_virus(data=test_dat_hepatitis,reference_model=reference_model_4,markov_order=c(3,3,3,3),sample_fraction = 0.1)
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


result_sars_1_percent_05 <- classify_virus(data=test_dat_sars_set_1,reference_model=reference_model_4,markov_order=c(3,3,3,3),sample_fraction = 0.05)
result_sars_2_percent_05 <- classify_virus(data=test_dat_sars_set_2,reference_model=reference_model_4,markov_order=c(3,3,3,3),sample_fraction = 0.05)
result_sars_3_percent_05 <- classify_virus(data=test_dat_sars_set_3,reference_model=reference_model_4,markov_order=c(3,3,3,3),sample_fraction = 0.05)
result_sars_4_percent_05 <- classify_virus(data=test_dat_sars_set_4,reference_model=reference_model_4,markov_order=c(3,3,3,3),sample_fraction = 0.05)
result_mers_percent_05 <- classify_virus(data=test_dat_mers,reference_model=reference_model_4,markov_order=c(3,3,3,3),sample_fraction = 0.05)
result_dengue_percent_05 <- classify_virus(data=test_dat_dengue,reference_model=reference_model_4,markov_order=c(3,3,3,3),sample_fraction = 0.05)
result_hepatitis_percent_05 <- classify_virus(data=test_dat_hepatitis,reference_model=reference_model_4,markov_order=c(3,3,3,3),sample_fraction = 0.05)
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

saveRDS(Confusion_MX_05,"Confusion_Matrix_5_Percent_model_4.rds")
saveRDS(Confusion_MX_10,"Confusion_Matrix_10_Percent_model_4.rds")
saveRDS(Confusion_MX_25,"Confusion_Matrix_25_Percent_model_4.rds")



####### Result: Model 5 ########


virus <- c("SARS-Cov-2","MERS","Dengue","Hepatitis")
true_virus_label <- c(rep("SARS-Cov-2",200), rep("MERS",50),rep("Dengue",100),rep("Hepatitis",150))



result_sars_1_percent_25 <- classify_virus(data=test_dat_sars_set_1,reference_model=reference_model_5,markov_order=c(4,4,4,3),sample_fraction = 0.25)
result_sars_2_percent_25 <- classify_virus(data=test_dat_sars_set_2,reference_model=reference_model_5,markov_order=c(4,4,4,3),sample_fraction = 0.25)
result_sars_3_percent_25 <- classify_virus(data=test_dat_sars_set_3,reference_model=reference_model_5,markov_order=c(4,4,4,3),sample_fraction = 0.25)
result_sars_4_percent_25 <- classify_virus(data=test_dat_sars_set_4,reference_model=reference_model_5,markov_order=c(4,4,4,3),sample_fraction = 0.25)
result_mers_percent_25 <- classify_virus(data=test_dat_mers,reference_model=reference_model_5,markov_order=c(4,4,4,3),sample_fraction = 0.25)
result_dengue_percent_25 <- classify_virus(data=test_dat_dengue,reference_model=reference_model_5,markov_order=c(4,4,4,3),sample_fraction = 0.25)
result_hepatitis_percent_25 <- classify_virus(data=test_dat_hepatitis,reference_model=reference_model_5,markov_order=c(4,4,4,3),sample_fraction = 0.25)

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


result_sars_1_percent_10 <- classify_virus(data=test_dat_sars_set_1,reference_model=reference_model_5,markov_order=c(4,4,4,3),sample_fraction = 0.1)
result_sars_2_percent_10 <- classify_virus(data=test_dat_sars_set_2,reference_model=reference_model_5,markov_order=c(4,4,4,3),sample_fraction = 0.1)
result_sars_3_percent_10 <- classify_virus(data=test_dat_sars_set_3,reference_model=reference_model_5,markov_order=c(4,4,4,3),sample_fraction = 0.1)
result_sars_4_percent_10 <- classify_virus(data=test_dat_sars_set_4,reference_model=reference_model_5,markov_order=c(4,4,4,3),sample_fraction = 0.1)
result_mers_percent_10 <- classify_virus(data=test_dat_mers,reference_model=reference_model_5,markov_order=c(4,4,4,3),sample_fraction = 0.1)
result_dengue_percent_10 <- classify_virus(data=test_dat_dengue,reference_model=reference_model_5,markov_order=c(4,4,4,3),sample_fraction = 0.1)
result_hepatitis_percent_10 <- classify_virus(data=test_dat_hepatitis,reference_model=reference_model_5,markov_order=c(4,4,4,3),sample_fraction = 0.1)
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


result_sars_1_percent_05 <- classify_virus(data=test_dat_sars_set_1,reference_model=reference_model_5,markov_order=c(4,4,4,3),sample_fraction = 0.05)
result_sars_2_percent_05 <- classify_virus(data=test_dat_sars_set_2,reference_model=reference_model_5,markov_order=c(4,4,4,3),sample_fraction = 0.05)
result_sars_3_percent_05 <- classify_virus(data=test_dat_sars_set_3,reference_model=reference_model_5,markov_order=c(4,4,4,3),sample_fraction = 0.05)
result_sars_4_percent_05 <- classify_virus(data=test_dat_sars_set_4,reference_model=reference_model_5,markov_order=c(4,4,4,3),sample_fraction = 0.05)
result_mers_percent_05 <- classify_virus(data=test_dat_mers,reference_model=reference_model_5,markov_order=c(4,4,4,3),sample_fraction = 0.05)
result_dengue_percent_05 <- classify_virus(data=test_dat_dengue,reference_model=reference_model_5,markov_order=c(4,4,4,3),sample_fraction = 0.05)
result_hepatitis_percent_05 <- classify_virus(data=test_dat_hepatitis,reference_model=reference_model_5,markov_order=c(4,4,4,3),sample_fraction = 0.05)
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

saveRDS(Confusion_MX_05,"Confusion_Matrix_5_Percent_model_5.rds")
saveRDS(Confusion_MX_10,"Confusion_Matrix_10_Percent_model_5.rds")
saveRDS(Confusion_MX_25,"Confusion_Matrix_25_Percent_model_5.rds")

