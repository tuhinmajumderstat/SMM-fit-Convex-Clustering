
library(cvxclustr)
library(gtools)
library(dplyr)
library(Matrix)
library(parallel)

#### Some Functions ####

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

create_partition <- function(p,g,min_size=2)
{
  min_allot <- min(floor(p/g),min_size)
  final_allot <- rep(min_allot,g)
  excess_n <- p-g*min_allot
  excess_assn <- sample(c(1:g),excess_n,replace=TRUE)
  final_allot <- final_allot+ sapply(c(1:g),function(x)length(which(excess_assn==x)))
  final_allot
}





#### Model 2: n=1000, 2000, m=3 ####

setwd("C://Users/tm389/OneDrive - Duke University/Desktop/Submitted Papers/SMM Fit/SMM-fit-Convex-Clustering")

source("SMM_Fit_Other_Methods_Functions.R")

d <- 4
m <- 3
p <- d^m

all_comb <- expand.grid(rep(list(1:d), m))

g <- 4
groups <- c(1:g)
set.seed(23456)
nums <- create_partition(p=p,g=g,min_size = 10)
nums
assignments <- rep(NA, p)
assignments[sample(p)] <- rep(c(1:g), nums)
group_id <- assignments

true_adjacency <- Matrix(0, nrow = p, ncol = p, sparse = TRUE)
for(i in 1:g)
{
  cluster_indices <- which(group_id==groups[i])
  n_g <- length(cluster_indices )
  if(n_g>1)
  {
    for(j in 1:(n_g-1))
    {
      for(l in (j+1):n_g)
      {
        true_adjacency[cluster_indices[j],cluster_indices[l]]=1
      }
      
    }
  }
  
}


R <- matrix(0.1,ncol=d,nrow=g)
diag(R) <- 0.7



Repl <- 1000


for(n in 1000*c(1:2))
{
  model2_dat_1K <- matrix(0,nrow=Repl,ncol=n)
  for(iter in 1:Repl)
  {
    d <- ncol(R)
    p <- d^m
    set.seed(12345+iter)
    init <- sample(c(1:p),1)
    chain <- array(dim=n)
    chain[1:m] <- as.numeric(all_comb[init,])
    for(i in (m+1):n)
    {
      prev <- chain[(i-m):(i-1)]
      prev_ind <- find_index(prev,d)
      grp_ind <- which(groups==group_id[prev_ind])
      chain[i] <- sample(c(1:d),1,prob=R[grp_ind,])
    }
    model2_dat_1K[iter,]=chain
  }
  
  model2_dat_1K <- apply(model2_dat_1K, 1, function(x) paste(x, collapse = ""))
  
  
  # Step 1: Convert triplets in all_comb to strings
  triplet_strings <- apply(all_comb, 1, paste0, collapse = "")
  
  # Step 2: Assign corresponding group_id values as a named list
  subset_belong_new <- as.list(group_id)
  names(subset_belong_new) <- triplet_strings
  
  
  model2_JCGS <- SMM$new(subset_belong= subset_belong_new,
                         subset_probs = R,
                         order = m, sigma = c("1","2","3","4"))
  
  print(paste0("Start DPMM:",n))
  model2_dpmmassigned_1K = run_dpmm(model2_dat_1K, 1000, 3, c("1","2","3","4"), 500, .1)
  saveRDS(model2_dpmmassigned_1K,paste0("JCGS_Model_2_dpmmassigned_n_",n,".rds"))
  
  print(paste0("Start Xiong:",n))
  model2_xiongassigned_1K = run_xiong(model2_dat_1K, 1000, 3, c("1","2","3","4"))
  saveRDS(model2_xiongassigned_1K,paste0("JCGS_Model_2_xiongassigned_n_",n,".rds"))
  
  print(paste0("Start GGL:",n))
  model2_gglassigned_1K = run_ggl(model2_dat_1K, 1000, 3, c("1","2","3","4"))
  saveRDS(model2_gglassigned_1K ,paste0("JCGS_Model_2_gglassigned_n_",n,".rds"))
  
  print(paste0("Start Kmeans:",n))
  model2_kmeansassigned_1K = run_kmeans(model2_dat_1K, 1000, 3, c("1","2","3","4"))
  saveRDS(model2_kmeansassigned_1K,paste0("JCGS_Model_2_kmeansassigned_n_",n,".rds"))
  
  model2_correct_1K = correct_matrix(model2_dat_1K, model2_JCGS)
  saveRDS(model2_correct_1K,paste0("JCGS_Model_2_correct_model_n_",n,".rds"))
  
  print(paste0("Completed All:",n))
}
