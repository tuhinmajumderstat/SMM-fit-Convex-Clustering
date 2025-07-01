library(clustRviz)
library(cvxclustr)
library(gtools)
library(dplyr)
library(Matrix)
library(parallel)

#### Some initial Functions ####

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


vector_to_symmetric_zero_diag <- function(vec, n) {
  mat <- matrix(0, n, n)
  mat[upper.tri(mat)] <- vec
  mat <- mat + t(mat)  # Reflect to lower triangle
  return(mat)
}


#### Data Generation ####


##   Number of states
d <- 4   
## Order of SMM
m <- 2   
## Number of m-tuples
p <- d^m 

all_comb <- expand.grid(rep(list(1:d), m))

## Number of partitions
g <- sqrt(p)
groups <- c(1:g)
set.seed(12345)
nums <- create_partition(p=p,g=g,min_size = g)
nums

## True Cluster assignment
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


## Generate True Transition Probabilities
set.seed(23456)
Z <- matrix(runif(g*d),nrow=g,ncol=d)
weights <- exp(Z)

R <- matrix(0,g,d,byrow=TRUE)
for(i in 1:g)
{
  R[i,] <- rdirichlet (1,weights[i,])
}
R

#### Modelfit by Weylandt ####

modelfit_weylandt <- function(iter,R,group_id,m,n )
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
  
  N <- matrix(0,p,d,byrow=TRUE)
  for(i in (m+1):n)
  {
    prev <- chain[(i-m):(i-1)]
    j <- chain[i]
    tuple_ind <- find_index(prev,d)
    N[tuple_ind,j] <- N[tuple_ind,j]+1
  }
  
  
  P <- N/rowSums(N)
  X <- t(P)
  keep_rows <- which(rowSums(N)>0)
  X <- X[,keep_rows]
  nu <- 1.8*(1/p)
  
  
  
  
  ## Pathwise Method: Weylandt ##
  
  X_weylandt <- t(X)
  
  model_all <- CARP(X_weylandt)
  all_result <- model_all$cluster_membership
  
  N_data <- as.data.frame(N[keep_rows,])
  names(N_data) <- paste0("state",c(1:d))
  
  max_iter <- max(model_all$cluster_membership$Iter)
  BIC <- matrix(0,max_iter,1,byrow=TRUE)
  
  
  
  for(j in 1:max_iter)
  {
    assignment <- all_result %>% filter(Iter==j)
    cluster_id <- unique(assignment$Cluster)
    N_clustered <- N_data %>% mutate(cluster=assignment$Cluster)
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
  }
  
  plot(1:max_iter,BIC,type="l")
  iter_min <- which.min(BIC)
  
  opt_assignment <- all_result %>% filter(Iter==iter_min)
  opt_assignment$cluster <- opt_assignment$Cluster
  unique_cluster_id <- unique(opt_assignment$cluster)
  
  estimated_adjacency <- Matrix(0, nrow = ncol(X), ncol = ncol(X), sparse = TRUE)
  for(i in 1:length(unique_cluster_id))
  {
    cluster_indices <- which(opt_assignment$cluster==unique_cluster_id[i])
    n_g <- length(cluster_indices )
    if(n_g>1)
    {
      for(j in 1:(n_g-1))
      {
        for(l in (j+1):n_g)
        {
          estimated_adjacency[cluster_indices[j],cluster_indices[l]]=1
        }
        
      }
    }
    
  }
  
  true_adjacency_new <- true_adjacency[keep_rows,keep_rows]
  
  mismatch <- sum(abs(estimated_adjacency-true_adjacency_new))
  RI <- 1-mismatch/(choose(ncol(X),2))
  
  a=length(which(which(as.matrix(true_adjacency_new==1))%in% which(as.matrix(estimated_adjacency==1))))
  d1=length(which(which(as.matrix(true_adjacency_new==0))%in% which(as.matrix(estimated_adjacency==0))))-(p*(p+1))/2
  b=length(which(which(as.matrix(true_adjacency_new==1))%in% which(as.matrix(estimated_adjacency==0))))
  c=length(which(which(as.matrix(true_adjacency_new==0))%in% which(as.matrix(estimated_adjacency==1))))
  
  ARI= (choose(ncol(X),2)*(a+d1)- ((a+b)*(a+c)+(c+d1)*(b+d1)))/((choose(ncol(X),2))^2-((a+b)*(a+c)+(c+d1)*(b+d1)))
  
  return(list("RI"=RI, "ARI"=ARI, "non_zero_row_id"=keep_rows,"cluster"=opt_assignment$cluster,
              "estimated_adjacency"=estimated_adjacency))
}


#### Fit SMM for each Replication and for each n ####

library(foreach)
library(doParallel)

# Number of replications
Repl <- 1000

## Fix number of cores depending on your machine
num_cores <- 20
num_cores <- min(num_cores, detectCores()-2)

for(n in 5000*c(1:5))
{
  modelfit1 <- function(iter)modelfit_weylandt(iter,R=R,n=n,m=2,group_id=group_id) 
  
  time=proc.time()
  # Set number of cores (you can adjust this based on your machine or HPC cluster)
  
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  # Run in parallel using foreach
  results_list <- foreach(iter = 1:Repl, .packages = c("Matrix", "dplyr", "cvxclustr", "clustRviz")) %dopar% {
    tryCatch({
      modelfit1(iter)
    }, error = function(e) {
      message(paste("Error in iteration", iter, ":", e$message))
      return(NULL)
    })
  }
  
  # Stop the cluster after use
  stopCluster(cl)
  
  time=proc.time()-time
  
  saveRDS(results_list,paste0("JCGS_Model1_weylandt_n_",n,".rds"))
  
  ari=sapply(results_list, function(x) x$ARI)
  print(paste0("Sample: ",n))
  print(mean(ari))
  print(length(which(ari==1))/1000)
  print(time)
  
}









