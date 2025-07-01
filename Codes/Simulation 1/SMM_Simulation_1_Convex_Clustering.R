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


#### Modelfit Function ####

modelfit <- function(iter,R,group_id,m=2,n=10000,lambda_seq = seq(0,0.05,by=0.0001),phi=100,knn=5,type="kernel")
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
  
  if(type=="kernel")
  {
    w <- kernel_weights(X,phi) 
    w <- knn_weights(w,k=knn,p)
  }else{
    if(type=="l_inf")
    {
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
          w[start] <- exp(-phi*wt[i,j]^2)
          start <- start+1
        }
      }
      w <- knn_weights(w,k=knn,p)
    }else{
      w <- rep(1,choose(ncol(X),2))
    }
    
    
  }
  
  
  
  ## BIC tuning param select ##
  
  
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
  }
  
  plot(lambda_seq,BIC,type="l")
  lambda_min <- lambda_seq[which.min(BIC)]
  A <- create_adjacency(model_all$V[[which.min(BIC)]],w,ncol(X))
  opt_assignment <- find_clusters(A)
  opt_assignment
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
  
  return(list("RI"=RI, "ARI"=ARI, "Lambda"=lambda_min,"non_zero_row_id"=keep_rows,"cluster"=opt_assignment$cluster,
              "size"=opt_assignment$size,"estimated_adjacency"=estimated_adjacency))
}



#### Fit SMM for each Replication ####

## Sample size, varies from 5000 to 25000
n <- 10000

## Number of replications
Repl <- 1000



lambda_seq <- c(seq(0,0.4,by=0.0002),seq(0.5,5,by=0.1))
modelfit1 <- function(iter)modelfit(iter,R=R,n=n,group_id=group_id,lambda_seq =lambda_seq,phi=100,knn=5,type="l_inf")
modelfit2 <- function(iter)modelfit(iter,R=R,n=n,group_id=group_id,lambda_seq =lambda_seq,phi=100,knn=5,type="kernel")
modelfit3 <- function(iter)modelfit(iter,R=R,n=n,group_id=group_id,lambda_seq =lambda_seq,phi=100,knn=3,type="l_inf")
modelfit4 <- function(iter)modelfit(iter,R=R,n=n,group_id=group_id,lambda_seq =lambda_seq,phi=100,knn=3,type="kernel")
modelfit5 <- function(iter)modelfit(iter,R=R,n=n,m=m,group_id=group_id,lambda_seq =lambda_seq,type="uniform")


## Fix number of cores depending on your machine
ncores <- 16
ncores <- min(ncores, detectCores()-2)

cl=makeCluster(ncores)
clusterExport(cl,c("modelfit","p","n","all_comb","m","find_index","d","groups",
                   "group_id","R","cvxclust","create_adjacency","find_clusters","%>%","mutate",
                   "group_by","summarise_all","select","Matrix","true_adjacency","lambda_seq",
                   "kernel_weights","knn_weights"))
Output_wt_1 <- parSapply(cl=cl,as.matrix(c(1:Repl)),modelfit1)
stopCluster(cl)


cl=makeCluster(ncores)
clusterExport(cl,c("modelfit","p","n","all_comb","m","find_index","d","groups",
                   "group_id","R","cvxclust","create_adjacency","find_clusters","%>%","mutate",
                   "group_by","summarise_all","select","Matrix","true_adjacency","lambda_seq",
                   "kernel_weights","knn_weights"))
Output_wt_2 <- parSapply(cl=cl,as.matrix(c(1:Repl)),modelfit2)
stopCluster(cl)


cl=makeCluster(ncores)
clusterExport(cl,c("modelfit","p","n","all_comb","m","find_index","d","groups",
                   "group_id","R","cvxclust","create_adjacency","find_clusters","%>%","mutate",
                   "group_by","summarise_all","select","Matrix","true_adjacency","lambda_seq",
                   "kernel_weights","knn_weights"))
Output_wt_3 <- parSapply(cl=cl,as.matrix(c(1:Repl)),modelfit3)
stopCluster(cl)


cl=makeCluster(ncores)
clusterExport(cl,c("modelfit","p","n","all_comb","m","find_index","d","groups",
                   "group_id","R","cvxclust","create_adjacency","find_clusters","%>%","mutate",
                   "group_by","summarise_all","select","Matrix","true_adjacency","lambda_seq",
                   "kernel_weights","knn_weights"))
Output_wt_4 <- parSapply(cl=cl,as.matrix(c(1:Repl)),modelfit4)
stopCluster(cl)



cl=makeCluster(ncores)
clusterExport(cl,c("modelfit","p","n","all_comb","m","find_index","d","groups",
                   "group_id","R","cvxclust","create_adjacency","find_clusters","%>%","mutate",
                   "group_by","summarise_all","select","Matrix","true_adjacency","lambda_seq",
                   "kernel_weights","knn_weights"))
Output_wt_5 <- parSapply(cl=cl,as.matrix(c(1:Repl)),modelfit5)
stopCluster(cl)

saveRDS(Output_wt_1,file=paste0("Result_Weighted_phi_100_l_inf_knn_5_Model_1_n_",n,"_d_",d,"_m_",m,"_g_",g,".rds"))
saveRDS(Output_wt_2,file=paste0("Result_Weighted_phi_100_l_2_knn_5_Model_1_n_",n,"_d_",d,"_m_",m,"_g_",g,".rds"))
saveRDS(Output_wt_3,file=paste0("Result_Weighted_phi_100_l_inf_knn_3_Model_1_n_",n,"_d_",d,"_m_",m,"_g_",g,".rds"))
saveRDS(Output_wt_4,file=paste0("Result_Weighted_phi_100_l_2_knn_3_Model_1_n_",n,"_d_",d,"_m_",m,"_g_",g,".rds"))
saveRDS(Output_wt_5,file=paste0("Result_Weighted_Uniform_Model_1_n_",n,"_d_",d,"_m_",m,"_g_",g,".rds"))
