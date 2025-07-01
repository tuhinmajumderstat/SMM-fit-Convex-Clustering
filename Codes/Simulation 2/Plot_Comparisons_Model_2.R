library(Matrix)
library(clues)

#### Functions ####

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

compute_rand = function(assigned_matrix, correct_matrix){
  output = rep(0, nrow(assigned_matrix))
  for(i in 1:nrow(assigned_matrix)){
    #Remove histories that did not appear in simulation
    assigned = assigned_matrix[i,which(assigned_matrix[i,] > 0)]
    correct = correct_matrix[i,which(correct_matrix[i,] > 0)]
    if(length(correct) != length(assigned)){
      print(i)
    }
    output[i] = adjustedRand(assigned, correct, "HA")
  }
  return(output)
}


#### Data Generation: ####

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


setwd("C://Users/tm389/OneDrive - Duke University/Desktop/Submitted Papers/SMM Fit/Results/Simulation 2")



n_list <- c(1000, 2000)
Repl = 1000



list_method <- c("l_2","l_inf","l_1")
phi_list <- c(100,10,10)

#### Aggregrate: m=3 ####




aggregate_df <- NULL

for(T1 in 1:length(list_method))
{
  phi <- phi_list[T1]
  method <- list_method[T1]
  knn <- 15
  Summary_df  <- data.frame("n"=n_list)
  
  
  for(i in 1:length(n_list))
  {
    n <- n_list[i]
    
    d <- 4
    dat <- readRDS(paste0("Result_Weighted_phi_",phi,"_",method,"_knn_",knn,"_Model_2_n_",n,"_d_",d,"_m_",m,"_g_",g,".rds"))
    Rand_Index <- unlist(dat[1,])
    Adj_Rand_Index <- unlist(dat[2,])
    
    Summary_df$RI[i] <- mean(Rand_Index)
    Summary_df$ARI[i] <- mean(Adj_Rand_Index)
    Summary_df$ARI_sd[i] <- sd(Adj_Rand_Index)
    Summary_df$True_Recovery[i] <- length(which(Rand_Index==1))/Repl 
    Summary_df$method[i] <- paste0("CVX_knn_",knn,"_weight_",method)
    
  }
  aggregate_df <- rbind(aggregate_df,Summary_df)
}

    


Summary_df  <- data.frame("n"=n_list)
for(i in 1:length(n_list))
{
  n=n_list[i]
  
  dat <- readRDS(paste0("C://Users/tm389/OneDrive - Duke University/Desktop/Submitted Papers/SMM Fit/Results/Simulation 2/Compared Methods/JCGS_Model2_weylandt_n_",n,".rds"))
  Rand_Index <- unlist(lapply(dat,function(x)x$RI))
  Adj_Rand_Index <- unlist(lapply(dat,function(x)x$ARI))
  
  Summary_df$RI[i] <- mean(Rand_Index)
  Summary_df$ARI[i] <- mean(Adj_Rand_Index)
  Summary_df$ARI_sd[i] <- sd(Adj_Rand_Index)
  Summary_df$True_Recovery[i] <- length(which(Rand_Index==1))/Repl 
  Summary_df$method[i] <- paste0("Weylandt")
}

aggregate_df <- rbind(aggregate_df,Summary_df)



compared_methods <- c("dpmm","xiong","ggl","kmeans")
for(method in compared_methods)
{
  Summary_df  <- data.frame("n"=n_list)
  for(i in 1:length(n_list))
  {
    n=n_list[i]

    dat <- readRDS(paste0("C://Users/tm389/OneDrive - Duke University/Desktop/Submitted Papers/SMM Fit/Results/Simulation 2/Compared Methods/JCGS_Model_2_",method,"assigned_n_",n,".rds"))
    if(method=="dpmm")
    {
      dat <- dat$output_matrix
    }
    model2_correct <- readRDS(paste0("./Compared Methods/JCGS_Model_2_correct_model_n_",n,".rds"))
    Adj_Rand_Index <- compute_rand(dat, model2_correct) 
    Summary_df$RI[i] <- NA
    Summary_df$ARI[i] <- mean(Adj_Rand_Index)
    Summary_df$ARI_sd[i] <- sd(Adj_Rand_Index)
    Summary_df$True_Recovery[i] <- length(which(Adj_Rand_Index==1))/Repl 
    Summary_df$method[i] <- method
  }
  aggregate_df <- rbind(aggregate_df,Summary_df)
}

aggregate_df$ARI_sd <- round(aggregate_df$ARI_sd,3)


#### Plot: ARI ####


library(ggplot2)
library(grid)
library(dplyr)

# Color: distinct, colorblind-friendly palette
method_colors <- c(
  "CVX_knn_15_weight_l_2"   = "#1b9e77",
  "CVX_knn_15_weight_l_1"   = "#d95f02",
  "Weylandt" = "#e7298a",
  "CVX_knn_15_weight_l_inf" = "#7570b3",
  "dpmm"               = "#66a61e",
  "xiong"              = "#e6ab02",
  "ggl"                = "#a6761d",
  "kmeans"             = "#666666"
)

# Linetypes: consistent, 8-style
method_linetypes <- c(
  "CVX_knn_15_weight_l_2" ="dotdash",
  "CVX_knn_15_weight_l_1"  = "dotted",
  "Weylandt" = "dashed",
  "CVX_knn_15_weight_l_inf"   = "solid",
  "dpmm"               = "longdash",
  "xiong"              = "twodash",
  "ggl"                = "dotted",
  "kmeans"             = "dotdash"
)

plot1_df <- aggregate_df %>% filter(method %in% names(method_colors)[1:4])

ggplot(plot1_df, aes(x = n, y = ARI, color = method, linetype = method)) +
  geom_line(size = 1.5) +
  labs(x = "Sample Size (n)", y = "Mean Adjusted Rand Index (ARI)", title = "Mean ARI – knn methods") +
  scale_color_manual(values = method_colors) +
  scale_linetype_manual(values = method_linetypes) +
  guides(color = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2)) +
  theme_bw() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    plot.title = element_text(size = 35, face = 'bold'),
    legend.title = element_text(size = 20, face = 'bold'),
    legend.text = element_text(size = 20),
    legend.position = "bottom",
    legend.key.width = unit(2, "cm"),
    legend.spacing.x = unit(1.0, 'cm'),
    axis.title = element_text(face = "bold", size = 30),
    axis.text = element_text(size = 25)
  )


plot2_df <- aggregate_df %>% filter(method %in% names(method_colors)[3:8])

ggplot(plot2_df, aes(x = n, y = ARI, color = method, linetype = method)) +
  geom_line(size = 1.5) +
  labs(x = "Sample Size (n)", y = "Mean Adjusted Rand Index (ARI)", title = "Mean ARI – other methods") +
  scale_color_manual(values = method_colors) +
  scale_linetype_manual(values = method_linetypes) +
  guides(color = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2)) +
  theme_bw() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    plot.title = element_text(size = 35, face = 'bold'),
    legend.title = element_text(size = 20, face = 'bold'),
    legend.text = element_text(size = 20),
    legend.position = "bottom",
    legend.key.width = unit(2, "cm"),
    legend.spacing.x = unit(1.0, 'cm'),
    axis.title = element_text(face = "bold", size = 30),
    axis.text = element_text(size = 25)
  )


#### Plot: True Recovery ####


plot1_df <- aggregate_df %>% filter(method %in% names(method_colors)[1:4])

ggplot(plot1_df, aes(x = n, y = True_Recovery, color = method, linetype = method)) +
  geom_line(size = 1.5) +
  labs(x = "Sample Size (n)", y = "True Recovery Proportion", title = "True Recovery – knn methods") +
  scale_color_manual(values = method_colors) +
  scale_linetype_manual(values = method_linetypes) +
  guides(color = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2)) +
  theme_bw() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    plot.title = element_text(size = 35, face = 'bold'),
    legend.title = element_text(size = 20, face = 'bold'),
    legend.text = element_text(size = 20),
    legend.position = "bottom",
    legend.key.width = unit(2, "cm"),
    legend.spacing.x = unit(1.0, 'cm'),
    axis.title = element_text(face = "bold", size = 30),
    axis.text = element_text(size = 25)
  )


plot2_df <- aggregate_df %>% filter(method %in% names(method_colors)[3:8])

ggplot(plot2_df, aes(x = n, y = True_Recovery, color = method, linetype = method)) +
  geom_line(size = 1.5) +
  labs(x = "Sample Size (n)", y = "True Recovery Proportion", title = "True Recovery – other methods") +
  scale_color_manual(values = method_colors) +
  scale_linetype_manual(values = method_linetypes) +
  guides(color = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2)) +
  theme_bw() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    plot.title = element_text(size = 35, face = 'bold'),
    legend.title = element_text(size = 20, face = 'bold'),
    legend.text = element_text(size = 20),
    legend.position = "bottom",
    legend.key.width = unit(2, "cm"),
    legend.spacing.x = unit(1.0, 'cm'),
    axis.title = element_text(face = "bold", size = 30),
    axis.text = element_text(size = 25)
  )


