setwd("C://Users/tm389/OneDrive - Duke University/Desktop/Submitted Papers/SMM Fit/Results/Data Analysis")


#### Model 1 ####

## CVX Confusion Matrix ##

percents <- c(5,10,25)
aggregate_df <- data.frame("percent"=percents)
aggregate_df$missclassification <- 0 
aggregate_df$method <- "CVX"
for(i in 1:3)
{
  err <- 1-sum(diag(readRDS(paste0("./Confusion Matrices/Convex Clustering Approach/Confusion_Matrix_",percents[i],"_Percent_model_1.rds"))))/500
  aggregate_df$missclassification[i] = err
}

## Other methods ##

methods <- c("dpmm","xiong","ggl","kmeans")
methods_renamed <- setNames(c("GSDPMM","Xiong","GGL","K-means"),compared_methods )
for(method in methods)
{
  summary_df <- data.frame("percent"=percents)
  summary_df$missclassification <- 0 
  summary_df$method <- compared_methods_renamed[[method]]
  
  for(i in 1:3)
  {
    err <- 1-sum(diag(readRDS(paste0("./Confusion Matrices/Compared Methods/Confusion_MX_",percents[i],"_model_1_method_",method,".rds"))))/500
    summary_df$missclassification[i]=err
  }
  aggregate_df <- rbind(aggregate_df,summary_df)
}


library(ggplot2)
library(grid)
library(dplyr)


# Color: distinct, colorblind-friendly palette
method_colors <- c(
  "CVX" = "#e7298a",
  "GSDPMM"                 = "#66a61e",
  "Xiong"                  = "#e6ab02",
  "GGL"                    = "#a6761d",
  "K-means"                = "#666666"
)

# Linetypes: consistent, 8-style
method_linetypes <- c(
  "CVX"="twodash",
  "GSDPMM" = "solid",
  "Xiong"  = "dashed",
  "GGL"    = "dotted",
  "K-means"= "longdash"
)

plot1_df <- aggregate_df %>% filter(method %in% names(method_colors))

p_model1 <- ggplot(plot1_df, aes(x = percent, y = missclassification, color = method, linetype = method)) +
  geom_line(size = 1.5) +
  labs(x = "Percentage of Available DNA String", y = "Misclassification Rate",
       title = "Missclassification in DNA Data: Model 1",color="Method",linetype="Method") +
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
    legend.title = element_text(size = 30, face = 'bold'),
    legend.text = element_text(size = 25),
    legend.position = "bottom",
    legend.key.width = unit(2, "cm"),
    legend.spacing.x = unit(1.0, 'cm'),
    axis.title = element_text(face = "bold", size = 30),
    axis.text = element_text(size = 25)
  )

p_model1
ggsave("plot_missclassification_model_1.jpeg",p_model1,height=8,width=12,dpi=600)
#### Model 2 ####

## CVX Confusion Matrix ##

percents <- c(5,10,25)
aggregate_df <- data.frame("percent"=percents)
aggregate_df$missclassification <- 0 
aggregate_df$method <- "CVX"
for(i in 1:3)
{
  err <- 1-sum(diag(readRDS(paste0("./Confusion Matrices/Convex Clustering Approach/Confusion_Matrix_",percents[i],"_Percent_model_4.rds"))))/500
  aggregate_df$missclassification[i] = err
}

## Other methods ##

methods <- c("dpmm","xiong","ggl","kmeans")
methods_renamed <- setNames(c("GSDPMM","Xiong","GGL","K-means"),compared_methods )
for(method in methods)
{
  summary_df <- data.frame("percent"=percents)
  summary_df$missclassification <- 0 
  summary_df$method <- compared_methods_renamed[[method]]
  
  for(i in 1:3)
  {
    err <- 1-sum(diag(readRDS(paste0("./Confusion Matrices/Compared Methods/Confusion_MX_",percents[i],"_model_2_method_",method,".rds"))))/500
    summary_df$missclassification[i]=err
  }
  aggregate_df <- rbind(aggregate_df,summary_df)
}


library(ggplot2)
library(grid)
library(dplyr)
# Color: distinct, colorblind-friendly palette
method_colors <- c(
  "CVX" = "#e7298a",
  "GSDPMM"                 = "#66a61e",
  "Xiong"                  = "#e6ab02",
  "GGL"                    = "#a6761d",
  "K-means"                = "#666666"
)

# Linetypes: consistent, 8-style
method_linetypes <- c(
  "CVX"="twodash",
  "GSDPMM" = "solid",
  "Xiong"  = "dashed",
  "GGL"    = "dotted",
  "K-means"= "longdash"
)

plot1_df <- aggregate_df %>% filter(method %in% names(method_colors))

p_model2 <- ggplot(plot1_df, aes(x = percent, y = missclassification, color = method, linetype = method)) +
  geom_line(size = 1.5) +
  labs(x = "Percentage of Available DNA String", y = "Misclassification Rate",
       title = "Missclassification in DNA Data: Model 2",color="Method",linetype="Method") +
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
    legend.title = element_text(size = 30, face = 'bold'),
    legend.text = element_text(size = 25),
    legend.position = "bottom",
    legend.key.width = unit(2, "cm"),
    legend.spacing.x = unit(1.0, 'cm'),
    axis.title = element_text(face = "bold", size = 30),
    axis.text = element_text(size = 25)
  )
p_model2
ggsave("plot_missclassification_model_2.jpeg",p_model2,height=8,width=12,dpi=600)
