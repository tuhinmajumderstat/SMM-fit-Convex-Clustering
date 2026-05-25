library(seqinr)
library(cvxclustr)
library(gtools)
library(dplyr)
library(Matrix)
library(parallel)
library(purrr)

## Change the location to the locations of the full data sets

#### SARS ####

full_dat_sars_set_1 <- read.fasta(file="C:/Users/tmajumd/Desktop/Graduate Studies/PhD/Data Analysis/Gene Classification SMM/Complete Test Genes/SARS_Cov_2_Sample_Set_1.fasta")
set.seed(12345)
samp_name_sars_1 <- sample(names(full_dat_sars_set_1),50)
test_dat_sars_set_1 <- full_dat_sars_set_1 %>% keep(names(.) %in% samp_name_sars_1)
rm(full_dat_sars_set_1)


full_dat_sars_set_2 <- read.fasta(file="C:/Users/tmajumd/Desktop/Graduate Studies/PhD/Data Analysis/Gene Classification SMM/Complete Test Genes/SARS_Cov_2_Sample_Set_2.fasta")
set.seed(12345)
samp_name_sars_2 <- sample(names(full_dat_sars_set_2),50)
test_dat_sars_set_2 <- full_dat_sars_set_2 %>% keep(names(.) %in% samp_name_sars_2)
rm(full_dat_sars_set_2)


full_dat_sars_set_3 <- read.fasta(file="C:/Users/tmajumd/Desktop/Graduate Studies/PhD/Data Analysis/Gene Classification SMM/Complete Test Genes/SARS_Cov_2_Sample_Set_3.fasta")
set.seed(12345)
samp_name_sars_3 <- sample(names(full_dat_sars_set_3),50)
test_dat_sars_set_3 <- full_dat_sars_set_3 %>% keep(names(.) %in% samp_name_sars_3)
rm(full_dat_sars_set_3)

full_dat_sars_set_4 <- read.fasta(file="C:/Users/tmajumd/Desktop/Graduate Studies/PhD/Data Analysis/Gene Classification SMM/Complete Test Genes/SARS_Cov_2_Sample_Set_4.fasta")
set.seed(12345)
samp_name_sars_4 <- sample(names(full_dat_sars_set_4),50)
test_dat_sars_set_4 <- full_dat_sars_set_4 %>% keep(names(.) %in% samp_name_sars_4)
rm(full_dat_sars_set_4)


setwd("C:/Users/tmajumd/Desktop/Graduate Studies/PhD/Data Analysis/Gene Classification SMM/Complete Test Genes")
saveRDS(test_dat_sars_set_1,"SARS-Cov-2_Test_Sample_Set_1.rds")
saveRDS(test_dat_sars_set_2,"SARS-Cov-2_Test_Sample_Set_2.rds")
saveRDS(test_dat_sars_set_3,"SARS-Cov-2_Test_Sample_Set_3.rds")
saveRDS(test_dat_sars_set_4,"SARS-Cov-2_Test_Sample_Set_4.rds")




#### MERS ####

full_dat_mers_set_1 <- read.fasta(file="C:/Users/tmajumd/Desktop/Graduate Studies/PhD/Data Analysis/Gene Classification SMM/Complete Test Genes/MERS_Sample_Set_1.fasta")
set.seed(12345)
samp_name_mers <- sample(names(full_dat_mers_set_1),50)
test_dat_mers <- full_dat_mers_set_1 %>% keep(names(.) %in% samp_name_mers)
rm(full_dat_mers_set_1)

saveRDS(test_dat_mers,"MERS_Test_Sample_Set.rds")



#### Dengue ####

full_dat_dengue_set_1 <- read.fasta(file="C:/Users/tmajumd/Desktop/Graduate Studies/PhD/Data Analysis/Gene Classification SMM/Complete Test Genes/Dengue_Sample_Set_1.fasta")
set.seed(12345)
samp_name_dengue <- sample(names(full_dat_dengue_set_1),100)
test_dat_dengue <- full_dat_dengue_set_1 %>% keep(names(.) %in% samp_name_dengue)
rm(full_dat_dengue_set_1)

saveRDS(test_dat_dengue,"Dengue_Test_Sample_Set.rds")



#### Hepatits ####

full_dat_hepatitis_set_1 <- read.fasta(file="C:/Users/tmajumd/Desktop/Graduate Studies/PhD/Data Analysis/Gene Classification SMM/Complete Test Genes/Hepatitis_Sample_Set_1.fasta")
set.seed(12345)
samp_name_hepatitis <- sample(names(full_dat_hepatitis_set_1),150)
test_dat_hepatitis <- full_dat_hepatitis_set_1 %>% keep(names(.) %in% samp_name_hepatitis)
rm(full_dat_hepatitis_set_1)

saveRDS(test_dat_hepatitis,"Hepatitis_Test_Sample_Set.rds")
