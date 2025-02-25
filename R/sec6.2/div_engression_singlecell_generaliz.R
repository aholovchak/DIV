library(ggplot2)
library(MASS)
library(splines)

library(engression)
library(tidyr)
library(dplyr)

setwd("~/Documents/ethz/DIV/r_package/DistributionIV")
library(devtools)
load_all()

dataset_rpe1 <- read.csv("data/sec6.2/dataset_rpe1.csv", stringsAsFactors=TRUE)
test_envs <- read.csv("data/sec6.2/single_cell_test_envs.csv", stringsAsFactors=TRUE)


#########################
# identified case
#########################
interv_genes <- colnames(dataset_rpe1)[1:9]

train_data <- dataset_rpe1 %>% filter(interventions %in% c(interv_genes, "non-targeting"))
#train_data$interventions <- droplevels(train_data$interventions)

test_data <- dataset_rpe1 %>% filter(interventions %in% test_envs$genes)
#test_data$interventions <- droplevels(test_data$interventions)

Xtr <- train_data[,1:9]
Ytr <- train_data[,10]
Ztr <- droplevels(train_data[,11])


set.seed(2901)

num_repeats <- 10
DIV_singlecell_10runs <- data.frame(matrix(nrow=length(test_data$interventions), ncol=num_repeats))
colnames(DIV_singlecell_10runs) <- paste0("Run", 1:num_repeats)

###################
# DIV
###################
for (rep in 1:num_repeats) {
  cat("DIV Running iteration", rep, "out of", num_repeats, "\n")
  div_mod <- div(Z = Ztr, X = Xtr, Y = Ytr, num_layer = 4, num_epochs = 20000, lr = 1e-3, joint = TRUE,
                 epsx_dim = 50, epsh_dim = 50, epsy_dim = 50)
  mean_div <- predict(div_mod, Xtest = test_data[,1:9], type = "mean", nsample = 1000)
  DIV_singlecell_10runs[,rep] <- mean_div
}

################################################################################

set.seed(2901)

engr_singlecell_10runs <- data.frame(matrix(nrow=length(test_data$interventions), ncol=num_repeats))
colnames(engr_singlecell_10runs) <- paste0("Run", 1:num_repeats)

###################
# Engression
###################
for (rep in 1:num_repeats) {
  cat("Engression Running iteration", rep, "out of", num_repeats, "\n")
  engr_mod <- engression(X = Xtr, Y = Ytr, num_layer = 4, num_epochs = 1, 
                        lr = 1e-3, hidden_dim = 100, noise_dim = 100)
  mean_engr <- predict(engr_mod, Xtest = test_data[,1:9], type = "mean", nsample = 1)
  engr_singlecell_10runs[,rep] <- mean_engr
}

write.csv(DIV_singlecell_10runs, "results/sec6.2/div_singlecell_10runs.csv", row.names = FALSE)
write.csv(engr_singlecell_10runs, "results/sec6.2/engr_singlecell_10runs.csv", row.names = FALSE)

