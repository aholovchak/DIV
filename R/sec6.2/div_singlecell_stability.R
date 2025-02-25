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
interv_genes <- colnames(dataset_rpe1)[1:9]  # Select first 9 columns as intervention genes

# Filter training data based on interventions
train_data <- dataset_rpe1 %>% filter(interventions %in% c(interv_genes, "non-targeting"))
train_data$interventions <- droplevels(train_data$interventions)

# No need to manually encode categorical variables since the model does this internally

Xtr <- train_data[, 1:9]   # First 9 columns
Ytr <- train_data[, 10]    # 10th column (target)
Ztr <- train_data[, 11] # Columns 12-20 (environmental variables)
summary(Ztr)

test_data <- dataset_rpe1 %>% filter(interventions %in% test_envs$genes)
Xtest <- test_data[,1:9]

# Matrix to store results
mean_div <- matrix(NA, nrow = nrow(Xtest), ncol = length(interv_genes))

# Loop through each intervention gene and exclude one per iteration
for (i in seq_along(interv_genes)) {
  gene_to_remove <- interv_genes[i]
  cat("\nIteration", i, "- Removed:", gene_to_remove, "\n")
  
  valid_rows <- train_data$interventions %in% c("non-targeting", setdiff(interv_genes, gene_to_remove))
  cat("Valid rows count:", sum(valid_rows), "out of", nrow(train_data), "\n")
  
  Xtr_excl <- Xtr[valid_rows, ]
  Ytr_excl <- Ytr[valid_rows]  
  Ztr_excl <- Ztr[valid_rows, drop = TRUE]
  
  # Train the model without the removed gene
  div_mod <- div(Z = Ztr_excl, X = Xtr_excl, Y = Ytr_excl, 
                 num_layer = 4, num_epochs = 20000, lr = 1e-3, joint = TRUE,
                 epsx_dim = 50, epsh_dim = 50, epsy_dim = 50)
  
  # Make predictions
  mean_div[, i] <- predict(div_mod, Xtest = Xtest, type = "mean", nsample = 1000)
}

colnames(mean_div) <- interv_genes

# Save results
write.csv(mean_div, "results/sec6.2/div_singlecell_pairwise.csv", row.names = FALSE)