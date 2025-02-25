library(ggplot2)
library(MASS)
library(engression)
library(tidyr)
library(dplyr)

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

test_data <- dataset_rpe1 %>% filter(interventions %in% test_envs$genes)
Xtest <- test_data[,1:9]

# Matrix to store results
mean_engression <- matrix(NA, nrow = nrow(Xtest), ncol = length(interv_genes))

# Loop through each intervention gene and exclude one per iteration
for (i in seq_along(interv_genes)) {
  gene_to_remove <- interv_genes[i]
  cat("\nIteration", i, "- Removed:", gene_to_remove, "\n")
  
  valid_rows <- train_data$interventions %in% c("non-targeting", setdiff(interv_genes, gene_to_remove))
  cat("Valid rows count:", sum(valid_rows), "out of", nrow(train_data), "\n")
  
  Xtr_excl <- Xtr[valid_rows, , drop = FALSE]
  Ytr_excl <- Ytr[valid_rows]  

  # Train the model without the removed gene
  engr_mod <- engression(X = Xtr_excl, Y = Ytr_excl, 
                 num_layer = 4, num_epochs = 20000, lr = 1e-3, noise_dim = 100)
  
  # Make predictions
  mean_engression[, i] <- predict(engr_mod, Xtest = Xtest, type = "mean", nsample = 1000)
}

colnames(mean_engression) <- interv_genes

# Save results
write.csv(mean_engression, "results/sec6.2/engression_singlecell.csv", row.names = FALSE)