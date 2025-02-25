library(ggplot2)
library(MASS)
library(splines)
library(engression)
library(tidyr)
library(dplyr)

# Load datasets
dataset_rpe1 <- read.csv("data/sec6.2/dataset_rpe1.csv", stringsAsFactors=TRUE)
test_envs <- read.csv("data/sec6.2/single_cell_test_envs.csv", stringsAsFactors=TRUE)

#########################
# Identified Case - Control Function NONLINEAR (ns) and LINEAR
#########################

interv_genes <- colnames(dataset_rpe1)[1:9]  # Select first 9 columns as intervention genes

# Filter training data
train_data <- dataset_rpe1 %>% filter(interventions %in% c(interv_genes, "non-targeting"))
train_data$interventions <- droplevels(train_data$interventions)

Xtr_all <- train_data[, 1:9]   # First 9 columns
Ytr_all <- train_data[, 10]    # 10th column (target)
Ztr_all <- train_data[, 11]    # 11th column (environmental variables)

test_data <- dataset_rpe1 %>% filter(interventions %in% test_envs$genes)
Xtest <- test_data[,1:9]

# Matrices to store predictions
mean_cf_nonlin <- matrix(NA, nrow = nrow(test_data), ncol = length(interv_genes))
mean_cf_lin <- matrix(NA, nrow = nrow(test_data), ncol = length(interv_genes))

# Loop through each intervention gene and exclude one per iteration
for (i in seq_along(interv_genes)) {
  gene_to_remove <- interv_genes[i]
  cat("\nIteration", i, "- Removed:", gene_to_remove, "\n")
  
  valid_rows <- train_data$interventions %in% c("non-targeting", setdiff(interv_genes, gene_to_remove))
  cat("Valid rows count:", sum(valid_rows), "out of", nrow(train_data), "\n")
  
  Xtr <- Xtr_all[valid_rows, , drop = FALSE]
  Ytr <- Ytr_all[valid_rows]
  Ztr <- Ztr_all[valid_rows]
  
  ##############
  ### CF LINEAR
  ##############
  
  # Step 1: X ~ I
  step1_X1 <- lm(Xtr[,1] ~ Ztr)
  step1_X2 <- lm(Xtr[,2] ~ Ztr)
  step1_X3 <- lm(Xtr[,3] ~ Ztr)
  step1_X4 <- lm(Xtr[,4] ~ Ztr)
  step1_X5 <- lm(Xtr[,5] ~ Ztr)
  step1_X6 <- lm(Xtr[,6] ~ Ztr)
  step1_X7 <- lm(Xtr[,7] ~ Ztr)
  step1_X8 <- lm(Xtr[,8] ~ Ztr)
  step1_X9 <- lm(Xtr[,9] ~ Ztr)
  
  resid_X1 <- Xtr[,1] - predict(step1_X1)
  resid_X2 <- Xtr[,2] - predict(step1_X2)
  resid_X3 <- Xtr[,3] - predict(step1_X3)
  resid_X4 <- Xtr[,4] - predict(step1_X4)
  resid_X5 <- Xtr[,5] - predict(step1_X5)
  resid_X6 <- Xtr[,6] - predict(step1_X6)
  resid_X7 <- Xtr[,7] - predict(step1_X7)
  resid_X8 <- Xtr[,8] - predict(step1_X8)
  resid_X9 <- Xtr[,9] - predict(step1_X9)
  
  data_step2 <- cbind(Xtr, resid_X1, resid_X2,
                      resid_X3, resid_X4,
                      resid_X5, resid_X6,
                      resid_X7, resid_X8, resid_X9)
  
  # Step 2: Y ~ X + resid
  step2_lin <- lm(Ytr ~ ENSG00000187514 +
                    ENSG00000075624 +
                    ENSG00000147604 +
                    ENSG00000110700 +
                    ENSG00000172757 +
                    ENSG00000133112 +
                    ENSG00000067225 +
                    ENSG00000108518 +
                    ENSG00000125691 +
                    resid_X1 + resid_X2 + resid_X3 + resid_X4 +
                    resid_X5 + resid_X6 + resid_X7 + resid_X8 +
                    resid_X9, data = data_step2)
  
  data_pred <- data.frame(cbind(test_data[,1:9],
                                resid_X1 = rep(0, length(test_data[,1])),
                                resid_X2 = rep(0, length(test_data[,1])),
                                resid_X3 = rep(0, length(test_data[,1])),
                                resid_X4 = rep(0, length(test_data[,1])),
                                resid_X5 = rep(0, length(test_data[,1])),
                                resid_X6 = rep(0, length(test_data[,1])),
                                resid_X7 = rep(0, length(test_data[,1])),
                                resid_X8 = rep(0, length(test_data[,1])),
                                resid_X9 = rep(0, length(test_data[,1])),
                                interventions = test_data$interventions))

  mean_cf_lin[, i] <- predict(step2_lin, newdata = data_pred)
  
  
  ##############
  ### CF nonLINEAR
  ##############
  
  # Step 1: same as CF linear
  
  # Step 2: Y ~ X + resid
  step2_ns <- lm(Ytr ~ ns(ENSG00000187514, df = 5) +
                   ns(ENSG00000075624, df = 5) +
                   ns(ENSG00000147604, df = 5) +
                   ns(ENSG00000110700, df = 5) +
                   ns(ENSG00000172757, df = 5) +
                   ns(ENSG00000133112, df = 5) +
                   ns(ENSG00000067225, df = 5) +
                   ns(ENSG00000108518, df = 5) +
                   ns(ENSG00000125691, df = 5) +
                   ns(resid_X1, df = 5) + ns(resid_X2, df = 5) + ns(resid_X3, df = 5) + ns(resid_X4, df = 5) +
                   ns(resid_X5, df = 5) + ns(resid_X6, df = 5) + ns(resid_X7, df = 5) + ns(resid_X8, df = 5) +
                   ns(resid_X9, df = 5), data = data_step2)
  
  # data_pred: same as CF linear
  
  mean_cf_nonlin[, i] <- predict(step2_ns, newdata = data_pred)
  
}

# Assign column names to results
colnames(mean_cf_nonlin) <- interv_genes
colnames(mean_cf_lin) <- interv_genes

# Save results
write.csv(mean_cf_nonlin, "results/sec6.2/cfnonlin_singlecell_pairwise.csv", row.names = FALSE)
write.csv(mean_cf_lin, "results/sec6.2/cflin_singlecell_pairwise.csv", row.names = FALSE)