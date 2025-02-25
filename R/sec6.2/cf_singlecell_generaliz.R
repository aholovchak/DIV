library(ggplot2)
library(MASS)
library(splines)
library(tidyr)
library(dplyr)


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
Ztr <- train_data[,11]

###################
# CF NONLINEAR (ns)
###################

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

mse_envs_cf_nonlin <- data.frame(env = factor(), mse = numeric())

for(interv in test_envs$genes){
  test_data_interv <- data_pred %>% filter(interventions == interv)
  names(test_data_interv) <- c(names(data_step2), "interventions")
  mean_cf_nonlin <- predict(step2_ns, newdata = test_data_interv[1,-19])
  mse <- mean((test_data_interv[,10] - mean_cf_nonlin)^2)
  mse_envs_cf_nonlin <- rbind(mse_envs_cf_nonlin, c(interv, mse))
}

names(mse_envs_cf_nonlin) <- c("Env", "mse")

cfnonlin_mean_quantiles <- quantile(as.numeric(mse_envs_cf_nonlin$mse), c(0.0, 0.05, 0.25, 0.5, 0.75, 0.95, 1))
print(cfnonlin_mean_quantiles)

write.csv(cfnonlin_mean_quantiles, "results/sec6.2/cfnonlin_singlecell_10runs.csv", row.names = FALSE)

###################
# CF LINEAR
###################

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

mse_envs_cf_lin <- data.frame(env = factor(), mse = numeric())

for(interv in test_envs$genes){
  test_data_interv <- data_pred %>% filter(interventions == interv)
  names(test_data_interv) <- c(names(data_step2), "interventions")
  mean_cf_lin <- predict(step2_lin, newdata = test_data_interv[1,-19])
  mse <- mean((test_data_interv[,10] - mean_cf_lin)^2)
  mse_envs_cf_lin <- rbind(mse_envs_cf_lin, c(interv, mse))
}

names(mse_envs_cf_lin) <- c("Env", "mse")

cflin_mean_quantiles <- quantile(as.numeric(mse_envs_cf_lin$mse), c(0.0, 0.05, 0.25, 0.5, 0.75, 0.95, 1))
print(cflin_mean_quantiles)

write.csv(cflin_mean_quantiles, "results/sec6.2/cflin_singlecell_10runs.csv", row.names = FALSE)
