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

base_path <- "results/sec6.2/"

###################
# DIV
###################
div <- read.csv(paste0(base_path, "div_singlecell_10runs.csv"))

div_quantiles_list <- list()

# Iterate over each simulation run (each column in div)
for (run in 1:10) {
  mse_envs_div <- data.frame(Env = character(), mse = numeric(), stringsAsFactors = FALSE)
  for (interv in test_envs$genes) {
    test_data_interv <- test_data %>% filter(interventions == interv)
    mean_div <- div[, run]  # Take the corresponding column of div
    mse <- mean((test_data_interv[, 10] - mean_div)^2)
    mse_envs_div <- rbind(mse_envs_div, c(interv, mse))
  }
  names(mse_envs_div) <- c("Env", "mse")
  mse_envs_div$mse <- as.numeric(mse_envs_div$mse)
  div_quantiles_list[[run]] <- quantile(mse_envs_div$mse, c(0.0, 0.05, 0.25, 0.5, 0.75, 0.95, 1))
}
div_mean_quantiles <- Reduce("+", div_quantiles_list) / length(div_quantiles_list)
print(div_mean_quantiles)

###################
# Engression
###################
engr <- read.csv(paste0(base_path, "engr_singlecell_10runs.csv"))

engr_quantiles_list <- list()

# Iterate over each simulation run (each column in div)
for (run in 1:10) {
  mse_envs_engr <- data.frame(Env = character(), mse = numeric(), stringsAsFactors = FALSE)
  for (interv in test_envs$genes) {
    test_data_interv <- test_data %>% filter(interventions == interv)
    mean_engr <- engr[, run]  # Take the corresponding column of div
    mse <- mean((test_data_interv[, 10] - mean_engr)^2)
    mse_envs_engr <- rbind(mse_envs_engr, c(interv, mse))
  }
  names(mse_envs_engr) <- c("Env", "mse")
  mse_envs_engr$mse <- as.numeric(mse_envs_engr$mse)
  engr_quantiles_list[[run]] <- quantile(mse_envs_engr$mse, c(0.0, 0.05, 0.25, 0.5, 0.75, 0.95, 1))
}
engr_mean_quantiles <- Reduce("+", engr_quantiles_list) / length(engr_quantiles_list)
print(engr_mean_quantiles)

###################
# CF NONLINEAR
###################

cfnonlin_mean_quantiles <- as.vector(read.csv(paste0(base_path, "cfnonlin_singlecell_10runs.csv")))

###################
# CF LINEAR
###################

cflin_mean_quantiles <- as.vector(read.csv(paste0(base_path, "cflin_singlecell_10runs.csv")))

###################
# DeepGMM
###################
deepgmm <- read.csv(paste0(base_path, "deepgmm_singlecell_10runs.csv"))

deepgmm_quantiles_list <- list()

# Iterate over each simulation run (each column in div)
for (run in 1:10) {
  mse_envs_deepgmm <- data.frame(Env = character(), mse = numeric(), stringsAsFactors = FALSE)
  for (interv in test_envs$genes) {
    test_data_interv <- test_data %>% filter(interventions == interv)
    mean_deepgmm <- deepgmm[, run]  # Take the corresponding column of div
    mse <- mean((test_data_interv[, 10] - mean_deepgmm)^2)
    mse_envs_deepgmm <- rbind(mse_envs_deepgmm, c(interv, mse))
  }
  names(mse_envs_deepgmm) <- c("Env", "mse")
  mse_envs_deepgmm$mse <- as.numeric(mse_envs_deepgmm$mse)
  deepgmm_quantiles_list[[run]] <- quantile(mse_envs_deepgmm$mse, c(0.0, 0.05, 0.25, 0.5, 0.75, 0.95, 1))
}
deepgmm_mean_quantiles <- Reduce("+", deepgmm_quantiles_list) / length(deepgmm_quantiles_list)
print(deepgmm_mean_quantiles)

###################
# DeepIV
# not considered due to NAs in several iterations
###################

# deepiv <- read.csv(paste0(base_path, "deepIV_singlecell_10runs.csv"))
# 
# deepiv_quantiles_list <- list()
# 
# # Iterate over each simulation run (each column in div)
# for (run in 1:10) {
#   mse_envs_deepiv <- data.frame(Env = character(), mse = numeric(), stringsAsFactors = FALSE)
#   for (interv in test_envs$genes) {
#     test_data_interv <- test_data %>% filter(interventions == interv)
#     mean_deepiv <- deepiv[, run]  # Take the corresponding column of div
#     mse <- mean((test_data_interv[, 10] - mean_deepiv)^2, na.rm = TRUE)
#     mse_envs_deepiv <- rbind(mse_envs_deepiv, c(interv, mse))
#   }
#   names(mse_envs_deepiv) <- c("Env", "mse")
#   mse_envs_deepiv$mse <- as.numeric(mse_envs_deepiv$mse)
#   deepiv_quantiles_list[[run]] <- quantile(mse_envs_deepiv$mse, c(0.0, 0.05, 0.25, 0.5, 0.75, 0.95, 1), na.rm = TRUE)
# }
# deepiv_quantiles_list <- deepiv_quantiles_list[!sapply(deepiv_quantiles_list, function(x) all(is.na(x)))]
# deepiv_mean_quantiles <- Reduce("+", deepiv_quantiles_list) / length(deepiv_quantiles_list)
# print(deepiv_mean_quantiles)


###################
# HSIC-X
###################
hsicx <- read.csv(paste0(base_path, "hsicx_singlecell_10runs.csv"))

hsicx_quantiles_list <- list()

# Iterate over each simulation run (each column in div)
for (run in 1:10) {
  mse_envs_hsicx <- data.frame(Env = character(), mse = numeric(), stringsAsFactors = FALSE)
  for (interv in test_envs$genes) {
    test_data_interv <- test_data %>% filter(interventions == interv)
    mean_hsicx <- hsicx[, run]  # Take the corresponding column of div
    mse <- mean((test_data_interv[, 10] - mean_hsicx)^2, na.rm = TRUE)
    mse_envs_hsicx <- rbind(mse_envs_hsicx, c(interv, mse))
  }
  names(mse_envs_hsicx) <- c("Env", "mse")
  mse_envs_hsicx$mse <- as.numeric(mse_envs_hsicx$mse)
  hsicx_quantiles_list[[run]] <- quantile(mse_envs_hsicx$mse, c(0.0, 0.05, 0.25, 0.5, 0.75, 0.95, 1), na.rm = TRUE)
}
hsicx_quantiles_list <- hsicx_quantiles_list[!sapply(hsicx_quantiles_list, function(x) all(is.na(x)))]
hsicx_mean_quantiles <- Reduce("+", hsicx_quantiles_list) / length(hsicx_quantiles_list)
print(hsicx_mean_quantiles)


####
# RESULTS
####

# Create a data frame with method names and mean quantiles
quantiles_df <- data.frame(
  Method = c("DIV", "Engression", "CF linear", "CF nonlinear", "DeepGMM", "HSIC-X"),
  rbind(
    div_mean_quantiles,
    engr_mean_quantiles,
    cflin_mean_quantiles$x,
    cfnonlin_mean_quantiles$x,
    deepgmm_mean_quantiles,
    #deepiv_mean_quantiles,
    hsicx_mean_quantiles
  )
)

colnames(quantiles_df) <- c("method", "Q00", "Q05", "Q25", "Q50", "Q75", "Q95", "Q100")

print(quantiles_df, digits = 4, row.names = FALSE)

# library(xtable)
# print(xtable(quantiles_df, digits = 4), include.rownames=FALSE)
