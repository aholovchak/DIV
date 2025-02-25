library(foreach)
library(doParallel)
library(ggplot2)
library(splines)
library(MASS)
library(readr)

setwd("~/Documents/ethz/DIV/r_package/DistributionIV")
library(devtools)
load_all()

div_erm <- function(data_train, data_test, num_layer, nsample, num_epochs){
  # DIV interventional mean
  div_mod <- div(X = data_train$X, Z = data_train$Z, Y = data_train$Y,
                 epsx_dim = 50, epsy_dim = 50, epsh_dim = 50,
                 num_epochs = num_epochs, num_layer = num_layer, lr = 1e-3)
  div_mean <- predict.DIV(div_mod, Xtest = data_test$Xint, type = "mean", nsample = 1000)
  #
  mean((div_mean - data_test$mean_int)^2)
  
  div_mean_grid <- predict.DIV(div_mod, Xtest = data_test$Xtest_grid, type = "mean", nsample = 1000)
  div_sample <- predict.DIV(div_mod, Xtest = data_test$Xint, type = "sample", nsample = 1)
  # Clean up
  div_mod <- NULL
  gc()  # Run garbage collection to free memory
  
  return(list(div_mean = div_mean, div_mean_grid = div_mean_grid, div_sample = div_sample))
}

# List of train and test data files
train_files <- c(
  "data/sec5.1/train_Zbin_g_lin_f_lin.csv",
  "data/sec5.1/train_Zcont_g_lin_f_lin.csv",
  "data/sec5.1//train_Zbin_g_lin_f_log_case.csv",
  "data/sec5.1/train_Zcont_g_lin_f_log_case.csv",
  "data/sec5.1/train_Zbin_g_lin_f_sin_lin.csv",
  "data/sec5.1/train_Zcont_g_lin_f_sin_lin.csv"
)

test_files <- c(
  "data/sec5.1/test_Zbin_g_lin_f_lin.csv",
  "data/sec5.1/test_Zcont_g_lin_f_lin.csv",
  "data/sec5.1/test_Zbin_g_lin_f_log_case.csv",
  "data/sec5.1/test_Zcont_g_lin_f_log_case.csv",
  "data/sec5.1/test_Zbin_g_lin_f_sin_lin.csv",
  "data/sec5.1/test_Zcont_g_lin_f_sin_lin.csv"
)

# Path for saving output
base_path <- "results/sec5.1"

# Number of runs for each setting
runs <- 10

# Loop through each train-test pair
for (i in 1:length(train_files)) {
# for (i in 2:2) {
  # Read the training and test data for the current setting
  data_train <- read_csv(train_files[i])
  data_test <- read_csv(test_files[i])
  
  # Run the simulation 'runs' times for the current setting
  sim_results <- foreach(j = 1:runs, .packages = c("ggplot2", "splines", "MASS", "engression", "readr"), .options.RNG = 123) %do% {
    div_erm(data_train, data_test, num_layer = 4, nsample = 5000, num_epochs = 10000)
  }
  
  div_mean_list <- lapply(sim_results, function(x) x$div_mean)
  div_mean_df <- do.call(cbind, div_mean_list)
  
  div_mean_grid_list <- lapply(sim_results, function(x) x$div_mean_grid)
  div_mean_grid_df <- do.call(cbind, div_mean_grid_list)
  
  # Name the columns
  col_names <- paste("Run", 1:runs, sep="_")
  names(div_mean_df) <- col_names
  names(div_mean_grid_df) <- col_names
  
  write.csv(div_mean_df, paste0(base_path, "div2step_result_mse_dgp", i, ".csv"), row.names = FALSE)
  write.csv(div_mean_grid_df, paste0(base_path, "div2step_result_plot_dgp", i, ".csv"), row.names = FALSE)
}

