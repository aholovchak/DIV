library(foreach)
library(ggplot2)
library(splines)
library(MASS)
library(readr)

setwd("~/Documents/ethz/DIV/r_package/DistributionIV")
library(devtools)
load_all()

div_erm <- function(data_train, data_test, num_layer, nsample, num_epochs){
  # DIV ATE
  div_mod <- div(X = data_train$X, Z = data_train$Z, Y = data_train$Y,
                 epsx_dim = 50, epsh_dim = 50, epsy_dim = 50,
                 num_epochs = num_epochs, num_layer = num_layer)
  div_ate <- predict.DIV(div_mod, Xtest = data_test$Xint, type = "mean", nsample = 1000)
  div_ate_grid <- predict.DIV(div_mod, Xtest = data_test$Xtest_grid, type = "mean", nsample = 1000)
  div_sample <- predict.DIV(div_mod, Xtest = data_test$Xint, type = "sample", nsample = 1)
  return(list(div_ate = div_ate, div_ate_grid = div_ate_grid, div_sample = div_sample))
}

# List of train and test data files
train_files <- c(
  "data/sec5.4/train_Zcont_g_fct_f_softplusalpha0.csv",
  "data/sec5.4/train_Zcont_g_fct_f_softplusalpha1.csv",
  "data/sec5.4/train_Zcont_g_fct_f_softplusalpha5.csv"
)

test_files <- c(
  "data/sec5.4/test_Zcont_g_fct_f_softplusalpha0.csv",
  "data/sec5.4/test_Zcont_g_fct_f_softplusalpha1.csv",
  "data/sec5.4/test_Zcont_g_fct_f_softplusalpha5.csv"
)

# Path for saving output
base_path <- "results/sec5.4/"

# Number of runs for each setting
runs <- 10
alpha <- c(0, 1, 5)
# Loop through each train-test pair
for (i in 1:length(train_files)) {
  # Read the training and test data for the current setting
  data_train <- read_csv(train_files[i])
  data_test <- read_csv(test_files[i])
  
  # Run the simulation 'runs' times for the current setting
  sim_results <- foreach(j = 1:runs, .packages = c("MASS", "readr"), .options.RNG = 123) %do% {
    div_erm(data_train, data_test, num_layer = 4, nsample = 1000, num_epochs = 10000)
  }
  
  # Aggregate the results for the current setting
  div_ate_list <- lapply(sim_results, function(x) x$div_ate)
  div_ate_df <- do.call(cbind, div_ate_list)
  
  div_ate_grid_list <- lapply(sim_results, function(x) x$div_ate_grid)
  div_ate_grid_df <- do.call(cbind, div_ate_grid_list)
  
  # Name the columns
  col_names <- paste("Run", 1:runs, sep="_")
  names(div_ate_df) <- col_names
  names(div_ate_grid_df) <- col_names
  
  # Save the results
  write.csv(div_ate_df, paste0(base_path, "div_result_mse_alpha", alpha[i], ".csv"), row.names = FALSE)
  #write.csv(div_ate_grid_df, paste0(base_path, "div_result_plot_alpha", alpha[i], ".csv"), row.names = FALSE)
}

