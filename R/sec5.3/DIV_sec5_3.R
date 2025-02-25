# Load necessary libraries
library(readr)

setwd("R/sec5.3/DIV_flin")
library(devtools)
load_all()

setwd("../../..")
base_path <- "data/sec5.3/"
file_sizes <- c(1000, 10000)  
nruns <- 10
target_beta <- c(1, 2)

for (file_size in file_sizes) {
  train_file <- paste0(base_path, "train_Zbin_g_mult_g_log_f_lin", file_size, ".csv")
  test_file <- paste0(base_path, "test_Zbin_g_mult_g_log_f_lin", file_size, ".csv")
  
  gen_f_params_list <- vector("list", nruns)
  gen_f_beta_list <- vector("list", nruns)
  beta_norm <- numeric(nruns)
  
  data_train <- read_csv(train_file)
  data_test <- read_csv(test_file)
  
  for (i in 1:nruns) {
    set.seed(1511 + i)
    
    div_mod <- div(X = data.frame(data_train$X1, data_train$X2), Z = data_train$Z, Y = data_train$Y,
                   num_epochs = 5000, num_layer = 4, lr = 1e-3,
                   epsx_dim = 50, epsh_dim = 50, epsy_dim = 50, mod_f_lin = TRUE
                   )
    div_ate <- predict.DIV(div_mod, Xtest = data.frame(data_test$Xint1, data_test$Xint2), 
                           type = "mean", nsample = 1000)
    estimated_beta <- as.numeric(div_mod$lin_coef)[1:2]
    beta_diff <- estimated_beta - target_beta
    beta_norm[i] <- sqrt(sum(beta_diff^2))
  }
  
  # Save results to a file
  output_file <- paste0("ui_DIV_flin_n1e", log10(file_size), ".csv")
  setwd("results/sec5.3")
  write.csv(data.frame(beta_norm = beta_norm), output_file, row.names = FALSE)
}
