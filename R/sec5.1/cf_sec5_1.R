library(splines)
library(MASS)
library(readr)

cf_lin_nonlin <- function(data_train, data_test, df = 5){
  
  # CF (linear) -----------------------------------------------------------
  # Step 1: X ~ Z
  step1_lin <- lm(X ~ Z, data = data_train)
  resid_lin <- data_train$X - predict(step1_lin)
  
  # Step 2: Y ~ X + resid
  step2_lin <- lm(Y ~ X + resid_lin, data = data_train)
  
  # Control functions (nonlinear) -------------------------------------------
  # Step 1: X ~ Z
  step1_nonlin <- lm(X ~ ns(Z, df = df), data = data_train)
  resid_nonlin <- data_train$X - predict(step1_nonlin)
  
  # Step 2: Y ~ X + resid
  step2_nonlin <- lm(Y ~ ns(X, df = df) + resid_nonlin, data = data_train)
  
  # interventional mean --------------------------------------------------------
  Y_hat_cf_lin <- predict(step2_lin, newdata = data.frame(X = data_test$Xint, resid_lin = rep(0, length(data_test$Xint))))
  Y_hat_cf_nonlin <- predict(step2_nonlin, newdata = data.frame(X = data_test$Xint, resid_nonlin = rep(0, length(data_test$Xint))))
  
  Y_hat_cf_lin_grid <- predict(step2_lin, newdata = data.frame(X = data_test$Xtest_grid, resid_lin = rep(0, length(data_test$Xtest_grid))))
  Y_hat_cf_nonlin_grid <- predict(step2_nonlin, newdata = data.frame(X = data_test$Xtest_grid, resid_nonlin = rep(0, length(data_test$Xtest_grid))))
  
  return(list(Y_hat_cf_lin = Y_hat_cf_lin, Y_hat_cf_lin_grid = Y_hat_cf_lin_grid,
              Y_hat_cf_nonlin = Y_hat_cf_nonlin, Y_hat_cf_nonlin_grid = Y_hat_cf_nonlin_grid))
}

# List of train and test data files
train_files <- c(
  "data/sec5.1/train_Zbin_g_lin_f_lin.csv",
  "data/sec5.1/train_Zcont_g_lin_f_lin.csv",
  "data/sec5.1/train_Zbin_g_lin_f_log_case.csv",
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
base_path <- "results/sec5.1/"

# Number of runs for each setting
runs <- 10

# The function `cf_lin_nonlin()` is deterministic for the same data,
# meaning that every run produces identical results. However, we
# repeat it `runs` times to maintain consistency with other methods.

for (i in 1:length(train_files)) {
  data_train <- read_csv(train_files[i])
  data_test <- read_csv(test_files[i])
  
  cflin_df_mse <- data.frame(matrix(nrow = 5000, ncol = runs))
  cflin_df_plot <- data.frame(matrix(nrow = 5000, ncol = runs))
  cfnonlin_df_mse <- data.frame(matrix(nrow = 5000, ncol = runs))
  cfnonlin_df_plot <- data.frame(matrix(nrow = 5000, ncol = runs))
  
  set.seed(1998)
  
  for(j in 1:runs){
    temp_res <- cf_lin_nonlin(data_train, data_test, df = 5)
    cflin_df_mse[,j] <- temp_res$Y_hat_cf_lin
    cflin_df_plot[,j] <- temp_res$Y_hat_cf_lin_grid
    
    cfnonlin_df_mse[,j] <- temp_res$Y_hat_cf_nonlin
    cfnonlin_df_plot[,j] <- temp_res$Y_hat_cf_nonlin_grid
  }
  
  # Plot (optional)
  # plot(data_test$Xtest_grid, data_test$ate_grid)
  # lines(data_test$Xtest_grid, cflin_df_plot$X1, col = "red", lwd = 2)
  # lines(data_test$Xtest_grid, cfnonlin_df_plot$X1, col = "green", lwd = 2)
  
  write.csv(cflin_df_mse, paste0(base_path, "cflin_result_mse_dgp", i, ".csv"), row.names = FALSE)
  write.csv(cflin_df_plot, paste0(base_path, "cflin_result_plot_dgp", i, ".csv"), row.names = FALSE)
  write.csv(cfnonlin_df_mse, paste0(base_path, "cfnonlin_result_mse_dgp", i, ".csv"), row.names = FALSE)
  write.csv(cfnonlin_df_plot, paste0(base_path, "cfnonlin_result_plot_dgp", i, ".csv"), row.names = FALSE)
}