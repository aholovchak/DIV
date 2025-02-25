library(tidyverse)

# Load test data
data_test_alpha0 <- read_csv("data/test_Zbin_g_fct_f_softplusalpha0.csv")
data_test_alpha1 <- read_csv("data/test_Zbin_g_fct_f_softplusalpha1.csv")
data_test_alpha5 <- read_csv("data/test_Zbin_g_fct_f_softplusalpha5.csv")

# Define methods and base path
methods <- c("div", "cflin", "cfnonlin", "deepiv", "deepgmm", "hsic")
base_path <- "results/"

# Function to compute mean MSE
compute_mean_mse <- function(alpha, data_test) {
  mse_data <- lapply(methods, function(method) {
    mse_file <- paste0(base_path, method, "_result_mse_alpha", alpha, ".csv")
    mse_df <- read_csv(mse_file) %>% select(starts_with("V"))
    
    # Compute mean MSE over 10 runs
    mean_mse <- mean(rowMeans((mse_df - data_test$ate_int)^2))
    
    tibble(method = method, mean_mse = mean_mse)
  })
  
  bind_rows(mse_data) %>% mutate(dgp = paste0("DGP", dgp))
}

# Compute mean MSE for alpha \in {0,1,5}
results_mse <- bind_rows(
  compute_mean_mse(0, data_test_dgp1),
  compute_mean_mse(1, data_test_dgp2),
  compute_mean_mse(5, data_test_dgp3)
)

# Display results
print(results_mse)