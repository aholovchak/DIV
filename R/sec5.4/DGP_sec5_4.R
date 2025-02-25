library(ggplot2)

simulation <- function(g, f, nsample, instrument = c("cont", "bin"), alpha){
  # Observational data
  n_train <- nsample
  
  eps_X <- runif(n_train, -1, 1)
  eps_Y <- runif(n_train, -1, 1)
  
  H <- runif(n_train, -1, 1)
  
  if(instrument == "bin"){
    Z <- rbinom(1, n = n_train, prob = 0.5)
    Zint <- rbinom(1, n = n_train, prob = 0.5)
  }
  else if(instrument == "cont"){
    Z <- runif(n_train, -3, 3)
    Zint <- runif(n_train, -3, 3)
  }
  
  X <- g(Z = Z, H = H, alpha = alpha, eps_X = eps_X)
  Y <- f(X = X, H = H, eps_Y = eps_Y)
  
  # Interventional data
  eps_Xint <- runif(n_train, -1, 1)
  eps_Yint <- runif(n_train, -1, 1)
  
  H1int <- runif(n_train, -1, 1)
  H2int <- runif(n_train, -1, 1)
  
  Xint <- g(Z = Zint, H = H1int, alpha = alpha, eps_X = eps_Xint)
  Yint <- f(X = Xint, H = H2int, eps_Y = eps_Yint)
  
  data <- data.frame(
    X = c(X, Xint),
    Y = c(Y, Yint),
    Group = factor(rep(c("Observational Data", "Interventional Data"), each = length(X)))
  )
  
  # Create the ggplot
  p <- ggplot(data, aes(x = X, y = Y, color = Group)) +
    geom_point(alpha = 0.2) +
    labs(x = "X", y = "Y") +
    theme_minimal() +
    guides(color = guide_legend(title = "")) +
    scale_color_manual(values = c("Observational Data" = "deepskyblue", "Interventional Data" = "darkgoldenrod1"))
  
  print(p)
  
  # True mean (interventional data, for MSE)
  
  # Initialize an empty vector to store the mean values
  mean_int <- numeric(length(Xint))
  num_rep <- 5000
  
  # Loop over each element in Xtest
  for (i in 1:length(Xint)) {
    # Sample new H and eps_Y vectors for each Xtest value
    
    eps_Ymean <- runif(num_rep, -1, 1)
    Hmean <- runif(num_rep, -1, 1)
    
    # Compute the mean for the current Xtest value
    Ymean <- f(X = rep(Xint[i], num_rep), H = Hmean, eps_Y = eps_Ymean)
    mean_value <- mean(Ymean, trim = 0.05)
    
    # Store the mean value in the mean_vector
    mean_int[i] <- mean_value
  }
  
  # Grid X (for plots)
  Xtest_grid <- seq(min(X), max(X), length.out = 5000)
  
  # True mean (interventional data, for plot)
  
  # Initialize an empty vector to store the mean values
  mean_grid <- numeric(length(Xtest_grid))
  num_rep <- 5000
  
  # Loop over each element in Xtest
  for (i in 1:length(Xtest_grid)) {
    # Sample new H and eps_Y vectors for each Xtest value
    eps_Ymean <- runif(num_rep, -1, 1)
    Hmean <- runif(num_rep, -1, 1)
    
    # Compute the mean for the current Xtest value
    Ymean <- f(X = rep(Xtest_grid[i], num_rep), H = Hmean, eps_Y = eps_Ymean)
    mean_value <- mean(Ymean, trim = 0.05)
    
    # Store the mean value in the mean_vector
    mean_grid[i] <- mean_value
  }
  Ysample_grid <- f(X = Xtest_grid, H = runif(num_rep, -1, 1), eps_Y = runif(num_rep, -1, 1))
  
  g_name <- deparse(substitute(g))
  f_name <- deparse(substitute(f))
  train_df_name <- paste("data/sec5.4/", "train_Z", instrument, "_", g_name, "_", f_name, "alpha", alpha, ".csv", sep = "")
  write.csv(data.frame(Z, X, Y), train_df_name, row.names = FALSE)
  
  test_df_name <- paste("data/sec5.4/", "test_Z", instrument, "_", g_name, "_", f_name, "alpha", alpha, ".csv", sep = "")
  write.csv(data.frame(Xint, mean_int, Xtest_grid, mean_grid), test_df_name, row.names = FALSE)
  return(data = data.frame(Xint, mean_int, Xtest_grid, mean_grid))
}

#########################################
## E(X|Z) = \alpha*Z
#########################################

g_fct <- function(Z, H, eps_X, alpha){
  X <- 0.5 * Z * (3 * H + eps_X + alpha)
}

f_softplus <- function(X, H, eps_Y){
  return(0.5 * log(1 + exp(2 * X - 3 * H + eps_Y)))
}

######
######

# \alpha \in {0,1,5}
set.seed(1998)
result <- simulation(g = g_fct, f = f_softplus, alpha = 0, instrument = "cont",
                     nsample = 10000)
