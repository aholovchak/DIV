library(ggplot2)

setwd("data/sec5.3")

simulation <- function(g1, g2, f, nsample, instrument = c("cont", "bin")){
  # Observational data
  n_train <- nsample
  n_test <- 5000
  eps_X1 <- rnorm(n_train)
  eps_X2 <- rnorm(n_train)

  eps_Y <- rnorm(n_train)

  H <- rnorm(n_train)

  if(instrument == "bin"){
    Z <- rbinom(1, n = n_train, prob = 0.5)
    Zint <- rbinom(1, n = n_test, prob = 0.5)
  }
  else if(instrument == "cont"){
    Z <- runif(n_train, 0, 3)
    Zint <- runif(n_test, 0, 3)
  }

  X1 <- g1(Z = Z, H = H, eps_X = eps_X1)
  X2 <- g2(Z = Z, H = H, eps_X = eps_X2)
  Y <- f(X1 = X1, X2 = X2, H = H, eps_Y = eps_Y)

  # Interventional data
  n_test <- 5000

  eps_X1int <- rnorm(n_test)
  eps_X2int <- rnorm(n_test)

  eps_Yint <- rnorm(n_test)

  H11int <- rnorm(n_test)
  H12int <- rnorm(n_test)
  
  H2int <- rnorm(n_test)

  Xint1 <- g1(Z = Zint, H = H11int, eps_X = eps_X1int)
  Xint2 <- g2(Z = Zint, H = H12int, eps_X = eps_X2int)

  Yint <- f(X1 = Xint1, X2 = Xint2, H = H2int, eps_Y = eps_Yint)

  # Assuming you have a data frame with your data
  # Create a data frame combining both X and Xint along with their corresponding Y and Yint
  data <- data.frame(
    X1 = c(X1, Xint1),
    X2 = c(X2, Xint2),
    Y = c(Y, Yint),
    Group = factor(rep(c("Observational Data", "Interventional Data"), times = c(length(X1), length(Xint1))))
  )

  # Create the ggplot
  p1 <- ggplot(data, aes(x = X1, y = Y, color = Group)) +
    geom_point(alpha = 0.5) +
    labs(x = "X1", y = "Y") +
    theme_minimal() +
    guides(color = guide_legend(title = ""))

  print(p1)

  p2 <- ggplot(data, aes(x = X2, y = Y, color = Group)) +
    geom_point(alpha = 0.5) +
    labs(x = "X2", y = "Y") +
    theme_minimal() +
    guides(color = guide_legend(title = ""))

  print(p2)

  # True mean (interventional data, for MSE)

  # Number of points to sample
  num_points <- 5000
  num_rep <- 5000
  
  # Define the range for x1 and x2
  x1_range <- seq(min(Xint1), max(Xint1), length.out = sqrt(num_points))
  x2_range <- seq(min(Xint2), max(Xint2), length.out = sqrt(num_points))
  
  # Create a grid of x1 and x2
  grid <- expand.grid(x1 = x1_range, x2 = x2_range)
  
  # Limit the grid to exactly num_points if needed
  grid <- grid[1:num_points, ]
  
  # Initialize a vector to store the results
  mean_int <- numeric(nrow(grid))
  
  # Loop over each row of the grid
  for (i in 1:nrow(grid)) {
    # Sample new H and eps_Y vectors
    eps_Ymean <- rnorm(num_rep)
    Hmean <- rnorm(num_rep)
    
    # Compute Ymean for the current x1 and x2 values
    Ymean <- f(X1 = grid$x1[i],
              X2 = grid$x2[i],
              H = Hmean, eps_Y = eps_Ymean)
    
    # Compute the trimmed mean and store it
    mean_value <- mean(Ymean, trim = 0.05)
    mean_int[i] <- mean_value
  }

  # Grid X (for plots)
  X1test_grid <- seq(min(X1), max(X1), length.out = 5000)
  X2test_fix <- median(X2)

  # True mean (interventional data, for plot)

  # Initialize an empty vector to store the mean values
  mean_grid <- numeric(length(X1test_grid))
  num_rep <- 5000

  # Loop over each element in Xtest
  for (i in 1:length(X1test_grid)) {
    # Sample new H and eps_Y vectors for each Xtest value
    eps_Ymean <- rnorm(num_rep)
    Hmean <- rnorm(num_rep)

    # Compute the mean for the current Xtest value
    Ymean_grid <- f(X1 = X1test_grid[i], X2 = X2test_fix,
              H = Hmean, eps_Y = eps_Ymean)
    mean_grid <- mean(Ymean_grid, trim = 0.05)

    # Store the mean value in the mean_vector
    mean_grid[i] <- mean_grid
  }

  g1_name <- deparse(substitute(g1))
  g2_name <- deparse(substitute(g2))

  f_name <- deparse(substitute(f))
  train_df_name <- paste("train_Z", instrument, "_", g1_name, "_", g2_name, "_", f_name, nsample, ".csv", sep = "")
  write.csv(data.frame(Z, X1, X2, Y), train_df_name, row.names = FALSE)

  test_df_name <- paste("test_Z", instrument, "_", g1_name, "_", g2_name, "_", f_name, nsample, ".csv", sep = "")
  write.csv(data.frame(Xint1 = grid$x1, Xint2 = grid$x2, mean_int, X1test_grid, X2test_fix, mean_grid), test_df_name, row.names = FALSE)
  return(data = data.frame(Xint1, Xint2, mean_int, X1test_grid, X2test_fix, mean_grid))
}

######
######
# Define functions g and f

g_log <- function(Z, H, eps_X){
  return(log(2 * Z + H + eps_X + 7))
}

g_mult <- function(Z, H, eps_X){
  return(Z * (2 * H - 0.5 * eps_X))
}

f_lin <- function(X1, X2, H, eps_Y){
  return(1*X1 + 2*X2 + 2*H + eps_Y)
}


######
######

set.seed(1998)
result <- simulation(g1 = g_mult, g2 = g_log, f = f_lin, instrument = "bin",
                     nsample = 1000)
# result <- simulation(g1 = g_mult, g2 = g_log, f = f_lin, instrument = "bin",
#                      nsample = 10000)

