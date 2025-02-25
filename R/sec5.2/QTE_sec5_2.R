### Simulation DIVE
### LK 2024

library(reticulate)
use_virtualenv("r-tensorflow", required = TRUE)

library(tidyverse)
library(tram)
library(IVQR)
library(dare)

set.seed(12)

# choose scenario '1' or '2'
scenario <- 1

# Settings ----------------------------------------------------------------

dgp <- switch(
  scenario,
  "1" = function(n = 1e2, do = FALSE) {
    Z <- rlogis(n)
    H <- rlogis(n)
    D <- as.numeric(4 * Z + (1 - do) * 4 * H > rlogis(n))
    Y <- log(1 + exp(18 + 8 * D + 6 * H))
    data.frame(Y = Y, H = H, D = D, Z = Z)
  },
  "2" = function(n = 1e2, do = FALSE) {
    Z <- rlogis(n)
    H <- rlogis(n)
    D <- as.numeric(4 * Z + (1 - do) * 4 * H > rlogis(n))
    Y <- 2 + (D + 1)^2 + 3 * (D + 1) + (2 * H + rlogis(n))
    data.frame(Y = Y, H = H, D = D, Z = Z)
  }
)


# Number of repetitions
Nrep <- 10

# Quantiles for tau
taus <- c(0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 0.99)

f <- function(y, cdf, newx) {
  s <- spline(x = y, y = cdf, method = "natural")
  # Ensure newx is within the range of s$y
  newx <- pmin(pmax(newx, min(s$y)), max(s$y))
  approx(x = s$y, y = s$x, xout = newx, yleft = min(y), yright = max(y))$y
}

# ORACLE ------------------------------------------------------------------

dint <- dgp(1e4, do = TRUE)
F0 <- ecdf(dint$Y[dint$D == 0])
F1 <- ecdf(dint$Y[dint$D == 1])
oracle <- Vectorize(\(y, d) d * F1(y) + (1 - d) * F0(y))

# Empty list to store QTE results for each method across 50 samples
pd_list <- list()


# Params DIVE ------------------------------------------------------------------

nep <- 3e3
wep <- 3e3
rep <- 1
ns <- 10^(2)

tlr <- 0.1
tord <- 50
iou <- 1
alp <- 0.1
tss <- 10

# Run the simulation 10 times for each method
for (iter in 1:Nrep) {
  
  # Generate training data for each sample
  set.seed(12 + iter)
  data_train <- dgp(1000)  # Assuming dgp is defined in the context
  
  # ### DIVE Method ###

  # Fit the DIVE model (assuming BoxCoxNN and fit_adaptive are defined)
  mtmp <- BoxCoxNN(Y | D ~ 1, data = data_train, order = 50,
                   optimizer = optimizer_adam(1e-2), tf_seed = iter)
  fit(mtmp, epochs = 3000, validation_split = 0, callbacks = list(
    callback_reduce_lr_on_plateau("loss", factor = 0.9, patience = 20, min_delta = 1e-3),
    callback_early_stopping("loss", patience = 60, min_delta = 1e-3)), verbose = FALSE)
  tmp <- get_weights(mtmp$model)

  # Fit DIVE model
  setwd("~/Documents/ethz/DIV/r_package/dive-master")
  devtools::load_all()

  args <- list(formula = Y | D ~ 1, data = data_train, anchor = ~ Z,
               loss = "indep", xi = 1, trafo_options = trafo_control(order_bsp = 50, support = range(data_train$Y)))
  cb <- \() list(callback_reduce_lr_on_plateau(
    "loss", patience = 20, factor = 0.9, min_delta = 1e-4),
    callback_early_stopping("loss", patience = 60, min_delta = 1e-4))
  m <- fit_adaptive(args, 3000, max_iter = 10, ws = tmp, modFUN = "BoxCoxDA", verbose = FALSE, lr = 0.1,
                    cb = list(callback_reduce_lr_on_plateau("loss", patience = 20, factor = 0.9, min_delta = 1e-4),
                              callback_early_stopping("loss", patience = 60, min_delta = 1e-4)), start_xi = TRUE, stepsize = 10,
                    indep_over_unif = 1, alpha = 0.1)

  # Predict CDF using DIVE
  DIVE_vals <- c(predict(m, type = "cdf"))

  # Generate a grid of D and Y for QTE calculation
  nd <- expand_grid(D = sort(unique(data_train$D)), Y = unique(data_train$Y))

  # Store DIVE QTE
  nd$DIVE <- c(predict(m, type = "cdf", newdata = data.frame(Y = nd$Y, D = nd$D)))
  nd$ORACLE <- oracle(nd$Y, nd$D)  # Assuming oracle function is defined

  # Reshape data for plotting
  pd_dive <- nd %>%
    pivot_longer(DIVE:ORACLE, names_to = "method", values_to = "tau") %>%
    group_by(D, method) %>%
    mutate(ntau = seq(0, 1, length.out = n())) %>%
    mutate(ny = f(Y, tau, ntau)) %>%
    select(D, ntau, ny, method) %>%
    pivot_wider(names_from = "D", values_from = "ny") %>%
    mutate(qte = `1` - `0`) %>%
    select(ntau, method, qte)
  pd_dive$iter <- iter

  ### DIV ###

  setwd("DistributionIV")
  devtools::load_all()

  div_mod <- div(X = data_train$D, Z = data_train$Z, Y = data_train$Y,
                       epsx_dim = 50, epsh_dim = 50, epsy_dim = 50,
                       num_epochs = 20000, lr=1e-4, num_layer = 4)
  q0_pred <- predict.DIV(div_mod, Xtest = rep(0, 1), type = "quantile", quantiles = taus, nsample = 1000)
  q1_pred <- predict.DIV(div_mod, Xtest = rep(1, 1), type = "quantile", quantiles = taus, nsample = 1000)

  # Store DIV QTE
  div_df <- data.frame(
    ntau = taus,
    qte = c(q1_pred - q0_pred),
    method = "DIV",
    iter = iter
  )

  
### IVQR Method ###
  grd <- seq(0, 15, length.out = 1e2)
  fit <- ivqr(Y ~ D | Z | 1, taus = taus, grid = grd, data = data_train)

  # Store IVQR QTE
  ivqr_df <- data.frame(
    ntau = fit$taus,
    qte = c(fit$coef$endg_var),
    method = "IVQR",
    iter = iter
  )

  # Combine results for this iteration
  pd_list[[iter]] <- bind_rows(pd_dive,
                               div_df,
                               ivqr_df)
}

# Combine all results into one data frame
pd_combined <- bind_rows(pd_list)

# Plot 10 lines for each method
# ggplot(pd_combined, aes(x = ntau, y = qte, color = method, group = interaction(method, iter))) +
#   geom_line(alpha = 0.3) +
#   labs(title = "QTE Comparison Across Methods, 10 Samples each",
#        x = bquote(alpha),
#        y = "QTE",
#        color = "Method") +
#   theme_minimal() +
#   theme(legend.position = "top")

###

setwd("results/sec5.2/")
file_name <- paste0("qte_scenario", scenario, ".csv")
write.csv(pd_combined, file_name, row.names = FALSE)

