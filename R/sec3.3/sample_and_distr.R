library(ggplot2)
library(MASS)
library(readr)
library(tidyr)
library(dplyr)

setwd("/Users/aholovchak/Documents/ethz/DIV/r_package/DistributionIV")
library(devtools)
load_all()

data_train <- read_csv("data/sec3.3/train_Zcont_g_softplus_f_softplus.csv")
data_test <- read_csv("data/sec3.3/test_Zcont_g_softplus_f_softplus.csv")

# Train DIV model & predict
div_mod <- div(X = data_train$X, Z = data_train$Z, Y = data_train$Y, num_epochs = 5000, num_layer = 4)
div_sample <- predict.DIV(div_mod, Xtest = data_test$Xint, type = "sample", nsample = 1)

# TRUE function
f_softplus <- function(X, H, eps_Y){
  return(0.5 * log(1 + exp(2 * X + 2 * H + 4 * eps_Y)))
}
f <- f_softplus
Ytest_Xint <- f(X = data_test$Xint, H = rnorm(1000), eps_Y = rnorm(1000))


# Sample from interventional distribution

true_sample_q25 <- f(X = quantile(data_test$Xint, 0.25), H = rnorm(1000),
                     eps_Y = rnorm(1000))
true_sample_q50 <- f(X = quantile(data_test$Xint, 0.50), H = rnorm(1000),
                     eps_Y = rnorm(1000))
true_sample_q75 <- f(X = quantile(data_test$Xint, 0.75), H = rnorm(1000),
                     eps_Y = rnorm(1000))

div_sample_q25 <- predict.DIV(div_mod, Xtest = quantile(data_test$Xint, 0.25), type = "sample", 
                              nsample = 1000)
div_sample_q50 <- predict.DIV(div_mod, Xtest = quantile(data_test$Xint, 0.50), type = "sample", 
                              nsample = 1000)
div_sample_q75 <- predict.DIV(div_mod, Xtest = quantile(data_test$Xint, 0.75), type = "sample", 
                              nsample = 1000)

# Plot sample
plot_data_1 <- data.frame(Xint = data_test$Xint, Ytest_Xint = data_test$Yint, Div_Sample = div_sample)

ggplot(plot_data_1) +
  geom_point(data = plot_data_1, aes(x = Xint, y = Ytest_Xint, color = "Ytest_Xint"), alpha = 0.5) +
  geom_point(data = plot_data_1, aes(x = Xint, y = Div_Sample, color = "Div_Sample"), alpha = 0.5) +
  theme_minimal() +
  labs(x = "X", y = "Y") +
  scale_color_manual(values = c("Ytest_Xint" = "deepskyblue", "Div_Sample" = "darkgoldenrod1"),
                     name = "",
                     labels = c(Ytest_Xint = "True", Div_Sample = "Estimated")) +
  theme(
    legend.position = c(0.97, 0.03),
    legend.justification = c(1, 0), 
    legend.box.just = "left",
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "mm"),
    legend.background = element_rect(fill = "white", linetype = "solid", colour = "white"),
    legend.text = element_text(size = 18),
    axis.text.x = element_text(size = 14), 
    axis.text.y = element_text(size = 14), 
    axis.title.x = element_text(size = 16), 
    axis.title.y = element_text(size = 16)
  )


###

# Interventional quantile function

data_frame_q25 <- data.frame(Value = c(true_sample_q25, div_sample_q25), 
                             Type = rep(c("True", "DIV"), each = 1000), 
                             Quantile = "Q25")
data_frame_q50 <- data.frame(Value = c(true_sample_q50, div_sample_q50), 
                             Type = rep(c("True", "DIV"), each = 1000),
                             Quantile = "Q50")
data_frame_q75 <- data.frame(Value = c(true_sample_q75, div_sample_q75), 
                             Type = rep(c("True", "DIV"), each = 1000),
                             Quantile = "Q75")

# Combine data frames
plot_data_2 <- rbind(data_frame_q25, data_frame_q50, data_frame_q75)
plot_data_2$Type <- factor(plot_data_2$Type, levels = c("DIV", "True"))
plot_data_2$Quantile <- factor(plot_data_2$Quantile, 
                               levels = c("Q25", "Q50", "Q75"),
                               labels = c("x[Q25]", "x[Q50]", "x[Q75]"))
ggplot(plot_data_2, aes(x = Value, fill = Type)) + 
  geom_density(alpha = 0.7, aes(group = rev(Type))) +   
  facet_grid(rows = vars(Quantile), labeller = label_parsed) + 
  theme_minimal() + 
  labs(x = "Y", y = expression(P[Y]^{do(X==x)})) + 
  xlim(-3, 12) + 
  scale_fill_manual(values = c("True" = "deepskyblue", "DIV" = "darkgoldenrod1"), 
                    name = "", labels = c("DIV" = "Estimated", "True" = "True")) + 
  theme(legend.position = c(0.97, 0.87), 
        legend.justification = c(1, 0), 
        legend.box.just = "left", 
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "mm"), 
        legend.background = element_rect(fill = "white", linetype = "solid", colour = "white"), 
        legend.text = element_text(size = 18), 
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16)
  )


div_quant <- predict.DIV(div_mod, Xtest = data_test$Xtest_grid, type = "quantile", 
                         quantiles = c(0.1, 0.5, 0.9), nsample = 1000)

plot_data_quant <- data.frame(
  Xtest_grid = data_test$Xtest_grid,
  Q10 = div_quant[, 1],
  Q50 = div_quant[, 2],
  Q90 = div_quant[, 3]
)

q10_grid <- numeric(length(data_test$Xtest_grid))
q50_grid <- numeric(length(data_test$Xtest_grid))
q90_grid <- numeric(length(data_test$Xtest_grid))

num_rep <- 5000

# Loop over each element in Xtest
for (i in 1:length(data_test$Xtest_grid)) {
  # Sample new H and eps_Y vectors for each Xtest value
  eps_Yate <- rnorm(num_rep)
  Hate <- rnorm(num_rep)
  
  # Compute the mean for the current Xtest value
  Yvec <- f(X = rep(data_test$Xtest_grid[i], num_rep), H = Hate, eps_Y = eps_Yate)
  true_q10 <- quantile(Yvec, 0.10)
  true_q50 <- quantile(Yvec, 0.5)
  true_q90 <- quantile(Yvec, 0.90)
  
  # Store the mean value in the mean_vector
  q10_grid[i] <- true_q10
  q50_grid[i] <- true_q50
  q90_grid[i] <- true_q90
}

library(ggplot2)
library(dplyr)
library(tidyr)

# Assuming `q10_grid`, `q50_grid`, `q90_grid` contain the true quantiles
# And `div_quant` contains the estimated quantiles
# Create data frames for both true and DIV quantiles

# True quantiles
plot_data_true <- data.frame(
  Xtest_grid = data_test$Xtest_grid,
  Q10 = q10_grid,
  Q50 = q50_grid,
  Q90 = q90_grid,
  Source = "True"
)

# DIV quantiles
plot_data_div <- data.frame(
  Xtest_grid = data_test$Xtest_grid,
  Q10 = div_quant[, 1],
  Q50 = div_quant[, 2],
  Q90 = div_quant[, 3],
  Source = "DIV"
)

# Combine both data frames
plot_data_combined <- bind_rows(plot_data_true, plot_data_div)

# Convert to long format for ggplot2
plot_data_long_combined <- plot_data_combined %>%
  pivot_longer(cols = c(Q10, Q50, Q90), names_to = "Quantile", values_to = "Value") %>%
  mutate(Source = ifelse(Source == "True", "Q_True", "Q_Estimated"))

plot_data_long_combined$Source <- as.factor(plot_data_long_combined$Source)
plot_data_long_combined$Quantile <- as.factor(plot_data_long_combined$Quantile)

plot_data_1 <- pivot_longer(plot_data_1, cols = c(Ytest_Xint, Div_Sample), 
                          names_to = "Type", values_to = "Value")

# Replace Type with "True" for Ytest_Xint and "Estimated" for Div_Sample
plot_data_1$Type <- recode(plot_data_1$Type, Ytest_Xint = "True", Div_Sample = "Estimated")


ggplot() + 
  # Scatter plot for both true and estimated samples, distinguished by color
  geom_point(data = plot_data_1, aes(x = Xint, y = Value, color = Type), alpha = 0.7) +

  # Line plot for quantiles with updated source labels
  geom_line(data = plot_data_long_combined, 
            aes(x = Xtest_grid, y = Value, group = interaction(Quantile, Source), 
                color = Source), size = 1) +
  
  # Customize plot labels and legend
  labs(x = "X", y = "Y") +
  xlim(min(plot_data_1$Xint), 6) + 
  ylim(0, 15) + 
  
  
  # Customize color for scatter plot points and line colors for quantiles
  scale_color_manual(
    values = c("True" = "deepskyblue", "Estimated" = "darkgoldenrod1", 
               "Q_Estimated" = "grey", "Q_True" = "black"),
    breaks = c("True", "Estimated", "Q_True", "Q_Estimated"), # Specify legend order
    labels = c(
      bquote("Sample from" ~ P[Y]^{do(X==x)}),  
      bquote("Sample from" ~ hat(P)[Y]^{do(X==x)}), 
      bquote(q[alpha]^"*" * (x)), 
      bquote(hat(q)[alpha]^"*" * (x))
    )
  ) +
  
  # Manually create the two legends with proper overrides for shape and linetype
  guides(
    color = guide_legend(
      override.aes = list(
        shape = c(16, 16, NA, NA),   # Dots for scatter points, no shape for lines
        linetype = c(0, 0, 1, 1)     # Solid lines for quantiles, no lines for points
      ),
      order = 1  # Order the guide to make sure it shows up first
    )
  ) +
  
  # Custom theme to manage layout
  theme_classic() + 
  theme(
    legend.position = c(0.18, 0.84),           
    legend.box = "vertical",                
    legend.box.just = "left",                 
    legend.text = element_text(size = 16),
    legend.title = element_blank(),     
    axis.text.x = element_text(size = 14), 
    axis.text.y = element_text(size = 14),  
    axis.title.x = element_text(size = 16), 
    axis.title.y = element_text(size = 16)
  )

