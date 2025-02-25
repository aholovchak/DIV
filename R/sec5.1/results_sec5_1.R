library(tidyverse)
library(ggplot2)
library(patchwork)
library(stringr)


# load all data

# Test data
data_test_dgp1 <- read_csv("data/sec5.1/test_Zbin_g_lin_f_lin.csv")
data_test_dgp2 <- read_csv("data/sec5.1/test_Zcont_g_lin_f_lin.csv")
data_test_dgp3 <- read_csv("data/sec5.1/test_Zbin_g_lin_f_log_case.csv")
data_test_dgp4 <- read_csv("data/sec5.1/test_Zcont_g_lin_f_log_case.csv")
data_test_dgp5 <- read_csv("data/sec5.1/test_Zbin_g_lin_f_sin_lin.csv")
data_test_dgp6 <- read_csv("data/sec5.1/test_Zcont_g_lin_f_sin_lin.csv")

# Train data
data_train_dgp1 <- read_csv("data/sec5.1/train_Zbin_g_lin_f_lin.csv")
data_train_dgp2 <- read_csv("data/sec5.1/train_Zcont_g_lin_f_lin.csv")
data_train_dgp3 <- read_csv("data/sec5.1/train_Zbin_g_lin_f_log_case.csv")
data_train_dgp4 <- read_csv("data/sec5.1/train_Zcont_g_lin_f_log_case.csv")
data_train_dgp5 <- read_csv("data/sec5.1/train_Zbin_g_lin_f_sin_lin.csv")
data_train_dgp6 <- read_csv("data/sec5.1/train_Zcont_g_lin_f_sin_lin.csv")


base_path <- "results/sec5.1/"

methods <- c("div", "cflin", "cfnonlin", "deepiv", "deepgmm", "hsic")
result_types <- c("mse", "plot")
results_list <- list()

for (method in methods) {
  for (result_type in result_types) {
    file_paths <- paste0(base_path, method, "_result_", result_type, "_dgp", 1:6, ".csv")
    
    # Read files into a list
    results_list[[paste0(method, "_", result_type)]] <- map(file_paths, read_csv)
    names(results_list[[paste0(method, "_", result_type)]]) <- paste0("dgp", 1:6)
  }
}


mse_plot_fct <- function(dgp, data_test, data_train) {
  # List of data frames to be renamed
  data_frames_mse <- list(results_list[["div_mse"]][[dgp]][1:10], 
                          results_list[["cflin_mse"]][[dgp]],
                          results_list[["cfnonlin_mse"]][[dgp]],
                          results_list[["deepiv_mse"]][[dgp]],
                          results_list[["deepgmm_mse"]][[dgp]],
                          results_list[["hsic_mse"]][[dgp]])
  
  div_colnames_mse <- colnames(results_list[["div_mse"]][[dgp]])[1:10]
  
  # Loop through each data frame in MSE and update column names
  for (i in seq_along(data_frames_mse)) {
    colnames(data_frames_mse[[i]]) <- div_colnames_mse
  }
  
  
  # Bind all MSE data frames
  data_all_mse <- do.call(rbind, data_frames_mse)
  data_all_mse$X_test <- rep(data_test$Xtest_grid, 6)
  data_all_mse$method <- rep(c("DIV", "CF linear", "CF nonlinear", "DeepIV", "DeepGMM", "HSIC-X"), each = length(as.data.frame(data_test)[,1]))
  
  # Pivot and mutate MSE data
  long_data_all_mse <- data_all_mse %>%
    pivot_longer(
      cols = starts_with("V"),
      names_to = "num_run", 
      values_to = "Y_hat"
    ) %>%
    mutate(num_run = as.integer(gsub("V", "", num_run)))
  
  # Calculate mean MSE
  mean_mse_data <- long_data_all_mse %>%
    group_by(method, num_run) %>%
    summarize(MSE = mean((Y_hat - data_test$mean_int)^2), .groups = 'drop') %>%
    group_by(method) %>%
    summarize(Mean_MSE = mean(MSE))
  
  mean_maxse_data <- long_data_all_mse %>%
    group_by(method, num_run) %>%
    summarize(maxSE = max((Y_hat - data_test$mean_int)^2), .groups = 'drop') %>%
    group_by(method) %>%
    summarize(Mean_maxSE = mean(maxSE))
  
  # List of data frames for plotting
  data_frames_plot <- list(results_list[["div_plot"]][[dgp]][1:10], 
                          results_list[["cflin_plot"]][[dgp]],
                          results_list[["cfnonlin_plot"]][[dgp]],
                          results_list[["deepiv_plot"]][[dgp]],
                          results_list[["deepgmm_plot"]][[dgp]],
                          results_list[["hsic_plot"]][[dgp]])
  
  div_colnames_plot <- colnames(results_list[["div_plot"]][[dgp]])[1:10]
  
  # Loop through each data frame in Plot and update column names
  for (i in seq_along(data_frames_plot)) {
    colnames(data_frames_plot[[i]]) <- div_colnames_plot
  }
  
  # Bind all plot data frames
  data_all_plot <- do.call(rbind, data_frames_plot)
  data_all_plot$X_test <- rep(data_test$Xtest_grid, 6)
  data_all_plot$method <- rep(c("DIV", "CF linear", "CF nonlinear", 
                                "DeepIV", "DeepGMM", "HSIC-X"), 
                              each = length(as.data.frame(data_test)[,1]))
  
  # Pivot and mutate plot data
  long_data_all_plot <- data_all_plot %>%
    pivot_longer(
      cols = starts_with("V"), # Selecting columns V1 to V10
      names_to = "num_run", 
      values_to = "Y_hat"
    ) %>%
    mutate(num_run = as.integer(gsub("V", "", num_run)))
  
  # Add method and run information to data_test and data_train data frames
  data_test$method <- rep("data_test", length(as.data.frame(data_test)[,1]))
  data_test$num_run <- rep("0", length(as.data.frame(data_test)[,1]))
  data_train$method <- rep("", length(as.data.frame(data_train)[,1]))
  data_train$num_run <- rep("", length(as.data.frame(data_train)[,1]))
  
  # Generate summary data
  summary_data <- long_data_all_plot %>%
    group_by(X_test, method) %>%
    summarize(
      Y_median = median(Y_hat),
      Y_min = min(Y_hat),
      Y_max = max(Y_hat),
      .groups = 'drop'
    )
  
  # Create base plot
  base_plot <- ggplot() +
    geom_point(data = data_train, aes(x = X, y = Y), color = "grey", alpha = 0.1) +
    labs(x = "X", y = "Y", color = "method") +
    theme_minimal() +
    scale_fill_brewer(palette = "Set2")
  
  summary_data$method <- factor(summary_data$method, levels = c("DIV",
                                                                "HSIC-X", "CF linear", "CF nonlinear",
                                                                "DeepGMM", "DeepIV", "Engression"))
  
  # Create final plot
  final_plot <- base_plot + 
    geom_line(data = summary_data, aes(x = X_test, y = Y_median, color = method), linewidth = 0.8,) +
    geom_ribbon(data = summary_data, aes(x = X_test, ymin = Y_min, ymax = Y_max, fill = method), alpha = 0.1) +
    geom_line(data = data_test, aes(x = Xtest_grid, y = mean_grid), linetype = "dotted", linewidth = 0.65, color = "black") +
    scale_color_brewer(palette = "Set2") +
    scale_fill_brewer(palette = "Set2") +
    theme(legend.position = "bottom")
  
  return(list(final_plot, mean_mse_data, mean_maxse_data))
}

# Create plots for each DGP
plot_list <- lapply(1:6, function(dgp) {
  mse_plot_fct(dgp = paste0("dgp", dgp), data_test = get(paste0("data_test_dgp", dgp)), data_train = get(paste0("data_train_dgp", dgp)))[[1]]
})

# Arrange the plots
combined_plot <- (plot_list[[1]] | plot_list[[3]] | plot_list[[5]]) /
  (plot_list[[2]] | plot_list[[4]] | plot_list[[6]]) +
  plot_layout(guides = 'collect') &
  theme(
    legend.position = "top",
    legend.background = element_rect(color = "black", fill = NA)  # Adds a rectangle around the legend
  )

combined_plot


######
# MSE
######


# Function to extract the legend from a ggplot object
g_legend <- function(a.gplot) {
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# Assuming `long_df_maxse` is your data frame
# Split the data into subsets for different y-axes
df_1_2 <- long_df_mse %>% filter(DGP %in% c("1", "2"))
df_3_4 <- long_df_mse %>% filter(DGP %in% c("3", "4"))
df_5_6 <- long_df_mse %>% filter(DGP %in% c("5", "6"))

p1 <- ggplot(df_1_2, aes(x = DGP, y = Value, color = method, shape = method)) +
  geom_point(size = 4, stroke = 1.2, position = position_dodge(width = 0.4)) +
  scale_color_brewer(palette = "Set2") +
  scale_shape_manual(values = c(16, 17, 18, 19, 15, 8, 7)) +
  theme_minimal() +
  labs(y = "mean MSE", x = NULL) +  # Remove x-axis label
  theme(axis.text.x = element_text(angle = 0, hjust = 1), legend.position = "none") + 
  ylim(0.0, 0.4)

p2 <- ggplot(df_3_4, aes(x = DGP, y = Value, color = method, shape = method)) +
  geom_point(size = 4, stroke = 1.2, position = position_dodge(width = 0.4)) +
  scale_color_brewer(palette = "Set2") +
  scale_shape_manual(values = c(16, 17, 18, 19, 15, 8, 7)) +
  theme_minimal() +
  labs(y = NULL, x = NULL) +  # Remove x-axis label
  theme(axis.text.x = element_text(angle = 0, hjust = 1), legend.position = "none") + 
  ylim(0.0, 0.4)

p3 <- ggplot(df_5_6, aes(x = DGP, y = Value, color = method, shape = method)) +
  geom_point(size = 4, stroke = 1.2, position = position_dodge(width = 0.4)) +
  scale_color_brewer(palette = "Set2") +
  scale_shape_manual(values = c(16, 17, 18, 19, 15, 8, 7)) +
  theme_minimal() +
  labs(y = NULL, x = NULL) +  # Remove x-axis label
  theme(axis.text.x = element_text(angle = 0, hjust = 1), legend.position = "none")

# Extract the legend from one of the plots
legend <- g_legend(
  ggplot(df_5_6, aes(x = DGP, y = log(Value + 1), color = method, shape = method)) +
    geom_point(size = 4, stroke = 1.2, position = position_dodge(width = 0.5)) +
    scale_color_brewer(palette = "Set2") +
    scale_shape_manual(values = c(16, 17, 18, 19, 15, 8, 7)) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      legend.justification = "center",
      legend.background = element_rect(color = "black", fill = NA)
    )
)

# Create a common x-axis label grob
x_axis_label <- textGrob("Scenario", gp = gpar(fontsize = 12))

# Arrange plots side by side and add the common x-axis label and shared legend
combined_plot <- grid.arrange(
  legend,
  arrangeGrob(p1, p2, p3, ncol = 3),
  x_axis_label,
  nrow = 3,
  heights = c(1, 6, 1)  # Adjust height proportions to fit the label
)

# Display the combined plot
combined_plot

