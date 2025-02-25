library(ggplot2)
library(dplyr)
library(gridExtra)
library(patchwork)
library(readr)
library(tidyr)

create_rmse_plot <- function(data, desired_taus, show_y_axis = TRUE) {
  # Convert method to factor and update method names
  data$method <- as.factor(data$method)
  
  method_labels <- c("DIV" = "DIV", 
                     "DIVE" = "DIVE", 
                     "IVQR" = "IVQR", 
                     "ORACLE" = "ORACLE")
  
  data$method <- recode(data$method, !!!method_labels)
  
  filtered_data <- data %>%
    group_by(method) %>%  # Group by method
    mutate(rounded_tau = sapply(ntau, function(x) {
      desired_taus[which.min(abs(desired_taus - x))]
    })) %>%
    group_by(method, rounded_tau) %>%  # Group by method and rounded_tau
    slice_min(abs(ntau - rounded_tau), n = 1) %>%  # Keep the closest `ntau` to `rounded_tau`
    ungroup()  # Ungroup to remove grouping after processing
  
  desired_order <- c("DIV", "DIVE", "IVQR", "ORACLE")
  
  # Change the order of the method factor levels
  filtered_data$method <- factor(filtered_data$method, levels = desired_order)
  
  # Calculate RMSE for each method relative to ORACLE for each quantile and iteration
  rmse_data <- filtered_data %>%
    left_join(filtered_data %>% filter(method == "ORACLE") %>% 
                select(rounded_tau, iter, oracle_qte = qte), 
              by = c("rounded_tau", "iter")) %>%
    filter(method != "ORACLE") %>%
    # Compute RMSE with oracle QTE values
    mutate(rmse = sqrt((qte - oracle_qte)^2)) %>%
    select(rounded_tau, method, iter, rmse)
  
  # Create a boxplot for RMSE per method
  plot <- ggplot(rmse_data, aes(x = as.factor(rounded_tau), y = rmse, fill = method)) +
    geom_boxplot() +
    coord_cartesian(ylim = c(0, 10)) +
    labs(x = expression(alpha), y = if (show_y_axis) "RMSE" else NULL, fill = "method")
  theme_minimal() +
    theme(
      axis.title.y = element_text(size = 12),
      axis.text.y = if (show_y_axis) element_text() else element_blank(),
      axis.ticks.y = if (show_y_axis) element_line() else element_blank()
    )
  
  return(plot)
}


setwd("results/sec5.2/")

# Load data for scenario 1 and 2
pd_combined1 <- read.csv("qte_scenario1.csv")
pd_combined2 <- read.csv("qte_scenario2.csv")

# Define the desired taus
desired_taus <- c(0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 0.99)

# Function to extract the legend from a ggplot object
g_legend <- function(a.gplot) {
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# Extract the legend from one of the plots (e.g., plot1)
mylegend <- g_legend(
  create_rmse_plot(pd_combined1, desired_taus) +
    theme_minimal() +
    theme(legend.position = "bottom", legend.justification = "center")
)

plot1 <- create_rmse_plot(pd_combined1, desired_taus, show_y_axis = TRUE) +
  coord_cartesian(ylim = c(0, 15)) +  
  labs(y = "RMSE") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position = "none")

plot2 <- create_rmse_plot(pd_combined2, desired_taus, show_y_axis = FALSE) +
  coord_cartesian(ylim = c(0, 7)) +
  labs(y = "RMSE") +
  theme_minimal() +  
  theme(legend.position = "none") + 
  scale_fill_brewer(palette = "Set2")

mylegend <- g_legend(
  create_rmse_plot(pd_combined2, desired_taus) +
    theme_minimal() + 
    theme(
      legend.position = "bottom", 
      legend.justification = "center",
      legend.background = element_rect(color = "black", fill = NA),
      legend.box.background = element_rect(color = "black")
    ) +
    scale_fill_brewer(palette = "Set2")
  
)

combined_plot <- grid.arrange(
  arrangeGrob(plot1, plot2, ncol = 1),  
  mylegend,       
  nrow = 2,
  heights = c(10, 1)
)