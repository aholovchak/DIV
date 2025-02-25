setwd("results/sec5.3")
library(readr)

# DIV

ui_DIV_flin_n1e3 <- read_csv("ui_DIV_flin_n1e3.csv")
ui_DIV_flin_n1e4 <- read_csv("ui_DIV_flin_n1e4.csv")

(mean_DIV_flin_n1e3 <- round(colMeans(ui_DIV_flin_n1e3, na.rm = TRUE), 3))
(mean_DIV_flin_n1e4 <- round(colMeans(ui_DIV_flin_n1e4, na.rm = TRUE), 3))

### 

# HSIC-X

ui_hsicx_flin_n1e3 <- read_csv("ui_hsicx_flin_n1000.csv")
ui_hsicx_flin_n1e4 <- read_csv("ui_hsicx_flin_n10000.csv")

beta_norm_hsicx_1000 <- apply(ui_hsicx_flin_n1e3 - c(1, 2), 2, function(x) sqrt(sum(x^2)))
beta_norm_hsicx_10000 <- apply(ui_hsicx_flin_n1e4 - c(1, 2), 2, function(x) sqrt(sum(x^2)))

mean(beta_norm_hsicx_1000)
mean(beta_norm_hsicx_10000)
