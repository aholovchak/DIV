library(httr)
library(readstata13)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ivreg)

# Data description: https://ivmodels.readthedocs.io/en/latest/examples/acemoglu2001colonial.html

url4 <- "https://www.dropbox.com/scl/fi/3yuv9j514zuajzjfluoc1/maketable4.zip?rlkey=pq9l7bxktw1iqxe6fmoh26g79&e=1&dl=1"
temp_file <- tempfile(fileext = ".zip")
GET(url4, write_disk(temp_file, overwrite = TRUE))

# Extract the DTA file from the ZIP archive
temp_dir <- tempdir()
unzip(temp_file, exdir = temp_dir)
dta_file <- list.files(temp_dir, pattern = "\\.dta$", full.names = TRUE)

df4 <- read.dta13(dta_file)
df4 <- df4 %>% filter(baseco == 1)

df4 <- df4 %>%
  mutate(other_continent = as.integer(shortnam %in% c("AUS", "MLT", "NZL")))

head(df4)

# Specify the endogenous and instrumental variables
X <- df4$avexpr  # Average protection against expropriation risk
Z <- df4$logem4 # Log European settler mortality
Y <- df4$logpgp95

setwd("DistributionIV")
library(devtools)
load_all()

div_mod <- div(X = X, Z = Z, Y = Y,
                     epsx_dim = 50, epsh_dim = 50, epsy_dim = 50, lr = 1e-4,
                     num_epochs = 20000, num_layer = 4)

ols_mod <-lm(Y ~ X)
tsls_mod <- ivreg(Y ~ X | Z)
 
avexpr_seq <- seq(min(X), max(X), length.out = 100)


# DIV
mean_pred_div <- predict.DIV(div_mod,
                               Xtest = avexpr_seq, 
                               nsample = 1000)

data <- data.frame(
  avexpr_seq = avexpr_seq,
  mean_pred_div = mean_pred_div,
  mean_pred_2sls = tsls_mod$coefficients[1] + tsls_mod$coefficients[1] * avexpr_seq,
  mean_pred_ols = ols_mod$coefficients[1] + ols_mod$coefficients[2] * avexpr_seq
)


# Combine data into a long-format data frame for ggplot
data_long <- data %>%
  pivot_longer(
    cols = c(mean_pred_div,
             mean_pred_2sls,
             mean_pred_ols),
    names_to = "method",
    values_to = "value"
  )

# Rename methods for better legend labels
data_long$method <- recode(data_long$method,
                           "mean_pred_div" = "DIV",
                           "mean_pred_2sls" = "2SLS",
                           "mean_pred_ols" = "OLS")

ggplot() +
  geom_point(data = df4, aes(x = avexpr, y = logpgp95), 
             color = "grey", alpha = 0.5) + 
  geom_line(data = data_long, aes(x = avexpr_seq, y = value, color = method), 
            size = 1) + 
  scale_color_brewer(
    palette = "Set2",
    name = "method",
    limits = c("DIV", "2SLS", "OLS")
  ) +
  # Add labels and title
  labs(
    title = "",
    x = "Average protection against expropriation risk, 1985-95",
    y = "Log GDP per capita, 1995"
  ) +
  # Style adjustments
  theme_minimal() +
  theme(legend.position = "top",
        legend.justification = "center",
        legend.background = element_rect(color = "black", fill = NA))