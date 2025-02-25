library(ggplot2)
library(MASS)
library(splines)
library(tidyr)
library(dplyr)

engr <- read.csv("results/sec6.2/engression_singlecell_pairwise.csv")
div <- read.csv("results/sec6.2/div_singlecell_pairwise.csv")
hsicx <- read.csv("results/sec6.2/hsic_singlecell_pairwise.csv")
deepgmm <- read.csv("results/sec6.2/deepgmm_singlecell_pairwise.csv")
# deepiv <- read.csv("results/sec6.2/deepiv_singlecell_pairwise.csv")
cflin <- read.csv("results/sec6.2/cflin_singlecell_pairwise.csv")
cfnonlin <- read.csv("results/sec6.2/cfnonlin_singlecell_pairwise.csv")


num_vars <- ncol(div)

##############
## Engression
ssd_results_engr <- apply(engr, 1, function(row) {
  pairwise_diffs <- combn(row, 2, FUN = function(x) (x[1] - x[2])^2)  # Compute squared differences for all pairs
  sum(pairwise_diffs)
})
engr$SSD <- ssd_results_engr
mean(engr$SSD)


##############
## DIV
ssd_results_div <- apply(div, 1, function(row) {
  pairwise_diffs <- combn(row, 2, FUN = function(x) (x[1] - x[2])^2)
  sum(pairwise_diffs)
})
div$SSD <- ssd_results_div
mean(div$SSD)

##############
## HSIC-X
ssd_results_hsic <- apply(hsicx, 1, function(row) {
  pairwise_diffs <- combn(row, 2, FUN = function(x) (x[1] - x[2])^2)
  sum(pairwise_diffs)
})
hsicx$SSD <- ssd_results_hsic
mean(hsicx$SSD)

##############
## DeepGMM
ssd_results_deepgmm <- apply(deepgmm, 1, function(row) {
  pairwise_diffs <- combn(row, 2, FUN = function(x) (x[1] - x[2])^2)
  sum(pairwise_diffs)
})
deepgmm$SSD <- ssd_results_deepgmm
mean(deepgmm$SSD)

##############
# ## DeepIV
# 
# # remove NAs columns
# deepiv <- deepiv[, colSums(!is.na(deepiv)) > 0]
# 
# 
# ssd_results_deepiv <- apply(deepiv, 1, function(row) {
#   pairwise_diffs <- combn(row, 2, FUN = function(x) (x[1] - x[2])^2)
#   sum(pairwise_diffs)
# })
# deepiv$SSD <- ssd_results_deepiv
# mean(deepiv$SSD)

##############
## CF linear
ssd_results_cflin <- apply(cflin, 1, function(row) {
  pairwise_diffs <- combn(row, 2, FUN = function(x) (x[1] - x[2])^2)
  sum(pairwise_diffs)
})
cflin$SSD <- ssd_results_cflin
mean(cflin$SSD)

##############
## CF nonlinear
ssd_results_cfnonlin <- apply(cfnonlin, 1, function(row) {
  pairwise_diffs <- combn(row, 2, FUN = function(x) (x[1] - x[2])^2)
  sum(pairwise_diffs)
})
cfnonlin$SSD <- ssd_results_cfnonlin
mean(cfnonlin$SSD)


##############
# RESULTS

# Compute mean SSD for each method
ssd_summary <- data.frame(
  method = c("DIV", "Engression", "CF linear", "CF nonlinear", "DeepGMM", "HSIC-X"),
  mean_SSD = c(mean(div$SSD, na.rm = TRUE),
               mean(engr$SSD, na.rm = TRUE),
               mean(cflin$SSD, na.rm = TRUE),
               mean(cfnonlin$SSD, na.rm = TRUE),
               mean(deepgmm$SSD, na.rm = TRUE),
               # mean(deepiv$SSD, na.rm = TRUE),
               mean(hsicx$SSD, na.rm = TRUE))
)

# Print the results
print(ssd_summary, row.names = FALSE)

