#' Determine Device for Torch Computations
#'
#' This function selects the most appropriate device (e.g., CUDA, MPS, or CPU)
#' for Torch computations based on system availability.
#'
#' @keywords internal
#' @import torch
#' @noRd
use_device <- function() {
  # Check if torch is installed and available
  if (!requireNamespace("torch", quietly = TRUE)) {
    warning("Torch is not installed. Defaulting to CPU.")
    return("cpu")
  }

  # Ensure torch dependencies are installed
  if (!torch::torch_is_installed()) {
    warning("Torch dependencies are missing. Installing now...")
    torch::install_torch()
  }

  # Check for available device
  if (torch::cuda_is_available()) {
    return(eval(parse(text = "torch::cuda_device()")))
  } else if (torch::backends_mps_is_available()) {
    return(torch::torch_device("mps"))
  } else {
    return(torch::torch_device("cpu"))
  }
}
