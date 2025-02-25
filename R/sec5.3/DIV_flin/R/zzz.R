.onLoad <- function(libname, pkgname) {
  if (!requireNamespace("torch", quietly = TRUE)) {
    warning("Torch is not installed. Some functionality may not work.")
  } else if (!torch::torch_is_installed()) {
    torch::install_torch()
  }
}
