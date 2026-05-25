# ============================================================
# install_packages.R
# Installs all required packages for reproducing the analysis
# ============================================================

cat("Starting package installation...\n")

# -------------------------------
# 1. CRAN packages
# -------------------------------

cran_packages <- c(
  "Matrix",
  "dplyr",
  "gtools",
  "parallel",
  "foreach",
  "doParallel",
  "ggplot2",
  "grid",
  "mclust",
  "plyr",
  "purrr",
  "stringr",
  "geometry",
  "combinat",
  "fossil",
  "gstat",
  "seqinr",
  "sets",
  "spatstat",
  "rlang",
  "clues"   # will be overridden if local version is installed
)

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing CRAN package:", pkg, "\n")
    install.packages(pkg)
  } else {
    cat("Already installed:", pkg, "\n")
  }
}

# -------------------------------
# 2. Install local packages (.tar.gz)
# -------------------------------

local_pkg_dir <- "Codes/Utils/External_Packages"

if (dir.exists(local_pkg_dir)) {
  tar_files <- list.files(local_pkg_dir, pattern = "\\.tar\\.gz$", full.names = TRUE)
  
  if (length(tar_files) > 0) {
    cat("Installing local source packages...\n")
    for (pkg_file in tar_files) {
      cat("Installing from:", pkg_file, "\n")
      tryCatch({
        install.packages(pkg_file, repos = NULL, type = "source")
      }, error = function(e) {
        cat("Failed to install:", pkg_file, "\n")
      })
    }
  } else {
    cat("No .tar.gz files found\n")
  }
} else {
  cat("No local package directory found\n")
}

# -------------------------------
# 3. GitHub package: clustRviz
# -------------------------------

# install.packages("devtools")
devtools::install_github("DataSlingers/clustRviz")

# -------------------------------
# 4. Final check
# -------------------------------

all_packages <- c(cran_packages, "cvxclustr","clustRviz")

cat("\nChecking installed packages:\n")

for (pkg in all_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat("✔", pkg, "\n")
  } else {
    cat("✘", pkg, "NOT INSTALLED\n")
  }
}

cat("\nPackage installation complete.\n")
