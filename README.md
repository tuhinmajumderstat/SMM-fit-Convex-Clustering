# Fitting Sparse Markov Model with Convex Clustering

## Dependencies
Two major packages are installed from the CRAN archiv or other github repos - `cvxclustr` for convex clustering and `clustRviz` for applying Weylandt's pathwise technique.
- Install `cvxclustr` using the following command
  ```{r}
  install.packages("https://cran.r-project.org/src/contrib/Archive/cvxclustr/cvxclustr_1.1.1.tar.gz", repos = NULL, type = "source")
  ```
  If there is an error, you can manually download the `cvxclustr_1.1.1.tar.gz` file locally from this [link](https://cran.r-project.org/src/contrib/Archive/cvxclustr/) and install in R. `Rtools` is needed to run this package.
- Install `clustRviz` using
  ```{r}
  # install.packages("devtools")
  devtools::install_github("DataSlingers/clustRviz")
  ```
- The other packages that need to be installed before running the scripts
 ```{r}
  packages <- c(
  "gtools", 
  "dplyr", 
  "Matrix", 
  "parallel", 
  "ggplot2", 
  "clues", 
  "combinat", 
  "fossil", 
  "geometry", 
  "gstat", 
  "LaplacesDemon", 
  "plyr", 
  "sets", 
  "spatstat", 
  "stringr", 
  "mclust"
)

# Install missing packages
installed_packages <- rownames(installed.packages())
for (pkg in packages) {
  if (!(pkg %in% installed_packages)) {
    install.packages(pkg)
  }
  require(pkg, character.only = TRUE)
}

  ```

## File Description
- The codes are saved in `/SMM-fit-Convex-Clustering/Codes`
  - The simulation 1 codes are in `/SMM-fit-Convex-Clustering/Codes/Simulation 1`, including the convex clustering approach, other smm fitting approaches and the codes for generating the plots. To run `SMM_Simulation_1_Other_Methods.R`, please import the functions from `SMM_Fit_Other_Methods_Functions.R`.
  - Similar notations go for Simulation 2 and Real Data Analysis.   
- The raw datasets are saved in `/SMM-fit-Convex-Clustering/Data`.
- Results for Simulation 1, Simulation 2 and Real Data analysis are stored in separate sub directories in `/SMM-fit-Convex-Clustering/Results`.
