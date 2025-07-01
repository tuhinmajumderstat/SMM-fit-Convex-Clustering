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
