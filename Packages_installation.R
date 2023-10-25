# package installation and loading code using the install_and_load_package function.
install_and_load_package <- function(package) {
  if (!(package %in% installed.packages()[,"Package"])) {
    install.packages(package)
  }
  library(package, character.only = TRUE)
}

required_packages <- c(
  "readxl", "gbm", "DiagrammeR", "DiagrammeRsvg", "Rcpp", "rpart",
  "rpart.plot", "randomForest", "corrplot", "ModelMetrics", "ROCR", "caret",
  "mlbench", "kernlab", "e1071", "MASS", "klaR", "dplyr", "naivebayes",
  "DescTools", "pROC", "rsample", "titanic", "funModeling", "GA", "DALEX", "optparse",
  "plotly", "ggplot2", "smotefamily", "ingredients","r2d3","ROSE", "patchwork", "devtools", "ggpubr",
  "catboost", "devEMF"
)

for (package in required_packages) {
  install_and_load_package(package)
