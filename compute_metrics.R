<<<<<<< HEAD
source('F1-score.R')
# Define function to evaluate models (accuracy, precision,...) with confidence intervals
compute_metrics <- function(confusion_matrix) {
  accuracy <- BinomCI(sum(diag(confusion_matrix)), sum(confusion_matrix), conf.level = 0.95, method = "clopper-pearson")
  precision <- BinomCI(confusion_matrix[2, 2], sum(confusion_matrix[, 2]), conf.level = 0.95, method = "clopper-pearson")
  sensitivity <- BinomCI(confusion_matrix[2, 2], sum(confusion_matrix[2, ]), conf.level = 0.95, method = "clopper-pearson")
  specificity <- BinomCI(confusion_matrix[1, 1], sum(confusion_matrix[1, ]), conf.level = 0.95, method = "clopper-pearson")
  f1_score <- f1scores(confusion_matrix)$Confidence.Interval[2,]  
  
  result <- list(
    Accuracy = accuracy,
    Precision = precision,
    Sensitivity = sensitivity,
    Specificity = specificity,
    F1_Score = f1_score
  )
  
  return(result)
=======
source('F1-score.R')
# Define function to evaluate models (accuracy, precision,...) with confidence intervals
compute_metrics <- function(confusion_matrix) {
  accuracy <- BinomCI(sum(diag(confusion_matrix)), sum(confusion_matrix), conf.level = 0.95, method = "clopper-pearson")
  precision <- BinomCI(confusion_matrix[2, 2], sum(confusion_matrix[, 2]), conf.level = 0.95, method = "clopper-pearson")
  sensitivity <- BinomCI(confusion_matrix[2, 2], sum(confusion_matrix[2, ]), conf.level = 0.95, method = "clopper-pearson")
  specificity <- BinomCI(confusion_matrix[1, 1], sum(confusion_matrix[1, ]), conf.level = 0.95, method = "clopper-pearson")
  f1_score <- f1scores(confusion_matrix)$Confidence.Interval[2,]  
  
  result <- list(
    Accuracy = accuracy,
    Precision = precision,
    Sensitivity = sensitivity,
    Specificity = specificity,
    F1_Score = f1_score
  )
  
  return(result)
>>>>>>> c18c5be1391a6c1982e85029db52adc20aff265f
}