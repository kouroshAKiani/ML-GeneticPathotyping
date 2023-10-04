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
}

source('F1-score.R')
source('compute_metrics.R')

parser <- OptionParser(option_list = list(
  make_option(c("-f", "--file_path"), type="character", help="The file path of the data file."),
  make_option(c("-t", "--task"), type="character", help="Specifying task"),
  make_option(c("-o", "--model"), type="character", help="Specifying an ML model to train and evaluate"),
  make_option(c("-m", "--method"), type="character", help="Specifying feature selection method"),
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE, help="Print detailed information")
))

args = parse_args(parser)

file_path = args$file_path
task = args$task
model = args$model
method = args$method
verbose = args$verbose

# Read data from the file
Bio_Data <- read_xlsx(file_path, sheet = "Sheet1")

# Define conversion mappings
## Genes      :      1 = Gene Presence    0 = Gene Absence
## Antibiotics :      1 =  isolate's Resistance    0 = isolate's Susceptiblity

# Convert categorical variables to factors
categorical_columns <- colnames(Bio_Data)
Bio_Data[, categorical_columns[-1]] <- lapply(Bio_Data[, categorical_columns[-1]], factor)

# Add a new level for Ampicillin variable
levels(Bio_Data$Ampicillin) <- c(levels(Bio_Data$Ampicillin), "0")

# Display summary and dimensions of the data
cat("Summary:\n")
summary(Bio_Data)
cat("Dimensions:", dim(Bio_Data), "\n")

# Calculate the proportion of OXA48_gene
proportion_OXA48_gene <- table(Bio_Data$OXA48_gene) %>% prop.table()
cat("Proportion of OXA48 gene:\n")
print(proportion_OXA48_gene)


# Divide the dataset into train and test sets
set.seed(123)
split_strat  <- initial_split(Bio_Data, prop = 0.8, strata = "OXA48_gene")
train_strat_OXA48  <- training(split_strat)
test_strat_OXA48 <- testing(split_strat)

# Display the dimensions of the train and test sets
cat("train set dimensions:", dim(train_strat_OXA48), "\n")
cat("test set dimensions:", dim(test_strat_OXA48), "\n")


if (!is.null(task)) {
  if (task == "model_evaluation") {
    # Task: Model Evaluation
    # You can add flags for specific model evaluation methods to use, e.g., --model_evaluation=confusion_matrix or --model_evaluation=roc_auc.
    # Perform the selected model evaluation method based on the provided flags.
    # Example: If user passes --model_evaluation=confusion_matrix, evaluate the model using confusion matrix.
    # If user passes --model_evaluation=roc_auc, evaluate the model using ROC AUC.
    # ...
    if (model == "LR-Chi Squared test"){
      # Logistic regression (Model 1)--------------------------------------------------
      model_Log1_OXA48 <- glm(OXA48_gene ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                                Cefepime + Ceftazidime + Ciprofloxacin +
                                Meropenem + Ceftriaxone + Aztreonam 
                              , family = "binomial", data = train_strat_OXA48)
      
      cat("Summary of LR model based on Chi Squared test including Null deviance, Residual deviance, and AIC, etc.:\n")
      print(summary(model_Log1_OXA48))
      
      test_strat_OXA48$probs1_OXA48 <- predict(model_Log1_OXA48, test_strat_OXA48, type = "response")
      
      test_strat_OXA48$pred_logreg1_OXA48 <- ifelse(test_strat_OXA48$probs1_OXA48 >= 0.5, 1, 0)
      
      # Model 1 evaluation
      confm_logreg1_OXA48 <- table(actual = test_strat_OXA48$OXA48_gene, prediction = test_strat_OXA48$pred_logreg1_OXA48)
      
      metrics_logreg1_OXA48 <- compute_metrics(confm_logreg1_OXA48)
      
      cat("Performance of LR model based on Chi-Squared test results:\n")
      metrics_logreg1_OXA48
    }
    else if (model == "LR-model agnostic/Wald test"){
      model_Log2_OXA48 <- glm(OXA48_gene ~ Meropenem + Imipenem
                              , family = "binomial", data = train_strat_OXA48)
      
      cat("Summary of LR model based on Wald test including Null deviance, Residual deviance, and AIC, etc.:\n")
      print(summary(model_Log2_OXA48))
      
      ##Prediction on test_strat_OXA48 (Model 2)--------------------------------------------------
      test_strat_OXA48$probs2_OXA48 <- predict(model_Log2_OXA48, test_strat_OXA48, type = "response")
      
      test_strat_OXA48$pred_logreg2_OXA48 <- ifelse(test_strat_OXA48$probs2_OXA48 >= 0.5, 1, 0)
      
      # Model 2 evaluation
      confm_logreg2_OXA48 <- table(actual = test_strat_OXA48$OXA48_gene, prediction = test_strat_OXA48$pred_logreg2_OXA48)
      
      metrics_logreg2_OXA48 <- compute_metrics(confm_logreg2_OXA48)
      
      cat("Performance of LR model based on Wald test and model-agnostic approach:\n")
      metrics_logreg2_OXA48
    }
    
    else if (model == "LR-Simplicity Principle_MEM"){
      model_Log3_OXA48 <- glm(OXA48_gene ~ Meropenem
                              , family = "binomial", data = train_strat_OXA48)
      
      cat("Summary of LR model based on simplicity principle (MEM) including Null deviance, Residual deviance, and AIC, etc.:\n")
      print(summary(model_Log3_OXA48))
      
      ##Prediction on test_strat_OXA48
      test_strat_OXA48$probs3_OXA48 <- predict(model_Log3_OXA48, test_strat_OXA48, type = "response")
      
      test_strat_OXA48$pred_logreg3_OXA48 <- ifelse(test_strat_OXA48$probs3_OXA48 >= 0.5, 1, 0)
      
      confm_logreg3_OXA48 <- table(actual = test_strat_OXA48$OXA48_gene, prediction = test_strat_OXA48$pred_logreg3_OXA48)
      
      metrics_logreg3_OXA48 <- compute_metrics(confm_logreg3_OXA48)
      
      cat("Performance of LR model based on simplicity principle (MEM):\n")
      metrics_logreg3_OXA48
    }
    
    else if (model == "NBC-Chi Squared test"){
      train_control <- trainControl(
        method = "cv",
        number = 10
      )
      
      search_grid_OXA48 <- expand.grid(
        usekernel = c(FALSE,TRUE),
        fL = seq(0,5,0.5),
        adjust = seq(0,5,0.5)
      )
      
      model_nb_OXA48 <- train(OXA48_gene ~  Amikacin + Cefotaxime + Gentamicin + Imipenem +
                                Cefepime + Ceftazidime + Ciprofloxacin +
                                Meropenem + Ceftriaxone + Aztreonam,
                              data = train_strat_OXA48,
                              method = "nb",
                              metric = "Accuracy",
                              trControl = train_control,
                              tuneGrid = search_grid_OXA48)
      
      
      cat("Tuned hyper-parameters through 10-fold Cross-Validation:\n")
      print(model_nb_OXA48$bestTune)
      
      ##Prediction on test_strat_OXA48 (Naive Bayes Model)----------------------------------
      test_strat_OXA48$pred_nb_OXA48 <- predict(model_nb_OXA48, test_strat_OXA48)
      
      # Model evaluation
      confm_nb_OXA_48 <- table(actual = test_strat_OXA48$OXA48_gene, prediction = test_strat_OXA48$pred_nb_OXA48)
      metrics_nb_OXA48<- compute_metrics(confm_nb_OXA_48)
      cat("Performance of NBC model based on Chi-Squared test results:\n")
      metrics_nb_OXA48
    }
    else if (model == "NBC-model agnostic/Wald test"){
      train_control <- trainControl(
        method = "cv",
        number = 10
      )
      
      search_grid_OXA48 <- expand.grid(
        usekernel = c(FALSE,TRUE),
        fL = seq(0,5,0.5),
        adjust = seq(0,5,0.5)
      )
      
      model_nb1_OXA48 <- train(OXA48_gene ~  Meropenem + Imipenem,
                               data = train_strat_OXA48,
                               method = "nb",
                               metric = "Accuracy",
                               trControl = train_control,
                               tuneGrid = search_grid_OXA48)
      
      cat("Tuned hyper-parameters through 10-fold Cross-Validation:\n")
      print(model_nb1_OXA48$bestTune)
      
      ##Prediction on test_strat_OXA48 (Naive Bayes Model)----------------------------------
      test_strat_OXA48$pred_nb1_OXA48 <- predict(model_nb1_OXA48, test_strat_OXA48)
      
      # Model evaluation
      confm_nb1_OXA_48 <- table(actual = test_strat_OXA48$OXA48_gene, prediction = test_strat_OXA48$pred_nb1_OXA48)
      
      metrics_nb1_OXA48<- compute_metrics(confm_nb1_OXA_48)
      cat("Performance of NBC model based on model-agnostic approach and Wald test results:\n")
      metrics_nb1_OXA48
    }
    
    else if (model == "NBC-Simplicity Principle_MEM"){
      
      train_control <- trainControl(
        method = "cv",
        number = 10
      )
      
      search_grid_OXA48 <- expand.grid(
        usekernel = c(FALSE,TRUE),
        fL = seq(0,5,0.5),
        adjust = seq(0,5,0.5)
      )
      
      model_nb2_OXA48 <- train(OXA48_gene ~ Meropenem,
                               data = train_strat_OXA48,
                               method = "nb",
                               metric = "Accuracy",
                               trControl = train_control,
                               tuneGrid = search_grid_OXA48)
      
      cat("Tuned hyper-parameters through 10-fold Cross-Validation:\n")
      print(model_nb2_OXA48$bestTune)
      
      ##Prediction on test_strat_OXA48 (Naive Bayes Model)
      test_strat_OXA48$pred_nb2_OXA48 <- predict(model_nb2_OXA48, test_strat_OXA48)
      
      ##confusion matrix(Naive Bayes Model)
      confm_nb2_OXA_48 <- table(actual = test_strat_OXA48$OXA48_gene, prediction = test_strat_OXA48$pred_nb2_OXA48)
      
      metrics_nb2_OXA48<- compute_metrics(confm_nb2_OXA_48)
      
      cat("Performance of NBC model based on Simplicity Principle (MEM):\n")
      metrics_nb2_OXA48
    }
    
    else if (model == "LDA-Chi Squared test"){
      
      model_ldadalex_OXA48 <- train(OXA48_gene ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                                      Cefepime + Ceftazidime + Ciprofloxacin +Aztreonam+
                                      Meropenem + Ceftriaxone  ,
                                    data = train_strat_OXA48,
                                    method = "lda",
                                    metric = "Accuracy")
      
      cat("LDA model based on Chi-Squared test results:\n")
      print(model_ldadalex_OXA48)
      
      #Prediction on test_strat_OXA48 (LDA Model)-----------------------------------------------
      test_strat_OXA48$pred_lda_OXA48 <- predict(model_ldadalex_OXA48, test_strat_OXA48)
      
      # Model Evaluation-------------------------------------------------
      confm_lda_OXA_48 = table(actual = test_strat_OXA48$OXA48_gene, prediction = test_strat_OXA48$pred_lda_OXA48)
      
      metrics_lda_OXA48<- compute_metrics(confm_lda_OXA_48)
      
      cat("Performance of LDA model based on Chi-Squared test results:\n")
      metrics_lda_OXA48
    }
    else if (model == "LDA-model agnostic"){
      model_lda_OXA48 <- train(OXA48_gene ~ Meropenem,
                               data = train_strat_OXA48,
                               method = "lda",
                               metric = "Accuracy")
      
      cat("LDA model based on model-agnostic approach:\n")
      print(model_lda_OXA48)
      
      #Prediction on test_strat_OXA48 (LDA Model)-----------------------------------------------
      test_strat_OXA48$pred_lda1_OXA48 <- predict(model_lda_OXA48, test_strat_OXA48)
      
      # Model evaluation
      confm_lda1_OXA_48 = table(actual = test_strat_OXA48$OXA48_gene, prediction = test_strat_OXA48$pred_lda1_OXA48)
      metrics_lda1_OXA48 <- compute_metrics(confm_lda1_OXA_48)
      cat("Performance of LDA model based on model-agnostic approach:\n")
      metrics_lda1_OXA48
    }
    
    else if (model == "LDA-Wald test"){
      model_lda2_OXA48 <- train(OXA48_gene ~ Imipenem + Meropenem,
                               data = train_strat_HVKP,
                               method = "lda",
                               metric = "Accuracy",
                               trControl = train_control)
      
      cat("LDA model based on the Wald test results:\n")
      print(model_lda2_OXA48)
      
      test_strat_OXA48$pred_lda2_OXA48 <- predict(model_lda2_OXA48, test_strat_OXA48)
      
      confm_lda2_OXA_48 = table(actual = test_strat_OXA48$OXA48_gene, prediction = test_strat_OXA48$pred_lda2_OXA48)
      
      metrics_lda2_OXA48<- compute_metrics(confm_lda2_OXA_48)
      
      cat("Performance of LDA model based on the Wald test results:\n")
      metrics_lda2_OXA48
    }
    
    else if (model == "SVM-Sig-Chi Squared test"){
      #Fit Support Vector Machine model to train data set
      #This code performs a grid search for tuning the hyper-parameters of the SVM model using the tune function
      # from the e1071 package. It then fits the SVM model with the specified hyper-parameters using the svm function.
      # Model 1
      set.seed(1234)
      tune_out_OXA48 <- e1071::tune("svm", OXA48_gene ~ Amikacin + Cefotaxime + Gentamicin + 
                                      Imipenem + Cefepime + Ceftazidime + Ciprofloxacin + Meropenem + 
                                      Ceftriaxone + Aztreonam
                                    , data = train_strat_OXA48, kernel = "sigmoid", tunecontrol=tune.control(cross=10),
                                    ranges = list(cost = seq(0.1, 5, length = 20)
                                                  , gamma = seq(0.1, 5, length = 20)))
      
      cat("Tuned hyper-parameters throgh 10-Fold Cross Validation:\n")
      print(summary(tune_out_OXA48))
      
      model_svm_OXA48 = svm(OXA48_gene ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                              Cefepime + Ceftazidime + Ciprofloxacin +
                              Meropenem + Ceftriaxone + Aztreonam, data = train_strat_OXA48, 
                             probability = T, kernel = "sigmoid", cost = tune_out_OXA48$best.parameters$cost,
                             gamma = tune_out_OXA48$best.parameters$gamma )
      
      cat("Sigmoid SVM based on Chi-Squared test results:\n")
      print(model_svm_OXA48)
      
      #Prediction on test_strat_OXA48 (SVM)----------------------------------------------------
      test_strat_OXA48$pred_svm0 <- predict(model_svm_OXA48, test_strat_OXA48)
      
      # Model evaluation
      confm0_svm <- table(actual = test_strat_OXA48$OXA48_gene, prediction = test_strat_OXA48$pred_svm0)
      
      metrics_svm1_OXA48_sig<- compute_metrics(confm0_svm)
      
      cat("Performance of SVM-Sigmoid model based on Chi-Squared test results:\n")
      metrics_svm1_OXA48_sig
    }
    else if (model=="SVM-Sig-model agnostic/Wald test"){
      
      set.seed(1234)
      tune_out2_OXA48 <- e1071::tune("svm", OXA48_gene ~ Meropenem + Imipenem
                                     , data = train_strat_OXA48, kernel = "sigmoid", tunecontrol=tune.control(cross=10),
                                     ranges = list(cost = seq(0.1, 5, length = 20)
                                                   , gamma = seq(0.1, 5, length = 20)))
      cat("Tuned hyper-parameters throgh 10-Fold Cross Validation:\n")
      print(summary(tune_out2_OXA48))
      
      model_svm2_OXA48 = svm(OXA48_gene ~ Meropenem + Imipenem, data = train_strat_OXA48, 
                             probability = T, kernel = "sigmoid", cost = tune_out2_OXA48$best.parameters$cost,
                             gamma = tune_out2_OXA48$best.parameters$gamma)
      
      cat("Sigmoid SVM based on model-agnostic approach and Wald test results:\n")
      print(model_svm2_OXA48)
      
      test_strat_OXA48$pred_svm <- predict(model_svm2_OXA48, test_strat_OXA48)
      
      confm_svm <- table(actual = test_strat_OXA48$OXA48_gene, prediction = test_strat_OXA48$pred_svm)
      
      metrics_svm2_OXA48_sig<- compute_metrics(confm_svm)
      
      cat("Performance of SVM-Sigmoid model based on model-agnostic approach and Wald test results:\n")
      metrics_svm2_OXA48_sig
    }
    else if (model=="SVM-Simplicity Principle_MEM"){
      set.seed(1234)
      tune_out3_OXA48 <- e1071::tune("svm", OXA48_gene ~ Meropenem
                                     , data = train_strat_OXA48, kernel = "sigmoid", tunecontrol=tune.control(cross=10),
                                     ranges = list(cost = seq(0.1, 5, length = 20)
                                                   , gamma = seq(0.1, 5, length = 20)))
      
      cat("Tuned hyper-parameters throgh 10-Fold Cross Validation:\n")
      print(summary(tune_out3_OXA48))
      
      model_svm3_OXA48 = svm(OXA48_gene ~ Meropenem, data = train_strat_OXA48, 
                             probability = T, kernel = "sigmoid", cost = tune_out3_OXA48$best.parameters$cost,
                             gamma = tune_out3_OXA48$best.parameters$gamma )
      
      cat("Sigmoid SVM based on the Simplicity principle (MEM):\n")
      print(model_svm3_OXA48)
      
      #Prediction on test_strat_OXA48 (SVM)
      test_strat_OXA48$pred_svm2 <- predict(model_svm3_OXA48, test_strat_OXA48)
      
      # Model evaluation
      confm_svm2 <- table(actual = test_strat_OXA48$OXA48_gene, prediction = test_strat_OXA48$pred_svm2)
      
      metrics_svm3_OXA48_sig<- compute_metrics(confm_svm2)
      
      cat("Performance of SVM-Sigmoid model based on the Simplicity principle (MEM) results:\n")
      metrics_svm3_OXA48_sig
    } 
    
    else if (model=="SVM-Poly-Chi Squared test"){
      set.seed(1234)
      tune_out4_OXA48 <- e1071::tune("svm", OXA48_gene ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                                       Cefepime + Ceftazidime + Ciprofloxacin +
                                       Meropenem + Ceftriaxone + Aztreonam  
                                     , data = train_strat_OXA48, kernel = "polynomial", tunecontrol=tune.control(cross=10),
                                     ranges = list(degree = c(1, 2, 3), cost = seq(0.1, 5, length = 20)
                                                   , gamma = seq(0.1, 5, length = 20)))
      
      cat("Tuned hyper-parameters throgh 10-Fold Cross Validation:\n")
      print(summary(tune_out4_OXA48))
      
      # Tuned hyper-parameters are placed in the model.
      model_svm4_OXA48 = svm(OXA48_gene ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                               Cefepime + Ceftazidime + Ciprofloxacin +
                               Meropenem + Ceftriaxone + Aztreonam
                             , data = train_strat_OXA48, probability = T, 
                             kernel = "polynomial", cost = tune_out4_OXA48$best.parameters$cost, 
                             degree = tune_out4_OXA48$best.parameters$degree , gamma = tune_out4_OXA48$best.parameters$gamma )
      
      cat("Polynomial SVM based on Chi-Squared test results:\n")
      print(model_svm4_OXA48)
      
      #Prediction on test (SVM 2)----------------------------------------------------
      test_strat_OXA48$pred_svm4 <- predict(model_svm4_OXA48, test_strat_OXA48)
      
      # Model evaluation
      confm_svm4 <- table(actual = test_strat_OXA48$OXA48_gene, prediction = test_strat_OXA48$pred_svm4)
      
      metrics_svm_OXA48<- compute_metrics(confm_svm4)
      
      cat("Performance of SVM-Polynomial model based on Chi-Squared test results:\n")
      metrics_svm_OXA48
    }
    else if (model=="SVM-Poly-model agnostic/Wald test"){
      
      set.seed(1234)
      tune_out4dalex_OXA48 <- e1071::tune("svm", OXA48_gene ~ Meropenem + Imipenem
                                          , data = train_strat_OXA48, kernel = "polynomial", tunecontrol=tune.control(cross=10),
                                          ranges = list(cost = seq(0.1, 5, length = 20), degree = c(1, 2, 3)
                                                        , gamma = seq(0.1, 5, length = 20)))
      
      cat("Tuned hyper-parameters throgh 10-Fold Cross Validation:\n")
      print(summary(tune_out4dalex_OXA48))
      
      model_svm4dalex_OXA48 = svm(OXA48_gene ~ Meropenem + Imipenem, data = train_strat_OXA48, 
                             probability = T, kernel = "polynomial", cost = tune_out4dalex_OXA48$best.parameters$cost,
                             degree = tune_out4dalex_OXA48$best.parameters$degree , gamma = tune_out4dalex_OXA48$best.parameters$gamma )
      
      cat("Polynomial SVM based on model-agnostic approach / Wald test results:\n")
      print(model_svm4dalex_OXA48)
      
      #Prediction on test_strat_OXA48 (SVM)----------------------------------------------------
      test_strat_OXA48$pred_svm4dalex <- predict(model_svm4dalex_OXA48, test_strat_OXA48)
      
      # Model evaluation
      confm_svm4dalex <- table(actual = test_strat_OXA48$OXA48_gene, prediction = test_strat_OXA48$pred_svm4dalex)
      
      metrics_svm1_OXA48<- compute_metrics(confm_svm4dalex)
      
      cat("Performance of SVM-Polynomial model based on model-agnostic approach / Wald test results:\n")
      metrics_svm1_OXA48
      
    }
    
    else if (model=="SVM-Poly-Simplicity principle_MEM"){
      
      set.seed(1234)
      tune_out5_OXA48 <- e1071::tune("svm", OXA48_gene ~ Meropenem
                                     , data = train_strat_OXA48, kernel = "polynomial", tunecontrol=tune.control(cross=10),
                                     ranges = list(cost = seq(0.1, 5, length = 20), degree = c(1, 2, 3)
                                                   , gamma = seq(0.1, 5, length = 20)))
      
      cat("Tuned hyper-parameters throgh 10-Fold Cross Validation:\n")
      print(summary(tune_out5_OXA48))
      
      model_svm5_OXA48 = svm(OXA48_gene ~ Meropenem , data = train_strat_OXA48, 
                             probability = T, kernel = "polynomial", cost = tune_out5_OXA48$best.parameters$cost,
                             degree = tune_out5_OXA48$best.parameters$degree , gamma = tune_out5_OXA48$best.parameters$gamma )
      
      cat("Polynomial SVM based on the Simplicity Principle (MEM):\n")
      print(model_svm5_OXA48)
      
      #Prediction on test_strat_OXA48 (SVM)----------------------------------------------------
      test_strat_OXA48$pred_svm5 <- predict(model_svm5_OXA48, test_strat_OXA48)
      
      # Model evaluation
      confm_svm5 <- table(actual = test_strat_OXA48$OXA48_gene, prediction = test_strat_OXA48$pred_svm5)
      
      metrics_svm3_OXA48<- compute_metrics(confm_svm5)
      
      cat("Performance of SVM-Polynomial model based on the Simplicity Principle (MEM):\n")
      metrics_svm3_OXA48
    }
    
    else if (model=="CatBoost-Chi Squared test"){
      # CatBoost Models
      Bio_Data <- as.data.frame(Bio_Data)
      
      #Convert categorical variables to factor
      categorical_var <- c(colnames(Bio_Data))
      
      Bio_Data[, categorical_var[-1]] <- lapply(Bio_Data[, categorical_var[-1]], factor)
      
      levels(Bio_Data$Ampicillin) <- c(levels(Bio_Data$Ampicillin),0)
      
      #Divide Data set into train and test
      set.seed(123)
      
      split_strat  <- initial_split(Bio_Data, prop = 0.8, strata = "OXA48_gene")
      
      train_strat_OXA48  <- training(split_strat)
      
      test_strat_OXA48   <- testing(split_strat)
      
      
      grid <- expand.grid(
        depth = c(2,3,4,5),
        learning_rate = c(0.1,0.2,0.3),
        l2_leaf_reg = c(2,3,4),
        rsm = 1,
        border_count = 1,
        iterations = seq(10,200,10)
      )
      fitControl <- trainControl(method = "cv",
                                       number = 10,
                                       classProbs = TRUE)
      
      set.seed(1234)
      model_catboost_OXA48 <- train(x = train_strat_OXA48[,c("Amikacin", "Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                             "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")],
                                    y = make.names(train_strat_OXA48[,"OXA48_gene"]),
                                    maximize = TRUE,
                                    method = catboost.caret, metric = "Accuracy", 
                                    tuneGrid =  grid,
                                    trControl = fitControl)
      
      cat("Tuned hyper-parameters throgh 10-Fold Cross Validation:\n")
      print(model_catboost_OXA48$bestTune)  
      
      #Prediction on test_strat_OXA48 (CatBoost)----------------------------------------------------
      test_strat_OXA48$catOXA48 = predict(model_catboost_OXA48, test_strat_OXA48)
      
      # Model evaluation
      cat_test_strat_OXA48 <- table(actual = test_strat_OXA48$OXA48_gene, prediction = test_strat_OXA48$catOXA48)
      
      metrics_cat_OXA48<- compute_metrics(cat_test_strat_OXA48)
      
      cat("Performance of CatBoost model based on Chi-Squared test results:\n")
      metrics_cat_OXA48
    }
    
    else if (model=="CatBoost-model agnostic/Wald test"){
      # CatBoost Models
      Bio_Data <- as.data.frame(Bio_Data)
      
      #Convert categorical variables to factor
      categorical_var <- c(colnames(Bio_Data))
      
      Bio_Data[, categorical_var[-1]] <- lapply(Bio_Data[, categorical_var[-1]], factor)
      
      levels(Bio_Data$Ampicillin) <- c(levels(Bio_Data$Ampicillin),0)
      
      #Divide Data set into train and test
      set.seed(123)
      
      split_strat  <- initial_split(Bio_Data, prop = 0.8, strata = "OXA48_gene")
      
      train_strat_OXA48  <- training(split_strat)
      
      test_strat_OXA48   <- testing(split_strat)
      
      grid2 <- expand.grid(
        depth = c(2,3,4),
        learning_rate = c(0.1,0.2,0.3),
        l2_leaf_reg = c(1,2,3),
        rsm = 1,
        border_count = 1,
        iterations = seq(10,100,5)
      )
      fitControl2 <- trainControl(method = "cv",
                                  number = 10,
                                  classProbs = TRUE)
      set.seed(1234)
      model_catboost2_OXA48 <- train(x = train_strat_OXA48[,c("Imipenem", "Meropenem")],
                                     y = make.names(train_strat_OXA48[,"OXA48_gene"]),
                                     maximize = TRUE,
                                     method = catboost.caret, metric = "Accuracy", 
                                     tuneGrid =  grid2,
                                     trControl = fitControl2)
      
      cat("Tuned hyper-parameters throgh 10-Fold Cross Validation:\n")
      print(model_catboost2_OXA48$bestTune)  
      
      test_strat_OXA48$cat1 = predict(model_catboost2_OXA48, test_strat_OXA48)
      
      # Model Evaluation
      confm_cat <- table(actual = test_strat_OXA48$OXA48_gene, prediction = test_strat_OXA48$cat1)
      
      metrics_cat1_OXA48<- compute_metrics(confm_cat)
      
      cat("Performance of CatBoost model based on model-agnostic approach results:\n")
      metrics_cat1_OXA48
    }
    
    else if (model=="CatBoost-Simplicity principle_MEM"){
      # CatBoost Models
      Bio_Data <- as.data.frame(Bio_Data)
      
      #Convert categorical variables to factor
      categorical_var <- c(colnames(Bio_Data))
      
      Bio_Data[, categorical_var[-1]] <- lapply(Bio_Data[, categorical_var[-1]], factor)
      
      levels(Bio_Data$Ampicillin) <- c(levels(Bio_Data$Ampicillin),0)
      
      #Divide Data set into train and test
      set.seed(123)
      
      split_strat  <- initial_split(Bio_Data, prop = 0.8, strata = "OXA48_gene")
      
      train_strat_OXA48  <- training(split_strat)
      
      test_strat_OXA48   <- testing(split_strat)
      
      grid3 <- expand.grid(
        depth = c(2,3,4),
        learning_rate = c(0.1,0.2,0.3),
        l2_leaf_reg = c(1,2,3),
        rsm = 1,
        border_count = 1,
        iterations = seq(10,100,10)
      )
      
      fitControl3 <- trainControl(method = "cv",
                                  number = 10,
                                  classProbs = TRUE)
      
      model_catboost3 <- train(x = train_strat_OXA48["Meropenem"],
                               y = make.names(train_strat_OXA48[["OXA48_gene"]]),
                               maximize = TRUE,
                               method = catboost.caret, metric = "Accuracy", 
                               tuneGrid =  grid3,
                               trControl = fitControl3)
      
      
      cat("Tuned hyper-parameters throgh 10-Fold Cross Validation:\n")
      print(model_catboost2_OXA48$bestTune)  
      
      test_strat_OXA48$cat = predict(model_catboost3, test_strat_OXA48)
      # Model evaluation
      confm_cat2 <- table(actual = test_strat_OXA48$OXA48_gene, prediction = test_strat_OXA48$cat)
      
      metrics_cat3_OXA48<- compute_metrics(confm_cat2)
      
      cat("Performance of CatBoost model based on the simplicity principle (MEM):\n")
      metrics_cat3_OXA48
    }
  }
  
  else if (task == 'feature_selectio') {
    # Task: Feature Selection
    # You can add flags for specific feature selection methods to use, e.g., --feature_selection=dalex or --feature_selection=da.
    # Perform the selected feature selection method based on the provided flags.
    # Example: If user passes --feature_selection=dalex, perform feature selection using DALEX.
    # If user passes --feature_selection=da, perform feature selection using model-agnostic approach.
    # ...
    if (method == "Chi-Squared test"){
      # Function to perform chi-square test and print p-value
      perform_chi_square <- function(data, variable) {
        result <- chisq.test(table(data[[variable]], data$OXA48_gene))
        cat("p-value for", variable, ":", result$p.value, "\n")
      }
      
      # Perform chi-square test for multiple variables
      variables <- c("Amikacin", "Cefotaxime", "Gentamicin", "Imipenem", "Cefepime",
                     "Ceftazidime", "Ciprofloxacin", "Meropenem", "Ampicillin",
                     "Ceftriaxone", "Aztreonam")
      
      for (variable in variables) {
        perform_chi_square(Bio_Data, variable)
      } 
    }
  }
  
  else if (method == "Wald test"){
    model_Log1_OXA48 <- glm(OXA48_gene ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                              Cefepime + Ceftazidime + Ciprofloxacin +
                              Meropenem + Ceftriaxone + Aztreonam 
                            , family = "binomial", data = train_strat_OXA48)
    
    cat("P-values show Wald test results:\n")
    print(summary(model_Log1_OXA48))
  }
  else if (method == "model agnostic-LR"){
    # Feature Selection by DALEX Package--------------------------------------
    # the 'explained_glm_TEM' object is created using the 'explain' function from the DALEX package. 
    # The variable_importance function is used to calculate the permutation-based 
    # feature importance ('fi_glm_TEM') with 50 permutations. Finally, the feature importance is displayed 
    # using 'fi_glm_TEM' and plotted using' plot(fi_glm_TEM)'.
    model_Log1_OXA48 <- glm(OXA48_gene ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                              Cefepime + Ceftazidime + Ciprofloxacin +
                              Meropenem + Ceftriaxone + Aztreonam 
                            , family = "binomial", data = train_strat_OXA48)
    
    explained_glm_OXA48 <- explain(model = model_Log1_OXA48, data=train_strat_OXA48[,c(2:9,11:12)], variables = colnames(train_strat_OXA48[,c(2:9,11:12)]),
                                   y=as.vector(as.numeric(train_strat_OXA48$OXA48_gene))-1, label = "LR", 
                                   type = "classification")
    #50 permuatation
    fi_glm_OXA48 = variable_importance(explained_glm_OXA48, B=50 ,variables = colnames(train_strat_OXA48[,c(2:9,11:12)]), loss_function = loss_root_mean_square, type = "raw")
    
    cat("Root mean squared dropout loss due to a feature elimination through model-agnostic approach in LR:\n")
    print(fi_glm_OXA48)
  }
  else if (method == "model agnostic-NBC"){
    train_control <- trainControl(
      method = "cv",
      number = 10
    )
    
    search_grid_OXA48 <- expand.grid(
      usekernel = c(FALSE,TRUE),
      fL = seq(0,5,0.5),
      adjust = seq(0,5,0.5)
    )
    
    model_nb_OXA48 <- train(OXA48_gene ~  Amikacin + Cefotaxime + Gentamicin + Imipenem +
                              Cefepime + Ceftazidime + Ciprofloxacin +
                              Meropenem + Ceftriaxone + Aztreonam,
                            data = train_strat_OXA48,
                            method = "nb",
                            metric = "Accuracy",
                            trControl = train_control,
                            tuneGrid = search_grid_OXA48)
    # Feature selection through model-agnostic approach (DALEX)
    # the 'explain' function from the DALEX package is used to explain the Naive Bayes classifier model 
    # 'model_nbdalex_TEM'. The relevant variables are selected and provided as data along with their 
    # corresponding column names. The target variable 'OXA48_gene' is transformed to numeric values 
    # (y = as.vector(as.numeric(train_strat_OXA48$OXA48_gene)) - 1). The label is set to "NBC" for Naive Bayes classifier,
    # and the type is set to "classification".
    explained_nb_OXA48 <- explain(model = model_nb_OXA48, data=train_strat_OXA48[,c(2:9,11:12)], variables = colnames(train_strat_OXA48[,c(2:9,11:12)]),
                                  y=as.vector(as.numeric(train_strat_OXA48$OXA48_gene))-1, label = "NBC", 
                                  type = "classification")
    
    # The 'variable_importance' function is then used to calculate the variable importance based on the
    # explained model. The 'B' parameter is set to 50 for the number of permutations, and the loss function
    # is set to "loss_root_mean_square". The resulting variable importance is stored in 'fitnb_TEM'
    fi_nb_OXA48 = variable_importance(explained_nb_OXA48, B=50 ,variables = colnames(train_strat_OXA48[,c(2:9,11:12)]), loss_function = loss_root_mean_square, type = "raw")
    
    cat("Root mean squared dropout loss due to a feature elimination through model-agnostic approach in NBC:\n")
    print(fi_nb_OXA48)
  }
  else if (method == "model agnostic-LDA"){
    model_ldadalex_OXA48 <- train(OXA48_gene ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                                    Cefepime + Ceftazidime + Ciprofloxacin +Aztreonam+
                                    Meropenem + Ceftriaxone  ,
                                  data = train_strat_OXA48,
                                  method = "lda",
                                  metric = "Accuracy")
    
    explda_OXA48 = explain(model = model_ldadalex_OXA48, data=train_strat_OXA48[,c("Amikacin", "Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                                   "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")], variables = colnames(train_strat_OXA48[,c("Amikacin", "Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                                                                                                                                                       "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")]),
                           y=as.vector(as.numeric(train_strat_OXA48$OXA48_gene))-1, label = "LDA", 
                           type = "classification")
    
    fitlda_OXA48 = variable_importance(explda_OXA48, B=50 ,variables = colnames(train_strat_OXA48[,c("Amikacin", "Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                                                     "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")]), loss_function = loss_root_mean_square, type = "raw" )
    cat("Root mean squared dropout loss due to a feature elimination through model-agnostic approach in LDA:\n")
    print(fitlda_OXA48)
  }
  else if (method == "model agnostic-Sig SVM"){
    set.seed(1234)
    tune_out_OXA48 <- e1071::tune("svm", OXA48_gene ~ Amikacin + Cefotaxime + Gentamicin + 
                                    Imipenem + Cefepime + Ceftazidime + Ciprofloxacin + Meropenem + 
                                    Ceftriaxone + Aztreonam
                                  , data = train_strat_OXA48, kernel = "sigmoid", tunecontrol=tune.control(cross=10),
                                  ranges = list(cost = seq(0.1, 5, length = 20)
                                                , gamma = seq(0.1, 5, length = 20)))
    

    model_svm_OXA48 = svm(OXA48_gene ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                            Cefepime + Ceftazidime + Ciprofloxacin +
                            Meropenem + Ceftriaxone + Aztreonam, data = train_strat_OXA48, 
                          probability = T, kernel = "sigmoid", cost = tune_out_OXA48$best.parameters$cost,
                          gamma = tune_out_OXA48$best.parameters$gamma )
    
    explained_SVM_OXA48 <- explain(model = model_svm_OXA48, data=train_strat_OXA48[,c(2:9,11:12)], variables = colnames(train_strat_OXA48[,c(2:9,11:12)]),
                                   y=as.vector(as.numeric(train_strat_OXA48$OXA48_gene))-1, label = "SVM-Sigmoid",
                                   type = "classification")
    
    #50 permuatation
    fi_SVM_OXA48 = variable_importance(explained_SVM_OXA48, B=50,variables = colnames(train_strat_OXA48[,c(2:9,11:12)]), loss_function = loss_root_mean_square, type = "raw" )
    
    cat("Root mean squared dropout loss due to a feature elimination through model-agnostic approach in sigmoid kernel of SVM:\n")
    print(fi_SVM_OXA48)
  }
  
  else if (method == "model agnostic-Poly SVM"){
    # it is highly recommeded that first train the model, then specify following hyper-parameters according to the tuned hyper-parameters.
    # Following values of the hyper-parameters have been specified in this way.
    set.seed(1234)
    tune_out4_OXA48 <- e1071::tune("svm", OXA48_gene ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                                     Cefepime + Ceftazidime + Ciprofloxacin +
                                     Meropenem + Ceftriaxone + Aztreonam  
                                   , data = train_strat_OXA48, kernel = "polynomial", tunecontrol=tune.control(cross=10),
                                   ranges = list(degree = c(1, 2, 3), cost = seq(0.1, 5, length = 20)
                                                 , gamma = seq(0.1, 5, length = 20)))
    

    # Tuned hyper-parameters are placed in the model.
    model_svm4_OXA48 = svm(OXA48_gene ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                             Cefepime + Ceftazidime + Ciprofloxacin +
                             Meropenem + Ceftriaxone + Aztreonam
                           , data = train_strat_OXA48, probability = T, 
                           kernel = "polynomial", cost = tune_out4_OXA48$best.parameters$cost, 
                           degree = tune_out4_OXA48$best.parameters$degree , gamma = tune_out4_OXA48$best.parameters$gamma )
    
    explained_SVM2_OXA48 <- explain(model = model_svm4_OXA48, data=train_strat_OXA48[,c(2:9,11:12)], variables = colnames(train_strat_OXA48[,c(2:9,11:12)]),
                                    y=as.vector(as.numeric(train_strat_OXA48$OXA48_gene))-1, label = "SVM-Polynomial",
                                    type = "classification")
    
    #50 permuatation
    fi_SVM_OXA482 = variable_importance(explained_SVM2_OXA48, B=50,variables = colnames(train_strat_OXA48[,c(2:9,11:12)]), loss_function = loss_root_mean_square, type = "raw" )
    
    cat("Root mean squared dropout loss due to a feature elimination through model-agnostic approach in polynomial kernel of SVM:\n")
    print(fi_SVM_OXA482)
  }
  else if (method == "model agnostic-CatBoost"){
    
    # CatBoost Models
    Bio_Data <- as.data.frame(Bio_Data)
    
    #Convert categorical variables to factor
    categorical_var <- c(colnames(Bio_Data))
    
    Bio_Data[, categorical_var[-1]] <- lapply(Bio_Data[, categorical_var[-1]], factor)
    
    levels(Bio_Data$Ampicillin) <- c(levels(Bio_Data$Ampicillin),0)
    
    #Divide Data set into train and test
    set.seed(123)
    
    split_strat  <- initial_split(Bio_Data, prop = 0.8, strata = "OXA48_gene")
    
    train_strat_OXA48  <- training(split_strat)
    
    test_strat_OXA48   <- testing(split_strat)
    
    
    grid <- expand.grid(
      depth = c(2,3,4,5),
      learning_rate = c(0.1,0.2,0.3),
      l2_leaf_reg = c(2,3,4),
      rsm = 1,
      border_count = 1,
      iterations = seq(10,200,10)
    )
    fitControl <- trainControl(method = "cv",
                               number = 10,
                               classProbs = TRUE)
    
    set.seed(1234)
    model_catboost_OXA48 <- train(x = train_strat_OXA48[,c("Amikacin", "Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                           "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")],
                                  y = make.names(train_strat_OXA48[,"OXA48_gene"]),
                                  maximize = TRUE,
                                  method = catboost.caret, metric = "Accuracy", 
                                  tuneGrid =  grid,
                                  trControl = fitControl)
    
    expsda_OXA48 = explain(model = model_catboost_OXA48, data=train_strat_OXA48[,c("Amikacin", "Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                                   "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")], variables = colnames(train_strat_OXA48[,c("Amikacin", "Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                                                                                                                                                       "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")]),
                           y=as.vector(as.numeric(train_strat_OXA48$OXA48_gene))-1, label = "CatBoost", 
                           type = "classification")
    
    fitcat_OXA48 = variable_importance(expsda_OXA48, B=50 ,variables = colnames(train_strat_OXA48[,c("Amikacin", "Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                                                     "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")]), loss_function = loss_root_mean_square, type = "raw" )
    
    cat("Root mean squared dropout loss due to a feature elimination through model-agnostic approach in CatBoost:\n")
    print(fitcat_OXA48)
  }
  else {
    cat("Invalid task. Please specify one of the following tasks: feature_selectio, model_evaluation\n")
    
  } 
}
if (is.null(task)){
  cat("You should specify a task\n")
}

