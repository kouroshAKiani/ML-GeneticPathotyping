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

# Calculate the proportion of CTX_M15_gene
proportion_CTX_M15_gene <- prop.table(table(Bio_Data$CTX_M15_gene))
cat("Proportion of CTX_M15_gene:\n")
print(proportion_CTX_M15_gene)


# Divide the dataset into train and test sets
set.seed(123)
train_cases <- sample(1:nrow(Bio_Data), nrow(Bio_Data) * 0.8)
train <- Bio_Data[train_cases,]
test <- Bio_Data[-train_cases,]

# Display the dimensions of the train and test sets
cat("Train set dimensions:", dim(train), "\n")
cat("Test set dimensions:", dim(test), "\n")


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
      model_Log1_CTX_M15 <- glm(CTX_M15_gene ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                                  Cefepime + Ceftazidime + Ciprofloxacin +
                                  Meropenem + Ceftriaxone + Aztreonam 
                                , family = "binomial", data = train)
      
      cat("Summary of LR model based on Chi Squared test including Null deviance, Residual deviance, and AIC, etc.:\n")
      print(summary(model_Log1_CTX_M15))
      
      test$probs1_CTX_M15 <- predict(model_Log1_CTX_M15, test, type = "response")
      
      test$pred_logreg1_CTX_M15 <- ifelse(test$probs1_CTX_M15 >= 0.5, 1, 0)
      
      # Model 1 evaluation
      confm_logreg1_CTX_M15 <- table(actual = test$CTX_M15_gene, prediction = test$pred_logreg1_CTX_M15)
      
      metrics_logreg1_CTX_M15 <- compute_metrics(confm_logreg1_CTX_M15)
      
      cat("Performance of LR model based on Chi-Squared test results:\n")
      metrics_logreg1_CTX_M15
    }
    else if (model == "LR-Wald test"){
      model_Log2_CTX_M15 <- glm(CTX_M15_gene ~ Cefepime + Cefotaxime
                            , family = "binomial", data = train)
      cat("Summary of LR model based on Wald test including Null deviance, Residual deviance, and AIC, etc.:\n")
      print(summary(model_Log2_CTX_M15))
      
      ##Prediction on test (Model 2)--------------------------------------------------
      test$probs2_CTX_M15 <- predict(model_Log2_CTX_M15, test, type = "response")
      
      test$pred_logreg2_CTX_M15 <- ifelse(test$probs2_CTX_M15 >= 0.5, 1, 0)
      
      # Model 2 evaluation
      confm_logreg2_CTX_M15 <- table(actual = test$CTX_M15_gene, prediction = test$pred_logreg2_CTX_M15)
      
      metrics_logreg2_CTX_M15 <- compute_metrics(confm_logreg2_CTX_M15)
      
      cat("Performance of LR model based on Wald test results:\n")
      metrics_logreg2_CTX_M15
    }
    else if (model == "LR-model agnostic"){
      model_Log3_CTX_M15 <- glm(CTX_M15_gene ~ Cefotaxime
                                , family = "binomial", data = train)
      
      cat("Summary of LR model based on model-agnostic approach including Null deviance, Residual deviance, and AIC, etc.:\n")
      print(summary(model_Log3_CTX_M15))
      
      test$probs3_CTX_M15 <- predict(model_Log3_CTX_M15, test, type = "response")

      test$pred_logreg3_CTX_M15 <- ifelse(test$probs3_CTX_M15 >= 0.5, 1, 0)
      
      confm_logreg3_CTX_M15 <- table(actual = test$CTX_M15_gene, prediction = test$pred_logreg3_CTX_M15)
      
      metrics_logreg3_CTX_M15 <- compute_metrics(confm_logreg3_CTX_M15)
      
      cat("Performance of LR model based on model-agnostic approach result:\n")
      metrics_logreg3_CTX_M15
    }
    else if (model == "LR-Simplicity Principle_FEP"){
      model_Log4_CTX_M15 <- glm(CTX_M15_gene ~ Cefepime
                                , family = "binomial", data = train)
      
      cat("Summary of LR model based on simplicity principle (FEP) including Null deviance, Residual deviance, and AIC, etc.:\n")
      print(summary(model_Log4_CTX_M15))
      
      ##Prediction on test
      test$probs4_CTX_M15 <- predict(model_Log4_CTX_M15, test, type = "response")

      test$pred_logreg4_CTX_M15 <- ifelse(test$probs4_CTX_M15 >= 0.5, 1, 0)
      
      confm_logreg4_CTX_M15 <- table(actual = test$CTX_M15_gene, prediction = test$pred_logreg4_CTX_M15)
      
      metrics_logreg4_CTX_M15 <- compute_metrics(confm_logreg4_CTX_M15)
      
      cat("Performance of LR model based on simplicity principle (FEP):\n")
      metrics_logreg4_CTX_M15
    }
    else if (model == "LR-Simplicity Principle_CIP"){
      model_Log5_CTX_M15 <- glm(CTX_M15_gene ~ Ciprofloxacin
                                , family = "binomial", data = train)
      
      cat("Summary of LR model based on simplicity principle (CIP) including Null deviance, Residual deviance, and AIC, etc.:\n")
      print(summary(model_Log5_CTX_M15))
      
      ##Prediction on test
      test$probs5_CTX_M15 <- predict(model_Log5_CTX_M15, test, type = "response")

      test$pred_logreg5_CTX_M15 <- ifelse(test$probs5_CTX_M15 >= 0.5, 1, 0)
      
      confm_logreg5_CTX_M15 <- table(actual = test$CTX_M15_gene, prediction = test$pred_logreg5_CTX_M15)
      
      metrics_logreg5_CTX_M15 <- compute_metrics(confm_logreg5_CTX_M15)
      
      cat("Performance of LR model based on simplicity principle (CIP):\n")
      metrics_logreg5_CTX_M15
    }
    else if (model == "NBC-Chi Squared test"){
      train_control <- trainControl(
        method = "cv", 
        number = 10
      )
      
      search_grid_CTX_M15 <- expand.grid(
        usekernel = c(FALSE,TRUE),
        fL = seq(0,5,0.5),
        adjust = seq(0,5,0.5)
      )
      
      
      model_nbdalex_CTX_M15 <- train(CTX_M15_gene ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                                       Cefepime + Ceftazidime + Ciprofloxacin +
                                       Meropenem + Ceftriaxone + Aztreonam,
                                     data = train,
                                     method = "nb",
                                     metric = "Accuracy",
                                     trControl = train_control,
                                     tuneGrid = search_grid_CTX_M15)
      
      cat("Tuned hyper-parameters through 10-fold Cross-Validation:\n")
      print(model_nbdalex_CTX_M15$bestTune)
      
      ##Prediction on test (Naive Bayes Model)----------------------------------
      test$pred_nb_CTX_M15 <- predict(model_nbdalex_CTX_M15, test)
      
      # Model evaluation
      confm_nb_CTX_M15 <- table(actual = test$CTX_M15_gene, prediction = test$pred_nb_CTX_M15)
      metrics_nb_CTX_M15<- compute_metrics(confm_nb_CTX_M15)
      cat("Performance of NBC model based on Chi-Squared test results:\n")
      metrics_nb_CTX_M15
    }
    else if (model == "NBC-model agnostic"){
      train_control <- trainControl(
        method = "cv",
        number = 10
      )
      
      search_grid_CTX_M15 <- expand.grid(
        usekernel = c(FALSE,TRUE),
        fL = seq(0,5,0.5),
        adjust = seq(0,5,0.5)
      )
      
      model_nbdalex1_CTX_M15 <- train(CTX_M15_gene ~ Cefepime + Ciprofloxacin,
                                      data = train,
                                      method = "nb",
                                      metric = "Accuracy",
                                      trControl = train_control,
                                      tuneGrid = search_grid_CTX_M15)
      
      cat("Tuned hyper-parameters through 10-fold Cross-Validation:\n")
      print(model_nbdalex1_CTX_M15$bestTune)
      
      ##Prediction on test (Naive Bayes Model)----------------------------------
      test$pred_nb1_CTX_M15 <- predict(model_nbdalex1_CTX_M15, test)
      
      # Model evaluation
      confm_nb1_CTX_M15 <- table(actual = test$CTX_M15_gene, prediction = test$pred_nb1_CTX_M15)
      metrics_nb1_CTX_M15<- compute_metrics(confm_nb1_CTX_M15)
      cat("Performance of NBC model based on model-agnostic approach results:\n")
      metrics_nb1_CTX_M15
    }
    
    else if (model == "NBC-Wald test"){
      
      train_control <- trainControl(
        method = "cv",
        number = 10
      )
      
      search_grid_CTX_M15 <- expand.grid(
        usekernel = c(FALSE,TRUE),
        fL = seq(0,5,0.5),
        adjust = seq(0,5,0.5)
      )
      
      model_nb2_CTX_M15 <- train(CTX_M15_gene ~ Cefepime + Cefotaxime,
                                 data = train,
                                 method = "nb",
                                 metric = "Accuracy",
                                 trControl = train_control,
                                 tuneGrid = search_grid_CTX_M15)
      
      cat("Tuned hyper-parameters through 10-fold Cross-Validation:\n")
      print(model_nb2_CTX_M15$bestTune)
      
      ##Prediction on test (Naive Bayes Model)
      test$pred_nb2_CTX_M15 <- predict(model_nb2_CTX_M15, test)
      
      ##confusion matrix(Naive Bayes Model)
      confm_nb2_CTX_M15 <- table(actual = test$CTX_M15_gene, prediction = test$pred_nb2_CTX_M15)
      
      metrics_nb2_CTX_M15<- compute_metrics(confm_nb2_CTX_M15)
      
      cat("Performance of NBC model based on Wald test results:\n")
      metrics_nb2_CTX_M15
    }
    
    else if (model == "NBC-Simplicity Principle_CTX"){
      
      train_control <- trainControl(
        method = "cv",
        number = 10
      )
      
      search_grid_CTX_M15 <- expand.grid(
        usekernel = c(FALSE,TRUE),
        fL = seq(0,5,0.5),
        adjust = seq(0,5,0.5)
      )
      
      model_nb3_CTX_M15 <- train(CTX_M15_gene ~ Cefotaxime,
                                 data = train,
                                 method = "nb",
                                 metric = "Accuracy",
                                 trControl = train_control,
                                 tuneGrid = search_grid_CTX_M15)
      
      cat("Tuned hyper-parameters through 10-fold Cross-Validation:\n")
      print(model_nb3_CTX_M15$bestTune)
      
      ##Prediction on test (Naive Bayes Model)
      test$pred_nb3_CTX_M15 <- predict(model_nb3_CTX_M15, test)
      
      ##confusion matrix(Naive Bayes Model)
      confm_nb3_CTX_M15 <- table(actual = test$CTX_M15_gene, prediction = test$pred_nb3_CTX_M15)
      
      metrics_nb3_CTX_M15<- compute_metrics(confm_nb3_CTX_M15)
      
      cat("Performance of NBC model based on Simplicity Principle (CTX):\n")
      metrics_nb3_CTX_M15
    }
    
    else if (model == "NBC-Simplicity Principle_FEP"){
      
      train_control <- trainControl(
        method = "cv",
        number = 10
      )
      
      search_grid_CTX_M15 <- expand.grid(
        usekernel = c(FALSE,TRUE),
        fL = seq(0,5,0.5),
        adjust = seq(0,5,0.5)
      )
      
      model_nb4_CTX_M15 <- train(CTX_M15_gene ~ Cefepime,
                                 data = train,
                                 method = "nb",
                                 metric = "Accuracy",
                                 trControl = train_control,
                                 tuneGrid = search_grid_CTX_M15)
      
      cat("Tuned hyper-parameters through 10-fold Cross-Validation:\n")
      print(model_nb4_CTX_M15$bestTune)
      
      ##Prediction on test (Naive Bayes Model)
      test$pred_nb4_CTX_M15 <- predict(model_nb4_CTX_M15, test)
      
      ##confusion matrix(Naive Bayes Model)
      confm_nb4_CTX_M15 <- table(actual = test$CTX_M15_gene, prediction = test$pred_nb4_CTX_M15)
      
      metrics_nb4_CTX_M15<- compute_metrics(confm_nb4_CTX_M15)
      
      cat("Performance of NBC model based on Simplicity Principle (FEP):\n")
      metrics_nb4_CTX_M15
    }
    
    else if (model == "NBC-Simplicity Principle_CIP"){
      
      train_control <- trainControl(
        method = "cv",
        number = 10
      )
      
      search_grid_CTX_M15 <- expand.grid(
        usekernel = c(FALSE,TRUE),
        fL = seq(0,5,0.5),
        adjust = seq(0,5,0.5)
      )
      
      model_nb5_CTX_M15 <- train(CTX_M15_gene ~ Ciprofloxacin,
                                 data = train,
                                 method = "nb",
                                 metric = "Accuracy",
                                 trControl = train_control,
                                 tuneGrid = search_grid_CTX_M15)
      
      cat("Tuned hyper-parameters through 10-fold Cross-Validation:\n")
      print(model_nb5_CTX_M15$bestTune)
      
      ##Prediction on test (Naive Bayes Model)
      test$pred_nb5_CTX_M15 <- predict(model_nb5_CTX_M15, test)
      
      ##confusion matrix(Naive Bayes Model)
      confm_nb5_CTX_M15 <- table(actual = test$CTX_M15_gene, prediction = test$pred_nb5_CTX_M15)
      
      metrics_nb5_CTX_M15<- compute_metrics(confm_nb5_CTX_M15)
      
      cat("Performance of NBC model based on Simplicity Principle (CIP):\n")
      metrics_nb5_CTX_M15
    }
    
    else if (model == "LDA-Chi Squared test"){
      
      model_ldadalex_CTX_M15 <- train(CTX_M15_gene ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                                        Cefepime + Ceftazidime + Ciprofloxacin +
                                        Meropenem + Ceftriaxone + Aztreonam ,
                                      data = train,
                                      method = "lda",
                                      metric = "Accuracy")
      
      cat("LDA model based on Chi-Squared test results:\n")
      print(model_ldadalex_CTX_M15)
      
      #Prediction on test (LDA Model)-----------------------------------------------
      test$pred_lda_CTX_M15 <- predict(model_ldadalex_CTX_M15, test)
      
      # Model Evaluation-------------------------------------------------
      confm_lda_CTX_M15 <- table(actual = test$CTX_M15_gene, prediction = test$pred_lda_CTX_M15)
      metrics_lda_CTX_M15<- compute_metrics(confm_lda_CTX_M15)
      cat("Performance of LDA model based on Chi-Squared test results:\n")
      metrics_lda_CTX_M15
    }
    else if (model == "LDA-model agnostic/Wald test"){
      model_ldadalex1_CTX_M15 <- train(CTX_M15_gene ~ Cefotaxime + Cefepime,
                                       data = train,
                                       method = "lda",
                                       metric = "Accuracy",
                                       trControl = train_control)
      
      cat("LDA model based on model-agnostic approach/Wald test results:\n")
      print(model_ldadalex1_CTX_M15)
      
      #Prediction on test (LDA Model)-----------------------------------------------
      test$pred_lda1_CTX_M15 <- predict(model_ldadalex1_CTX_M15, test)
      
      # Model evaluation
      confm_lda1_CTX_M15 <- table(actual = test$CTX_M15_gene, prediction = test$pred_lda1_CTX_M15)
      metrics_lda1_CTX_M15<- compute_metrics(confm_lda1_CTX_M15)
      cat("Performance of LDA model based on model-agnostic approach/Wald test results:\n")
      metrics_lda1_CTX_M15
    }
    
    else if (model == "LDA-Simplicity Principle_FEP"){
      model_lda2_CTX_M15 <- train(CTX_M15_gene ~ Cefepime,
                                  data = train,
                                  method = "lda",
                                  metric = "Accuracy",
                                  trControl = train_control)
      
      cat("LDA model based on Simplicity Principle (FEP):\n")
      print(model_lda2_CTX_M15)
      
      test$pred_lda2_CTX_M15 <- predict(model_lda2_CTX_M15, test)
      
      confm_lda2_CTX_M15 <- table(actual = test$CTX_M15_gene, prediction = test$pred_lda2_CTX_M15)
      
      metrics_lda2_CTX_M15<- compute_metrics(confm_lda2_CTX_M15)
      
      cat("Performance of LDA model based on Simplicity Principle (FEP):\n")
      metrics_lda2_CTX_M15
    }
    
    else if (model == "LDA-Simplicity Principle_CTX"){
      
      model_lda3_CTX_M15 <- train(CTX_M15_gene ~ Cefotaxime,
                                  data = train,
                                  method = "lda",
                                  metric = "Accuracy",
                                  trControl = train_control)
      
      cat("LDA model based on Simplicity Principle (CTX):\n")
      print(model_lda3_CTX_M15)
      
      test$pred_lda3_CTX_M15 <- predict(model_lda3_CTX_M15, test)
      
      confm_lda3_CTX_M15 <- table(actual = test$CTX_M15_gene, prediction = test$pred_lda3_CTX_M15)
      
      metrics_lda3_CTX_M15<- compute_metrics(confm_lda3_CTX_M15)
      
      cat("Performance of LDA model based on Simplicity Principle (CTX):\n")
      metrics_lda3_CTX_M15
    } 
    
    else if (model == "LDA-Simplicity Principle_CIP"){
      
      model_lda4_CTX_M15 <- train(CTX_M15_gene ~ Ciprofloxacin,
                                  data = train,
                                  method = "lda",
                                  metric = "Accuracy")
      
      cat("LDA model based on Simplicity Principle (CIP):\n")
      print(model_lda4_CTX_M15)
      
      test$pred_lda4_CTX_M15 <- predict(model_lda4_CTX_M15, test)
      
      confm_lda4_CTX_M15 <- table(actual = test$CTX_M15_gene, prediction = test$pred_lda4_CTX_M15)
      
      metrics_lda4_CTX_M15<- compute_metrics(confm_lda4_CTX_M15)
      
      cat("Performance of LDA model based on Simplicity Principle (CIP):\n")
      metrics_lda4_CTX_M15
    }
    else if (model == "SVM-Sig-Chi Squared test"){
      #Fit Support Vector Machine model to train data set
      #This code performs a grid search for tuning the hyper-parameters of the SVM model using the tune function
      # from the e1071 package. It then fits the SVM model with the specified hyper-parameters using the svm function.
      # Model 1
      set.seed(1234)
      tune_out_CTX_M15_sig <- e1071::tune("svm", CTX_M15_gene ~ Amikacin + Cefotaxime + Gentamicin + 
                                            Imipenem + Cefepime + Ceftazidime + Ciprofloxacin + Meropenem + 
                                            Ceftriaxone + Aztreonam 
                                          , data = train, kernel = "sigmoid",tunecontrol=tune.control(cross=10),
                                          ranges = list(cost = seq(0.1, 5, length = 20)
                                                        , gamma = seq(0.5, 5, length = 20)))
      
      cat("Tuned hyper-parameters throgh 10-Fold Cross Validation:\n")
      print(summary(tune_out_CTX_M15_sig))
      
      model_svm_CTX_M15_sig = svm(CTX_M15_gene ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                            Cefepime + Ceftazidime + Ciprofloxacin +
                            Meropenem + Ceftriaxone + Aztreonam, data = train, 
                          probability = T, kernel = "sigmoid", cost = tune_out_CTX_M15_sig$best.parameters$cost,
                          gamma = tune_out_CTX_M15_sig$best.parameters$gamma )
      
      cat("Sigmoid SVM based on Chi-Squared test results:\n")
      print(model_svm_CTX_M15_sig)
      
      #Prediction on test (SVM)----------------------------------------------------
      test$pred_svm_CTX_M15_sig <- predict(model_svm_CTX_M15_sig, test)
      
      # Model evaluation
      confm_svm_CTX_M15_sig <- table(actual = test$CTX_M15_gene, prediction = test$pred_svm_CTX_M15_sig)
      
      metrics_svm1_CTX_M15_sig<- compute_metrics(confm_svm_CTX_M15_sig)
      
      cat("Performance of SVM-Sigmoid model based on Chi-Squared test results:\n")
      metrics_svm1_CTX_M15_sig
    }
    else if (model=="SVM-Sig-model agnostic"){
      
      set.seed(1234)
      tune_out1_CTX_M15_sig <- e1071::tune("svm", CTX_M15_gene ~ Imipenem + Aztreonam + Meropenem + Amikacin 
                                           , data = train, kernel = "sigmoid",tunecontrol=tune.control(cross=10),
                                           ranges = list(cost = seq(0.1, 5, length = 20)
                                                         , gamma = seq(0.1, 5, length = 20)))
      cat("Tuned hyper-parameters throgh 10-Fold Cross Validation:\n")
      print(summary(tune_out1_CTX_M15_sig))
      
      model_svm1_CTX_M15_sig = svm(CTX_M15_gene ~ Imipenem + Aztreonam + Meropenem + Amikacin, data = train, 
                           probability = T, kernel = "sigmoid", cost = tune_out1_CTX_M15_sig$best.parameters$cost,
                           gamma = tune_out1_CTX_M15_sig$best.parameters$gamma)
      
      cat("Sigmoid SVM based on model-agnostic approach results:\n")
      print(model_svm1_CTX_M15_sig)
      
      test$pred_svm1_CTX_M15_sig <- predict(model_svm1_CTX_M15_sig, test)
      
      confm_svm1_CTX_M15_sig <- table(actual = test$CTX_M15_gene, prediction = test$pred_svm1_CTX_M15_sig)
      
      metrics_svm2_CTX_M15_sig<- compute_metrics(confm_svm1_CTX_M15_sig)
      
      cat("Performance of SVM-Sigmoid model based on model-agnostic approach results:\n")
      metrics_svm2_CTX_M15_sig
    }
    else if (model=="SVM-Sig-Wald test"){
      set.seed(1234)
      tune_out_CTX_M15_2_sig <- e1071::tune("svm", CTX_M15_gene ~ Cefotaxime + Cefepime
                                            , data = train, kernel = "sigmoid", tunecontrol=tune.control(cross=10),
                                            ranges = list(cost = seq(0.1, 5, length = 20)
                                                          , gamma = seq(0.1, 5, length = 20)))
      
      cat("Tuned hyper-parameters throgh 10-Fold Cross Validation:\n")
      print(summary(tune_out_CTX_M15_2_sig))
      
      model_svm_CTX_M15_2_sig = svm(CTX_M15_gene ~ Cefotaxime + Cefepime, data = train, 
                           probability = T, kernel = "sigmoid", cost = tune_out_CTX_M15_2_sig$best.parameters$cost,
                           gamma = tune_out_CTX_M15_2_sig$best.parameters$gamma )
      
      cat("Sigmoid SVM based on the Wald test result:\n")
      print(model_svm_CTX_M15_2_sig)
      
      #Prediction on test (SVM)
      test$pred_svm_CTX_M15_2_sig <- predict(model_svm_CTX_M15_2_sig, test)
      
      # Model evaluation
      confm_svm2_CTX_M15_sig <- table(actual = test$CTX_M15_gene, prediction = test$pred_svm_CTX_M15_2_sig)
      
      metrics_svm3_CTX_M15_sig<- compute_metrics(confm_svm2_CTX_M15_sig)
      
      cat("Performance of SVM-Sigmoid model based on the Wald test results:\n")
      metrics_svm3_CTX_M15_sig
      } 
    
    else if (model=="SVM-Sig-Simplicity Principle_CTX"){
      set.seed(1234)
      tune_out_CTX_M15_3_sig <- e1071::tune("svm", CTX_M15_gene ~ Cefotaxime
                                            , data = train, kernel = "sigmoid", tunecontrol=tune.control(cross=10),
                                            ranges = list(cost = seq(0.1, 5, length = 20)
                                                          , gamma = seq(0.1, 5, length = 20)))
      
      cat("Tuned hyper-parameters throgh 10-Fold Cross Validation:\n")
      print(summary(tune_out_CTX_M15_3_sig))
      
      model_svm_CTX_M15_3_sig = svm(CTX_M15_gene ~ Cefotaxime, data = train, 
                           probability = T, kernel = "sigmoid", cost = tune_out_CTX_M15_3_sig$best.parameters$cost,
                           gamma = tune_out_CTX_M15_3_sig$best.parameters$gamma )
      
      cat("Sigmoid SVM based on Simplicity Principle (CTX):\n")
      print(model_svm_CTX_M15_3_sig)
      
      #Prediction on test (SVM)
      test$pred_svm_CTX_M15_3_sig <- predict(model_svm_CTX_M15_3_sig, test)
      
      # Model evaluation
      confm_svm3_CTX_M15_sig <- table(actual = test$CTX_M15_gene, prediction = test$pred_svm_CTX_M15_3_sig)
      
      metrics_svm4_CTX_M15_sig<- compute_metrics(confm_svm3_CTX_M15_sig)
      
      cat("Performance of SVM-Sigmoid model based on Simplicity Principle (CTX):\n")
      metrics_svm4_CTX_M15_sig
    } 
    
    else if (model=="SVM-Sig-Simplicity Principle_FEP"){
      set.seed(1234)
      tune_out_CTX_M15_4_sig <- e1071::tune("svm", CTX_M15_gene ~ Cefepime
                                            , data = train, kernel = "sigmoid", tunecontrol=tune.control(cross=10),
                                            ranges = list(cost = seq(0.1, 5, length = 20)
                                                          , gamma = seq(0.1, 5, length = 20)))
      
      cat("Tuned hyper-parameters throgh 10-Fold Cross Validation:\n")
      print(summary(tune_out_CTX_M15_4_sig))
      
      model_svm_CTX_M15_4_sig = svm(CTX_M15_gene ~ Cefepime, data = train, 
                                    probability = T, kernel = "sigmoid", cost = tune_out_CTX_M15_4_sig$best.parameters$cost,
                                    gamma = tune_out_CTX_M15_4_sig$best.parameters$gamma )
      
      cat("Sigmoid SVM based on Simplicity Principle (FEP):\n")
      print(model_svm_CTX_M15_4_sig)
      
      #Prediction on test (SVM)
      test$pred_svm1_CTX_M15_sig <- predict(model_svm1_CTX_M15_sig, test)
      
      # Model evaluation
      confm_svm4_CTX_M15_sig <- table(actual = test$CTX_M15_gene, prediction = test$pred_svm_CTX_M15_4_sig)
      
      metrics_svm5_CTX_M15_sig<- compute_metrics(confm_svm4_CTX_M15_sig)
      
      cat("Performance of SVM-Sigmoid model based on Simplicity Principle (FEP):\n")
      metrics_svm5_CTX_M15_sig
    }
    
    else if (model=="SVM-Sig-Simplicity Principle_CIP"){
      
      set.seed(1234)
      tune_out_CTX_M15_5_sig <- e1071::tune("svm", CTX_M15_gene ~ Ciprofloxacin
                                            , data = train, kernel = "sigmoid", tunecontrol=tune.control(cross=10),
                                            ranges = list(cost = seq(0.1, 5, length = 20)
                                                          , gamma = seq(0.1, 5, length = 20)))
      
      cat("Tuned hyper-parameters throgh 10-Fold Cross Validation:\n")
      print(summary(tune_out_CTX_M15_5_sig))
      
      model_svm_CTX_M15_5_sig = svm(CTX_M15_gene ~ Ciprofloxacin, data = train, 
                           probability = T, kernel = "polynomial", cost = tune_out_CTX_M15_5_sig$best.parameters$cost,
                           gamma = tune_out_CTX_M15_5_sig$best.parameters$gamma )
      
      cat("Sigmoid SVM based on Simplicity Principle (CIP):\n")
      print(model_svm_CTX_M15_5_sig)
      
      #Prediction on test (SVM)----------------------------------------------------
      test$pred_svm_CTX_M15_5_sig <- predict(model_svm_CTX_M15_5_sig, test)
      
      # Model evaluation
      confm_svm5_CTX_M15_sig <- table(actual = test$CTX_M15_gene, prediction = test$pred_svm_CTX_M15_5_sig)
      
      metrics_svm6_CTX_M15_sig<- compute_metrics(confm_svm5_CTX_M15_sig)
      
      cat("Performance of SVM-Sigmoid model based on Simplicity Principle (CIP):\n")
      metrics_svm6_CTX_M15_sig
    }
    else if (model=="SVM-Poly-Chi Squared test"){
      set.seed(1234)
      tune_out_CTX_M15 <- e1071::tune("svm", CTX_M15_gene ~ Amikacin + Cefotaxime + Gentamicin + 
                                        Imipenem + Cefepime + Ceftazidime + Ciprofloxacin + Meropenem + 
                                        Ceftriaxone + Aztreonam 
                                      , data = train, kernel = "polynomial",tunecontrol=tune.control(cross=10),
                                      ranges = list(degree = c(1, 2, 3), cost = seq(0.1, 5, length = 20)
                                                    , gamma = seq(0.5, 5, length = 20)))
      
      cat("Tuned hyper-parameters throgh 10-Fold Cross Validation:\n")
      print(summary(tune_out_CTX_M15))
      
      # Tuned hyper-parameters are placed in the model.
      model_svm_CTX_M15 = svm(CTX_M15_gene ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                              Cefepime + Ceftazidime + Ciprofloxacin +
                              Meropenem + Ceftriaxone + Aztreonam, data = train, probability = T, 
                            kernel = "polynomial", cost = tune_out_CTX_M15$best.parameters$cost, 
                            degree = tune_out_CTX_M15$best.parameters$degree , gamma = tune_out_CTX_M15$best.parameters$gamma )
      
      cat("Polynomial SVM based on Chi-Squared test results:\n")
      print(model_svm_CTX_M15)
      
      #Prediction on test (SVM)----------------------------------------------------
      test$pred_svm_CTX_M15 <- predict(model_svm_CTX_M15, test)
      
      # Model evaluation
      confm_svm_CTX_M15 <- table(actual = test$CTX_M15_gene, prediction = test$pred_svm_CTX_M15)
      
      metrics_svm_CTX_M15<- compute_metrics(confm_svm_CTX_M15)
      
      cat("Performance of SVM-Polynomial model based on Chi-Squared test results:\n")
      metrics_svm_CTX_M15
    }
    else if (model=="SVM-Poly-model agnostic"){
      
      set.seed(1234)
      tune_out1_CTX_M15 <- e1071::tune("svm", CTX_M15_gene ~ Imipenem + Gentamicin + Meropenem + Amikacin 
                                       , data = train, kernel = "polynomial",tunecontrol=tune.control(cross=10),
                                       ranges = list(degree = c(1, 2, 3), cost = seq(0.1, 5, length = 20)
                                                     , gamma = seq(0.1, 5, length = 20)))
      
      cat("Tuned hyper-parameters throgh 10-Fold Cross Validation:\n")
      print(summary(tune_out1_CTX_M15))
      
      model_svm1_CTX_M15 = svm(CTX_M15_gene ~ Imipenem + Gentamicin + Meropenem + Amikacin, data = train, 
                           probability = T, kernel = "polynomial", cost = tune_out1_CTX_M15$best.parameters$cost,
                           degree = tune_out1_CTX_M15$best.parameters$degree , gamma = tune_out1_CTX_M15$best.parameters$gamma )
      
      cat("Polynomial SVM based on model-agnostic approach results:\n")
      print(model_svm1_CTX_M15)
      
      #Prediction on test (SVM)----------------------------------------------------
      test$pred_svm1_CTX_M15 <- predict(model_svm1_CTX_M15, test)
      
      # Model evaluation
      confm_svm1_CTX_M15 <- table(actual = test$CTX_M15_gene, prediction = test$pred_svm1_CTX_M15)
      
      metrics_svm1_CTX_M15<- compute_metrics(confm_svm1_CTX_M15)
      
      cat("Performance of SVM-Polynomial model based on model-agnostic approach results:\n")
      metrics_svm1_CTX_M15
      
    }
    
    else if (model=="SVM-Poly-Wald test"){
      
      set.seed(1234)
      tune_out_CTX_M15_2 <- e1071::tune("svm", CTX_M15_gene ~ Cefotaxime + Cefepime
                                        , data = train, kernel = "polynomial", tunecontrol=tune.control(cross=10),
                                        ranges = list(degree = c(1, 2, 3), cost = seq(0.1, 10, length = 20)
                                                      , gamma = seq(0.1, 5, length = 20)))
      
      cat("Tuned hyper-parameters throgh 10-Fold Cross Validation:\n")
      print(summary(tune_out_CTX_M15_2))
      
      model_svm_CTX_M15_2 = svm(CTX_M15_gene ~ Cefotaxime + Cefepime , data = train, 
                           probability = T, kernel = "polynomial", cost = tune_out_CTX_M15_2$best.parameters$cost,
                           degree = tune_out_CTX_M15_2$best.parameters$degree , gamma = tune_out_CTX_M15_2$best.parameters$gamma )
      
      cat("Polynomial SVM based on the Wald test result\n")
      print(model_svm_CTX_M15_2)
      
      #Prediction on test (SVM)----------------------------------------------------
      test$pred_svm_CTX_M15_2 <- predict(model_svm_CTX_M15_2, test)
      
      # Model evaluation
      confm_svm2_CTX_M15 <- table(actual = test$CTX_M15_gene, prediction = test$pred_svm_CTX_M15_2)
      
      metrics_svm2_CTX_M15<- compute_metrics(confm_svm2_CTX_M15)
      
      cat("Performance of SVM-Polynomial model based on the Wald test results:\n")
      metrics_svm2_CTX_M15
    }
    
    else if (model=="SVM-Poly-Simplicity principle_CTX"){
      
      set.seed(1234)
      tune_out_CTX_M15_3 <- e1071::tune("svm", CTX_M15_gene ~ Cefotaxime
                                        , data = train, kernel = "polynomial", tunecontrol=tune.control(cross=10),
                                        ranges = list(degree = c(1, 2, 3), cost = seq(0.1, 10, length = 20)
                                                      , gamma = seq(0.1, 5, length = 20)))
      
      cat("Tuned hyper-parameters throgh 10-Fold Cross Validation:\n")
      print(summary(tune_out_CTX_M15_3))
      
      model_svm_CTX_M15_3 = svm(CTX_M15_gene ~ Cefotaxime , data = train, 
                           probability = T, kernel = "polynomial", cost = tune_out_CTX_M15_3$best.parameters$cost,
                           degree = tune_out_CTX_M15_3$best.parameters$degree , gamma = tune_out_CTX_M15_3$best.parameters$gamma )
      
      cat("Polynomial SVM based on the Simplicity Principle (CTX):\n")
      print(model_svm_CTX_M15_3)
      
      #Prediction on test (SVM)----------------------------------------------------
      test$pred_svm_CTX_M15_3 <- predict(model_svm_CTX_M15_3, test)
      
      # Model evaluation
      confm_svm3_CTX_M15 <- table(actual = test$CTX_M15_gene, prediction = test$pred_svm_CTX_M15_3)
      
      metrics_svm3_CTX_M15<- compute_metrics(confm_svm3_CTX_M15)
      
      cat("Performance of SVM-Polynomial model based on the Simplicity Principle (CTX):\n")
      metrics_svm3_CTX_M15
    }
    
    else if (model=="SVM-Poly-Simplicity principle_FEP"){
      
      set.seed(1234)
      tune_out_CTX_M15_4 <- e1071::tune("svm", CTX_M15_gene ~ Cefepime
                                        , data = train, kernel = "polynomial", tunecontrol=tune.control(cross=10),
                                        ranges = list(degree = c(1, 2, 3), cost = seq(0.1, 5, length = 20)
                                                      , gamma = seq(0.1, 5, length = 20)))
      
      cat("Tuned hyper-parameters throgh 10-Fold Cross Validation:\n")
      print(summary(tune_out_CTX_M15_4))
      
      model_svm_CTX_M15_4 = svm(CTX_M15_gene ~ Cefepime , data = train, 
                                probability = T, kernel = "polynomial", cost = tune_out_CTX_M15_4$best.parameters$cost,
                                degree = tune_out_CTX_M15_4$best.parameters$degree , gamma = tune_out_CTX_M15_4$best.parameters$gamma )
      
      cat("Polynomial SVM based on the Simplicity Principle (FEP):\n")
      print(model_svm_CTX_M15_4)
      
      #Prediction on test (SVM)----------------------------------------------------
      test$pred_svm_CTX_M15_4 <- predict(model_svm_CTX_M15_4, test)
      
      # Model evaluation
      confm_svm4_CTX_M15 <- table(actual = test$CTX_M15_gene, prediction = test$pred_svm_CTX_M15_4)
      
      metrics_svm4_CTX_M15<- compute_metrics(confm_svm4_CTX_M15)
      
      cat("Performance of SVM-Polynomial model based on the Simplicity Principle (FEP):\n")
      metrics_svm4_CTX_M15
    }
    
    else if (model=="SVM-Poly-Simplicity principle_CIP"){
      
      set.seed(1234)
      tune_out_CTX_M15_5 <- e1071::tune("svm", CTX_M15_gene ~ Ciprofloxacin
                                        , data = train, kernel = "polynomial", tunecontrol=tune.control(cross=10),
                                        ranges = list(degree = c(1, 2, 3), cost = seq(0.1, 5, length = 20)
                                                      , gamma = seq(0.1, 5, length = 20)))
      
      cat("Tuned hyper-parameters throgh 10-Fold Cross Validation:\n")
      print(summary(tune_out_CTX_M15_5))
      
      model_svm_CTX_M15_5 = svm(CTX_M15_gene ~ Ciprofloxacin , data = train, 
                                probability = T, kernel = "polynomial", cost = tune_out_CTX_M15_5$best.parameters$cost,
                                degree = tune_out_CTX_M15_5$best.parameters$degree , gamma = tune_out_CTX_M15_5$best.parameters$gamma )
      
      cat("Polynomial SVM based on the Simplicity Principle (CIP):\n")
      print(model_svm_CTX_M15_5)
      
      #Prediction on test (SVM)----------------------------------------------------
      test$pred_svm_CTX_M15_5 <- predict(model_svm_CTX_M15_5, test)
      
      # Model evaluation
      confm_svm5_CTX_M15 <- table(actual = test$CTX_M15_gene, prediction = test$pred_svm_CTX_M15_5)
      
      metrics_svm5_CTX_M15<- compute_metrics(confm_svm5_CTX_M15)
      
      cat("Performance of SVM-Polynomial model based on the Simplicity Principle (CIP):\n")
      metrics_svm5_CTX_M15
    }
    
    else if (model=="CatBoost-Chi Squared test"){
      # CatBoost Models
      Bio_Data <- as.data.frame(Bio_Data)
      
      #Convert categorical variables to factor
      categorical_var <- c(colnames(Bio_Data))
      
      Bio_Data[, categorical_var[-1]] <- lapply(Bio_Data[, categorical_var[-1]], factor)
      
      levels(Bio_Data$Ampicillin) <- c(levels(Bio_Data$Ampicillin),0)
      
      #Divide Data set into Train and Test
      set.seed(123)
      
      train_cases <- sample(1:nrow(Bio_Data), nrow(Bio_Data) * 0.8)
      
      train <- Bio_Data[train_cases,]
      
      test  <- Bio_Data[- train_cases,]
      
      grid <- expand.grid(
        depth = c(2,3,4,5),
        learning_rate = c(0.1,0.2,0.3),
        l2_leaf_reg = c(2,3,4),
        rsm = 1,
        border_count = 1,
        iterations = seq(10,100,5)
      )
      fitControl <- trainControl(method = "cv",
                                 number = 10,
                                 classProbs = TRUE)
      
      set.seed(1234)
      model_catboost_CTX_M15 <- train(x = train[,c("Amikacin", "Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                   "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")],
                                      y = make.names(train[,"CTX_M15_gene"]),
                                      maximize = TRUE,
                                      method = catboost.caret, metric = "Accuracy", 
                                      tuneGrid =  grid,
                                      trControl = fitControl)
      
      cat("Tuned hyper-parameters throgh 10-Fold Cross Validation:\n")
      print(model_catboost_CTX_M15$bestTune)  
      
      #Prediction on test (CatBoost)----------------------------------------------------
      test$cat_CTX_M15 = predict(model_catboost_CTX_M15, test)
      
      # Model evaluation
      cat_test_CTX_M15 <- table(actual = test$CTX_M15_gene, prediction = test$cat_CTX_M15)
      
      metrics_cat_CTX_M15<- compute_metrics(cat_test_CTX_M15)
      
      cat("Performance of CatBoost model based on Chi-Squared test results:\n")
      metrics_cat_CTX_M15
    }
    
    else if (model=="CatBoost-model agnostic"){
      # CatBoost Models--------------------------------------------------------------------
      Bio_Data <- as.data.frame(Bio_Data)
      
      #Convert categorical variables to factor----------------------------------
      categorical_var <- c(colnames(Bio_Data))
      
      Bio_Data[, categorical_var[-1]] <- lapply(Bio_Data[, categorical_var[-1]], factor)
      
      levels(Bio_Data$Ampicillin) <- c(levels(Bio_Data$Ampicillin),0)
      
      #Divide Data set into Train and Test---------------------------------------
      set.seed(123)
      
      train_cases <- sample(1:nrow(Bio_Data), nrow(Bio_Data) * 0.8)
      
      train <- Bio_Data[train_cases,]
      
      test  <- Bio_Data[- train_cases,]
      
      grid2_CTX_M15 <- expand.grid(
        depth = c(2,3,4),
        learning_rate = c(0.1,0.2,0.3),
        l2_leaf_reg = c(1,2,3),
        rsm = 1,
        border_count = 1,
        iterations = seq(10,50,1)
      )
      fitControl2_CTX_M15 <- trainControl(method = "cv",
                                  number = 10,
                                  classProbs = TRUE)
      
      set.seed(1234)
      model_catboost2_CTX_M15 <- train(x = train[,c("Cefotaxime",  "Cefepime", "Ciprofloxacin")],
                                       y = make.names(train[,"CTX_M15_gene"]),
                                       maximize = TRUE,
                                       method = catboost.caret, metric = "Accuracy", 
                                       tuneGrid =  grid2_CTX_M15,
                                       trControl = fitControl2_CTX_M15)
      
      cat("Tuned hyper-parameters throgh 10-Fold Cross Validation:\n")
      print(model_catboost2_CTX_M15$bestTune)  
      
      test$cat_CTX_M15 = predict(model_catboost2_CTX_M15, test)
      
      # Model Evaluation
      confm_cat_CTX_M15 <- table(actual = test$CTX_M15_gene, prediction = test$cat_CTX_M15)
      
      metrics_cat1_CTX_M15<- compute_metrics(confm_cat_CTX_M15)
      
      cat("Performance of CatBoost model based on model-agnostic approach results:\n")
      metrics_cat1_CTX_M15
    }
    
    else if (model=="CatBoost-Wald test"){
      # CatBoost Models--------------------------------------------------------------------
      Bio_Data <- as.data.frame(Bio_Data)
      
      #Convert categorical variables to factor----------------------------------
      categorical_var <- c(colnames(Bio_Data))
      
      Bio_Data[, categorical_var[-1]] <- lapply(Bio_Data[, categorical_var[-1]], factor)
      
      levels(Bio_Data$Ampicillin) <- c(levels(Bio_Data$Ampicillin),0)
      
      #Divide Data set into Train and Test---------------------------------------
      set.seed(123)
      
      train_cases <- sample(1:nrow(Bio_Data), nrow(Bio_Data) * 0.8)
      
      train <- Bio_Data[train_cases,]
      
      test  <- Bio_Data[- train_cases,]
      
      grid3_CTX_M15 <- expand.grid(
        depth = c(2,3,4),
        learning_rate = c(0.1,0.2,0.3),
        l2_leaf_reg = c(1,2,3,4),
        rsm = 1,
        border_count = 1,
        iterations = seq(10,20,1)
      )
      fitControl3_CTX_M15 <- trainControl(method = "cv",
                                  number = 10,
                                  classProbs = TRUE)
      
      model_catboost3_CTX_M15 <- train(x = train[,c("Cefotaxime",  "Cefepime")],
                                       y = make.names(train[,"CTX_M15_gene"]),
                                       maximize = TRUE,
                                       method = catboost.caret, metric = "Accuracy", 
                                       tuneGrid =  grid3_CTX_M15,
                                       trControl = fitControl3_CTX_M15)
      
      cat("Tuned hyper-parameters throgh 10-Fold Cross Validation:\n")
      print(model_catboost3_CTX_M15$bestTune)  
      
      test$cat2_CTX_M15 = predict(model_catboost3_CTX_M15, test)
      # Model evaluation
      confm_cat2_CTX_M15 <- table(actual = test$CTX_M15_gene, prediction = test$cat2_CTX_M15)
      
      metrics_cat2_CTX_M15<- compute_metrics(confm_cat2_CTX_M15)
      
      cat("Performance of CatBoost model based on the Wald test results:\n")
      metrics_cat2_CTX_M15
    }
    
    else if (model=="CatBoost-Simplicity principle_CTX"){
      # CatBoost Models--------------------------------------------------------------------
      Bio_Data <- as.data.frame(Bio_Data)
      
      #Convert categorical variables to factor----------------------------------
      categorical_var <- c(colnames(Bio_Data))
      
      Bio_Data[, categorical_var[-1]] <- lapply(Bio_Data[, categorical_var[-1]], factor)
      
      levels(Bio_Data$Ampicillin) <- c(levels(Bio_Data$Ampicillin),0)
      
      #Divide Data set into Train and Test---------------------------------------
      set.seed(123)
      
      train_cases <- sample(1:nrow(Bio_Data), nrow(Bio_Data) * 0.8)
      
      train <- Bio_Data[train_cases,]
      
      test  <- Bio_Data[- train_cases,]
      
      grid4_CTX_M15 <- expand.grid(
        depth = c(2,3),
        learning_rate = c(0.1,0.2,0.3),
        l2_leaf_reg = c(1,2,3),
        rsm = 1,
        border_count = 1,
        iterations = seq(10,30,1)
      )
      fitControl4_CTX_M15 <- trainControl(method = "cv",
                                          number = 10,
                                          classProbs = TRUE)
      model_catboost4_CTX_M15 <- train(x = train["Cefotaxime"],
                                       y = make.names(train[["CTX_M15_gene"]]),
                                       maximize = TRUE,
                                       method = catboost.caret, metric = "Accuracy", 
                                       tuneGrid =  grid4_CTX_M15,
                                       trControl = fitControl4_CTX_M15)
      
      cat("Tuned hyper-parameters throgh 10-Fold Cross Validation:\n")
      print(model_catboost4_CTX_M15$bestTune)  
      
      test$cat3_CTX_M15 = predict(model_catboost4_CTX_M15, test)
      # Model evaluation
      confm_cat3_CTX_M15 <- table(actual = test$CTX_M15_gene, prediction = test$cat3_CTX_M15)
      
      metrics_cat3_CTX_M15<- compute_metrics(confm_cat3_CTX_M15)
      
      cat("Performance of CatBoost model based on the simplicity principle (CTX):\n")
      metrics_cat3_CTX_M15
    }
  
  
  else if (model=="CatBoost-Simplicity principle_FEP"){
    # CatBoost Models--------------------------------------------------------------------
    Bio_Data <- as.data.frame(Bio_Data)
    
    #Convert categorical variables to factor----------------------------------
    categorical_var <- c(colnames(Bio_Data))
    
    Bio_Data[, categorical_var[-1]] <- lapply(Bio_Data[, categorical_var[-1]], factor)
    
    levels(Bio_Data$Ampicillin) <- c(levels(Bio_Data$Ampicillin),0)
    
    #Divide Data set into Train and Test---------------------------------------
    set.seed(123)
    
    train_cases <- sample(1:nrow(Bio_Data), nrow(Bio_Data) * 0.8)
    
    train <- Bio_Data[train_cases,]
    
    test  <- Bio_Data[- train_cases,]
    
    grid5_CTX_M15 <- expand.grid(
      depth = c(2,3,4),
      learning_rate = c(0.1,0.2,0.3),
      l2_leaf_reg = c(1,2,3),
      rsm = 1,
      border_count = 1,
      iterations = seq(10,500,50)
    )
    fitControl5_CTX_M15 <- trainControl(method = "cv",
                                        number = 10,
                                        classProbs = TRUE)
    model_catboost5_CTX_M15 <- train(x = train["Cefepime"],
                                     y = make.names(train[["CTX_M15_gene"]]),
                                     maximize = TRUE,
                                     method = catboost.caret, metric = "Accuracy", 
                                     tuneGrid =  grid5_CTX_M15,
                                     trControl = fitControl5_CTX_M15)
    
    cat("Tuned hyper-parameters throgh 10-Fold Cross Validation:\n")
    print(model_catboost5_CTX_M15$bestTune)  
    
    test$cat4_CTX_M15 = predict(model_catboost5_CTX_M15, test)
    # Model evaluation
    confm_cat4_CTX_M15 <- table(actual = test$CTX_M15_gene, prediction = test$cat4_CTX_M15)
    
    metrics_cat4_CTX_M15<- compute_metrics(confm_cat4_CTX_M15)
    
    cat("Performance of CatBoost model based on the simplicity principle (FEP):\n")
    metrics_cat4_CTX_M15
  }


  else if (model=="CatBoost-Simplicity principle_CIP"){
  # CatBoost Models--------------------------------------------------------------------
  Bio_Data <- as.data.frame(Bio_Data)
  
  #Convert categorical variables to factor----------------------------------
  categorical_var <- c(colnames(Bio_Data))
  
  Bio_Data[, categorical_var[-1]] <- lapply(Bio_Data[, categorical_var[-1]], factor)
  
  levels(Bio_Data$Ampicillin) <- c(levels(Bio_Data$Ampicillin),0)
  
  #Divide Data set into Train and Test---------------------------------------
  set.seed(123)
  
  train_cases <- sample(1:nrow(Bio_Data), nrow(Bio_Data) * 0.8)
  
  train <- Bio_Data[train_cases,]
  
  test  <- Bio_Data[- train_cases,]
  
  grid6_CTX_M15 <- expand.grid(
    depth = c(2,3,4),
    learning_rate = c(0.1,0.2,0.3),
    l2_leaf_reg = c(1,2,3),
    rsm = 1,
    border_count = 1,
    iterations = seq(10,200,40)
  )
  fitControl6_CTX_M15 <- trainControl(method = "cv",
                                      number = 10,
                                      classProbs = TRUE)
  model_catboost6_CTX_M15 <- train(x = train["Ciprofloxacin"],
                                   y = make.names(train[["CTX_M15_gene"]]),
                                   maximize = TRUE,
                                   method = catboost.caret, metric = "Accuracy", 
                                   tuneGrid =  grid5_CTX_M15,
                                   trControl = fitControl5_CTX_M15)
  
  cat("Tuned hyper-parameters throgh 10-Fold Cross Validation:\n")
  print(model_catboost6_CTX_M15$bestTune)  
  
  test$cat5_CTX_M15 = predict(model_catboost6_CTX_M15, test)
  # Model evaluation
  confm_cat5_CTX_M15 <- table(actual = test$CTX_M15_gene, prediction = test$cat5_CTX_M15)
  
  metrics_cat5_CTX_M15<- compute_metrics(confm_cat5_CTX_M15)
  
  cat("Performance of CatBoost model based on the simplicity principle (CIP):\n")
  metrics_cat5_CTX_M15
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
        result <- chisq.test(table(data[[variable]], data$CTX_M15_gene))
        cat("p-value for", variable, ":", result$p.value, "\n")
      }
      
      # Perform chi-square tests for multiple variables
      variables <- c("Amikacin", "Cefotaxime", "Gentamicin", "Imipenem", "Cefepime",
                     "Ceftazidime", "Ciprofloxacin", "Meropenem", "Ampicillin",
                     "Ceftriaxone", "Aztreonam")
      
      for (variable in variables) {
        perform_chi_square(Bio_Data, variable)
      } 
    }
  }
   else if (method == "Wald test"){
     model_Log1_CTX_M15 <- glm(CTX_M15_gene ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                                 Cefepime + Ceftazidime + Ciprofloxacin +
                                 Meropenem + Ceftriaxone + Aztreonam 
                               , family = "binomial", data = train)
    
    cat("P-values show Wald test results:\n")
    print(summary(model_Log1_CTX_M15))
  }
   else if (method == "model agnostic-LR"){
    # Feature Selection by DALEX Package--------------------------------------
    # the 'explained_glm_TEM' object is created using the 'explain' function from the DALEX package. 
    # The variable_importance function is used to calculate the permutation-based 
    # feature importance ('fi_glm_TEM') with 50 permutations. Finally, the feature importance is displayed 
    # using 'fi_glm_TEM' and plotted using' plot(fi_glm_TEM)'.
     model_Log1_CTX_M15 <- glm(CTX_M15_gene ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                                 Cefepime + Ceftazidime + Ciprofloxacin +
                                 Meropenem + Ceftriaxone + Aztreonam 
                               , family = "binomial", data = train)
    
     explained_glm_CTX_M15 <- explain(model = model_Log1_CTX_M15, data=train[,c(2:9,11:12)], variables = colnames(train[,c(2:9,11:12)]),
                                      y=as.vector(as.numeric(train$CTX_M15_gene))-1, label = "LR", 
                                      type = "classification")
     #50 permuatation
     fi_glm_CTX_M15 = variable_importance(explained_glm_CTX_M15, B=50 ,variables = colnames(train[,c(2:9,11:12)]),loss_function = loss_root_mean_square, type = "raw")
    
    cat("Root mean squared dropout loss due to a feature elimination through model-agnostic approach in LR:\n")
    print(fi_glm_CTX_M15)
  }
  else if (method == "model agnostic-NBC"){
    train_control <- trainControl(
      method = "cv", 
      number = 10
    )
    
    search_grid_CTX_M15 <- expand.grid(
      usekernel = c(FALSE,TRUE),
      fL = seq(0,5,0.5),
      adjust = seq(0,5,0.5)
    )
    
    
    model_nbdalex_CTX_M15 <- train(CTX_M15_gene ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                                     Cefepime + Ceftazidime + Ciprofloxacin +
                                     Meropenem + Ceftriaxone + Aztreonam,
                                   data = train,
                                   method = "nb",
                                   metric = "Accuracy",
                                   trControl = train_control,
                                   tuneGrid = search_grid_CTX_M15)
    # Feature selection through model-agnostic approach (DALEX)
    # the 'explain' function from the DALEX package is used to explain the Naive Bayes classifier model 
    # 'model_nbdalex_TEM'. The relevant variables are selected and provided as data along with their 
    # corresponding column names. The target variable 'CTX_M15_gene' is transformed to numeric values 
    # (y = as.vector(as.numeric(train$CTX_M15_gene)) - 1). The label is set to "NBC" for Naive Bayes classifier,
    # and the type is set to "classification".
    expnb_CTX_M15 = explain(model = model_nbdalex_CTX_M15, data=train[,c("Amikacin", "Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                         "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")], variables = colnames(train[,c("Amikacin", "Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                                                                                                                                 "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")]),
                            y=as.vector(as.numeric(train$CTX_M15_gene))-1, label = "NBC", 
                            type = "classification")
    
    # The 'variable_importance' function is then used to calculate the variable importance based on the
    # explained model. The 'B' parameter is set to 50 for the number of permutations, and the loss function
    # is set to "loss_root_mean_square". The resulting variable importance is stored in 'fitnb_TEM'
    fitnb_CTX_M15 = variable_importance(expnb_CTX_M15, B=50 ,variables = colnames(train[,c("Amikacin", "Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                                           "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")]), loss_function = loss_root_mean_square, type = "raw" )
    
    cat("Root mean squared dropout loss due to a feature elimination through model-agnostic approach in NBC:\n")
    print(fitnb_CTX_M15)
  }
  else if (method == "model agnostic-LDA"){
    model_ldadalex_CTX_M15 <- train(CTX_M15_gene ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                                      Cefepime + Ceftazidime + Ciprofloxacin +
                                      Meropenem + Ceftriaxone + Aztreonam ,
                                    data = train,
                                    method = "lda",
                                    metric = "Accuracy")
    
    explda_CTX_M15 = explain(model = model_ldadalex_CTX_M15, data=train[,c("Amikacin", "Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                           "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")], variables = colnames(train[,c("Amikacin", "Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                                                                                                                                   "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")]),
                             y=as.vector(as.numeric(train$CTX_M15_gene))-1, label = "LDA", 
                             type = "classification")
    
    fitlda_CTX_M15 = variable_importance(explda_CTX_M15, B=50 ,variables = colnames(train[,c("Amikacin", "Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                                             "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")]), loss_function = loss_root_mean_square, type = "raw" )
    cat("Root mean squared dropout loss due to a feature elimination through model-agnostic approach in LDA:\n")
    print(fitlda_CTX_M15)
  }
  else if (method == "model agnostic-Sig SVM"){
    set.seed(1234)
    tune_out_CTX_M15_sig <- e1071::tune("svm", CTX_M15_gene ~ Amikacin + Cefotaxime + Gentamicin + 
                                          Imipenem + Cefepime + Ceftazidime + Ciprofloxacin + Meropenem + 
                                          Ceftriaxone + Aztreonam 
                                        , data = train, kernel = "sigmoid",tunecontrol=tune.control(cross=10),
                                        ranges = list(cost = seq(0.1, 5, length = 20)
                                                      , gamma = seq(0.5, 5, length = 20)))
    
    model_svm_CTX_M15_sig = svm(CTX_M15_gene ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                                  Cefepime + Ceftazidime + Ciprofloxacin +
                                  Meropenem + Ceftriaxone + Aztreonam, data = train, 
                                probability = T, kernel = "sigmoid", cost = tune_out_CTX_M15_sig$best.parameters$cost,
                                gamma = tune_out_CTX_M15_sig$best.parameters$gamma)
    
    explained_SVM_CTX_M15 <- explain(model = model_svm_CTX_M15, data=train[,c(2:9,11:12)], variables = colnames(train[,c(2:9,11:12)]),
                                     y=as.vector(as.numeric(train$CTX_M15_gene))-1, label = "SVM-Polynomial",
                                     type = "classification")
    
    #50 permuatation
    fi_SVM_CTX_M15 = variable_importance(explained_SVM_CTX_M15, B=50,variables = colnames(train[,c(2:9,11:12)]), loss_function = loss_root_mean_square, type = "raw" )
    
    cat("Root mean squared dropout loss due to a feature elimination through model-agnostic approach in sigmoid kernel of SVM:\n")
    print(fi_SVM_CTX_M15)
  }
  else if (method == "model agnostic-Poly SVM"){
    # it is highly recommeded that first train the model, then specify following hyper-parameters according to the tuned hyper-parameters.
    # Following values of the hyper-parameters have been specified in this way.
    set.seed(1234)
    tune_out_CTX_M15 <- e1071::tune("svm", CTX_M15_gene ~ Amikacin + Cefotaxime + Gentamicin + 
                                      Imipenem + Cefepime + Ceftazidime + Ciprofloxacin + Meropenem + 
                                      Ceftriaxone + Aztreonam 
                                    , data = train, kernel = "polynomial",tunecontrol=tune.control(cross=10),
                                    ranges = list(degree = c(1, 2, 3), cost = seq(0.1, 5, length = 20)
                                                  , gamma = seq(0.5, 5, length = 20)))
    
    model_svm_CTX_M15 = svm(CTX_M15_gene ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                              Cefepime + Ceftazidime + Ciprofloxacin +
                              Meropenem + Ceftriaxone + Aztreonam, data = train, probability = T, 
                            kernel = "polynomial", cost = tune_out_CTX_M15$best.parameters$cost, 
                            degree = tune_out_CTX_M15$best.parameters$degree , gamma = tune_out_CTX_M15$best.parameters$gamma)
    
    explained_SVM_CTX_M15 <- explain(model = model_svm_CTX_M15, data=train[,c(2:9,11:12)], variables = colnames(train[,c(2:9,11:12)]),
                                     y=as.vector(as.numeric(train$CTX_M15_gene))-1, label = "SVM-Polynomial",
                                     type = "classification")
    
    #50 permuatation
    fi_SVM_CTX_M15 = variable_importance(explained_SVM_CTX_M15, B=50,variables = colnames(train[,c(2:9,11:12)]), loss_function = loss_root_mean_square, type = "raw" )
   
    cat("Root mean squared dropout loss due to a feature elimination through model-agnostic approach in polynomial kernel of SVM:\n")
    print(fi_SVM_CTX_M15)
  }
  else if (method == "model agnostic-CatBoost"){
    Bio_Data <- as.data.frame(Bio_Data)
    
    #Convert categorical variables to factor----------------------------------
    categorical_var <- c(colnames(Bio_Data))
    
    Bio_Data[, categorical_var[-1]] <- lapply(Bio_Data[, categorical_var[-1]], factor)
    
    levels(Bio_Data$Ampicillin) <- c(levels(Bio_Data$Ampicillin),0)
    
    #Divide Data set into Train and Test---------------------------------------
    set.seed(123)
    
    train_cases <- sample(1:nrow(Bio_Data), nrow(Bio_Data) * 0.8)
    
    train <- Bio_Data[train_cases,]
    
    test  <- Bio_Data[- train_cases,]
    
    grid <- expand.grid(
      depth = c(2,3,4,5),
      learning_rate = c(0.1,0.2,0.3),
      l2_leaf_reg = c(2,3,4),
      rsm = 1,
      border_count = 1,
      iterations = seq(10,100,5)
    )
    fitControl <- trainControl(method = "cv",
                               number = 10,
                               classProbs = TRUE)
    
    set.seed(1234)
    model_catboost_CTX_M15 <- train(x = train[,c("Amikacin", "Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                 "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")],
                                    y = make.names(train[,"CTX_M15_gene"]),
                                    maximize = TRUE,
                                    method = catboost.caret, metric = "Accuracy", 
                                    tuneGrid =  grid,
                                    trControl = fitControl)
    
    expsda_CTX_M15 = explain(model = model_catboost_CTX_M15, data=train[,c("Amikacin", "Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                           "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")], variables = colnames(train[,c("Amikacin", "Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                                                                                                                                   "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")]),
                             y=as.vector(as.numeric(train$CTX_M15_gene))-1, label = "CatBoost", 
                             type = "classification")
    
    fitcat_CTX_M15 = variable_importance(expsda_CTX_M15, B=50 ,variables = colnames(train[,c("Amikacin", "Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                                             "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")]), loss_function = loss_root_mean_square, type = "raw" )
    
    cat("Root mean squared dropout loss due to a feature elimination through model-agnostic approach in CatBoost:\n")
    print(fitcat_CTX_M15)
  }
  else {
    cat("Invalid task. Please specify one of the following tasks: feature_selection, model_evaluation\n")
    
  } 
}
if (is.null(task)){
  cat("You should specify a task\n")
}

