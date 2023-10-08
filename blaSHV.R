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

# Calculate the proportion of SHV_gene
proportion_SHV_gene <- prop.table(table(Bio_Data$SHV_gene))
cat("Proportion of SHV_gene:\n")
print(proportion_SHV_gene)


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
      model_Log1_SHV <- glm(SHV_gene ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                              Cefepime + Ceftazidime + Ciprofloxacin +
                              Meropenem + Ceftriaxone + Aztreonam,
                            family = "binomial", data = train)
      cat("Summary of LR model based on Chi Squared test including Null deviance, Residual deviance, and AIC, etc.:\n")
      print(summary(model_Log1_SHV))
      
      test$probs1 <- predict(model_Log1_SHV, test, type = "response")
      test$pred_logreg1 <- ifelse(test$probs1 >= 0.5, 1, 0)
      # Model 1 evaluation
      confm_logreg1 <- table(actual = test$SHV_gene, prediction = test$pred_logreg1)
      metrics_logreg1 <- compute_metrics(confm_logreg1)
      cat("Performance of LR model based on Chi-Squared test results:\n")
      metrics_logreg1
    }
    else if (model == "LR-model agnostic/Wald test"){
      model_Log2 <- glm(SHV_gene ~ Cefotaxime
                        , family = "binomial", data = train)
      cat("Summary of LR model based on model agnostic approach / Wald test including Null deviance, Residual deviance, and AIC, etc.:\n")
      print(summary(model_Log2))
      
      ##Prediction on test (Model 2)--------------------------------------------------
      test$probs2 <- predict(model_Log2, test, type = "response")
      test$pred_logreg2 <- ifelse(test$probs2 >= 0.5, 1, 0)
      
      # Model 2 evaluation
      confm_logreg2 <- table(actual = test$SHV_gene, prediction = test$pred_logreg2)
      metrics_logreg2 <- compute_metrics(confm_logreg2)
      cat("Performance of LR model based on model agnostic approach / Wald test results:\n")
      metrics_logreg2
    }
    else if (model == "NBC-Chi Squared test"){
      train_control <- trainControl(
        method = "cv", 
        number = 10
      )
      
      search_grid <- expand.grid(
        usekernel = c(FALSE,TRUE),
        fL = seq(0,5,0.5),
        adjust = seq(0,5,0.5)
      )
      
      model_nbdalex_SHV <- train(SHV_gene ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                                   Cefepime + Ceftazidime + Ciprofloxacin +
                                   Meropenem + Ceftriaxone + Aztreonam,
                                 data = train,
                                 method = "nb",
                                 metric = "Accuracy",
                                 trControl = train_control,
                                 tuneGrid = search_grid,
      )
      cat("Tuned hyper-parameters through 10-fold Cross-Validation:\n")
      print(model_nbdalex_SHV$bestTune)
      
      ##Prediction on test (Naive Bayes Model)----------------------------------
      test$pred_nb <- predict(model_nbdalex_SHV, test)
      
      # Model evaluation
      confm_nb <- table(actual = test$SHV_gene, prediction = test$pred_nb)
      metrics_nb<- compute_metrics(confm_nb)
      cat("Performance of NBC model based on Chi-Squared test results:\n")
      metrics_nb
    }
    else if (model == "NBC-model agnostic/Wald test"){
      train_control <- trainControl(
        method = "cv",
        number = 10
      )
      
      search_grid <- expand.grid(
        usekernel = c(FALSE,TRUE),
        fL = seq(0,5,0.5),
        adjust = seq(0,5,0.5)
      )
      
      model_nbdalex1_SHV <- train(SHV_gene ~ Cefotaxime,
                                  data = train,
                                  method = "nb",
                                  metric = "Accuracy",
                                  trControl = train_control,
                                  tuneGrid = search_grid)
      
      cat("Tuned hyper-parameters through 10-fold Cross-Validation:\n")
      print(model_nbdalex1_SHV$bestTune)
      
      ##Prediction on test (Naive Bayes Model)----------------------------------
      test$pred_nb1 <- predict(model_nbdalex1_SHV, test)
      
      # Model evaluation
      confm_nb1 <- table(actual = test$SHV_gene, prediction = test$pred_nb1)
      metrics_nb1<- compute_metrics(confm_nb1)
      cat("Performance of NBC model based on model-agnostic approach/Wald test results:\n")
      metrics_nb1
    }
    else if (model == "LDA-Chi Squared test"){
      #Linear Discriminant Analysis-------------------------------------------------
      model_lda_dalex_SHV <- train(SHV_gene ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                                     Cefepime + Ceftazidime + Ciprofloxacin +
                                     Meropenem + Ceftriaxone + Aztreonam,
                                   data = train,
                                   method = "lda",
                                   metric = "Accuracy")
      
      cat("LDA model based on Chi-Squared test results:\n")
      print(model_lda_dalex_SHV)
      
      #Prediction on test (LDA Model)-----------------------------------------------
      test$pred_lda <- predict(model_lda_dalex_SHV, test)
      
      # Model Evaluation-------------------------------------------------
      confm_lda <- table(actual = test$SHV_gene, prediction = test$pred_lda)
      metrics_lda<- compute_metrics(confm_lda)
      cat("Performance of LDA model based on Chi-Squared test results:\n")
      metrics_lda
    }
    else if (model == "LDA-model agnostic/Wald test"){
      model_lda_dalex1_SHV <- train(SHV_gene ~ Cefotaxime,
                                    data = train,
                                    method = "lda",
                                    metric = "Accuracy")
      
      cat("LDA model based on model-agnostic approach/Wald test results:\n")
      print(model_lda_dalex1_SHV)
      
      #Prediction on test (LDA Model)-----------------------------------------------
      test$pred_lda1 <- predict(model_lda_dalex1_SHV, test)  
      
      # Model evaluation
      confm_lda1 <- table(actual = test$SHV_gene, prediction = test$pred_lda1)
      metrics_lda1<- compute_metrics(confm_lda1)
      cat("Performance of LDA model based on model-agnostic approach/Wald test results:\n")
      metrics_lda1
    }
    else if (model == "SVM-Sig-Chi Squared test"){
      #Fit Support Vector Machine model to train data set
      #This code performs a grid search for tuning the hyper-parameters of the SVM model using the tune function
      # from the e1071 package. It then fits the SVM model with the specified hyper-parameters using the svm function.
      # Model 1
      set.seed(123)
      tune_out <- e1071::tune("svm", SHV_gene ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                                Cefepime + Ceftazidime + Ciprofloxacin +
                                Meropenem + Ceftriaxone + Aztreonam
                              , data = train, kernel = "sigmoid", tunecontrol=tune.control(cross=10),
                              ranges = list(cost = seq(0.001, 10, length = 20)
                                            , gamma = seq(0, 5, length = 20)))
      
      cat("Tuned hyper-parameters through 10-Fold Cross Validation:\n")
      print(summary(tune_out))
      
      model_svm_SHV = svm(SHV_gene ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                            Cefepime + Ceftazidime + Ciprofloxacin +
                            Meropenem + Ceftriaxone + Aztreonam, data = train, 
                          probability = T, kernel = "sigmoid", cost = tune_out$best.parameters$cost ,
                          gamma = tune_out$best.parameters$gamma)
      
      cat("Sigmoid SVM based on Chi-Squared test results:\n")
      print(model_svm_SHV)
      
      #Prediction on test (SVM)----------------------------------------------------
      test$pred_svm_SHV <- predict(model_svm_SHV, test)
      
      # Model evaluation
      confm_svm_SHV <- table(actual = test$SHV_gene, prediction = test$pred_svm_SHV)
      
      metrics_svm1<- compute_metrics(confm_svm_SHV)
      
      cat("Performance of SVM-Sigmoid model based on Chi-Squared test results:\n")
      metrics_svm1
    }
    else if (model=="SVM-Sig-model agnostic"){
      # Cefepime, Aztreonam and Ceftazidime has been selected based on DALEX Package for the 2nd model-------------
      set.seed(1234)
      tune_out2 <- e1071::tune("svm", SHV_gene ~ Cefepime + Aztreonam + Ceftazidime
                               , data = train, kernel = "sigmoid", tunecontrol=tune.control(cross=10),
                               ranges = list(cost = seq(0.001, 5, length = 20)
                                             , gamma = seq(0, 5, length = 20)))
      cat("Tuned hyper-parameters through 10-Fold Cross Validation:\n")
      print(summary(tune_out2))
      
      model_svm_SHV2 = svm(SHV_gene ~ Cefepime + Aztreonam + Ceftazidime, data = train, 
                           probability = T, kernel = "sigmoid", cost =tune_out2$best.parameters$cost  ,
                           gamma = tune_out2$best.parameters$gamma )
      
      cat("Sigmoid SVM based on model-agnostic approach results:\n")
      print(model_svm_SHV2)
      
      #Prediction on test (SVM)----------------------------------------------------
      test$pred_svm1 <- predict(model_svm_SHV2, test)
      
      # Model evaluation
      confm_svm1 <- table(actual = test$SHV_gene, prediction = test$pred_svm1)
      
      metrics_svm2<- compute_metrics(confm_svm1)
      
      cat("Performance of SVM-Sigmoid model based on model-agnostic approach results:\n")
      metrics_svm2
    }
    else if (model=="SVM-Sig-Wald test"){
      set.seed(1234)
      tune_out3 <- e1071::tune("svm", SHV_gene ~ Cefotaxime
                               , data = train, kernel = "sigmoid", tunecontrol=tune.control(cross=10),
                               ranges = list(cost = seq(0.001, 5, length = 20)
                                             , gamma = seq(0, 5, length = 20)))
      
      cat("Tuned hyper-parameters through 10-Fold Cross Validation:\n")
      print(summary(tune_out3))
      
      model_svm_SHV3 = svm(SHV_gene ~ Cefotaxime, data = train, 
                           probability = T, kernel = "sigmoid", cost = tune_out3$best.parameters$cost ,
                           gamma = tune_out3$best.parameters$gamma )
      
      cat("Sigmoid SVM based on Wald test results:\n")
      print(model_svm_SHV3)
      
      #Prediction on test (SVM)----------------------------------------------------
      test$pred_svm2 <- predict(model_svm_SHV3, test)
      
      # Model evaluation
      confm_svm2 <- table(actual = test$SHV_gene, prediction = test$pred_svm2)
      
      metrics_svm3<- compute_metrics(confm_svm2)
      
      cat("Performance of SVM-Sigmoid model based on Wald test results:\n")
      metrics_svm3
      
    } 
    
    else if (model=="SVM-Poly-Chi Squared test"){
      set.seed(96875)
      tune_out4 <- e1071::tune("svm", SHV_gene ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                                 Cefepime + Ceftazidime + Ciprofloxacin +
                                 Meropenem + Ceftriaxone + Aztreonam  
                               , data = train, kernel = "polynomial", tunecontrol=tune.control(cross=10),
                               ranges = list(degree = c(1, 2, 3), cost = seq(0.1, 10, length = 20)
                                             , gamma = seq(0.1, 5, length = 20)))
      
      cat("Tuned hyper-parameters through 10-Fold Cross Validation:\n")
      print(summary(tune_out4))
      
      model_svm_SHV4 = svm(SHV_gene ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                             Cefepime + Ceftazidime + Ciprofloxacin +
                             Meropenem + Ceftriaxone + Aztreonam, data = train, 
                           probability = T, kernel = "polynomial", cost = tune_out4$best.parameters$cost,
                           degree = tune_out4$best.parameters$degree , gamma = tune_out4$best.parameters$gamma )
      
      cat("Polynomial SVM based on Chi-Squared test results:\n")
      print(model_svm_SHV4)
      
      #Prediction on test (SVM)----------------------------------------------------
      test$pred_svm3 <- predict(model_svm_SHV4, test)
      
      # Model evaluation
      confm_svm3 <- table(actual = test$SHV_gene, prediction = test$pred_svm3)
      
      metrics_svm4<- compute_metrics(confm_svm3)
      
      cat("Performance of SVM-Polynomial model based on Chi-Squared test results:\n")
      metrics_svm4
    }
    else if (model=="SVM-Poly-model agnostic"){
      set.seed(1234)
      tune_out5 <- e1071::tune("svm", SHV_gene ~ Gentamicin + Amikacin + Ceftriaxone + Cefepime + Ciprofloxacin
                               , data = train, kernel = "polynomial", tunecontrol=tune.control(cross=10),
                               ranges = list(cost = seq(0.1, 5, length = 20), degree = c(1, 2, 3)
                                             , gamma = seq(0.1, 5, length = 20)))
      
      cat("Tuned hyper-parameters through 10-Fold Cross Validation:\n")
      print(summary(tune_out5))
      
      # Tuned hyper-parameters are placed in the model.
      model_svm_SHV5 = svm(SHV_gene ~ Gentamicin + Amikacin + Ceftriaxone + Cefepime + Ciprofloxacin,
                           data = train, probability = T, kernel = "polynomial", cost = tune_out5$best.parameters$cost, 
                           degree = tune_out5$best.parameters$degree , gamma = tune_out5$best.parameters$gamma )
      
      cat("Polynomial SVM based on model agnostic approach results:\n")
      print(model_svm_SHV5)
      
      #Prediction on test (SVM)----------------------------------------------------
      test$pred_svm4 <- predict(model_svm_SHV5, test)
      
      # Model evaluation
      confm_svm4 <- table(actual = test$SHV_gene, prediction = test$pred_svm4)
      
      metrics_svm5<- compute_metrics(confm_svm4)
      
      cat("Performance of SVM-Polynomial model based on model-agnostic approach results:\n")
      metrics_svm5
    }
    else if (model=="SVM-Poly-Wald test"){
      set.seed(1234)
      
      tune_out6 <- e1071::tune("svm", SHV_gene ~ Cefotaxime
                               , data = train, kernel = "polynomial", tunecontrol=tune.control(cross=10),
                               ranges = list(cost = seq(0.001, 5, length = 20), degree = c(1, 2, 3)
                                             , gamma = seq(0, 5, length = 20)))
      
      cat("Tuned hyper-parameters through 10-Fold Cross Validation:\n")
      print(summary(tune_out6))
      
      model_svm_SHV6 = svm(SHV_gene ~ Cefotaxime, data = train, 
                           probability = T, kernel = "polynomial", cost = tune_out6$best.parameters$cost, 
                           degree = tune_out6$best.parameters$degree , gamma = tune_out6$best.parameters$gamma )
      
      cat("Polynomial SVM based on Wald test results:\n")
      print(model_svm_SHV6)
      
      #Prediction on test (SVM)----------------------------------------------------
      test$pred_svm5 <- predict(model_svm_SHV6, test)
      
      # Model evaluation
      confm_svm5 <- table(actual = test$SHV_gene, prediction = test$pred_svm5)
      
      metrics_svm6<- compute_metrics(confm_svm5)
      
      cat("Performance of SVM-Polynomial model based on Wald test results:\n")
      metrics_svm6
      
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
        depth = c(1,2,3),
        learning_rate = c(0.1,0.2,0.3),
        l2_leaf_reg = c(1,2,3),
        rsm = 1,
        border_count = 1,
        iterations = seq(10,150,10)
      )
      
      fitControl <- trainControl(method = "cv",
                                 number = 10,
                                 classProbs = TRUE)
      
      set.seed(1234)
      
      model_catboost_SHV <- train(x = train[,c("Amikacin", "Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                               "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")],
                                  y = make.names(train[,"SHV_gene"]),
                                  maximize = TRUE,
                                  method = catboost.caret, metric = "Accuracy",
                                  tuneGrid =  grid,
                                  trControl = fitControl)
      
      cat("Tuned hyper-parameters through 10-Fold Cross Validation:\n")
      print(model_catboost_SHV$bestTune)  
      
      #Prediction on test (CatBoost)----------------------------------------------------
      test$cat_SHV = predict(model_catboost_SHV, test)
      
      # Model evaluation
      confm_cat_SHV <- table(actual = test$SHV_gene, prediction = test$cat_SHV)
      
      metrics_cat<- compute_metrics(confm_cat_SHV)
      
      cat("Performance of CatBoost model based on Chi-Squared test results:\n")
      metrics_cat
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
      
      grid2 <- expand.grid(
        depth = c(2,3),
        learning_rate = c(0.1,0.2),
        l2_leaf_reg = c(1,2),
        rsm = 1,
        border_count = 1,
        iterations = seq(10,100,10)
      )
      
      fitControl2 <- trainControl(method = "cv",
                                  number = 10,
                                  classProbs = TRUE)
      
      set.seed(1234)
      
      model_catboost2_SHV <- train(x = train[,c("Cefotaxime",  "Aztreonam", "Imipenem")],
                                   y = make.names(train[,"SHV_gene"]),
                                   maximize = TRUE,
                                   method = catboost.caret, metric = "Accuracy", 
                                   tuneGrid =  grid2,
                                   trControl = fitControl2)
      
      cat("Tuned hyper-parameters through 10-Fold Cross Validation:\n")
      print(model_catboost2_SHV$bestTune)  
      
      test$cat_SHV1 = predict(model_catboost2_SHV, test)
      
      # Model Evaluation
      confm_cat_SHV1 <- table(actual = test$SHV_gene, prediction = test$cat_SHV1)
      
      metrics_cat1<- compute_metrics(confm_cat_SHV1)
      
      cat("Performance of CatBoost model based on model-agnostic approach results:\n")
      metrics_cat1
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
      
      grid3 <- expand.grid(
        depth = c(2,3),
        learning_rate = c(0.1,0.2),
        l2_leaf_reg = c(1,2),
        rsm = 1,
        border_count = 1,
        iterations = seq(10,100,10)
      )
      
      fitControl3 <- trainControl(method = "cv",
                                  number = 10,
                                  classProbs = TRUE)
      
      model_catboost3_SHV <- train(x = train["Cefotaxime"],
                                   y = make.names(train[["SHV_gene"]]),
                                   maximize = TRUE,
                                   method = catboost.caret, metric = "Accuracy", 
                                   tuneGrid =  grid3,
                                   trControl = fitControl3)
      
      cat("Tuned hyper-parameters through 10-Fold Cross Validation:\n")
      print(model_catboost3_SHV$bestTune)  
      
      test$cat_SHV2 = predict(model_catboost3_SHV, test)
      # Model evaluation
      confm_cat2_SHV <- table(actual = test$SHV_gene, prediction = test$cat_SHV2)  
      
      metrics_cat2<- compute_metrics(confm_cat2_SHV)
      
      cat("Performance of CatBoost model based on Wald test results:\n")
      metrics_cat2
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
        result <- chisq.test(table(data[[variable]], data$SHV_gene))
        cat("p-value for", variable, ":", result$p.value, "\n")
      }
      
      # Perform chi-square tests for multiple variables
      variables <- c("Amikacin", "Cefotaxime", "Gentamicin", "Imipenem", "Cefepime",
                     "Ceftazidime", "Ciprofloxacin", "Meropenem", "Ampicillin",
                     "Ceftriaxone", "Aztreonam")
      
      for (variable in variables) {
        perform_chi_square(Bio_Data, variable)
      } 
      # p-value for Ampicillin is not available due to resistance of all the isolates to this antibiotic
    }
  }
  else if (method == "Wald test"){
    model_Log1_SHV <- glm(SHV_gene ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                            Cefepime + Ceftazidime + Ciprofloxacin +
                            Meropenem + Ceftriaxone + Aztreonam,
                          family = "binomial", data = train)
    
    cat("P-values show Wald test results:\n")
    print(summary(model_Log1_SHV))  # p-values show result of wald test for each antibiotic
  }
  else if (method == "model agnostic-LR"){
    # Feature Selection by DALEX Package--------------------------------------
    # the 'explained_glm_SHV' object is created using the 'explain' function from the DALEX package. 
    # The variable_importance function is used to calculate the permutation-based 
    # feature importance ('fi_glm_SHV') with 50 permutations. Finally, the feature importance is displayed 
    # using 'fi_glm_SHV' and plotted using' plot(fi_glm_SHV)'.
    model_Log1_SHV <- glm(SHV_gene ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                            Cefepime + Ceftazidime + Ciprofloxacin +
                            Meropenem + Ceftriaxone + Aztreonam,
                          family = "binomial", data = train)
    
    explained_glm_SHV <- explain(model = model_Log1_SHV, data=train[,c(2:9,11:12)], variables = colnames(train[,c(2:9,11:12)]),
                                 y=as.vector(as.numeric(train$SHV_gene))-1, label = "LR", 
                                 type = "classification")
    #50 permuatations
    fi_glm_SHV = variable_importance(explained_glm_SHV, B=50,variables = colnames(train[,c(2:9,11:12)]),loss_function = loss_root_mean_square, type = "raw")
    
    cat("Root mean squared dropout loss due to a feature elimination:\n")
    print(fi_glm_SHV)
  }
  else if (method == "model agnostic-NBC"){
    train_control <- trainControl(
      method = "cv", 
      number = 10
    )
    
    search_grid <- expand.grid(
      usekernel = c(FALSE,TRUE),
      fL = seq(0,5,0.5),
      adjust = seq(0,5,0.5)
    )
    
    model_nbdalex_SHV <- train(SHV_gene ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                                 Cefepime + Ceftazidime + Ciprofloxacin +
                                 Meropenem + Ceftriaxone + Aztreonam,
                               data = train,
                               method = "nb",
                               metric = "Accuracy",
                               trControl = train_control,
                               tuneGrid = search_grid,
    )
    # Feature selection through model-agnostic approach (DALEX)
    # the 'explain' function from the DALEX package is used to explain the Naive Bayes classifier model 
    # 'model_nbdalex_SHV'. The relevant variables are selected and provided as data along with their 
    # corresponding column names. The target variable 'SHV_gene' is transformed to numeric values 
    # (y = as.vector(as.numeric(train$SHV_gene)) - 1). The label is set to "NBC" for Naive Bayes classifier,
    # and the type is set to "classification".
    expnb_SHV = explain(model = model_nbdalex_SHV, data=train[,c("Amikacin", "Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                 "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")], variables = colnames(train[,c("Amikacin", "Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                                                                                                                         "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")]),
                        y=as.vector(as.numeric(train$SHV_gene))-1, label = "NBC", 
                        type = "classification")
    
    # The 'variable_importance' function is then used to calculate the variable importance based on the
    # explained model. The 'B' parameter is set to 50 for the number of permutations, and the loss function
    # is set to "loss_root_mean_square". The resulting variable importance is stored in 'fitnb_SHV'
    fitnb_SHV = variable_importance(expnb_SHV, B=50 ,variables = colnames(train[,c("Amikacin", "Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                                   "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")]), loss_function = loss_root_mean_square, type = "raw" )
    
    cat("Root mean squared dropout loss due to a feature elimination:\n")
    print(fitnb_SHV)
  }
  else if (method == "model agnostic-LDA"){
    model_lda_dalex_SHV <- train(SHV_gene ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                                   Cefepime + Ceftazidime + Ciprofloxacin +
                                   Meropenem + Ceftriaxone + Aztreonam,
                                 data = train,
                                 method = "lda",
                                 metric = "Accuracy")
    
    explda_SHV = explain(model = model_lda_dalex_SHV, data=train[,c("Amikacin", "Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                    "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")], variables = colnames(train[,c("Amikacin", "Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                                                                                                                            "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")]),
                         y=as.vector(as.numeric(train$SHV_gene))-1, label = "LDA", 
                         type = "classification")
    
    fitlda_SHV = variable_importance(explda_SHV, B=50 ,variables = colnames(train[,c("Amikacin", "Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                                     "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")]), loss_function = loss_root_mean_square, type = "raw" )
    cat("Root mean squared dropout loss due to a feature elimination:\n")
    print(fitlda_SHV)
  }
  else if (method == "model agnostic-Sig SVM"){
    set.seed(123)
    tune_out <- e1071::tune("svm", SHV_gene ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                              Cefepime + Ceftazidime + Ciprofloxacin +
                              Meropenem + Ceftriaxone + Aztreonam
                            , data = train, kernel = "sigmoid", tunecontrol=tune.control(cross=10),
                            ranges = list(cost = seq(0.001, 10, length = 20)
                                          , gamma = seq(0, 5, length = 20)))
    
    cat("Tuned hyper-parameters through 10-Fold Cross Validation:\n")
    print(summary(tune_out))
    
    model_svm_SHV = svm(SHV_gene ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                          Cefepime + Ceftazidime + Ciprofloxacin +
                          Meropenem + Ceftriaxone + Aztreonam, data = train, 
                        probability = T, kernel = "sigmoid", cost = tune_out$best.parameters$cost ,
                        gamma = tune_out$best.parameters$gamma)
    
    explained_SVM_SHV <- explain(model = model_svm_SHV, data=train[,c(2:9,11:12)], variables = colnames(train[,c(2:9,11:12)]),
                                 y=as.vector(as.numeric(train$SHV_gene))-1, label = "SVM-Sigmoid",
                                 type = "classification")
    
    #50 permuatation
    fi_SVM_SHV = variable_importance(explained_SVM_SHV, B=50,variables = colnames(train[,c(2:9,11:12)]), loss_function = loss_root_mean_square, type = "raw" )
    
    cat("Root mean squared dropout loss due to a feature elimination:\n")
    print(fi_SVM_SHV)
  }
  else if (method == "model agnostic-Poly SVM"){
    # it is highly recommeded that first train the model, then specify following hyper-parameters according to the tuned hyper-parameters.
    # Following values of the hyper-parameters have been specified in this way.
    set.seed(96875)
    tune_out4 <- e1071::tune("svm", SHV_gene ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                               Cefepime + Ceftazidime + Ciprofloxacin +
                               Meropenem + Ceftriaxone + Aztreonam  
                             , data = train, kernel = "polynomial", tunecontrol=tune.control(cross=10),
                             ranges = list(degree = c(1, 2, 3), cost = seq(0.1, 10, length = 20)
                                           , gamma = seq(0.1, 5, length = 20)))
    
    model_svm_SHV4 = svm(SHV_gene ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                           Cefepime + Ceftazidime + Ciprofloxacin +
                           Meropenem + Ceftriaxone + Aztreonam, data = train, 
                         probability = T, kernel = "polynomial", cost = tune_out4$best.parameters$cost,
                         degree = tune_out4$best.parameters$degree , gamma = tune_out4$best.parameters$gamma)
    
    explained_SVM_SHV2 <- explain(model = model_svm_SHV4, data=train[,c(2:9,11:12)], variables = colnames(train[,c(2:9,11:12)]),
                                  y=as.vector(as.numeric(train$SHV_gene))-1, label = "SVM-Polynomial",
                                  type = "classification")
    
    #50 permuatation
    fi_SVM_SHV2 = variable_importance(explained_SVM_SHV2, B=50,variables = colnames(train[,c(2:9,11:12)]), loss_function = loss_root_mean_square, type = "raw" )
    
    print(fi_SVM_SHV2)
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
      depth = c(1,2,3),
      learning_rate = c(0.1,0.2,0.3),
      l2_leaf_reg = c(1,2,3),
      rsm = 1,
      border_count = 1,
      iterations = seq(10,150,10)# seq(40,100,1)
    )
    
    fitControl <- trainControl(method = "cv",
                               number = 10,
                               classProbs = TRUE)
    
    set.seed(1234)
    
    model_catboost_SHV <- train(x = train[,c("Amikacin", "Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                             "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")],
                                y = make.names(train[,"SHV_gene"]),
                                maximize = TRUE,
                                method = catboost.caret, metric = "Accuracy",
                                tuneGrid =  grid,
                                trControl = fitControl)
    expsda_SHV = explain(model = model_catboost_SHV, data=train[,c("Amikacin", "Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                   "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")], variables = colnames(train[,c("Amikacin", "Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                                                                                                                           "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")]),
                         y=as.vector(as.numeric(train$SHV_gene))-1, label = "CatBoost", 
                         type = "classification")
    
    fitcat_SHV = variable_importance(expsda_SHV, B=50 ,variables = colnames(train[,c("Amikacin", "Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                                     "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")]), loss_function = loss_root_mean_square, type = "raw" )
    print(fitcat_SHV)
  }
  else {
    cat("Invalid task. Please specify one of the following tasks: feature_selection, model_evaluation\n")
    
  } 
}
if (is.null(task)){
  cat("You should specify a task\n")
}
