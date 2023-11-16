source('F1-score.R') # Call F1-Score function for computing point estimation and its 95% confidence interval

source('compute_metrics.R') # Call compute_metrics function for computing all classification metrics and their confidence interval.

source('Packages_installation.R') # Installation of all required packages


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
verbose = args$verbose
method = args$method

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

# Calculate the proportion of MDR
proportion_MDR <- prop.table(table(Bio_Data$MDR))
cat("Proportion of MDR:\n")
print(proportion_MDR)


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
    
    if (model == "LR-Chi Squared test/Wald test"){
      # Logistic regression (Model 1)--------------------------------------------------
      model_Log1_MDR <- glm(MDR ~ Imipenem + Gentamicin + Ciprofloxacin + Meropenem + Aztreonam +
                              Amikacin + Cefotaxime + Cefepime + Ceftazidime + Ceftriaxone
                            , family = "binomial", data = train)
      cat("Summary of LR model based on Chi Squared test including Null deviance, Residual deviance, and AIC, etc.:\n")
      print(summary(model_Log1_MDR))
      
      test$probs1_MDR <- predict(model_Log1_MDR, test, type = "response")
      
      test$pred_logreg1_MDR <- ifelse(test$probs1_MDR >= 0.5, 1, 0)
      # Model 1 evaluation
      confm_logreg1_MDR <- table(actual = test$MDR, prediction = test$pred_logreg1_MDR)
      metrics_logreg1_MDR <- compute_metrics(confm_logreg1_MDR)
      cat("Performance of LR model based on Chi-Squared test results:\n")
      metrics_logreg1_MDR
    }
    else if (model == "LR-model agnostic"){
      model_Log2_MDR <- glm(MDR ~ Gentamicin + Ciprofloxacin + Amikacin + Meropenem + Cefotaxime
                            , family = "binomial", data = train)
      cat("Summary of LR model based on model agnostic approach including Null deviance, Residual deviance, and AIC, etc.:\n")
      print(summary(model_Log2_MDR))
      
      ##Prediction on test (Model 2)--------------------------------------------------
      test$probs2_MDR <- predict(model_Log2_MDR, test, type = "response")
      
      test$pred_logreg2_MDR <- ifelse(test$probs2_MDR >= 0.5, 1, 0)
      
      # Model 2 evaluation
      confm_logreg2_MDR <- table(actual = test$MDR, prediction = test$pred_logreg2_MDR)
      metrics_logreg2_MDR <- compute_metrics(confm_logreg2_MDR)
      cat("Performance of LR model based on model agnostic approach results:\n")
      metrics_logreg2_MDR
    }
    
    else if (model == "LR-Simplicity Principle_GN"){
      model_Log3_MDR <- glm(MDR ~ Gentamicin 
                            , family = "binomial", data = train)
      cat("Summary of LR model based on the simplicity principle (GN) including Null deviance, Residual deviance, and AIC, etc.:\n")
      print(summary(model_Log3_MDR))
      
      ##Prediction on test (Model 2)--------------------------------------------------
      test$probs3_MDR <- predict(model_Log3_MDR, test, type = "response")
      
      test$pred_logreg3_MDR <- ifelse(test$probs3_MDR >= 0.5, 1, 0)
      
      # Model 2 evaluation
      confm_logreg3_MDR <- table(actual = test$MDR, prediction = test$pred_logreg3_MDR)
      
      metrics_logreg3_MDR <- compute_metrics(confm_logreg3_MDR)
      
      cat("Performance of LR model based on the simplicity principle (GN) results:\n")
      
      metrics_logreg3_MDR
    }
    
    else if (model == "NBC-Chi Squared test/Wald test"){
      train_control <- trainControl(
        method = "cv", 
        number = 10
      )
      
      search_grid_MDR <- expand.grid(
        usekernel = c(FALSE,TRUE),
        fL = seq(0,5,0.5),
        adjust = seq(0,5,0.5)
      )
      
      model_nbdalex_MDR <- train(MDR ~ Imipenem + Gentamicin + Ciprofloxacin + Meropenem + Aztreonam +
                                   Amikacin + Cefotaxime + Cefepime + Ceftazidime + Ceftriaxone,
                                 data = train,
                                 method = "nb",
                                 metric = "Accuracy",
                                 trControl = train_control,
                                 tuneGrid = search_grid_MDR,
      )
      
      cat("Tuned hyper-parameters through 10-fold Cross-Validation:\n")
      print(model_nbdalex_MDR$bestTune)
      
      ##Prediction on test (Naive Bayes Model)----------------------------------
      test$pred_nb_MDR <- predict(model_nbdalex_MDR, test)
      
      # Model evaluation
      confm_nb_MDR <- table(actual = test$MDR, prediction = test$pred_nb_MDR)
      
      metrics_nb_MDR<- compute_metrics(confm_nb_MDR)
      
      cat("Performance of NBC model based on Chi-Squared test and Wald test results:\n")
      metrics_nb_MDR
    }
    else if (model == "NBC-model agnostic"){
      train_control <- trainControl(
        method = "cv",
        number = 10
      )
      
      search_grid_MDR <- expand.grid(
        usekernel = c(FALSE,TRUE),
        fL = seq(0,5,0.5),
        adjust = seq(0,5,0.5)
      )
      
      model_nb_MDR <- train(MDR ~ Gentamicin,
                            data = train,
                            method = "nb",
                            metric = "Accuracy",
                            trControl = train_control,
                            tuneGrid = search_grid_MDR,
      )
      
      cat("Tuned hyper-parameters through 10-fold Cross-Validation:\n")
      print(model_nb_MDR$bestTune)
      
      ##Prediction on test (Naive Bayes Model)----------------------------------
      test$pred_nb1_MDR <- predict(model_nb_MDR, test)
      
      # Model evaluation
      confm_nb1_MDR <- table(actual = test$MDR, prediction = test$pred_nb1_MDR)
      
      metrics_nb1_MDR<- compute_metrics(confm_nb1_MDR)
      cat("Performance of NBC model based on model-agnostic approach results:\n")
      metrics_nb1_MDR
    }
    else if (model == "LDA-Chi Squared test/Wald test"){
      #Linear Discriminant Analysis-------------------------------------------------
      model_ldadalex_MDR <- train(MDR ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                                    Cefepime + Ceftazidime + Ciprofloxacin +
                                    Meropenem + Ceftriaxone + Aztreonam,
                                  data = train,
                                  method = "lda",
                                  metric = "Accuracy")
      
      cat("LDA model based on Chi-Squared test results:\n")
      print(model_ldadalex_MDR)
      
      #Prediction on test (LDA Model)-----------------------------------------------
      test$pred_lda_MDR <- predict(model_ldadalex_MDR, test)
      
      # Model Evaluation-------------------------------------------------
      confm_lda_MDR <- table(actual = test$MDR, prediction = test$pred_lda_MDR)
      
      metrics_lda_MDR<- compute_metrics(confm_lda_MDR)
      cat("Performance of LDA model based on Chi-Squared test and Wald test results:\n")
      metrics_lda_MDR
    }
    else if (model == "LDA-model agnostic"){
      model_lda_MDR <- train(MDR ~ Gentamicin,
                             data = train,
                             method = "lda",
                             metric = "Accuracy")
      
      cat("LDA model based on model-agnostic approach results:\n")
      print(model_lda_MDR)
      
      #Prediction on test (LDA Model)-----------------------------------------------
      test$pred_lda1_MDR <- predict(model_lda_MDR, test) 
      
      # Model evaluation
      confm_lda1_MDR <- table(actual = test$MDR, prediction = test$pred_lda1_MDR)
      
      metrics_lda1_MDR<- compute_metrics(confm_lda1_MDR)
      
      cat("Performance of LDA model based on model-agnostic approach results:\n")
      metrics_lda1_MDR
    }
    else if (model == "SVM-Sig-Chi Squared test/Wald test"){
      #Fit Support Vector Machine model to train data set
      #This code performs a grid search for tuning the hyper-parameters of the SVM model using the tune function
      # from the e1071 package. It then fits the SVM model with the specified hyper-parameters using the svm function.
      # Model 1
      set.seed(1234)
      tune_out_MDR <- e1071::tune("svm", MDR ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                                    Cefepime + Ceftazidime + Ciprofloxacin +
                                    Meropenem + Ceftriaxone + Aztreonam
                                  , data = train, kernel = "sigmoid", tunecontrol=tune.control(cross=10),
                                  ranges = list(cost = seq(0.1, 10, length = 20)
                                                , gamma = seq(0.1, 5, length = 20)))
      
      cat("Tuned hyper-parameters through 10-Fold Cross Validation:\n")
      print(summary(tune_out_MDR))
      
      model_svm_MDR = svm(MDR ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                            Cefepime + Ceftazidime + Ciprofloxacin +
                            Meropenem + Ceftriaxone + Aztreonam, data = train, 
                          probability = T, kernel = "sigmoid", cost = tune_out_MDR$best.parameters$cost ,
                          gamma = tune_out_MDR$best.parameters$gamma)
      
      cat("Sigmoid SVM based on Chi-Squared test and Wald test results:\n")
      print(model_svm_MDR)
      
      #Prediction on test (SVM)----------------------------------------------------
      test$pred_svm_MDR <- predict(model_svm_MDR, test)
      
      # Model evaluation
      confm_svm_MDR <- table(actual = test$MDR, prediction = test$pred_svm_MDR)
      
      metrics_svm1_MDR<- compute_metrics(confm_svm_MDR)
      
      cat("Performance of SVM-Sigmoid model based on Chi-Squared test and Wald test results:\n")
      metrics_svm1_MDR
    }
    else if (model=="SVM-Sig-model agnostic"){
      # Cefepime, Aztreonam and Ceftazidime has been selected based on DALEX Package for the 2nd model-------------
      set.seed(1234)
      tune_out2_MDR <- e1071::tune("svm", MDR ~ Cefepime + Cefotaxime + Aztreonam + Imipenem + Ceftriaxone + Ceftazidime
                                   , data = train, kernel = "sigmoid", tunecontrol=tune.control(cross=10),
                                   ranges = list(cost = seq(0.1, 5, length = 20)
                                                 , gamma = seq(0.1, 5, length = 20)))
      
      cat("Tuned hyper-parameters through 10-Fold Cross Validation:\n")
      print(summary(tune_out2_MDR))
      
      model_svm_MDR2 = svm(MDR ~ Cefepime + Cefotaxime + Aztreonam + Imipenem + Ceftriaxone + Ceftazidime,
                           data = train, probability = T, kernel = "sigmoid", cost =tune_out2_MDR$best.parameters$cost  ,
                           gamma = tune_out2_MDR$best.parameters$gamma )
      
      cat("Sigmoid SVM based on model-agnostic approach results:\n")
      print(model_svm_MDR2)
      
      #Prediction on test (SVM)----------------------------------------------------
      test$pred_svm1_MDR <- predict(model_svm_MDR2, test)
      
      # Model evaluation
      confm_svm1_MDR <- table(actual = test$MDR, prediction = test$pred_svm1_MDR)
      
      metrics_svm2_MDR<- compute_metrics(confm_svm1_MDR)
      
      cat("Performance of SVM-Sigmoid model based on model-agnostic approach results:\n")
      metrics_svm2_MDR
    }
    else if (model=="SVM-Sig-Simplicity Principle_GN"){
      set.seed(1234)
      tune_out4_MDR <- e1071::tune("svm", MDR ~ Gentamicin
                                   , data = train, kernel = "sigmoid", tunecontrol=tune.control(cross=10),
                                   ranges = list(cost = seq(0.1, 5, length = 20)
                                                 , gamma = seq(0.1, 5, length = 20)))
      
      cat("Tuned hyper-parameters through 10-Fold Cross Validation:\n")
      print(summary(tune_out4_MDR))
      
      model_svm4_MDR = svm(MDR ~ Gentamicin, data = train, 
                           probability = T, kernel = "sigmoid", cost = tune_out4_MDR$best.parameters$cost ,
                           gamma = tune_out4_MDR$best.parameters$gamma )
      
      cat("Sigmoid SVM based on the simplicity principle (GN) results:\n")
      print(model_svm4_MDR)
      
      #Prediction on test (SVM)----------------------------------------------------
      test$pred_svm3_MDR <- predict(model_svm4_MDR, test)
      
      # Model evaluation
      confm_svm3_MDR <- table(actual = test$MDR, prediction = test$pred_svm3_MDR)
      
      metrics_svm3_MDR<- compute_metrics(confm_svm3_MDR)
      
      cat("Performance of SVM-Sigmoid model based on the simplicity principle (GN) results:\n")
      metrics_svm3_MDR
      
    } 
    
    else if (model=="SVM-Poly-Chi Squared test/Wald test"){
      set.seed(1234)
      tune_out5_MDR <- e1071::tune("svm", MDR ~ Imipenem + Gentamicin + Ciprofloxacin + Meropenem + 
                                     Aztreonam + Amikacin + Cefotaxime + Cefepime + Ceftazidime + 
                                     Ceftriaxone  
                                   , data = train, kernel = "polynomial", tunecontrol=tune.control(cross=10),
                                   ranges = list(degree = c(1, 2, 3), cost = seq(0.1, 10, length = 20)
                                                 , gamma = seq(0.1, 5, length = 20)))
      
      cat("Tuned hyper-parameters through 10-Fold Cross Validation:\n")
      print(summary(tune_out5_MDR))
      
      model_svm5_MDR = svm(MDR ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                             Cefepime + Ceftazidime + Ciprofloxacin +
                             Meropenem + Ceftriaxone + Aztreonam, data = train, 
                           probability = T, kernel = "polynomial", cost = tune_out5_MDR$best.parameters$cost,
                           degree = tune_out5_MDR$best.parameters$degree , gamma = tune_out5_MDR$best.parameters$gamma )
      
      cat("Polynomial SVM based on Chi-Squared test and Wald test results:\n")
      print(model_svm5_MDR)
      
      #Prediction on test (SVM)----------------------------------------------------
      test$pred_svm4_MDR <- predict(model_svm5_MDR, test)
      
      # Model evaluation
      confm_svm4_MDR <- table(actual = test$MDR, prediction = test$pred_svm4_MDR)
      
      metrics_svm4_MDR<- compute_metrics(confm_svm4_MDR)
      
      cat("Performance of SVM-Polynomial model based on Chi-Squared test and Wald test results:\n")
      metrics_svm4_MDR
    }
    else if (model=="SVM-Poly-model agnostic"){
      set.seed(1234)
      tune_out4_MDR <- e1071::tune("svm", MDR ~ Imipenem + Aztreonam + Meropenem + Ceftriaxone + Ceftazidime
                                   , data = train, kernel = "polynomial", tunecontrol=tune.control(cross=10),
                                   ranges = list(cost = seq(0.1, 5, length = 20), degree = c(1, 2, 3)
                                                 , gamma = seq(0.1, 5, length = 20)))
      
      cat("Tuned hyper-parameters through 10-Fold Cross Validation:\n")
      print(summary(tune_out4_MDR))
      
      # Tuned hyper-parameters are placed in the model.
      model_svm51_MDR = svm(MDR ~ Imipenem + Aztreonam + Meropenem + Ceftriaxone + Ceftazidime,
                           data = train, probability = T, kernel = "polynomial", cost = tune_out4_MDR$best.parameters$cost, 
                           degree = tune_out4_MDR$best.parameters$degree , gamma = tune_out4_MDR$best.parameters$gamma )
      
      cat("Polynomial SVM based on model agnostic approach results:\n")
      print(model_svm51_MDR)
      
      #Prediction on test (SVM)----------------------------------------------------
      test$pred_svm41_MDR <- predict(model_svm51_MDR, test)
      
      # Model evaluation
      confm_svm41_MDR <- table(actual = test$MDR, prediction = test$pred_svm41_MDR)
      
      metrics_svm5_MDR<- compute_metrics(confm_svm41_MDR)
      
      cat("Performance of SVM-Polynomial model based on model-agnostic approach results:\n")
      metrics_svm5_MDR
    }
    else if (model=="SVM-Poly-Simplicity Principle_GN"){
      set.seed(1234)
      
      tune_out6_MDR <- e1071::tune("svm", MDR ~ Gentamicin
                                   , data = train, kernel = "polynomial", tunecontrol=tune.control(cross=10),
                                   ranges = list(cost = seq(0.1, 5, length = 20), degree = c(1, 2, 3)
                                                 , gamma = seq(0.1, 5, length = 20)))
      
      cat("Tuned hyper-parameters through 10-Fold Cross Validation:\n")
      print(summary(tune_out6_MDR))
      
      model_svm7_MDR = svm(MDR ~ Gentamicin, data = train, 
                           probability = T, kernel = "polynomial", cost = tune_out6_MDR$best.parameters$cost, 
                           degree = tune_out6_MDR$best.parameters$degree , gamma = tune_out6_MDR$best.parameters$gamma )
      
      cat("Polynomial SVM based on the Simplicity Principle (GN) results:\n")
      print(model_svm7_MDR)
      
      #Prediction on test (SVM)----------------------------------------------------
      test$pred_svm6_MDR <- predict(model_svm7_MDR, test)
      
      # Model evaluation
      confm_svm6_MDR <- table(actual = test$MDR, prediction = test$pred_svm6_MDR)
      
      metrics_svm6_MDR<- compute_metrics(confm_svm6_MDR)
      
      cat("Performance of SVM-Polynomial model based on the Simplicity Principle (GN) results:\n")
      metrics_svm6_MDR
      
    }
    else if (model=="CatBoost-Chi Squared test/Wald test"){
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
        depth = c(2,3,4),
        learning_rate = c(0.1,0.2,0.3),
        l2_leaf_reg = c(1,2,3),
        rsm = 1,
        border_count = 1,
        iterations = seq(20,40,1)
      )
      fitControl <- trainControl(method = "cv",
                                 number = 10,
                                 classProbs = TRUE)
      set.seed(1234)
      model_catboost_MDR <- train(x = train[,c("Amikacin", "Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                               "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")],
                                  y = make.names(train[,"MDR"]),
                                  maximize = TRUE,
                                  method = catboost.caret, metric = "Accuracy", 
                                  tuneGrid =  grid,
                                  trControl = fitControl)
      
      cat("Tuned hyper-parameters through 10-Fold Cross Validation:\n")
      print(model_catboost_MDR$bestTune)  
      
      #Prediction on test (CatBoost)----------------------------------------------------
      test$cat_MDR = predict(model_catboost_MDR, test)
      
      # Model evaluation
      cat_test_MDR <- table(actual = test$MDR, prediction = test$cat_MDR)
      
      metrics_cat_MDR<- compute_metrics(cat_test_MDR)
      
      cat("Performance of CatBoost model based on Chi-Squared test and Wald test results:\n")
      metrics_cat_MDR
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
        depth = c(2,3,4),
        learning_rate = c(0.1,0.2,0.3),
        l2_leaf_reg = c(1,2,3),
        rsm = 1,
        border_count = 1,
        iterations = seq(40,60,1)
      )
      fitControl2 <- trainControl(method = "cv",
                                  number = 10,
                                  classProbs = TRUE)
      set.seed(1234)
      model_catboost2_MDR <- train(x = train[,c("Ciprofloxacin","Gentamicin","Amikacin")],
                                   y = make.names(train[,"MDR"]),
                                   maximize = TRUE,
                                   method = catboost.caret, metric = "Accuracy", 
                                   tuneGrid =  grid2,
                                   trControl = fitControl2)
      
      cat("Tuned hyper-parameters through 10-Fold Cross Validation:\n")
      print(model_catboost2_MDR$bestTune)  
      
      test$cat2_MDR = predict(model_catboost2_MDR, test)
      
      # Model Evaluation
      confm_cat_MDR <- table(actual = test$MDR, prediction = test$cat2_MDR)
      
      metrics_cat1_MDR<- compute_metrics(confm_cat_MDR)
      
      cat("Performance of CatBoost model based on model-agnostic approach results:\n")
      metrics_cat1_MDR
    }
    
    else if (model=="CatBoost-Simplicity Principle_GN"){
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
      
      grid4 <- expand.grid(
        depth = c(2,3,4),
        learning_rate = c(0.1,0.2,0.3),
        l2_leaf_reg = c(1,2,3,4),
        rsm = 1,
        border_count = 1,
        iterations = seq(10,20,1)
      )
      fitControl4 <- trainControl(method = "cv",
                                  number = 5,
                                  classProbs = TRUE)
      model_catboost4_MDR <- train(x = train["Gentamicin"],
                                   y = make.names(train[["MDR"]]),
                                   maximize = TRUE,
                                   method = catboost.caret, metric = "Accuracy", 
                                   tuneGrid =  grid4,
                                   trControl = fitControl4)
      
      cat("Tuned hyper-parameters through 10-Fold Cross Validation:\n")
      print(model_catboost4_MDR$bestTune)  
      
      test$cat4_MDR = predict(model_catboost4_MDR, test)
      # Model evaluation
      confm_cat4_MDR <- table(actual = test$MDR, prediction = test$cat4_MDR)  
      
      metrics_cat2_MDR<- compute_metrics(confm_cat4_MDR)
      
      cat("Performance of CatBoost model based on the Simplicity Principle (GN) results:\n")
      metrics_cat2_MDR
    }
  }
  else if (task == 'feature_selectio') {

        if (method == "Chi-Squared test"){
      # Function to perform chi-square test and print p-value
      perform_chi_square <- function(data, variable) {
        result <- chisq.test(table(data[[variable]], data$MDR))
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
    model_Log1_MDR <- glm(MDR ~ Imipenem + Gentamicin + Ciprofloxacin + Meropenem + Aztreonam +
                            Amikacin + Cefotaxime + Cefepime + Ceftazidime + Ceftriaxone
                          , family = "binomial", data = train)
    
    cat("P-values show Wald test results:\n")
    print(summary(model_Log1_MDR))  # p-values show result of wald test for each antibiotic
  }
  else if (method == "model agnostic-LR"){

    model_Log1_MDR <- glm(MDR ~ Imipenem + Gentamicin + Ciprofloxacin + Meropenem + Aztreonam +
                            Amikacin + Cefotaxime + Cefepime + Ceftazidime + Ceftriaxone
                          , family = "binomial", data = train)
    
    explained_glm_MDR <- explain(model = model_Log1_MDR, data=train[,c(2:9,11:12)], variables = colnames(train[,c(2:9,11:12)]),
                                 y=as.vector(as.numeric(train$MDR))-1, label = "LR", 
                                 type = "classification")
    #50 permuatations
    fi_glm_MDR = variable_importance(explained_glm_MDR, B=50 ,variables = colnames(train[,c(2:9,11:12)]),loss_function = loss_root_mean_square, type = "raw")
    
    cat("Root mean squared dropout loss due to a feature elimination:\n")
    print(fi_glm_MDR)
  }
  else if (method == "model agnostic-NBC"){
    train_control <- trainControl(
      method = "cv", 
      number = 10
    )
    
    search_grid_MDR <- expand.grid(
      usekernel = c(FALSE,TRUE),
      fL = seq(0,5,0.5),
      adjust = seq(0,5,0.5)
    )
    
    model_nbdalex_MDR <- train(MDR ~ Imipenem + Gentamicin + Ciprofloxacin + Meropenem + Aztreonam +
                                 Amikacin + Cefotaxime + Cefepime + Ceftazidime + Ceftriaxone,
                               data = train,
                               method = "nb",
                               metric = "Accuracy",
                               trControl = train_control,
                               tuneGrid = search_grid_MDR,
    )
    
    expnb_MDR = explain(model = model_nbdalex_MDR, data=train[,c("Amikacin", "Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                 "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")], variables = colnames(train[,c("Amikacin", "Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                                                                                                                         "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")]),
                        y=as.vector(as.numeric(train$MDR))-1, label = "NBC", 
                        type = "classification")
    
    fitnb_MDR = variable_importance(expnb_MDR, B=50 ,variables = colnames(train[,c("Amikacin", "Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                                   "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")]), loss_function = loss_root_mean_square, type = "raw" )
    
    cat("Root mean squared dropout loss due to a feature elimination:\n")
    print(fitnb_MDR)
  }
  else if (method == "model agnostic-LDA"){
    model_ldadalex_MDR <- train(MDR ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                                  Cefepime + Ceftazidime + Ciprofloxacin +
                                  Meropenem + Ceftriaxone + Aztreonam,
                                data = train,
                                method = "lda",
                                metric = "Accuracy")
    
    explda_MDR = explain(model = model_ldadalex_MDR, data=train[,c("Amikacin", "Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                   "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")], variables = colnames(train[,c("Amikacin", "Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                                                                                                                           "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")]),
                         y=as.vector(as.numeric(train$MDR))-1, label = "LDA", 
                         type = "classification")
    
    fitlda_MDR = variable_importance(explda_MDR, B=50 ,variables = colnames(train[,c("Amikacin", "Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                                     "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")]), loss_function = loss_root_mean_square, type = "raw" )
    cat("Root mean squared dropout loss due to a feature elimination:\n")
    print(fitlda_MDR)
  }
  else if (method == "model agnostic-Sig SVM"){
    set.seed(1234)
    tune_out_MDR <- e1071::tune("svm", MDR ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                                  Cefepime + Ceftazidime + Ciprofloxacin +
                                  Meropenem + Ceftriaxone + Aztreonam
                                , data = train, kernel = "sigmoid", tunecontrol=tune.control(cross=10),
                                ranges = list(cost = seq(0.1, 10, length = 20)
                                              , gamma = seq(0.1, 5, length = 20)))
    

    model_svm_MDR = svm(MDR ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                          Cefepime + Ceftazidime + Ciprofloxacin +
                          Meropenem + Ceftriaxone + Aztreonam, data = train, 
                        probability = T, kernel = "sigmoid", cost = tune_out_MDR$best.parameters$cost ,
                        gamma = tune_out_MDR$best.parameters$gamma)
    
    explained_SVM_MDR <- explain(model = model_svm_MDR, data=train[,c(2:9,11:12)], variables = colnames(train[,c(2:9,11:12)]),
                                 y=as.vector(as.numeric(train$MDR))-1, label = "SVM-Sigmoid",
                                 type = "classification")
    
    #50 permuatation
    fi_SVM_MDR = variable_importance(explained_SVM_MDR, B=50,variables = colnames(train[,c(2:9,11:12)]), loss_function = loss_root_mean_square, type = "raw" )
    
    cat("Root mean squared dropout loss due to a feature elimination:\n")
    print(fi_SVM_MDR)
  }
  else if (method == "model agnostic-Poly SVM"){
    
    set.seed(1234)
    tune_out5_MDR <- e1071::tune("svm", MDR ~ Imipenem + Gentamicin + Ciprofloxacin + Meropenem + 
                                   Aztreonam + Amikacin + Cefotaxime + Cefepime + Ceftazidime + 
                                   Ceftriaxone  
                                 , data = train, kernel = "polynomial", tunecontrol=tune.control(cross=10),
                                 ranges = list(degree = c(1, 2, 3), cost = seq(0.1, 10, length = 20)
                                               , gamma = seq(0.1, 5, length = 20)))
    

    model_svm5_MDR = svm(MDR ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                           Cefepime + Ceftazidime + Ciprofloxacin +
                           Meropenem + Ceftriaxone + Aztreonam, data = train, 
                         probability = T, kernel = "polynomial", cost = tune_out5_MDR$best.parameters$cost,
                         degree = tune_out5_MDR$best.parameters$degree , gamma = tune_out5_MDR$best.parameters$gamma )
    
    explained_SVM2_MDR <- explain(model = model_svm5_MDR, data=train[,c(2:9,11:12)], variables = colnames(train[,c(2:9,11:12)]),
                                  y=as.vector(as.numeric(train$MDR))-1, label = "SVM-Polynomial",
                                  type = "classification")
    
    #50 permuatation
    fi_SVM_MDR2 = variable_importance(explained_SVM2_MDR, B=50,variables = colnames(train[,c(2:9,11:12)]), loss_function = loss_root_mean_square, type = "raw" )
    
    print(fi_SVM_MDR2)
  }
  else if (method == "model agnostic-CatBoost"){
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
      depth = c(2,3,4),
      learning_rate = c(0.1,0.2,0.3),
      l2_leaf_reg = c(1,2,3),
      rsm = 1,
      border_count = 1,
      iterations = seq(20,40,1)
    )
    fitControl <- trainControl(method = "cv",
                               number = 10,
                               classProbs = TRUE)
    set.seed(1234)
    model_catboost_MDR <- train(x = train[,c("Amikacin", "Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                             "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")],
                                y = make.names(train[,"MDR"]),
                                maximize = TRUE,
                                method = catboost.caret, metric = "Accuracy", 
                                tuneGrid =  grid,
                                trControl = fitControl)
    
    expsda_MDR = explain(model = model_catboost_MDR, data=train[,c("Amikacin", "Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                   "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")], variables = colnames(train[,c("Amikacin", "Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                                                                                                                           "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")]),
                         y=as.vector(as.numeric(train$MDR))-1, label = "CatBoost", 
                         type = "classification")
    
    fitcat_MDR = variable_importance(expsda_MDR, B=50 ,variables = colnames(train[,c("Amikacin", "Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                                     "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")]), loss_function = loss_root_mean_square, type = "raw" )
    print(fitcat_MDR)
  }
  else {
    cat("Invalid task. Please specify one of the following tasks: feature_selection, model_evaluation\n")
    
  } 
}
if (is.null(task)){
  cat("You should specify a task\n")
}
