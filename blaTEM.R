<<<<<<< HEAD
source('F1-score.R') # Call F1-Score function for computing point estimation and its 95% confidence interval

source('compute_metrics.R') # Call compute_metrics function for computing all classification metrics and their confidence interval.

source('Packages_installation.R') # Installation of all required packages
=======
source('F1-score.R')
source('compute_metrics.R')
>>>>>>> c18c5be1391a6c1982e85029db52adc20aff265f

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

# Calculate the proportion of TEM_gene
proportion_TEM_gene <- prop.table(table(Bio_Data$TEM_gene))
cat("Proportion of TEM_gene:\n")
print(proportion_TEM_gene)


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

    if (model == "LR-Chi Squared test"){
      # Logistic regression (Model 1)--------------------------------------------------
      model_Log1_TEM <- glm(TEM_gene ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                              Cefepime + Ceftazidime + Ciprofloxacin +
                              Meropenem + Ceftriaxone + Aztreonam 
                            , family = "binomial", data = train)
      cat("Summary of LR model based on Chi Squared test including Null deviance, Residual deviance, and AIC, etc.:\n")
      print(summary(model_Log1_TEM))
      
      test$probs1_TEM <- predict(model_Log1_TEM, test, type = "response")
      test$pred_logreg1_TEM <- ifelse(test$probs1_TEM >= 0.5, 1, 0)
      # Model 1 evaluation
      confm_logreg1_TEM <- table(actual = test$TEM_gene, prediction = test$pred_logreg1_TEM)
      metrics_logreg1_TEM <- compute_metrics(confm_logreg1_TEM)
      cat("Performance of LR model based on Chi-Squared test results:\n")
      metrics_logreg1_TEM
    }
    else if (model == "LR-model agnostic/Wald test"){
      model_Log2_TEM <- glm(TEM_gene ~ Ciprofloxacin + Cefotaxime
                            , family = "binomial", data = train)
      cat("Summary of LR model based on model agnostic approach/Wald test including Null deviance, Residual deviance, and AIC, etc.:\n")
      print(summary(model_Log2_TEM))
      
      ##Prediction on test (Model 2)--------------------------------------------------
      test$probs2_TEM <- predict(model_Log2_TEM, test, type = "response")
      
      test$pred_logreg2_TEM <- ifelse(test$probs2_TEM >= 0.5, 1, 0)
      
      # Model 2 evaluation
      confm_logreg2_TEM <- table(actual = test$TEM_gene, prediction = test$pred_logreg2_TEM)
      
      metrics_logreg2_TEM <- compute_metrics(confm_logreg2_TEM)
      
      cat("Performance of LR model based on model agnostic approach / Wald test results:\n")
      metrics_logreg2_TEM
    }
    else if (model == "LR-Simplicity Principle_CTX"){
      model_Log3_TEM <- glm(TEM_gene ~ Cefotaxime
                            , family = "binomial", data = train)
      
      cat("Summary of LR model based on simplicity principle (CTX) including Null deviance, Residual deviance, and AIC, etc.:\n")
      print(summary(model_Log3_TEM))
      
      test$probs3_TEM <- predict(model_Log3_TEM, test, type = "response")
      
      test$pred_logreg3_TEM <- ifelse(test$probs3_TEM >= 0.5, 1, 0)
      
      confm_logreg3_TEM <- table(actual = test$TEM_gene, prediction = test$pred_logreg3_TEM)
      
      metrics_logreg3_TEM <- compute_metrics(confm_logreg3_TEM)
      
      cat("Performance of LR model based on simplicity principle (CTX):\n")
      metrics_logreg3_TEM
    }
    else if (model == "LR-Simplicity Principle_CIP"){
      model_Log4_TEM <- glm(TEM_gene ~ Ciprofloxacin
                            , family = "binomial", data = train)
      
      cat("Summary of LR model based on simplicity principle (CIP) including Null deviance, Residual deviance, and AIC, etc.:\n")
      print(summary(model_Log4_TEM))
      
      ##Prediction on test
      test$probs4_TEM <- predict(model_Log4_TEM, test, type = "response")
      
      test$pred_logreg4_TEM <- ifelse(test$probs4_TEM >= 0.5, 1, 0)
      
      confm_logreg4_TEM <- table(actual = test$TEM_gene, prediction = test$pred_logreg4_TEM)
      
      metrics_logreg4_TEM <- compute_metrics(confm_logreg4_TEM)
      
      cat("Performance of LR model based on simplicity principle (CIP):\n")
      metrics_logreg4_TEM
    }
    else if (model == "NBC-Chi Squared test"){
      train_control <- trainControl(
        method = "cv", 
        number = 10
      )
      
      search_grid_TEM <- expand.grid(
        usekernel = c(FALSE,TRUE),
        fL = seq(0,5,0.5),
        adjust = seq(0,5,0.5)
      )
      
      model_nbdalex_TEM <- train(TEM_gene ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                                   Cefepime + Ceftazidime + Ciprofloxacin +
                                   Meropenem + Ceftriaxone + Aztreonam,
                                 data = train,
                                 method = "nb",
                                 metric = "Accuracy",
                                 trControl = train_control,
                                 tuneGrid = search_grid_TEM)
      
      cat("Tuned hyper-parameters through 10-fold Cross-Validation:\n")
      print(model_nbdalex_TEM$bestTune)
      
      ##Prediction on test (Naive Bayes Model)----------------------------------
      test$pred_nb_TEM <- predict(model_nbdalex_TEM, test)
      
      # Model evaluation
      confm_nb_TEM <- table(actual = test$TEM_gene, prediction = test$pred_nb_TEM)
      metrics_nb_TEM<- compute_metrics(confm_nb_TEM)
      cat("Performance of NBC model based on Chi-Squared test results:\n")
      metrics_nb_TEM
    }
    else if (model == "NBC-model agnostic"){
      train_control <- trainControl(
        method = "cv",
        number = 10
      )
      
      search_grid_TEM <- expand.grid(
        usekernel = c(FALSE,TRUE),
        fL = seq(0,5,0.5),
        adjust = seq(0,5,0.5)
      )
      
      model_nbdalex1_TEM <- train(TEM_gene ~ Ciprofloxacin,
                                  data = train,
                                  method = "nb",
                                  metric = "Accuracy",
                                  trControl = train_control,
                                  tuneGrid = search_grid_TEM)
      
      cat("Tuned hyper-parameters through 10-fold Cross-Validation:\n")
      print(model_nbdalex1_TEM$bestTune)
      
      ##Prediction on test (Naive Bayes Model)----------------------------------
      test$pred_nb1_TEM <- predict(model_nbdalex1_TEM, test)
      
      # Model evaluation
      confm_nb1_TEM <- table(actual = test$TEM_gene, prediction = test$pred_nb1_TEM)
      metrics_nb1_TEM<- compute_metrics(confm_nb1_TEM)
      cat("Performance of NBC model based on model-agnostic approach results:\n")
      metrics_nb1_TEM
    }
    
    else if (model == "NBC-Wald test"){
      
      train_control <- trainControl(
        method = "cv",
        number = 10
      )
      
      search_grid_TEM <- expand.grid(
        usekernel = c(FALSE,TRUE),
        fL = seq(0,5,0.5),
        adjust = seq(0,5,0.5)
      )
      
      model_nb2_TEM <- train(TEM_gene ~ (Cefotaxime+ Ciprofloxacin),
                             data = train,
                             method = "nb",
                             metric = "Accuracy",
                             trControl = train_control,
                             tuneGrid = search_grid_TEM)
      
      cat("Tuned hyper-parameters through 10-fold Cross-Validation:\n")
      print(model_nb2_TEM$bestTune)
      
      ##Prediction on test (Naive Bayes Model)
      test$pred_nb2_TEM <- predict(model_nb2_TEM, test)
      
      ##confusion matrix(Naive Bayes Model)
      confm_nb2_TEM <- table(actual = test$TEM_gene, prediction = test$pred_nb2_TEM)
      
      metrics_nb2_TEM<- compute_metrics(confm_nb2_TEM)
      
      cat("Performance of NBC model based on Wald test results:\n")
      metrics_nb2_TEM
    }
    
    else if (model == "NBC-Simplicity Principle_CTX"){
      
      train_control <- trainControl(
        method = "repeatedcv",
        repeats = 5,
        number = 10
      )
      
      search_grid_TEM <- expand.grid(
        usekernel = c(TRUE,FALSE),
        fL = seq(0,0.25,3),
        adjust = seq(0,0.25,3)
      )
      
      model_nb3_TEM <- train(TEM_gene ~ Cefotaxime,
                             data = train,
                             method = "nb",
                             metric = "Accuracy",
                             trControl = train_control,
                             tuneGrid = search_grid_TEM)
      
      cat("Tuned hyper-parameters through 10-fold Cross-Validation:\n")
      print(model_nb3_TEM$bestTune)
      
      ##Prediction on test (Naive Bayes Model)
      test$pred_nb3_TEM <- predict(model_nb3_TEM, test)
      
      ##confusion matrix(Naive Bayes Model)
      confm_nb3_TEM <- table(actual = test$TEM_gene, prediction = test$pred_nb3_TEM)
      
      metrics_nb3_TEM<- compute_metrics(confm_nb3_TEM)
      
      cat("Performance of NBC model based on Simplicity Principle (CTX):\n")
      metrics_nb3_TEM
    }
      
    else if (model == "LDA-Chi Squared test"){
      #Linear Discriminant Analysis-------------------------------------------------
      model_ldadalex_TEM <- train(TEM_gene ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                                    Cefepime + Ceftazidime + Ciprofloxacin +
                                    Meropenem + Ceftriaxone + Aztreonam,
                                  data = train,
                                  method = "lda",
                                  metric = "Accuracy")
      
      cat("LDA model based on Chi-Squared test results:\n")
      print(model_ldadalex_TEM)
      
      #Prediction on test (LDA Model)-----------------------------------------------
      test$pred_lda_TEM <- predict(model_ldadalex_TEM, test)
      
      # Model Evaluation-------------------------------------------------
      confm_lda_TEM <- table(actual = test$TEM_gene, prediction = test$pred_lda_TEM)
      metrics_lda_TEM<- compute_metrics(confm_lda_TEM)
      cat("Performance of LDA model based on Chi-Squared test results:\n")
      metrics_lda_TEM
    }
    else if (model == "LDA-model agnostic/Wald test"){
      model_ldadalex1_TEM <- train(TEM_gene ~ Cefotaxime + Ciprofloxacin,
                                   data = train,
                                   method = "lda",
                                   metric = "Accuracy")
      
      cat("LDA model based on model-agnostic approach/Wald test results:\n")
      print(model_ldadalex1_TEM)
      
      #Prediction on test (LDA Model)-----------------------------------------------
      test$pred_lda1_TEM <- predict(model_ldadalex1_TEM, test)  
      
      # Model evaluation
      confm_lda1_TEM <- table(actual = test$TEM_gene, prediction = test$pred_lda1_TEM)
      metrics_lda1_TEM<- compute_metrics(confm_lda1_TEM)
      cat("Performance of LDA model based on model-agnostic approach/Wald test results:\n")
      metrics_lda1_TEM
    }
    
    else if (model == "LDA-Simplicity Principle_CIP"){
      model_lda2_TEM <- train(TEM_gene ~ Ciprofloxacin,
                              data = train,
                              method = "lda",
                              metric = "Accuracy",
                              trControl = train_control)
      
      cat("LDA model based on Simplicity Principle (CIP):\n")
      print(model_lda2_TEM)
      
      test$pred_lda2_TEM <- predict(model_lda2_TEM, test)
      
      confm_lda2_TEM <- table(actual = test$TEM_gene, prediction = test$pred_lda2_TEM)
      
      metrics_lda2_TEM<- compute_metrics(confm_lda2_TEM)
      
      cat("Performance of LDA model based on Simplicity Principle (CIP):\n")
      metrics_lda2_TEM
    }
     
    else if (model == "LDA-Simplicity Principle_CTX"){
      
      model_lda3_TEM <- train(TEM_gene ~ Cefotaxime,
                              data = train,
                              method = "lda",
                              metric = "Accuracy",
                              trControl = train_control)
      
      cat("LDA model based on Simplicity Principle (CTX):\n")
      print(model_lda3_TEM)
      
      test$pred_lda3_TEM <- predict(model_lda3_TEM, test)
      
      confm_lda3_TEM <- table(actual = test$TEM_gene, prediction = test$pred_lda3_TEM)
      
      metrics_lda3_TEM<- compute_metrics(confm_lda3_TEM)
      
      cat("Performance of LDA model based on Simplicity Principle (CTX):\n")
      metrics_lda3_TEM
    } 
    else if (model == "SVM-Sig-Chi Squared test"){
      #Fit Support Vector Machine model to train data set
      #This code performs a grid search for tuning the hyper-parameters of the SVM model using the tune function
      # from the e1071 package. It then fits the SVM model with the specified hyper-parameters using the svm function.
      # Model 1
      set.seed(1234)
      tune_out_TEM <- e1071::tune("svm", TEM_gene ~ Amikacin + Cefotaxime + Gentamicin + 
                                    Imipenem + Cefepime + Ceftazidime + Ciprofloxacin + Meropenem + 
                                    Ceftriaxone + Aztreonam
                                  , data = train, kernel = "sigmoid", tunecontrol=tune.control(cross=10),
                                  ranges = list(cost = seq(0.1, 5, length = 25)
                                                , gamma = seq(0.1, 5, length = 25)))
      
      cat("Tuned hyper-parameters throgh 10-Fold Cross Validation:\n")
      print(summary(tune_out_TEM))
      
      model_svm_TEM = svm(TEM_gene ~ Amikacin + Cefotaxime + Gentamicin + 
                            Imipenem + Cefepime + Ceftazidime + Ciprofloxacin + Meropenem + 
                            Ceftriaxone + Aztreonam, data = train, 
                          probability = T, kernel = "sigmoid", cost = tune_out_TEM$best.parameters$cost,
                          gamma = tune_out_TEM$best.parameters$gamma )
      
      cat("Sigmoid SVM based on Chi-Squared test results:\n")
      print(model_svm_TEM)
      
      #Prediction on test (SVM)----------------------------------------------------
      test$pred_svm_TEM <- predict(model_svm_TEM, test)
      
      # Model evaluation
      confm_svm_TEM <- table(actual = test$TEM_gene, prediction = test$pred_svm_TEM)
      
      metrics_svm1_TEM<- compute_metrics(confm_svm_TEM)
      
      cat("Performance of SVM-Sigmoid model based on Chi-Squared test results:\n")
      metrics_svm1_TEM
    }
    else if (model=="SVM-Sig-model agnostic"){
      
      set.seed(1234)
      tune_out1_TEM <- e1071::tune("svm", TEM_gene ~ Amikacin +  Meropenem + Aztreonam + Ceftriaxone
                                   , data = train, kernel = "sigmoid", tunecontrol=tune.control(cross=10),
                                   ranges = list(cost = seq(0.1, 5, length = 20)
                                                 , gamma = seq(0, 5, length = 20)))
      cat("Tuned hyper-parameters throgh 10-Fold Cross Validation:\n")
      print(summary(tune_out1_TEM))
      
      model_svm_TEM1 = svm(TEM_gene ~ Amikacin +  Meropenem + Aztreonam + Ceftriaxone, data = train, 
                           probability = T, kernel = "sigmoid", cost = tune_out1_TEM$best.parameters$cost,
                           gamma = tune_out1_TEM$best.parameters$gamma)
      
      cat("Sigmoid SVM based on model-agnostic approach results:\n")
      print(model_svm_TEM1)
      
      test$pred_svm_TEM1 <- predict(model_svm_TEM1, test)
      
      confm_svm_TEM1 <- table(actual = test$TEM_gene, prediction = test$pred_svm_TEM1)
      
      metrics_svm2_TEM<- compute_metrics(confm_svm_TEM1)
      
      cat("Performance of SVM-Sigmoid model based on model-agnostic approach results:\n")
      metrics_svm2_TEM
    }
    else if (model=="SVM-Sig-Simplicity Principle_CTX"){
      set.seed(1234)
      tune_out2_TEM <- e1071::tune("svm", TEM_gene ~ Cefotaxime
                                   , data = train, kernel = "sigmoid", tunecontrol=tune.control(cross=10),
                                   ranges = list(cost = seq(0.01, 5, length = 20)
                                                 , gamma = seq(0.01, 5, length = 20)))
      
      cat("Tuned hyper-parameters throgh 10-Fold Cross Validation:\n")
      print(summary(tune_out2_TEM))
      
      model_svm_TEM2 = svm(TEM_gene ~ Cefotaxime, data = train, 
                           probability = T, kernel = "sigmoid", cost = tune_out2_TEM$best.parameters$cost,
                           gamma = tune_out2_TEM$best.parameters$gamma )
      
      cat("Sigmoid SVM based on Simplicity Principle (CTX):\n")
      print(model_svm_TEM2)
      
      #Prediction on test (SVM)
      test$pred_svm_TEM2 <- predict(model_svm_TEM2, test)
      
      # Model evaluation
      confm_svm2_TEM <- table(actual = test$TEM_gene, prediction = test$pred_svm_TEM2)
      
      metrics_svm3_TEM<- compute_metrics(confm_svm2_TEM)
      
      cat("Performance of SVM-Sigmoid model based on Simplicity Principle (CTX):\n")
      metrics_svm3_TEM
      
    } 
    
    else if (model=="SVM-Sig-Simplicity Principle_CIP"){
      set.seed(1234)
      tune_out3_TEM <- e1071::tune("svm", TEM_gene ~ Ciprofloxacin
                                   , data = train, kernel = "sigmoid", tunecontrol=tune.control(cross=10),
                                   ranges = list(cost = seq(0.1, 5, length = 20)
                                                 , gamma = seq(0.01, 5, length = 20)))
      
      cat("Tuned hyper-parameters throgh 10-Fold Cross Validation:\n")
      print(summary(tune_out3_TEM))
      
      model_svm_TEM3 = svm(TEM_gene ~ Ciprofloxacin, data = train, 
                           probability = T, kernel = "sigmoid", cost = tune_out3_TEM$best.parameters$cost,
                           gamma = tune_out3_TEM$best.parameters$gamma )
      
      cat("Sigmoid SVM based on Simplicity Principle (CIP):\n")
      print(model_svm_TEM3)
      
      #Prediction on test (SVM)
      test$pred_svm_TEM3 <- predict(model_svm_TEM3, test)
      
      # Model evaluation
      confm_svm3_TEM <- table(actual = test$TEM_gene, prediction = test$pred_svm_TEM3)
      
      metrics_svm4_TEM<- compute_metrics(confm_svm3_TEM)
      
      cat("Performance of SVM-Sigmoid model based on Simplicity Principle (CIP):\n")
      metrics_svm4_TEM
    } 
    
    else if (model=="SVM-Sig-Wald test"){
      
      set.seed(1234)
      tune_out4_TEM <- e1071::tune("svm", TEM_gene ~ Cefotaxime + Ciprofloxacin
                                   , data = train, kernel = "sigmoid", tunecontrol=tune.control(cross=10),
                                   ranges = list(cost = seq(0.01, 5, length = 20)
                                                 , gamma = seq(0.01, 5, length = 20)))
      
      cat("Tuned hyper-parameters throgh 10-Fold Cross Validation:\n")
      print(summary(tune_out4_TEM))
      
      model_svm_TEM4 = svm(TEM_gene ~ Cefotaxime + Ciprofloxacin, data = train, 
                           probability = T, kernel = "sigmoid", cost = tune_out4_TEM$best.parameters$cost,
                           gamma = tune_out4_TEM$best.parameters$gamma )
      
      cat("Sigmoid SVM based on the Wald test results:\n")
      print(model_svm_TEM4)
      
      #Prediction on test (SVM)
      test$pred_svm_TEM4 <- predict(model_svm_TEM4, test)
      
      # Model evaluation
      confm_svm4_TEM <- table(actual = test$TEM_gene, prediction = test$pred_svm_TEM4)
      
      metrics_svm5_TEM<- compute_metrics(confm_svm4_TEM)
      
      cat("Performance of SVM-Sigmoid model based on the Wald test results:\n")
      metrics_svm5_TEM
    }
    
    else if (model=="SVM-Poly-Chi Squared test"){
      
      set.seed(1234)
      tune_out5_TEM <- e1071::tune("svm", TEM_gene ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                                     Cefepime + Ceftazidime + Ciprofloxacin +
                                     Meropenem + Ceftriaxone + Aztreonam  
                                   , data = train, kernel = "polynomial", tunecontrol=tune.control(cross=10),
                                   ranges = list(degree = c(1, 2, 3), cost = seq(0.1, 5, length = 20)
                                                 , gamma = seq(0.1, 5, length = 20)))
      
      cat("Tuned hyper-parameters throgh 10-Fold Cross Validation:\n")
      print(summary(tune_out5_TEM))
      
      model_svm_TEM5 = svm(TEM_gene ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                             Cefepime + Ceftazidime + Ciprofloxacin +
                             Meropenem + Ceftriaxone + Aztreonam, data = train, 
                           probability = T, kernel = "polynomial", cost = tune_out5_TEM$best.parameters$cost,
                           degree = tune_out5_TEM$best.parameters$degree , gamma = tune_out5_TEM$best.parameters$gamma )
      
      cat("Polynomial SVM based on Chi-Squared test results:\n")
      print(model_svm_TEM5)
      
      #Prediction on test (SVM)----------------------------------------------------
      test$pred_svm5_TEM <- predict(model_svm_TEM5, test)
      
      # Model evaluation
      confm_svm5_TEM <- table(actual = test$TEM_gene, prediction = test$pred_svm5_TEM)
      
      metrics_svm6_TEM<- compute_metrics(confm_svm5_TEM)
      
      cat("Performance of SVM-Polynomial model based on Chi-Squared test results:\n")
      metrics_svm6_TEM
    }
    else if (model=="SVM-Poly-model agnostic"){
      set.seed(1234)
      tune_out6_TEM <- e1071::tune("svm", TEM_gene ~ Meropenem + Gentamicin + Aztreonam + Amikacin
                               + Imipenem + Ceftriaxone
                               , data = train, kernel = "polynomial", tunecontrol=tune.control(cross=10),
                               ranges = list(cost = seq(0.1, 5, length = 20), degree = c(1, 2, 3)
                                             , gamma = seq(0.1, 5, length = 20)))
      
      cat("Tuned hyper-parameters throgh 10-Fold Cross Validation:\n")
      print(summary(tune_out6_TEM))
      
      # Tuned hyper-parameters are placed in the model.
      model_svm_TEM51 = svm(TEM_gene ~ Gentamicin + Amikacin + Ceftriaxone + Cefepime + Ciprofloxacin,
                           data = train, probability = T, kernel = "polynomial", cost = tune_out6_TEM$best.parameters$cost, 
                           degree = tune_out6_TEM$best.parameters$degree , gamma = tune_out6_TEM$best.parameters$gamma )
      
      cat("Polynomial SVM based on model agnostic approach results:\n")
      print(model_svm_TEM51)
      
      #Prediction on test (SVM)----------------------------------------------------
      test$pred_svm51_TEM <- predict(model_svm_TEM51, test)
      
      # Model evaluation
      confm_svm51_TEM <- table(actual = test$TEM_gene, prediction = test$pred_svm51_TEM)
      
      metrics_svm7_TEM<- compute_metrics(confm_svm51_TEM)
      
      cat("Performance of SVM-Polynomial model based on model-agnostic approach results:\n")
      metrics_svm7_TEM
    }
    else if (model=="SVM-Poly-Wald test"){
      
      set.seed(1234)
      tune_out7_TEM <- e1071::tune("svm", TEM_gene ~ Cefotaxime + Ciprofloxacin
                               , data = train, kernel = "polynomial", tunecontrol=tune.control(cross=10),
                               ranges = list(cost = seq(0.01, 5, length = 20), degree = c(1, 2, 3)
                                             , gamma = seq(0.01, 5, length = 20)))
      
      cat("Tuned hyper-parameters throgh 10-Fold Cross Validation:\n")
      print(summary(tune_out7_TEM))
      
      model_svm_TEM6 = svm(TEM_gene ~ Cefotaxime + Ciprofloxacin, data = train, 
                           probability = T, kernel = "polynomial", cost = tune_out7_TEM$best.parameters$cost,
                           degree = tune_out7_TEM$best.parameters$degree , gamma = tune_out7_TEM$best.parameters$gamma )
      
      cat("Polynomial SVM based on Wald test results:\n")
      print(model_svm_TEM6)
      
      #Prediction on test (SVM)----------------------------------------------------
      test$pred_svm6 <- predict(model_svm_TEM6, test)
      
      # Model evaluation
      confm_svm6_TEM <- table(actual = test$TEM_gene, prediction = test$pred_svm6)
      
      metrics_svm8_TEM<- compute_metrics(confm_svm6_TEM)
      
      cat("Performance of SVM-Polynomial model based on Wald test results:\n")
      metrics_svm8_TEM
      
    }
    
    else if (model=="SVM-Poly-Simplicity principle_CTX"){
      
      set.seed(1234)
      tune_out8_TEM <- e1071::tune("svm", TEM_gene ~ Cefotaxime
                                   , data = train, kernel = "polynomial", tunecontrol=tune.control(cross=10),
                                   ranges = list(cost = seq(0.01, 5, length = 20), degree = c(1, 2, 3)
                                                 , gamma = seq(0.01, 5, length = 20)))
      
      cat("Tuned hyper-parameters throgh 10-Fold Cross Validation:\n")
      print(summary(tune_out8_TEM))
      
      model_svm_TEM7 = svm(TEM_gene ~ Cefotaxime , data = train, 
                           probability = T, kernel = "polynomial", cost = tune_out8_TEM$best.parameters$cost,
                           degree = tune_out8_TEM$best.parameters$degree , gamma = tune_out8_TEM$best.parameters$gamma )
      
      cat("Polynomial SVM based on the Simplicity Principle (CTX):\n")
      print(model_svm_TEM7)
      
      #Prediction on test (SVM)----------------------------------------------------
      test$pred_svm7_TEM <- predict(model_svm_TEM7, test)
      
      # Model evaluation
      confm_svm7_TEM <- table(actual = test$TEM_gene, prediction = test$pred_svm7_TEM)
      
      metrics_svm9_TEM<- compute_metrics(confm_svm7_TEM)
      
      cat("Performance of SVM-Polynomial model based on the Simplicity Principle (CTX):\n")
      metrics_svm9_TEM
    }
    
    else if (model=="SVM-Poly-Simplicity principle_CIP"){
      
      set.seed(1234)
      tune_out9_TEM <- e1071::tune("svm", TEM_gene ~ Ciprofloxacin
                                   , data = train, kernel = "polynomial", tunecontrol=tune.control(cross=10),
                                   ranges = list(cost = seq(0.01, 5, length = 20), degree = c(1, 2, 3)
                                                 , gamma = seq(0.01, 5, length = 20)))
      
      cat("Tuned hyper-parameters throgh 10-Fold Cross Validation:\n")
      print(summary(tune_out9_TEM))
      
      model_svm_TEM8 = svm(TEM_gene ~ Ciprofloxacin , data = train, 
                           probability = T, kernel = "polynomial", cost = tune_out9_TEM$best.parameters$cost,
                           degree = tune_out9_TEM$best.parameters$degree , gamma = tune_out9_TEM$best.parameters$gamma )
      
      cat("Polynomial SVM based on the Simplicity Principle (CIP):\n")
      print(model_svm_TEM8)
      
      #Prediction on test (SVM)----------------------------------------------------
      test$pred_svm8_TEM <- predict(model_svm_TEM8, test)
      
      # Model evaluation
      confm_svm8_TEM <- table(actual = test$TEM_gene, prediction = test$pred_svm8_TEM)
      
      metrics_svm10_TEM<- compute_metrics(confm_svm8_TEM)
      
      cat("Performance of SVM-Polynomial model based on the Simplicity Principle (CIP):\n")
      metrics_svm10_TEM
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
        depth = c(2,3,4),
        learning_rate = c(0.1,0.2,0.3),
        l2_leaf_reg = c(2,3,4),
        rsm = 1,
        border_count = 1,
        iterations = seq(40,60,1)
      )
      fitControl <- trainControl(method = "cv",
                                 number = 10,
                                 classProbs = TRUE)
      
      set.seed(1234)
      
      model_catboost_TEM <- train(x = train[,c("Amikacin", "Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                               "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")],
                                  y = make.names(train[,"TEM_gene"]),
                                  maximize = TRUE,
                                  method = catboost.caret, metric = "Accuracy", 
                                  tuneGrid =  grid,
                                  trControl = fitControl)
      
      cat("Tuned hyper-parameters throgh 10-Fold Cross Validation:\n")
      print(model_catboost_TEM$bestTune)  
      
      #Prediction on test (CatBoost)----------------------------------------------------
      test$cat_TEM = predict(model_catboost_TEM, test)
      
      # Model evaluation
      cat_test_TEM <- table(actual = test$TEM_gene, prediction = test$cat_TEM)
      
      metrics_cat_TEM<- compute_metrics(cat_test_TEM)
      
      cat("Performance of CatBoost model based on Chi-Squared test results:\n")
      metrics_cat_TEM
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
        iterations = seq(10,20,1)
      )
      fitControl2 <- trainControl(method = "cv",
                                  number = 10,
                                  classProbs = TRUE)
      
      set.seed(1234)
      
      model_catboost2_TEM <- train(x = train["Ciprofloxacin"],
                                   y = make.names(train[["TEM_gene"]]),
                                   maximize = TRUE,
                                   method = catboost.caret, metric = "Accuracy", 
                                   tuneGrid =  grid2,
                                   trControl = fitControl2)
      
      cat("Tuned hyper-parameters throgh 10-Fold Cross Validation:\n")
      print(model_catboost2_TEM$bestTune)  
      
      test$cat_TEM2 = predict(model_catboost2_TEM, test)
      
      # Model Evaluation
      confm_cat_TEM <- table(actual = test$TEM_gene, prediction = test$cat_TEM)
      
      metrics_cat1_TEM<- compute_metrics(confm_cat_TEM)
      
      cat("Performance of CatBoost model based on model-agnostic approach results:\n")
      metrics_cat1_TEM
    }
    
    else if (model=="CatBoost-Simplicity Principle_CTX"){
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
        depth = c(2,3,4),
        learning_rate = c(0.1,0.2,0.3),
        l2_leaf_reg = c(1,2,3,4),
        rsm = 1,
        border_count = 1,
        iterations = seq(10,700,50)
      )
      fitControl3 <- trainControl(method = "cv",
                                  number = 10,
                                  classProbs = TRUE)
      
      model_catboost3_TEM <- train(x = train["Cefotaxime"],
                                   y = make.names(train[["TEM_gene"]]),
                                   maximize = TRUE,
                                   method = catboost.caret, metric = "Accuracy", 
                                   tuneGrid =  grid3,
                                   trControl = fitControl3)
      
      cat("Tuned hyper-parameters throgh 10-Fold Cross Validation:\n")
      print(model_catboost3_TEM$bestTune)  
      
      test$cat_TEM3 = predict(model_catboost3_TEM, test)
      # Model evaluation
      confm_cat2_TEM <- table(actual = test$TEM_gene, prediction = test$cat_TEM3)  
      
      metrics_cat2_TEM<- compute_metrics(confm_cat2_TEM)
      
      cat("Performance of CatBoost model based on Wald test results:\n")
      metrics_cat2_TEM
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
      
      grid4 <- expand.grid(
        depth = c(2,3),
        learning_rate = c(0.1,0.2,0.3),
        l2_leaf_reg = c(1,2,3),
        rsm = 1,
        border_count = 1,
        iterations = seq(10,100,10)
      )
      fitControl4 <- trainControl(method = "cv",
                                  number = 10,
                                  classProbs = TRUE)
      model_catboost4_TEM <- train(x = train[,c("Cefotaxime", "Ciprofloxacin")],
                               y = make.names(train[,"TEM_gene"]),
                               maximize = TRUE,
                               method = catboost.caret, metric = "Accuracy", 
                               tuneGrid =  grid4,
                               trControl = fitControl4)
      
      cat("Tuned hyper-parameters throgh 10-Fold Cross Validation:\n")
      print(model_catboost4_TEM$bestTune)  
      
      test$cat_TEM4 = predict(model_catboost4_TEM, test)
      # Model evaluation
      confm_cat3_TEM <- table(actual = test$TEM_gene, prediction = test$cat_TEM4)  
      
      metrics_cat3_TEM<- compute_metrics(confm_cat3_TEM)
      
      cat("Performance of CatBoost model based on Wald test results:\n")
      metrics_cat3_TEM
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
        result <- chisq.test(table(data[[variable]], data$TEM_gene))
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
    model_Log1_TEM <- glm(TEM_gene ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                            Cefepime + Ceftazidime + Ciprofloxacin +
                            Meropenem + Ceftriaxone + Aztreonam 
                          , family = "binomial", data = train)
    
    cat("P-values show Wald test results:\n")
    print(summary(model_Log1_TEM))
  }
  else if (method == "model agnostic-LR"){
    # Feature Selection by DALEX Package--------------------------------------
    # the 'explained_glm_TEM' object is created using the 'explain' function from the DALEX package. 
    # The variable_importance function is used to calculate the permutation-based 
    # feature importance ('fi_glm_TEM') with 50 permutations. Finally, the feature importance is displayed 
    # using 'fi_glm_TEM' and plotted using' plot(fi_glm_TEM)'.
    model_Log1_TEM <- glm(TEM_gene ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                            Cefepime + Ceftazidime + Ciprofloxacin +
                            Meropenem + Ceftriaxone + Aztreonam 
                          , family = "binomial", data = train)
    
    explained_glm_TEM <- explain(model = model_Log1_TEM, data=train[,c(2:9,11:12)], variables = colnames(train[,c(2:9,11:12)]),
                                 y=as.vector(as.numeric(train$TEM_gene))-1, label = "LR", 
                                 type = "classification")
    #50 permuatation
    fi_glm_TEM = variable_importance(explained_glm_TEM, B=50 ,variables = colnames(train[,c(2:9,11:12)]),loss_function = loss_root_mean_square, type = "raw")

    cat("Root mean squared dropout loss due to a feature elimination through model-agnostic approach:\n")
    print(fi_glm_TEM)
  }
  else if (method == "model agnostic-NBC"){
    train_control <- trainControl(
      method = "cv", 
      number = 10
    )
    
    search_grid_TEM <- expand.grid(
      usekernel = c(FALSE,TRUE),
      fL = seq(0,5,0.5),
      adjust = seq(0,5,0.5)
    )
    
    model_nbdalex_TEM <- train(TEM_gene ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                                 Cefepime + Ceftazidime + Ciprofloxacin +
                                 Meropenem + Ceftriaxone + Aztreonam,
                               data = train,
                               method = "nb",
                               metric = "Accuracy",
                               trControl = train_control,
                               tuneGrid = search_grid_TEM)
    
    # Feature selection through model-agnostic approach (DALEX)
    # the 'explain' function from the DALEX package is used to explain the Naive Bayes classifier model 
    # 'model_nbdalex_TEM'. The relevant variables are selected and provided as data along with their 
    # corresponding column names. The target variable 'TEM_gene' is transformed to numeric values 
    # (y = as.vector(as.numeric(train$TEM_gene)) - 1). The label is set to "NBC" for Naive Bayes classifier,
    # and the type is set to "classification".
    expnb_TEM = explain(model = model_nbdalex_TEM, data=train[,c("Amikacin", "Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                 "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")], variables = colnames(train[,c("Amikacin", "Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                                                                                                                         "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")]),
                        y=as.vector(as.numeric(train$TEM_gene))-1, label = "NBC", 
                        type = "classification")
    
    # The 'variable_importance' function is then used to calculate the variable importance based on the
    # explained model. The 'B' parameter is set to 50 for the number of permutations, and the loss function
    # is set to "loss_root_mean_square". The resulting variable importance is stored in 'fitnb_TEM'
    fitnb_TEM = variable_importance(expnb_TEM, B=50 ,variables = colnames(train[,c("Amikacin", "Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                                   "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")]), loss_function = loss_root_mean_square, type = "raw" )
    
    cat("Root mean squared dropout loss due to a feature elimination through model-agnostic approach:\n")
    print(fitnb_TEM)
  }
  else if (method == "model agnostic-LDA"){
    model_ldadalex_TEM <- train(TEM_gene ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                                  Cefepime + Ceftazidime + Ciprofloxacin +
                                  Meropenem + Ceftriaxone + Aztreonam,
                                data = train,
                                method = "lda",
                                metric = "Accuracy")
    
    explda_TEM = explain(model = model_ldadalex_TEM, data=train[,c("Amikacin", "Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                   "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")], variables = colnames(train[,c("Amikacin", "Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                                                                                                                           "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")]),
                         y=as.vector(as.numeric(train$TEM_gene))-1, label = "LDA", 
                         type = "classification")
    
    fitlda_TEM = variable_importance(explda_TEM, B=50 ,variables = colnames(train[,c("Amikacin", "Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                                     "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")]), loss_function = loss_root_mean_square, type = "raw" )
    cat("Root mean squared dropout loss due to a feature elimination through model-agnostic approach:\n")
    print(fitlda_TEM)
  }
  else if (method == "model agnostic-Sig SVM"){
    set.seed(1234)
    tune_out_TEM <- e1071::tune("svm", TEM_gene ~ Amikacin + Cefotaxime + Gentamicin + 
                                  Imipenem + Cefepime + Ceftazidime + Ciprofloxacin + Meropenem + 
                                  Ceftriaxone + Aztreonam
                                , data = train, kernel = "sigmoid", tunecontrol=tune.control(cross=10),
                                ranges = list(cost = seq(0.1, 5, length = 25)
                                              , gamma = seq(0.1, 5, length = 25)))
    
    model_svm_TEM = svm(TEM_gene ~ Amikacin + Cefotaxime + Gentamicin + 
                          Imipenem + Cefepime + Ceftazidime + Ciprofloxacin + Meropenem + 
                          Ceftriaxone + Aztreonam, data = train, 
                        probability = T, kernel = "sigmoid", cost = tune_out_TEM$best.parameters$cost ,
                        gamma = tune_out_TEM$best.parameters$gamma)
    
    explained_SVM_TEM <- explain(model = model_svm_TEM, data=train[,c(2:9,11:12)], variables = colnames(train[,c(2:9,11:12)]),
                                 y=as.vector(as.numeric(train$TEM_gene))-1, label = "SVM-Sigmoid",
                                 type = "classification")
    
    #50 permuatation
    fi_SVM_TEM = variable_importance(explained_SVM_TEM, B=50,variables = colnames(train[,c(2:9,11:12)]), loss_function = loss_root_mean_square, type = "raw" )
    
    cat("Root mean squared dropout loss due to a feature elimination through model-agnostic approach:\n")
    print(fi_SVM_TEM)
  }
  else if (method == "model agnostic-Poly SVM"){
    # it is highly recommeded that first train the model, then specify following hyper-parameters according to the tuned hyper-parameters.
    # Following values of the hyper-parameters have been specified in this way.
    set.seed(1234)
    tune_out5_TEM <- e1071::tune("svm", TEM_gene ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                                   Cefepime + Ceftazidime + Ciprofloxacin +
                                   Meropenem + Ceftriaxone + Aztreonam  
                                 , data = train, kernel = "polynomial", tunecontrol=tune.control(cross=10),
                                 ranges = list(degree = c(1, 2, 3), cost = seq(0.1, 5, length = 20)
                                               , gamma = seq(0.1, 5, length = 20)))
    
    model_svm_TEM4 = svm(TEM_gene ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                           Cefepime + Ceftazidime + Ciprofloxacin +
                           Meropenem + Ceftriaxone + Aztreonam, data = train, 
                         probability = T, kernel = "polynomial", cost = tune_out5_TEM$best.parameters$cost,
                         degree = tune_out5_TEM$best.parameters$degree , gamma = tune_out5_TEM$best.parameters$gamma)
    
    explained_SVM_TEM2 <- explain(model = model_svm_TEM5, data=train[,c(2:9,11:12)], variables = colnames(train[,c(2:9,11:12)]),
                                  y=as.vector(as.numeric(train$TEM_gene))-1, label = "SVM-Polynomial",
                                  type = "classification")
    
    #50 permuatation
    fi_SVM_TEM2 = variable_importance(explained_SVM_TEM2, B=50,variables = colnames(train[,c(2:9,11:12)]), loss_function = loss_root_mean_square, type = "raw" )
    
    print(fi_SVM_TEM2)
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
      depth = c(2,3,4),
      learning_rate = c(0.1,0.2,0.3),
      l2_leaf_reg = c(2,3,4),
      rsm = 1,
      border_count = 1,
      iterations = seq(40,60,1)
    )
    fitControl <- trainControl(method = "cv",
                               number = 10,
                               classProbs = TRUE)
    
    set.seed(1234)
    model_catboost_TEM <- train(x = train[,c("Amikacin", "Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                             "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")],
                                y = make.names(train[,"TEM_gene"]),
                                maximize = TRUE,
                                method = catboost.caret, metric = "Accuracy", 
                                tuneGrid =  grid,
                                trControl = fitControl)
    
    expsda_TEM = explain(model = model_catboost_TEM, data=train[,c("Amikacin", "Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                   "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")], variables = colnames(train[,c("Amikacin", "Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                                                                                                                           "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")]),
                         y=as.vector(as.numeric(train$TEM_gene))-1, label = "CatBoost", 
                         type = "classification")
    
    fitcat_TEM = variable_importance(expsda_TEM, B=50 ,variables = colnames(train[,c("Amikacin", "Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                                     "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")]), loss_function = loss_root_mean_square, type = "raw" )
    print(fitcat_TEM)
  }
  else {
    cat("Invalid task. Please specify one of the following tasks: feature_selection, model_evaluation\n")
    
  } 
}
if (is.null(task)){
  cat("You should specify a task\n")
}
