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

# Calculate the proportion of NDM_1_gene
proportion_NDM_1_gene <- table(Bio_Data$NDM_1_gene) %>% prop.table()
cat("Proportion of NDM_1_gene:\n")
print(proportion_NDM_1_gene)


# Divide the dataset into train and test sets
set.seed(123)
split_strat  <- initial_split(Bio_Data, prop = 0.8, strata = "NDM_1_gene")
train_strat_NDM_1  <- training(split_strat)
test_strat_NDM_1 <- testing(split_strat)

# Display the dimensions of the train and test sets
cat("train set dimensions:", dim(train_strat_NDM_1), "\n")
cat("test set dimensions:", dim(test_strat_NDM_1), "\n")


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
      model_Log1_NDM_1 <- glm(NDM_1_gene ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                                Cefepime + Ceftazidime + Ciprofloxacin + Meropenem  
                              , family = "binomial", data = train_strat_NDM_1)
      
      cat("Summary of LR model based on Chi Squared test including Null deviance, Residual deviance, and AIC, etc.:\n")
      print(summary(model_Log1_NDM_1))
      
      test_strat_NDM_1$probs1_NDM_1 <- predict(model_Log1_NDM_1, test_strat_NDM_1, type = "response")
      
      test_strat_NDM_1$pred_logreg1_NDM_1 <- ifelse(test_strat_NDM_1$probs1_NDM_1 >= 0.5, 1, 0)
      
      # Model 1 evaluation
      confm_logreg1_NDM_1 <- table(actual = test_strat_NDM_1$NDM_1_gene, prediction = test_strat_NDM_1$pred_logreg1_NDM_1)
      
      metrics_logreg1_NDM_1 <- compute_metrics(confm_logreg1_NDM_1)
      
      cat("Performance of LR model based on Chi-Squared test results:\n")
      metrics_logreg1_NDM_1
    }
    else if (model == "LR-Wald test"){
      model_Log2_NDM_1 <- glm(NDM_1_gene ~ Imipenem + Amikacin
                              , family = "binomial", data = train_strat_NDM_1)
      cat("Summary of LR model based on Wald test including Null deviance, Residual deviance, and AIC, etc.:\n")
      print(summary(model_Log2_NDM_1))
      
      ##Prediction on test_strat_NDM_1 (Model 2)--------------------------------------------------
      test_strat_NDM_1$probs2_NDM_1 <- predict(model_Log2_NDM_1, test_strat_NDM_1, type = "response")
      
      test_strat_NDM_1$pred_logreg2_NDM_1 <- ifelse(test_strat_NDM_1$probs2_NDM_1 >= 0.5, 1, 0)
      
      # Model 2 evaluation
      confm_logreg2_NDM_1 <- table(actual = test_strat_NDM_1$NDM_1_gene, prediction = test_strat_NDM_1$pred_logreg2_NDM_1)
      
      metrics_logreg2_NDM_1 <- compute_metrics(confm_logreg2_NDM_1)
      
      cat("Performance of LR model based on Wald test results:\n")
      metrics_logreg2_NDM_1
    }
    else if (model == "LR-model agnostic"){
      model_Log3_NDM_1 <- glm(NDM_1_gene ~ Imipenem 
                              , family = "binomial", data = train_strat_NDM_1)
      
      cat("Summary of LR model based on model-agnostic approach including Null deviance, Residual deviance, and AIC, etc.:\n")
      print(summary(model_Log3_NDM_1))
      
      test_strat_NDM_1$probs3_NDM_1 <- predict(model_Log3_NDM_1, test_strat_NDM_1, type = "response")
      
      test_strat_NDM_1$pred_logreg3_NDM_1 <- ifelse(test_strat_NDM_1$probs3_NDM_1 >= 0.5, 1, 0)
      
      confm_logreg3_NDM_1 <- table(actual = test_strat_NDM_1$NDM_1_gene, prediction = test_strat_NDM_1$pred_logreg3_NDM_1)
      
      metrics_logreg3_NDM_1 <- compute_metrics(confm_logreg3_NDM_1)
      
      cat("Performance of LR model based on model-agnostic approach result:\n")
      metrics_logreg3_NDM_1
    }
    else if (model == "LR-Simplicity Principle_AK"){
      model_Log4_NDM_1 <- glm(NDM_1_gene ~ Amikacin 
                              , family = "binomial", data = train_strat_NDM_1)
      
      cat("Summary of LR model based on simplicity principle (AK) including Null deviance, Residual deviance, and AIC, etc.:\n")
      print(summary(model_Log4_NDM_1))
      
      ##Prediction on test_strat_NDM_1
      test_strat_NDM_1$probs4_NDM_1 <- predict(model_Log4_NDM_1, test_strat_NDM_1, type = "response")
      
      test_strat_NDM_1$pred_logreg4_NDM_1 <- ifelse(test_strat_NDM_1$probs4_NDM_1 >= 0.5, 1, 0)
      
      confm_logreg4_NDM_1 <- table(actual = test_strat_NDM_1$NDM_1_gene, prediction = test_strat_NDM_1$pred_logreg4_NDM_1)
      
      metrics_logreg4_NDM_1 <- compute_metrics(confm_logreg4_NDM_1)
      
      cat("Performance of LR model based on simplicity principle (AK):\n")
      metrics_logreg4_NDM_1
    }
    
    else if (model == "NBC-Chi Squared test"){
      train_control <- trainControl(
        method = "cv", 
        number = 10
      )
      
      search_grid_NDM_1 <- expand.grid(
        usekernel = c(FALSE,TRUE),
        fL = seq(0,5,0.5),
        adjust = seq(0,5,0.5)
      )
      
      
      model_nbdalex_NDM_1 <- train(NDM_1_gene ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                                     Cefepime + Ceftazidime + Ciprofloxacin + Meropenem,
                                   data = train_strat_NDM_1,
                                   method = "nb",
                                   metric = "Accuracy",
                                   trControl = train_control,
                                   tuneGrid = search_grid_NDM_1)
      
      
      cat("Tuned hyper-parameters through 10-fold Cross-Validation:\n")
      print(model_nbdalex_NDM_1$bestTune)
      
      ##Prediction on test_strat_NDM_1 (Naive Bayes Model)----------------------------------
      test_strat_NDM_1$pred_nb_NDM_1 <- predict(model_nbdalex_NDM_1, test_strat_NDM_1)
      
      # Model evaluation
      confm_nb_NDM_1 <- table(actual = test_strat_NDM_1$NDM_1_gene, prediction = test_strat_NDM_1$pred_nb_NDM_1)
      metrics_nb_NDM_1<- compute_metrics(confm_nb_NDM_1)
      cat("Performance of NBC model based on Chi-Squared test results:\n")
      metrics_nb_NDM_1
    }
    else if (model == "NBC-model agnostic"){
      train_control <- trainControl(
        method = "cv",
        number = 10
      )
      
      search_grid_NDM_1 <- expand.grid(
        usekernel = c(FALSE,TRUE),
        fL = seq(0,5,0.5),
        adjust = seq(0,5,0.5)
      )
      
      model_nbdalex1_NDM_1 <- train(NDM_1_gene ~ Imipenem + Cefepime + Ciprofloxacin,
                                    data = train_strat_NDM_1,
                                    method = "nb",
                                    metric = "Accuracy",
                                    trControl = train_control,
                                    tuneGrid = search_grid_NDM_1)
      
      cat("Tuned hyper-parameters through 10-fold Cross-Validation:\n")
      print(model_nbdalex1_NDM_1$bestTune)
      
      ##Prediction on test_strat_NDM_1 (Naive Bayes Model)----------------------------------
      test_strat_NDM_1$pred_nb1_NDM_1 <- predict(model_nbdalex1_NDM_1, test_strat_NDM_1)
      
      # Model evaluation
      confm_nb1_NDM_1 <- table(actual = test_strat_NDM_1$NDM_1_gene, prediction = test_strat_NDM_1$pred_nb1_NDM_1)
      
      metrics_nb1_NDM_1<- compute_metrics(confm_nb1_NDM_1)
      cat("Performance of NBC model based on model-agnostic approach results:\n")
      metrics_nb1_NDM_1
    }
    
    else if (model == "NBC-Wald test"){
      
      train_control <- trainControl(
        method = "cv",
        number = 10
      )
      
      search_grid_NDM_1 <- expand.grid(
        usekernel = c(FALSE,TRUE),
        fL = seq(0,5,0.5),
        adjust = seq(0,5,0.5)
      )
      
      model_nb2_NDM_1 <- train(NDM_1_gene ~ (Imipenem + Amikacin),
                               data = train_strat_NDM_1,
                               method = "nb",
                               metric = "Accuracy",
                               trControl = train_control,
                               tuneGrid = search_grid_NDM_1)
      
      cat("Tuned hyper-parameters through 10-fold Cross-Validation:\n")
      print(model_nb2_NDM_1$bestTune)
      
      ##Prediction on test_strat_NDM_1 (Naive Bayes Model)
      test_strat_NDM_1$pred_nb2_NDM_1 <- predict(model_nb2_NDM_1, test_strat_NDM_1)
      
      ##confusion matrix(Naive Bayes Model)
      confm_nb2_NDM_1 <- table(actual = test_strat_NDM_1$NDM_1_gene, prediction = test_strat_NDM_1$pred_nb2_NDM_1)
      
      metrics_nb2_NDM_1<- compute_metrics(confm_nb2_NDM_1)
      
      cat("Performance of NBC model based on Wald test results:\n")
      metrics_nb2_NDM_1
    }
    
    else if (model == "NBC-Simplicity Principle_IMI"){
      
      train_control <- trainControl(
        method = "cv",
        number = 10
      )
      
      search_grid_NDM_1 <- expand.grid(
        usekernel = c(FALSE,TRUE),
        fL = seq(0,5,0.5),
        adjust = seq(0,5,0.5)
      )
      
      model_nb3_NDM_1 <- train(NDM_1_gene ~ (Imipenem),
                               data = train_strat_NDM_1,
                               method = "nb",
                               metric = "Accuracy",
                               trControl = train_control,
                               tuneGrid = search_grid_NDM_1)
      
      cat("Tuned hyper-parameters through 10-fold Cross-Validation:\n")
      print(model_nb3_NDM_1$bestTune)
      
      ##Prediction on test_strat_NDM_1 (Naive Bayes Model)
      test_strat_NDM_1$pred_nb3_NDM_1 <- predict(model_nb3_NDM_1, test_strat_NDM_1)
      
      ##confusion matrix(Naive Bayes Model)
      confm_nb3_NDM_1 <- table(actual = test_strat_NDM_1$NDM_1_gene, prediction = test_strat_NDM_1$pred_nb3_NDM_1)
      
      metrics_nb3_NDM_1<- compute_metrics(confm_nb3_NDM_1)
      
      cat("Performance of NBC model based on Simplicity Principle (IMI):\n")
      metrics_nb3_NDM_1
    }
    
    else if (model == "NBC-Simplicity Principle_AK"){
      
      train_control <- trainControl(
        method = "cv",
        number = 10
      )
      
      search_grid_NDM_1 <- expand.grid(
        usekernel = c(FALSE,TRUE),
        fL = seq(0,5,0.5),
        adjust = seq(0,5,0.5)
      )
      
      model_nb4_NDM_1 <- train(NDM_1_gene ~ (Amikacin),
                               data = train_strat_NDM_1,
                               method = "nb",
                               metric = "Accuracy",
                               trControl = train_control,
                               tuneGrid = search_grid_NDM_1)
      
      cat("Tuned hyper-parameters through 10-fold Cross-Validation:\n")
      print(model_nb4_NDM_1$bestTune)
      
      ##Prediction on test_strat_NDM_1 (Naive Bayes Model)
      test_strat_NDM_1$pred_nb4_NDM_1 <- predict(model_nb4_NDM_1, test_strat_NDM_1)
      
      ##confusion matrix(Naive Bayes Model)
      confm_nb4_NDM_1 <- table(actual = test_strat_NDM_1$NDM_1_gene, prediction = test_strat_NDM_1$pred_nb4_NDM_1)
      
      metrics_nb4_NDM_1<- compute_metrics(confm_nb4_NDM_1)
      
      cat("Performance of NBC model based on Simplicity Principle (AK):\n")
      metrics_nb4_NDM_1
    }
    
    else if (model == "LDA-Chi Squared test"){
      
      model_ldadalex_NDM_1 <- train(NDM_1_gene ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                                      Cefepime + Ceftazidime + Ciprofloxacin + Meropenem  ,
                                    data = train_strat_NDM_1,
                                    method = "lda",
                                    metric = "Accuracy")
      
      cat("LDA model based on Chi-Squared test results:\n")
      print(model_ldadalex_NDM_1)
      
      #Prediction on test_strat_NDM_1 (LDA Model)-----------------------------------------------
      test_strat_NDM_1$pred_lda_NDM_1 <- predict(model_ldadalex_NDM_1, test_strat_NDM_1)
      
      # Model Evaluation-------------------------------------------------
      confm_lda_NDM_1 = table(actual = test_strat_NDM_1$NDM_1_gene, prediction = test_strat_NDM_1$pred_lda_NDM_1)
      metrics_lda_NDM_1<- compute_metrics(confm_lda_NDM_1)
      cat("Performance of LDA model based on Chi-Squared test results:\n")
      metrics_lda_NDM_1
    }
    else if (model == "LDA-model agnostic"){
      model_lda1_NDM_1 <- train(NDM_1_gene ~ Imipenem,
                                data = train_strat_NDM_1,
                                method = "lda",
                                metric = "Accuracy")
      
      cat("LDA model based on model-agnostic approach:\n")
      print(model_lda1_NDM_1)
      
      #Prediction on test_strat_NDM_1 (LDA Model)-----------------------------------------------
      test_strat_NDM_1$pred_lda1_NDM_1 <- predict(model_lda1_NDM_1, test_strat_NDM_1)
      
      # Model evaluation
      confm_lda1_NDM_1 = table(actual = test_strat_NDM_1$NDM_1_gene, prediction = test_strat_NDM_1$pred_lda1_NDM_1)
      metrics_lda1_NDM_1 <- compute_metrics(confm_lda1_NDM_1)
      cat("Performance of LDA model based on model-agnostic approach:\n")
      metrics_lda1_NDM_1
    }
    
    else if (model == "LDA-Wald test"){
      model_lda2_NDM_1 <- train(NDM_1_gene ~ Imipenem + Amikacin  ,
                                data = train_strat_NDM_1,
                                method = "lda",
                                metric = "Accuracy")
      
      cat("LDA model based on the Wald test results:\n")
      print(model_lda2_NDM_1)
      
      test_strat_NDM_1$pred_lda2_NDM_1 <- predict(model_lda2_NDM_1, test_strat_NDM_1)
      
      confm_lda2_NDM_1 = table(actual = test_strat_NDM_1$NDM_1_gene, prediction = test_strat_NDM_1$pred_lda2_NDM_1)
      
      metrics_lda2_NDM_1<- compute_metrics(confm_lda2_NDM_1)
      
      cat("Performance of LDA model based on the Wald test results:\n")
      metrics_lda2_NDM_1
    }
    
    else if (model == "LDA-Simplicity Principle_AK"){
      
      model_lda3_NDM_1 <- train(NDM_1_gene ~ Amikacin  ,
                                data = train_strat_NDM_1,
                                method = "lda",
                                metric = "Accuracy")
      
      cat("LDA model based on Simplicity Principle (AK):\n")
      print(model_lda3_NDM_1)
      
      test_strat_NDM_1$pred_lda3_NDM_1 <- predict(model_lda3_NDM_1, test_strat_NDM_1)
      
      confm_lda3_NDM_1 = table(actual = test_strat_NDM_1$NDM_1_gene, prediction = test_strat_NDM_1$pred_lda3_NDM_1)
      
      metrics_lda3_NDM_1<- compute_metrics(confm_lda3_NDM_1)
      
      cat("Performance of LDA model based on Simplicity Principle (AK):\n")
      metrics_lda3_NDM_1
    } 
    
    else if (model == "SVM-Sig-Chi Squared test"){
      #Fit Support Vector Machine model to train data set
      #This code performs a grid search for tuning the hyper-parameters of the SVM model using the tune function
      # from the e1071 package. It then fits the SVM model with the specified hyper-parameters using the svm function.
      # Model 1
      set.seed(1234)
      tune_out1_NDM_1 <- e1071::tune("svm", NDM_1_gene ~ Amikacin + Cefotaxime + Gentamicin + 
                                       Imipenem + Cefepime + Ceftazidime + Ciprofloxacin + Meropenem
                                     , data = train_strat_NDM_1, kernel = "sigmoid", tunecontrol=tune.control(cross=10),
                                     ranges = list(cost = seq(0.1, 5, length = 20)
                                                   , gamma = seq(0.1,5, length = 20)))
      
      cat("Tuned hyper-parameters throgh 10-Fold Cross Validation:\n")
      print(summary(tune_out1_NDM_1))
      
      model_svm1_NDM_1 = svm(NDM_1_gene ~ Amikacin + Cefotaxime + Gentamicin + 
                                  Imipenem + Cefepime + Ceftazidime + Ciprofloxacin + Meropenem, data = train_strat_NDM_1, 
                                  probability = T, kernel = "sigmoid", cost = tune_out1_NDM_1$best.parameters$cost,
                                  gamma = tune_out1_NDM_1$best.parameters$gamma )
      
      cat("Sigmoid SVM based on Chi-Squared test results:\n")
      print(model_svm1_NDM_1)
      
      #Prediction on test_strat_NDM_1 (SVM)----------------------------------------------------
      test_strat_NDM_1$pred_svm_NDM_1 <- predict(model_svm1_NDM_1, test_strat_NDM_1)
      
      # Model evaluation
      confm_svm_NDM_1 <- table(actual = test_strat_NDM_1$NDM_1_gene, prediction = test_strat_NDM_1$pred_svm_NDM_1)
      
      metrics_svm1_NDM_1_sig<- compute_metrics(confm_svm_NDM_1)
      
      cat("Performance of SVM-Sigmoid model based on Chi-Squared test results:\n")
      metrics_svm1_NDM_1_sig
    }
    else if (model=="SVM-Sig-model agnostic"){
      
      set.seed(1234)
      tune_out2_NDM_1 <- e1071::tune("svm", NDM_1_gene ~ Amikacin,
                                     data = train_strat_NDM_1,
                                     kernel = "sigmoid", tunecontrol=tune.control(cross=10),
                                     ranges = list(cost = seq(0.1, 5, length = 25)
                                                   , gamma = seq(0, 5, length = 25)))
      cat("Tuned hyper-parameters throgh 10-Fold Cross Validation:\n")
      print(summary(tune_out2_NDM_1))
      
      model_svm2_NDM_1 = svm(NDM_1_gene ~ Amikacin, data = train_strat_NDM_1, 
                                   probability = T, kernel = "sigmoid", cost = tune_out2_NDM_1$best.parameters$cost,
                                   gamma = tune_out2_NDM_1$best.parameters$gamma)
      
      cat("Sigmoid SVM based on model-agnostic approach results:\n")
      print(model_svm2_NDM_1)
      
      test_strat_NDM_1$pred_svm2_NDM_1 <- predict(model_svm2_NDM_1, test_strat_NDM_1)
      
      confm_svm2_NDM_1 <- table(actual = test_strat_NDM_1$NDM_1_gene, prediction = test_strat_NDM_1$pred_svm2_NDM_1)
      
      metrics_svm2_NDM_1_sig<- compute_metrics(confm_svm2_NDM_1)
      
      cat("Performance of SVM-Sigmoid model based on model-agnostic approach results:\n")
      metrics_svm2_NDM_1_sig
    }
    else if (model=="SVM-Sig-Wald test"){
      set.seed(1234)
      tune_out3_NDM_1 <- e1071::tune("svm", NDM_1_gene ~ Imipenem + Amikacin
                                     , data = train_strat_NDM_1, kernel = "sigmoid", tunecontrol=tune.control(cross=10),
                                     ranges = list(cost = seq(0.1, 5, length = 25)
                                                   , gamma = seq(0.1, 1, length = 25)))
      
      cat("Tuned hyper-parameters throgh 10-Fold Cross Validation:\n")
      print(summary(tune_out3_NDM_1))
      
      model_svm3_NDM_1 = svm(NDM_1_gene ~ Imipenem + Amikacin, data = train_strat_NDM_1, 
                                    probability = T, kernel = "sigmoid", cost = tune_out3_NDM_1$best.parameters$cost,
                                    gamma = tune_out3_NDM_1$best.parameters$gamma )
      
      cat("Sigmoid SVM based on the Wald test result:\n")
      print(model_svm3_NDM_1)
      
      #Prediction on test_strat_NDM_1 (SVM)
      test_strat_NDM_1$pred_svm3_NDM_1 <- predict(model_svm3_NDM_1, test_strat_NDM_1)
      
      # Model evaluation
      confm_svm3_NDM_1 <- table(actual = test_strat_NDM_1$NDM_1_gene, prediction = test_strat_NDM_1$pred_svm3_NDM_1)
      
      metrics_svm3_NDM_1_sig<- compute_metrics(confm_svm3_NDM_1)
      
      cat("Performance of SVM-Sigmoid model based on the Wald test results:\n")
      metrics_svm3_NDM_1_sig
    } 
    
    else if (model=="SVM-Sig-Simplicity Principle_IMI"){
      set.seed(1234)
      tune_out4_NDM_1 <- e1071::tune("svm", NDM_1_gene ~ Imipenem
                                     , data = train_strat_NDM_1, kernel = "sigmoid", tunecontrol=tune.control(cross=10),
                                     ranges = list(cost = seq(0.1, 5, length = 25)
                                                   , gamma = seq(0.1, 5, length = 25)))
      
      cat("Tuned hyper-parameters throgh 10-Fold Cross Validation:\n")
      print(summary(tune_out4_NDM_1))
      
      model_svm4_NDM_1 = svm(NDM_1_gene ~ Imipenem, data = train_strat_NDM_1, 
                                    probability = T, kernel = "sigmoid", cost = tune_out4_NDM_1$best.parameters$cost,
                                    gamma = tune_out4_NDM_1$best.parameters$gamma )
      
      cat("Sigmoid SVM based on Simplicity Principle (IMI):\n")
      print(model_svm4_NDM_1)
      
      #Prediction on test_strat_NDM_1 (SVM)
      test_strat_NDM_1$pred_svm4_NDM_1 <- predict(model_svm4_NDM_1, test_strat_NDM_1)
      
      # Model evaluation
      confm_svm4_NDM_1 <- table(actual = test_strat_NDM_1$NDM_1_gene, prediction = test_strat_NDM_1$pred_svm4_NDM_1)
      
      metrics_svm4_NDM_1_sig<- compute_metrics(confm_svm4_NDM_1)
      
      cat("Performance of SVM-Sigmoid model based on Simplicity Principle (IMI):\n")
      metrics_svm4_NDM_1_sig
    } 
    
    else if (model=="SVM-Poly-Chi Squared test"){
      set.seed(1234)
      tune_out5_NDM_1 <- e1071::tune("svm", NDM_1_gene ~ Amikacin + Cefotaxime + Gentamicin + 
                                       Imipenem + Cefepime + Ceftazidime + Ciprofloxacin + Meropenem  
                                     , data = train_strat_NDM_1, kernel = "polynomial", tunecontrol=tune.control(cross=10),
                                     ranges = list(cost = seq(0.1, 5, length = 20), degree = seq(1,3,1)
                                                   , gamma = seq(0.1,5, length = 20)))
      
      cat("Tuned hyper-parameters throgh 10-Fold Cross Validation:\n")
      print(summary(tune_out5_NDM_1))
      
      # Tuned hyper-parameters are placed in the model.
      model_svm5_NDM_1 = svm(NDM_1_gene ~ Amikacin + Cefotaxime + Gentamicin + 
                              Imipenem + Cefepime + Ceftazidime + Ciprofloxacin + Meropenem
                            , data = train_strat_NDM_1, probability = T, 
                              kernel = "polynomial", cost = tune_out5_NDM_1$best.parameters$cost, 
                              degree = tune_out5_NDM_1$best.parameters$degree , gamma = tune_out5_NDM_1$best.parameters$gamma )
      
      cat("Polynomial SVM based on Chi-Squared test results:\n")
      print(model_svm5_NDM_1)
      
      #Prediction on test (SVM 2)----------------------------------------------------
      test_strat_NDM_1$pred_svm5_NDM_1 <- predict(model_svm5_NDM_1, test_strat_NDM_1)
      
      # Model evaluation
      confm_svm5_NDM_1 <- table(actual = test_strat_NDM_1$NDM_1_gene, prediction = test_strat_NDM_1$pred_svm5_NDM_1)
      
      metrics_svm_NDM_1<- compute_metrics(confm_svm5_NDM_1)
      
      cat("Performance of SVM-Polynomial model based on Chi-Squared test results:\n")
      metrics_svm_NDM_1
    }
    else if (model=="SVM-Poly-model agnostic"){
      
      set.seed(1234)
      tune_out6_NDM_1 <- e1071::tune("svm", NDM_1_gene ~ Imipenem + Amikacin,
                                     data = train_strat_NDM_1,
                                     kernel = "polynomial", tunecontrol=tune.control(cross=10),
                                     ranges = list(cost = seq(0.1, 5, length = 20), degree = seq(1,3,1)
                                                   , gamma = seq(0, 5, length = 20)))
      
      cat("Tuned hyper-parameters throgh 10-Fold Cross Validation:\n")
      print(summary(tune_out6_NDM_1))
      
      model_svm6_NDM_1 = svm(NDM_1_gene ~ Imipenem + Amikacin, data = train_strat_NDM_1, 
                               probability = T, kernel = "polynomial", cost = tune_out6_NDM_1$best.parameters$cost,
                               degree = tune_out6_NDM_1$best.parameters$degree , gamma = tune_out6_NDM_1$best.parameters$gamma )
      
      cat("Polynomial SVM based on model-agnostic approach / Wald test results:\n")
      print(model_svm6_NDM_1)
      
      #Prediction on test_strat_NDM_1 (SVM)----------------------------------------------------
      test_strat_NDM_1$pred_svm6_NDM_1 <- predict(model_svm6_NDM_1, test_strat_NDM_1)
      
      # Model evaluation
      confm_svm6_NDM_1 <- table(actual = test_strat_NDM_1$NDM_1_gene, prediction = test_strat_NDM_1$pred_svm6_NDM_1)
      
      metrics_svm1_NDM_1<- compute_metrics(confm_svm6_NDM_1)
      
      cat("Performance of SVM-Polynomial model based on model-agnostic approach / Wald test results:\n")
      metrics_svm1_NDM_1
      
    }
    
    else if (model=="SVM-Poly-Simplicity principle_IMI"){
      
      set.seed(1234)
      tune_out7_NDM_1 <- e1071::tune("svm", NDM_1_gene ~ Imipenem,
                                     data = train_strat_NDM_1,
                                     kernel = "polynomial", tunecontrol=tune.control(cross=10),
                                     ranges = list(cost = seq(0.1, 5, length = 20), degree = seq(1,3,1)
                                                   , gamma = seq(0.1, 5, length = 20)))
      
      cat("Tuned hyper-parameters throgh 10-Fold Cross Validation:\n")
      print(summary(tune_out7_NDM_1))
      
      model_svm7_NDM_1 = svm(NDM_1_gene ~ Imipenem , data = train_strat_NDM_1, 
                                probability = T, kernel = "polynomial", cost = tune_out7_NDM_1$best.parameters$cost,
                                degree = tune_out7_NDM_1$best.parameters$degree , gamma = tune_out7_NDM_1$best.parameters$gamma )
      
      cat("Polynomial SVM based on the Simplicity Principle (IMI):\n")
      print(model_svm7_NDM_1)
      
      #Prediction on test_strat_NDM_1 (SVM)----------------------------------------------------
      test_strat_NDM_1$pred_svm7_NDM_1 <- predict(model_svm7_NDM_1, test_strat_NDM_1)
      
      # Model evaluation
      confm_svm7_NDM_1 <- table(actual = test_strat_NDM_1$NDM_1_gene, prediction = test_strat_NDM_1$pred_svm7_NDM_1)
      
      metrics_svm3_NDM_1<- compute_metrics(confm_svm7_NDM_1)
      
      cat("Performance of SVM-Polynomial model based on the Simplicity Principle (IMI):\n")
      metrics_svm3_NDM_1
    }
    
    else if (model=="SVM-Poly-Simplicity principle_AK"){
      
      set.seed(1234)
      tune_out8_NDM_1 <- e1071::tune("svm", NDM_1_gene ~ Amikacin,
                                     data = train_strat_NDM_1,
                                     kernel = "polynomial", tunecontrol=tune.control(cross=10),
                                     ranges = list(cost = seq(0.1, 5, length = 20), degree = seq(1,3,1)
                                                   , gamma = seq(0.1, 5, length = 20)))
      
      cat("Tuned hyper-parameters throgh 10-Fold Cross Validation:\n")
      print(summary(tune_out8_NDM_1))
      
      model_svm8_NDM_1 = svm(NDM_1_gene ~ Amikacin , data = train_strat_NDM_1, 
                                probability = T, kernel = "polynomial", cost = tune_out8_NDM_1$best.parameters$cost,
                                degree = tune_out8_NDM_1$best.parameters$degree , gamma = tune_out8_NDM_1$best.parameters$gamma )
      
      cat("Polynomial SVM based on the Simplicity Principle (AK):\n")
      print(model_svm8_NDM_1)
      
      #Prediction on test_strat_NDM_1 (SVM)----------------------------------------------------
      test_strat_NDM_1$pred_svm8_NDM_1 <- predict(model_svm8_NDM_1, test_strat_NDM_1)
      
      # Model evaluation
      confm_svm8_NDM_1 <- table(actual = test_strat_NDM_1$NDM_1_gene, prediction = test_strat_NDM_1$pred_svm8_NDM_1)
      
      metrics_svm4_NDM_1<- compute_metrics(confm_svm8_NDM_1)
      
      cat("Performance of SVM-Polynomial model based on the Simplicity Principle (AK):\n")
      metrics_svm4_NDM_1
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
      
      split_strat_NDM_1  <- initial_split(Bio_Data, prop = 0.8, strata = "NDM_1_gene")
      
      train_strat_NDM_1  <- training(split_strat_NDM_1)
      
      test_strat_NDM_1   <- testing(split_strat_NDM_1)
      
      grid_NDM_1 <- expand.grid(
        depth = c(2,3,4,5),
        learning_rate = c(0.1,0.2,0.3),
        l2_leaf_reg = c(2,3,4),
        rsm = 1,
        border_count = 1,
        iterations = seq(10,100,5)
      )
      fitControl_NDM_1 <- trainControl(method = "cv",
                                 number = 10,
                                 classProbs = TRUE)
      
      set.seed(1234)
      model_catboost_NDM_1 <- train(x = train_strat_NDM_1[,c("Amikacin", "Cefotaxime","Imipenem", "Gentamicin",
                                                             "Cefepime", "Ceftazidime", "Ciprofloxacin", "Meropenem")],
                                    y = make.names(train_strat_NDM_1[,"NDM_1_gene"]),
                                    maximize = TRUE,
                                    method = catboost.caret, metric = "Accuracy", 
                                    tuneGrid =  grid_NDM_1,
                                    trControl = fitControl_NDM_1)
      
      cat("Tuned hyper-parameters throgh 10-Fold Cross Validation:\n")
      print(model_catboost_NDM_1$bestTune)  
      
      #Prediction on test_strat_NDM_1 (CatBoost)----------------------------------------------------
      test_strat_NDM_1$cat_NDM0_1 = predict(model_catboost_NDM_1, test_strat_NDM_1)
      
      # Model evaluation
      confm_cat0_NDM_1 <- table(actual = test_strat_NDM_1$NDM_1_gene, prediction = test_strat_NDM_1$cat_NDM0_1)
      
      metrics_cat_NDM_1<- compute_metrics(confm_cat0_NDM_1)
      
      cat("Performance of CatBoost model based on Chi-Squared test results:\n")
      metrics_cat_NDM_1
    }
    
    else if (model=="CatBoost-model agnostic"){
      # CatBoost Models
      Bio_Data <- as.data.frame(Bio_Data)
      
      #Convert categorical variables to factor
      categorical_var <- c(colnames(Bio_Data))
      
      Bio_Data[, categorical_var[-1]] <- lapply(Bio_Data[, categorical_var[-1]], factor)
      
      levels(Bio_Data$Ampicillin) <- c(levels(Bio_Data$Ampicillin),0)
      
      #Divide Data set into train and test_strat_NDM_1
      set.seed(123)
      
      split_strat_NDM_1  <- initial_split(Bio_Data, prop = 0.8, strata = "NDM_1_gene")
      
      train_strat_NDM_1  <- training(split_strat_NDM_1)
      
      test_strat_NDM_1   <- testing(split_strat_NDM_1)
      
      grid_NDM_1 <- expand.grid(
        depth = c(2,3,4),
        learning_rate = c(0.1,0.2,0.3),
        l2_leaf_reg = c(1,2,3),
        rsm = 1,
        border_count = 1,
        iterations = seq(10,100,10)
      )
      fitControl_NDM_1 <- trainControl(method = "cv",
                                       number = 10,
                                       classProbs = TRUE)
      
      model_catboost1_NDM_1 <- train(x = train_strat_NDM_1["Imipenem"],
                                     y = make.names(train_strat_NDM_1[["NDM_1_gene"]]),
                                     maximize = TRUE,
                                     method = catboost.caret, metric = "Accuracy", 
                                     tuneGrid =  grid_NDM_1,
                                     trControl = fitControl_NDM_1)
      
      cat("Tuned hyper-parameters throgh 10-Fold Cross Validation:\n")
      print(model_catboost1_NDM_1$bestTune)  
      
      test_strat_NDM_1$cat_NDM_1 = predict(model_catboost1_NDM_1, test_strat_NDM_1)
      
      # Model Evaluation
      confm_cat_NDM_1 <- table(actual = test_strat_NDM_1$NDM_1_gene, prediction = test_strat_NDM_1$cat_NDM_1)
      
      metrics_cat1_NDM_1<- compute_metrics(confm_cat_NDM_1)
      
      cat("Performance of CatBoost model based on model-agnostic approach results:\n")
      metrics_cat1_NDM_1
    }
    
    else if (model=="CatBoost-Wald test"){
      # CatBoost Models
      Bio_Data <- as.data.frame(Bio_Data)
      
      #Convert categorical variables to factor
      categorical_var <- c(colnames(Bio_Data))
      
      Bio_Data[, categorical_var[-1]] <- lapply(Bio_Data[, categorical_var[-1]], factor)
      
      levels(Bio_Data$Ampicillin) <- c(levels(Bio_Data$Ampicillin),0)
      
      #Divide Data set into train and test_strat_NDM_1
      set.seed(123)
      
      split_strat_NDM_1  <- initial_split(Bio_Data, prop = 0.8, strata = "NDM_1_gene")
      
      train_strat_NDM_1  <- training(split_strat_NDM_1)
      
      test_strat_NDM_1   <- testing(split_strat_NDM_1)
      
      grid3_NDM_1 <- expand.grid(
        depth = c(2,3,4),
        learning_rate = c(0.1,0.2,0.3),
        l2_leaf_reg = c(1,2,3),
        rsm = 1,
        border_count = 1,
        iterations = seq(10,30,1)
      )
      
      fitControl3_NDM_1 <- trainControl(method = "cv",
                                        number = 10,
                                        classProbs = TRUE)
      
      model_catboost3_NDM_1 <- train(x = train_strat_NDM_1[,c("Amikacin","Imipenem")],
                                     y = make.names(train_strat_NDM_1[["NDM_1_gene"]]),
                                     maximize = TRUE,
                                     method = catboost.caret, metric = "Accuracy", 
                                     tuneGrid =  grid3_NDM_1,
                                     trControl = fitControl3_NDM_1)
      
      cat("Tuned hyper-parameters throgh 10-Fold Cross Validation:\n")
      print(model_catboost3_NDM_1$bestTune)  
      
      test_strat_NDM_1$cat_NDM3_1 = predict(model_catboost3_NDM_1, test_strat_NDM_1)
      # Model evaluation
      confm_cat3_NDM_1 <- table(actual = test_strat_NDM_1$NDM_1_gene, prediction = test_strat_NDM_1$cat_NDM3_1)
      
      metrics_cat2_NDM_1<- compute_metrics(confm_cat3_NDM_1)
      
      cat("Performance of CatBoost model based on the Wald test results:\n")
      metrics_cat2_NDM_1
    }
    
    else if (model=="CatBoost-Simplicity principle_AK"){
      # CatBoost Models
      Bio_Data <- as.data.frame(Bio_Data)
      
      #Convert categorical variables to factor
      categorical_var <- c(colnames(Bio_Data))
      
      Bio_Data[, categorical_var[-1]] <- lapply(Bio_Data[, categorical_var[-1]], factor)
      
      levels(Bio_Data$Ampicillin) <- c(levels(Bio_Data$Ampicillin),0)
      
      #Divide Data set into train and test_strat_NDM_1
      set.seed(123)
      
      split_strat_NDM_1  <- initial_split(Bio_Data, prop = 0.8, strata = "NDM_1_gene")
      
      train_strat_NDM_1  <- training(split_strat_NDM_1)
      
      test_strat_NDM_1   <- testing(split_strat_NDM_1)
      
      grid2_NDM_1 <- expand.grid(
        depth = c(2,3,4),
        learning_rate = c(0.1,0.2,0.3),
        l2_leaf_reg = c(1,2,3),
        rsm = 1,
        border_count = 1,
        iterations = seq(10,100,10)
      )
      
      fitControl2_NDM_1 <- trainControl(method = "cv",
                                        number = 10,
                                        classProbs = TRUE)
      
      model_catboost2_NDM_1 <- train(x = train_strat_NDM_1["Amikacin"],
                                     y = make.names(train_strat_NDM_1[["NDM_1_gene"]]),
                                     maximize = TRUE,
                                     method = catboost.caret, metric = "Accuracy", 
                                     tuneGrid =  grid2_NDM_1,
                                     trControl = fitControl2_NDM_1)
      
      cat("Tuned hyper-parameters throgh 10-Fold Cross Validation:\n")
      print(model_catboost2_NDM_1$bestTune)  
      
      test_strat_NDM_1$cat_NDM2_1 = predict(model_catboost2_NDM_1, test_strat_NDM_1)
      # Model evaluation
      confm_cat2_NDM_1 <- table(actual = test_strat_NDM_1$NDM_1_gene, prediction = test_strat_NDM_1$cat_NDM2_1)
      
      metrics_cat3_NDM_1<- compute_metrics(confm_cat2_NDM_1)
      
      cat("Performance of CatBoost model based on the simplicity principle (AK):\n")
      metrics_cat3_NDM_1
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
        result <- chisq.test(table(data[[variable]], data$NDM_1_gene))
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
    model_Log1_NDM_1 <- glm(NDM_1_gene ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                              Cefepime + Ceftazidime + Ciprofloxacin + Meropenem  
                            , family = "binomial", data = train_strat_NDM_1)
    
    cat("P-values show Wald test results:\n")
    print(summary(model_Log1_NDM_1))
  }
  else if (method == "model agnostic-LR"){
    # Feature Selection by DALEX Package--------------------------------------
    # the 'explained_glm_TEM' object is created using the 'explain' function from the DALEX package. 
    # The variable_importance function is used to calculate the permutation-based 
    # feature importance ('fi_glm_TEM') with 50 permutations. Finally, the feature importance is displayed 
    # using 'fi_glm_TEM' and plotted using' plot(fi_glm_TEM)'.
    model_Log1_NDM_1 <- glm(NDM_1_gene ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                                Cefepime + Ceftazidime + Ciprofloxacin +
                                Meropenem + Ceftriaxone + Aztreonam 
                              , family = "binomial", data = train_strat_NDM_1)
    
    explained_glm_NDM_1 <- explain(model = model_Log1_NDM_1, data=train_strat_NDM_1[,c(2:9)], variables = colnames(train_strat_NDM_1[,c(2:9)]),
                                   y=as.vector(as.numeric(train_strat_NDM_1$NDM_1_gene))-1, label = "LR", 
                                   type = "classification")
    #50 permuatation
    fi_glm_NDM_1 = variable_importance(explained_glm_NDM_1, B=50 ,variables = colnames(train_strat_NDM_1[,c(2:9)]),loss_function = loss_root_mean_square, type = "raw")
    
    cat("Root mean squared dropout loss due to a feature elimination through model-agnostic approach in LR:\n")
    print(fi_glm_NDM_1)
  }
  else if (method == "model agnostic-NBC"){
    train_control <- trainControl(
      method = "cv",
      number = 10
    )
    
    search_grid_NDM_1 <- expand.grid(
      usekernel = c(FALSE,TRUE),
      fL = seq(0,5,0.5),
      adjust = seq(0,5,0.5)
    )
    
    model_nbdalex_NDM_1 <- train(NDM_1_gene ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                                   Cefepime + Ceftazidime + Ciprofloxacin + Meropenem,
                                 data = train_strat_NDM_1,
                                 method = "nb",
                                 metric = "Accuracy",
                                 trControl = train_control,
                                 tuneGrid = search_grid_NDM_1)
    # Feature selection through model-agnostic approach (DALEX)
    # the 'explain' function from the DALEX package is used to explain the Naive Bayes classifier model 
    # 'model_nbdalex_TEM'. The relevant variables are selected and provided as data along with their 
    # corresponding column names. The target variable 'NDM_1_gene' is transformed to numeric values 
    # (y = as.vector(as.numeric(train_strat_NDM_1$NDM_1_gene)) - 1). The label is set to "NBC" for Naive Bayes classifier,
    # and the type is set to "classification".
    expnb_NDM_1 = explain(model = model_nbdalex_NDM_1, data=train_strat_NDM_1[,c("Amikacin", "Cefotaxime","Imipenem", "Gentamicin",
                                                                                 "Cefepime", "Ceftazidime", "Ciprofloxacin", "Meropenem")], variables = colnames(train_strat_NDM_1[,c("Amikacin", "Cefotaxime","Imipenem", "Gentamicin",
                                                                                                                                                                                      "Cefepime", "Ceftazidime", "Ciprofloxacin", "Meropenem")]),
                          y=as.vector(as.numeric(train_strat_NDM_1$NDM_1_gene))-1, label = "NBC", 
                          type = "classification")
    
    # The 'variable_importance' function is then used to calculate the variable importance based on the
    # explained model. The 'B' parameter is set to 50 for the number of permutations, and the loss function
    # is set to "loss_root_mean_square". The resulting variable importance is stored in 'fitnb_TEM'
    fitnb_NDM_1 = variable_importance(expnb_NDM_1, B=50 ,variables = colnames(train_strat_NDM_1[,c("Amikacin", "Cefotaxime","Imipenem", "Gentamicin",
                                                                                                   "Cefepime", "Ceftazidime", "Ciprofloxacin", "Meropenem")]), loss_function = loss_root_mean_square, type = "raw" )
    
    cat("Root mean squared dropout loss due to a feature elimination through model-agnostic approach in NBC:\n")
    print(fitnb_NDM_1)
  }
  else if (method == "model agnostic-LDA"){
    model_ldadalex_NDM_1 <- train(NDM_1_gene ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                                    Cefepime + Ceftazidime + Ciprofloxacin + Meropenem  ,
                                  data = train_strat_NDM_1,
                                  method = "lda",
                                  metric = "Accuracy")
    
    explda_NDM_1 = explain(model = model_ldadalex_NDM_1, data=train_strat_NDM_1[,c("Amikacin", "Cefotaxime","Imipenem", "Gentamicin",
                                                                                   "Cefepime", "Ceftazidime", "Ciprofloxacin", "Meropenem")], variables = colnames(train_strat_NDM_1[,c("Amikacin", "Cefotaxime","Imipenem", "Gentamicin",
                                                                                                                                                                                        "Cefepime", "Ceftazidime", "Ciprofloxacin", "Meropenem")]),
                           y=as.vector(as.numeric(train_strat_NDM_1$NDM_1_gene))-1, label = "LDA", 
                           type = "classification")
    
    fitlda_NDM_1 = variable_importance(explda_NDM_1, B=50 ,variables = colnames(train_strat_NDM_1[,c("Amikacin", "Cefotaxime","Imipenem", "Gentamicin",
                                                                                                     "Cefepime", "Ceftazidime", "Ciprofloxacin", "Meropenem")]), loss_function = loss_root_mean_square, type = "raw" )
    cat("Root mean squared dropout loss due to a feature elimination through model-agnostic approach in LDA:\n")
    print(fitlda_NDM_1)
  }
  else if (method == "model agnostic-Sig SVM"){
    set.seed(1234)
    tune_out1_NDM_1 <- e1071::tune("svm", NDM_1_gene ~ Amikacin + Cefotaxime + Gentamicin + 
                                     Imipenem + Cefepime + Ceftazidime + Ciprofloxacin + Meropenem
                                   , data = train_strat_NDM_1, kernel = "sigmoid", tunecontrol=tune.control(cross=10),
                                   ranges = list(cost = seq(0.1, 5, length = 20)
                                                 , gamma = seq(0.1,5, length = 20)))
    
    model_svm1_NDM_1 = svm(NDM_1_gene ~ Amikacin + Cefotaxime + Gentamicin + 
                             Imipenem + Cefepime + Ceftazidime + Ciprofloxacin + Meropenem, data = train_strat_NDM_1, 
                           probability = T, kernel = "sigmoid", cost = tune_out1_NDM_1$best.parameters$cost,
                           gamma = tune_out1_NDM_1$best.parameters$gamma )
    
    explained_SVM1_NDM_1 <- explain(model = model_svm1_NDM_1, data=train_strat_NDM_1[,c(2:9)], variables = colnames(train_strat_NDM_1[,c(2:9)]),
                                    y=as.vector(as.numeric(train_strat_NDM_1$NDM_1_gene))-1, label = "SVM-Sigmoid", 
                                    type = "classification")
    
    #50 permuatation
    fi_SVM1_NDM_1 = variable_importance(explained_SVM1_NDM_1, B=50 ,variables = colnames(train_strat_NDM_1[,c(2:9)]),loss_function = loss_root_mean_square, type = "raw")
    
    cat("Root mean squared dropout loss due to a feature elimination through model-agnostic approach in sigmoid kernel of SVM:\n")
    print(fi_SVM1_NDM_1)
  }
  else if (method == "model agnostic-Poly SVM"){
    # it is highly recommeded that first train the model, then specify following hyper-parameters according to the tuned hyper-parameters.
    # Following values of the hyper-parameters have been specified in this way.
    set.seed(1234)
    tune_out5_NDM_1 <- e1071::tune("svm", NDM_1_gene ~ Amikacin + Cefotaxime + Gentamicin + 
                                     Imipenem + Cefepime + Ceftazidime + Ciprofloxacin + Meropenem  
                                   , data = train_strat_NDM_1, kernel = "polynomial", tunecontrol=tune.control(cross=10),
                                   ranges = list(cost = seq(0.1, 5, length = 20), degree = seq(1,3,1)
                                                 , gamma = seq(0.1,5, length = 20)))
    
    # Tuned hyper-parameters are placed in the model.
    model_svm5_NDM_1 = svm(NDM_1_gene ~ Amikacin + Cefotaxime + Gentamicin + 
                             Imipenem + Cefepime + Ceftazidime + Ciprofloxacin + Meropenem
                           , data = train_strat_NDM_1, probability = T, 
                           kernel = "polynomial", cost = tune_out5_NDM_1$best.parameters$cost, 
                           degree = tune_out5_NDM_1$best.parameters$degree , gamma = tune_out5_NDM_1$best.parameters$gamma )
    
    explained_SVM5_NDM_1 <- explain(model = model_svm5_NDM_1, data=train_strat_NDM_1[,c(2:9)], variables = colnames(train_strat_NDM_1[,c(2:9)]),
                                    y=as.vector(as.numeric(train_strat_NDM_1$NDM_1_gene))-1, label = "SVM-Polynomial", 
                                    type = "classification")
    
    #50 permuatation
    fi_SVM5_NDM_1 = variable_importance(explained_SVM5_NDM_1, B=50 ,variables = colnames(train_strat_NDM_1[,c(2:9)]),loss_function = loss_root_mean_square, type = "raw")
    
    cat("Root mean squared dropout loss due to a feature elimination through model-agnostic approach in polynomial kernel of SVM:\n")
    print(fi_SVM5_NDM_1)
  }
  else if (method == "model agnostic-CatBoost"){

    Bio_Data <- as.data.frame(Bio_Data)
    
    #Convert categorical variables to factor
    categorical_var <- c(colnames(Bio_Data))
    
    Bio_Data[, categorical_var[-1]] <- lapply(Bio_Data[, categorical_var[-1]], factor)
    
    levels(Bio_Data$Ampicillin) <- c(levels(Bio_Data$Ampicillin),0)
    
    #Divide Data set into train and test
    set.seed(123)
    
    split_strat_NDM_1  <- initial_split(Bio_Data, prop = 0.8, strata = "NDM_1_gene")
    
    train_strat_NDM_1  <- training(split_strat_NDM_1)
    
    test_strat_NDM_1   <- testing(split_strat_NDM_1)
    
    grid_NDM_1 <- expand.grid(
      depth = c(2,3,4,5),
      learning_rate = c(0.1,0.2,0.3),
      l2_leaf_reg = c(2,3,4),
      rsm = 1,
      border_count = 1,
      iterations = seq(10,100,5)
    )
    fitControl_NDM_1 <- train_strat_NDM_1Control(method = "cv",
                                                 number = 10,
                                                 classProbs = TRUE)
    
    set.seed(1234)
    model_catboost_NDM_1 <- train(x = train_strat_NDM_1[,c("Amikacin", "Cefotaxime","Imipenem", "Gentamicin",
                                                           "Cefepime", "Ceftazidime", "Ciprofloxacin", "Meropenem")],
                                  y = make.names(train_strat_NDM_1[,"NDM_1_gene"]),
                                  maximize = TRUE,
                                  method = catboost.caret, metric = "Accuracy", 
                                  tuneGrid =  grid_NDM_1,
                                  trControl = fitControl_NDM_1)
    
    expsda_NDM_1 = explain(model = model_catboost_NDM_1, data=train_strat_NDM_1[,c("Amikacin", "Cefotaxime","Imipenem", "Gentamicin",
                                                                                   "Cefepime", "Ceftazidime", "Ciprofloxacin", "Meropenem")], variables = colnames(train_strat_NDM_1[,c("Amikacin", "Cefotaxime","Imipenem", "Gentamicin",
                                                                                                                                                                                        "Cefepime", "Ceftazidime", "Ciprofloxacin", "Meropenem")]),
                           y=as.vector(as.numeric(train_strat_NDM_1$NDM_1_gene))-1, label = "CatBoost", 
                           type = "classification")
    
    fitcat_NDM_1 = variable_importance(expsda_NDM_1, B=50 ,variables = colnames(train_strat_NDM_1[,c("Amikacin", "Cefotaxime","Imipenem", "Gentamicin",
                                                                                                     "Cefepime", "Ceftazidime", "Ciprofloxacin", "Meropenem")]), loss_function = loss_root_mean_square, type = "raw" )
    
    cat("Root mean squared dropout loss due to a feature elimination through model-agnostic approach in CatBoost:\n")
    print(fitcat_NDM_1)
  }
  else {
    cat("Invalid task. Please specify one of the following tasks: feature_selection, model_evaluation\n")
    
  } 
}
if (is.null(task)){
  cat("You should specify a task\n")
}

