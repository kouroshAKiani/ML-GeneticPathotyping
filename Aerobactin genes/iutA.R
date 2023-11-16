source('F1-score.R') # Call F1-Score function for computing point estimation and its 95% confidence interval

source('compute_metrics.R') # Call compute_metrics function for computing all classification metrics and their confidence interval.

source('Packages_installation.R') # Installation of all required packages

# Parse command argumnets 
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

# Calculate the proportion of IUT_gene

proportion_IUT <- table(Bio_Data$IUT) %>% prop.table()

cat("Proportion of IUT gene:\n")

print(proportion_IUT)

# We will use Over and under sampling combination for data splitting

train_strat_IUT <- ovun.sample(IUT~Imipenem + Gentamicin + Ciprofloxacin + Meropenem + Aztreonam +
                                 Cefotaxime + Cefepime + Ceftazidime + Ceftriaxone, data=train, 
                               method = "both", p = 0.5, seed = 1234, N = 272)$data

table(train_strat_IUT$IUT) %>% prop.table() 

set.seed(123)

train_cases <- sample(1:nrow(Bio_Data), nrow(Bio_Data) * 0.8)

test <- Bio_Data[-train_cases,]

# Display the dimensions of the train and test sets
cat("train set dimensions:", dim(train_strat_IUT), "\n")
cat("test set dimensions:", dim(test), "\n")


if (!is.null(task)) {
  if (task == "model_evaluation") {
    # Task: Model Evaluation
    # You can add flags for specific model training and evaluation to use, e.g., -o LR-Chi Squared test selects LR algorithm based on the 
    # selected antibiotics through Chi Squared test method
    if (model == "LR-Chi Squared test"){
      # Logistic regression (Model 1)--------------------------------------------------
      model_Log1_IUT <- glm(IUT ~ Imipenem + Gentamicin + Ciprofloxacin + Meropenem + Aztreonam +
                              Cefotaxime + Cefepime + Ceftazidime + Ceftriaxone
                            , family = "binomial", data = train_strat_IUT)
      
      cat("Summary of LR model based on Chi Squared test including Null deviance, Residual deviance, and AIC, etc.:\n")
      print(summary(model_Log1_IUT))
      
      test$probs1_IUT<- predict(model_Log1_IUT, test, type = "response")
      
      test$pred_logreg1_IUT <- ifelse(test$probs1_IUT >= 0.5, 1, 0)
      
      # Model 1 evaluation
      confm_logreg1_IUT <- table(actual = test$IUT, prediction = test$pred_logreg1_IUT)
      
      metrics_logreg1_IUT <- compute_metrics(confm_logreg1_IUT)
      
      cat("Performance of LR model based on Chi-Squared test results:\n")
      metrics_logreg1_IUT
    }
    else if (model == "LR-model agnostic"){
      model_Log2_IUT <- glm(IUT ~  Meropenem 
                            , family = "binomial", data = train_strat_IUT)
      
      cat("Summary of LR model based on the model-agnostic approach including Null deviance, Residual deviance, and AIC, etc.:\n")
      print(summary(model_Log2_IUT))
      
      ##Prediction on test (Model 2)--------------------------------------------------
      test$probs2_IUT <- predict(model_Log2_IUT, test, type = "response")
      
      test$pred_logreg2_IUT <- ifelse(test$probs2_IUT >= 0.5, 1, 0)
      
      # Model 2 evaluation
      confm_logreg2_IUT <- table(actual = test$IUT, prediction = test$pred_logreg2_IUT)
      
      metrics_logreg2_IUT <- compute_metrics(confm_logreg2_IUT)
      
      cat("Performance of LR model based on the model-agnostic approach:\n")
      metrics_logreg2_IUT
    }
    
    else if (model == "LR-Wald test"){
      model_Log3_IUT <- glm(IUT ~ Meropenem + Aztreonam + Cefotaxime + Cefepime + Ceftriaxone
                            , family = "binomial", data = train_strat_IUT)
      
      cat("Summary of LR model based on the Wald test including Null deviance, Residual deviance, and AIC, etc.:\n")
      print(summary(model_Log3_IUT))
      
      ##Prediction on test
      test$probs3_IUT <- predict(model_Log3_IUT, test, type = "response")
      
      test$pred_logreg3_IUT <- ifelse(test$probs3_IUT >= 0.5, 1, 0)
      
      confm_logreg3_IUT <- table(actual = test$IUT, prediction = test$pred_logreg3_IUT)
      
      metrics_logreg3_IUT <- compute_metrics(confm_logreg3_IUT)
      
      cat("Performance of LR model based on the Wald test:\n")
      metrics_logreg3_IUT
    }
    
    else if (model == "LR-Simplicity Principle_CRO"){
      model_Log4_IUT <- glm(IUT ~ Ceftriaxone
                            , family = "binomial", data = train_strat_IUT)
      
      cat("Summary of LR model based on the simplicity principle (CRO) including Null deviance, Residual deviance, and AIC, etc.:\n")
      print(summary(model_Log4_IUT))
      
      ##Prediction on test
      test$probs4_IUT <- predict(model_Log4_IUT, test, type = "response")
      
      test$pred_logreg4_IUT <- ifelse(test$probs4_IUT >= 0.5, 1, 0)
      
      confm_logreg4_IUT <- table(actual = test$IUT, prediction = test$pred_logreg4_IUT)
      
      metrics_logreg4_IUT <- compute_metrics(confm_logreg4_IUT)
      
      cat("Performance of LR model based on the simplicity principle (CRO):\n")
      metrics_logreg4_IUT
    }
    
    else if (model == "NBC-Chi Squared test"){
      train_control <- trainControl(
        method = "cv",
        number = 10
      )
      
      search_grid_IUT <- expand.grid(
        usekernel = c(FALSE,TRUE),
        fL = seq(0,5,0.5),
        adjust = seq(0,5,0.5)
      )
      
      model_nbdalex_IUT <- train(IUT ~ Imipenem + Gentamicin + Ciprofloxacin + Meropenem + Aztreonam +
                                   Cefotaxime + Cefepime + Ceftazidime + Ceftriaxone,
                                 data = train_strat_IUT,
                                 method = "nb",
                                 metric = "Accuracy",
                                 trControl = train_control,
                                 tuneGrid = search_grid_IUT,
      )
      
      
      cat("Tuned hyper-parameters through 10-fold Cross-Validation:\n")
      print(model_nbdalex_IUT$bestTune)
      
      ##Prediction on test (Naive Bayes Model)----------------------------------
      test$pred_nb_IUT <- predict(model_nbdalex_IUT, test)
      
      # Model evaluation
      confm_nb_IUT <- table(actual = test$IUT, prediction = test$pred_nb_IUT)
      metrics_nb_IUT<- compute_metrics(confm_nb_IUT)
      cat("Performance of NBC model based on Chi-Squared test results:\n")
      metrics_nb_IUT
    }
    else if (model == "NBC-model agnostic"){
      train_control <- trainControl(
        method = "cv",
        number = 10
      )
      
      search_grid_IUT <- expand.grid(
        usekernel = c(FALSE,TRUE),
        fL = seq(0,5,0.5),
        adjust = seq(0,5,0.5)
      )
      
      model_nb_IUT <- train(IUT ~ Ceftriaxone,
                            data = train_strat_IUT,
                            method = "nb",
                            metric = "Accuracy",
                            trControl = train_control,
                            tuneGrid = search_grid_IUT,
      )
      
      cat("Tuned hyper-parameters through 10-fold Cross-Validation:\n")
      print(model_nb_IUT$bestTune)
      
      ##Prediction on test (Naive Bayes Model)----------------------------------
      test$pred_nb1_IUT <- predict(model_nb_IUT, test)
      
      # Model evaluation
      confm_nb1_IUT <- table(actual = test$IUT, prediction = test$pred_nb1_IUT)
      
      metrics_nb1_IUT<- compute_metrics(confm_nb1_IUT)
      cat("Performance of NBC model based on model-agnostic approach results:\n")
      metrics_nb1_IUT
    }
    
    else if (model == "NBC-Wald test"){
      
      train_control <- trainControl(
        method = "cv",
        number = 10
      )
      
      search_grid_IUT <- expand.grid(
        usekernel = c(FALSE,TRUE),
        fL = seq(0,5,0.5),
        adjust = seq(0,5,0.5)
      )
      
      model_nb2_IUT <- train( IUT ~ Meropenem + Aztreonam + Cefotaxime + Cefepime + Ceftriaxone,
                              data = train_strat_IUT,
                              method = "nb",
                              metric = "Accuracy",
                              trControl = train_control,
                              tuneGrid = search_grid_IUT,
      )
      
      cat("Tuned hyper-parameters through 10-fold Cross-Validation:\n")
      print(model_nb2_IUT$bestTune)
      
      ##Prediction on test (Naive Bayes Model)
      test$pred_nb2_IUT <- predict(model_nb2_IUT, test)
      
      ##confusion matrix(Naive Bayes Model)
      confm_nb2_IUT <- table(actual = test$IUT, prediction = test$pred_nb2_IUT)
      
      metrics_nb2_IUT<- compute_metrics(confm_nb2_IUT)
      
      cat("Performance of NBC model based on the Wald test results:\n")
      metrics_nb2_IUT
    }
    
    else if (model == "NBC-Simplicity Principle_MEM"){
      
      train_control <- trainControl(
        method = "cv",
        number = 10
      )
      
      search_grid_IUT <- expand.grid(
        usekernel = c(FALSE,TRUE),
        fL = seq(0,5,0.5),
        adjust = seq(0,5,0.5)
      )
      
      model_nb3_IUT <- train( IUT ~ Meropenem,
                              data = train_strat_IUT,
                              method = "nb",
                              metric = "Accuracy",
                              trControl = train_control,
                              tuneGrid = search_grid_IUT,
      )
      
      cat("Tuned hyper-parameters through 10-fold Cross-Validation:\n")
      print(model_nb3_IUT$bestTune)
      
      ##Prediction on test (Naive Bayes Model)
      test$pred_nb3_IUT <- predict(model_nb3_IUT, test)
      
      ##confusion matrix(Naive Bayes Model)
      confm_nb3_IUT <- table(actual = test$IUT, prediction = test$pred_nb3_IUT)
      
      metrics_nb3_IUT<- compute_metrics(confm_nb3_IUT)
      
      cat("Performance of NBC model based on the simplicity principle (MEM) results:\n")
      metrics_nb3_IUT
    }
    
    else if (model == "LDA-Chi Squared test"){
      
      model_ldadalex_IUT <- train(IUT ~ Cefotaxime + Gentamicin + Imipenem +
                                    Cefepime + Ceftazidime + Ciprofloxacin +
                                    Meropenem + Ceftriaxone + Aztreonam,
                                  data = train_strat_IUT,
                                  method = "lda",
                                  metric = "Accuracy")
      
      cat("LDA model based on Chi-Squared test results:\n")
      print(model_ldadalex_IUT)
      
      #Prediction on test (LDA Model)-----------------------------------------------
      test$pred_lda_IUT <- predict(model_ldadalex_IUT, test)
      
      # Model Evaluation-------------------------------------------------
      confm_lda_IUT <- table(actual = test$IUT, prediction = test$pred_lda_IUT)
      
      metrics_lda_IUT<- compute_metrics(confm_lda_IUT)
      
      cat("Performance of LDA model based on Chi-Squared test results:\n")
      metrics_lda_IUT
    }
    else if (model == "LDA-model agnostic"){
      model_lda_IUT <- train(IUT ~ Meropenem,
                             data = train_strat_IUT,
                             method = "lda",
                             metric = "Accuracy")
      
      cat("LDA model based on model-agnostic approach:\n")
      print(model_lda_IUT)
      
      #Prediction on test (LDA Model)-----------------------------------------------
      test$pred_lda1_IUT <- predict(model_lda_IUT, test)
      
      # Model evaluation
      confm_lda1_IUT <- table(actual = test$IUT, prediction = test$pred_lda1_IUT)
      metrics_lda1_IUT <- compute_metrics(confm_lda1_IUT)
      cat("Performance of LDA model based on model-agnostic approach:\n")
      metrics_lda1_IUT
    }
    
    else if (model == "LDA-Wald test"){
      model_lda2_IUT <- train(IUT ~ Meropenem + Aztreonam + Cefotaxime + Cefepime + Ceftriaxone,
                              data = train_strat_IUT,
                              method = "lda",
                              metric = "Accuracy")
      
      cat("LDA model based on the Wald test results:\n")
      print(model_lda2_IUT)
      
      test$pred_lda2_IUT <- predict(model_lda2_IUT, test)
      
      confm_lda2_IUT <- table(actual = test$IUT, prediction = test$pred_lda2_IUT)
      
      metrics_lda2_IUT<- compute_metrics(confm_lda2_IUT)
      
      cat("Performance of LDA model based on the Wald test results:\n")
      metrics_lda2_IUT
    }
    
    else if (model == "LDA-Simplicity Principle_CRO"){
      model_lda3_IUT <- train(IUT ~ Ceftriaxone ,
                              data = train_strat_IUT,
                              method = "lda",
                              metric = "Accuracy",
                              trControl = train_control)
      
      cat("LDA model based on the Simplicity Principle (CRO) results:\n")
      print(model_lda3_IUT)
      
      test$pred_lda3_IUT <- predict(model_lda3_IUT, test)
      
      confm_lda3_IUT <- table(actual = test$IUT, prediction = test$pred_lda3_IUT)
      
      metrics_lda3_IUT<- compute_metrics(confm_lda3_IUT)
      
      cat("Performance of LDA model based on the Simplicity Principle (CRO) results:\n")
      metrics_lda3_IUT
    }
    
    else if (model == "SVM-Sig-Chi Squared test"){
      # Fit Support Vector Machine model to train data set
      # This code performs a grid search for tuning the hyper-parameters of the SVM model using the tune function
      # from the e1071 package. It then fits the SVM model with the specified hyper-parameters using the svm function.
      # Model 1
      set.seed(1234)
      tune_out2_IUT <- e1071::tune("svm", IUT ~ Imipenem + Gentamicin + Ciprofloxacin + Meropenem + 
                                     Aztreonam + Cefotaxime + Cefepime + Ceftazidime + Ceftriaxone  
                                   , data = train_strat_IUT, kernel = "sigmoid", tunecontrol=tune.control(cross=10),
                                   ranges = list(cost = seq(0.1, 10, length = 25)
                                                 , gamma = seq(0.1,5, length = 25)))
      
      cat("Tuned hyper-parameters through 10-Fold Cross Validation:\n")
      print(summary(tune_out2_IUT))
      
      model_svm2_IUT = svm(IUT ~ Imipenem + Gentamicin + Ciprofloxacin + Meropenem + Aztreonam +
                             Cefotaxime + Cefepime + Ceftazidime + Ceftriaxone, data = train_strat_IUT, 
                           probability = T, kernel = "sigmoid", cost = tune_out2_IUT$best.parameters$cost,
                           gamma = tune_out2_IUT$best.parameters$gamma )
      
      cat("Sigmoid SVM based on Chi-Squared test results:\n")
      print(model_svm2_IUT)
      
      #Prediction on test (SVM)----------------------------------------------------
      test$pred_svm_IUT <- predict(model_svm2_IUT, test)
      
      # Model evaluation
      confm_svm_IUT <- table(actual = test$IUT, prediction = test$pred_svm_IUT)
      
      metrics_svm1_IUT_sig<- compute_metrics(confm_svm_IUT)
      
      cat("Performance of SVM-Sigmoid model based on Chi-Squared test results:\n")
      metrics_svm1_IUT_sig
    }
    
    else if (model=="SVM-Sig-model agnostic"){
      
      set.seed(1234)
      tune_out21_IUT <- e1071::tune("svm", IUT ~ Meropenem ,
                                    data = train_strat_IUT,
                                    kernel = "sigmoid", tunecontrol=tune.control(cross=10),
                                    ranges = list(cost = seq(0.1, 5, length = 25)
                                                  , gamma = seq(0, 5, length = 25)))
      
      cat("Tuned hyper-parameters through 10-Fold Cross Validation:\n")
      print(summary(tune_out21_IUT))
      
      model_svm21_IUT = svm(IUT ~ Meropenem, data = train_strat_IUT, 
                            probability = T, kernel = "sigmoid", cost = tune_out21_IUT$best.parameters$cost,
                            gamma = tune_out21_IUT$best.parameters$gamma)
      
      cat("Sigmoid SVM based on model-agnostic approach:\n")
      print(model_svm21_IUT)
      
      test$pred_svm21_IUT <- predict(model_svm21_IUT, test)
      
      confm_svm21_IUT <- table(actual = test$IUT, prediction = test$pred_svm21_IUT)
      
      metrics_svm21_IUT_sig<- compute_metrics(confm_svm21_IUT)
      
      cat("Performance of SVM-Sigmoid model based on model-agnostic approach:\n")
      metrics_svm21_IUT_sig
    }
    
    else if (model=="SVM-Wald test"){
      set.seed(1234)
      tune_out3_IUT <- e1071::tune("svm", IUT ~ Meropenem + Aztreonam + Cefotaxime + Ceftriaxone + Cefepime
                                   , data = train_strat_IUT, kernel = "sigmoid", tunecontrol=tune.control(cross=10),
                                   ranges = list(cost = seq(0.1, 5, length = 25)
                                                 , gamma = seq(0.1, 1, length = 20)))
      
      cat("Tuned hyper-parameters through 10-Fold Cross Validation:\n")
      print(summary(tune_out3_IUT))
      
      model_svm3_IUT = svm(IUT ~ Meropenem + Aztreonam + Cefotaxime + Ceftriaxone + Cefepime, data = train_strat_IUT, 
                           probability = T, kernel = "sigmoid", cost = tune_out3_IUT$best.parameters$cost,
                           gamma = tune_out3_IUT$best.parameters$gamma )
      
      cat("Sigmoid SVM based on the Wald test:\n")
      print(model_svm3_IUT)
      
      #Prediction on test (SVM)
      test$pred_svm3_IUT <- predict(model_svm3_IUT, test)
      
      # Model evaluation
      confm_svm3_IUT <- table(actual = test$IUT, prediction = test$pred_svm3_IUT)
      
      metrics_svm3_IUT_sig<- compute_metrics(confm_svm3_IUT)
      
      cat("Performance of SVM-Sigmoid model based on the Wald test results:\n")
      metrics_svm3_IUT_sig
    } 
    
    else if (model=="SVM-Simplicity Principle_CRO"){
      set.seed(1234)
      tune_out4_IUT <- e1071::tune("svm", IUT ~ Ceftriaxone ,
                                   data = train_strat_IUT,
                                   kernel = "sigmoid", tunecontrol=tune.control(cross=10),
                                   ranges = list(cost = seq(0.1, 5, length = 25)
                                                 , gamma = seq(0, 5, length = 25)))
      
      cat("Tuned hyper-parameters through 10-Fold Cross Validation:\n")
      print(summary(tune_out4_IUT))
      
      model_svm4_IUT = svm(IUT ~ Ceftriaxone, data = train_strat_IUT, 
                           probability = T, kernel = "sigmoid", cost = tune_out4_IUT$best.parameters$cost,
                           gamma = tune_out4_IUT$best.parameters$gamma )
      
      cat("Sigmoid SVM based on the Simplicity Principle (CRO):\n")
      print(model_svm4_IUT)
      
      #Prediction on test (SVM)
      test$pred_svm4_IUT <- predict(model_svm4_IUT, test)
      
      # Model evaluation
      confm_svm4_IUT <- table(actual = test$IUT, prediction = test$pred_svm4_IUT)
      
      metrics_svm4_IUT_sig<- compute_metrics(confm_svm4_IUT)
      
      cat("Performance of SVM-Sigmoid model based on the Simplicity Principle (CRO) results:\n")
      metrics_svm4_IUT_sig
    } 
    
    else if (model=="SVM-Poly-Chi Squared test"){
      set.seed(1234)
      tune_out5_IUT <- e1071::tune("svm", IUT ~ Imipenem + Gentamicin + Ciprofloxacin + Meropenem + 
                                     Aztreonam + Cefotaxime + Cefepime + Ceftazidime + Ceftriaxone  
                                   , data = train_strat_IUT, kernel = "polynomial", tunecontrol=tune.control(cross=10),
                                   ranges = list(cost = seq(0.1, 5, length = 20), degree = seq(1,4,1)
                                                 , gamma = seq(0.1,5, length = 20)))
      
      cat("Tuned hyper-parameters through 10-Fold Cross Validation:\n")
      print(summary(tune_out5_IUT))
      
      # Tuned hyper-parameters are placed in the model.
      model_svm5_IUT = svm(IUT ~ Imipenem + Gentamicin + Ciprofloxacin + Meropenem + Aztreonam +
                             Cefotaxime + Cefepime + Ceftazidime + Ceftriaxone
                           , data = train_strat_IUT, probability = T, 
                           kernel = "polynomial", cost = tune_out5_IUT$best.parameters$cost, 
                           degree = tune_out5_IUT$best.parameters$degree , gamma = tune_out5_IUT$best.parameters$gamma )
      
      cat("Polynomial SVM based on Chi-Squared test results:\n")
      print(model_svm5_IUT)
      
      #Prediction on test (SVM 2)----------------------------------------------------
      test$pred_svm5_IUT <- predict(model_svm5_IUT, test)
      
      # Model evaluation
      confm_svm5_IUT <- table(actual = test$IUT, prediction = test$pred_svm5_IUT)
      
      metrics_svm_IUT<- compute_metrics(confm_svm5_IUT)
      
      cat("Performance of SVM-Polynomial model based on Chi-Squared test results:\n")
      metrics_svm_IUT
    }
    else if (model=="SVM-Poly-model agnostic"){
      
      set.seed(1234)
      tune_out51_IUT <- e1071::tune("svm", IUT ~ Meropenem ,
                                    data = train_strat_IUT,
                                    kernel = "polynomial", tunecontrol=tune.control(cross=10),
                                    ranges = list(cost = seq(0.1, 5, length = 20), degree = seq(1,3,1)
                                                  , gamma = seq(0, 5, length = 20)))
      
      cat("Tuned hyper-parameters through 10-Fold Cross Validation:\n")
      print(summary(tune_out51_IUT))
      
      model_svm51_IUT = svm(IUT ~ Meropenem, data = train_strat_IUT, 
                            probability = T, kernel = "polynomial", cost = tune_out51_IUT$best.parameters$cost,
                            degree = tune_out51_IUT$best.parameters$degree , gamma = tune_out51_IUT$best.parameters$gamma )
      
      cat("Polynomial SVM based on model-agnostic approach:\n")
      print(model_svm51_IUT)
      
      #Prediction on test (SVM)----------------------------------------------------
      test$pred_svm51_IUT <- predict(model_svm51_IUT, test)
      
      # Model evaluation
      confm_svm51_IUT <- table(actual = test$IUT, prediction = test$pred_svm51_IUT)
      
      metrics_svm1_IUT<- compute_metrics(confm_svm51_IUT)
      
      cat("Performance of SVM-Polynomial model based on model-agnostic approach:\n")
      metrics_svm1_IUT
      
    }
    
    else if (model=="SVM-Poly-Wald test"){
      
      set.seed(1234)
      tune_out6_IUT <- e1071::tune("svm", IUT ~ Meropenem + Cefepime +
                                     Aztreonam + Cefotaxime + Ceftriaxone,
                                   data = train_strat_IUT,
                                   kernel = "polynomial", tunecontrol=tune.control(cross=10),
                                   ranges = list(cost = seq(0.1, 5, length = 20), degree = seq(1,4,1)
                                                 , gamma = seq(0, 5, length = 20)))
      
      cat("Tuned hyper-parameters through 10-Fold Cross Validation:\n")
      print(summary(tune_out6_IUT))
      
      model_svm6_IUT = svm(IUT ~ Meropenem + Cefepime +
                             Aztreonam + Cefotaxime + Ceftriaxone , data = train_strat_IUT, 
                           probability = T, kernel = "polynomial", cost = tune_out6_IUT$best.parameters$cost,
                           degree = tune_out6_IUT$best.parameters$degree , gamma = tune_out6_IUT$best.parameters$gamma )
      
      cat("Polynomial SVM based on the Wald test:\n")
      print(model_svm6_IUT)
      
      #Prediction on test (SVM)----------------------------------------------------
      test$pred_svm6_IUT <- predict(model_svm6_IUT, test)
      
      # Model evaluation
      confm_svm6_IUT <- table(actual = test$IUT, prediction = test$pred_svm6_IUT)
      
      metrics_svm2_IUT<- compute_metrics(confm_svm6_IUT)
      
      cat("Performance of SVM-Polynomial model based on the Wald test:\n")
      metrics_svm2_IUT
    }
    
    else if (model=="SVM-Poly-Simplicity Principle_CRO"){
      
      set.seed(1234)
      tune_out7_IUT <- e1071::tune("svm", IUT ~ Ceftriaxone,
                                   data = train_strat_IUT,
                                   kernel = "polynomial", tunecontrol=tune.control(cross=10),
                                   ranges = list(cost = seq(0.1, 5, length = 20), degree = seq(1,3,1)
                                                 , gamma = seq(0, 5, length = 20)))
      
      cat("Tuned hyper-parameters through 10-Fold Cross Validation:\n")
      print(summary(tune_out7_IUT))
      
      model_svm7_IUT = svm(IUT ~ Ceftriaxone, data = train_strat_IUT, 
                           probability = T, kernel = "polynomial", cost = tune_out7_IUT$best.parameters$cost,
                           degree = tune_out7_IUT$best.parameters$degree , gamma = tune_out7_IUT$best.parameters$gamma )
      
      cat("Polynomial SVM based on the Simplicity Principle (CRO):\n")
      print(model_svm7_IUT)
      
      #Prediction on test (SVM)----------------------------------------------------
      test$pred_svm7_IUT <- predict(model_svm7_IUT, test)
      
      # Model evaluation
      confm_svm7_IUT <- table(actual = test$IUT, prediction = test$pred_svm7_IUT)
      
      metrics_svm3_IUT<- compute_metrics(confm_svm7_IUT)
      
      cat("Performance of SVM-Polynomial model based on the Simplicity Principle (CRO):\n")
      metrics_svm3_IUT
    }
    
    else if (model=="CatBoost-Chi Squared test"){
      # CatBoost Models
      Bio_Data <- as.data.frame(Bio_Data)
      
      #Convert categorical variables to factor
      categorical_var <- c(colnames(Bio_Data))
      
      Bio_Data[, categorical_var[-1]] <- lapply(Bio_Data[, categorical_var[-1]], factor)
      
      levels(Bio_Data$Ampicillin) <- c(levels(Bio_Data$Ampicillin),0)
      
      train_strat_IUT <- ovun.sample(IUT~Imipenem + Gentamicin + Ciprofloxacin + Meropenem + Aztreonam +
                                       Cefotaxime + Cefepime + Ceftazidime + Ceftriaxone, data=train, 
                                     method = "both", p = 0.5, seed = 1234, N = 272)$data
      
      set.seed(123)
      
      train_cases <- sample(1:nrow(Bio_Data), nrow(Bio_Data) * 0.8)
      
      test <- Bio_Data[-train_cases,]
      
      # Training CatBoost model based on the selected antibiotics through chi squared test
      grid_IUT <- expand.grid(
        depth = c(3,4,5),
        learning_rate = c(0.2,0.3,0.4),
        l2_leaf_reg = c(1,2,3),
        rsm = 1,
        border_count = 1,
        iterations = seq(10,30,1)
      )
      fitControl_IUT <- trainControl(method = "cv",
                                     number = 10,
                                     classProbs = TRUE)
      set.seed(1234)
      model_catboost_IUT <- train(x = train_strat_IUT[,c("Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                         "Ceftazidime", "Ciprofloxacin", "Meropenem", "Ceftriaxone", "Aztreonam")],
                                  y = make.names(train_strat_IUT[,"IUT"]),
                                  maximize = TRUE,
                                  method = catboost.caret, metric = "Accuracy", 
                                  tuneGrid =  grid_IUT,
                                  trControl = fitControl_IUT)
      
      cat("Tuned hyper-parameters through 10-Fold Cross Validation:\n")
      print(model_catboost_IUT$bestTune)  
      
      #Prediction on test (CatBoost)----------------------------------------------------
      test$cat_IUT = predict(model_catboost_IUT, test)
      
      # Model evaluation
      cat_test_IUT <- table(actual = test$IUT, prediction = test$cat_IUT)
      
      metrics_cat_IUT<- compute_metrics(cat_test_IUT)
      
      cat("Performance of CatBoost model based on Chi-Squared test results:\n")
      metrics_cat_IUT
    }
    
    else if (model=="CatBoost-model agnostic"){
      # CatBoost Models
      Bio_Data <- as.data.frame(Bio_Data)
      
      #Convert categorical variables to factor
      categorical_var <- c(colnames(Bio_Data))
      
      Bio_Data[, categorical_var[-1]] <- lapply(Bio_Data[, categorical_var[-1]], factor)
      
      levels(Bio_Data$Ampicillin) <- c(levels(Bio_Data$Ampicillin),0)
      
      train_strat_IUT <- ovun.sample(IUT~Imipenem + Gentamicin + Ciprofloxacin + Meropenem + Aztreonam +
                                       Cefotaxime + Cefepime + Ceftazidime + Ceftriaxone, data=train, 
                                     method = "both", p = 0.5, seed = 1234, N = 272)$data
      
      set.seed(123)
      
      train_cases <- sample(1:nrow(Bio_Data), nrow(Bio_Data) * 0.8)
      
      test <- Bio_Data[-train_cases,]
      
      grid_IUT <- expand.grid(
        depth = c(2,3),
        learning_rate = c(0.1,0.2),
        l2_leaf_reg = c(1,2),
        rsm = 1,
        border_count = 1,
        iterations = seq(30,60,1)
      )
      fitControl_IUT <- trainControl(method = "cv",
                                     number = 10,
                                     classProbs = TRUE)
      model_catboost1_IUT <- train(x = train_strat_IUT[,c("Ceftriaxone", "Aztreonam", "Cefepime", "Meropenem")],
                                   y = make.names(train_strat_IUT[,"IUT"]),
                                   maximize = TRUE,
                                   method = catboost.caret, metric = "Accuracy", 
                                   tuneGrid =  grid_IUT,
                                   trControl = fitControl_IUT)

      cat("Tuned hyper-parameters through 10-Fold Cross Validation:\n")
      print(model_catboost1_IUT$bestTune)  
      
      test$cat_IUT1 = predict(model_catboost1_IUT, test)
      
      # Model Evaluation
      confm_cat1_IUT <- table(actual = test$IUT, prediction = test$cat_IUT1)
      
      metrics_cat1_IUT<- compute_metrics(confm_cat1_IUT)
      
      cat("Performance of CatBoost model based on model-agnostic approach results:\n")
      metrics_cat1_IUT
    }
    
    else if (model=="CatBoost-Wald test"){
      # CatBoost Models
      Bio_Data <- as.data.frame(Bio_Data)
      
      #Convert categorical variables to factor
      categorical_var <- c(colnames(Bio_Data))
      
      Bio_Data[, categorical_var[-1]] <- lapply(Bio_Data[, categorical_var[-1]], factor)
      
      levels(Bio_Data$Ampicillin) <- c(levels(Bio_Data$Ampicillin),0)
      
      train_strat_IUT <- ovun.sample(IUT~Imipenem + Gentamicin + Ciprofloxacin + Meropenem + Aztreonam +
                                       Cefotaxime + Cefepime + Ceftazidime + Ceftriaxone, data=train, 
                                     method = "both", p = 0.5, seed = 1234, N = 272)$data
      
      set.seed(123)
      
      train_cases <- sample(1:nrow(Bio_Data), nrow(Bio_Data) * 0.8)
      
      test <- Bio_Data[-train_cases,]
      
      grid2_IUT <- expand.grid(
        depth = c(2,3,4),
        learning_rate = c(0.1,0.2,0.3),
        l2_leaf_reg = c(1,2,3),
        rsm = 1,
        border_count = 1,
        iterations = seq(20,40,1)
      )
      
      fitControl2_IUT <- trainControl(method = "cv",
                                      number = 10,
                                      classProbs = TRUE)
      model_catboost2_IUT <- train(x = train_strat_IUT[,c("Aztreonam","Ceftriaxone","Meropenem","Cefotaxime","Cefepime")],
                                   y = make.names(train_strat_IUT[["IUT"]]),
                                   maximize = TRUE,
                                   method = catboost.caret, metric = "Accuracy", 
                                   tuneGrid =  grid2_IUT,
                                   trControl = fitControl2_IUT)
      
      
      cat("Tuned hyper-parameters through 10-Fold Cross Validation:\n")
      print(model_catboost2_IUT$bestTune)  
      
      test$cat_IUT2 = predict(model_catboost2_IUT, test)
      # Model evaluation
      confm_cat2_IUT <- table(actual = test$IUT, prediction = test$cat_IUT2)
      
      metrics_cat3_IUT<- compute_metrics(confm_cat2_IUT)
      
      cat("Performance of CatBoost model based on the Wald test results:\n")
      metrics_cat3_IUT
    }
    
    else if (model=="CatBoost-Simplicity Principle_MEM"){
      # CatBoost Models
      Bio_Data <- as.data.frame(Bio_Data)
      
      #Convert categorical variables to factor
      categorical_var <- c(colnames(Bio_Data))
      
      Bio_Data[, categorical_var[-1]] <- lapply(Bio_Data[, categorical_var[-1]], factor)
      
      levels(Bio_Data$Ampicillin) <- c(levels(Bio_Data$Ampicillin),0)
      
      train_strat_IUT <- ovun.sample(IUT~Imipenem + Gentamicin + Ciprofloxacin + Meropenem + Aztreonam +
                                       Cefotaxime + Cefepime + Ceftazidime + Ceftriaxone, data=train, 
                                     method = "both", p = 0.5, seed = 1234, N = 272)$data
      
      set.seed(123)
      
      train_cases <- sample(1:nrow(Bio_Data), nrow(Bio_Data) * 0.8)
      
      test <- Bio_Data[-train_cases,]
      
      grid3_IUT <- expand.grid(
        depth = c(2,3,4),
        learning_rate = c(0.1,0.2,0.3),
        l2_leaf_reg = c(1,2,3),
        rsm = 1,
        border_count = 1,
        iterations = seq(10,100,1)
      )
      
      fitControl3_IUT <- trainControl(method = "cv",
                                      number = 10,
                                      classProbs = TRUE)
      model_catboost3_IUT <- train(x = train_strat_IUT["Meropenem"],
                                   y = make.names(train_strat_IUT[["IUT"]]),
                                   maximize = TRUE,
                                   method = catboost.caret, metric = "Accuracy", 
                                   tuneGrid =  grid3_IUT,
                                   trControl = fitControl3_IUT)
      
      
      cat("Tuned hyper-parameters through 10-Fold Cross Validation:\n")
      print(model_catboost3_IUT$bestTune)  
      
      test$cat_IUT3 = predict(model_catboost3_IUT, test)
      # Model evaluation
      confm_cat3_IUT <- table(actual = test$IUT, prediction = test$cat_IUT3)
      
      metrics_cat4_IUT<- compute_metrics(confm_cat3_IUT)
      
      cat("Performance of CatBoost model based on the Simplicity Principle (MEM):\n")
      metrics_cat4_IUT
    }
    
    else if (model=="CatBoost-Simplicity Principle_CRO"){
      # CatBoost Models
      Bio_Data <- as.data.frame(Bio_Data)
      
      #Convert categorical variables to factor
      categorical_var <- c(colnames(Bio_Data))
      
      Bio_Data[, categorical_var[-1]] <- lapply(Bio_Data[, categorical_var[-1]], factor)
      
      levels(Bio_Data$Ampicillin) <- c(levels(Bio_Data$Ampicillin),0)
      
      train_strat_IUT <- ovun.sample(IUT~Imipenem + Gentamicin + Ciprofloxacin + Meropenem + Aztreonam +
                                       Cefotaxime + Cefepime + Ceftazidime + Ceftriaxone, data=train, 
                                     method = "both", p = 0.5, seed = 1234, N = 272)$data
      
      set.seed(123)
      
      train_cases <- sample(1:nrow(Bio_Data), nrow(Bio_Data) * 0.8)
      
      test <- Bio_Data[-train_cases,]
      
      grid4_IUT <- expand.grid(
        depth = c(2,3,4),
        learning_rate = c(0.1,0.2,0.3),
        l2_leaf_reg = c(1,2,3),
        rsm = 1,
        border_count = 1,
        iterations = seq(10,100,10)
      )
      
      fitControl4_IUT <- trainControl(method = "cv",
                                      number = 5,
                                      classProbs = TRUE)
      model_catboost4_IUT <- train(x = train_strat_IUT["Ceftriaxone"],
                                   y = make.names(train_strat_IUT[["IUT"]]),
                                   maximize = TRUE,
                                   method = catboost.caret, metric = "Accuracy", 
                                   tuneGrid =  grid4_IUT,
                                   trControl = fitControl4_IUT)
      
      
      cat("Tuned hyper-parameters through 10-Fold Cross Validation:\n")
      print(model_catboost4_IUT$bestTune)  
      
      test$cat_IUT4 = predict(model_catboost4_IUT, test)
      # Model evaluation
      confm_cat4_IUT <- table(actual = test$IUT, prediction = test$cat_IUT4)
      
      metrics_cat5_IUT<- compute_metrics(confm_cat4_IUT)
      
      cat("Performance of CatBoost model based on the Simplicity Principle (CRO):\n")
      metrics_cat5_IUT
    }
  }
  
  else if (task == 'feature_selectio') {

    if (method == "Chi-Squared test"){
      # Function to perform chi-square test and print p-value
      perform_chi_square <- function(data, variable) {
        result <- chisq.test(table(data[[variable]], data$IUT))
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
    model_Log1_IUT <- glm(IUT ~ Imipenem + Gentamicin + Ciprofloxacin + Meropenem + Aztreonam +
                            Cefotaxime + Cefepime + Ceftazidime + Ceftriaxone
                          , family = "binomial", data = train_strat_IUT)
    
    cat("P-values show Wald test results:\n")
    print(summary(model_Log1_IUT))
  }
  else if (method == "model agnostic-LR"){
    # Feature Selection by DALEX Package--------------------------------------
    model_Log1_IUT <- glm(IUT ~ Imipenem + Gentamicin + Ciprofloxacin + Meropenem + Aztreonam +
                            Cefotaxime + Cefepime + Ceftazidime + Ceftriaxone
                          , family = "binomial", data = train_strat_IUT)
    
    explained_glm_IUT <- explain(model = model_Log1_IUT, data=train_strat_IUT[,c(1:9)], variables = colnames(train_strat_IUT[,c(1:9)]),
                                 y=as.vector(as.numeric(train_strat_IUT$IUT))-1, label = "LR", 
                                 type = "classification")
    #50 permuatation
    fi_glm_IUT = variable_importance(explained_glm_IUT, B=50 ,variables = colnames(train_strat_IUT[,c(1:9)]), loss_function = loss_root_mean_square, type = "raw")
    
    cat("Root mean squared dropout loss due to a feature elimination through model-agnostic approach in LR:\n")
    print(fi_glm_IUT)
  }
  else if (method == "model agnostic-NBC"){
    train_control <- trainControl(
      method = "repeatedcv", 
      repeats = 5,
      number = 10
    )
    
    search_grid_IUT <- expand.grid(
      usekernel = c(FALSE,TRUE),
      fL = seq(0,5,0.5),
      adjust = seq(0,5,0.5)
    )
    
    model_nbdalex_IUT <- train(IUT ~ Imipenem + Gentamicin + Ciprofloxacin + Meropenem + Aztreonam +
                                 Cefotaxime + Cefepime + Ceftazidime + Ceftriaxone,
                               data = train_strat_IUT,
                               method = "nb",
                               metric = "Accuracy",
                               trControl = train_control,
                               tuneGrid = search_grid_IUT,
    )
    
    expnb_IUT = explain(model = model_nbdalex_IUT, data=train_strat_IUT[,c("Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                           "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")], variables = colnames(train_strat_IUT[,c("Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                                                                                                                                             "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")]),
                        y=as.vector(as.numeric(train_strat_IUT$IUT))-1, label = "NBC", 
                        type = "classification")
    
    # The 'variable_importance' function is then used to calculate the variable importance based on the
    # explained model. The 'B' parameter is set to 50 for the number of permutations, and the loss function
    # is set to "loss_root_mean_square". The resulting variable importance is stored in 'fitnb_TEM'
    fitnb_IUT = variable_importance(expnb_IUT, B=50 ,variables = colnames(train_strat_IUT[,c("Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                                             "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")]), loss_function = loss_root_mean_square, type = "raw" )
    
    cat("Root mean squared dropout loss due to a feature elimination through model-agnostic approach in NBC:\n")
    print(fitnb_IUT)
  }
  
  else if (method == "model agnostic-LDA"){
    model_ldadalex_IUT <- train(IUT ~ Cefotaxime + Gentamicin + Imipenem +
                                  Cefepime + Ceftazidime + Ciprofloxacin +
                                  Meropenem + Ceftriaxone + Aztreonam,
                                data = train_strat_IUT,
                                method = "lda",
                                metric = "Accuracy")
    
    explda_IUT = explain(model = model_ldadalex_IUT, data=train_strat_IUT[,c("Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                             "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")], variables = colnames(train_strat_IUT[,c("Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                                                                                                                                               "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")]),
                         y=as.vector(as.numeric(train_strat_IUT$IUT))-1, label = "LDA", 
                         type = "classification")
    
    fitlda_IUT = variable_importance(explda_IUT, B=50 ,variables = colnames(train_strat_IUT[,c("Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                                               "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")]), loss_function = loss_root_mean_square, type = "raw" )
    
    cat("Root mean squared dropout loss due to a feature elimination through model-agnostic approach in LDA:\n")
    print(fitlda_IUT)
  }
  else if (method == "model agnostic-Sig SVM"){
    set.seed(1234)
    tune_out2_IUT <- e1071::tune("svm", IUT ~ Imipenem + Gentamicin + Ciprofloxacin + Meropenem + 
                                   Aztreonam + Cefotaxime + Cefepime + Ceftazidime + Ceftriaxone  
                                 , data = train_strat_IUT, kernel = "sigmoid", tunecontrol=tune.control(cross=10),
                                 ranges = list(cost = seq(0.1, 10, length = 25)
                                               , gamma = seq(0.1,5, length = 25)))
    
    model_svm2_IUT = svm(IUT ~ Imipenem + Gentamicin + Ciprofloxacin + Meropenem + Aztreonam +
                           Cefotaxime + Cefepime + Ceftazidime + Ceftriaxone, data = train_strat_IUT, 
                         probability = T, kernel = "sigmoid", cost = tune_out2_IUT$best.parameters$cost,
                         gamma = tune_out2_IUT$best.parameters$gamma )
    
    explained_SVM2_IUT <- explain(model = model_svm2_IUT, data=train_strat_IUT[,c(1:9)], variables = colnames(train_strat_IUT[,c(1:9)]),
                                  y=as.vector(as.numeric(train_strat_IUT$IUT))-1, label = "SVM-Sigmoid", 
                                  type = "classification")
    
    #50 permuatation
    fi_SVM2_IUT = variable_importance(explained_SVM2_IUT, B=50 ,variables = colnames(train_strat_IUT[,c(1:9)]),loss_function = loss_root_mean_square, type = "raw")
    
    cat("Root mean squared dropout loss due to a feature elimination through model-agnostic approach in sigmoid kernel of SVM:\n")
    print(fi_SVM2_IUT)
  }
  
  else if (method == "model agnostic-Poly SVM"){
    
    set.seed(1234)
    tune_out5_IUT <- e1071::tune("svm", IUT ~ Imipenem + Gentamicin + Ciprofloxacin + Meropenem + 
                                   Aztreonam + Cefotaxime + Cefepime + Ceftazidime + Ceftriaxone  
                                 , data = train_strat_IUT, kernel = "polynomial", tunecontrol=tune.control(cross=10),
                                 ranges = list(cost = seq(0.1, 5, length = 20), degree = seq(1,4,1)
                                               , gamma = seq(0.1,5, length = 20)))
    
    model_svm5_IUT = svm(IUT ~ Imipenem + Gentamicin + Ciprofloxacin + Meropenem + Aztreonam +
                           Cefotaxime + Cefepime + Ceftazidime + Ceftriaxone
                         , data = train_strat_IUT, probability = T, 
                         kernel = "polynomial", cost = tune_out5_IUT$best.parameters$cost, 
                         degree = tune_out5_IUT$best.parameters$degree , gamma = tune_out5_IUT$best.parameters$gamma )
    
    explained_SVM5_IUT <- explain(model = model_svm5_IUT, data=train_strat_IUT[,c(1:9)], variables = colnames(train_strat_IUT[,c(1:9)]),
                                  y=as.vector(as.numeric(train_strat_IUT$IUT))-1, label = "SVM-Polynomial", 
                                  type = "classification")
    
    #50 permuatation
    fi_SVM5 = variable_importance(explained_SVM5_IUT, B=50 ,variables = colnames(train_strat_IUT[,c(1:9)]),loss_function = loss_root_mean_square, type = "raw")
    
    cat("Root mean squared dropout loss due to a feature elimination through model-agnostic approach in polynomial kernel of SVM:\n")
    print(fi_SVM5)
  }
  else if (method == "model agnostic-CatBoost"){
    # CatBoost Models
    Bio_Data <- as.data.frame(Bio_Data)
    
    #Convert categorical variables to factor
    categorical_var <- c(colnames(Bio_Data))
    
    Bio_Data[, categorical_var[-1]] <- lapply(Bio_Data[, categorical_var[-1]], factor)
    
    levels(Bio_Data$Ampicillin) <- c(levels(Bio_Data$Ampicillin),0)
    
    train_strat_IUT <- ovun.sample(IUT~Imipenem + Gentamicin + Ciprofloxacin + Meropenem + Aztreonam +
                                     Cefotaxime + Cefepime + Ceftazidime + Ceftriaxone, data=train, 
                                   method = "both", p = 0.5, seed = 1234, N = 272)$data
    
    set.seed(123)
    
    train_cases <- sample(1:nrow(Bio_Data), nrow(Bio_Data) * 0.8)
    
    test <- Bio_Data[-train_cases,]
    
    # Training CatBoost model based on the selected antibiotics through chi squared test
    grid_IUT <- expand.grid(
      depth = c(3,4,5),
      learning_rate = c(0.2,0.3,0.4),
      l2_leaf_reg = c(1,2,3),
      rsm = 1,
      border_count = 1,
      iterations = seq(10,30,1)
    )
    fitControl_IUT <- trainControl(method = "cv",
                                   number = 10,
                                   classProbs = TRUE)
    set.seed(1234)
    model_catboost_IUT <- train(x = train_strat_IUT[,c("Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                       "Ceftazidime", "Ciprofloxacin", "Meropenem", "Ceftriaxone", "Aztreonam")],
                                y = make.names(train_strat_IUT[,"IUT"]),
                                maximize = TRUE,
                                method = catboost.caret, metric = "Accuracy", 
                                tuneGrid =  grid_IUT,
                                trControl = fitControl_IUT)
    
    expsda_IUT = explain(model = model_catboost_IUT, data=train_strat_IUT[,c("Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                             "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")], variables = colnames(train_strat_IUT[,c("Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                                                                                                                                               "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")]),
                         y=as.vector(as.numeric(train_strat_IUT$IUT))-1, label = "CatBoost", 
                         type = "classification")
    
    fitcat_IUT = variable_importance(expsda_IUT, B=50 ,variables = colnames(train_strat_IUT[,c("Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                                               "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")]), loss_function = loss_root_mean_square, type = "raw" )
    
    cat("Root mean squared dropout loss due to a feature elimination through model-agnostic approach in CatBoost:\n")
    print(fitcat_IUT)
  }
  else {
    cat("Invalid task. Please specify one of the following tasks: feature_selectio, model_evaluation\n")
    
  } 
}
if (is.null(task)){
  cat("You should specify a task\n")
}

