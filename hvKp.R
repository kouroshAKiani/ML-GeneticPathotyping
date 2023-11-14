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

# Calculate the proportion of HVKP_gene

proportion_HVKP <- table(Bio_Data$HVKP) %>% prop.table()

cat("Proportion of HVKP gene:\n")

print(proportion_HVKP)
# Response variable is imbalanced (0 : 70% , 1 : 30%) so we're using stratified sampling

# Divide the data set into train and test sets through stratified sampling method

set.seed(123)

split_strat  <- initial_split(Bio_Data, prop = 0.8, strata = "HVKP")

train_strat_HVKP  <- training(split_strat)

test_strat_HVKP   <- testing(split_strat)

# Display the dimensions of the train and test sets
cat("train set dimensions:", dim(train_strat_HVKP), "\n")
cat("test set dimensions:", dim(test_strat_HVKP), "\n")


if (!is.null(task)) {
  if (task == "model_evaluation") {
    
    if (model == "LR-Chi Squared test"){
      # Logistic regression (Model 1)--------------------------------------------------
      model_Log1_HVKP <- glm(HVKP ~ Imipenem + Gentamicin + Ciprofloxacin + Meropenem + Aztreonam +
                               Cefotaxime + Cefepime + Ceftazidime + Ceftriaxone
                             , family = "binomial", data = train_strat_HVKP)
      
      cat("Summary of LR model based on Chi Squared test including Null deviance, Residual deviance, and AIC, etc.:\n")
      print(summary(model_Log1_HVKP))
      
      test_strat_HVKP$probs1_HVKP <- predict(model_Log1_HVKP, test_strat_HVKP, type = "response")
      
      test_strat_HVKP$pred_logreg1_HVKP <- ifelse(test_strat_HVKP$probs1_HVKP >= 0.5, 1, 0)
      
      # Model 1 evaluation
      confm_logreg1_HVKP <- table(actual = test_strat_HVKP$HVKP, prediction = test_strat_HVKP$pred_logreg1_HVKP)
      
      metrics_logreg1_HVKP <- compute_metrics(confm_logreg1_HVKP)
      
      cat("Performance of LR model based on Chi-Squared test results:\n")
      metrics_logreg1_HVKP
    }
    
    else if (model == "LR-Wald test"){
      model_Log2_HVKP <- glm( HVKP ~  Imipenem + Meropenem + Aztreonam + Ceftriaxone + Ceftazidime
                              , family = "binomial", data = train_strat_HVKP)
      
      cat("Summary of LR model based on the Wald test including Null deviance, Residual deviance, and AIC, etc.:\n")
      print(summary(model_Log2_HVKP))
      
      ##Prediction on test_strat_HVKP
      test_strat_HVKP$probs2_HVKP <- predict(model_Log2_HVKP, test_strat_HVKP, type = "response")
      
      test_strat_HVKP$pred_logreg2_HVKP <- ifelse(test_strat_HVKP$probs2_HVKP >= 0.5, 1, 0)
      
      confm_logreg2_HVKP <- table(actual = test_strat_HVKP$HVKP, prediction = test_strat_HVKP$pred_logreg2_HVKP)
      
      metrics_logreg2_HVKP <- compute_metrics(confm_logreg2_HVKP)
      
      cat("Performance of LR model based on the Wald test:\n")
      metrics_logreg2_HVKP
    }
    
    else if (model == "LR-model agnostic"){
      model_Log3_HVKP <- glm(HVKP ~  Meropenem 
                             , family = "binomial", data = train_strat_HVKP)
      
      cat("Summary of LR model based on the model-agnostic approach including Null deviance, Residual deviance, and AIC, etc.:\n")
      print(summary(model_Log3_HVKP))
      
      ##Prediction on test_strat_HVKP (Model 2)--------------------------------------------------
      test_strat_HVKP$probs3_HVKP <- predict(model_Log3_HVKP, test_strat_HVKP, type = "response")
      
      test_strat_HVKP$pred_logreg3_HVKP <- ifelse(test_strat_HVKP$probs3_HVKP >= 0.5, 1, 0)
      
      # Model 2 evaluation
      confm_logreg3_HVKP <- table(actual = test_strat_HVKP$HVKP, prediction = test_strat_HVKP$pred_logreg3_HVKP)
      
      metrics_logreg3_HVKP <- compute_metrics(confm_logreg3_HVKP)
      
      cat("Performance of LR model based on the model-agnostic approach:\n")
      metrics_logreg3_HVKP
    }
    
    else if (model == "NBC-Chi Squared test"){
      train_control <- trainControl(
        method = "cv",
        number = 10
      )
      
      search_grid_HVKP <- expand.grid(
        usekernel = c(FALSE,TRUE),
        fL = seq(0,5,0.5),
        adjust = seq(0,5,0.5)
      )
      
      model_nbdalex_HVKP <- train(HVKP ~ Imipenem + Gentamicin + Ciprofloxacin + Meropenem + Aztreonam +
                                    Cefotaxime + Cefepime + Ceftazidime + Ceftriaxone,
                                  data = train_strat_HVKP,
                                  method = "nb",
                                  metric = "Accuracy",
                                  trControl = train_control,
                                  tuneGrid = search_grid_HVKP
      )
      
      
      cat("Tuned hyper-parameters through 10-fold Cross-Validation:\n")
      print(model_nbdalex_HVKP$bestTune)
      
      ##Prediction on test_strat_HVKP (Naive Bayes Model)----------------------------------
      test_strat_HVKP$pred_nb_HVKP <- predict(model_nbdalex_HVKP, test_strat_HVKP)
      
      # Model evaluation
      confm_nb_HVKP <- table(actual = test_strat_HVKP$HVKP, prediction = test_strat_HVKP$pred_nb_HVKP)
      metrics_nb_HVKP<- compute_metrics(confm_nb_HVKP)
      cat("Performance of NBC model based on Chi-Squared test results:\n")
      metrics_nb_HVKP
    }
    else if (model == "NBC-model agnostic"){
      train_control <- trainControl(
        method = "cv",
        number = 10
      )
      
      search_grid_HVKP <- expand.grid(
        usekernel = c(FALSE,TRUE),
        fL = seq(0,5,0.5),
        adjust = seq(0,5,0.5)
      )
      
      model_nb_HVKP <- train(HVKP ~ Meropenem + Ceftazidime,
                             data = train_strat_HVKP,
                             method = "nb",
                             metric = "Accuracy",
                             trControl = train_control,
                             tuneGrid = search_grid_HVKP
      )
      
      
      cat("Tuned hyper-parameters through 10-fold Cross-Validation:\n")
      print(model_nb_HVKP$bestTune)
      
      ##Prediction on test_strat_HVKP (Naive Bayes Model)----------------------------------
      test_strat_HVKP$pred_nb1_HVKP <- predict(model_nb_HVKP, test_strat_HVKP)
      
      # Model evaluation
      confm_nb1_HVKP <- table(actual = test_strat_HVKP$HVKP, prediction = test_strat_HVKP$pred_nb1_HVKP)
      
      metrics_nb1_HVKP<- compute_metrics(confm_nb1_HVKP)
      cat("Performance of NBC model based on model-agnostic approach results:\n")
      metrics_nb1_HVKP
    }
    
    else if (model == "NBC-Wald test"){
      
      train_control <- trainControl(
        method = "cv",
        number = 10
      )
      
      search_grid_HVKP <- expand.grid(
        usekernel = c(FALSE,TRUE),
        fL = seq(0,5,0.5),
        adjust = seq(0,5,0.5)
      )
      
      model_nb2_HVKP <- train(HVKP ~ Imipenem + Meropenem + Aztreonam + Ceftriaxone + Ceftazidime,
                              data = train_strat_HVKP,
                              method = "nb",
                              metric = "Accuracy",
                              trControl = train_control,
                              tuneGrid = search_grid_HVKP,
      )
      
      cat("Tuned hyper-parameters through 10-fold Cross-Validation:\n")
      print(model_nb2_HVKP$bestTune)
      
      ##Prediction on test_strat_HVKP (Naive Bayes Model)
      test_strat_HVKP$pred_nb2_HVKP <- predict(model_nb2_HVKP, test_strat_HVKP)
      
      ##confusion matrix(Naive Bayes Model)
      confm_nb2_HVKP <- table(actual = test_strat_HVKP$HVKP, prediction = test_strat_HVKP$pred_nb2_HVKP)
      
      metrics_nb2_HVKP<- compute_metrics(confm_nb2_HVKP)
      
      cat("Performance of NBC model based on the Wald test results:\n")
      metrics_nb2_HVKP
    }
    
    else if (model == "NBC-Simplicity Principle_MEM"){
      
      train_control <- trainControl(
        method = "cv",
        number = 10
      )
      
      search_grid_HVKP <- expand.grid(
        usekernel = c(FALSE,TRUE),
        fL = seq(0,5,0.5),
        adjust = seq(0,5,0.5)
      )
      
      model_nb3_HVKP <- train(HVKP ~ Meropenem,
                              data = train_strat_HVKP,
                              method = "nb",
                              metric = "Accuracy",
                              trControl = train_control,
                              tuneGrid = search_grid_HVKP,
      )
      
      cat("Tuned hyper-parameters through 10-fold Cross-Validation:\n")
      print(model_nb3_HVKP$bestTune)
      
      ##Prediction on test_strat_HVKP (Naive Bayes Model)
      test_strat_HVKP$pred_nb3_HVKP <- predict(model_nb3_HVKP, test_strat_HVKP)
      
      ##confusion matrix(Naive Bayes Model)
      confm_nb3_HVKP <- table(actual = test_strat_HVKP$HVKP, prediction = test_strat_HVKP$pred_nb3_HVKP)
      
      metrics_nb3_HVKP<- compute_metrics(confm_nb3_HVKP)
      
      cat("Performance of NBC model based on the Simplicity Principle (MEM) results:\n")
      metrics_nb3_HVKP
    }
    
    else if (model == "LDA-Chi Squared test"){
      
      model_ldadalex_HVKP <- train(HVKP ~ Imipenem + Gentamicin + Ciprofloxacin + Meropenem + Aztreonam +
                                     Cefotaxime + Cefepime + Ceftazidime + Ceftriaxone,
                                   data = train_strat_HVKP,
                                   method = "lda",
                                   metric = "Accuracy")
      
      cat("LDA model based on Chi-Squared test results:\n")
      print(model_ldadalex_HVKP)
      
      #Prediction on test_strat_HVKP (LDA Model)-----------------------------------------------
      test_strat_HVKP$pred_lda_HVKP <- predict(model_ldadalex_HVKP, test_strat_HVKP)
      
      # Model Evaluation-------------------------------------------------
      confm_lda_HVKP <- table(actual = test_strat_HVKP$HVKP, prediction = test_strat_HVKP$pred_lda_HVKP)
      
      metrics_lda_HVKP<- compute_metrics(confm_lda_HVKP)
      
      cat("Performance of LDA model based on Chi-Squared test results:\n")
      metrics_lda_HVKP
    }
    
    else if (model == "LDA-model agnostic"){
      model_lda_HVKP <- train(HVKP ~ Meropenem,
                              data = train_strat_HVKP,
                              method = "lda",
                              metric = "Accuracy")
      
      cat("LDA model based on model-agnostic approach:\n")
      print(model_lda_HVKP)
      
      #Prediction on test_strat_HVKP (LDA Model)-----------------------------------------------
      test_strat_HVKP$pred_lda1_HVKP <- predict(model_lda_HVKP, test_strat_HVKP)
      
      # Model evaluation
      confm_lda1_HVKP <- table(actual = test_strat_HVKP$HVKP, prediction = test_strat_HVKP$pred_lda1_HVKP)
      metrics_lda1_HVKP <- compute_metrics(confm_lda1_HVKP)
      cat("Performance of LDA model based on model-agnostic approach:\n")
      metrics_lda1_HVKP
    }
    
    else if (model == "LDA-Wald test"){
      model_lda2_HVKP <- train(HVKP ~ Imipenem + Meropenem + Aztreonam + Ceftriaxone + Ceftazidime,
                               data = train_strat_HVKP,
                               method = "lda",
                               metric = "Accuracy")
                            
      cat("LDA model based on the Wald test results:\n")
      print(model_lda2_HVKP)
      
      test_strat_HVKP$pred_lda2_HVKP <- predict(model_lda2_HVKP, test_strat_HVKP)
      
      confm_lda2_HVKP <- table(actual = test_strat_HVKP$HVKP, prediction = test_strat_HVKP$pred_lda2_HVKP)
      
      metrics_lda2_HVKP<- compute_metrics(confm_lda2_HVKP)
      
      cat("Performance of LDA model based on the Wald test results:\n")
      metrics_lda2_HVKP
    }
    
    else if (model == "SVM-Sig-Chi Squared test"){
      # Fit Support Vector Machine model to train data set
      # This code performs a grid search for tuning the hyper-parameters of the SVM model using the tune function
      # from the e1071 package. It then fits the SVM model with the specified hyper-parameters using the svm function.
      # Model 1
      set.seed(1234)
      tune_out2_HVKP <- e1071::tune("svm", HVKP ~ Imipenem + Gentamicin + Ciprofloxacin + Meropenem + Aztreonam +
                                      Cefotaxime + Cefepime + Ceftazidime + Ceftriaxone 
                                    , data = train_strat_HVKP, kernel = "sigmoid", tunecontrol=tune.control(cross=10),
                                    ranges = list(cost = seq(0.1, 5, length = 20)
                                                  , gamma = seq(0.1,5, length = 20)))
      
      cat("Tuned hyper-parameters through 10-Fold Cross Validation:\n")
      print(summary(tune_out2_HVKP))
      
      model_svm2_HVKP = svm(HVKP ~ Imipenem + Gentamicin + Ciprofloxacin + Meropenem + Aztreonam +
                              Cefotaxime + Cefepime + Ceftazidime + Ceftriaxone, data = train_strat_HVKP, 
                           probability = T, kernel = "sigmoid", cost = tune_out2_HVKP$best.parameters$cost,
                           gamma = tune_out2_HVKP$best.parameters$gamma )
      
      cat("Sigmoid SVM based on Chi-Squared test results:\n")
      print(model_svm2_HVKP)
      
      #Prediction on test_strat_HVKP (SVM)----------------------------------------------------
      test_strat_HVKP$pred_svm_HVKP <- predict(model_svm2_HVKP, test_strat_HVKP)
      
      # Model evaluation
      confm_svm_HVKP <- table(actual = test_strat_HVKP$HVKP, prediction = test_strat_HVKP$pred_svm_HVKP)
      
      metrics_svm1_HVKP_sig<- compute_metrics(confm_svm_HVKP)
      
      cat("Performance of SVM-Sigmoid model based on Chi-Squared test results:\n")
      metrics_svm1_HVKP_sig
    }
    
    else if (model=="SVM-Sig-model agnostic"){
      
      set.seed(1234)
      tune_out21_HVKP <- e1071::tune("svm", HVKP ~ Meropenem,
                                     data = train_strat_HVKP,
                                     kernel = "sigmoid", tunecontrol=tune.control(cross=10),
                                     ranges = list(cost = seq(0.1, 5, length = 25)
                                                   , gamma = seq(0, 5, length = 25)))
      
      cat("Tuned hyper-parameters through 10-Fold Cross Validation:\n")
      print(summary(tune_out21_HVKP))
      
      model_svm21_HVKP = svm(HVKP ~ Meropenem, data = train_strat_HVKP, 
                            probability = T, kernel = "sigmoid", cost = tune_out21_HVKP$best.parameters$cost,
                            gamma = tune_out21_HVKP$best.parameters$gamma)
      
      cat("Sigmoid SVM based on model-agnostic approach:\n")
      print(model_svm21_HVKP)
      
      test_strat_HVKP$pred_svm21_HVKP <- predict(model_svm21_HVKP, test_strat_HVKP)
      
      confm_svm21_HVKP <- table(actual = test_strat_HVKP$HVKP, prediction = test_strat_HVKP$pred_svm21_HVKP)
      
      metrics_svm2_HVKP_sig<- compute_metrics(confm_svm21_HVKP)
      
      cat("Performance of SVM-Sigmoid model based on model-agnostic approach:\n")
      metrics_svm2_HVKP_sig
    }
    
    else if (model=="SVM-Wald test"){
      set.seed(1234)
      tune_out3_HVKP <- e1071::tune("svm", HVKP ~ Imipenem + Meropenem + Aztreonam + Ceftriaxone + Ceftazidime
                                    , data = train_strat_HVKP, kernel = "sigmoid", tunecontrol=tune.control(cross=10),
                                    ranges = list(cost = seq(0.1, 5, length = 25)
                                                  , gamma = seq(0.1, 1, length = 25)))
      
      cat("Tuned hyper-parameters through 10-Fold Cross Validation:\n")
      print(summary(tune_out3_HVKP))
      
      model_svm3_HVKP = svm(HVKP ~ Imipenem + Meropenem + Aztreonam + Ceftriaxone + Ceftazidime, data = train_strat_HVKP, 
                           probability = T, kernel = "sigmoid", cost = tune_out3_HVKP$best.parameters$cost,
                           gamma = tune_out3_HVKP$best.parameters$gamma )
      
      cat("Sigmoid SVM based on the Wald test:\n")
      print(model_svm3_HVKP)
      
      #Prediction on test_strat_HVKP (SVM)
      test_strat_HVKP$pred_svm3_HVKP <- predict(model_svm3_HVKP, test_strat_HVKP)
      
      # Model evaluation
      confm_svm3_HVKP <- table(actual = test_strat_HVKP$HVKP, prediction = test_strat_HVKP$pred_svm3_HVKP)
      
      metrics_svm3_HVKP_sig<- compute_metrics(confm_svm3_HVKP)
      
      cat("Performance of SVM-Sigmoid model based on the Wald test results:\n")
      metrics_svm3_HVKP_sig
    } 
    
    else if (model=="SVM-Poly-Chi Squared test"){
      set.seed(1234)
      tune_out4_HVKP <- e1071::tune("svm", HVKP ~ Imipenem + Gentamicin + Ciprofloxacin + Meropenem + 
                                      Aztreonam + Cefotaxime + Cefepime + Ceftazidime + Ceftriaxone  
                                    , data = train_strat_HVKP, kernel = "polynomial", tunecontrol=tune.control(cross=10),
                                    ranges = list(cost = seq(0.1, 5, length = 20), degree = seq(1,3,1)
                                                  , gamma = seq(0.1,5, length = 20)))
      
      cat("Tuned hyper-parameters through 10-Fold Cross Validation:\n")
      print(summary(tune_out4_HVKP))
      
      # Tuned hyper-parameters are placed in the model.
      model_svm4_HVKP = svm(HVKP ~ Imipenem + Gentamicin + Ciprofloxacin + Meropenem + 
                              Aztreonam + Cefotaxime + Cefepime + Ceftazidime + Ceftriaxone
                           , data = train_strat_HVKP, probability = T, 
                           kernel = "polynomial", cost = tune_out4_HVKP$best.parameters$cost, 
                           degree = tune_out4_HVKP$best.parameters$degree , gamma = tune_out4_HVKP$best.parameters$gamma )
      
      cat("Polynomial SVM based on Chi-Squared test results:\n")
      print(model_svm4_HVKP)
      
      #Prediction on test (SVM 2)----------------------------------------------------
      test_strat_HVKP$pred_svm4_HVKP <- predict(model_svm4_HVKP, test_strat_HVKP)
      
      # Model evaluation
      confm_svm4_HVKP <- table(actual = test_strat_HVKP$HVKP, prediction = test_strat_HVKP$pred_svm4_HVKP)
      
      metrics_svm_HVKP<- compute_metrics(confm_svm4_HVKP)
      
      cat("Performance of SVM-Polynomial model based on Chi-Squared test results:\n")
      metrics_svm_HVKP
    }
    else if (model=="SVM-Poly-model agnostic"){
      
      set.seed(1234)
      tune_out41_HVKP <- e1071::tune("svm", HVKP ~ Meropenem  + Ceftriaxone,
                                     data = train_strat_HVKP,
                                     kernel = "polynomial", tunecontrol=tune.control(cross=10),
                                     ranges = list(cost = seq(0.1, 5, length = 20), degree = seq(1,3,1)
                                                   , gamma = seq(0, 5, length = 20)))
      
      cat("Tuned hyper-parameters through 10-Fold Cross Validation:\n")
      print(summary(tune_out41_HVKP))
      
      model_svm41_HVKP = svm(HVKP ~ Meropenem  + Ceftriaxone, data = train_strat_HVKP, 
                            probability = T, kernel = "polynomial", cost = tune_out41_HVKP$best.parameters$cost,
                            degree = tune_out41_HVKP$best.parameters$degree , gamma = tune_out41_HVKP$best.parameters$gamma )
      
      cat("Polynomial SVM based on model-agnostic approach:\n")
      print(model_svm41_HVKP)
      
      #Prediction on test_strat_HVKP (SVM)----------------------------------------------------
      test_strat_HVKP$pred_svm41_HVKP <- predict(model_svm41_HVKP, test_strat_HVKP)
      
      # Model evaluation
      confm_svm41_HVKP <- table(actual = test_strat_HVKP$HVKP, prediction = test_strat_HVKP$pred_svm41_HVKP)
      
      metrics_svm41_HVKP<- compute_metrics(confm_svm41_HVKP)
      
      cat("Performance of SVM-Polynomial model based on model-agnostic approach:\n")
      metrics_svm41_HVKP
      
    }
    
    else if (model=="SVM-Poly-Simplicity Principle_MEM"){
      
      set.seed(1234)
      tune_out5_HVKP <- e1071::tune("svm", HVKP ~ Meropenem,
                                    data = train_strat_HVKP,
                                    kernel = "polynomial", tunecontrol=tune.control(cross=10),
                                    ranges = list(cost = seq(0.1, 5, length = 20), degree = seq(1,3,1)
                                                  , gamma = seq(0, 5, length = 20)))
      
      cat("Tuned hyper-parameters through 10-Fold Cross Validation:\n")
      print(summary(tune_out5_HVKP))
      
      model_svm5_HVKP = svm(HVKP ~ Meropenem, data = train_strat_HVKP, 
                             probability = T, kernel = "polynomial", cost = tune_out5_HVKP$best.parameters$cost,
                             degree = tune_out5_HVKP$best.parameters$degree , gamma = tune_out5_HVKP$best.parameters$gamma )
      
      cat("Polynomial SVM based on Simplicity Principle (MEM):\n")
      print(model_svm5_HVKP)
      
      #Prediction on test_strat_HVKP (SVM)----------------------------------------------------
      test_strat_HVKP$pred_svm5_HVKP <- predict(model_svm5_HVKP, test_strat_HVKP)
      
      # Model evaluation
      confm_svm5_HVKP <- table(actual = test_strat_HVKP$HVKP, prediction = test_strat_HVKP$pred_svm5_HVKP)
      
      metrics_svm5_HVKP<- compute_metrics(confm_svm5_HVKP)
      
      cat("Performance of SVM-Polynomial model based on Simplicity Principle (MEM):\n")
      metrics_svm5_HVKP
      
    }
    
    else if (model=="SVM-Poly-Wald test"){
      
      set.seed(1234)
      tune_out6_HVKP <- e1071::tune("svm", HVKP ~ Imipenem + Meropenem + Aztreonam + Ceftriaxone + Ceftazidime,
                                    data = train_strat_HVKP,
                                    kernel = "polynomial", tunecontrol=tune.control(cross=10),
                                    ranges = list(cost = seq(0.1, 5, length = 20), degree = seq(1,3,1)
                                                  , gamma = seq(0, 5, length = 20)))
      
      cat("Tuned hyper-parameters through 10-Fold Cross Validation:\n")
      print(summary(tune_out6_HVKP))
      
      model_svm6_HVKP = svm(HVKP ~ Imipenem + Meropenem + Aztreonam + Ceftriaxone + Ceftazidime, data = train_strat_HVKP, 
                           probability = T, kernel = "polynomial", cost = tune_out6_HVKP$best.parameters$cost,
                           degree = tune_out6_HVKP$best.parameters$degree , gamma = tune_out6_HVKP$best.parameters$gamma )
      
      cat("Polynomial SVM based on the Wald test:\n")
      print(model_svm6_HVKP)
      
      #Prediction on test_strat_HVKP (SVM)----------------------------------------------------
      test_strat_HVKP$pred_svm6_HVKP <- predict(model_svm6_HVKP, test_strat_HVKP)
      
      # Model evaluation
      confm_svm6_HVKP <- table(actual = test_strat_HVKP$HVKP, prediction = test_strat_HVKP$pred_svm6_HVKP)
      
      metrics_svm6_HVKP<- compute_metrics(confm_svm6_HVKP)
      
      cat("Performance of SVM-Polynomial model based on the Wald test:\n")
      metrics_svm6_HVKP
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
      
      split_strat_HVKP  <- initial_split(Bio_Data, prop = 0.8, strata = "HVKP")
      
      train_strat_HVKP  <- training(split_strat_HVKP)
      
      test_strat_HVKP   <- testing(split_strat_HVKP)
      
      # Training CatBoost model based on the selected antibiotics through chi squared test
      grid_HVKP <- expand.grid(
        depth = c(2,3,4,5),
        learning_rate = c(0.1,0.2,0.3),
        l2_leaf_reg = c(1,2,3,4),
        rsm = 1,
        border_count = 1,
        iterations = seq(10,130,1)
      )
      fitControl_HVKP <- trainControl(method = "cv",
                                      number = 10,
                                      classProbs = TRUE)
      set.seed(1234)
      model_catboost_HVKP <- train(x = train_strat_HVKP[,c("Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                           "Ceftazidime", "Ciprofloxacin", "Meropenem", "Ceftriaxone", "Aztreonam")],
                                   y = make.names(train_strat_HVKP[,"HVKP"]),
                                   maximize = TRUE,
                                   method = catboost.caret, metric = "Accuracy", 
                                   tuneGrid =  grid_HVKP,
                                   trControl = fitControl_HVKP)
      
      cat("Tuned hyper-parameters through 10-Fold Cross Validation:\n")
      print(model_catboost_HVKP$bestTune)  
      
      #Prediction on test_strat_HVKP (CatBoost)----------------------------------------------------
      test_strat_HVKP$cat_HVKP = predict(model_catboost_HVKP, test_strat_HVKP)
      
      # Model evaluation
      cat_test_HVKP <- table(actual = test_strat_HVKP$HVKP, prediction = test_strat_HVKP$cat_HVKP)
      
      metrics_cat_HVKP<- compute_metrics(cat_test_HVKP)
      
      cat("Performance of CatBoost model based on Chi-Squared test results:\n")
      metrics_cat_HVKP
    }
    
    else if (model=="CatBoost-model agnostic"){
      # CatBoost Models
      Bio_Data <- as.data.frame(Bio_Data)
      
      #Convert categorical variables to factor
      categorical_var <- c(colnames(Bio_Data))
      
      Bio_Data[, categorical_var[-1]] <- lapply(Bio_Data[, categorical_var[-1]], factor)
      
      levels(Bio_Data$Ampicillin) <- c(levels(Bio_Data$Ampicillin),0)
      
      #Divide Data set into train and test
      set.seed(123)
      
      split_strat_HVKP  <- initial_split(Bio_Data, prop = 0.8, strata = "HVKP")
      
      train_strat_HVKP  <- training(split_strat_HVKP)
      
      test_strat_HVKP   <- testing(split_strat_HVKP)
      
      grid_HVKP <- expand.grid(
        depth = c(2,3),
        learning_rate = c(0.1,0.2,0.3),
        l2_leaf_reg = c(1,2,3),
        rsm = 1,
        border_count = 1,
        iterations = seq(10,60,1)
      )
      fitControl_HVKP <- trainControl(method = "cv",
                                      number = 10,
                                      classProbs = TRUE)
      model_catboost1_HVKP <- train(x = train_strat_HVKP[,c("Cefepime", "Ceftazidime", "Meropenem")],
                                    y = make.names(train_strat_HVKP[,"HVKP"]),
                                    maximize = TRUE,
                                    method = catboost.caret, metric = "Accuracy", 
                                    tuneGrid =  grid_HVKP,
                                    trControl = fitControl_HVKP)
      
      cat("Tuned hyper-parameters through 10-Fold Cross Validation:\n")
      print(model_catboost1_HVKP$bestTune)  
      
      test_strat_HVKP$cat1_HVKP = predict(model_catboost1_HVKP, test_strat_HVKP)
      
      # Model Evaluation
      confm_cat1_HVKP <- table(actual = test_strat_HVKP$HVKP, prediction = test_strat_HVKP$cat1_HVKP)
      
      metrics_cat1_HVKP<- compute_metrics(confm_cat1_HVKP)
      
      cat("Performance of CatBoost model based on model-agnostic approach results:\n")
      metrics_cat1_HVKP
    }
    
    else if (model=="CatBoost-Wald test"){
      # CatBoost Models
      Bio_Data <- as.data.frame(Bio_Data)
      
      #Convert categorical variables to factor
      categorical_var <- c(colnames(Bio_Data))
      
      Bio_Data[, categorical_var[-1]] <- lapply(Bio_Data[, categorical_var[-1]], factor)
      
      levels(Bio_Data$Ampicillin) <- c(levels(Bio_Data$Ampicillin),0)
      
      #Divide Data set into train and test
      set.seed(123)
      
      split_strat_HVKP  <- initial_split(Bio_Data, prop = 0.8, strata = "HVKP")
      
      train_strat_HVKP  <- training(split_strat_HVKP)
      
      test_strat_HVKP   <- testing(split_strat_HVKP)
      
      grid2_HVKP <- expand.grid(
        depth = c(3,4,5,6,7),
        learning_rate = c(0.1,0.2,0.3,0.4,0.5),
        l2_leaf_reg = c(2,3,4),
        rsm = 1,
        border_count = 1,
        iterations = seq(90,120,1)
      )
      fitControl2_HVKP <- trainControl(method = "cv",
                                       number = 10,
                                       classProbs = TRUE)
      model_catboost2_HVKP <- train(x = train_strat_HVKP[,c("Aztreonam","Imipenem","Meropenem","Ceftazidime","Ceftriaxone")],
                                    y = make.names(train_strat_HVKP[["HVKP"]]),
                                    maximize = TRUE,
                                    method = catboost.caret, metric = "Accuracy", 
                                    tuneGrid =  grid2_HVKP,
                                    trControl = fitControl2_HVKP)
      
      
      cat("Tuned hyper-parameters through 10-Fold Cross Validation:\n")
      print(model_catboost2_HVKP$bestTune)  
      
      test_strat_HVKP$cat_HVKP2 = predict(model_catboost2_HVKP, test_strat_HVKP)
      # Model evaluation
      confm_cat2_HVKP <- table(actual = test_strat_HVKP$HVKP, prediction = test_strat_HVKP$cat_HVKP2)
      
      metrics_cat2_HVKP<- compute_metrics(confm_cat2_HVKP)
      
      cat("Performance of CatBoost model based on the Wald test results:\n")
      metrics_cat2_HVKP
    }
    
    else if (model=="CatBoost-Simplicity Principle_MEM"){
      # CatBoost Models
      Bio_Data <- as.data.frame(Bio_Data)
      
      #Convert categorical variables to factor
      categorical_var <- c(colnames(Bio_Data))
      
      Bio_Data[, categorical_var[-1]] <- lapply(Bio_Data[, categorical_var[-1]], factor)
      
      levels(Bio_Data$Ampicillin) <- c(levels(Bio_Data$Ampicillin),0)
      
      #Divide Data set into train and test
      set.seed(123)
      
      split_strat_HVKP  <- initial_split(Bio_Data, prop = 0.8, strata = "HVKP")
      
      train_strat_HVKP  <- training(split_strat_HVKP)
      
      test_strat_HVKP   <- testing(split_strat_HVKP)
      
      grid3_HVKP <- expand.grid(
        depth = c(2,3,4),
        learning_rate = c(0.1,0.2,0.3),
        l2_leaf_reg = c(1,2,3),
        rsm = 1,
        border_count = 1,
        iterations = seq(10,200,10)
      )
      
      fitControl3_HVKP <- trainControl(method = "cv",
                                       number = 10,
                                       classProbs = TRUE)
      
      model_catboost3_HVKP <- train(x = train_strat_HVKP["Meropenem"],
                                    y = make.names(train_strat_HVKP[["HVKP"]]),
                                    maximize = TRUE,
                                    method = catboost.caret, metric = "Accuracy", 
                                    tuneGrid =  grid3_HVKP,
                                    trControl = fitControl3_HVKP)
      
      
      cat("Tuned hyper-parameters through 10-Fold Cross Validation:\n")
      print(model_catboost3_HVKP$bestTune)  
      
      test_strat_HVKP$cat_HVKP3 = predict(model_catboost3_HVKP, test_strat_HVKP)
      # Model evaluation
      confm_cat3_HVKP <- table(actual = test_strat_HVKP$HVKP, prediction = test_strat_HVKP$cat_HVKP3)
      
      metrics_cat3_HVKP<- compute_metrics(confm_cat3_HVKP)
      
      cat("Performance of CatBoost model based on the Simplicity Principle (MEM):\n")
      metrics_cat3_HVKP
    }
  }
  
  else if (task == 'feature_selectio') {
    # Task: Feature Selection
    # You can add flags for specific feature selection methods to use, e.g., -m "Chi-Squared test" or -m "Wald test".
    # In order to investigate results of feature selection methods, you should specify "task" to "feature_selection"
    # then you can choose a specific method by -m, e.g., if you want to see results of model-agnostic approach in LR 
    # algorithm, use -t "feature_selection" -m "model-agnostic-LR"
    
    if (method == "Chi-Squared test"){
      # Function to perform chi-square test and print p-value
      perform_chi_square <- function(data, variable) {
        result <- chisq.test(table(data[[variable]], data$HVKP))
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
    model_Log1_HVKP <- glm(HVKP ~ Imipenem + Gentamicin + Ciprofloxacin + Meropenem + Aztreonam +
                             Cefotaxime + Cefepime + Ceftazidime + Ceftriaxone
                           , family = "binomial", data = train_strat_HVKP)
    
    cat("P-values show Wald test results:\n")
    print(summary(model_Log1_HVKP))
  }
  else if (method == "model agnostic-LR"){
    # Feature Selection by DALEX Package--------------------------------------
    model_Log1_HVKP <- glm(HVKP ~ Imipenem + Gentamicin + Ciprofloxacin + Meropenem + Aztreonam +
                             Cefotaxime + Cefepime + Ceftazidime + Ceftriaxone
                           , family = "binomial", data = train_strat_HVKP)
    
    explained_glm_HVKP <- explain(model = model_Log1_HVKP, data=train_strat_HVKP[,c(3:9,11:12)], variables = colnames(train_strat_HVKP[,c(3:9,11:12)]),
                                  y=as.vector(as.numeric(train_strat_HVKP$HVKP)), label = "LR", 
                                  type = "classification")
    #50 permuatation
    fi_glm_HVKP = variable_importance(explained_glm_HVKP, B=50 ,variables = colnames(train_strat_HVKP[,c(3:9,11:12)]), loss_function = loss_root_mean_square, type = "raw")
    
    cat("Root mean squared dropout loss due to a feature elimination through model-agnostic approach in LR:\n")
    print(fi_glm_HVKP)
  }
  
  else if (method == "model agnostic-NBC"){
    train_control <- trainControl(
      method = "cv",
      number = 10
    )
    
    search_grid_HVKP <- expand.grid(
      usekernel = c(FALSE,TRUE),
      fL = seq(0,5,0.5),
      adjust = seq(0,5,0.5)
    )
    
    model_nbdalex_HVKP <- train(HVKP ~ Imipenem + Gentamicin + Ciprofloxacin + Meropenem + Aztreonam +
                                  Cefotaxime + Cefepime + Ceftazidime + Ceftriaxone,
                                data = train_strat_HVKP,
                                method = "nb",
                                metric = "Accuracy",
                                trControl = train_control,
                                tuneGrid = search_grid_HVKP
    )
    # Feature selection through model-agnostic approach (DALEX)
    explained_nb_HVKP <- explain(model = model_nbdalex_HVKP, data=train_strat_HVKP[,c(3:9,11:12)], variables = colnames(train_strat_HVKP[,c(3:9,11:12)]),
                                 y=as.vector(as.numeric(train_strat_HVKP$HVKP)), label = "NBC", 
                                 type = "classification")
    
    # The 'variable_importance' function is then used to calculate the variable importance based on the
    # explained model. The 'B' parameter is set to 50 for the number of permutations, and the loss function
    # is set to "loss_root_mean_square". The resulting variable importance is stored in 'fitnb_TEM'
    fi_nb_HVKP = variable_importance(explained_nb_HVKP, B=50 ,variables = colnames(train_strat_HVKP[,c(3:9,11:12)]), loss_function = loss_root_mean_square, type = "raw")
    
    cat("Root mean squared dropout loss due to a feature elimination through model-agnostic approach in NBC:\n")
    print(fi_nb_HVKP)
  }
  
  else if (method == "model agnostic-LDA"){
    
    model_ldadalex_HVKP <- train(HVKP ~ Imipenem + Gentamicin + Ciprofloxacin + Meropenem + Aztreonam +
                                   Cefotaxime + Cefepime + Ceftazidime + Ceftriaxone,
                                 data = train_strat_HVKP,
                                 method = "lda",
                                 metric = "Accuracy")
    
    explained_lda_HVKP <- explain(model = model_ldadalex_HVKP, data=train_strat_HVKP[,c(3:9,11:12)], variables = colnames(train_strat_HVKP[,c(3:9,11:12)]),
                                  y=as.vector(as.numeric(train_strat_HVKP$HVKP)), label = "LDA", 
                                  type = "classification")
    
    fi_lda_HVKP = variable_importance(explained_lda_HVKP, B=50 ,variables = colnames(train_strat_HVKP[,c(3:9,11:12)]), loss_function = loss_root_mean_square, type = "raw")
    
    cat("Root mean squared dropout loss due to a feature elimination through model-agnostic approach in LDA:\n")
    print(fi_lda_HVKP)
  }
  
  else if (method == "model agnostic-Sig SVM"){
    set.seed(1234)
    tune_out2_HVKP <- e1071::tune("svm", HVKP ~ Imipenem + Gentamicin + Ciprofloxacin + Meropenem + Aztreonam +
                                    Cefotaxime + Cefepime + Ceftazidime + Ceftriaxone 
                                  , data = train_strat_HVKP, kernel = "sigmoid", tunecontrol=tune.control(cross=10),
                                  ranges = list(cost = seq(0.1, 5, length = 20)
                                                , gamma = seq(0.1,5, length = 20)))
    
    model_svm2_HVKP = svm(HVKP ~ Imipenem + Gentamicin + Ciprofloxacin + Meropenem + Aztreonam +
                            Cefotaxime + Cefepime + Ceftazidime + Ceftriaxone, data = train_strat_HVKP, 
                          probability = T, kernel = "sigmoid", cost = tune_out2_HVKP$best.parameters$cost,
                          gamma = tune_out2_HVKP$best.parameters$gamma )
    
    explained_SVM2_HVKP <- explain(model = model_svm2_HVKP, data=train_strat_HVKP[,c(3:9,11:12)], variables = colnames(train_strat_HVKP[,c(3:9,11:12)]),
                                   y=as.vector(as.numeric(train_strat_HVKP$HVKP)), label = "SVM-Sigmoid", 
                                   type = "classification")
    #50 permuatation
    fi_SVM2_HVKP = variable_importance(explained_SVM2_HVKP, B=50 ,variables = colnames(train_strat_HVKP[,c(3:9,11:12)]),loss_function = loss_root_mean_square, type = "raw")
    
    cat("Root mean squared dropout loss due to a feature elimination through model-agnostic approach in sigmoid kernel of SVM:\n")
    print(fi_SVM2_HVKP)
  }
  
  else if (method == "model agnostic-Poly SVM"){
    
    set.seed(1234)
    tune_out4_HVKP <- e1071::tune("svm", HVKP ~ Imipenem + Gentamicin + Ciprofloxacin + Meropenem + 
                                    Aztreonam + Cefotaxime + Cefepime + Ceftazidime + Ceftriaxone  
                                  , data = train_strat_HVKP, kernel = "polynomial", tunecontrol=tune.control(cross=10),
                                  ranges = list(cost = seq(0.1, 5, length = 20), degree = seq(1,3,1)
                                                , gamma = seq(0.1,5, length = 20)))
    
    model_svm4_HVKP = svm(HVKP ~ Imipenem + Gentamicin + Ciprofloxacin + Meropenem + 
                            Aztreonam + Cefotaxime + Cefepime + Ceftazidime + Ceftriaxone
                          , data = train_strat_HVKP, probability = T, 
                          kernel = "polynomial", cost = tune_out4_HVKP$best.parameters$cost, 
                          degree = tune_out4_HVKP$best.parameters$degree , gamma = tune_out4_HVKP$best.parameters$gamma )
    
    explained_SVM4_HVKP <- explain(model = model_svm4_HVKP, data=train_strat_HVKP[,c(3:9,11:12)], variables = colnames(train_strat_HVKP[,c(3:9,11:12)]),
                                   y=as.vector(as.numeric(train_strat_HVKP$HVKP))-1, label = "SVM-Polynomial", 
                                   type = "classification")
    
    #50 permuatation
    fi_SVM4_HVKP = variable_importance(explained_SVM4_HVKP, B=50 ,variables = colnames(train_strat_HVKP[,c(3:9,11:12)]),loss_function = loss_root_mean_square, type = "raw")
    
    cat("Root mean squared dropout loss due to a feature elimination through model-agnostic approach in polynomial kernel of SVM:\n")
    print(fi_SVM4_HVKP)
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
    
    split_strat_HVKP  <- initial_split(Bio_Data, prop = 0.8, strata = "HVKP")
    
    train_strat_HVKP  <- training(split_strat_HVKP)
    
    test_strat_HVKP   <- testing(split_strat_HVKP)
    
    # Training CatBoost model based on the selected antibiotics through chi squared test
    grid_HVKP <- expand.grid(
      depth = c(2,3,4,5),
      learning_rate = c(0.1,0.2,0.3),
      l2_leaf_reg = c(1,2,3,4),
      rsm = 1,
      border_count = 1,
      iterations = seq(10,130,1)
    )
    fitControl_HVKP <- trainControl(method = "cv",
                                    number = 10,
                                    classProbs = TRUE)
    set.seed(1234)
    model_catboost_HVKP <- train(x = train_strat_HVKP[,c("Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                         "Ceftazidime", "Ciprofloxacin", "Meropenem", "Ceftriaxone", "Aztreonam")],
                                 y = make.names(train_strat_HVKP[,"HVKP"]),
                                 maximize = TRUE,
                                 method = catboost.caret, metric = "Accuracy", 
                                 tuneGrid =  grid_HVKP,
                                 trControl = fitControl_HVKP)
    
    expsda_HVKP = explain(model = model_catboost_HVKP, data=train_strat_HVKP[,c("Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                                "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")], variables = colnames(train_strat_HVKP[,c("Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                                                                                                                                                   "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")]),
                          y=as.vector(as.numeric(train_strat_HVKP$HVKP))-1, label = "CatBoost", 
                          type = "classification")
    
    fitcat_HVKP = variable_importance(expsda_HVKP, B=50 ,variables = colnames(train_strat_HVKP[,c("Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                                                  "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")]), loss_function = loss_root_mean_square, type = "raw" )
    
    cat("Root mean squared dropout loss due to a feature elimination through model-agnostic approach in CatBoost:\n")
    print(fitcat_HVKP)
  }
  else {
    cat("Invalid task. Please specify one of the following tasks: feature_selectio, model_evaluation\n")
    
  } 
}
if (is.null(task)){
  cat("You should specify a task\n")
}
