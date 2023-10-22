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

# Calculate the proportion of IUC_gene

proportion_IUC <- table(Bio_Data$IUC) %>% prop.table()

cat("Proportion of IUC gene:\n")

print(proportion_IUC)
# Response variable is imbalanced (0 : 72% , 1 : 28%) so we're using stratified sampling

# Divide the data set into train and test sets through stratified sampling method

set.seed(6546)

split_strat_IUC  <- initial_split(Bio_Data, prop = 0.8, strata = "IUC")

train_strat_IUC  <- training(split_strat_IUC)

test_strat_IUC   <- testing(split_strat_IUC)

# Display the dimensions of the train and test sets
cat("train set dimensions:", dim(train_strat_IUC), "\n")
cat("test set dimensions:", dim(test_strat_IUC), "\n")


if (!is.null(task)) {
  if (task == "model_evaluation") {

    if (model == "LR-Chi Squared test"){
      # Logistic regression (Model 1)--------------------------------------------------
      model_Log1_IUC <- glm(IUC ~ Imipenem + Gentamicin + Ciprofloxacin + Meropenem + Aztreonam +
                              Cefotaxime + Cefepime + Ceftazidime + Ceftriaxone
                            , family = "binomial", data = train_strat_IUC)
      
      cat("Summary of LR model based on Chi Squared test including Null deviance, Residual deviance, and AIC, etc.:\n")
      print(summary(model_Log1_IUC))
      
      test_strat_IUC$probs1_IUC<- predict(model_Log1_IUC, test_strat_IUC, type = "response")
      
      test_strat_IUC$pred_logreg1_IUC <- ifelse(test_strat_IUC$probs1_IUC >= 0.5, 1, 0)
      
      # Model 1 evaluation
      confm_logreg1_IUC <- table(actual = test_strat_IUC$IUC, prediction = test_strat_IUC$pred_logreg1_IUC)
      
      metrics_logreg1_IUC <- compute_metrics(confm_logreg1_IUC)
      
      cat("Performance of LR model based on Chi-Squared test results:\n")
      metrics_logreg1_IUC
    }
    else if (model == "LR-model agnostic"){
      model_Log2_IUC <- glm(IUC ~  Meropenem 
                            , family = "binomial", data = train_strat_IUC)
      
      cat("Summary of LR model based on the model-agnostic approach including Null deviance, Residual deviance, and AIC, etc.:\n")
      print(summary(model_Log2_IUC))
      
      ##Prediction on test_strat_IUC (Model 2)--------------------------------------------------
      test_strat_IUC$probs2_IUC <- predict(model_Log2_IUC, test_strat_IUC, type = "response")
      
      test_strat_IUC$pred_logreg2_IUC <- ifelse(test_strat_IUC$probs2_IUC >= 0.5, 1, 0)
      
      # Model 2 evaluation
      confm_logreg2_IUC <- table(actual = test_strat_IUC$IUC, prediction = test_strat_IUC$pred_logreg2_IUC)
      
      metrics_logreg2_IUC <- compute_metrics(confm_logreg2_IUC)
      
      cat("Performance of LR model based on the model-agnostic approach:\n")
      metrics_logreg2_IUC
    }
    
    else if (model == "LR-Wald test"){
      model_Log3_IUC <- glm(IUC ~ Imipenem + Meropenem + Aztreonam + Ceftazidime + Cefepime
                            , family = "binomial", data = train_strat_IUC)
      
      cat("Summary of LR model based on the Wald test including Null deviance, Residual deviance, and AIC, etc.:\n")
      print(summary(model_Log3_IUC))
      
      ##Prediction on test_strat_IUC
      test_strat_IUC$probs3_IUC <- predict(model_Log3_IUC, test_strat_IUC, type = "response")
      
      test_strat_IUC$pred_logreg3_IUC <- ifelse(test_strat_IUC$probs3_IUC >= 0.5, 1, 0)
      
      confm_logreg3_IUC <- table(actual = test_strat_IUC$IUC, prediction = test_strat_IUC$pred_logreg3_IUC)
      
      metrics_logreg3_IUC <- compute_metrics(confm_logreg3_IUC)
      
      cat("Performance of LR model based on the Wald test:\n")
      metrics_logreg3_IUC
    }
    
    else if (model == "NBC-Chi Squared test"){
      train_control <- trainControl(
        method = "cv",
        number = 10
      )
      
      search_grid_IUC <- expand.grid(
        usekernel = c(FALSE,TRUE),
        fL = seq(0,5,0.5),
        adjust = seq(0,5,0.5)
      )
      
      model_nbdalex_IUC <- train(IUC ~ Imipenem + Gentamicin + Ciprofloxacin + Meropenem + Aztreonam +
                                   Cefotaxime + Cefepime + Ceftazidime + Ceftriaxone,
                                 data = train_strat_IUC,
                                 method = "nb",
                                 metric = "Accuracy",
                                 trControl = train_control,
                                 tuneGrid = search_grid_IUC
      )
      
      
      cat("Tuned hyper-parameters through 10-fold Cross-Validation:\n")
      print(model_nbdalex_IUC$bestTune)
      
      ##Prediction on test_strat_IUC (Naive Bayes Model)----------------------------------
      test_strat_IUC$pred_nb_IUC <- predict(model_nbdalex_IUC, test_strat_IUC)
      
      # Model evaluation
      confm_nb_IUC <- table(actual = test_strat_IUC$IUC, prediction = test_strat_IUC$pred_nb_IUC)
      metrics_nb_IUC<- compute_metrics(confm_nb_IUC)
      cat("Performance of NBC model based on Chi-Squared test results:\n")
      metrics_nb_IUC
    }
    else if (model == "NBC-model agnostic"){
      train_control <- trainControl(
        method = "cv",
        number = 10
      )
      
      search_grid_IUC <- expand.grid(
        usekernel = c(FALSE,TRUE),
        fL = seq(0,5,0.5),
        adjust = seq(0,5,0.5)
      )
      
      model_nbdalex1_IUC <- train(IUC ~ Meropenem,
                                  data = train_strat_IUC,
                                  method = "nb",
                                  metric = "Accuracy",
                                  trControl = train_control,
                                  tuneGrid = search_grid_IUC,
      )
      
      cat("Tuned hyper-parameters through 10-fold Cross-Validation:\n")
      print(model_nbdalex1_IUC$bestTune)
      
      ##Prediction on test_strat_IUC (Naive Bayes Model)----------------------------------
      test_strat_IUC$pred_nb1_IUC <- predict(model_nbdalex1_IUC, test_strat_IUC)
      
      # Model evaluation
      confm_nb1_IUC <- table(actual = test_strat_IUC$IUC, prediction = test_strat_IUC$pred_nb1_IUC)
      
      metrics_nb1_IUC<- compute_metrics(confm_nb1_IUC)
      cat("Performance of NBC model based on model-agnostic approach results:\n")
      metrics_nb1_IUC
    }
    
    else if (model == "NBC-Wald test"){
      
      train_control <- trainControl(
        method = "cv",
        number = 10
      )
      
      search_grid_IUC <- expand.grid(
        usekernel = c(FALSE,TRUE),
        fL = seq(0,5,0.5),
        adjust = seq(0,5,0.5)
      )
      
      model_nb2_IUC <- train( IUC ~ Meropenem + Aztreonam + Ceftazidime + Cefepime + Imipenem,
                              data = train_strat_IUC,
                              method = "nb",
                              metric = "Accuracy",
                              trControl = train_control,
                              tuneGrid = search_grid_IUC,
      )
      
      cat("Tuned hyper-parameters through 10-fold Cross-Validation:\n")
      print(model_nb2_IUC$bestTune)
      
      ##Prediction on test_strat_IUC (Naive Bayes Model)
      test_strat_IUC$pred_nb2_IUC <- predict(model_nb2_IUC, test_strat_IUC)
      
      ##confusion matrix(Naive Bayes Model)
      confm_nb2_IUC <- table(actual = test_strat_IUC$IUC, prediction = test_strat_IUC$pred_nb2_IUC)
      
      metrics_nb2_IUC<- compute_metrics(confm_nb2_IUC)
      
      cat("Performance of NBC model based on the Wald test results:\n")
      metrics_nb2_IUC
    }
    
    else if (model == "LDA-Chi Squared test"){
      
      model_ldadalex_IUC <- train(IUC ~ Cefotaxime + Gentamicin + Imipenem +
                                    Cefepime + Ceftazidime + Ciprofloxacin +
                                    Meropenem + Ceftriaxone + Aztreonam,
                                  data = train_strat_IUC,
                                  method = "lda",
                                  metric = "Accuracy")
      
      cat("LDA model based on Chi-Squared test results:\n")
      print(model_ldadalex_IUC)
      
      #Prediction on test_strat_IUC (LDA Model)-----------------------------------------------
      test_strat_IUC$pred_lda_IUC <- predict(model_ldadalex_IUC, test_strat_IUC)
      
      # Model Evaluation-------------------------------------------------
      confm_lda_IUC <- table(actual = test_strat_IUC$IUC, prediction = test_strat_IUC$pred_lda_IUC)
      
      metrics_lda_IUC<- compute_metrics(confm_lda_IUC)
      
      cat("Performance of LDA model based on Chi-Squared test results:\n")
      metrics_lda_IUC
    }
    else if (model == "LDA-model agnostic"){
      model_ldadalex1_IUC <- train(IUC ~ Meropenem,
                                   data = train_strat_IUC,
                                   method = "lda",
                                   metric = "Accuracy")
      
      cat("LDA model based on model-agnostic approach:\n")
      print(model_ldadalex1_IUC)
      
      #Prediction on test_strat_IUC (LDA Model)-----------------------------------------------
      test_strat_IUC$pred_lda1_IUC <- predict(model_ldadalex1_IUC, test_strat_IUC)
      
      # Model evaluation
      confm_lda1_IUC <- table(actual = test_strat_IUC$IUC, prediction = test_strat_IUC$pred_lda1_IUC)
      metrics_lda1_IUC <- compute_metrics(confm_lda1_IUC)
      cat("Performance of LDA model based on model-agnostic approach:\n")
      metrics_lda1_IUC
    }
    
    else if (model == "LDA-Wald test"){
      model_lda2_IUC <- train(IUC ~ Imipenem + Meropenem + Aztreonam + Ceftazidime + Cefepime,
                              data = train_strat_IUC,
                              method = "lda",
                              metric = "Accuracy")
      
      cat("LDA model based on the Wald test results:\n")
      print(model_lda2_IUC)
      
      test_strat_IUC$pred_lda2_IUC <- predict(model_lda2_IUC, test_strat_IUC)
      
      confm_lda2_IUC <- table(actual = test_strat_IUC$IUC, prediction = test_strat_IUC$pred_lda2_IUC)
      
      metrics_lda2_IUC<- compute_metrics(confm_lda2_IUC)
      
      cat("Performance of LDA model based on the Wald test results:\n")
      metrics_lda2_IUC
    }
    
    else if (model == "SVM-Sig-Chi Squared test"){
      # Fit Support Vector Machine model to train data set
      # This code performs a grid search for tuning the hyper-parameters of the SVM model using the tune function
      # from the e1071 package. It then fits the SVM model with the specified hyper-parameters using the svm function.
      # Model 1
      set.seed(1234)
      tune_out2_IUC <- e1071::tune("svm", IUC ~ Imipenem + Gentamicin + Ciprofloxacin + Meropenem + 
                                     Aztreonam + Cefotaxime + Cefepime + Ceftazidime + Ceftriaxone  
                                   , data = train_strat_IUC, kernel = "sigmoid", tunecontrol=tune.control(cross=10),
                                   ranges = list(cost = seq(0.1, 10, length = 30)
                                                 , gamma = seq(0.1,1, length = 10)))
      
      cat("Tuned hyper-parameters through 10-Fold Cross Validation:\n")
      print(summary(tune_out_IUC))
      
      model_svm2_IUC = svm(IUC ~ Imipenem + Gentamicin + Ciprofloxacin + Meropenem + 
                            Aztreonam + Cefotaxime + Cefepime + Ceftazidime + Ceftriaxone  , data = train_strat_IUC, 
                            probability = T, kernel = "sigmoid", cost = tune_out2_IUC$best.parameters$cost,
                            gamma = tune_out2_IUC$best.parameters$gamma )
      
      cat("Sigmoid SVM based on Chi-Squared test results:\n")
      print(model_svm2_IUC)
      
      #Prediction on test_strat_IUC (SVM)----------------------------------------------------
      test_strat_IUC$pred_svm_IUC <- predict(model_svm2_IUC, test_strat_IUC)
      
      # Model evaluation
      confm_svm_IUC <- table(actual = test_strat_IUC$IUC, prediction = test_strat_IUC$pred_svm_IUC)
      
      metrics_svm1_IUC_sig<- compute_metrics(confm_svm_IUC)
      
      cat("Performance of SVM-Sigmoid model based on Chi-Squared test results:\n")
      metrics_svm1_IUC_sig
    }
    else if (model=="SVM-Sig-model agnostic"){
      
      set.seed(1234)
      tune_out2_IUC <- e1071::tune("svm", IUC ~ Meropenem,
                                   data = train_strat_IUC,
                                   kernel = "sigmoid", tunecontrol=tune.control(cross=10),
                                   ranges = list(cost = seq(0.1, 5, length = 30)
                                                 , gamma = seq(0, 5, length = 25)))
      
      cat("Tuned hyper-parameters through 10-Fold Cross Validation:\n")
      print(summary(tune_out2_IUC))
      
      model_svm21_IUC = svm(IUC ~ Meropenem, data = train_strat_IUC, 
                             probability = T, kernel = "sigmoid", cost = tune_out2_IUC$best.parameters$cost,
                             gamma = tune_out2_IUC$best.parameters$gamma)
      
      cat("Sigmoid SVM based on model-agnostic approach:\n")
      print(model_svm21_IUC)
      
      test_strat_IUC$pred_svm2_IUC <- predict(model_svm21_IUC, test_strat_IUC)
      
      confm_svm2_IUC <- table(actual = test_strat_IUC$IUC, prediction = test_strat_IUC$pred_svm2_IUC)
      
      metrics_svm2_IUC_sig<- compute_metrics(confm_svm2_IUC)
      
      cat("Performance of SVM-Sigmoid model based on model-agnostic approach:\n")
      metrics_svm2_IUC_sig
    }
    
    else if (model=="SVM-Wald test"){
      set.seed(1234)
      tune_out3_IUC <- e1071::tune("svm", IUC ~ Imipenem + Meropenem + Aztreonam + Ceftazidime + Cefepime
                                   , data = train_strat_IUC, kernel = "sigmoid", tunecontrol=tune.control(cross=10),
                                   ranges = list(cost = seq(0.1, 5, length = 25)
                                                 , gamma = seq(0.1, 1, length = 20)))
      
      cat("Tuned hyper-parameters through 10-Fold Cross Validation:\n")
      print(summary(tune_out3_IUC))
      
      model_svm3_IUC = svm(IUC ~ Imipenem + Meropenem + Aztreonam + Ceftazidime + Cefepime, data = train_strat_IUC, 
                             probability = T, kernel = "sigmoid", cost = tune_out3_IUC$best.parameters$cost,
                             gamma = tune_out3_IUC$best.parameters$gamma )
      
      cat("Sigmoid SVM based on the Wald test:\n")
      print(model_svm3_IUC)
      
      #Prediction on test_strat_IUC (SVM)
      test_strat_IUC$pred_svm3_IUC <- predict(model_svm3_IUC, test_strat_IUC)
      
      # Model evaluation
      confm_svm3_IUC <- table(actual = test_strat_IUC$IUC, prediction = test_strat_IUC$pred_svm3_IUC)
      
      metrics_svm3_IUC_sig<- compute_metrics(confm_svm3_IUC)
      
      cat("Performance of SVM-Sigmoid model based on the Wald test results:\n")
      metrics_svm3_IUC_sig
    } 
    
    else if (model=="SVM-Poly-Chi Squared test"){
      set.seed(1234)
      tune_out4_IUC <- e1071::tune("svm", IUC ~ Imipenem + Gentamicin + Ciprofloxacin + Meropenem + 
                                     Aztreonam + Cefotaxime + Cefepime + Ceftazidime + Ceftriaxone  
                                   , data = train_strat_IUC, kernel = "polynomial", tunecontrol=tune.control(cross=10),
                                   ranges = list(cost = seq(0.1, 5, length = 20), degree = seq(1,4,1)
                                                 , gamma = seq(0.1,5, length = 20)))
      
      cat("Tuned hyper-parameters through 10-Fold Cross Validation:\n")
      print(summary(tune_out4_IUC))
      
      # Tuned hyper-parameters are placed in the model.
      model_svm4_IUC = svm(IUC ~ Imipenem + Gentamicin + Ciprofloxacin + Meropenem + Aztreonam +
                             Cefotaxime + Cefepime + Ceftazidime + Ceftriaxone
                             , data = train_strat_IUC, probability = T, 
                             kernel = "polynomial", cost = tune_out4_IUC$best.parameters$cost, 
                             degree = tune_out4_IUC$best.parameters$degree , gamma = tune_out4_IUC$best.parameters$gamma )
      
      cat("Polynomial SVM based on Chi-Squared test results:\n")
      print(model_svm4_IUC)
      
      #Prediction on test (SVM 2)----------------------------------------------------
      test_strat_IUC$pred_svm4_IUC <- predict(model_svm4_IUC, test_strat_IUC)
      
      # Model evaluation
      confm_svm4_IUC <- table(actual = test_strat_IUC$IUC, prediction = test_strat_IUC$pred_svm4_IUC)
      
      metrics_svm_IUC<- compute_metrics(confm_svm4_IUC)
      
      cat("Performance of SVM-Polynomial model based on Chi-Squared test results:\n")
      metrics_svm_IUC
    }
    else if (model=="SVM-Poly-model agnostic"){
      
      set.seed(1234)
      tune_out41_IUC <- e1071::tune("svm", IUC ~ Meropenem ,
                                    data = train_strat_IUC,
                                    kernel = "polynomial", tunecontrol=tune.control(cross=10),
                                    ranges = list(cost = seq(0.1, 5, length = 20), degree = seq(1,3,1)
                                                  , gamma = seq(0, 5, length = 20)))
      
      cat("Tuned hyper-parameters through 10-Fold Cross Validation:\n")
      print(summary(tune_out41_IUC))
      
      model_svm41_IUC = svm(IUC ~ Meropenem, data = train_strat_IUC, 
                                  probability = T, kernel = "polynomial", cost = tune_out41_IUC$best.parameters$cost,
                                  degree = tune_out41_IUC$best.parameters$degree , gamma = tune_out41_IUC$best.parameters$gamma )
      
      cat("Polynomial SVM based on model-agnostic approach:\n")
      print(model_svm41_IUC)
      
      #Prediction on test_strat_IUC (SVM)----------------------------------------------------
      test_strat_IUC$pred_svm41_IUC <- predict(model_svm41_IUC, test_strat_IUC)
      
      # Model evaluation
      confm_svm41_IUC <- table(actual = test_strat_IUC$IUC, prediction = test_strat_IUC$pred_svm41_IUC)
      
      metrics_svm41_IUC<- compute_metrics(confm_svm41_IUC)
      
      cat("Performance of SVM-Polynomial model based on model-agnostic approach:\n")
      metrics_svm41_IUC
      
    }
    
    else if (model=="SVM-Poly-Wald test"){
      
      set.seed(1234)
      tune_out5_IUC <- e1071::tune("svm", IUC ~ Meropenem + Ceftazidime + Cefepime +
                                     Aztreonam + Imipenem,
                                   data = train_strat_IUC,
                                   kernel = "polynomial", tunecontrol=tune.control(cross=10),
                                   ranges = list(cost = seq(0.1, 5, length = 20), degree = seq(1,3,1)
                                                 , gamma = seq(0, 5, length = 20)))
      
      cat("Tuned hyper-parameters through 10-Fold Cross Validation:\n")
      print(summary(tune_out5_IUC))
      
      model_svm5_IUC = svm(IUC ~ Meropenem + Ceftazidime + Cefepime +
                             Aztreonam + Imipenem , data = train_strat_IUC, 
                             probability = T, kernel = "polynomial", cost = tune_out5_IUC$best.parameters$cost,
                             degree = tune_out5_IUC$best.parameters$degree , gamma = tune_out5_IUC$best.parameters$gamma )
      
      cat("Polynomial SVM based on the Wald test:\n")
      print(model_svm5_IUC)
      
      #Prediction on test_strat_IUC (SVM)----------------------------------------------------
      test_strat_IUC$pred_svm5_IUC <- predict(model_svm5_IUC, test_strat_IUC)
      
      # Model evaluation
      confm_svm5_IUC <- table(actual = test_strat_IUC$IUC, prediction = test_strat_IUC$pred_svm5_IUC)
      
      metrics_svm3_IUC<- compute_metrics(confm_svm5_IUC)
      
      cat("Performance of SVM-Polynomial model based on the Wald test:\n")
      metrics_svm3_IUC
    }
    
    else if (model=="CatBoost-Chi Squared test"){
      # CatBoost Models
      Bio_Data <- as.data.frame(Bio_Data)
      
      #Convert categorical variables to factor
      categorical_var <- c(colnames(Bio_Data))
      
      Bio_Data[, categorical_var[-1]] <- lapply(Bio_Data[, categorical_var[-1]], factor)
      
      levels(Bio_Data$Ampicillin) <- c(levels(Bio_Data$Ampicillin),0)
      
      #Divide Data set into train and test
      set.seed(6546)
      
      split_strat_IUC  <- initial_split(Bio_Data, prop = 0.8, strata = "IUC")
      
      train_strat_IUC  <- training(split_strat_IUC)
      
      test_strat_IUC   <- testing(split_strat_IUC)
      
      # Training CatBoost model based on the selected antibiotics through chi squared test
      grid_IUC <- expand.grid(
        depth = c(2,3,4,5),
        learning_rate = c(0.1,0.2,0.3),
        l2_leaf_reg = c(1,2,3),
        rsm = 1,
        border_count = 1,
        iterations = seq(10,150,10)
      )
      
      fitControl_IUC <- trainControl(method = "cv",
                                 number = 10,
                                 classProbs = TRUE)
      
      set.seed(1234)
      model_catboost_IUC <- train(x = train_strat_IUC[,c("Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                         "Ceftazidime", "Ciprofloxacin", "Meropenem", "Ceftriaxone", "Aztreonam")],
                                  y = make.names(train_strat_IUC[,"IUC"]),
                                  maximize = TRUE,
                                  method = catboost.caret, metric = "Accuracy", 
                                  tuneGrid =  grid_IUC,
                                  trControl = fitControl_IUC)
      
      cat("Tuned hyper-parameters through 10-Fold Cross Validation:\n")
      print(model_catboost_IUC$bestTune)  
      
      #Prediction on test_strat_IUC (CatBoost)----------------------------------------------------
      test_strat_IUC$cat_IUC = predict(model_catboost_IUC, test_strat_IUC)
      
      # Model evaluation
      cat_test_IUC <- table(actual = test_strat_IUC$IUC, prediction = test_strat_IUC$cat_IUC)
      
      metrics_cat_IUC<- compute_metrics(cat_test_IUC)
      
      cat("Performance of CatBoost model based on Chi-Squared test results:\n")
      metrics_cat_IUC
    }
    
    else if (model=="CatBoost-model agnostic"){
      # CatBoost Models
      Bio_Data <- as.data.frame(Bio_Data)
      
      #Convert categorical variables to factor
      categorical_var <- c(colnames(Bio_Data))
      
      Bio_Data[, categorical_var[-1]] <- lapply(Bio_Data[, categorical_var[-1]], factor)
      
      levels(Bio_Data$Ampicillin) <- c(levels(Bio_Data$Ampicillin),0)
      
      #Divide Data set into train and test
      set.seed(6546)
      
      split_strat_IUC  <- initial_split(Bio_Data, prop = 0.8, strata = "IUC")
      
      train_strat_IUC  <- training(split_strat_IUC)
      
      test_strat_IUC   <- testing(split_strat_IUC)
      
      grid_IUC <- expand.grid(
        depth = c(2,3,4),
        learning_rate = c(0.2,0.3,0.4),
        l2_leaf_reg = c(2,3,4),
        rsm = 1,
        border_count = 1,
        iterations = seq(10,35,1)
      )
      
      fitControl_IUC <- trainControl(method = "cv",
                                     number = 10,
                                     classProbs = TRUE)
      
      model_catboost1_IUC <- train(x = train_strat_IUC[,c("Imipenem",  "Meropenem", "Aztreonam")],
                                   y = make.names(train_strat_IUC[,"IUC"]),
                                   maximize = TRUE,
                                   method = catboost.caret, metric = "Accuracy", 
                                   tuneGrid =  grid_IUC,
                                   trControl = fitControl_IUC)
      
      cat("Tuned hyper-parameters through 10-Fold Cross Validation:\n")
      print(model_catboost1_IUC$bestTune)  
      
      test_strat_IUC$cat_IUC1 = predict(model_catboost1_IUC, test_strat_IUC)
      
      # Model Evaluation
      confm_cat_IUC1 <- table(actual = test_strat_IUC$IUC, prediction = test_strat_IUC$cat_IUC1)
      
      metrics_cat1_IUC<- compute_metrics(confm_cat_IUC1)
      
      cat("Performance of CatBoost model based on model-agnostic approach results:\n")
      metrics_cat1_IUC
    }
    
    else if (model=="CatBoost-Wald test"){
      # CatBoost Models
      Bio_Data <- as.data.frame(Bio_Data)
      
      #Convert categorical variables to factor
      categorical_var <- c(colnames(Bio_Data))
      
      Bio_Data[, categorical_var[-1]] <- lapply(Bio_Data[, categorical_var[-1]], factor)
      
      levels(Bio_Data$Ampicillin) <- c(levels(Bio_Data$Ampicillin),0)
      
      #Divide Data set into train and test
      set.seed(6546)
      
      split_strat_IUC  <- initial_split(Bio_Data, prop = 0.8, strata = "IUC")
      
      train_strat_IUC  <- training(split_strat_IUC)
      
      test_strat_IUC   <- testing(split_strat_IUC)
      
      grid2_IUC <- expand.grid(
        depth = c(2,3,4,5,6),
        learning_rate = c(0.1,0.2,0.3,0.4,0.5),
        l2_leaf_reg = c(1,2,3,4,5),
        rsm = 1,
        border_count = 1,
        iterations = seq(10,30,1)
      )
      
      fitControl2_IUC <- trainControl(method = "cv",
                                  number = 10,
                                  classProbs = TRUE)
      
      model_catboost2_IUC <- train(x = train_strat_IUC[,c("Aztreonam","Imipenem","Meropenem","Ceftazidime","Cefepime")],
                                   y = make.names(train_strat_IUC[["IUC"]]),
                                   maximize = TRUE,
                                   method = catboost.caret, metric = "Accuracy", 
                                   tuneGrid =  grid2_IUC,
                                   trControl = fitControl2_IUC)
      
      
      cat("Tuned hyper-parameters through 10-Fold Cross Validation:\n")
      print(model_catboost2_IUC$bestTune)  
      
      test_strat_IUC$cat_IUC2 = predict(model_catboost2_IUC, test_strat_IUC)
      # Model evaluation
      confm_cat2_IUC <- table(actual = test_strat_IUC$IUC, prediction = test_strat_IUC$cat_IUC2)
      
      metrics_cat3_IUC<- compute_metrics(confm_cat2_IUC)
      
      cat("Performance of CatBoost model based on the Wald test results:\n")
      metrics_cat3_IUC
    }
    
    else if (model=="CatBoost-Simplicity Principle_MEM"){
      # CatBoost Models
      Bio_Data <- as.data.frame(Bio_Data)
      
      #Convert categorical variables to factor
      categorical_var <- c(colnames(Bio_Data))
      
      Bio_Data[, categorical_var[-1]] <- lapply(Bio_Data[, categorical_var[-1]], factor)
      
      levels(Bio_Data$Ampicillin) <- c(levels(Bio_Data$Ampicillin),0)
      
      #Divide Data set into train and test
      set.seed(6546)
      
      split_strat_IUC  <- initial_split(Bio_Data, prop = 0.8, strata = "IUC")
      
      train_strat_IUC  <- training(split_strat_IUC)
      
      test_strat_IUC   <- testing(split_strat_IUC)
      
      grid3_IUC <- expand.grid(
        depth = c(2,3,4),
        learning_rate = c(0.1,0.2),
        l2_leaf_reg = c(1,2,3),
        rsm = 1,
        border_count = 1,
        iterations = seq(10,30,1)
      )
      
      fitControl3_IUC <- trainControl(method = "cv",
                                      number = 10,
                                      classProbs = TRUE)
      
      model_catboost3_IUC <- train(x = train_strat_IUC["Meropenem"],
                                   y = make.names(train_strat_IUC[["IUC"]]),
                                   maximize = TRUE,
                                   method = catboost.caret, metric = "Accuracy", 
                                   tuneGrid =  grid3_IUC,
                                   trControl = fitControl3_IUC)
      
      
      cat("Tuned hyper-parameters through 10-Fold Cross Validation:\n")
      print(model_catboost3_IUC$bestTune)  
      
      test_strat_IUC$cat_IUC3 = predict(model_catboost3_IUC, test_strat_IUC)
      # Model evaluation
      confm_cat3_IUC <- table(actual = test_strat_IUC$IUC, prediction = test_strat_IUC$cat_IUC3)
      
      metrics_cat4_IUC<- compute_metrics(confm_cat3_IUC)
      
      cat("Performance of CatBoost model based on the Simplicity Principle (MEM):\n")
      metrics_cat4_IUC
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
        result <- chisq.test(table(data[[variable]], data$IUC))
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
    model_Log1_IUC <- glm(IUC ~ Imipenem + Gentamicin + Ciprofloxacin + Meropenem + Aztreonam +
                            Cefotaxime + Cefepime + Ceftazidime + Ceftriaxone
                          , family = "binomial", data = train_strat_IUC)
    
    cat("P-values show Wald test results:\n")
    print(summary(model_Log1_IUC))
  }
  else if (method == "model agnostic-LR"){
    # Feature Selection by DALEX Package--------------------------------------
    # the 'explained_glm_TEM' object is created using the 'explain' function from the DALEX package. 
    # The variable_importance function is used to calculate the permutation-based 
    # feature importance ('fi_glm_TEM') with 50 permutations. Finally, the feature importance is displayed 
    # using 'fi_glm_TEM' and plotted using' plot(fi_glm_TEM)'.
    model_Log1_IUC <- glm(IUC ~ Imipenem + Gentamicin + Ciprofloxacin + Meropenem + Aztreonam +
                            Cefotaxime + Cefepime + Ceftazidime + Ceftriaxone
                          , family = "binomial", data = train_strat_IUC)
    
    explained_glm_IUC <- explain(model = model_Log1_IUC, data=train_strat_IUC[,c(3:9,11:12)], variables = colnames(train_strat_IUC[,c(3:9,11:12)]),
                                 y=as.vector(as.numeric(train_strat_IUC$IUC))-1, label = "LR", 
                                 type = "classification")
    #50 permuatation
    fi_glm_IUC = variable_importance(explained_glm_IUC, B=50 ,variables = colnames(train_strat_IUC[,c(3:9,11:12)]), loss_function = loss_root_mean_square, type = "raw")
    
    cat("Root mean squared dropout loss due to a feature elimination through model-agnostic approach in LR:\n")
    print(fi_glm_IUC)
  }
  else if (method == "model agnostic-NBC"){
    train_control <- trainControl(
      method = "cv",
      number = 10
    )
    
    search_grid_IUC <- expand.grid(
      usekernel = c(FALSE,TRUE),
      fL = seq(0,5,0.5),
      adjust = seq(0,5,0.5)
    )
    
    model_nbdalex_IUC <- train(IUC ~ Imipenem + Gentamicin + Ciprofloxacin + Meropenem + Aztreonam +
                                 Cefotaxime + Cefepime + Ceftazidime + Ceftriaxone,
                               data = train_strat_IUC,
                               method = "nb",
                               metric = "Accuracy",
                               trControl = train_control,
                               tuneGrid = search_grid_IUC
    )
    # Feature selection through model-agnostic approach (DALEX)
    # the 'explain' function from the DALEX package is used to explain the Naive Bayes classifier model 
    # 'model_nbdalex_TEM'. The relevant variables are selected and provided as data along with their 
    # corresponding column names. The target variable 'IUC_gene' is transformed to numeric values 
    # (y = as.vector(as.numeric(train_strat_IUC$IUC_gene)) - 1). The label is set to "NBC" for Naive Bayes classifier,
    # and the type is set to "classification".
    expnb_IUC = explain(model = model_nbdalex_IUC, data=train_strat_IUC[,c("Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                           "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")], variables = colnames(train_strat_IUC[,c("Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                                                                                                                                             "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")]),
                        y=as.vector(as.numeric(train_strat_IUC$IUC))-1, label = "NBC", 
                        type = "classification")
    
    # The 'variable_importance' function is then used to calculate the variable importance based on the
    # explained model. The 'B' parameter is set to 50 for the number of permutations, and the loss function
    # is set to "loss_root_mean_square". The resulting variable importance is stored in 'fitnb_TEM'
    fitnb_IUC = variable_importance(expnb_IUC, B=50 ,variables = colnames(train_strat_IUC[,c("Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                                             "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")]), loss_function = loss_root_mean_square, type = "raw" )
    
    cat("Root mean squared dropout loss due to a feature elimination through model-agnostic approach in NBC:\n")
    print(fitnb_IUC)
  }
  
  else if (method == "model agnostic-LDA"){
    model_ldadalex_IUC <- train(IUC ~ Cefotaxime + Gentamicin + Imipenem +
                                  Cefepime + Ceftazidime + Ciprofloxacin +
                                  Meropenem + Ceftriaxone + Aztreonam,
                                data = train_strat_IUC,
                                method = "lda",
                                metric = "Accuracy")
    
    explda_IUC = explain(model = model_ldadalex_IUC, data=train_strat_IUC[,c("Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                             "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")], variables = colnames(train_strat_IUC[,c("Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                                                                                                                                               "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")]),
                         y=as.vector(as.numeric(train_strat_IUC$IUC))-1, label = "LDA", 
                         type = "classification")
    
    fitlda_IUC = variable_importance(explda_IUC, B=50 ,variables = colnames(train_strat_IUC[,c("Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                                               "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")]), loss_function = loss_root_mean_square, type = "raw" )
    
    cat("Root mean squared dropout loss due to a feature elimination through model-agnostic approach in LDA:\n")
    print(fitlda_IUC)
  }
  else if (method == "model agnostic-Sig SVM"){
    set.seed(1234)
    tune_out2_IUC <- e1071::tune("svm", IUC ~ Imipenem + Gentamicin + Ciprofloxacin + Meropenem + 
                                   Aztreonam + Cefotaxime + Cefepime + Ceftazidime + Ceftriaxone  
                                 , data = train_strat_IUC, kernel = "sigmoid", tunecontrol=tune.control(cross=10),
                                 ranges = list(cost = seq(0.1, 10, length = 30)
                                               , gamma = seq(0.1,1, length = 10)))
    
    
    model_svm2_IUC = svm(IUC ~ Imipenem + Gentamicin + Ciprofloxacin + Meropenem + Aztreonam +
                           Cefotaxime + Cefepime + Ceftazidime + Ceftriaxone, data = train_strat_IUC, 
                          probability = T, kernel = "sigmoid", cost = tune_out2_IUC$best.parameters$cost,
                          gamma = tune_out2_IUC$best.parameters$gamma )
    
    explained_SVM2_IUC <- explain(model = model_svm2_IUC, data=train_strat_IUC[,c(3:9,11:12)], variables = colnames(train_strat_IUC[,c(3:9,11:12)]),
                                  y=as.vector(as.numeric(train_strat_IUC$IUC))-1, label = "SVM-Sigmoid", 
                                  type = "classification")
    
    #50 permuatation
    fi_SVM2_IUC = variable_importance(explained_SVM2_IUC, B=50 ,variables = colnames(train_strat_IUC[,c(3:9,11:12)]),loss_function = loss_root_mean_square, type = "raw")
    
    cat("Root mean squared dropout loss due to a feature elimination through model-agnostic approach in sigmoid kernel of SVM:\n")
    print(fi_SVM2_IUC)
  }
  
  else if (method == "model agnostic-Poly SVM"){
    
    set.seed(1234)
    tune_out4_IUC <- e1071::tune("svm", IUC ~ Imipenem + Gentamicin + Ciprofloxacin + Meropenem + 
                                   Aztreonam + Cefotaxime + Cefepime + Ceftazidime + Ceftriaxone  
                                 , data = train_strat_IUC, kernel = "polynomial", tunecontrol=tune.control(cross=10),
                                 ranges = list(cost = seq(0.1, 5, length = 20), degree = seq(1,4,1)
                                               , gamma = seq(0.1,5, length = 20)))
    
    
    # Tuned hyper-parameters are placed in the model.
    model_svm4_IUC = svm(IUC ~ Imipenem + Gentamicin + Ciprofloxacin + Meropenem + 
                           Aztreonam + Cefotaxime + Cefepime + Ceftazidime + Ceftriaxone
                           , data = train_strat_IUC, probability = T, 
                           kernel = "polynomial", cost = tune_out4_IUC$best.parameters$cost, 
                           degree = tune_out4_IUC$best.parameters$degree , gamma = tune_out4_IUC$best.parameters$gamma )
    
    explained_SVM4_IUC <- explain(model = model_svm4_IUC, data=train_strat_IUC[,c(3:9,11:12)], variables = colnames(train_strat_IUC[,c(3:9,11:12)]),
                                  y=as.vector(as.numeric(train_strat_IUC$IUC))-1, label = "SVM-Polynomial", 
                                  type = "classification")
    
    #50 permuatation
    fi_SVM4_IUC = variable_importance(explained_SVM4_IUC, B=50 ,variables = colnames(train_strat_IUC[,c(3:9,11:12)]),loss_function = loss_root_mean_square, type = "raw")
    
    cat("Root mean squared dropout loss due to a feature elimination through model-agnostic approach in polynomial kernel of SVM:\n")
    print(fi_SVM4_IUC)
  }
  else if (method == "model agnostic-CatBoost"){
    
    # CatBoost Models
    Bio_Data <- as.data.frame(Bio_Data)
    
    #Convert categorical variables to factor
    categorical_var <- c(colnames(Bio_Data))
    
    Bio_Data[, categorical_var[-1]] <- lapply(Bio_Data[, categorical_var[-1]], factor)
    
    levels(Bio_Data$Ampicillin) <- c(levels(Bio_Data$Ampicillin),0)
    
    #Divide Data set into train and test
    set.seed(6546)
    
    split_strat_IUC  <- initial_split(Bio_Data, prop = 0.8, strata = "IUC")
    
    train_strat_IUC  <- training(split_strat_IUC)
    
    test_strat_IUC   <- testing(split_strat_IUC)
    
    
    grid_IUC <- expand.grid(
      depth = c(3,4,5),
      learning_rate = c(0.1,0.2,0.3),
      l2_leaf_reg = c(1,2),
      rsm = 1,
      border_count = 1,
      iterations = seq(10,61,10)
    )
    
    fitControl_IUC <- trainControl(method = "cv",
                               number = 10,
                               classProbs = TRUE)
    
    set.seed(1234)
    model_catboost_IUC <- train(x = train_strat_IUC[,c("Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                       "Ceftazidime", "Ciprofloxacin", "Meropenem", "Ceftriaxone", "Aztreonam")],
                                y = make.names(train_strat_IUC[,"IUC"]),
                                maximize = TRUE,
                                method = catboost.caret, metric = "Accuracy", 
                                tuneGrid =  grid_IUC,
                                trControl = fitControl_IUC)
    
    expsda_IUC = explain(model = model_catboost_IUC, data=train_strat_IUC[,c("Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                             "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")], variables = colnames(train_strat_IUC[,c("Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                                                                                                                                               "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")]),
                         y=as.vector(as.numeric(train_strat_IUC$IUC))-1, label = "CatBoost", 
                         type = "classification")
    
    fitcat_IUC = variable_importance(expsda_IUC, B=50 ,variables = colnames(train_strat_IUC[,c("Cefotaxime", "Gentamicin","Imipenem", "Cefepime",
                                                                                               "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")]), loss_function = loss_root_mean_square, type = "raw" )
    
    cat("Root mean squared dropout loss due to a feature elimination through model-agnostic approach in CatBoost:\n")
    print(fitcat_IUC)
  }
  else {
    cat("Invalid task. Please specify one of the following tasks: feature_selectio, model_evaluation\n")
    
  } 
}
if (is.null(task)){
  cat("You should specify a task\n")
}

