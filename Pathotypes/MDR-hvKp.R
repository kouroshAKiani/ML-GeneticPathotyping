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

proportion_HVKPMDR <- table(Bio_Data$`HVKP&MDR`) %>% prop.table()

cat("Proportion of IUC gene:\n")

print(proportion_HVKPMDR)
# Response variable is imbalanced (0 : 78% , 1 : 22%) so we're using stratified sampling

# Divide the data set into train and test sets through stratified sampling method

set.seed(123)

split_strat_HVKPMDR  <- initial_split(Bio_Data, prop = 0.8, strata = "HVKP&MDR")

train_strat_HVKPMDR  <- training(split_strat_HVKPMDR)

test_strat_HVKPMDR   <- testing(split_strat_HVKPMDR)

# Display the dimensions of the train and test sets
cat("train set dimensions:", dim(train_strat_HVKPMDR), "\n")
cat("test set dimensions:", dim(test_strat_HVKPMDR), "\n")


if (!is.null(task)) {
  if (task == "model_evaluation") {
    
    if (model == "LR-Chi Squared test/Wald test"){
      
      model_Log1_HVKPMDR <- glm(`HVKP&MDR` ~ Imipenem + Gentamicin + Ciprofloxacin + Meropenem + Aztreonam +
                                  Cefotaxime + Cefepime + Ceftazidime + Ceftriaxone + Amikacin
                                , family = "binomial", data = train_strat_HVKPMDR)
      
      cat("Summary of LR model based on Chi Squared test and Wald test including Null deviance, Residual deviance, and AIC, etc.:\n")
      print(summary(model_Log1_HVKPMDR))
      
      test_strat_HVKPMDR$probs1_HVKPMDR <- predict(model_Log1_HVKPMDR, test_strat_HVKPMDR, type = "response")
      
      test_strat_HVKPMDR$pred_logreg1_HVKPMDR <- ifelse(test_strat_HVKPMDR$probs1_HVKPMDR >= 0.5, 1, 0)
      
      # Model 1 evaluation
      confm_logreg1_HVKPMDR <- table(actual = test_strat_HVKPMDR$`HVKP&MDR`, prediction = test_strat_HVKPMDR$pred_logreg1_HVKPMDR)
      
      metrics_logreg1_HVKPMDR <- compute_metrics(confm_logreg1_HVKPMDR)
      
      cat("Performance of LR model based on Chi-Squared test and Wald test results:\n")
      metrics_logreg1_HVKPMDR
    }
    else if (model == "LR-model agnostic"){
      model_Log2_HVKPMDR <- glm( `HVKP&MDR` ~ Meropenem
                                 , family = "binomial", data = train_strat_HVKPMDR)
      
      cat("Summary of LR model based on the model-agnostic approach including Null deviance, Residual deviance, and AIC, etc.:\n")
      print(summary(model_Log2_HVKPMDR))
      
      ##Prediction on test_strat_IUC (Model 2)--------------------------------------------------
      test_strat_HVKPMDR$probs2_HVKPMDR <- predict(model_Log2_HVKPMDR, test_strat_HVKPMDR, type = "response")
      
      test_strat_HVKPMDR$pred_logreg2_HVKPMDR <- ifelse(test_strat_HVKPMDR$probs2_HVKPMDR >= 0.5, 1, 0)
      
      confm_logreg2_HVKPMDR <- table(actual = test_strat_HVKPMDR$`HVKP&MDR`, prediction = test_strat_HVKPMDR$pred_logreg2_HVKPMDR)
      
      metrics_logreg2_HVKPMDR <- compute_metrics(confm_logreg2_HVKPMDR)
      
      cat("Performance of LR model based on the model-agnostic approach:\n")
      metrics_logreg2_HVKPMDR
    }
    
    else if (model == "NBC-Chi Squared test/Wald test"){
      train_control <- trainControl(
        method = "cv",
        number = 10
      )
      
      search_grid_HVKPMDR <- expand.grid(
        usekernel = c(FALSE,TRUE),
        fL = seq(0,5,0.5),
        adjust = seq(0,5,0.5)
      )
      
      model_nbdalex_HVKPMDR <- train(`HVKP&MDR` ~ Imipenem + Gentamicin + Ciprofloxacin + Meropenem + Aztreonam +
                                       Cefotaxime + Cefepime + Ceftazidime + Ceftriaxone + Amikacin,
                                     data = train_strat_HVKPMDR,
                                     method = "nb",
                                     metric = "Accuracy",
                                     trControl = train_control,
                                     tuneGrid = search_grid_HVKPMDR,
      )
      
      
      cat("Tuned hyper-parameters through 10-fold Cross-Validation:\n")
      print(model_nbdalex_HVKPMDR$bestTune)
      
      ##Prediction on test_strat_IUC (Naive Bayes Model)----------------------------------
      test_strat_HVKPMDR$pred_nb_HVKPMDR <- predict(model_nbdalex_HVKPMDR, test_strat_HVKPMDR)
      
      # Model evaluation
      confm_nb_HVKPMDR <- table(actual = test_strat_HVKPMDR$`HVKP&MDR`, prediction = test_strat_HVKPMDR$pred_nb_HVKPMDR)
      metrics_nb_HVKPMDR<- compute_metrics(confm_nb_HVKPMDR)
      cat("Performance of NBC model based on Chi-Squared test and Wald test results:\n")
      metrics_nb_HVKPMDR
    }
    else if (model == "NBC-model agnostic"){
      train_control <- trainControl(
        method = "cv",
        number = 10
      )
      
      search_grid_HVKPMDR <- expand.grid(
        usekernel = c(FALSE,TRUE),
        fL = seq(0,5,0.5),
        adjust = seq(0,5,0.5)
      )
      
      model_nb_HVKPMDR <- train(`HVKP&MDR` ~ Meropenem,
                                data = train_strat_HVKPMDR,
                                method = "nb",
                                metric = "Accuracy",
                                trControl = train_control,
                                tuneGrid = search_grid_HVKPMDR,
      )
      
      cat("Tuned hyper-parameters through 10-fold Cross-Validation:\n")
      print(model_nb_HVKPMDR$bestTune)
      
      ##Prediction on test_strat_IUC (Naive Bayes Model)----------------------------------
      test_strat_HVKPMDR$pred_nb1_HVKPMDR <- predict(model_nb_HVKPMDR, test_strat_HVKPMDR)
      
      # Model evaluation
      confm_nb1_HVKPMDR <- table(actual = test_strat_HVKPMDR$`HVKP&MDR`, prediction = test_strat_HVKPMDR$pred_nb1_HVKPMDR)
      
      metrics_nb1_HVKPMDR<- compute_metrics(confm_nb1_HVKPMDR)
      cat("Performance of NBC model based on model-agnostic approach results:\n")
      metrics_nb1_HVKPMDR
    }
    
    else if (model == "LDA-Chi Squared test/Wald test"){
      
      model_ldadalex_HVKPMDR <- train(`HVKP&MDR` ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                                        Cefepime + Ceftazidime + Ciprofloxacin +
                                        Meropenem + Ceftriaxone + Aztreonam,
                                      data = train_strat_HVKPMDR,
                                      method = "lda",
                                      metric = "Accuracy")
      
      cat("LDA model based on Chi-Squared test and Wald test results:\n")
      print(model_ldadalex_HVKPMDR)
      
      #Prediction on test_strat_IUC (LDA Model)-----------------------------------------------
      test_strat_HVKPMDR$pred_lda_HVKPMDR <- predict(model_ldadalex_HVKPMDR, test_strat_HVKPMDR)
      
      # Model Evaluation-------------------------------------------------
      confm_lda_HVKPMDR <- table(actual = test_strat_HVKPMDR$`HVKP&MDR`, prediction = test_strat_HVKPMDR$pred_lda_HVKPMDR)
      
      metrics_lda_HVKPMDR<- compute_metrics(confm_lda_HVKPMDR)
      
      cat("Performance of LDA model based on Chi-Squared test and Wald test results:\n")
      metrics_lda_HVKPMDR
    }
    else if (model == "LDA-model agnostic"){
      model_lda_HVKPMDR <- train(`HVKP&MDR` ~ Meropenem,
                                 data = train_strat_HVKPMDR,
                                 method = "lda",
                                 metric = "Accuracy")
      
      cat("LDA model based on model-agnostic approach:\n")
      print(model_lda_HVKPMDR)
      
      #Prediction on test_strat_IUC (LDA Model)-----------------------------------------------
      test_strat_HVKPMDR$pred_lda1_HVKPMDR <- predict(model_lda_HVKPMDR, test_strat_HVKPMDR)
      
      # Model evaluation
      confm_lda1_HVKPMDR <- table(actual = test_strat_HVKPMDR$`HVKP&MDR`, prediction = test_strat_HVKPMDR$pred_lda1_HVKPMDR)
      
      metrics_lda1_HVKPMDR <- compute_metrics(confm_lda1_HVKPMDR)
      
      cat("Performance of LDA model based on model-agnostic approach:\n")
      metrics_lda1_HVKPMDR
    }
    
    else if (model == "SVM-Sig-Chi Squared test/Wald test"){
      # Fit Support Vector Machine model to train data set
      # This code performs a grid search for tuning the hyper-parameters of the SVM model using the tune function
      # from the e1071 package. It then fits the SVM model with the specified hyper-parameters using the svm function.
      set.seed(1234)
      tune_out_HVKPMDR <- e1071::tune("svm", `HVKP&MDR` ~ Imipenem + Gentamicin + Ciprofloxacin + 
                                        Meropenem + Aztreonam + Cefotaxime + Cefepime + Ceftazidime + 
                                        Ceftriaxone + Amikacin
                                      , data = train_strat_HVKPMDR, kernel = "sigmoid", tunecontrol=tune.control(cross=10),
                                      ranges = list( cost = seq(0.1, 5, length = 20)
                                                     , gamma = seq(0.1, 5, length = 20)))
      
      cat("Tuned hyper-parameters through 10-Fold Cross Validation:\n")
      print(summary(tune_out_HVKPMDR))
      
      model_svm_HVKPMDR = svm(`HVKP&MDR` ~ Imipenem + Gentamicin + Ciprofloxacin + 
                                Meropenem + Aztreonam + Cefotaxime + Cefepime + Ceftazidime + 
                                Ceftriaxone + Amikacin  , data = train_strat_HVKPMDR, 
                           probability = T, kernel = "sigmoid", cost = tune_out_HVKPMDR$best.parameters$cost,
                           gamma = tune_out_HVKPMDR$best.parameters$gamma )
      
      cat("Sigmoid SVM based on Chi-Squared test and Wald test results:\n")
      print(model_svm_HVKPMDR)
      
      #Prediction on test_strat_IUC (SVM)----------------------------------------------------
      test_strat_HVKPMDR$pred_svm_HVKPMDR <- predict(model_svm_HVKPMDR, test_strat_HVKPMDR)
      
      # Model evaluation
      confm_svm_HVKPMDR <- table(actual = test_strat_HVKPMDR$`HVKP&MDR`, prediction = test_strat_HVKPMDR$pred_svm_HVKPMDR)
      
      metrics_svm1_HVKPMDR_sig<- compute_metrics(confm_svm_HVKPMDR)
      
      cat("Performance of SVM-Sigmoid model based on Chi-Squared test and Wald test results:\n")
      metrics_svm1_HVKPMDR_sig
    }
    else if (model=="SVM-Sig-model agnostic"){
      
      set.seed(1234)
      tune_out2_HVKPMDR <- e1071::tune("svm", `HVKP&MDR` ~ Meropenem,
                                       data = train_strat_HVKPMDR,
                                       kernel = "sigmoid", tunecontrol=tune.control(cross=10),
                                       ranges = list(cost = seq(0.1, 5, length = 25)
                                                     , gamma = seq(0.1, 5, length = 25)))
      
      cat("Tuned hyper-parameters through 10-Fold Cross Validation:\n")
      print(summary(tune_out2_HVKPMDR))
      
      model_svm1_HVKPMDR = svm(`HVKP&MDR` ~ Meropenem, data = train_strat_HVKPMDR, 
                            probability = T, kernel = "sigmoid", cost = tune_out2_HVKPMDR$best.parameters$cost,
                            gamma = tune_out2_HVKPMDR$best.parameters$gamma)
      
      cat("Sigmoid SVM based on model-agnostic approach:\n")
      print(model_svm1_HVKPMDR)
      
      test_strat_HVKPMDR$pred_svm1_HVKPMDR <- predict(model_svm1_HVKPMDR, test_strat_HVKPMDR)
      
      confm_svm1_HVKPMDR <- table(actual = test_strat_HVKPMDR$`HVKP&MDR`, prediction = test_strat_HVKPMDR$pred_svm1_HVKPMDR)
      
      metrics_svm2_HVKPMDR_sig<- compute_metrics(confm_svm1_HVKPMDR)
      
      cat("Performance of SVM-Sigmoid model based on model-agnostic approach:\n")
      metrics_svm2_HVKPMDR_sig
    }
   
    else if (model=="SVM-Poly-Chi Squared test/Wald test"){
      set.seed(1234)
      tune_out3_HVKPMDR <- e1071::tune("svm", `HVKP&MDR` ~ Imipenem + Gentamicin + Ciprofloxacin + 
                                         Meropenem + Aztreonam + Cefotaxime + Cefepime + Ceftazidime + 
                                         Ceftriaxone + Amikacin
                                       , data = train_strat_HVKPMDR, kernel = "polynomial", tunecontrol=tune.control(cross=10),
                                       ranges = list( cost = seq(0.1, 5, length = 20), degree = c(1,2,3)
                                                      , gamma = seq(0.1, 5, length = 20)))
      
      cat("Tuned hyper-parameters through 10-Fold Cross Validation:\n")
      print(summary(tune_out3_HVKPMDR))
      
      # Tuned hyper-parameters are placed in the model.
      model_svm2_HVKPMDR = svm(`HVKP&MDR` ~ Imipenem + Gentamicin + Ciprofloxacin + 
                                 Meropenem + Aztreonam + Cefotaxime + Cefepime + Ceftazidime + 
                                 Ceftriaxone + Amikacin
                           , data = train_strat_HVKPMDR, probability = T, 
                           kernel = "polynomial", cost = tune_out3_HVKPMDR$best.parameters$cost, 
                           degree = tune_out3_HVKPMDR$best.parameters$degree , gamma = tune_out3_HVKPMDR$best.parameters$gamma )
      
      cat("Polynomial SVM based on Chi-Squared test and Wald test results:\n")
      print(model_svm2_HVKPMDR)
      
      #Prediction on test (SVM 2)----------------------------------------------------
      test_strat_HVKPMDR$pred_svm2_HVKPMDR <- predict(model_svm2_HVKPMDR, test_strat_HVKPMDR)
      
      # Model evaluation
      confm_svm2_HVKPMDR <- table(actual = test_strat_HVKPMDR$`HVKP&MDR`, prediction = test_strat_HVKPMDR$pred_svm2_HVKPMDR)
      
      metrics_svm_HVKPMDR<- compute_metrics(confm_svm2_HVKPMDR)
      
      cat("Performance of SVM-Polynomial model based on Chi-Squared test and Wald test results:\n")
      metrics_svm_HVKPMDR
    }
    else if (model=="SVM-Poly-model agnostic"){
      
      set.seed(1234)
      tune_out4_HVKPMDR <- e1071::tune("svm", `HVKP&MDR` ~ Meropenem,
                                       data = train_strat_HVKPMDR,
                                       kernel = "polynomial", tunecontrol=tune.control(cross=10),
                                       ranges = list(cost = seq(0.1, 5, length = 20), degree = c(1,2,3)
                                                     , gamma = seq(0.1, 5, length = 20)))
      
      cat("Tuned hyper-parameters through 10-Fold Cross Validation:\n")
      print(summary(tune_out4_HVKPMDR))
      
      model_svm3_HVKPMDR = svm(`HVKP&MDR` ~ Meropenem, data = train_strat_HVKPMDR, 
                            probability = T, kernel = "polynomial", cost = tune_out4_HVKPMDR$best.parameters$cost,
                            degree = tune_out4_HVKPMDR$best.parameters$degree , gamma = tune_out4_HVKPMDR$best.parameters$gamma )
      
      cat("Polynomial SVM based on model-agnostic approach:\n")
      print(model_svm3_HVKPMDR)
      
      #Prediction on test_strat_IUC (SVM)----------------------------------------------------
      test_strat_HVKPMDR$pred_svm3_HVKPMDR <- predict(model_svm3_HVKPMDR, test_strat_HVKPMDR)
      
      # Model evaluation
      confm_svm3_HVKPMDR <- table(actual = test_strat_HVKPMDR$`HVKP&MDR`, prediction = test_strat_HVKPMDR$pred_svm3_HVKPMDR)
      
      metrics_svm3_HVKPMDR<- compute_metrics(confm_svm3_HVKPMDR)
      
      cat("Performance of SVM-Polynomial model based on model-agnostic approach:\n")
      metrics_svm3_HVKPMDR
      
    }
    
    else if (model=="CatBoost-Chi Squared test/Wald test"){
      # CatBoost Models--------------------------------------------------------------------
      Bio_Data <- as.data.frame(Bio_Data)
      
      #Convert categorical variables to factor----------------------------------
      categorical_var <- c(colnames(Bio_Data))
      
      Bio_Data[, categorical_var[-1]] <- lapply(Bio_Data[, categorical_var[-1]], factor)
      
      levels(Bio_Data$Ampicillin) <- c(levels(Bio_Data$Ampicillin),0)
      
      #Divide Data set into Train and Test---------------------------------------
      set.seed(123)
     
       split_strat_HVKPMDR  <- initial_split(Bio_Data, prop = 0.8, strata = "HVKP&MDR")
      
       train_strat_HVKPMDR  <- training(split_strat_HVKPMDR)
      
       test_strat_HVKPMDR   <- testing(split_strat_HVKPMDR)
      
      grid_HVKPMDDR <- expand.grid(
        depth = c(2,3),
        learning_rate = c(0.1,0.2),
        l2_leaf_reg = c(1,2),
        rsm = 1,
        border_count = 1,
        iterations = seq(130,151,1)
      )
      fitControl_HVKPMDR <- trainControl(method = "cv",
                                         number = 10,
                                         classProbs = TRUE)
      set.seed(1234)
      model_catboost_HVKPMDR <- train(x = train_strat_HVKPMDR[,c("Cefotaxime", "Gentamicin","Imipenem", "Cefepime", "Amikacin",
                                                                 "Ceftazidime", "Ciprofloxacin", "Meropenem", "Ceftriaxone", "Aztreonam")],
                                      y = make.names(train_strat_HVKPMDR[,"HVKP&MDR"]),
                                      maximize = TRUE,
                                      method = catboost.caret, metric = "Accuracy", 
                                      tuneGrid =  grid_HVKPMDDR,
                                      trControl = fitControl_HVKPMDR)
      
      cat("Tuned hyper-parameters through 10-Fold Cross Validation:\n")
      print(model_catboost_HVKPMDR$bestTune)  
      
      test_strat_HVKPMDR$cat_HVKPMDR = predict(model_catboost_HVKPMDR, test_strat_HVKPMDR)
      
      cat_test_HVKPMDR <- table(actual = test_strat_HVKPMDR$`HVKP&MDR`, prediction = test_strat_HVKPMDR$cat_HVKPMDR)
      
      metrics_cat_HVKPMDR<- compute_metrics(cat_test_HVKPMDR)
      
      cat("Performance of CatBoost model based on Chi-Squared test and Wald test results:\n")
      metrics_cat_HVKPMDR
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
      
      split_strat_HVKPMDR  <- initial_split(Bio_Data, prop = 0.8, strata = "HVKP&MDR")
      
      train_strat_HVKPMDR  <- training(split_strat_HVKPMDR)
      
      test_strat_HVKPMDR   <- testing(split_strat_HVKPMDR)
      
      grid_HVKPMDR <- expand.grid(
        depth = c(2,3,4),
        learning_rate = c(0.1,0.2,0.3),
        l2_leaf_reg = c(1,2,3),
        rsm = 1,
        border_count = 1,
        iterations = seq(10,30,1)
      )
      
      fitControl_HVKPMDR <- trainControl(method = "cv",
                                     number = 10,
                                     classProbs = TRUE)
      
      model_catboost2_HVKPMDR <- train(x = train_strat_HVKPMDR["Meropenem"],
                                       y = make.names(train_strat_HVKPMDR[["HVKP"]]),
                                       maximize = TRUE,
                                       method = catboost.caret, metric = "Accuracy", 
                                       tuneGrid =  grid_HVKPMDR,
                                       trControl = fitControl_HVKPMDR)
      
      cat("Tuned hyper-parameters through 10-Fold Cross Validation:\n")
      print(model_catboost2_HVKPMDR$bestTune)  
      
      test_strat_HVKPMDR$cat_HVKPMDR = predict(model_catboost_HVKPMDR, test_strat_HVKPMDR)
      
      confm_cat_HVKPMDR <- table(actual = test_strat_HVKPMDR$`HVKP&MDR`, prediction = test_strat_HVKPMDR$cat_HVKPMDR)
      
      metrics_cat1_HVKPMDR<- compute_metrics(confm_cat_HVKPMDR)
      
      cat("Performance of CatBoost model based on model-agnostic approach results:\n")
      metrics_cat1_HVKPMDR
    }
  }
  
  else if (task == 'feature_selectio') {

    if (method == "Chi-Squared test"){
      # Function to perform chi-square test and print p-value
      perform_chi_square <- function(data, variable) {
        result <- chisq.test(table(data[[variable]], data$`HVKP&MDR`))
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
    model_Log1_HVKPMDR <- glm(`HVKP&MDR` ~ Imipenem + Gentamicin + Ciprofloxacin + Meropenem + Aztreonam +
                                Cefotaxime + Cefepime + Ceftazidime + Ceftriaxone + Amikacin
                              , family = "binomial", data = train_strat_HVKPMDR)
    
    cat("P-values show Wald test results:\n")
    print(summary(model_Log1_HVKPMDR))
  }
  else if (method == "model agnostic-LR"){
    
    model_Log1_HVKPMDR <- glm(`HVKP&MDR` ~ Imipenem + Gentamicin + Ciprofloxacin + Meropenem + Aztreonam +
                                Cefotaxime + Cefepime + Ceftazidime + Ceftriaxone + Amikacin
                              , family = "binomial", data = train_strat_HVKPMDR)
    
    explained_glm_HVKPMDR <- explain(model = model_Log1_HVKPMDR, data=train_strat_HVKPMDR[,c(2:9,11:12)], variables = colnames(train_strat_HVKPMDR[,c(2:9,11:12)]),
                                     y=as.vector(as.numeric(train_strat_HVKPMDR$`HVKP&MDR`))-1, label = "LR", 
                                     type = "classification")
    #50 permuatation
    fi_glm_HVKPMDR = variable_importance(explained_glm_HVKPMDR, B=50 ,variables = colnames(train_strat_HVKPMDR[,c(2:9,11:12)]), loss_function = loss_root_mean_square, type = "raw")
    
    cat("Root mean squared dropout loss due to a feature elimination through model-agnostic approach in LR:\n")
    print(fi_glm_HVKPMDR)
  }
  
  else if (method == "model agnostic-NBC"){
    train_control <- trainControl(
      method = "cv",
      number = 10
    )
    
    search_grid_HVKPMDR <- expand.grid(
      usekernel = c(FALSE,TRUE),
      fL = seq(0,5,0.5),
      adjust = seq(0,5,0.5)
    )
    
    model_nbdalex_HVKPMDR <- train(`HVKP&MDR` ~ Imipenem + Gentamicin + Ciprofloxacin + Meropenem + Aztreonam +
                                     Cefotaxime + Cefepime + Ceftazidime + Ceftriaxone + Amikacin,
                                   data = train_strat_HVKPMDR,
                                   method = "nb",
                                   metric = "Accuracy",
                                   trControl = train_control,
                                   tuneGrid = search_grid_HVKPMDR,
    )

    expnb_HVKPMDR = explain(model = model_nbdalex_HVKPMDR, data=train_strat_HVKPMDR[,c("Cefotaxime", "Gentamicin","Imipenem", "Cefepime", "Amikacin",
                                                                                       "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")], variables = colnames(train_strat_HVKPMDR[,c("Cefotaxime", "Gentamicin","Imipenem", "Cefepime", "Amikacin",
                                                                                                                                                                                                             "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")]),
                            y=as.vector(as.numeric(train_strat_HVKPMDR$`HVKP&MDR`))-1, label = "NBC", 
                            type = "classification")
    
    fitnb_HVKPMDR = variable_importance(expnb_HVKPMDR, B=50 ,variables = colnames(train_strat_HVKPMDR[,c("Cefotaxime", "Gentamicin","Imipenem", "Cefepime", "Amikacin",
                                                                                                         "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")]), loss_function = loss_root_mean_square, type = "raw" )
    
    cat("Root mean squared dropout loss due to a feature elimination through model-agnostic approach in NBC:\n")
    print(fitnb_HVKPMDR)
  }
  
  else if (method == "model agnostic-LDA"){
    model_ldadalex_HVKPMDR <- train(`HVKP&MDR` ~ Amikacin + Cefotaxime + Gentamicin + Imipenem +
                                      Cefepime + Ceftazidime + Ciprofloxacin +
                                      Meropenem + Ceftriaxone + Aztreonam,
                                    data = train_strat_HVKPMDR,
                                    method = "lda",
                                    metric = "Accuracy")
    
    explda_HVKPMDR = explain(model = model_ldadalex_HVKPMDR, data=train_strat_HVKPMDR[,c("Cefotaxime", "Gentamicin","Imipenem", "Cefepime", "Amikacin",
                                                                                         "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")], variables = colnames(train_strat_HVKPMDR[,c("Cefotaxime", "Gentamicin","Imipenem", "Cefepime", "Amikacin",
                                                                                                                                                                                                               "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")]),
                             y=as.vector(as.numeric(train_strat_HVKPMDR$`HVKP&MDR`))-1, label = "LDA", 
                             type = "classification")
    
    fitlda_HVKPMDR = variable_importance(explda_HVKPMDR, B=50 ,variables = colnames(train_strat_HVKPMDR[,c("Cefotaxime", "Gentamicin","Imipenem", "Cefepime", "Amikacin",
                                                                                                           "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")]), loss_function = loss_root_mean_square, type = "raw" )
    
    cat("Root mean squared dropout loss due to a feature elimination through model-agnostic approach in LDA:\n")
    print(fitlda_HVKPMDR)
  }
  
  else if (method == "model agnostic-Sig SVM"){
    set.seed(1234)
    tune_out_HVKPMDR <- e1071::tune("svm", `HVKP&MDR` ~ Imipenem + Gentamicin + Ciprofloxacin + 
                                      Meropenem + Aztreonam + Cefotaxime + Cefepime + Ceftazidime + 
                                      Ceftriaxone + Amikacin
                                    , data = train_strat_HVKPMDR, kernel = "sigmoid", tunecontrol=tune.control(cross=10),
                                    ranges = list( cost = seq(0.1, 5, length = 20)
                                                   , gamma = seq(0.1, 5, length = 20)))
    
    model_svm_HVKPMDR = svm(`HVKP&MDR` ~ Imipenem + Gentamicin + Ciprofloxacin + 
                              Meropenem + Aztreonam + Cefotaxime + Cefepime + Ceftazidime + 
                              Ceftriaxone + Amikacin  , data = train_strat_HVKPMDR, 
                            probability = T, kernel = "sigmoid", cost = tune_out_HVKPMDR$best.parameters$cost,
                            gamma = tune_out_HVKPMDR$best.parameters$gamma )
    
    explained_SVM_HVKPMDR <- explain(model = model_svm_HVKPMDR, data=train_strat_HVKPMDR[,c(2:9,11:12)], variables = colnames(train_strat_HVKPMDR[,c(2:9,11:12)]),
                                     y=as.vector(as.numeric(train_strat_HVKPMDR$`HVKP&MDR`))-1, label = "SVM-Sigmoid", 
                                     type = "classification")
    
    #50 permuatation
    fi_SVM_HVKPMDR = variable_importance(explained_SVM_HVKPMDR, B=50 ,variables = colnames(train_strat_HVKPMDR[,c(2:9,11:12)]),loss_function = loss_root_mean_square, type = "raw")
    
    cat("Root mean squared dropout loss due to a feature elimination through model-agnostic approach in sigmoid kernel of SVM:\n")
    print(fi_SVM_HVKPMDR)
  }
  
  else if (method == "model agnostic-Poly SVM"){
    
    set.seed(1234)
    tune_out3_HVKPMDR <- e1071::tune("svm", `HVKP&MDR` ~ Imipenem + Gentamicin + Ciprofloxacin + 
                                       Meropenem + Aztreonam + Cefotaxime + Cefepime + Ceftazidime + 
                                       Ceftriaxone + Amikacin
                                     , data = train_strat_HVKPMDR, kernel = "polynomial", tunecontrol=tune.control(cross=10),
                                     ranges = list( cost = seq(0.1, 5, length = 20), degree = c(1,2,3)
                                                    , gamma = seq(0.1, 5, length = 20)))
    
    model_svm2_HVKPMDR = svm(`HVKP&MDR` ~ Imipenem + Gentamicin + Ciprofloxacin + 
                               Meropenem + Aztreonam + Cefotaxime + Cefepime + Ceftazidime + 
                               Ceftriaxone + Amikacin
                             , data = train_strat_HVKPMDR, probability = T, 
                             kernel = "polynomial", cost = tune_out3_HVKPMDR$best.parameters$cost, 
                             degree = tune_out3_HVKPMDR$best.parameters$degree , gamma = tune_out3_HVKPMDR$best.parameters$gamma )
    
    explained_SVM2_HVKPMDR <- explain(model = model_svm2_HVKPMDR, data=train_strat_HVKPMDR[,c(2:9,11:12)], variables = colnames(train_strat_HVKPMDR[,c(2:9,11:12)]),
                                      y=as.vector(as.numeric(train_strat_HVKPMDR$`HVKP&MDR`))-1, label = "SVM-Polynomial", 
                                      type = "classification")
    
    #50 permuatation
    fi_SVM2_HVKPMDR = variable_importance(explained_SVM2_HVKPMDR, B=50 ,variables = colnames(train_strat_HVKPMDR[,c(2:9,11:12)]),loss_function = loss_root_mean_square, type = "raw")
    
    cat("Root mean squared dropout loss due to a feature elimination through model-agnostic approach in polynomial kernel of SVM:\n")
    print(fi_SVM2_HVKPMDR)
  }
  
  else if (method == "model agnostic-CatBoost"){
    
    # CatBoost Models--------------------------------------------------------------------
    Bio_Data <- as.data.frame(Bio_Data)
    
    #Convert categorical variables to factor----------------------------------
    categorical_var <- c(colnames(Bio_Data))
    
    Bio_Data[, categorical_var[-1]] <- lapply(Bio_Data[, categorical_var[-1]], factor)
    
    levels(Bio_Data$Ampicillin) <- c(levels(Bio_Data$Ampicillin),0)
    
    #Divide Data set into Train and Test---------------------------------------
    set.seed(123)
    
    split_strat_HVKPMDR  <- initial_split(Bio_Data, prop = 0.8, strata = "HVKP&MDR")
    
    train_strat_HVKPMDR  <- training(split_strat_HVKPMDR)
    
    test_strat_HVKPMDR   <- testing(split_strat_HVKPMDR)
    
    grid_HVKPMDDR <- expand.grid(
      depth = c(2,3),
      learning_rate = c(0.1,0.2),
      l2_leaf_reg = c(1,2),
      rsm = 1,
      border_count = 1,
      iterations = seq(130,151,1)
    )
    fitControl_HVKPMDR <- trainControl(method = "cv",
                                       number = 10,
                                       classProbs = TRUE)
    set.seed(1234)
    model_catboost_HVKPMDR <- train(x = train_strat_HVKPMDR[,c("Cefotaxime", "Gentamicin","Imipenem", "Cefepime", "Amikacin",
                                                               "Ceftazidime", "Ciprofloxacin", "Meropenem", "Ceftriaxone", "Aztreonam")],
                                    y = make.names(train_strat_HVKPMDR[,"HVKP&MDR"]),
                                    maximize = TRUE,
                                    method = catboost.caret, metric = "Accuracy", 
                                    tuneGrid =  grid_HVKPMDDR,
                                    trControl = fitControl_HVKPMDR)
    
    expsda_HVKPMDR = explain(model = model_catboost_HVKPMDR, data=train_strat_HVKPMDR[,c("Cefotaxime", "Gentamicin","Imipenem", "Cefepime", "Amikacin",
                                                                                         "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")], variables = colnames(train_strat_HVKPMDR[,c("Cefotaxime", "Gentamicin","Imipenem", "Cefepime", "Amikacin",
                                                                                                                                                                                                               "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")]),
                             y=as.vector(as.numeric(train_strat_HVKPMDR$`HVKP&MDR`))-1, label = "CatBoost", 
                             type = "classification")
    
    fitcat_HVKPMDR = variable_importance(expsda_HVKPMDR, B=50 ,variables = colnames(train_strat_HVKPMDR[,c("Cefotaxime", "Gentamicin","Imipenem", "Cefepime", "Amikacin",
                                                                                                           "Ceftazidime", "Ciprofloxacin","Meropenem", "Ceftriaxone", "Aztreonam")]), loss_function = loss_root_mean_square, type = "raw" )
    
    cat("Root mean squared dropout loss due to a feature elimination through model-agnostic approach in CatBoost:\n")
    print(fitcat_HVKPMDR)
  }
  else {
    cat("Invalid task. Please specify one of the following tasks: feature_selectio, model_evaluation\n")
    
  } 
}
if (is.null(task)){
  cat("You should specify a task\n")
}
