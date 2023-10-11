# Antimicrobial Resistance Gene and Pathotype Prediction (R)

This repository contains code and resources for predicting and classifying antimicrobial resistance (AMR) genes and pathotypes based on antibiotic susceptibility data through R language. We employ multiple machine learning algorithms, including Logistic Regression (LR), Naive Bayes Classifier (NBC), Linear Discriminant Analysis (LDA), Support Vector Machines (SVM), and CatBoost, for the prediction and classification tasks. Additionally, we apply three feature selection methods, including Pearson's Chi-Squared test, Wald test, and a model-agnostic approach using the DALEX package, to identify the most important antibiotics for each gene and pathotype prediction

## Table of Contents

- [Introduction](#introduction)
- [Features](#features)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Results](#results)
- [Contributing](#contributing)
- [License](#license)

## Introduction

Antimicrobial resistance is a pressing global health concern. This project leverages R and a range of machine learning techniques to predict and classify AMR genes and pathotypes using antibiotic susceptibility data. Our objective is to provide a robust tool that helps identify important antibiotics associated with specific genes and pathotypes.

## Features

- **ML Algorithms**: We utilize various machine learning algorithms, such as LR, NBC, LDA, SVM, and CatBoost, to perform gene and pathotype prediction and classification tasks.

- **Feature Selection**: We offer three feature selection methods, including Pearson's Chi-Squared test, Wald test, and a model-agnostic approach through the DALEX package, to identify crucial antibiotics for gene and pathotype detection.

- **Evaluation Metrics**: Due to our goal for predicting presence of genes or pathotypes, thus, Sensitivity, precision, and F1-Score are included to assess the performance of our models, ensuring reliable detections. However, other metrcis like Accuracy and Specificity are reported.

## Requirements

The method is developed in R (4.2.2) and the following packages can be installed as follow:

```bash
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
```

## Installation

To get started with this repository, clone the repository to your local machine:

   ```bash
    git clone https://github.com/kouroshAKiani/ML-GeneticPathotyping.git
   ```

## Usage
### 1. Preparation:
 Download the dataset (mach-learn.xlsx) and the code files in your directory. Change your work directory to the path that the files and dataset have been saved:

```bash
    cd /path/to/your_directory
   ```

### 2. Feature Selection:
 Choose one of the feature selection methods provided in the repository and run it on your dataset to select the most important antibiotics as follow:

```bash
   Rscript gene/pathotype name.R -f mach-learn.xlsx -t "feature_selection" -m "feature selection method"
```

For example if you want to inspect model-agnostic approach in LR for prediction of blaOXA-48, you should execute following command in your command-line (terminal):
 
```bash
   Rscript blaOXA-48.R -f mach-learn.xlsx -t "feature_selection" -m "model agnostic-LR"
```
You can examine other feature selection methods by replacing one of the following items with "model agnostic-LR" in the above-mentioned command.

Chi-Squared test, Wald test, model agnostic-NBC, model agnostic-LDA, model agnostic-Sig SVM, model agnostic-Poly SVM, model agnostic-CatBoost

for example if you want to review the results of Pearson's Chi-Squared test in the case of detection of blaSHV gene, you should execute following command in your terminal:

```bash
   Rscript blaSHV.R -f mach-learn.xlsx -t "feature_selection" -m "Chi-Squared test"
```

### 2. Model training and evaluation:
In order to investigate ML models performance that are reported in tables 2 to 11 in the article, you should choose the algorithm and the feature selection method that results in intended antibiotic combination in the algorithm. For this purpose you should consider following command in your command-line:


   
   

