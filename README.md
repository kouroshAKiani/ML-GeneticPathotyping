# Antimicrobial Resistance Gene and Pathotype Prediction

This repository contains code and resources for predicting and classifying antimicrobial resistance (AMR) genes and pathotypes based on antibiotic susceptibility data through R. Moreover, the data distribution is visualized through Python. We employ multiple machine learning algorithms, including Logistic Regression (LR), Naive Bayes Classifier (NBC), Linear Discriminant Analysis (LDA), Support Vector Machines (SVM), and CatBoost, for the prediction and classification tasks. Additionally, we apply three feature selection methods, including Pearson's Chi-Squared test, Wald test, and a model-agnostic approach using the DALEX package, to identify the most important antibiotics for each gene and pathotype detection.

## Table of Contents

- [Introduction](#introduction)
- [Features](#features)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Results](#results)
- [Authors](#authors)
- [License](#license)
- [Reference](#reference)

## Introduction

Antimicrobial resistance is a pressing global health concern. This project leverages a range of machine learning techniques to predict and classify AMR genes and pathotypes using antibiotic susceptibility data. Our objective is to provide a robust tool that helps identify important antibiotics associated with specific genes and pathotypes.

## Features

- **ML Algorithms**: We utilize various machine learning algorithms, such as LR, NBC, LDA, SVM, and CatBoost, to perform gene and pathotype prediction and classification tasks.

- **Feature Selection**: We offer three feature selection methods, including Pearson's Chi-Squared test, Wald test, and a model-agnostic approach through the DALEX package, to identify crucial antibiotics for gene and pathotype detection.

- **Evaluation Metrics**: Due to our goal for predicting presence of genes or pathotypes, thus, Sensitivity, precision, and F1-Score are included to assess the performance of our models, ensuring reliable detections. However, other metrcis like Accuracy and Specificity are reported but they are not considered for models evaluations and decision making.

## Requirements

The method is developed in R (4.2.2) and Python (3.9). You should save Packages_installation.R (for required packages installations), compute_metrics.R (for model evaluation metrics) , and F1-score.R (for confidence interval computation of F1-Score) files in the same directory that other genes/pathotypes files that are going to be saved. 

## Installation

To get started with this repository, clone the repository to your local machine:

   ```bash
    git clone https://github.com/kouroshAKiani/ML-GeneticPathotyping.git
   ```

## Usage
### 1. Preparation:
 Download the dataset (mach-learn.xlsx) and save it in the same directory that the other files have been saved. Change your work directory to the path that the files and dataset have been saved:

```bash
    cd /path/to/your_directory
   ```

### 2. Feature Selection:
 Choose one of the feature selection methods provided in the repository and run it on your dataset to select the most important antibiotics as follow:

```bash
   Rscript gene/pathotype name.R -f mach-learn.xlsx -t "feature_selection" -m "feature selection method"
```

For example if you want to examine model-agnostic approach in LR for prediction of blaOXA-48, you should execute following command in your command-line (terminal):
 
```bash
   Rscript blaOXA-48.R -f mach-learn.xlsx -t "feature_selection" -m "model agnostic-LR"
```
You can examine other feature selection methods by replacing one of the following items with "model agnostic-LR" in the above-mentioned command.

Chi-Squared test, Wald test, model agnostic-NBC, model agnostic-LDA, model agnostic-Sig SVM, model agnostic-Poly SVM, model agnostic-CatBoost

for example if you want to review the results of Pearson's Chi-Squared test in the case of detection of blaSHV gene, you should execute following command in your terminal:

```bash
   Rscript blaSHV.R -f mach-learn.xlsx -t "feature_selection" -m "Chi-Squared test"
```

### 3. Model training and evaluation:
In order to investigate ML models performance that are reported in tables 2 to 11 in the article, you should choose the algorithm and the feature selection methods based on the intended antibiotic(s) combination. For this purpose you should follow below command in your command-line:

```bash
   Rscript blaSHV.R -f mach-learn.xlsx -t "model_evaluation" -o "algorithm name-feature selection method name"
```

It should be noted that feature selection methods names have been introduced before, in part 2 (feature selection). Moreoever, if there are two methods resulting in the same antibiotic(s) combination (according to tables 2 to 11), then you can consider following objects to replace with feature selection method name:

if model-agnostic approach and Wald result in the same antibiotic(s) combination, you should replace "feature selection method name" with "model agnostic/Wald test". Moreover, if Chi-Squared test and the Wald test results are the same, then you should replace "feature selection method name" with "Chi Squared test/Wald test".

For example if the antibiotic(s) combination selected through the model-agnostic approach in LR and the Wald test are the same, you can follow below command:

```bash
   Rscript blaSHV.R -f mach-learn.xlsx -t "model_evaluation" -o "LR-model agnostic/Wald test"
```
We should emphasize that computation of confidence interval for F1-Score for multi-class classification has been computed by [1]. By getting inspiration from this work and applying some changes to their computations, we can extract a function for computation of F1-Score and its confidence interval in a binary classification that is provided in F1-score.R file.

## Results

Performance of all the models based on the sensitivity, precision, and F1-Score are represented in the tables 2 to 11. Additionally, feature importance analysis based on the model-agnostic approach (through DALEX package), residual distribution, precision-recall curves, and the data distribution figures, reported in the article, are available in feature-importance.R, 3-plots.R, and data-disturb.py files, respectively.

## Authors

This project was made possible thanks to the contributions of the following individuals:

- [Kourosh Alizadeh Kiani](https://github.com/kouroshAKiani)
- [Mehdi Soroush](https://github.com/MehdiSoroush)
- [Seyed Ahmad Sanikhani](https://github.com/AhmadSanikhani)
- [Sajad Tavakoli](https://github.com/sajadtavakoli) 

Special thanks to all the contributors for their hard work and dedication to this project.

## License

This project is licensed under the MIT License.
   
## Reference
1. Takahashi K, Yamamoto K, Kuchiba A, Koyama T. Confidence interval for micro-averaged F1 and macro-averaged F1 scores.  Appl Intell (Dordr). 2022 Mar;52(5):4961-4972. doi: 10.1007/s10489-021-02635-5. Epub 2021 Jul 31.


