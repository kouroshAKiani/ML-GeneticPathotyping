# Antimicrobial Resistance Gene and Pathotype Prediction

This repository contains code and resources for predicting and classifying antimicrobial resistance (AMR) genes and pathotypes based on antibiotic susceptibility data. We employ multiple machine learning algorithms, including Logistic Regression (LR), Naive Bayes Classifier (NBC), Linear Discriminant Analysis (LDA), Support Vector Machines (SVM), and CatBoost, for the prediction and classification tasks. Additionally, we apply three feature selection methods, including Pearson's Chi-Squared test, Wald test, and a model-agnostic approach using the DALEX package, to identify the most important antibiotics for each gene and pathotype prediction

## Table of Contents

- [Introduction](#introduction)
- [Features](#features)
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
   Rscript gene / pathotype name.R -f mach-learn.xlsx -t "feature_selection" -m "feature selection method"
```

For example if you want to inspect model-agnostic approach in LR for prediction of blaOXA-48, you should execute following command in your command-line (terminal):
 
```bash
   Rscript blaOXA-48.R -f mach-learn.xlsx -t "feature_selection" -m "model agnostic-LR"
```



   
   

