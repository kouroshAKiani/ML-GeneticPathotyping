# Antimicrobial Resistance Gene and Pathotype Prediction

This repository contains code and resources for predicting and classifying antimicrobial resistance (AMR) genes and pathotypes based on antibiotic susceptibility data. We employ multiple machine learning algorithms, including Logistic Regression (LR), Naive Bayes Classifier (NBC), Linear Discriminative Analysis, Support Vector Machines (SVM), and CatBoost, for the prediction and classification tasks. Additionally, we apply three feature selection methods, including Pearson's Chi-Squared test, Wald test, and a model-agnostic approach using the DALEX package, to identify the most important antibiotics for each gene and pathotype prediction

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

- **ML Algorithms**: We utilize various machine learning algorithms, such as LR, NBC, Linear Discriminative Analysis, SVM, and CatBoost, to perform gene and pathotype prediction and classification tasks.

- **Feature Selection**: We offer three feature selection methods, including Pearson's Chi-Squared test, Wald test, and a model-agnostic approach through the DALEX package, to identify crucial antibiotics for gene and pathotype prediction.

- **Evaluation Metrics**: Due to our goal for predicting presence of genes or pathotypes, thus, Sensitivity, precision, and F1-Score are included to assess the performance of our models, ensuring reliable predictions. However, other metrcis like Accuracy and Specificity are reported.

## Installation

To get started with this repository, follow these installation steps:

1. Clone the repository to your local machine:

   ```bash
   git clone https://github.com/kouroshAKiani/ML-GeneticPathotyping.git

