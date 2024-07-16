# The code implement of siAGE Prediction Model

## 1. Introduction

This is the code implementation of siAGE prediction model. The model is based on the randomForest, which is implemented by R language. The model is trained on the scRNA genes profiles dataset of 56 healthy individuals. The model can predict the immune age of the individual based on the scRNA data.  

### 2. Dependencies

This project requires the following packages:

- **R**: 4.3.0
- **Seurat**: 5.0.3
- **randomForest**: 1.0.10
- **ggplot2**: 3.5.1
- **data.table**: 1.15.1
- **tibble**: 3.2.1

Except for the Seurat package, the annotated reference data is also required. The reference data could be downloaded at: [syn61609846](https://www.synapse.org/Synapse:syn61609846).

## 3. Usage

### 3.1 CellType Prediction

To predict the cell of a new datasets, the siAge prediction model requires the same input features as the trained model. The new dataset should be preprocessed in the same way as the training dataset (shown in the manuscript). Then, predict the cell type of the new dataset using [predicteCellType.R](./1.predicteCellType.R).  

### 3.2 CP10K Calculation  

Calculate the CP10K of the new dataset using [calculateCP10K.R](./2.calculateCP10K.R).  

### 3.3 siAGE Prediction

Predict the siAGE of the new dataset using [predictsiAGE.R](./3.predictsiAGE.R).

