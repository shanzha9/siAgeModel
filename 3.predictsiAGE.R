library(randomForest)
library(data.table)
library(tibble)

# load model
randomforest_model <- readRDS("./all_selectedModel.rds")

# load data
data <- fread("./cp10k.csv")
data <- column_to_rownames(data, var = "V1")
data_nm <- as.data.frame(t(data_nm))

# feature selection
validation_data <- data[, rownames(randomforest_model$importance)]
validation_data[] <- apply(validation_data, 2, as.numeric)

# predict
predict_siage <- predict(randomforest_model, validation_data)
predict_siage <- as.data.frame(predict_siage)
predict_siage$samplename <- rownames(predict_siage)

# save
write.csv(predict_siage, "./siAGE_prediction.csv", row.names = F)

