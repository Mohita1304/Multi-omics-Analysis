# Insatll and load the following libraries

library("survival")
library("survminer")
library('pROC')
library('rpart')

# Survival prediction using the risk score for each sample

# Read the data
# Training datasets
df <- read.csv("C:/Users/RS_c.csv") 
 # external test datasets 1
data <- read.csv("C:/User/RS_m.csv")
# external test dataset 2
data1 <- read.csv("C:/User/RS_m2.csv") 
df1 <- df[, 2:4]
df2 <- data[, 2:4]
df3 <- data1[, 2:4]

# Define the proportion for the split
train_proportion <- 0.8
n_train <- round(nrow(df1) * train_proportion)

# Specify the best seed found from the previous analysis
best_seed <- 89 # Replace this with the seed found from the previous analysis

# Set the seed
set.seed(best_seed)

# Shuffle the data
df1 <- df1[sample(nrow(df1)), ]

# Split the data into training and testing sets
train_data <- df1[1:n_train, ]
test_data <- df1[(n_train + 1):nrow(df1), ]

# Fit Cox proportional hazards model on training data
surv_model <- coxph(Surv(OS_MONTH, OS_STATUS) ~ ., data = train_data)

# Predict survival probabilities for test data
surv_prob1 <- predict(surv_model, newdata = test_data, type = "risk")

# Predict survival probabilities for df2
surv_prob2 <- predict(surv_model, newdata = df2, type = "risk")

# Predict survival probabilities for df3
surv_prob3 <- predict(surv_model, newdata = df3, type = "risk")

# Define time points for ROC curve calculation (1 year = 12 months)
time_point <- 60

# Calculate survival status at the specified time point for test data
status_at_time_test <- ifelse(test_data$OS_MONTH <= time_point & test_data$OS_STATUS == 1, 1, 0)

# Calculate survival status at the specified time point for df2
status_at_time_df2 <- ifelse(df2$OS_MONTH <= time_point & df2$OS_STATUS == 1, 1, 0)

# Calculate survival status at the specified time point for df3
status_at_time_df3 <- ifelse(df3$OS_MONTH <= time_point & df3$OS_STATUS == 1, 1, 0)

# Calculate ROC curve for test data
roc_test <- roc(status_at_time_test, 1 - surv_prob1, levels = c(0, 1))
auc_test <- auc(roc_test)

# Calculate ROC curve for df2
roc_df2 <- roc(status_at_time_df2, 1 - surv_prob2, levels = c(0, 1))
auc_df2 <- auc(roc_df2)

# Calculate ROC curve for df3
roc_df3 <- roc(status_at_time_df3, 1 - surv_prob3, levels = c(0, 1))
auc_df3 <- auc(roc_df3)

# Plot ROC curves
plot(roc_test, main = "ROC Curves 5 year", col = "red", lwd = 2)
lines(roc_df2, col = "blue", lwd = 2)
lines(roc_df3, col = "green", lwd = 2)
legend("bottomright", legend = c(paste("cBioPortal AUC =", round(auc_test, 2)), 
                                 paste("GSE29623 AUC =", round(auc_df2, 2)), 
                                 paste("GSE17536 AUC =", round(auc_df3, 2))), 
       col = c("red", "blue", "green"), lwd = 2)

# Print the AUC values
print(paste("Test data AUC:", auc_test))
print(paste("DF2 AUC:", auc_df2))
print(paste("DF3 AUC:", auc_df3))

