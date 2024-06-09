library("survival")
library("survminer")
library('pROC')
library('rpart')
setwd("/home/sukanta-mondal/Documents/cBioPortal/Multi_Omics_Analysis/CNA_mRNA/OS_SA/DEG_CNA_OS/ROC/")
df<-read.csv("Four_Genes.csv", header = TRUE, sep=",", row.names = 1)
#View(df)
AUC_Value=c()
i=0
for (i in i:1001) {
  #df<-read.csv("Four_Genes.csv", header = TRUE, sep=",", row.names = 1)
  
  #df<-read.csv("Two_Genes.csv", header = TRUE, sep=",", row.names = 1)
  
  # Separate features (X) and target variable (y)
  X <- df[, !(names(df) %in% c('OS_MONTHS', 'OS_STATUS'))]
  y <- df[, c('OS_STATUS', 'OS_MONTHS')]
  
  # Split the data into training and testing sets
  #set.seed(42)  # For reproducibility
  train_index <- sample(nrow(df), 0.8 * nrow(df))
  X_train <- X[train_index, ]
  X_test <- X[-train_index, ]
  y_train <- y[train_index, ]
  y_test <- y[-train_index, ]
  
  # Example: Training your survival model
  surv_model <- coxph(Surv(OS_MONTHS, OS_STATUS) ~ ., data = data.frame(X_train, y_train))
  #surv_model <- rpart(Surv(OS_MONTHS, OS_STATUS) ~ ., data = data.frame(X_train, y_train))
  
  # Example: Getting predictions for the test data
  surv_prob <- predict(surv_model, newdata = data.frame(X_test), se.fit = TRUE)
  
  # Example: Calculating time-dependent ROC curves for 1, 3, and 5 years
  time_points <- c(12)  # Time points in months
  for (time_point in time_points) {
    # Calculate survival status at the specified time point
    status_at_time <- ifelse(y_test$OS_MONTHS <= time_point & y_test$OS_STATUS == 1, 1, 0)
    
    # Calculate sensitivity and specificity for different thresholds
    #roc_data <- roc(status_at_time, 1 - surv_prob$fit, levels=c(0,1))
    #roc_data <- roc(status_at_time, 1 - surv_prob$, levels = c(0, 1))
    if (sum(status_at_time) > 0) {
      roc_data <- roc(status_at_time, 1 - surv_prob$fit, levels=c(0,1))
    }
    else
    {
      next
    }
    # Plot ROC curve
    # z<-sprintf('/home/sukanta-mondal/Documents/cBioPortal/Multi_Omics_Analysis/CNA_mRNA/OS_SA/DEG_CNA_OS/ROC/COX/OneYear/% s.png', i)
    #png(z, width = 800, height = 600)
    # plot(roc_data, main=paste("ROC Curve at", time_point, "months"))
    
    
    
    # Calculate AUC and display it on the plot
    auc_val <- auc(roc_data)
    
    #print(auc_val)
    AUC_Value= c(AUC_Value, auc_val)
    #text(0.8, 0.2, paste("AUC =", round(auc_val, 2)))
    #dev.off()
  }
}
print(AUC_Value)
RPART_AUC<-as.data.frame(AUC_Value)
write.csv(RPART_AUC, file="/home/sukanta-mondal/Documents/cBioPortal/Multi_Omics_Analysis/CNA_mRNA/OS_SA/DEG_CNA_OS/ROC/COX/COX_AUC_1years_4.csv")