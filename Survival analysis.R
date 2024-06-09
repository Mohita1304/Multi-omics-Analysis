# univariate Cox proportional hazards regression analysis was used to screen for genes significantly associated with CRC prognosis

library("survival")
library("survminer")
library('ggplot2')
data_folder="/Path/" # the gene expression files with clinical inforamtion, column are genes and row are the sampels
library(readr)
my_files<-list.files(data_folder, pattern = "\\.csv$")
t=0
for (csv_file in my_files) {
  file_path <- file.path(data_folder, csv_file)
  data <- read.csv(file_path, header = TRUE, sep = ",")
  res.cox2 <- coxph(Surv(data$OS_DAYS, data$OS_STATUS) ~ RSEM, data =  data)
  A <- summary(res.cox2)
  P_value=A$logtest["pvalue"]
  if (P_value< 0.05){
    t=t+1
    p<-substr(csv_file, 1, nchar(csv_file) - 4)
    coefficients_df <- data.frame(
      Variable = rownames(A$coefficients),
      Coefficient = A$coefficients[, "coef"], HR = A$coefficients[, "exp(coef)"],
      Standard_Error = A$coefficients[, "se(coef)"],
      Z_Value = A$coefficients[, "z"],
      P_Value_coef = A$coefficients[, "Pr(>|z|)"],
      P_value = A$logtest["pvalue"]
    )
    z<-sprintf("/folder path/% s", csv_file) # 
    write.csv(coefficients_df, file = z)} save the gene expression and clinical infroamtion for individual genes
}

# Least absolute shrinkage and selection operator (LASSO) Cox analysis

setwd("/Path/OS/")
library(glmnet)
d<-read.csv("Gene expression_survialInfo.csv", row.names = 1) # the file containg time, status and gene expression
View(d)
d1=as.matrix(d)
x=d[, -1:-2]
View(x)
x=as.matrix(x)
y=d[, 1:2]
y=as.matrix(y)

fit <- glmnet(x, y, family = "cox")
cv_model <- cv.glmnet(x, y, family = "cox", alpha = 1)
best_lambda <- cv_model$lambda.min
best_lambda
selected_predictors <- which(coef(cv_model, s = best_lambda) != 0)
selected_predictor_names <- colnames(x)[selected_predictors]

# Extract coefficients and selected predictor names
selected_coefficients <- coef(cv_model, s = best_lambda)[selected_predictors]
selected_predictor_names <- colnames(x)[selected_predictors]

# Create dataframe
selected_predictors_df <- data.frame(
  Predictor = selected_predictor_names,
  Coefficient = selected_coefficients
)
write.csv(selected_predictors_df, file="Final_selected_LASSO_Gene_All_Samples.csv")


# Multivariate Cox proportional hazards regression analysis
d<-read.csv("Gene expression_survialInfo.csv", row.names = 1)  ## the file containg time, status and gene expression
res.cox <- coxph(Surv(d$time, d$status) ~ Gene1+Gene2+..., data =  d)
A <- summary(res.cox3)
ggsurvplot(survfit(Surv(d$time, d$status) ~ Gene1+Gene2+..., data =  d), pval=TRUE, risk.table = TRUE, legend="right")

