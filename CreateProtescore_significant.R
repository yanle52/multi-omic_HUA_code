library(dplyr)
library(readxl)
library(survival)
library(glmnet)
library(parallel)
library(doParallel)
library(caret)
library(openxlsx)
library(CsChange)
library(survIDINRI)
#导入数据及合并数据集
HUAanalysis <- read.csv("/dssg/home/acct-sunkun/duxihao-1/UKB/HUA/HUA_datasetimputed.csv")
protein <-read.csv("/dssg/home/acct-sunkun/duxihao-1/UKB/HUA/protein_UKB_filled.csv")
covariates <- read.csv("/dssg/home/acct-sunkun/duxihao-1/UKB/HUA/covariates.csv")
ProSign <- read.csv("/dssg/home/acct-sunkun/duxihao-1/UKB/HUA/Signatures/ProteinSign.csv")

#蛋白组
HUAprotein <- merge(HUAanalysis, protein, by.x="n_eid",by.y="eid")
HUAprotein[, 38:2957] <- as.data.frame(scale(HUAprotein[, 38:2957]))

diseases <- c("cad","stroke","hf","arrhyth","hyperten","t2d","ckd","gout","death")
HUAprotein_cad <- subset(HUAprotein,time_cad>0)
HUAprotein_stroke <- subset(HUAprotein,time_stroke>0)
HUAprotein_hf <- subset(HUAprotein,time_hf>0)
HUAprotein_hyperten <- subset(HUAprotein,time_hyperten>0)
HUAprotein_arrhyth <- subset(HUAprotein,time_arrhyth>0)
HUAprotein_t2d <- subset(HUAprotein,time_t2d>0)
HUAprotein_ckd <- subset(HUAprotein,time_ckd>0)
HUAprotein_gout <- subset(HUAprotein,time_gout>0)
HUAprotein_death<- subset(HUAprotein,time_death>0)

num_cores <- 20
cl <- makeCluster(num_cores)
registerDoParallel(cl)
all_disease_results <- list()

for (disease in diseases) {
  data_name <- paste0("HUAprotein_", disease)
  current_data <- get(data_name)
  protein_cols <- ProSign %>%
    filter(Diseases == disease) %>%
    pull(Variables)
  outer_folds <- createFolds(current_data$n_eid, k = 10, list = TRUE, returnTrain = FALSE)
  results <- foreach(i = 1:10, .packages = c("glmnet", "survival", "caret", "dplyr")) %dopar% {
    tryCatch({
      test_indices <- outer_folds[[i]]
      train_indices <- setdiff(1:nrow(current_data), test_indices)
      
      train.data <- current_data[train_indices, ]
      test.data <- current_data[test_indices, ]
      
      scores <- data.frame(ID = current_data$n_eid[test_indices])
      
      x <- as.matrix(train.data[, protein_cols])
      y <- Surv(train.data[[paste0("time_", disease)]], 
                train.data[[paste0("event_",disease)]])
      
      set.seed(123)
      cv <- cv.glmnet(x, y, family = 'cox', alpha = 1, nfolds = 10)
      model <- glmnet(x, y, family = 'cox', alpha = 1, lambda = cv$lambda.min)
      
      fold_betas <- coef(model)[, 1]
      
      test_x <- as.matrix(test.data[, protein_cols])
      scores$Protescore <- as.vector(predict(model, newx = test_x, type = "link"))
      
      list(scores = scores, betas = fold_betas)
      
    }, error = function(e) {
      cat("Error in", disease, "fold", i, ":", conditionMessage(e), "\n")
      return(NULL)
    })
  }
  
  all_predictions <- list()
  beta_list <- list()
  for (i in 1:10) {
    if (!is.null(results[[i]])) {
      all_predictions[[i]] <- results[[i]]$scores
      beta_list[[i]] <- results[[i]]$betas
    }
  }
  
  all_predictions_df <- do.call(rbind, all_predictions)
  
  disease_betas <- do.call(cbind, beta_list)
  colnames(disease_betas) <- paste0("Fold", 1:10)
  rownames(disease_betas) <- protein_cols
  
  all_disease_results[[disease]] <- list(
    predictions = all_predictions_df,
    betas = disease_betas
  )
}

for (disease in diseases) {
  write.csv(all_disease_results[[disease]]$predictions, 
            paste0("/dssg/home/acct-sunkun/duxihao-1/UKB/HUA/Signatures/Protecore_", disease, ".csv"), 
            row.names = FALSE)
}

wb <- createWorkbook()
for (disease in diseases) {
  addWorksheet(wb, sheetName = disease)
  writeData(wb, sheet = disease, x = all_disease_results[[disease]]$betas, rowNames = TRUE)
  setColWidths(wb, sheet = disease, cols = 1:(ncol(all_disease_results[[disease]]$betas)+1), widths = "auto")
}
saveWorkbook(wb, file = "/dssg/home/acct-sunkun/duxihao-1/UKB/HUA/Signatures/Protescore_beta.xlsx", overwrite = TRUE)


# #读取蛋白质数据文件
# merge_proteins <- function(protein_data, prote_file, eid_col = "n_eid", id_col = "ID") {
#   prote_data <- read.csv(prote_file)
#   merged_data <- merge(protein_data, prote_data, by.x = eid_col, by.y = id_col)
#   return(merged_data)
# }
# HUAprotein_cad1 <- merge_proteins(HUAprotein_cad, "/dssg/home/acct-sunkun/duxihao-1/UKB/HUA/Signatures/Protecore_cad.csv")
# HUAprotein_stroke1 <- merge_proteins(HUAprotein_stroke, "/dssg/home/acct-sunkun/duxihao-1/UKB/HUA/Signatures/Protecore_stroke.csv")
# HUAprotein_hf1 <- merge_proteins(HUAprotein_hf, "/dssg/home/acct-sunkun/duxihao-1/UKB/HUA/Signatures/Protecore_hf.csv")
# HUAprotein_hyperten1 <- merge_proteins(HUAprotein_hyperten, "/dssg/home/acct-sunkun/duxihao-1/UKB/HUA/Signatures/Protecore_hyperten.csv")
# HUAprotein_arrhyth1 <- merge_proteins(HUAprotein_arrhyth, "/dssg/home/acct-sunkun/duxihao-1/UKB/HUA/Signatures/Protecore_arrhyth.csv")
# HUAprotein_t2d1 <- merge_proteins(HUAprotein_t2d, "/dssg/home/acct-sunkun/duxihao-1/UKB/HUA/Signatures/Protecore_t2d.csv")
# HUAprotein_ckd1 <- merge_proteins(HUAprotein_ckd, "/dssg/home/acct-sunkun/duxihao-1/UKB/HUA/Signatures/Protecore_ckd.csv")
# HUAprotein_gout1 <- merge_proteins(HUAprotein_gout, "/dssg/home/acct-sunkun/duxihao-1/UKB/HUA/Signatures/Protecore_gout.csv")
# HUAprotein_death1 <- merge_proteins(HUAprotein_death, "/dssg/home/acct-sunkun/duxihao-1/UKB/HUA/Signatures/Protecore_death.csv")
# 
# 
# FUNC <- function(model) {
#   c <- round(summary(model)$concordance[1], 3)
#   c.lwr <- round(summary(model)$concordance[1] - 1.96 * summary(model)$concordance[2], 3)
#   c.upr <- round(summary(model)$concordance[1] + 1.96 * summary(model)$concordance[2], 3)
#   paste0(c, " (", c.lwr, ", ", c.upr, ")")
# }
# 
# FUNChange <- function(model0, model1,data) {
#   modelcompare <- CsChange(model0, model1, data=data, nb=200)
#   change <- round(modelcompare[[1]]$change, 3)
#   change.lwr <- round(modelcompare[[1]]$low, 3)
#   change.upr <- round(modelcompare[[1]]$up, 3)
#   paste0(change, " (", change.lwr, ", ", change.upr, ")")
# }
# 
# wb <- createWorkbook()
# 
# for(disease in diseases) {
#   data_name <- paste0("HUAprotein_", disease, "1")
#   current_data <- get(data_name)
#   
#   model0 <- coxph(as.formula(paste0("Surv(time_", disease, ",event_", disease,") ~ Age+Sex+Ethnicity")), data = current_data)
#   model1 <- coxph(as.formula(paste0("Surv(time_", disease, ",event_", disease,") ~ Age+Sex+Ethnicity+Tdi+Edu")), data = current_data)
#   model2 <- coxph(as.formula(paste0("Surv(time_", disease, ",event_", disease,") ~ Age+Sex+Ethnicity+Tdi+Edu+BMI+Smoke+Drink+Dietscore+Mets+BMI")), data = current_data)
#   model3 <- coxph(as.formula(paste0("Surv(time_", disease, ",event_", disease,") ~ Age+Sex+Ethnicity+Tdi+Edu+BMI+Smoke+Drink+Dietscore+Mets+BMI+Uriacid")), data = current_data)
#   model4 <- coxph(as.formula(paste0("Surv(time_", disease, ",event_", disease,") ~ Age+Sex+Ethnicity+Tdi+Edu+BMI+Smoke+Drink+Dietscore+Mets+BMI+Uriacid+Protescore")), data = current_data)
#   
#   c_stats <- data.frame(
#     Model = c("Base", "SES", "Lifestyle", "UA", "Full"),
#     C_index = c(
#       FUNC(model0),
#       FUNC(model1),
#       FUNC(model2),
#       FUNC(model3),
#       FUNC(model4)
#     )
#   )
#   
#   changes <- data.frame(
#     Comparison = c("Base vs Full", "SES vs Full", "Lifestyle vs Full", "UA vs Full"),
#     Change = c(
#       FUNChange(model0, model4, current_data),
#       FUNChange(model1, model4, current_data),
#       FUNChange(model2, model4, current_data),
#       FUNChange(model3, model4, current_data)
#     )
#   )
#   
#   addWorksheet(wb, disease)
#   writeData(wb, disease, "C-index Results", startCol = 1, startRow = 1)
#   writeData(wb, disease, c_stats, startCol = 1, startRow = 2)
#   writeData(wb, disease, "C-index Changes", startCol = 1, startRow = nrow(c_stats) + 4)
#   writeData(wb, disease, changes, startCol = 1, startRow = nrow(c_stats) + 5)
#   
#   setColWidths(wb, disease, cols = 1:3, widths = "auto")
# }
# 
# saveWorkbook(wb, "/dssg/home/acct-sunkun/duxihao-1/UKB/HUA/Signatures/Model_Evaluation_Results.xlsx", overwrite = TRUE)
