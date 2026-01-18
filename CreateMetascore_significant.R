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
meta <- read.csv("/dssg/home/acct-sunkun/duxihao-1/UKB/HUA/metabolism_base_filled.csv")
MetaSign <- read.csv("/dssg/home/acct-sunkun/duxihao-1/UKB/HUA/Signatures/MetaSign.csv")

#蛋白组
HUAmeta <- merge(HUAanalysis, meta, by.x="n_eid",by.y="eid")
HUAmeta[, 38:288] <- as.data.frame(scale(HUAmeta[, 38:288]))

diseases <- c("cad","stroke","hf","arrhyth","hyperten","t2d","ckd","gout","death")
HUAmeta_cad <- subset(HUAmeta,time_cad>0)
HUAmeta_stroke <- subset(HUAmeta,time_stroke>0)
HUAmeta_hf <- subset(HUAmeta,time_hf>0)
HUAmeta_hyperten <- subset(HUAmeta,time_hyperten>0)
HUAmeta_arrhyth <- subset(HUAmeta,time_arrhyth>0)
HUAmeta_t2d <- subset(HUAmeta,time_t2d>0)
HUAmeta_ckd <- subset(HUAmeta,time_ckd>0)
HUAmeta_gout <- subset(HUAmeta,time_gout>0)
HUAmeta_death<- subset(HUAmeta,time_death>0)

num_cores <- 20
cl <- makeCluster(num_cores)
registerDoParallel(cl)
all_disease_results <- list()

for (disease in diseases) {
  data_name <- paste0("HUAmeta_", disease)
  current_data <- get(data_name)
  protein_cols <- MetaSign %>%
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
      scores$Metascore <- as.vector(predict(model, newx = test_x, type = "link"))
      
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
            paste0("/dssg/home/acct-sunkun/duxihao-1/UKB/HUA/Signatures/Metacore_", disease, ".csv"), 
            row.names = FALSE)
}

wb <- createWorkbook()
for (disease in diseases) {
  addWorksheet(wb, sheetName = disease)
  writeData(wb, sheet = disease, x = all_disease_results[[disease]]$betas, rowNames = TRUE)
  setColWidths(wb, sheet = disease, cols = 1:(ncol(all_disease_results[[disease]]$betas)+1), widths = "auto")
}
saveWorkbook(wb, file = "/dssg/home/acct-sunkun/duxihao-1/UKB/HUA/Signatures/Metascore_beta.xlsx", overwrite = TRUE)


