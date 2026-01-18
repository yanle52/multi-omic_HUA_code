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

HUAanalysis <- read.csv("/dssg/home/acct-sunkun/duxihao-1/UKB/HUA/HUA_datasetimputed.csv")
gene <- read.csv("/dssg/home/acct-sunkun/duxihao-1/UKB/HUA/PRS_HUA.csv")
#蛋白组
HUAgene <- merge(HUAanalysis, gene, by="n_eid")
HUAgene[, 38:46] <- as.data.frame(scale(HUAgene[, 38:46]))

diseases <- c("cad","stroke","hf","arrhyth","hyperten","t2d","ckd","gout","death")
HUAgene_cad <- subset(HUAgene,time_cad>0)
HUAgene_stroke <- subset(HUAgene,time_stroke>0)
HUAgene_hf <- subset(HUAgene,time_hf>0)
HUAgene_hyperten <- subset(HUAgene,time_hyperten>0)
HUAgene_arrhyth <- subset(HUAgene,time_arrhyth>0)
HUAgene_t2d <- subset(HUAgene,time_t2d>0)
HUAgene_ckd <- subset(HUAgene,time_ckd>0)
HUAgene_gout <- subset(HUAgene,time_gout>0)
HUAgene_death <- subset(HUAgene,time_death>0)

FUNC <- function(model) {
  c <- round(summary(model)$concordance[1], 3)
  c.lwr <- round(summary(model)$concordance[1] - 1.96 * summary(model)$concordance[2], 3)
  c.upr <- round(summary(model)$concordance[1] + 1.96 * summary(model)$concordance[2], 3)
  paste0(c, " (", c.lwr, ", ", c.upr, ")")
}

FUNChange <- function(model0, model1,data) {
  modelcompare <- CsChange(model0, model1, data=data, nb=200)
  change <- round(modelcompare[[1]]$change, 3)
  change.lwr <- round(modelcompare[[1]]$low, 3)
  change.upr <- round(modelcompare[[1]]$up, 3)
  paste0(change, " (", change.lwr, ", ", change.upr, ")")
}

wb <- createWorkbook()

for(disease in diseases) {
  data_name <- paste0("HUAgene_", disease)
  current_data <- get(data_name)
  
  model0 <- coxph(as.formula(paste0("Surv(time_", disease, ",event_", disease,") ~ Age+Sex+Ethnicity")), data = current_data)
  model1 <- coxph(as.formula(paste0("Surv(time_", disease, ",event_", disease,") ~ Age+Sex+Ethnicity+Tdi+Edu")), data = current_data)
  model2 <- coxph(as.formula(paste0("Surv(time_", disease, ",event_", disease,") ~ Age+Sex+Ethnicity+Tdi+Edu+BMI+Smoke+Drink+Dietscore+Mets+BMI")), data = current_data)
  model3 <- coxph(as.formula(paste0("Surv(time_", disease, ",event_", disease,") ~ Age+Sex+Ethnicity+Tdi+Edu+BMI+Smoke+Drink+Dietscore+Mets+BMI+Uriacid")), data = current_data)
  model4 <- coxph(as.formula(paste0("Surv(time_", disease, ",event_", disease,") ~ Age+Sex+Ethnicity+Tdi+Edu+BMI+Smoke+Drink+Dietscore+Mets+BMI+Uriacid+PRS_",disease)), data = current_data)
  
  c_stats <- data.frame(
    Model = c("Base", "SES", "Lifestyle", "UA", "Full"),
    C_index = c(
      FUNC(model0),
      FUNC(model1),
      FUNC(model2),
      FUNC(model3),
      FUNC(model4)
    )
  )
  
  changes <- data.frame(
    Comparison = c("Base vs Full", "SES vs Full", "Lifestyle vs Full", "UA vs Full"),
    Change = c(
      FUNChange(model0, model4, current_data),
      FUNChange(model1, model4, current_data),
      FUNChange(model2, model4, current_data),
      FUNChange(model3, model4, current_data)
    )
  )
  
  addWorksheet(wb, disease)
  writeData(wb, disease, "C-index Results", startCol = 1, startRow = 1)
  writeData(wb, disease, c_stats, startCol = 1, startRow = 2)
  writeData(wb, disease, "C-index Changes", startCol = 1, startRow = nrow(c_stats) + 4)
  writeData(wb, disease, changes, startCol = 1, startRow = nrow(c_stats) + 5)
  
  setColWidths(wb, disease, cols = 1:3, widths = "auto")
}

saveWorkbook(wb, "/dssg/home/acct-sunkun/duxihao-1/UKB/HUA/Signatures/Top20/Gene_Model_Evaluation_Results.xlsx", overwrite = TRUE)
