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
diseases <- c("cad","stroke","hf","arrhyth","hyperten","t2d","ckd","gout","death")
gene <- read.csv("/dssg/home/acct-sunkun/duxihao-1/UKB/HUA/PRS_HUA.csv")

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
  Protein <- read.csv(paste0("/dssg/home/acct-sunkun/duxihao-1/UKB/HUA/Signatures/Top20/Protescore_top20_", disease, ".csv"))
  Bio <- read.csv(paste0("/dssg/home/acct-sunkun/duxihao-1/UKB/HUA/Signatures/Top20/biocore_top20_", disease, ".csv"))
  Meta <- read.csv(paste0("/dssg/home/acct-sunkun/duxihao-1/UKB/HUA/Signatures/Top20/Metacore_top20_", disease, ".csv"))
  
  HUApro <- merge(HUAanalysis, Protein, by.x="n_eid",by.y="ID")
  HUAbiopro <- merge(HUApro, Bio, by.x="n_eid",by.y="ID")
  HUAbioprometa <- merge(HUAbiopro, Meta, by.x="n_eid",by.y="ID")
  HUAbioprometagene <- merge(HUAbioprometa, gene, by="n_eid")
  HUAbioprometagene[, 41:49] <- as.data.frame(scale(HUAbioprometagene[, 41:49]))
  
  model0 <- coxph(as.formula(paste0("Surv(time_", disease, ",event_", disease,") ~ Age+Sex+Ethnicity")), data = HUAbioprometagene)
  model1 <- coxph(as.formula(paste0("Surv(time_", disease, ",event_", disease,") ~ Age+Sex+Ethnicity+Tdi+Edu")), data = HUAbioprometagene)
  model2 <- coxph(as.formula(paste0("Surv(time_", disease, ",event_", disease,") ~ Age+Sex+Ethnicity+Tdi+Edu+BMI+Smoke+Drink+Dietscore+Mets+BMI")), data = HUAbioprometagene)
  model3 <- coxph(as.formula(paste0("Surv(time_", disease, ",event_", disease,") ~ Age+Sex+Ethnicity+Tdi+Edu+BMI+Smoke+Drink+Dietscore+Mets+BMI+Uriacid")), data = HUAbioprometagene)
  model4 <- coxph(as.formula(paste0("Surv(time_", disease, ",event_", disease,") ~ Age+Sex+Ethnicity+Tdi+Edu+BMI+Smoke+Drink+Dietscore+Mets+BMI+Uriacid+PRS_",disease,"+bioscore+Metascore+Protescore")), data = HUAbioprometagene)
  
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
      FUNChange(model0, model4, HUAbioprometagene),
      FUNChange(model1, model4, HUAbioprometagene),
      FUNChange(model2, model4, HUAbioprometagene),
      FUNChange(model3, model4, HUAbioprometagene)
    )
  )
  
  addWorksheet(wb, disease)
  writeData(wb, disease, "C-index Results", startCol = 1, startRow = 1)
  writeData(wb, disease, c_stats, startCol = 1, startRow = 2)
  writeData(wb, disease, "C-index Changes", startCol = 1, startRow = nrow(c_stats) + 4)
  writeData(wb, disease, changes, startCol = 1, startRow = nrow(c_stats) + 5)
  
  setColWidths(wb, disease, cols = 1:3, widths = "auto")
}

saveWorkbook(wb, "/dssg/home/acct-sunkun/duxihao-1/UKB/HUA/Signatures/Top20/Multiomics_Model_Evaluation_Results.xlsx", overwrite = TRUE)
