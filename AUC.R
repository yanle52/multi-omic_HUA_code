library(openxlsx)
library(survival)
library(riskRegression)
library(ggplot2)

HUAanalysis <- read.csv("/dssg/home/acct-sunkun/duxihao-1/UKB/HUA/HUA_datasetimputed.csv")
diseases <- c("cad","stroke","hf","arrhyth","hyperten","t2d","ckd","gout","death")
colors <- c('#7F7F80', '#55719D', '#5C996C', "#E2C7BA", '#FF9045', '#BF1E20', '#8B5E3C')
gene <- read.csv("/dssg/home/acct-sunkun/duxihao-1/UKB/HUA/PRS_HUA.csv")
results_auc <- list()
for(disease in diseases) {
  Protein <- read.csv(paste0("/dssg/home/acct-sunkun/duxihao-1/UKB/HUA/Signatures/Top20/Protescore_top20_", disease, ".csv"))
  Bio <- read.csv(paste0("/dssg/home/acct-sunkun/duxihao-1/UKB/HUA/Signatures/Top20/biocore_top20_", disease, ".csv"))
  Meta <- read.csv(paste0("/dssg/home/acct-sunkun/duxihao-1/UKB/HUA/Signatures/Top20/Metacore_top20_", disease, ".csv"))
  
  HUApro <- merge(HUAanalysis, Protein, by.x="n_eid",by.y="ID")
  HUAbiopro <- merge(HUApro, Bio, by.x="n_eid",by.y="ID")
  HUAbioprometa <- merge(HUAbiopro, Meta, by.x="n_eid",by.y="ID")
  HUAbioprometagene <- merge(HUAbioprometa, gene, by="n_eid")
  HUAbioprometagene[, 41:49] <- as.data.frame(scale(HUAbioprometagene[, 41:49]))
  
  is_zero_bioscore <- all(HUAbioprometagene$bioscore == 0)
  is_zero_metascore <- all(HUAbioprometagene$Metascore == 0)
  is_zero_protescore <- all(HUAbioprometagene$Protescore == 0)
  
  base_vars <- c("Age", "Sex", "Ethnicity", "Tdi", "Edu")
  lifestyle_vars <- c(base_vars, "BMI", "Smoke", "Drink", "Dietscore", "Mets")
  clinical_vars <- c(lifestyle_vars, "Uriacid")
  
  formula_ses <- reformulate(termlabels = base_vars, 
                             response = paste0("Surv(time_", disease, ",event_", disease, ")"))
  formula_lifestyle <- reformulate(termlabels = lifestyle_vars,
                                   response = paste0("Surv(time_", disease, ",event_", disease, ")"))
  formula_bio <- if(is_zero_bioscore) {
    reformulate(termlabels = clinical_vars,
                response = paste0("Surv(time_", disease, ",event_", disease, ")"))
  } else {
    reformulate(termlabels = c(clinical_vars, "bioscore"),
                response = paste0("Surv(time_", disease, ",event_", disease, ")"))
  }
  
  formula_gene <- reformulate(termlabels = c(clinical_vars, paste0("PRS_", disease)),
                              response = paste0("Surv(time_", disease, ",event_", disease, ")"))
  
  formula_meta <- if(is_zero_metascore) {
    reformulate(termlabels = c(clinical_vars, "bioscore"),
                response = paste0("Surv(time_", disease, ",event_", disease, ")"))
  } else {
    reformulate(termlabels = c(clinical_vars, "Metascore"),
                response = paste0("Surv(time_", disease, ",event_", disease, ")"))
  }
  
  formula_prote <- if(is_zero_protescore) {
    reformulate(termlabels = clinical_vars,
                response = paste0("Surv(time_", disease, ",event_", disease, ")"))
  } else {
    reformulate(termlabels = c(clinical_vars, "Protescore"),
                response = paste0("Surv(time_", disease, ",event_", disease, ")"))
  }
  
  full_vars <- clinical_vars
  if(!is_zero_protescore) full_vars <- c(full_vars, "Protescore")
  if(!is_zero_bioscore) full_vars <- c(full_vars, "bioscore")
  if(!is_zero_metascore) full_vars <- c(full_vars, "Metascore")
  full_vars <- c(full_vars, paste0("PRS_", disease))
  formula_full <- reformulate(termlabels = full_vars,
                              response = paste0("Surv(time_", disease, ",event_", disease, ")"))
  
  model_list <- list(
    SES_model = coxph(formula_ses, data = HUAbioprometagene, x = TRUE),
    Lifestyle_model = coxph(formula_lifestyle, data = HUAbioprometagene, x = TRUE),
    Clinical_model = coxph(formula_bio, data = HUAbioprometagene, x = TRUE),
    Genetics_model = coxph(formula_gene, data = HUAbioprometagene, x = TRUE),
    Metabolomics_model = coxph(formula_meta, data = HUAbioprometagene, x = TRUE),
    Proteomics_model = coxph(formula_prote, data = HUAbioprometagene, x = TRUE),
    Multiomics_model = coxph(formula_full, data = HUAbioprometagene, x = TRUE)
  )
  
  xs <- Score(model_list,
              formula = as.formula(paste0("Hist(time_", disease, ",event_", disease, ") ~ 1")),
              data = HUAbioprometagene, times = 15, plots = "roc", metrics = "auc")
  
  results_auc[[disease]] <- list(auc = xs[["AUC"]][["score"]][["AUC"]], 
                                 p_values = xs[["AUC"]][["contrasts"]][["p"]],
                                 roc_data = xs)
  
  pdf(file = paste0("/dssg/home/acct-sunkun/duxihao-1/UKB/HUA/Signatures/Top20/AUC/", disease, "_ROC.pdf"), 
      width = 8, height = 6)
  plotROC(xs, col = colors)
  dev.off()
}


wb <- createWorkbook()
addWorksheet(wb, "AUC Results")

# 修改 AUC 结果格式
format_auc_data <- function(results_auc) {
  data.frame(
    Disease = diseases,
    SES_model = sapply(results_auc, function(x) {
      auc_val <- x$roc_data$AUC$score$AUC[1]
      cilow <- x$roc_data$AUC$score$lower[1]
      cihigh <- x$roc_data$AUC$score$upper[1]
      sprintf("%.3f (%.3f to %.3f)", auc_val,cilow, cihigh)
    }),
    Lifestyle_model = sapply(results_auc, function(x) {
      auc_val <- x$roc_data$AUC$score$AUC[2]
      cilow <- x$roc_data$AUC$score$lower[2]
      cihigh <- x$roc_data$AUC$score$upper[2]
      sprintf("%.3f (%.3f to %.3f)", auc_val,cilow, cihigh)
    }),
    Clinical_model = sapply(results_auc, function(x) {
      auc_val <- x$roc_data$AUC$score$AUC[3]
      cilow <- x$roc_data$AUC$score$lower[3]
      cihigh <- x$roc_data$AUC$score$upper[3]
      sprintf("%.3f (%.3f to %.3f)", auc_val,cilow, cihigh)
    }),
    Genetics_model = sapply(results_auc, function(x) {
      auc_val <- x$roc_data$AUC$score$AUC[4]
      cilow <- x$roc_data$AUC$score$lower[4]
      cihigh <- x$roc_data$AUC$score$upper[4]
      sprintf("%.3f (%.3f to %.3f)", auc_val,cilow, cihigh)
    }),
    Metabolomics_model = sapply(results_auc, function(x) {
      auc_val <- x$roc_data$AUC$score$AUC[5]
      cilow <- x$roc_data$AUC$score$lower[5]
      cihigh <- x$roc_data$AUC$score$upper[5]
      sprintf("%.3f (%.3f to %.3f)", auc_val,cilow, cihigh)
    }),
    Proteomics_model = sapply(results_auc, function(x) {
      auc_val <- x$roc_data$AUC$score$AUC[6]
      cilow <- x$roc_data$AUC$score$lower[6]
      cihigh <- x$roc_data$AUC$score$upper[6]
      sprintf("%.3f (%.3f to %.3f)", auc_val,cilow, cihigh)
    }),
    Multiomics_model = sapply(results_auc, function(x) {
      auc_val <- x$roc_data$AUC$score$AUC[7]
      cilow <- x$roc_data$AUC$score$lower[7]
      cihigh <- x$roc_data$AUC$score$upper[7]
      sprintf("%.3f (%.3f to %.3f)", auc_val,cilow, cihigh)
    })
  )
}

addWorksheet(wb, "P Values")
format_p_values <- function(results_auc) {
  p_value_matrix <- matrix(nrow = length(diseases), ncol = 21)
  colnames(p_value_matrix) <- c(
    "SES_vs_Lifestyle", "SES_vs_Clinical", "SES_vs_Genetics",
    "SES_vs_Metabolomics", "SES_vs_Proteomics", "SES_vs_Multiomics",
    "Lifestyle_vs_Clinical", "Lifestyle_vs_Genetics", "Lifestyle_vs_Metabolomics",
    "Lifestyle_vs_Proteomics", "Lifestyle_vs_Multiomics", "Clinical_vs_Genetics",
    "Clinical_vs_Metabolomics", "Clinical_vs_Proteomics", "Clinical_vs_Multiomics",
    "Genetics_vs_Metabolomics", "Genetics_vs_Proteomics", "Genetics_vs_Multiomics",
    "Metabolomics_vs_Proteomics", "Metabolomics_vs_Multiomics", "Proteomics_vs_Multiomics"
  )
  rownames(p_value_matrix) <- diseases
  
  for(i in seq_along(diseases)) {
    pvals <- results_auc[[diseases[i]]]$roc_data$AUC$contrasts$p
    if(length(pvals) == 21) {
      p_value_matrix[i,] <- round(pvals, 3)
    } else {
      p_value_matrix[i,] <- NA
    }
  }
  
  data.frame(
    Disease = diseases,
    p_value_matrix
  )
}

auc_data <- format_auc_data(results_auc)
p_value_data <- format_p_values(results_auc)

colnames(p_value_data)[2:22] <- c(
  "SES_vs_Lifestyle", "SES_vs_Clinical", "SES_vs_Genetics",
  "SES_vs_Metabolomics", "SES_vs_Proteomics", "SES_vs_Multiomics",
  "Lifestyle_vs_Clinical", "Lifestyle_vs_Genetics", "Lifestyle_vs_Metabolomics",
  "Lifestyle_vs_Proteomics", "Lifestyle_vs_Multiomics", "Clinical_vs_Genetics",
  "Clinical_vs_Metabolomics", "Clinical_vs_Proteomics", "Clinical_vs_Multiomics",
  "Genetics_vs_Metabolomics", "Genetics_vs_Proteomics", "Genetics_vs_Multiomics",
  "Metabolomics_vs_Proteomics", "Metabolomics_vs_Multiomics", "Proteomics_vs_Multiomics"
)

writeData(wb, "AUC Results", auc_data)
writeData(wb, "P Values", p_value_data)
setColWidths(wb, 1, cols = 1:ncol(auc_data), widths = "auto")
setColWidths(wb, 2, cols = 1:ncol(p_value_data), widths = "auto")

saveWorkbook(wb, "/dssg/home/acct-sunkun/duxihao-1/UKB/HUA/Signatures/Top20/AUC/AUC_Analysis_Results.xlsx", overwrite = TRUE)



# 死亡
pdf(file = paste0("/dssg/home/acct-sunkun/duxihao-1/UKB/HUA/Signatures/Top20/AUC/Death_ROC.pdf"), onefile = FALSE)

Protein <- read.csv("/dssg/home/acct-sunkun/duxihao-1/UKB/HUA/Signatures/Top20/Protescore_top20_death.csv")
Bio <- read.csv("/dssg/home/acct-sunkun/duxihao-1/UKB/HUA/Signatures/Top20/biocore_top20_death.csv")
Meta <- read.csv("/dssg/home/acct-sunkun/duxihao-1/UKB/HUA/Signatures/Top20/Metacore_top20_death.csv")

HUApro <- merge(HUAanalysis, Protein, by.x="n_eid",by.y="ID")
HUAbiopro <- merge(HUApro, Bio, by.x="n_eid",by.y="ID")
HUAbioprometa <- merge(HUAbiopro, Meta, by.x="n_eid",by.y="ID")
HUAbioprometagene <- merge(HUAbioprometa, gene, by="n_eid")
HUAbioprometagene[, 41:48] <- as.data.frame(scale(HUAbioprometagene[, 41:48]))

# 构建不同模型公式
disease <- "death"
formula_ses <- reformulate(termlabels = c("Age", "Sex", "Ethnicity", "Tdi", "Edu"), 
                           response = paste0("Surv(time_", disease, ",event_", disease, ")"))
formula_lifestyle <- reformulate(termlabels = c("Age", "Sex", "Ethnicity", "Tdi", "Edu", "BMI", "Smoke","Drink","Dietscore","Mets"), 
                                 response = paste0("Surv(time_", disease, ",event_", disease, ")"))
formula_bio <- reformulate(termlabels = c("Age", "Sex", "Ethnicity", "Tdi", "Edu", "BMI", "Smoke","Drink","Dietscore","Mets","Uriacid", "bioscore"), 
                           response = paste0("Surv(time_", disease, ",event_", disease, ")"))
formula_gene <- reformulate(termlabels = c("Age", "Sex", "Ethnicity", "Tdi", "Edu", "BMI", "Smoke","Drink","Dietscore","Mets","Uriacid",paste0("PRS_",disease)), 
                            response = paste0("Surv(time_", disease, ",event_", disease, ")"))
formula_meta <- reformulate(termlabels = c("Age", "Sex", "Ethnicity", "Tdi", "Edu", "BMI", "Smoke","Drink","Dietscore","Mets","Uriacid","Metascore"), 
                            response = paste0("Surv(time_", disease, ",event_", disease, ")"))
formula_prote <- reformulate(termlabels = c("Age", "Sex", "Ethnicity", "Tdi", "Edu", "BMI", "Smoke","Drink","Dietscore","Mets","Uriacid","Protescore"), 
                             response = paste0("Surv(time_", disease, ",event_", disease, ")"))
formula_full <- reformulate(termlabels = c("Age", "Sex", "Ethnicity", "Tdi", "Edu", "BMI", "Smoke","Drink","Dietscore","Mets","Uriacid","Protescore", "bioscore", 
                                           "Metascore",paste0("PRS_",disease)), 
                            response = paste0("Surv(time_", disease, ",event_", disease, ")"))

# 拟合模型
fs1 <- coxph(formula_ses, data = HUAbioprometagene, x = TRUE)
fs2 <- coxph(formula_lifestyle, data = HUAbioprometagene, x = TRUE)
fs3 <- coxph(formula_bio, data = HUAbioprometagene, x = TRUE)
fs4 <- coxph(formula_gene, data = HUAbioprometagene, x = TRUE)
fs5 <- coxph(formula_meta, data = HUAbioprometagene, x = TRUE)
fs6 <- coxph(formula_prote, data = HUAbioprometagene, x = TRUE)
fs7 <- coxph(formula_full, data = HUAbioprometagene, x = TRUE)


# 计算 AUC
xs <- Score(list(Base_model = fs1, SES_model = fs2, Lifestyle_Score = fs3, MRI_Score = fs4),
            formula = as.formula(paste("Hist(time,status) ~ 1", sep = "")),
            data = final_cmd_dep_3, times = 15, plots = "roc", metrics = "auc")
resultAUC <- list()
resultp <- list()

# 存储 AUC 结果
resultAUC <- cbind(resultAUC, xs[["AUC"]][["score"]][["AUC"]])
resultp <- cbind(resultp, xs[["AUC"]][["contrasts"]][["p"]][1:11])

# 绘制 ROC 曲线
plotROC(xs, col = colors)

dev.off()


# 保存结果
write.csv(resultAUC, "/lustre/home/acct-wenze.zhong/jingxuanw/UKB/dataset/resultAUC.csv")
write.csv(resultp, "/lustre/home/acct-wenze.zhong/jingxuanw/UKB/dataset/resultp.csv")
