library(dplyr)
library(mice)
library(readxl)
library(survival)
library(parallel)
library(doParallel)
library(caret)
library(openxlsx)

ICD10 <- read.csv("/dssg/home/acct-sunkun/duxihao-1/UKB/HUA/ICD10.csv")
# incident depression
create_disease_event <- function(data, icd_pattern, disease_name) {
  detect_event <- function(i) {
    disease_time <- numeric(0)
    for (j in 0:258) {
      column_name <- paste0("s_41270_0_", j)
      if (!is.na(data[[column_name]][i]) && grepl(icd_pattern, data[[column_name]][i])) {
        disease_time <- c(disease_time, data[[paste0("s_41280_0_", j)]][i])
      }
    }
    
    if (length(disease_time) > 0) {
      time <- (min(disease_time, na.rm = TRUE) - data$s_53_0_0[i]) / 365
      event <- 1
    } else if (!is.na(data$s_40000_0_0[i])) {
      time <- (data$s_40000_0_0[i] - data$s_53_0_0[i]) / 365
      event <- 0
    } else {
      time <- (23303 - data$s_53_0_0[i]) / 365
      event <- 0
    }
    
    return(c(event, time))
  }
  
  results <- mclapply(1:nrow(data), detect_event, mc.cores = 5)
  results_df <- do.call(rbind, results)
  
  data[[paste0("event_", disease_name)]] <- results_df[, 1]
  data[[paste0("time_", disease_name)]] <- results_df[, 2]
  
  return(data)
}
#CAD
CAD_pattern <- "^I20[0-9]|^I21[0-9]|^I22[0-9]|^I23[0-9]|^I24[0-9]|^I25[0-9]"
ICD10 <- create_disease_event(ICD10, CAD_pattern, "cad")

#stroke
stroke_pattern <- "^I60[0-9]|^I61[0-9]|^I62[0-9]|^I63[0-9]|^I64|^I69[0-9]"
ICD10 <- create_disease_event(ICD10, stroke_pattern, "stroke")

#HF
HF_pattern <- "^I50[0-9]"
ICD10 <- create_disease_event(ICD10, HF_pattern, "hf")

#Hypertension
Hypertension_pattern <- "^I10|^I11[0-9]|^I12[0-9]|^I13[0-9]|^I14[0-9]|^I15[0-9]"
ICD10 <- create_disease_event(ICD10, Hypertension_pattern, "hyperten")

#arrhythmias
Arrhythmias_pattern <- "^I44[0-9]|^I45[0-9]|^I46[0-9]|^I47[0-9]|^I48[0-9]|^I49[0-9]"
ICD10 <- create_disease_event(ICD10, Arrhythmias_pattern, "arrhyth")

#T2D
T2D_pattern <- "^E11[0-9]"
ICD10 <- create_disease_event(ICD10, T2D_pattern, "t2d")

# Renal failure
Renalfail_pattern <- "^N17[0-9]|^N18[0-9]|^N19[0-9]"
ICD10 <- create_disease_event(ICD10, Renalfail_pattern, "ckd")

# gout
Gout_pattern <- "^M10[0-9]"
ICD10 <- create_disease_event(ICD10, Gout_pattern, "gout")

# death
event_death <- numeric(nrow(ICD10))
time_death <- numeric(nrow(ICD10))

for (i in 1:nrow(ICD10)) {
  if (!is.na(ICD10$s_40000_0_0[i])) {
    event_death[i] <- 1
    time_death[i] <- (ICD10$s_40000_0_0[i] - ICD10$s_53_0_0[i]) / 365
  } else {
    event_death[i] <- 0
    time_death[i] <- (23315 - ICD10$s_53_0_0[i]) / 365
  }
}
ICD10$time_death <- time_death
ICD10$event_death <- event_death

data <- ICD10[,c(1,744:761)]
Urate <- read.csv("/dssg/home/acct-sunkun/duxihao-1/UKB/HUA/Urate.csv")
data1 <- merge(data,Urate,by="n_eid")
write.csv(data1,"/dssg/home/acct-sunkun/duxihao-1/UKB/HUA/HUAdiseases.csv",row.names = F)
#导入数据及合并数据集
HUAanalysis <- read.csv("/dssg/home/acct-sunkun/duxihao-1/UKB/HUA/HUAdiseases.csv")
protein <-read.csv("/dssg/home/acct-sunkun/duxihao-1/UKB/HUA/protein_UKB_filled.csv")
covariates <- read.csv("/dssg/home/acct-sunkun/duxihao-1/UKB/HUA/covariates.csv")
HUAanalysis <- merge(HUAanalysis, covariates, by="n_eid")

HUAanalysis$HUA <- ifelse(HUAanalysis$Uriacid > 420,1,0)
HUAanalysis$HUA[is.na(HUAanalysis$HUA)] = 0  
HUAanalysis <- subset(HUAanalysis,HUA==1)
HUAanalysis <- HUAanalysis %>%
  select(-c(HUA))
##多重插补协变量
covariates1 <- c("Mets", "Ethnicity", "Sleeptime", "Edu", "Employed", "Smoke", "Drink", "Tdi", "Dietscore", 
                 "Age", "Sex", "BMI") 
data_to_impute <- HUAanalysis[, covariates1]
imputed_data <- mice(data_to_impute, method = 'pmm', m = 10, seed = 123)
completed_data <- complete(imputed_data)
HUAanalysis[, covariates1] <- completed_data
factvar <- c("Ethnicity", "Edu", "Employed", "Smoke", "Drink", "Sex")
HUAanalysis[factvar] <- lapply(HUAanalysis[factvar], as.factor)
write.csv(HUAanalysis,"/dssg/home/acct-sunkun/duxihao-1/UKB/HUA/HUA_datasetimputed.csv",row.names = F)

dput(names(HUAanalysis))
library(tableone)
myVars <-c("Age","Sex", "Ethnicity", "Edu", "Employed", "Tdi", "Dietscore","Mets", 
           "Sleeptime","Smoke", "Drink", "BMI","event_cad",
           "event_stroke","event_hf","event_hyperten","event_arrhyth",
           "event_t2d","event_ckd","event_gout","event_death")
chavariables <- c("Sex", "Ethnicity", "Edu", "Employed","Smoke","Drink","event_cad",
                  "event_stroke","event_hf","event_hyperten","event_arrhyth",
                  "event_t2d","event_ckd","event_gout","event_death")
tab2 <- CreateTableOne(vars = myVars, data = HUAanalysis,strata = "event_gout", factorVars = chavariables,addOverall = TRUE)
tabone <- print(tab2, contDigits = 2,showAllLevels = TRUE, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
tabone
write.csv(tabone,"/dssg/home/acct-sunkun/duxihao-1/UKB/HUA/baseline.csv")

#发病人数
total_n <- nrow(HUAanalysis)
results <- data.frame(
  Group = c("CAD", "Stroke", "HF", "Arrhythmia", "hypertension","T2D", "CKD", "Gout","Death"),
  N = c(
    sum(HUAanalysis$event_cad == 1),
    sum(HUAanalysis$event_stroke == 1),
    sum(HUAanalysis$event_hf == 1),
    sum(HUAanalysis$event_arrhyth == 1),
    sum(HUAanalysis$event_hyperten == 1),
    sum(HUAanalysis$event_t2d == 1),
    sum(HUAanalysis$event_ckd == 1),
    sum(HUAanalysis$event_gout == 1),
    sum(HUAanalysis$event_death == 1)
  ),
  Percentage = c(
    round(sum(HUAanalysis$event_cad == 1) / total_n * 100, 1),
    round(sum(HUAanalysis$event_stroke == 1) / total_n * 100, 1),
    round(sum(HUAanalysis$event_hf == 1) / total_n * 100, 1),
    round(sum(HUAanalysis$event_arrhyth == 1) / total_n * 100, 1),
    round(sum(HUAanalysis$event_hyperten == 1) / total_n * 100, 1),
    round(sum(HUAanalysis$event_t2d == 1) / total_n * 100, 1),
    round(sum(HUAanalysis$event_ckd == 1) / total_n * 100, 1),
    round(sum(HUAanalysis$event_gout == 1) / total_n * 100, 1),
    round(sum(HUAanalysis$event_death == 1) / total_n * 100, 1)
  )
)
results$Output <- paste(results$N, " (", results$Percentage, ")", sep = "")
print(results[, c("Group", "Output")])


#蛋白组关联分析
HUAprotein <- merge(HUAanalysis, protein, by.x="n_eid",by.y="eid")
HUAprotein[, 38:2957] <- as.data.frame(scale(HUAprotein[, 38:2957]))

wb <- createWorkbook()

analyze_disease <- function(disease, data, protein_cols) {
  result <- data.frame(Disease = character(), Variable = character(), Estimate = character(), Pvalue = numeric(), stringsAsFactors = FALSE)
  for(protein in protein_cols) {
    formula <- as.formula(paste0("Surv(time_", disease, ", event_", disease,")~ ", protein, " + Age + Sex + Mets + Ethnicity + Sleeptime + Edu + Employed + Smoke + Drink + Tdi + Dietscore + BMI"))
    model <- coxph(formula, data = data)
    summary_table <- summary(model)
    esti <- summary_table$conf.int[1, 1]
    upper <- summary_table$conf.int[1, 4]
    down <- summary_table$conf.int[1, 3]
    pvalue <- summary_table$coefficients[1, 5]
    ci <- sprintf("%.2f (%.2f to %.2f)", esti, down, upper)
    result <- rbind(result, data.frame(Disease = disease, Variable = protein, Estimate = ci, Pvalue = pvalue))
  }
  result$Fdr_pvalue <- p.adjust(result$Pvalue, method = "fdr")
  result$Bonferroni_pvalue <- p.adjust(result$Pvalue, method = "bonferroni")
  return(result)
}

diseases <- c("cad", "stroke", "hf", "arrhyth", "hyperten", "t2d", "ckd", "gout", "death")
protein_cols <- colnames(HUAprotein)[38:2957]

for(disease in diseases) {
  data_subset <- subset(HUAprotein, get(paste0("time_", disease)) > 0)
  results <- analyze_disease(disease, data_subset, protein_cols)
  addWorksheet(wb, disease)
  writeData(wb, disease, results)
  setColWidths(wb, disease, cols = 1:ncol(results), widths = "auto")
  addFilter(wb, disease, row = 1, cols = 1:ncol(results))
}

saveWorkbook(wb, "/dssg/home/acct-sunkun/duxihao-1/UKB/HUA/Signatures/Protein_Disease_Associations.xlsx", overwrite = TRUE)



#代谢组学
meta <- read.csv("/dssg/home/acct-sunkun/duxihao-1/UKB/HUA/metabolism_base_filled.csv")
HUAmeta <- merge(HUAanalysis, meta, by.x="n_eid",by.y="eid")

HUAmeta[, 38:288] <- as.data.frame(scale(HUAmeta[, 38:288]))
wb <- createWorkbook()
analyze_disease <- function(disease, data, meta_cols) {
  result <- data.frame(Disease = character(), Variable = character(), Estimate = character(), Pvalue = numeric(), stringsAsFactors = FALSE)
  for(meta in meta_cols) {
    formula <- as.formula(paste0("Surv(time_", disease, ", event_", disease,")~ ", meta, " + Age + Sex + Mets + Ethnicity + Sleeptime + Edu + Employed + Smoke + Drink + Tdi + Dietscore + BMI"))
    model <- coxph(formula, data = data)
    summary_table <- summary(model)
    esti <- summary_table$conf.int[1, 1]
    upper <- summary_table$conf.int[1, 4]
    down <- summary_table$conf.int[1, 3]
    pvalue <- summary_table$coefficients[1, 5]
    ci <- sprintf("%.2f (%.2f to %.2f)", esti, down, upper)
    result <- rbind(result, data.frame(Disease = disease, Variable = meta, Estimate = ci, Pvalue = pvalue))
  }
  result$Fdr_pvalue <- p.adjust(result$Pvalue, method = "fdr")
  result$Bonferroni_pvalue <- p.adjust(result$Pvalue, method = "bonferroni")
  return(result)
}

diseases <- c("cad", "stroke", "hf", "arrhyth", "hyperten", "t2d", "ckd", "gout", "death")
meta_cols <- colnames(HUAmeta)[38:288]

for(disease in diseases) {
  data_subset <- subset(HUAmeta, get(paste0("time_", disease)) > 0)
  results <- analyze_disease(disease, data_subset, meta_cols)
  addWorksheet(wb, disease)
  writeData(wb, disease, results)
  setColWidths(wb, disease, cols = 1:ncol(results), widths = "auto")
  addFilter(wb, disease, row = 1, cols = 1:ncol(results))
}

saveWorkbook(wb, "/dssg/home/acct-sunkun/duxihao-1/UKB/HUA/Signatures/Metabolomics_Disease_Associations.xlsx", overwrite = TRUE)


#血生化
bio <- read.csv("/dssg/home/acct-sunkun/duxihao-1/UKB/HUA/biochemistry.csv")
####缺失情况
complete_cases <- colSums(!is.na(bio))
complete_cases_df <- data.frame(
  Variable = names(complete_cases),
  N = complete_cases,
  Percentage = round(complete_cases/nrow(bio) * 100, 2)
)
complete_cases_df <- complete_cases_df[order(-complete_cases_df$N),]
bio <- bio %>%
  select(-c(n_30660_0_0, n_30790_0_0, n_30800_0_0, n_30820_0_0, n_30880_0_0))

biovars <- colnames(bio[,2:55])
data_to_impute <- bio[, biovars]
imputed_data <- mice(data_to_impute, method = 'pmm', m = 10, seed = 1234)
completed_data <- complete(imputed_data)
bio[, biovars] <- completed_data

HUAbio <- merge(HUAanalysis, bio, by="n_eid")
HUAbio[, 38:91] <- as.data.frame(scale(HUAbio[, 38:91]))
write.csv(HUAbio,"/dssg/home/acct-sunkun/duxihao-1/UKB/HUA/HUAbiodataset.csv",row.names = F)
bio_cols <- colnames(HUAbio)[c(38:91)]

wb <- createWorkbook()

analyze_disease <- function(disease, data, bio_cols) {
  result <- data.frame(Disease = character(), Variable = character(), Estimate = character(), Pvalue = numeric(), stringsAsFactors = FALSE)
  for(bio in bio_cols) {
    formula <- as.formula(paste0("Surv(time_", disease, ", event_", disease,")~ ", bio, " + Age + Sex + Mets + Ethnicity + Sleeptime + Edu + Employed + Smoke + Drink + Tdi + Dietscore + BMI"))
    model <- coxph(formula, data = data)
    summary_table <- summary(model)
    esti <- summary_table$conf.int[1, 1]
    upper <- summary_table$conf.int[1, 4]
    down <- summary_table$conf.int[1, 3]
    pvalue <- summary_table$coefficients[1, 5]
    ci <- sprintf("%.2f (%.2f to %.2f)", esti, down, upper)
    result <- rbind(result, data.frame(Disease = disease, Variable = bio, Estimate = ci, Pvalue = pvalue))
  }
  result$Fdr_pvalue <- p.adjust(result$Pvalue, method = "fdr")
  result$Bonferroni_pvalue <- p.adjust(result$Pvalue, method = "bonferroni")
  return(result)
}

diseases <- c("cad", "stroke", "hf", "arrhyth", "hyperten", "t2d", "ckd", "gout", "death")

for(disease in diseases) {
  data_subset <- subset(HUAbio, get(paste0("time_", disease)) > 0)
  results <- analyze_disease(disease, data_subset, bio_cols)
  addWorksheet(wb, disease)
  writeData(wb, disease, results)
  setColWidths(wb, disease, cols = 1:ncol(results), widths = "auto")
  addFilter(wb, disease, row = 1, cols = 1:ncol(results))
}

saveWorkbook(wb, "/dssg/home/acct-sunkun/duxihao-1/UKB/HUA/Signatures/Biochemical_Disease_Associations.xlsx", overwrite = TRUE)

#基因
HUAanalysis <- read.csv("/dssg/home/acct-sunkun/duxihao-1/UKB/HUA/HUA_datasetimputed.csv")
prs <- read.csv("/dssg/home/acct-sunkun/duxihao-1/UKB/HUA/PRS_HUA.csv")
HUAprs <- merge(HUAanalysis, prs, by="n_eid")
HUAprs[, 38:46] <- as.data.frame(scale(HUAprs[, 38:46]))

prs_cols <- colnames(HUAprs)[38:46]
wb <- createWorkbook()

analyze_disease <- function(disease, data, prs_cols) {
  result <- data.frame(Disease = character(), Variable = character(), Estimate = character(), Pvalue = numeric(), stringsAsFactors = FALSE)
  for(prs in prs_cols) {
    formula <- as.formula(paste0("Surv(time_", disease, ", event_", disease,")~ ", prs, " + Age + Sex + Ethnicity"))
    model <- coxph(formula, data = data)
    summary_table <- summary(model)
    esti <- summary_table$conf.int[1, 1]
    upper <- summary_table$conf.int[1, 4]
    down <- summary_table$conf.int[1, 3]
    pvalue <- summary_table$coefficients[1, 5]
    ci <- sprintf("%.2f (%.2f to %.2f)", esti, down, upper)
    result <- rbind(result, data.frame(Disease = disease, Variable = prs, Estimate = ci, Pvalue = pvalue))
  }
  result$Fdr_pvalue <- p.adjust(result$Pvalue, method = "fdr")
  result$Bonferroni_pvalue <- p.adjust(result$Pvalue, method = "bonferroni")
  return(result)
}

diseases <- c("cad", "stroke", "hf", "arrhyth", "hyperten", "t2d", "ckd", "gout", "death")

for(disease in diseases) {
  data_subset <- subset(HUAprs, get(paste0("time_", disease)) > 0)
  results <- analyze_disease(disease, data_subset, prs_cols)
  addWorksheet(wb, disease)
  writeData(wb, disease, results)
  setColWidths(wb, disease, cols = 1:ncol(results), widths = "auto")
  addFilter(wb, disease, row = 1, cols = 1:ncol(results))
}

saveWorkbook(wb, "/dssg/home/acct-sunkun/duxihao-1/UKB/HUA/Signatures/PRS_Disease_Associations.xlsx", overwrite = TRUE)




wb <- loadWorkbook("/dssg/home/acct-sunkun/duxihao-1/UKB/HUA/Signatures/Biochemical_Disease_Associations.xlsx")
sheets <- getSheetNames("/dssg/home/acct-sunkun/duxihao-1/UKB/HUA/Signatures/Biochemical_Disease_Associations.xlsx")

significant_associations <- data.frame(Variables = character(), Diseases = character())

for(sheet in sheets) {
  data <- read.xlsx(wb, sheet)
  sig_vars <- data %>%
    filter(Fdr_pvalue < 0.05) %>%
    select(Variable)
  
  if(nrow(sig_vars) > 0) {
    temp_df <- data.frame(
      Variables = sig_vars$Variable,
      Diseases = rep(sheet, nrow(sig_vars))
    )
    significant_associations <- rbind(significant_associations, temp_df)
  }
}

write.csv(significant_associations, "/dssg/home/acct-sunkun/duxihao-1/UKB/HUA/Signatures/BioSign.csv", row.names = FALSE)
