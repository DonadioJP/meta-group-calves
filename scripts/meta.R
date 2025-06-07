# =============================================
# META-ANALYSIS: CALF HOUSING SYSTEMS
# AUTHOR: JP DONADIO
# =============================================

library(metafor)
library(dplyr)
library(ggplot2)
library(tidyr)
library(readxl)
library(purrr)
library(stringr)
library(tidyverse)
library(flextable)

# ----------------------------
# DATA PREPARATION
# ----------------------------

meta <- read_excel("~/meta-group-calves/data/02_meta.xlsx")

outcome_data <- read_excel("~/meta-group-calves/data/03_outcomes.xlsx")
names(outcome_data)
str(outcome_data)
outcome_data$Variable <- as.factor(outcome_data$Variable)
outcome_data <- outcome_data %>%
  mutate(across(
    where(~ is.character(.)),
    as.numeric
  ))

full_data <- outcome_data %>%
  left_join(meta, by = "Study") %>%
  mutate(across(
    where(~ is.character(.)),
    as.factor
  ))

View(full_data)
names(full_data)
str(full_data)

# ----------------------------
# PRIMARY META-ANALYSES
# ----------------------------
analyze_variable <- function(data, var_name, smd = FALSE, min_studies = 3) {
  df <- data %>% 
    filter(Variable == var_name) %>%
    mutate(
      Treatment_SD = Treatment_SEM * sqrt(Treatment_N),
      Control_SD = Control_SEM * sqrt(Control_N),
      valid = !is.na(Treatment_Mean) & !is.na(Control_Mean) &
        !is.na(Treatment_SD) & !is.na(Control_SD) &
        Treatment_SD > 0 & Control_SD > 0 &
        Treatment_N > 1 & Control_N > 1
    ) %>%
    filter(valid)
  
  if(nrow(df) < min_studies) {
    warning(var_name, ": Only ", nrow(df), " valid studies")
    return(NULL)
  }
  
  tryCatch({
    esc <- escalc(
      measure = ifelse(smd, "SMD", "MD"),
      m1i = Treatment_Mean,
      sd1i = Treatment_SD,   
      n1i = Treatment_N,
      m2i = Control_Mean,
      sd2i = Control_SD,    
      n2i = Control_N,
      data = df
    )
    
    model <- rma(yi, vi, data = esc, method = "DL")
    funnel_test <- regtest(model)
    
    list(
      n_studies = nrow(df),
      model = model,
      funnel = funnel_test,
      data = esc
    )
  }, error = function(e) {
    message("Error in ", var_name, ": ", e$message)
    print(df %>% select(contains(c("Mean", "SD", "N"))))
    NULL
  })
}

# WMD Analyses (Performance/Health)
vars_wmd <- c("WW", "ADG", "DMI", "CI", 
              "Blood_glucose", "BUN", "TNF")

results_wmd <- lapply(vars_wmd, function(x) analyze_variable(full_data, x))
names(results_wmd) <- vars_wmd

# SMD Analyses (Behavior)
vars_smd <- c("Lying", "Standing", "Self-grooming", "Feeding", 
              "Playing", "Interacting_with_pen", "Consuming_milk")

results_smd <- lapply(vars_smd, function(x) analyze_variable(full_data, x, smd = TRUE))
names(results_smd) <- vars_smd

# ----------------------------
# REMOVING OUTLIERS
# ----------------------------
funnel_test <- function(results, smd = FALSE) {
  map_dfr(results, function(res) {
    data.frame(
      Variable = res$model$slab[1],
      Control_Mean = mean(res$data$Control_Mean),
      N_Studies = res$n_studies,
      Funnel_P = ifelse(res$funnel$pval < 0.001, "<0.001", 
                        sprintf("%.3f", res$funnel$pval))
    )
  }, .id = "Outcome")
}

#WMD variables

funnel_test(results_wmd)
#ADG, Blood Glucose

## ADG
# Detecting outliers in funnel plot
funnel(results_wmd$ADG$model, 
       xlab = "Effect Size (WMD)", 
       ylab = "Standard Error",
       main = "Funnel Plot with Study IDs")

text(x = results_wmd$ADG$model$yi, 
     y = sqrt(results_wmd$ADG$model$vi), 
     labels = paste0("(", 1:nrow(results_wmd$ADG$data), ")"), 
     pos = 3,  
     cex = 0.7, 
     col = "blue")

# Removing outliers
studies_to_remove <- c(11)

results_wmd$ADG$data <- results_wmd$ADG$data[-studies_to_remove, ]

#Update the data
results_wmd$ADG$model <- rma(yi, vi, 
                       data = results_wmd$ADG$data, 
                       method = "DL")

results_wmd$ADG$n_studies <- nrow(results_wmd$ADG$data)

results_wmd$ADG$funnel <- regtest(results_wmd$ADG$model)

# Checking the new funnel plot
regtest(results_wmd$ADG$model)

funnel(results_wmd$ADG$model, 
       xlab = "Effect Size (WMD)", 
       ylab = "Standard Error",
       main = "Funnel Plot with Study IDs")

## Blood_glucose
# Detecting outliers in funnel plot
funnel(results_wmd$Blood_glucose$model, 
       xlab = "Effect Size (WMD)", 
       ylab = "Standard Error",
       main = "Funnel Plot with Study IDs")

text(x = results_wmd$Blood_glucose$model$yi, 
     y = sqrt(results_wmd$Blood_glucose$model$vi), 
     labels = paste0("(", 1:nrow(results_wmd$Blood_glucose$data), ")"), 
     pos = 3,  
     cex = 0.7, 
     col = "blue")

# Removing outliers
studies_to_remove <- c(4)

results_wmd$Blood_glucose$data <- results_wmd$Blood_glucose$data[-studies_to_remove, ]

#Update the data
results_wmd$Blood_glucose$model <- rma(yi, vi, 
                             data = results_wmd$Blood_glucose$data, 
                             method = "DL")

results_wmd$Blood_glucose$n_studies <- nrow(results_wmd$Blood_glucose$data)

results_wmd$Blood_glucose$funnel <- regtest(results_wmd$Blood_glucose$model)

# Checking the new funnel plot
regtest(results_wmd$Blood_glucose$model)

funnel(results_wmd$Blood_glucose$model, 
       xlab = "Effect Size (WMD)", 
       ylab = "Standard Error",
       main = "Funnel Plot with Study IDs")

#New values
funnel_test(results_wmd)


#SMD variables

funnel_test(results_smd)
#Lying, Standing, Self-grooming, Feeding, Playing, Consuming milk

## Lying
# Detecting outliers in funnel plot
funnel(results_smd$Lying$model, 
       xlab = "Effect Size (SMD)", 
       ylab = "Standard Error",
       main = "Funnel Plot with Study IDs")

text(x = results_smd$Lying$model$yi, 
     y = sqrt(results_smd$Lying$model$vi), 
     labels = paste0("(", 1:nrow(results_smd$Lying$data), ")"), 
     pos = 3,
     cex = 0.7, 
     col = "blue")

# Removing outliers
studies_to_remove <- c(10)

results_smd$Lying$data <- results_smd$Lying$data[-studies_to_remove, ]

#Update the data
results_smd$Lying$model <- rma(yi, vi, 
                                       data = results_smd$Lying$data, 
                                       method = "REML")


results_smd$Lying$n_studies <- nrow(results_smd$Lying$data)

results_smd$Lying$funnel <- regtest(results_smd$Lying$model)

# Checking the new funnel plot
regtest(results_smd$Lying$model)

funnel(results_smd$Lying$model, 
       xlab = "Effect Size (smd)", 
       ylab = "Standard Error",
       main = "Funnel Plot with Study IDs")

## Standing
# Detecting outliers in funnel plot
funnel(results_smd$Standing$model, 
       xlab = "Effect Size (SMD)", 
       ylab = "Standard Error",
       main = "Funnel Plot with Study IDs")

text(x = results_smd$Standing$model$yi, 
     y = sqrt(results_smd$Standing$model$vi), 
     labels = paste0("(", 1:nrow(results_smd$Standing$data), ")"), 
     pos = 3,
     cex = 0.7, 
     col = "blue")

# Removing outliers
studies_to_remove <- c(10, 17)

results_smd$Standing$data <- results_smd$Standing$data[-studies_to_remove, ]

#Update the data
results_smd$Standing$model <- rma(yi, vi, 
                               data = results_smd$Standing$data, 
                               method = "REML")

results_smd$Standing$n_studies <- nrow(results_smd$Standing$data)

results_smd$Standing$funnel <- regtest(results_smd$Standing$model)

# Checking the new funnel plot
regtest(results_smd$Standing$model)

funnel(results_smd$Standing$model, 
       xlab = "Effect Size (smd)", 
       ylab = "Standard Error",
       main = "Funnel Plot with Study IDs")

## Self-grooming
# Detecting outliers in funnel plot
funnel(results_smd$'Self-grooming'$model, 
       xlab = "Effect Size (SMD)", 
       ylab = "Standard Error",
       main = "Funnel Plot with Study IDs")

text(x = results_smd$'Self-grooming'$model$yi, 
     y = sqrt(results_smd$'Self-grooming'$model$vi), 
     labels = paste0("(", 1:nrow(results_smd$'Self-grooming'$data), ")"), 
     pos = 3,
     cex = 0.7, 
     col = "blue")

# Removing outliers
studies_to_remove <- c(5)

results_smd$'Self-grooming'$data <- results_smd$'Self-grooming'$data[-studies_to_remove, ]

#Update the data
results_smd$'Self-grooming'$model <- rma(yi, vi, 
                                  data = results_smd$'Self-grooming'$data, 
                                  method = "REML")

results_smd$'Self-grooming'$n_studies <- nrow(results_smd$'Self-grooming'$data)

results_smd$'Self-grooming'$funnel <- regtest(results_smd$'Self-grooming'$model)

# Checking the new funnel plot
regtest(results_smd$'Self-grooming'$model)

funnel(results_smd$'Self-grooming'$model, 
       xlab = "Effect Size (smd)", 
       ylab = "Standard Error",
       main = "Funnel Plot with Study IDs")

## Feeding
# Detecting outliers in funnel plot
funnel(results_smd$Feeding$model, 
       xlab = "Effect Size (SMD)", 
       ylab = "Standard Error",
       main = "Funnel Plot with Study IDs")

text(x = results_smd$Feeding$model$yi, 
     y = sqrt(results_smd$Feeding$model$vi), 
     labels = paste0("(", 1:nrow(results_smd$Feeding$data), ")"), 
     pos = 3,
     cex = 0.7, 
     col = "blue")

# Removing outliers
studies_to_remove <- c(4)

results_smd$Feeding$data <- results_smd$Feeding$data[-studies_to_remove, ]

#Update the data
results_smd$Feeding$model <- rma(yi, vi, 
                                         data = results_smd$Feeding$data, 
                                         method = "REML")

results_smd$Feeding$n_studies <- nrow(results_smd$Feeding$data)

results_smd$Feeding$funnel <- regtest(results_smd$Feeding$model)

# Checking the new funnel plot
regtest(results_smd$Feeding$model)

funnel(results_smd$Feeding$model, 
       xlab = "Effect Size (smd)", 
       ylab = "Standard Error",
       main = "Funnel Plot with Study IDs")

## Playing
# Detecting outliers in funnel plot
funnel(results_smd$Playing$model, 
       xlab = "Effect Size (SMD)", 
       ylab = "Standard Error",
       main = "Funnel Plot with Study IDs")

text(x = results_smd$Playing$model$yi, 
     y = sqrt(results_smd$Playing$model$vi), 
     labels = paste0("(", 1:nrow(results_smd$Playing$data), ")"), 
     pos = 3,
     cex = 0.7, 
     col = "blue")

# Removing outliers
studies_to_remove <- c(7)

results_smd$Playing$data <- results_smd$Playing$data[-studies_to_remove, ]

#Update the data
results_smd$Playing$model <- rma(yi, vi, 
                                 data = results_smd$Playing$data, 
                                 method = "REML")

results_smd$Playing$n_studies <- nrow(results_smd$Playing$data)

results_smd$Playing$funnel <- regtest(results_smd$Playing$model)

# Checking the new funnel plot
regtest(results_smd$Playing$model)

funnel(results_smd$Playing$model, 
       xlab = "Effect Size (smd)", 
       ylab = "Standard Error",
       main = "Funnel Plot with Study IDs")

## Consuming_milk
# Detecting outliers in funnel plot
funnel(results_smd$Consuming_milk$model, 
       xlab = "Effect Size (SMD)", 
       ylab = "Standard Error",
       main = "Funnel Plot with Study IDs")

text(x = results_smd$Consuming_milk$model$yi, 
     y = sqrt(results_smd$Consuming_milk$model$vi), 
     labels = paste0("(", 1:nrow(results_smd$Consuming_milk$data), ")"), 
     pos = 3,
     cex = 0.7, 
     col = "blue")

# Removing outliers
studies_to_remove <- c(3)

results_smd$Consuming_milk$data <- results_smd$Consuming_milk$data[-studies_to_remove, ]

#Update the data
results_smd$Consuming_milk$model <- rma(yi, vi, 
                                 data = results_smd$Consuming_milk$data, 
                                 method = "REML")

results_smd$Consuming_milk$n_studies <- nrow(results_smd$Consuming_milk$data)

results_smd$Consuming_milk$funnel <- regtest(results_smd$Consuming_milk$model)

# Checking the new funnel plot
regtest(results_smd$Consuming_milk$model)

funnel(results_smd$Consuming_milk$model, 
       xlab = "Effect Size (smd)", 
       ylab = "Standard Error",
       main = "Funnel Plot with Study IDs")

# ----------------------------
# META-ANALYSIS TABLES
# ----------------------------
# Function to create summary tables
create_table <- function(results, smd = FALSE) {
  map_dfr(results, function(res) {
    data.frame(
      Control_Mean = mean(res$data$Control_Mean),
      N_Studies = res$n_studies,
      Effect_Size = ifelse(smd, "SMD", "WMD"),
      Estimate = sprintf("%.2f (%.2f to %.2f)", 
                         res$model$b, res$model$ci.lb, res$model$ci.ub),
      P_Value = ifelse(res$model$pval < 0.001, "<0.001", 
                       sprintf("%.3f", res$model$pval)),
      I2 = sprintf("%.1f%%", res$model$I2),
      I2_P = ifelse(res$model$QEp < 0.001, "<0.001", 
                    sprintf("%.3f", res$model$QEp)),
      Funnel_P = ifelse(res$funnel$pval < 0.001, "<0.001", 
                        sprintf("%.3f", res$funnel$pval))
    )
  }, .id = "Outcome")
}

table_wmd <- create_table(results_wmd)
View(table_wmd)


table_smd <- create_table(results_smd, smd = TRUE) %>% select(-Control_Mean)
View(table_smd)

# ----------------------------
# META-REGRESSION
# ----------------------------

run_meta_reg <- function(results_list, var_name) {
 
  tryCatch({
    mod <- rma(yi, vi, 
               mods = ~ Weaning_age_class + Sex + Milk_allowance_class +
                 Housing_type + Total_area_per_calf + Group_size_class,
               data = results_list[[var_name]]$data,
               method = "DL")
    list(
      variable = var_name,
      n_studies = nrow(results_list[[var_name]]$data),
      adj_r2 = mod$R2,
      model = mod,
      data = results_list[[var_name]]$data
    )
  }, error = function(e) {
    message("Error in meta-regression for ", var_name, ": ", e$message)
    NULL
  })
}

# Criterias meet
# 1 - p<0.05 for heterogeneity = WW, ADG, CI, Lying, Feeding, Playing
# 2 - p>0.05 for funnel plot = WW, ADG, CI, Lying, Feeding, Playing
# 3 - no observations with values for studentized residuals 
# out of the range −2.5 to 2.5 (outliers) = ADG, CI, Lying, Playing
# 4 - high heterogeneity (I² statistic >50%) = ADG, CI, Lying, Playing

# ADG, CI, Lying, Playing
meta_reg_vars <- c("ADG", "Lying", "Playing")
mod_base_names <- c("Weaning_age_class", "Sex", "Milk_allowance_class",
                    "Housing_type", "Total_area_per_calf", "Group_size_class")

meta_reg_results <- lapply(meta_reg_vars, function(x) {
  if(x %in% vars_wmd) {
    run_meta_reg(results_wmd, x)
  } else {
    run_meta_reg(results_smd, x)
  }
})
names(meta_reg_results) <- meta_reg_vars


for(var in meta_reg_vars) {
  if(!is.null(meta_reg_results[[var]])) {
    cat("\n=== Detailed results for", var, "===\n")
    print(summary(meta_reg_results[[var]]$model))
  }
}

extract_meta_reg <- function(meta_reg_results) {
  mod_base_names <- c("Weaning_age_class", "Sex", "Milk_allowance_class",
                      "Housing_type", "Total_area_per_calf", "Group_size_class")
  
  map_dfr(meta_reg_results, function(res) {
    if(is.null(res)) return(data.frame())
    
    coefs <- coef(res$model)[-1]
    pvals <- res$model$pval[-1]
    mod_groups <- names(coefs)
    
    mod_summary <- map_dfr(unique(mod_groups), function(mod) {
        idx <- which(mod_groups == mod)
      
      most_sig_idx <- idx[which.min(pvals[idx])]
      
      data.frame(
        Moderator = mod,
        Estimate = sprintf("%.2f", coefs[most_sig_idx]),
        P_value = ifelse(pvals[most_sig_idx] < 0.001, "<0.001", 
                         sprintf("%.3f", pvals[most_sig_idx])),
        stringsAsFactors = FALSE
      )
    })
    
    data.frame(
      Variable = res$variable,
      N_Studies = res$n_studies,
      Adj_R2 = sprintf("%.1f%%", res$adj_r2),
      mod_summary,
      stringsAsFactors = FALSE
    )
  })
}

meta_reg_table <- extract_meta_reg(meta_reg_results)
View(meta_reg_table)

list_original_moderators <- function(meta_reg_results) {
  unique_moderators <- map_df(meta_reg_results, function(res) {
    if (!is.null(res)) {
      tibble(
        Variable = res$variable,
        Moderator = names(coef(res$model))[-1]
      )}}) %>% 
    distinct(Moderator) %>%
    arrange(Moderator)
  
  return(unique_moderators)
}
original_names <- list_original_moderators(meta_reg_results)
print(original_names)

# Reference categories
reference_info <- tribble(
  ~Factor,            ~Reference_Category,
  "Weaning Age",      "≥8 weeks (baseline)",
  "Sex",              "Female (baseline)",
  "Milk Allowance",   "High (baseline)", 
  "Housing Type",     "Indoor (baseline)",
  "Space Allowance",  "<2.0 m² (baseline)",
  "Group Size",       "Individual (baseline)"
)

# Factor-level mapping
factor_level_mapping <- tribble(
  ~original_name,                           ~Factor,           ~Level,
  "Weaning_age_classless than 8 wk",        "Weaning Age",     "<8 weeks",
  "Weaning_age_classmore than 8 wk",        "Weaning Age",     "≥8 weeks",
  "SexMale",                                "Sex",             "Male",
  "SexMale and Female",                     "Sex",             "Mixed",
  "Milk_allowance_classLow",                "Milk Allowance",  "Low",
  "Milk_allowance_classMedium",             "Milk Allowance",  "Medium",
  "Milk_allowance_classNA",                 "Milk Allowance",  "Not reported",
  "Housing_typeoutdoor",                    "Housing Type",    "Outdoor",
  "Total_area_per_calf> 4.0",               "Space Allowance", ">4.0 m²",
  "Total_area_per_calf2.0 to 4.0",          "Space Allowance", "2.0-4.0 m²",
  "Total_area_per_calfNA",                  "Space Allowance", "Not reported",
  "Group_size_classgroup of 7 to 15 calves","Group Size",      "7-15 calves",
  "Group_size_classpair",                   "Group Size",      "Pair"
)

meta_reg_table_clean <- meta_reg_table %>%
  left_join(factor_level_mapping, by = c("Moderator" = "original_name")) %>%
  mutate(
    Factor = ifelse(is.na(Factor), 
                    str_to_title(gsub("_", " ", gsub("class|type", "", Moderator))), 
                    Factor),
    Level = ifelse(is.na(Level), "", Level),
    P_value = case_when(
      P_value == "<0.001" ~ "<0.001",
      is.na(as.numeric(P_value)) ~ NA_character_,
      as.numeric(P_value) < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", as.numeric(P_value))
    )
  ) %>%
  left_join(reference_info, by = "Factor")

intercept_data <- meta_reg_results %>%
  map_dfr(function(res) {
    if (!is.null(res)) {
      tibble(
        Variable = res$variable,
        N_Studies = res$n_studies,
        Adj_R2 = sprintf("%.1f%%", res$adj_r2),
        Factor = "Intercept",
        Level = "Baseline",
        Estimate = sprintf("%.2f", coef(res$model)[1]),
        P_value = ifelse(res$model$pval[1] < 0.001, "<0.001",
                         sprintf("%.3f", res$model$pval[1])),
        Reference_Category = "All reference levels"
      )
    }
  })

complete_table <- bind_rows(intercept_data, meta_reg_table_clean) %>%
  arrange(Variable, Factor == "Intercept", Factor)

complete_table %>%
  flextable(col_keys = c("Variable", "N_Studies", "Adj_R2", "Factor", "Level", 
                         "Estimate", "P_value", "Reference_Category")) %>%
  set_header_labels(
    Variable = "Outcome",
    N_Studies = "Studies",
    Adj_R2 = "Adj. R²",
    Factor = "Factor",
    Level = "Level",
    Estimate = "Estimate",
    P_value = "P-value",
    Reference_Category = "Reference"
  ) %>%
  add_header_row(
    values = c("", "Model Summary", "Moderator Terms", "Effects"),
    colwidths = c(1, 2, 2, 3)  
  ) %>%
  theme_booktabs() %>%
  autofit() %>%
  align(align = "center", part = "all") %>%
  merge_v(j = c("Variable", "N_Studies", "Adj_R2")) %>%
  bg(i = ~ Factor == "Intercept", bg = "#f7f7f7") %>%
  bold(i = ~ Factor == "Intercept") %>%
  bold(~ P_value == "<0.001" | (!is.na(P_value) & as.numeric(P_value) < 0.05), 
       ~ P_value + Estimate) %>%
  footnote(
    i = 1, j = 6:8, part = "header",
    value = as_paragraph(
      "Intercept represents estimate at reference levels.",
      "Bold values indicate p < 0.05"
    )
  ) %>%
  padding(padding = 3, part = "all")

# ----------------------------
# FOREST PLOTS
# ----------------------------

create_final_forest <- function(data, var_name, moderator, smd = FALSE) {
  df <- data %>% 
    filter(Variable == var_name) %>%
    mutate(
      moderator_factor = factor(get(moderator)),
      Study_Label = paste(Authors, Year, sep = ", "),
      Treatment_SD = Treatment_SEM * sqrt(Treatment_N),
      Control_SD = Control_SEM * sqrt(Control_N)
    )
  
  measure <- ifelse(smd, "SMD", "MD")
  df <- escalc(measure = measure,
               m1i = Treatment_Mean, sd1i = Treatment_SD, n1i = Treatment_N,
               m2i = Control_Mean, sd2i = Control_SD, n2i = Control_N,
               data = df)
  
  df$weight <- weights(rma(yi, vi, data = df, method = "DL"))
  
  subgroups <- levels(df$moderator_factor)
  n_studies <- nrow(df)
  n_subgroups <- length(subgroups)
  
  rows_needed <- n_studies + (n_subgroups * 2) + 5  
  
  par(mar = c(4, 7, 2, 3), xpd = TRUE)
  plot(NA, xlim = c(-6, 4), ylim = c(0.5, rows_needed + 0.5),
       xlab = ifelse(smd, "Standardized Mean Difference", "Weighted Mean Difference"),
       ylab = "", yaxt = "n", bty = "n", xaxt = "n")
  
  axis(1, at = seq(-2, 2, by = 0.5), cex.axis = 0.9)
  
  abline(v = 0, lty = "dashed", col = "gray", lwd = 1.2)
  
  text(-5.5, rows_needed, "Study", pos = 4, font = 2, cex = 1.0)
  text(1.5, rows_needed, "Effect (95% CI)", font = 2, cex = 1.0)
  text(3.2, rows_needed, "% Weight", font = 2, cex = 1.0)
  
  current_row <- rows_needed - 1  
  
  study_spacing <- 1.0
  subgroup_spacing <- 1.5
  stats_spacing <- 0.8
  diamond_height <- 0.35
  
  subgroup_results <- list()
  for (subgroup in subgroups) {
    subgroup_studies <- df[df$moderator_factor == subgroup, ]
    n_in_subgroup <- nrow(subgroup_studies)
    
    text(-5.5, current_row, subgroup, pos = 4, font = 2, cex = 1.0)
    current_row <- current_row - study_spacing
    
    for (i in 1:n_in_subgroup) {
      study <- subgroup_studies[i, ]
      
      text(-5.5, current_row, study$Study_Label, pos = 4, cex = 0.9)
      
      points(study$yi, current_row, pch = 18, cex = 1.3, col = "#1f78b4")
      lines(c(study$yi - 1.96*sqrt(study$vi), 
              study$yi + 1.96*sqrt(study$vi)), 
            rep(current_row, 2), col = "#1f78b4", lwd = 1.8)
      
      est <- sprintf("%.2f", study$yi)
      ci <- sprintf("(%.2f, %.2f)", 
                    study$yi - 1.96*sqrt(study$vi),
                    study$yi + 1.96*sqrt(study$vi))
      text(1.5, current_row, paste(est, ci), cex = 0.9)
      
      text(3.2, current_row, sprintf("%.1f%%", study$weight), cex = 0.9)
      
      current_row <- current_row - study_spacing
    }
    
    subgroup_res <- rma(yi, vi, data = subgroup_studies, method = "DL")
    subgroup_results[[subgroup]] <- subgroup_res
    
    polygon(c(subgroup_res$ci.lb, subgroup_res$beta, subgroup_res$ci.ub, subgroup_res$beta),
            c(current_row - diamond_height, current_row, 
              current_row - diamond_height, current_row - (2*diamond_height)),
            col = "#33a02c", border = NA)
    
    subgroup_label <- sprintf("Subgroup, DL (I² = %.1f%%, p = %.3f)", 
                              subgroup_res$I2, subgroup_res$QEp)
    text(-5.5, current_row - (diamond_height*1.5), subgroup_label, font = 2, cex = 0.9, pos = 4)
    
    est <- sprintf("%.2f", subgroup_res$beta)
    ci <- sprintf("(%.2f, %.2f)", subgroup_res$ci.lb, subgroup_res$ci.ub)
    text(1.5, current_row - (diamond_height*1.5), paste(est, ci), font = 2, cex = 1.0)
    
    text(3.2, current_row - (diamond_height*1.5), 
         sprintf("%.1f%%", sum(subgroup_studies$weight)), font = 2, cex = 1.0)
    
    current_row <- current_row - subgroup_spacing
  }
  
  res_overall <- rma(yi, vi, data = df, method = "DL")
  
  if (n_subgroups > 1) {
    bg_test <- anova.rma(res_overall, btt = 2:n_subgroups)
    bg_p <- bg_test$pval
  } else {
    bg_p <- NA
  }
  
  polygon(c(res_overall$ci.lb, res_overall$beta, res_overall$ci.ub, res_overall$beta),
          c(current_row - diamond_height, current_row, 
            current_row - diamond_height, current_row - (2*diamond_height)),
          col = "#e31a1c", border = NA)
  
  abline(v = res_overall$beta, col = "#e31a1c", lty = "dotted", lwd = 2)
  
  text(-5.5, current_row - (diamond_height*1.5), "Overall Effect", pos = 4, font = 2, cex = 1.0)
  
  est <- sprintf("%.2f", res_overall$beta)
  ci <- sprintf("(%.2f, %.2f)", res_overall$ci.lb, res_overall$ci.ub)
  text(1.5, current_row - (diamond_height*1.5), paste(est, ci), font = 2, cex = 1.0)
  
  overall_dl_label <- sprintf("Overall, DL (I² = %.1f%%, p = %.3f)", 
                              res_overall$I2, res_overall$QEp)
  text(-5.5, current_row - (diamond_height*3), overall_dl_label, font = 2, cex = 0.9, pos = 4)
  
  text(3.2, current_row - (diamond_height*1.5), "100.0%", font = 2, cex = 1.0)
  
  if (n_subgroups > 1 && !is.na(bg_p)) {
    bg_text <- sprintf("Heterogeneity between groups: p = %.3f", bg_p)
    text(-5.5, current_row - (diamond_height*4.5), bg_text, font = 2, cex = 0.9, pos = 4, col = "#e31a1c")
  }
  
  segments(x0 = seq(-2, 2, by = 0.5), y0 = 0.4, y1 = 0.6, lwd = 1.2)
  text(x = seq(-2, 2, by = 0.5), y = 0.2, labels = seq(-2, 2, by = 0.5), cex = 0.8)
}

create_final_forest(full_data, "ADG", "Sex")
create_final_forest(full_data, "ADG", "Housing_type")
create_final_forest(full_data, "ADG", "Total_area_per_calf")

create_final_forest(full_data, "Lying", "Housing_type", smd = TRUE)
create_final_forest(full_data, "Lying", "Milk_allowance_class", smd = TRUE)
create_final_forest(full_data, "Lying", "Sex", smd = TRUE)

create_final_forest(full_data, "Playing", "Weaning_age_class", smd = TRUE)
create_final_forest(full_data, "Playing", "Sex", smd = TRUE)

# ----------------------------
# SUBGROUP BAR PLOTS
# ----------------------------

create_subgroup_barplot <- function(data, var_name, moderator, smd = FALSE) {
  df <- data %>% 
    filter(Variable == var_name) %>%
    mutate(
      moderator_factor = factor(get(moderator)),
      Study_Label = paste(Authors, Year, sep = ", "),
      Treatment_SD = Treatment_SEM * sqrt(Treatment_N),
      Control_SD = Control_SEM * sqrt(Control_N)
    )
  
  measure <- ifelse(smd, "SMD", "MD")
  df <- escalc(measure = measure,
               m1i = Treatment_Mean, sd1i = Treatment_SD, n1i = Treatment_N,
               m2i = Control_Mean, sd2i = Control_SD, n2i = Control_N,
               data = df)
  
  subgroups <- levels(df$moderator_factor)
  
  subgroup_effects <- list()
  for (subgroup in subgroups) {
    subgroup_data <- df[df$moderator_factor == subgroup, ]
    res <- rma(yi, vi, data = subgroup_data, method = "DL")
    subgroup_effects[[subgroup]] <- data.frame(
      subgroup = subgroup,
      estimate = res$beta[1],
      ci.lb = res$ci.lb,
      ci.ub = res$ci.ub,
      n_studies = nrow(subgroup_data),
      stringsAsFactors = FALSE
    )
  }
  
  plot_data <- do.call(rbind, subgroup_effects)
  
  overall <- rma(yi, vi, data = df, method = "DL")
  
  par(mar = c(5, 8, 4, 2))
  y_pos <- barplot_height <- length(subgroups):1
  xlim <- range(c(plot_data$ci.lb, plot_data$ci.ub, overall$ci.lb, overall$ci.ub))
  xlim <- c(min(xlim[1], -0.5), max(xlim[2], 0.5))
  
  plot(NA, xlim = xlim, ylim = c(0.5, length(subgroups) + 1.5),
       xlab = ifelse(smd, "Standardized Mean Difference", "Mean Difference"),
       ylab = "", yaxt = "n", main = paste("Subgroup Analysis:", var_name))
  
  abline(v = 0, lty = 2, col = "gray")
  
  for (i in 1:nrow(plot_data)) {
    rect(0, y_pos[i] - 0.3, plot_data$estimate[i], y_pos[i] + 0.3,
         col = ifelse(plot_data$estimate[i] > 0, "#1f78b4", "#e31a1c"), border = NA)
    
    arrows(plot_data$ci.lb[i], y_pos[i], plot_data$ci.ub[i], y_pos[i],
           length = 0.05, angle = 90, code = 3, lwd = 1.5)
    
    axis(2, at = y_pos[i], labels = paste0(plot_data$subgroup[i], " (n=", plot_data$n_studies[i], ")"), 
         las = 1, tick = FALSE)
  }
  
  polygon(c(overall$ci.lb, overall$beta, overall$ci.ub, overall$beta),
          c(max(y_pos) + 1 - 0.15, max(y_pos) + 1, max(y_pos) + 1 - 0.15, max(y_pos) + 1 - 0.3),
          col = "#33a02c", border = NA)
  
  text(xlim[1], max(y_pos) + 1, "Overall Effect", pos = 4, font = 2)
  
  het_text <- sprintf("I² = %.1f%%, p = %.3f", overall$I2, overall$QEp)
  text(mean(xlim), max(y_pos) + 1.5, het_text, font = 2)
  
  if (length(subgroups) > 1) {
    bg_test <- anova.rma(overall, btt = 2:length(subgroups))
    bg_text <- sprintf("Between-group p = %.3f", bg_test$pval)
    text(mean(xlim), max(y_pos) + 1.3, bg_text, font = 2, col = "#e31a1c")
  }
}

create_subgroup_barplot(full_data, "ADG", "Sex")
create_subgroup_barplot(full_data, "ADG", "Housing_type")
create_subgroup_barplot(full_data, "ADG", "Total_area_per_calf")

create_subgroup_barplot(full_data, "Lying", "Housing_type", smd = TRUE)
create_subgroup_barplot(full_data, "Lying", "Milk_allowance_class", smd = TRUE)
create_subgroup_barplot(full_data, "Lying", "Sex", smd = TRUE)

create_subgroup_barplot(full_data, "Playing", "Weaning_age_class", smd = TRUE)
create_subgroup_barplot(full_data, "Playing", "Sex", smd = TRUE)

