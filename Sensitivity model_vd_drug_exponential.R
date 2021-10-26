
options(digits = 3, scipen = 999)

# Crude
drug_crude <- survreg(Surv(fu_yr, hz) ~vd_prescription, 
                        data = bd_i, dist = "exponential") %>% 
      tidy(conf.int = T, exponentiate = T) %>% 
      filter(term == "vd_prescriptionHad vitamin D prescription")
drug_crude[,-1] <- map(drug_crude[,-1], exp)
drug_crude$model <- "crude"

# Partially adjusted model
drug_partial <- survreg(Surv(fu_yr, hz) ~vd_prescription + sex + age_c, 
                     data = bd_i, dist = "exponential") %>% 
      tidy(conf.int = T, exponentiate = T) %>% 
      filter(term == "vd_prescriptionHad vitamin D prescription")
drug_partial[,-1] <- map(drug_partial[,-1], exp)

drug_partial$model <- "Partially adjusted"

# full model:
drug_full <- survreg(Surv(fu_yr, hz) ~vd_prescription + sex + age_c + ethnic + 
                           bmi_group + drink_freq_c + smoke_stat + imd_bd_q + regions +
                           season_c + asthma + ckd+ copd + depress + dm + ibd + ra + sle + immunosuppression, 
                     data = bd_i, dist = "exponential") %>% 
      tidy(conf.int = T, exponentiate = T) %>% 
      filter(term == "vd_prescriptionHad vitamin D prescription")
drug_full[,-1] <- map(drug_full[,-1], exp)
drug_full$model <- "fully adjusted"


exp_result <- bind_rows(drug_crude,drug_partial,drug_full)


exp_result$statistic <- NULL
exp_result$std.error <- NULL


exp_result %>% write_csv("exponential_reg_result.csv")
