# Model check : make sure the model fits the dataset 
# 1. check the Deviance residuals is close to 0
# 2. check goodness of fit using chi-square. H0: model works well; HA: model doesn't fit#    


## 1: Primary exposure vitamin D and hz ----
# Create a survival object for person years-----
hz_surv <- Surv(time = as.numeric(bd_i$time_in)/365.25, 
                time2 = as.numeric(bd_i$time_out)/365.25, 
                event = bd_i$hz)

## crude
vd_m1 <- glm(hz~ vitd_s + offset(log(fu_yr)), 
             data = bd_i, family = "poisson")
summary(vd_m1)
# Deviance Residuals: 
#       Min      1Q  Median      3Q     Max  
#     -0.281  -0.147  -0.066   0.000   5.610  
# Goodness of fit
pchisq(vd_m1$deviance, df=vd_m1$df.residual, lower.tail=FALSE)
# p = 1

## partially adjusted 
vd_m2<- glm(hz~ vitd_s + sex + age_c + offset(log(fu_yr)), 
            data = bd_i, family = "poisson")
summary(vd_m2)
# Deviance Residuals: 
#       Min      1Q  Median      3Q     Max  
#     -0.289  -0.141  -0.070   0.007   5.597

pchisq(vd_m2$deviance, df=vd_m2$df.residual, lower.tail=FALSE)
# p = 1

## Fully adjusted
vd_m3<- glm(hz~ vitd_s + offset(log(fu_yr)) + sex + age_c + ethnic + 
                  bmi_group + drink_freq_c + smoke_stat + imd_bd_q + regions +
                  season_c + asthma + ckd+ copd + depress + dm + ibd + ra + sle + 
                  immunosuppression, 
            data = bd_i, family = "poisson")

summary(vd_m3)
# Deviance Residuals: 
#       Min      1Q  Median      3Q     Max  
#     -0.436  -0.120  -0.068  -0.018   5.597

pchisq(vd_m3$deviance, df=vd_m3$df.residual, lower.tail=FALSE)
# p = 1


## 2nd exposure: self-reported supplementation ----
supp_m1 <- glm(hz~all_supp + offset(log(fu_yr)), 
               data = bd_i, family = "poisson") 

summary(supp_m1)
# Deviance Residuals: 
#       Min      1Q  Median      3Q     Max  
#     -0.257  -0.152  -0.073  -0.006   5.087  

pchisq(supp_m1$deviance, df = supp_m1$df.residual, lower.tail = F)
# p = 1

# Partially adjusted
supp_m2 <- glm(hz~all_supp + offset(log(fu_yr)) + sex + age_c, 
               data = bd_i, family = "poisson")

summary(supp_m2)
# Deviance Residuals: 
#       Min      1Q  Median      3Q     Max  
#     -0.283  -0.147  -0.074  -0.001   5.076 

pchisq(supp_m2$deviance, df = supp_m2$df.residual, lower.tail = F)
# p = 1

# Fully adjusted
supp_m3 <- vd_supp_full <- glm(hz~ all_supp + offset(log(fu_yr)) + sex + age_c + ethnic + 
                                     bmi_group + drink_freq_c + smoke_stat + imd_bd_q + regions +
                                     season_c + asthma + ckd+ copd + depress + dm + ibd + ra + sle + immunosuppression, 
                               data = bd_i, family = "poisson")

summary(supp_m3)
# Deviance Residuals: 
#       Min      1Q  Median      3Q     Max  
#     -0.417  -0.126  -0.073  -0.022   5.078 

pchisq(supp_m3$deviance, df=supp_m3$df.residual, lower.tail = F)
# p = 1


## 2nd exposure: prescribed supplementation ----
# crude
drug_m1 <- glm(hz~vd_prescription + offset(log(fu_yr)), 
               data = bd_i, family = "poisson") 

summary(drug_m1)
# Deviance Residuals: 
#       Min      1Q  Median      3Q     Max  
#     -0.297  -0.147  -0.069   0.000   5.610 

pchisq(drug_m1$deviance, df = drug_m1$df.residual, lower.tail = F)
# p = 1

# Partially adjusted 
drug_m2 <- glm(hz~vd_prescription + offset(log(fu_yr)) + sex + age_c, 
               data = bd_i, family = "poisson")

summary(drug_m2)
# Deviance Residuals: 
#       Min      1Q  Median      3Q     Max  
#     -0.328  -0.141  -0.070   0.007   5.597  

pchisq(drug_m2$deviance, df=drug_m2$df.residual, lower.tail = F)
# p = 1

# Fully adjusted
drug_m3<- glm(hz~ vd_prescription + offset(log(fu_yr)) + sex + age_c + ethnic + 
                    bmi_group + drink_freq_c + smoke_stat + imd_bd_q + regions +
                    season_c + asthma + ckd+ copd + depress + dm + ibd + ra + sle + immunosuppression, 
              data = bd_i, family = "poisson")

summary(drug_m3)
# Deviance Residuals: 
#       Min      1Q  Median      3Q     Max  
#     -0.434  -0.120  -0.068  -0.018   5.597 

pchisq(drug_m3$deviance, df=drug_m3$df.residual, lower.tail = F)
# p = 1