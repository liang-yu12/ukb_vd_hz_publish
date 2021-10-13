# Model check : make sure the model fits the dataset 
# 1. check the Deviance residuals is close to 0
# 2. check goodness of fit using chi-square. H0: model works well; HA: model doesn't fit#    
library(survival)
library(lmtest)


## 1: Primary exposure vitamin D and hz ----

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


## Model check using Weibull regression ----
# primary
hz_vd.weibull <- survreg( Surv(fu_yr, hz) ~ vitd_s + sex + age_c + ethnic + 
                                bmi_group + drink_freq_c + smoke_stat + imd_bd_q + regions +
                                season_c + asthma + ckd+ copd + depress + dm + ibd + ra + sle + 
                                immunosuppression, 
                          data = bd_i, dist = "weibull")
hz_vd.weibull %>% summary()
# # Call:
# survreg(formula = Surv(fu_yr, hz) ~ vitd_s + sex + age_c + ethnic + 
#               bmi_group + drink_freq_c + smoke_stat + imd_bd_q + regions + 
#               season_c + asthma + ckd + copd + depress + dm + ibd + ra + 
#               sle + immunosuppression, data = bd_i, dist = "weibull")
# Value Std. Error      z                    p
# (Intercept)                         8.61630    0.17029  50.60 < 0.0000000000000002
# vitd_s0_deficiency                  0.00613    0.05207   0.12              0.90627
# vitd_s1_insufficiency              -0.02747    0.03397  -0.81              0.41880
# sexMale                             0.27155    0.03257   8.34 < 0.0000000000000002
# age_c                              -0.03864    0.00214 -18.04 < 0.0000000000000002
# ethnic2_mixed                      -0.07715    0.22674  -0.34              0.73366
# ethnic3_asian                       0.06467    0.13550   0.48              0.63320
# ethnic4_black                       0.16980    0.18426   0.92              0.35678
# ethnic5_chinese                    -0.10422    0.33813  -0.31              0.75790
# ethnic6_others                      0.15419    0.17896   0.86              0.38892
# bmi_group0_underweight             -0.08530    0.22011  -0.39              0.69838
# bmi_group2_over weight             -0.02448    0.03692  -0.66              0.50726
# bmi_group3_obese                    0.02711    0.04262   0.64              0.52479
# drink_freq_csometimes              -0.00823    0.06154  -0.13              0.89364
# drink_freq_cweekly                 -0.03673    0.05855  -0.63              0.53040
# drink_freq_cdaily                   0.04390    0.06540   0.67              0.50205
# smoke_statex-smoker                -0.13907    0.03340  -4.16    0.000031251318833
# smoke_statcurrent-smoker           -0.01874    0.05506  -0.34              0.73363
# imd_bd_q2_deprived                  0.00377    0.05010   0.08              0.94009
# imd_bd_q3_deprived                 -0.08798    0.04889  -1.80              0.07194
# imd_bd_q4_deprived                 -0.06209    0.04972  -1.25              0.21174
# imd_bd_qmost_deprived              -0.01793    0.05189  -0.35              0.72972
# regionsLondon                       0.12010    0.07447   1.61              0.10682
# regionsNorth East                  -0.03881    0.06068  -0.64              0.52242
# regionsNorth West                  -0.20045    0.06860  -2.92              0.00348
# regionsSouth East                   0.00615    0.12094   0.05              0.95943
# regionsSouth West                  -0.09172    0.09745  -0.94              0.34660
# regionsWest Midlands               -0.00544    0.07921  -0.07              0.94525
# regionsYorkshire and The Humber     0.02055    0.05457   0.38              0.70653
# regionsWales                       -0.59373    0.06375  -9.31 < 0.0000000000000002
# regionsScotland                     0.31631    0.06643   4.76    0.000001922119857
# season_c1_Spring                    0.04343    0.04453   0.98              0.32940
# season_c3_Autumn                    0.05353    0.04488   1.19              0.23302
# season_c4_Winter                    0.01082    0.04773   0.23              0.82069
# asthmaAsthma                       -0.14131    0.04363  -3.24              0.00120
# ckdCKD                             -0.41553    0.12993  -3.20              0.00138
# copdCOPD                           -0.29275    0.07830  -3.74              0.00019
# depressdepression                  -0.08067    0.04705  -1.71              0.08642
# dmHave DM                          -0.03264    0.06181  -0.53              0.59748
# ibdInflammatory bowel disease      -0.30727    0.08129  -3.78              0.00016
# raRA                               -0.26472    0.10398  -2.55              0.01090
# sleSLE                             -0.32239    0.28083  -1.15              0.25097
# immunosuppressionImmunosuppression -0.54955    0.07130  -7.71    0.000000000000013
# Log(scale)                          0.19382    0.01238  15.66 < 0.0000000000000002
# 
# Scale= 1.21 
# 
# Weibull distribution
# Loglik(model)= -41442   Loglik(intercept only)= -41913
# Chisq= 941 on 42 degrees of freedom, p= 0.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000062 
# Number of Newton-Raphson Iterations: 10 
# n=171349 (6223 observations deleted due to missingness)

hz_vd.exp <- survreg( Surv(fu_yr, hz) ~ vitd_s + sex + age_c + ethnic + 
                            bmi_group + drink_freq_c + smoke_stat + imd_bd_q + regions +
                            season_c + asthma + ckd+ copd + depress + dm + ibd + ra + sle + 
                            immunosuppression, 
                      data = bd_i, dist = "exponential")
hz_vd.exp %>% summary

# Call:
#       survreg(formula = Surv(fu_yr, hz) ~ vitd_s + sex + age_c + ethnic + 
#                     bmi_group + drink_freq_c + smoke_stat + imd_bd_q + regions + 
#                     season_c + asthma + ckd + copd + depress + dm + ibd + ra + 
#                     sle + immunosuppression, data = bd_i, dist = "exponential")
# Value Std. Error      z                    p
# (Intercept)                         7.52274    0.12525  60.06 < 0.0000000000000002
# vitd_s0_deficiency                  0.00430    0.04290   0.10              0.92014
# vitd_s1_insufficiency              -0.02263    0.02799  -0.81              0.41885
# sexMale                             0.22253    0.02668   8.34 < 0.0000000000000002
# age_c                              -0.03211    0.00172 -18.62 < 0.0000000000000002
# ethnic2_mixed                      -0.06233    0.18679  -0.33              0.73860
# ethnic3_asian                       0.05639    0.11163   0.51              0.61343
# ethnic4_black                       0.14094    0.15178   0.93              0.35311
# ethnic5_chinese                    -0.08477    0.27855  -0.30              0.76088
# ethnic6_others                      0.12688    0.14742   0.86              0.38943
# bmi_group0_underweight             -0.07260    0.18133  -0.40              0.68888
# bmi_group2_over weight             -0.02000    0.03041  -0.66              0.51070
# bmi_group3_obese                    0.02241    0.03511   0.64              0.52336
# drink_freq_csometimes              -0.00635    0.05070  -0.13              0.90029
# drink_freq_cweekly                 -0.02977    0.04823  -0.62              0.53702
# drink_freq_cdaily                   0.03688    0.05388   0.68              0.49359
# smoke_statex-smoker                -0.11530    0.02748  -4.20   0.0000271066130494
# smoke_statcurrent-smoker           -0.01715    0.04535  -0.38              0.70525
# imd_bd_q2_deprived                  0.00342    0.04127   0.08              0.93398
# imd_bd_q3_deprived                 -0.07222    0.04027  -1.79              0.07292
# imd_bd_q4_deprived                 -0.05088    0.04096  -1.24              0.21413
# imd_bd_qmost_deprived              -0.01490    0.04275  -0.35              0.72748
# regionsLondon                       0.08769    0.06132   1.43              0.15272
# regionsNorth East                  -0.03204    0.04998  -0.64              0.52149
# regionsNorth West                  -0.15798    0.05646  -2.80              0.00514
# regionsSouth East                   0.01547    0.09963   0.16              0.87662
# regionsSouth West                  -0.07932    0.08028  -0.99              0.32309
# regionsWest Midlands               -0.01615    0.06525  -0.25              0.80457
# regionsYorkshire and The Humber     0.01431    0.04495   0.32              0.75020
# regionsWales                       -0.47916    0.05210  -9.20 < 0.0000000000000002
# regionsScotland                     0.27966    0.05468   5.11   0.0000003148206959
# season_c1_Spring                    0.03513    0.03668   0.96              0.33823
# season_c3_Autumn                    0.04549    0.03698   1.23              0.21858
# season_c4_Winter                    0.01064    0.03933   0.27              0.78683
# asthmaAsthma                       -0.11645    0.03592  -3.24              0.00119
# ckdCKD                             -0.34616    0.10696  -3.24              0.00121
# copdCOPD                           -0.24554    0.06444  -3.81              0.00014
# depressdepression                  -0.06737    0.03875  -1.74              0.08213
# dmHave DM                          -0.02884    0.05092  -0.57              0.57113
# ibdInflammatory bowel disease      -0.25498    0.06690  -3.81              0.00014
# raRA                               -0.21705    0.08563  -2.53              0.01125
# sleSLE                             -0.26310    0.23133  -1.14              0.25539
# immunosuppressionImmunosuppression -0.46016    0.05852  -7.86   0.0000000000000037
# 
# Scale fixed at 1 
# 
# Exponential distribution
# Loglik(model)= -41573   Loglik(intercept only)= -42049
# Chisq= 952 on 42 degrees of freedom, p= 0.00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000032 
# Number of Newton-Raphson Iterations: 7 
# n=171349 (6223 observations deleted due to missingness)
# examined with lrtest 

lrtest(hz_vd.exp,hz_vd.weibull)
# Likelihood ratio test
# 
# Model 1: Surv(fu_yr, hz) ~ vitd_s + sex + age_c + ethnic + bmi_group + 
#       drink_freq_c + smoke_stat + imd_bd_q + regions + season_c + 
#       asthma + ckd + copd + depress + dm + ibd + ra + sle + immunosuppression
# Model 2: Surv(fu_yr, hz) ~ vitd_s + sex + age_c + ethnic + bmi_group + 
#       drink_freq_c + smoke_stat + imd_bd_q + regions + season_c + 
#       asthma + ckd + copd + depress + dm + ibd + ra + sle + immunosuppression
# #Df LogLik Df Chisq          Pr(>Chisq)    
# 1  43 -41573                                 
# 2  44 -41442  1   262 <0.0000000000000002 ***
#       ---
#       Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# In exponential distribution, the baseline risk is constant over time, 
# which is similar to Poisson; but the lrtest doesn't support this assumption


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


# Model checking using Weibull regression -----

hz_vdsupp.w <- survreg(Surv(fu_yr, hz) ~all_supp + sex + age_c + ethnic + 
                             bmi_group + drink_freq_c + smoke_stat + imd_bd_q + regions +
                             season_c + asthma + ckd+ copd + depress + dm + ibd + ra + sle + immunosuppression, 
                       data = bd_i, dist = "weibull")
hz_vdsupp.w %>% summary
# Call:
#       survreg(formula = Surv(fu_yr, hz) ~ all_supp + sex + age_c + 
#                     ethnic + bmi_group + drink_freq_c + smoke_stat + imd_bd_q + 
#                     regions + season_c + asthma + ckd + copd + depress + dm + 
#                     ibd + ra + sle + immunosuppression, data = bd_i, dist = "weibull")
# Value Std. Error     z                    p
# (Intercept)                            8.02528    0.27765 28.90 < 0.0000000000000002
# all_supp1_vitD and mineral supplement  0.03585    0.07443  0.48              0.63008
# sexMale                                0.19165    0.05146  3.72              0.00020
# age_c                                 -0.02909    0.00360 -8.08  0.00000000000000067
# ethnic2_mixed                         -0.06171    0.34915 -0.18              0.85972
# ethnic3_asian                         -0.05252    0.19939 -0.26              0.79224
# ethnic4_black                          0.08509    0.26800  0.32              0.75086
# ethnic5_chinese                        0.22099    0.60287  0.37              0.71394
# ethnic6_others                        -0.13013    0.25363 -0.51              0.60790
# bmi_group0_underweight                 0.15048    0.38310  0.39              0.69446
# bmi_group2_over weight                 0.01389    0.05727  0.24              0.80837
# bmi_group3_obese                       0.01281    0.06572  0.19              0.84544
# drink_freq_csometimes                  0.03191    0.09684  0.33              0.74180
# drink_freq_cweekly                    -0.09416    0.09143 -1.03              0.30306
# drink_freq_cdaily                      0.02492    0.10228  0.24              0.80752
# smoke_statex-smoker                   -0.12836    0.05118 -2.51              0.01215
# smoke_statcurrent-smoker               0.00164    0.09619  0.02              0.98637
# imd_bd_q2_deprived                    -0.03761    0.07824 -0.48              0.63073
# imd_bd_q3_deprived                    -0.20825    0.07574 -2.75              0.00596
# imd_bd_q4_deprived                    -0.13814    0.07778 -1.78              0.07572
# imd_bd_qmost_deprived                 -0.06305    0.08283 -0.76              0.44656
# regionsLondon                          0.09062    0.11418  0.79              0.42737
# regionsNorth East                     -0.07481    0.09359 -0.80              0.42408
# regionsNorth West                     -0.25841    0.10199 -2.53              0.01129
# regionsSouth East                     -0.15703    0.18004 -0.87              0.38308
# regionsSouth West                     -0.01401    0.15327 -0.09              0.92719
# regionsWest Midlands                  -0.04016    0.12348 -0.33              0.74498
# regionsYorkshire and The Humber        0.10661    0.08538  1.25              0.21177
# regionsWales                          -0.63618    0.09876 -6.44  0.00000000011798275
# regionsScotland                        0.36456    0.10813  3.37              0.00075
# season_c1_Spring                       0.03862    0.06706  0.58              0.56471
# season_c3_Autumn                       0.11954    0.07034  1.70              0.08923
# season_c4_Winter                      -0.00194    0.07129 -0.03              0.97834
# asthmaAsthma                          -0.06934    0.06419 -1.08              0.27999
# ckdCKD                                -0.33936    0.19633 -1.73              0.08390
# copdCOPD                              -0.35083    0.11399 -3.08              0.00209
# depressdepression                      0.01002    0.07270  0.14              0.89034
# dmHave DM                              0.11479    0.09878  1.16              0.24521
# ibdInflammatory bowel disease         -0.46468    0.11296 -4.11  0.00003893566456552
# raRA                                  -0.23132    0.14391 -1.61              0.10796
# sleSLE                                -0.25233    0.36693 -0.69              0.49164
# immunosuppressionImmunosuppression    -0.54768    0.10450 -5.24  0.00000015980365335
# Log(scale)                             0.18378    0.01947  9.44 < 0.0000000000000002
# 
# Scale= 1.2 
# 
# Weibull distribution
# Loglik(model)= -16493   Loglik(intercept only)= -16654
# Chisq= 323 on 41 degrees of freedom, p= 0.0000000000000000000000000000000000000000000021 
# Number of Newton-Raphson Iterations: 9 
# n=62416 (115156 observations deleted due to missingness)

# Exponential distribution

hz_vdsupp.exp <- survreg(Surv(fu_yr, hz) ~all_supp + sex + age_c + ethnic + 
                               bmi_group + drink_freq_c + smoke_stat + imd_bd_q + regions +
                               season_c + asthma + ckd+ copd + depress + dm + ibd + ra + sle + immunosuppression, 
                         data = bd_i, dist = "exponential")
hz_vdsupp.exp %>% summary
# # Call:
# survreg(formula = Surv(fu_yr, hz) ~ all_supp + sex + age_c + 
#               ethnic + bmi_group + drink_freq_c + smoke_stat + imd_bd_q + 
#               regions + season_c + asthma + ckd + copd + depress + dm + 
#               ibd + ra + sle + immunosuppression, data = bd_i, dist = "exponential")
# Value Std. Error     z                    p
# (Intercept)                            7.082414   0.212437 33.34 < 0.0000000000000002
# all_supp1_vitD and mineral supplement  0.030180   0.061934  0.49              0.62605
# sexMale                                0.158264   0.042701  3.71              0.00021
# age_c                                 -0.024470   0.002965 -8.25 < 0.0000000000000002
# ethnic2_mixed                         -0.050632   0.290536 -0.17              0.86165
# ethnic3_asian                         -0.041204   0.165916 -0.25              0.80387
# ethnic4_black                          0.071643   0.222994  0.32              0.74800
# ethnic5_chinese                        0.186400   0.501644  0.37              0.71021
# ethnic6_others                        -0.109272   0.211044 -0.52              0.60462
# bmi_group0_underweight                 0.125841   0.318782  0.39              0.69302
# bmi_group2_over weight                 0.011889   0.047652  0.25              0.80298
# bmi_group3_obese                       0.010556   0.054690  0.19              0.84695
# drink_freq_csometimes                  0.027222   0.080581  0.34              0.73550
# drink_freq_cweekly                    -0.077903   0.076060 -1.02              0.30573
# drink_freq_cdaily                      0.021626   0.085105  0.25              0.79941
# smoke_statex-smoker                   -0.107287   0.042539 -2.52              0.01167
# smoke_statcurrent-smoker              -0.000104   0.080033  0.00              0.99897
# imd_bd_q2_deprived                    -0.031087   0.065102 -0.48              0.63300
# imd_bd_q3_deprived                    -0.173100   0.062934 -2.75              0.00595
# imd_bd_q4_deprived                    -0.114170   0.064682 -1.77              0.07755
# imd_bd_qmost_deprived                 -0.052474   0.068917 -0.76              0.44642
# regionsLondon                          0.064506   0.094977  0.68              0.49703
# regionsNorth East                     -0.061827   0.077864 -0.79              0.42717
# regionsNorth West                     -0.208525   0.084734 -2.46              0.01386
# regionsSouth East                     -0.120972   0.149774 -0.81              0.41927
# regionsSouth West                     -0.015077   0.127538 -0.12              0.90589
# regionsWest Midlands                  -0.043405   0.102753 -0.42              0.67272
# regionsYorkshire and The Humber        0.086268   0.071018  1.21              0.22447
# regionsWales                          -0.519969   0.081414 -6.39        0.00000000017
# regionsScotland                        0.321630   0.089880  3.58              0.00035
# season_c1_Spring                       0.031142   0.055796  0.56              0.57675
# season_c3_Autumn                       0.100845   0.058506  1.72              0.08477
# season_c4_Winter                      -0.000317   0.059337 -0.01              0.99573
# asthmaAsthma                          -0.057860   0.053401 -1.08              0.27858
# ckdCKD                                -0.286556   0.163282 -1.75              0.07926
# copdCOPD                              -0.295690   0.094697 -3.12              0.00179
# depressdepression                      0.007674   0.060496  0.13              0.89906
# dmHave DM                              0.094592   0.082175  1.15              0.24969
# ibdInflammatory bowel disease         -0.388563   0.093710 -4.15        0.00003376733
# raRA                                  -0.192185   0.119713 -1.61              0.10841
# sleSLE                                -0.208853   0.305304 -0.68              0.49392
# immunosuppressionImmunosuppression    -0.461713   0.086582 -5.33        0.00000009678
# 
# Scale fixed at 1 
# 
# Exponential distribution
# Loglik(model)= -16541   Loglik(intercept only)= -16703
# Chisq= 325 on 41 degrees of freedom, p= 0.00000000000000000000000000000000000000000000067 
# Number of Newton-Raphson Iterations: 7 
# n=62416 (115156 observations deleted due to missingness)

# lrtest
lrtest(hz_vdsupp.exp, hz_vdsupp.w)
# Likelihood ratio test
# 
# Model 1: Surv(fu_yr, hz) ~ all_supp + sex + age_c + ethnic + bmi_group + 
#       drink_freq_c + smoke_stat + imd_bd_q + regions + season_c + 
#       asthma + ckd + copd + depress + dm + ibd + ra + sle + immunosuppression
# Model 2: Surv(fu_yr, hz) ~ all_supp + sex + age_c + ethnic + bmi_group + 
#       drink_freq_c + smoke_stat + imd_bd_q + regions + season_c + 
#       asthma + ckd + copd + depress + dm + ibd + ra + sle + immunosuppression
# #Df LogLik Df Chisq          Pr(>Chisq)    
# 1  42 -16541                                 
# 2  43 -16493  1  94.9 <0.0000000000000002 ***
#       ---
#       Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


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

# Model testing using Weibull regression ----
vd_drug.w <- survreg(Surv(fu_yr, hz) ~vd_prescription + sex + age_c + ethnic + 
                           bmi_group + drink_freq_c + smoke_stat + imd_bd_q + regions +
                           season_c + asthma + ckd+ copd + depress + dm + ibd + ra + sle + immunosuppression, 
                     data = bd_i, dist = "weibull")
vd_drug.w %>% summary
# Call:
#       survreg(formula = Surv(fu_yr, hz) ~ vd_prescription + sex + age_c + 
#                     ethnic + bmi_group + drink_freq_c + smoke_stat + imd_bd_q + 
#                     regions + season_c + asthma + ckd + copd + depress + dm + 
#                     ibd + ra + sle + immunosuppression, data = bd_i, dist = "weibull")
# Value Std. Error      z                    p
# (Intercept)                                8.60318    0.16907  50.89 < 0.0000000000000002
# vd_prescriptionHad vitamin D prescription -0.10005    0.09675  -1.03              0.30110
# sexMale                                    0.26912    0.03266   8.24 < 0.0000000000000002
# age_c                                     -0.03847    0.00214 -18.00 < 0.0000000000000002
# ethnic2_mixed                             -0.07881    0.22664  -0.35              0.72803
# ethnic3_asian                              0.07541    0.13456   0.56              0.57520
# ethnic4_black                              0.16938    0.18401   0.92              0.35730
# ethnic5_chinese                           -0.10603    0.33800  -0.31              0.75375
# ethnic6_others                             0.15625    0.17880   0.87              0.38218
# bmi_group0_underweight                    -0.08143    0.22013  -0.37              0.71144
# bmi_group2_over weight                    -0.02646    0.03689  -0.72              0.47312
# bmi_group3_obese                           0.02360    0.04229   0.56              0.57679
# drink_freq_csometimes                     -0.00964    0.06153  -0.16              0.87556
# drink_freq_cweekly                        -0.03784    0.05845  -0.65              0.51732
# drink_freq_cdaily                          0.04288    0.06536   0.66              0.51184
# smoke_statex-smoker                       -0.13914    0.03340  -4.17    0.000030927565903
# smoke_statcurrent-smoker                  -0.01878    0.05491  -0.34              0.73232
# imd_bd_q2_deprived                         0.00323    0.05010   0.06              0.94854
# imd_bd_q3_deprived                        -0.08820    0.04889  -1.80              0.07125
# imd_bd_q4_deprived                        -0.06216    0.04972  -1.25              0.21115
# imd_bd_qmost_deprived                     -0.01820    0.05185  -0.35              0.72556
# regionsLondon                              0.12021    0.07447   1.61              0.10646
# regionsNorth East                         -0.03950    0.06067  -0.65              0.51502
# regionsNorth West                         -0.20206    0.06861  -2.95              0.00323
# regionsSouth East                          0.00358    0.12096   0.03              0.97640
# regionsSouth West                         -0.09077    0.09744  -0.93              0.35160
# regionsWest Midlands                      -0.00666    0.07921  -0.08              0.93302
# regionsYorkshire and The Humber            0.02135    0.05457   0.39              0.69558
# regionsWales                              -0.59501    0.06375  -9.33 < 0.0000000000000002
# regionsScotland                            0.31414    0.06628   4.74    0.000002138454921
# season_c1_Spring                           0.04053    0.04326   0.94              0.34880
# season_c3_Autumn                           0.05236    0.04478   1.17              0.24232
# season_c4_Winter                           0.00780    0.04619   0.17              0.86600
# asthmaAsthma                              -0.14082    0.04363  -3.23              0.00125
# ckdCKD                                    -0.41029    0.13000  -3.16              0.00160
# copdCOPD                                  -0.29176    0.07829  -3.73              0.00019
# depressdepression                         -0.08007    0.04704  -1.70              0.08873
# dmHave DM                                 -0.03272    0.06177  -0.53              0.59631
# ibdInflammatory bowel disease             -0.30446    0.08134  -3.74              0.00018
# raRA                                      -0.26246    0.10402  -2.52              0.01163
# sleSLE                                    -0.31309    0.28094  -1.11              0.26509
# immunosuppressionImmunosuppression        -0.53976    0.07198  -7.50    0.000000000000065
# Log(scale)                                 0.19380    0.01238  15.66 < 0.0000000000000002
# 
# Scale= 1.21 
# 
# Weibull distribution
# Loglik(model)= -41442   Loglik(intercept only)= -41913
# Chisq= 941 on 41 degrees of freedom, p= 0.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000012 
# Number of Newton-Raphson Iterations: 10 
# n=171349 (6223 observations deleted due to missingness)

# Exponential distribution
vd_drug.exp <- survreg(Surv(fu_yr, hz) ~vd_prescription + sex + age_c + ethnic + 
                             bmi_group + drink_freq_c + smoke_stat + imd_bd_q + regions +
                             season_c + asthma + ckd+ copd + depress + dm + ibd + ra + sle + immunosuppression, 
                       data = bd_i, dist = "exponential")

vd_drug.exp %>% summary

# lrtest
lrtest(vd_drug.exp,vd_drug.w)
# Likelihood ratio test
# 
# Model 1: Surv(fu_yr, hz) ~ vd_prescription + sex + age_c + ethnic + bmi_group + 
#       drink_freq_c + smoke_stat + imd_bd_q + regions + season_c + 
#       asthma + ckd + copd + depress + dm + ibd + ra + sle + immunosuppression
# Model 2: Surv(fu_yr, hz) ~ vd_prescription + sex + age_c + ethnic + bmi_group + 
#       drink_freq_c + smoke_stat + imd_bd_q + regions + season_c + 
#       asthma + ckd + copd + depress + dm + ibd + ra + sle + immunosuppression
# #Df LogLik Df Chisq          Pr(>Chisq)    
# 1  42 -41573                                 
# 2  43 -41442  1   262 <0.0000000000000002 ***
#       ---
#       Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1