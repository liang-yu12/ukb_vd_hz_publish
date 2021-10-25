# Script for visualizing the regression model
# Supplementary Fig 8: vitamin D supplement/prescription and herpes zoster using Poisson regression

# data:
source("./master_script.R")

# Load additional library: 
lapply(c("ggpubr","grid","gridExtra","forcats"), require, character.only=T)


# survival objects for 2nd exposure

# hz_surv <- survival::Surv(time = as.numeric(bd_i$time_in)/365.25, 
#                           time2 = as.numeric(bd_i$time_out)/365.25, 
#                           event = bd_i$hz)


# Vitamin D supplementation ----
# crude----
vd_supp_crude_glm <- glm(hz~all_supp + offset(log(fu_yr)), 
                         data = bd_i, family = "poisson") 

vd_supp_crude_glm %>% idr.display()  # check output

# Poisson regression predicting hz with offset = NULL 
# 
# IDR(95%CI)        P(Wald's test) P(LR-test)
# all_supp: 1_vitD and mineral supplement vs 0_no vitD supplement 0.91 (0.81,1.02)  0.121          0.116     
#                                                                                                            
# Log-likelihood = -14378.3598
# No. of observations = 64630
# AIC value = 28760.7197

vd_supp_crude <- vd_supp_crude_glm %>% 
      tidy(conf.int = T, conf.level = 0.95, exponentiate = T) %>% 
      filter(term != "(Intercept)") %>% 
      dplyr::select(term, estimate, conf.low, conf.high)

vd_supp_crude %<>% add_row(
      term = "No_supp",
      estimate = 1,
      conf.low=1,
      conf.high=1) # add reference

vd_supp_crude$model <- "Crude"

# add n
vd_supp_crude_n <- vd_supp_crude_glm %>% model.matrix() %>% 
      colSums()%>%                           # show the total number in each level
      as.data.frame() %>% t %>% as.data.frame() # transpose it make it easier to select


vd_supp_crude %<>% mutate(N=case_when(
      term == "No_supp" ~ (vd_supp_crude_n$`(Intercept)` - vd_supp_crude_n$`all_supp1_vitD and mineral supplement`),
      term == "all_supp1_vitD and mineral supplement" ~ vd_supp_crude_n$`all_supp1_vitD and mineral supplement`
))

rm(vd_supp_crude_glm, vd_supp_crude_n)# house keeping

# partial ----
vd_supp_partial_glm <- glm(hz~all_supp + offset(log(fu_yr)) + sex + age_c, 
                           data = bd_i, family = "poisson") 

vd_supp_partial_glm %>% idr.display()

# Poisson regression predicting hz with offset = NULL 
# crude IDR(95%CI)  adj. IDR(95%CI)   P(Wald's test) P(LR-test)
# all_supp: 1_vitD and mineral supplement vs 0_no vitD supplement 0.91 (0.81,1.02)  0.97 (0.86,1.09)  0.625          0.623     
#                                                                                                                              
# sex                                                             0.89 (0.82,0.96)  0.87 (0.8,0.94)   < 0.001        < 0.001   
#                                                                                                                              
# age_c (cont. var.)                                              1.03 (1.02,1.03)  1.03 (1.02,1.03)  < 0.001        < 0.001   
#                                                                                                                              
# Log-likelihood = -14331.1419
# No. of observations = 64630
# AIC value = 28670.2837

vd_supp_partial <- vd_supp_partial_glm %>% 
      tidy(conf.int = T, conf.level = 0.95, exponentiate = T) %>% 
      filter(term == "all_supp1_vitD and mineral supplement") %>% 
      dplyr::select(term, estimate, conf.low, conf.high)

vd_supp_partial %<>% add_row(
      term = "No_supp",
      estimate = 1,
      conf.low=1,
      conf.high=1) # add reference

vd_supp_partial$model <- "Partially adjusted" # marked the model

# add N
vd_supp_partial_n <- vd_supp_partial_glm %>% model.matrix() %>% 
      colSums()%>%                           # show the total number in each level
      as.data.frame() %>% t %>% as.data.frame() # transpose it make it easier to select

vd_supp_partial %<>% mutate(N=case_when(
      term == "No_supp" ~ (vd_supp_partial_n$`(Intercept)` - vd_supp_partial_n$`all_supp1_vitD and mineral supplement`),
      term == "all_supp1_vitD and mineral supplement" ~ vd_supp_partial_n$`all_supp1_vitD and mineral supplement`
))

rm(vd_supp_partial_glm, vd_supp_partial_n) # house keeping


# Full adjusted model ----
vd_supp_full_glm <- glm(hz~ all_supp + offset(log(fu_yr)) + sex + age_c + ethnic + 
                              bmi_group + drink_freq_c + smoke_stat + imd_bd_q + regions +
                              season_c + asthma + ckd+ copd + depress + dm + ibd + ra + sle + immunosuppression, 
                        data = bd_i, family = "poisson") 

vd_supp_full_glm %>% idr.display()

# Poisson regression predicting hz with offset = NULL 
# 
# crude IDR(95%CI)        adj. IDR(95%CI)         P(Wald's test)
# all_supp: 1_vitD and mineral supplement vs 0_no vitD supplement   0.9 (0.8,1.01)          0.97 (0.86,1.1)         0.626         
#                                                                                                                                 
# sex                                                               0.89 (0.82,0.96)        0.85 (0.79,0.93)        < 0.001       
#                                                                                                                                 
# age_c (cont. var.)                                                1.03 (1.02,1.03)        1.02 (1.02,1.03)        < 0.001       
#                                                                                                                                 
# ethnic: ref.=1_white                                                                                                            
#    2_mixed                                                        0.95 (0.54,1.68)        1.05 (0.6,1.86)         0.862         
#    3_asian                                                        0.9 (0.66,1.23)         1.04 (0.75,1.44)        0.804         
#    4_black                                                        0.78 (0.5,1.19)         0.93 (0.6,1.44)         0.748         
#    5_chinese                                                      0.67 (0.25,1.78)        0.83 (0.31,2.22)        0.71          
#    6_others                                                       1 (0.67,1.51)           1.12 (0.74,1.69)        0.605         
#                                                                                                                                 
# bmi_group: ref.=1_healthy weight                                                                                                
#    0_underweight                                                  0.87 (0.46,1.62)        0.88 (0.47,1.65)        0.693         
#    2_over weight                                                  1 (0.92,1.1)            0.99 (0.9,1.08)         0.803         
#    3_obese                                                        1.02 (0.92,1.13)        0.99 (0.89,1.1)         0.847         
#                                                                                                                                 
# drink_freq_c: ref.=never                                                                                                        
#    sometimes                                                      0.95 (0.81,1.11)        0.97 (0.83,1.14)        0.735         
#    weekly                                                         1.02 (0.88,1.18)        1.08 (0.93,1.25)        0.306         
#    daily                                                          0.95 (0.81,1.12)        0.98 (0.83,1.16)        0.799         
#                                                                                                                                 
# smoke_stat: ref.=non-smoker                                                                                                     
#    ex-smoker                                                      1.15 (1.06,1.25)        1.11 (1.02,1.21)        0.012         
#    current-smoker                                                 0.9651 (0.8276,1.1253)  1.0001 (0.8549,1.17)    0.999         
#                                                                                                                                 
# imd_bd_q: ref.=least_deprived                                                                                                   
#    2_deprived                                                     1.03 (0.91,1.17)        1.03 (0.91,1.17)        0.633         
#    3_deprived                                                     1.17 (1.04,1.33)        1.19 (1.05,1.35)        0.006         
#    4_deprived                                                     1.1 (0.98,1.25)         1.12 (0.99,1.27)        0.078         
#    most_deprived                                                  1.03 (0.91,1.17)        1.05 (0.92,1.21)        0.446         
#                                                                                                                                 
# regions: ref.=East Midlands                                                                                                     
#    London                                                         0.93 (0.78,1.12)        0.94 (0.78,1.13)        0.497         
#    North East                                                     1.06 (0.91,1.23)        1.06 (0.91,1.24)        0.427         
#    North West                                                     1.24 (1.05,1.46)        1.23 (1.04,1.45)        0.014         
#    South East                                                     1.07 (0.8,1.43)         1.13 (0.84,1.51)        0.419         
#    South West                                                     0.99 (0.77,1.27)        1.02 (0.79,1.3)         0.906         
#    West Midlands                                                  1.06 (0.87,1.29)        1.04 (0.85,1.28)        0.673         
#    Yorkshire and The Humber                                       0.93 (0.81,1.07)        0.92 (0.8,1.05)         0.224         
#    Wales                                                          1.68 (1.43,1.96)        1.68 (1.43,1.97)        < 0.001       
#    Scotland                                                       0.71 (0.6,0.85)         0.72 (0.61,0.86)        < 0.001       
#                                                                                                                                 
# season_c: ref.=2_Summer                                                                                                         
#    1_Spring                                                       1.01 (0.91,1.12)        0.97 (0.87,1.08)        0.577         
#    3_Autumn                                                       0.91 (0.81,1.01)        0.9 (0.81,1.01)         0.085         
#    4_Winter                                                       1.041 (0.931,1.1639)    1.0003 (0.8905,1.1237)  0.996         
#                                                                                                                                 
# asthma                                                            1.12 (1.01,1.24)        1.06 (0.95,1.18)        0.279         
#                                                                                                                                 
# ckd: CKD vs No CKD                                                1.63 (1.19,2.24)        1.33 (0.97,1.83)        0.079         
#                                                                                                                                 
# copd: COPD vs No COPD                                             1.57 (1.32,1.87)        1.34 (1.12,1.62)        0.002         
#                                                                                                                                 
# depress: depression vs No depression                              0.9861 (0.8774,1.1082)  0.9924 (0.8814,1.1173)  0.899         
#                                                                                                                                 
# dm: Have DM vs No DM                                              0.99 (0.85,1.16)        0.91 (0.77,1.07)        0.25          
#                                                                                                                                 
# ibd: Inflammatory bowel disease vs No inflammatory bowel diseaase 1.58 (1.32,1.89)        1.47 (1.23,1.77)        < 0.001       
#                                                                                                                                 
# ra: RA vs No RA                                                   1.61 (1.29,2.01)        1.21 (0.96,1.53)        0.108         
#                                                                                                                                 
# sle: SLE vs No SLE                                                1.54 (0.85,2.79)        1.23 (0.68,2.24)        0.494         
#                                                                                                                                 
# immunosuppression: Immunosuppression vs Not immunosuppressive     1.89 (1.61,2.21)        1.59 (1.34,1.88)        < 0.001       
#                                                                                                                                 
#                                                                   P(LR-test)
# all_supp: 1_vitD and mineral supplement vs 0_no vitD supplement   0.625     
#                                                                             
# sex                                                               < 0.001   
#                                                                             
# age_c (cont. var.)                                                < 0.001   
#                                                                             
# ethnic: ref.=1_white                                              0.987     
#    2_mixed                                                                  
#    3_asian                                                                  
#    4_black                                                                  
#    5_chinese                                                                
#    6_others                                                                 
#                                                                             
# bmi_group: ref.=1_healthy weight                                  0.976     
#    0_underweight                                                            
#    2_over weight                                                            
#    3_obese                                                                  
#                                                                             
# drink_freq_c: ref.=never                                          0.103     
#    sometimes                                                                
#    weekly                                                                   
#    daily                                                                    
#                                                                             
# smoke_stat: ref.=non-smoker                                       0.034     
#    ex-smoker                                                                
#    current-smoker                                                           
#                                                                             
# imd_bd_q: ref.=least_deprived                                     0.039     
#    2_deprived                                                               
#    3_deprived                                                               
#    4_deprived                                                               
#    most_deprived                                                            
#                                                                             
# regions: ref.=East Midlands                                       < 0.001   
#    London                                                                   
#    North East                                                               
#    North West                                                               
#    South East                                                               
#    South West                                                               
#    West Midlands                                                            
#    Yorkshire and The Humber                                                 
#    Wales                                                                    
#    Scotland                                                                 
#                                                                             
# season_c: ref.=2_Summer                                           0.265     
#    1_Spring                                                                 
#    3_Autumn                                                                 
#    4_Winter                                                                 
#                                                                             
# asthma                                                            0.281     
#                                                                             
# ckd: CKD vs No CKD                                                0.093     
#                                                                             
# copd: COPD vs No COPD                                             0.003     
#                                                                             
# depress: depression vs No depression                              0.899     
#                                                                             
# dm: Have DM vs No DM                                              0.244     
#                                                                             
# ibd: Inflammatory bowel disease vs No inflammatory bowel diseaase < 0.001   
#                                                                             
# ra: RA vs No RA                                                   0.117     
#                                                                             
# sle: SLE vs No SLE                                                0.508     
#                                                                             
# immunosuppression: Immunosuppression vs Not immunosuppressive     < 0.001   
#                                                                             
# Log-likelihood = -13664.3283
# No. of observations = 62416
# AIC value = 27412.6565


vd_supp_full <- vd_supp_full_glm %>% 
      tidy(conf.int = T, conf.level = 0.95, exponentiate = T, n_digits = 2) %>% 
      filter(term == "all_supp1_vitD and mineral supplement") %>% 
      dplyr::select(term, estimate, conf.low, conf.high)

vd_supp_full %<>% add_row(
      term = "No_supp",
      estimate = 1,
      conf.low=1,
      conf.high=1) # add reference

vd_supp_full$model <- "Fully adjusted" # marked the model


# add N
vd_supp_full_n <- vd_supp_full_glm %>% model.matrix() %>% 
      colSums()%>%                           # show the total number in each level
      as.data.frame() %>% t %>% as.data.frame() # transpose it make it easier to select


vd_supp_full %<>% mutate(N=case_when(
      term == "No_supp" ~ (vd_supp_full_n$`(Intercept)` - vd_supp_full_n$`all_supp1_vitD and mineral supplement`),
      term == "all_supp1_vitD and mineral supplement" ~ vd_supp_full_n$`all_supp1_vitD and mineral supplement`
))

rm(vd_supp_full_glm, vd_supp_full_n) # clean


fig3a_data <- bind_rows(vd_supp_crude, vd_supp_partial,vd_supp_full)

rm(list = ls(pattern = "supp"))# hosue keeping

fig3a_data %>% write_csv("dataset_wd/fig3a_data.csv")
# fig3a_data  <-  read_csv("dataset_wd/fig3a_data.csv")

# Fig3a: data management for plotting: ----
fig3a_data<- fig3a_data %>%
      mutate(conf.high2 = format(conf.high, nmall=2)) %>%             # 2 digits strings
      mutate(conf.low2 = format(conf.low, nmall=2)) %>%               # 2 digits strings
      mutate(conf.rr = format(estimate, nmall=2)) %>%                 # 2 digits strings
      mutate(ci = paste(conf.low2, conf.high2, sep = ",")) %>%        # lci, hci
      mutate(ci_l = paste0("(", ci)) %>%                              # (lci, hci
      mutate(ci_r = paste0(ci_l, ")")) %>%                            # (lci, hci)
      mutate(rr = paste(conf.rr, ci_r, sep = " ") ) %>%               # rr (lci, hci)
      mutate(RR = ifelse(estimate==1.00, "1 (Reference)",rr)) %>%     # set 1.0 into 1 (reference)
      mutate(vd_intake=case_when(
            term == "all_supp1_vitD and mineral supplement" ~ "Taking supplements",
            term == "No_supp" ~ "No supplements")) %>%     # rename for labelling
      mutate(model_order=case_when(
            model == "Crude" ~ 1,
            model == "Partially adjusted" ~ 2,
            model == "Fully adjusted" ~3)) %>%                      # for facet order
      mutate(Model = ifelse(vd_intake =="No supplements", 
                            model, NA)) %>%                           # for table labelling
      arrange( model_order,                                    # arrange order by vd status
               match(vd_intake, c("No supplements", "Taking supplements"))) %>% 
      dplyr::select(model_order, Model, vd_intake, N, estimate, conf.low, conf.high, RR)

# deal with factors issues:
fig3a_data[,1:3] <- map(fig3a_data[,1:3], as.factor)                      # factor coercion
fig3a_data$model_order %<>% factor(
      labels = c("Crude", "Partially adjusted", "Fully adjusted"))    # labelled for facet

fig3a_data$vd_intake %<>% relevel(ref = "Taking supplements")
fig3a_data$order <- 6:1

# Fig 3a: main foreat plot ----
fig3a <- fig3a_data %>% 
      ggplot(aes(y = order, x = estimate, 
                 xmin=conf.low, xmax=conf.high, color=vd_intake)) + 
      xlim(0.8,1.8)+
      geom_point(size=5, shape=19) + 
      geom_errorbarh(height=.3) +
      geom_vline(xintercept=1, lty=2) +
      ggtitle("        ")+
      theme_classic2() +
      theme(plot.title = element_text(size = 20 ,lineheight=.001),
            text = element_text(size = 20),
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_blank(),
            strip.text.y = element_blank(),
            legend.position = "none",
            axis.title.x = element_text(colour = "black", size = 20),
            axis.title.y = element_blank(),
            axis.text.x = element_text(colour = "black", size = 20),
            axis.text.y = element_blank())+ xlab("Rate ratio")
fig3a

# Fig 3a: left table ----

fig3a_table_l <- fig3a_data %>% 
      ggplot(aes(y = order)) + xlim(0,2) +
      geom_text(aes(x = 0, label = Model ),lineheight = 0.001, 
                hjust = 0, size = 6, colour = "black") + 
      geom_text(aes(x = 0.8, label = vd_intake),lineheight = 0.001, 
                hjust = 0 ,size = 6, colour = "black") + 
      geom_text(aes(x = 1.8, label = N),lineheight = 0.001, 
                hjust = 0 ,size = 6, colour = "black") + 
      theme_classic2() +
      ggtitle("   Model               Vitamin D intake         N")+    # use title for header
      theme(plot.title = element_text(size = 20 ,lineheight=.001, face="bold"),
            strip.text.x = element_blank(),
            strip.text.y = element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.border =  element_blank(),
            axis.title.x = element_text(colour = "white", size = 22), # make the text invisible
            axis.title.y = element_blank(),
            axis.line = element_line(colour = "white"),               # make the lines invisible
            axis.text.x = element_text(colour = "white", size = 22),  # make the text invisible
            axis.text.y = element_blank(),
            axis.ticks = element_blank()) + xlab(" ") 

fig3a_table_l

# Fig 3a: right table ----

fig3a_table_r <- fig3a_data %>% 
      ggplot(aes(y = order)) + xlim(0,0.5) +
      geom_text(aes(x = 0, label = RR),lineheight = 0.01, hjust = 0, size = 6, colour = "black") + 
      theme_classic2( ) +
      ggtitle("   RR(95%CI)")+
      theme(plot.title = element_text(size = 20 ,lineheight=.001, face="bold"),
            strip.text.x = element_blank(),
            strip.text.y = element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.border =  element_blank(),
            axis.title.x = element_text(colour = "white", size = 22),# make the text invisible
            axis.title.y = element_blank(),
            axis.text.x = element_text(colour = "white", size = 22), # make the text invisible
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_line(colour = "white"),              # make the lines invisible
            panel.background=element_blank(),
            plot.background=element_blank()) + xlab(" ") 

fig3a_table_r


# Fig3a: combine tables and forest plot ----

fig3a_combine <- grid.arrange(fig3a_table_l, fig3a, fig3a_table_r, ncol = 3)
# W:1800  H:400

# Fig 3b Vitamin D prescription and HZ ----

# crude ----
vd_drug_crude_glm <- glm(hz~vd_prescription + offset(log(fu_yr)), 
                         data = bd_i, family = "poisson")

vd_drug_crude_glm %>% idr.display()

# Poisson regression predicting hz with offset = NULL 
# 
# IDR(95%CI)        P(Wald's test) P(LR-test)
# vd_prescription: Had vitamin D prescription vs No vitamin D prescription 1.48 (1.27,1.71)  < 0.001        < 0.001   
#                                                                                                                     
# Log-likelihood = -36341.7258
# No. of observations = 177572
# AIC value = 72687.4515

vd_drug_crude <- vd_drug_crude_glm %>% 
      tidy(conf.int = T, conf.level = 0.95, exponentiate = T) %>% 
      filter(term != "(Intercept)") %>% 
      dplyr::select(term, estimate, conf.low, conf.high)

vd_drug_crude %<>% add_row(
      term = "No_supp",
      estimate = 1,
      conf.low=1,
      conf.high=1) # add reference

vd_drug_crude$model <- "Crude"

# add N
vd_drug_crude_n <- vd_drug_crude_glm %>% model.matrix() %>% 
      colSums()%>%                           # show the total number in each level
      as.data.frame() %>% t %>% as.data.frame() # transpose it make it easier to select

vd_drug_crude %<>%  mutate(N=case_when(
      term == "No_supp" ~ (vd_drug_crude_n$`(Intercept)` -vd_drug_crude_n$`vd_prescriptionHad vitamin D prescription`),
      term == "vd_prescriptionHad vitamin D prescription" ~ vd_drug_crude_n$`vd_prescriptionHad vitamin D prescription`
))

rm(vd_drug_crude_glm,vd_drug_crude_n)

# Partially adjusted for sex and age
vd_drug_partial_glm <- glm(hz~vd_prescription + offset(log(fu_yr)) + sex + age_c, 
                           data = bd_i, family = "poisson") 
vd_drug_partial_glm %>% idr.display()

# Poisson regression predicting hz with offset = NULL 
# 
# crude IDR(95%CI)  adj. IDR(95%CI)   P(Wald's test)
# vd_prescription: Had vitamin D prescription vs No vitamin D prescription 1.48 (1.27,1.71)  1.22 (1.05,1.42)  0.009         
#                                                                                                                            
# sex                                                                      0.82 (0.78,0.87)  0.81 (0.77,0.85)  < 0.001       
#                                                                                                                            
# age_c (cont. var.)                                                       1.03 (1.03,1.04)  1.03 (1.03,1.04)  < 0.001       
#                                                                                                                            
#                                                                          P(LR-test)
# vd_prescription: Had vitamin D prescription vs No vitamin D prescription 0.011     
#                                                                                    
# sex                                                                      < 0.001   
#                                                                                    
# age_c (cont. var.)                                                       < 0.001   
#                                                                                    
# Log-likelihood = -36086.8578
# No. of observations = 177572
# AIC value = 72181.7157

vd_drug_partial <- vd_drug_partial_glm %>% 
      tidy(conf.int = T, conf.level = 0.95, exponentiate = T) %>% 
      filter(term == "vd_prescriptionHad vitamin D prescription") %>% 
      dplyr::select(term, estimate, conf.low, conf.high)

vd_drug_partial %<>% add_row(
      term = "No_supp",
      estimate = 1,
      conf.low=1,
      conf.high=1) # add reference

vd_drug_partial$model <- "Partially adjusted"

# add N

vd_drug_partial_n <- vd_drug_partial_glm %>% model.matrix() %>% 
      colSums()%>%                           # show the total number in each level
      as.data.frame() %>% t %>% as.data.frame() # transpose it make it easier to select

vd_drug_partial %<>%  mutate(N=case_when(
      term == "No_supp" ~ (vd_drug_partial_n$`(Intercept)`-vd_drug_partial_n$`vd_prescriptionHad vitamin D prescription`),
      term == "vd_prescriptionHad vitamin D prescription" ~ vd_drug_partial_n$`vd_prescriptionHad vitamin D prescription`
))

rm(vd_drug_partial_glm,vd_drug_partial_n)  # clean


# Fully adjusted 
vd_drug_full_glm <- glm(hz~ vd_prescription + offset(log(fu_yr)) + sex + age_c + ethnic + 
                              bmi_group + drink_freq_c + smoke_stat + imd_bd_q + regions +
                              season_c + asthma + ckd+ copd + depress + dm + ibd + ra + sle + immunosuppression, 
                        data = bd_i, family = "poisson") 

vd_drug_full_glm %>% idr.display()

# Poisson regression predicting hz with offset = NULL 
# 
# crude IDR(95%CI)        adj. IDR(95%CI)        
# vd_prescription: Had vitamin D prescription vs No vitamin D prescription 1.49 (1.28,1.73)        1.09 (0.93,1.27)       
# 
# sex                                                                      0.82 (0.78,0.87)        0.8 (0.76,0.85)        
# 
# age_c (cont. var.)                                                       1.03 (1.03,1.04)        1.03 (1.03,1.04)       
# 
# ethnic: ref.=1_white                                                                                                    
# 2_mixed                                                               0.93 (0.65,1.34)        1.07 (0.74,1.54)       
# 3_asian                                                               0.79 (0.64,0.97)        0.94 (0.75,1.16)       
# 4_black                                                               0.71 (0.53,0.95)        0.87 (0.65,1.17)       
# 5_chinese                                                             0.9 (0.52,1.55)         1.09 (0.63,1.88)       
# 6_others                                                              0.78 (0.59,1.04)        0.88 (0.66,1.17)       
# 
# bmi_group: ref.=1_healthy weight                                                                                        
# 0_underweight                                                         1.08 (0.76,1.54)        1.07 (0.75,1.53)       
# 2_over weight                                                         1.04 (0.98,1.11)        1.02 (0.96,1.08)       
# 3_obese                                                               1.04 (0.97,1.11)        0.98 (0.92,1.05)       
# 
# drink_freq_c: ref.=never                                                                                                
# sometimes                                                             0.9725 (0.8817,1.0727)  1.0075 (0.9123,1.1128) 
# weekly                                                                0.94 (0.86,1.03)        1.03 (0.94,1.13)       
# daily                                                                 0.93 (0.84,1.03)        0.96 (0.87,1.07)       
# 
# smoke_stat: ref.=non-smoker                                                                                             
# ex-smoker                                                             1.19 (1.13,1.26)        1.12 (1.06,1.18)       
# current-smoker                                                        0.97 (0.89,1.06)        1.02 (0.93,1.11)       
# 
# imd_bd_q: ref.=least_deprived                                                                                           
# 2_deprived                                                            0.9966 (0.9197,1.0801)  0.997 (0.9196,1.081)   
# 3_deprived                                                            1.07 (0.99,1.15)        1.08 (0.99,1.16)       
# 4_deprived                                                            1.03 (0.96,1.12)        1.05 (0.97,1.14)       
# most_deprived                                                         0.99 (0.91,1.07)        1.02 (0.93,1.1)        
# 
# regions: ref.=East Midlands                                                                                             
# London                                                                0.89 (0.79,1)           0.92 (0.81,1.03)       
# North East                                                            1.02 (0.92,1.12)        1.03 (0.94,1.14)       
# North West                                                            1.18 (1.06,1.32)        1.17 (1.05,1.31)       
# South East                                                            0.94 (0.77,1.14)        0.99 (0.81,1.2)        
# South West                                                            1.05 (0.9,1.23)         1.08 (0.92,1.27)       
# West Midlands                                                         1 (0.88,1.13)           1.02 (0.9,1.16)        
# Yorkshire and The Humber                                              0.98 (0.89,1.07)        0.99 (0.9,1.08)        
# Wales                                                                 1.59 (1.44,1.76)        1.62 (1.46,1.79)       
# Scotland                                                              0.73 (0.66,0.81)        0.76 (0.68,0.84)       
# 
# season_c: ref.=2_Summer                                                                                                 
# 1_Spring                                                              1 (0.94,1.07)           0.97 (0.9,1.04)        
# 3_Autumn                                                              0.96 (0.89,1.03)        0.96 (0.89,1.03)       
# 4_Winter                                                              1.023 (0.9523,1.099)    0.992 (0.9207,1.0689)  
# 
# asthma                                                                   1.21 (1.13,1.29)        1.12 (1.05,1.2)        
# 
# ckd: CKD vs No CKD                                                       1.85 (1.51,2.28)        1.41 (1.14,1.74)       
# 
# copd: COPD vs No COPD                                                    1.64 (1.45,1.84)        1.28 (1.13,1.45)       
# 
# depress: depression vs No depression                                     1.09 (1.01,1.17)        1.07 (0.99,1.15)       
# 
# dm: Have DM vs No DM                                                     1.16 (1.05,1.28)        1.03 (0.93,1.14)       
# 
# ibd: Inflammatory bowel disease vs No inflammatory bowel diseaase        1.44 (1.27,1.64)        1.29 (1.13,1.47)       
# 
# ra: RA vs No RA                                                          1.74 (1.48,2.04)        1.24 (1.05,1.47)       
# 
# sle: SLE vs No SLE                                                       1.71 (1.09,2.68)        1.29 (0.82,2.03)       
# 
# immunosuppression: Immunosuppression vs Not immunosuppressive            1.96 (1.76,2.18)        1.57 (1.4,1.76)        
# 
# P(Wald's test) P(LR-test)
# vd_prescription: Had vitamin D prescription vs No vitamin D prescription 0.289          0.294     
#                                                                                                   
# sex                                                                      < 0.001        < 0.001   
#                                                                                                   
# age_c (cont. var.)                                                       < 0.001        < 0.001   
#                                                                                                   
# ethnic: ref.=1_white                                                                    0.825     
#    2_mixed                                                               0.732                    
#    3_asian                                                               0.557                    
#    4_black                                                               0.354                    
#    5_chinese                                                             0.756                    
#    6_others                                                              0.383                    
#                                                                                                   
# bmi_group: ref.=1_healthy weight                                                        0.597     
#    0_underweight                                                         0.702                    
#    2_over weight                                                         0.476                    
#    3_obese                                                               0.577                    
#                                                                                                   
# drink_freq_c: ref.=never                                                                0.273     
#    sometimes                                                             0.882                    
#    weekly                                                                0.524                    
#    daily                                                                 0.503                    
#                                                                                                   
# smoke_stat: ref.=non-smoker                                                             < 0.001   
#    ex-smoker                                                             < 0.001                  
#    current-smoker                                                        0.703                    
#                                                                                                   
# imd_bd_q: ref.=least_deprived                                                           0.227     
#    2_deprived                                                            0.942                    
#    3_deprived                                                            0.072                    
#    4_deprived                                                            0.213                    
#    most_deprived                                                         0.723                    
#                                                                                                   
# regions: ref.=East Midlands                                                             < 0.001   
#    London                                                                0.152                    
#    North East                                                            0.514                    
#    North West                                                            0.005                    
#    South East                                                            0.894                    
#    South West                                                            0.328                    
#    West Midlands                                                         0.792                    
#    Yorkshire and The Humber                                              0.739                    
#    Wales                                                                 < 0.001                  
#    Scotland                                                              < 0.001                  
#                                                                                                   
# season_c: ref.=2_Summer                                                                 0.578     
#    1_Spring                                                              0.36                     
#    3_Autumn                                                              0.228                    
#    4_Winter                                                              0.833                    
#                                                                                                   
# asthma                                                                   0.001          0.001     
#                                                                                                   
# ckd: CKD vs No CKD                                                       0.001          0.002     
#                                                                                                   
# copd: COPD vs No COPD                                                    < 0.001        < 0.001   
#                                                                                                   
# depress: depression vs No depression                                     0.084          0.087     
#                                                                                                   
# dm: Have DM vs No DM                                                     0.57           0.571     
#                                                                                                   
# ibd: Inflammatory bowel disease vs No inflammatory bowel diseaase        < 0.001        < 0.001   
#                                                                                                   
# ra: RA vs No RA                                                          0.012          0.015     
#                                                                                                   
# sle: SLE vs No SLE                                                       0.27           0.289     
#                                                                                                   
# immunosuppression: Immunosuppression vs Not immunosuppressive            < 0.001        < 0.001   
#                                                                                                   
# Log-likelihood = -34519.0292
# No. of observations = 171349
# AIC value = 69122.0583



vd_drug_full <- vd_drug_full_glm %>% 
      tidy(conf.int = T, conf.level = 0.95, exponentiate = T, n_digits = 2) %>% 
      filter(term == "vd_prescriptionHad vitamin D prescription") %>% 
      dplyr::select(term, estimate, conf.low, conf.high)

vd_drug_full %<>% add_row(
      term = "No_supp",
      estimate = 1,
      conf.low=1,
      conf.high=1) # add reference

vd_drug_full$model <- "Fully adjusted"


# add N

vd_drug_full_n <- vd_drug_full_glm %>% model.matrix() %>% 
      colSums()%>%                           # show the total number in each level
      as.data.frame() %>% t %>% as.data.frame() # transpose it make it easier to select

vd_drug_full %<>% mutate(N=case_when(
      term == "No_supp" ~ (vd_drug_full_n$`(Intercept)`-vd_drug_full_n$`vd_prescriptionHad vitamin D prescription`),
      term == "vd_prescriptionHad vitamin D prescription" ~ vd_drug_full_n$`vd_prescriptionHad vitamin D prescription`
))

rm(vd_drug_full_n, vd_drug_full_glm) # clean


# combine data
fig3b_data <- bind_rows(vd_drug_crude, vd_drug_partial, vd_drug_full) # combine output tables

fig3b_data %>% write_csv("dataset_wd/fig3b_data.csv")
# fig3b_data <- read_csv("dataset_wd/fig3b_data.csv")

rm(list = ls(pattern = "drug"))

# data management for plotting: ----

fig3b_data <- fig3b_data %>% 
      mutate(conf.high2 = format(conf.high, nmall=2)) %>%             # 2 digits strings
      mutate(conf.low2 = format(conf.low, nmall=2)) %>%               # 2 digits strings
      mutate(conf.rr = format(estimate, nmall=2)) %>%                 # 2 digits strings
      mutate(ci = paste(conf.low2, conf.high2, sep = ",")) %>%        # lci, hci
      mutate(ci_l = paste0("(", ci)) %>%                              # (lci, hci
      mutate(ci_r = paste0(ci_l, ")")) %>%                            # (lci, hci)
      mutate(rr = paste(conf.rr, ci_r, sep = " ") ) %>%               # rr (lci, hci)
      mutate(RR = ifelse(estimate==1.00, "1 (Reference)",rr)) %>%     # set 1.0 into 1 (reference)
      mutate(vd_intake=case_when(
            term == "vd_prescriptionHad vitamin D prescription" ~ "Receiving prescriptions",
            term == "No_supp" ~ "No prescriptions")) %>%     # rename for labelling
      mutate(model_order=case_when(
            model == "Crude" ~ 1,
            model == "Partially adjusted" ~ 2,
            model == "Fully adjusted" ~3)) %>%                      # for facet order
      mutate(Model = ifelse(vd_intake =="No prescriptions", 
                            model, NA)) %>%                           # for table labelling
      arrange( model_order,                                    # arrange order by vd status
               match(vd_intake, c("No prescriptions", "Receiving prescriptions"))) %>% 
      dplyr::select(model_order, Model, vd_intake,N , estimate, conf.low, conf.high, RR)


# deal with factors issues:
fig3b_data[,1:3] <- map(fig3b_data[,1:3], as.factor)                      # factor coercion
fig3b_data$model_order %<>% factor(
      labels = c("Crude", "Partially adjusted", "Fully adjusted"))    # labelled for facet

fig3b_data$vd_intake %<>% relevel(ref = "Receiving prescriptions")
fig3b_data$order <- 6:1

# Fig3b: forest plot ----

fig3b <- fig3b_data %>% 
      ggplot(aes(y = order, x = estimate, 
                 xmin=conf.low, xmax=conf.high, color=vd_intake)) + 
      xlim(0.8,1.8)+
      geom_point(size=5, shape=19) + 
      geom_errorbarh(height=.3) +
      geom_vline(xintercept=1, lty=2) +
      ggtitle("        ")+
      theme_classic2() +
      theme(plot.title = element_text(size = 20 ,lineheight=.001),
            text = element_text(size = 20),
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_blank(),
            strip.text.y = element_blank(),
            legend.position = "none",
            axis.title.x = element_text(colour = "black", size = 20),
            axis.title.y = element_blank(),
            axis.text.x = element_text(colour = "black", size = 20),
            axis.text.y = element_blank())+ xlab("Rate ratio")
fig3b


# fig3b: left table ----

fig3b_table_l <- fig3b_data %>% 
      ggplot(aes(y =order)) + xlim(0,2) +
      geom_text(aes(x = 0, label = Model ),lineheight = 0.001, 
                hjust = 0, size = 6, colour = "black") + 
      geom_text(aes(x = 0.8, label = vd_intake),lineheight = 0.001, 
                hjust = 0 ,size = 6, colour = "black") + 
      geom_text(aes(x = 1.8, label = N),lineheight = 0.001, 
                hjust = 0 ,size = 6, colour = "black") + 
      theme_classic2() +
      ggtitle("   Model               Vitamin D intake         N")+      # use title for header
      theme(plot.title = element_text(size = 20 ,lineheight=.001, face="bold"),
            strip.text.x = element_blank(),
            strip.text.y = element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.border =  element_blank(),
            axis.title.x = element_text(colour = "white", size = 22), # make the text invisible
            axis.title.y = element_blank(),
            axis.line = element_line(colour = "white"),               # make the lines invisible
            axis.text.x = element_text(colour = "white", size = 22),  # make the text invisible
            axis.text.y = element_blank(),
            axis.ticks = element_blank()) + xlab(" ") 
fig3b_table_l 

# fig3b: right table ----
fig3b_table_r <- fig3b_data %>% 
      ggplot(aes(y = order)) + xlim(0,0.5) +
      geom_text(aes(x = 0, label = RR),lineheight = 0.01, hjust = 0, size = 6, colour = "black") + 
      theme_classic2( ) +
      ggtitle("   RR(95%CI)")+
      theme(plot.title = element_text(size = 20 ,lineheight=.001, face="bold"),
            strip.text.x = element_blank(),
            strip.text.y = element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.border =  element_blank(),
            axis.title.x = element_text(colour = "white", size = 22),# make the text invisible
            axis.title.y = element_blank(),
            axis.text.x = element_text(colour = "white", size = 22), # make the text invisible
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_line(colour = "white"),              # make the lines invisible
            panel.background=element_blank(),
            plot.background=element_blank()) + xlab(" ") 

fig3b_table_r

# Fig3b: combine tables and forestplot
fig3b_combine <- grid.arrange(fig3b_table_l, fig3b, fig3b_table_r, ncol = 3)
# W:1800 H:400

# Combine two Figures using ggpubr ---

ggarrange(fig3a_combine, fig3b_combine,
          labels = c("a", "b"),
          ncol = 1, nrow = 2,
          align = "v",
          widths = 1)

# output size: W: 1500 H:1000