# Script for visualizing regression model 
# the association between vitamin D status and herpes zoster using poission regression

# data: 
source("./master_script.R")

# load additional packages
lapply(c("ggpubr","grid","gridExtra","forcats"), require, character.only=T)

# option
options(digits = 2, scipen = 999)

# Create a survival object for person years-----
hz_surv <- Surv(time = as.numeric(bd_i$time_in)/365.25, 
                time2 = as.numeric(bd_i$time_out)/365.25, 
                event = bd_i$hz)

# Number, person-year, rate by variables ----

# Vitamin D status
pyears(hz_surv ~ vitd_s, data = bd_i) %>% summary(n = T, rate = T, ci.r = T)

# Other coavariates
covars <- c("sex" ,"ethnic" ,"bmi_group" ,"drink_freq_c" ,"smoke_stat" ,
            "imd_bd_q" ,"regions" ,"season_c" ,"asthma" ,"ckd" ,"copd" ,"depress" ,"dm" ,
            "ibd" ,"ra" ,"sle" ,"immunosuppression")

map(bd_i[covars], function(x){pyears(hz_surv ~ x, data = bd_i) %>% summary(n = T, rate = T, ci.r = T)})
# need to find a better way to automatically save the results into a data.frame


## Poisson regression ---
# crude analysis: 

crude <- glm(hz~ vitd_s + offset(log(fu_yr)), 
             data = bd_i, family = "poisson")

crude %>% idr.display() # take a look of the results
# Poisson regression predicting hz with offset = NULL 
# 
# IDR(95%CI)        P(Wald's test) P(LR-test)
# vitd_s: ref.=2_sufficiency                                  0.006     
#    0_deficiency            0.89 (0.82,0.96)  0.002                    
#    1_insufficiency         0.97 (0.93,1.03)  0.332                    
#                                                                       
# Log-likelihood = -36348.2134
# No. of observations = 177572
# AIC value = 72702.4268

# use model matrix to show the number in each factor level
crude_n <- crude %>% model.matrix() %>% 
      colSums()%>%                           # show the total number in each level
      as.data.frame() %>% t %>% as.data.frame() # transpose it make it easier to select

crude_reg <-crude %>% 
      tidy(conf.int = T, conf.level = 0.95, exponentiate = T) %>% 
      filter(term != "(Intercept)") %>% 
      dplyr::select(term, estimate, conf.low, conf.high)

crude_reg %<>% add_row(
      term = "vitd_s2_sufficiency",
      estimate = 1,
      conf.low=1,
      conf.high=1)

crude_reg$model <- "Crude"

crude_reg %<>% mutate(N=case_when(
      term == "vitd_s2_sufficiency" ~ (crude_n$`(Intercept)` - crude_n$vitd_s0_deficiency - crude_n$vitd_s1_insufficiency),
      term == "vitd_s1_insufficiency" ~ crude_n$vitd_s1_insufficiency,
      term == "vitd_s0_deficiency" ~ crude_n$vitd_s0_deficiency
))

# house keeping
rm(crude_n)

# Partially adjusted for sex and age_c

partial <- glm(hz~ vitd_s + sex + age_c + offset(log(fu_yr)), 
               data = bd_i, family = "poisson")

partial %>% idr.display()

# Poisson regression predicting hz with offset = NULL 
# 
# crude IDR(95%CI)  adj. IDR(95%CI)   P(Wald's test) P(LR-test)
# vitd_s: ref.=2_sufficiency                                                    0.46      
#    0_deficiency            0.89 (0.82,0.96)  0.96 (0.89,1.04)  0.325                    
#    1_insufficiency         0.97 (0.93,1.03)  1.01 (0.96,1.06)  0.703                    
#                                                                                         
# sex: Male vs Female        0.82 (0.78,0.87)  0.8 (0.76,0.84)   < 0.001        < 0.001   
#                                                                                         
# age_c (cont. var.)         1.03 (1.03,1.04)  1.03 (1.03,1.04)  < 0.001        < 0.001   
#                                                                                         
# Log-likelihood = -36089.3248
# No. of observations = 177572
# AIC value = 72188.6495 

partial_n <- partial %>% model.matrix() %>% 
      colSums()%>%                           # show the total number in each level
      as.data.frame() %>% t %>% as.data.frame() %>% # transpose it make it easier to select
      dplyr::select(c("(Intercept)", "vitd_s0_deficiency", "vitd_s1_insufficiency"))

partial_reg <- partial %>% 
      tidy(conf.int = T, conf.level = 0.95, exponentiate = T) %>% 
      filter(term == "vitd_s0_deficiency" | term == "vitd_s1_insufficiency") %>% 
      dplyr::select(term, estimate, conf.low, conf.high)

partial_reg %<>% add_row(
      term = "vitd_s2_sufficiency",
      estimate = 1,
      conf.low=1,
      conf.high=1)

partial_reg$model <- "Partially adjusted" 

# add n to the table

partial_reg %<>% mutate(N=case_when(
      term == "vitd_s2_sufficiency" ~ (partial_n$`(Intercept)` - partial_n$vitd_s0_deficiency - partial_n$vitd_s1_insufficiency),
      term == "vitd_s1_insufficiency" ~ partial_n$vitd_s1_insufficiency,
      term == "vitd_s0_deficiency" ~ partial_n$vitd_s0_deficiency
))

rm(partial_n)# house keeping


# Full adjusted for all covariates
full <- glm(hz~ vitd_s + offset(log(fu_yr)) + sex + age_c + ethnic + 
                  bmi_group + drink_freq_c + smoke_stat + imd_bd_q + regions +
                  season_c + asthma + ckd+ copd + depress + dm + ibd + ra + sle + 
                  immunosuppression, 
            data = bd_i, family = "poisson")

full %>% idr.display()

# Poisson regression predicting hz with offset = NULL 
# 
# crude IDR(95%CI)        adj. IDR(95%CI)         P(Wald's test)
# vitd_s: ref.=2_sufficiency                                                                                                      
#    0_deficiency                                                   0.8916 (0.8252,0.9634)  0.9957 (0.9154,1.083)   0.92          
#    1_insufficiency                                                0.98 (0.93,1.03)        1.02 (0.97,1.08)        0.419         
#                                                                                                                                 
# sex                                                               0.82 (0.78,0.87)        0.8 (0.76,0.84)         < 0.001       
#                                                                                                                                 
# age_c (cont. var.)                                                1.03 (1.03,1.04)        1.03 (1.03,1.04)        < 0.001       
#                                                                                                                                 
# ethnic: ref.=1_white                                                                                                            
#    2_mixed                                                        0.93 (0.65,1.34)        1.06 (0.74,1.53)        0.739         
#    3_asian                                                        0.79 (0.64,0.97)        0.95 (0.76,1.18)        0.613         
#    4_black                                                        0.71 (0.53,0.95)        0.87 (0.65,1.17)        0.353         
#    5_chinese                                                      0.9 (0.52,1.55)         1.09 (0.63,1.88)        0.761         
#    6_others                                                       0.78 (0.59,1.04)        0.88 (0.66,1.18)        0.389         
#                                                                                                                                 
# bmi_group: ref.=1_healthy weight                                                                                                
#    0_underweight                                                  1.08 (0.76,1.54)        1.08 (0.75,1.53)        0.689         
#    2_over weight                                                  1.04 (0.98,1.11)        1.02 (0.96,1.08)        0.511         
#    3_obese                                                        1.04 (0.97,1.11)        0.98 (0.91,1.05)        0.523         
#                                                                                                                                 
# drink_freq_c: ref.=never                                                                                                        
#    sometimes                                                      0.9725 (0.8817,1.0727)  1.0064 (0.9112,1.1115)  0.9           
#    weekly                                                         0.94 (0.86,1.03)        1.03 (0.94,1.13)        0.537         
#    daily                                                          0.93 (0.84,1.03)        0.96 (0.87,1.07)        0.494         
#                                                                                                                                 
# smoke_stat: ref.=non-smoker                                                                                                     
#    ex-smoker                                                      1.19 (1.13,1.26)        1.12 (1.06,1.18)        < 0.001       
#    current-smoker                                                 0.97 (0.89,1.06)        1.02 (0.93,1.11)        0.705         
#                                                                                                                                 
# imd_bd_q: ref.=least_deprived                                                                                                   
#    2_deprived                                                     0.9966 (0.9197,1.0801)  0.9966 (0.9191,1.0806)  0.934         
#    3_deprived                                                     1.07 (0.99,1.15)        1.07 (0.99,1.16)        0.073         
#    4_deprived                                                     1.03 (0.96,1.12)        1.05 (0.97,1.14)        0.214         
#    most_deprived                                                  0.99 (0.91,1.07)        1.02 (0.93,1.1)         0.727         
#                                                                                                                                 
# regions: ref.=East Midlands                                                                                                     
#    London                                                         0.89 (0.79,1)           0.92 (0.81,1.03)        0.153         
#    North East                                                     1.02 (0.92,1.12)        1.03 (0.94,1.14)        0.521         
#    North West                                                     1.18 (1.06,1.32)        1.17 (1.05,1.31)        0.005         
#    South East                                                     0.94 (0.77,1.14)        0.98 (0.81,1.2)         0.877         
#    South West                                                     1.05 (0.9,1.23)         1.08 (0.92,1.27)        0.323         
#    West Midlands                                                  1 (0.88,1.13)           1.02 (0.89,1.15)        0.805         
#    Yorkshire and The Humber                                       0.98 (0.89,1.07)        0.99 (0.9,1.08)         0.75          
#    Wales                                                          1.59 (1.44,1.76)        1.61 (1.46,1.79)        < 0.001       
#    Scotland                                                       0.73 (0.66,0.81)        0.76 (0.68,0.84)        < 0.001       
#                                                                                                                                 
# season_c: ref.=2_Summer                                                                                                         
#    1_Spring                                                       1 (0.94,1.07)           0.97 (0.9,1.04)         0.338         
#    3_Autumn                                                       0.96 (0.89,1.03)        0.96 (0.89,1.03)        0.219         
#    4_Winter                                                       1.02 (0.95,1.1)         0.99 (0.92,1.07)        0.787         
#                                                                                                                                 
# asthma                                                            1.21 (1.13,1.29)        1.12 (1.05,1.21)        0.001         
#                                                                                                                                 
# ckd: CKD vs No CKD                                                1.85 (1.51,2.28)        1.41 (1.15,1.74)        0.001         
#                                                                                                                                 
# copd: COPD vs No COPD                                             1.64 (1.45,1.84)        1.28 (1.13,1.45)        < 0.001       
#                                                                                                                                 
# depress: depression vs No depression                              1.09 (1.01,1.17)        1.07 (0.99,1.15)        0.082         
#                                                                                                                                 
# dm: Have DM vs No DM                                              1.16 (1.05,1.28)        1.03 (0.93,1.14)        0.571         
#                                                                                                                                 
# ibd: Inflammatory bowel disease vs No inflammatory bowel diseaase 1.44 (1.27,1.64)        1.29 (1.13,1.47)        < 0.001       
#                                                                                                                                 
# ra: RA vs No RA                                                   1.74 (1.48,2.04)        1.24 (1.05,1.47)        0.011         
#                                                                                                                                 
# sle: SLE vs No SLE                                                1.71 (1.09,2.68)        1.3 (0.83,2.05)         0.255         
#                                                                                                                                 
# immunosuppression: Immunosuppression vs Not immunosuppressive     1.96 (1.76,2.18)        1.58 (1.41,1.78)        < 0.001       
#                                                                                                                                 
#                                                                   P(LR-test)
# vitd_s: ref.=2_sufficiency                                        0.649     
#    0_deficiency                                                             
#    1_insufficiency                                                          
#                                                                             
# sex                                                               < 0.001   
#                                                                             
# age_c (cont. var.)                                                < 0.001   
#                                                                             
# ethnic: ref.=1_white                                              0.842     
#    2_mixed                                                                  
#    3_asian                                                                  
#    4_black                                                                  
#    5_chinese                                                                
#    6_others                                                                 
#                                                                             
# bmi_group: ref.=1_healthy weight                                  0.579     
#    0_underweight                                                            
#    2_over weight                                                            
#    3_obese                                                                  
#                                                                             
# drink_freq_c: ref.=never                                          0.275     
#    sometimes                                                                
#    weekly                                                                   
#    daily                                                                    
#                                                                             
# smoke_stat: ref.=non-smoker                                       < 0.001   
#    ex-smoker                                                                
#    current-smoker                                                           
#                                                                             
# imd_bd_q: ref.=least_deprived                                     0.226     
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
# season_c: ref.=2_Summer                                           0.566     
#    1_Spring                                                                 
#    3_Autumn                                                                 
#    4_Winter                                                                 
#                                                                             
# asthma                                                            0.001     
#                                                                             
# ckd: CKD vs No CKD                                                0.002     
#                                                                             
# copd: COPD vs No COPD                                             < 0.001   
#                                                                             
# depress: depression vs No depression                              0.085     
#                                                                             
# dm: Have DM vs No DM                                              0.573     
#                                                                             
# ibd: Inflammatory bowel disease vs No inflammatory bowel diseaase < 0.001   
#                                                                             
# ra: RA vs No RA                                                   0.014     
#                                                                             
# sle: SLE vs No SLE                                                0.275     
#                                                                             
# immunosuppression: Immunosuppression vs Not immunosuppressive     < 0.001   
#                                                                             
# Log-likelihood = -34519.1463
# No. of observations = 171349
# AIC value = 69124.2925


full_reg <- full %>% 
      tidy(conf.int = T, conf.level = 0.95, exponentiate = T) %>% 
      filter(term == "vitd_s0_deficiency" | term == "vitd_s1_insufficiency") %>% 
      dplyr::select(term, estimate, conf.low, conf.high)

full_reg %<>% add_row(
      term = "vitd_s2_sufficiency",
      estimate = 1,
      conf.low=1,
      conf.high=1)

full_reg$model <- "Fully adjusted"

# add N
full_n <-  full %>% model.matrix() %>% 
      colSums()%>%                           # show the total number in each level
      as.data.frame() %>% t %>% as.data.frame() %>% # transpose it make it easier to select
      dplyr::select(c("(Intercept)", "vitd_s0_deficiency", "vitd_s1_insufficiency"))

full_reg %<>% mutate(N=case_when(
      term == "vitd_s2_sufficiency" ~ (full_n$`(Intercept)` - full_n$vitd_s0_deficiency - full_n$vitd_s1_insufficiency),
      term == "vitd_s1_insufficiency" ~ full_n$vitd_s1_insufficiency,
      term == "vitd_s0_deficiency" ~ full_n$vitd_s0_deficiency
))

rm(full_n)

# combine the output tables

fig2_data <- bind_rows(crude_reg, partial_reg, full_reg)
# term                  estimate conf.low conf.high model                  N
# <chr>                    <dbl>    <dbl>     <dbl> <chr>              <dbl>
# 1 vitd_s0_deficiency       0.885    0.820     0.955 Crude              25274
# 2 vitd_s1_insufficiency    0.975    0.926     1.03  Crude              74963
# 3 vitd_s2_sufficiency      1        1         1     Crude              77335
# 4 vitd_s0_deficiency       0.962    0.891     1.04  Partially adjusted 25274
# 5 vitd_s1_insufficiency    1.01     0.959     1.06  Partially adjusted 74963
# 6 vitd_s2_sufficiency      1        1         1     Partially adjusted 77335
# 7 vitd_s0_deficiency       0.996    0.915     1.08  Fully adjusted     24248
# 8 vitd_s1_insufficiency    1.02     0.968     1.08  Fully adjusted     72464
# 9 vitd_s2_sufficiency      1        1         1     Fully adjusted     74637

fig2_data %>% write_csv("dataset_wd/fig2_data.csv")
# fig2_data <- read_csv("dataset_wd/fig2_data.csv")


# Step 1: data manipulation of tidyr output table -----

fig2_data<- fig2_data %>%
      mutate(conf.high2 = format(conf.high, nsmall=2)) %>%             # 2 digits strings
      mutate(conf.low2 = format(conf.low, nsmall=2)) %>%               # 2 digits strings
      mutate(conf.rr = format(estimate, nsmall=2)) %>%                 # 2 digits strings
      mutate(ci = paste(conf.low2, conf.high2, sep = ",")) %>%        # lci, hci
      mutate(ci_l = paste0("(", ci)) %>%                              # (lci, hci
      mutate(ci_r = paste0(ci_l, ")")) %>%                            # (lci, hci)
      mutate(rr = paste(conf.rr, ci_r, sep = " ") ) %>%               # rr (lci, hci)
      mutate(RR = ifelse(estimate==1.00, "1 (Reference)",rr)) %>%     # set 1.0 into 1 (reference)
      mutate(vd_status=case_when(
            term == "vitd_s0_deficiency" ~ "Deficiency",
            term == "vitd_s1_insufficiency" ~ "Insufficiency",
            term == "vitd_s2_sufficiency" ~ "Sufficiency")) %>%     # rename for labelling
      mutate(model_order=case_when(
            model == "Crude" ~ 1,
            model == "Partially adjusted" ~ 2,
            model == "Fully adjusted" ~3)) %>%                      # for facet order
      mutate(Model = ifelse(vd_status =="Sufficiency", 
                            model, NA)) %>%                           # for table labelling
      arrange( model_order,                                    # arrange order by vd status
               match(vd_status, c("Sufficiency", "Insufficiency", "Deficiency"))) %>% 
      dplyr::select(model_order, Model, vd_status, N, estimate, conf.low, conf.high, RR)


# deal with factors issues:
fig2_data[,1:3] <- map(fig2_data[,1:3], as.factor)                      # factor coercion
fig2_data$model_order %<>% factor(
      labels = c("Crude", "Partially adjusted", "Fully adjusted"))    # labelled for facet
fig2_data

fig2_data$order <- 9:1
fig2_data$order %<>% as.factor
fig2_data$order %>% levels()

# Step 2: main forest plot -----
fig2_plot <- fig2_data %>% 
      ggplot(aes(y = order, x = estimate, xmin=conf.low, xmax=conf.high, color=vd_status)) + 
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
            axis.text.y = element_blank())+
      # ggtitle("Vitamin D status and the risk of herpes zoster") +
      xlab("Rate ratio")

fig2_plot

# Step 3: Plot the table on the left -----

fig2_table <- 
      ggplot(data = fig2_data, aes(y = order)) + xlim(0,0.8) +
      geom_text(aes(x = 0, label = Model ),lineheight = 0.001, hjust = 0, size = 6, colour = "black") + 
      geom_text(aes(x = 0.4, label = vd_status),lineheight = 0.001, hjust = 0 ,size = 6, colour = "black") +
      geom_text(aes(x = 0.6, label = N),lineheight = 0.001, hjust = 0 ,size = 6, colour = "black") + 
      theme_classic2() + 
      ggtitle("   Model          Vitamin D status     N")+   
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
fig2_table


# Step 4: Plot the table on the right -----

fig2_table2 <- 
      ggplot(data = fig2_data,  aes(y = order)) + xlim(0,0.5) +
      geom_text(aes(x = 0, label = RR),lineheight = 0.01, hjust = 0, size = 6, colour = "black") + 
      theme_classic2( ) +
      ggtitle("    RR(95%CI)")+
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

fig2_table2

# Step 5: combine the three figure together

grid.arrange(fig2_table, fig2_plot, fig2_table2, ncol=3)

# size: W:1800 H:600
