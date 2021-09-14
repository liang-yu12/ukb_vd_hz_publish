# sensitivity analysis: exclude records after 2013 (vaccination programe)
# Supplementary figure 2 vitamin D supp/drug

# data: 
source("./master_script.R")

bd_i$hz_se %<>% as.numeric()

# additional library: 
lapply(c("ggpubr","grid","gridExtra","forcats"), require, character.only=T)

# Sup Fig 2a Vitamin D supplementation: regression ----
### crude----
vd_supp_crude_se <- glm(hz_se~all_supp + offset(log(fu_yr_se)), 
                        data = bd_i, family = "poisson") %>% 
   tidy(conf.int = T, conf.level = 0.95, exponentiate = T) %>% 
   filter(term != "(Intercept)") %>% 
   dplyr::select(term, estimate, conf.low, conf.high)

vd_supp_crude_se %<>% add_row(
   term = "No_supp",
   estimate = 1,
   conf.low=1,
   conf.high=1) # add reference

vd_supp_crude_se$model <- "Crude"

### partial ----
vd_supp_partial_se <- glm(hz_se~all_supp + offset(log(fu_yr_se)) + sex + age_c, 
                          data = bd_i, family = "poisson") %>% 
   tidy(conf.int = T, conf.level = 0.95, exponentiate = T) %>% 
   filter(term == "all_supp1_vitD and mineral supplement") %>% 
   dplyr::select(term, estimate, conf.low, conf.high)

vd_supp_partial_se %<>% add_row(
   term = "No_supp",
   estimate = 1,
   conf.low=1,
   conf.high=1) # add reference

vd_supp_partial_se$model <- "Partially adjusted" # marked the model

### Full adjusted model ----
vd_supp_full_se <- glm(hz_se~ all_supp + offset(log(fu_yr_se)) + sex + age_c + ethnic + 
                          bmi_group + drink_freq_c + smoke_stat + imd_bd_q + regions +
                          season_c + asthma + ckd+ copd + depress + dm + ibd + ra + sle + immunosuppression, 
                       data = bd_i, family = "poisson") %>% 
   tidy(conf.int = T, conf.level = 0.95, exponentiate = T, n_digits = 2) %>% 
   filter(term == "all_supp1_vitD and mineral supplement") %>% 
   dplyr::select(term, estimate, conf.low, conf.high)

vd_supp_full_se %<>% add_row(
   term = "No_supp",
   estimate = 1,
   conf.low=1,
   conf.high=1) # add reference

vd_supp_full_se$model <- "Fully adjusted" # marked the model


### Combine regression tables
sup_fig2a_data <- bind_rows(vd_supp_crude_se, vd_supp_partial_se,vd_supp_full_se)
sup_fig2a_data %>% write_csv("dataset_wd/sup_fig2a_data.csv")
# sup_fig2a_data <-  read_csv("dataset_wd/sup_fig2a_data.csv")

# Sup Fig2a: data management for plotting: ----
sup_fig2a_data <- sup_fig2a_data %>%    
   mutate(conf.high2 = format(conf.high, nsmall=2)) %>%             # 2 digits strings
   mutate(conf.low2 = format(conf.low, nsmall=2)) %>%               # 2 digits strings
   mutate(conf.rr = format(estimate, nsmall=2)) %>%                 # 2 digits strings
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
   dplyr::select(model_order, Model, vd_intake, estimate, conf.low, conf.high, RR)

# deal with factors issues:
sup_fig2a_data[,1:3] <- map(sup_fig2a_data[,1:3], as.factor)                      # factor coercion
sup_fig2a_data$model_order %<>% factor(
   labels = c("Crude", "Partially adjusted", "Fully adjusted"))    # labelled for facet

sup_fig2a_data$vd_intake %<>% relevel(ref = "Taking supplements")
sup_fig2a_data$order <- 6:1

### Sup Fig2a: forest plot ----
sup_fig2a <- sup_fig2a_data %>% 
   ggplot(aes(y = order, x = estimate, 
              xmin=conf.low, xmax=conf.high, color=vd_intake)) + 
   xlim(0.8,1.8)+
   geom_point(size=4, shape=15) + 
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
sup_fig2a

### Sup Fig2a: left table ----
sup_fig2a_table_l <- sup_fig2a_data %>% 
   ggplot(aes(y = order)) + xlim(0,2) +
   geom_text(aes(x = 0, label = Model ),lineheight = 0.001, 
             hjust = 0, size = 6, colour = "black") + 
   geom_text(aes(x = 0.8, label = vd_intake),lineheight = 0.001, 
             hjust = 0 ,size = 6, colour = "black") + 
   theme_classic2() +
   ggtitle("   Model            Vitamin D intake ")+      # use title for header
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

sup_fig2a_table_l



### Sup Fig2b: right table ----
sup_fig2a_table_r <- sup_fig2a_data %>% 
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

sup_fig2a_table_r


# Sup Fig2a: combine tables and plots ----
sup_fig2a_combine <- grid.arrange(sup_fig2a_table_l, sup_fig2a, sup_fig2a_table_r, ncol = 3)


# Sup Fig 2b Vitamin D prescription and HZ: regression ----

# crude ----
vd_drug_crude_se <- glm(hz_se~vd_prescription + offset(log(fu_yr_se)), 
                        data = bd_i, family = "poisson") %>% 
   tidy(conf.int = T, conf.level = 0.95, exponentiate = T) %>% 
   filter(term != "(Intercept)") %>% 
   dplyr::select(term, estimate, conf.low, conf.high)

vd_drug_crude_se %<>% add_row(
   term = "No_supp",
   estimate = 1,
   conf.low=1,
   conf.high=1) # add reference

vd_drug_crude_se$model <- "Crude"

# Partially adjusted for sex and age
vd_drug_partial_se <- glm(hz_se~vd_prescription + offset(log(fu_yr_se)) + sex + age_c, 
                          data = bd_i, family = "poisson") %>% 
   tidy(conf.int = T, conf.level = 0.95, exponentiate = T) %>% 
   filter(term == "vd_prescriptionHad vitamin D prescription") %>% 
   dplyr::select(term, estimate, conf.low, conf.high)

vd_drug_partial_se %<>% add_row(
   term = "No_supp",
   estimate = 1,
   conf.low=1,
   conf.high=1) # add reference

vd_drug_partial_se$model <- "Partially adjusted"


# Fully adjusted 
vd_drug_full_se <- glm(hz_se~ vd_prescription + offset(log(fu_yr_se)) + sex + age_c + ethnic + 
                          bmi_group + drink_freq_c + smoke_stat + imd_bd_q + regions +
                          season_c + asthma + ckd+ copd + depress + dm + ibd + ra + sle + immunosuppression, 
                       data = bd_i, family = "poisson") %>% 
   tidy(conf.int = T, conf.level = 0.95, exponentiate = T, n_digits = 2) %>% 
   filter(term == "vd_prescriptionHad vitamin D prescription") %>% 
   dplyr::select(term, estimate, conf.low, conf.high)

vd_drug_full_se %<>% add_row(
   term = "No_supp",
   estimate = 1,
   conf.low=1,
   conf.high=1) # add reference

vd_drug_full_se$model <- "Fully adjusted"

sup_fig2b_data <- bind_rows(vd_drug_crude_se, vd_drug_partial_se, vd_drug_full_se) # combine output tables
sup_fig2b_data %>% write_csv("dataset_wd/sup_fig2b_data.csv")
# sup_fig2b_data <- read_csv("dataset_wd/sup_fig2b_data.csv")

# Sup Fig2b: data management for plotting: ----

sup_fig2b_data %<>% 
   mutate(conf.high2 = format(conf.high, nsmall=2)) %>%             # 2 digits strings
   mutate(conf.low2 = format(conf.low, nsmall=2)) %>%               # 2 digits strings
   mutate(conf.rr = format(estimate, nsmall=2)) %>%                 # 2 digits strings
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
   dplyr::select(model_order, Model, vd_intake, estimate, conf.low, conf.high, RR)


# deal with factors issues:
sup_fig2b_data[,1:3] <- map(sup_fig2b_data[,1:3], as.factor)                      # factor coercion
sup_fig2b_data$model_order %<>% factor(
   labels = c("Crude", "Partially adjusted", "Fully adjusted"))    # labelled for facet
sup_fig2b_data$vd_intake %<>% relevel(ref = "Receiving prescriptions")
sup_fig2b_data$order <- 6:1

### Sup Fig 2b: forest plot ----

sup_fig2b <- sup_fig2b_data %>% 
   ggplot(aes(y = order, x = estimate, 
              xmin=conf.low, xmax=conf.high, color=vd_intake)) + 
   xlim(0.8,1.8)+
   geom_point(size=4, shape=15) + 
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
sup_fig2b


### Sup Fig 2b: left table ----

sup_fig2b_table_l <- sup_fig2b_data %>% 
   ggplot(aes(y = order)) + xlim(0,2) +
   geom_text(aes(x = 0, label = Model ),lineheight = 0.001, 
             hjust = 0, size = 6, colour = "black") + 
   geom_text(aes(x = 0.8, label = vd_intake),lineheight = 0.001, 
             hjust = 0 ,size = 6, colour = "black") + 
   theme_classic2() +
   ggtitle("   Model            Vitamin D intake ")+      # use title for header
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

sup_fig2b_table_l   

### Sup fig2b: right table
sup_fig2b_table_r <- sup_fig2b_data %>% 
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

sup_fig2b_table_r

# Sup fig2b: combine tables and plots ----
sup_fig2b_combine <- grid.arrange(sup_fig2b_table_l, sup_fig2b, sup_fig2b_table_r, ncol = 3)


# Combine two Figures using ggpubr ---

ggarrange(sup_fig2a_combine, sup_fig2b_combine,
          labels = c("a", "b"),
          ncol = 1, nrow = 2,
          align = "v",
          widths = 0.5) 

# size: W: 1500  H: 1000