# Packages, dataset and setting -----

# packages needed
.libPaths("V:/VOLUMEO/uk_biobank/2.software/r_packages")
lapply(c("tidyverse", "markdown", "rmarkdown", "survival", "survminer",
         "broom", "magrittr", "ggpubr","grid","gridExtra","forcats"), 
       require, character.only=T)

# dataset:
bd_i <- read_rds("./dataset_wd/bd_i.RDS")

# option
options(digits = 2, scipen = 999)

# Function for organizing regression output
tidy_output <- function(x){
   x %>% tidy(conf.int = T, conf.level = 0.95,exponentiate = TRUE) %>% 
      dplyr::select(term, estimate,conf.low, conf.high)
}


# 0. Ristset 
bd_i$hz %<>% as.numeric() # prepare for coxph

risk_vd <- Surv(time = as.numeric(bd_i$time_in), 
                time2 = as.numeric(bd_i$time_out), 
                origin = as.numeric(bd_i$ymob), 
                event = bd_i$hz)

## f5. Vitamin D supplementation
# Crude model 
supp_crude_cox <- coxph(
   risk_vd ~ all_supp,
   data = bd_i,
   id = f.eid
)

supp_crude_cox %>% tidy_output()



# Partially adjusted for sex
supp_partial_cox <- coxph(
   risk_vd ~ all_supp + sex,
   data = bd_i,
   id = f.eid
)
supp_partial_cox %>% tidy_output()


# Fully adjusted
supp_full_cox <- coxph(
   risk_vd ~ all_supp + sex + ethnic + bmi_group + drink_freq_c + 
      smoke_stat + imd_bd_q + regions + season_c + asthma + ckd + 
      copd + depress + dm + ibd + ra + sle + immunosuppression,
   data = bd_i,
   id = f.eid
)

supp_full_cox %>% tidy_output()


# model check # not violate
supp_crude_cox %>% cox.zph() %>% print()
supp_partial_cox %>% cox.zph() %>% print()
supp_full_cox %>% cox.zph() %>% print() 

# f5a. Data management for plotting ----

sup_cox_data <- bind_rows((supp_crude_cox %>% tidy_output() %>% 
                              add_row( # add a reference
                                 term = "No_supp",
                                 estimate = 1,
                                 conf.low=1,
                                 conf.high=1) %>% mutate(model="Crude")),
                          (supp_partial_cox %>% tidy_output() %>% 
                              add_row(
                                 term = "No_supp",
                                 estimate = 1,
                                 conf.low=1,
                                 conf.high=1) %>% mutate(model="Partially adjusted")),
                          (supp_full_cox %>% tidy_output() %>% 
                              add_row(
                                 term = "No_supp",
                                 estimate = 1,
                                 conf.low=1,
                                 conf.high=1) %>% mutate(model="Fully adjusted"))) %>% 
   mutate(conf.high2 = format(conf.high, nsmall=2)) %>%             # 2 digits strings
   mutate(conf.low2 = format(conf.low, nsmall=2)) %>%               # 2 digits strings
   mutate(conf.rr = format(estimate, nsmall=2)) %>%                 # 2 digits strings
   mutate(ci = paste(conf.low2, conf.high2, sep = ",")) %>%        # lci, hci
   mutate(ci_l = paste0("(", ci)) %>%                              # (lci, hci
   mutate(ci_r = paste0(ci_l, ")")) %>%                            # (lci, hci)
   mutate(rr = paste(conf.rr, ci_r, sep = " ") ) %>%               # rr (lci, hci)
   mutate(HR = ifelse(estimate==1.00, "1 (Reference)",rr)) %>%     # set 1.0 into 1 (reference)
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
   dplyr::select(model_order, Model, vd_intake, estimate, conf.low, conf.high, HR) %>% 
   filter(!is.na(vd_intake))

# deal with factors issues:

sup_cox_data[,1:3] <- map(sup_cox_data[,1:3], as.factor)                      # factor coercion
sup_cox_data$model_order %<>% factor(
   labels = c("Crude", "Partially adjusted", "Fully adjusted"))    # labelled for facet

sup_cox_data$vd_intake %<>% relevel(ref = "Taking supplements")

sup_cox_data$order <- 6:1
sup_cox_data %>% write_csv("dataset_wd/supp_fig5a_data.csv")
# sup_cox_data <- read_csv("dataset_wd/supp_fig5a_data.csv")

# f5a: forestplot
f5a_forest <- sup_cox_data %>% 
   ggplot(aes(y = order, x = estimate, 
              xmin=conf.low, xmax=conf.high, color=vd_intake)) + 
   xlim(0,3)+
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
         axis.text.y = element_blank())+ xlab("Hazard ratio")
f5a_forest      

# f5a: left table
f5a_table_left <- sup_cox_data %>% 
   ggplot(aes(y = order)) + xlim(0,2) +
   geom_text(aes(x = 0, label = Model ),lineheight = 0.001, 
             hjust = 0, size = 6, colour = "black") + 
   geom_text(aes(x = 0.8, label = vd_intake),lineheight = 0.001, 
             hjust = 0 ,size = 6, colour = "black") + 
   theme_classic2() +
   ggtitle("   Model           Vitamin D intake ")+      # use title for header
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
f5a_table_left

# f5a: right table
f5a_table_right <- sup_cox_data %>% 
   ggplot(aes(y = order)) + xlim(0,0.5) +
   geom_text(aes(x = 0, label = HR),lineheight = 0.01, hjust = 0, size = 6, colour = "black") + 
   theme_classic2( ) +
   ggtitle("   HR(95%CI)")+
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
f5a_table_right

sup_f5a <- grid.arrange(f5a_table_left,f5a_forest,f5a_table_right, ncol = 3)


# f5b: vitamin D prescription

## Vitamin D prescriptions
# crude
vddrug_crude_cox <- coxph(
   risk_vd ~ vd_prescription,
   data = bd_i,
   id = f.eid
)
vddrug_crude_cox %>% tidy_output()


# Partially adjusted
vddrug_partial_cox <- coxph(
   risk_vd ~ vd_prescription + sex,
   data = bd_i,
   id = f.eid
)
vddrug_partial_cox %>% tidy_output()

# Full model
vddrug_full_cox <- coxph(
   risk_vd ~ vd_prescription + sex + ethnic + bmi_group + drink_freq_c + 
      smoke_stat + imd_bd_q + regions + season_c + asthma + ckd + 
      copd + depress + dm + ibd + ra + sle + immunosuppression,
   data = bd_i,
   id = f.eid
)
vddrug_full_cox %>% tidy_output()

# model check
vddrug_crude_cox %>% cox.zph() %>% print() # not violate
vddrug_partial_cox %>% cox.zph() %>% print() # not violate
vddrug_full_cox %>% cox.zph() %>% print() # not violate

# f5b: data management for plotting -----
drug_cox_data <- bind_rows((vddrug_crude_cox %>% tidy_output() %>% 
                               add_row(
                                  term = "No_drug",
                                  estimate = 1,
                                  conf.low=1,
                                  conf.high=1) %>% mutate(model="Crude")),
                           (vddrug_partial_cox %>% tidy_output()  %>% 
                               add_row(
                                  term = "No_drug",
                                  estimate = 1,
                                  conf.low=1,
                                  conf.high=1) %>% mutate(model="Partially adjusted")),
                           (vddrug_full_cox %>% tidy_output()  %>% 
                               add_row(
                                  term = "No_drug",
                                  estimate = 1,
                                  conf.low=1,
                                  conf.high=1) %>% mutate(model="Fully adjusted"))) %>% 
   mutate(conf.high2 = format(conf.high, nsmall=2)) %>%             # 2 digits strings
   mutate(conf.low2 = format(conf.low, nsmall=2)) %>%               # 2 digits strings
   mutate(conf.rr = format(estimate, nsmall=2)) %>%                 # 2 digits strings
   mutate(ci = paste(conf.low2, conf.high2, sep = ",")) %>%        # lci, hci
   mutate(ci_l = paste0("(", ci)) %>%                              # (lci, hci
   mutate(ci_r = paste0(ci_l, ")")) %>%                            # (lci, hci)
   mutate(rr = paste(conf.rr, ci_r, sep = " ") ) %>%               # rr (lci, hci)
   mutate(HR = ifelse(term == "No_drug", "1 (Reference)",rr)) %>%     # set 1.0 into 1 (reference)
   mutate(vd_intake=case_when(
      term == "vd_prescriptionHad vitamin D prescription" ~ "Receiving prescriptions",
      term == "No_drug" ~ "No prescriptions")) %>%     # rename for labelling
   mutate(model_order=case_when(
      model == "Crude" ~ 1,
      model == "Partially adjusted" ~ 2,
      model == "Fully adjusted" ~3)) %>%                      # for facet order
   mutate(Model = ifelse(vd_intake =="No prescriptions", 
                         model, NA)) %>%                           # for table labelling
   arrange( model_order,                                    # arrange order by vd status
            match(vd_intake, c("No prescriptions", "Receiving prescriptions"))) %>% 
   dplyr::select(model_order, Model, vd_intake, estimate, conf.low, conf.high, HR) %>% 
   filter(!is.na(vd_intake))

drug_cox_data[,1:3] <- map(drug_cox_data[,1:3], as.factor)
drug_cox_data$model_order%<>% factor(
   labels = c("Crude", "Partially adjusted", "Fully adjusted")) 

drug_cox_data$vd_intake %<>% relevel(ref = "Receiving prescriptions")

drug_cox_data$order <- 6:1

# save for later use
drug_cox_data %>% write_csv("dataset_wd/supp_fig5b_data.csv")
# drug_cox_data <- read_csv("dataset_wd/supp_fig5b_data.csv")


# f5b: forest plot
f5b_forest <- drug_cox_data %>% 
   ggplot(aes(y = order, x = estimate, 
              xmin=conf.low, xmax=conf.high, color=vd_intake)) + 
   xlim(0,3)+
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
         axis.text.y = element_blank())+ xlab("Hazard ratio")
f5b_forest


# f5b: table on the left
f5b_table_left <- drug_cox_data %>% 
   ggplot(aes(y = order)) + xlim(0,2) +
   geom_text(aes(x = 0, label = Model ),lineheight = 0.001, 
             hjust = 0, size = 6, colour = "black") + 
   geom_text(aes(x = 0.8, label = vd_intake),lineheight = 0.001, 
             hjust = 0 ,size = 6, colour = "black") + 
   theme_classic2() +
   ggtitle("   Model           Vitamin D intake ")+     # use title for header
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
f5b_table_left


# f5b: right table
f5b_table_right <- drug_cox_data %>% 
   ggplot(aes(y = order)) + xlim(0,0.5) +
   geom_text(aes(x = 0, label = HR),lineheight = 0.01, hjust = 0, size = 6, colour = "black") + 
   theme_classic2( ) +
   ggtitle("   HR(95%CI)")+
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
f5b_table_right

# combine
sup_f5b <- grid.arrange(f5b_table_left, f5b_forest, f5b_table_right, ncol = 3)


# Combine two Figures using ggpubr ---

ggarrange(sup_f5a, sup_f5b,
          labels = c("a", "b"),
          ncol = 1, nrow = 2,
          align = "v",
          widths = 1)

# size: W:2000  H: 1000