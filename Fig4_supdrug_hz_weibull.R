## 2nd exposure using Weibull regression ----

# data: 
source("./master_script.R")

# load additional packages
lapply(c("ggpubr","grid","gridExtra","forcats", "survival"), require, character.only=T)

options(digits = 3, scipen = 999)

stats_to_convert<- c("estimate",  "conf.low",  "conf.high")

# Data management of the regression model : 
# 1. weibull regression -> save the model
# 2. use tidy() to organize the regression model; 
# 2.5 Important: fixed the output number (reverse +/-)
# 3. exponentiation of the results
# 4. add reference row
# 5. use model.frame() to see how many people were included in the model
# 6. create new var N from 5.
# 7. rename the table and clean the global environment

# Functions to use:

# tidy output
tidy_reg <- function(x){
      tidy(x, conf.int = T, conf.level = 0.95, exponentiate = TRUE) %>% 
            filter(term != "(Intercept)") %>% 
            dplyr::select(term, estimate, conf.low, conf.high) %>% 
            filter(term == "all_supp1_vitD and mineral supplement")
}
tidy_reg_drug <- function(x){
      tidy(x, conf.int = T, conf.level = 0.95, exponentiate = TRUE) %>% 
            filter(term != "(Intercept)") %>% 
            dplyr::select(term, estimate, conf.low, conf.high) %>% 
            filter(term == "vd_prescriptionHad vitamin D prescription")
}

# rename the confidence interval
rename_ci <- function(x){
      x %>% 
            rename("conf.high2"="conf.low") %>% 
            rename("conf.low"="conf.high") %>% 
            rename("conf.high"="conf.high2")
}

# exponentiate and round
exp_round <- function(x){
      round(exp(x),digits = 2)
}

# add reference 
ref_row <- function(x){
      add_row(x,
              term = "No supp",
              estimate = 1,
              conf.low=1,
              conf.high=1)  # add the reference
}

# for extracting N
n_supp <- function(x){
      table(model.frame(x)$all_supp) %>% as.data.frame() %>% 
            pivot_wider(names_from = Var1, values_from =  Freq)
}
n_drug  <- function(x){
      table(model.frame(x)$vd_prescription) %>% as.data.frame() %>% 
            pivot_wider(names_from = Var1, values_from =  Freq)
}

# a. all_supp ----
# crude model -----
hz_supp_crude_w <- survreg(Surv(fu_yr, hz) ~all_supp, data = bd_i, dist = "weibull")

# tidy output
hz_vdsupp_crude_tidy <- hz_supp_crude_w %>% tidy_reg()

# reverse +/- 
hz_vdsupp_crude_tidy[,-1] <- hz_vdsupp_crude_tidy[,-1]*-1

# rename
hz_vdsupp_crude_tidy %<>% rename_ci()

# exp and round
hz_vdsupp_crude_tidy[stats_to_convert] <- map(hz_vdsupp_crude_tidy[stats_to_convert], exp_round)

# add ref row
hz_vdsupp_crude_tidy %<>% ref_row()

# get N
hz_supp_crude_n <- hz_supp_crude_w %>% n_supp()

# add N
hz_vdsupp_crude_tidy %<>% mutate(N=case_when(
      term == "all_supp1_vitD and mineral supplement" ~ hz_supp_crude_n$`1_vitD and mineral supplement`,
      term == "No supp" ~ hz_supp_crude_n$`0_no vitD supplement`
))
# rename and clean
supp_crude <- hz_vdsupp_crude_tidy
rm(list = ls(pattern = "hz"))


# Partial adjusted model -----
hz_vdsupp.partial <- survreg(Surv(fu_yr, hz) ~all_supp + sex + age_c, 
                             data = bd_i, dist = "weibull")

# tidy output
hz_supp_tidy <- hz_vdsupp.partial %>% tidy_reg()

# fix +/-
hz_supp_tidy[,-1] <- hz_supp_tidy[,-1]*-1

# rename
hz_supp_tidy %<>% rename_ci()

# exp and round
hz_supp_tidy[stats_to_convert] <- map(hz_supp_tidy[stats_to_convert], exp_round)

# add reference 
hz_supp_tidy %<>% ref_row()

# get N
hz_supp_partial_n <- hz_vdsupp.partial %>% n_supp()

# add N
hz_supp_tidy %<>% mutate(N=case_when(
      term == "all_supp1_vitD and mineral supplement" ~ hz_supp_partial_n$`1_vitD and mineral supplement`,
      term == "No supp" ~ hz_supp_partial_n$`0_no vitD supplement`
))

# rename and clean
supp_partial <- hz_supp_tidy
rm(list = ls(pattern = "hz"))

# Fully adjusted model -----
hz_vdsupp_full.w <- survreg(Surv(fu_yr, hz) ~all_supp + sex + age_c + ethnic + 
                                  bmi_group + drink_freq_c + smoke_stat + imd_bd_q + regions +
                                  season_c + asthma + ckd+ copd + depress + dm + ibd + ra + sle + immunosuppression, 
                            data = bd_i, dist = "weibull")

# tidy output
hz_supp_full_w <- hz_vdsupp_full.w %>% tidy_reg()

# fix +/-
hz_supp_full_w[,-1] <- hz_supp_full_w[,-1]*-1

# rename
hz_supp_full_w %<>% rename_ci

# exp and round
hz_supp_full_w[stats_to_convert] <- map(hz_supp_full_w[stats_to_convert], exp_round)

# add ref
hz_supp_full_w %<>% ref_row()

# get N
hz_supp_full_n <- hz_vdsupp_full.w %>% n_supp()

# add N
hz_supp_full_w %<>% mutate(N=case_when(
      term == "all_supp1_vitD and mineral supplement" ~ hz_supp_full_n$`1_vitD and mineral supplement`,
      term == "No supp" ~ hz_supp_full_n$`0_no vitD supplement`
))

# rename and clean
supp_full <- hz_supp_full_w
rm(list = ls(pattern = "hz"))

# add models -----
supp_crude$model <- "Crude"
supp_partial$model <- "Partially adjusted"
supp_full$model <- "Fully adjusted"


fig3a_data <- bind_rows(supp_crude, supp_partial, supp_full)

# save file
fig3a_data %>% write_csv("dataset_wd/fig3a_data_weibull.csv")
# fig3a_data <- read_csv("dataset_wd/fig3a_data_weibull.csv")

# data manipulation for figure ------
fig3a_data_w <- fig3a_data %>% 
      mutate(conf.high2 = format(conf.high, nmall=2)) %>%             # 2 digits strings
      mutate(conf.low2 = format(conf.low, nmall=2)) %>%               # 2 digits strings
      mutate(conf.rr = format(estimate, nmall=2)) %>%                 # 2 digits strings
      mutate(ci = paste(conf.low2, conf.high2, sep = ",")) %>%        # lci, hci
      mutate(ci_l = paste0("(", ci)) %>%                              # (lci, hci
      mutate(ci_r = paste0(ci_l, ")")) %>%                            # (lci, hci)
      mutate(rr = paste(conf.rr, ci_r, sep = " ") ) %>%               # rr (lci, hci)
      mutate(HR = ifelse(estimate==1.00, "1 (Reference)",rr)) %>%     # set 1.0 into 1 (reference)
      mutate(vd_intake=case_when(
            term == "all_supp1_vitD and mineral supplement" ~ "Taking supplements",
            term == "No supp" ~ "No supplements")) %>%     # rename for labelling
      mutate(model_order=case_when(
            model == "Crude" ~ 1,
            model == "Partially adjusted" ~ 2,
            model == "Fully adjusted" ~3)) %>%                      # for group order
      mutate(supp_order=ifelse(vd_intake=="No supplements", 1,2)) %>% 
      mutate(Model = ifelse(vd_intake =="No supplements", 
                            model, NA)) %>%                           # for table labelling
      arrange( model_order, supp_order) %>% 
      dplyr::select(model_order, Model, vd_intake, N, estimate, conf.low, conf.high, HR)

fig3a_data_w$order <- 6:1

# Plotting
# forest plot: ----
fig3a_plot_w <- fig3a_data_w %>% 
      ggplot(aes(y = order, x = estimate, 
                 xmin=conf.low, xmax=conf.high, color=vd_intake)) + 
      xlim(0.6,2.0)+
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
            axis.text.y = element_blank())+ xlab("Hazard ratio")
fig3a_plot_w

# Left table ----
fig3a_left_w <-  fig3a_data_w %>% 
      ggplot(aes(y = order)) + xlim(0,2) +
      geom_text(aes(x = 0, label = Model ),lineheight = 0.001, 
                hjust = 0, size = 6, colour = "black") + 
      geom_text(aes(x = 0.8, label = vd_intake),lineheight = 0.001, 
                hjust = 0 ,size = 6, colour = "black") + 
      geom_text(aes(x = 1.8, label = N),lineheight = 0.001, 
                hjust = 0 ,size = 6, colour = "black") + 
      theme_classic2() +
      ggtitle("   Model           Vitamin D supplement       N")+    # use title for header
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
fig3a_left_w 

# right table -----
fig3a_right_w <- fig3a_data_w %>% 
      ggplot(aes(y = order)) + xlim(0,0.5) +
      geom_text(aes(x = 0, label = HR),lineheight = 0.01, hjust = 0, size = 6, colour = "black") + 
      theme_classic2( ) +
      ggtitle("    HR(95%CI)")+
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
fig3a_right_w 

# combine 3a ----
fig3a_w_combine <- grid.arrange(fig3a_left_w,
                                fig3a_plot_w,
                                fig3a_right_w, ncol = 3)


# b: vd_prescription -------

# crude ----
vd_drugcrude.w <- survreg(Surv(fu_yr, hz) ~vd_prescription, 
                          data = bd_i, dist = "weibull")

# tidy output
vd_drug_crude_tidy <- vd_drugcrude.w %>% tidy_reg_drug()

# fix +/-
vd_drug_crude_tidy[,-1] <- vd_drug_crude_tidy[,-1]*-1

# rename
vd_drug_crude_tidy %<>% rename_ci

# exp and round
vd_drug_crude_tidy[stats_to_convert] <- map(vd_drug_crude_tidy[stats_to_convert], exp_round)

# add refernce row
vd_drug_crude_tidy %<>% ref_row()

# get N
vd_drug_crude_n <- vd_drugcrude.w %>% n_drug()

# add N
vd_drug_crude_tidy %<>% mutate(N=case_when(
      term == "vd_prescriptionHad vitamin D prescription" ~ vd_drug_crude_n$`Had vitamin D prescription`,
      term == "No supp" ~ vd_drug_crude_n$`No vitamin D prescription`
))

# rename and clean
drug_crude <- vd_drug_crude_tidy

rm(list = ls(pattern = "vd"))


# partially adjusted -----

# model:
vd_drug_partial.w <- survreg(Surv(fu_yr, hz) ~vd_prescription + sex + age_c, 
                             data = bd_i, dist = "weibull")
# tidy 
vd_drug_partial_tidy <- vd_drug_partial.w %>% tidy_reg_drug()

# fix +/-
vd_drug_partial_tidy[,-1] <- vd_drug_partial_tidy[,-1]*-1

# rename
vd_drug_partial_tidy %<>% rename_ci()

# exp and round
vd_drug_partial_tidy[stats_to_convert] <- map(vd_drug_partial_tidy[stats_to_convert],exp_round)

# add ref
vd_drug_partial_tidy %<>% ref_row

# get N
vd_drug_partial_n <- vd_drug_partial.w %>% n_drug()

# add N
vd_drug_partial_tidy %<>% mutate(N=case_when(
      term == "vd_prescriptionHad vitamin D prescription" ~ vd_drug_partial_n$`Had vitamin D prescription`,
      term == "No supp" ~ vd_drug_partial_n$`No vitamin D prescription`
))

# rename and clean
drug_partial <- vd_drug_partial_tidy
rm(list = ls(pattern = "vd"))

# full ----
# model
vd_drug_full.w <- survreg(Surv(fu_yr, hz) ~vd_prescription + sex + age_c + ethnic + 
                                bmi_group + drink_freq_c + smoke_stat + imd_bd_q + regions +
                                season_c + asthma + ckd+ copd + depress + dm + ibd + ra + sle + immunosuppression, 
                          data = bd_i, dist = "weibull")
# tidy 
vd_drug_full_tidy <- vd_drug_full.w %>% tidy_reg_drug()

# fix +/-
vd_drug_full_tidy[,-1] <- vd_drug_full_tidy[,-1]*-1

# rename
vd_drug_full_tidy %<>% rename_ci()

# exp round
vd_drug_full_tidy[stats_to_convert] <- map(vd_drug_full_tidy[stats_to_convert], exp_round)

# add ref
vd_drug_full_tidy %<>% ref_row()

# get N
vd_drug_full_n <- vd_drug_full.w %>% n_drug()

# add N
vd_drug_full_tidy %<>%  mutate(N=case_when(
      term == "vd_prescriptionHad vitamin D prescription" ~ vd_drug_full_n$`Had vitamin D prescription`,
      term == "No supp" ~ vd_drug_full_n$`No vitamin D prescription`
))

# rename and clean
drug_full <- vd_drug_full_tidy
rm(list = ls(pattern = "vd"))

# add model 
drug_crude$model <- "Crude"
drug_partial$model <- "Partially adjusted"
drug_full$model <- "Fully adjusted"

# combine data -----
fig3b_data <- bind_rows(drug_crude, drug_partial, drug_full)
# save files
fig3b_data %>% write_csv("dataset_wd/fig3b_data_weibull.csv")
# fig3b_data <- read_csv("dataset_wd/fig3b_data_weibull.csv")

# manipulation for plotting -----
fig3b_data_w <- fig3b_data %>% 
      mutate(conf.high2 = format(conf.high, nsmall=2)) %>%             # 2 digits strings
      mutate(conf.low2 = format(conf.low, nsmall=2)) %>%               # 2 digits strings
      mutate(conf.rr = format(estimate, nsmall=2)) %>%                 # 2 digits strings
      mutate(ci = paste(conf.low2, conf.high2, sep = ",")) %>%        # lci, hci
      mutate(ci_l = paste0("(", ci)) %>%                              # (lci, hci
      mutate(ci_r = paste0(ci_l, ")")) %>%                            # (lci, hci)
      mutate(rr = paste(conf.rr, ci_r, sep = " ") ) %>%               # rr (lci, hci)
      mutate(HR = ifelse(estimate==1.00, "1 (Reference)",rr)) %>%     # set 1.0 into 1 (reference)
      mutate(vd_intake=case_when(
            term == "vd_prescriptionHad vitamin D prescription" ~ "Receiving prescriptions",
            term == "No supp" ~ "No prescriptions")) %>%     # rename for labelling
      mutate(model_order=case_when(
            model == "Crude" ~ 1,
            model == "Partially adjusted" ~ 2,
            model == "Fully adjusted" ~3)) %>%                      # for model order
      mutate(supp_order=ifelse(vd_intake=="No prescriptions", 1,2)) %>% # exposure order
      mutate(Model = ifelse(vd_intake =="No prescriptions", 
                            model, NA)) %>%                           # for table labelling
      arrange( model_order, supp_order) %>% 
      dplyr::select(model_order, Model, vd_intake, N, estimate, conf.low, conf.high, HR)

fig3b_data_w$order <- 6:1

# 3b plotting ------
fig3b_plot_w <- fig3b_data_w %>% 
      ggplot(aes(y = order, x = estimate, 
                 xmin=conf.low, xmax=conf.high, color=vd_intake)) + 
      xlim(0.6,2.0)+
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
            axis.text.y = element_blank())+ xlab("Hazard ratio")
fig3b_plot_w

# 3b left table ----
fig3b_left_w <- fig3b_data_w %>% 
      ggplot(aes(y = order)) + xlim(0,2) +
      geom_text(aes(x = 0, label = Model ),lineheight = 0.001, 
                hjust = 0, size = 6, colour = "black") + 
      geom_text(aes(x = 0.8, label = vd_intake),lineheight = 0.001, 
                hjust = 0 ,size = 6, colour = "black") + 
      geom_text(aes(x = 1.8, label = N),lineheight = 0.001, 
                hjust = 0 ,size = 6, colour = "black") + 
      theme_classic2() +
      ggtitle("   Model      Vitamin D prescription       N")+    # use title for header
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
fig3b_left_w 


# 3b right table ----
fig3b_right_w <- fig3b_data_w %>% 
      ggplot(aes(y = order)) + xlim(0,0.5) +
      geom_text(aes(x = 0, label = HR),lineheight = 0.01, hjust = 0, size = 6, colour = "black") + 
      theme_classic2( ) +
      ggtitle("    HR(95%CI)")+
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
fig3b_right_w

# combine 3b ----
fig3b_w_combine <- grid.arrange(fig3b_left_w,
                                fig3b_plot_w,
                                fig3b_right_w, ncol = 3)



# Combine two Figures using ggpubr ---

ggarrange(fig3a_w_combine, fig3b_w_combine,
          labels = c("a", "b"),
          ncol = 1, nrow = 2,
          align = "v",
          widths = 1)

# output size: W: 1800 H:1000