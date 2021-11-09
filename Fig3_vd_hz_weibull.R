## Model check using Weibull regression ----

# data: 
source("./master_script.R")

# load additional packages
lapply(c("ggpubr","grid","gridExtra","forcats", "survival", "survminer"), require, character.only=T)

options(digits = 3, scipen = 999)

stats_to_convert<- c("estimate",  "conf.low",  "conf.high")

# Data management of the regression model : 
# 1. weibull regression -> save the model
# 2. use tidy() to organize the regression model
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
            filter(term == "vitd_s0_deficiency" | term == "vitd_s1_insufficiency")
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
              term = "vitd_s2_sufficiency",
              estimate = 1,
              conf.low=1,
              conf.high=1)  # add the reference
}

# for extracting N
n_included <- function(x){
      table(model.frame(x)$vitd_s) %>% as.data.frame() %>% 
            pivot_wider(names_from = Var1, values_from =  Freq)
}


# primary exposure -----

# Crude model -----
hz_vdw_crude <- survreg( Surv(fu_yr, hz) ~ vitd_s, 
                         data = bd_i, dist = "weibull")

# tidy regression output
hz_vdw_crude_tidy <- hz_vdw_crude %>% tidy_reg()

# reverse the positive and negative value
hz_vdw_crude_tidy[,2:4] <- hz_vdw_crude_tidy[,2:4]*-1

# rename
hz_vdw_crude_tidy %<>% rename_ci()

# exponential and round 
hz_vdw_crude_tidy[stats_to_convert] <- map(hz_vdw_crude_tidy[stats_to_convert], exp_round)

# add reference row
hz_vdw_crude_tidy <- hz_vdw_crude_tidy %>% ref_row()

# get N
hz_vdw_crude_n <- hz_vdw_crude %>% n_included()

# add N
hz_vdw_crude_tidy$term %>% tab1
hz_vdw_crude_tidy %<>% mutate(N = case_when(
      term == "vitd_s0_deficiency" ~ hz_vdw_crude_n$`0_deficiency`,
      term == "vitd_s1_insufficiency" ~ hz_vdw_crude_n$`1_insufficiency`,
      term == "vitd_s2_sufficiency" ~ hz_vdw_crude_n$`2_sufficiency`
))

# rename and clean
crude_w<- hz_vdw_crude_tidy #rename for output organization
rm(list = ls(pattern = "vd"))

# Partially adjusted model ----
# Model
hz_vdw_partial <-  survreg( Surv(fu_yr, hz) ~ vitd_s + sex + age_c, 
                            data = bd_i, dist = "weibull")

# tidy output
hz_vdw_partial_w <- hz_vdw_partial %>% tidy_reg()

# reverse negative and positive
hz_vdw_partial_w[,-1] <- hz_vdw_partial_w[,-1]*-1

# rename
hz_vdw_partial_w %<>% rename_ci()

# exponential and round 
hz_vdw_partial_w[stats_to_convert] <- map(hz_vdw_partial_w[stats_to_convert], exp_round)

# add reference
hz_vdw_partial_w <- hz_vdw_partial_w %>% ref_row()

# get N
hz_vdw_partial_n <- hz_vdw_partial %>% n_included()

# add N
hz_vdw_partial_w %<>% mutate(N=case_when(
      term == "vitd_s0_deficiency" ~ hz_vdw_partial_n$`0_deficiency`,
      term == "vitd_s1_insufficiency" ~ hz_vdw_partial_n$`1_insufficiency`,
      term == "vitd_s2_sufficiency" ~ hz_vdw_partial_n$`2_sufficiency`
))

# rename and clean
partial_w <- hz_vdw_partial_w
rm(list = ls(pattern = "vd"))

# Fully adjusted model -----
# model
hz_vd.weibull <- survreg( Surv(fu_yr, hz) ~ vitd_s + sex + age_c + ethnic + 
                                bmi_group + drink_freq_c + smoke_stat + imd_bd_q + regions +
                                season_c + asthma + ckd+ copd + depress + dm + ibd + ra + sle + 
                                immunosuppression, 
                          data = bd_i, dist = "weibull")
# tidy output
hz_vdw_full_w <- hz_vd.weibull %>% tidy_reg()

# reverse + / - 
hz_vdw_full_w[,-1]<- hz_vdw_full_w[,-1]*-1

# rename
hz_vdw_full_w %<>% rename_ci()

# exp and round
hz_vdw_full_w[stats_to_convert] <- map(hz_vdw_full_w[stats_to_convert], exp_round)

# add ref row
hz_vdw_full_w <- hz_vdw_full_w %>% ref_row()

# get N
hz_vdw_full_n <- hz_vd.weibull %>% n_included()

# add N
hz_vdw_full_w %<>% mutate(N=case_when(
      term == "vitd_s0_deficiency" ~ hz_vdw_full_n$`0_deficiency`,
      term == "vitd_s2_sufficiency" ~ hz_vdw_full_n$`2_sufficiency`,
      term == "vitd_s1_insufficiency" ~ hz_vdw_full_n$`1_insufficiency`
))


# rename and clean
full_w <- hz_vdw_full_w
rm(list = ls(pattern = "vd"))

# add model names
crude_w$model <- "Crude"
partial_w$model <- "Partially adjusted"
full_w$model <- "Fully adjusted"

# combine ----
fig2_data <- bind_rows(crude_w, partial_w, full_w)
fig2_data %>% write_csv("dataset_wd/fig2_data_weibull.csv")
# fig2_data <- read_csv("dataset_wd/fig2_data_weibull.csv")

# Further data management for forestplot ----

# a: vd status:
fig2_data_w <- fig2_data %>% 
      mutate(conf.high2 = format(conf.high, nsmall=2)) %>%             # 2 digits strings
      mutate(conf.low2 = format(conf.low, nsmall=2)) %>%               # 2 digits strings
      mutate(conf.rr = format(estimate, nsmall=2)) %>%                 # 2 digits strings
      mutate(ci = paste(conf.low2, conf.high2, sep = ",")) %>%        # lci, hci
      mutate(ci_l = paste0("(", ci)) %>%                              # (lci, hci
      mutate(ci_r = paste0(ci_l, ")")) %>%                            # (lci, hci)
      mutate(rr = paste(conf.rr, ci_r, sep = " ") ) %>%               # rr (lci, hci)
      mutate(HR = ifelse(estimate==1.00, "1 (Reference)",rr)) %>%     # set 1.0 into 1 (reference)
      mutate(vd_status=case_when(
            term == "vitd_s0_deficiency" ~ "Deficiency",
            term == "vitd_s1_insufficiency" ~ "Insufficiency",
            term == "vitd_s2_sufficiency" ~ "Sufficiency")) %>%     # rename for labelling
      mutate(vd_status=case_when(
            term == "vitd_s0_deficiency" ~ "Deficiency",
            term == "vitd_s1_insufficiency" ~ "Insufficiency",
            term == "vitd_s2_sufficiency" ~ "Sufficiency")) %>%     # rename for labelling
      mutate(Model = ifelse(vd_status =="Sufficiency", 
                            model, NA)) %>%  
      mutate(model_order=case_when(
            model == "Crude" ~ 1,
            model == "Partially adjusted" ~ 2,
            model == "Fully adjusted" ~3)) %>%   # group order
      mutate(exp_order=case_when(
            vd_status == "Sufficiency" ~ 1,
            vd_status == "Insufficiency" ~ 2,
            vd_status ==  "Deficiency" ~ 3)) %>%  # order within a group 
      arrange(model_order, exp_order) %>% # arrange the order by model then exposure
      dplyr::select(Model, vd_status, N, estimate, conf.low, conf.high, HR)
fig2_data_w$order <- 9:1

# Plotting: -----
# central figure 
fig2_plot_w <- fig2_data_w %>% 
      ggplot(aes(y = order, x = estimate, xmin=conf.low, xmax=conf.high, color=vd_status)) + 
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
            axis.text.y = element_blank())+
      xlab("Hazard Ratio")
fig2_plot_w

# left table
fig2_left_w <- fig2_data_w %>% 
      ggplot(aes(y = order)) + xlim(0,0.8) +
      geom_text(aes(x = 0, label = Model ),lineheight = 0.001, hjust = 0, size = 6, colour = "black") + 
      geom_text(aes(x = 0.3, label = vd_status),lineheight = 0.001, hjust = 0 ,size = 6, colour = "black") +
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
fig2_left_w

# Right table
fig2_right_w <-  fig2_data_w %>% 
      ggplot(aes(y = order)) + xlim(0,0.5) +
      geom_text(aes(x = 0, label = HR),lineheight = 0.01, hjust = 0, size = 6, colour = "black") + 
      theme_classic2( ) +
      ggtitle("  HR(95%CI)")+
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
fig2_right_w

# Combine 
Fig2_w <- grid.arrange(fig2_left_w,
                       fig2_plot_w, 
                       fig2_right_w, ncol=3)

# size: W:1800 H:600
