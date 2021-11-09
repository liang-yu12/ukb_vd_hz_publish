## Sensitivity analysis using Weibull regression ----

# data: 
source("./master_script.R")

# load additional packages
lapply(c("ggpubr","grid","gridExtra","forcats", "survival", "survminer"), require, character.only=T)

options(digits = 3, scipen = 999)

stats_to_convert<- c("estimate",  "conf.low",  "conf.high")

# Data management of the regression model : 
# 1. weibull regression -> save the model
# 2. use tidy() to organize the regression model
# 2.5 fixed the reverse +-
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


# Original full model ------
# Fully adjusted model -----
# model
hz_vd.weibull <- survreg( Surv(fu_yr, hz) ~ vitd_s + sex + age_c + ethnic + 
                                bmi_group + drink_freq_c + smoke_stat + imd_bd_q + regions +
                                season_c + asthma + ckd+ copd + depress + dm + ibd + ra + sle + 
                                immunosuppression, 
                          data = bd_i, dist = "weibull")
# tidy output
hz_vdw_full_w <- hz_vd.weibull %>% tidy_reg()

# fixed +-
hz_vdw_full_w[,-1] <- hz_vdw_full_w[,-1]*-1

# rename
hz_vdw_full_w %<>% rename_ci

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
original_full_w <- hz_vdw_full_w
rm(list = ls(pattern = "vd"))


# sensitivity model: without self-reported health disease -----
full_reg_se1.w<- survreg(Surv(fu_yr, hz) ~ vitd_s + sex + age_c + ethnic + 
                               bmi_group + drink_freq_c + smoke_stat + imd_bd_q + regions +
                               season_c + asthma_se + ckd_se+ copd_se + depress_se + dm_se + 
                               ibd_se + ra_se + sle_se + 
                               immunosuppression, 
                         data = bd_i, dist = "weibull")
# tidy
full_se1 <- full_reg_se1.w %>% tidy_reg()

# fixed +-
full_se1[,-1] <- full_se1[,-1]*-1

# rename
full_se1 %<>% rename_ci

# exp round
full_se1[stats_to_convert] <- map(full_se1[stats_to_convert], exp_round)

# add ref
full_se1 %<>% ref_row()

# get N
full_se1_n <- full_reg_se1.w %>% n_included()

# add N
full_se1 %<>% mutate(N=case_when(
      term == "vitd_s0_deficiency" ~ full_se1_n$`0_deficiency`,
      term == "vitd_s2_sufficiency" ~ full_se1_n$`2_sufficiency`,
      term == "vitd_s1_insufficiency" ~ full_se1_n$`1_insufficiency` 
))

# rename and clean
se_1_full_w <- full_se1
rm(list = ls(pattern = "se1"))


# SE2: # Only include people with high-dose oral steroid, with self-reported disease----

full_reg_se2.w <- survreg( Surv(fu_yr, hz) ~ vitd_s + sex + age_c + ethnic + 
                                 bmi_group + drink_freq_c + smoke_stat + imd_bd_q + regions +
                                 season_c + asthma + ckd+ copd + depress + dm + ibd + ra + sle + 
                                 immunosuppression_2, 
                           data = bd_i, dist = "weibull")

# tidy
full_reg_se2 <- full_reg_se2.w %>% tidy_reg

# fixed +-
full_reg_se2[,-1] <- full_reg_se2[,-1]*-1

# rename
full_reg_se2 %<>% rename_ci

# exp_round
full_reg_se2[stats_to_convert] <- map(full_reg_se2[stats_to_convert], exp_round)

# add ref
full_reg_se2 %<>% ref_row()

# get N
full_reg_se2_n <- full_reg_se2.w %>% n_included()

# add N 
full_reg_se2 %<>% mutate(N=case_when(
      term == "vitd_s0_deficiency" ~ full_reg_se2_n$`0_deficiency`,
      term == "vitd_s2_sufficiency" ~ full_reg_se2_n$`2_sufficiency`,
      term == "vitd_s1_insufficiency" ~ full_reg_se2_n$`1_insufficiency` 
))

# clean and rename 
se_2_full_w <- full_reg_se2
rm(list = ls(pattern = "se2"))


# SE3 Without self-reported outcome, and only high dose oral steroid -----
full_reg_se3.w<- survreg(Surv(fu_yr, hz) ~ vitd_s + sex + age_c + ethnic + 
                               bmi_group + drink_freq_c + smoke_stat + imd_bd_q + regions +
                               season_c + asthma_se + ckd_se+ copd_se + depress_se + dm_se + 
                               ibd_se + ra_se + sle_se + 
                               immunosuppression_2, 
                         data = bd_i, dist = "weibull")

# tidy output
full_se3 <- full_reg_se3.w %>% tidy_reg()

# fixed +-
full_se3[,-1] <- full_se3[,-1]*-1

# rename
full_se3 %<>% rename_ci

# exp round
full_se3[stats_to_convert] <- map(full_se3[stats_to_convert], exp_round)

# add ref
full_se3 %<>% ref_row()

# get M
full_se3_n <- full_reg_se3.w %>% n_included()

# add N
full_se3 %<>% mutate(N=case_when(
      term == "vitd_s0_deficiency" ~ full_se3_n$`0_deficiency`,
      term == "vitd_s2_sufficiency" ~ full_se3_n$`2_sufficiency`,
      term == "vitd_s1_insufficiency" ~ full_se3_n$`1_insufficiency` 
))

# rename and clean 
se_3_full_w <- full_se3
rm(list = ls(pattern = "se3"))


# mark and combine ------
original_full_w$model <- "Original model"
se_1_full_w$model <- "Sensitivity analysis model 1"
se_2_full_w$model <- "Sensitivity analysis model 2"
se_3_full_w$model <- "Sensitivity analysis model 3"

supp_fig3_data_w <- bind_rows(original_full_w,se_1_full_w,se_2_full_w,se_3_full_w)

supp_fig3_data_w %>% write_csv("dataset_wd/supp_fig3_data_weibull.csv")
# supp_fig3_data_w <- read_csv("dataset_wd/supp_fig3_data_weibull.csv")


# Data manipulation for plotting -------

supp_fig3_data  <- supp_fig3_data_w %>%
      mutate(conf.high2 = format(conf.high, nsmall=2)) %>%             # 2 digits strings
      mutate(conf.low2 = format(conf.low, nsmall=2)) %>%               # 2 digits strings
      mutate(conf.rr = format(estimate, nsmall=2)) %>%                 # 2 digits strings
      mutate(ci = paste(conf.low2, conf.high2, sep = ",")) %>%        # lci, hci
      mutate(ci_l = paste0("(", ci)) %>%                              # (lci, hci
      mutate(ci_r = paste0(ci_l, ")")) %>%                            # (lci, hci)
      mutate(rr = paste(conf.rr, ci_r, sep = " ") ) %>%               # rr (lci, hci)
      mutate(RR = ifelse(term == "vitd_s2_sufficiency", "1 (Reference)",rr)) %>%     # set 1.0 into 1 (reference)
      mutate(vd_status=case_when(
            term == "vitd_s0_deficiency" ~ "Deficiency",
            term == "vitd_s1_insufficiency" ~ "Insufficiency",
            term == "vitd_s2_sufficiency" ~ "Sufficiency")) %>%          # rename for labelling
      mutate(model_order=case_when(
            model == "Original model" ~ 1,
            model == "Sensitivity analysis model 1" ~ 2,
            model == "Sensitivity analysis model 2" ~ 3,
            model == "Sensitivity analysis model 3" ~ 4)) %>%            # for facet order
      mutate(Model = ifelse(vd_status =="Sufficiency", 
                            model, NA)) %>%                           # for table labelling
      arrange( model_order,                                           # arrange order by vd status
               match(vd_status, c("Sufficiency", "Insufficiency", "Deficiency"))) %>% 
      dplyr::select(model_order, Model, vd_status,N , estimate, conf.low, conf.high, RR)

# deal with factors issues:
supp_fig3_data[,1:3] <- map(supp_fig3_data[,1:3], as.factor)       # factor coercion
supp_fig3_data$model_order %<>% factor(
      labels = c("Original model", "Sensitivity analysis model 1",    # labelled for facet
                 "Sensitivity analysis model 2", "Sensitivity analysis model 3"))    

supp_fig3_data$order <- 12:1


# Step 2: main forest plot -----
supp_fig3 <- supp_fig3_data %>% 
      ggplot(aes(y =order, x = estimate, xmin=conf.low, xmax=conf.high, 
                 color=vd_status)) + 
      xlim(0.5,1.5)+
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
            axis.text.y = element_blank())+
      xlab("Hazard ratio")

supp_fig3


# Step 3: Plot the table on the left -----

supp_fig3_table <- supp_fig3_data %>% 
      ggplot(aes(y = order)) + xlim(0,2) +
      geom_text(aes(x = 0, label = Model ),lineheight = 0.001, hjust = 0, size = 6, colour = "black") + 
      geom_text(aes(x = 1.2, label = vd_status),lineheight = 0.001, hjust = 0 ,size = 6, colour = "black") + 
      geom_text(aes(x = 1.8, label = N),lineheight = 0.001, hjust = 0 ,size = 6, colour = "black") + 
      theme_classic2() +
      ggtitle("   Model                        Vitamin D status       N")+    # use title for header
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
supp_fig3_table

# Step 4: Plot the table on the right -----

supp_fig3_table2 <- supp_fig3_data %>% 
      ggplot(aes(y = order)) + xlim(0,0.5) +
      geom_text(aes(x = 0, label = RR),lineheight = 0.01, hjust = 0, size = 6, colour = "black") + 
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

supp_fig3_table2 

# Step 5: combine the three figure together

grid.arrange(supp_fig3_table, supp_fig3, supp_fig3_table2, ncol=3)

# W:2000 H:600

