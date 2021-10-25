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


# 1. Create a survival object: set origin as date of birth ----
sens_ob <- Surv(time = as.numeric(bd_i$time_in), 
                time2 = as.numeric(bd_i$time_out), 
                origin = as.numeric(bd_i$ymob), 
                event = bd_i$hz)

# because the model violate the proportional hazard assumption, 
# use stratified models

# 2. Stratify the follow-up time by vaccination time -----
bd_i$hz %<>% as.numeric # the outcome needs to be numeric or you'll get error
bd_split <- 
   survSplit(
      Surv(time = as.numeric(bd_i$time_in), 
           time2 = as.numeric(bd_i$time_out), 
           origin = as.numeric(bd_i$ymob), 
           event = bd_i$hz)~.,
      data = bd_i,
      cut = as.numeric(as.Date("2013-09-01")),
      start = "start_d",
      end = "end_d",
      episode = "vaccine"
   )

# A:Stratified model ----
## a1: stratified cox model ----
primary_split <- coxph(Surv(start_d, end_d, event) ~ vitd_s*strata(vaccine),
                       data = bd_split)

primary_split %>%  tidy_output()

## a2. Data management of output table ----
primary_split_data <- primary_split %>%  tidy_output() %>% mutate(period=ifelse(
   term == "vitd_s0_deficiency" | term == "vitd_s1_insufficiency", 1, 2))%>% 
   add_row(term = "vitd_s2_sufficiency", estimate = 1, conf.low=1, conf.high=1, period = 1) %>% 
   add_row(term = "vitd_s2_sufficiency", estimate = 1, conf.low=1, conf.high=1, period = 2) %>% 
   mutate(model="Crude") %>% 
   mutate(conf.high2 = format(conf.high, nsmall=2)) %>%             # 2 digits strings
   mutate(conf.low2 = format(conf.low, nsmall=2)) %>%               # 2 digits strings
   mutate(conf.rr = format(estimate, nsmall=2)) %>%                 # 2 digits strings
   mutate(ci = paste(conf.low2, conf.high2, sep = ",")) %>%        # lci, hci
   mutate(ci_l = paste0("(", ci)) %>%                              # (lci, hci
   mutate(ci_r = paste0(ci_l, ")")) %>%                            # (lci, hci)
   mutate(hr = paste(conf.rr, ci_r, sep = " ") ) %>%               # rr (lci, hci)
   mutate(HR = ifelse(term=="vitd_s2_sufficiency", "1 (Reference)",hr)) %>%     # set 1.0 into 1 (reference)
   mutate(vd_status=case_when(
      term == "vitd_s0_deficiency" ~ "Deficiency",
      term == "vitd_s0_deficiency:strata(vaccine)vaccine=2" ~ "Deficiency",
      term == "vitd_s1_insufficiency" ~ "Insufficiency",
      term == "vitd_s1_insufficiency:strata(vaccine)vaccine=2" ~ "Insufficiency",
      term == "vitd_s2_sufficiency" ~ "Sufficiency")) %>%     # rename for labelling
   mutate(Model = ifelse(vd_status =="Sufficiency" & period == 1, 
                         model, NA)) %>%                           # for table labelling
   arrange( period,                                    # arrange order by vd status
            match(vd_status, c("Sufficiency", "Insufficiency", "Deficiency"))) %>% 
   mutate(Period = case_when(
      vd_status =="Sufficiency" & period == 1 ~ "Before vaccination program",
      vd_status =="Sufficiency" & period == 2 ~ "After vaccination program"
   )) %>% 
   dplyr::select(Model, Period, period, vd_status, estimate, conf.low, conf.high, HR) 

primary_split_data$vd_status %<>% as.factor()
primary_split_data$vd_status %<>% factor(levels = c(c("Deficiency","Insufficiency","Sufficiency")))

primary_split_data %>% write_csv("dataset_wd/supfig4_cox_primary_split_data.csv")
# primary_split_data  <- read_csv("dataset_wd/supfig4_cox_primary_split_data.csv")

primary_split_data$order <- 6:1

## a3. main forestplot
a3_forest <- primary_split_data %>% 
   ggplot(aes(y = order, x = estimate, xmin=conf.low, xmax=conf.high, color=vd_status)) + 
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
         axis.text.y = element_blank())+
   xlab("Hazard ratio")

a3_forest

## a4. Table on the left ----
a3_table_left <- primary_split_data %>% 
   ggplot(aes(y = order)) + xlim(0,10) +
   geom_text(aes(x = 0, label = Model ),lineheight = 0.001, hjust = 0, size = 6, colour = "black") + 
   geom_text(aes(x = 3, label = Period ),lineheight = 0.001, hjust = 0, size = 6, colour = "black") + 
   geom_text(aes(x = 8, label = vd_status),lineheight = 0.001, hjust = 0 ,size = 6, colour = "black") + 
   theme_classic2() +
   ggtitle("    Model             Period                Vitamin D status")+       # use title for header
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
a3_table_left

## a4. Table on the right ----
a3_table_right <- primary_split_data %>% 
   ggplot(aes(y = order)) + xlim(0,0.5) +
   geom_text(aes(x = 0, label = HR),lineheight = 0.01, hjust = 0, size = 6, colour = "black") + 
   theme_classic2( ) +
   ggtitle("     HR(95%CI)")+
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
a3_table_right 

supp_fig4a <- grid.arrange(a3_table_left, a3_forest, a3_table_right, ncol=3)
# width: 2000  height: 400


# B: partially adjusted for sex  ----
## b1: stratified model ----
partial_split <- coxph(Surv(start_d, end_d, event) ~ vitd_s*strata(vaccine) + sex,
                       data = bd_split)
partial_split %>%  tidy_output()

## b2: data management for plotting -----
partial_split_data <- partial_split %>%  tidy_output() %>% mutate(period=ifelse(
   term == "vitd_s0_deficiency" | term == "vitd_s1_insufficiency", 1, 2))%>% 
   add_row(term = "vitd_s2_sufficiency", estimate = 1, conf.low=1, conf.high=1, period = 1) %>% 
   add_row(term = "vitd_s2_sufficiency", estimate = 1, conf.low=1, conf.high=1, period = 2) %>% 
   mutate(model="Partially adjusted") %>% 
   mutate(conf.high2 = format(conf.high, nsmall=2)) %>%             # 2 digits strings
   mutate(conf.low2 = format(conf.low, nsmall=2)) %>%               # 2 digits strings
   mutate(conf.rr = format(estimate, nsmall=2)) %>%                 # 2 digits strings
   mutate(ci = paste(conf.low2, conf.high2, sep = ",")) %>%        # lci, hci
   mutate(ci_l = paste0("(", ci)) %>%                              # (lci, hci
   mutate(ci_r = paste0(ci_l, ")")) %>%                            # (lci, hci)
   mutate(hr = paste(conf.rr, ci_r, sep = " ") ) %>%               # rr (lci, hci)
   mutate(HR = ifelse(term=="vitd_s2_sufficiency", "1 (Reference)",hr)) %>%     # set 1.0 into 1 (reference)
   mutate(vd_status=case_when(
      term == "vitd_s0_deficiency" ~ "Deficiency",
      term == "vitd_s0_deficiency:strata(vaccine)vaccine=2" ~ "Deficiency",
      term == "vitd_s1_insufficiency" ~ "Insufficiency",
      term == "vitd_s1_insufficiency:strata(vaccine)vaccine=2" ~ "Insufficiency",
      term == "vitd_s2_sufficiency" ~ "Sufficiency")) %>%     # rename for labelling
   mutate(Model = ifelse(vd_status =="Sufficiency" & period == 1, 
                         model, NA)) %>%                           # for table labelling
   arrange( period,                                    # arrange order by vd status
            match(vd_status, c("Sufficiency", "Insufficiency", "Deficiency"))) %>% 
   mutate(Period = case_when(
      vd_status =="Sufficiency" & period == 1 ~ "Before vaccination program",
      vd_status =="Sufficiency" & period == 2 ~ "After vaccination program"
   )) %>% 
   dplyr::select(Model, Period, period, vd_status, estimate, conf.low, conf.high, HR) %>% 
   filter(!is.na(vd_status)) # only keep key exposure

partial_split_data$vd_status %<>% as.factor()
partial_split_data$vd_status %<>% factor(levels = c(c("Deficiency","Insufficiency","Sufficiency")))

partial_split_data %>% write_csv("dataset_wd/supfig4_cox_partial_split_data.csv")
partial_split_data <-  read_csv("dataset_wd/supfig4_cox_partial_split_data.csv")

partial_split_data$order <- 6:1

## b3: forest plot ----

b3_forest <- partial_split_data %>% 
   ggplot(aes(y = order, x = estimate, xmin=conf.low, xmax=conf.high, color=vd_status)) + 
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
         axis.text.y = element_blank())+
   xlab("Hazard ratio")
b3_forest

## b3: left table ----
b3_table_left <- partial_split_data %>% 
   ggplot(aes(y = order)) + xlim(0,10) +
   geom_text(aes(x = 0, label = Model ),lineheight = 0.001, hjust = 0, size = 6, colour = "black") + 
   geom_text(aes(x = 3, label = Period ),lineheight = 0.001, hjust = 0, size = 6, colour = "black") + 
   geom_text(aes(x = 8, label = vd_status),lineheight = 0.001, hjust = 0 ,size = 6, colour = "black") + 
   theme_classic2() +
   ggtitle("    Model             Period                Vitamin D status")+       # use title for header
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
b3_table_left

## b3: right table ----
b3_table_right <- partial_split_data %>% 
   ggplot(aes(y = order)) + xlim(0,0.5) +
   geom_text(aes(x = 0, label = HR),lineheight = 0.01, hjust = 0, size = 6, colour = "black") + 
   theme_classic2( ) +
   ggtitle("     HR(95%CI)")+
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
b3_table_right

# combine
supp_fig4b <- grid.arrange(b3_table_left, b3_forest, b3_table_right, ncol=3)


# C: fully adjusted model -----
## c1: split model -----
full_split <- coxph(Surv(start_d, end_d, event) ~ vitd_s*strata(vaccine)  + sex + ethnic + 
                       bmi_group + drink_freq_c + smoke_stat + imd_bd_q + regions +
                       season_c + asthma + ckd+ copd + depress + dm + ibd + ra + sle + 
                       immunosuppression,
                    data = bd_split)
full_split %>%  tidy_output() 

## c2: data management  -----
full_split_data <- full_split %>%  tidy_output() %>% mutate(period=ifelse(
   term == "vitd_s0_deficiency" | term == "vitd_s1_insufficiency", 1, 2))%>% 
   add_row(term = "vitd_s2_sufficiency", estimate = 1, conf.low=1, conf.high=1, period = 1) %>% 
   add_row(term = "vitd_s2_sufficiency", estimate = 1, conf.low=1, conf.high=1, period = 2) %>% 
   mutate(model="Fully adjusted") %>% 
   mutate(conf.high2 = format(conf.high, nsmall=2)) %>%             # 2 digits strings
   mutate(conf.low2 = format(conf.low, nsmall=2)) %>%               # 2 digits strings
   mutate(conf.rr = format(estimate, nsmall=2)) %>%                 # 2 digits strings
   mutate(ci = paste(conf.low2, conf.high2, sep = ",")) %>%        # lci, hci
   mutate(ci_l = paste0("(", ci)) %>%                              # (lci, hci
   mutate(ci_r = paste0(ci_l, ")")) %>%                            # (lci, hci)
   mutate(hr = paste(conf.rr, ci_r, sep = " ") ) %>%               # rr (lci, hci)
   mutate(HR = ifelse(term=="vitd_s2_sufficiency", "1 (Reference)",hr)) %>%     # set 1.0 into 1 (reference)
   mutate(vd_status=case_when(
      term == "vitd_s0_deficiency" ~ "Deficiency",
      term == "vitd_s0_deficiency:strata(vaccine)vaccine=2" ~ "Deficiency",
      term == "vitd_s1_insufficiency" ~ "Insufficiency",
      term == "vitd_s1_insufficiency:strata(vaccine)vaccine=2" ~ "Insufficiency",
      term == "vitd_s2_sufficiency" ~ "Sufficiency")) %>%     # rename for labelling
   mutate(Model = ifelse(vd_status =="Sufficiency" & period == 1, 
                         model, NA)) %>%                           # for table labelling
   arrange( period,                                    # arrange order by vd status
            match(vd_status, c("Sufficiency", "Insufficiency", "Deficiency"))) %>% 
   mutate(Period = case_when(
      vd_status =="Sufficiency" & period == 1 ~ "Before vaccination program",
      vd_status =="Sufficiency" & period == 2 ~ "After vaccination program"
   )) %>% 
   dplyr::select(Model, Period, period, vd_status, estimate, conf.low, conf.high, HR) %>% 
   filter(!is.na(vd_status)) # only keep key exposure


full_split_data$vd_status %<>% as.factor()
full_split_data$vd_status %<>% factor(levels = c(c("Deficiency","Insufficiency","Sufficiency")))

full_split_data %>% write_csv("dataset_wd/supfig4_cox_full_split_data.csv")
# full_split_data <-  read_csv("dataset_wd/supfig4_cox_full_split_data.csv")
full_split_data$order <- 6:1

# c3: main forest plot -----
c3_forest <- full_split_data %>% 
   ggplot(aes(y = order, x = estimate, xmin=conf.low, xmax=conf.high, color=vd_status)) + 
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
         axis.text.y = element_blank())+
   xlab("Hazard ratio")
c3_forest

# c4: table on the left
c3_table_left <- full_split_data %>% 
   ggplot(aes(y = order)) + xlim(0,10) +
   geom_text(aes(x = 0, label = Model ),lineheight = 0.001, hjust = 0, size = 6, colour = "black") + 
   geom_text(aes(x = 3, label = Period ),lineheight = 0.001, hjust = 0, size = 6, colour = "black") + 
   geom_text(aes(x = 8, label = vd_status),lineheight = 0.001, hjust = 0 ,size = 6, colour = "black") + 
   theme_classic2() +
   ggtitle("    Model             Period                Vitamin D status")+     # use title for header
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
c3_table_left

# c3: table on the right-----
c3_table_right <- full_split_data %>% 
   ggplot(aes(y = order)) + xlim(0,0.5) +
   geom_text(aes(x = 0, label = HR),lineheight = 0.01, hjust = 0, size = 6, colour = "black") + 
   theme_classic2( ) +
   ggtitle("     HR(95%CI)")+
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
c3_table_right 

supp_fig4c <- grid.arrange(c3_table_left, c3_forest, c3_table_right, ncol=3)



# 3. Combine a, b, c
ggarrange(supp_fig4a, supp_fig4b, supp_fig4c,
          labels = c("a", "b", "c"),
          ncol = 1, nrow = 3,
          align = "v",
          widths = 1)

# output size: W:2000 H:1600



