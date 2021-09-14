# Script for visualizing the regression model
# Fig 3: vitamin D supplement/prescription and herpes zoster

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
vd_supp_crude <- glm(hz~all_supp + offset(log(fu_yr)), 
                     data = bd_i, family = "poisson") %>% 
   tidy(conf.int = T, conf.level = 0.95, exponentiate = T) %>% 
   filter(term != "(Intercept)") %>% 
   dplyr::select(term, estimate, conf.low, conf.high)

vd_supp_crude %<>% add_row(
   term = "No_supp",
   estimate = 1,
   conf.low=1,
   conf.high=1) # add reference

vd_supp_crude$model <- "Crude"

# partial ----
vd_supp_partial <- glm(hz~all_supp + offset(log(fu_yr)) + sex + age_c, 
                       data = bd_i, family = "poisson") %>% 
   tidy(conf.int = T, conf.level = 0.95, exponentiate = T) %>% 
   filter(term == "all_supp1_vitD and mineral supplement") %>% 
   dplyr::select(term, estimate, conf.low, conf.high)

vd_supp_partial %<>% add_row(
   term = "No_supp",
   estimate = 1,
   conf.low=1,
   conf.high=1) # add reference

vd_supp_partial$model <- "Partially adjusted" # marked the model

# Full adjusted model ----
vd_supp_full <- glm(hz~ all_supp + offset(log(fu_yr)) + sex + age_c + ethnic + 
                       bmi_group + drink_freq_c + smoke_stat + imd_bd_q + regions +
                       season_c + asthma + ckd+ copd + depress + dm + ibd + ra + sle + immunosuppression, 
                    data = bd_i, family = "poisson") %>% 
   tidy(conf.int = T, conf.level = 0.95, exponentiate = T, n_digits = 2) %>% 
   filter(term == "all_supp1_vitD and mineral supplement") %>% 
   dplyr::select(term, estimate, conf.low, conf.high)

vd_supp_full %<>% add_row(
   term = "No_supp",
   estimate = 1,
   conf.low=1,
   conf.high=1) # add reference

vd_supp_full$model <- "Fully adjusted" # marked the model


fig3a_data <- bind_rows(vd_supp_crude, vd_supp_partial,vd_supp_full)



fig3a_data %>% write_csv("dataset_wd/fig3a_data.csv")
# fig3a_data  <-  read_csv("dataset_wd/fig3a_data.csv")

# Fig3a: data management for plotting: ----
fig3a_data<- fig3a_data %>%
   mutate(conf.high2 = format(conf.high, nsmall=2)) %>%             # 2 digits strings
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
   dplyr::select(model_order, Model, vd_intake, estimate, conf.low, conf.high, RR)

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

# Fig 3b Vitamin D prescription and HZ ----

# crude ----
vd_drug_crude <- glm(hz~vd_prescription + offset(log(fu_yr)), 
                     data = bd_i, family = "poisson") %>% 
   tidy(conf.int = T, conf.level = 0.95, exponentiate = T) %>% 
   filter(term != "(Intercept)") %>% 
   dplyr::select(term, estimate, conf.low, conf.high)

vd_drug_crude %<>% add_row(
   term = "No_supp",
   estimate = 1,
   conf.low=1,
   conf.high=1) # add reference

vd_drug_crude$model <- "Crude"

# Partially adjusted for sex and age
vd_drug_partial <- glm(hz~vd_prescription + offset(log(fu_yr)) + sex + age_c, 
                       data = bd_i, family = "poisson") %>% 
   tidy(conf.int = T, conf.level = 0.95, exponentiate = T) %>% 
   filter(term == "vd_prescriptionHad vitamin D prescription") %>% 
   dplyr::select(term, estimate, conf.low, conf.high)

vd_drug_partial %<>% add_row(
   term = "No_supp",
   estimate = 1,
   conf.low=1,
   conf.high=1) # add reference

vd_drug_partial$model <- "Partially adjusted"


# Fully adjusted 
vd_drug_full <- glm(hz~ vd_prescription + offset(log(fu_yr)) + sex + age_c + ethnic + 
                       bmi_group + drink_freq_c + smoke_stat + imd_bd_q + regions +
                       season_c + asthma + ckd+ copd + depress + dm + ibd + ra + sle + immunosuppression, 
                    data = bd_i, family = "poisson") %>% 
   tidy(conf.int = T, conf.level = 0.95, exponentiate = T, n_digits = 2) %>% 
   filter(term == "vd_prescriptionHad vitamin D prescription") %>% 
   dplyr::select(term, estimate, conf.low, conf.high)

vd_drug_full %<>% add_row(
   term = "No_supp",
   estimate = 1,
   conf.low=1,
   conf.high=1) # add reference

vd_drug_full$model <- "Fully adjusted"


fig3b_data <- bind_rows(vd_drug_crude, vd_drug_partial, vd_drug_full) # combine output tables

fig3b_data %>% write_csv("dataset_wd/fig3b_data.csv")
# fig3b_data <- read_csv("dataset_wd/fig3b_data.csv")
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
   dplyr::select(model_order, Model, vd_intake, estimate, conf.low, conf.high, RR)


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

# Combine two Figures using ggpubr ---

ggarrange(fig3a_combine, fig3b_combine,
          labels = c("a", "b"),
          ncol = 1, nrow = 2,
          align = "v",
          widths = 1)

# output size: W: 1500 H:1000