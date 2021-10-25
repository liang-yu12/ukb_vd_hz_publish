# sensitivity analysis: exclude records after 2013 (vaccination programe)
# Supplementary figure 2 vitamin D status and HZ

# data: 
source("./master_script.R")

bd_i$hz_se %<>% as.numeric()

# load additional packages
lapply(c("ggpubr","grid","gridExtra","forcats"), require, character.only=T)


## Poisson regression ---
# crude analysis: 
crude_reg_se <- glm(hz_se~ vitd_s + offset(log(fu_yr_se)), 
                    data = bd_i, family = "poisson") %>% 
   tidy(conf.int = T, conf.level = 0.95, exponentiate = T) %>% 
   filter(term != "(Intercept)") %>% 
   dplyr::select(term, estimate, conf.low, conf.high)

crude_reg_se %<>% add_row(
   term = "vitd_s2_sufficiency",
   estimate = 1,
   conf.low=1,
   conf.high=1)

crude_reg_se$model <- "Crude"


# Partially adjusted for sex and age_c

partial_reg_se <- glm(hz_se~ vitd_s + sex + age_c + offset(log(fu_yr_se)), 
                      data = bd_i, family = "poisson") %>% 
   tidy(conf.int = T, conf.level = 0.95, exponentiate = T) %>% 
   filter(term == "vitd_s0_deficiency" | term == "vitd_s1_insufficiency") %>% 
   dplyr::select(term, estimate, conf.low, conf.high)

partial_reg_se %<>% add_row(
   term = "vitd_s2_sufficiency",
   estimate = 1,
   conf.low=1,
   conf.high=1)

partial_reg_se$model <- "Partially adjusted" 


# Full adjusted for all covariates
full_reg_se <- glm(hz_se~ vitd_s + offset(log(fu_yr)) + sex + age_c + ethnic + 
                      bmi_group + drink_freq_c + smoke_stat + imd_bd_q + regions +
                      season_c + asthma + ckd+ copd + depress + dm + ibd + ra + sle + 
                      immunosuppression, 
                   data = bd_i, family = "poisson") %>% 
   tidy(conf.int = T, conf.level = 0.95, exponentiate = T) %>% 
   filter(term == "vitd_s0_deficiency" | term == "vitd_s1_insufficiency") %>% 
   dplyr::select(term, estimate, conf.low, conf.high)

full_reg_se %<>% add_row(
   term = "vitd_s2_sufficiency",
   estimate = 1,
   conf.low=1,
   conf.high=1)

full_reg_se$model <- "Fully adjusted"


# Data management of output table
supp_fig1_data <- bind_rows(crude_reg_se, partial_reg_se, full_reg_se)

supp_fig1_data %>% write_csv("dataset_wd/supp_fig1_data.csv")
# supp_fig1_data <- read_csv("dataset_wd/supp_fig1_data.csv")

# rename the data for output
supp_fig1_data<- supp_fig1_data %>%
   mutate(conf.high2 = format(conf.high, nsmall=2)) %>%             # 2 digits strings
   mutate(conf.low2 = format(conf.low, nmall=2)) %>%               # 2 digits strings
   mutate(conf.rr = format(estimate, nmall=2)) %>%                 # 2 digits strings
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
   dplyr::select(model_order, Model, vd_status, estimate, conf.low, conf.high, RR)

# deal with factors issues:
supp_fig1_data[,1:3] <- map(supp_fig1_data[,1:3], as.factor)                      # factor coercion
supp_fig1_data$model_order %<>% factor(
   labels = c("Crude", "Partially adjusted", "Fully adjusted"))    # labelled for facet

supp_fig1_data$order <- 9:1

# Step 2: main forest plot -----

supp_fig1 <- supp_fig1_data %>% 
   ggplot(aes(y = order, x = estimate, xmin=conf.low, xmax=conf.high, color=vd_status)) + 
   xlim(0.8,1.8)+
   geom_point(size=3, shape=15) + 
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
   xlab("Rate ratio")

supp_fig1

# Step 3: Plot the table on the left -----

supp_fig1_table <- supp_fig1_data %>% 
   ggplot(aes(y =order)) + xlim(0,0.8) +
   geom_text(aes(x = 0, label = Model ),lineheight = 0.001, hjust = 0, size = 6, colour = "black") + 
   geom_text(aes(x = 0.5, label = vd_status),lineheight = 0.001, hjust = 0 ,size = 6, colour = "black") + 
   theme_classic2() +
   ggtitle("   Model                     Vitamin D status ")+      # use title for header
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
supp_fig1_table

# Step 4: Plot the table on the right -----

supp_fig1_table2 <- supp_fig1_data %>% 
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
supp_fig1_table2

# Step 5: combine the three figure together

grid.arrange(supp_fig1_table, supp_fig1, supp_fig1_table2, ncol=3)
# size: W: 1500 H:500