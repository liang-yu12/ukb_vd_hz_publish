# Script for visualizing regression model 
# Figure 2: the association between vitamin D status and herpes zoster

# data: 
source("./master_script.R")

# load additional packages
lapply(c("tidyverse","magrittr","ggpubr","grid","gridExtra","forcats"), require, character.only=T)

# option
options(digits = 2, scipen = 999)


## Poisson regression ---
# crude analysis: 
crude_reg <- glm(hz~ vitd_s + offset(log(fu_yr)), 
                 data = bd_i, family = "poisson") %>% 
   tidy(conf.int = T, conf.level = 0.95, exponentiate = T) %>%    # Use tidyr from broom package to organise the output
   filter(term != "(Intercept)") %>%                              # Only keep the coefficient of the exposures
   dplyr::select(term, estimate, conf.low, conf.high)

crude_reg %<>% add_row(
   term = "vitd_s2_sufficiency",
   estimate = 1,
   conf.low=1,
   conf.high=1)                                                   # Add the reference (sufficient vitamin D)

crude_reg$model <- "Crude"                                        # Label the model


# Partially adjusted for sex and age_c

partial_reg <- glm(hz~ vitd_s + sex + age_c + offset(log(fu_yr)), 
                   data = bd_i, family = "poisson") %>% 
   tidy(conf.int = T, conf.level = 0.95, exponentiate = T) %>%                # tidy the output
   filter(term == "vitd_s0_deficiency" | term == "vitd_s1_insufficiency") %>% # Only keep the coefficient of the exposures
   dplyr::select(term, estimate, conf.low, conf.high)

partial_reg %<>% add_row(
   term = "vitd_s2_sufficiency",
   estimate = 1,
   conf.low=1,
   conf.high=1)                                                   # Add the reference (sufficient vitamin D)

partial_reg$model <- "Partially adjusted"                         # Label the model


# Full adjusted for all covariates
full_reg <- glm(hz~ vitd_s + offset(log(fu_yr)) + sex + age_c + ethnic + 
                   bmi_group + drink_freq_c + smoke_stat + imd_bd_q + regions +
                   season_c + asthma + ckd+ copd + depress + dm + ibd + ra + sle + 
                   immunosuppression, 
                data = bd_i, family = "poisson") %>% 
   tidy(conf.int = T, conf.level = 0.95, exponentiate = T) %>%    # tidy the output
   filter(term == "vitd_s0_deficiency" | term == "vitd_s1_insufficiency") %>% # Only keep the coefficient of the exposures
   dplyr::select(term, estimate, conf.low, conf.high)

full_reg %<>% add_row(
   term = "vitd_s2_sufficiency",
   estimate = 1,
   conf.low=1,
   conf.high=1)                                                   # Add the reference (sufficient vitamin D)

full_reg$model <- "Fully adjusted"                                # Label the model



fig2_data <- bind_rows(crude_reg, partial_reg, full_reg)          # Combine the outputs into a tables

fig2_data %>% write_csv("dataset_wd/fig2_data.csv")               # Save the file for plotting
# fig2_data <- read_csv("dataset_wd/fig2_data.csv")


# Forest plot using ggplot2 -----

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
   dplyr::select(model_order, Model, vd_status, estimate, conf.low, conf.high, RR)


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
   theme_classic2() + 
   ggtitle("   Model                 Vitamin D status")+   
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

# size: W:1500 H:500
