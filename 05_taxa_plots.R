rm(list = ls())
#install.packages("librarian")
librarian::shelf(here, tidyverse, ggplot2, ggpubr)

#### PREP DATA ####
# read tables
rec_table <- read.csv(here::here("summary_tables", "recruit_glmm_effects.csv"))
juv_table <- read.csv(here::here("summary_tables", "juv_glmm_effects.csv"))

# filter to anomalies, add significance columns
rec_sub <- rec_table %>%
  filter(str_detect(term, "_dev"),
         taxa != "OCTO") %>%
  mutate(significance = case_when(p.value < 0.05 ~ "P < 0.05",
                                  T ~ "P > 0.05"))

juv_sub <- juv_table %>%
  filter(str_detect(term, "_dev"),
         taxa != "OCTO") %>%
  mutate(significance = case_when(p.value < 0.05 ~ "P < 0.05",
                                  T ~ "P > 0.05"))

#### PLOT ####
okabe_ito <- c(
  "adult_density_log_dev" = "#000080", 
  "recruit_density_log_dev" = "#33a8c7", 
  "macroalgae_sqrt_dev" = "#a0e426",  
  "depth_dev" = "yellow4",  
  "sst_mean_dev" = "#ffab00",  
  "dhw_log_dev" = "#f050ae",  
  "kd_log_dev" = "#d883ff"
)

legend_order <- c(
  "adult_density_log_dev",
  "recruit_density_log_dev",
  "macroalgae_sqrt_dev",
  "depth_dev",
  "sst_mean_dev",
  "dhw_log_dev",
  "kd_log_dev"
)

(rec_plot <- 
   ggplot() +
   geom_point(data = rec_sub, aes(x = estimate, y = taxa, colour = factor(term, levels = legend_order), shape = significance),
              position = position_dodge(0.5), size = 4) +
   geom_errorbar(data = rec_sub, aes(xmin = conf.low, xmax = conf.high,
                                     y = taxa, colour = factor(term, levels = legend_order),
                                     linetype = significance),
                 width = 0.1, position = position_dodge(0.5)) +
   scale_y_discrete(limits = c("SIDE", "PORI", "FAVI", "AGAR")) +
   geom_vline(xintercept = 0, linetype = 2) +
   scale_shape_manual(values = c("P < 0.05" = 19,
                                 "P > 0.05" = 1)) + 
   scale_colour_manual(values = okabe_ito,
                       labels = c(
                         "adult_density_log_dev" = "Adult corals",  
                         "macroalgae_sqrt_dev" = "Macroalgal cover",  
                         "depth_dev" = "Depth",  
                         "sst_mean_dev" = "Mean SST",  
                         "dhw_log_dev" = "DHW",  
                         "kd_log_dev" = "Turbidity"
                       )
   ) +
   labs(title = "a.", x = "Effect estimate (± 95% CI)", y = "Taxa", colour = "Predictor",
        shape = "Significance", linetype = "Significance")+ 
    theme_classic() +
    theme(axis.title.x = element_text(size=25), 
          axis.title.y = element_text(size=25), 
          axis.text.x = element_text(color = "black", size = 24, hjust = .5, vjust = .5, face = "plain", margin = margin(r = 15)),
          axis.text.y = element_text(color = "black", size = 24, hjust = .5, vjust = .5, face = "plain", margin = margin(r = 15)),
          strip.text.x = element_text(size = 25),
          plot.title = element_text(color = "black", size = 25, hjust = 0, vjust = 0, face = "plain"),
          panel.grid.major = element_blank(),  #remove major-grid labels
          panel.grid.minor = element_blank(),  #remove minor-grid labels
          legend.title = element_text(colour = "black", size = 25),
          legend.text = element_text(colour = "black", size = 24))
   )



(juv_plot <- 
    ggplot() +
    geom_point(data = juv_sub, aes(x = estimate, y = taxa, colour = factor(term, levels = legend_order),
                                   shape = significance),
               position = position_dodge(0.5), size = 4) +
    geom_errorbar(data = juv_sub, aes(xmin = conf.low, xmax = conf.high,
                                      y = taxa, colour = factor(term, levels = legend_order),
                                      linetype = significance),
                  width = 0.1, position = position_dodge(0.5)) +
    scale_y_discrete(limits = c("SIDE", "PORI", "FAVI", "AGAR")) +
    geom_vline(xintercept = 0, linetype = 2) +
    scale_shape_manual(values = c("P < 0.05" = 19,
                                  "P > 0.05" = 1)) + 
    scale_linetype_manual(values = c("P < 0.05" = 1,
                                    "P > 0.05" = 2)) +
    scale_colour_manual(values = okabe_ito,
                        breaks = c(
                          "adult_density_log_dev",
                          "recruit_density_log_dev",
                          "macroalgae_sqrt_dev",
                          "depth_dev",
                          "sst_mean_dev",
                          "dhw_log_dev",
                          "kd_log_dev"
                        ),
                        labels = c(
                          "adult_density_log_dev" = "Adult corals", 
                          "recruit_density_log_dev" = "Coral recruits",  
                          "macroalgae_sqrt_dev" = "Macroalgal cover",  
                          "depth_dev" = "Depth",  
                          "sst_mean_dev" = "Mean SST",  
                          "dhw_log_dev" = "DHW",  
                          "kd_log_dev" = "Turbidity"
                        )
    ) + 
    labs(title = "b.", x = "Effect estimate (± 95% CI)", y = "Taxa", colour = "Predictor",
         shape = "Significance", linetype = "Significance")+ 
    theme_classic() +
    theme(axis.title.x = element_text(size=25), 
          axis.title.y = element_text(size=25), 
          axis.text.x = element_text(color = "black", size = 24, hjust = .5, vjust = .5, face = "plain", margin = margin(r = 15)),
          axis.text.y = element_text(color = "black", size = 24, hjust = .5, vjust = .5, face = "plain", margin = margin(r = 15)),
          strip.text.x = element_text(size = 25),
          plot.title = element_text(color = "black", size = 25, hjust = 0, vjust = 0, face = "plain"),
          panel.grid.major = element_blank(),  #remove major-grid labels
          panel.grid.minor = element_blank(),  #remove minor-grid labels
          legend.title = element_text(colour = "black", size = 25),
          legend.text = element_text(colour = "black", size = 24))
)

# combine
(taxa_panel <- ggarrange(rec_plot, juv_plot, nrow = 1, ncol = 2, widths = c(1,1),
                         common.legend = T, legend.grob = get_legend(juv_plot),
                         legend = "right")
)

# ggsave(filename = "output/taxa_panel.png", taxa_panel,
#        width = 16, height = 12, dpi = "retina")

