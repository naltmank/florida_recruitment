rm(list = ls())
#install.packages("librarian")
librarian::shelf(here, janitor, lubridate, tidyverse, ggplot2, glmmTMB,  performance, DHARMa,
                 car, broom.mixed)

#### READ AND MERGE DATA ####
###### READ #####
env <- read.csv(here::here("clean_data", "env_data_clean.csv"))
macro <- read.csv(here::here("clean_data", "macro_cover_clean.csv")) %>%
  mutate(study_year = as.factor(recruit_year),
         macroalgae_sqrt = sqrt(macroalgae)) %>%
  select(-habitat)

# now including adults as predictors of juveniles
adult <- read.csv(here::here("clean_data", "lta_clean.csv")) %>%
  select(-lta)
adult_octo <- read.csv(here::here("clean_data", "adult_octo_density_clean.csv")) 
# reformat octo to be in "long" format for merge
adult_octo_wide <- adult_octo %>%
  pivot_longer(cols = OCTO,
               names_to = "taxa",
               values_to = "density")
adult_full <- rbind(adult, adult_octo_wide)

recruit <- read.csv(here::here("clean_data", "recruits_clean.csv"))

juv <- read.csv(here::here("clean_data", "juv_clean.csv"))
octo_juv <- read.csv(here::here("clean_data", "octo_juv_clean.csv")) %>%
  mutate(taxa = "OCTO")



###### MERGE DATASETS #####
adult_rec_sub <- adult_full %>%
  filter(!taxa %in% c("Other", "UNKS")) %>%
  rename(adult_density = density) %>%
  mutate(study_year = as.factor(recruit_year),
         adult_density = case_when(taxa != "OCTO" ~  adult_density*10E-04, # convert LTA to m2
                                   T ~ adult_density) # leave octocorals the same
  ) %>% 
  select(-c(cover_year, recruit_year, habitat)) %>%
  mutate(adult_density_log = log(adult_density + 1))

recruit_sub <- recruit %>%
  mutate(study_year = as.factor(recruit_year)) %>%
  filter(!taxa %in% c("Other", "UNKS", "TotalCor", "TotalSto")) %>%
  mutate(taxa = case_when(taxa == "TotalOct" ~ "OCTO",
                          T ~ taxa)) %>% 
  select(-c(abundance, total_tiles, tile_area)) %>%
  rename(recruit_density = density) %>%
  mutate(recruit_density_log = log(recruit_density + 1) )

# filter out UNKS, Other, and NA (Millepora)
juv_sub <- juv %>%
  filter(!taxa %in% c("Other", "UNKS", NA)) %>%
  mutate(study_year = as.factor(recruit_year)) %>%
  select(-c(cover_year, recruit_year, habitat)) %>%
  mutate(density_log = log(density + 1))

octo_juv_sub <- octo_juv %>%
  mutate(study_year = as.factor(recruit_year)) %>%
  select(-c(cover_year, recruit_year, habitat)) %>%
  mutate(density_log = log(density + 1))

juv_full <- rbind(juv_sub, octo_juv_sub)

# subset env to years relevant for recruits - 2015-2016
env_rec_sub <- env %>%
  filter(year %in% c(2015:2016)) %>%
  mutate(study_year = as.factor(year+1),
         chl_log = log(chl_mean + 1),
         rrs_log = log(rrs667_mean + 1),
         dhw_log = log(dhw + 1),
         kd_log = log(kd490_mean)
  ) %>%
  select(-c(year, habitat, site_code))

# not all species present in all sites for juvs - create new df where those exist as zeroes
# all site-year combos
site_year <- juv_full %>%
  distinct(study_year, site_name, region)

# all taxa
taxa <- juv_full %>%
  distinct(taxa)

# expand grid to create df with all taxa for all sites/years then merge with the juvenile df
juv_expand <- expand_grid(site_year, taxa) %>%
  left_join(
    juv_full,
    by = c("study_year", "site_name", "taxa", "region")
  ) %>%
  select(-c(region, quadrat_area, total_quadrats, abundance)) %>%
  mutate(
    density = replace_na(density, 0),
    density_log = replace_na(density_log, 0)
  )


# merge predictor datasets
juv_pred <- recruit_sub %>%
  left_join(env_rec_sub, by = c("study_year", "site_name")) %>%
  left_join(macro, by = c("study_year", "site_name", "region", "cover_year", "recruit_year")) %>%
  left_join(adult_rec_sub, by = c("study_year", "site_name", "region", "taxa"))

juv_complete <- juv_pred %>%
  left_join(juv_expand, by = c("study_year", "site_name", "taxa")) %>%
  mutate(survivorship = density/recruit_density) %>%
  filter(!is.na(survivorship)) %>%
  filter(is.finite(survivorship))

#### MODEL ####
surv_mod <- glmmTMB(survivorship ~ taxa + (1|region/site_name), family = ziGamma(link = "log"),
                    ziformula = ~taxa, data = juv_complete)
summary(surv_mod)
plot(simulateResiduals(surv_mod)) # no issues
car::Anova(surv_mod) # X2 = 288.9, P < 2.2E-16
emmeans(surv_mod, pairwise ~ taxa)

#### PLOT ####
# add connecting letters to juv_complete
juv_complete <- juv_complete %>%
  mutate(letters = case_when(
    taxa == "OCTO" ~ "d",
    taxa == "AGAR" ~ "a",
    taxa == "FAVI" ~ "b",
    taxa == "PORI" ~ "ab",
    taxa == "SIDE" ~ "c",
  ))

(surv_plot <- ggplot() +
  geom_boxplot(data = juv_complete, aes(x = taxa, y = survivorship)) +
  geom_jitter(data = juv_complete, aes(x = taxa, y = survivorship)) +
  geom_text(data = juv_complete, aes(x = taxa, y = 60, label = letters), size = 8) +
  scale_x_discrete(limits = c("AGAR", "FAVI", "PORI", "SIDE", "OCTO")) +
  scale_y_log10() +
  labs(title = "", x = "", y = "Juvenile density\n────────────\nRecruit density")+ 
  theme_classic() +
  theme(axis.title.x = element_text(size=25), 
        axis.title.y = element_text(size=23), 
        axis.text.x = element_text(color = "black", size = 24, hjust = .5, vjust = .5, face = "plain", margin = margin(r = 15)),
        axis.text.y = element_text(color = "black", size = 21, hjust = .5, vjust = .5, face = "plain", margin = margin(r = 15)),
        strip.text.x = element_text(size = 25),
        plot.title = element_text(color = "black", size = 25, hjust = 0, vjust = 0, face = "plain"),
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        legend.title = element_text(colour = "black", size = 25),
        legend.text = element_text(colour = "black", size = 24))
)

# ggsave(filename = "output/survivorship_metric_plot.png", surv_plot,
#        width = 10, height = 8, dpi = "retina")
       