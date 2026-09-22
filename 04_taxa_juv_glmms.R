rm(list = ls())
#install.packages("librarian")
librarian::shelf(here, janitor, lubridate, tidyverse, ggplot2, glmmTMB, performance, DHARMa,
                 car, broom.mixed)

#### FUNCTIONS ####
gmc_data <- function(df) {
  df %>%
    # in each year, take the regional mean
    group_by(study_year, region) %>%
    mutate(
      across(
        # do this for all the variables that are treated as explanatory variables
        all_of(explanatory_vars),
        list(
          # calculate the mundlak device (VAR_region) and anomaly (VAR_dev)
          region = ~ mean(.x, na.rm = TRUE),
          dev = ~ .x - mean(.x, na.rm = TRUE)
        ),
        .names = "{.col}_{.fn}"
      )
    ) %>%
    ungroup()
}

#### READ AND MERGE DATA ####
###### READ #####
env <- read.csv(here::here("clean_data", "env_data_clean.csv"))
macro <- read.csv(here::here("clean_data", "macro_cover_clean.csv")) %>%
  mutate(study_year = as.factor(recruit_year),
         macroalgae_sqrt = sqrt(macroalgae)) %>%
  select(-habitat)

recruit <- read.csv(here::here("clean_data", "recruits_clean.csv"))
juv <- read.csv(here::here("clean_data", "juv_clean.csv"))

###### MERGE DATASETS #####

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
site_year <- juv_sub %>%
  distinct(study_year, site_name, region)

# all taxa
taxa <- recruit_sub %>%
  distinct(taxa)

# expand grid to create df with all taxa for all sites/years then merge with the juvenile df
juv_expand <- expand_grid(site_year, taxa) %>%
  left_join(
    juv_sub,
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
  left_join(macro, by = c("study_year", "site_name", "region", "cover_year", "recruit_year")) 

juv_full <- juv_pred %>%
  left_join(juv_expand, by = c("study_year", "site_name", "taxa")) 

##### SUBSET AND CENTER #####
# define predictors
explanatory_vars <- c(
  "sst_mean", "dhw_log", "kd_log", "depth", "recruit_density_log", "macroalgae_sqrt"
)
# taxa are "AGAR" "FAVI" "PORI" "SIDE" "OCTO"

# AGAR
agar_juv <- subset(juv_full, taxa == "AGAR")
agar_juv_center <- gmc_data(agar_juv)

# FAVI
favi_juv <- subset(juv_full, taxa == "FAVI")
favi_juv_center <- gmc_data(favi_juv)

# PORI
pori_juv <- subset(juv_full, taxa == "PORI")
pori_juv_center <- gmc_data(pori_juv)

# SIDE
side_juv <- subset(juv_full, taxa == "SIDE")
side_juv_center <- gmc_data(side_juv)

# OCTO
octo_juv <- subset(juv_full, taxa == "OCTO")
octo_juv_center <- gmc_data(octo_juv)

#### JUV MODELS ####
###### AGAR ######
agar_juv_mod <- nlme::lme(density_log ~ 
                            sst_mean_region + sst_mean_dev +
                            dhw_log_region + dhw_log_dev +
                            kd_log_region + kd_log_dev +
                            depth_region + depth_dev +
                            recruit_density_log_region + recruit_density_log_dev +
                            macroalgae_sqrt_region + macroalgae_sqrt_dev,
                          random = ~ 1 | region/site_name,
                          na.action = na.omit,
                          data = agar_juv_center)
summary(agar_juv_mod) # nothing sig
performance::check_model(agar_juv_mod) # all look fine
Anova(agar_juv_mod) # depth region now sig
# extract effects
agar_effects <- broom.mixed::tidy(agar_juv_mod, effects = "fixed", conf.int = TRUE) %>%
  dplyr::filter(term != "(Intercept)") %>%
  select(term, estimate, conf.low, conf.high, p.value) %>%
  mutate(taxa = "AGAR")

###### FAVI ######
favi_juv_mod <- nlme::lme(density_log ~ 
                            sst_mean_region + sst_mean_dev +
                            dhw_log_region + dhw_log_dev +
                            kd_log_region + kd_log_dev +
                            depth_region + depth_dev +
                            recruit_density_log_region + recruit_density_log_dev +
                            macroalgae_sqrt_region + macroalgae_sqrt_dev,
                          random = ~ 1 | region/site_name,
                          na.action = na.omit,
                          data = favi_juv_center)
summary(favi_juv_mod) # nothing sig
performance::check_model(favi_juv_mod) # all look fine
Anova(favi_juv_mod) # reflect summary
# extract effects
favi_effects <- broom.mixed::tidy(favi_juv_mod, effects = "fixed", conf.int = TRUE) %>%
  dplyr::filter(term != "(Intercept)") %>%
  select(term, estimate, conf.low, conf.high, p.value) %>%
  mutate(taxa = "FAVI")

###### PORI ######
pori_juv_mod <- nlme::lme(density_log ~ 
                            sst_mean_region + sst_mean_dev +
                            dhw_log_region + dhw_log_dev +
                            kd_log_region + kd_log_dev +
                            depth_region + depth_dev +
                            recruit_density_log_region + recruit_density_log_dev +
                            macroalgae_sqrt_region + macroalgae_sqrt_dev,
                          random = ~ 1 | region/site_name,
                          na.action = na.omit,
                          data = pori_juv_center)
summary(pori_juv_mod) # depth region sig
performance::check_model(pori_juv_mod) # all look fine
Anova(pori_juv_mod) # reflect summary
# extract effects
pori_effects <- broom.mixed::tidy(pori_juv_mod, effects = "fixed", conf.int = TRUE) %>%
  dplyr::filter(term != "(Intercept)") %>%
  select(term, estimate, conf.low, conf.high, p.value) %>%
  mutate(taxa = "PORI")

###### SIDE ######
side_juv_mod <- nlme::lme(density_log ~ 
                            sst_mean_region + sst_mean_dev +
                            dhw_log_region + dhw_log_dev +
                            kd_log_region + kd_log_dev +
                            depth_region + depth_dev +
                            recruit_density_log_region + recruit_density_log_dev +
                            macroalgae_sqrt_region + macroalgae_sqrt_dev,
                          random = ~ 1 | region/site_name,
                          na.action = na.omit,
                          data = side_juv_center)
summary(side_juv_mod) # no sig
performance::check_model(side_juv_mod) # all look fine
Anova(side_juv_mod) # sst region now sig
# extract effects
side_effects <- broom.mixed::tidy(side_juv_mod, effects = "fixed", conf.int = TRUE) %>%
  dplyr::filter(term != "(Intercept)") %>%
  select(term, estimate, conf.low, conf.high, p.value) %>%
  mutate(taxa = "SIDE")

###### OCTO ######
octo_juv_mod <- nlme::lme(density_log ~ 
                            sst_mean_region + sst_mean_dev +
                            dhw_log_region + dhw_log_dev +
                            kd_log_region + kd_log_dev +
                            depth_region + depth_dev +
                            recruit_density_log_region + recruit_density_log_dev +
                            macroalgae_sqrt_region + macroalgae_sqrt_dev,
                          random = ~ 1 | region/site_name,
                          na.action = na.omit,
                          data = octo_juv_center)
summary(octo_juv_mod) # depth anomaly,dhw region, sst region sig
performance::check_model(octo_juv_mod) # all look fine
Anova(octo_juv_mod) # reflect summary
# extract effects
octo_effects <- broom.mixed::tidy(octo_juv_mod, effects = "fixed", conf.int = TRUE) %>%
  dplyr::filter(term != "(Intercept)") %>%
  select(term, estimate, conf.low, conf.high, p.value) %>%
  mutate(taxa = "OCTO")


# write summary tables
juv_summary <- rbind(agar_effects, favi_effects, pori_effects, side_effects) %>%
  mutate(response = "Juveniles")

# write.csv(juv_summary, "summary_tables/juv_glmm_effects.csv", row.names = F)
