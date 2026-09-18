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
adult <- read.csv(here::here("clean_data", "lta_clean.csv")) %>%
  select(-lta)
adult_octo <- read.csv(here::here("clean_data", "adult_octo_density_clean.csv")) 
# reformat octo to be in "long" format for merge
adult_octo_wide <- adult_octo %>%
  pivot_longer(cols = OCTO,
               names_to = "taxa",
               values_to = "density")
adult_full <- rbind(adult, adult_octo_wide)

env <- read.csv(here::here("clean_data", "env_data_clean.csv"))
macro <- read.csv(here::here("clean_data", "macro_cover_clean.csv")) %>%
  mutate(study_year = as.factor(recruit_year)) %>%
  select(-habitat)

recruit <- read.csv(here::here("clean_data", "recruits_clean.csv"))

###### MERGE RECRUITS #####
# filter out UNKS and Other
adult_rec_sub <- adult_full %>%
  filter(!taxa %in% c("Other", "UNKS")) %>%
  rename(adult_density = density) %>%
  mutate(study_year = as.factor(recruit_year)) %>%
  select(-c(cover_year, recruit_year, juv_year))

# subset env to years relevant for recruits - 2015-2017
env_rec_sub <- env %>%
  filter(year %in% c(2015:2017)) %>%
  mutate(study_year = as.factor(year+1)) %>%
  select(-c(year, habitat, site_code))

# merge predictor datasets
recruit_pred <- adult_rec_sub %>%
  left_join(env_rec_sub, by = c("study_year", "site_name", "region")) %>%
  left_join(macro, by = c("study_year", "site_name", "region")) 


recruit_sub <- recruit %>%
  mutate(study_year = as.factor(recruit_year)) %>%
  filter(!taxa %in% c("Other", "UNKS", "TotalCor", "TotalSto")) %>%
  mutate(taxa = case_when(taxa == "TotalOct" ~ "OCTO",
                          T ~ taxa)) %>%
  select(-c(abundance, total_tiles, tile_area))

recruit_full <- recruit_sub %>%
  left_join(recruit_pred, by = c("study_year", "site_name", "taxa",
                                 "recruit_year", "cover_year", "juv_year")) %>%
  mutate(density_log = log(density + 1) )

##### SUBSET AND CENTER #####
# define predictors
explanatory_vars <- c(
  "sst_mean", "dhw", "chl_mean", "rrs667_mean", "depth", "adult_density", "macroalgae"
)
# taxa are "AGAR" "FAVI" "PORI" "SIDE" "OCTO"

# AGAR
agar_rec <- subset(recruit_full, taxa == "AGAR")
agar_rec_center <- gmc_data(agar_rec)

# FAVI
favi_rec <- subset(recruit_full, taxa == "FAVI")
favi_rec_center <- gmc_data(favi_rec)

# PORI
pori_rec <- subset(recruit_full, taxa == "PORI")
pori_rec_center <- gmc_data(pori_rec)

# SIDE
side_rec <- subset(recruit_full, taxa == "SIDE")
side_rec_center <- gmc_data(side_rec)

# OCTO
octo_rec <- subset(recruit_full, taxa == "OCTO")
octo_rec_center <- gmc_data(octo_rec)

#### RECRUIT MODELS ####
###### AGAR ######
agar_rec_mod <- nlme::lme(density_log ~ 
                            sst_mean_region + sst_mean_dev +
                            dhw_region + dhw_dev +
                            chl_mean_region + chl_mean_dev +
                            rrs667_mean_region + rrs667_mean_dev +
                            depth_region + depth_dev +
                            adult_density_region + adult_density_dev +
                            macroalgae_region + macroalgae_dev,
                          random = ~ 1 | region/site_name,
                          na.action = na.omit,
                          data = agar_rec_center)
summary(agar_rec_mod) # adult density anomaly, sst_mean anomaly sig
performance::check_model(agar_rec_mod) # all look fine
Anova(agar_rec_mod) # reflect summary
# extract effects
agar_effects <- broom.mixed::tidy(agar_rec_mod, effects = "fixed", conf.int = TRUE) %>%
  dplyr::filter(term != "(Intercept)") %>%
  select(term, estimate, conf.low, conf.high, p.value) %>%
  mutate(taxa = "AGAR")

###### FAVI ######
favi_rec_mod <- nlme::lme(density_log ~ 
                          sst_mean_region + sst_mean_dev +
                          dhw_region + dhw_dev +
                          chl_mean_region + chl_mean_dev +
                          rrs667_mean_region + rrs667_mean_dev +
                          depth_region + depth_dev +
                          adult_density_region + adult_density_dev +
                          macroalgae_region + macroalgae_dev,
                          random = ~ 1 | region/site_name,
                        na.action = na.omit,
                        data = favi_rec_center)
summary(favi_rec_mod) # rrs anomaly
performance::check_model(favi_rec_mod) # all look fine
Anova(favi_rec_mod) # reflect summary
# extract effects
favi_effects <- broom.mixed::tidy(favi_rec_mod, effects = "fixed", conf.int = TRUE) %>%
  dplyr::filter(term != "(Intercept)") %>%
  select(term, estimate, conf.low, conf.high, p.value) %>%
  mutate(taxa = "FAVI")

###### pori ######
pori_rec_mod <- nlme::lme(density_log ~ 
                            sst_mean_region + sst_mean_dev +
                            dhw_region + dhw_dev +
                            chl_mean_region + chl_mean_dev +
                            rrs667_mean_region + rrs667_mean_dev +
                            depth_region + depth_dev +
                            adult_density_region + adult_density_dev +
                            macroalgae_region + macroalgae_dev,
                          random = ~ 1 | region/site_name,
                          na.action = na.omit,
                          data = pori_rec_center)
summary(pori_rec_mod) # adult density anomaly, sst_mean region sig
performance::check_model(pori_rec_mod) # all look fine
Anova(pori_rec_mod) # reflect summary
# extract effects
pori_effects <- broom.mixed::tidy(pori_rec_mod, effects = "fixed", conf.int = TRUE) %>%
  dplyr::filter(term != "(Intercept)") %>%
  select(term, estimate, conf.low, conf.high, p.value) %>%
  mutate(taxa = "PORI")

###### side ######
side_rec_mod <- nlme::lme(density_log ~ 
                            sst_mean_region + sst_mean_dev +
                            dhw_region + dhw_dev +
                            chl_mean_region + chl_mean_dev +
                            rrs667_mean_region + rrs667_mean_dev +
                            depth_region + depth_dev +
                            adult_density_region + adult_density_dev +
                            macroalgae_region + macroalgae_dev,
                          random = ~ 1 | region/site_name,
                          na.action = na.omit,
                          data = side_rec_center)
summary(side_rec_mod) # depth anomaly and region, rrs667 anomaly, dhw region, sst region sig
performance::check_model(side_rec_mod) # all look fine
Anova(side_rec_mod) # reflect summary
# extract effects
side_effects <- broom.mixed::tidy(side_rec_mod, effects = "fixed", conf.int = TRUE) %>%
  dplyr::filter(term != "(Intercept)") %>%
  select(term, estimate, conf.low, conf.high, p.value) %>%
  mutate(taxa = "SIDE")

###### octo ######
octo_rec_mod <- nlme::lme(density_log ~ 
                            sst_mean_region + sst_mean_dev +
                            dhw_region + dhw_dev +
                            chl_mean_region + chl_mean_dev +
                            rrs667_mean_region + rrs667_mean_dev +
                            depth_region + depth_dev +
                            adult_density_region + adult_density_dev +
                            macroalgae_region + macroalgae_dev,
                          random = ~ 1 | region/site_name,
                          na.action = na.omit,
                          data = octo_rec_center)
summary(octo_rec_mod) # depth anomaly, chl anomaly, dhw region, sst anomaly sig
performance::check_model(octo_rec_mod) # all look fine
Anova(octo_rec_mod) # reflect summary
# extract effects
octo_effects <- broom.mixed::tidy(octo_rec_mod, effects = "fixed", conf.int = TRUE) %>%
  dplyr::filter(term != "(Intercept)") %>%
  select(term, estimate, conf.low, conf.high, p.value) %>%
  mutate(taxa = "OCTO")
