rm(list = ls())
#install.packages("librarian")
librarian::shelf(here, janitor, lubridate, tidyverse, ggplot2, nlme, performance, piecewiseSEM)


#### FUNCTIONS ####
gmc_data <- function(df, # dataframe to pass the function to
                     timevar = study_year, # for if sites were surveyed more than one year
                     spatial_group = region # the spatial grouping for setting regional means
                     ) {  
  df %>%
    # in each year, take the regional mean
    group_by({{ timevar }}, {{ spatial_group }}) %>%
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

#### READ AND COMBINE DATA ####
adult <- read.csv(here::here("clean_data", "lta_clean.csv")) %>%
  select(-lta)

env <- read.csv(here::here("clean_data", "env_data_clean.csv"))
macro <- read.csv(here::here("clean_data", "macro_cover_clean.csv")) %>%
  mutate(study_year = as.factor(recruit_year),
         macroalgae_sqrt = sqrt(macroalgae)) %>%
  select(-habitat)

recruit <- read.csv(here::here("clean_data", "recruits_clean.csv"))
juv <- read.csv(here::here("clean_data", "juv_clean.csv"))
live_recruit <- read.csv(here::here("clean_data", "live_recruits_clean.csv"))

##### MERGE DATA #####
# aggregate LTA to total coral LTA
adult_agg <- adult %>%
  rename(adult_density = density) %>%
  mutate(study_year = as.factor(recruit_year)) %>%
  select(-c(cover_year, recruit_year)) %>%
  group_by(study_year, region, habitat, site_name) %>%
  summarise(adult_density_sum = sum(adult_density, na.rm = T)) %>%
  ungroup() %>%
  mutate(adult_density_sum = adult_density_sum*10E-04) %>% # convert LTA to m2
  mutate(adult_density_log = log(adult_density_sum + 1))

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

# merge recruit predictor datasets
recruit_pred <- adult_agg %>%
  left_join(env_rec_sub, by = c("study_year", "site_name", "region")) %>%
  left_join(macro, by = c("study_year", "site_name", "region")) 

# aggregate recruit data - TotalSto already serves as the sum of all stony coral recruits
recruit_agg <- recruit %>%
  mutate(study_year = as.factor(recruit_year)) %>%
  filter(taxa %in% c("TotalSto")) %>%
  select(-c(abundance, total_tiles, tile_area)) %>%
  select(-taxa) %>% # only have one taxa - recruit density
  rename(recruit_density = density) %>%
  mutate(recruit_density_log = log(recruit_density + 1))

# aggregate live recruit data - Stony already serves as the sum of all stony coral recruits
live_recruit_agg <- live_recruit %>%
  mutate(study_year = as.factor(recruit_year)) %>%
  filter(taxa %in% c("Stony")) %>%
  select(-c(abundance, total_tiles, tile_area)) %>%
  select(-taxa) %>% # only have one taxa - recruit density
  mutate(live_recruit_density_log = log(live_recruit_density + 1))

# merge recruits with predictor data
recruit_full <- recruit_pred %>%
  left_join(recruit_agg, c("study_year", "site_name",
                           "recruit_year", "cover_year"))  %>%
  left_join(live_recruit_agg, c("study_year", "site_name",
                           "recruit_year", "cover_year"))  %>%
  select(-c(cover_year, recruit_year, habitat)) 

# aggregate juvenile data
juv_agg <- juv %>%
  mutate(study_year = as.factor(recruit_year)) %>%
  select(-c(cover_year, recruit_year, habitat)) %>%
  group_by(study_year, site_name) %>%
  summarise(juv_density_sum = sum(density, na.rm = T)) %>%
  ungroup() %>%
  mutate(juv_density_log = log(juv_density_sum + 1))

# merge with the rest of the data
coral_full <- recruit_full %>%
  left_join(juv_agg, c("study_year", "site_name"))

#### CREATE GROUP MEAN CENTERED DF ####
# define predictors
explanatory_vars <- c(
  "sst_mean", "dhw_log", "kd_log", "depth", "adult_density_log", "macroalgae_sqrt", "recruit_density_log",
  "live_recruit_density_log"
)

# center data
coral_center <- gmc_data(coral_full)

#### COMPONENT MODELS ####
coral_rec_mod <- nlme::lme(recruit_density_log_dev ~ 
                            recruit_density_log_region +
                            sst_mean_region + sst_mean_dev +
                            dhw_log_region + dhw_log_dev +
                            kd_log_region + kd_log_dev +
                            depth_region + depth_dev +
                            adult_density_log_region + adult_density_log_dev +
                            macroalgae_sqrt_region + macroalgae_sqrt_dev,
                          random = ~ 1 | region/site_name,
                          na.action = na.omit,
                          data = coral_center)
summary(coral_rec_mod)
performance::check_model(coral_rec_mod) # looks fine - high VIF from some mundlak devices

coral_juv_mod <- nlme::lme(juv_density_log ~ 
                            sst_mean_region + sst_mean_dev +
                            dhw_log_region + dhw_log_dev +
                            kd_log_region + kd_log_dev +
                            depth_region + depth_dev +
                            recruit_density_log_region + recruit_density_log_dev +
                            adult_density_log_region + adult_density_log_dev +
                            macroalgae_sqrt_region + macroalgae_sqrt_dev,
                          random = ~ 1 | region/site_name,
                          na.action = na.omit,
                          data = coral_center)
summary(coral_juv_mod)
performance::check_model(coral_juv_mod)  # looks fine

#### SEM ####
sem <- psem(coral_rec_mod, coral_juv_mod)
(sem_summary <- summary(sem) )

#### WRITE SUMMARY TABLES ####
# R2 values
sem_r2 <- sem_summary$R2
# write.csv(sem_r2, "summary_tables/scor_sem_r2.csv", row.names = F)

# Standardized estimates (all)
sem_coefs_full <- sem_summary$coefficients %>%
  select(Response, Predictor, Std.Estimate, DF, Crit.Value, P.Value)
# write.csv(sem_coefs_full, "summary_tables/scor_sem_coefs_full.csv", row.names = F)

# Standardized estimates (just anomalies)
sem_anomaly <- sem_coefs_full %>%
  filter(str_detect(Predictor, "_dev"))
# write.csv(sem_anomaly, "summary_tables/scor_sem_coefs_anomaly.csv", row.names = F)

# tests of directed separation
sem_dsep <- sem_summary$dTable
# write.csv(sem_dsep, "summary_tables/scor_sem_dsep.csv", row.names = F)



#### COMPONENT MODELS ####
coral_liverec_mod <- nlme::lme(live_recruit_density_log_dev ~ 
                             live_recruit_density_log_region +
                             sst_mean_region + sst_mean_dev +
                             dhw_log_region + dhw_log_dev +
                             kd_log_region + kd_log_dev +
                             depth_region + depth_dev +
                             adult_density_log_region + adult_density_log_dev +
                             macroalgae_sqrt_region + macroalgae_sqrt_dev,
                           random = ~ 1 | region/site_name,
                           na.action = na.omit,
                           data = coral_center)
summary(coral_liverec_mod)
performance::check_model(coral_liverec_mod) # looks fine - high VIF from some mundlak devices

coral_livejuv_mod <- nlme::lme(juv_density_log ~ 
                             sst_mean_region + sst_mean_dev +
                             dhw_log_region + dhw_log_dev +
                             kd_log_region + kd_log_dev +
                             depth_region + depth_dev +
                             live_recruit_density_log_region + live_recruit_density_log_dev +
                             adult_density_log_region + adult_density_log_dev +
                             macroalgae_sqrt_region + macroalgae_sqrt_dev,
                           random = ~ 1 | region/site_name,
                           na.action = na.omit,
                           data = coral_center)
summary(coral_livejuv_mod)
performance::check_model(coral_livejuv_mod)  # looks fine

#### SEM ####
live_sem <- psem(coral_liverec_mod, coral_livejuv_mod)
(live_sem_summary <- summary(live_sem) )

#### WRITE SUMMARY TABLES ####
# R2 values
live_sem_r2 <- live_sem_summary$R2
# write.csv(live_sem_r2, "summary_tables/scor_live_sem_r2.csv", row.names = F)

# Standardized estimates (all)
live_sem_coefs_full <- live_sem_summary$coefficients %>%
  select(Response, Predictor, Std.Estimate, DF, Crit.Value, P.Value)
# write.csv(live_sem_coefs_full, "summary_tables/scor_live_sem_coefs_full.csv", row.names = F)

# Standardized estimates (just anomalies)
live_sem_anomaly <- live_sem_coefs_full %>%
  filter(str_detect(Predictor, "_dev"))
# write.csv(live_sem_anomaly, "summary_tables/scor_live_sem_coefs_anomaly.csv", row.names = F)

# tests of directed separation
live_sem_dsep <- live_sem_summary$dTable
# write.csv(live_sem_dsep, "summary_tables/scor_live_sem_dsep.csv", row.names = F)




