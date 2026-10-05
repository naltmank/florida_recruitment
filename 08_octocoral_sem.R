rm(list = ls())
#install.packages("librarian")
librarian::shelf(here, janitor, lubridate, tidyverse, ggplot2, nlme, performance, piecewiseSEM)


#### FUNCTIONS ####
gmc_data <- function(df, # dataframe to pass the function to
                     timevar = study_year, # for if sites were surveyed more than one year
                     spatial_group = region, # the spatial grouping for setting regional means
                     explanatory_vars# a vector of your explanatory variables
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
env <- read.csv(here::here("clean_data", "env_data_clean.csv"))
macro <- read.csv(here::here("clean_data", "macro_cover_clean.csv")) %>%
  mutate(study_year = as.factor(recruit_year),
         macroalgae_sqrt = sqrt(macroalgae)) %>%
  select(-habitat)

adult_octo <- read.csv(here::here("clean_data", "adult_octo_density_clean.csv")) 
# reformat octo to be in "long" format for merge
adult_octo_long <- adult_octo %>%
  pivot_longer(cols = OCTO,
               names_to = "taxa",
               values_to = "density") %>%
  rename(adult_density = density) %>%
  mutate(adult_density_log = log(adult_density + 1),
         study_year = as.factor(recruit_year))

recruit <- read.csv(here::here("clean_data", "recruits_clean.csv"))
recruit_sub <- recruit %>%
  mutate(study_year = as.factor(recruit_year)) %>%
  filter(taxa %in% c("TotalOct")) %>%
  mutate(taxa = case_when(taxa == "TotalOct" ~ "OCTO",
                          T ~ taxa)) %>% 
  select(-c(abundance, total_tiles, tile_area)) %>%
  rename(recruit_density = density) %>%
  mutate(recruit_density_log = log(recruit_density + 1) )

octo_juv <- read.csv(here::here("clean_data", "octo_juv_clean.csv")) %>%
  mutate(taxa = "OCTO")

octo_juv_sub <- octo_juv %>%
  mutate(study_year = as.factor(recruit_year)) %>%
  select(-c(cover_year, recruit_year, habitat, region, total_quadrats, quadrat_area)) %>%
  mutate(juv_density_log = log(density + 1))

##### MERGE #####
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
recruit_pred <- adult_octo_long %>%
  left_join(env_rec_sub, by = c("study_year", "site_name", "region")) %>%
  left_join(macro, by = c("study_year", "site_name", "region", "cover_year", "recruit_year"))

# add in recruit data
recruit_full <- recruit_pred %>%
  left_join(recruit_sub, c("study_year", "site_name", 
                           "recruit_year", "cover_year", "taxa")) %>%
  select(-c(cover_year, recruit_year))

# merge octocoral data
octocoral_full <- recruit_full %>%
  left_join(octo_juv_sub, c("study_year", "site_name", "taxa"))


#### CREATE GROUP MEAN CENTERED DF ####
# define predictors
explanatory_vars <- c(
  "sst_mean", "dhw_log", "kd_log", "depth", "adult_density_log", "macroalgae_sqrt", "recruit_density_log")

# center data
octocoral_center <- gmc_data(octocoral_full, explanatory_vars = explanatory_vars)

#### COMPONENT MODELS ####
octocoral_rec_mod <- nlme::lme(recruit_density_log_dev ~ 
                             recruit_density_log_region +
                             sst_mean_region + sst_mean_dev +
                             dhw_log_region + dhw_log_dev +
                             kd_log_region + kd_log_dev +
                             depth_region + depth_dev +
                             adult_density_log_region + adult_density_log_dev +
                             macroalgae_sqrt_region + macroalgae_sqrt_dev,
                           random = ~ 1 | region/site_name,
                           na.action = na.omit,
                           data = octocoral_center)
summary(octocoral_rec_mod)
performance::check_model(octocoral_rec_mod) # looks fine - high VIF from some mundlak devices

octocoral_juv_mod <- nlme::lme(juv_density_log ~ 
                             sst_mean_region + sst_mean_dev +
                             dhw_log_region + dhw_log_dev +
                             kd_log_region + kd_log_dev +
                             depth_region + depth_dev +
                             recruit_density_log_region + recruit_density_log_dev +
                             adult_density_log_region + adult_density_log_dev +
                             macroalgae_sqrt_region + macroalgae_sqrt_dev,
                           random = ~ 1 | region/site_name,
                           na.action = na.omit,
                           data = octocoral_center)
summary(octocoral_juv_mod)
performance::check_model(octocoral_juv_mod)  # looks fine

#### SEM ####
sem <- psem(octocoral_rec_mod, octocoral_juv_mod)
(sem_summary <- summary(sem) )

#### WRITE SUMMARY TABLES ####
# R2 values
sem_r2 <- sem_summary$R2
# write.csv(sem_r2, "summary_tables/octo_sem_r2.csv", row.names = F)

# Standardized estimates (all)
sem_coefs_full <- sem_summary$coefficients %>%
  select(Response, Predictor, Std.Estimate, DF, Crit.Value, P.Value)
# write.csv(sem_coefs_full, "summary_tables/octo_sem_coefs_full.csv", row.names = F)

# Standardized estimates (just anomalies)
sem_anomaly <- sem_coefs_full %>%
  filter(str_detect(Predictor, "_dev"))
# write.csv(sem_anomaly, "summary_tables/octo_sem_coefs_anomaly.csv", row.names = F)







