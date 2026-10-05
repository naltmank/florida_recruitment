rm(list = ls())
#install.packages("librarian")
librarian::shelf(here, janitor, lubridate, tidyverse, ggplot2, nlme,  performance, DHARMa,
                 car, broom.mixed)

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
  left_join(juv_expand, by = c("study_year", "site_name", "taxa"))

# quick plot - looks like porites, agaricia, and octos should be correlated
ggplot() +
  geom_point(data = juv_complete, aes(x = recruit_density_log, y = density_log)) +
  facet_wrap(~taxa)

##### SUBSET AND CENTER #####
# define predictors
explanatory_vars <- c(
  "sst_mean", "dhw_log", "kd_log", "depth", "recruit_density_log", "macroalgae_sqrt", "adult_density_log"
)
# taxa are "AGAR" "FAVI" "PORI" "SIDE" "OCTO"

# AGAR
agar_juv <- subset(juv_complete, taxa == "AGAR")
agar_juv_center <- gmc_data(agar_juv,
                            timevar = study_year,
                            spatial_group = region,
                            explanatory_vars = explanatory_vars)
ggplot() +
  geom_point(data = agar_juv_center, aes(x = recruit_density_log_dev, y = density_log)) +
  geom_smooth(data = agar_juv_center, aes(x = recruit_density_log_dev, y = density_log), method = "lm")


# FAVI
favi_juv <- subset(juv_complete, taxa == "FAVI")
favi_juv_center <- gmc_data(favi_juv,
                            timevar = study_year,
                            spatial_group = region,
                            explanatory_vars = explanatory_vars)

# PORI
pori_juv <- subset(juv_complete, taxa == "PORI")
pori_juv_center <- gmc_data(pori_juv,
                            timevar = study_year,
                            spatial_group = region,
                            explanatory_vars = explanatory_vars)

ggplot() +
  geom_point(data = pori_juv_center, aes(x = recruit_density_log_dev, y = density_log)) +
  geom_smooth(data = pori_juv_center, aes(x = recruit_density_log_dev, y = density_log), method = "lm")


# SIDE
side_juv <- subset(juv_complete, taxa == "SIDE")
side_juv_center <- gmc_data(side_juv,
                            timevar = study_year,
                            spatial_group = region,
                            explanatory_vars = explanatory_vars)

# OCTO
octo_juv <- subset(juv_complete, taxa == "OCTO")
octo_juv_center <- gmc_data(octo_juv,
                            timevar = study_year,
                            spatial_group = region,
                            explanatory_vars = explanatory_vars)

#### JUV MODELS ####
###### AGAR ######
agar_juv_mod <- nlme::lme(density_log ~ 
                            sst_mean_region + sst_mean_dev +
                            dhw_log_region + dhw_log_dev +
                            kd_log_region + kd_log_dev +
                            depth_region + depth_dev +
                            recruit_density_log_region + recruit_density_log_dev +
                            adult_density_log_region + adult_density_log_dev +
                            macroalgae_sqrt_region + macroalgae_sqrt_dev,
                          random = ~ 1 | region/site_name,
                          na.action = na.omit,
                          data = agar_juv_center)
summary(agar_juv_mod) # adults sig
performance::check_model(agar_juv_mod) # all look fine
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
                            adult_density_log_region + adult_density_log_dev +
                            macroalgae_sqrt_region + macroalgae_sqrt_dev,
                          random = ~ 1 | region/site_name,
                          na.action = na.omit,
                          data = favi_juv_center)
summary(favi_juv_mod) # nothing sig
performance::check_model(favi_juv_mod) # all look fine
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
                            adult_density_log_region + adult_density_log_dev +
                            macroalgae_sqrt_region + macroalgae_sqrt_dev,
                          random = ~ 1 | region/site_name,
                          na.action = na.omit,
                          data = pori_juv_center)
summary(pori_juv_mod) # adults sig
performance::check_model(pori_juv_mod) # all look fine
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
                            adult_density_log_region + adult_density_log_dev +
                            macroalgae_sqrt_region + macroalgae_sqrt_dev,
                          random = ~ 1 | region/site_name,
                          na.action = na.omit,
                          data = side_juv_center)
summary(side_juv_mod) # no sig
performance::check_model(side_juv_mod) # all look fine
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
                            adult_density_log_region + adult_density_log_dev +
                            macroalgae_sqrt_region + macroalgae_sqrt_dev,
                          random = ~ 1 | region/site_name,
                          na.action = na.omit,
                          data = octo_juv_center)
summary(octo_juv_mod) # adluts, recruits sig positive, depth, kd, dhw sig negative 
performance::check_model(octo_juv_mod) # all look fine
# extract effects
octo_effects <- broom.mixed::tidy(octo_juv_mod, effects = "fixed", conf.int = TRUE) %>%
  dplyr::filter(term != "(Intercept)") %>%
  select(term, estimate, conf.low, conf.high, p.value) %>%
  mutate(taxa = "OCTO")


# write summary tables
juv_summary <- rbind(agar_effects, favi_effects, pori_effects, side_effects, octo_effects) %>%
  mutate(response = "Juveniles")

# write.csv(juv_summary, "summary_tables/juv_glmm_effects.csv", row.names = F)
