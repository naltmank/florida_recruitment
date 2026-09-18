rm(list = ls())

#install.packages("librarian")
librarian::shelf(here, janitor, lubridate, tidyverse, ggplot2, vegan)

#### ADULT DENSITIES ####
adult <- read.csv(here::here("clean_data", "coral_density_clean.csv"))
octo <- read.csv(here::here("clean_data", "adult_octo_density_clean.csv")) 

# subset octocoral data to only relevant data for merge
octo_sub <- octo %>%
  mutate(study_year = as.factor(recruit_year)) %>%
  select(study_year, site_name, OCTO) 

# filter out UNKS and Other
adult_sub <- adult %>%
  filter(!taxa %in% c("Other", "UNKS")) 

# pivot wide
adult_wide <- adult_sub %>%
  pivot_wider(names_from = taxa,
              values_from = density) %>%
  # set year to factors
  mutate(study_year = as.factor(recruit_year), # setting recruit_year as the primary year of interest
         unique_id = paste0(site_name, "-", study_year),
         design = case_when(
           region == "BC" | region == "DC" ~ "SE FL",
           T ~ "FL Keys"
           )) %>%
  left_join(octo_sub, by = c("study_year", "site_name"))

# set metadata
adult_meta <- c("cover_year", "juv_year", "recruit_year", "region", "habitat", "site_name",
                "study_year", "unique_id", "design")

# subset species matrix
adult_matrix <- as.matrix(adult_wide[,-which(names(adult_wide) %in% adult_meta)])
rownames(adult_matrix) <- adult_wide$unique_id

# log transform for skewness
adult_matrix_log <- log(adult_matrix + 1)

# conduct PERMANOVA
adonis2(adult_matrix_log ~ study_year, method = "bray", by = "terms",
        data = adult_wide)
# R2 = 0.004, F = 0.17, P = 0.97


# run ordination
set.seed(1032)
adult_ord2 <- metaMDS(adult_matrix_log, distance = "bray", k = 2, maxit = 999, trymax = 999) # 0.11

adult_site_scores <- as.data.frame(scores(adult_ord2, "sites"))
adult_site_scores$unique_id <- rownames(adult_site_scores)
adult_site_scores_full <- adult_site_scores %>%
  left_join(adult_wide, by = "unique_id")

adult_species_scores <- as.data.frame(scores(adult_ord2, "species"))
adult_species_scores$species <- rownames(adult_species_scores)

# plot
(adult_nmds_plot <- 
    ggplot() + 
    # add sites
    geom_point(data= adult_site_scores_full, aes(x=NMDS1,y=NMDS2, colour = study_year, shape = design)) + 
    # add species
    geom_text(data = adult_species_scores, aes(x = NMDS1, y = NMDS2, label=species)) +
    stat_ellipse(data= adult_site_scores_full, aes(x=NMDS1,y=NMDS2, colour = study_year, linetype = design),
                 level = 0.95) +
    labs(title = "a.", colour = "Year", shape = "Region") +
    theme_bw() + 
    theme(axis.title.x = element_text(size=18), 
          axis.title.y = element_text(size=18), 
          strip.text.x = element_text(size = 25),
          panel.grid.major = element_blank(),  #remove major-grid labels
          panel.grid.minor = element_blank(),  #remove minor-grid labels
          legend.position = 'bottom', # move legend to bottom
          legend.box.background = element_rect(colour = "black")) 
)

#### RECRUITS ####
recruit <- read.csv(here::here("clean_data", "recruits_clean.csv"))




# filter out UNKS, Other, TotalCor, and TotalSto, mutate TotalOct to OCTO
recruit_sub <- recruit %>%
  filter(!taxa %in% c("Other", "UNKS", "TotalCor", "TotalSto")) %>%
  mutate(taxa = case_when(taxa == "TotalOct" ~ "OCTO",
                          T ~ taxa)) %>%
  select(-c(abundance, total_tiles, tile_area))

# pull out region data from adults df
region_data <- adult_wide %>%
  select(region, site_name, study_year)

# pivot wide
recruit_wide <- recruit_sub %>%
  pivot_wider(names_from = taxa,
              values_from = density) %>%
  # set year to factors
  mutate(study_year = as.factor(recruit_year), # setting recruit_year as the primary year of interest
         unique_id = paste0(site_name, "-", study_year)
         ) %>%
  left_join(region_data, by = c("site_name", "study_year") ) %>%
  mutate(
         design = case_when(
           region == "BC" | region == "DC" ~ "SE FL",
           T ~ "FL Keys"
         )
  )

# set metadata
recruit_meta <- c("cover_year", "juv_year", "recruit_year", "region", "habitat", "site_name",
                "study_year", "unique_id", "design")

# subset species matrix
recruit_matrix <- as.matrix(recruit_wide[,-which(names(recruit_wide) %in% recruit_meta)])
rownames(recruit_matrix) <- recruit_wide$unique_id

# log transform for skewness
recruit_matrix_log <- log(recruit_matrix + 1)

# conduct PERMANOVA
adonis2(recruit_matrix_log ~ study_year*design, method = "bray", by = "terms",
        data = recruit_wide)

adonis2(formula = recruit_matrix_log ~ study_year , data = recruit_wide, method = "bray", by = "terms")
# R2 = 0.17, F = 8.6, P = 0.001


# run ordination
set.seed(1032)
recruit_ord2 <- metaMDS(recruit_matrix_log, distance = "bray", k = 2, maxit = 999, trymax = 999) # 0.19
recruit_ord3 <- metaMDS(recruit_matrix_log, distance = "bray", k = 3, maxit = 999, trymax = 999) # 0.13


recruit_site_scores <- as.data.frame(scores(recruit_ord3, "sites"))
recruit_site_scores$unique_id <- rownames(recruit_site_scores)
recruit_site_scores_full <- recruit_site_scores %>%
  left_join(recruit_wide, by = "unique_id")

recruit_species_scores <- as.data.frame(scores(recruit_ord3, "species"))
recruit_species_scores$species <- rownames(recruit_species_scores)

# plot
(recruit_nmds_plot <- 
    ggplot() + 
    # add sites
    geom_point(data= recruit_site_scores_full, aes(x=NMDS1,y=NMDS2, colour = study_year, shape = design)) + 
    # add species
    geom_text(data = recruit_species_scores, aes(x = NMDS1, y = NMDS2, label=species)) +
    stat_ellipse(data= recruit_site_scores_full, aes(x=NMDS1,y=NMDS2, colour = study_year, linetype = design),
                 level = 0.95) +
    labs(title = "b.", colour = "Year", shape = "Region") +
    theme_bw() + 
    theme(axis.title.x = element_text(size=18), 
          axis.title.y = element_text(size=18), 
          strip.text.x = element_text(size = 25),
          panel.grid.major = element_blank(),  #remove major-grid labels
          panel.grid.minor = element_blank(),  #remove minor-grid labels
          legend.position = 'bottom', # move legend to bottom
          legend.box.background = element_rect(colour = "black")) 
)

#### JUVENILES ####

