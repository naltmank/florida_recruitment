rm(list = ls())

#install.packages("librarian")
librarian::shelf(here, janitor, lubridate, tidyverse, ggplot2, vegan, ggpubr)

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
adult_meta <- c("cover_year", "recruit_year", "region", "habitat", "site_name",
                "study_year", "unique_id", "design")

# subset species matrix
adult_matrix <- as.matrix(adult_wide[,-which(names(adult_wide) %in% adult_meta)])
rownames(adult_matrix) <- adult_wide$unique_id

# log transform for skewness
adult_matrix_log <- log(adult_matrix + 1)

# conduct PERMANOVA
adonis2(adult_matrix_log ~ study_year*design, method = "bray", by = "terms",
        data = adult_wide)
# Df SumOfSqs      R2       F Pr(>F)    
# study_year         1  0.00468 0.00203  0.2088  0.873    
# design             1  1.04278 0.45221 46.5094  0.001 ***
#   study_year:design  1  0.00293 0.00127  0.1305  0.911    
# Residual          56  1.25557 0.54449                   
# Total             59  2.30595 1.00000  


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
    labs(title = "a.", colour = "Year", shape = "Region", linetype = "Region") +
    theme_bw() + 
    scale_colour_manual(values = c("2016" = "#f77976",
                                   "2017" = "#9336fd") ) +
    guides(colour = guide_legend(title.position = "top"),
           shape = guide_legend(title.position = "top")) +
    theme(axis.title.x = element_text(size=18), 
          axis.title.y = element_text(size=18), 
          plot.title = element_text(size=20), 
          strip.text.x = element_text(size = 25),
          panel.background = element_blank(), 
          panel.grid.major = element_blank(),  #remove major-grid labels
          panel.grid.minor = element_blank(),  #remove minor-grid labels
          legend.position = 'bottom', # move legend to bottom
          legend.title = element_text(size = 25, hjust = 0.5),
          legend.text = element_text(size = 24)
    ) 
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
recruit_meta <- c("cover_year",  "recruit_year", "region", "habitat", "site_name",
                "study_year", "unique_id", "design")

# subset species matrix
recruit_matrix <- as.matrix(recruit_wide[,-which(names(recruit_wide) %in% recruit_meta)])
rownames(recruit_matrix) <- recruit_wide$unique_id

# log transform for skewness
recruit_matrix_log <- log(recruit_matrix + 1)

# conduct PERMANOVA
set.seed(1032)
adonis2(recruit_matrix_log ~ study_year*design, method = "bray", by = "terms",
        data = recruit_wide)
# Df SumOfSqs      R2      F Pr(>F)   
# study_year         1   0.4295 0.06496 4.2900  0.008 **
#   design             1   0.4200 0.06353 4.1954  0.008 **
#   study_year:design  1   0.1556 0.02353 1.5540  0.178   
# Residual          56   5.6067 0.84798                 
# Total             59   6.6118 1.00000 


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
    labs(title = "b.", colour = "Year", shape = "Region", linetype = "Region") +
    theme_bw() + 
    scale_colour_manual(values = c("2016" = "#f77976",
                                   "2017" = "#9336fd") ) +
    guides(colour = guide_legend(title.position = "top"),
           shape = guide_legend(title.position = "top")) +
    theme(axis.title.x = element_text(size=18), 
          axis.title.y = element_text(size=18), 
          plot.title = element_text(size=20), 
          strip.text.x = element_text(size = 25),
          panel.background = element_blank(), 
          panel.grid.major = element_blank(),  #remove major-grid labels
          panel.grid.minor = element_blank(),  #remove minor-grid labels
          legend.position = 'bottom', # move legend to bottom
          legend.title = element_text(size = 25, hjust = 0.5),
          legend.text = element_text(size = 24)
    ) 
)

#### JUVENILES ####
juv <- read.csv(here::here("clean_data", "juv_clean.csv"))

juv_sub <- juv %>%
  filter(!taxa %in% c("Other", "UNKS", NA)) %>%
  select(-c(abundance, total_quadrats, quadrat_area))


# pivot wide
juv_wide <- juv_sub %>%
  pivot_wider(names_from = taxa,
              values_from = density,
              values_fill = 0) %>%
  # set year to factors
  mutate(study_year = as.factor(recruit_year), # setting juv_year as the primary year of interest
         unique_id = paste0(site_name, "-", study_year)
  ) %>%
  mutate(
    design = case_when(
      region == "SEFL" ~ "SE FL",
      T ~ "FL Keys"
    )
  )

# set metadata
juv_meta <- c("cover_year",  "recruit_year", "region", "habitat", "site_name",
                  "study_year", "unique_id", "design")

# subset species matrix
juv_matrix <- as.matrix(juv_wide[,-which(names(juv_wide) %in% juv_meta)])
rownames(juv_matrix) <- juv_wide$unique_id

# log transform for skewness
juv_matrix_log <- log(juv_matrix + 1)

# conduct PERMANOVA
adonis2(juv_matrix_log ~ study_year*design, method = "bray", by = "terms",
        data = juv_wide)
# Df SumOfSqs      R2       F Pr(>F)    
# study_year         1   0.0265 0.00547  0.4338  0.743    
# design             1   1.3725 0.28315 22.4621  0.001 ***
#   study_year:design  1   0.0264 0.00545  0.4323  0.748    
# Residual          56   3.4218 0.70593                   
# Total             59   4.8472 1.00000 


# run ordination
set.seed(1032)
juv_ord2 <- metaMDS(juv_matrix_log, distance = "bray", k = 2, maxit = 999, trymax = 999) # 0.13
juv_ord3 <- metaMDS(juv_matrix_log, distance = "bray", k = 3, maxit = 999, trymax = 999) # 0.08


juv_site_scores <- as.data.frame(scores(juv_ord3, "sites"))
juv_site_scores$unique_id <- rownames(juv_site_scores)
juv_site_scores_full <- juv_site_scores %>%
  left_join(juv_wide, by = "unique_id")

juv_species_scores <- as.data.frame(scores(juv_ord3, "species"))
juv_species_scores$species <- rownames(juv_species_scores)

# plot
(juv_nmds_plot <- 
    ggplot() + 
    # add sites
    geom_point(data= juv_site_scores_full, aes(x=NMDS1,y=NMDS2, colour = study_year, shape = design)) + 
    # add species
    geom_text(data = juv_species_scores, aes(x = NMDS1, y = NMDS2, label=species)) +
    stat_ellipse(data= juv_site_scores_full, aes(x=NMDS1,y=NMDS2, colour = study_year, linetype = design),
                 level = 0.95) +
    labs(title = "c.", colour = "Year", shape = "Region", linetype = "Region") +
    theme_bw() + 
    scale_colour_manual(values = c("2016" = "#f77976",
                                   "2017" = "#9336fd") ) +
    guides(colour = guide_legend(title.position = "top"),
           shape = guide_legend(title.position = "top")) +
    theme(axis.title.x = element_text(size=18), 
          axis.title.y = element_text(size=18), 
          plot.title = element_text(size=20), 
          strip.text.x = element_text(size = 25),
          panel.background = element_blank(), 
          panel.grid.major = element_blank(),  #remove major-grid labels
          panel.grid.minor = element_blank(),  #remove minor-grid labels
          legend.position = 'bottom', # move legend to bottom
          legend.title = element_text(size = 25, hjust = 0.5),
          legend.text = element_text(size = 24)
    ) 
)

#### PANEL GRAPH ####
(nmds_panel <- ggarrange(adult_nmds_plot, recruit_nmds_plot, juv_nmds_plot,
                         nrow = 1, ncol = 3, widths = c(1,1,1),
                         common.legend = T, legend = "bottom")
)

# ggsave(filename = "output/nmds_panel.png", nmds_panel,
#        width = 16, height = 6, dpi = "retina")


