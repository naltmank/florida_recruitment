rm(list = ls())

#install.packages("librarian")
librarian::shelf(here, janitor, lubridate, tidyverse, ggplot2, vegan, ggpubr)

#### ADULT DENSITIES ####
adult <- read.csv(here::here("clean_data", "coral_density_clean_full.csv"))
octo <- read.csv(here::here("clean_data", "adult_octo_density_clean_full.csv")) 

# subset octocoral data to only relevant data for merge
octo_sub <- octo %>%
  mutate(year = as.factor(cover_year)) %>%
  select(year, site_name, OCTO) 

# filter out UNKS and Other
adult_sub <- adult %>%
  filter(!taxa %in% c("Other", "UNKS")) 

# pivot wide
adult_wide <- adult_sub %>%
  pivot_wider(names_from = taxa,
              values_from = density) %>%
  # set year to factors
  mutate(year = as.factor(cover_year), # setting recruit_year as the primary year of interest
         unique_id = paste0(site_name, "-", year),
         design = case_when(
           region == "BC" | region == "DC" ~ "SE FL",
           T ~ "FL Keys"
           )) %>%
  left_join(octo_sub, by = c("year", "site_name"))

# set metadata
adult_meta <- c("cover_year", "recruit_year", "region", "habitat", "site_name",
                "year", "unique_id", "design")

# subset species matrix and remove missing value from red dunn 2014
adult_matrix <- as.matrix(na.omit(adult_wide[,-which(names(adult_wide) %in% adult_meta)]))
rownames(adult_matrix) <- na.omit(adult_wide)$unique_id



# log transform for skewness
adult_matrix_log <- log(adult_matrix + 1)

# conduct PERMANOVA
adonis2(adult_matrix_log ~ year*design, method = "bray", by = "terms",
        data = na.omit(adult_wide))
# Df SumOfSqs      R2       F Pr(>F)    
# year          3   0.0276 0.00597  0.4072  0.882    
# design        1   2.0644 0.44647 91.3164  0.001 ***
#   year:design   3   0.0224 0.00484  0.3300  0.945    
# Residual    111   2.5094 0.54271                   
# Total       118   4.6238 1.00000       


# run ordination
set.seed(1032)
adult_ord2 <- metaMDS(adult_matrix_log, distance = "bray", k = 2, maxit = 999, trymax = 999) # 0.12

adult_site_scores <- as.data.frame(scores(adult_ord2, "sites"))
adult_site_scores$unique_id <- rownames(adult_site_scores)
adult_site_scores_full <- adult_site_scores %>%
  left_join(adult_wide, by = "unique_id")

adult_species_scores <- as.data.frame(scores(adult_ord2, "species"))
adult_species_scores$species <- rownames(adult_species_scores)
adult_species_scores <- adult_species_scores %>%
  mutate(species_fullname = case_when(species == "AGAR" ~ "Agaricidae",
                                      species == "FAVI" ~ "Faviidae",
                                      species == "PORI" ~ "Poritidae",
                                      species == "SIDE" ~ "Siderastridae",
                                      species == "OCTO" ~ "Octocorallia",
                                      T ~ species))

# plot
(adult_nmds_plot <- 
    ggplot() + 
    # add sites
    geom_point(data= adult_site_scores_full, aes(x=NMDS1,y=NMDS2, colour = year, shape = design)) + 
    # add species
    geom_text(data = adult_species_scores, aes(x = NMDS1, y = NMDS2, label=species_fullname)) +
    stat_ellipse(data= adult_site_scores_full, aes(x=NMDS1,y=NMDS2, colour = year, linetype = design),
                 level = 0.95) +
    labs(title = "a.", colour = "Year", shape = "Region", linetype = "Region") +
    theme_bw() + 
    scale_colour_manual(
      values = c(
        "2014" = "#E69F00",
        "2015" = "#56B4E9",
        "2016" = "#009E73",
        "2017" = "#9336fd",
        "2018" = "#CC79A7") ) +
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
recruit <- read.csv(here::here("clean_data", "recruits_clean_full.csv"))

# filter out UNKS, Other, TotalCor, and TotalSto, mutate TotalOct to OCTO
recruit_sub <- recruit %>%
  filter(!taxa %in% c("Other", "UNKS", "TotalCor", "TotalSto")) %>%
  mutate(taxa = case_when(taxa == "TotalOct" ~ "OCTO",
                          T ~ taxa)) %>%
  select(-c(abundance, total_tiles, tile_area))

# pull out region data from adults df
region_data <- adult_wide %>%
  distinct(region, site_name) 

# pivot wide
recruit_wide <- recruit_sub %>%
  
  pivot_wider(names_from = taxa,
              values_from = density) %>%
  left_join(region_data, by = c("site_name") ) %>%
  # set year to factors
  mutate(year = as.factor(recruit_year), # setting recruit_year as the primary year of interest
         unique_id = paste0(site_name, "-", year)
         ) %>%
  
  mutate(
         design = case_when(
           region == "BC" | region == "DC" ~ "SE FL",
           T ~ "FL Keys"
         )
  )

# set metadata
recruit_meta <- c("cover_year",  "recruit_year", "region", "habitat", "site_name",
                "year", "unique_id", "design")

# subset species matrix
recruit_matrix <- as.matrix(recruit_wide[,-which(names(recruit_wide) %in% recruit_meta)])
rownames(recruit_matrix) <- recruit_wide$unique_id

# log transform for skewness
recruit_matrix_log <- log(recruit_matrix + 1)

# conduct PERMANOVA
set.seed(1032)
adonis2(recruit_matrix_log ~ year*design, method = "bray", by = "terms",
        data = recruit_wide)
# Df SumOfSqs      R2      F Pr(>F)    
# year         2   1.6271 0.16527 9.2454  0.001 ***
#   design       1   0.6944 0.07053 7.8911  0.001 ***
#   year:design  2   0.1320 0.01341 0.7500  0.645    
# Residual    84   7.3916 0.75079                  
# Total       89   9.8451 1.00000     

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
recruit_species_scores <- recruit_species_scores %>%
  mutate(species_fullname = case_when(species == "AGAR" ~ "Agaricidae",
                                      species == "FAVI" ~ "Faviidae",
                                      species == "PORI" ~ "Poritidae",
                                      species == "SIDE" ~ "Siderastridae",
                                      species == "OCTO" ~ "Octocorallia",
                                      T ~ species))

# plot
(recruit_nmds_plot <- 
    ggplot() + 
    # add sites
    geom_point(data= recruit_site_scores_full, aes(x=NMDS1,y=NMDS2, colour = year, shape = design)) + 
    # add species
    geom_text(data = recruit_species_scores, aes(x = NMDS1, y = NMDS2, label=species_fullname)) +
    stat_ellipse(data= recruit_site_scores_full, aes(x=NMDS1,y=NMDS2, colour = year, linetype = design),
                 level = 0.95) +
    labs(title = "b.", colour = "Year", shape = "Region", linetype = "Region") +
    theme_bw() + 
    scale_colour_manual(
      values = c(
        "2014" = "#E69F00",
        "2015" = "#56B4E9",
        "2016" = "#009E73",
        "2017" = "#9336fd",
        "2018" = "#CC79A7") ) +
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
juv <- read.csv(here::here("clean_data", "juv_clean_full.csv"))

juv_sub <- juv %>%
  filter(!taxa %in% c("Other", "UNKS", NA)) %>%
  select(-c(abundance, total_quadrats, quadrat_area))

# add in octocorals
octo_juv <- read.csv(here::here("clean_data", "octo_juv_clean_full.csv"))
octo_juv_sub <- octo_juv %>%
  mutate(taxa = "OCTO") %>% # add octo taxa ID for merge
  select(-c(abundance, total_quadrats, quadrat_area))

# combine
juv_full <- rbind(juv_sub, octo_juv_sub)

# pivot wide
juv_wide <- juv_full %>%
  pivot_wider(names_from = taxa,
              values_from = density,
              values_fill = 0) %>%
  # set year to factors
  mutate(year = as.factor(recruit_year), # setting juv_year as the primary year of interest
         unique_id = paste0(site_name, "-", year)
  ) %>%
  mutate(
    design = case_when(
      region == "SEFL" ~ "SE FL",
      T ~ "FL Keys"
    )
  )

# set metadata
juv_meta <- c("cover_year",  "recruit_year", "region", "habitat", "site_name",
                  "year", "unique_id", "design")

# subset species matrix
juv_matrix <- as.matrix(juv_wide[,-which(names(juv_wide) %in% juv_meta)])
rownames(juv_matrix) <- juv_wide$unique_id

# log transform for skewness
juv_matrix_log <- log(juv_matrix + 1)

# conduct PERMANOVA
adonis2(juv_matrix_log ~ year*design, method = "bray", by = "terms",
        data = juv_wide)
# Df SumOfSqs      R2       F Pr(>F)    
# year         2   0.0492 0.01183  0.7378  0.589    
# design       1   1.2731 0.30626 38.1845  0.001 ***
#   year:design  2   0.0341 0.00819  0.5108  0.812    
# Residual    84   2.8005 0.67372                   
# Total       89   4.1568 1.00000                        


# run ordination
set.seed(1032)
juv_ord2 <- metaMDS(juv_matrix_log, distance = "bray", k = 2, maxit = 999, trymax = 999) # 0.17
juv_ord3 <- metaMDS(juv_matrix_log, distance = "bray", k = 3, maxit = 999, trymax = 999) # 0.10


juv_site_scores <- as.data.frame(scores(juv_ord3, "sites"))
juv_site_scores$unique_id <- rownames(juv_site_scores)
juv_site_scores_full <- juv_site_scores %>%
  left_join(juv_wide, by = "unique_id")

juv_species_scores <- as.data.frame(scores(juv_ord3, "species"))
juv_species_scores$species <- rownames(juv_species_scores)
juv_species_scores <- juv_species_scores %>%
  mutate(species_fullname = case_when(species == "AGAR" ~ "Agaricidae",
                                      species == "FAVI" ~ "Faviidae",
                                      species == "PORI" ~ "Poritidae",
                                      species == "SIDE" ~ "Siderastridae",
                                      species == "OCTO" ~ "Octocorallia",
                                      T ~ species))

# plot
(juv_nmds_plot <- 
    ggplot() + 
    # add sites
    geom_point(data= juv_site_scores_full, aes(x=NMDS1,y=NMDS2, colour = year, shape = design)) + 
    # add species
    geom_text(data = juv_species_scores, aes(x = NMDS1, y = NMDS2, label=species_fullname)) +
    stat_ellipse(data= juv_site_scores_full, aes(x=NMDS1,y=NMDS2, colour = year, linetype = design),
                 level = 0.95) +
    labs(title = "c.", colour = "Year", shape = "Region", linetype = "Region") +
    theme_bw() + 
   scale_colour_manual(
     values = c(
       "2014" = "#E69F00",
       "2015" = "#56B4E9",
       "2016" = "#009E73",
       "2017" = "#9336fd",
       "2018" = "#CC79A7") ) +
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

# ggsave(filename = "output/nmds_panel_v2.png", nmds_panel,
#        width = 16, height = 6, dpi = "retina")


