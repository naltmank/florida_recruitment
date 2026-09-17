rm(list = ls())

#install.packages("librarian")
librarian::shelf(here, janitor, lubridate, tidyverse, ggplot2)

#### ENVIRONMENTAL DATA ####
# clean and combine erddap and depth data

#### RECRUIT DATA ####
# read data
recruit_raw <- read.csv(here::here("raw_data", "recruit_data.csv"))

#  "Red Dun" should be "Red Dun Reef" 
#  "Admiral Patch" should be "Admiral" 
# "West Turtle Shoals" should be "West Turtle Shoal"

recruit_corr <- recruit_raw %>%
  # rename year to recruit_year for future merges
  rename(recruit_year = Year) %>%
  mutate(Site = case_when(
    Site == "Red Dun" ~ "Red Dun Reef",
    Site == "Admiral Patch" ~ "Admiral",
    Site == "West Turtle Shoals" ~ "West Turtle Shoal",
    T ~ Site
  ),
  # add year that links to benthic cover data - previous year's cover influences recruits
  cover_year = recruit_year - 1, 
  # add the year that links to juvenile data - recruits influence next year's juvs
  juv_year = recruit_year + 1 
    ) %>%
  # some tiles were scored even though they were partially destroyed in a hurricane
  group_by(recruit_year, Site, Tile) %>%
  # if a site does not have all sides, mark it as complete ~ FALSE
  mutate(
    complete = all(c("Bottom", "Side", "Top") %in% Loc)
  ) %>%
  ungroup() %>%
  # create a new column indicating the number of complete tiles in a given site/year
  group_by(recruit_year, Site) %>%
  mutate(
    total_tiles = n_distinct(Tile[complete])
  ) %>%
  ungroup() %>%
  # filter out incomplete tiles
  filter(complete == TRUE) %>%
  # add a total area column
  mutate(tile_area = total_tiles * # for all tiles in a site
           (2*(0.15^2 - pi*(0.064/2)^2) + # 2 15x15 surfaces each with a 0.64 cm diameter screw in it
              0.02*0.15*4 )) # 4 sides, each 2*15 cm surfaces

# extract site data
sites <- unique(recruit_corr$Site)

# aggregate to site level
recruit_sum <- recruit_corr %>%
  group_by(recruit_year, cover_year, juv_year, Site, tile_area, total_tiles) %>%
  summarise(across(all_of(c("AGAR", "FAVI", "PORI", "SIDE", "Other", "UNKS",
                            "TotalSto", "TotalOct", "TotalCor")), sum)
  )

# pivot long
recruit_long <- recruit_sum %>%
  pivot_longer(cols = c("AGAR", "FAVI", "PORI", "SIDE", "Other", "UNKS",
                         "TotalSto", "TotalOct", "TotalCor"),
               names_to = "taxa",
               values_to = "abundance") %>%
  mutate(density = abundance/tile_area) %>%
  rename(site_name = Site)

#### JUVENILE DATA ####
juv_raw <- read.csv(here::here("raw_data", "juv_data.csv"))
juv_map <- read.csv(here::here("raw_data", "juv_id_map.csv")) %>%
  # rename juv_id to ID for merge with the raw data
  rename(ID = juv_ids)

# subset according to years and sites
juv_sub <- juv_raw %>%
  filter(Project_Year %in% c(2017, 2018, 2019), # juveniles linked to one year after recruit data
         Site_Name %in% sites) %>%
  mutate(recruit_year = Project_Year - 1, # add column that indicates what recruit year it links to - year prior
         cover_year = Project_Year - 1) %>% # for linking to cover data
  rename(juv_year = Project_Year) %>%
  group_by(juv_year, Site_Name) %>%
  mutate(
    total_quadrats = n_distinct(Quadrat),
    complete = total_quadrats == 32,
    # one quadrat missing from BC1 in 2017 - need to account for this in area (0.25 m2 quads)
    quadrat_area = total_quadrats * 0.25
  ) %>%
  ungroup() %>%
  # merge in the id map
  left_join(juv_map, by = "ID") 

# check all sites accounted for
setdiff(unique(juv_sub$Site_Name), sites)
setdiff(sites, unique(juv_sub$Site_Name)) 

# aggregate to site level
juv_long <- juv_sub %>%
  filter(Colony_Type %in% "Juv") %>%
  group_by(recruit_year, cover_year, juv_year, Site_Name, Region, Habitat,
           recruit_map, total_quadrats, quadrat_area) %>%
  summarise(abundance = n(),
            .groups = "drop") %>%
  mutate(density = abundance/quadrat_area) %>%
  rename(site_name = Site_Name,
         region = Region,
         habitat = Habitat,
         taxa = recruit_map)

#### LTA DATA ####
lta_map <- read.csv(here::here("raw_data", "lta_id_map.csv")) %>%
  rename(taxa = lta_ids)
cremp_lta_raw <- read.csv(here::here("raw_data", "cremp_coral_lta.csv"))
cremp_lta_sub <- cremp_lta_raw %>%
  mutate(Site_name = case_when(
    Site_name == "West Washer Women" ~ "West Washerwoman",
    T ~ Site_name
  )) %>%
  filter(Site_name %in% sites,
         Year %in% c(2015, 2016, 2017)) %>%
  mutate(recruit_year = Year + 1, # cover affects recruits the following year
         juv_year = Year + 2
  ) %>%
  rename(cover_year = Year)

secremp_lta_raw <- read.csv(here::here("raw_data", "secremp_coral_lta.csv"))
secremp_lta_sub <- secremp_lta_raw %>%
  # site names are full name instead of abbrev
  mutate(Site_name = case_when(
    Site_name == "Broward County 1" ~ "BC1",
    Site_name == "Broward County 2" ~ "BC2",
    Site_name == "Broward County 3" ~ "BC3",
    Site_name == "Broward County 4" ~ "BC4",
    Site_name == "Dade County 1" ~ "DC1",
    Site_name == "Dade County 2" ~ "DC2",
    Site_name == "Dade County 3" ~ "DC3",
    Site_name == "Dade County 4" ~ "DC4",
    Site_name == "Dade County 5" ~ "DC5",
    Site_name == "Dade County 6" ~ "DC6",
    Site_name == "Dade County 7" ~ "DC7",
    Site_name == "Dade County 8" ~ "DC8",
    T ~ Site_name
  ) ) %>%
  filter(Site_name %in% sites,
         Year %in% c(2015, 2016, 2017)) %>%
  mutate(recruit_year = Year + 1, # cover affects recruits the following year
         juv_year = Year + 2
  ) %>%
  rename(cover_year = Year) %>%
  select(-TransectArea)

# check all sites are accounted for 
setdiff(c(unique(cremp_lta_sub$Site_name), unique(secremp_lta_sub$Site_name)), sites)
setdiff(sites, c(unique(cremp_lta_sub$Site_name), unique(secremp_lta_sub$Site_name)))

# pivot longer and merge datasets
meta <- c("Year", "Date", "Subregion", "Habitat", "SiteID", "Site_name",                  
          "StationID", "TransectArea", "cover_year", "juv_year", "recruit_year")
cremp_lta_long <- cremp_lta_sub %>%
  pivot_longer(cols = colnames(cremp_lta_sub)[-which(colnames(cremp_lta_sub) %in% meta)],
               names_to = "taxa",
               values_to = "lta") %>%
  left_join(lta_map, by = "taxa") %>%
  group_by(cover_year, juv_year, recruit_year, Subregion, Habitat, Site_name, recruit_map) %>%
  summarise(lta = sum(lta, na.rm = T)) %>%
  mutate(density = lta/10) %>% # cremp lta data on 10m2 transects
  ungroup() %>%
  rename(site_name = Site_name,
         region = Subregion,
         habitat = Habitat,
         taxa = recruit_map) %>%
  filter(!is.na(taxa)) # millepora listed as NA - filter those out

secremp_lta_long <- secremp_lta_sub %>%
  pivot_longer(cols = colnames(secremp_lta_sub)[-which(colnames(secremp_lta_sub) %in% meta)],
               names_to = "taxa",
               values_to = "lta") %>%
  left_join(lta_map, by = "taxa") %>%
  group_by(cover_year, juv_year, recruit_year, Subregion, Habitat, Site_name, recruit_map) %>%
  summarise(lta = sum(lta, na.rm = T)) %>%
  mutate(density = lta/22) %>% # secremp lta data on 22m2 transects
  ungroup() %>%
  rename(site_name = Site_name,
         region = Subregion,
         habitat = Habitat,
         taxa = recruit_map) %>%
  filter(!is.na(taxa)) # millepora listed as NA - filter those out

comb_lta_long <- rbind(cremp_lta_long, secremp_lta_long)

#### COVER DATA ####
# really only care about macroalgae
cremp_cover <- read.csv(here::here("raw_data", "cremp_allbenthos_pcover.csv"))
cremp_macro <- cremp_cover %>%
  mutate(site_name = case_when(
    site_name == "West Washer Women" ~ "West Washerwoman",
    T ~ site_name
  )) %>%
  filter(group_fine == "Macroalgae") %>%
  filter(site_name %in% sites,
         year %in% c(2015, 2016, 2017)) %>%
  mutate(recruit_year = year + 1, # cover affects recruits the following year
         juv_year = year + 2
  ) %>%
  rename(cover_year = year)

secremp_cover <- read.csv(here::here("raw_data", "secremp_allbenthos_pcover.csv"))
secremp_macro <- secremp_cover %>%
  mutate(site_name = case_when(
    site_name == "Broward County 1" ~ "BC1",
    site_name == "Broward County 2" ~ "BC2",
    site_name == "Broward County 3" ~ "BC3",
    site_name == "Broward County 4" ~ "BC4",
    site_name == "Dade County 1" ~ "DC1",
    site_name == "Dade County 2" ~ "DC2",
    site_name == "Dade County 3" ~ "DC3",
    site_name == "Dade County 4" ~ "DC4",
    site_name == "Dade County 5" ~ "DC5",
    site_name == "Dade County 6" ~ "DC6",
    site_name == "Dade County 7" ~ "DC7",
    site_name == "Dade County 8" ~ "DC8",
    T ~ site_name
  ) )  %>%
  filter(group_fine == "Macroalgae") %>%
  filter(site_name %in% sites,
         year %in% c(2015, 2016, 2017)) %>%
  mutate(recruit_year = year + 1, # cover affects recruits the following year
         juv_year = year + 2
  ) %>%
  rename(cover_year = year) %>%
  select(-id)

comb_macro <- rbind(cremp_macro, secremp_macro)

#### CORAL DENSITY ####



#### OCTOCORAL DENSITY ####
# octo data are already in terms of density
cremp_octo <- read.csv(here::here("raw_data", "cremp_octo.csv"))
cremp_octo_corr <- cremp_octo %>%
  mutate(site_name = case_when(
    site_name == "West Washer Women" ~ "West Washerwoman",
    T ~ site_name
  ))

secremp_octo <- read.csv(here::here("raw_data", "secremp_octo.csv"))
secremp_octo_corr <- secremp_octo %>%
  mutate(site_name = case_when(
    site_name == "Broward County 1" ~ "BC1",
    site_name == "Broward County 2" ~ "BC2",
    site_name == "Broward County 3" ~ "BC3",
    site_name == "Broward County 4" ~ "BC4",
    site_name == "Dade County 1" ~ "DC1",
    site_name == "Dade County 2" ~ "DC2",
    site_name == "Dade County 3" ~ "DC3",
    site_name == "Dade County 4" ~ "DC4",
    site_name == "Dade County 5" ~ "DC5",
    site_name == "Dade County 6" ~ "DC6",
    site_name == "Dade County 7" ~ "DC7",
    site_name == "Dade County 8" ~ "DC8",
    T ~ site_name
  ) ) 

comb_oct <- rbind(cremp_octo_corr, secremp_octo_corr) %>%
  filter(site_name %in% sites,
         year %in% c(2016, 2017, 2018)) %>%
  mutate(recruit_year = year + 1, 
         juv_year = year + 2) %>%
  rename(cover_year = year)

# check all sites are present
setdiff( sites, unique(comb_oct$site_name))
