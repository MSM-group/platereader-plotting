# Install packages
pacman::p_load("tidyverse", "readxl", "randomcoloR", "lubridate", "ggrepel",
               "RColorBrewer", "ggplot2", "broom", "maditr", "ggpubr")

# Read in the slopes
tmafils <- list.files(paste0("output/"), pattern = "64DHH", full.names = T)
tmafil_slopes <- tmafils[grep("_calculated_slopes.csv", tmafils)]
tmafil_slopes

### First analyze butyrate
tma_but <- tmafil_slopes[grep("butyrate", tmafil_slopes)]
tma_but_raw <- rawdat <- tibble(filename = tma_but) %>%
  mutate(file_contents = purrr::map(filename,          # read files into
                                    ~ read_csv(file.path(.))) # a new data column
  ) %>%
  unnest(.) 
summary(tma_but_raw$max_slope)
tma_but


### Then analyze trimethylacetate
tma_tri <- tmafil_slopes[grep("trimethylacetate", tmafil_slopes)]
tma_tri_raw <- rawdat <- tibble(filename = tma_tri) %>%
  mutate(file_contents = purrr::map(filename,          # read files into
                                    ~ read_csv(file.path(.))) # a new data column
  ) %>%
  unnest(.) 
tma_tri_raw
summary(tma_tri_raw$max_slope)

# Make a boxplot of the two
tma_comb <- bind_rows(tma_but_raw, tma_tri_raw) %>%
  dplyr::filter(max_slope < 0.75) %>%
  dplyr::mutate(substrate = ifelse(grepl("butyrate", filename), "butyrate", "trimethylacetate"))
table(tma_comb$substrate)

ggplot(tma_comb, aes(x = substrate, y = max_slope)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16,
             outlier.size=2, notch=FALSE) +
  theme_pubr() +
  geom_text_repel(data = subset(tma_comb, max_slope > 0.2), aes(label = org))
write_csv(tma_comb, "output/all_slopes_combined.csv")
