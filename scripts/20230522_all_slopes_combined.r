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
  unnest(cols = c(file_contents)) 
summary(tma_but_raw$max_slope)
tma_but # now includes 20230421

### Then analyze trimethylacetate
tma_tri <- tmafil_slopes[grep("trimethylacetate", tmafil_slopes)]
tma_tri_raw <- rawdat <- tibble(filename = tma_tri) %>%
  mutate(file_contents = purrr::map(filename,          # read files into
                                    ~ read_csv(file.path(.))) # a new data column
  ) %>%
  unnest(cols = c(file_contents)) 
tma_tri # now includes 20230421

# Make a boxplot of the two
tma_comb <- bind_rows(tma_but_raw, tma_tri_raw) %>%
  dplyr::mutate(substrate = ifelse(grepl("butyrate", filename), "butyrate", "trimethylacetate"))
table(tma_comb$substrate)

ggplot(tma_comb, aes(x = substrate, y = max_slope)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16,
             outlier.size=2, notch=FALSE) +
  theme_pubr() +
  geom_text_repel(data = subset(tma_comb, max_slope > 0.2), aes(label = org))

# Now prepare data for writing
colnames(tma_comb)
tma_comb <- tma_comb %>%
  dplyr::mutate(substrate_parse = stringr::word(string = filename, sep = "4NP_", -1)) %>%
  dplyr::mutate(substrate = gsub("_[[:digit:]]_all_data_calculated_slopes.csv|_all_data_calculated_slopes.csv", "", substrate_parse)) %>%
  dplyr::mutate(date_parse = stringr::word(string = filename, sep = "_", 1)) %>%
  dplyr::mutate(date = gsub("output//", "", date_parse)) %>%
  dplyr::mutate(max_slope = round(max_slope, 4)) %>%
  dplyr::mutate(r2 = round(r2, 4)) %>%
  dplyr::select(date, filename, org, max_slope, r2, substrate)

colnames(tma_comb)

write_csv(tma_comb, "output/20230522_all_slopes_combined.csv")
