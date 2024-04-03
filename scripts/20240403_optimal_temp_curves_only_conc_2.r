library("tidyverse")
library("readxl")
library("lubridate")
library("reshape2")
library("RColorBrewer")
library("randomcoloR")
library("hms")
library("ggpubr")

# Set the substrate (compound) and the date for enzyme activity screening
enzym <- "lacticaseibacillus" 
yyyymmdd <- "20230623"

# Read in the results
files <- list.files("output/", pattern = paste0("_slopes_normalized"), full.names = T)
files <- files[!grepl("_5_|_1_", files)] # remove old files

readin <- tibble(filename = files) %>%
  mutate(file_contents = purrr::map(filename,          
                                    ~ read_csv(file.path(.), col_names = T)) # a new data column
  ) %>%
  unnest(.,cols = c(file_contents))

merg_summ <- readin %>%
  dplyr::mutate(filesplit = gsub(yyyymmdd, "", filename)) %>%
  dplyr::mutate(concentration_factor = parse_number(filesplit)) %>%
  dplyr::mutate(temp_num = parse_number(pH)) %>%
  group_by(concentration_factor, temp_num) %>%
  summarise_each(funs(mean, sd), nm_slope)

merg_summ$concentration_factor <- as.factor(merg_summ$concentration_factor)


pdf(paste0("output/", yyyymmdd, "_", enzym, ID, "_combined_final_graph.pdf"))
pl <- ggplot(merg_summ,  aes(x = temp_num, y = mean, color= concentration_factor)) + 
  geom_point() +
  geom_line(aes(y = mean, group = concentration_factor), size = 1) +
  geom_errorbar(aes(ymax=mean + sd, ymin = mean - sd), width=0.3, size=0.6) +
  theme_pubr() +
  xlab("Temperature (°C)") +
  ylab("Activity \n nmol pNP/nmol enzyme/minute") +
  scale_color_manual(values = c("#bfbfff", "cornflowerblue", "navy"))
pl
dev.off()

pdf(paste0("output/", yyyymmdd, "_", enzym, ID, "_combined_final_graph.pdf"))
library(ggplot2)
library(dplyr)
library(ggpubr)  # Assuming theme_pubr() comes from the ggpubr package

# Filter the dataframe to only include rows where concentration_factor equals 2
filtered_data <- merg_summ %>% filter(concentration_factor == 2)

# Now plot using the filtered_data
pl <- ggplot(filtered_data,  aes(x = temp_num, y = mean, color = concentration_factor)) + 
  geom_point() +
  geom_line(aes(y = mean, group = concentration_factor), size = 1) +
  #geom_errorbar(aes(ymax = mean + sd, ymin = mean - sd), width = 0.3, size = 0.6) +
  theme_pubr() +  # Apply publication ready theme
  xlab("Temperature (°C)") +
  ylab("Activity \n nmol pNP/nmol enzyme/minute") +
  scale_color_manual(values = "black")+
  theme(legend.position = "none")

ggsave("20240304_temp_optimum_only_conc_2.png", pl, width = 5, height = 3, dpi=900)


ggsave(paste0("output/", yyyymmdd, "_", enzym, "_temperature_optimum_graph.png"), pl, width = 5, height = 5)

