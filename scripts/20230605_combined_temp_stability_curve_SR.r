# Install packages
library("tidyverse")
library("readxl")
library("lubridate")
library("reshape2")
library("RColorBrewer")
library("randomcoloR")
library("hms")
library("ggpubr")

# Set the substrate (compound) and the date for enzyme activity screening
rm(list = ls())
enzym <- "degrees" 
yyyymmdd <- "20230605"

# Read in the results
files <- list.files("output/", pattern = enzym, full.names = T)

readin <- tibble(filename = files) %>%
  mutate(file_contents = purrr::map(filename,          
                                    ~ read_csv(file.path(.), col_names = T)) # a new data column
  ) %>%
  unnest(.,cols = c(file_contents))

colnames(readin) <- c("filename",
                      "org",
                      "max_slope",
                      "r2",
                      "intercept")


readclean <- readin %>%
  dplyr::select(-1) %>%
  dplyr::filter(grepl("Lacticaseib", org)) %>%
  filter(!grepl("max_slope", max_slope)) %>% # remove the column names
  distinct() %>%
  dplyr::mutate(max_slope = as.numeric(max_slope))
readclean

merg_summ <- readclean %>%
  dplyr::mutate(temperature = parse_number(org)) %>%
  group_by(temperature) %>%
  summarise_each(funs(mean, sd), max_slope) 


pdf(paste0("output/", yyyymmdd, "_", enzym, "temp_stability_combined_final_graph.pdf"))
pl <- ggplot(merg_summ,  aes(x = temperature, y = mean)) + 
  geom_point() +
  geom_errorbar(aes(ymax=mean + sd, ymin = mean - sd), width=0.3,size=0.6) +
  theme_pubr()
pl
dev.off()
