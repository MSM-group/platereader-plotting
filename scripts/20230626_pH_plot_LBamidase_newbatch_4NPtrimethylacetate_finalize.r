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
enzym <- "paracetamol_amidase" 
yyyymmdd <- "20230626"

# Read in the results
files <- list.files("output/", pattern = paste0(yyyymmdd, "_", enzym, "_all_data_calculated_slopes.csv"), full.names = T)

readin <- tibble(filename = files) %>%
  mutate(file_contents = purrr::map(filename,          
                                    ~ read_csv(file.path(.), col_names = F)) # a new data column
  ) %>%
  unnest(.,cols = c(file_contents))

colnames(readin) <- c("filename",
"org",
"max_slope",
"r2",
"intercept",
"pH")


readclean <- readin %>%
  dplyr::select(-1) %>%
  filter(!grepl("max_slope", max_slope)) %>% # remove the column names
  distinct() %>%
  dplyr::mutate(max_slope = as.numeric(max_slope))
readclean

readclean$max_slope

merg_summ <- readclean %>%
  dplyr::mutate(pH_num = parse_number(pH)) %>%
  group_by(pH_num) %>%
  summarise_each(funs(mean, sd), max_slope)

# Fix the point5
grepl(".5", merg_summ$pH_num)

merg_summ$pH_num[grepl(".5", merg_summ$pH_num)] <- gsub("5", "\\.5", merg_summ$pH_num[grepl(".5", merg_summ$pH_num)])
merg_summ$pH_num <- as.numeric(merg_summ$pH_num)
merg_summ$pH_num
pdf(paste0("output/", yyyymmdd, "_", enzym, "_combined_final_graph.pdf"))
pl <- ggplot(merg_summ,  aes(x = pH_num, y = mean)) + 
  geom_point() +
  geom_errorbar(aes(ymax=mean + sd, ymin = mean - sd), width=0.3,size=0.6) +
  theme_pubr()
pl
dev.off()
ggsave(paste0("output/", yyyymmdd, "_", enzym, "_combined_final_graph.png"),pl)
pl
