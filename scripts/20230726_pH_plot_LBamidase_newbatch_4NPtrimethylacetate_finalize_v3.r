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
volume <- 20 # enzyme volume 

# Read in the results
files <- list.files("output/", pattern = paste0(yyyymmdd, "_", enzym, "_all_data_calculated_slopes.csv"), full.names = T)
files
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

# we add 20 µL of enzyme in 200 µL total volume
# The enzyme started out 0.22mg/mL
# That is 52652.53 g/mol
#(20 µL) * (1 mL/1000 µL) * (0.22 mg/1 mL) * (1 g/1000 mg) * (1 mol/52652.53 g) * (1e9 nmoles / 1 mol)
nmolesenzyme <- volume * (1/1000) * (0.22) * (1/1000) * (1/52652.53) * (1e9) 
#0.08356673 nmoles of enzyme per 20uL
# Conversion per Liter
# 0.08356673 * (1e6/200) = nm OR  µM
nmolesenzyme

readclean <- readin %>%
  dplyr::select(-1) %>%
  filter(!grepl("max_slope", max_slope)) %>% # remove the column names
  dplyr::mutate(sec_slope = as.numeric(max_slope)) %>%
  #dplyr::mutate(sec_slope = max_slope /120) %>% # REMOVE THIS CONVERSION
  dplyr::mutate(nm_slope = sec_slope /nmolesenzyme) %>%
  distinct() %>%
  dplyr::mutate(nm_slope = as.numeric(nm_slope))
readclean

readclean$nm_slope

merg_summ <- readclean %>%
  dplyr::mutate(pH_num = parse_number(pH)) %>%
  group_by(pH_num) %>%
  summarise_each(funs(mean, sd), nm_slope)
merg_summ
# Fix the point5
grepl(".5", merg_summ$pH_num)

merg_summ$pH_num[grepl(".5", merg_summ$pH_num)] <- gsub("5", "\\.5", merg_summ$pH_num[grepl(".5", merg_summ$pH_num)])
merg_summ$pH_num <- as.numeric(merg_summ$pH_num)
merg_summ$pH_num
pdf(paste0("output/", yyyymmdd, "_", enzym, "_combined_final_graph.pdf"))
pl <- ggplot(merg_summ,  aes(x = pH_num, y = mean)) + 
  geom_point() +
  geom_errorbar(aes(ymax=mean + sd, ymin = mean - sd), width=0.3,size=0.6) +
  ylab("Activity \n nmol pNP/nmol enzyme/minute") +
  xlab("pH") +
  theme_pubr()
pl
dev.off()
ggsave(paste0("output/", yyyymmdd, "_", enzym, "_combined_final_graph.png"), pl, width = 4, height = 4)
pl
