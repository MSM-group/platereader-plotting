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
files <- list.files("output/", pattern = paste0( "__all_data_calculated_slopes.csv"), full.names = T)
files <- files[grepl(yyyymmdd, files)]
files
readin <- tibble(filename = files) %>%
  mutate(file_contents = purrr::map(filename,          
                                    ~ read_csv(file.path(.), col_names = T)) # a new data column
  ) %>%
  unnest(.,cols = c(file_contents))

# we add 5 µL of enzyme in 200 µL total volume
# The enzyme started out 0.5 mg/ml
# That is 66004.46 g/mol
# (5 µL) * (1 mL/1000 µL) * (0.5 mg/1 mL) * (1 g/1000 mg) * (1 mol/66004.46 g) * (1e9 nmoles / 1 mol)
# 5 * (1/1000) * (0.5) * (1/1000) * (1/66004.46) * (1e9) = 0.03787623 nmoles in 200 µL
# Conversion per Liter
# 0.03787623 * (1e6/200) = 189 nm OR 0.189 µM

readclean <- readin %>%
  dplyr::select(-1) %>%
  filter(!grepl("max_slope", max_slope)) %>%
  dplyr::mutate(hr_slope = max_slope * 2 * 60) %>% # per hour 
  dplyr::mutate(nm_slope = max_slope * 27.02703) %>% # to make one nanomole of enzyme
  dplyr::mutate(log_slope = log10(nm_slope)) %>% # remove the column names
  distinct() %>%
  dplyr::mutate(log_slope = as.numeric(log_slope))

merg_summ <- readclean %>%
  dplyr::mutate(temp_num = parse_number(pH)) %>%
  group_by(temp_num) %>%
  summarise_each(funs(mean, sd), log_slope)

pdf(paste0("output/", yyyymmdd, "_", enzym, "_combined_final_graph.pdf"))
pl <- ggplot(merg_summ,  aes(x = temp_num, y = mean)) + 
  geom_point() +
  geom_errorbar(aes(ymax=mean + sd, ymin = mean - sd), width=0.3,size=0.6) +
  theme_pubr() +
  xlab("Temperature (°C)") +
  ylab("Activity \n log(nmol pNP/nmol enzyme/ hr)")
pl
dev.off()


pdf(paste0("output/", yyyymmdd, "_", enzym, "_combined_final_graph_std_dev.pdf"))
pl <- ggplot(merg_summ,  aes(x = temp_num, y = mean)) + 
  geom_point() +
  geom_line(aes(y = mean), size = 1) + 
  geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd), alpha = 0.5) +
  geom_errorbar(aes(ymax=mean + sd, ymin = mean - sd), width=0.3,size=0.6) +
  theme_pubr()+
  theme(text = element_text(size = 20)) +
  xlab("Temperature (°C)") +
  ylab("Activity \n log(nmol pNP/nmol enzyme/ hr)")   ### TO CALCULATE ###
pl
dev.off()
ggsave(paste0("output/", yyyymmdd, "_", enzym, "_temperature_optimum_graph.png"), pl, width = 8, height = 5)

       