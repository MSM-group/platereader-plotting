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
enzym <- "degrees_all_data" 
yyyymmdd <- "20230605"
volume <- 10

# Read in the results
files <- list.files("output/", pattern = enzym, full.names = T)
files <- files[!grepl("pdf|png|4degrees", files)]
files
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


# we add 10 µL of enzyme in 200 µL total volume
# The enzyme started out 0.48mg/mL
# That is 52652.53 g/mol
#(20 µL) * (1 mL/1000 µL) * (0.22 mg/1 mL) * (1 g/1000 mg) * (1 mol/52652.53 g) * (1e9 nmoles / 1 mol)
nmolesenzyme <- volume * (1/1000) * (0.48) * (1/1000) * (1/52652.53) * (1e9) 
nmolesenzyme
#0.08356673 nmoles of enzyme per 20uL

readclean <- readin %>%
  dplyr::select(-1) %>%
  dplyr::filter(grepl("Lacticaseib", org)) %>%
  filter(!grepl("max_slope", max_slope)) %>% # remove the column names
  dplyr::mutate(sec_slope = as.numeric(max_slope)) %>%
  #dplyr::mutate(sec_slope = max_slope /120) %>% # REMOVE THIS CONVERSION
  dplyr::mutate(nm_slope = sec_slope /nmolesenzyme) %>%
  distinct() %>%
readclean

merg_summ <- readclean %>%
  dplyr::mutate(temperature = parse_number(org)) %>%
  group_by(temperature) %>%
  summarise_each(funs(mean, sd), nm_slope) 

pdf(paste0("output/", yyyymmdd, "_", enzym, "temp_stability_combined_final_graph.pdf"))
pl <- ggplot(merg_summ,  aes(x = temperature, y = mean)) + 
  geom_point() +
  geom_errorbar(aes(ymax=mean + sd, ymin = mean - sd), width=0.3, size=0.6) +
  theme_pubr() +
  xlab("Temperature (°C)") +
  ylab("Activity \n nmol pNP/nmol enzyme/minute") 
pl
dev.off()

ggsave(plot = pl, filename = paste0("output/", yyyymmdd, "_", enzym, "temp_stability_combined_final_graph.png"), width = 8, height = 3)
