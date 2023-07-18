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
volume <- 1
ID <- "mean_2"
# Read in the results
files <- list.files("output/", pattern = paste0( "__all_data_calculated_slopes.csv"), full.names = T)
files <- files[grepl(yyyymmdd, files)]
files
readin <- tibble(filename = files) %>%
  mutate(file_contents = purrr::map(filename,          
                                    ~ read_csv(file.path(.), col_names = T)) # a new data column
  ) %>%
  unnest(.,cols = c(file_contents))

# we add 20 µL of enzyme in 200 µL total volume
# The enzyme started out 0.22mg/mL
# That is 52652.53 g/mol
#(20 µL) * (1 mL/1000 µL) * (0.22 mg/1 mL) * (1 g/1000 mg) * (1 mol/52652.53 g) * (1e9 nmoles / 1 mol)
nmolesenzyme <- volume * (1/1000) * (0.22) * (1/1000) * (1/52652.53) * (1e9) 
nmolesenzyme
#0.08356673 nmoles of enzyme per 20uL
# Conversion per Liter
# 0.08356673 * (1e6/200) = nm OR  µM

readclean <- readin %>%
  dplyr::select(-1) %>%
  filter(!grepl("max_slope", max_slope)) %>% # remove the column names
  dplyr::mutate(max_slope = as.numeric(max_slope))%>%
  dplyr::mutate(sec_slope = max_slope /120) %>% # 
  dplyr::mutate(nm_slope = sec_slope /nmolesenzyme) %>% # to make one nanomole of enzyme
  distinct() %>%
  dplyr::mutate(nm_slope = as.numeric(nm_slope)) %>%
  filter(org==ID)
readclean
write_csv(readclean, paste0("output/", yyyymmdd, "_", enzym, "_", volume, "_slopes_normalized.csv"))

merg_summ <- readclean %>%
  dplyr::mutate(temp_num = parse_number(pH)) %>%
  group_by(temp_num) %>%
  summarise_each(funs(mean, sd), nm_slope)

pdf(paste0("output/", yyyymmdd, "_", enzym, ID, "_combined_final_graph.pdf"))
pl <- ggplot(merg_summ,  aes(x = temp_num, y = mean)) + 
  geom_point() +
  geom_errorbar(aes(ymax=mean + sd, ymin = mean - sd), width=0.3,size=0.6) +
  theme_pubr() +
  xlab("Temperature (°C)") +
  ylab("Activity \n log(nmol pNP/nmol enzyme/ hr)")
pl
dev.off()

pl
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
pl
       