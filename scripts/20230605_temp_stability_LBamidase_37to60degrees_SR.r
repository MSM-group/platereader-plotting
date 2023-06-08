# Script for plotting platereader output for 4-nitrophenyl* substrates
# Comments use a hash symbol
# You can run the script line-by-line using Ctrl + Enter

# Check your working directory is correct e.g., should be github/platereader-plotting
getwd()
# If not correct then open the .RProj file File -> Open RProject

# Install packages
#install.packages("tidyverse")
#install.packages("readxl")
#install.packages("lubridate")
#install.packages("reshape2")
#install.packages("RColorBrewer")
#install.packages("randomcoloR")
#install.packages("hms")
library("tidyverse")
library("readxl")
library("lubridate")
library("reshape2")
library("RColorBrewer")
library("randomcoloR")
library("hms")

# Set the substrate (compound) and the date for enzyme activity screening
enzym <- "LBamidase" # change this for p055 etc.
ddmmyy <- "20230605"
temprange <-  "37to60degrees" #"70to90degrees"
  
# Read in the plate template
folder_path <- file.path(paste0("data/", ddmmyy, "/"))

file_path <- list.files(folder_path, pattern = paste0("Template_LBamidase_", temprange), full.names = T)

temp <- read_excel(file_path) %>%
  janitor::clean_names() %>%
  dplyr::select(-x1) %>% 
  as.matrix(bycol = T) %>%
  t %>%
  as.vector()
temp # lots of NAs are ok

# Read in the raw platereader data
tmafils <- list.files(paste0(folder_path, ""), pattern = temprange, full.names = T) 
tmafils <- tmafils[!grepl("~|setup|Bradford|screenshot|split|Template|Tris", tmafils)] # remove any temporary files
tmafils # check the file name is right

# Split out the files into three separate files for each triplicate
# NOTE: this is specific only for files where triplicates were done on the same plate
#check if lines correspond to the excel file
tma <- read_excel(tmafils, range = "B28:CU89", col_types = c("date", rep("numeric", 97)))
oldnam <- as.vector(t(colnames(tma))) 
oldnam
newnam <- c("time", "temperature_c", paste0(temp, "_rep", 1:3))
newnam

template_check <- bind_cols(oldnam, newnam)
template_check # compares template to platereader column names
# write_csv(template_check, "data/template_check.csv")
# Open the file and look at it to double check it matches
colnames(tma) <- make.unique(newnam) # set the column names
tma

#Clean up the data
tma1 <- tma %>%
  dplyr::select(-temperature_c) %>%   # remove the column for temperature (constant at 37C)
  dplyr::filter(!grepl("Time", time)) %>% # removes row if column names are duplicated
  dplyr::filter(complete.cases(.)) %>% # removes rows containing NAs
  dplyr::mutate(time = round(as.numeric(lubridate::hms(stringr::word(time, sep = " ", start = 2)))/60), 0) %>%
  dplyr::select(-contains("NA_"))# minutes 

# Convert from wide to long format
resbind <- tma1 %>%
  reshape2::melt(., id = 'time') %>%
  dplyr::mutate(value = as.numeric(value)) %>%
  dplyr::mutate(variable = gsub("\\.[[:digit:]]", "", variable))

# Plot the standard curve
pNPs <- resbind %>% 
  dplyr::filter(grepl("stdcurve", variable)) %>%
  dplyr::mutate(µL = as.numeric(word(variable, sep = "_", -2))) %>%
  dplyr::mutate(mM = µL * (8/200)) %>% # 8 mM stock solution, 200 µL final well volume
  dplyr::mutate(nM = mM * 1e6) %>%
  dplyr::mutate(nmol = nM * (1/1000) * (1/1000) * 200)  # nmoles = nmoles/L * 1L/1000 mL * 1ml/1000µL * 200 µL (final volume)

# Linear regression
pNP_fit <- lm(value ~ nmol, data = pNPs)
summary(pNP_fit)
pdf(paste0("output/", ddmmyy, "_", enzym, "_", temprange, "_standard_curve.pdf"))
pl <- ggplot(pNPs,  aes(x = nmol, y = value, color = time)) + 
  geom_point() +
  geom_smooth(method = "lm") +
  labs(y="Absorbance (410 nm)", x="nmol pNP") +
  theme(axis.line=element_line(color="black"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        panel.background=element_blank(),
        text = element_text(size = 8),
        legend.key= element_rect(fill=NA, color=NA),
        legend.position="right",
        legend.title = element_text(size = 8)) + 
  guides(shape = guide_legend(override.aes = list(size = 10)))
pl
dev.off()
pl

# Calculate slope and intercept of pNP standard curve
pNP_fit$coefficients
b <- pNP_fit$coefficients[1]
m <- pNP_fit$coefficients[2]
summary(pNP_fit)

# Now look at data
dat2 <- resbind %>%
  dplyr::filter(!grepl("stdcurve|Tris|^0$|emptvec|buffer", variable)) # Filter out variables you don't want

dat3 <- dat2 %>%
  dplyr::mutate(nmols_pNP = (value - b)/m) %>%
  group_by(variable, time)

# Calculate slopes
a <- dat3 %>%
  ungroup() %>%
  dplyr::mutate(variable = as.character(variable)) %>%
  dplyr::filter(!is.na(variable)) %>% 
  dplyr::filter(!grepl("stdcurve", variable)) %>% 
  dplyr::mutate(minutes = as.numeric(time)) %>%
  dplyr::mutate(minutes = case_when(grepl(max(time), time) ~ 60,
                                    TRUE ~ minutes)) %>%
  dplyr::group_by(variable) %>%
  dplyr::slice(0:30) %>%
  dplyr::ungroup() %>%
  dplyr::select(variable, minutes, nmols_pNP)
a

## Function that calculates slopes
slopes <- function(d) { 
  m <- lm(nmols_pNP ~ minutes, as.data.frame(d, stringsAsFactors = F))
  summ <- summary(m)
  r2 <- summ$r.squared
  intercept <- coef(m)[1]
  slope <- coef(m)[2]
  return(list(org = d$variable[1], r2 = r2, slope = slope, intercept = intercept))
}

## Take the slope from the first 5 minutes
windowsize <- 5 

# Calculate slopes for each enzyme
orgs <- unique(a$variable)

res <- list()
for(i in 1:length(orgs)) {
  tmp <- a[a$variable == orgs[i],]
  res[[i]] <- do.call(rbind.data.frame,lapply(seq(dim(tmp)[1] - windowsize),
                                              function(x) slopes(tmp[x:(x + windowsize),])))
}

names(res) <- orgs
resl <- plyr::ldply(res, data.frame)
resll <- do.call(rbind.data.frame, res)

# Find max slope for each enzyme and replicate
resmax <- resl %>%
  dplyr::filter(r2 >= 0.9) %>% # filters for R^2 is above or equal to 0.9
  group_by(org) %>%
  summarise_each(funs(max_slope = max), slope) %>%
  dplyr::filter(max_slope > 0.0 )
resmax

# Merge with the original dataset
slope_merg <- resmax %>%
  inner_join(., resl, by = "org") %>%
  dplyr::filter(r2 >= 0.9) %>%
  group_by(org) %>%
  dplyr::filter(slope == max(slope)) %>%
  dplyr::select(org, max_slope, r2, intercept)

# Plot slopes
merg_all <- slope_merg %>% 
  left_join(., a, by = c("org" = "variable")) %>% # to exclude inactive ones
  dplyr::mutate(winners = case_when(is.na(max_slope) ~ " inactive",
                                    TRUE ~ org))

# Visually assess results and remove any that look strange
slope_final <- slope_merg[order(slope_merg$max_slope, decreasing = T),]
slope_final
write_csv(slope_final, paste0("output/", ddmmyy, "_", enzym, "_", temprange, "_all_data_calculated_slopes.csv"))

