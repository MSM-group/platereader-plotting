# Script for plotting pH curves
# Install packages
library("tidyverse")
library("readxl")
library("lubridate")
library("reshape2")
library("RColorBrewer")
library("randomcoloR")
library("hms")

# Set the substrate (compound) and the date for enzyme activity screening
rm(list = ls())
enzym <- "paracetamol_amidase" 
yyyymmdd <- "20230626"
pH_row <- "pH9" # new entry, change for each pH

# Read in the plate template
folder_path <- file.path(paste0("data/", yyyymmdd, "/"))
file_path <- list.files(folder_path, pattern = "20230626_lacticaseibacillus_paracetamol_amidase_pH_optimum_Template", full.names = T)
file_path <- file_path[!grepl(pattern="~",file_path)]

temp <- read_excel(file_path) %>%
  janitor::clean_names() %>%
  dplyr::select(-x1) %>% 
  as.matrix(bycol = T) %>%
  t %>%
  as.vector()
temp[grepl("mutant",temp)]<-paste0(temp[grepl("mutant", temp)], "_", 1:3)

# Read in the raw platereader data
tmafils <- list.files(paste0(folder_path, "/"), pattern = enzym, full.names = T) 
tmafils <- tmafils[!grepl("~|setup|Bradford|screenshot|split|Template|Tris|pH4", tmafils)] # remove any temporary files
tmafils # check the file name is right# 

tma <- read_excel(tmafils, range = "B63:CU94", col_types = c("date", rep("numeric", 97)))
oldnam <- as.vector(t(colnames(tma))) 
oldnam
newnam <- c("time", "temperature_c", paste0(temp))
newnam
template_check <- bind_cols(oldnam, newnam)
template_check # compares template to platereader column names

# If the names match
colnames(tma) <- make.unique(newnam) # set the column names
colnames(tma)

# Clean up the data
tma1 <- tma %>%
  dplyr::select(-temperature_c) %>%   # remove the column for temperature (constant at 37C)
  dplyr::filter(!grepl("Time", time)) %>% # removes row if column names are duplicated
  dplyr::filter(complete.cases(.)) %>% # removes rows containing NAs
  dplyr::mutate(time = round(as.numeric(lubridate::hms(stringr::word(time, sep = " ", start = 2)))/60), 0)  %>%
  dplyr::select(-contains("NA"))# minutes 

# Convert from wide to long format
resbind <- tma1 %>%
  reshape2::melt(., id = 'time') %>%
  dplyr::mutate(value = as.numeric(value)) %>%
  dplyr::mutate(variable = gsub("\\.[[:digit:]]", "", variable)) %>%
  dplyr::filter(grepl(pH_row, variable)) # double check only one pH is selected

# Plot the standard curve
pNPs <- resbind %>% 
  dplyr::filter(grepl("stdcurve", variable)) %>% 
  dplyr::filter(!grepl("x", variable)) %>% # stdcurve_1 
  dplyr::mutate(µL = as.numeric(word(variable, sep = "_", -1))) %>%
  dplyr::mutate(mM = µL * (8/200)) %>% # 8 mM stock solution, 200 µL final well volume
  dplyr::mutate(nM = mM * 1e6) %>% 
  dplyr::mutate(nmol = nM * (1/1000) * (1/1000) * 200) %>% # to convert from nM to nanomoles -> first volume to µL and then * 200 µL
  dplyr::mutate(pH = word(variable, sep = "_", 1)) # nmoles = nmoles/L * 1L/1000 mL * 1ml/1000µL * 200 µL (final volume)
pNPs

# Linear regression
pNP_fit <- lm(value ~ nmol, data = pNPs)
summary(pNP_fit) 
b <- pNP_fit$coefficients[1]
m <- pNP_fit$coefficients[2]

pdf(paste0("output/", yyyymmdd, "_", enzym, "_", pH_row, "_standard_curve.pdf"))
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

# Now look at data
dat2 <- resbind %>%
  dplyr::filter(!grepl("stdcurve|Tris|^0$|emptvec", variable)) # Filter out variables you don't want

dat3 <- dat2 %>%
  dplyr::mutate(nmols_pNP = (value - b)/m) %>%
  group_by(variable, time) %>%
  summarise_each(funs(mean, sd), nmols_pNP) # calculate mean activity

# Custom color palette
pal <- colorRampPalette(brewer.pal(8,"Set1"))(8) 
pal2 <- c("gray80", "black", "dodgerblue", "goldenrod",  pal[c(1, 3:5, 8)], "blue", "gold1", distinctColorPalette(60))
pal2

pdf(paste0("output/", yyyymmdd, "_", enzym, "_", pH_row, "_without_errorbars.pdf"), width = 13, height = 8)
pl <- ggplot(dat3, aes(x=time, y=mean, color=variable)) +
  geom_point() +
  labs(y = "nmol 4NP produced", x = "Time (minutes)") +
  # geom_errorbar(aes(ymax=mean + sd, ymin = mean - sd), width=0.3,size=0.6)+
  theme(legend.title=element_blank(), axis.line=element_line(color="black"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        panel.background=element_blank(),
        text = element_text(size = 20),
        legend.key= element_rect(fill=NA, color=NA),
        legend.position="right") + 
  guides(shape = guide_legend(override.aes = list(size = 10))) +
  scale_color_manual(values=pal2) 
  
pl
dev.off()
pl

dat4 <- dat2 %>%
  dplyr::mutate(nmols_pNP = (value - b)/m) %>%
  dplyr::mutate(reps_combined = gsub("_1|_2|_3", "", variable)) %>%
  group_by(variable, time) %>%
  summarise_each(funs(mean), nmols_pNP) %>%
  tidyr::spread(., variable, nmols_pNP) # calculate mean activity

dat_normalized <- dat4 %>%
  dplyr::mutate(mean_1 = cur_data()[[2]]-cur_data()[[5]]) %>%
  dplyr::mutate(mean_2 = cur_data()[[3]]-cur_data()[[6]]) %>%
  dplyr::mutate(mean_3 = cur_data()[[4]]-cur_data()[[7]]) %>%
  dplyr::mutate(id = pH_row) 

dat_normalized_long <- dat_normalized %>%
  select(contains(c("time", "mean_", "id"))) %>%
  reshape2::melt(., id.vars = c("time", "id")) %>%
  dplyr::rename(.,  mean = value)
#colnames(dat4)

pdf(paste0("output/", yyyymmdd, "_", enzym,"_", pH_row, "_with_errorbars.pdf"), width = 14, height = 8)
pl <- ggplot(dat_normalized_long, aes(x=time, y=mean, color = variable)) +
  geom_point() +
  labs(y = "nmol 4NP produced", x = "Time (minutes)") +
  #  geom_errorbar(aes(ymax=mean + sd, ymin = mean - sd), width=0.3,size=0.6)+
  theme(legend.title=element_blank(), axis.line=element_line(color="black"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        panel.background=element_blank(),
        text = element_text(size = 20),
        legend.key= element_rect(fill=NA, color=NA),
        legend.position="right") + 
  guides(shape = guide_legend(override.aes = list(size = 10))) +
  scale_color_manual(values=pal2) +
  ylim(-10, 70)
pl
dev.off()
pl

# Calculate slopes
dat_normalized
a <- dat_normalized_long %>%
  ungroup() %>%
  dplyr::mutate(variable = as.character(variable)) %>%
  dplyr::filter(!is.na(variable)) %>% 
  dplyr::filter(!grepl("stdcurve", variable)) %>% # stdcurve|emptyvec
  dplyr::mutate(minutes = as.numeric(time)) %>%
  dplyr::mutate(minutes = case_when(grepl(max(time), time) ~ 60,
                                    TRUE ~ minutes)) %>%
  dplyr::group_by(variable) %>%
  dplyr::slice(0:30) %>%
  dplyr::ungroup() %>%
  dplyr::select(variable, minutes, mean)
a

## Function that calculates slopes
slopes <- function(d) { 
  m <- lm(mean ~ minutes, as.data.frame(d, stringsAsFactors = F))
  summ <- summary(m)
  r2 <- summ$r.squared
  intercept <- coef(m)[1]
  slope <- coef(m)[2]
  return(list(org = d$variable[1], r2 = r2, slope = slope, intercept = intercept))
}

## Take the slope from the first 4 minutes
windowsize <- 13 #windowsize 20 for pH4, pH5, pH6; windowsize 10 for pH7,pH8,pH85

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

# Find max slope for each organism
resmax <- resl %>%
  #dplyr::filter(r2 >= 0.8) %>% # make sure R^2 is above or equal to 0.8
  group_by(org) %>%
  summarise_each(funs(max_slope = max), slope) #%>%
  #dplyr::filter(max_slope > 0.25 ) # SET activity threshold of 0.1

# Merge with the original dataset
slope_merg <- resmax %>%
  inner_join(., resl, by = "org", multiple = "all") %>%
  #dplyr::filter(r2 >= 0.8) %>%
  group_by(org) %>%
  dplyr::filter(slope == max(slope)) %>%
  dplyr::mutate(pH = pH_row) %>%
  dplyr::select(org, max_slope, r2, intercept, pH)

merg_all <- slope_merg %>% 
  left_join(., a, by = c("org" = "variable"), multiple = "all") 

pdf(paste0("output/", yyyymmdd, "_", enzym, "_", pH_row, "_slopes_plotted.pdf"), width = 13, height = 8)
pl <- ggplot(merg_all,  aes(x = minutes, y = mean, color = org)) + 
  geom_point() +
  geom_abline(slope = merg_all$max_slope, intercept = merg_all$intercept) +
  labs(y = "nmol 4NP produced", x = "Time (minutes)") +
  theme(legend.title=element_blank(), axis.line=element_line(color="black"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        panel.background=element_blank(),
        text = element_text(size = 16),
        legend.position = "right",
        legend.key= element_rect(fill=NA, color=NA)) +
  # ylim(-10, 70) +
  guides(shape = guide_legend(override.aes = list(size = 10))) +
  scale_color_manual(values=pal2) 
# legend.position = "none")
pl
dev.off() 
pl

# Visually assess results and remove any that look strange
slope_final <- slope_merg[order(slope_merg$max_slope, decreasing = T),]
slope_final
write_csv(slope_final, paste0("output/", pH_row, yyyymmdd, "_", enzym, "_all_data_calculated_slopes.csv"))



