# Script for plotting platereader output for 4-nitrophenyl* substrates

# Install packages
library("tidyverse")
library("readxl")
library("lubridate")
library("reshape2")
library("RColorBrewer")
library("randomcoloR")
library("hms")

# Set the substrate (compound) and the date for enzyme activity screening
enzym <- "4NP_butyrate_p55_pH" 
yyyymmdd <- "20230503"

# Read in the plate template
folder_path <- file.path(paste0("data/", yyyymmdd, "/"))
file_path <- list.files(folder_path, pattern = "4NP_butyrate_Platelayout", full.names = T)
file_path <- file_path[!grepl(pattern="~",file_path)]
file_path
temp <- read_excel(file_path) %>%
  janitor::clean_names() %>%
  dplyr::select(-x1) %>% 
  as.matrix(bycol = T) %>%
  t %>%
  as.vector()
temp # lots of NAs are ok
temp[1:12] <- paste0("pH6_", temp[1:12])
temp[13:24] <- paste0("pH7_", temp[13:24])
temp[25:36] <- paste0("pH8_", temp[25:36])
temp[37:48] <- paste0("pH9_", temp[37:48])

# Read in the raw platereader data
tmafils <- list.files(paste0(folder_path, "/"), pattern = enzym, full.names = T) 
tmafils <- tmafils[!grepl("~|setup|Bradford|screenshot|split|template|Tris", tmafils)] # remove any temporary files
tmafils # check the file name is right# 

tma <- read_excel(tmafils, range = "B40:CU71", col_types = c("date", rep("numeric", 97))) # NOTE Changed CU60 -> CU90
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
  dplyr::mutate(variable = gsub("\\.[[:digit:]]", "", variable))

dat2 <- resbind %>%
  dplyr::filter(!grepl("stdcurve|Tris|^0$|emptvec", variable)) # Filter out variables you don't want

dat3 <- dat2 %>%
  group_by(variable, time) %>%
  summarise_each(funs(mean, sd), value)  # calculate mean activity

# Custom color palette
pal <- colorRampPalette(brewer.pal(8,"Set1"))(8) 
pal2 <- c("gray80", "black", "dodgerblue", "goldenrod",  pal[c(1, 3:5, 8)], "blue", "gold1", distinctColorPalette(60))
pal2

# pH6
dat4 <- tibble(dat3) %>%
  dplyr::filter(grepl("pH6", variable))
pdf(paste0("output/", yyyymmdd, "_", enzym, "_pH6.pdf"),width = 9, height = 5)
pl <- ggplot(dat4, aes(x=time, y=mean, color=variable)) +
  geom_point() +
  labs(y = "Absorbance (410 nm)", x = "Time (minutes)") +
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

# pH7
pdf(paste0("output/", yyyymmdd, "_", enzym, "_pH7.pdf"), width = 9, height = 5)
dat4 <- tibble(dat3) %>%
  dplyr::filter(grepl("pH7", variable))
pl <- ggplot(dat4, aes(x=time, y=mean, color=variable)) +
  geom_point() +
  labs(y = "Absorbance (410 nm)", x = "Time (minutes)") +
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

# pH8
dat4 <- tibble(dat3) %>%
  dplyr::filter(grepl("pH8", variable))
pdf(paste0("output/", yyyymmdd, "_", enzym, "_pH8.pdf"), width = 9, height = 5)
pl <- ggplot(dat4, aes(x=time, y=mean, color=variable)) +
  geom_point() +
  labs(y = "Absorbance (410 nm)", x = "Time (minutes)") +
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
# pH 9
dat4 <- tibble(dat3) %>%
  dplyr::filter(grepl("pH9", variable))
pdf(paste0("output/", yyyymmdd, "_", enzym, "_pH9.pdf"), width = 9, height = 5)
pl <- ggplot(dat4, aes(x=time, y=mean, color=variable)) +
  geom_point() +
  labs(y = "Absorbance (410 nm)", x = "Time (minutes)") +
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



pdf(paste0("output/", yyyymmdd, "_", enzym, "_without_errorbars.pdf"), width = 13, height = 8)
pl <- ggplot(dat3, aes(x=time, y=mean, color=variable)) +
  geom_point() +
  labs(y = "Absorbance (410 nm)", x = "Time (minutes)") +
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
  scale_color_manual(values=pal2) +
  ylim(-5, 5)
pl
dev.off()
pl


# Calculate slopes
a <- dat3 %>%
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

# Find max slope for each organism
resmax <- resl %>%
  dplyr::filter(r2 >= 0.5) %>% # make sure R^2 is above or equal to 0.7
  group_by(org) %>%
  summarise_each(funs(max_slope = max), slope) #%>%
  #dplyr::filter(max_slope > 0.0 ) # SET activity threshold of 0.1
resmax

write_csv(resmax, paste0("output/", yyyymmdd, "_", enzym, "_slope_differences.csv"))

# Merge with the original dataset
slope_merg <- resmax %>%
  inner_join(., resl, by = "org") %>%
  dplyr::filter(r2 >= 0.9) %>%
  group_by(org) %>%
  dplyr::filter(slope == max(slope)) %>%
  dplyr::select(org, max_slope, r2, intercept)
slope_merg
