# Script for plotting platereader output for 4-nitrophenyl* substrates
library("tidyverse")
library("readxl")
library("lubridate")
library("reshape2")
library("RColorBrewer")
library("randomcoloR")
library("hms")
library("ggpubr")

# Read in the data
tmafils <- "data/20251017/20251008_colorimetric_data_labeled.xlsx"
tma <- read_excel(tmafils, col_types = c("date", rep("numeric", 384))) %>%
  janitor::clean_names()

tma1 <- tma %>%  # remove the column for temperature (constant at 37C)
  dplyr::filter(complete.cases(.)) %>% # removes rows containing NAs
  dplyr::mutate(time = round(as.numeric(lubridate::hms(stringr::word(time, sep = " ", start = 2)))/60), 0) %>%
  dplyr::select(-contains("NA.")) %>%
  dplyr::select(-contains("trimethylacet")) %>% 
  dplyr::select(-contains("paracet"))

# Get all column names except the 'Time' column
all_cols <- names(tma1)[-1]
all_cols

# Extract the unique substrate names. The substrate name is everything after the first underscore.
# E.g., from "C8_paracetamol" or "RFP_paracetamol", we extract "paracetamol".
substrate_groups <- word(all_cols, sep = "_", 2)
unique_substrates <- tolower(unique(substrate_groups)[!is.na(unique(substrate_groups))])
# Initialize a new data frame with the 'Time' column
corrected_data <- data.frame(time = tma1$time)

# Loop through each unique substrate to perform background subtraction
for (substrate in unique_substrates) {
  
  # 1. Construct the expected RFP (blank) column name
  rfp_col_name <- paste0("rfp_", substrate)
  rfp_col_name
  
  # Get the RFP (blank) values
  rfp_values <- rowMeans(tma1[, grep(rfp_col_name, colnames(tma1))])
  print(rfp_values)
  colnames(rfp_values)
  
  # 2. Identify all concentration columns (C*) for this specific substrate
  # This regex finds columns starting with 'C' followed by digits, then '_', and ending with the substrate name.
  cx_col_names <- grep(substrate, names(tma1), value = TRUE)
  cx_col_names <- cx_col_names[!grepl("rfp", cx_col_names)]
  if (length(cx_col_names) == 0) {
    # If only the RFP column exists for a substrate, we skip it
    next 
  }
 
  # 3. Perform subtraction for each concentration column
  for (col_name in cx_col_names) {
    # Define the new corrected column name
    new_col_name <- paste0(col_name, "_corrected")
    
    # Calculate the corrected absorbance: Value - RFP Blank
    corrected_values <- tma1[[col_name]] - rfp_values
    
    # Add the new column to the output data frame
    corrected_data[[new_col_name]] <- corrected_values
  }
}

# --- Step 5: Save Results ---
write.csv(corrected_data, row.names = F,
          file = "data/20251017/datcorrect.csv")

resbind <- corrected_data %>%
  reshape2::melt(., id = 'time') %>%
  dplyr::mutate(value = as.numeric(value)) %>%
  dplyr::filter(variable != "0")
colnames(corrected_data)  

# Read in the Boccomero standard curve
stdcurve <- read_excel("data/20251017/pnitroaniline_stdcurve.xlsx")

# 200 µM of stock solution -> 70 µL added to a final volume 80µL
# (200)(70) = (x)(80)

pNP_fit <- lm(abs ~ conc, data = stdcurve)
summary(pNP_fit)
pdf(paste0("data/20251017/4nitroaniline_standard_curve.pdf"))
pl <- ggplot(stdcurve,  aes(x = conc, y = abs)) + 
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

# Calculate slope and intercept of pNP standard curve
pNP_fit$coefficients
b <- pNP_fit$coefficients[1] #y = mx +b
m <- pNP_fit$coefficients[2]

dat3 <- resbind %>%
  dplyr::mutate(nmols_pNP = (value - b)/m) 

# Calculate slopes
a <- dat3 %>%
  dplyr::mutate(variable = as.character(variable)) %>%
  dplyr::filter(!is.na(variable)) %>% 
  dplyr::filter(!grepl("stdcurve", variable)) %>% # stdcurve|emptyvec
  dplyr::mutate(minutes = as.numeric(time)) %>%
  dplyr::mutate(minutes = case_when(grepl(max(time), time) ~ 60,
                                    TRUE ~ minutes)) %>%
  dplyr::group_by(variable) %>%
  dplyr::slice(30:327) %>% # 1000 data points
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
  return(list(enz = d$variable[1], r2 = r2, slope = slope, intercept = intercept))
}

## Take the slope from rolling window of 30 minutes
windowsize <- 30 

# Calculate slopes for each enzyme
enz <- unique(a$variable)

res <- list()
for(i in 1:length(enz)) {
  tmp <- a[a$variable == enz[i],]
  res[[i]] <- do.call(rbind.data.frame,lapply(seq(dim(tmp)[1] - windowsize),
                                              function(x) slopes(tmp[x:(x + windowsize),])))
}

names(res) <- enz
resl <- plyr::ldply(res, data.frame)
resll <- do.call(rbind.data.frame, res)

# Find max slope for each organism
resmax <- resl %>%
  #dplyr::filter(r2 >= 0.2) %>% # make sure R^2 is above a threshold
  group_by(enz) %>%
  summarise_each(funs(max_slope = max), slope) %>%
  dplyr::filter(max_slope > 0) 
resmax

# Merge with the original dataset
slope_merg <- resmax %>%
  inner_join(., resl, by = "enz") %>%
  group_by(enz) %>%
  dplyr::filter(slope == max(slope)) %>%
  dplyr::select(enz, max_slope, r2, intercept)

# Plot winners on graph
merg_all <- slope_merg %>% 
  left_join(., a, by = c("enz" = "variable")) %>% # to exclude inactive ones
  dplyr::mutate(winners = case_when(is.na(max_slope) ~ " inactive",
                                    TRUE ~ enz)) %>%
  dplyr::mutate(substrate = stringr::word(enz, sep = "_", 2)) %>%
  dplyr::mutate(colony = stringr::word(enz, sep = "_", 1))

# Custom color palette
pal <- colorRampPalette(brewer.pal(8,"Set1"))(8) 
pal2 <- c("gray80", "gray80", "black", "dodgerblue", "goldenrod",  pal[c(1, 3:5, 8)], "blue", "gold1", distinctColorPalette(60))
pal2

pl <- ggplot(merg_all,  aes(x = minutes, y = nmols_pNP)) + 
  geom_point(alpha = 0.5, aes(color = colony), size = 1) +
  geom_abline(aes(slope = max_slope, intercept = intercept)) +
  labs(y = "nmol 4NP produced", x = "Time (minutes)") +
  theme(legend.title=element_blank(), axis.line=element_line(color="black"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        panel.background=element_blank(),
        text = element_text(size = 16),
        legend.position = "right",
        legend.key= element_rect(fill=NA, color=NA)) +
  facet_wrap(~substrate, scales = "free") +
  guides(shape = guide_legend(override.aes = list(size = 10))) +
  scale_color_manual(values=pal2) 
pl

pdf("output/prot_norm_split_by_substrate.pdf", width = 20)
pl
dev.off()

###### 

slope_final0 <- slope_merg[order(slope_merg$max_slope, decreasing = T),] %>%
  dplyr::mutate(substrate = word(enz, sep = "_", 2)) %>%
  dplyr::mutate(colony = word(enz, sep = "_", 1))

# Now read in the protein concentrations
prot2 <- read_excel("data/20251017/20251007_directed_evo_protein_concentration_expressers.xlsx")
prot2$Sample <- tolower(prot2$Sample)

slope_final <- slope_final0 %>%
  left_join(., prot2, by = c("colony" = "Sample"))

slope_final <- slope_final %>%
  mutate(max_slope = max_slope / Protein_Conc_nM_per_OD)
# write_csv(slope_final, paste0("data/20251017/calculated_slopes.csv"))

# Make the histogram
target_means <- slope_final %>%
  dplyr::filter(colony %in% c("wt", "rfp")) %>%
  group_by(substrate, colony) %>%
  summarise(
    mean_slope = mean(max_slope, na.rm = TRUE),
    .groups = 'drop'
  )
#data <- read_csv("data/20251017/prot_norm_calculated_slopes.csv")

# --- 2. Calculate WT Reference Means and Fold Change ---

# Calculate the mean 'wt' slope for each 'substrate'
# This serves as the baseline (Fold Change = 1.0) for each substrate.
wt_means <- slope_final %>%
  filter(colony == "wt") %>%
  group_by(substrate) %>%
  summarise(
    wt_mean_slope = mean(max_slope, na.rm = TRUE),
    .groups = 'drop'
  )

# Join the WT mean back to the original data and calculate the Fold Change
data_fc <- slope_final %>%
  left_join(wt_means, by = "substrate") %>%
  # Calculate fold change relative to the mean WT slope for that substrate
  mutate(fold_change = max_slope / wt_mean_slope)

# --- 3. Calculate Reference Means for V-Lines (in FC) ---
# Calculate the mean 'fold_change' for the 'wt' and 'rfp' colonies
# for each unique 'substrate' to use as vertical reference lines.
target_means_fc <- data_fc %>%
  filter(colony %in% c("wt", "rfp")) %>%
  group_by(substrate, colony) %>%
  summarise(
    mean_fc = mean(fold_change, na.rm = TRUE),
    .groups = 'drop'
  )

# --- 4. Generate Plot ---

slopes_plot <- data_fc %>% # Use the data_fc dataframe with the new column
  ggplot(aes(x = fold_change)) + # Set x-axis to the new 'fold_change'
  
  # A. Create the histogram (using density instead of count for normalization)
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 30,
    fill = "gray50", # Changed bar fill to gray
    color = "white", # Changed bar outline to white
    alpha = 0.8
  ) +
  
  # B. Add the vertical reference lines
  geom_vline(
    data = target_means_fc, # Use the new reference data (in fold change)
    aes(xintercept = mean_fc, color = colony),
    linetype = "solid",
    linewidth = 1.2
  ) +
  
  # C. Facet the plot by 'substrate'
  # scales = "free_y" ensures the y-axis (density) scales independently for each histogram.
  # UPDATED: Set to 3 columns wide and 2 rows tall.
  facet_wrap(~substrate, scales = "free_y", ncol = 3, nrow = 2) +
  
  # D. Manually assign colors as requested (green for 'wt', red for 'rfp')
  scale_color_manual(
    values = c("rfp" = "red", "wt" = "green4"),
    labels = c("rfp" = "RFP Mean FC", "wt" = "WT kcat min^-1"),
    name = "Reference Line"
  ) +
  
  # E. Labels and Theme
  labs(
    title = "Fold Change kcat (min^-1) over WT",
    x = "kcat min^-1 (Fold Change over WT)", # Updated x-axis label
    y = "Density"
  ) +
  # Using theme_minimal as a base, then customizing for dark background
  theme_minimal(base_size = 14) +
  theme(
    # Set entire plot background to transparent
    plot.background = element_rect(fill = "transparent", color = NA),
    # Set panel background (area behind data) to transparent
    panel.background = element_rect(fill = "transparent", color = NA),
    
    # Text colors to white
    text = element_text(color = "white"),
    plot.title = element_text(hjust = 0.5, face = "bold", margin = margin(b = 10), color = "white"),
    
    # Axes, ticks, and grid lines to white/light gray
    axis.text = element_text(color = "white"),
    axis.title = element_text(color = "white"),
    axis.line = element_line(color = "white"),
    panel.grid.major = element_line(color = "gray30"),
    panel.grid.minor = element_line(color = "gray30"),
    
    # Facet strip background to transparent, text to white
    strip.background = element_rect(fill = "transparent", color = "white"),
    strip.text = element_text(face = "bold", color = "white"),
    
    # Legend text and box to white/transparent
    legend.position = "bottom",
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.key = element_rect(fill = "transparent", color = NA)
  )

# Print the plot
pdf("data/20251017/histogram_fc_protein_norm.pdf", width = 12, height = 8, bg = "black")
print(slopes_plot)
invisible(dev.off())

# Print the plot to the console for Canvas preview (optional, but helpful)
print(slopes_plot)
