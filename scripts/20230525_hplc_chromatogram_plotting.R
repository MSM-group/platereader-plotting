# Install packages
rm(list = ls())
pacman::p_load("tidyverse", "readxl", "RColorBrewer", "ggpubr", "plotly","viridis" )

# Read in the data (finds all files to overlay)
fils <- list.files("data/diana_HPLC_data/", pattern = "*.xlsx", full.names = T)
fils
filrem <- fils[!grepl("~|Wasser|MiliQ|parac|LBIIpb|mgL|t4", fils)] # this is where you select patterns of the files you want or don't want
filrem
filcomb <- c(filrem,fils[grepl("\\/LBIIpbt6",fils)])
readin <- tibble(filename = filcomb) %>%
  mutate(file_contents = purrr::map(filename,          
                                    ~ read_excel(file.path(.), skip = 76)) # skips the first 76 line header
  ) %>%
  unnest(.,cols = c(file_contents)) 
readin

dat1 <- readin %>%
  janitor::clean_names() %>%
  dplyr::mutate(sampleraw = word(filename, sep = "diana_HPLC_data\\/", -1))%>%
  mutate(sample=case_when(str_detect(sampleraw, "t1") ~ " enzyme_0.25h",
                          str_detect(sampleraw, "t2") ~ " enzyme_1.10h",
                          str_detect(sampleraw, "t3") ~ " enzyme_17.6h",
                          str_detect(sampleraw, "t5") ~ " enzyme_24.0h",
                          str_detect(sampleraw, "LBIIpt6") ~ " enzyme_40.4h",
                          str_detect(sampleraw, "LBIIpbt6") ~ "boiled_enzyme_40.4h"))%>%
  filter(time_min<10)


pal1 <-  colorRampPalette(brewer.pal(12,"Paired"))(12)
pal2 <- c(magma(6)[1:5], "grey80")

# Plot the data 2D
pl1 <- ggplot(data = dat1,  aes(x = time_min, y = value_m_au, colour = sample)) +
  theme_pubr() +
  xlab("Retention time (min)") +
  ylab("Absorbance (246 nm)") +
 # facet_grid(rows = vars(grp)) +
  geom_point(size = 0.0000000001, alpha = 0.01) +
  geom_path(linewidth = 0.75) +
  scale_color_manual(values = pal2) +
  facet_grid(rows = vars(sample)) +
  xlim(3,10) +
  theme(strip.text.y = element_blank()) 
  #geom_label(data=subset(dat1,time_min>9.95),
            #aes(label=sample))
pl1
ggsave("output/LBamidase_paracetamol_degradation_for_poster.png", plot = pl1)

# Try out a 3D plot
plot_ly(dat1, x = ~sample, y = ~time_min, z = ~value_m_au, type = 'scatter3d', 
        mode = 'lines', color = ~sample)

