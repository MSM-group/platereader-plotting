# Install packages
rm(list = ls())
pacman::p_load("tidyverse", "readxl", "RColorBrewer", "ggpubr", "plotly")

# Read in the data (finds all files to overlay)
fils <- list.files("data/diana_HPLC_data/", pattern = "*.xlsx", full.names = T)
filrem <- fils[!grepl("~", fils)] # this is where you select patterns of the files you want or don't want

readin <- tibble(filename = filrem) %>%
  mutate(file_contents = purrr::map(filename,          
                                    ~ read_excel(file.path(.), skip = 76)) # skips the first 76 line header
  ) %>%
  unnest(.,cols = c(file_contents)) 
readin

dat1 <- readin %>%
  janitor::clean_names() %>%
  dplyr::mutate(conc = word(filename, sep = "\\/\\/", -1))


pal1 <-  colorRampPalette(brewer.pal(12,"Paired"))(6)

# Plot the data statically
pl1 <- ggplot(data = dat1,  aes(x = time_min, y = value_m_au, colour = conc)) +
  theme_pubr() +
  xlab("Retention time (min)") +
  ylab("Absorbance (246 nm)") +
 # facet_grid(rows = vars(grp)) +
  geom_point(size = 0.0000000001, alpha = 0.01) +
  geom_path(linewidth = 0.75) +
  scale_color_manual(values = pal1) +
  facet_grid(rows = vars(conc))
pl1


# Try out a 3D plot
plot_ly(dat1, x = ~conc, y = ~time_min, z = ~value_m_au, type = 'scatter3d', 
        mode = 'lines', color = ~conc)
