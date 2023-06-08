# Install packages
pacman::p_load("tidyverse", "ggpubr", "ggdist", "readxl", "viridis")

# Read in the data

z1 <- read_excel("data/hplc/HPLC_results_paper.xlsx", col_names = F) %>%
  janitor::clean_names() %>%
  dplyr::select(-x8, -x9, -x10) %>%
  dplyr::filter(!is.na(x2))

DF1 <- z1 %>%
  dplyr::slice((grep("DH", z1$x2)[1]+1):(grep("EV", z1$x2)[1]-1)) 
colnames(DF1) <- c(z1[grep("DH", z1$x2)[1],]) 
colnames(DF1)[1] <- "time"
DF1_l <- reshape2::melt(DF1, id.vars = "time")


DF2 <- z1 %>%
  dplyr::slice((grep("EV", z1$x2)[1]+1):(grep("p006", z1$x2)[1]-1))
colnames(DF2) <- z1[grep("EV", z1$x2)[1],]
colnames(DF2)[1] <- "time"
DF2_l <- reshape2::melt(DF2, id.vars = "time")


DF3 <- z1  %>%
  dplyr::slice((grep("DH", z1$x2)[2]+1):(grep("EV", z1$x2)[2]-1))
colnames(DF3) <- z1[grep("DH", z1$x2)[2],]
colnames(DF3)
colnames(DF3)[1] <- "time"
DF3_l <- reshape2::melt(DF3, id.vars = "time")


DF4 <-  z1  %>%
  dplyr::slice((grep("EV", z1$x2)[2]+1):(grep("p006", z1$x2)[2]-1)) 
colnames(DF4) <- z1[grep("EV", z1$x2)[2],]
colnames(DF4)[1] <- "time"
DF4_l <- reshape2::melt(DF4, id.vars = "time")

DF5 <-  z1  %>%
  dplyr::slice((grep("p006", z1$x2)[1]+1):(grep("DH", z1$x2)[2]-1))
colnames(DF5) <- z1[grep("p006", z1$x2)[1],]
colnames(DF5)[1] <- "time"
DF5_l <- reshape2::melt(DF5, id.vars = "time")

DF6 <-  z1  %>%
  dplyr::slice((grep("p006", z1$x2)[2]+1):nrow(z1))
colnames(DF6) <- z1[grep("p006", z1$x2)[2],]
colnames(DF6)[1] <- "time"
DF6_l <- reshape2::melt(DF6, id.vars = "time")


bindall <- bind_rows( DF5_l, DF6_l, DF1_l, DF3_l) #DF2_l, DF4_l,)

# Make the reactants and products
binddf <- bindall %>%
  dplyr::mutate(biorep = as.factor(word(variable, sep = "_", 2))) %>%
  dplyr::mutate(compound = as.factor(word(variable, sep = "_", 3))) %>%
  dplyr::mutate(time_hrs = as.numeric(word(time, sep = "h", 1))) %>%
  dplyr::mutate(conc = as.numeric(trimws(word(value, sep = " ", 1)))) %>%
  dplyr::mutate(enzyme = as.factor(word(variable, sep = "_", 1))) %>%
  dplyr::mutate(idvar = paste0(enzyme, "_", time, "_", compound)) %>%
  dplyr::filter(enzyme != "EV") %>%
  dplyr::filter(variable != "p006_2_MB")
 
summdf <- binddf %>%
  dplyr::group_by(idvar) %>%
  dplyr::mutate(mean = mean(conc, na.rm = T)) %>%
  dplyr::mutate(sd = sd(conc, na.rm = T)) %>%
  ungroup()
summdf$conc
summdf$mean_conc

#  Plot
pal3 <- c(magma(6)[c(1, 3:5)])

ggplot(summdf, aes(x = time_hrs, y = mean)) + 
  geom_point(aes(color = compound)) +
  geom_errorbar(aes(color = compound, ymax=mean + sd, ymin = mean - sd), width=0.3,size=0.6) +
  geom_line(aes(color = compound)) +
  facet_grid(rows = vars(enzyme)) +
  theme_pubr() +
  scale_color_manual(values = pal3) +
  ylab("Concentration (µM)") +
  xlab("Time (hours)")

ggplot(summdf, aes(x = time_hrs, y = mean)) + 
  geom_point(aes(color = compound)) +
  geom_errorbar(aes(color = compound, ymax=mean + sd, ymin = mean - sd), width=0.3,size=0.6) +
  geom_line(aes(color = compound)) +
  facet_grid(rows = vars(enzyme)) +
  theme_pubr() +
  scale_color_manual(values = pal3) +
  ylab("Concentration (µM)") +
  xlab("Time (hours)")

pl3 <- ggplot(summdf, aes(x = time_hrs, y = mean)) + 
  geom_point(aes(color = compound)) +
  geom_errorbar(aes(color = compound, ymax=mean + sd, ymin = mean - sd), width=0.3,size=0.6) +
  geom_line(aes(color = compound)) +
  facet_grid(cols = vars(enzyme)) +
  theme_pubr() +
  scale_color_manual(values = pal3) +
  ylab("Concentration (µM)") +
  xlab("Time (hours)")
pl3
ggsave(pl3, filename = "output/hplc_sample_plot.png", width = 7.5, heigh = 5)

binddf <- bindall %>%
  #dplyr::filter(str_detect("EV" variable)) %>%
  dplyr::mutate(biorep = as.factor(word(variable, sep = "_", 2))) %>%
  dplyr::mutate(compound = as.factor(word(variable, sep = "_", 3))) %>%
  dplyr::mutate(time_hrs = as.numeric(word(time, sep = "h", 1))) %>%
  dplyr::mutate(conc = as.numeric(trimws(word(value, sep = " ", 1)))) %>%
  dplyr::mutate(enzyme = as.factor(word(variable, sep = "_", 1))) %>%
  dplyr::mutate(idvar = paste0(enzyme, "_", time, "_", compound)) %>%
  dplyr::filter(variable != "p006_2_MB") %>%
  dplyr::filter(!grepl("DH_24h", idvar))
binddf$enzyme
binddf$idvar

summdf <- binddf %>%
  dplyr::group_by(idvar) %>%
  dplyr::mutate(mean = mean(conc, na.rm = T)) %>%
  dplyr::mutate(sd = sd(conc, na.rm = T)) %>%
  ungroup()

pl4 <- ggplot(summdf, aes(x = time_hrs, y = mean)) + 
  geom_point(aes(color = compound)) +
  geom_errorbar(aes(color = compound, ymax=mean + sd, ymin = mean - sd), width=0.3,size=0.6) +
  geom_line(aes(color = compound)) +
  facet_grid(cols = vars(enzyme), scales = 'free') +
  theme_pubr() +
  scale_color_manual(values = pal3) +
  ylab("Concentration (µM)") +
  xlab("Time (hours)")
pl4
ggsave(pl4, filename = "output/DEET_zoom_hplc.png", width = 6.5, height = 4)
