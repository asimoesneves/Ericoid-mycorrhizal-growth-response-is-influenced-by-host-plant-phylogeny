library(tidyverse)

################## RESPONSE SPECIFICITY FUNGAL SPECIES ########
#plots for response specificity fungi

create_stacked_barplot <- function(df) {
  # Get the order of species based on positive appearances
  species_order <- df %>%
    arrange(desc(MGR_positive)) %>%
    pull(fungal_sp)
  
  # Reshape data from wide to long format
  df_long <- df %>%
    mutate(fungal_sp = factor(fungal_sp, levels = species_order)) %>%
    pivot_longer(
      cols = starts_with("MGR_"),
      names_to = "appearance_type",
      values_to = "count"
    ) %>%
    mutate(
      appearance_type = str_remove(appearance_type, "MGR_")
    )
  
  # Create the plot
  ggplot(df_long, aes(x = fungal_sp, y = count, fill = appearance_type)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(
      values = c(
        "positive" = "#2ECC71",
        "negative" = "#E74C3C",
        "ns" = "#95A5A6"
      ),
      name = "",
      labels = c(
        "positive" = "MGR > 0",
        "negative" = "MGR < 0",
        "ns" = "MGR Non-significant"
      )
    ) +
    scale_y_continuous(breaks = function(x) seq(from = 0, to = ceiling(max(x)), by = 1)) +  # Integer breaks
    theme_classic(base_size = 16) +  # Increase base font size
    labs(
      title = "Fungal species response specificity",
      x = "",
      y = "Number of plant species"
    ) +
    theme(
      axis.text.x = element_text(angle = 65, hjust = 1, size = 12),
      axis.text.y = element_text(size = 16),
      axis.title = element_text(size = 16),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 1),
      title = element_text(size = 16),
      panel.grid.major.x = element_blank(),
      panel.border = element_blank(),
      legend.position = "right"
    )
}

##### LOAD DATA
### total MGR
dat_total_fungi = read.csv2("response_specificity_total_fungi.csv")
#str(dat_total_fungi)
#dat_total_fungi

### abvg MGR
dat_abvg_fungi = read.csv2("response_specificity_abvg_fungi.csv")
#str(dat_abvg_fungi)
#dat_abvg_fungi

### blwg MGR
dat_blwg_fungi = read.csv2("response_specificity_blwg_fungi.csv")
#str(dat_blwg_fungi)
#dat_blwg_fungi

### CREATE PLOTS
## total MGR
plot_total_fungi <- create_stacked_barplot(dat_total_fungi)
#print(plot_total_fungi)

png("figures/response_specificity/response_specificity_total_fungi.jpg", width = 7, height = 7, units ='in', res = 300)
plot_total_fungi
dev.off()

## abvg MGR
plot_abvg_fungi <- create_stacked_barplot(dat_abvg_fungi)
#print(plot_abvg_fungi)

png("figures/response_specificity/response_specificity_abvg_fungi.jpg", width = 7, height = 7, units ='in', res = 300)
plot_abvg_fungi
dev.off()

## blwg MGR
plot_blwg_fungi <- create_stacked_barplot(dat_blwg_fungi)
#print(plot_blwg_fungi)

png("figures/response_specificity/response_specificity_blwg_fungi.jpg", width = 7, height = 7, units ='in', res = 300)
plot_blwg_fungi
dev.off()

################################################################################
##################### RESPONSE SPECIFICITY FUNGAL SPECIES ######################
################################################################################
#plots for response specificity plant

create_stacked_barplot <- function(df) {
  # Get the order of species based on positive appearances
  species_order <- df %>%
    arrange(desc(MGR_positive)) %>%
    pull(plant_sp)
  
  # Reshape data from wide to long format
  df_long <- df %>%
    mutate(plant_sp = factor(plant_sp, levels = species_order)) %>%
    pivot_longer(
      cols = starts_with("MGR_"),
      names_to = "appearance_type",
      values_to = "count"
    ) %>%
    mutate(
      appearance_type = str_remove(appearance_type, "MGR_")
    )
  
  # Create the plot
  ggplot(df_long, aes(x = plant_sp, y = count, fill = appearance_type)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(
      values = c(
        "positive" = "#2ECC71",
        "negative" = "#E74C3C",
        "ns" = "#95A5A6"
      ),
      name = "",
      labels = c(
        "positive" = "MGR > 0",
        "negative" = "MGR < 0",
        "ns" = "MGR Non-significant"
      )
    ) +
    scale_y_continuous(breaks = function(x) seq(from = 0, to = ceiling(max(x)), by = 1)) +  # Integer breaks
    theme_classic(base_size = 16) +  # Increase base font size
    labs(
      title = "Plant species response specificity",
      x = "",
      y = "Number of fungal species"
    ) +
    theme(
      axis.text.x = element_text(angle = 65, hjust = 1, size = 12),
      axis.text.y = element_text(size = 16),
      axis.title = element_text(size = 16),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 1),
      title = element_text(size = 16),
      panel.grid.major.x = element_blank(),
      panel.border = element_blank(),
      legend.position = "right"
    )
}

##### LOAD DATA
### total MGR
dat_total_plant = read.csv2("response_specificity_total_plant.csv")
#str(dat_total_plant)
#dat_total_plant

### abvg MGR
dat_abvg_plant = read.csv2("response_specificity_abvg_plant.csv")
#str(dat_abvg_plant)
#dat_abvg_plant

### blwg MGR
dat_blwg_plant = read.csv2("response_specificity_blwg_plant.csv")
#str(dat_blwg_plant)
#dat_blwg_plant

### CREATE PLOTS
## total MGR
plot_total_plant <- create_stacked_barplot(dat_total_plant)
#print(plot_total_plant)

png("figures/response_specificity/response_specificity_total_plant.jpg", width = 7, height = 7, units ='in', res = 300)
plot_total_plant
dev.off()

## abvg MGR
plot_abvg_plant <- create_stacked_barplot(dat_abvg_plant)
#print(plot_abvg_plant)

png("figures/response_specificity/response_specificity_abvg_plant.jpg", width = 7, height = 7, units ='in', res = 300)
plot_abvg_plant
dev.off()

## blwg MGR
plot_blwg_plant <- create_stacked_barplot(dat_blwg_plant)
#print(plot_blwg_plant)

png("figures/response_specificity/response_specificity_blwg_plant.jpg", width = 7, height = 7, units ='in', res = 300)
plot_blwg_plant
dev.off()
