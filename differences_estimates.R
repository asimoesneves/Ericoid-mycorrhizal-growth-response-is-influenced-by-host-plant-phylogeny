library(tidyverse)
library(ggplot2)
library(ggsci)
library(RColorBrewer)
library(lme4)
library(lmtest)
library(lmerTest)
library(sjPlot)


dat_og = read.csv2("measurements_21stFebruary_Version.csv")
str(dat_og)

dat_og = dat_og[1:714, 1:9]

dat_og =dat_og %>%select(-c("Height.t1"))

dat2 = read.csv2("Harvested_Plants_june.csv")


str(dat2)

dat2 = dat2[1:643, 1:19]

dat2 = dat2 %>% select(-c("plant_sp", "fungal_sp", "replicate", "treatment",
                          "block", "Height.t0", "t1.t0", "Dead.Sick", "X.cassette"))

str(dat2)

dat = dat_og %>%
  left_join(dat2, by = "unique_id")

str(dat)

##### CHANGE PLANT SP NAMES FOR LONG NAME

dat <- transform(dat,plant_sp=gsub(pattern="CVU", replacement="C_vulgaris", plant_sp))
dat <- transform(dat,plant_sp=gsub(pattern="GSH", replacement="G_shallon", plant_sp))
dat <- transform(dat,plant_sp=gsub(pattern="KLA", replacement="K_latifolia", plant_sp))
dat <- transform(dat,plant_sp=gsub(pattern="PJA", replacement="P_japonica", plant_sp))
dat <- transform(dat,plant_sp=gsub(pattern="RAR", replacement="R_arboreum", plant_sp))
dat <- transform(dat,plant_sp=gsub(pattern="RFE", replacement="R_ferrugineum", plant_sp))
dat <- transform(dat,plant_sp=gsub(pattern="VAN", replacement="V_angustifolium", plant_sp))
dat <- transform(dat,plant_sp=gsub(pattern="VMY", replacement="V_myrtillus", plant_sp))
dat <- transform(dat,plant_sp=gsub(pattern="VVI", replacement="V_vitis-idaea", plant_sp))


##### CHANGE FUNGAL SP NAMES FOR LONG NAME

dat <- transform(dat,fungal_sp=gsub(pattern="JPK-132", replacement="Serendipitaceae_sp", fungal_sp))
dat <- transform(dat,fungal_sp=gsub(pattern="JPK-87", replacement="K_argillacea", fungal_sp))
dat <- transform(dat,fungal_sp=gsub(pattern="MEB", replacement="H_bicolor", fungal_sp))
dat <- transform(dat,fungal_sp=gsub(pattern="MEV", replacement="H_variabilis", fungal_sp))
dat <- transform(dat,fungal_sp=gsub(pattern="MGK", replacement="H_gryndleri", fungal_sp))
dat <- transform(dat,fungal_sp=gsub(pattern="OMA", replacement="O_maius", fungal_sp))
dat <- transform(dat,fungal_sp=gsub(pattern="RER", replacement="H_hepaticicola_1", fungal_sp))
dat <- transform(dat,fungal_sp=gsub(pattern="RHE", replacement="H_hepaticicola_2", fungal_sp))
dat <- transform(dat,fungal_sp=gsub(pattern="STR", replacement="Sterile", fungal_sp))


dat = dat %>% select(-c("replicate", "treatment", "block", "Height.t0",
                          "Dead.Sick", "unique_id", "Height.t1", "blockx", "shoot_weight"))
####

##calculate average wet/dry biomass ratio for each plant_sp

dat1 <- data.frame()

# add real values of wet and dry weight from experiment to the table
#need to do in 2 steps: first until row 243 - where largesubsample was done

ratio_real_values <- filter(dat, random_label<=243)
ratio_real_values = ratio_real_values %>% select(-c("fungal_sp", "total_root_wet", "random_label",
                                                    "small_subsample_wet", "total_root_dry"))

colnames(ratio_real_values)[colnames(ratio_real_values) == c("plant_sp", "large_subsample_wet","large_subsample_dry") ] <- c("Plant.sp", "weight.wet.root", "weight.dry.root")
dat1 <- rbind(dat1, ratio_real_values)

# and then from 244 until the end with total_root_wet and total_root_dry
ratio_real_values1 <- filter(dat, random_label>243)
ratio_real_values1 = ratio_real_values1 %>% select(-c("fungal_sp", "random_label",
                                                      "large_subsample_wet", "small_subsample_wet",
                                                      "large_subsample_dry"))

colnames(ratio_real_values1)[colnames(ratio_real_values1) == c("plant_sp", "total_root_wet","total_root_dry") ] <- c("Plant.sp", "weight.wet.root", "weight.dry.root")
dat1 <- rbind(dat1, ratio_real_values1)

#group by
#mean (in tidyverse)
dat1 = dat1 %>% drop_na()



dat1$ratio = dat1$weight.dry.root/dat1$weight.wet.root

xx=dat1%>%
  group_by(Plant.sp)%>%
  summarize(Mean = mean(ratio))

ratio_table = data.frame(xx)
colnames(ratio_table) <- c("plant_sp", "ratio")

dat1 = dat1 %>% left_join(xx)
dat1 = dat1 %>% select(-c("ratio"))


dat1$weight.dry.estimated = dat1$weight.wet.root*dat1$Mean


dat1$dif.est = dat1$weight.dry.root-dat1$weight.dry.estimated

df <- dat1 %>% slice(-146)


library(ggpubr)

dat1$Plant.sp <- factor(dat1$Plant.sp,levels = c("G_shallon", "V_myrtillus", "P_japonica",
                                             "V_vitis-idaea", "C_vulgaris", "V_angustifolium",
                                             "K_latifolia", "R_arboreum", "R_ferrugineum"))

png("figures/difference_real_values_and_estimates_updated_all rows.jpg", width = 17, height = 3.5, units ='in', res = 300)
ggplot(data=dat1, mapping=aes(x=weight.dry.root, y=weight.dry.estimated)) +
  facet_grid(cols = vars(Plant.sp)) +
  scale_color_brewer(palette = "Paired") +
  geom_point(alpha=0.3) +
  geom_smooth(method="lm", linewidth = 0.5, alpha=0.1) +
  theme_classic(base_size = 15) +
  theme(
    panel.background = element_rect(fill = NA, color = "grey", size = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "transparent", colour = NA_character_),
    strip.background = element_rect(fill = "grey90", color = "grey", size = 0.5),
    legend.background = element_rect(fill = "transparent"),
    legend.box.background = element_rect(fill = "transparent"),
    legend.key = element_rect(fill = "transparent"),
    axis.title.x = element_text(size = 16),
    axis.line = element_line(size = 0.5),
    strip.text = element_text(size = 15),
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title.y = element_text(size = 16),
    legend.title = element_text(size = 17),
    legend.text = element_text(size = 16)
  ) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = NA, color = "grey"),
        panel.spacing = unit(1, "lines"))+
  stat_cor(aes(label = after_stat(rr.label)), color = "black", geom = "label") +
  labs(y= "estimated dry weight", x = "real dry weight") + 
  theme(aspect.ratio = 1.2)
dev.off()

## without outlier

df$Plant.sp <- factor(df$Plant.sp,levels = c("G_shallon", "V_myrtillus", "P_japonica",
                                                     "V_vitis-idaea", "C_vulgaris", "V_angustifolium",
                                                     "K_latifolia", "R_arboreum", "R_ferrugineum"))

png("figures/difference_real_values_and_estimates_updated.jpg", width = 17, height = 3.5, units ='in', res = 300)
ggplot(data=df, mapping=aes(x=weight.dry.root, y=weight.dry.estimated)) +
  facet_grid(cols = vars(Plant.sp)) +
  scale_color_brewer(palette = "Paired") +
  geom_point(alpha=0.3) +
  geom_smooth(method="lm", linewidth = 0.5, alpha=0.1) +
  theme_classic(base_size = 15) +
  theme(
    panel.background = element_rect(fill = NA, color = "grey", size = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "transparent", colour = NA_character_),
    strip.background = element_rect(fill = "grey90", color = "grey", size = 0.5),
    legend.background = element_rect(fill = "transparent"),
    legend.box.background = element_rect(fill = "transparent"),
    legend.key = element_rect(fill = "transparent"),
    axis.title.x = element_text(size = 16),
    axis.line = element_line(size = 0.5),
    strip.text = element_text(size = 15),
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title.y = element_text(size = 16),
    legend.title = element_text(size = 17),
    legend.text = element_text(size = 16)
  ) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = NA, color = "grey"),
        panel.spacing = unit(1, "lines"))+
  stat_cor(aes(label = after_stat(rr.label)), color = "black", geom = "label") +
  labs(y= "estimated dry weight", x = "real dry weight") + 
  theme(aspect.ratio = 1.2)
dev.off()

