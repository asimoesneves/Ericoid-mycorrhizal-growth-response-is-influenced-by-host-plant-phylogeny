###### TEST FOR MODELS

library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(lme4)
library(nmixgof)
library(DHARMa)
library(emmeans)
library(glmmTMB)
library(lmtest)
library(Matrix)
library(emmeans)
library(lsmeans)
library(tidyr)
library(lmerTest)
library(ggsignif)
library(multcompView)
library(pbkrtest)
library(forcats)
library(sjPlot)
library(sjmisc)
library(dplyr)
library(stringr)



####
#homogeneity function
resid_plot_fit <- function(mod,dat) {
  E1 <- resid(mod, type = "pearson")
  N <- nrow(dat)
  p <- length(coef(dat))
  sum(E1^2)/(N-p)
  F1 <- fitted(mod)
  ggplot() +
    aes(F1, E1) +
    geom_point() +
    ggtitle("Residual plot") +
    xlab("Fitted Values") +
    ylab("Pearson residuals") 
}

### COLONIZATION ASSESSMENT for coils!

# import csv with intersections + sum all intersections of each sample
col_ass = read.csv2("Colonization_Assessment1.csv")
col_ass = col_ass[,1:8]
str(col_ass)
summary(col_ass)

col_ass = col_ass %>% drop_na()

# dataframe for coils
coils_ass = col_ass %>% 
  group_by(random_label) %>% 
  summarize(colonized = sum(Coils),
            intersections_c = sum(roots_crossed))

coils_ass$non_colonized = coils_ass$intersections_c-coils_ass$colonized
coils_ass$percent_colonization = coils_ass$colonized/coils_ass$intersections_c
str(coils_ass)
summary(coils_ass)

#dataframe for hyphae
hyphae_ass = col_ass %>% 
  group_by(random_label) %>% 
  summarize(colonized_h = sum(Hyphae),
            intersections_h = sum(roots_crossed))

hyphae_ass$non_colonized_h = hyphae_ass$intersections_h-hyphae_ass$colonized_h
hyphae_ass$percent_colonization_h = hyphae_ass$colonized_h/hyphae_ass$intersections_h
str(hyphae_ass)
summary(hyphae_ass)

#dataframe for conidiophores
spikes_ass = col_ass %>% 
  group_by(random_label) %>% 
  summarize(colonized_s = sum(Spikes),
            intersections_s = sum(roots_crossed))

spikes_ass$non_colonized_s = spikes_ass$intersections_s-spikes_ass$colonized_s
spikes_ass$percent_colonization_s = spikes_ass$colonized_s/spikes_ass$intersections_s
str(spikes_ass)
summary(spikes_ass)

# import master spread sheet and see just colonization with coils
dat = read.csv2("Harvested_Plants_Processed_try1_updated.csv")

dat = dat %>% select(-c("replicate", "treatment", "blockx", "above_below_ratio",
                        "Dead.Sick", "Height.t1", "abvg_biomass_STR", "blwg_biomass_STR",
                        "total_root_wet", "small_subsample_wet", "large_subsample_wet",
                        "shoot_weight", "ratio", "biomass_STR",
                        "total_biomass","large_subsample_dry",
                        "total_root_dry", "small_subsample_dry", "ratio_STR", "MGR_ratio"))

dat = dat%>%left_join(coils_ass)
dat = dat%>%left_join(hyphae_ass)
dat = dat%>%left_join(spikes_ass)
dat = dat %>% drop_na(percent_colonization)

#percentage of samples colonized
m = dat %>% count(percent_colonization!=0)
percentage = m[2, 2]/(m[1,2]+m[2, 2])*100

m1 = dat %>% count(percent_colonization_h!=0)
percentage1 = m1[2, 2]/(m1[1,2]+m1[2, 2])*100

m2 = dat %>% count(percent_colonization_s!=0)
percentage2 = m2[2, 2]/(m2[1,2]+m2[2, 2])*100

#################################################################

## 2. str vs plant_sp*fungal_sp

mod = glmmTMB(cbind(colonized, non_colonized) ~fungal_sp*plant_sp,
             family = binomial, data=dat)

# check assumptions --> ??
resid_plot_fit(mod,dat)
testZeroInflation(mod)
testDispersion(mod)

#save model output
# Extract fixed effects coefficients and statistics
fixed_effects <- as.data.frame(summary(mod)$coefficients$cond)


# Write to CSV
write.csv2(fixed_effects, "C:/Users/Dpao/Desktop/Master's Biology/Master Thesis/#2 semester project/Ericoid Project/colonization_model_output_5thfeb.csv", row.names=TRUE)


ref1 = emmeans::emmeans(mod, ~fungal_sp*plant_sp, type="response")
ref1
summary(ref1)

#Add p-value and statistical tests
ref1 = update(ref1, infer = c(TRUE, TRUE), null = log(35), type="response",
               calc = c(n = ".wgt."))
summary(ref1)

ref.table1 = as.data.frame(ref1)

# Write to CSV
write.csv2(ref.table1, "C:/Users/Dpao/Desktop/Master's Biology/Master Thesis/#2 semester project/Ericoid Project/colonization_emmeans_output_5thfeb.csv", row.names=TRUE)

ref.table1$fungal_sp <- gsub(fixed("Sterile"), "Non-inoculated", ref.table1$fungal_sp)

ref.table1$plant_sp <- factor(ref.table1$plant_sp,levels = c("G. shallon", "V. myrtillus", "P. japonica",
                                               "V. vitis-idaea", "C. vulgaris", "V. angustifolium",
                                               "K. latifolia", "R. arboreum", "R. ferrugineum"))


# Create color palette with ordered species
library(RColorBrewer)

# colors_definition.R
fungal_colors <- c(
  "Non-inoculated" = "#66C2A5",
  "H. gryndleri" = "#FC8D62",
  "H. hepaticicola_PK 135-3" = "#8DA0CB",
  "H. hepaticicola_UAMH7357/ICMP18553" = "#E78AC3",
  "K. argillacea" = "#A6D854",
  "H. bicolor" = "#FFD92F",
  "H. variabilis" = "#E5C494",
  "O. maius" = "#B3B3B3",
  "Serendipitaceae_sp" = "#7570B3"
)

saveRDS(fungal_colors, "fungal_colors.rds")

ref.table1$fungal_sp <- factor(ref.table1$fungal_sp,
                               levels = c("Non-inoculated", "H. gryndleri", "H. hepaticicola_PK 135-3", "H. hepaticicola_UAMH7357/ICMP18553",
                                          "K. argillacea", "H. bicolor", "H. variabilis",
                                          "O. maius", "Serendipitaceae_sp"))

png("figures/colonization/col_str_treatment_update.jpg", width = 18, height = 5.4, units ='in', res = 300)
  ggplot(ref.table1, aes(fungal_sp, prob, color = fungal_sp, fill = fungal_sp)) +  
  facet_grid(cols = vars(plant_sp)) +
  geom_bar(stat="identity", alpha = 0.08, width=0.8) +
  geom_point(size = 1.2) +
  geom_errorbar(aes(ymin=prob-asymp.LCL, ymax=prob+asymp.UCL), width=.65, linewidth=0.6, position=position_dodge(1)) +
  scale_fill_manual(values = fungal_colors) +
  scale_color_manual(values = fungal_colors) +
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
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line = element_line(size = 0.5),
      strip.text = element_text(size = 14),
      plot.title = element_text(hjust = 0.5, size = 16),
      axis.title.y = element_text(size = 16),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 13)
    ) +
  ylab("Root colonization rate") +
  geom_hline(yintercept=0, linetype='dotted', col = "black", size=0.1)
dev.off()



tab_results = pairs(ref1)

makeStars <- function(x){
  stars <- c("****", "***", "**", "*", "ns")
  vec <- c(0, 0.0001, 0.001, 0.01, 0.05, 1.01)
  i <- findInterval(x, vec)
  stars[i]
}


tab_results = as.data.frame(tab_results)

tab_results = separate(tab_results, col=contrast, into=c('comb1', 'comb2'), sep=' / ')
tab_results$p_value_stars <- makeStars(tab_results$p.value)

#remove "(" and ")" from some of the cells
tab_results$comb1 <- gsub("\\(", "", tab_results$comb1)
tab_results$comb2 <- gsub("\\(", "", tab_results$comb2)
tab_results$comb1 <- gsub("\\)", "", tab_results$comb1)
tab_results$comb2 <- gsub("\\)", "", tab_results$comb2)


#separate comb1 into fungi1 and plant1 and comb2 into fungi2 plant2

# split the combination based on the plant name pattern
split_combination <- function(combination) {
  # Create a pattern with all possible plant names
  # Escape special characters like "-" and "."
  plant_patterns <- c(
    "G\\. shallon", "V\\. myrtillus", "P\\. japonica",
    "V\\. vitis-idaea", "C\\. vulgaris", "V\\. angustifolium",
    "K\\. latifolia", "R\\. arboreum", "R\\. ferrugineum"
  )
  
  # Create a regex pattern that matches any of the plant names
  plant_pattern <- paste0("(", paste(plant_patterns, collapse="|"), ")$")
  
  # Extract the plant name
  plant_name <- str_extract(combination, plant_pattern)
  
  # Extract the fungi name by removing the plant name
  fungi_name <- str_replace(combination, paste0(" ", plant_pattern), "")
  
  return(list(fungi = fungi_name, plant = plant_name))
}


tab_results = tab_results %>%
  mutate(
    comb1_split = map(comb1, split_combination),
    comb2_split = map(comb2, split_combination),
    
    fungi1 = map_chr(comb1_split, "fungi"),
    plant1 = map_chr(comb1_split, "plant"),
    fungi2 = map_chr(comb2_split, "fungi"),
    plant2 = map_chr(comb2_split, "plant")
  ) %>%
  select(-comb1_split, -comb2_split)  # Remove the intermediate columns

tab_results = tab_results %>%
  select(-comb1, -comb2)



# keep rows with "sterile" and plant1=plant2
filtered_results <- tab_results[
  (tab_results$fungi1 == "Sterile" | tab_results$fungi2 == "Sterile") & 
    (tab_results$plant1 == tab_results$plant2),
]


filtered_results <- filtered_results[order(filtered_results$plant1), ]
summary_table <- table(filtered_results$plant1)
print(summary_table)

write.csv2(filtered_results, "C:/Users/Dpao/Desktop/Master's Biology/Master Thesis/#2 semester project/Ericoid Project/colonization_pairs_3rdmar.csv", row.names=FALSE)


##### colonization and MGR

# make categorical
dat$fungal_sp <- to_factor(dat$fungal_sp)
dat$plant_sp <- to_factor(dat$plant_sp)

dat1 = subset(dat, dat[,"fungal_sp"]!="Sterile")
dat1 = dat1 %>% drop_na(MGR)
  
m = min(dat1$MGR)
n = m*-1+0.01


#tried with:
colm = lmer(log(MGR+n)~plant_sp*percent_colonization*fungal_sp +   Height.t0 + (1|block), data=dat1)
colm1 = lmer(log(MGR+n)~plant_sp*percent_colonization + fungal_sp*percent_colonization +
              plant_sp*fungal_sp + Height.t0 + (1|block), data=dat1)

colm2 = lmer(log(MGR+n)~fungal_sp*percent_colonization + plant_sp + Height.t0 + (1|block), data=dat1)
colm3 = lmer(log(MGR+n)~plant_sp*percent_colonization + fungal_sp + Height.t0 + (1|block), data=dat1)

colm4 = lmer(log(MGR+n)~plant_sp + percent_colonization + fungal_sp + Height.t0 + (1|block), data=dat1)

summary(colm4)
a = anova(colm4)
plot(colm4)

#aboveground
o = min(dat1$MGR_abvg)
p = o*-1+0.01
colm = lmer(log(MGR_abvg+p)~plant_sp + percent_colonization + fungal_sp + Height.t0 + (1|block), data=dat1)
summary(colm)
b = anova(colm)
plot(colm)


#belowground
q = min(dat1$MGR_blwg)
r = q*-1+0.01
colmx = lmer(log(MGR_blwg+r)~plant_sp + percent_colonization + fungal_sp + Height.t0 + (1|block), data=dat1)
summary(colmx)
c = anova(colmx)
plot(colmx)


write.csv2(a, "C:/Users/Dpao/Desktop/Master's Biology/Master Thesis/#2 semester project/Ericoid Project/MGRtotal_col_updated.csv", row.names=TRUE)
write.csv2(b, "C:/Users/Dpao/Desktop/Master's Biology/Master Thesis/#2 semester project/Ericoid Project/MGRabvg_col_updated.csv", row.names=TRUE)
write.csv2(c, "C:/Users/Dpao/Desktop/Master's Biology/Master Thesis/#2 semester project/Ericoid Project/MGRblwg_col_updated.csv", row.names=TRUE)

#### colonization and MGR (hyphal coils) - subset of plant-fungi combinations
## that had significantly higher colonization than control
#from filtered_results take combinations that are not significantly different
combinations_contaminated = filtered_results %>%
  select (-c(odds.ratio, SE, df, null, z.ratio, p.value, fungi2, plant2))

combinations_contaminated <- combinations_contaminated %>%
  filter(!(p_value_stars!= "ns")) %>%
  rename(fungal_sp = fungi1, plant_sp = plant1) %>%
  select(-c(p_value_stars))

###from dat1 remove the plant-fungi combinations that are present in combinations_contaminated
dat_subset <- dat1 %>%
  anti_join(combinations_contaminated, by = c("plant_sp", "fungal_sp"))

ms = min(dat_subset$MGR)
ns = ms*-1+0.01


#tried with:
colms = lmer(log(MGR+ns)~plant_sp*percent_colonization*fungal_sp +   Height.t0 + (1|block), data=dat_subset)
colm1s = lmer(log(MGR+ns)~plant_sp*percent_colonization + fungal_sp*percent_colonization +
               plant_sp*fungal_sp + Height.t0 + (1|block), data=dat_subset)

colm2s = lmer(log(MGR+ns)~fungal_sp*percent_colonization + plant_sp + Height.t0 + (1|block), data=dat_subset)
colm3s = lmer(log(MGR+ns)~plant_sp*percent_colonization + fungal_sp + Height.t0 + (1|block), data=dat_subset)

colm4s = lmer(log(MGR+ns)~plant_sp + percent_colonization + fungal_sp + Height.t0 + (1|block), data=dat_subset)

summary(colm4s)
as = anova(colm4s)
plot(colm4s)

#aboveground
os = min(dat_subset$MGR_abvg)
ps = os*-1+0.01
colms = lmer(log(MGR_abvg+ps)~plant_sp + percent_colonization + fungal_sp + Height.t0 + (1|block), data=dat_subset)
summary(colms)
bs = anova(colms)
plot(colms)


#belowground
qs = min(dat_subset$MGR_blwg)
rs = qs*-1+0.01
colmxs = lmer(log(MGR_blwg+rs)~plant_sp + percent_colonization + fungal_sp + Height.t0 + (1|block), data=dat_subset)
summary(colmxs)
cs = anova(colmxs)
plot(colmxs)



write.csv2(as, "C:/Users/Dpao/Desktop/Master's Biology/Master Thesis/#2 semester project/Ericoid Project/MGRtotal_col_subset.csv", row.names=TRUE)
write.csv2(bs, "C:/Users/Dpao/Desktop/Master's Biology/Master Thesis/#2 semester project/Ericoid Project/MGRabvg_col_subset.csv", row.names=TRUE)
write.csv2(cs, "C:/Users/Dpao/Desktop/Master's Biology/Master Thesis/#2 semester project/Ericoid Project/MGRblwg_col_subset.csv", row.names=TRUE)



##### colonization and MGR - hyphae

#total biomass
m = min(dat1$MGR)
n = m*-1+0.01


colm_h1 = lmer(log(MGR+n)~plant_sp + percent_colonization_h + fungal_sp + Height.t0 + (1|block), data=dat1)

summary(colm_h1)
d = anova(colm_h1)
plot(colm_h1)

#aboveground
o = min(dat1$MGR_abvg)
p = o*-1+0.01
colm_h2 = lmer(log(MGR_abvg+p)~plant_sp + percent_colonization_h + fungal_sp + Height.t0 + (1|block), data=dat1)
summary(colm_h2)
e = anova(colm_h2)
plot(colm_h2)


#belowground
q = min(dat1$MGR_blwg)
r = q*-1+0.01
colm_h3 = lmer(log(MGR_blwg+r)~plant_sp + percent_colonization_h + fungal_sp + Height.t0 + (1|block), data=dat1)
summary(colm_h3)
f = anova(colm_h3)
plot(colm_h3)


write.csv2(d, "C:/Users/Dpao/Desktop/Master's Biology/Master Thesis/#2 semester project/Ericoid Project/MGRtotal_hyphae_updated.csv", row.names=TRUE)
write.csv2(e, "C:/Users/Dpao/Desktop/Master's Biology/Master Thesis/#2 semester project/Ericoid Project/MGRabvg_hyphae_updated.csv", row.names=TRUE)
write.csv2(f, "C:/Users/Dpao/Desktop/Master's Biology/Master Thesis/#2 semester project/Ericoid Project/MGRblwg_hyphae_updated.csv", row.names=TRUE)


##### colonization and MGR - conidiophores

#total biomass
m = min(dat1$MGR)
n = m*-1+0.01


colm_c1 = lmer(log(MGR+n)~plant_sp + percent_colonization_s + fungal_sp + Height.t0 + (1|block), data=dat1)

summary(colm_c1)
g = anova(colm_c1)
plot(colm_c1)

#aboveground
o = min(dat1$MGR_abvg)
p = o*-1+0.01
colm_c2 = lmer(log(MGR_abvg+p)~plant_sp + percent_colonization_s + fungal_sp + Height.t0 + (1|block), data=dat1)
summary(colm_c2)
h = anova(colm_c2)
plot(colm_c2)


#belowground
q = min(dat1$MGR_blwg)
r = q*-1+0.01
colm_c3 = lmer(log(MGR_blwg+r)~plant_sp + percent_colonization_s + fungal_sp + Height.t0 + (1|block), data=dat1)
summary(colm_c3)
i = anova(colm_c3)
plot(colm_c3)


write.csv2(g, "C:/Users/Dpao/Desktop/Master's Biology/Master Thesis/#2 semester project/Ericoid Project/MGRtotal_conidiophores_updated.csv", row.names=TRUE)
write.csv2(h, "C:/Users/Dpao/Desktop/Master's Biology/Master Thesis/#2 semester project/Ericoid Project/MGRabvg_conidiophores_updated.csv", row.names=TRUE)
write.csv2(i, "C:/Users/Dpao/Desktop/Master's Biology/Master Thesis/#2 semester project/Ericoid Project/MGRblwg_conidiophores_updated.csv", row.names=TRUE)




#######################################################################
######################### HYPHAE & SPIKES #############################
#######################################################################

# models didn't work, so plot real values in boxplot
#################################################################
dat$fungal_sp <- gsub(fixed("Sterile"), "Non-inoculated", dat$fungal_sp)


dat$plant_sp <- factor(dat$plant_sp,levels = c("G. shallon", "V. myrtillus", "P. japonica",
                                               "V. vitis-idaea", "C. vulgaris", "V. angustifolium",
                                               "K. latifolia", "R. arboreum", "R. ferrugineum"))


dat$fungal_sp <- factor(dat$fungal_sp,levels = c("Non-inoculated", "H. gryndleri", "H. hepaticicola_PK 135-3", "H. hepaticicola_UAMH7357/ICMP18553",
                                                 "K. argillacea", "H. bicolor", "H. variabilis",
                                                 "O. maius", "Serendipitaceae_sp"))

#plot hyphal colonization
png("figures/colonization/hyphae_col_assessment_updated.jpg", width = 18, height = 5.4, units ='in', res = 300)
ggplot(dat, aes(fungal_sp, percent_colonization_h, color = fungal_sp, fill = fungal_sp)) +  
  facet_grid(cols = vars(plant_sp)) +
  geom_boxplot(alpha = 0.1, width=0.4) +
  scale_fill_manual(values = fungal_colors) +
  scale_color_manual(values = fungal_colors) +
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
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line = element_line(size = 0.5),
    strip.text = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title.y = element_text(size = 16),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 13)
  ) +
  ylab("Root colonization rate - Hyphae") +
  theme(axis.title.x=element_blank())
dev.off()

#plot conidiophores
png("figures/colonization/conidiophores_col_assessment_updated.jpg", width = 18, height = 5.4, units ='in', res = 300)
ggplot(dat, aes(fungal_sp, percent_colonization_s, color = fungal_sp, fill = fungal_sp)) +  
  facet_grid(cols = vars(plant_sp)) +
  geom_boxplot(alpha = 0.1, width=0.4) +
  scale_fill_manual(values = fungal_colors) +
  scale_color_manual(values = fungal_colors) +
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
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line = element_line(size = 0.5),
    strip.text = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title.y = element_text(size = 16),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 13)
  ) +
  ylab("Root colonization rate - Conidiophores") +
  theme(axis.title.x=element_blank())
dev.off()
