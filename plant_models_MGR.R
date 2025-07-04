library(tidyverse)
library(lme4)
library(lmtest)
library(Matrix)
library(emmeans)
library(ggplot2)
library(lsmeans)
library(tidyr)
library(lmerTest)
library(ggsignif)
library(multcompView)
library(pbkrtest)
library(forcats)
library(egg)

################################################################################
########### FUNCTIONS NEEDED IN THIS SCRIPT

# Load color data
fungal_colors <- readRDS("fungal_colors.rds")

#make stars for p-value
makeStars <- function(x){
  stars <- c("****", "***", "**", "*", "ns")
  vec <- c(0, 0.0001, 0.001, 0.01, 0.05, 1.01)
  i <- findInterval(x, vec)
  stars[i]
}

################################################################################
##using measurements.csv - height of sick-looking plants was deleted
######to give the correct number of samples in the plots
#read data

dat = read.csv2("Harvested_Plants_Processed_try1_updated.csv")
str(dat)

dat$plant_sp=as.factor(dat$plant_sp)
dat$fungal_sp=as.factor(dat$fungal_sp)

summary(dat)
str(dat)

#removing columns from master spreadsheet that are not important here
# keeping: 
dat = dat%>%
  select(-c("Dead.Sick", "Height.t1", "blockx", "total_root_wet", "small_subsample_wet",
            "large_subsample_wet", "large_subsample_dry", "shoot_weight",
            "ratio", "small_subsample_dry", "total_root_dry", "abvg_biomass_STR",
            "blwg_biomass_STR", "total_biomass", "biomass_STR", "ratio_STR"))
str(dat)
summary(dat)
dat = dat %>% drop_na()

summary(dat)
str(dat)


#table w/out STR treatment
dat1 = subset(dat, dat[,"fungal_sp"]!="Sterile")

##############################################################################
#############################ABOVEGROUND######################################
##############################################################################

#1. model
#2. emmeans with plant*fungi interaction and only plant
#3. plot MGR
#4. pairwise comparison

################ model ################
m = min(dat1$MGR_abvg)
n = m*-1+0.01

bm1 = lmer(log(MGR_abvg+n)~fungal_sp * plant_sp + Height.t0 + (1|block), data=dat1)
summary(bm1)
anova(bm1)
plot(bm1)


#save model output, extract fixed effects coefficients and statistics
abvg_model_output <- as.data.frame(summary(bm1)$coefficients)

# Write to CSV - model output
write.csv2(abvg_model_output, "C:/Users/Dpao/Desktop/Master's Biology/Master Thesis/#2 semester project/Ericoid Project/MGR_abvg_model_output_5thjfeb.csv", row.names=TRUE)

############# emmeans ###############
means = emmeans::emmeans(bm1, ~fungal_sp*plant_sp, type="response")
means
summary(means)

#Add p-value and statistical tests
means = update(means, infer = c(TRUE, TRUE), null = log(35), type="response",
               calc = c(n = ".wgt."))
summary(means)

meas_abvg = summary(means)
meas_abvg$response

# Write to CSV - emmeans output
write.csv2(meas_abvg, "C:/Users/Dpao/Desktop/Master's Biology/Master Thesis/#2 semester project/Ericoid Project/MGR_abvg_emmeans_output_5thjfeb.csv", row.names=TRUE)

## emmeans for plant only
means_abvg_p = emmeans::emmeans(bm1, ~plant_sp, type="response")

means_abvg_p
summary(means_abvg_p)

#Add p-value and statistical tests
means_abvg_p = update(means_abvg_p, infer = c(TRUE, TRUE), null = log(35), type="response",
                       calc = c(n = ".wgt."))
summary(means_abvg_p)

meas_abvg_p = summary(means_abvg_p)

# Write to CSV - emmeans plants across fungal treatment
write.csv2(meas_abvg_p, "C:/Users/Dpao/Desktop/Master's Biology/Master Thesis/#2 semester project/Ericoid Project/MGR_abvg_emmeans_plants_output_5thfeb.csv", row.names=TRUE)

## emmeans for fungi only
means_abvg_f = emmeans::emmeans(bm1, ~fungal_sp, type="response")

means_abvg_f
summary(means_abvg_f)

#Add p-value and statistical tests
means_abvg_f = update(means_abvg_f, infer = c(TRUE, TRUE), null = log(35), type="response",
                       calc = c(n = ".wgt."))
summary(means_abvg_f)

means_abvg_f = summary(means_abvg_f)

# Write to CSV
write.csv2(means_abvg_f, "C:/Users/Dpao/Desktop/Master's Biology/Master Thesis/#2 semester project/Ericoid Project/MGR_abvg_emmeans_fungi_output_5thfeb.csv", row.names=TRUE)


####################### Plot ABOVEGROUND ##########################

#add column to specify if error bars are going to be dashed or not
meas_abvg$lines=ifelse(sign(meas_abvg$lower.CL)==sign(meas_abvg$upper.CL),"1","2")


meas_abvg$plant_sp <- factor(meas_abvg$plant_sp,levels = c("G. shallon", "V. myrtillus", "P. japonica",
                                                     "V. vitis-idaea", "C. vulgaris", "V. angustifolium",
                                                     "K. latifolia", "R. arboreum", "R. ferrugineum"))


#Add Sterile Treatment to Dataset to have the empty space to make colonization and MGR plots match

# add "Sterile" to fungal species factor levels
meas_abvg$fungal_sp <- factor(meas_abvg$fungal_sp,
                               levels = c(levels(meas_abvg$fungal_sp), "Non-inoculated"))

# data frame with all plant species and the Sterile treatment
sterile_rows <- data.frame(
  plant_sp = levels(meas_abvg$plant_sp),
  fungal_sp = rep("Non-inoculated", length(levels(meas_abvg$plant_sp)))
)

# Add any other necessary columns with NA or appropriate default values
# Assuming other columns exist in meas_total, add them to sterile_rows with NA
other_cols <- setdiff(names(meas_abvg), c("plant_sp", "fungal_sp"))
for(col in other_cols) {
  sterile_rows[[col]] <- 0
}

# Combine the original data with the new sterile rows
meas_abvg_with_sterile <- rbind(meas_abvg, sterile_rows)

meas_abvg_with_sterile$fungal_sp <- factor(meas_abvg_with_sterile$fungal_sp,
                                            levels = c("Non-inoculated","H. gryndleri", "H. hepaticicola_PK 135-3", "H. hepaticicola_UAMH7357/ICMP18553",
                                                       "K. argillacea", "H. bicolor", "H. variabilis",
                                                       "O. maius", "Serendipitaceae_sp"))


png("figures/abvg_MGR/MGR_abvg_emmeans_updatedtryout.jpg", width = 18, height = 5.4, units ='in', res = 300)
ggplot(data=meas_abvg_with_sterile, mapping=aes(x=fungal_sp, color=fungal_sp)) +
  facet_grid(cols = vars(plant_sp)) +
  # Split geom_point into two layers - one for Sterile and one for other treatments
  geom_point(data = subset(meas_abvg_with_sterile, fungal_sp == "Non-inoculated"),
             aes(x= fungal_sp, y=response), 
             position=position_dodge(width =1), 
             size=1.5,
             alpha=0) +  
  geom_point(data = subset(meas_abvg_with_sterile, fungal_sp != "Non-inoculated"),
             aes(x= fungal_sp, y=response), 
             position=position_dodge(width =1), 
             size=1.5) +
  geom_errorbar(data = subset(meas_abvg_with_sterile, fungal_sp == "Non-inoculated"),
                aes(ymin=lower.CL, ymax=upper.CL), 
                position=position_dodge(width =1), 
                width=.75, 
                linewidth=0.7,
                alpha=0) +
  geom_errorbar(data = subset(meas_abvg_with_sterile, fungal_sp != "Non-inoculated"),
                aes(ymin=lower.CL, ymax=upper.CL, linetype=lines), 
                position=position_dodge(width =1), 
                width=.75, 
                linewidth=0.7) +
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
  ylab("Mycorrhizal Growth Response") +
  ylim (-1.5, 25) +
  guides(linetype = "none") +
  labs(color = "Fungal Isolates") +
  geom_hline(yintercept=0, linetype='dotted', col = "black")
dev.off()


################## pairwise comparison #######################
tab_results = pairs(means)
tab_results = as.data.frame(tab_results)


tab_results = separate(tab_results, col=contrast, into=c('comb1', 'comb2'), sep=' - ')

tab_results$p_value_stars <- makeStars(tab_results$p.value)

#remove "(" and ")" from some of the cells
tab_results$comb1 <- gsub("\\(", "", tab_results$comb1)
tab_results$comb2 <- gsub("\\(", "", tab_results$comb2)
tab_results$comb1 <- gsub("\\)", "", tab_results$comb1)
tab_results$comb2 <- gsub("\\)", "", tab_results$comb2)


#separate comb1 into fungi1 and plant1 and comb2 into fungi2 plant2
# function that will split the combination based on the plant name pattern
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


# keep only rows plant1=plant2
filtered_results_abvg <- tab_results[tab_results$plant1 == tab_results$plant2, ]

write.csv2(filtered_results_abvg, "C:/Users/Dpao/Desktop/Master's Biology/Master Thesis/#2 semester project/Ericoid Project/MGR_abvg_pairs_3rdmar.csv", row.names=FALSE)




################################################################################
############################# BELOWGROUND ######################################
################################################################################

#1. model
#2. emmeans with plant*fungi interaction and only plant
#3. plot MGR
#4. pairwise comparison

######################## model #########################
o = min(dat1$MGR_blwg)
p = o*-1+0.01

bm2 = lmer(log(MGR_blwg+p)~fungal_sp * plant_sp + Height.t0 + (1|block), data=dat1)
summary(bm2)
anova(bm2)
plot(bm2)


#save model output
# Extract fixed effects coefficients and statistics
blwg_model_output <- as.data.frame(summary(bm2)$coefficients)

# Write to CSV - model output
write.csv2(blwg_model_output, "C:/Users/Dpao/Desktop/Master's Biology/Master Thesis/#2 semester project/Ericoid Project/MGR_blwg_model_output_5thfeb.csv", row.names=TRUE)

############################ emmeans ####################################
means_blwg = emmeans::emmeans(bm2, ~fungal_sp*plant_sp, type="response")

means_blwg
summary(means_blwg)

#Add p-value and statistical tests
means_blwg = update(means_blwg, infer = c(TRUE, TRUE), null = log(35), type="response",
                 calc = c(n = ".wgt."))
summary(means_blwg)

meas_blwg = summary(means_blwg)
meas_blwg$response


# Write to CSV - emmeans output
write.csv2(meas_blwg, "C:/Users/Dpao/Desktop/Master's Biology/Master Thesis/#2 semester project/Ericoid Project/MGR_blwg_emmeans_output_5thfeb.csv", row.names=TRUE)

## emmeans for plant only
means_blwg_p = emmeans::emmeans(bm2, ~plant_sp, type="response")

means_blwg_p
summary(means_blwg_p)

#Add p-value and statistical tests
means_blwg_p = update(means_blwg_p, infer = c(TRUE, TRUE), null = log(35), type="response",
                      calc = c(n = ".wgt."))
summary(means_blwg_p)

meas_blwg_p = summary(means_blwg_p)

# Write to CSV - emmeans plants across fungal treatment
write.csv2(meas_blwg_p, "C:/Users/Dpao/Desktop/Master's Biology/Master Thesis/#2 semester project/Ericoid Project/MGR_blwg_emmeans_plants_output_5thfeb.csv", row.names=TRUE)



## emmeans for fungi only
means_blwg_f = emmeans::emmeans(bm2, ~fungal_sp, type="response")

means_blwg_f
summary(means_blwg_f)

#Add p-value and statistical tests
means_blwg_f = update(means_blwg_f, infer = c(TRUE, TRUE), null = log(35), type="response",
                      calc = c(n = ".wgt."))
summary(means_blwg_f)

means_blwg_f = summary(means_blwg_f)

# Write to CSV - emmeans fungi across plant species
write.csv2(means_blwg_f, "C:/Users/Dpao/Desktop/Master's Biology/Master Thesis/#2 semester project/Ericoid Project/MGR_blwg_emmeans_fungi_output_5thfeb.csv", row.names=TRUE)



################################ Plot BELOWGROUND #########################################
#add column to specify if error bars are going to be dashed or not
meas_blwg$lines=ifelse(sign(meas_blwg$lower.CL)==sign(meas_blwg$upper.CL),"1","2")


meas_blwg$plant_sp <- factor(meas_blwg$plant_sp,levels = c("G. shallon", "V. myrtillus", "P. japonica",
                                                           "V. vitis-idaea", "C. vulgaris", "V. angustifolium",
                                                           "K. latifolia", "R. arboreum", "R. ferrugineum"))


#Add Sterile Treatment to Dataset - same reason as aboveground

meas_blwg$fungal_sp <- factor(meas_blwg$fungal_sp,
                              levels = c(levels(meas_blwg$fungal_sp), "Non-inoculated"))

# Create a data frame with all plant species and the Sterile treatment
sterile_rows <- data.frame(
  plant_sp = levels(meas_blwg$plant_sp),
  fungal_sp = rep("Non-inoculated", length(levels(meas_blwg$plant_sp)))
)

# Add any other necessary columns with NA or appropriate default values
# Assuming other columns exist in meas_total, add them to sterile_rows with NA
other_cols <- setdiff(names(meas_blwg), c("plant_sp", "fungal_sp"))
for(col in other_cols) {
  sterile_rows[[col]] <- 0
}

# Combine the original data with the new sterile rows
meas_blwg_with_sterile <- rbind(meas_blwg, sterile_rows)

meas_blwg_with_sterile$fungal_sp <- factor(meas_blwg_with_sterile$fungal_sp,
                                           levels = c("Non-inoculated","H. gryndleri", "H. hepaticicola_PK 135-3", "H. hepaticicola_UAMH7357/ICMP18553",
                                                      "K. argillacea", "H. bicolor", "H. variabilis",
                                                      "O. maius", "Serendipitaceae_sp"))


png("figures/blwg_MGR/MGR_blwg_emmeans_updatedtryout.jpg", width = 18, height = 5.4, units ='in', res = 300)
ggplot(data=meas_blwg_with_sterile, mapping=aes(x=fungal_sp, color=fungal_sp)) +
  facet_grid(cols = vars(plant_sp)) +
  # Split geom_point into two layers - one for Sterile and one for other treatments
  geom_point(data = subset(meas_blwg_with_sterile, fungal_sp == "Non-inoculated"),
             aes(x= fungal_sp, y=response), 
             position=position_dodge(width =1), 
             size=1.5,
             alpha=0) +  
  geom_point(data = subset(meas_blwg_with_sterile, fungal_sp != "Non-inoculated"),
             aes(x= fungal_sp, y=response), 
             position=position_dodge(width =1), 
             size=1.5) +
  geom_errorbar(data = subset(meas_blwg_with_sterile, fungal_sp == "Non-inoculated"),
                aes(ymin=lower.CL, ymax=upper.CL), 
                position=position_dodge(width =1), 
                width=.75, 
                linewidth=0.7,
                alpha=0) +
  geom_errorbar(data = subset(meas_blwg_with_sterile, fungal_sp != "Non-inoculated"),
                aes(ymin=lower.CL, ymax=upper.CL, linetype=lines), 
                position=position_dodge(width =1), 
                width=.75, 
                linewidth=0.7) +
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
  ylab("Mycorrhizal Growth Response") +
  ylim (-1.5, 35.5) +
  guides(linetype = "none") +
  labs(color = "Fungal Isolates") +
  geom_hline(yintercept=0, linetype='dotted', col = "black")
dev.off()

######################### pairwise comparisons #############################3

tab_results_b = pairs(means_blwg)
tab_results_b = as.data.frame(tab_results_b)

tab_results_b = separate(tab_results_b, col=contrast, into=c('comb1', 'comb2'), sep=' - ')

tab_results_b$p_value_stars <- makeStars(tab_results_b$p.value)


#remove "(" and ")" from some of the cells
tab_results_b$comb1 <- gsub("\\(", "", tab_results_b$comb1)
tab_results_b$comb2 <- gsub("\\(", "", tab_results_b$comb2)
tab_results_b$comb1 <- gsub("\\)", "", tab_results_b$comb1)
tab_results_b$comb2 <- gsub("\\)", "", tab_results_b$comb2)


# Apply the function to both columns
tab_results_b = tab_results_b %>%
  mutate(
    comb1_split = map(comb1, split_combination),
    comb2_split = map(comb2, split_combination),
    
    fungi1 = map_chr(comb1_split, "fungi"),
    plant1 = map_chr(comb1_split, "plant"),
    fungi2 = map_chr(comb2_split, "fungi"),
    plant2 = map_chr(comb2_split, "plant")
  ) %>%
  select(-comb1_split, -comb2_split)  # Remove the intermediate columns

tab_results_b = tab_results_b %>%
  select(-comb1, -comb2)


# keep rows plant1=plant2
filtered_results_blwg <- tab_results_b[tab_results_b$plant1 == tab_results_b$plant2, ]

write.csv2(filtered_results_blwg, "C:/Users/Dpao/Desktop/Master's Biology/Master Thesis/#2 semester project/Ericoid Project/MGR_blwg_pairs_3rdmar.csv", row.names=FALSE)


############################################################################
########################### TOTAL BIOMASS ##################################
############################################################################

#1. model
#2. emmeans with plant*fungi interaction and only plant
#3. plot MGR
#4. pairwise comparison

########################## model ##############################
n = min(dat1$MGR)
h = n*-1+0.01

bm5 = lmer(log(MGR+h) ~fungal_sp * plant_sp + Height.t0 + (1|block), data=dat1)
summary(bm5)
anova(bm5)
plot(bm5)


#save model output
# Extract fixed effects coefficients and statistics
total_model_output <- as.data.frame(summary(bm5)$coefficients)

# Write to CSV
write.csv2(total_model_output, "C:/Users/Dpao/Desktop/Master's Biology/Master Thesis/#2 semester project/Ericoid Project/MGR_total_model_output_5thfeb.csv", row.names=TRUE)



################################## emmeans ####################################
means_total = emmeans::emmeans(bm5, ~fungal_sp*plant_sp, type="response")

means_total
summary(means_total)

#Add p-value and statistical tests
means_total = update(means_total, infer = c(TRUE, TRUE), null = log(35), type="response",
                     calc = c(n = ".wgt."))
summary(means_total)

meas_total = summary(means_total)
meas_total$response

# Write to CSV
write.csv2(meas_total, "C:/Users/Dpao/Desktop/Master's Biology/Master Thesis/#2 semester project/Ericoid Project/MGR_total_emmeans_output_5thfeb.csv", row.names=TRUE)

## emmeans for plant only
means_total_p = emmeans::emmeans(bm5, ~plant_sp, type="response")

means_total_p
summary(means_total_p)

#Add p-value and statistical tests
means_total_p = update(means_total_p, infer = c(TRUE, TRUE), null = log(35), type="response",
                     calc = c(n = ".wgt."))
summary(means_total_p)

meas_total_p = summary(means_total_p)

# Write to CSV
write.csv2(meas_total_p, "C:/Users/Dpao/Desktop/Master's Biology/Master Thesis/#2 semester project/Ericoid Project/MGR_total_emmeans_plants_output_5thfeb.csv", row.names=TRUE)

## emmeans for fungi only
means_total_f = emmeans::emmeans(bm5, ~fungal_sp, type="response")

means_total_f
summary(means_total_f)

#Add p-value and statistical tests
means_total_f = update(means_total_f, infer = c(TRUE, TRUE), null = log(35), type="response",
                       calc = c(n = ".wgt."))
summary(means_total_f)

meas_total_f = summary(means_total_f)

# Write to CSV
write.csv2(meas_total_f, "C:/Users/Dpao/Desktop/Master's Biology/Master Thesis/#2 semester project/Ericoid Project/MGR_total_emmeans_fungi_output_5thfeb.csv", row.names=TRUE)


################################ Plot TOTAL BIOMASS ####################################
#add column to specify if error bars are going to be dashed or not
meas_total$lines=ifelse(sign(meas_total$lower.CL)==sign(meas_total$upper.CL),"1","2")


meas_total$plant_sp <- factor(meas_total$plant_sp,levels = c("G. shallon", "V. myrtillus", "P. japonica",
                                                             "V. vitis-idaea", "C. vulgaris", "V. angustifolium",
                                                             "K. latifolia", "R. arboreum", "R. ferrugineum"))


#Add Sterile Treatment to Dataset - same as aboveground

# add "Sterile" to fungal species factor levels
meas_total$fungal_sp <- factor(meas_total$fungal_sp,
                               levels = c(levels(meas_total$fungal_sp), "Non-inoculated"))

# Create a data frame with all plant species and the Sterile treatment
sterile_rows <- data.frame(
  plant_sp = levels(meas_total$plant_sp),
  fungal_sp = rep("Non-inoculated", length(levels(meas_total$plant_sp)))
)

# Add any other necessary columns with NA or appropriate default values
# Assuming other columns exist in meas_total, add them to sterile_rows with NA
other_cols <- setdiff(names(meas_total), c("plant_sp", "fungal_sp"))
for(col in other_cols) {
  sterile_rows[[col]] <- 0
}

# Combine the original data with the new sterile rows
meas_total_with_sterile <- rbind(meas_total, sterile_rows)

meas_total_with_sterile$fungal_sp <- factor(meas_total_with_sterile$fungal_sp,
                               levels = c("Non-inoculated","H. gryndleri", "H. hepaticicola_PK 135-3", "H. hepaticicola_UAMH7357/ICMP18553",
                                          "K. argillacea", "H. bicolor", "H. variabilis",
                                          "O. maius", "Serendipitaceae_sp"))


png("figures/total_MGR/MGR_total_emmeans_updatedtryout.jpg", width = 18, height = 5.4, units ='in', res = 300)
ggplot(data=meas_total_with_sterile, mapping=aes(x=fungal_sp, color=fungal_sp)) +
  facet_grid(cols = vars(plant_sp)) +
  # Split geom_point into two layers - one for Sterile and one for other treatments
  geom_point(data = subset(meas_total_with_sterile, fungal_sp == "Non-inoculated"),
             aes(x= fungal_sp, y=response), 
             position=position_dodge(width =1), 
             size=1.5,
             alpha=0) +  
  geom_point(data = subset(meas_total_with_sterile, fungal_sp != "Non-inoculated"),
             aes(x= fungal_sp, y=response), 
             position=position_dodge(width =1), 
             size=1.5) +
  geom_errorbar(data = subset(meas_total_with_sterile, fungal_sp == "Non-inoculated"),
                aes(ymin=lower.CL, ymax=upper.CL), 
                position=position_dodge(width =1), 
                width=.75, 
                linewidth=0.7,
                alpha=0) +
  geom_errorbar(data = subset(meas_total_with_sterile, fungal_sp != "Non-inoculated"),
                aes(ymin=lower.CL, ymax=upper.CL, linetype=lines), 
                position=position_dodge(width =1), 
                width=.75, 
                linewidth=0.7) +
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
  ylab("Mycorrhizal Growth Response") +
  ylim (-1.5, 25) +
  guides(linetype = "none") +
  labs(color = "Fungal Isolates") +
  geom_hline(yintercept=0, linetype='dotted', col = "black")
dev.off()

################################ pairwise comparisons ############################

tab_results_t = pairs(means_total)
tab_results_t = as.data.frame(tab_results_t)

tab_results_t = separate(tab_results_t, col=contrast, into=c('comb1', 'comb2'), sep=' - ')
tab_results_t$p_value_stars <- makeStars(tab_results_t$p.value)

#remove "(" and ")" from some of the cells
tab_results_t$comb1 <- gsub("\\(", "", tab_results_t$comb1)
tab_results_t$comb2 <- gsub("\\(", "", tab_results_t$comb2)
tab_results_t$comb1 <- gsub("\\)", "", tab_results_t$comb1)
tab_results_t$comb2 <- gsub("\\)", "", tab_results_t$comb2)


# Apply the function to both columns
tab_results_t = tab_results_t %>%
  mutate(
    comb1_split = map(comb1, split_combination),
    comb2_split = map(comb2, split_combination),
    
    fungi1 = map_chr(comb1_split, "fungi"),
    plant1 = map_chr(comb1_split, "plant"),
    fungi2 = map_chr(comb2_split, "fungi"),
    plant2 = map_chr(comb2_split, "plant")
  ) %>%
  select(-comb1_split, -comb2_split)  # Remove the intermediate columns

tab_results_t = tab_results_t %>%
  select(-comb1, -comb2)

# Filter to keep only rows where plant species match
filtered_results <- tab_results_t[tab_results_t$plant1 == tab_results_t$plant2, ]

write.csv2(filtered_results, "C:/Users/Dpao/Desktop/Master's Biology/Master Thesis/#2 semester project/Ericoid Project/MGR_total_pairs_3rdmar.csv", row.names=FALSE)


############################################################################
########################### RATIO ABVG/BLWG ##################################

a = min(dat1$MGR_ratio)
b = a*-1+0.01

bm6 = lmer(log(MGR_ratio+b)~fungal_sp * plant_sp + Height.t0 + (1|block), data=dat1)
summary(bm6)
ratio_anova = anova(bm6) #interaction is not significant, so plot much simpler and not go to emmeans
plot(bm6)

#save model output
# Extract fixed effects coefficients and statistics
ratio_model_output <- as.data.frame(summary(bm6)$coefficients)

# Write to CSV
write.csv2(ratio_model_output, "C:/Users/Dpao/Desktop/Master's Biology/Master Thesis/#2 semester project/Ericoid Project/root_shoot_ratio_model_output_5thfeb.csv", row.names=TRUE)
write.csv2(ratio_anova, "C:/Users/Dpao/Desktop/Master's Biology/Master Thesis/#2 semester project/Ericoid Project/root_shoot_ratio_model_anova_5thfeb.csv", row.names=TRUE)



#for fungi
means_ratio_fungi = emmeans::emmeans(bm6, ~fungal_sp, type="response")
means_ratio_fungi
summary(means_ratio_fungi)
#Add p-value and statistical tests
means_ratio_fungi = update(means_ratio_fungi, infer = c(TRUE, TRUE), null = log(35), type="response",
                     calc = c(n = ".wgt."))
summary(means_ratio_fungi)
meas_ratio_fungi = summary(means_ratio_fungi)
meas_ratio_fungi$response

# Write to CSV for fungi
write.csv2(meas_ratio_fungi, "C:/Users/Dpao/Desktop/Master's Biology/Master Thesis/#2 semester project/Ericoid Project/root_shoot_ratio_emmeans_fungi_output_5thfeb.csv", row.names=TRUE)


#for plant
means_ratio_plant = emmeans::emmeans(bm6, ~plant_sp, type="response")
means_ratio_plant
summary(means_ratio_plant)

#Add p-value and statistical tests
means_ratio_plant = update(means_ratio_plant, infer = c(TRUE, TRUE), null = log(35), type="response",
                           calc = c(n = ".wgt."))
summary(means_ratio_plant)

meas_ratio_plant = summary(means_ratio_plant)
meas_ratio_plant$response

# Write to CSV for fungi
write.csv2(meas_ratio_plant, "C:/Users/Dpao/Desktop/Master's Biology/Master Thesis/#2 semester project/Ericoid Project/root_shoot_ratio_emmeans_plant_output_5thfeb.csv", row.names=TRUE)



####################### Plots ratio fungi ##########################

#add column to specify if error bars are going to be dashed or not
meas_ratio_fungi$lines=ifelse(sign(meas_ratio_fungi$lower.CL)==sign(meas_ratio_fungi$upper.CL),"1","2")


meas_ratio_fungi$fungal_sp <- factor(meas_ratio_fungi$fungal_sp,levels = c("Sterile","H. gryndleri", "H. hepaticicola_PK 135-3", "H. hepaticicola_UAMH7357/ICMP18553",
                                                                           "K. argillacea", "H. bicolor", "H. variabilis",
                                                                           "O. maius", "Serendipitaceae_sp"))

# code to make plots with ore than one plant species to compare it
ratio_fungi = ggplot(data=meas_ratio_fungi, mapping=aes(x= fungal_sp, color=fungal_sp)) +
  geom_point(aes(x= fungal_sp, y=response), position=position_dodge(width =1),size=1.5, alpha=0.3) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL, linetype=lines), position=position_dodge(width =1), width=0.4, linewidth=0.7) +
  theme_classic(base_size = 16) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent",
                                       colour = NA_character_),
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_color_manual(values = fungal_colors) +
  theme(panel.background = element_rect(fill = NA, color = "grey"))+
  theme(plot.title = element_text(hjust = 0.5))+
  ylab("Difference in shoot:root\nratio due to inoculation") +
  guides(linetype = "none") +
  labs(color = "Fungal Isolates") +
  geom_hline(yintercept=0, linetype='dotted', col = "black")

####################### Plots ratio plant ##########################
#add column to specify if error bars are going to be dashed or not
meas_ratio_plant$lines=ifelse(sign(meas_ratio_plant$lower.CL)==sign(meas_ratio_plant$upper.CL),"1","2")


meas_ratio_plant$plant_sp <- factor(meas_ratio_plant$plant_sp,levels = c("G. shallon", "V. myrtillus", "P. japonica",
                                                                         "V. vitis-idaea", "C. vulgaris", "V. angustifolium",
                                                                         "K. latifolia", "R. arboreum", "R. ferrugineum"))

# code to make plots with ore than one plant species to compare it
ratio_plants = ggplot(data=meas_ratio_plant, mapping=aes(x= plant_sp, color=plant_sp)) +
  geom_point(aes(x= plant_sp, y=response), position=position_dodge(width =1),size=1.5, alpha=0.3) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL, linetype=lines), position=position_dodge(width =1), width=.4, linewidth=0.7) +
  theme_classic(base_size = 16) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent",
                                       colour = NA_character_),
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(panel.background = element_rect(fill = NA, color = "grey"))+
  theme(plot.title = element_text(hjust = 0.5))+
  ylab("Difference in shoot:root\nratio due to inoculation") +
  guides(linetype = "none") +
  labs(color = "Plant species") +
  geom_hline(yintercept=0, linetype='dotted', col = "black")



plot_ratio = ggarrange(ratio_plants, ratio_fungi, nrow = 1)
png("figures/ratio_both_emmeans.jpg", width = 12, height = 4.5, units ='in', res = 300)
plot_ratio
dev.off()
