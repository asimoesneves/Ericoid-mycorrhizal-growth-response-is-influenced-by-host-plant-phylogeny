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

########################################################################
############################# ASSUMPTIONS ##############################
########################################################################
E1 <- resid(bm3, type = "pearson")
N <- nrow(dat1)
p <-length(coef(bm3))
sum(E1^2)/(N-p)
F1<-fitted(bm3)

plot(x = F1, y = E1, xlab = "Fitted values", ylab ="Pearson residuals")
abline(h = 0, lty = 2)

#homogeneity function
resid_plot_fit <- function(mod,dat) {
  E1 <- resid(mod, type = "pearson")
  N <- nrow(dat1)
  p <- length(coef(dat1))
  sum(E1^2)/(N-p)
  F1 <- fitted(mod)
  ggplot() +
    aes(F1, E1) +
    geom_point() +
    ggtitle("Residual plot") +
    xlab("Fitted Values") +
    ylab("Pearson residuals") 
}
resid_plot_fit(bm3,dat1[dat1[,"fungal_sp"]==MEB,])


plot(E1~dat1$Height.t0)

hist(resid(bm3))
qqnorm(resid(bm3))

plot(cooks.distance(bm3), type = "h", xlab ="Observation", ylab ="Cook distance")

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

# Function to extract plant species from a combination string
##### used after pairs() to extract only the interest combinations
extract_plant <- function(combo_string) {
  # Split string and take the last part (plant species)
  parts <- strsplit(combo_string, " ")[[1]]
  return(parts[2])  # Returns the plant species part
}


################################################################################
##using measurements.csv - height of sick-looking plants was deleted
######to give the correct number of samples in the plots
#read data


dat = read.csv2("Harvested_Plants_Processed_try1.csv")
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
#4. pairwaise comparison

m = min(dat1$MGR_abvg)
n = m*-1+0.01

bm1 = lmer(log(MGR_abvg+n)~fungal_sp * plant_sp + Height.t0 + (1|block), data=dat1)
summary(bm1)
anova(bm1)
plot(bm1)


#save model output
# Extract fixed effects coefficients and statistics
abvg_model_output <- as.data.frame(summary(bm1)$coefficients)

# Write to CSV - model output
write.csv2(abvg_model_output, "C:/Users/Dpao/Desktop/Master's Biology/Master Thesis/#2 semester project/Ericoid Project/MGR_abvg_model_output_5thjan.csv", row.names=TRUE)

##emmeans
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
write.csv2(meas_abvg, "C:/Users/Dpao/Desktop/Master's Biology/Master Thesis/#2 semester project/Ericoid Project/MGR_abvg_emmeans_output_5thjan.csv", row.names=TRUE)

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
write.csv2(meas_abvg_p, "C:/Users/Dpao/Desktop/Master's Biology/Master Thesis/#2 semester project/Ericoid Project/MGR_abvg_emmeans_plants_output_5thjan.csv", row.names=TRUE)

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
write.csv2(means_abvg_f, "C:/Users/Dpao/Desktop/Master's Biology/Master Thesis/#2 semester project/Ericoid Project/MGR_abvg_emmeans_fungi_output_6thjan.csv", row.names=TRUE)


####################### Plot ABOVEGROUND ##########################

#add column to specify if error bars are going to be dashed or not
meas_abvg$lines=ifelse(sign(meas_abvg$lower.CL)==sign(meas_abvg$upper.CL),"1","2")


meas_abvg$plant_sp <- factor(meas_abvg$plant_sp,levels = c("G_shallon", "V_myrtillus", "P_japonica",
                                                     "V_vitis-idaea", "C_vulgaris", "V_angustifolium",
                                                     "K_latifolia", "R_arboreum", "R_ferrugineum"))


#Add Sterile Treatment to Dataset

# First, add "Sterile" to fungal species factor levels
meas_abvg$fungal_sp <- factor(meas_abvg$fungal_sp,
                               levels = c(levels(meas_abvg$fungal_sp), "Sterile"))

# Create a data frame with all plant species and the Sterile treatment
sterile_rows <- data.frame(
  plant_sp = levels(meas_abvg$plant_sp),
  fungal_sp = rep("Sterile", length(levels(meas_abvg$plant_sp)))
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
                                            levels = c("Sterile","H_gryndleri", "H_hepaticicola_1", "H_hepaticicola_2",
                                                       "K_argillacea", "H_bicolor", "H_variabilis",
                                                       "O_maius", "Serendipitaceae_sp"))


png("figures/abvg_MGR/MGR_abvg_emmeans_updatedtryout.jpg", width = 18, height = 5.4, units ='in', res = 300)
ggplot(data=meas_abvg_with_sterile, mapping=aes(x=fungal_sp, color=fungal_sp)) +
  facet_grid(cols = vars(plant_sp)) +
  # Split geom_point into two layers - one for Sterile and one for other treatments
  geom_point(data = subset(meas_abvg_with_sterile, fungal_sp == "Sterile"),
             aes(x= fungal_sp, y=response), 
             position=position_dodge(width =1), 
             size=1.5,
             alpha=0) +  
  geom_point(data = subset(meas_abvg_with_sterile, fungal_sp != "Sterile"),
             aes(x= fungal_sp, y=response), 
             position=position_dodge(width =1), 
             size=1.5) +
  geom_errorbar(data = subset(meas_abvg_with_sterile, fungal_sp == "Sterile"),
                aes(ymin=lower.CL, ymax=upper.CL), 
                position=position_dodge(width =1), 
                width=.75, 
                linewidth=0.7,
                alpha=0) +
  geom_errorbar(data = subset(meas_abvg_with_sterile, fungal_sp != "Sterile"),
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
    strip.text = element_text(size = 15),
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title.y = element_text(size = 16),
    legend.title = element_text(size = 17),
    legend.text = element_text(size = 16)
  ) +
  ggtitle(label = "Aboveground Biomass") +
  ylab("Mycorrhizal Growth Response") +
  ylim (-1.5, 25) +
  guides(linetype = "none") +
  labs(color = "Fungal Isolates") +
  geom_hline(yintercept=0, linetype='dotted', col = "black")
dev.off()



################################################################################
############################# BELOWGROUND ######################################
################################################################################

#1. model
#2. emmeans with plant*fungi interaction and only plant
#3. plot MGR
#4. pairwaise comparison

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
write.csv2(blwg_model_output, "C:/Users/Dpao/Desktop/Master's Biology/Master Thesis/#2 semester project/Ericoid Project/MGR_blwg_model_output_5thjan.csv", row.names=TRUE)

##emmeans
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
write.csv2(meas_blwg, "C:/Users/Dpao/Desktop/Master's Biology/Master Thesis/#2 semester project/Ericoid Project/MGR_blwg_emmeans_output_5thjan.csv", row.names=TRUE)

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
write.csv2(meas_blwg_p, "C:/Users/Dpao/Desktop/Master's Biology/Master Thesis/#2 semester project/Ericoid Project/MGR_blwg_emmeans_plants_output_5thjan.csv", row.names=TRUE)



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
write.csv2(means_blwg_f, "C:/Users/Dpao/Desktop/Master's Biology/Master Thesis/#2 semester project/Ericoid Project/MGR_blwg_emmeans_fungi_output_6thjan.csv", row.names=TRUE)



############ BELOWGROUND plots
#add column to specify if error bars are going to be dashed or not
meas_blwg$lines=ifelse(sign(meas_blwg$lower.CL)==sign(meas_blwg$upper.CL),"1","2")


meas_blwg$plant_sp <- factor(meas_blwg$plant_sp,levels = c("G_shallon", "V_myrtillus", "P_japonica",
                                                     "V_vitis-idaea", "C_vulgaris", "V_angustifolium",
                                                     "K_latifolia", "R_arboreum", "R_ferrugineum"))


#Add Sterile Treatment to Dataset

# First, add "Sterile" to fungal species factor levels
meas_blwg$fungal_sp <- factor(meas_blwg$fungal_sp,
                              levels = c(levels(meas_blwg$fungal_sp), "Sterile"))

# Create a data frame with all plant species and the Sterile treatment
sterile_rows <- data.frame(
  plant_sp = levels(meas_blwg$plant_sp),
  fungal_sp = rep("Sterile", length(levels(meas_blwg$plant_sp)))
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
                                           levels = c("Sterile","H_gryndleri", "H_hepaticicola_1", "H_hepaticicola_2",
                                                      "K_argillacea", "H_bicolor", "H_variabilis",
                                                      "O_maius", "Serendipitaceae_sp"))


png("figures/blwg_MGR/MGR_blwg_emmeans_updatedtryout.jpg", width = 18, height = 5.4, units ='in', res = 300)
ggplot(data=meas_blwg_with_sterile, mapping=aes(x=fungal_sp, color=fungal_sp)) +
  facet_grid(cols = vars(plant_sp)) +
  # Split geom_point into two layers - one for Sterile and one for other treatments
  geom_point(data = subset(meas_blwg_with_sterile, fungal_sp == "Sterile"),
             aes(x= fungal_sp, y=response), 
             position=position_dodge(width =1), 
             size=1.5,
             alpha=0) +  
  geom_point(data = subset(meas_blwg_with_sterile, fungal_sp != "Sterile"),
             aes(x= fungal_sp, y=response), 
             position=position_dodge(width =1), 
             size=1.5) +
  geom_errorbar(data = subset(meas_blwg_with_sterile, fungal_sp == "Sterile"),
                aes(ymin=lower.CL, ymax=upper.CL), 
                position=position_dodge(width =1), 
                width=.75, 
                linewidth=0.7,
                alpha=0) +
  geom_errorbar(data = subset(meas_blwg_with_sterile, fungal_sp != "Sterile"),
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
    strip.text = element_text(size = 15),
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title.y = element_text(size = 16),
    legend.title = element_text(size = 17),
    legend.text = element_text(size = 16)
  ) +
  ggtitle(label = "Belowground Biomass") +
  ylab("Mycorrhizal Growth Response") +
  ylim (-1.5, 35.5) +
  guides(linetype = "none") +
  labs(color = "Fungal Isolates") +
  geom_hline(yintercept=0, linetype='dotted', col = "black")
dev.off()




############################################################################
########################### TOTAL BIOMASS ##################################
############################################################################

#1. model
#2. emmeans with plant*fungi interaction and only plant
#3. plot MGR
#4. pairwaise comparison

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
write.csv2(total_model_output, "C:/Users/Dpao/Desktop/Master's Biology/Master Thesis/#2 semester project/Ericoid Project/MGR_total_model_output_4thjan.csv", row.names=TRUE)



## emmeans
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
write.csv2(meas_total, "C:/Users/Dpao/Desktop/Master's Biology/Master Thesis/#2 semester project/Ericoid Project/MGR_total_emmeans_output_4thjan.csv", row.names=TRUE)

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
write.csv2(meas_total_p, "C:/Users/Dpao/Desktop/Master's Biology/Master Thesis/#2 semester project/Ericoid Project/MGR_total_emmeans_plants_output_4thjan.csv", row.names=TRUE)

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
write.csv2(meas_total_f, "C:/Users/Dpao/Desktop/Master's Biology/Master Thesis/#2 semester project/Ericoid Project/MGR_total_emmeans_fungi_output_6thjan.csv", row.names=TRUE)


############ TOTAL BIOMASS plots
#add column to specify if error bars are going to be dashed or not
meas_total$lines=ifelse(sign(meas_total$lower.CL)==sign(meas_total$upper.CL),"1","2")


meas_total$plant_sp <- factor(meas_total$plant_sp,levels = c("G_shallon", "V_myrtillus", "P_japonica",
                                                     "V_vitis-idaea", "C_vulgaris", "V_angustifolium",
                                                     "K_latifolia", "R_arboreum", "R_ferrugineum"))


#Add Sterile Treatment to Dataset

# First, add "Sterile" to fungal species factor levels
meas_total$fungal_sp <- factor(meas_total$fungal_sp,
                               levels = c(levels(meas_total$fungal_sp), "Sterile"))

# Create a data frame with all plant species and the Sterile treatment
sterile_rows <- data.frame(
  plant_sp = levels(meas_total$plant_sp),
  fungal_sp = rep("Sterile", length(levels(meas_total$plant_sp)))
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
                               levels = c("Sterile","H_gryndleri", "H_hepaticicola_1", "H_hepaticicola_2",
                                          "K_argillacea", "H_bicolor", "H_variabilis",
                                          "O_maius", "Serendipitaceae_sp"))


png("figures/total_MGR/MGR_total_emmeans_updatedtryout.jpg", width = 18, height = 5.4, units ='in', res = 300)
ggplot(data=meas_total_with_sterile, mapping=aes(x=fungal_sp, color=fungal_sp)) +
  facet_grid(cols = vars(plant_sp)) +
  # Split geom_point into two layers - one for Sterile and one for other treatments
  geom_point(data = subset(meas_total_with_sterile, fungal_sp == "Sterile"),
             aes(x= fungal_sp, y=response), 
             position=position_dodge(width =1), 
             size=1.5,
             alpha=0) +  
  geom_point(data = subset(meas_total_with_sterile, fungal_sp != "Sterile"),
             aes(x= fungal_sp, y=response), 
             position=position_dodge(width =1), 
             size=1.5) +
  geom_errorbar(data = subset(meas_total_with_sterile, fungal_sp == "Sterile"),
                aes(ymin=lower.CL, ymax=upper.CL), 
                position=position_dodge(width =1), 
                width=.75, 
                linewidth=0.7,
                alpha=0) +
  geom_errorbar(data = subset(meas_total_with_sterile, fungal_sp != "Sterile"),
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
    strip.text = element_text(size = 15),
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title.y = element_text(size = 16),
    legend.title = element_text(size = 17),
    legend.text = element_text(size = 16)
  ) +
  ggtitle(label = "Total Biomass") +
  ylab("Mycorrhizal Growth Response") +
  ylim (-1.5, 25) +
  guides(linetype = "none") +
  labs(color = "Fungal Isolates") +
  geom_hline(yintercept=0, linetype='dotted', col = "black")
dev.off()



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
write.csv2(ratio_model_output, "C:/Users/Dpao/Desktop/Master's Biology/Master Thesis/#2 semester project/Ericoid Project/root_shoot_ratio_model_output_5thjan.csv", row.names=TRUE)
write.csv2(ratio_anova, "C:/Users/Dpao/Desktop/Master's Biology/Master Thesis/#2 semester project/Ericoid Project/root_shoot_ratio_model_anova_6thjan.csv", row.names=TRUE)



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
write.csv2(meas_ratio_fungi, "C:/Users/Dpao/Desktop/Master's Biology/Master Thesis/#2 semester project/Ericoid Project/root_shoot_ratio_emmeans_fungi_output_5thjan.csv", row.names=TRUE)


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
write.csv2(meas_ratio_plant, "C:/Users/Dpao/Desktop/Master's Biology/Master Thesis/#2 semester project/Ericoid Project/root_shoot_ratio_emmeans_plant_output_5thjan.csv", row.names=TRUE)



####################### Plots ratio fungi ##########################

#add column to specify if error bars are going to be dashed or not
meas_ratio_fungi$lines=ifelse(sign(meas_ratio_fungi$lower.CL)==sign(meas_ratio_fungi$upper.CL),"1","2")


meas_ratio_fungi$fungal_sp <- factor(meas_ratio_fungi$fungal_sp,levels = c("H_gryndleri", "H_hepaticicola_1", "H_hepaticicola_2",
                                                       "K_argillacea", "H_bicolor", "H_variabilis",
                                                       "O_maius", "Serendipitaceae_sp"))

# code to make plots with ore than one plant species to compare it
png("figures/ratio_fungi_emmeans.jpg", width = 7, height = 5, units ='in', res = 300)
ggplot(data=meas_ratio_fungi, mapping=aes(x= fungal_sp, color=fungal_sp)) +
  geom_point(aes(x= fungal_sp, y=response), position=position_dodge(width =1),size=1.5, alpha=0.3) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL, linetype=lines), position=position_dodge(width =1), width=0.4, linewidth=0.7) +
  theme_classic(base_size = 20) +
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
dev.off()

####################### Plots ratio plant ##########################
#add column to specify if error bars are going to be dashed or not
meas_ratio_plant$lines=ifelse(sign(meas_ratio_plant$lower.CL)==sign(meas_ratio_plant$upper.CL),"1","2")


meas_ratio_plant$plant_sp <- factor(meas_ratio_plant$plant_sp,levels = c("G_shallon", "V_myrtillus", "P_japonica",
                                                             "V_vitis-idaea", "C_vulgaris", "V_angustifolium",
                                                             "K_latifolia", "R_arboreum", "R_ferrugineum"))

# code to make plots with ore than one plant species to compare it
png("figures/ratio_plant_emmeans.jpg", width = 7, height = 5, units ='in', res = 300)
ggplot(data=meas_ratio_plant, mapping=aes(x= plant_sp, color=plant_sp)) +
  geom_point(aes(x= plant_sp, y=response), position=position_dodge(width =1),size=1.5, alpha=0.3) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL, linetype=lines), position=position_dodge(width =1), width=.4, linewidth=0.7) +
  theme_classic(base_size = 20) +
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
dev.off()






##########################################################
###Pairs() to make pair contrasts:

### ABOVEGROUND
tab_results = pairs(means)
tab_results = as.data.frame(tab_results)


tab_results = separate(tab_results, col=contrast, into=c('comb1', 'comb2'), sep=' - ')

tab_results$p_value_stars <- makeStars(tab_results$p.value)


# Add plant species columns
tab_results$plant1 <- sapply(tab_results$comb1, extract_plant)
tab_results$plant2 <- sapply(tab_results$comb2, extract_plant)

# Filter to keep only rows where plant species match
filtered_results_abvg <- tab_results[tab_results$plant1 == tab_results$plant2, ]

# Optional: remove the helper columns if you don't need them
filtered_results_abvg$plant1 <- NULL
filtered_results_abvg$plant2 <- NULL

# Sort by the plant species (if you want to keep them grouped)
filtered_results_abvg <- filtered_results_abvg[order(sapply(filtered_results_abvg$comb1, extract_plant)), ]

write.csv2(filtered_results_abvg, "C:/Users/Dpao/Desktop/Master's Biology/Master Thesis/#2 semester project/Ericoid Project/MGR_abvg_pairs_5thjan.csv", row.names=FALSE)



### BELOWGROUND
tab_results_b = pairs(means_blwg)
tab_results_b = as.data.frame(tab_results_b)


tab_results_b = separate(tab_results_b, col=contrast, into=c('comb1', 'comb2'), sep=' - ')

tab_results_b$p_value_stars <- makeStars(tab_results_b$p.value)


# Add plant species columns
tab_results_b$plant1 <- sapply(tab_results_b$comb1, extract_plant)
tab_results_b$plant2 <- sapply(tab_results_b$comb2, extract_plant)

# Filter to keep only rows where plant species match
filtered_results_blwg <- tab_results[tab_results_b$plant1 == tab_results_b$plant2, ]

# Optional: remove the helper columns if you don't need them
filtered_results_blwg$plant1 <- NULL
filtered_results_blwg$plant2 <- NULL

# Sort by the plant species (if you want to keep them grouped)
filtered_results_blwg <- filtered_results_blwg[order(sapply(filtered_results_blwg$comb1, extract_plant)), ]

write.csv2(filtered_results_blwg, "C:/Users/Dpao/Desktop/Master's Biology/Master Thesis/#2 semester project/Ericoid Project/MGR_blwg_pairs_5thjan.csv", row.names=FALSE)

### TOTAL BIOMASS
tab_results_t = pairs(means_total)
tab_results_t = as.data.frame(tab_results_t)


tab_results_t = separate(tab_results_t, col=contrast, into=c('comb1', 'comb2'), sep=' - ')

tab_results_t$p_value_stars <- makeStars(tab_results_t$p.value)

#write.csv2(tab_results_t, "C:/Users/Dpao/Desktop/Master's Biology/Master Thesis/#2 semester project/Ericoid Project/MGR_total_contrasts.csv", row.names=FALSE)


# Function to extract plant species from a combination string
extract_plant <- function(combo_string) {
  # Split string and take the last part (plant species)
  parts <- strsplit(combo_string, " ")[[1]]
  return(parts[2])  # Returns the plant species part
}

# Add plant species columns
tab_results_t$plant1 <- sapply(tab_results_t$comb1, extract_plant)
tab_results_t$plant2 <- sapply(tab_results_t$comb2, extract_plant)

# Filter to keep only rows where plant species match
filtered_results <- tab_results_t[tab_results_t$plant1 == tab_results_t$plant2, ]

# Optional: remove the helper columns if you don't need them
filtered_results$plant1 <- NULL
filtered_results$plant2 <- NULL

# Sort by the plant species (if you want to keep them grouped)
filtered_results <- filtered_results[order(sapply(filtered_results$comb1, extract_plant)), ]

write.csv2(filtered_results, "C:/Users/Dpao/Desktop/Master's Biology/Master Thesis/#2 semester project/Ericoid Project/MGR_total_pairs_4thjan.csv", row.names=FALSE)



### RATIO
tab_results_r = pairs(means_r)
tab_results_r = as.data.frame(tab_results_r)


tab_results_r = separate(tab_results_r, col=contrast, into=c('comb1', 'comb2'), sep=' - ')

tab_results_r$p_value_stars <- makeStars(tab_results_r$p.value)

write.csv2(tab_results_r, "C:/Users/Dpao/Desktop/Master's Biology/Master Thesis/#2 semester project/Ericoid Project/MGR_ratio_contrasts.csv", row.names=FALSE)
