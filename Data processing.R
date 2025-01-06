library(tidyverse)
library(ggplot2)
library(ggsci)
library(RColorBrewer)
library(lme4)
library(lmtest)
library(lmerTest)
library(sjPlot)

###### In this R script, the dry root biomass is calculated using the conversion ratios;
########### the total biomass and MGR are also calculated and everything is added to the
########### dataframe that is exported as a csv table ("Harvested_Plants_Processed.csv")
########### and is ready to be used for further modelling

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

###########################################################################################
########################################  RATIO  ##########################################
###########################################################################################

##calculate average wet/dry biomass ratio for each plant_sp

dat1 <- data.frame()

# add real values of wet and dry weight from experiment to the table
#need to do in 2 steps: first until row 243 - where largesubsample was done

ratio_real_values <- filter(dat, random_label<=243)
ratio_real_values = ratio_real_values %>% select(-c("unique_id", "fungal_sp", "replicate", "treatment",
                          "block", "Height.t0", "Dead.Sick", "random_label", "Height.t1",
                          "blockx", "total_root_wet", "small_subsample_wet", "total_root_dry",
                          "shoot_weight"))

colnames(ratio_real_values)[colnames(ratio_real_values) == c("plant_sp", "large_subsample_wet","large_subsample_dry") ] <- c("Plant.sp", "weight.wet.root", "weight.dry.root")
dat1 <- rbind(dat1, ratio_real_values)

# and then from 244 until the end with total_root_wet and total_root_dry
ratio_real_values1 <- filter(dat, random_label>243)
ratio_real_values1 = ratio_real_values1 %>% select(-c("unique_id", "fungal_sp", "replicate", "treatment",
                                                    "block", "Height.t0", "Dead.Sick", "random_label", "Height.t1",
                                                    "blockx", "large_subsample_wet", "small_subsample_wet", "large_subsample_dry",
                                                    "shoot_weight"))

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

#ratio_table$percentage_ratio = ratio_table$ratio*100
write.csv2(ratio_table, "ratio_table.csv", row.names = FALSE)

png("figures/biomass_ratio/Distribution_ratio_dots_updated1.jpg", width = 8, height = 3.5, units ='in', res = 300)
ggplot(dat1, aes(x=Plant.sp, y=ratio)) +
  geom_boxplot(aes(fill=Plant.sp)) +
  scale_fill_brewer(palette = "Paired") +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, binwidth = 0.011) +
  ggtitle(label = "Distribution of ratio dry/wet per plant sp") +
  theme_minimal(base_size = 15) +
  theme(
    panel.background = element_rect(fill = NA, color = "grey", size = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "transparent", colour = NA_character_),
    strip.background = element_rect(fill = "grey90", color = "grey", size = 0.5),
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
  labs(y= "wet weight:dry weight ratio") + 
  theme(axis.title.x=element_blank())
dev.off()


###########################################################################################
############################ root biomass calculations ####################################
###########################################################################################

#calculate small subsample wet weight
dat$small_subsample_wet = dat$total_root_wet-dat$large_subsample_wet

# calculate dry biomass for large subsample - those that the scale didn't detect
### use average ratio!


#join ratio in dataframe of main excel spreadsheet and
ratio_avrg = read.csv2("ratio_table.csv")
dat = dat %>% left_join(ratio_avrg)
#dat = dat %>% select(-c("percentage_ratio"))

i1 <- is.na(dat$large_subsample_dry)
dat$large_subsample_dry[i1] <- with(dat, large_subsample_wet[i1] *ratio[i1])


# calculate dry biomass for small subsample
## 1. if don't have large subsample -> use average ratio
## 2. if have large subsample -> use equation
dat = dat %>% drop_na(large_subsample_wet)
for (i in 1:630){
  if (dat$large_subsample_wet[i] == 0) {
    dat$small_subsample_dry[i] = dat$small_subsample_wet[i]*dat$ratio[i]
  }
  else {dat$small_subsample_dry[i] =
    (dat$large_subsample_dry[i]/dat$large_subsample_wet[i])*dat$small_subsample_wet[i]}
}


i2 <- is.na(dat$total_root_dry)
dat$total_root_dry[i2] <- with(dat, dat$small_subsample_dry[i2]+dat$large_subsample_dry[i2])



#################################################################################
################################  BIOMASS  #####################################
#################################################################################
#######################TOTAL BIOMASS, MGR CALCULATION############################
#################################################################################
#################################################################################

##MGR FOR ABOVEGROUND
STR=dat[dat[,"fungal_sp"]=="Sterile",c("plant_sp","shoot_weight")]

avr_STR = rep(0, times=9)
plant_nam = c("C_vulgaris", "G_shallon", "K_latifolia", "P_japonica", "R_arboreum",
              "R_ferrugineum", "V_angustifolium", "V_myrtillus", "V_vitis-idaea")

for (i in 1:9){
  x = STR[STR[,"plant_sp"]==plant_nam[i],]
  avr_STR[i] = mean(x[,"shoot_weight"], na.rm=TRUE)
}

abvg_biomass_STR_table = data.frame(plant_sp=plant_nam, abvg_biomass_STR=avr_STR)

###left join in dat
dat = dat%>%left_join(abvg_biomass_STR_table)

###calculate MGR = (Biomass EMF - Biomass STR)/Biomass STR in new column in the dataframe

dat$MGR_abvg = (dat$shoot_weight-dat$abvg_biomass_STR)/dat$abvg_biomass_STR


##MGR FOR BELOWGROUND
STR_b=dat[dat[,"fungal_sp"]=="Sterile",c("plant_sp","total_root_dry")]

avr_STR_b = rep(0, times=9)
plant_nam = c("C_vulgaris", "G_shallon", "K_latifolia",
              "P_japonica", "R_arboreum", "R_ferrugineum", "V_angustifolium", "V_myrtillus", "V_vitis-idaea")

for (i in 1:9){
  x = STR_b[STR_b[,"plant_sp"]==plant_nam[i],]
  avr_STR_b[i] = mean(x[,"total_root_dry"], na.rm=TRUE)
}

blwg_biomass_STR_table = data.frame(plant_sp=plant_nam, blwg_biomass_STR=avr_STR_b)

###left join in dat
dat = dat%>%left_join(blwg_biomass_STR_table)

###calculate MGR = (Biomass EMF - Biomass STR)/Biomass STR in new column in the dataframe

dat$MGR_blwg = (dat$total_root_dry-dat$blwg_biomass_STR)/dat$blwg_biomass_STR



#calculate total biomass = root + shoot (below and aboveground)

dat$total_biomass = dat$total_root_dry + dat$shoot_weight

#calculate average total biomass for sterile treatments
###subsample of only STR, avr of the ones that are the same plant sp and create new df

STR_total =dat[dat[,"fungal_sp"]=="Sterile",c("plant_sp","total_biomass")]
avr_STR_total = rep(0, times=9)
plant_nam = c("C_vulgaris", "G_shallon", "K_latifolia", "P_japonica", "R_arboreum",
              "R_ferrugineum", "V_angustifolium", "V_myrtillus", "V_vitis-idaea")

for (i in 1:9){
  x = STR_total[STR_total[,"plant_sp"]==plant_nam[i],]
  avr_STR_total[i] = mean(x[,"total_biomass"], na.rm=TRUE)
}


biomass_STR_table = data.frame(plant_sp=plant_nam, biomass_STR=avr_STR_total)

###left join in dat
dat = dat%>%left_join(biomass_STR_table)

###calculate MGR = (Biomass EMF - Biomass STR)/Biomass STR in new column in the dataframe

dat$MGR = (dat$total_biomass-dat$biomass_STR)/dat$biomass_STR

################################################################################
################# RATIO ABVG/BLWG - SHIFTS IN ALLOCATION #######################
################################################################################
dat$above_below_ratio = dat$shoot_weight/dat$total_root_dry

STR_ratio =dat[dat[,"fungal_sp"]=="Sterile",c("plant_sp","above_below_ratio")]
avr_STR_ratio = rep(0, times=9)
plant_nam = c("C_vulgaris", "G_shallon", "K_latifolia", "P_japonica", "R_arboreum",
              "R_ferrugineum", "V_angustifolium", "V_myrtillus", "V_vitis-idaea")

for (i in 1:9){
  x = STR_ratio[STR_total[,"plant_sp"]==plant_nam[i],]
  avr_STR_ratio[i] = mean(x[,"above_below_ratio"], na.rm=TRUE)
}

ratio_STR_table = data.frame(plant_sp=plant_nam, ratio_STR=avr_STR_ratio)
dat = dat%>%left_join(ratio_STR_table)

###calculate MGR = (Biomass EMF - Biomass STR)/Biomass STR in new column in the dataframe

dat$MGR_ratio = (dat$above_below_ratio-dat$ratio_STR)/dat$ratio_STR



png("figures/abvg_blwg_ratio.jpg", width = 30, height = 8, units ='in', res = 300)
ggplot(data=dat, mapping=aes(x=fungal_sp, fill=fungal_sp, y=above_below_ratio)) +
  facet_grid(cols = vars(plant_sp)) +
  scale_color_brewer(palette = "Paired") +
  geom_boxplot()+
  #geom_text(aes(label=round(total_percent_colonization, digits = 2)), position=position_dodge(width =1.2), vjust=1.6, color="black", size=3.5)+
  theme(panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        plot.background = element_rect(fill = "transparent",
                                       colour = NA_character_), # necessary to avoid drawing plot outline
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent")) +
  ggtitle(label = "Ratio abvg/blwg") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept=0, linetype='dotted', col = "black")
dev.off()


write.csv2(dat, "Harvested_Plants_Processed_try1.csv", row.names = FALSE)
