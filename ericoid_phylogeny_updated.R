library(ape)
packageVersion("ape")
library(tidyverse)
packageVersion("tidyverse")
library(abdiv)
packageVersion("abdiv")
library(doParallel)
packageVersion("doParallel")
library(foreach)
packageVersion("foreach")
library(pez)
packageVersion("pez")
library(phylosignal)
library(phytools)
library(phylobase)
library(reshape2)
library(xts)


#####################################
############# make dat ##############
#####################################

plants = c("Vaccinium_angustifolium", "Calluna_vulgaris",
           "Rhododendron_ferrugineum", "Rhododendron_arboreum",
           "Vaccinium_myrtillus", "Kalmia_latifolia", "Pieris_japonica",
           "Vaccinium_vitis-idaea", "Gaultheria_shallon")                                    
plant_sp = as.data.frame(plants)

#####################################
######## Build tree per plot ########
#####################################

load("GBOTB.extended.rda")                                                                                 
GBOTB.extended.congen <- congeneric.merge(GBOTB.extended,plants, split = "_", cite = FALSE)    
alltreesp <- c(GBOTB.extended.congen$tip.label) %>% as.data.frame()                                             
colnames(alltreesp) <- "plants"                                                                 

registerDoParallel(48)                                                                                         
filename = paste("trees/tree_", "plants_ericoid", sep = "")                                              
nodrop <- plant_sp%>% select(plants)                                                                 
todrop <- setdiff(alltreesp, nodrop) %>% as.vector()                                                   
todrop <- c(t(todrop))                                                                                   
p.tree <- drop.tip(GBOTB.extended.congen, unlist(todrop))                                                                                                                                                           
filepath = paste(filename, ".tre", sep = "")                                                               
write.tree(p.tree, filepath)                                                                                

##############################

p.tree = read.tree("trees/tree_plants_ericoid.tre")
plot.phylo(p.tree)

#############################################################
###################### ABOVEGROUND ##########################
#############################################################

dat_abvg = read.csv2("MGR_response_abvg.csv")
str(dat_abvg)

dat_abvg = dat_abvg%>%
  select(-c("df", "n","lower.CL", "upper.CL", "null", "t.ratio", "p.value", "lines"))
str(dat_abvg)
summary(dat_abvg)

dat_abvg <- transform(dat_abvg,plant_sp=gsub(pattern="C_vulgaris", replacement="Calluna_vulgaris", plant_sp))
dat_abvg <- transform(dat_abvg,plant_sp=gsub(pattern="G_shallon", replacement="Gaultheria_shallon", plant_sp))
dat_abvg <- transform(dat_abvg,plant_sp=gsub(pattern="K_latifolia", replacement="Kalmia_latifolia", plant_sp))
dat_abvg <- transform(dat_abvg,plant_sp=gsub(pattern="P_japonica", replacement="Pieris_japonica", plant_sp))
dat_abvg <- transform(dat_abvg,plant_sp=gsub(pattern="R_arboreum", replacement="Rhododendron_arboreum", plant_sp))
dat_abvg <- transform(dat_abvg,plant_sp=gsub(pattern="R_ferrugineum", replacement="Rhododendron_ferrugineum", plant_sp))
dat_abvg <- transform(dat_abvg,plant_sp=gsub(pattern="V_angustifolium", replacement="Vaccinium_angustifolium", plant_sp))
dat_abvg <- transform(dat_abvg,plant_sp=gsub(pattern="V_myrtillus", replacement="Vaccinium_myrtillus", plant_sp))
dat_abvg <- transform(dat_abvg,plant_sp=gsub(pattern="V_vitis-idaea", replacement="Vaccinium_vitis-idaea", plant_sp))

##### CREATE HEAT MAP OF MGR

dat_abvg1 = dat_abvg%>% select(-c("SE"))
response_abvg <- dcast(data = dat_abvg1, formula = plant_sp ~ fungal_sp)
rownames(response_abvg) <- response_abvg$plant_sp
response_abvg = response_abvg%>% select(-c("plant_sp"))

#min(response_abvg)
#max(response_abvg)

par(mar=c(0,0,0,0))

png("figures/heatmap_abvg.jpg", width = 10, height = 10, units ='in', res = 300)
phylo.heatmap(p.tree,response_abvg,
              split=c(0.8,0.7),fsize=c(1.5,0.5,1),
              standardize=FALSE,pts=FALSE, legend=TRUE, labels=TRUE)
dev.off()

par(mar=c(5.1,4.1,4.1,2.1)) ## reset margins to default


###### PHYLOGENETIC SIGNAL FOR EACH FUNGAL SPECIES

############## H. bicolor
dat_bicolor_abvg = dat_abvg %>% subset(dat_abvg$fungal_sp=="H_bicolor") %>%
  remove_rownames %>% column_to_rownames(var="plant_sp")

response_bicolor_abvg = setNames(dat_bicolor_abvg$response, rownames(dat_bicolor_abvg))
se_bicolor_abvg = setNames(dat_bicolor_abvg$SE, rownames(dat_bicolor_abvg))

# K-test
set.seed(123)
K_bicolor_abvg <-phylosig(p.tree, response_bicolor_abvg, se=se_bicolor_abvg,method="K", test=TRUE)
#print(K_bicolor_abvg); plot(K_bicolor_abvg)
# Pagel's Lambda
set.seed(123)
lambda_bicolor_abvg <-phylosig(p.tree, response_bicolor_abvg, se=se_bicolor_abvg,method="lambda", test=TRUE)
#print(lambda_bicolor_abvg); plot(lambda_bicolor_abvg)


############## H_gryndleri
dat_gryndleri_abvg = dat_abvg %>% subset(dat_abvg$fungal_sp=="H_gryndleri") %>%
  remove_rownames %>% column_to_rownames(var="plant_sp")

response_gryndleri_abvg = setNames(dat_gryndleri_abvg$response, rownames(dat_gryndleri_abvg))
se_gryndleri_abvg = setNames(dat_gryndleri_abvg$SE, rownames(dat_gryndleri_abvg))

# K-test
set.seed(123)
K_gryndleri_abvg <-phylosig(p.tree, response_gryndleri_abvg, se=se_gryndleri_abvg, method="K", test=TRUE)
#print(K_gryndleri_abvg); plot(K_gryndleri_abvg)
# Pagel's Lambda
set.seed(123)
lambda_gryndleri_abvg <-phylosig(p.tree, response_gryndleri_abvg, se=se_gryndleri_abvg, method="lambda", test=TRUE)
#print(lambda_gryndleri_abvg); plot(lambda_gryndleri_abvg)


############## H_hepaticicola_1
dat_hep1_abvg = dat_abvg %>% subset(dat_abvg$fungal_sp=="H_hepaticicola_1") %>%
  remove_rownames %>% column_to_rownames(var="plant_sp")

response_hep1_abvg = setNames(dat_hep1_abvg$response, rownames(dat_hep1_abvg))
se_hep1_abvg = setNames(dat_hep1_abvg$SE, rownames(dat_hep1_abvg))

# K-test
set.seed(123)
K_hep1_abvg <-phylosig(p.tree, response_hep1_abvg, se=se_hep1_abvg, method="K", test=TRUE)
#print(K_hep1_abvg); plot(K_hep1_abvg)
# Pagel's Lambda
set.seed(123)
lambda_hep1_abvg <-phylosig(p.tree, response_hep1_abvg, se=se_hep1_abvg, method="lambda", test=TRUE)
#print(lambda_hep1_abvg); plot(lambda_hep1_abvg)


############## H_hepaticicola_2
dat_hep2_abvg = dat_abvg %>% subset(dat_abvg$fungal_sp=="H_hepaticicola_2") %>%
  remove_rownames %>% column_to_rownames(var="plant_sp")

response_hep2_abvg = setNames(dat_hep2_abvg$response, rownames(dat_hep2_abvg))
se_hep2_abvg = setNames(dat_hep2_abvg$SE, rownames(dat_hep2_abvg))

# K-test
set.seed(123)
K_hep2_abvg <-phylosig(p.tree, response_hep2_abvg, se=se_hep2_abvg,method="K", test=TRUE)
#print(K_hep2_abvg); plot(K_hep2_abvg)
# Pagel's Lambda
set.seed(123)
lambda_hep2_abvg <-phylosig(p.tree, response_hep2_abvg, se=se_hep2_abvg, method="lambda", test=TRUE)
#print(lambda_hep2_abvg); plot(lambda_hep2_abvg)


############## H_variabilis
dat_variabilis_abvg = dat_abvg %>% subset(dat_abvg$fungal_sp=="H_variabilis") %>%
  remove_rownames %>% column_to_rownames(var="plant_sp")

response_variabilis_abvg = setNames(dat_variabilis_abvg$response, rownames(dat_variabilis_abvg))
se_variabilis_abvg = setNames(dat_variabilis_abvg$SE, rownames(dat_variabilis_abvg))

# K-test
set.seed(123)
K_variabilis_abvg <-phylosig(p.tree, response_variabilis_abvg, se=se_variabilis_abvg, method="K", test=TRUE)
#print(K_variabilis_abvg); plot(K_variabilis_abvg)
# Pagel's Lambda
set.seed(123)
lambda_variabilis_abvg <-phylosig(p.tree, response_variabilis_abvg, se=se_variabilis_abvg, method="lambda", test=TRUE)
#print(lambda_variabilis_abvg); plot(lambda_variabilis_abvg)


############## K_argillacea
dat_argillacea_abvg = dat_abvg %>% subset(dat_abvg$fungal_sp=="K_argillacea") %>%
  remove_rownames %>% column_to_rownames(var="plant_sp")

response_argillacea_abvg = setNames(dat_argillacea_abvg$response, rownames(dat_argillacea_abvg))
se_argillacea_abvg = setNames(dat_argillacea_abvg$SE, rownames(dat_argillacea_abvg))

# K-test
set.seed(123)
K_argillacea_abvg <-phylosig(p.tree, response_argillacea_abvg, se=se_argillacea_abvg, method="K", test=TRUE)
#print(K_argillacea_abvg); plot(K_argillacea_abvg)
# Pagel's Lambda
set.seed(123)
lambda_argillacea_abvg <-phylosig(p.tree, response_argillacea_abvg, se=se_argillacea_abvg,method="lambda", test=TRUE)
#print(lambda_argillacea_abvg); plot(lambda_argillacea_abvg)


############## O_maius
dat_maius_abvg = dat_abvg %>% subset(dat_abvg$fungal_sp=="O_maius") %>%
  remove_rownames %>% column_to_rownames(var="plant_sp")

response_maius_abvg = setNames(dat_maius_abvg$response, rownames(dat_maius_abvg))
se_maius_abvg = setNames(dat_maius_abvg$SE, rownames(dat_maius_abvg))

# K-test
set.seed(123)
K_maius_abvg <-phylosig(p.tree, response_maius_abvg, se=se_maius_abvg, method="K", test=TRUE)
#print(K_maius_abvg); plot(K_maius_abvg)
# Pagel's Lambda
set.seed(123)
lambda_maius_abvg <-phylosig(p.tree, response_maius_abvg, se=se_maius_abvg, method="lambda", test=TRUE)
#print(lambda_maius_abvg); plot(lambda_maius_abvg)


############## Serendipitaceae_sp
dat_Serendipitaceae_abvg = dat_abvg %>% subset(dat_abvg$fungal_sp=="Serendipitaceae_sp") %>%
  remove_rownames %>% column_to_rownames(var="plant_sp")

response_Serendipitaceae_abvg = setNames(dat_Serendipitaceae_abvg$response, rownames(dat_Serendipitaceae_abvg))
se_Serendipitaceae_abvg = setNames(dat_Serendipitaceae_abvg$SE, rownames(dat_Serendipitaceae_abvg))

# K-test
set.seed(123)
K_Serendipitaceae_abvg <-phylosig(p.tree, response_Serendipitaceae_abvg, se=se_Serendipitaceae_abvg, method="K", test=TRUE)
#print(K_Serendipitaceae_abvg); plot(K_Serendipitaceae_abvg)
# Pagel's Lambda
set.seed(123)
lambda_Serendipitaceae_abvg <-phylosig(p.tree, response_Serendipitaceae_abvg, se=se_Serendipitaceae_abvg, method="lambda", test=TRUE)
#print(lambda_Serendipitaceae_abvg); plot(lambda_Serendipitaceae_abvg)



### DATAFRAME WITH RESULTS FROM K AND LAMBDA TESTS FOR ALL FUNGI
column_names = c("K_test", "K_p_value", "Lambda_test", "Lambda_p_value")
H_bicolor_abvg = c(K_bicolor_abvg$K, K_bicolor_abvg$P, lambda_bicolor_abvg$lambda, lambda_bicolor_abvg$P)
H_gryndleri_abvg = c(K_gryndleri_abvg$K, K_gryndleri_abvg$P, lambda_gryndleri_abvg$lambda, lambda_gryndleri_abvg$P)
H_hepaticicola_1_abvg = c(K_hep1_abvg$K, K_hep1_abvg$P, lambda_hep1_abvg$lambda, lambda_hep1_abvg$P)
H_hepaticicola_2_abvg = c(K_hep2_abvg$K, K_hep2_abvg$P, lambda_hep2_abvg$lambda, lambda_hep2_abvg$P)
H_variabilis_abvg = c(K_variabilis_abvg$K, K_variabilis_abvg$P, lambda_variabilis_abvg$lambda, lambda_variabilis_abvg$P)
K_argillacea_abvg = c(K_argillacea_abvg$K, K_argillacea_abvg$P, lambda_argillacea_abvg$lambda, lambda_argillacea_abvg$P)
O_maius_abvg = c(K_maius_abvg$K, K_maius_abvg$P, lambda_maius_abvg$lambda, lambda_maius_abvg$P)
Serendipitaceae_sp_abvg = c(K_Serendipitaceae_abvg$K, K_Serendipitaceae_abvg$P, lambda_Serendipitaceae_abvg$lambda, lambda_Serendipitaceae_abvg$P)

df_abvg <- as.data.frame(rbind(H_bicolor_abvg, H_gryndleri_abvg, H_hepaticicola_1_abvg,
                          H_hepaticicola_2_abvg, H_variabilis_abvg, K_argillacea_abvg,
                          O_maius_abvg, Serendipitaceae_sp_abvg))
colnames(df_abvg)=column_names
df_abvg
write.csv2(df_abvg, "C:/Users/Dpao/Desktop/Master's Biology/Master Thesis/#2 semester project/Ericoid Project/phylogenetic_analysis_abvg.csv", row.names=TRUE)



#############################################################
###################### BELOWGROUND ##########################
#############################################################

dat_blwg = read.csv2("MGR_response_blwg.csv")
str(dat_blwg)

dat_blwg = dat_blwg%>%
  select(-c("df", "n","lower.CL", "upper.CL", "null", "t.ratio", "p.value", "lines"))
str(dat_blwg)
summary(dat_blwg)

dat_blwg <- transform(dat_blwg,plant_sp=gsub(pattern="C_vulgaris", replacement="Calluna_vulgaris", plant_sp))
dat_blwg <- transform(dat_blwg,plant_sp=gsub(pattern="G_shallon", replacement="Gaultheria_shallon", plant_sp))
dat_blwg <- transform(dat_blwg,plant_sp=gsub(pattern="K_latifolia", replacement="Kalmia_latifolia", plant_sp))
dat_blwg <- transform(dat_blwg,plant_sp=gsub(pattern="P_japonica", replacement="Pieris_japonica", plant_sp))
dat_blwg <- transform(dat_blwg,plant_sp=gsub(pattern="R_arboreum", replacement="Rhododendron_arboreum", plant_sp))
dat_blwg <- transform(dat_blwg,plant_sp=gsub(pattern="R_ferrugineum", replacement="Rhododendron_ferrugineum", plant_sp))
dat_blwg <- transform(dat_blwg,plant_sp=gsub(pattern="V_angustifolium", replacement="Vaccinium_angustifolium", plant_sp))
dat_blwg <- transform(dat_blwg,plant_sp=gsub(pattern="V_myrtillus", replacement="Vaccinium_myrtillus", plant_sp))
dat_blwg <- transform(dat_blwg,plant_sp=gsub(pattern="V_vitis-idaea", replacement="Vaccinium_vitis-idaea", plant_sp))

##### CREATE HEAT MAP OF MGR

dat_blwg1 = dat_blwg%>% select(-c("SE"))
response_blwg <- dcast(data = dat_blwg1, formula = plant_sp ~ fungal_sp)
rownames(response_blwg) <- response_blwg$plant_sp
response_blwg = response_blwg%>% select(-c("plant_sp"))

#min(response_blwg)
#max(response_blwg)

par(mar=c(0,0,0,0))

png("figures/heatmap_blwg.jpg", width = 10, height = 10, units ='in', res = 300)
phylo.heatmap(p.tree,response_blwg,
              split=c(0.8,0.7),fsize=c(1.5,0.5,1),
              standardize=FALSE,pts=FALSE, legend=TRUE, labels=TRUE)
dev.off()

par(mar=c(5.1,4.1,4.1,2.1)) ## reset margins to default


###### PHYLOGENETIC SIGNAL FOR EACH FUNGAL SPECIES

############## H. bicolor
dat_bicolor_blwg = dat_blwg %>% subset(dat_blwg$fungal_sp=="H_bicolor") %>%
  remove_rownames %>% column_to_rownames(var="plant_sp")

response_bicolor_blwg = setNames(dat_bicolor_blwg$response, rownames(dat_bicolor_blwg))
se_bicolor_blwg = setNames(dat_bicolor_blwg$SE, rownames(dat_bicolor_blwg))

# K-test
set.seed(123)
K_bicolor_blwg <-phylosig(p.tree, response_bicolor_blwg, se=se_bicolor_blwg, method="K", test=TRUE)
#print(K_bicolor_blwg); plot(K_bicolor_blwg)
# Pagel's Lambda
set.seed(123)
lambda_bicolor_blwg <-phylosig(p.tree, response_bicolor_blwg, se=se_bicolor_blwg, method="lambda", test=TRUE)
#print(lambda_bicolor_blwg); plot(lambda_bicolor_blwg)


############## H_gryndleri
dat_gryndleri_blwg = dat_blwg %>% subset(dat_blwg$fungal_sp=="H_gryndleri") %>%
  remove_rownames %>% column_to_rownames(var="plant_sp")

response_gryndleri_blwg = setNames(dat_gryndleri_blwg$response, rownames(dat_gryndleri_blwg))
se_gryndleri_blwg = setNames(dat_gryndleri_blwg$SE, rownames(dat_gryndleri_blwg))

# K-test
set.seed(123)
K_gryndleri_blwg <-phylosig(p.tree, response_gryndleri_blwg, se=se_gryndleri_blwg, method="K", test=TRUE)
#print(K_gryndleri_blwg); plot(K_gryndleri_blwg)
# Pagel's Lambda
set.seed(123)
lambda_gryndleri_blwg <-phylosig(p.tree, response_gryndleri_blwg, se=se_gryndleri_blwg, method="lambda", test=TRUE)
#print(lambda_gryndleri_blwg); plot(lambda_gryndleri_blwg)


############## H_hepaticicola_1
dat_hep1_blwg = dat_blwg %>% subset(dat_blwg$fungal_sp=="H_hepaticicola_1") %>%
  remove_rownames %>% column_to_rownames(var="plant_sp")

response_hep1_blwg = setNames(dat_hep1_blwg$response, rownames(dat_hep1_blwg))
se_hep1_blwg = setNames(dat_hep1_blwg$SE, rownames(dat_hep1_blwg))

# K-test
set.seed(123)
K_hep1_blwg <-phylosig(p.tree, response_hep1_blwg, se=se_hep1_blwg, method="K", test=TRUE)
#print(K_hep1_blwg); plot(K_hep1_blwg)
# Pagel's Lambda
set.seed(123)
lambda_hep1_blwg <-phylosig(p.tree, response_hep1_blwg, se=se_hep1_blwg, method="lambda", test=TRUE)
#print(lambda_hep1_blwg); plot(lambda_hep1_blwg)


############## H_hepaticicola_2
dat_hep2_blwg = dat_blwg %>% subset(dat_blwg$fungal_sp=="H_hepaticicola_2") %>%
  remove_rownames %>% column_to_rownames(var="plant_sp")

response_hep2_blwg = setNames(dat_hep2_blwg$response, rownames(dat_hep2_blwg))
se_hep2_blwg = setNames(dat_hep2_blwg$SE, rownames(dat_hep2_blwg))

# K-test
set.seed(123)
K_hep2_blwg <-phylosig(p.tree, response_hep2_blwg, se=se_hep2_blwg, method="K", test=TRUE)
#print(K_hep2_blwg); plot(K_hep2_blwg)
# Pagel's Lambda
set.seed(123)
lambda_hep2_blwg <-phylosig(p.tree, response_hep2_blwg, se=se_hep2_blwg, method="lambda", test=TRUE)
#print(lambda_hep2_blwg); plot(lambda_hep2_blwg)


############## H_variabilis
dat_variabilis_blwg = dat_blwg %>% subset(dat_blwg$fungal_sp=="H_variabilis") %>%
  remove_rownames %>% column_to_rownames(var="plant_sp")

response_variabilis_blwg = setNames(dat_variabilis_blwg$response, rownames(dat_variabilis_blwg))
se_variabilis_blwg = setNames(dat_variabilis_blwg$SE, rownames(dat_variabilis_blwg))

# K-test
set.seed(123)
K_variabilis_blwg <-phylosig(p.tree, response_variabilis_blwg, se=se_variabilis_blwg, method="K", test=TRUE)
#print(K_variabilis_blwg); plot(K_variabilis_blwg)
# Pagel's Lambda
set.seed(123)
lambda_variabilis_blwg <-phylosig(p.tree, response_variabilis_blwg, se=se_variabilis_blwg, method="lambda", test=TRUE)
#print(lambda_variabilis_blwg); plot(lambda_variabilis_blwg)


############## K_argillacea
dat_argillacea_blwg = dat_blwg %>% subset(dat_blwg$fungal_sp=="K_argillacea") %>%
  remove_rownames %>% column_to_rownames(var="plant_sp")

response_argillacea_blwg = setNames(dat_argillacea_blwg$response, rownames(dat_argillacea_blwg))
se_argillacea_blwg = setNames(dat_argillacea_blwg$SE, rownames(dat_argillacea_blwg))

# K-test
set.seed(123)
K_argillacea_blwg <-phylosig(p.tree, response_argillacea_blwg, se=se_argillacea_blwg, method="K", test=TRUE)
#print(K_argillacea_blwg); plot(K_argillacea_blwg)
# Pagel's Lambda
set.seed(123)
lambda_argillacea_blwg <-phylosig(p.tree, response_argillacea_blwg, se=se_argillacea_blwg, method="lambda", test=TRUE)
#print(lambda_argillacea_blwg); plot(lambda_argillacea_blwg)


############## O_maius
dat_maius_blwg = dat_blwg %>% subset(dat_blwg$fungal_sp=="O_maius") %>%
  remove_rownames %>% column_to_rownames(var="plant_sp")

response_maius_blwg = setNames(dat_maius_blwg$response, rownames(dat_maius_blwg))
se_maius_blwg = setNames(dat_maius_blwg$SE, rownames(dat_maius_blwg))

# K-test
set.seed(123)
K_maius_blwg <-phylosig(p.tree, response_maius_blwg, se=se_maius_blwg, method="K", test=TRUE)
#print(K_maius_blwg); plot(K_maius_blwg)
# Pagel's Lambda
set.seed(123)
lambda_maius_blwg <-phylosig(p.tree, response_maius_blwg, se=se_maius_blwg, method="lambda", test=TRUE)
#print(lambda_maius_blwg); plot(lambda_maius_blwg)


############## Serendipitaceae_sp
dat_Serendipitaceae_blwg = dat_blwg %>% subset(dat_blwg$fungal_sp=="Serendipitaceae_sp") %>%
  remove_rownames %>% column_to_rownames(var="plant_sp")

response_Serendipitaceae_blwg = setNames(dat_Serendipitaceae_blwg$response, rownames(dat_Serendipitaceae_blwg))
se_Serendipitaceae_blwg = setNames(dat_Serendipitaceae_blwg$SE, rownames(dat_Serendipitaceae_blwg))

# K-test
set.seed(123)
K_Serendipitaceae_blwg <-phylosig(p.tree, response_Serendipitaceae_blwg, se=se_Serendipitaceae_blwg, method="K", test=TRUE)
#print(K_Serendipitaceae_blwg); plot(K_Serendipitaceae_blwg)
# Pagel's Lambda
set.seed(123)
lambda_Serendipitaceae_blwg <-phylosig(p.tree, response_Serendipitaceae_blwg, se=se_Serendipitaceae_blwg, method="lambda", test=TRUE)
#print(lambda_Serendipitaceae_blwg); plot(lambda_Serendipitaceae_blwg)



### DATAFRAME WITH RESULTS FROM K AND LAMBDA TESTS FOR ALL FUNGI
column_names = c("K_test", "K_p_value", "Lambda_test", "Lambda_p_value")
H_bicolor_blwg = c(K_bicolor_blwg$K, K_bicolor_blwg$P, lambda_bicolor_blwg$lambda, lambda_bicolor_blwg$P)
H_gryndleri_blwg = c(K_gryndleri_blwg$K, K_gryndleri_blwg$P, lambda_gryndleri_blwg$lambda, lambda_gryndleri_blwg$P)
H_hepaticicola_1_blwg = c(K_hep1_blwg$K, K_hep1_blwg$P, lambda_hep1_blwg$lambda, lambda_hep1_blwg$P)
H_hepaticicola_2_blwg = c(K_hep2_blwg$K, K_hep2_blwg$P, lambda_hep2_blwg$lambda, lambda_hep2_blwg$P)
H_variabilis_blwg = c(K_variabilis_blwg$K, K_variabilis_blwg$P, lambda_variabilis_blwg$lambda, lambda_variabilis_blwg$P)
K_argillacea_blwg = c(K_argillacea_blwg$K, K_argillacea_blwg$P, lambda_argillacea_blwg$lambda, lambda_argillacea_blwg$P)
O_maius_blwg = c(K_maius_blwg$K, K_maius_blwg$P, lambda_maius_blwg$lambda, lambda_maius_blwg$P)
Serendipitaceae_sp_blwg = c(K_Serendipitaceae_blwg$K, K_Serendipitaceae_blwg$P, lambda_Serendipitaceae_blwg$lambda, lambda_Serendipitaceae_blwg$P)

df_blwg <- as.data.frame(rbind(H_bicolor_blwg, H_gryndleri_blwg, H_hepaticicola_1_blwg,
                               H_hepaticicola_2_blwg, H_variabilis_blwg, K_argillacea_blwg,
                               O_maius_blwg, Serendipitaceae_sp_blwg))
colnames(df_blwg)=column_names
df_blwg
write.csv2(df_blwg, "C:/Users/Dpao/Desktop/Master's Biology/Master Thesis/#2 semester project/Ericoid Project/phylogenetic_analysis_blwg.csv", row.names=TRUE)



#############################################################
##################### TOTAL BIOMASS #########################
#############################################################

dat_t = read.csv2("MGR_response_total.csv")
str(dat_t)

dat_t = dat_t%>%
  select(-c("df", "n","lower.CL", "upper.CL", "null", "t.ratio", "p.value", "lines"))
str(dat_t)
summary(dat_t)

dat_t <- transform(dat_t,plant_sp=gsub(pattern="C_vulgaris", replacement="Calluna_vulgaris", plant_sp))
dat_t <- transform(dat_t,plant_sp=gsub(pattern="G_shallon", replacement="Gaultheria_shallon", plant_sp))
dat_t <- transform(dat_t,plant_sp=gsub(pattern="K_latifolia", replacement="Kalmia_latifolia", plant_sp))
dat_t <- transform(dat_t,plant_sp=gsub(pattern="P_japonica", replacement="Pieris_japonica", plant_sp))
dat_t <- transform(dat_t,plant_sp=gsub(pattern="R_arboreum", replacement="Rhododendron_arboreum", plant_sp))
dat_t <- transform(dat_t,plant_sp=gsub(pattern="R_ferrugineum", replacement="Rhododendron_ferrugineum", plant_sp))
dat_t <- transform(dat_t,plant_sp=gsub(pattern="V_angustifolium", replacement="Vaccinium_angustifolium", plant_sp))
dat_t <- transform(dat_t,plant_sp=gsub(pattern="V_myrtillus", replacement="Vaccinium_myrtillus", plant_sp))
dat_t <- transform(dat_t,plant_sp=gsub(pattern="V_vitis-idaea", replacement="Vaccinium_vitis-idaea", plant_sp))

##### CREATE HEAT MAP OF MGR

dat_t1 = dat_t%>% select(-c("SE"))
response_t <- dcast(data = dat_t1, formula = plant_sp ~ fungal_sp)
rownames(response_t) <- response_t$plant_sp
response_t = response_t%>% select(-c("plant_sp"))

#min(response_t)
#max(response_t)

par(mar=c(0,0,0,0))

png("figures/heatmap_total.jpg", width = 10, height = 10, units ='in', res = 300)
phylo.heatmap(p.tree,response_t,
              split=c(0.8,0.7),fsize=c(1.5,0.5,1),
              standardize=FALSE,pts=FALSE, legend=TRUE, labels=TRUE)
dev.off()

par(mar=c(5.1,4.1,4.1,2.1)) ## reset margins to default


###### PHYLOGENETIC SIGNAL FOR EACH FUNGAL SPECIES

############## H. bicolor
dat_bicolor_t = dat_t %>% subset(dat_t$fungal_sp=="H_bicolor") %>%
  remove_rownames %>% column_to_rownames(var="plant_sp")

response_bicolor_t = setNames(dat_bicolor_t$response, rownames(dat_bicolor_t))
se_bicolor_t = setNames(dat_bicolor_t$SE, rownames(dat_bicolor_t))

# K-test
set.seed(123)
K_bicolor_t <-phylosig(p.tree, response_bicolor_t, se=se_bicolor_t, method="K", test=TRUE)
#print(K_bicolor_t); plot(K_bicolor_t)
# Pagel's Lambda
set.seed(123)
lambda_bicolor_t <-phylosig(p.tree, response_bicolor_t, se=se_bicolor_t, method="lambda", test=TRUE)
#print(lambda_bicolor_t); plot(lambda_bicolor_t)


############## H_gryndleri
dat_gryndleri_t = dat_t %>% subset(dat_t$fungal_sp=="H_gryndleri") %>%
  remove_rownames %>% column_to_rownames(var="plant_sp")

response_gryndleri_t = setNames(dat_gryndleri_t$response, rownames(dat_gryndleri_t))
se_gryndleri_t = setNames(dat_gryndleri_t$SE, rownames(dat_gryndleri_t))

# K-test
set.seed(123)
K_gryndleri_t <-phylosig(p.tree, response_gryndleri_t, se=se_gryndleri_t, method="K", test=TRUE)
#print(K_gryndleri_t); plot(K_gryndleri_t)
# Pagel's Lambda
set.seed(123)
lambda_gryndleri_t <-phylosig(p.tree, response_gryndleri_t, se=se_gryndleri_t, method="lambda", test=TRUE)
#print(lambda_gryndleri_t); plot(lambda_gryndleri_t)


############## H_hepaticicola_1
dat_hep1_t = dat_t %>% subset(dat_t$fungal_sp=="H_hepaticicola_1") %>%
  remove_rownames %>% column_to_rownames(var="plant_sp")

response_hep1_t = setNames(dat_hep1_t$response, rownames(dat_hep1_t))
se_hep1_t = setNames(dat_hep1_t$SE, rownames(dat_hep1_t))

# K-test
set.seed(123)
K_hep1_t <-phylosig(p.tree, response_hep1_t, se=se_hep1_t, method="K", test=TRUE)
#print(K_hep1_t); plot(K_hep1_t)
# Pagel's Lambda
set.seed(123)
lambda_hep1_t <-phylosig(p.tree, response_hep1_t, se=se_hep1_t, method="lambda", test=TRUE)
#print(lambda_hep1_t); plot(lambda_hep1_t)


############## H_hepaticicola_2
dat_hep2_t = dat_t %>% subset(dat_t$fungal_sp=="H_hepaticicola_2") %>%
  remove_rownames %>% column_to_rownames(var="plant_sp")

response_hep2_t = setNames(dat_hep2_t$response, rownames(dat_hep2_t))
se_hep2_t = setNames(dat_hep2_t$SE, rownames(dat_hep2_t))

# K-test
set.seed(123)
K_hep2_t <-phylosig(p.tree, response_hep2_t, se=se_hep2_t, method="K", test=TRUE)
#print(K_hep2_t); plot(K_hep2_t)
# Pagel's Lambda
set.seed(123)
lambda_hep2_t <-phylosig(p.tree, response_hep2_t, se=se_hep2_t, method="lambda", test=TRUE)
#print(lambda_hep2_t); plot(lambda_hep2_t)


############## H_variabilis
dat_variabilis_t = dat_t %>% subset(dat_t$fungal_sp=="H_variabilis") %>%
  remove_rownames %>% column_to_rownames(var="plant_sp")

response_variabilis_t = setNames(dat_variabilis_t$response, rownames(dat_variabilis_t))
se_variabilis_t = setNames(dat_variabilis_t$SE, rownames(dat_variabilis_t))

# K-test
set.seed(123)
K_variabilis_t <-phylosig(p.tree, response_variabilis_t, se=se_variabilis_t, method="K", test=TRUE)
#print(K_variabilis_t); plot(K_variabilis_t)
# Pagel's Lambda
set.seed(123)
lambda_variabilis_t <-phylosig(p.tree, response_variabilis_t, se=se_variabilis_t, method="lambda", test=TRUE)
#print(lambda_variabilis_t); plot(lambda_variabilis_t)


############## K_argillacea
dat_argillacea_t = dat_t %>% subset(dat_t$fungal_sp=="K_argillacea") %>%
  remove_rownames %>% column_to_rownames(var="plant_sp")

response_argillacea_t = setNames(dat_argillacea_t$response, rownames(dat_argillacea_t))
se_argillacea_t = setNames(dat_argillacea_t$SE, rownames(dat_argillacea_t))

# K-test
set.seed(123)
K_argillacea_t <-phylosig(p.tree, response_argillacea_t, se=se_argillacea_t, method="K", test=TRUE)
#print(K_argillacea_t); plot(K_argillacea_t)
# Pagel's Lambda
set.seed(123)
lambda_argillacea_t <-phylosig(p.tree, response_argillacea_t, se=se_argillacea_t, method="lambda", test=TRUE)
#print(lambda_argillacea_t); plot(lambda_argillacea_t)


############## O_maius
dat_maius_t = dat_t %>% subset(dat_t$fungal_sp=="O_maius") %>%
  remove_rownames %>% column_to_rownames(var="plant_sp")

response_maius_t = setNames(dat_maius_t$response, rownames(dat_maius_t))
se_maius_t = setNames(dat_maius_t$SE, rownames(dat_maius_t))

# K-test
set.seed(123)
K_maius_t <-phylosig(p.tree, response_maius_t, se=se_maius_t, method="K", test=TRUE)
#print(K_maius_t); plot(K_maius_t)
# Pagel's Lambda
set.seed(123)
lambda_maius_t <-phylosig(p.tree, response_maius_t, se=se_maius_t, method="lambda", test=TRUE)
#print(lambda_maius_t); plot(lambda_maius_t)


############## Serendipitaceae_sp
dat_Serendipitaceae_t = dat_t %>% subset(dat_t$fungal_sp=="Serendipitaceae_sp") %>%
  remove_rownames %>% column_to_rownames(var="plant_sp")

response_Serendipitaceae_t = setNames(dat_Serendipitaceae_t$response, rownames(dat_Serendipitaceae_t))
se_Serendipitaceae_t = setNames(dat_Serendipitaceae_t$SE, rownames(dat_Serendipitaceae_t))

# K-test
set.seed(123)
K_Serendipitaceae_t <-phylosig(p.tree, response_Serendipitaceae_t, se=se_Serendipitaceae_t, method="K", test=TRUE)
#print(K_Serendipitaceae_t); plot(K_Serendipitaceae_t)
# Pagel's Lambda
set.seed(123)
lambda_Serendipitaceae_t <-phylosig(p.tree, response_Serendipitaceae_t, se=se_Serendipitaceae_t, method="lambda", test=TRUE)
#print(lambda_Serendipitaceae_t); plot(lambda_Serendipitaceae_t)



### DATAFRAME WITH RESULTS FROM K AND LAMBDA TESTS FOR ALL FUNGI
column_names = c("K_test", "K_p_value", "Lambda_test", "Lambda_p_value")
H_bicolor_t = c(K_bicolor_t$K, K_bicolor_t$P, lambda_bicolor_t$lambda, lambda_bicolor_t$P)
H_gryndleri_t = c(K_gryndleri_t$K, K_gryndleri_t$P, lambda_gryndleri_t$lambda, lambda_gryndleri_t$P)
H_hepaticicola_1_t = c(K_hep1_t$K, K_hep1_t$P, lambda_hep1_t$lambda, lambda_hep1_t$P)
H_hepaticicola_2_t = c(K_hep2_t$K, K_hep2_t$P, lambda_hep2_t$lambda, lambda_hep2_t$P)
H_variabilis_t = c(K_variabilis_t$K, K_variabilis_t$P, lambda_variabilis_t$lambda, lambda_variabilis_t$P)
K_argillacea_t = c(K_argillacea_t$K, K_argillacea_t$P, lambda_argillacea_t$lambda, lambda_argillacea_t$P)
O_maius_t = c(K_maius_t$K, K_maius_t$P, lambda_maius_t$lambda, lambda_maius_t$P)
Serendipitaceae_sp_t = c(K_Serendipitaceae_t$K, K_Serendipitaceae_t$P, lambda_Serendipitaceae_t$lambda, lambda_Serendipitaceae_t$P)

df_t <- as.data.frame(rbind(H_bicolor_t, H_gryndleri_t, H_hepaticicola_1_t,
                               H_hepaticicola_2_t, H_variabilis_t, K_argillacea_t,
                               O_maius_t, Serendipitaceae_sp_t))
colnames(df_t)=column_names
df_t
write.csv2(df_t, "C:/Users/Dpao/Desktop/Master's Biology/Master Thesis/#2 semester project/Ericoid Project/phylogenetic_analysis_total.csv", row.names=TRUE)
