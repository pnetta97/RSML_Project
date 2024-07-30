# FERNERKUNDUNG PROJEKT

# Packages laden
suppressPackageStartupMessages({
  library(sp)
  library(terra)
  library(sf)
  library(lidR)
  library(future)
  library(dplyr)
  library(mapview)
  library(CAST)
  library(caret)
  library(lattice)
  library(ggplot2)
})


## Kronen Basteln

predlid <- readLAS("D:/prnetta/cloud9071b553a468cdae.las")
crop_liso <- clip_rectangle(predlid,349650, 5775400, 350650, 5776400)

#digital terrain model um anschschließen cutoff bei Höhe zu machen
dtm <- rasterize_terrain(crop_liso, 1, knnidw())
crop_liso <- crop_liso - dtm 

cropy_liso <- filter_poi(crop_liso, Z >= 13)

#Einzelne Bäume identifizieren
predcanopy <- rasterize_canopy(cropy_liso, 0.5, p2r(subcircle = 0.2), pkg = "terra")
predttopss <- locate_trees(predcanopy, lmf(ws = 5))
algo <- dalponte2016(predcanopy, predttopss)
predseg <- segment_trees(predlid, algo)
#Kronen (Polygone) der Bäume und 3D metrics errechnen
crown_met <- crown_metrics(predseg, geom = "convex", func =  ~list(z_max = max(Z), z_mean = mean(Z)))

# RGB Stack
LntnStack <-rast("D:/prnetta/LntnStack.tif")
LntnStack <- crop(LntnStack, c(349650, 5775400, 350650, 5776400))
LntnStackklein <- aggregate(LntnStack, fact = 2)

names(LntnStack) <- c("RGB_1","RGB_2","RGB_3","alpha") #RGB

#RGB mit Crown metrics und polygonen verbinden
ras_cm <- extract(LntnStackklein, crown_met, bind = T, fun = "mean")
#umwandeln und na's raus werfen
ras_cm <- st_as_sf(ras_cm)
ras_cm <- na.omit(ras_cm)

##Model trainieren
lntn_crowns <- st_read("D:/prnetta/TrainingCrowns.gpkg") #trainingsdaten aus dem Feld
lntn_crowns <- st_sample(lntn_crowns,1000) # mehr samples aus den polygonen ziehen

lntn_joined <- st_join(lntn_crowns, ras_cm, join = st_nearest_feature, left = T) # <- trainingsdaten Tabelle
#check
table(is.na(lntn_joined$RGB_1)) # mal sehen wieviele na's dabei sind
#NA raus schmeißen
#lntn_train <- na.omit(lntn_joined) #gibbet nit
lntn_train <- lntn_joined
predictors <- c("z_max","z_mean","RGB_1","RGB_2", "RGB_3", "alpha")
#extr <- extract(lntn_joined,lntn_crowns, bind=TRUE, fun = "mean")
#lntn_train <- st_as_sf(extr)



lntn_train <- st_drop_geometry(lntn_train)

model <- train(lntn_train[,predictors],
               lntn_train$art,
               # method = "ranger",
               #  importance = "permutation",
               #num.trees = 100,
               #importance=T,
               method="rf")
               # tuneLength=5,
               #ntree=50
               
model

#prediction
pred = predict(object = model, newdata = st_drop_geometry(ras_cm))
#species den polygonen zuweisen
ras_cm$species = pred

##plotten und visualisieren
plot(ras_cm[,"species"], col = cols)
mapview::mapview(ras_cm, zcol="species", col.regions = c("red","blue", "olivedrab3","sienna","pink",
                                                         "darkgreen","yellow2", "lightblue", "orange"))

st_write(ras_cm,"D:/prnetta/FinaleTabelle.shp")

##NDVI und VSDI in Tabelle den Polygonen zuorden
#Sentinel-2 daten handling
ixstack_10m <- rast(c("D:/prnetta/IMG_DATA/R10m/T32ULC_20240625T103631_B02_10m.jp2",
                      "D:/prnetta/IMG_DATA/R10m/T32ULC_20240625T103631_B04_10m.jp2",
                      "D:/prnetta/IMG_DATA/R10m/T32ULC_20240625T103631_B08_10m.jp2"))

ixstack_20m <- rast(c("D:/prnetta/IMG_DATA/R20m/T32ULC_20240625T103631_B12_20m.jp2",
                      "D:/prnetta/IMG_DATA/R20m/T32ULC_20240625T103631_B11_20m.jp2"))

ixstack_10re <- resample(ixstack_20m, ixstack_10m)
ixstack <- c(ixstack_10m, ixstack_10re)
ixstack <- project(ixstack,"EPSG:32632")

spatex <-c(349650, 350650, 5775400, 5776400)
ixstack <- crop(ixstack,spatex)

names(ixstack) <- substr(names(ixstack),
                         nchar(names(ixstack))-6, 
                         nchar(names(ixstack))-4) 

ixstack

#VSDI
ixstack$vsdi <- 1- (+(ixstack$B12 - ixstack$B02) + (ixstack$B04 - ixstack$B02))

#NDVI
ixstack$ndvi <- (ixstack$B08 - ixstack$B04)/(ixstack$B08 + ixstack$B04)

indizes <- c(ixstack$vsdi,ixstack$ndvi)

indizes <- st_as_sf(indizes)
#ab in die Tabelle damit
Luenten_Table <- extract(indizes, ras_cm, bind= T, fun = "mean")
Luenten_Table <- st_as_sf(Luenten_Table)
Luenten_Table

st_write(Luenten_Table,"D:/prnetta/Luenten_Tab_fin.shp") #<- Finale Tabelle


lntn_tab <- st_read("C:/Users/Paul/Documents/RData/Luenten_Tab_fin.shp")


#lntn_tab <- lntn_tab[,-c("z_max","z_mean","RGB_1","RGB_2", "RGB_3", "alpha")]

print(lntn_tab)

#Ergebnisse plotten
ggplot(lntn_tab, aes(x = vsdi, y = ndvi, color = species)) +
  scale_color_manual(values = 
                       c("red","blue", "olivedrab3","sienna","pink",
                                "darkgreen","yellow2", "lightblue", "orange"))+
  labs( x = "VSDI", y = "NDVI") +
  geom_point(size =0.5)+
  geom_smooth(method="lm", se=FALSE) +
  theme_light()

ggplot(lntn_tab, aes(x = vsdi, y = ndvi, color = species)) +
  facet_wrap(~ species) +
  scale_color_manual(values = c("red","blue", "olivedrab3","sienna","pink",
                                "darkgreen","yellow2", "lightblue", "orange"))+
  labs( x = "VSDI", y = "NDVI") +
  geom_point(size =0.3)+
  geom_smooth(method="lm", se=FALSE) +
  theme_light()

#Art individuen zählen
df_counts <- lntn_tab %>%
  group_by(lntn_tab$species) %>%
  summarise(count = n())

df_counts

#mapview species identification
mapview::mapview(lntn_tab, zcol="species", col.regions = c("red","blue", "olivedrab3","sienna","pink",
                                                         "darkgreen","yellow2", "lightblue", "orange"))
#<- subset(lntn_crowns,lntn_crowns$convhull_area >10)

#tmap species identification
library(tmap)
tm_shape(lntn_tab) +
  tm_polygons(
    col = "species",  # Die Variable für die Farbzuweisung
    palette = c("red", "blue", "olivedrab3", "sienna", "pink",
                "darkgreen", "yellow2", "lightblue", "orange"),
    title = "Species",  # Titel für die Legende
    size = 0.1  # Größe der Punkte
  ) +
  tm_legend(legend.position = c("left", "bottom")) +
  tm_layout(legend.outside = TRUE,  # Legende außerhalb des Plots
            legend.outside.position = "right")  # Position der Legende       


