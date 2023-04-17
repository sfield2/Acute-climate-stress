############## ACUTE CLIMATE STRESS ################


#### ENVIRONMENTS & DATA #######################################################
packages <-c('tidyverse','sf','ggplot2','rgdal','rgeos','magrittr','plyr','FedData',
             'raster','dplyr','leastcostpath','gdistance', 'ggpubr','zoo', 'reshape2')
for(p in packages) if(p %in% rownames(installed.packages()) == F) { install.packages(p) }
for(p in packages) suppressPackageStartupMessages(library(p,quietly=T,character.only=T))

setwd("...")
theme_set(theme_bw())

######### STEP 1.0: IMPORT FUNCTIONS & DATA #####################################
### functions from Reese (2020)
base::source('./FUNCTIONS/Reese2020paleocar.R')

### CRS
longlat.projection <- sp::CRS('+proj=longlat +datum=WGS84 +ellps=WGS84')

### data
dem <- raster('./DATA/Climate/VEPIIN_NED_1.tif',na.rm=T)
fv <- readOGR(dsn="./DATA/Climate", layer="Farview")
mvl <- readOGR(dsn="./DATA/Climate", layer="Mesa_verde_landform")
fv_rm <- readOGR(dsn="./DATA/settlement",layer="fv_rm")
fv_kiva <- readOGR(dsn="./DATA/settlement",layer="fv_kivas")


### re-project data to ensure accurate
dem <- raster::projectRaster(dem,crs=longlat.projection)
fv <- spTransform(fv,crs(longlat.projection))
mvl <- spTransform(mvl,crs(longlat.projection))
fv_rm <- spTransform(fv_rm,crs(longlat.projection))
fv_kiva <- spTransform(fv_kiva,crs(longlat.projection))

### check data
plot(dem)
points(fv,col="black")
lines(mvl,col="blue")


### ceramic data
# calibration
cal_full <- read.csv("C:/users/seanp/Documents/R/DATING_COMPARISON/DATA/VEP_calibration_groupcounts.csv",header=T,fileEncoding = 'UTF-8-BOM')
cal <- cal_full[,4:18]

# settlement data
fv_ceramics <- read.csv("./DATA/settlement/Farview_ceramics.csv",header=T,fileEncoding = 'UTF-8-BOM')

### Create Brainerd-Robinson function (Peeples 2011)
BR <- function(x) {
  rd <- dim(x)[1]
  results <- matrix(0,rd,rd)
  for (s1 in 1:rd) {
    for (s2 in 1:rd) {
      x1Temp <- as.numeric(x[s1, ])
      x2Temp <- as.numeric(x[s2, ])
      br.temp <- 0
      results[s1,s2] <- 200 - (sum(abs(x1Temp - x2Temp)))}}
  row.names(results) <- row.names(x)
  colnames(results) <- row.names(x)
  return(results)}

### Set appropriate time ranges for occupation periods
ranges<- c("0700-0799","0800-0879","0880-0919","0920-0979","0980-1019", "1020-1059",
           "1060-1099","1100-1179","1180-1249","1250-1280")
period <- 1:10
ranges <- data.frame(ranges,period)

######### STEP 2.0: DEFINE STUDAY AREAS ########################################
#### define areas within 2 hr round-trip travel from center of village

# set limits for computational savings and create accumulated cost surfaces
neigh <- 16
dem.aggregate <- 5
dem.aggregated <- aggregate(dem, fact=dem.aggregate, fun=mean)
tobler_cs <- leastcostpath::create_slope_cs(dem = dem.aggregated, cost_function = "tobler", neighbours = neigh)
acc_fv_cs<- accCost(tobler_cs,fv)

# limit accumulated cost surfaces
limit_fv<- acc_fv_cs<3600%>%
limit_fv_l<-rasterToContour(limit_fv)
limit_fv_l <- limit_fv_l[1,]
limit_fv_l_sf <- st_as_sf(limit_fv_l)

# check the limits
plot(limit_fv_l_sf,col="white",add=T)

# export catchment areas example
fv_catchment <- as(limit_fv_l_sf, 'Spatial')
writeOGR(fv_catchment,"./output","Farview_catchment",driver="ESRI Shapefile")

# clean environment
# rm(list=ls())

######### STEP 3.0: CLIMATE AND MAIZE NICHE RECON ##############################
# incorporate catchment
fv_catchment <- readOGR(dsn="./DATA/Climate", layer="Farview_catchment")
fv_catchment <- spTransform(fv_catchment,crs(longlat.projection))

######### STEP 3.1: LOCAL CLIMATE CONDITIONS ###################################
growing_niche_fv <- paleocar(template = extent(fv_catchment),
                          label = 'paleoproductivity',
                          raw.dir = paste0('./output/FARVIEW/spatial-products/paleocar/raw/'),
                          extraction.dir = paste0('./output/FARVIEW/spatial-products/paleocar/extractions/'),
                          prcp_threshold = 350,
                          gdd_threshold = 1800,
                          years = 1:2000,
                          force.redo = F)

### extract reconstruction for study and take average
niche_fv <- growing_niche_fv$niche%>%
  as.data.frame()%>%
  set_colnames(1:2000)%>%
  subset(select=c(600:1300))%>%
  t%>%
  as.data.frame()%>%
  mutate(sum = rowSums(across(where(is.numeric))))%>%
  mutate(pct_in_niche = (sum/110)*100)%>%
  mutate(year = 600:1300)

precip_fv <- growing_niche_fv$PPT%>%
  as.data.frame()%>%
  set_colnames(1:2000)%>%
  select(c(600:1300))%>%
  t%>%
  as.data.frame()%>%
  mutate(average_precip = rowMeans(across(where(is.numeric))))%>%
  mutate(average_precip_cm = average_precip/10)%>%
  mutate(year = 600:1300)

gdd_fv <- growing_niche_fv$GDD%>%
  as.data.frame()%>%
  set_colnames(1:2000)%>%
  select(c(600:1300))%>%
  t%>%
  as.data.frame()%>%
  mutate(average_gdd = rowMeans(across(where(is.numeric))))%>%
  mutate(year = 600:1300)


### create data frame and export
fv_recon <- as.data.frame(cbind(niche_fv$year,niche_fv$pct_in_niche,precip_fv$average_precip_cm,gdd_fv$average_gdd))
colnames(fv_recon) <- c("Year", "Percent_in_niche", "Average_precipitation_in_cm", "Average_growing_degree_days")

# write.csv(fv_recon,"./output/farview_reconstruction.csv", row.names = FALSE)

######### STEP 3.2: REGIONAL CLIMATE CONDITIONS ################################
growing_niche_mvl <- paleocar(template = extent(mvl),
                             label = 'paleoproductivity',
                             raw.dir = paste0('./output/MESAVERDELANDFORM/spatial-products/paleocar/raw/'),
                             extraction.dir = paste0('./output/MESAVERDELANDFORM/spatial-products/paleocar/extractions/'),
                             prcp_threshold = 350,
                             gdd_threshold = 1800,
                             years = 1:2000,
                             force.redo = F)

niche_mvl <- growing_niche_mvl$niche%>%
  as.data.frame()%>%
  set_colnames(1:2000)%>%
  select(c(600:1300))%>%
  t%>%
  as.data.frame()%>%
  mutate(sum = rowSums(across(where(is.numeric))))%>%
  mutate(pct_in_niche = (sum/980)*100)%>%
  mutate(year = 600:1300)

precip_mvl <- growing_niche_mvl$PPT%>%
  as.data.frame()%>%
  set_colnames(1:2000)%>%
  select(c(600:1300))%>%
  t%>%
  as.data.frame()%>%
  mutate(average_precip = rowMeans(across(where(is.numeric))))%>%
  mutate(average_precip_cm = average_precip/10)%>%
  mutate(year = 600:1300)

gdd_mvl <- growing_niche_mvl$GDD%>%
  as.data.frame()%>%
  set_colnames(1:2000)%>%
  select(c(600:1300))%>%
  t%>%
  as.data.frame()%>%
  mutate(average_gdd = rowMeans(across(where(is.numeric))))%>%
  mutate(year = 600:1300)

### create data frame and export
mvl_recon <- as.data.frame(cbind(niche_mvl$year,niche_mvl$pct_in_niche,precip_mvl$average_precip_cm,gdd_mvl$average_gdd))
colnames(mvl_recon) <- c("Year", "Percent_in_niche", "Average_precipitation_in_cm", "Average_growing_degree_days")

# write.csv(mvl_recon,"./output/mesaverdelandform_reconstruction.csv", row.names = FALSE)

######### STEP 4.0: ANALYZE CLIMATE RECONSTRUCTIONS ############################
# read in select outputs 
fv_recon <- read.csv("./DATA/Climate/farview_reconstruction.csv",header=T)
mvl_recon <- read.csv("./DATA/Climate/mesaverdelandform_reconstruction.csv",header=T)

######### STEP 4.1: IDENTIFY ACCUTE STRESS PERIODS #############################
### calculate mean conditions for each local area
## MIGHT DELETE???
fv_av_precip <- mean(fv_recon$Average_precipitation_in_cm)
fv_av_temp <- mean(fv_recon$Average_growing_degree_days)

### Run following calculations for each village
### 4.1.1 calculate generational averages in both conditions
### 4.1.2 calculate percent difference between generational av and annual conditions in each year
### 4.1.3 transform into quartiles
### 4.1.4 incorporate percentage of region in niche

fv_variab<- fv_recon %>%
  # 4.1.1
  mutate(Gen_Precipitation = rollmean(Average_precipitation_in_cm,20,fill=NA,align="center"))%>%
  mutate(Gen_Temperature = rollmean(Average_growing_degree_days,20,fill=NA,align="center"))%>%
  # 4.1.2
  mutate(Precipitation_diff = abs(1 - (Average_precipitation_in_cm/Gen_Precipitation)))%>%
  mutate(Temperature_diff = abs(1 - (Average_growing_degree_days/Gen_Temperature)))%>%
  # 4.1.3
  mutate(Precipitation_diff_quartile = ntile(Precipitation_diff,4))%>%
  mutate(Temperature_diff_quartile = ntile(Temperature_diff,4))%>%
  # 4.1.4 
  mutate(Region_low_niche = mvl_recon$Percent_in_niche)

# identify years with any of the three measures
fv_stress <- list(fv_variab$Year[fv_variab$Precipitation_diff_quartile == "4"],
             fv_variab$Year[fv_variab$Temperature_diff_quartile == "4"],
             fv_variab$Year[fv_variab$Percent_in_niche < 50],
             fv_variab$Year[fv_variab$Region_low_niche < 50 & fv_variab$Percent_in_niche <50])

# count how many stress measures in each year 
fv_stress <- as.data.frame(unlist(fv_stress))
stress_ct <- count(fv_stress$`unlist(fv_stress)`)

# build complete stress regime
fv_stress_regime <- fv_variab$Year%>%
  data.frame(Year=.)%>%
  mutate(stress = stress_ct$freq[match(Year,stress_ct$x)])%>%
  mutate_all(funs(ifelse(is.na(.), 0, .)))%>%
  mutate(gen_av = rollmean(stress, 9, align="right", fill=NA))

# determine periods with 5 or more consecutive years where more than one measure is met
potential_stress_periods <- subset(fv_stress_regime, gen_av > "1")

sp <- NULL 
sapply(seq(nrow(potential_stress_periods)), function(x)
{
  ifelse((sum(diff(potential_stress_periods[x:(x+4), "Year"], 1)) == 4 &
            sum(diff(potential_stress_periods[x:(x+4), "Year"], 1) == 1) == 4),
         sp <<- rbind(sp, potential_stress_periods[x:(x+4),]),"")
})

fv_stress_periods <- unique(sp)

######### STEP 4.2: PLOT ACUTE STRESS REGIME ###############################
# build metric for plotting
fv_stress_regime$gen_av_quant <- ntile(fv_stress_regime$gen_av,5)
fv_stress_regime$gen_av_quant[is.na(fv_stress_regime$gen_av_quant)]<-1

#build plot
p1<-ggplot()+
  geom_col(data=fv_stress_regime,aes(Year,1,fill=factor(gen_av_quant)),width=2)+
  scale_fill_manual(values=c("#fee5d9","#fcae91","#fb6a4a","#de2d26","#a50f15"),
                    name="Stress",labels=c("Low Stress","","","","Acute Stress"))+
  scale_x_continuous(limits=c(595,1305),breaks=seq(600,1300,by=50))+
  theme(legend.position = 'top',
        axis.ticks.length.x=unit(.15, "cm"),
        legend.text = element_text(size=18,color="black"),
        axis.ticks.y = element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.title=element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        title=element_text(size=25,color="black"))+ 
  guides(fill = guide_legend(title.position = "top", 
                             title.hjust = 0.5,
                             label.position = "bottom")) 




p2<-ggplot()+
  geom_col(data = fv_stress_periods, aes(x = Year, y = 1), fill = "coral1",width=1,alpha=0.5)+
  scale_x_continuous(limits=c(595,1305),breaks=seq(600,1300,by=50))+
  annotate("text", x = 600, y = .75, col="grey30",size=6,alpha=0.5, fontface =2,hjust=0,
           label = "Acute Stress
Period")+
  annotate("segment", x = 700, xend = 800, y = .8, yend = .8,size = 1,
           colour = "grey30", alpha=0.5)+
  
  
  annotate("segment", x = 700, xend = 800, y = .2, yend = .2,size = 1,
           colour = "grey30", alpha=0.5)+
  annotate("text", x = 600, y = .2, col="grey30",size=6,alpha=0.5, fontface =2,hjust=0,
           label = "No. of Years")+
  
  
  annotate("text", x = 820, y = .15, col="coral3",size=6,alpha=0.85, fontface =2,hjust=0,
           label = "8")+
  annotate("text", x = 985, y = .15, col="coral3",size=6,alpha=0.85, fontface =2,hjust=0,
           label = "7")+
  annotate("text", x = 1017, y = .15, col="coral3",size=6,alpha=0.85, fontface =2,hjust=0,
           label = "6")+
  annotate("text", x = 1057, y = .15, col="coral3",size=6,alpha=0.85, fontface =2,hjust=0,
           label = "6")+
  annotate("text", x = 1104, y = .15, col="coral3",size=6,alpha=0.85, fontface =2,hjust=0,
           label = "9")+
  annotate("text", x = 1147, y = .15, col="coral3",size=6,alpha=0.85, fontface =2,hjust=0,
           label = "5")+
  annotate("text", x = 1208, y = .15, col="coral3",size=6,alpha=0.85, fontface =2,hjust=0,
           label = "9")+
  annotate("text", x = 1275, y = .15, col="coral3",size=6,alpha=0.85, fontface =2,hjust=0,
           label = "9")+
  
  
  theme(legend.position="none",
        axis.ticks.length.x=unit(.15, "cm"),
        axis.ticks.y = element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.title=element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())


plot1<-ggarrange(p1,p2,
                 heights = c(1,.35),
                 ncol=1, nrow=2,align="v")

png('./FIGURES/Acutestressregime_example.png',height=800,width=1000)
plot1
dev.off()


######### STEP 5.0: CALCULATE TOTAL NUMBER OF SHERDS IN SAMPLE DATA ###########
df <- fv_ceramics
df[is.na(df)] <- 0

### Limit survey data to diagnostic types (Wetherill considered as Mancos BW)
df_diag <- subset(df, select=c(Site.Number,Chapin.Gray,Moccasin.Gray,Mancos.Gray,
                               Mancos.Corrugated,Dolores.Corrugated,Mesa.Verde.Corrugated,
                               Chapin.B.W,Piedra.B.W,Cortez.B.W,
                               McElmo.B.W,Mesa.Verde.B.W,Abajo.R.O,
                               Bluff.B.R,Deadman.s.B.R))%>%
  mutate(Mancos.B.W = df$Mancos.B.W+df$Wetherill.B.W)%>%
  replace(is.na(.),0)%>%
  group_by(Site.Number) %>% 
  summarise_each(funs(sum))


### reorder so types are ordered chronologically by ware
df_diag <- df_diag[, c(1,2,3,4,5,6,7,8,9,10,16,11,12,13,14,15)]

######### STEP 5.1: CALCULATE PROPORTION OF SHERDS IN SAMPLE DATA #############
rownames(df_diag) <-NULL
df_diag <- column_to_rownames(df_diag,"Site.Number")
df_allprop <-prop.table(as.matrix(df_diag),1)*100

######### STEP 6.0: MEASURE DIFFERENCE BETWEEN SURVEY DATA AND CALIBRATION #### 
### Merge calibration and sample data
test <-rbind(cal_allprop,df_allprop)

### Run BR and extract relevant data
test_br_out <- as.data.frame(BR(test))
test_br_out <- test_br_out[11:nrow(test),1:10]

### normalize
test_br_zscore<- as.data.frame(matrix(NA,nrow(df_diag),nrow(cal)))
colnames(test_br_zscore) <- 1:10
rownames(test_br_zscore) <- rownames(df_diag)

n <- nrow(test_br_zscore)
n2 <- ncol(test_br_zscore)

for(i in 1:n){
  for(j in 1:n2){
    test_br_zscore[i,j] <- (test_br_out[i,j]-rowMeans(test_br_out[i,]))/sd(test_br_out[i,])
  }
}

######### STEP 7.0: TURN INTO TIME SERIES  ####################################
test_br_zscore <- rownames_to_column(test_br_zscore, "Site.Number")

test_br_zscore_ts<- test_br_zscore%>%
  gather(key="Timeperiod", value = "Value",`1`:`10`)%>%
  mutate(range = ranges$ranges[match(Timeperiod,ranges$period)])%>%
  mutate(start=substr(range,1L,4L))%>%
  mutate(end=substr(range,6L,9L))%>%
  rowwise() %>%
  do(data.frame(Site = .$Site.Number, phase = .$Timeperiod, range = .$range, start =.$start, end= .$end,
                Value=.$Value,Year=seq(.$start,.$end,by=1)))

######### STEP 8.0: CALCULATE NUMBER OF ROOMS #################################
fv_rmct <- as.data.frame(fv_rm)
fv_rmct[is.na(fv_rmct)] <- 0

### calculate number of rooms based on Pueblo Decomposition Model 
fv_rmct_agg <- fv_rmct%>%
  mutate(rb_vol = (Shape_Area*Height)*.75)%>%
  mutate(rm_ct = ceiling(rb_vol/9.68))%>%
  mutate(rm_ct_down = ceiling(rb_vol/(9.68/.65)))%>%
  subset(select=c(2,9:11))%>%
  group_by(newsiten_1) %>% 
  summarise_each(funs(sum))

######### STEP 8.1: INCORPORATE NUMBER OF KIVAS ###############################
fv_kv <- as.data.frame(fv_kiva)
kv_tbl <- table(fv_kv['newsitenum'])

######### STEP 8.2: LINK ROOM & KIVA COUNT WITH OCCUPATION EST. ############### 
### Join kiva count with room count
fv_rmct_agg$kv_ct <- kv_tbl[match(fv_rmct_agg$newsiten_1, names(kv_tbl))]
colnames(fv_rmct_agg)[1]<-"Site.Number"

### Join occupation estimates with room and kiva counts
fv_occ <- merge(test_br_zscore,fv_rmct_agg,by="Site.Number")

### Build into time series
fv_comp_stand_ts<- fv_occ%>%
  gather(key="Timeperiod", value = "Value",`1`:`10`)%>%
  mutate(range = ranges$ranges[match(Timeperiod,ranges$period)])%>%
  mutate(start=substr(range,1L,4L))%>%
  mutate(end=substr(range,6L,9L))

######### STEP 9.0: ASSIGN RM & KIVA CT. & OCC. BASED ON SURVEY/EXC. DATA #####
### Far View Tower 
fv_comp_stand_ts$rm_ct_down[fv_comp_stand_ts$Site.Number == "5MV00810"] <- 16

### Megalithic House
fv_comp_stand_ts$rm_ct_down[fv_comp_stand_ts$Site.Number == "5MV00790"] <- 7

### One Clan House
fv_comp_stand_ts$rm_ct_down[fv_comp_stand_ts$Site.Number == "5MV00835"] <- 13

### 5MV864
fv_comp_stand_ts$rm_ct_down[fv_comp_stand_ts$Site.Number == "5MV00864"] <- 10 

### 5MV836
fv_comp_stand_ts$rm_ct_down[fv_comp_stand_ts$Site.Number == "5MV00836"] <- 4

### 5MV837
fv_comp_stand_ts$rm_ct_down[fv_comp_stand_ts$Site.Number == "5MV00837"] <- 4

### 5MV807
fv_comp_stand_ts$rm_ct_down[fv_comp_stand_ts$Site.Number == "5MV00807"] <- 3

### Far View House
# room count
fv_comp_stand_ts$rm_ct_down[fv_comp_stand_ts$Site.Number == "5MV00808" & fv_comp_stand_ts$Timeperiod <= 2] <- 0
fv_comp_stand_ts$rm_ct_down[fv_comp_stand_ts$Site.Number == "5MV00808" & fv_comp_stand_ts$Timeperiod == 3] <- 18
fv_comp_stand_ts$rm_ct_down[fv_comp_stand_ts$Site.Number == "5MV00808" & fv_comp_stand_ts$Timeperiod == 4] <- 32
fv_comp_stand_ts$rm_ct_down[fv_comp_stand_ts$Site.Number == "5MV00808" & fv_comp_stand_ts$Timeperiod == 5] <- 32
fv_comp_stand_ts$rm_ct_down[fv_comp_stand_ts$Site.Number == "5MV00808" & fv_comp_stand_ts$Timeperiod == 6] <- 47
fv_comp_stand_ts$rm_ct_down[fv_comp_stand_ts$Site.Number == "5MV00808" & fv_comp_stand_ts$Timeperiod == 7] <- 47
fv_comp_stand_ts$rm_ct_down[fv_comp_stand_ts$Site.Number == "5MV00808" & fv_comp_stand_ts$Timeperiod == 8] <- 57
fv_comp_stand_ts$rm_ct_down[fv_comp_stand_ts$Site.Number == "5MV00808" & fv_comp_stand_ts$Timeperiod == 9] <- 57
fv_comp_stand_ts$rm_ct_down[fv_comp_stand_ts$Site.Number == "5MV00808" & fv_comp_stand_ts$Timeperiod == 10] <- 0
# kiva count
fv_comp_stand_ts$kv_ct[fv_comp_stand_ts$Site.Number == "5MV00808" & fv_comp_stand_ts$Timeperiod <= 2] <- 0
fv_comp_stand_ts$kv_ct[fv_comp_stand_ts$Site.Number == "5MV00808" & fv_comp_stand_ts$Timeperiod == 3] <- 2
fv_comp_stand_ts$kv_ct[fv_comp_stand_ts$Site.Number == "5MV00808" & fv_comp_stand_ts$Timeperiod == 4] <- 2
fv_comp_stand_ts$kv_ct[fv_comp_stand_ts$Site.Number == "5MV00808" & fv_comp_stand_ts$Timeperiod == 5] <- 2
fv_comp_stand_ts$kv_ct[fv_comp_stand_ts$Site.Number == "5MV00808" & fv_comp_stand_ts$Timeperiod == 6] <- 2
fv_comp_stand_ts$kv_ct[fv_comp_stand_ts$Site.Number == "5MV00808" & fv_comp_stand_ts$Timeperiod == 7] <- 2
fv_comp_stand_ts$kv_ct[fv_comp_stand_ts$Site.Number == "5MV00808" & fv_comp_stand_ts$Timeperiod == 8] <- 3
fv_comp_stand_ts$kv_ct[fv_comp_stand_ts$Site.Number == "5MV00808" & fv_comp_stand_ts$Timeperiod == 9] <- 4
fv_comp_stand_ts$kv_ct[fv_comp_stand_ts$Site.Number == "5MV00808" & fv_comp_stand_ts$Timeperiod == 10] <- 0
# occupation
fv_comp_stand_ts$Value[fv_comp_stand_ts$Site.Number == "5MV00808" & fv_comp_stand_ts$Timeperiod == 3] <- 1
fv_comp_stand_ts$Value[fv_comp_stand_ts$Site.Number == "5MV00808" & fv_comp_stand_ts$Timeperiod == 4] <- 1
fv_comp_stand_ts$Value[fv_comp_stand_ts$Site.Number == "5MV00808" & fv_comp_stand_ts$Timeperiod == 8] <- 1
fv_comp_stand_ts$Value[fv_comp_stand_ts$Site.Number == "5MV00808" & fv_comp_stand_ts$Timeperiod == 9] <- 1

### Pipe Shrine House
# room count
fv_comp_stand_ts$rm_ct_down[fv_comp_stand_ts$Site.Number == "5MV00809" & fv_comp_stand_ts$Timeperiod <= 3] <- 0
fv_comp_stand_ts$rm_ct_down[fv_comp_stand_ts$Site.Number == "5MV00809" & fv_comp_stand_ts$Timeperiod == 4] <- 10
fv_comp_stand_ts$rm_ct_down[fv_comp_stand_ts$Site.Number == "5MV00809" & fv_comp_stand_ts$Timeperiod == 5] <- 10
fv_comp_stand_ts$rm_ct_down[fv_comp_stand_ts$Site.Number == "5MV00809" & fv_comp_stand_ts$Timeperiod == 6] <- 21
fv_comp_stand_ts$rm_ct_down[fv_comp_stand_ts$Site.Number == "5MV00809" & fv_comp_stand_ts$Timeperiod == 7] <- 21
fv_comp_stand_ts$rm_ct_down[fv_comp_stand_ts$Site.Number == "5MV00809" & fv_comp_stand_ts$Timeperiod >= 8] <- 0
# kiva count
fv_comp_stand_ts$kv_ct[fv_comp_stand_ts$Site.Number == "5MV00809" & fv_comp_stand_ts$Timeperiod <= 3] <- 0
fv_comp_stand_ts$kv_ct[fv_comp_stand_ts$Site.Number == "5MV00809" & fv_comp_stand_ts$Timeperiod == 4] <- 1
fv_comp_stand_ts$kv_ct[fv_comp_stand_ts$Site.Number == "5MV00809" & fv_comp_stand_ts$Timeperiod == 5] <- 1
fv_comp_stand_ts$kv_ct[fv_comp_stand_ts$Site.Number == "5MV00809" & fv_comp_stand_ts$Timeperiod == 6] <- 1
fv_comp_stand_ts$kv_ct[fv_comp_stand_ts$Site.Number == "5MV00809" & fv_comp_stand_ts$Timeperiod == 7] <- 1
fv_comp_stand_ts$kv_ct[fv_comp_stand_ts$Site.Number == "5MV00809" & fv_comp_stand_ts$Timeperiod >= 8] <- 0

### Coyote Village
# room count
fv_comp_stand_ts$rm_ct_down[fv_comp_stand_ts$Site.Number == "5MV00820" & fv_comp_stand_ts$Timeperiod <= 2] <- 0
fv_comp_stand_ts$rm_ct_down[fv_comp_stand_ts$Site.Number == "5MV00820" & fv_comp_stand_ts$Timeperiod == 3] <- 19
fv_comp_stand_ts$rm_ct_down[fv_comp_stand_ts$Site.Number == "5MV00820" & fv_comp_stand_ts$Timeperiod == 4] <- 21
fv_comp_stand_ts$rm_ct_down[fv_comp_stand_ts$Site.Number == "5MV00820" & fv_comp_stand_ts$Timeperiod == 5] <- 27
fv_comp_stand_ts$rm_ct_down[fv_comp_stand_ts$Site.Number == "5MV00820" & fv_comp_stand_ts$Timeperiod == 6] <- 27
fv_comp_stand_ts$rm_ct_down[fv_comp_stand_ts$Site.Number == "5MV00820" & fv_comp_stand_ts$Timeperiod == 7] <- 35
fv_comp_stand_ts$rm_ct_down[fv_comp_stand_ts$Site.Number == "5MV00809" & fv_comp_stand_ts$Timeperiod >= 8] <- 0
# kiva count
fv_comp_stand_ts$kv_ct[fv_comp_stand_ts$Site.Number == "5MV00820" & fv_comp_stand_ts$Timeperiod <= 2] <- 0
fv_comp_stand_ts$kv_ct[fv_comp_stand_ts$Site.Number == "5MV00820" & fv_comp_stand_ts$Timeperiod == 3] <- 3
fv_comp_stand_ts$kv_ct[fv_comp_stand_ts$Site.Number == "5MV00820" & fv_comp_stand_ts$Timeperiod == 4] <- 3
fv_comp_stand_ts$kv_ct[fv_comp_stand_ts$Site.Number == "5MV00820" & fv_comp_stand_ts$Timeperiod == 5] <- 3
fv_comp_stand_ts$kv_ct[fv_comp_stand_ts$Site.Number == "5MV00820" & fv_comp_stand_ts$Timeperiod == 6] <- 3
fv_comp_stand_ts$kv_ct[fv_comp_stand_ts$Site.Number == "5MV00820" & fv_comp_stand_ts$Timeperiod == 7] <- 5
fv_comp_stand_ts$kv_ct[fv_comp_stand_ts$Site.Number == "5MV00809" & fv_comp_stand_ts$Timeperiod >= 8] <- 0

### 5MV866
# room count
fv_comp_stand_ts$rm_ct_down[fv_comp_stand_ts$Site.Number == "5MV00866" & fv_comp_stand_ts$Timeperiod == 5] <- 10
fv_comp_stand_ts$rm_ct_down[fv_comp_stand_ts$Site.Number == "5MV00866" & fv_comp_stand_ts$Timeperiod == 6] <- 10
fv_comp_stand_ts$rm_ct_down[fv_comp_stand_ts$Site.Number == "5MV00866" & fv_comp_stand_ts$Timeperiod == 7] <- 10
fv_comp_stand_ts$rm_ct_down[fv_comp_stand_ts$Site.Number == "5MV00866" & fv_comp_stand_ts$Timeperiod == 8] <- 1
# kiva count
fv_comp_stand_ts$kv_ct[fv_comp_stand_ts$Site.Number == "5MV00866" & fv_comp_stand_ts$Timeperiod == 5] <- 3
fv_comp_stand_ts$kv_ct[fv_comp_stand_ts$Site.Number == "5MV00866" & fv_comp_stand_ts$Timeperiod == 6] <- 3
fv_comp_stand_ts$kv_ct[fv_comp_stand_ts$Site.Number == "5MV00866" & fv_comp_stand_ts$Timeperiod == 7] <- 3
fv_comp_stand_ts$kv_ct[fv_comp_stand_ts$Site.Number == "5MV00866" & fv_comp_stand_ts$Timeperiod == 8] <- 0

### 5MV875
# room count
fv_comp_stand_ts$rm_ct_down[fv_comp_stand_ts$Site.Number == "5MV00875" & fv_comp_stand_ts$Timeperiod == 5] <- 17
fv_comp_stand_ts$rm_ct_down[fv_comp_stand_ts$Site.Number == "5MV00875" & fv_comp_stand_ts$Timeperiod == 6] <- 17
fv_comp_stand_ts$rm_ct_down[fv_comp_stand_ts$Site.Number == "5MV00875" & fv_comp_stand_ts$Timeperiod == 7] <- 0
fv_comp_stand_ts$rm_ct_down[fv_comp_stand_ts$Site.Number == "5MV00875" & fv_comp_stand_ts$Timeperiod == 8] <- 15
# kiva count
fv_comp_stand_ts$kv_ct[fv_comp_stand_ts$Site.Number == "5MV00875" & fv_comp_stand_ts$Timeperiod == 5] <- 3
fv_comp_stand_ts$kv_ct[fv_comp_stand_ts$Site.Number == "5MV00875" & fv_comp_stand_ts$Timeperiod == 6] <- 3
fv_comp_stand_ts$kv_ct[fv_comp_stand_ts$Site.Number == "5MV00875" & fv_comp_stand_ts$Timeperiod == 7] <- 0
fv_comp_stand_ts$kv_ct[fv_comp_stand_ts$Site.Number == "5MV00875" & fv_comp_stand_ts$Timeperiod == 8] <- 1

### 5MV499
# room count
fv_comp_stand_ts$rm_ct_down[fv_comp_stand_ts$Site.Number == "5MV00499" & fv_comp_stand_ts$Timeperiod == 1] <- 1
fv_comp_stand_ts$rm_ct_down[fv_comp_stand_ts$Site.Number == "5MV00499" & fv_comp_stand_ts$Timeperiod >= 5] <- 12
# kiva count
fv_comp_stand_ts$kv_ct[fv_comp_stand_ts$Site.Number == "5MV00499" & fv_comp_stand_ts$Timeperiod == 1] <- 1
fv_comp_stand_ts$kv_ct[fv_comp_stand_ts$Site.Number == "5MV00499" & fv_comp_stand_ts$Timeperiod >= 5] <- 2


######### STEP 9.1: DETERMINE POPULATION BASED ON ROOM & KIVA COUNTS ##########
### population estimate based on room counts and kiva counts
fv_comp_stand_ts <- fv_comp_stand_ts%>%
  mutate(occupied = ifelse(Value > 0,1,
                           ifelse(Value < 0, 0,NA)))%>%
  mutate(rf_area = rm_ct_down*6.787)%>%
  mutate(pop_min = (ceiling(rf_area/10))*occupied)%>%
  mutate(kv_pop = (kv_ct*6)*occupied)

fv_comp_stand_ts[is.na(fv_comp_stand_ts)] <- 0

### Built into time series
fv_occ_ts <- fv_comp_stand_ts %>%
  rowwise() %>%
  do(data.frame(Value=.$Value,Year=seq(.$start,.$end,by=1),sitenum=.$Site.Number,pop_min_pdm=.$pop_min,pop_min_kvct=.$kv_pop))

######### STEP 9.2: AGGREGATE TO DETERMINE COMMUNITY OCCUPATION PATTERN #######
# community wide 
fv_occ_ts_c <- subset(fv_occ_ts,select=c("Year","pop_min_pdm","pop_min_kvct"))
fv_occ_ts_c <- aggregate(. ~ Year, data=fv_occ_ts_c, FUN=sum)
fv_occ_ts_c$pop_min_av <- (fv_occ_ts_c$pop_min_pdm+fv_occ_ts_c$pop_min_kvct)/2

######### STEP 10.0: PLOT COMMUNITY PATTERNS ###################################
settlement <- read.csv("./OUTPUT/all_community_occupations.csv",header=T,fileEncoding = 'UTF-8-BOM')
acutestress <- read.csv("./OUTPUT/acute_stress_periods.csv",header=T,fileEncoding = 'UTF-8-BOM')

### fill in some plotting values
settlement[is.na(settlement)] <- 0
settlement[1,2:10] <- 0

### extract for plotting easier
acutestress_fv <- subset(acutestress, comm == "fv")
acutestress_more <- subset(acutestress, comm == "more")
acutestress_yjp <- subset(acutestress, comm == "yjp")

### plot
p1<-ggplot()+
  geom_col(data = acutestress_fv, aes(x = Year, y =400), fill = "coral1",width=1,alpha=0.4)+
  geom_line(data = settlement, aes(x = Year, y = fv_pop_rmblk), col = "grey70",size=.75)+
  geom_line(data = settlement, aes(x = Year, y = fv_pop_kv), col = "grey70",size=.75)+
  geom_line(data = settlement, aes(x = Year, y = fv_pop_av), col = "grey30",size=1.75)+
  scale_x_continuous(limits=c(595,1305),breaks=seq(600,1300,by=50))+
  scale_y_continuous(limits=c(0,410),breaks=seq(0,400,by=100))+
  labs(y = "Population",
       title= "Far View Community")+
  theme(legend.position="none",
        axis.ticks.length.x=unit(.25, "cm"),
        axis.text.y=element_text(size=25,color="black"),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y = element_text(size=32,color="black"),
        panel.grid.minor.y = element_blank(),
        title=element_text(size=32,color="black"))

p2<-ggplot()+
  geom_col(data = acutestress_more, aes(x = Year, y =400), fill = "coral1",width=1,alpha=0.4)+
  geom_line(data = settlement, aes(x = Year, y = m_pop_rmblk), col = "grey70",size=.75)+
  geom_line(data = settlement, aes(x = Year, y = m_pop_kv), col = "grey70",size=.75)+
  geom_line(data = settlement, aes(x = Year, y = m_pop_av), col = "grey30",size=1.75)+
  scale_x_continuous(limits=c(595,1305),breaks=seq(600,1300,by=50))+
  scale_y_continuous(limits=c(0,410),breaks=seq(0,400,by=100))+
  labs(y = "Population",
       title= "Morefield Canyon Great House Village")+
  theme(legend.position="none",
        axis.ticks.length.x=unit(.25, "cm"),
        axis.text.y=element_text(size=25,color="black"),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y = element_text(size=32,color="black"),
        panel.grid.minor.y = element_blank(),
        title=element_text(size=32,color="black"))


p3<-ggplot()+
  geom_col(data = acutestress_yjp, aes(x = Year, y = 1000), fill = "coral1",width=1,alpha=0.4)+
  geom_line(data = settlement, aes(x = Year, y = yjp_pop_rmblk), col = "grey70",size=.75)+
  geom_line(data = settlement, aes(x = Year, y = yjp_pop_kv), col = "grey70",size=.75)+
  geom_line(data = settlement, aes(x = Year, y = yjp_pop_av), col = "grey30",size=1.75)+
  scale_x_continuous(limits=c(595,1305),breaks=seq(600,1300,by=50))+
  scale_y_continuous(limits=c(0,1010),breaks=seq(0,1000,by=200))+
  labs(y="Population",
       x= "",
       title="Yellow Jacket Pueblo")+
  theme(legend.position="none",
        axis.ticks.length.x=unit(.25, "cm"),
        axis.text=element_text(size=25,color="black"),
        axis.title = element_text(size=32,color="black"),
        title=element_text(size=32,color="black"),
        panel.grid.minor.y = element_blank())



plot2<-ggarrange(p1,p2,p3,
                 ncol=1, nrow=3,align="v")

png('./FIGURES/Acutestress&occupation_example.png',height=1400,width=1200)
plot2
dev.off()

