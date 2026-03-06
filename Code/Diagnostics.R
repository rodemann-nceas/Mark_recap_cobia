#########################################################
# Diagnostic plots for Cobia AT data
# Date Last Updated: 3/6/2026
# Author: J. Rodemann
# Purpose: Create plots to visualize movement data after download to determine:
#   - Efficiency of array
#   - Quality of Data
#   - Initial movement and habitat use patterns
##########################################################

## Filtering the raw data and combining it ##
library(stringr)
library(tidyverse)
#library(VTrack)
ComputeDistance <- function(x,y,x1,y1){
  dist <- sqrt((x-x1)^2+(y-y1)^2)
  return(dist)
}

# The following code pulls in downloads from 2 different dates
# If code does not work, there is a chance that one or more of the csvs do not have the same heading for Date
adata_fwcrepro <- read.csv('Data/COBREPRO_OBIS_20260226.csv')
adata_fwcgom <- read.csv('Data/COBGOM_AcousticDetections_20260406.csv')
str(adata_fwcgom)
head(adata_fwcgom)

#change Date_time into posixct
adata_fwcgom$datecollected <- as.POSIXct(adata_fwcgom$dateCollectedUTC, format='%Y-%m-%d %H:%M:%S')
odata$Date.and.Time..UTC. <- as.POSIXct(odata$Date.and.Time..UTC., format='%Y-%m-%d %H:%M:%S')


#add columns to show lag lat and long and calculate distance btw points
dat_gom <- adata_fwcgom %>% 
  arrange(datecollected) %>% #sorts it chronologically
  group_by(catalognumber) %>% #group by tag number
  mutate(llon=lag(longitude), llat=lag(latitude)) %>% #last position as a column - where you were before and where you are now
  rowwise() %>% 
  filter(!is.na(llon)) %>% #get rid of na rows for llon
  mutate(dist=ComputeDistance(llat, latitude, llon, longitude)) %>% 
  filter(!is.na(dist)) %>% 
  filter(dist > 0.01) #filter out times when individual did not change position - think critically!!!
 
#calculate total distance moved by each fish in Km according to acoustics
dist_fish <- dat_fish1 %>% 
  group_by(Transmitter) %>% 
  dplyr::summarise(distance = sum(dist))

head(dist_fish)
## Determine efficiency of array ##

##bar plot of detections per receiver (use dataset from before)
#add latitude and longitude to dataset and organize - just for Rankin
meta_rec <- rec_meta %>% dplyr::select(-Station)
meta_rec <- meta_rec %>% dplyr::rename(Station=Station_number)
sta <- merge(meta_rec, Trouts, by='Station')
str(sta)
sta <- sta %>% dplyr::select(Station, Latitude, Longitude, n)
sta <- sta %>% dplyr::distinct()
str(sta)
sta <- sta %>% 
  mutate(Order = as.numeric(gsub("RBr", "", Station))) %>% 
  arrange(Order) %>% 
  mutate(Station = factor(Station, levels = Station))
str(sta)

#bar plot of detections per receiver
bar_rec <- ggplot(sta, aes(x = Station, y = n))+
  geom_bar(stat = 'identity')+
  theme_classic()+
  ylab("Number of Detections")+
  ggtitle("Detections per receiver")+
  xlab("Station")+
  theme(axis.text.x=element_text(angle =- 90, vjust = 0.5))
bar_rec
ggsave(filename = 'E:/FIU/Project/Downloads/June_2020/Diagnostic/barplot_receiver.jpeg', bar_rec)

##put on map
require(rgdal)
require(rgeos)
require(adehabitatHR)
require(mapproj)
require(sf)
require(rnaturalearth)
require(rnaturalearthdata)
require(igraph)
require(ggraph)
require(tidygraph)
require(RColorBrewer)
install.packages(c('mapproj', 'igraph', 'ggraph', 'tidygraph'))


world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

data <- data %>% group_by(decimalLatitude, decimalLongitude) %>% mutate(n = n()) %>% ungroup()
data$n <- as.numeric(data$n)
data <- data %>% mutate(size = n/1000)
data <- data %>% mutate(month = month(datecollected))
data$month <- as.factor(data$month)
data$season <- as.character(NA)

for(i in 1:nrow(data)){
  if(data$month[i] %in% c(1,2,12)){
    data$season[i] <- 'Winter'
  }else if(data$month[i] %in% c(3,4,5)){
    data$season[i] <- 'Spring'
  }else if(data$month[i] %in% c(6,7,8)){
    data$season[i] <- "Summer"
  }else{
    data$season[i] <- 'Fall'
  }
}
map <- ggplot(data=world)+
  geom_sf()+
  coord_sf(xlim = c(-88, -77), ylim = c(24, 31), expand = FALSE)+
  geom_point(data=data, 
             aes(decimalLongitude, decimalLatitude, fill = month), inherit.aes=F, pch=21, size = log(data$n)+2)+
  scale_fill_brewer(type = 'div', palette = 'RdBu')+
  facet_wrap(~catalogNumber)+
  guides(legend = 'none', color = guide_legend(override.aes = list(size = 5)))
map

# 1. Extract a list of unique catalog numbers
individuals <- unique(data$catalogNumber)

# 2. Open the PDF graphics device
# This creates the file and gets it ready to receive pages
pdf(file = "Plots/All_Individuals_Maps_season.pdf", width = 8, height = 6)

# 3. Loop through each individual
for (ind in individuals) {
  
  # Subset the data
  ind_data <- data %>% filter(catalogNumber == ind)
  
  # Build the plot
  p <- ggplot(data = world) +
    geom_sf() +
    coord_sf(xlim = c(-88, -77), ylim = c(24, 31), expand = FALSE) +
    geom_point(data = ind_data, 
               aes(x = decimalLongitude, y = decimalLatitude, fill = season, size = log(n) + 2), 
               inherit.aes = FALSE, pch = 21) +
    scale_fill_brewer(type = 'div', palette = 'RdBu') +
    scale_size_identity() +
    guides(fill = guide_legend(override.aes = list(size = 5))) +
    # Added a title so you know which individual is on which page!
    ggtitle(paste("Individual:", ind)) 
  
  # Print the plot to the open PDF device (this creates a new page)
  print(p)
}

# 4. Close the graphics device to finish and save the file
# If you forget this step, the PDF might be corrupted or empty!
dev.off()
ggsave(filename = 'E:/FIU/Project/Downloads/June_2020/Diagnostic/map_receiver.jpeg', map)

##bar plot of detections per transmitter (dataset from before)
bar_trans <- ggplot(trouts, aes(x = Transmitter, y = n))+
  geom_bar(stat = 'identity')+
  theme_classic()+
  ylab("Number of Detections")+
  ggtitle("Detections per Transmitter")+
  xlab("Transmitter")+
  theme(axis.text.x=element_text(angle =- 90, vjust = 0.5))
bar_trans
ggsave(filename = 'E:/FIU/Project/Downloads/June_2020/Diagnostic/bar_transmitter.jpeg', bar_trans)

##overall abacus plot
#load in RDS file
data <- readRDS(file = "E:/FIU/Project/Downloads/June_2020/Data/Acoustic/trout_June_2020.rds")
str(data)

library(splitstackshape) ###package to use cSplit
library(chron)

# Break up the new datetime column into  a date, and a time
#Jon's Data
data <- adata_fwcgom
data$Date_Time2 <- data$datecollected

data <- cSplit(data, "Date_Time2", sep = " ", type.convert = F)

head(data)


# Rename columns so they make sense, it'll just be your last 2 column numbers, in this case the 15th and 16th column s
colnames(data)[48:49]<-c("Date", "Time")
head(data)

# Then I repeat this and parse the date into year, month, and day (then hour, minute, second), so I can easily subset by year
# Copy date to parse out, maintain original date
data$Date2 <- data$Date


data<-cSplit(data, "Date2", sep = "-", type.convert = FALSE)
head(data)
colnames(data)[50:52]<-c("Year", "Month", "Day")
head(data)

data$Time2 <- data$Time
head(data)

data<-cSplit(data, "Time2", sep = ":", type.convert = FALSE)
head(data)
colnames(data)[53:55]<-c("Hour", "Minute", "Second")
head(data)

# Assign classes so data plays nice for plotting and summaries
data$Date <- as.Date(data$Date)
data$Time <- as.times(data$Time)
data$Year <- as.numeric(data$Year)
data$Month <- as.numeric(data$Month)
data$Day <- as.numeric(data$Day)
data$Hour <- as.numeric(data$Hour)
data$Minute <- as.numeric(data$Minute)
data$Second <- as.numeric(data$Second)
data$catalogNumber <- as.factor(data$catalogNumber)
summary(data)
head(data)
max(data$datecollected)
unique(data$catalognumber)

#graph
cobia_fwcgom <- ggplot(data,aes(datecollected, catalogNumber, group = 1)) + 
  geom_point(aes(color = decimalLatitude)) +
  scale_color_gradient(low="darkblue", high="lightgreen") +
  scale_x_date(date_breaks = "1 year", limits = as.Date(c("2017-04-01", "2025-06-01")))+
  labs(x = "Date",  y = "Transmitter") +
  # Change the angle and text size of x-axis tick labels to make more readable in final plots
  theme_bw()+
  theme(axis.text.x=element_text(angle= 50, size=10, vjust = 1, hjust=1),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.y=element_text(size=8))+
  labs(title = "FWC GOM Cobia detections")
cobia_fwcgom
ggsave(filename = 'E:/FIU/Project/Downloads/June_2020/Diagnostic/Overall_abacus.jpeg', trout_det)

#abacus plot with colors for sav
#bring in SAV data
SAV <- read.csv('E:/FIU/Project/Downloads/June_2020/SAV_Stations.csv')

SAV_dat <- merge(SAV, data, by='Station')
str(SAV_dat)

trout_det_SAV <- ggplot(SAV_dat,aes(Date, Transmitter))+ 
  geom_vline(aes(xintercept = (as.Date("2020-06-15"))), linetype=4, size=0.75, color="red")+    
  geom_vline(aes(xintercept = (as.Date("2020-02-29"))), linetype=4, size=0.75, color="red")+    
  geom_point(aes(colour = SAV)) +
  scale_color_gradient(low="darkblue", high="lightgreen") +
  scale_x_date(date_breaks = "1 month", limits = as.Date(c("2020-02-01", "2020-07-01")))+
  labs(x = "Date",  y = "Transmitter") +
  # Change the angle and text size of x-axis tick labels to make more readable in final plots
  theme_bw()+
  theme(axis.text.x=element_text(angle= 50, size=10, vjust = 1, hjust=1),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.y=element_text(size=8))+
  scale_color_manual(values = c('black', 'red', 'blue', 'green'))+
  labs(title = "Trout Detections")
trout_det_SAV
ggsave(filename = 'E:/FIU/Project/Downloads/June_2020/Diagnostic/SAV_abacus.jpeg', trout_det_SAV)

##Kernel Densities
#organize data to be read by adehabitatHR
data2 <- data
data3 <- data2 %>% dplyr::select(long, lat, Transmitter, Station)
data4 <- data3 %>% group_by(Transmitter) %>% summarize(n=n_distinct(Station))
data4
data5 <- data4 %>% subset(n>1)
data6 <- merge(data5, data2, by="Transmitter")
str(data6)
data7 <- data6 %>% filter(!is.na(n))
data_ker2 <- data7 %>% dplyr::select(lat, long, Transmitter)

data_ker <- data7 %>% dplyr::select(lat, long, Transmitter) %>% droplevels()
str(data_ker)
data_ker <- data_ker %>% dplyr::rename(id=Transmitter)
coordinates(data_ker) <- ~long+lat
proj4string(data_ker)<-CRS("+proj=longlat +datum=WGS84")
data_ker2<-spTransform(data_ker, CRS("+proj=utm +zone=17 +datum=NAD83 +units=m +no_defs"))
head(as.data.frame(data_ker2))
str(data_ker2)


#run kernels
kernel_fish <- kernelUD(data_ker[,1], h = 'href', grid=100, extent=1)
kernel_fish
image(kernel_fish)
kernel_vol <- getvolumeUD(kernel_fish)
kernel_poly <- getverticeshr(kernel_fish, 25)
print(kernel_poly)
image(kernel_vol)
plot(kernel_poly)
kernel_poly@data$id <- as.factor(kernel_poly@data$id)

pl_names <- names(kernel_fish)
ud <- lapply(kernel_fish, function(x) try(getverticeshr(x, 25)))
sapply(1:length(ud), function(i) 
  row.names(ud[[i]]) <<- pl_names[i]
)
plf_poly <- Reduce(rbind, ud)
plot(plf_poly)
df <- fortify(plf_poly)
str(df)

g <- ggplot(df, aes(x=long, y=lat, fill=id, group=group))+
  geom_polygon(alpha=.4)
g

ggsave(filename = 'E:/FIU/Project/Downloads/June_2020/Diagnostic/Kernel_density_individual.jpeg', g)

#population kernel density
dat_total <- data2 %>% dplyr::select(long, lat)
coordinates(dat_total) <- ~long+lat
proj4string(dat_total)<-CRS("+proj=longlat +datum=WGS84")


Kernel_total <- kernelUD(dat_total)
Kernel_total
image(Kernel_total)
ker_tot_vol <- getvolumeUD(Kernel_total)
image(ker_tot_vol)
ker_tot_poly <- getverticeshr(Kernel_total, 95)
plot(ker_tot_poly)
df2 <- fortify(ker_tot_poly)
g2 <- ggplot(df2, aes(x=long, y=lat))+
  geom_polygon(alpha=.4)
g2

ggsave(filename = 'E:/FIU/Project/Downloads/June_2020/Diagnostic/Kernel_density_pop.jpeg', g2)
