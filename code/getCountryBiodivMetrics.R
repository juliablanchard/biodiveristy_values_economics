### BII data form natural history museum



library(tidyverse)
library(readr)
library(dplyr)

# Load the data
bii_data <- readRDS("data/long_data.rds")
bii_data <- read.csv("data/resource.csv")
# View the first few rows
head(bii_data)

# Check the structure of the dataset
str(bii_data)
bii_data <- bii_data %>%
  mutate(area_code = gsub("-", ".", area_code))

library(jsonlite)


# Read area code JSON, countries only
area_mapping<- fromJSON("data/bii.json")[c(1:128,153:279),]
# Check structure of the JSON
str(area_mapping)
names(area_mapping)<-c("area_code","area_name","parent_code")

# add area_names
bii_data_named <- bii_data %>%
  left_join(area_mapping, by = "area_code")


# Choose a few example countries
selected_countries <- c("Brazil", "Kenya", "India", "Indonesia", "Australia")

bii_subset <- bii_data_named %>%
  filter(area_name %in% selected_countries)

ggplot(bii_subset, aes(x = year, y = value, color = scenario)) +
  geom_point() +
  facet_wrap(~ area_name, scales = "free_y") +
  labs(title = "BII Time Series by Country",
       x = "Year", y = "Biodiversity Intactness Index (%)") +
  theme_minimal()




#### Code to aggregate biodiversity metrics (from Pereira et al. 2023) to country level
### Using data and extraction scripts from Pereira et al. (2024). Global trends and scenarios for terrestrial biodiversity and ecosystem services from 1900-2050. Science https://doi.org/science.adn3441
### Zonal statistics ----
### Calculate zonal statistics per Country for all netCDFs
#### Output 1: Values_IPBES-regions_raw-changes.xlsx - all models: mean, sum, median, area-weighted sum and area-weighted mean
#### Output 2: Values_IPBES-regions_delta-changes.xlsx - ESS models only: sum and area-weighted sum are calculated from the absolute values
####     for those years (when they exist), by first summing all cell values for each year, and then calculating the change
### Project UTAS Biodiversity & Economics Workshop (for Natalie Stoeckl et al.)
### Created May 2025, Julia Blanchard


### 1 - Initializations ----
# clear workspace
rm(list=ls())

#libraries
library(terra)
library(dplyr)
library(ebvcube)
#library(xlsx)
library(stringr)
#Base maps
library(rnaturalearth)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("rhdf5")

#setting working directory to the current file source location 
root<-setwd("/Users/juliab6/Library/CloudStorage/OneDrive-UniversityofTasmania/biodiversity_economics_workshop/biodiveristy_values_economics/")  #only works R studio


root_output<-file.path(root,"output")

#list of nc files in directory with biodiversity model outputs as ebvcubes (they are downloaded below)

netcdfs_path <- list.files(file.path(root,"data"), patter='*nc', full.names = T)


# Load NetCDF file
ebv_raster <- rast(netcdfs_path[1])  # Multi-layer SpatRaster

# Check layer names (they should reflect years, e.g., "X2001", "X2002", etc.)
names(ebv_raster)

ebv_raster[[2]]





#functions
areaweight_mean<-function(raster1)
{
  #produce a raster that has the area of each degree cell
  a<-cellSize(raster1)
  a <- mask(a,raster1)
  #plot(a)
  #calculate area weighted mean
  asum<-global(a, "sum", na.rm=TRUE)$sum
  w<-a/asum
  #plot(raster1*w)
  areaw_mean<-global(raster1*w, "sum", na.rm=TRUE)$sum
  return(list(Area_W_Mean=areaw_mean, Area_W_Sum=areaw_mean*asum))
}



calc_delta_lq<-function(df,year0,year1,stats)
{
  (df[[stats]][df$Years==year1][[1]]-df[[stats]][df$Years==year0][[1]])/df[[stats]][df$Years==year0][[1]]*100
}

calc_delta_cube_lq<-function(es)
{
  es2<-es
  es2$Metric <- str_replace(es2$Metric," \\(.*", "")
  metrics <- unique(es2$Metric)
  scenarios <- unique(es2$Scenario)
  y1900 <- setdiff(unique(es2$Years),c("2015","2050"))
  entities <- unique(es2$Entity)
  for (m in metrics){
    if(m %in% es$Metric){
      #print(m)
      for (s in scenarios){
        absol<-subset(es2, Metric==m & Scenario==s & Entity==e & ! str_detect(Units,'%'))
        index<-which(es2$Metric==m & es2$Scenario==s & es2$Entity==e & str_detect(es2$Units,'%')) #=="Percentage (%)"
        for (i in index)
        {
          if ((!es2$Years[i] %in% c("2015","2050")) | (length(y1900)==0 & (es2$Years[i]=="2015")))
          {
            es2$Sum[i] <- 0
            es2$Area_W_Sum[i] <- 0
          } else if ((es2$Years[i] =="2015") & (!length(y1900)==0))
          {
            es2$Sum[i] <- calc_delta_lq(absol,y1900,"2015","Sum")
            es2$Area_W_Sum[i] <- calc_delta_lq(absol,y1900,"2015","Area_W_Sum")
          } else
          {
            es2$Sum[i] <- calc_delta_lq(absol,"2015","2050","Sum")
            es2$Area_W_Sum[i] <- calc_delta_lq(absol,"2015","2050","Area_W_Sum")
          }
        }
      }
      
    }else{
      print(paste0('No change calculated because no absolute values available.'))
    }
    
  }
  
  return(es2)
}

### 2 - Import data ----

# load country shapefile

world_vec<-terra::vect(ne_countries(returnclass = "sf"))


#ipbes_vec <- terra::vect(file.path(root_countries, 'IPBES_subregions_simp_v6.shp'))

### 3 - Calculations ----
#prepare result data.frame
cnames <- c('Scenario', 'LUCC', 'Model', 'Entity', 'Region','Years', 'Metric','Units', 'Sum','Mean', 'Median','Area_W_Mean','Area_W_Sum')
result <- data.frame(matrix(ncol = length(cnames), nrow = 0))
colnames(result) <- cnames

# loop through netCDFs
for (nc.i in 1:length(netcdfs_path)){ #
  nc <- netcdfs_path[nc.i]
  #get relevant netCDF info 
  cubes <- ebv_datacubepaths(nc)
  dates <- ebv_properties(nc)@temporal$dates
  entities <- ebv_properties(nc)@general$entity_names
  model <- ebv_properties(nc)@general$title
  pos<-last(str_locate_all(model,"BES-SIM ")[[1]])
  model <- str_sub(model,pos+1,-2) 
  
  cat('---- processing file ', model, ' (', nc.i, '/', length(netcdfs_path), ')\n', sep = '')
  
  #loop through cubes
  for (c.i in 1:dim(cubes)[1]){
    c.path <- cubes[c.i, 1]
    units <- ebv_properties(nc, c.path)@ebv_cube$units
    cat('-> cube ', c.path, ' (', c.i, '/', dim(cubes)[1], ')\n', sep = '')
    
    #loop through entities
    for(e in entities){
      
      #loop through ipbes regions
      for(region in unique(ipbes_vec$IPBES_sub)){
        
        #subset the vector data
        region.vec <- ipbes_vec[ipbes_vec$IPBES_sub == region,]
        
        #build result data.frame
        part <- data.frame(matrix(ncol = length(cnames), nrow = length(dates)))
        colnames(part) <- cnames
        r.i <- 1
        
        #loop through dates
        for(ts.i in 1:length(dates)){
          
          #read data
          data <- ebv_read(
            filepath = nc,
            datacubepath = c.path, 
            entity = e, 
            timestep = ts.i,
            verbose = F,
            ignore_RAM = TRUE
          )
          
          
          #extract values -> no weights + decide on touches!
          temp.raster <- terra::rasterize(region.vec, data, touches=T) #touches true to not loose the small pixels (e.g. Oceania)
          
          #mask the subset
          region.raster <- terra::mask(data, temp.raster)
          
          #calculate the metrics
          mean <- terra::global(region.raster, 'mean', na.rm=TRUE)
          sum <-global(region.raster, "sum", na.rm=TRUE)
          #median <-global(region.raster, 'median', na.rm=TRUE)
          median <- median(as.array(region.raster),na.rm=T)
          stats <-areaweight_mean(region.raster)
          Area_W_Mean <- stats$Area_W_Mean
          Area_W_Sum <- stats$Area_W_Sum
          
          #create part data.frame
          part$Scenario[r.i] <-  substr(cubes[c.i, 2],1,12) 
          part$LUCC[r.i] <- substr(cubes[c.i, 2],13,16)
          part$Model[r.i] <- model
          part$Entity[r.i] <- e
          part$Region[r.i] <- region
          part$Years[r.i] <- str_sub(as.character(dates[ts.i]), 1,4)
          part$Metric[r.i] <- cubes[c.i, 3]
          part$Units[r.i] <- units
          part$Sum[r.i] <- sum
          part$Mean[r.i] <- mean
          part$Median[r.i] <- median
          part$Area_W_Mean[r.i] <- Area_W_Mean
          part$Area_W_Sum[r.i] <- Area_W_Sum
          
          r.i <- r.i +1
          
        }
        
        
        #add data to result data.frame
        result <- rbind(result, part)
      }
    }
    
  }
  
}

