#This code helps us check and format the raw data to be used for the trend analyses
library(ggplot2)

#For activity trend analyses only using the output from Kaleidoscope\
bat.data <- readxl::read_excel("Data/Raw/ActivityModels/AllYearActivitybyNight_2016-2025_SM2_dual_only_withcov.xlsx", sheet = "OccupancyKSP",.name_repair = ~ gsub("\\s+", "_", trimws(.x)))
bat.data$Year <- lubridate::year(bat.data$Night)
#remove data from 2015 and prior since it is not consistent or clean enough
bat.data <- bat.data[bat.data$Year>2015,]
#create column for jNight
bat.data$jNight <- as.numeric(format(as.Date(bat.data$Night), "%j"))
#update water nearby based on new water distance
bat.data$Water_Nearby <- ifelse(bat.data$DIST_WATER_M<100,0,1)

#uplaod transect data from manually verified results only
bat.transect.data <- readxl::read_excel("Data/Raw/AllYearActivitybyNight_2016-2025_BC_withCovariates.xlsx", sheet = "OccupancyManual", .name_repair = ~ gsub("\\s+", "_", trimws(.x)))
#---Make sure it is properly formatted and remove stationary
bat.transect.data <- bat.transect.data[bat.transect.data$Quadrant =="Transects",]
bat.transect.data$Year <- lubridate::year(bat.transect.data$Night)
#create column for jNight
bat.transect.data$jNight <- as.numeric(format(as.Date(bat.transect.data$Night), "%j"))
#update water nearby based on new water distance
bat.transect.data$Water_Nearby <- ifelse(bat.transect.data$DIST_WATER_M<100,0,1)
#make sure Tlength is numeric
bat.transect.data$Tlength <- as.numeric(bat.transect.data$Tlength)

#------CHECK OF THE DATA-----------
# Check how many sites per year-----------------------------------------------------------------
xtabs(~Region+Year, data=bat.data, exclude=NULL, na.action=na.pass)

#check the grid cells--------------------------------------------------------------------------
length(unique(bat.data$GRTS_Cell_ID))

#check which grid were sampled which years-----------------------------------------------------
xtabs(~GRTS_Cell_ID+Year, data=bat.data, exclude=NULL, na.action=na.pass)

# check quadrant sampling-----------------------------------------------------------------------
xtabs(~Year+Quadrant, data=bat.data, exclude=NULL, na.action=na.pass)
bat.data$Quadrant <- paste0(toupper(substr(bat.data$Quadrant,1,1)), substring(bat.data$Quadrant,2))
xtabs(~Year+Quadrant, data=bat.data, exclude=NULL, na.action=na.pass)

xtabs(~GRTS_Cell_ID+Quadrant, data=bat.data, exclude=NULL, na.action=na.pass)

#check which dates data is collected in--------------------------------------------------------------
bat.data$jNight <- lubridate::yday(bat.data$Night)
myplot <- ggplot(data=bat.data[ bat.data$Year > 2010,], aes(x=jNight, y=Year))+
  ggtitle("Julian nights of data collection")+
  geom_point(position=position_jitter(h=.1, w=.1))
myplot
#truncate any data from after julian day 240 (End of August)
bat.data <- bat.data[ bat.data$jNight < 240,]
bat.data$jNightM15 <- bat.data$jNight - lubridate::yday("2020-05-14") #get the number of days since May 14. 

# Get general idea of where sampling occurred-------------------------------------------------------------
ggplot(data=bat.data, aes(x=Longitude, y=Latitude,color=as.factor(Year)))+
  ggtitle("Where collected")+
  geom_point( position=position_jitter(h=.1, w=.1), shape=1)

# Check covariates:Distance Clutter-------------------------------------------------------------------------------
xtabs(~`Distance_to_Clutter_(m)`,data=bat.data, exclude=NULL, na.action=na.pass)
#cap the maximum distance to 75m for modeling purposes
bat.data$`Distance_to_Clutter_(m)`[bat.data$`Distance_to_Clutter_(m)` > 75] <- 75
xtabs(~`Distance_to_Clutter_(m)`,data=bat.data, exclude=NULL, na.action=na.pass)
#get the percentage of rows with no distance to clutter
mean(is.na(bat.data$`Distance_to_Clutter_(m)`))

# Check covariates:Percent Clutter-------------------------------------------------------------------------------
xtabs(~Percent_Clutter,data=bat.data, exclude=NULL, na.action=na.pass)
#get the percentage of rows with no percent clutter
mean(is.na(bat.data$Percent_Clutter))

# Check covariates:Distance to water-------------------------------------------------------------------------------
mean(is.na(bat.data$DIST_WATER_M))
#update the water nearby column to make sure it is correct (if water within 100m then yes = 1)
bat.data$Water_Nearby[bat.data$DIST_WATER_M<=100] <- 1
bat.data$Water_Nearby[bat.data$DIST_WATER_M>100] <- 0

# Check covariate: Nightly Precipitation-------------------------------------------------------------------------------
mean(is.na(bat.data$Nightly_Precipitation))
  #over 70% of rows for nightly precipitation are empty. This means that we can't reliably use this metric

# Check covariate: Nightly Temperature-------------------------------------------------------------------------------
mean(is.na(bat.data$Nightly_Min_Temp))
mean(is.na(bat.data$Nightly_Max_Temp))
mean(is.na(bat.data$Nightly_Mean_Temp))

hist(bat.data$Nightly_Max_Temp)
hist(bat.data$Nightly_Min_Temp)
hist(bat.data$Nightly_Mean_Temp)

# Check covariate: Nightly RH-------------------------------------------------------------------------------
mean(is.na(bat.data$Nightly_Min_RH))
mean(is.na(bat.data$Nightly_Max_RH))
mean(is.na(bat.data$Nightly_Mean_RH))

hist(bat.data$Nightly_Min_RH)
hist(bat.data$Nightly_Max_RH)
hist(bat.data$Nightly_Mean_RH)

# Check covariate: Nightly Wind speed-------------------------------------------------------------------------------
mean(is.na(bat.data$Nightly_Mean_Windsp))
hist(bat.data$Nightly_Mean_Windsp)

# Check covariate: Moon Illumination and fraction-------------------------------------------------------------------------------
mean(is.na(bat.data$moon_illumination))
mean(is.na(bat.data$moon_fraction_above_horizon))
#create a new column that shows the product of illumination and time above horizon. This is our proxy for how much moon light is present
bat.data$moon <- bat.data$moon_illumination*bat.data$moon_fraction_above_horizon




#----------Now that you've checked to make sure all the data looks good set it up so it's ready for running models----------------------------------

bat.data.species <- c("ANPA", "COTO", "EPFU", "EUMA", "LABO", "LACI", "LANO", "MYCA", 
                      "MYCI", "MYEV", "MYLU", "MYTH", "MYVO", "MYSE", "MYYU") #, "PAHE", "TABR")

#separate transect from stationary data
bat.data <- bat.data[ bat.data$Quadrant != "Transects",]

#create summarized detection counts for all nights for each species. 

#This code has different columns for single, all couplets included, and couplets include only when there is a single
#list all species that you are interested in
species.id <- bat.data.species

##For Stationary
bat.data.long <- plyr::ldply(species.id, function(species, bat.data){
  bat.data2 <- bat.data
  bat.data$SpeciesGroup    <- species
  SpeciesGroupCountVars    <- bat.data.species[ grepl(species, bat.data.species)] 
  bat.data$SpeciesFullPool <- apply(bat.data[, SpeciesGroupCountVars],1,sum,na.rm=TRUE)
  bat.data$SpeciesSingleton<- as.vector(bat.data[, species, drop=TRUE])
  bat.data2[,SpeciesGroupCountVars] <- (diag( as.numeric(bat.data2[,species]>0))) %*% as.matrix(bat.data2[,SpeciesGroupCountVars] )
  bat.data$SpeciesPartPool <- apply(bat.data2[, SpeciesGroupCountVars],1,sum,na.rm=TRUE)
  bat.data
}, bat.data=bat.data)


#For transects
bat.transect.data.long <- plyr::ldply(species.id, function(species, bat.transect.data){
  bat.transect.data2 <- bat.transect.data
  bat.transect.data$SpeciesGroup    <- species
  SpeciesGroupCountVars    <- bat.data.species[ grepl(species, bat.data.species)]
  bat.transect.data$SpeciesFullPool <- apply(bat.transect.data[, SpeciesGroupCountVars],1,sum,na.rm=TRUE)
  bat.transect.data$SpeciesSingleton<- as.vector(bat.transect.data[, species, drop=TRUE])
  bat.transect.data2[,SpeciesGroupCountVars] <- (diag( as.numeric(bat.transect.data2[,species]>0))) %*% as.matrix(bat.transect.data2[,SpeciesGroupCountVars] )
  bat.transect.data$SpeciesPartPool <- apply(bat.transect.data2[, SpeciesGroupCountVars],1,sum,na.rm=TRUE)
  bat.transect.data
}, bat.transect.data=bat.transect.data)


#look at raw detentions summarized as a quick check

#---Get plot of the change in average detection per night at each region
# Get unique regions
unique_regions <- unique(bat.data.long$Region)
print(paste("Creating average per night plots for regions:", paste(unique_regions, collapse = ", ")))

# Create a list to store plots
regional_avg_plots <- list()

# Loop through each region to create individual plots
for(region in unique_regions) {
  
  # Filter data for current region
  regional_data <- bat.data.long[bat.data.long$Region == region, ]
  
  # Calculate average detections per night by species and year
  avg_detections <- regional_data %>%
    group_by(SpeciesGroup, Year) %>%
    summarise(
      total_detections = sum(SpeciesSingleton, na.rm = TRUE),
      total_nights = n_distinct(Night),  # Count unique nights
      avg_per_night = total_detections / total_nights,
      .groups = 'drop'
    )
  
  # Create faceted plot for this region
  regional_avg_plot <- ggplot(avg_detections, aes(x = Year, y = avg_per_night)) +
    geom_line(color = "steelblue", size = 1) +
    geom_point(color = "steelblue", size = 1.5) +
    geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "dashed", alpha = 0.7) +
    facet_wrap(~SpeciesGroup, scales = "free_y", ncol = 4) +
    labs(title = paste("Average Species Detections Per Night -", region, "Region"),
         subtitle = "Blue = actual data, Red dashed = linear trend line with confidence interval",
         x = "Year",
         y = "Average Detections per Night") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          strip.text = element_text(face = "bold"))
  
  # Store the plot in the list
  regional_avg_plots[[region]] <- regional_avg_plot
  
  # Print the plot
  print(regional_avg_plot)
}

## Now use the bat.data.long and bat.transect.data.long to run the models in 02_ModelSelection