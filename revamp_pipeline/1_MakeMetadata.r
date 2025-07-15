#! /usr/bin/Rscript

library(dplyr)
library(tidyverse)

# Calculate weather data and correct full metadata file (has slightly incorrect weather info)

# Load data
full_metadata = read.delim("0_data_files/G2f_2019_sub_corrected_metadata.tsv")
weather_2019 <- read.csv("0_data_files/G2F_2019_Weather_ARK_cleaned.csv")
g2f_2019_metadata <- read.delim("0_data_files/sample_metadata.tsv") %>%
  filter(!is.na(location) ) # Remove blanks, which don't have location data
time_range <- g2f_2019_metadata[,c("location","Month","day","year")] %>%
  unique()

# Define columns to do operations on
mean_cols = c("Temperature..C.", "Dew.Point..C.", "Relative.Humidity....", 
              "Wind.Speed..m.s.", "Wind.Direction..degrees.", "Wind.Gust..m.s.",
              "Soil.Temperature..C.", "Soil.Moisture...VWC." )
sum_cols = c("Solar.Radiation..W.m2.", "Rainfall..mm." )


# Create data frame to add results into
weather_result <- list()

# Loop over locations  
for(i in 1:nrow(time_range)){
  sample_location = time_range$location[i]
  temp_location_envriomental_list = subset(weather_2019,Field.Location == sample_location)
  
  #print(temp_location_envriomental_list)
  
  mean_result <- temp_location_envriomental_list %>%
    select(all_of(mean_cols)) %>%
    colMeans(na.rm=TRUE) %>%
    as.data.frame() %>%
    t() %>%
    as.data.frame() %>%
    mutate(location = sample_location)
  Sum_result <- temp_location_envriomental_list %>%
    select(all_of(sum_cols)) %>%
    colSums(na.rm = TRUE)%>%
    as.data.frame() %>%
    t() %>%
    as.data.frame() %>%
    mutate(location = sample_location)
  
  weather_result[[i]] <- dplyr::left_join(mean_result,Sum_result,by="location")
#  weather_result <-  dplyr::bind_rows(weather_result,weather_result_per_location)
}

G2F_weather_result <- weather_result %>% 
  bind_rows() %>%
  column_to_rownames("location")
write.table(G2F_weather_result,"1_parsed_files/1_G2F_2019_weather_average_and_sum.tsv",row.names = TRUE,sep = "\t")


##############
# Correct existing full metadata file
##############

cols_to_replace = names(G2F_weather_result)
weather = G2F_weather_result %>%
  rownames_to_column("location")
fixed_metadata = full_metadata %>%
  select(-all_of(cols_to_replace)) %>%
  left_join(weather, by="location")
write.table(fixed_metadata, file="1_parsed_files/1_G2f_2019_sub_corrected_metadata.fix_weather.tsv",
            sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)


