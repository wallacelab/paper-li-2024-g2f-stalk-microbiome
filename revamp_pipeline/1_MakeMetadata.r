library(dplyr)

calculate_weather_information_all_time_range <- function(weather_file,time_range_file){
  
  time_range_file_row <- nrow(time_range_file)
  time_range_file_col <- nrow(time_range_file)
  
  weather_result <- data.frame(matrix(NA, nrow=0, ncol=11))
  
  for(i in 1:time_range_file_row){
    #get location, plant time 
    sample_location = time_range_file[i,1]
    #sample_Month = time_range_file[i,2]
    #sample_Day = time_range_file[i,3]
    
    #print(sample_location)
    #print(sample_Month)
    #print(sample_Day)
    #select two time interval  One is from planting to harvest month  the other is the days in harvest month 
    #temp_location_encriomental_list_1 = subset(weather_file,Field.Location == sample_location & Month < sample_Month)
    #temp_location_encriomental_list_2 = subset(weather_file,Field.Location == sample_location & Month == sample_Month & Day <= sample_Day)
    #temp_location_encriomental_list = dplyr::bind_rows(temp_location_encriomental_list_1,temp_location_encriomental_list_2)
    temp_location_encriomental_list = subset(weather_file,Field.Location == sample_location)
    
    #print(temp_location_encriomental_list)
    
    mean_result <- as.data.frame(colMeans(temp_location_encriomental_list[,c(11:13,16:20)],na.rm = TRUE))
    mean_result <- as.data.frame(t(mean_result))
    mean_result$location <- temp_location_encriomental_list$Field.Location[1]
    Sum_result <- as.data.frame(colSums(temp_location_encriomental_list[,c(14:15)],na.rm = TRUE))
    Sum_result <- as.data.frame(t(Sum_result))
    Sum_result$location <- temp_location_encriomental_list$Field.Location[1]
    print(Sum_result)
    print(mean_result)
    
    weather_result_per_location <- dplyr::left_join(mean_result,Sum_result,by="location")
    #colnames(weather_result) <- colnames(weather_result_per_location)
    weather_result <-  dplyr::bind_rows(weather_result,weather_result_per_location)
    
  }
  
  return(weather_result)
}




##############################

#calculate the average of weather data 

weather_2019 <- read.csv("/home/hl46161/G2F_data_dada2/g2f_2019_enviromental_data/Weather/2019_Weather_ARK_cleaned.csv",sep=",")
g2f_2019_metadata <- read.csv("/home/hl46161/G2F_data_dada2/2018_G2F_metadata.tsv",sep="\t",stringsAsFactors = FALSE)
g2f_2019_metadata <-  g2f_2019_metadata[-c(1:12),]
time_range <- g2f_2019_metadata[,c("location","Month","day","year")]
time_range$year <- 2019
time_range <- unique(time_range)
time_range

class(weather_2019$Month[1]) 
class(time_range$Month[1]) 

G2F_weather_result <- calculate_weather_information(weather_2019,time_range)
G2F_weather_result <- G2F_weather_result[,12:ncol(G2F_weather_result)]
rownames(G2F_weather_result) <- G2F_weather_result$location
G2F_weather_result$location <- NULL
write.table(G2F_weather_result,"/home/hl46161/new_G2F_dada2/G2F_2019_weather_aervage_and_sum.tsv",row.names = TRUE,sep = "\t")
