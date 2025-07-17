#! /usr/bin/Rscript

library(dplyr)
library(tidyverse)

# Calculate weather data and correct full metadata file (has slightly incorrect weather info)
# NOTE: NEH2 lists 2 different planting dates a month apart, and 

# Load data
full_metadata = read.delim("0_data_files/G2f_2019_sub_corrected_metadata.tsv")
weather_2019 <- read.csv("0_data_files/G2F_2019_Weather_ARK_cleaned.csv")
sample_metadata <- read.delim("0_data_files/sample_metadata.tsv") %>%
  filter(!is.na(location) ) # Remove blanks, which don't have location data
plot_data = read.csv("0_data_files/g2f_2019_phenotypic_clean_data.csv")

# Format the weather data timestamp is in
timestamp_format="%m/%d/%y %H:%M %p"

# Define columns to do operations on
mean_cols = c("Temperature..C.", "Dew.Point..C.", "Relative.Humidity....", 
              "Wind.Speed..m.s.", "Wind.Direction..degrees.", "Wind.Gust..m.s.",
              "Soil.Temperature..C.", "Soil.Moisture...VWC." )
sum_cols = c("Solar.Radiation..W.m2.", "Rainfall..mm." )

# Create key of time range info
planting_dates = plot_data %>%
  filter(Year==2019) %>%
  select(Field.Location, Date.Plot.Planted..MM.DD.YY.) %>%
  rename(location=Field.Location, planting_date = Date.Plot.Planted..MM.DD.YY.) %>%
  filter(planting_date != "") %>% # Remove any empty ones
  mutate(planting_date = sub(planting_date, pattern="/19$", repl="/2019")) %>% # Standardize
  unique() %>%
  mutate(start = as.Date(planting_date, format="%m/%d/%Y"))

# Reduce to earliest planting date per location
counts = table(planting_dates$location)
if(any(counts != 1)){
  warning("Some locations don't have a single planting value; taking earliest for:",
          names(counts)[counts!=1])
  planting_dates = planting_dates %>%
    group_by(location) %>%
    arrange(start) %>%
    slice_head(n=1)
}

# Make a time range data frame with start and stop for each location
time_range <- sample_metadata[,c("location","Month","day","year")] %>%
  unique() %>% 
  mutate(end = paste(year, Month, day, sep="-")) %>%
  mutate(end = as.Date(end, format="%Y-%m-%d")) %>%
  select(location, end) %>%
  full_join(planting_dates, by="location") %>%
  filter(!is.na(end)) %>%
  select(location, start, end)




# Create data frame to add results into
weather_result <- list()
weather_subsets = list()

# Loop over locations and calculate weather
for(i in 1:nrow(time_range)){
  sample_location = time_range$location[i]
  
  # Pull out weather for just this location and time window
  myweather = subset(weather_2019,Field.Location == sample_location) %>%
    mutate(date=as.Date(DATE_Local, format=timestamp_format)) %>%
    filter(date >= time_range$start[i] & date <= time_range$end[i])

  mean_result <- myweather %>%
    select(all_of(mean_cols)) %>%
    colMeans(na.rm=TRUE) %>%
    as.data.frame() %>%
    t() %>%
    as.data.frame() %>%
    mutate(location = sample_location)
  Sum_result <- myweather %>%
    select(all_of(sum_cols)) %>%
    colSums(na.rm = TRUE)%>%
    as.data.frame() %>%
    t() %>%
    as.data.frame() %>%
    mutate(location = sample_location)
  
  weather_result[[i]] <- dplyr::left_join(mean_result,Sum_result,by="location")
  weather_subsets[[i]] = myweather
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


#####################
# Diagnostic plots
#####################

# Basic plots
weatherlong = weather %>%
  pivot_longer(cols=-location, names_to="parameter", values_to="value") %>%
  mutate(summary_type = ifelse(parameter %in% mean_cols, yes="Mean", no="None")) %>%
  mutate(summary_type = ifelse(parameter %in% sum_cols, yes="Sum", no=summary_type)) 
summaryplot = ggplot(weatherlong) +
  aes(x=location, y=value, fill=summary_type) +
  geom_col() + 
  facet_wrap(~ parameter, scales="free") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5)) 
ggsave(summaryplot, file="1_parsed_files/1_weather.summaries.png", width=12, height=8)


# Plot of timepoints to make sure average make sense
weatherfull = bind_rows(weather_subsets) %>%
  rename(location=Field.Location) %>%
  mutate(date=as.POSIXct(DATE_Local, format=timestamp_format)) %>%
  select(location, date, all_of(c(mean_cols, sum_cols))) %>%
  pivot_longer(cols=-c(location, date), names_to="parameter", values_to="value") %>%
  filter(is.finite(value))
bigplot = ggplot(weatherfull) + 
  aes(x=date, y=value, color=parameter) +
  geom_point(alpha=0.25) + 
  facet_grid(parameter ~ location, scales="free_y") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5)) +
  geom_hline(data=weatherlong %>% filter(parameter %in% mean_cols), mapping=aes(yintercept=value))
ggsave(bigplot, file="1_parsed_files/1_weather.fullplot.png", width=20, height=20)
  
  