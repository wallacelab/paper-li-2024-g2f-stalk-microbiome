#! /usr/bin/Rscript

# Plot G2F locations. Also used to develop a standard color key for locations

# TODO: more were sampled. Need to differentiate between "not sampled" and "filtered out"

library(tidyverse) 
library(phyloseq)
library(ggplot2)
library(ggmap)
library(pals)

# Load data
mapfile="1_parsed_files/1b_map_data.rds"  # File to store map data so only have to call once
metadata=read.csv("0_data_files/g2f_2019_field_metadata.csv")
raw_locs = read.delim("0_data_files/G2f_2019_sub_corrected_metadata.tsv")$location %>% 
  unique() %>% na.omit() %>% as.character()  # Pull out all original locations
filtered = readRDS("1_parsed_files/1a_asv_table_no_taxa_from_blanks.rarefy1500.phyloseq.rds")
filtered_locs = sample_data(filtered)$location %>% unique()


#########
# Cleaning and formatting
#########


# Streamline data table
plotdata=metadata %>%
  rename(location=Experiment_Code,
         latitude=Latitude_of_Field_Corner_.1..lower.left.,
         longitude=Longitude_of_Field_Corner_.1..lower.left.) %>%
  select(location, latitude, longitude)

# Fix error; NEH2 lat/long not recorded for anything, so copy from NEH1 since are both in Lincoln, NE
plotdata$latitude[plotdata$location=="NEH2"] = plotdata$latitude[plotdata$location=="NEH1"]
plotdata$longitude[plotdata$location=="NEH2"] = plotdata$longitude[plotdata$location=="NEH1"]

# Determine if sampled & if passed filters
plotdata$sampled = plotdata$location %in% raw_locs
plotdata$pass_filters = plotdata$location %in% filtered_locs

# Group into regions
plotdata = plotdata %>% mutate(region = case_when(
  location %in% c("IAH1", "IAH2", "IAH3", "IAH4", "ILH1", "MNH1", "MOH1", "NEH1", "NEH2") ~ "West",
  location %in% c("INH1", "MIH1", "OHH1", "WIH1", "WIH2") ~ "Central",
  location %in% c("GAH1", "GAH2", "NCH1", "SCH1", "TXH1", "TXH2", "TXH3", "TXH4") ~ "South",
  location %in% c("DEH1", "NYH1", "NYH2", "NYH3") ~ "Northeast",
  location %in% c("COH1", "GEH1", "ONH2") ~ "Other"
))

# Define location color key
ncolors = sum(plotdata$pass_filters)
plotdata$location_color=NA
plotdata$location_color[!plotdata$pass_filters]="black" # Black for ones sampled but filtered out
plotdata$location_color[!plotdata$sampled]="gray50" # Gray for unsampled locations
plotdata$location_color[plotdata$pass_filters]= pals::kelly()[1:ncolors +1]  # +1 offset to avoid light gray & make higher contrast among sites in same region
## Make a variable for use later
lockey = plotdata$location_color
names(lockey) = plotdata$location

# Define region color key; colorblind-friendly palette from https://davidmathlogic.com/colorblind/
regionkey = c("Central" = "#cd9a00",  # Slightly darker than original for easier visibility
              "West" = "#D81B60", 
              "Northeast" = "#1E88E5", 
              "South" = "#004D40", 
              "Other" = "gray50")
plotdata$region_color = regionkey[plotdata$region]

# write out location key for future use
write.csv(plotdata, file="1_parsed_files/1b_location_key.csv", row.names=FALSE)

#########
# Color test of within-region differentiation of colors 
#########

subcolor = subset(plotdata, plotdata$pass_filters)
colorkey = subcolor$location_color
names(colorkey)=subcolor$location
colortest = ggplot(subcolor) +
  aes(x=location, y=pass_filters, fill=location) +
  geom_col() +
  scale_fill_manual(values=colorkey) +
  facet_wrap(~region, scale="free_x") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5),
        legend.position="none")
ggsave(colortest, file="1_parsed_files/1b_location_colortest.png")



###########
# Plot locations
###########

# Define bounding box
lower_left = c(-130, 20)
upper_right=c(-50, 50)
window = c(lower_left, upper_right )

# Get map
if(file.exists(mapfile)){
  mymap=readRDS(mapfile)
}else{
  mymap = get_map(location=window, zoom=5, maptype="stamen_terrain_background")
  saveRDS(mymap, mapfile)
}
states=map_data("state")

# Helper function to make base map and standard stuff
basemap = function(m, xlim=c(-110, -70), ylim=c(25, 47)){
  ggmap(m) +
    geom_polygon(mapping=aes(x=long, y=lat, group=group), data=states, fill="#00000000", color="black") +
    coord_map(projection="mercator", xlim=xlim, ylim=ylim) +
    labs(x="Longitude", y="Latitude")
}

# Standard parameters
final = plotdata %>% filter(pass_filters==TRUE)
pointsize=3
pointshape=21
pointstroke=0.5

# Map colored by region (all)
regionplot = basemap(mymap) +
  geom_point(mapping=aes(x=longitude, y=latitude, fill=region), 
             data=plotdata, size=pointsize, shape=pointshape, stroke=pointstroke) +
  scale_fill_manual(values=regionkey)
ggsave(regionplot, file="1_parsed_files/1b_g2f_locations.by_region.png", width=5, height=4)  

# Map colored by region (Sampled)
regionplot = basemap(mymap, xlim=c(-100, -70), ylim=c(29, 47)) +
  geom_point(mapping=aes(x=longitude, y=latitude, fill=region), 
             data=final, size=pointsize, shape=pointshape, stroke=pointstroke) +
  scale_fill_manual(values=regionkey)
ggsave(regionplot, file="1_parsed_files/1b_g2f_locations.by_region.sampled.png", width=5, height=4)  


# Map colored by locations (all)
locationplot = basemap(mymap) +
  geom_point(mapping=aes(x=longitude, y=latitude, fill=location), 
             data=plotdata, size=pointsize, shape=pointshape, stroke=pointstroke) +
  scale_fill_manual(values=lockey)
ggsave(locationplot, file="1_parsed_files/1b_g2f_locations.by_location.png", width=5, height=4)  

# Subset to just the sampled locations
sampleplot = basemap(mymap, xlim=c(-100, -70), ylim=c(29, 47)) +
  geom_point(mapping=aes(x=longitude, y=latitude, fill=location), 
             data=final, size=pointsize, shape=pointshape, stroke=pointstroke) +
  scale_fill_manual(values=lockey)
ggsave(sampleplot, file="1_parsed_files/1b_g2f_locations.by_location.sampled.png", width=5, height=4)  


################
# Publication Graphic
################

# Adjustments
final$abbreviation = sub(final$location, pattern="H.$", repl="")
regions = data.frame(region = names(regionkey), 
                     region_color=regionkey,
                     latitude =c( 38.5,    40.5,    40.7,   34.0,   0),  # 0 is for "other"
                     longitude=c(-85,  -93.6,      -76,  -79.5,   0))
  

# Final graphic for publication
pubplot = basemap(mymap, xlim=c(-100, -70), ylim=c(29, 47)) +
  geom_point(mapping=aes(x=longitude, y=latitude, fill=region), 
             data=final, size=3.5, shape=21, stroke=0.5) +
  scale_fill_manual(values=regionkey) +
  geom_text(mapping=aes(x=longitude, y=latitude, label=abbreviation), 
            data=final, color="white", size=1.5, fontface="bold") +
  annotate(geom="label", x=regions$longitude, y=regions$latitude, fill=regions$region_color,
           label=regions$region, color = "white", fontface = "bold", size=3, alpha=0.8) +
  theme(legend.position="none")
ggsave(pubplot, file="1_parsed_files/1b_g2f_locations.publication.png", width=4, height=4)  

