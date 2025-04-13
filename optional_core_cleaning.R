#cleaning and unify format of 10 years data
#library
library(dplyr)
library(ggplot2)
library(ggmap)
library(tidyr)
library(ggrepel)




setwd("~/PERMA/R")
file_list <- c(2013:2023) %>%
  paste0("~/PERMA/",.,".csv")

data_list <- lapply(file_list, read.csv)
names(data_list) <- c(2013:2023)

#examine shared columns
column_names_list <- lapply(data_list, names)
shared_column_names <- Reduce(intersect, column_names_list)

#keep shared columns
filtered_data_list <- lapply(data_list, function(df) df[, shared_column_names, drop = FALSE])

for (i in 1:length(filtered_data_list)){
  filtered_data_list[[i]]$year <- names(filtered_data_list[i])
}

all.data <- do.call(rbind, filtered_data_list)
select_names <- shared_column_names[c(1,2,4,5,9,13,14,16,17,18,21,25)]
data_clean <- all.data[,select_names]
colnames(data_clean) <- c("sample", "core", "Tsoil", "pH", "GPS", "site", "habitat", "WTD", "ALD", "depth_interval", "depth", "year")
#ph measn"pH_porewater" 
# 
# #reclassify depth group
# lines <- strsplit(data_clean$depth_interval, "-")
# intervals <- unlist(strsplit(lines[2], "\\s+"))
# values <- as.numeric(unlist(strsplit(lines[4], "\\s+")))
# 
# classify_interval <- function(interval) {
#   if (grepl("-", interval)) {
#     parts <- as.numeric(unlist(strsplit(interval, "-")))
#     start <- parts[1]
#     end <- parts[2]
#   } else {
#     start <- as.numeric(interval)
#     end <- start
#   }
#   
#   group_start <- floor(start / 10) * 10
#   group_end <- group_start + 9
#   group <- paste0(group_start, "-", group_end)
#   
#   return(group)
# }
# data_clean$group <- sapply(data_clean$depth_interval, classify_interval)
# 
# write.csv(data_clean, "core_metadata.csv", row.names = T)



#correct format of GPS
convert <- function(dms_string) {
  # Check if the string starts with 'N' or 'E'
  if (startsWith(dms_string, "N ")) {
    # Split the string by space and comma
    parts <- strsplit(dms_string, "[ ,]+")[[1]]
    
    # Extract the degrees, minutes, and direction for latitude and longitude
    lat_deg <- as.numeric(parts[2])
    lat_min <- as.numeric(parts[3])
    long_deg <- as.numeric(parts[5])
    long_min <- as.numeric(parts[6])
    
    # Convert to decimal degrees
    lat_decimal <- lat_deg + (lat_min / 60)
    long_decimal <- long_deg + (long_min / 60)
    
    # Reattach the direction
    lat_decimal <- paste0(format(lat_decimal, digits = 10), parts[1])
    long_decimal <- paste0(format(long_decimal, digits = 10), parts[4])
    
    # Combine the latitude and longitude
    return(paste(lat_decimal, long_decimal, sep = ", "))
  } else {
    # If the string does not start with 'N' or 'E', return it as is
    return(dms_string)
  }
}

data_clean$GPS <- gsub("\\. ", ".", data_clean$GPS)

data_clean <- data_clean %>%
  mutate(GPS = sapply(GPS, convert)) 

#1. check coordinates 
GPS_check <- data_clean[,c("sample", "GPS", "habitat", "site", "core" )]
GPS_check$sample <- paste0(data_clean$year, "_", GPS_check$sample, "_", data_clean$core)

GPS_check_unique <- GPS_check %>%
  distinct(sample, GPS, .keep_all = TRUE) %>%
  separate(GPS, into = c("Latitude", "Longitude"), sep = ", ") %>%
  mutate(Latitude = as.numeric(sub("N|S", "", Latitude)),
         Longitude = as.numeric(sub("E|W", "", Longitude)))

GPS_check_unique<- GPS_check_unique[complete.cases(GPS_check_unique$Latitude, GPS_check_unique$Longitude), ]
GPS_check_unique$year <- substring(GPS_check_unique$sample, 1, 4)

color <- c("Palsa" = "#B8912F", "Bog" = "#5E813F", "Fen" = "#4273B0", "Collapsed Palsa" = "orange" , "Poor Fen" = "blue")

GPS <- GPS_check_unique[complete.cases(GPS_check_unique$Latitude, GPS_check_unique$Longitude), ]

ggplot(data = GPS, aes(x = Longitude, y = Latitude)) +
  geom_point(aes(color = habitat)) +
  labs(x = "Longitude", y = "Latitude") +
  scale_color_manual(values = color) +
  ggtitle("GPS Coordinates") +
  theme_minimal()

filtered_GPS <- GPS[grep("autochamber|Inc|Avni", GPS$site, ignore.case = TRUE), ]
filtered_GPS <- filtered_GPS[grep("\\bautochamber\\b|\\btransect C\\b|\\bInc-Sphagnum\\b|\\bInc-Eriophorum\\b", filtered_GPS$site, ignore.case = TRUE), ]

filtered_GPS_23 <- filtered_GPS[filtered_GPS$year %in% c(2022, 2023), ]

centroids <- filtered_GPS_23 %>%
  group_by(site) %>%
  summarize(Longitude = mean(Longitude), Latitude = mean(Latitude), habitat = first(habitat))  # Use the first habitat for color

location.23 <- ggplot(data = filtered_GPS_23, aes(x = Longitude, y = Latitude)) +
  geom_point(aes(color = habitat)) +
  labs(x = "Longitude", y = "Latitude") +
  scale_color_manual(values = color) +
  ggtitle("GPS Coordinates") +
  theme_minimal()+
  geom_text(data = centroids, aes(label = site),  size = 4, check_overlap = TRUE)  # Adjust nudge_y as needed

ggsave("GPS_23.pdf",location.23,  width = 5, height = 5)


center_lon <- (max(filtered_GPS_23$Longitude) - min(filtered_GPS_23$Longitude))/2 +  min(filtered_GPS_23$Longitude)
center_lat <- (max(filtered_GPS_23$Latitude) - min(filtered_GPS_23$Latitude))/2 +  min(filtered_GPS_23$Latitude)

register_google(key = #use your google API, write = TRUE)

map <- get_googlemap(center = c(lon = center_lon, lat = center_lat), zoom = 16, maptype = "satellite")

map.all<- ggmap(map) +
  geom_point(data = GPS, aes(x = Longitude, y = Latitude, color = habitat))+
  labs(x = "Longitude", y = "Latitude") +
  scale_color_manual(values = color) +
  ggtitle("GPS Coordinates all sites")

centroids <- centroids %>%
  rename(lon = Longitude, lat = Latitude)

map.filtered <- ggmap(map) +
  geom_point(data = filtered_GPS, aes(x = Longitude, y = Latitude, color = habitat)) +
  labs(x = "Longitude", y = "Latitude") +
  scale_color_manual(values = color) +
  ggtitle("GPS Coordinates Avni sites") +
  theme_minimal() 

ggsave("map.filtered.pdf", map.filtered, height = 7.5, width = 10, units = "in")

ggsave("map_all.png", map.all, height = 7.5, width = 10, units = "in", dpi = 300 )

GPS_avni <- subset(GPS, grepl("Avni", site))
ggmap(map) +
  geom_point(data = GPS_avni, aes(x = Longitude, y = Latitude, color = habitat))+
  labs(x = "Longitude", y = "Latitude") +
  scale_color_manual(values = color) +
  ggtitle("GPS Coordinates all sites")


GPS_ac <- subset(GPS, grepl("Autochamber", site))
GPS_ac$label <- paste0(substr(GPS_ac$habitat, 1, 1), GPS_ac$core)


center_lon <- (max(GPS_ac$Longitude) - 19.04654)/2 +  19.04654
center_lat <- (68.35357 - min(GPS_ac$Latitude))/2 +  min(GPS_ac$Latitude)
map <- get_googlemap(center = c(lon = center_lon, lat = center_lat), zoom = 19, maptype = "satellite")

years <- unique(GPS_ac$year)

for (year in years) {
  # Subset data for the current year
  subset_data <- GPS_ac[GPS_ac$year == year, ]
  
  # Create the plot for the current year
  p <- ggmap(map) +
    geom_point(data = subset_data, aes(x = Longitude, y = Latitude, color = habitat), size = 3) +
    geom_text(data = subset_data, aes(x = Longitude, y = Latitude, label = label), vjust = -1, hjust = 0.5, color = "orange") +             
    labs(x = NULL, y = NULL) + # Remove axis labels
    scale_color_manual(values = color) +
    ggtitle(paste("GPS Coordinates on Google Map - Year", year)) + 
    theme(legend.position = "none")+
    theme(axis.line = element_blank(),         # Remove axis lines
          axis.text = element_blank(),         # Remove axis text
          axis.ticks = element_blank(),        # Remove axis ticks
          axis.title = element_blank(),        # Remove axis titles
          plot.title = element_text(hjust = 0.5), # Center the plot title
          panel.grid = element_blank(),        # Remove grid lines
          panel.border = element_blank())      # Remove panel border
  
  ggsave(filename = paste0("GPS_Coordinates_Year_", year, ".png"), plot = p, width = 7.5, height = 7.5, units = "in" , dpi = 300)
}

#TEMPERATURES -----
T.ac <- subset(data_clean, grepl("Autochamber", site))
T.ac$Tsoil <-as.numeric(T.ac$Tsoil)
T.ac <- subset(T.ac, T.ac$group %in% c("0-9", "10-19", "20-29", "30-39", "40-49"))

Tsoil.base<-ggplot(data = T.ac, aes(x = depth, y = Tsoil, color = habitat)) +
  geom_point(alpha = 0.4) +
  geom_smooth(aes(fill = habitat), method = "loess", formula = y ~ x, se = TRUE) + # Adjusted formula to use 'y ~ x'
  labs(x = "Depth", y = "Soil Temperature (Tsoil)") +
  ggtitle("Scatter Plot of Depth vs Soil Temperature") +
  scale_color_manual(values = color) +
  scale_fill_manual(values = color) +
  scale_x_reverse() + # Reverse the x-axis to properly reflect depth
  coord_flip() +
  theme_minimal()

ggplot(data = T.ac, aes(x = year, y = Tsoil, color = habitat)) +
  geom_point(alpha = 0.4) +
  geom_smooth(aes(fill = habitat), method = "loess", formula = y ~ x, se = TRUE) + # Adjusted formula to use 'y ~ x'
  labs(x = "Year", y = "Soil Temperature (Tsoil)") +
  ggtitle("Soil Temperature yearly change") +
  scale_color_manual(values = color) +
  scale_fill_manual(values = color) +
  theme_minimal()+
  facet_wrap(~group)


years <- unique(T.ac$year)

for (year in years) {
  # Subset data for the current year
  subset_data <- T.ac[T.ac$year == year, ]
  
  subset_data <- T.ac[T.ac$year == 2022, ]
  # Create the plot for the current year
  p <- ggplot(data = subset_data, aes(x = depth, y = Tsoil, color = habitat)) +
    geom_point(alpha = 0.4) +
    geom_smooth(aes(fill = habitat),method = "loess", formula = y ~ x, se = TRUE) + # Adjusted formula to use 'y ~ x'
    labs(x = "Depth", y = "Soil Temperature (Tsoil)") +
    ggtitle(paste0(year,"_Depth vs Soil Temperature")) +
    scale_color_manual(values = color) +
    scale_fill_manual(values = color) +
    scale_y_continuous(limits = c(0, 20)) +  # Set x-axis limits
    scale_x_continuous(limits = c(0, 80))+  # Set y-axis limits
    scale_x_reverse() + # Reverse the x-axis to properly reflect depth
    coord_flip()+
    theme_bw()

  ggsave(filename = paste0("Tsoil_Year_", year, ".png"), plot = p, width = 5, height = 4, units = "in" , dpi = 300)
}


T.air <- read.csv("~/PERMA/WeatherHawk_Stordalen_MetData_2013-2019_10_min.csv")
# Calculate daily average temperature
daily_avg <- T.air %>%
  group_by(year, day) %>%
  summarize(Avg_TairC = mean(Tair_C, na.rm = TRUE))

# Calculate monthly average temperature
monthly_avg <- T.air %>%
  group_by(year, month) %>%
  summarize(Avg_TairC = mean(Tair_C, na.rm = TRUE))

# Convert year and day of year to a Date format
daily_avg$Date <- as.Date(ISOdate(daily_avg$year, 1, 1) + daily_avg$day - 1)
daily_avg$Date <- as.Date(paste(daily_avg$year, daily_avg$day), format = "%Y %j")
monthly_avg$Date <- as.Date(ISOdate(monthly_avg$year, monthly_avg$month, 15))

Tchange <- ggplot(daily_avg, aes(x = Date, y = Avg_TairC)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +  # Add linear regression line for all data
  geom_smooth(data = filter(monthly_avg, month == 7), aes(x = Date, y = Avg_TairC, group = 1), method = "lm", se = FALSE, color = "red") +  # Add linear regression line for July data
  geom_point(data = filter(monthly_avg, month == 7), aes(x = Date, y = Avg_TairC), color = "red", size = 3) +  # Highlight July points
  labs(x = "Date", y = "Average Temperature (TairC)", title = "Average Daily Temperature Change Over Time") +
  theme_minimal()

ggsave("airT.png", Tchange, width = 14, height = 3, unit = "in", dpi = 300)

# look for collapsed palsa
T.cp <- subset(data_clean, grepl("Collapsed Palsa", site) | grepl("Collapsed Palsa", habitat))
T.cp$Tsoil <-as.numeric(T.cp$Tsoil)
T.cp <- T.cp %>%
  mutate(CoreID = paste(sample, core, year, sep = "_"))

labels <- T.cp %>%
  filter(!is.na(Tsoil)) %>%
  group_by(CoreID) %>%
  slice_max(depth, n = 1, with_ties = FALSE) %>%
  ungroup()

cp <- Tsoil.base + 
  geom_point(data = T.cp, aes(x = depth, y = Tsoil),color = "red", size = 1.5)+
  geom_line(data = T.cp, aes(x = depth, y = Tsoil, group = CoreID), color = "red")+ 
  geom_text_repel(data = labels, aes(x = depth, y = Tsoil, label = CoreID), 
                  size = 3, color = "blue", nudge_y = 0.5, max.overlaps = 10)+
  labs(title = "Collapsed Palsa")


ggsave("Tsoil_all.png", Tsoil.base, width =5, height = 4, dpi = 300)
ggsave("Tsoil_collapsed palsa.png", cp, width =5, height = 4, dpi = 300)


T.pf <- subset(data_clean, grepl("Poor Fen", habitat))
T.pf$Tsoil <-as.numeric(T.pf$Tsoil)
T.pf <- T.pf %>%
  mutate(CoreID = paste(sample, core, year, sep = "_"))

labels <- T.pf %>%
  filter(!is.na(Tsoil)) %>%
  group_by(CoreID) %>%
  slice_max(depth, n = 1, with_ties = FALSE) %>%
  ungroup()

pf <- Tsoil.base + 
  geom_point(data = T.pf, aes(x = depth, y = Tsoil),color = "red", size = 1.5)+
  geom_line(data = T.pf, aes(x = depth, y = Tsoil, group = CoreID), color = "red")+ 
  geom_text_repel(data = labels, aes(x = depth, y = Tsoil, label = CoreID), 
                  size = 3, color = "blue", nudge_y = 0.5, max.overlaps = Inf)+
  labs(title = "Poor fen")

ggsave("Tsoil_poorfen.png", pf, width =5, height = 4, dpi = 300)



T.avni <- data_clean[grepl("Inc|Avni|Anvi", data_clean$site), ]
T.avni$Tsoil <-as.numeric(T.avni$Tsoil)
T.avni <- T.avni %>%
  mutate(CoreID = paste(sample, core, year, sep = "_"))

T.avni.p <- T.avni[grepl("Palsa", T.avni$habitat), ]
T.avni.f <- T.avni[grepl("Fen", T.avni$habitat), ]
T.avni.b <- T.avni[grepl("Bog", T.avni$habitat), ]


labels <- T.avni.p %>%
  filter(!is.na(Tsoil)) %>%
  group_by(CoreID) %>%
  slice_max(depth, n = 1, with_ties = FALSE) %>%
  ungroup()

avni.p <- Tsoil.base + 
  geom_point(data = T.avni.p, aes(x = depth, y = Tsoil),color = "red", size = 1.5)+
  geom_line(data = T.avni.p, aes(x = depth, y = Tsoil, group = CoreID), color = "red")+
  geom_text_repel(data = labels, aes(x = depth, y = Tsoil, label = CoreID), size = 3, color = "blue", nudge_y = 0.5, max.overlaps = 10)


labels <- T.avni.b %>%
  filter(!is.na(Tsoil)) %>%
  group_by(CoreID) %>%
  slice_max(depth, n = 1, with_ties = FALSE) %>%
  ungroup()

avni.b <- Tsoil.base + 
  geom_point(data = T.avni.b, aes(x = depth, y = Tsoil),color = "red", size = 1.5)+
  geom_line(data = T.avni.b, aes(x = depth, y = Tsoil, group = CoreID), color = "red")+
  geom_text_repel(data = labels, aes(x = depth, y = Tsoil, label = CoreID), size = 3, color = "blue", nudge_y = 0.5, max.overlaps = Inf)

labels <- T.avni.f %>%
  filter(!is.na(Tsoil)) %>%
  group_by(CoreID) %>%
  slice_max(depth, n = 1, with_ties = FALSE) %>%
  ungroup()

avni.f <- Tsoil.base + 
  geom_point(data = T.avni.f, aes(x = depth, y = Tsoil),color = "red", size = 1.5)+
  geom_line(data = T.avni.f, aes(x = depth, y = Tsoil, group = CoreID), color = "red")+
  geom_text_repel(data = labels, aes(x = depth, y = Tsoil, label = CoreID), size = 3, color = "blue", nudge_y = 0.5, max.overlaps = Inf)

ggsave("Tsoil_avni.b.png", avni.b + labs(title = "avni bog"), width =5, height = 4, dpi = 300)
ggsave("Tsoil_avni.p.png", avni.p + labs(title = "avni palsa"), width =5, height = 4, dpi = 300)
ggsave("Tsoil_avni.f.png", avni.f + labs(title = "avni fen"), width =5, height = 4, dpi = 300)


#testing hydro
ALD_df <- data_clean %>%
  dplyr::select(ALD) %>%
  mutate(ALD = as.numeric(ALD))
ALD_df$habitat <- data_clean$habitat

ggplot(ALD_df, aes(x = ALD, color = habitat)) +
  geom_density(alpha = 0.7, size = 1) +
  labs(title = "Density Plot of ALD for Each Habitat", x = "ALD", y = "Density") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = color) +
  scale_color_manual(values = color)


WTD_df <- data_clean %>%
  dplyr::select(WTD) %>%
  mutate(WTD = as.numeric(WTD))
WTD_df$habitat <- data_clean$habitat


ggplot(WTD_df, aes(x = WTD, color = habitat)) +
  geom_density(alpha = 0.7, size = 1) +
  labs(title = "Water table depth for Each Habitat", x = "WTD", y = "Density") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = color) +
  scale_color_manual(values = color)


#track site change by record.
reduced_df <- GPS %>%
  dplyr::select(site, year, habitat) %>%
  distinct(site, year, habitat)

changes_df <- reduced_df %>%
  group_by(site) %>%
  arrange(year) %>%
  mutate(habitat_change = habitat != lag(habitat, default = first(habitat)))

habitat_changes <- reduced_df %>%
  group_by(site) %>%
  summarise(habitats = n_distinct(habitat)) %>%
  filter(habitats > 1)

#track site change by soilT
T.track <- T.avni[which(T.avni$depth < 35.1 & T.avni$depth > 31.9),]

T.track$Tsoil_habitat <- ifelse(T.track$Tsoil < 2.5, "Palsa", 
                         ifelse(T.track$Tsoil > 9, "Fen", "Bog"))

#samples have inconsistant name over years, change them.
T.track$sample <- gsub("AvC-P1|AvnC-P", "AvC-P", T.track$sample)
T.track$sample <- gsub("AvC-C", "AvCe-C", T.track$sample)
T.track$sample <- gsub("Inc-S", "IncS", T.track$sample)
T.track$sample <- gsub("Inc-E Auto|Inc-E", "IncE", T.track$sample)
T.track$sample <- gsub("AvC-F1|AvCe_F1|AvCe-F", "AvCe-F1", T.track$sample)

metadata <- read.csv("metadata.csv", row.names = 1, header = T)
metadata$Year <- factor(metadata$Year, levels = c(2013:2023)) 

# Filter out the values "BOG", "PALSA", and "FEN"
Avni_metabolites_sample <- metadata$Habitat_plant[!(metadata$Habitat_plant %in% c("BOG", "PALSA", "FEN"))]
# Get the unique values from the filtered list
Avni_metabolites_sample <- unique(Avni_metabolites_sample)
Avni_metabolites_sample <- gsub("_", "-", Avni_metabolites_sample)

T.track.subset <- subset(T.track, T.track$sample %in%Avni_metabolites_sample)

T.track.subset.compare <- T.track.subset[,c(1,7,12,13,14)]
