library(vegan)
library(tidyverse)
library(gridExtra)
library(mixOmics)
library(ggrepel)
library(colorspace)

setwd("~/PERMA/R")

#load data and function -------
metadata <- read.csv("metadata.csv", row.names = 1, header = T)
metadata$Year <- factor(metadata$Year, levels = c(2013:2023)) 
all.row.names <- rownames(metadata)
# cleaning core data and merge with metadata ----
core.metadata <- read.csv("core_metadata.csv")
core.metadata <- core.metadata %>%
  mutate(sample = dplyr::recode(
    sample,
    "Palsa" = "P",
    "Sphagnum" = "S",
    "Eriophorum" = "E",
    "Bog" = "S",
    "Fen" = "E",
    "AvCe-F" = "AvCe-F1"
  ))

core.metadata$depth_interval <- gsub("-", "_", core.metadata$depth_interval)
core.metadata$Sample_code_by_core <- paste0(core.metadata$year, "_", core.metadata$sample, "_", core.metadata$core, "_",core.metadata$depth_interval )
core.metadata$Sample_code_by_core <- gsub("-", "_", core.metadata$Sample_code_by_core)

metadata["2019_P2_30_34", "Sample_code"] <- "core_2"
metadata$depth <- gsub("90_92", "90_94", metadata$depth)
metadata$label <-metadata$Habitat_plant
metadata <- metadata %>%
  mutate(label = dplyr::recode(label, 
                        "PALSA" = "P", 
                        "BOG" = "S", 
                        "FEN" = "E",
  ))
metadata$Sample_code_by_core <- paste0(metadata$Year, "_", metadata$label, "_", substr(metadata$Sample_code, nchar(metadata$Sample_code), nchar(metadata$Sample_code)), "_", metadata$depth)
metadata <- left_join(metadata, core.metadata[,c(4:10,14)])
rownames(metadata) <- all.row.names
write.csv(metadata, "metadata.all.csv", row.names =  T)

moi <- read.csv("soil_moi_all.csv", sep = ",")
colnames(moi)[colnames(moi) == "Moisture_contend..."] <- "Moi"

metadata <- left_join(metadata, moi)
metadata$Sample_code_by_core <- gsub("2022_IncE_2", "2022_IncE", metadata$Sample_code_by_core)
rownames(metadata) <- metadata$Sample_code_by_core
rownames(metadata) <- gsub("_", "", rownames(metadata))

HNEG <- read.csv("HNEG.csv", row.names = 1, header = T) %>%
  t(.) %>%
  as.data.frame(.) 

new_row_names <- rownames(HNEG) %>%
  gsub("^X", "", .) %>%  # Remove "X" at the beginning of each row name
  gsub("\\.raw\\.\\.F.*", "", .)%>%
  gsub("_", "", .) # Remove everything after ".raw..F" including ".raw..F"
rownames(HNEG) <- new_row_names

HNEG <- HNEG[rownames(metadata),]
#2014_S_2_20_24 was lost during processing for HILIC

RPPOS <- read.csv("RPPOS.csv", header = T, row.names = 1) %>%
  t(.) %>%
  as.data.frame(.)
new_row_names <- rownames(RPPOS) %>%
  gsub("^X", "", .) %>%  # Remove "X" at the beginning of each row name
  gsub("\\.raw\\.\\.F.*", "", .) %>%  # Remove everything after ".raw..F" including ".raw..F"
  gsub("_", "", .) 
rownames(RPPOS) <- new_row_names

#remove insource fragments
frag_RPPOS <- read.csv("~/PERMA/in_source_RP.csv")
RPPOS_remove <- frag_RPPOS %>%
  # Extract the fragments column
  pull(fragments) %>%
  # Split each entry by ';' and remove spaces
  str_split(';') %>%
  # Flatten the list into a single vector
  unlist() %>%
  # Trim whitespace from the fragments
  str_trim()

frag_HNEG <- read.csv("~/PERMA/in_source_HN.csv")
HNEG_remove <- frag_HNEG %>%
  # Extract the fragments column
  pull(fragments) %>%
  # Split each entry by ';' and remove spaces
  str_split(';') %>%
  # Flatten the list into a single vector
  unlist() %>%
  # Trim whitespace from the fragments
  str_trim()


RPPOS <- RPPOS%>%
  dplyr::select(-any_of(RPPOS_remove))
RPPOS <- RPPOS[rownames(metadata),]

HNEG <- HNEG%>%
  dplyr::select(-any_of(HNEG_remove))
HNEG <- HNEG[rownames(metadata),]


compounds <- read.csv("compound_link.csv", header = T, row.names = 1)

#Pareto Scaling:
PS_helper <- function(x) {
  (x - mean(x)) / sqrt(sd(x, na.rm = T))
}	

pareto_scale <- function(x){
  mtb_scaled <- data.frame(apply(x, 2, PS_helper))
  return(mtb_scaled)
}
#Log Transformation Functions:
log_helper <- function(x, min.val) {
  log2((x + sqrt(x ^ 2 + min.val ^ 2)) / 2)
}


#Log Scaling:
log_transform <- function(x){
  x_nz <- x[ ,which(apply(x, 2, sum) != 0)] # remove 0 features
  min.val <- min(abs(x_nz[x_nz!=0]))/10
  x_log_trans <- data.frame(apply(x_nz, 2, log_helper, min.val))
  return(x_log_trans)
}



#### ordination, limit depth to 34 

#HNEG_year changes,limite depth to 1-34
HNEG.pareto<- na.omit(HNEG) %>%
  log_transform() %>%
  pareto_scale()

RPPOS.pareto<- na.omit(RPPOS) %>%
  log_transform() %>%
  pareto_scale()

metadata$Tsoil <- as.numeric(metadata$Tsoil)
metadata$pH <- as.numeric(metadata$pH)
metadata$Moi  <- as.numeric(metadata$Moi)

metadata.HNEG.year <- subset(metadata, depth %in% c("1_5", "10_14", "20_24", "30_34"))
metadata.HNEG.year$depth <- factor(metadata.HNEG.year$depth, levels = c("1_5", "10_14", "20_24", "30_34"))
metadata.HNEG.year<- metadata.HNEG.year %>%
  arrange(Habitat_plant, depth)
#one sampe is missing in HNEG
metadata.HNEG.year <- metadata.HNEG.year[rownames(metadata.HNEG.year) != "2014S22024", ]
metadata.HNEG.ac <- metadata.HNEG.year[grepl("Autochamber", metadata.HNEG.year$site), ]

metadata.RPPOS.year <- subset(metadata, depth %in% c("1_5", "10_14", "20_24", "30_34"))
metadata.RPPOS.year$depth <- factor(metadata.RPPOS.year$depth, levels = c("1_5", "10_14", "20_24", "30_34"))

metadata.RPPOS.ac <- metadata.RPPOS.year[grepl("Autochamber", metadata.RPPOS.year$site), ]
metadata.RPPOS.ac$Habitat_plant <- factor(metadata.RPPOS.ac$Habitat_plant, levels = c("PALSA", "BOG", "FEN"))

# Convert Year to numeric (but keep as discrete values in the plot)
metadata.RPPOS.ac.t <- metadata.RPPOS.ac %>%
  mutate(Year = as.numeric(as.character(Year)))  # Ensure Year is numeric

# Create the plot
Tsoil_plot <- ggplot(metadata.RPPOS.ac.t, aes(x = Year, y = Tsoil, color = Habitat_plant, linetype = as.factor(depth))) +
  geom_point(aes(shape = as.factor(depth)), size = 3, alpha = 0.7, na.rm = TRUE) +  # Scatter plot with different shapes for depth
  geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = F, linewidth = 1, na.rm = TRUE) + 
  facet_wrap(~ Habitat_plant, scales = "free_y", ncol = 1) +  # Facet by Habitat_plant
  scale_color_manual(values = c("PALSA" = "#B8912F", "BOG" = "#5E813F", "FEN" = "#4273B0")) +  # Assign colors to habitats
  scale_shape_manual(values = c(16, 17, 3, 4)) +  # Different shapes for depth
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash")) +  # Different line types for depth
  scale_x_continuous(breaks = unique(metadata.RPPOS.ac.t$Year), labels = as.integer(unique(metadata.RPPOS.ac.t$Year))) +  # Ensure integer breaks for Year
  theme_classic() +
  labs(x = "Year", y = "Tsoil (°C)", color = "Habitat", shape = "Depth", linetype = "Depth") +
  theme(legend.position = "bottom")

ggsave("Tsoil_decadal_change.pdf",Tsoil_plot, height = 8, width = 3) 
### the follwing is just for MALAK to share with EMERGE people -----

site_color <- c("#B8912F","#5E813F","#4273B0", "#5E813F","#4273B0","purple4","red4","red1","#B8912F")
names(site_color) <-  c("PALSA","BOG","FEN","IncS","IncE", "AvCe_F1", "AvCe_C", "AvC_P", "AvCj_P")
habitat_order <- c("PALSA","AvCj_P","BOG", "IncS", "AvC_P","AvCe_C", "AvCe_F1","FEN","IncE")

# Convert Habitat_plant to a factor with the custom order
metadata.RPPOS.year$Habitat_plant <- factor(metadata.RPPOS.year$Habitat_plant, levels = habitat_order)

ggplot(metadata.RPPOS.year, aes(x = as.factor(Habitat_plant), y = Moi, color = Habitat_plant)) +
  geom_boxplot() +
  geom_point()+
  scale_color_manual(values = site_color) +
  facet_wrap(~ depth) +
  theme_minimal() +
  labs(x = "habitat", y = "Moisture Content", color = "Habitat Plant") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_nmds <- function(metadata, pareto_data, group_column, group_value, dataset_name) {
  
  # Subset data for the specified group
  subset_data <- metadata[metadata[[group_column]] == group_value,]
  a <- pareto_data[rownames(subset_data),]
  # Compute distance matrix and NMDS
  dist <- vegdist(pareto_data[rownames(subset_data),], method = "euclidean")
  nmds_result <- metaMDS(dist)
  
  adonis_result <- adonis2(dist ~ Habitat_plant+depth, data = subset_data)
  r2_habitat <- round(adonis_result$R2[1], 3)
  r2_depth <- round(adonis_result$R2[2], 3)
  stress_value <- round(nmds_result$stress, 3)
  
  subset_data$NMDS1 <- nmds_result$points[,1]
  subset_data$NMDS2 <- nmds_result$points[,2]
  # 
  # pareto.fit <- envfit(nmds_result, subset_data[,c("Tsoil","pH", "Moi")], perm = 9999, na.rm = TRUE)
  # pareto.fit.df <- as.data.frame(pareto.fit$vectors$arrows*pareto.fit$vectors$r*50)
  # pareto.fit.df$species <- rownames(pareto.fit.df)
  # pareto.fit.df$P <- pareto.fit$vectors$pvals
  # pareto.fit.df$R <- pareto.fit$vectors$r
  # pareto.fit.df <- pareto.fit.df[which(pareto.fit.df$P < 0.05),]
  # 
  # Compute centroids for each Habitat_plant group
  centroids <- subset_data %>%
    group_by(Habitat_plant) %>%
    summarise(NMDS1 = mean(NMDS1, na.rm = TRUE), 
              NMDS2 = mean(NMDS2, na.rm = TRUE))
  
  # Plot NMDS
  p <- ggplot(data = subset_data, aes(x = NMDS1, y = NMDS2)) +
    # Plot ellipses for each Habitat_plant group
    stat_ellipse(aes(color = Habitat_plant), level = 0.95, linewidth = 1, alpha = 0.5) +
    # Add centroids
    geom_point(data = centroids, aes(x = NMDS1, y = NMDS2, color = Habitat_plant), size = 4, shape = 15) +
    
    # Customize color scale
    scale_color_manual(values = site_color) +
    
    # Add envfit vectors
    #geom_segment(data = pareto.fit.df, aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2), 
    #             arrow = arrow(length = unit(0.3, "cm")), colour = "black") +
    
    #geom_text_repel(data = pareto.fit.df, aes(x = NMDS1, y = NMDS2, label = species), 
    #                size = 4, max.overlaps = Inf, alpha = 0.5) +
    
    # Add titles and labels
    ggtitle(paste(dataset_name, group_value)) +
    #annotate("text", x = min(subset_data$NMDS1), y = min(subset_data$NMDS2),
    #         label = paste("\nR²(Habitat):", r2_habitat, "\nR²(Depth):", r2_depth, "\nStress:", stress_value), 
    #         hjust = 0, vjust = 0, size = 4, color = "black") +
    
    theme_classic()
  
  summary_table <- data.frame(
    Year = group_value,
    R2_Habitat = r2_habitat,
    R2_Depth = r2_depth,
    Stress = stress_value
  )
  
  return(list(plot = p, summary = summary_table))
}

# Example usage for RPPOS dataset and different depths
depths <- c("1_5", "10_14", "20_24", "30_34")
plots_RPPOS <- list()
plots_HNEG <- list()
for (depth in depths) {
  plots_RPPOS[[depth]] <- print(plot_nmds(metadata.RPPOS.year, RPPOS.pareto, "depth", depth, "RPPOS"))
}

# Example usage for HNEG dataset and different depths
for (depth in depths) {
  plots_HNEG[[depth]] <- print(plot_nmds(metadata.HNEG.year, HNEG.pareto, "depth", depth, "HNEG"))
}

# plot by year.
plots_RPPOS_year <- list()
summary_tables_RPPOS_year <- list()

years <- c(2013:2020, 2022:2023)

# Loop over each year to generate plots and summary tables
for (year in years) {
  result <- plot_nmds(metadata.RPPOS.ac, RPPOS.pareto, "Year", year, "RPPOS")
  
  # Store the plot and summary table for each year
  plots_RPPOS_year[[as.character(year)]] <- result$plot
  summary_tables_RPPOS_year[[as.character(year)]] <- result$summary
}


plots_HNEG_year <- list()
summary_tables_HNEG_year <- list()

years <- c(2013:2020, 2022:2023)

# Loop over each year to generate plots and summary tables
for (year in years) {
  result <- plot_nmds(metadata.HNEG.ac, HNEG.pareto, "Year", year, "HNEG")
  
  # Store the plot and summary table for each year
  plots_HNEG_year[[as.character(year)]] <- result$plot
  summary_tables_HNEG_year[[as.character(year)]] <- result$summary
}
# Combine all summary tables into a single data frame
summary_table_all_years <- do.call(rbind, summary_tables_RPPOS_year)

library(ggpubr)

NMDS_contour_by_depth <- ggarrange(
  plots_HNEG[["1_5"]], plots_RPPOS[["1_5"]],
  plots_HNEG[["10_14"]], plots_RPPOS[["10_14"]],
  plots_HNEG[["20_24"]], plots_RPPOS[["20_24"]],
  plots_HNEG[["30_34"]], plots_RPPOS[["30_34"]],
  ncol = 2, nrow = 4,  # Define the number of columns and rows
  labels = c("A", "B", "C", "D", "E", "F", "G", "H"),  # Optional: Add labels for each subplot
  common.legend = TRUE,   # Use a common legend for all plots
  legend = "bottom"       # Position the legend at the bottom
)


NMDS_by_year <- ggarrange(
  # First row: Combine all plots from plots_HNEG_year
  plots_HNEG_year[[1]], plots_HNEG_year[[2]], plots_HNEG_year[[3]], plots_HNEG_year[[4]], 
  plots_HNEG_year[[5]], plots_HNEG_year[[6]], plots_HNEG_year[[7]], plots_HNEG_year[[8]], 
  plots_HNEG_year[[9]], plots_HNEG_year[[10]],
  
  # Second row: Combine all plots from plots_RPPOS_year
  plots_RPPOS_year[[1]], plots_RPPOS_year[[2]], plots_RPPOS_year[[3]], plots_RPPOS_year[[4]], 
  plots_RPPOS_year[[5]], plots_RPPOS_year[[6]], plots_RPPOS_year[[7]], plots_RPPOS_year[[8]], 
  plots_RPPOS_year[[9]], plots_RPPOS_year[[10]],
  # Define layout: 10 plots per row, 2 rows
  ncol = 10, nrow = 2,
  common.legend = TRUE,legend = "right"
)


#plot R2 plots for each year.
extract_r2_data <- function(summary_list, dataset_name) {
  bind_rows(lapply(names(summary_list), function(year) {
    data.frame(
      Year = as.numeric(year),
      R2_Habitat = summary_list[[year]]$R2_Habitat,
      R2_Depth = summary_list[[year]]$R2_Depth,
      Dataset = dataset_name
    )
  }))
}

# Extract data for both RPPOS and HNEG
rppos_data <- extract_r2_data(summary_tables_RPPOS_year, "RPPOS")
hneg_data <- extract_r2_data(summary_tables_HNEG_year, "HNEG")

# Combine and reshape data for stacking
plot_data <- bind_rows(rppos_data, hneg_data) %>%
  pivot_longer(cols = c(R2_Habitat, R2_Depth), names_to = "Factor", values_to = "R2")

# Define colors for Habitat and Depth
factor_colors <- c("R2_Habitat" = "#AB0520", "R2_Depth" = "#0C234B")

# Create the stacked bar plot
plot_NMDS_R2_by_year <- ggplot(plot_data, aes(x = as.factor(Year), y = R2, fill = Factor)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~Dataset, ncol = 2) +  # Separate RPPOS and HNEG
  scale_fill_manual(values = factor_colors) +
  labs(x = "Year", y = expression(R^2), fill = "Factor") +
  theme_classic() +
  theme(legend.position = "bottom")

plot_NMDS_R2_by_year <- ggplot(plot_data, aes(x = as.factor(Year), y = R2, fill = Factor)) +
  geom_bar(stat = "identity", position = "stack") +
  coord_flip() +  # Flip the coordinates for horizontal bars
  facet_wrap(~Dataset, ncol = 2) +  # Separate RPPOS and HNEG
  scale_fill_manual(values = factor_colors) +
  labs(x = "Year", y = expression(R^2), fill = "Factor") +
  theme_classic() +
  theme(legend.position = "bottom")

dist.rp.all <-  vegdist(RPPOS.pareto[rownames(metadata.RPPOS.year),], method = "euclidean")
nmds_result.rp.all <- metaMDS(dist.rp.all)
metadata.RPPOS.year$NMDS1 <- nmds_result.rp.all$points[,1]
metadata.RPPOS.year$NMDS2 <- nmds_result.rp.all$points[,2]

filtered_data <- metadata.RPPOS.year %>%
  mutate(Year = as.numeric(as.character(Year)))%>%
  filter(grepl("Autochamber", site))

nmds1_limits <- range(metadata.RPPOS.year$NMDS1, na.rm = TRUE)
nmds2_limits <- range(metadata.RPPOS.year$NMDS2, na.rm = TRUE)

single_year_data <- metadata.RPPOS.year %>%
  filter(Habitat_plant %in% c("AvC_P", "AvCe_C", "AvCe_F1"))

# Plot NMDS1 over time with smoothing
rp_ac_nmds1 <-ggplot(filtered_data, aes(x = Year, y = NMDS1, color = Habitat_plant, fill = Habitat_plant)) +
  geom_point(aes(shape = depth), alpha = 0.5, size = 2) +  # Individual points
  geom_smooth(method = "loess", se = TRUE, linewidth = 1) +  # Smoothed trend with confidence interval
  scale_color_manual(values = c("PALSA" = "#B8912F", "BOG" = "#5E813F", "FEN" = "#4273B0")) +
  scale_fill_manual(values = c("PALSA" = "#B8912F", "BOG" = "#5E813F", "FEN" = "#4273B0")) +  
  scale_shape_manual(values = c(16, 17, 3, 4)) +
  theme_classic() +
  scale_x_continuous(breaks = unique(filtered_data$Year), labels = as.integer(unique(filtered_data$Year))) +  # Ensure integer breaks for Year
  labs(x = "Year", y = "NMDS1", color = "Habitat") +
  ylim(nmds1_limits) +  
  theme(legend.position = "none")

rp_transition_nmds1 <- ggplot(single_year_data, aes(x = Habitat_plant, y = NMDS1, fill = Habitat_plant)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7, aes(color = Habitat_plant, shape = depth)) +  # Add jittered points
  scale_fill_manual(values = c("AvC_P" = "red1", "AvCe_C" = "red4", "AvCe_F1" = "purple4")) +
  scale_color_manual(values = c("AvC_P" = "red1", "AvCe_C" = "red4", "AvCe_F1" = "purple4")) +
  theme_classic() +
  labs(y = "NMDS1", x = "Habitat") +
  ylim(nmds1_limits) +  # Ensure same y-range as NMDS1 temporal plot
  theme(legend.position = "none")

rp_ac_nmds2 <- ggplot(filtered_data, aes(x = Year, y = NMDS2, color = Habitat_plant, fill = Habitat_plant)) +
  geom_point(aes(shape = depth), alpha = 0.5, size = 2) +  # Individual points
  geom_smooth(method = "loess", se = TRUE, linewidth = 1) +  # Smoothed trend with confidence interval
  scale_color_manual(values = c("PALSA" = "#B8912F", "BOG" = "#5E813F", "FEN" = "#4273B0")) +
  scale_fill_manual(values = c("PALSA" = "#B8912F", "BOG" = "#5E813F", "FEN" = "#4273B0")) +  
  theme_classic() +
  scale_x_continuous(breaks = unique(filtered_data$Year), labels = as.integer(unique(filtered_data$Year))) +  # Ensure integer breaks for Year
  labs(x = "Year", y = "NMDS2", color = "Habitat") +
  ylim(nmds2_limits) + 
  theme(legend.position = "none")

rp_transition_nmds2 <-ggplot(single_year_data, aes(x = Habitat_plant, y = NMDS2, fill = Habitat_plant)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7, aes(color = Habitat_plant, shape = depth)) +  # Add jittered points
  scale_fill_manual(values = c("AvC_P" = "red1", "AvCe_C" = "red4", "AvCe_F1" = "purple4")) +
  scale_color_manual(values = c("AvC_P" = "red1", "AvCe_C" = "red4", "AvCe_F1" = "purple4")) +
  theme_classic() +
  labs(y = "NMDS2", x = "Habitat") +
  ylim(nmds2_limits) +  # Ensure same y-range as NMDS1 temporal plot
  theme(legend.position = "none")

library(ggpubr)

rp_combined_plot <- ggarrange(
  rp_ac_nmds1, rp_transition_nmds1,        # Place p1 and p3 in the top row
  rp_ac_nmds2, rp_transition_nmds2,        # Place p2 and p4 in the bottom row
  ncol = 2,       # Number of columns
  nrow = 2,       # Number of rows
  labels = c("A", "B", "C", "D"), # Optional labels for each plot
  heights = c(1, 1), # Control the height of the rows
  widths = c(3, 1)   # Control the width of the columns
)

# Plot NMDS1 of HN over time with smoothing


dist.hn.all <-  vegdist(HNEG.pareto[rownames(metadata.HNEG.year),], method = "euclidean")
nmds_result.hn.all <- metaMDS(dist.hn.all)
metadata.HNEG.year$NMDS1 <- nmds_result.hn.all$points[,1]
metadata.HNEG.year$NMDS2 <- nmds_result.hn.all$points[,2]

filtered_data <- metadata.HNEG.year %>%
  mutate(Year = as.numeric(as.character(Year)))%>%
  filter(grepl("Autochamber", site))

nmds1_limits <- range(metadata.HNEG.year$NMDS1, na.rm = TRUE)
nmds2_limits <- range(metadata.HNEG.year$NMDS2, na.rm = TRUE)

single_year_data <- metadata.HNEG.year %>%
  filter(Habitat_plant %in% c("AvC_P", "AvCe_C", "AvCe_F1"))

hn_ac_nmds1 <-ggplot(filtered_data, aes(x = Year, y = NMDS1, color = Habitat_plant, fill = Habitat_plant)) +
  geom_point(aes(shape = depth), alpha = 0.5, size = 2) +  # Individual points
  geom_smooth(method = "loess", se = TRUE, linewidth = 1) +  # Smoothed trend with confidence interval
  scale_color_manual(values = c("PALSA" = "#B8912F", "BOG" = "#5E813F", "FEN" = "#4273B0")) +
  scale_fill_manual(values = c("PALSA" = "#B8912F", "BOG" = "#5E813F", "FEN" = "#4273B0")) +  
  scale_shape_manual(values = c(16, 17, 3, 4)) +
  theme_classic() +
  scale_x_continuous(breaks = unique(filtered_data$Year), labels = as.integer(unique(filtered_data$Year))) +  # Ensure integer breaks for Year
  labs(x = "Year", y = "NMDS1", color = "Habitat") +
  ylim(nmds1_limits) +  
  theme(legend.position = "none")

hn_transition_nmds1 <- ggplot(single_year_data, aes(x = Habitat_plant, y = NMDS1, fill = Habitat_plant)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7, aes(color = Habitat_plant, shape = depth)) +  # Add jittered points
  scale_fill_manual(values = c("AvC_P" = "red1", "AvCe_C" = "red4", "AvCe_F1" = "purple4")) +
  scale_color_manual(values = c("AvC_P" = "red1", "AvCe_C" = "red4", "AvCe_F1" = "purple4")) +
  theme_classic() +
  labs(y = "NMDS1", x = "Habitat") +
  ylim(nmds1_limits) +  # Ensure same y-range as NMDS1 temporal plot
  theme(legend.position = "none")

hn_ac_nmds2 <- ggplot(filtered_data, aes(x = Year, y = NMDS2, color = Habitat_plant, fill = Habitat_plant)) +
  geom_point(aes(shape = depth), alpha = 0.5, size = 2) +  # Individual points
  geom_smooth(method = "loess", se = TRUE, linewidth = 1) +  # Smoothed trend with confidence interval
  scale_color_manual(values = c("PALSA" = "#B8912F", "BOG" = "#5E813F", "FEN" = "#4273B0")) +
  scale_fill_manual(values = c("PALSA" = "#B8912F", "BOG" = "#5E813F", "FEN" = "#4273B0")) +  
  theme_classic() +
  scale_x_continuous(breaks = unique(filtered_data$Year), labels = as.integer(unique(filtered_data$Year))) +  # Ensure integer breaks for Year
  labs(x = "Year", y = "NMDS2", color = "Habitat") +
  ylim(nmds2_limits) + 
  theme(legend.position = "none")

hn_transition_nmds2 <-ggplot(single_year_data, aes(x = Habitat_plant, y = NMDS2, fill = Habitat_plant)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7, aes(color = Habitat_plant, shape = depth)) +  # Add jittered points
  scale_fill_manual(values = c("AvC_P" = "red1", "AvCe_C" = "red4", "AvCe_F1" = "purple4")) +
  scale_color_manual(values = c("AvC_P" = "red1", "AvCe_C" = "red4", "AvCe_F1" = "purple4")) +
  theme_classic() +
  labs(y = "NMDS2", x = "Habitat") +
  ylim(nmds2_limits) +  # Ensure same y-range as NMDS1 temporal plot
  theme(legend.position = "none")


hn_combined_plot <- ggarrange(
  hn_ac_nmds1, hn_transition_nmds1,        # Place p1 and p3 in the top row
  hn_ac_nmds2, hn_transition_nmds2,        # Place p2 and p4 in the bottom row
  ncol = 2,       # Number of columns
  nrow = 2,       # Number of rows
  labels = c("A", "B", "C", "D"), # Optional labels for each plot
  heights = c(1, 1), # Control the height of the rows
  widths = c(3, 1)   # Control the width of the columns
)

ggsave("hn_temporal_nmds.pdf",hn_combined_plot, width = 6, height = 4)
ggsave("rp_temporal_nmds.pdf",rp_combined_plot, width = 6.5, height = 4)

ggsave("R2_change_by_year.pdf", plot_NMDS_R2_by_year, width = 4, height = 4)


plot_nmds_all <- function(metadata, pareto_data, dataset_name) {
  # Subset data for the specified depth
  
  # Compute distance matrix and NMDS
  dist <- vegdist(pareto_data[rownames(metadata),], method = "euclidean")
  nmds_result <- metaMDS(dist)
  
  adonis_result <- adonis2(dist ~ Habitat_plant, data = metadata)
  r2_habitat <- round(adonis_result$R2[1], 3)
  stress_value <- round(nmds_result$stress, 3)
  
  metadata$NMDS1 <- nmds_result$points[,1]
  metadata$NMDS2 <- nmds_result$points[,2]
  
  
  # Plot NMDS
  p <- ggplot(data = metadata, aes(x = NMDS1, y = NMDS2)) +
    # Plot density contours
    geom_point(data = metadata[grepl("Autochamber", metadata$site), ],
               aes(color = Habitat_plant), size = 1, alpha = 0.5) +
    
    # Plot points
    geom_point(data = metadata[!grepl("Autochamber", metadata$site), ],
               aes(color = Habitat_plant), size = 2, shape = 17) +
    
    # Customize color scale
    scale_color_manual(values = site_color) +
    # Add titles and labels
    ggtitle(paste(dataset_name)) +
    annotate("text", x = min(metadata$NMDS1), y = min(metadata$NMDS2),
             label = paste("\nR²(Habitat):", r2_habitat, "\nStress:", stress_value), 
             hjust = 0, vjust = 0, size = 4, color = "black") +
    theme_minimal()
  return(p)
}

NMDS_RP<-plot_nmds_all(metadata.RPPOS.year, RPPOS.pareto, "RPPOS")
NMDS_HN<-plot_nmds_all(metadata.HNEG.year, HNEG.pareto, "HNEG")

NMDS.all <- ggarrange(NMDS_RP, NMDS_HN,ncol = 2,labels = c("A", "B"),common.legend = TRUE, legend = "bottom")

ggsave("NMDS_together_EMRGE.pdf", NMDS.all, width = 10, height = 5)

ggsave("NMDS_contour_by_depth_EMRGE.pdf", NMDS_contour_by_depth, width = 10, height = 16)

### all factors together but only stable sites.----


### focus on difference betwee avni fen, palsa and autochamber palsa and fen.----
site_color <- c("#5E813F","#4273B0","purple4","red4","red1", "#B8912F")
names(site_color) <-  c("IncS","IncE", "AvCe_F1", "AvCe_C", "AvC_P", "AvCj_P")



plot_nmds <- function(metadata, pareto_data, depth, dataset_name) {
  # Subset data for the specified depth
  subset_data <- metadata[metadata$depth == depth,] 
  # Compute distance matrix and NMDS
  dist <- vegdist(pareto_data[rownames(subset_data),], method = "euclidean")
  nmds_result <- metaMDS(dist)
  
  adonis_result <- adonis2(dist ~ Habitat_plant, data = subset_data)
  r2_habitat <- round(adonis_result$R2[1], 3)
  stress_value <- round(nmds_result$stress, 3)

  subset_data$NMDS1 <- nmds_result$points[,1]
  subset_data$NMDS2 <- nmds_result$points[,2]
  
   # compounds$ID <- rownames(compounds)
   # 
   # pareto.classy <- pareto_data %>%
   #   t(.)%>%
   #   as.data.frame(.)%>%
   #   rownames_to_column("ID")%>%
   #   left_join(., compounds) %>% 
   #   column_to_rownames("ID")
   # 
   # pareto.classy_mean <- pareto.classy %>%
   #   group_by(Superclass) %>%
   #   summarize(across(where(is.numeric), \(x) mean(x, na.rm = TRUE))) %>%
   #   t(.) %>%
   #   as.data.frame(.)
   # 
   # colnames(pareto.classy_mean) <- pareto.classy_mean[1, ]
   # pareto.classy_mean <- pareto.classy_mean[-1, ]
   # pareto.classy_mean <- pareto.classy_mean[, !(is.na(colnames(pareto.classy_mean)) | colnames(pareto.classy_mean) == "")]
   # 
   # pareto.classy_mean <- pareto.classy_mean[rownames(subset_data),]
   # pareto.classy_mean<- as.data.frame(apply(pareto.classy_mean, 2, as.numeric))
   
   pareto.fit <- envfit(nmds_result, subset_data[,c("Tsoil","pH", "Moi")], perm = 9999, na.rm = TRUE)
   pareto.fit.df <- as.data.frame(pareto.fit$vectors$arrows*pareto.fit$vectors$r*50)
   pareto.fit.df$species <- rownames(pareto.fit.df)
   pareto.fit.df$P <- pareto.fit$vectors$pvals
   pareto.fit.df$R <- pareto.fit$vectors$r
   pareto.fit.df <- pareto.fit.df[which(pareto.fit.df$P < 0.05),]
   
   
  
  
  # Plot NMDS
  p <- ggplot(data = subset_data, aes(x = NMDS1, y = NMDS2)) +
    # Plot density contours
    stat_density_2d(data = subset_data[subset_data$site == "Palsa Autochamber Site", ],
                    aes(color = Habitat_plant), linewidth = 0.5, color = "#B8912F" , contour = TRUE, alpha = 0.5) +
    
    stat_density_2d(data = subset_data[subset_data$site == "Eriophorum Autochamber Site", ],
                    aes(color = Habitat_plant), linewidth = 0.5, color = "#4273B0", contour = TRUE, alpha = 0.5) +
    
    stat_density_2d(data = subset_data[subset_data$site == "Sphagnum Autochamber Site", ],
                    aes(color = Habitat_plant), linewidth = 0.5, color = "#5E813F", contour = TRUE, alpha = 0.5) +
    
    # Plot points
    geom_point(data = subset_data[!grepl("Autochamber", subset_data$site), ],
               aes(color = Habitat_plant), size = 2) +
    
    # Customize color scale
    scale_color_manual(values = site_color) +
    
    #envfit
    geom_segment(data = pareto.fit.df, aes (x = 0, xend = NMDS1, y = 0, yend = NMDS2), arrow = arrow(length = unit(0.3, "cm")), colour = "black")+
    #   geom_text_repel(data = HNEG.fit.df, aes (x = NMDS1, y = NMDS2, label = species), size = 5, max.overlaps = Inf)
    geom_text_repel(data = pareto.fit.df, aes (x = NMDS1, y = NMDS2, label = species), size = 4, max.overlaps = Inf, alpha = 0.5)+
    # Add titles and labels
    ggtitle(paste(dataset_name, depth)) +
    annotate("text", x = min(subset_data$NMDS1), y = min(subset_data$NMDS2),
             label = paste("\nR²(Habitat):", r2_habitat, "\nStress:", stress_value), 
             hjust = 0, vjust = 0, size = 4, color = "black") +
    theme_classic()
  return(p)
}

# Example usage for RPPOS dataset and different depths
depths <- c("1_5", "10_14", "20_24", "30_34")
plots_RPPOS <- list()
plots_HNEG <- list()
for (depth in depths) {
  plots_RPPOS[[depth]] <- print(plot_nmds(metadata.RPPOS.year, RPPOS.pareto, depth, "RPPOS"))
}

# Example usage for HNEG dataset and different depths
for (depth in depths) {
  plots_HNEG[[depth]] <- print(plot_nmds(metadata.HNEG.year, HNEG.pareto, depth, "HNEG"))
}

plot_nmds(metadata.HNEG.year, HNEG.pareto, "20_24", "HNEG")

metadata.HNEG.mineral <-  metadata.HNEG.year[!metadata.HNEG.year$Habitat_plant %in% c("IncS", "IncE", "AvCj_P"), ]

# Filter for specific Habitat_plant and Sample_code conditions
metadata.HNEG.mineral <- metadata.HNEG.mineral[
  !(metadata.HNEG.mineral$Habitat_plant == "AvCe_F1" & metadata.HNEG.mineral$Sample_code != "core _2") &
    !(metadata.HNEG.mineral$Habitat_plant == "AvCe_C" & metadata.HNEG.mineral$Sample_code != "core _2") &
    !(metadata.HNEG.mineral$Habitat_plant == "AvC_P" & metadata.HNEG.mineral$Sample_code != "core _1"), 
]
HNEG.pareto.mineral <- HNEG.pareto[rownames(metadata.HNEG.mineral),]

metadata.RPPOS.mineral <-  metadata.RPPOS.year[!metadata.RPPOS.year$Habitat_plant %in% c("IncS", "IncE", "AvCj_P"), ]

# Filter for specific Habitat_plant and Sample_code conditions
metadata.RPPOS.mineral <- metadata.RPPOS.mineral[
  !(metadata.RPPOS.mineral$Habitat_plant == "AvCe_F1" & metadata.RPPOS.mineral$Sample_code != "core _2") &
    !(metadata.RPPOS.mineral$Habitat_plant == "AvCe_C" & metadata.RPPOS.mineral$Sample_code != "core _2") &
    !(metadata.RPPOS.mineral$Habitat_plant == "AvC_P" & metadata.RPPOS.mineral$Sample_code != "core _1"), 
]
RPPOS.pareto.mineral <- RPPOS.pareto[rownames(metadata.RPPOS.mineral),]

S4a <- plot_nmds(metadata.HNEG.mineral, HNEG.pareto.mineral, "20_24", "HNEG")
S4c <- plot_nmds(metadata.HNEG.mineral, HNEG.pareto.mineral, "30_34", "HNEG")

S4b <- plot_nmds(metadata.RPPOS.mineral, RPPOS.pareto.mineral, "20_24", "RPPOS")
S4d <- plot_nmds(metadata.RPPOS.mineral, RPPOS.pareto.mineral, "30_34", "RPPOS")
library(ggpubr)

NMDS_mineral <- ggarrange(
  S4a, S4b, S4c, S4d,
  ncol = 2, nrow = 2,  # Define the number of columns and rows
  labels = c("A", "B", "C", "D"),  # Optional: Add labels for each subplot
  common.legend = TRUE,   # Use a common legend for all plots
  legend = "bottom"       # Position the legend at the bottom
)
# Save the combined plot
ggsave("NMDS_mineral.pdf", NMDS_mineral, width = 8, height = 6)


NMDS_contour_by_depth <- ggarrange(
  plots_HNEG[["1_5"]], plots_RPPOS[["1_5"]],
  plots_HNEG[["10_14"]], plots_RPPOS[["10_14"]],
  plots_HNEG[["20_24"]], plots_RPPOS[["20_24"]],
  plots_HNEG[["30_34"]], plots_RPPOS[["30_34"]],
  ncol = 2, nrow = 4,  # Define the number of columns and rows
  labels = c("A", "B", "C", "D", "E", "F", "G", "H"),  # Optional: Add labels for each subplot
  common.legend = TRUE,   # Use a common legend for all plots
  legend = "bottom"       # Position the legend at the bottom
)
# Save the combined plot
ggsave("NMDS_contour_by_depth_envfit.pdf", NMDS_contour_by_depth, width = 10, height = 16)


important_features <- read_csv("important_features_module_1224.csv")
ID <- important_features$ID

#violin plot for mineral layers
HNEG.log<- na.omit(HNEG) %>%
  log_transform()

HNEG.log <- HNEG.log[, colnames(HNEG.log) %in% ID]

RPPOS.log<- na.omit(RPPOS) %>%
  log_transform()

RPPOS.log <- RPPOS.log[, colnames(RPPOS.log) %in% ID]

prepare_data <- function(log_data, metadata) {
  log_data$Sample <- rownames(log_data)
  data_long <- log_data %>%
    pivot_longer(cols = -Sample, names_to = "Feature", values_to = "Intensity") %>%
    right_join(metadata, by = "Sample")
  
  return(data_long)
}
metadata.HNEG.year$Sample <- rownames(metadata.HNEG.year)
metadata.RPPOS.year$Sample <- rownames(metadata.RPPOS.year)

# Prepare data
HNEG_data <- prepare_data(HNEG.log, metadata.HNEG.year)
RPPOS_data <- prepare_data(RPPOS.log, metadata.RPPOS.year)

depth_colors <- c("1_5" = "#654321",  # Dark brown
                  "10_14" = "#987554",  # Transitional brown
                  "20_24" = "#e5d3b3",  # Transitional grey
                  "30_34" = "#dddad3")  # Grey

facet_order <- c("PALSA", "BOG", "FEN", "AvC_P", "AvCe_C", "AvCe_F1", "IncE", "IncS", "AvCj_P")


plot_violin <- function(data, title) {
  data$Habitat_plant <- factor(data$Habitat_plant, levels = facet_order)  # Set facet order
  ggplot(data, aes(x = factor(depth), y = Intensity, fill = depth)) +
    geom_violin(trim = TRUE, alpha = 0.8) +
    scale_fill_manual(values = depth_colors) +
    stat_summary(fun = mean, geom = "crossbar", color = "red", width = 0.3) +
    facet_wrap(~Habitat_plant, scales = "free", nrow = 3) +
    labs(title = title, x = "Depth", y = "Log-transformed Intensity") +
    theme_minimal() +
    theme(legend.position = "none")
}

# Generate plots
HNEG_plot <- plot_violin(HNEG_data, "Distribution of log-transfomred Feature Intensities in HNEG")
RPPOS_plot <- plot_violin(RPPOS_data, "Distribution of log-transfomred Feature Intensities in RPPOS")

combined_plot <- ggarrange(HNEG_plot, RPPOS_plot, ncol = 1, nrow = 2)

ggsave("log_intensity.pdf", combined_plot, width = 8.5, height = 11, units = "in")



# autochamber only
metadata.RPPOS.ac <- subset(metadata.RPPOS.year, grepl("Autochamber", metadata.RPPOS.year$site))
RPPOS.ac <- RPPOS[rownames(metadata.RPPOS.ac),]
RPPOS.ac$habitat <- metadata.RPPOS.ac$Habitat_plant
RPPOS.pareto.ac <- RPPOS.pareto[rownames(metadata.RPPOS.ac),]

metadata.HNEG.ac <- subset(metadata.HNEG.year, grepl("Autochamber", metadata.HNEG.year$site))
HNEG.ac <- HNEG[rownames(metadata.HNEG.ac),]
HNEG.ac$habitat <- metadata.HNEG.ac$Habitat_plant
HNEG.pareto.ac <- HNEG.pareto[rownames(metadata.HNEG.ac),]

write.table(RPPOS.ac,"all.rp.csv",sep = ",")
write.table(HNEG.ac,"all.hn.csv", sep = ",")

set.seed(123)
### auto chamber RP-----
dist.rp.ac <- vegdist(RPPOS.pareto.ac[rownames(metadata.RPPOS.ac),], method = "euclidean")
rp.result.ac<- adonis2(dist.rp.ac ~ Habitat_plant + depth + Year, data = metadata.RPPOS.ac)
#all significant p < 0.001, but habitat R2 = 0.233; Year 0.105, depth 0.064
adonis_table_rp <- as.data.frame(rp.result.ac)
adonis_table_rp$Variable <- rownames(adonis_table_rp)
adonis_table_rp <- adonis_table_rp[, c("Variable", "R2", "Pr(>F)")]

nmds_result.rp <- metaMDS(dist.rp.ac)
metadata.RPPOS.ac$NMDS1 <- nmds_result.rp$points[, 1]
metadata.RPPOS.ac$NMDS2 <- nmds_result.rp$points[, 2]
metadata.RPPOS.ac$Habitat_plant <- factor(metadata.RPPOS.ac$Habitat_plant, levels = c("PALSA", "BOG", "FEN"))

rp.fit <- envfit(nmds_result.rp, metadata.RPPOS.ac[,c("Tsoil","pH")], perm = 9999, na.rm = TRUE)
rp.fit.df <- as.data.frame(rp.fit$vectors$arrows*rp.fit$vectors$r*50)
rp.fit.df$species <- rownames(rp.fit.df)
rp.fit.df$P <- rp.fit$vectors$pvals
rp.fit.df$R <- rp.fit$vectors$r
rp.fit.df <- rp.fit.df[which(rp.fit.df$P < 0.05 & rp.fit.df$R > 0.3),]

# auto chamber HN-----
dist.hn.ac <- vegdist(HNEG.pareto.ac[rownames(metadata.HNEG.ac),], method = "euclidean")
hn.result.ac<- adonis2(dist.hn.ac ~ Habitat_plant + depth + Year, data = metadata.HNEG.ac)
#all significant p < 0.001, but habitat R2 = 0.278; Year 0.119, depth 0.06
# divide by habitat
nmds_result.hn <- metaMDS(dist.hn.ac)
metadata.HNEG.ac$NMDS1 <- nmds_result.hn$points[, 1]
metadata.HNEG.ac$NMDS2 <- nmds_result.hn$points[, 2]
metadata.HNEG.ac$Habitat_plant <- factor(metadata.HNEG.ac$Habitat_plant, levels = c("PALSA", "BOG", "FEN"))
adonis_table_hn <- as.data.frame(hn.result.ac)
adonis_table_hn$Variable <- rownames(adonis_table_hn)
adonis_table_hn <- adonis_table_hn[, c("Variable", "R2", "Pr(>F)")]

hn.fit <- envfit(nmds_result.hn, metadata.HNEG.ac[,c("Tsoil","pH")], perm = 9999, na.rm = TRUE)
hn.fit.df <- as.data.frame(hn.fit$vectors$arrows*hn.fit$vectors$r*50)
hn.fit.df$species <- rownames(hn.fit.df)
hn.fit.df$P <- hn.fit$vectors$pvals
hn.fit.df$R <- hn.fit$vectors$r
hn.fit.df <- hn.fit.df[which(hn.fit.df$P < 0.05 & hn.fit.df$R > 0.3),]

color <- c("PALSA" = "#B8912F", "BOG" = "#5E813F", "FEN" = "#4273B0")
shape <- c("1_5" = 15, "10_14" = 16, "20_24" = 18, "30_34" =17)

# Create the plot for ppt
RPPOS.base <- ggplot(metadata.RPPOS.ac, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = Habitat_plant, shape = depth), 
             size = 2) +
  scale_color_manual(values = color) +
  scale_shape_manual(values = shape) +
  geom_segment(data = rp.fit.df, aes (x = 0, xend = NMDS1, y = 0, yend = NMDS2), arrow = arrow(length = unit(0.3, "cm")), colour = "black")+
  #   geom_text_repel(data = HNEG.fit.df, aes (x = NMDS1, y = NMDS2, label = species), size = 5, max.overlaps = Inf)
  geom_text_repel(data = rp.fit.df, aes (x = NMDS1, y = NMDS2, label = species), size = 4, max.overlaps = Inf, alpha = 0.5)+
  theme_minimal() +
  theme(
    legend.position = "bottom", 
    legend.box = "vertical", 
    legend.margin = margin(),
  ) 

HNEG.base <- ggplot(metadata.HNEG.ac, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = Habitat_plant, shape = depth), 
             size = 2) +
  scale_color_manual(values = color) +
  scale_shape_manual(values = shape) +
  geom_segment(data = hn.fit.df, aes (x = 0, xend = NMDS1, y = 0, yend = NMDS2), arrow = arrow(length = unit(0.3, "cm")), colour = "black")+
  #   geom_text_repel(data = HNEG.fit.df, aes (x = NMDS1, y = NMDS2, label = species), size = 5, max.overlaps = Inf)
  geom_text_repel(data = hn.fit.df, aes (x = NMDS1, y = NMDS2, label = species), size = 4, max.overlaps = Inf, alpha = 0.5)+
  theme_minimal() +
  theme(
    legend.position = "bottom", 
    legend.box = "vertical", 
    legend.margin = margin(),
  ) 

# Create tables as grobs
table_rp <- tableGrob(adonis_table_rp, rows = NULL)
table_hn <- tableGrob(adonis_table_hn, rows = NULL)

# Adjust plots for stacking
RPPOS.base <- RPPOS.base + labs(title = "NMDS plot for long-term stable sites - RP mode")
HNEG.base <- HNEG.base + labs(title = "NMDS plot for long-term stable sites - Hilic mode")

# Combine NMDS plots and tables
RPPOS_final <- plot_grid(RPPOS.base, ggdraw() + draw_grob(table_rp), ncol = 1, rel_heights = c(3, 1))
HNEG_final <- plot_grid(HNEG.base, ggdraw() + draw_grob(table_hn), ncol = 1, rel_heights = c(3, 1))
library(ggpubr)
final_plot <- ggarrange(
  RPPOS_final, HNEG_final,
  nrow = 1, labels = c("A", "B"),
  common.legend = TRUE, legend = "bottom"
)

ggsave(
  filename = "NMDS_combined_plot_ac.pdf", 
  plot = final_plot, 
  width = 8, 
  height = 6
)
#let's do a envfit on relative abundance of classes; basically average the Pareto scaling results of same classes
#1. add the function pareto and average.


compounds$ID <- rownames(compounds)

RPPOS.classy <- RPPOS.pareto %>%
  t(.)%>%
  as.data.frame(.)%>%
  rownames_to_column("ID")%>%
  left_join(., compounds)

RPPOS.classy_mean <- RPPOS.classy %>%
  group_by(Superclass) %>%
  summarize(across(where(is.numeric), \(x) mean(x, na.rm = TRUE))) %>%
  t(.) %>%
  as.data.frame(.)

colnames(RPPOS.classy_mean) <- RPPOS.classy_mean[1, ]
RPPOS.classy_mean <- RPPOS.classy_mean[-1, ]
RPPOS.classy_mean <- RPPOS.classy_mean[rownames(metadata.RPPOS.year),]
RPPOS.classy_mean<- as.data.frame(apply(RPPOS.classy_mean, 2, as.numeric))

RPPOS.fit <- envfit(nmds_result.rp, RPPOS.classy_mean, perm = 9999)
RPPOS.fit.df <- as.data.frame(RPPOS.fit$vectors$arrows*RPPOS.fit$vectors$r*50)
RPPOS.fit.df$species <- rownames(RPPOS.fit.df)
RPPOS.fit.df$P <- RPPOS.fit$vectors$pvals
RPPOS.fit.df$R <- RPPOS.fit$vectors$r
RPPOS.fit.df <- RPPOS.fit.df[which(RPPOS.fit.df$P < 0.05),]

library(ggrepel)
RPPOS.autochambar <- RPPOS.base + 
  geom_segment(data = RPPOS.fit.df, aes (x = 0, xend = NMDS1, y = 0, yend = NMDS2), arrow = arrow(length = unit(0.3, "cm")), colour = "black")+
  geom_text_repel(data = RPPOS.fit.df, aes (x = NMDS1, y = NMDS2, label = species), size = 5, max.overlaps = Inf)



#lets' label the points by sample!
# Function to create segment data for a specific Habitat_plant
create_segment_data <- function(data, habitat_plant) {
  # Filter data by the specified Habitat_plant
  filtered_data <- data %>% filter(Habitat_plant == habitat_plant)
  
  # Sort the filtered data by sample_code and depth
  sorted_data <- filtered_data %>%
    arrange(filtered_data$Sample_code, depth)
  
  # Create segments data frame
  segment_data <- sorted_data %>%
    group_by(Sample_code) %>%
    mutate(
      xend = lead(NMDS1),
      yend = lead(NMDS2)
    ) %>%
    filter(!is.na(xend) & !is.na(yend)) %>%
    ungroup()
  
  return(segment_data)
}

habitat_plant <- "AvCj_P"
labeled.NMDS_HNEG <- function(habitat_plant) {
  segment_data <- create_segment_data(metadata.HNEG.year, habitat_plant)
  labeled_plot <- HNEG.base +
    geom_segment(data = segment_data, aes(x = NMDS1, y = NMDS2, xend = xend, yend = yend),
                 arrow = arrow(type = "closed", length = unit(0.05, "inches")), color = "red") 
    #geom_text_repel(data = label[which(label$Habitat_plant == habitat_plant),] +
     #                 aes(x = NMDS1, y = NMDS2, label = label), size = 3, color = "blue", nudge_y = 0.5, max.overlaps = 10)
  return(labeled_plot)
}

# ignore the following for now.

ggsave("AvCj_P.png", labeled.NMDS_HNEG("AvCj_P"), width = 5, height = 3, dpi = 300)
ggsave("AvC_P.png", labeled.NMDS_HNEG("AvC_P"), width = 5, height = 3, dpi = 300)
ggsave("AvCe_C.png", labeled.NMDS_HNEG("AvCe_C"), width = 5, height = 3, dpi = 300)
ggsave("AvCe_F1.png", labeled.NMDS_HNEG("AvCe_F1"), width = 5, height = 3, dpi = 300)
ggsave("IncE.png",labeled.NMDS_HNEG("IncE"), width = 5, height = 3, dpi = 300)
ggsave("IncS.png",labeled.NMDS_HNEG("IncS"), width = 5, height = 3, dpi = 300)



# test the effect of year on metabolomes
#min(HNEG.pareto[rownames(metadata.HNEG.year),]) = -6.885563, please 7 to make everything positive

HNEG.final <- cbind(metadata.HNEG.year, HNEG.pareto[rownames(metadata.HNEG.year),]+7)
HNEG.final$habitat <- factor(HNEG.final$habitat, levels = c("Palsa", "Bog", "Fen"))
HNEG.final.ac <- subset(HNEG.final, grepl("Autochamber", site))
depths <- c("1_5", "10_14", "20_24", "30_34")




#individual year
plot_nmds <- function(metadata, habitat, d) {
  # Subset metadata by Habitat_plant and depth
  metadata_subset <- subset(metadata, Habitat_plant == habitat & depth == d)
  # Perform NMDS analysis
  dist <- vegdist(HNEG.pareto[rownames(metadata_subset),], method = "euclidean")
  nmds_result <- metaMDS(dist)
  
  metadata_subset$NMDS1 <- nmds_result$points[, 1]
  metadata_subset$NMDS2 <- nmds_result$points[, 2]
  
  # Perform PERMANOVA
  perm_anova <- adonis2(dist ~ Year, data = metadata_subset)
  p_value <- perm_anova$`Pr(>F)`[1]
  R2 <- perm_anova$R2[1]
  
  # Calculate centroids and confidence intervals
  centroid <- metadata_subset %>%
    group_by(Year) %>%
    summarize(X = mean(NMDS1),Y = mean(NMDS2), ci_NMDS1 = 1.96 * sd(NMDS1) / sqrt(n()), ci_NMDS2 = 1.96 * sd(NMDS2) / sqrt(n()))
  
  
  segments_df <- centroid %>%
    arrange(Year) %>%
    mutate(next_X = lead(X),
           next_Y = lead(Y),
           next_Year = lead(Year))
  
  # Plot NMDS results with centroids and confidence intervals
  ggplot(centroid, aes(x = X, y = Y, color = as.factor(Year))) +
    geom_errorbar(aes(ymin = Y - ci_NMDS2, ymax = Y + ci_NMDS2), width = 0) +
    geom_errorbarh(aes(xmin = X - ci_NMDS1, xmax = X + ci_NMDS1), height = 0) +
    geom_point(size = 3) +
    geom_segment(data = segments_df, aes(x = X, y = Y, xend = next_X, yend = next_Y), arrow = arrow(length = unit(0.03, "npc")), color = "black") +
    theme_bw() +
    theme(legend.position = "bottom")+
    labs(caption = paste("p-value:", round(p_value, 4), "R2:", round(R2, 4))) +
    ggtitle(paste("NMDS Plot for", habitat, d))
}

ggsave("test.pdf", plot_nmds(metadata.HNEG.year, "Erio", "10_14"), width = 4, height = 4)

plots_list <- list()

# Loop through each combination of plant habitat and depth
for (depth in c("1_5", "10_14", "20_24", "30_34")) {
  for (habitat in c("Erio", "Palsa", "Sphagnum")){
    # Generate the plot for the current combination
    plot <- plot_nmds(metadata.HNEG.year, habitat, depth)
    
    # Add the plot to the list
    plots_list[[paste(habitat, depth, sep = "_")]] <- plot
  }
}

# Arrange the plots into a grid
ggsave("HNEG_year.pdf",grid.arrange(grobs = plots_list, ncol = 3), width = 12, height = 16)
#NMDS won't give infomrative results.





# PLSDA

#plot_nmds <- function(metadata, habitat, d) {
  # Subset metadata by Habitat_plant and depth, #HNEG
data <- HNEG.pareto[rownames(metadata.HNEG.year),]
plant<-plsda(data, metadata.HNEG.year$Habitat_plant, ncomp = 10)

perf.plant <- perf(plant, validation = 'Mfold', folds = 3, 
                         progressBar = FALSE,  
                         nrepeat = 10)  

plot(perf.plant, sd = TRUE, legend.position = 'horizontal') # 2 is enough
plant<-plsda(data, metadata.HNEG.year$Habitat_plant,ncomp = 2)

plotIndiv(plant, ind.names = F, ellipse = TRUE, legend = TRUE, title = "HNEG_PLSDA", group = metadata.HNEG.year$Habitat_plant, col.per.group = c("#B8912F","#5E813F","#4273B0","orange","blue","purple","green") )
plant.load <- loadings(plant)[["X"]]

depth<-plsda(data, metadata$depth, ncomp = 10)
perf.depth <- perf(depth, validation = 'Mfold', folds = 3, 
                   progressBar = FALSE,  
                   nrepeat = 10)  
plot(perf.depth, sd = TRUE, legend.position = 'horizontal')

#RPPOS
data <- RPPOS.pareto[rownames(metadata.RPPOS.year),]
plant<-plsda(data, metadata.RPPOS.year$Habitat_plant, ncomp = 10)

perf.plant <- perf(plant, validation = 'Mfold', folds = 3, 
                   progressBar = FALSE,  
                   nrepeat = 10)  

plot(perf.plant, sd = TRUE, legend.position = 'horizontal') # 2 is enough
plant<-plsda(data, metadata.RPPOS.year$Habitat_plant, ncomp = 2)

plotIndiv(plant, ind.names = F, ellipse = TRUE, legend = TRUE, title = "RPPOS_PLSDA", group = metadata.RPPOS.year$Habitat_plant, col.per.group = c("#B8912F","#5E813F","#4273B0","orange","blue","purple","green") )
plant.load <- loadings(plant)[["X"]]

depth<-plsda(data, metadata$depth, ncomp = 10)
perf.depth <- perf(depth, validation = 'Mfold', folds = 3, 
                   progressBar = FALSE,  
                   nrepeat = 10)  
plot(perf.depth, sd = TRUE, legend.position = 'horizontal')

write.table(metadata.RPPOS.year[,c(1:5)], "metadata.year.txt", row.names = TRUE)


# 6 seems to be enough
depth<-plsda(data, metadata$depth, ncomp = 6)
loadings(depth)
plotIndiv(depth, c(1,2), ind.names = TRUE, ellipse = TRUE, legend = TRUE)

depth.load <- loadings(depth)[["X"]]
plotIndiv(depth, ind.names = TRUE, ellipse = TRUE, legend = TRUE)
