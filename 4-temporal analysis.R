#temporal analysis. 
library(performance)
library(vegan)
library(lme4)
library(lmerTest)
library(car)
library(tidyverse)
setwd("~/PERMA/R")

#Function ----
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
#important_features <- read.csv("top_5%_features.csv")

plot_avni <- function(feature_names, habitat_col = "habitat") {
  # Check if feature_names exist in the data
  if (length(feature_names) == 0) {
    stop("None of the specified feature names are found in the data.")
  }
  
  df_name <- ifelse(grepl("RPPOS", feature_names), "RPPOS.avni.final", "HNEG.avni.final")
  data <- get(df_name)
  
  # Select the valid features from the data
  selected_data <- data %>%
    select(all_of(feature_names)) %>%
    pivot_longer(cols = everything(), names_to = "Feature", values_to = "Value")
  
  metadata <- data %>%
    select(habitat, depth, Habitat_plant)
  
  # Merge with metadata to include habitat and depth
  plot_data <- cbind(selected_data, metadata)
  plot_data$Habitat_plant <- factor(plot_data$Habitat_plant, 
                                    levels = c("PALSA", "BOG", "AvC_P", "AvCe_C", "AvCe_F1", "FEN"))
  
  # Create the box plot
  ggplot(plot_data, aes(x = Habitat_plant, y = Value, fill = Habitat_plant, shape = depth)) +
    geom_jitter(position = position_jitterdodge(dodge.width = 0.6), aes(color = Habitat_plant), alpha = 0.7) +
    scale_fill_manual(values = habitat_colors) +
    scale_color_manual(values = site_colors) +
    labs(title = feature_names,
         x = "habitat",
         y = "Intensity") +
    theme_minimal()
}

site_colors <-  c("PALSA" = "#B8912F", "BOG" = "#5E813F", "FEN" = "#4273B0", "AvCe_F1"= "purple4", "AvCe_C" = "red4",  "AvC_P" = "red1")
habitat_colors <- c("Palsa" = "#B8912F", "Bog" = "#5E813F", "Fen" = "#4273B0")

plot_feature_temporal <- function(feature_names, habitat_col = "habitat") {
  
  # Check if feature_names exist in the data
  if (length(feature_names) == 0) {
    stop("None of the specified feature names are found in the data.")
  }
  
  df_name <- ifelse(grepl("RPPOS", feature_names), "RPPOS.ac.final", "HNEG.ac.final")
  data <- get(df_name)
  
  # Select the valid features from the data
  selected_data <- data %>%
    select(all_of(feature_names)) 
  
  metadata <- data %>%
    select(Year, habitat, depth) 
  
  # Merge with metadata to include habitat and depth
  plot_data <- cbind(selected_data[rownames(metadata),], metadata)
  colnames(plot_data)[1] <- "Value"
  plot_data$Year <- as.factor(plot_data$Year)
  
  # Load regression results
  regression_results <- get("all_results_linear")
  
  # Create a new column combining habitat and depth
  plot_data <- plot_data %>%
    mutate(Habitat_plant_Depth = paste(toupper(.data[[habitat_col]]), depth, sep = "_"))
  
  plot_data.t <- plot_data %>%
    left_join(regression_results %>% 
                select(Habitat_plant_Depth, all_of(feature_names))%>%
                rename_with(~ "Coefficient", all_of(feature_names))
                )%>% mutate(Label = ifelse(!is.na(Coefficient) & Coefficient != "ns", paste("Coefficient =", Coefficient), NA))
  

  coef_labels <- plot_data.t %>%
    filter(!is.na(Coefficient) & Coefficient != "ns") %>%
    group_by(.data[[habitat_col]], depth) %>%
    summarize(Coefficient = unique(Coefficient), .groups = "drop") %>%
    mutate(Label = paste("Coefficient =", Coefficient))
  
  # Plot
  p <- ggplot(plot_data.t, aes(x = Year, y = Value, fill = .data[[habitat_col]])) +
    geom_point(aes(color = .data[[habitat_col]]), size = 1, alpha = 0.5) +
     geom_smooth(data = subset(plot_data.t, !is.na(Coefficient) & Coefficient != "ns"), 
                 method = "lm", se = TRUE, 
                 aes(group = interaction(.data[[habitat_col]], depth), 
                     color = .data[[habitat_col]])) + 
    scale_fill_manual(values = habitat_colors) +
    scale_color_manual(values = habitat_colors) +
    scale_shape_manual(values = c(16, 17, 18, 15)) +  
    scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash")) +  
    labs(title = paste(feature_names),
         x = "Year",
         y = "Pareto-scaled Intensity",
         shape = "Depth",
         linetype = "Depth") +
    theme_classic() +
    facet_grid(rows = vars(depth), cols = vars(.data[[habitat_col]])) +
    geom_text(data = coef_labels, aes(x = Inf, y = Inf, label = Label), 
              hjust = 1.1, vjust = 1.5, inherit.aes = FALSE, size = 2.1, color = "black")+
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          text = element_text(size = 7))
  
  
  return(p)
}

plot_feature_depth <- function(feature_names, habitat_col = "habitat", depth_col = "depth") {
  # Check if feature_names exist in the data
  if (length(feature_names) == 0) {
    stop("None of the specified feature names are found in the data.")
  }
  
  df_name <- ifelse(grepl("RPPOS", feature_names), "RPPOS.ac.final", "HNEG.ac.final")
  data <- get(df_name)
  
  # Select the valid features from the data
  selected_data <- data %>%
    select(all_of(feature_names)) %>%
    pivot_longer(cols = everything(), names_to = "Feature", values_to = "Value")
  
  metadata <- data %>%
    select(Year, habitat, depth)
  
  # Merge with metadata to include habitat and depth
  plot_data <- cbind(selected_data, metadata)
  
  # Create the box plot
  ggplot(plot_data, aes(x = .data[[depth_col]], y = Value, fill = .data[[habitat_col]])) +
    geom_jitter(position = position_jitterdodge(dodge.width = 0.6), aes(color = .data[[habitat_col]]), alpha = 0.7) +
    geom_smooth(method = "loess", se = TRUE, aes(group = .data[[habitat_col]], color = .data[[habitat_col]])) + 
    scale_fill_manual(values = habitat_colors) +
    scale_color_manual(values = habitat_colors) +
    labs(title = feature_names,
         x = "Year",
         y = "Intensity") +
    theme_minimal()
}



#load data ----
metadata <- read.csv("metadata.all.csv", row.names = 1)
metadata$Sample_code <- gsub(" ", "", metadata$Sample_code)
metadata$habitat <- factor(metadata$habitat, levels = c("Palsa", "Bog", "Fen"))
depth_lookup <- data.frame(depth_category = c("1_5", "10_14", "20_24", "30_34", "40_44", "50_54", "60_64", "70_74", "80_84", "90_94"), 
                           depth_numeric = c(3, 12, 22, 32, 42, 52, 62, 72, 82, 92)
)
metadata$depth_n <- depth_lookup$depth_numeric[match(metadata$depth, depth_lookup$depth_category)]
metadata <- metadata[which(metadata$depth_n < 40),]
metadata.ac <- subset(metadata, grepl("Autochamber",metadata$site))



### annotate temporal features ----
#read HNEG
HNEG <- read.csv("HNEG.csv", row.names = 1, header = T) %>%
  t(.) %>%
  as.data.frame(.) 
#2014_S_2_20_24 was lost during processing for HILIC
new_row_names <- rownames(HNEG) %>%
  gsub("^X", "", .) %>%  # Remove "X" at the beginning of each row name
  gsub("\\.raw\\.\\.F.*", "", .)  # Remove everything after ".raw..F" including ".raw..F"

#cleaning out in-source fragmentation
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

HNEG <- HNEG %>%
  dplyr::select(-all_of(HNEG_remove))
#1076 --> 1046

rownames(HNEG) <- new_row_names
metadata.HNEG.ac <- metadata.ac[rownames(metadata.ac) != "2014_S_2_20_24", ]
HNEG.ac <- HNEG[rownames(metadata.HNEG.ac),]
HNEG.ac.pareto <- na.omit(HNEG.ac) %>%
  log_transform() %>%
  pareto_scale()

HNEG.ac.log <- na.omit(HNEG.ac) %>%
  log_transform() 

#read RPPOS
RPPOS <- read.csv("RPPOS.csv", row.names = 1, header = T) %>%
  t(.) %>%
  as.data.frame(.) 

new_row_names <- rownames(RPPOS) %>%
  gsub("^X", "", .) %>%  # Remove "X" at the beginning of each row name
  gsub("\\.raw\\.\\.F.*", "", .)  # Remove everything after ".raw..F" including ".raw..F"
rownames(RPPOS) <- new_row_names

#cleaning out in-source fragmentation
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

RPPOS_remove <- intersect(RPPOS_remove, colnames(RPPOS))


RPPOS <- RPPOS %>%
  dplyr::select(-all_of(RPPOS_remove))
#4101 --> 3886

metadata.RPPOS.ac <- metadata.ac
RPPOS.ac <- RPPOS[rownames(metadata.RPPOS.ac),]
RPPOS.ac.pareto <- na.omit(RPPOS.ac) %>%
  log_transform() %>%
  pareto_scale()

RPPOS.ac.log <- na.omit(RPPOS.ac) %>%
  log_transform()



###
HNEG.ac.final <- cbind(HNEG.ac.pareto, metadata.HNEG.ac[rownames(HNEG.ac.pareto),])
RPPOS.ac.final <- cbind(RPPOS.ac.pareto, metadata.RPPOS.ac[rownames(RPPOS.ac.pareto),])


metadata.avni <- subset(metadata, Year == "2023" & (depth == "1_5" | depth == "10_14"))%>%
  mutate(Habitat_plant = case_when(
    Habitat_plant == "AvCj_P" ~ "PALSA",
    Habitat_plant == "IncE" ~ "FEN",
    TRUE ~ Habitat_plant  # Keeps other values unchanged
  ))


HNEG.avni <- HNEG[rownames(metadata.avni),]

HNEG.avni.pareto <- na.omit(HNEG.avni) %>%
  log_transform() %>%
  pareto_scale()
HNEG.avni.final <- cbind(HNEG.avni.pareto, metadata.avni[rownames(HNEG.avni.pareto),])


RPPOS.avni <- RPPOS[rownames(metadata.avni),]

RPPOS.avni.pareto <- na.omit(RPPOS.avni) %>%
  log_transform() %>%
  pareto_scale()
RPPOS.avni.final <- cbind(RPPOS.avni.pareto, metadata.avni[rownames(RPPOS.avni.pareto),])


# HNEG.ac.final <- cbind(HNEG.ac.log, metadata.HNEG.ac[rownames(HNEG.ac.log),])
# RPPOS.ac.final <- cbind(RPPOS.ac.log, metadata.RPPOS.ac[rownames(RPPOS.ac.log),])
# 
# 
# HNEG.ac.final <- cbind(HNEG.ac, metadata.HNEG.ac[rownames(HNEG.ac),])
# RPPOS.ac.final <- cbind(RPPOS.ac, metadata.RPPOS.ac[rownames(RPPOS.ac),])

# Function to annotate temporal change for a single feature
annotate_temporal_change <- function(df, feature) {
  results <- data.frame(Habitat_plant = character(), Annotation = character(), stringsAsFactors = FALSE)
  
  # Subset data by unique Habitat_plant
  for (habitat in unique(df$Habitat_plant)) {
    subset_data <- df[df$Habitat_plant == habitat, ]
    
    # Initialize variables
    p_value <- NA
    slope <- NA
    r2 <- NA
    
    # Perform linear model
    tryCatch({
      lm_result <- lm(as.formula(paste(feature, "~ Year")), data = subset_data)
      p_value <- summary(lm_result)$coefficients[2, 4]
      slope <- round(summary(lm_result)$coefficients[2, 1], 3)
      r2 <- summary(lm_result)$r.squared
    }, error = function(e) {
      # Handle errors and keep NA values
      warning(paste("Error in feature:", feature, "for Habitat_plant:", habitat, "-", e$message))
    })
    
    # Annotate based on significance and slope, set p < 0.05 and R2 > 0.3
    if (!is.na(p_value) & p_value < 0.05 & r2 > 0.2) {
      annotation <- slope
    } else {
      annotation <- "ns"
    }
    
    results <- rbind(results, data.frame(Habitat_plant = habitat, Annotation = annotation))
  }
  
  return(results)
}



important_annotated <- read_csv("/Users/benyang/PERMA/R/important_features_module_1224.csv")
list_of_features <- important_annotated$ID
all_results_linear <- data.frame()  # Initialize an empty dataframe to store results

for (feature in list_of_features) {
  # Determine the appropriate dataframe
  df_name <- ifelse(grepl("RPPOS", feature), "RPPOS.ac.final", "HNEG.ac.final")
  df <- get(df_name)
  
  # Check if the feature exists in the dataframe
  if (feature %in% colnames(df)) {
    # Subset the relevant columns
    df_subset <- df[, c("Year", "Habitat_plant", "depth", feature)]
    
    # Initialize a temporary storage for this feature
    feature_results_all <- data.frame()
    
    # Iterate over each combination of Habitat_plant and Depth
    for (habitat in unique(df_subset$Habitat_plant)) {
      for (depth in c("1_5", "10_14", "20_24", "30_34")) {
        subset_data <- df_subset[df_subset$Habitat_plant == habitat & df_subset$depth == depth, ]
        
        if (nrow(subset_data) > 1) {  # Ensure enough data points for regression
          feature_results <- annotate_temporal_change(subset_data, feature)
          feature_results$Depth <- depth  # Add depth information
          feature_results$Feature <- feature  # Track which feature this corresponds to
          
          # Store results
          feature_results_all <- rbind(feature_results_all, feature_results)
        }
      }
    }
    
    # Reshape so that each feature has its own column
    if (nrow(feature_results_all) > 0) {
      feature_results_all <- feature_results_all %>%
        unite("Habitat_plant_Depth", Habitat_plant, Depth, sep = "_") %>%
        pivot_wider(names_from = Feature, values_from = Annotation)
      
      # Merge with main results dataframe
      if (nrow(all_results_linear) == 0) {
        all_results_linear <- feature_results_all
      } else {
        all_results_linear <- full_join(all_results_linear, feature_results_all, by = "Habitat_plant_Depth")
      }
    }
  } else {
    warning(paste("Feature", feature, "not found in dataframe", df_name))
  }
}

regression_stats <- function(feature) {
  # Determine the appropriate dataframe
  df_name <- ifelse(grepl("RPPOS", feature), "RPPOS.ac.final", "HNEG.ac.final")
  df <- get(df_name)
  
  # Check if the feature exists in the dataframe
  if (!(feature %in% colnames(df))) {
    warning(paste("Feature", feature, "not found in dataframe", df_name))
    return(NULL)
  }
  
  # Subset the relevant columns
  df_subset <- df[, c("Year", "Habitat_plant", "depth", feature)]
  
  # Initialize storage for regression results
  regression_results <- data.frame()
  
  # Iterate over each combination of Habitat_plant and Depth
  for (habitat in unique(df_subset$Habitat_plant)) {
    for (depth in c("1_5", "10_14", "20_24", "30_34")) {
      subset_data <- df_subset[df_subset$Habitat_plant == habitat & df_subset$depth == depth, ]
      
      if (nrow(subset_data) > 1) {  # Ensure enough data points for regression
        model <- lm(reformulate("Year", response = feature), data = subset_data)
        
        # Extract statistics
        p_value <- round(summary(model)$coefficients[2, 4], 3)
        slope <- round(summary(model)$coefficients[2, 1], 3)
        r2 <- round(summary(model)$r.squared, 3)  # Use regular R² instead of adjusted R²

# Store results
regression_results <- rbind(regression_results, 
                            data.frame(Habitat_plant = habitat,
                                       Depth = depth,
                                       Feature = feature,
                                       Coefficient = slope,
                                       R2 = r2,
                                       P_value = p_value))
      }
    }
  }
  
  return(regression_results)
}



important_annotated_linear <- all_results_linear%>%
  column_to_rownames("Habitat_plant_Depth")%>%
  t(.)%>%
  as.data.frame(.)%>%
  rownames_to_column("ID")%>%
   left_join(important_annotated,.)
  



regression_stats("FT_445.31565_19.506_RPPOS")
#no info but every thing increased in 10 years.
plot_feature_temporal("FT_370.20113_20.227_RPPOS")

plot_feature_temporal("FT_507.3174_1.293_HNEG")
colnames(all_results_linear)
plot_feature_temporal("FT_335.1753_16.365_RPPOS")
#(+)-Catechin, decreasing should include. 
f5a <- plot_feature_temporal("FT_289.07189_1.99_HNEG")
regression_stats("FT_289.07189_1.99_HNEG")

# both 2 water soluable decrease overtime.
# carnitine
f5b <- plot_feature_temporal("FT_204.1229_1.9_RPPOS")
regression_stats("FT_204.1229_1.9_RPPOS")


regression_stats("FT_329.06691_1.251_HNEG")
regression_stats("FT_335.09277_1.154_HNEG")


f5c <- plot_feature_temporal("FT_309.06037_13.879_RPPOS")
regression_stats("FT_309.06037_13.879_RPPOS")
# Pantothenicacid(VitaminB5)
plot_feature_temporal("FT_220.11782_1.91_RPPOS")

f5d <- plot_feature_temporal("FT_369.13296_16.375_RPPOS")
regression_stats("FT_369.13296_16.375_RPPOS")

plot_feature_temporal("FT_577.13365_11.236_RPPOS")

f5e <-plot_feature_temporal("FT_110.95803_2.204_HNEG")
regression_stats("FT_110.95803_2.204_HNEG")



f5f <-plot_feature_temporal("FT_575.43501_1.161_HNEG")
regression_stats("FT_575.43501_1.161_HNEG")

temporal <- ggarrange(f5a, f5b, f5c, f5d, f5e, f5f,
                      labels = c("A", "B", "C", "D", "E", "F"),
                           ncol = 2, nrow = 3,
                           common.legend = TRUE, legend = "bottom")

# Save the combined plot as a PDF (letter size: 8.5 x 11 inches)
ggsave("temporal.pdf", temporal, width = 8.5, height = 11, units = "in")

#Triterpenoids
#Cucurbitacin F
plot_feature_temporal("FT_517.31714_1.507_HNEG")
#2alpha,19alpha-Dihydroxy-3-oxo-12-ursen-28-oic acid
plot_feature_temporal("FT_485.32734_1.182_HNEG")

plot_feature_temporal("FT_331.08105_15.931_RPPOS")

plot_feature_temporal("FT_414.1604_7.104_RPPOS")
#Icaritin
plot_feature_temporal("FT_369.13296_16.375_RPPOS") # stronger 
plot_feature_temporal("FT_369.13301_19.662_RPPOS")

#Albanin A
plot_feature_temporal("FT_355.11735_15.939_RPPOS")

#if brown means decrease after transition stoped, increase in 10 years means unstable? 
#small sulfur reduced

plot_avni("FT_110.95803_2.204_HNEG")
#big sulfur increased
plot_feature_temporal("FT_575.43499_1.244_HNEG")
plot_avni("FT_575.43499_1.244_HNEG")

plot_feature_temporal("FT_575.43501_1.161_HNEG")
plot_avni("FT_575.43501_1.161_HNEG")

plot_feature_temporal("FT_616.46131_1.339_HNEG")


#brown-bog pos- O-methlated flavonoids
plot_feature_temporal("FT_283.06132_1.56_HNEG")
plot_avni("FT_283.06132_1.56_HNEG")


#bog-pos # cannot be group into any, but keep increase in 10 years.
plot_feature_temporal("FT_398.43535_20.583_RPPOS")
plot_avni("FT_398.43535_20.583_RPPOS")

plot_feature_temporal("FT_335.09277_1.154_HNEG")
plot_avni("FT_335.09277_1.154_HNEG")

plot_feature_temporal("FT_331.08105_15.931_RPPOS")
plot_avni("FT_761.18629_14.007_RPPOS")

plot_avni("FT_349.12942_1.119_HNEG")

plot_avni("FT_469.1294_1.576_HNEG")

write.table(important_annotated_linear,"linear_change_r20.csv", sep = ",")



nosc <- important_LR_final %>%
  #Use OrgMassSpecR package to parse the chemical formula - Returns a named list
  mutate(breakdown = purrr::map(Formula, function(x) data.frame(OrgMassSpecR::ListFormula(x)))) %>%
  #Unnest that list into multiple columns, each named after the element
  unnest('breakdown') %>%
  #Use those columns to calculate the NOSC value
  mutate(NOSC = 4 - ((4*C + H - 2*O - 3*N - 2*S + 5*P)/C)) 

important_LR_final$nosc <- nosc$NOSC

write.csv(important_LR_final, "important_features_LR_addinc_nofilt.csv", row.names = FALSE)






plot_feature_temporal("FT_289.07189_1.99_HNEG")

plot_feature_temporal("FT_455.14857_11.4_RPPOS")
plot_feature_temporal("FT_453.13448_1.656_HNEG")
plot_feature_temporal("FT_359.09266_1.656_HNEG")
plot_feature_temporal("FT_679.19677_1.635_HNEG")
plot_feature_temporal("FT_677.18109_1.448_HNEG")
plot_feature_temporal("FT_693.17603_1.632_HNEG")
plot_feature_temporal("FT_677.18108_1.506_HNEG")

plot_feature_depth("FT_455.14857_11.4_RPPOS")



#benzone
benzene <- c("FT_200.07048_13.022_RPPOS",
"FT_259.09636_12.413_RPPOS",
"FT_281.12832_15.836_RPPOS",
"FT_247.05993_16.263_RPPOS",
"FT_274.1072_12.385_RPPOS",
"FT_469.33097_17.922_RPPOS",
"FT_331.09229_15.037_RPPOS",
"FT_179.03494_1.18_HNEG",
"FT_575.43501_1.161_HNEG")

lapply(benzene, plot_feature_temporal)

benzopyrans <- c("FT_371.11237_14.944_RPPOS",
"FT_271.09635_14.147_RPPOS",
"FT_301.07044_17.529_RPPOS",
"FT_259.05998_11.626_RPPOS",
"FT_331.11748_16.397_RPPOS",
"FT_261.11197_16.485_RPPOS")
lapply(benzopyrans, plot_feature_temporal)


plot_feature_temporal("FT_375.08766_1.597_HNEG")

plot_feature_temporal("FT_179.03494_1.18_HNEG")
plot_feature_temporal("FT_369.13303_16.153_RPPOS")
plot_feature_temporal("FT_435.32545_19.845_RPPOS")

plot_feature_temporal("FT_499.15522_11.689_RPPOS")

plot_feature_temporal("FT_411.08765_1.046_HNEG")
#Avni plots






FT_455.14857_11.4_RPPOS
FT_453.13448_1.656_HNEG
FT_679.19677_1.635_HNEG
FT_677.18109_1.448_HNEG
FT_693.17603_1.632_HNEG
FT_677.18108_1.506_HNEG

plot_avni("FT_455.14857_11.4_RPPOS")
plot_avni("FT_453.13448_1.656_HNEG")
plot_avni("FT_679.19677_1.635_HNEG")
plot_avni("FT_677.18109_1.448_HNEG")
plot_avni("FT_693.17603_1.632_HNEG")
plot_avni("FT_677.18108_1.506_HNEG")




plot_feature_temporal("FT_289.07189_1.99_HNEG")
plot_feature_temporal(list_of_features[1])
plot_feature_temporal(list_of_features[5])
plot_feature_temporal("FT_167.1429_20.919_RPPOS")
plot_feature_temporal("FT_418.22211_16.867_RPPOS")

plot_feature_temporal("FT_451.12308_10.989_RPPOS")

plot_feature_temporal("FT_905.25863_13.872_RPPOS")
# such weird pattern!!!!


plot_feature_temporal("FT_329.10169_16.673_RPPOS")
plot_feature_temporal("FT_295.13273_14.791_RPPOS")

poly.temporal.result <- poly.temporal.result %>%
  mutate(across(everything(), as.numeric))
poly.temporal.result[is.na(poly.temporal.result)] <- 0
linear.temporal.result <- linear.temporal.result %>%
  mutate(across(everything(), as.numeric))
linear.temporal.result[is.na(linear.temporal.result)] <- 0

plot_feature_temporal("shannon_RPPOS")
plot_feature_temporal("shannon_HNEG")

kruskal.test(metadata.HNEG.ac$shannon_HNEG,metadata.HNEG.ac$depth) #depth not sig
kruskal.test(metadata.RPPOS.ac$shannon_RPPOS,metadata.RPPOS.ac$depth)


kruskal.test(metadata.HNEG.ac$shannon_HNEG,metadata.HNEG.ac$Habitat_plant) #Habitat sig
kruskal.test(metadata.RPPOS.ac$shannon_RPPOS,metadata.RPPOS.ac$Habitat_plant) #Habitat sig

#add module infor

plot_feature_temporal("FT_110.95803_2.204_HNEG")
plot_feature_temporal("FT_616.46131_1.339_HNEG")

plot_feature_temporal("FT_204.1229_1.9_RPPOS")
plot_feature_temporal("FT_289.07189_1.99_HNEG")
plot_feature_temporal("FT_447.15098_4.578_HNEG")


plot_feature_temporal("FT_289.07189_1.99_HNEG")
plot_feature_temporal("FT_297.13443_1.673_HNEG")

plot_feature_temporal("FT_329.06691_1.251_HNEG")
plot_feature_temporal("FT_437.08805_1.673_HNEG")
plot_feature_temporal("FT_577.13365_11.236_RPPOS")
plot_feature_temporal("FT_110.95803_2.204_HNEG")

plot_avni("FT_110.95803_2.204_HNEG")

plot_feature_temporal("FT_309.06037_13.879_RPPOS")

plot_feature_temporal("FT_163.08649_13.817_RPPOS")
plot_feature_temporal("FT_335.1753_16.365_RPPOS")
plot_feature_temporal("FT_575.43499_1.244_HNEG")
plot_feature_temporal("FT_355.11735_15.939_RPPOS")
plot_feature_temporal("FT_355.1189_1.32_HNEG")

plot_feature_temporal("FT_327.08617_19.178_RPPOS")

plot_feature_temporal("FT_355.11735_15.939_RPPOS")

plot_feature_temporal("FT_204.1229_1.9_RPPOS")
plot_feature_temporal("FT_289.07189_1.99_HNEG")

icartin <- c("FT_369.13296_16.375_RPPOS",
             "FT_369.13301_19.662_RPPOS",
             "FT_369.13303_16.153_RPPOS",
             "FT_369.13307_18.594_RPPOS")

lapply(icartin, plot_feature_temporal)


results_list <- list()

# Define a function to perform regression analysis for each dataframe
run_regression_analysis <- function(df) {
  # Filter only columns starting with "FT_" (metabolites) and required metadata
  metabolite_cols <- grep("^FT_", colnames(df), value = TRUE)
  habitats <- c("PALSA", "BOG", "FEN")
  
  # Initialize an empty dataframe to store results
  regression_results <- data.frame(
    Feature = character(),
    Habitat_plant = character(),
    P_value = numeric(),
    R2 = numeric(),
    Slope = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Loop over each habitat
  for (habitat in habitats) {
    habitat_data <- df %>% filter(Habitat_plant == !!habitat)
    # Loop over each metabolite column
    for (feature in metabolite_cols) {
      # Run linear regression
      model <- lm(habitat_data[[feature]] ~ habitat_data$Year)
      
      # Extract p-value, R², and slope
      p_value <- summary(model)$coefficients[2, 4]   # p-value for slope
      r_squared <- summary(model)$r.squared          # R² value
      slope <- coef(model)[2]                        # slope coefficient
      
      # Store results in the dataframe
      regression_results <- rbind(regression_results, data.frame(
        Feature = feature,
        Habitat_plant = habitat,
        P_value = p_value,
        R2 = r_squared,
        Slope = slope
      ))
    }
  }
  
  return(regression_results)
}

results_list$HNEG <- run_regression_analysis(HNEG.ac.final)
results_list$RPPOS <- run_regression_analysis(RPPOS.ac.final)

results_combined <- bind_rows(results_list)%>%
  mutate(P_value_corrected = p.adjust(P_value, method = "BH"))%>%
  mutate(Significance = ifelse(P_value_corrected > 0.05, "ns", P_value_corrected))

results_wide <- results_combined %>%
  select(Feature,Habitat_plant, Significance, R2, Slope)%>%
  pivot_wider(
    names_from = Habitat_plant, 
    values_from = c(R2, Slope, Significance),
    names_glue = "{Habitat_plant}_{.value}"
  )

colnames(results_wide)[1] <- "ID"

important.module <- read.csv("important_features_module.csv", row.names = 1)%>%
  left_join(.,results_wide)

write.csv(important.module,"important.module.linear.temporal.csv")

