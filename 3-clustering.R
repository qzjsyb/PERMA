library(tidyverse)

library(ggtree)
library(ggplot2)
library(ggnewscale)
library(ape)
library(scales)

#pattern clustering
#function
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

#loading
setwd("~/PERMA/R")
metadata <- read.csv("metadata.all.csv", row.names = 1)
metadata$habitat <- factor(metadata$habitat, levels = c("Palsa", "Bog", "Fen"))
depth_lookup <- data.frame(depth_category = c("1_5", "10_14", "20_24", "30_34", "40_44", "50_54", "60_64", "70_74", "80_84", "90_94"), 
                           depth_numeric = c(3, 12, 22, 32, 42, 52, 62, 72, 82, 92)
)
metadata$depth_n <- depth_lookup$depth_numeric[match(metadata$depth, depth_lookup$depth_category)]
metadata <- metadata[which(metadata$depth_n < 40),]

metadata$core <- paste0(metadata$Year,metadata$Habitat_plant,metadata$Sample_code)
#RPPOS----
RPPOS <- read.csv("RPPOS.csv", row.names = 1, header = T) %>%
  t(.) %>%
  as.data.frame(.) 


new_row_names <- rownames(RPPOS) %>%
  gsub("^X", "", .) %>%  # Remove "X" at the beginning of each row name
  gsub("\\.raw\\.\\.F.*", "", .)  # Remove everything after ".raw..F" including ".raw..F"
rownames(RPPOS) <- new_row_names

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




#HNEG----
HNEG <- read.csv("HNEG.csv", row.names = 1, header = T) %>%
  t(.) %>%
  as.data.frame(.) 
#2014_S_2_20_24 was lost during processing for HILIC
new_row_names <- rownames(HNEG) %>%
  gsub("^X", "", .) %>%  # Remove "X" at the beginning of each row name
  gsub("\\.raw\\.\\.F.*", "", .)  # Remove everything after ".raw..F" including ".raw..F"
rownames(HNEG) <- new_row_names

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



important_features<-read.csv("important_features_annotated_1015.csv", header = T) %>%
  mutate(across(where(is.character), ~na_if(., "")))
important_features <- important_features[,c(1:15)]
important_features <- important_features %>%
  mutate(Group = case_when(
    Palsa > 0 ~ "Palsa_pos",
    Bog > 0 ~ "Bog_pos",
    Fen > 0 ~ "Fen_pos"
  ))
important_features <- important_features[!is.na(important_features$Group),]

all_features<-read.csv("compound_link.csv", header = T) %>%
  mutate(across(where(is.character), ~na_if(., "")))

###try with removing all NAs

metadata.23 <- metadata %>%
  filter(depth %in% c("1_5", "10_14"))%>%
  filter(Year %in% c("2023"))%>%
  mutate(Habitat_plant = case_when(
    Habitat_plant == "AvCj_P" ~ "PALSA",
    Habitat_plant == "IncE" ~ "FEN",
    Habitat_plant == "IncS" ~ "BOG",
    TRUE ~ Habitat_plant  # Keep other values unchanged
  ))

data <- cbind(RPPOS[rownames(metadata.23),],HNEG[rownames(metadata.23),])
#data <- data %>%
#  select(all_of(important_features$ID))

data.log <- na.omit(data) %>%
  log_transform()

data.pareto <- na.omit(data) %>%
  log_transform() %>%
  pareto_scale()

library(WGCNA)  # Load the package
allowWGCNAThreads() 
input <- data.pareto

# Choose a soft-thresholding power
powers <- c(1:20)  # Set of powers to test
sft <- pickSoftThreshold(input, powerVector = powers, verbose = 5)

# Plot the scale-free topology fit index as a function of the power
par(mfrow = c(1, 2))  # Set up a 1x2 plot layout

powers <- sft$fitIndices[,1]  # Power values tested
scale_free_fit <- sft$fitIndices[,2]  # Scale-free topology fit (R²)
mean_connectivity <- sft$fitIndices[,5]  # Mean connectivity (k)

plot(powers, scale_free_fit)
text(powers, scale_free_fit, 
     labels = powers, cex = 0.9, col = "red")

# Plot 2: Mean connectivity as a function of the power
plot(powers, mean_connectivity)
text(powers, mean_connectivity, 
     labels = powers, cex = 0.9, col = "red")


softPower <- 7  # Example of chosen power based on previous step, R2 = 0.7710987198, k = 0.241016849
temp_cor <- cor  
cor <- WGCNA::cor 


netwk <- blockwiseModules(input, power = softPower, networkType = "signed", 
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          # detectCutHeight = 0.75,
                          minModuleSize = 30,
                          maxBlockSize = 4000,
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          saveTOMs = T,
                          saveTOMFileBase = "ER",
                          numericLabels = T,
                          verbose = 3)

column_indices <- netwk$blockGenes[[1]]
dendrogram_order <- netwk$dendrograms[[1]]$order

# Map tip labels to input column names
tip_to_column <- colnames(input)[column_indices[dendrogram_order]]


# Convert labels to colors for plotting
mergedColors = labels2colors(netwk$colors)
length(unique(mergedColors))
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )

# check modules
MEs <- moduleEigengenes(input, mergedColors)$eigengenes

#### small subsection do correlation----
moi <- read.csv("soil_moi_all.csv", sep = ",")
colnames(moi)[colnames(moi) == "Moisture_contend..."] <- "Moi"
moi$Sample_code_by_core <- gsub("_", "", moi$Sample_code_by_core)
metadata$Sample_code_by_core <- gsub("_", "", rownames(metadata))

moi_metadata <- left_join(moi, metadata)

MEs <- MEs%>%rownames_to_column("Sample_code_by_core")
MEs$Sample_code_by_core <- gsub("_", "", MEs$Sample_code_by_core)

# Left join with ME_moi
ME_moi <- left_join(MEs, moi_metadata, by = "Sample_code_by_core")

site_color <- c("#5E813F","#4273B0","#4273B0","purple4","red4","red1", "#B8912F", "#B8912F")
names(site_color) <-  c("BOG","FEN","IncE", "AvCe_F1", "AvCe_C", "AvC_P", "PALSA", "AvCj_P")
# Remove all underscores from Sample_code_by_core in metadata


cor.test(ME_moi$MEyellow, ME_moi$Tsoil, method = "pearson")
cor.test(ME_moi$MEbrown, ME_moi$Moi, method = "pearson")

library(lme4)
library(lmerTest)
predictors <- c("Tsoil", "pH", "Moi")
responses <- c("MEyellow", "MEturquoise", "MEbrown", "MEpink", "MEgreen")

# Fit models and store results
results <- data.frame(Response = character(), Predictor = character(),
                      Coefficient = numeric(), p_value = numeric(), Sig = character(), stringsAsFactors = FALSE)

for (x in predictors) {
  for (y in responses) {
    
    formula <- as.formula(paste(y, "~", x, "+ (1 | habitat)"))
    model <- lmerTest::lmer(formula, data = ME_moi)
    
    # Extract coefficient and p-value
    model_summary <- summary(model)
    coef_value <- model_summary$coefficients[2, 1]  # Fixed effect estimate
    p_value <- model_summary$coefficients[2, 5]     # p-value from lmerTest
    
    # Determine significance level
    sig_label <- ifelse(p_value > 0.1, "ns",
                        ifelse(p_value > 0.05, ".",
                               ifelse(p_value > 0.01, "*",
                                      ifelse(p_value > 0.001, "**", "***"))))
    
    # Store results (Include sig_label)
    results <- rbind(results, data.frame(Response = y, Predictor = x,
                                         Coefficient = coef_value, p_value = p_value, Sig = sig_label))
  }
}

# Print results
print(results)


p1 <- ggplot(ME_moi, aes(x = Tsoil, y = MEyellow, color = Habitat_plant)) +
  geom_point(alpha = 0.7, size = 3) +  # Scatter points
  geom_smooth(method = "lm", color = "yellow", se = TRUE) +  # Regression line
  scale_color_manual(values = site_color) +  # Custom colors
  labs(title = "MEyellow vs Tsoil", x = "Tsoil", y = "MEyellow") +
  theme_minimal()

ggplot(ME_moi, aes(x = Tsoil, y = MEbrown, color = Habitat_plant)) +
  geom_point(alpha = 0.7, size = 3) +  # Scatter points
  geom_smooth(method = "lm", color = "yellow", se = TRUE) +  # Regression line
  scale_color_manual(values = site_color) +  # Custom colors
  labs(title = "MEyellow vs Tsoil", x = "Tsoil", y = "MEyellow") +
  theme_minimal()

predictor_sets <- list(c("Moi", "pH"), c("Moi", "Tsoil"))
responses <- c("MEyellow", "MEturquoise", "MEbrown", "MEpink", "MEgreen")

# Fit models and store results
results <- data.frame(Response = character(), Predictors = character(),
                      Coefficient_Moi = numeric(), Coefficient_Other = numeric(), Coefficient_Interaction = numeric(),
                      p_value_Moi = numeric(), p_value_Other = numeric(), p_value_Interaction = numeric(),
                      Sig_Moi = character(), Sig_Other = character(), Sig_Interaction = character(),
                      stringsAsFactors = FALSE)

for (predictors in predictor_sets) {
  x1 <- predictors[1]
  x2 <- predictors[2]
  
  for (y in responses) {
    formula <- as.formula(paste(y, "~", x1, "*", x2))  # Include interaction
    model <- lm(formula, data = ME_moi)
    
    # Extract coefficient and p-values
    model_summary <- summary(model)
    coef_values <- model_summary$coefficients[2:4, 1]  # Fixed effect estimates
    p_values <- model_summary$coefficients[2:4, 4]     # p-values
    
    # Determine significance levels
    sig_labels <- sapply(p_values, function(p) {
      if (p > 0.1) "ns" else if (p > 0.05) "." else if (p > 0.01) "*" else if (p > 0.001) "**" else "***"
    })
    
    # Store results
    results <- rbind(results, data.frame(Response = y, Predictors = paste(x1, x2, sep = " * "),
                                         Coefficient_Moi = coef_values[1], Coefficient_Other = coef_values[2], Coefficient_Interaction = coef_values[3],
                                         p_value_Moi = p_values[1], p_value_Other = p_values[2], p_value_Interaction = p_values[3],
                                         Sig_Moi = sig_labels[1], Sig_Other = sig_labels[2], Sig_Interaction = sig_labels[3]))
  }
}

# Print results
print(results)

ggplot(ME_moi, aes(x = Moi, y = MEgreen, color = factor(pH))) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Interaction Effect of Moi and pH on MEyellow",
       x = "Moi", y = "MEyellow", color = "pH Level") +
  theme_minimal()



# Plot MEbrown vs Moi, colored by Habitat_plant
p2 <- ggplot(ME_moi, aes(x = Moi, y = MEbrown, color = Habitat_plant)) +
  geom_point(alpha = 0.7, size = 3) +
  geom_smooth(method = "lm", color = "brown", se = TRUE) +
  scale_color_manual(values = site_color) +
  labs(title = "MEbrown vs Moi", x = "Moi", y = "MEbrown") +
  theme_minimal()


# Extract list of features in each module-----
moduleFeatures <- split(colnames(input), mergedColors)
module_df <- stack(moduleFeatures)
colnames(module_df) <- c("ID", "module")

module_df <- module_df[match(tip_to_column, module_df$ID), ]

important_features_module <- left_join(important_features, module_df)
all_features_module <- left_join(all_features, module_df, by = c("CompoundID" = "ID")) %>%
  left_join(., important_features %>% select(ID, Group), by = c("CompoundID" = "ID"))

all_features_module <- all_features_module %>%
  filter(!CompoundID %in% c(HNEG_remove, RPPOS_remove))


### let me exact test see if any superclass is enriched in any module: ----
perform_fisher_test <- function(all_features_module, rare_num = 5) {
  # Generate count table for each module
  module_count_all <- as.data.frame(table(all_features_module$module)) %>%
    filter(Freq > 0)
  
  # Get total counts of each Superclass in the full dataset
  total_superclass_count <- as.data.frame(table(all_features_module$Superclass))
  colnames(total_superclass_count) <- c("Superclass", "total_count")
  
  # Initialize list to store results for each module
  module_results <- list()
  
  # Loop through each module
  for (mod in module_count_all$Var1) {
    # Filter data for the current module
    module_data <- all_features_module %>% filter(module == mod)
    
    # Get superclass count in this module
    superclass_count_module <- as.data.frame(table(module_data$Superclass))
    colnames(superclass_count_module) <- c("Superclass", "count_in_module")
    
    # Merge with total superclass count
    merged_data <- merge(superclass_count_module, total_superclass_count, by = "Superclass", all.x = TRUE)
    
    # Initialize vectors to store Fisher test results
    p_values <- numeric(nrow(merged_data))
    odds_ratios <- numeric(nrow(merged_data))
    
    # Loop through each Superclass
    for (i in 1:nrow(merged_data)) {
      class_in_module <- merged_data$count_in_module[i]
      non_class_in_module <- module_count_all$Freq[module_count_all$Var1 == mod] - class_in_module
      class_in_dataset <- merged_data$total_count[i] - class_in_module
      non_class_in_dataset <- sum(total_superclass_count$total_count) - sum(superclass_count_module$count_in_module) - class_in_dataset
      
      # Construct contingency table
      contingency_table <- matrix(c(class_in_module, non_class_in_module,
                                    class_in_dataset, non_class_in_dataset),
                                  nrow = 2, byrow = TRUE)
      
      # Perform Fisher's Exact Test
      fisher_test_result <- fisher.test(contingency_table)
      
      # Store results
      p_values[i] <- fisher_test_result$p.value
      odds_ratios[i] <- fisher_test_result$estimate
    }
    
    # Add results to dataframe
    merged_data$p_value <- p_values
    merged_data$p_value_adj <- p.adjust(p_values, method = "BH")
    merged_data$odds_ratio <- odds_ratios
    
    # Store in list
    module_results[[mod]] <- merged_data
  }
  
  return(module_results)
}

# Run the function
final_results <- perform_fisher_test(all_features_module, rare_num = 5)
filtered_results <- lapply(final_results, function(df) {
  df %>% filter(p_value_adj <= 0.05)
})


#for S-containing:
all_features_module$S_containing <- grepl(" S", all_features_module$Formula)

# Count S-containing in the Brown module
brown_module <- all_features_module %>% filter(module == "brown")
S_in_brown <- sum(brown_module$S_containing)
non_S_in_brown <- nrow(brown_module) - S_in_brown

# Count S-containing in the entire dataset
S_in_dataset <- sum(all_features_module$S_containing)
non_S_in_dataset <- nrow(all_features_module) - S_in_dataset

# Construct contingency table
contingency_table <- matrix(c(S_in_brown, non_S_in_brown,
                              S_in_dataset - S_in_brown, non_S_in_dataset - non_S_in_brown),
                            nrow = 2, byrow = TRUE)

# Perform Fisher's Exact Test
fisher_test_result <- fisher.test(contingency_table)

fisher_test_result

####----

module_counts <- table(important_features_module$Group, important_features_module$module)
module_counts

plot_df <- cbind(MEs, metadata.23[rownames(MEs),])
plot_df$Habitat_plant <- factor(plot_df$Habitat_plant, levels = c("PALSA", "BOG", "AvC_P", "AvCe_C","AvCe_F1", "FEN"))

plot_df$Habitat_plant <- fct_recode(
  plot_df$Habitat_plant,
  RT = "AvC_P",
  ET = "AvCe_C",
  PT = "AvCe_F1"
)


site_color <- c("#4273B0","purple4","red4","red1", "#B8912F","#B8912F","#4273B0", "#5E813F", "#5E813F")
names(site_color) <-  c("IncE", "PT", "ET", "RT", "AvCj_P","PALSA", "FEN", "BOG", "IncS")


plot_MEcolor <- function(color_input) {
  # Create the column name dynamically based on input color
  column_name <- paste0("ME", color_input)
  
  # Ensure that the column exists in plot_df
  if (!column_name %in% colnames(plot_df)) {
    stop(paste("Column", column_name, "not found in plot_df"))
  }
  
  # Plot using the dynamic color and column
  ggplot(plot_df, aes(x = Habitat_plant, y = .data[[column_name]], group = 1)) +
    stat_summary(fun = mean, geom = "line", size = 1.5, color = color_input) +  # Line plot in input color
    stat_summary(fun = mean, geom = "point", size = 4, color = color_input) +  # Points in input color
    stat_summary(fun.data = mean_cl_normal, geom = "ribbon", alpha = 0.1, fill = color_input) +  # Confidence region in input color
    labs(title = paste("ME", color_input), x = "Habitat", y = paste("Mean ME", color_input, "Value")) +
    theme_classic() +
    theme(legend.position = "none",
          text = element_text(size = 8))
}



plot_MEcolor("brown")
plot_MEcolor("green") 
plot_MEcolor("pink")

plot_MEcolor("turquoise")

plot_MEcolor("yellow")
plot_MEcolor("purple")
plot_MEcolor("black")



ggsave("brown.pdf", plot_MEcolor("brown"),width = 2.5, height = 2.5 )
ggsave("green.pdf", plot_MEcolor("green"), width = 2.5, height = 2.5)
ggsave("pink.pdf",plot_MEcolor("pink"), width = 2.5, height = 2.5)
ggsave("yellow.pdf",plot_MEcolor("yellow"), width = 2.5, height = 2.5)
ggsave("turquoise.pdf",plot_MEcolor("turquoise"), width = 2.5, height = 2.5)


# Match the Group colors to the data columns
feature_group_colors <- group_colors[important_features$Group[match(colnames(input), important_features$ID)]]


plotDendroAndColors(netwk$dendrograms[[1]], feature_group_colors[netwk$blockGenes[[1]]], "feature_group_colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
library(RColorBrewer)
all_features <- all_features %>%
  mutate(Superclass_merged = case_when(
    Superclass %in% c("Alkaloids and derivatives",
                      "Homogeneous non-metal compounds",
                      "Hydrocarbons",
                      "Lignans, neolignans and related compounds",
                      "Mixed metal/non-metal compounds",
                      "Nucleosides, nucleotides, and analogues",
                      "Organic 1,3-dipolar compounds",
                      "Organophosphorus compounds",
                      "Organosulfur compounds") ~ "Other",
    TRUE ~ Superclass  # Keep the other Superclasses unchanged
  ))
unique_superclasses <- unique(all_features$Superclass_merged)
superclass_colors <- brewer.pal(length(na.omit(unique_superclasses)), "Set1")
names(superclass_colors) <- na.omit(unique_superclasses)
feature_superclass_colors <- superclass_colors[all_features$Superclass_merged[match(colnames(input), all_features$CompoundID)]]


plotDendroAndColors(netwk$dendrograms[[1]], feature_superclass_colors[netwk$blockGenes[[1]]], "superclass",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)




tree_data <- data.frame(
  CompoundID = colnames(input)[netwk$blockGenes[[1]]]) %>%
  left_join(.,all_features_module%>%select(CompoundID, Superclass_merged, Group, module))

circ <- ggtree(netwk$dendrograms[[1]], 
               layout = "fan", 
               open.angle = 30, 
               xlim = c(-1, -0.5)) +  # Manually set x-axis limits
  theme_minimal()

tip_data <- circ$data[circ$data$isTip, ]%>%
  mutate(Group = tree_data$Group) %>%
  filter(!is.na(Group))


circ <- circ + geom_tippoint(aes(color = Group), data = tip_data, size = 3) +
  scale_color_manual(values = group_colors, name = "group")

p1 <- gheatmap(circ, tree_data["module"], 
           width = 0.1,
           colnames_angle = 90, 
           colnames_offset_y = 0.25) +
  scale_fill_identity(mergedColors)+
  theme(legend.position = "right")

library(ggnewscale)
p2 <- p1 + new_scale_fill() 
p2<-  gheatmap(p2, tree_data["Superclass_merged"], offset = 0.05, width = 0.1,
               colnames_angle = 90, colnames_offset_y = 0.25) +
  scale_fill_manual(values = superclass_colors, 
                    name = "Superclass_merged")+
  theme(legend.position = "right",
        text = element_text(size = 8))


ggsave("tree.pdf", p2, width = 15, heigh = 10, )


# export module table
nosc <- important_features%>%
  #Use OrgMassSpecR package to parse the chemical formula - Returns a named list
  mutate(breakdown = purrr::map(Formula, function(x) data.frame(OrgMassSpecR::ListFormula(x)))) %>%
  #Unnest that list into multiple columns, each named after the element
  unnest('breakdown') %>%
  #Use those columns to calculate the NOSC value
  mutate(NOSC = 4 - ((4*C + H - 2*O - 3*N - 2*S + 5*P)/C)) 
important_features$nosc <- nosc$NOSC

important_features %>%
  filter(nosc > -1.6 & nosc < 1.2) %>%
  ggplot(aes(x = module, y = nosc, fill = module, color = module)) +
  
  geom_violin(alpha = 0.5) +  # Violin plot for group-wise distribution
  
  geom_jitter(position = position_jitter(width = 0.2), size = 1.5, alpha = 0.7) +
  
  # Add horizontal line at y = -0.3
  geom_hline(yintercept = 0, linetype = "solid", color = "darkgreen") +
  annotate("text", x = 3.5, y = 0.05, label = "Average for Plant Metabolites", color = "darkgreen", hjust = 0) +
  
  geom_hline(yintercept = -0.3, linetype = "solid", color = "black") +
  annotate("text", x = 3.5, y = -0.25, label = "Average for Soil Organic Matter", color = "black", hjust = 0) +
  
  # Add titles and labels
  labs(title = "Violin Plot of NOSC across Groups of Predictive Features",
       x = "Group",
       y = "NOSC") +
  
  # Apply minimal theme and adjust legend
  theme_linedraw() +
  theme(legend.position = "bottom", legend.direction = "horizontal") 
  # geom_signif(comparisons = list(c("Palsa_pos", "Overall"), 
  #                                c("Bog_pos", "Overall"),
  #                                c("Fen_pos", "Overall")),
  #             annotations = annotations,  # Annotate with D and p-values
  #             y_position = c(1.6, 1.4, 1.2, 1),
  #             color = "black",# Adjust y positions for the lines
  #             tip_length = 0.1)  # Add some space for the lines


write.csv(important_features_module, "important_features_module_1224.csv")


# Example usage:
plot_MEcolor(unique(dynamicColors)[1]) 
plot_MEcolor(unique(dynamicColors)[2])
plot_MEcolor(unique(dynamicColors)[3])
plot_MEcolor(unique(dynamicColors)[4])
plot_MEcolor(unique(dynamicColors)[5])
plot_MEcolor(unique(dynamicColors)[6])


ggsave("tur.pdf", plot_MEcolor(unique(dynamicColors)[3]),width = 3.5, height = 3 )
ggsave("blue.pdf", plot_MEcolor(unique(dynamicColors)[2]), width = 3.5, height = 3)
ggsave("brown.pdf",plot_MEcolor(unique(dynamicColors)[1]), width = 3.5, height = 3)
ggsave("yellow.pdf",plot_MEcolor(unique(dynamicColors)[5]), width = 3.5, height = 3)


#optional, hypothetical temporal change plot
plot_avni_change <- function(feature_name) {
  
  site_color <- c("#4273B0","purple4","red4","red1", "#B8912F","#B8912F","#4273B0", "#5E813F", "#5E813F")
  names(site_color) <-  c("IncE", "AvCe_F1", "AvCe_C", "AvC_P", "AvCj_P","PALSA", "FEN", "BOG", "IncS")
  
 
  # Check if feature_names exist in the data
  if (grepl("_HNEG", feature_name)) {
    data_to_use <- HNEG   # Assuming HNEG.all is your dataset for HNEG features
    metadata_to_use <- metadata.23  # Assuming metadata.HNEG is your metadata for HNEG features
  } else if (grepl("_RPPOS", feature_name)) {
    data_to_use <- RPPOS # Assuming RPPOS.all is your dataset for RPPOS features
    metadata_to_use <- metadata.23  # Assuming metadata.RPPOS is your metadata for RPPOS features
  }
  
  # Merge with metadata to include habitat and depth
  plot_data <- cbind(metadata_to_use, data_to_use[rownames(metadata_to_use),])
  
  plot_data <- plot_data %>%
    dplyr::rename(feature_value = feature_name)
  
  plot_data$Habitat_plant <- factor(plot_data$Habitat_plant, levels = c("PALSA", "BOG", "AvC_P", "AvCe_C","AvCe_F1", "FEN"))
  avg_data <- aggregate(feature_value ~ Habitat_plant, data = plot_data, FUN = mean)
  
  # Create the box plot
  ggplot(plot_data, aes(x = Habitat_plant, y = feature_value, fill = Habitat_plant)) +
    geom_jitter(position = position_jitterdodge(dodge.width = 0.6), aes(color = Habitat_plant, alpha = 0.3),size = 1) +
    stat_summary(fun.data = mean_cl_normal, geom = "point", size = 4,position = position_dodge(width = 0.6), aes(color = Habitat_plant)) +
    stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.2,position = position_dodge(width = 0.6), aes(color = Habitat_plant)) +
    scale_fill_manual(values = site_color) +
    scale_color_manual(values = site_color) +
    labs(title = feature_name,
         x = "Habitat",
         y = "LC-MS/MS Intensity") +
    theme_minimal()+
    theme(legend.position = "none")+
    scale_y_log10()
}

lapply(moduleFeatures$blue,plot_avni_change) #increased p - b - transit - f.
lapply(moduleFeatures$brown,plot_avni_change)
lapply(moduleFeatures$turquoise,plot_avni_change)

plot_avni_change("FT_285.07561_15.408_RPPOS")

NAS <- c("FT_679.19677_1.635_HNEG",
         "FT_283.05996_15.736_RPPOS",
         "FT_315.04982_13.659_RPPOS",
         "FT_319.03663_15.992_RPPOS",
         "FT_349.12942_1.119_HNEG",
         "FT_407.07589_15.386_RPPOS",
         "FT_197.00916_1.942_HNEG",
         "FT_221.0443_11.687_RPPOS",
         "FT_363.07096_16.475_RPPOS",
         "FT_275.05485_16.884_RPPOS",
         "FT_309.06037_13.879_RPPOS",
         "FT_361.17967_17.267_RPPOS",
         "FT_389.13826_16.067_RPPOS",
         "FT_227.09253_1.671_HNEG",
         "FT_474.24829_17.609_RPPOS",
         "FT_499.15522_11.689_RPPOS",
         "FT_451.10217_15.269_RPPOS",
         "FT_297.13443_1.673_HNEG",
         "FT_435.07236_1.672_HNEG",
         "FT_403.08251_1.265_HNEG",
         "FT_301.07044_17.529_RPPOS",
         "FT_329.06691_1.251_HNEG",
         "FT_331.08105_15.931_RPPOS",
         "FT_391.1538_16.706_RPPOS",
         "FT_543.12838_16.806_RPPOS",
         "FT_453.13448_1.656_HNEG",
         "FT_455.14849_11.686_RPPOS",
         "FT_283.05986_18.535_RPPOS",
         "FT_335.09277_1.154_HNEG",
         "FT_349.10682_11.686_RPPOS",
         "FT_279.10138_17.555_RPPOS",
         "FT_469.1294_1.576_HNEG",
         "FT_202.97528_1.425_HNEG",
         "FT_236.93633_1.243_HNEG",
         "FT_269.08071_14.345_RPPOS",
         "FT_273.04054_2.882_HNEG",
         "FT_281.0443_14.974_RPPOS",
         "FT_281.04432_14.693_RPPOS",
         "FT_398.43535_20.583_RPPOS",
         "FT_403.08257_1.338_HNEG",
         "FT_405.09664_18.78_RPPOS",
         "FT_485.06314_15.603_RPPOS",
         "FT_521.12171_1.651_HNEG",
         "FT_521.31204_1.166_HNEG",
         "FT_525.13897_11.044_RPPOS",
         "FT_525.13903_11.885_RPPOS",
         "FT_541.11267_17.415_RPPOS",
         "FT_551.11831_11.992_RPPOS",
         "FT_569.31888_1.531_HNEG",
         "FT_727.18074_13.584_RPPOS",
         "FT_761.18629_14.007_RPPOS",
         "FT_905.25928_2.194_HNEG",
         "FT_905.25933_1.679_HNEG",
         "FT_273.00414_1.083_HNEG",
         "FT_437.08805_1.673_HNEG",
         "FT_485.12433_1.672_HNEG",
         "FT_273.0327_13.767_RPPOS",
         "FT_267.16024_1.664_HNEG",
         "FT_375.08766_1.597_HNEG",
         "FT_215.14292_14.773_RPPOS",
         "FT_573.10383_1.583_HNEG",
         "FT_349.21605_18.806_RPPOS",
         "FT_363.19527_17.556_RPPOS",
         "FT_365.21099_16.701_RPPOS",
         "FT_405.13305_15.694_RPPOS")

NAS_plot <- lapply(NAS,plot_avni_change)



