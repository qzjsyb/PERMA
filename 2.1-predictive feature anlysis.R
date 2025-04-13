library(tidyverse)




#loading ----
setwd("~/PERMA/R")
compounds<-read.csv("compound_link.csv", header = T, row.names = 1) %>%
  rownames_to_column("ID")

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

compounds <- subset(compounds,!compounds$ID%in%c(HNEG_remove, RPPOS_remove))
important_features<-read.csv("important_features_module_1224.csv", header = T) %>%
  mutate(across(where(is.character), ~na_if(., "")))


### plot stacking plot----
count_HNEG <- subset(compounds, grepl("HNEG", ID)) #1046
count_HNEG <- count_HNEG %>%
  mutate(across(everything(), ~na_if(., "")))
nrow(count_HNEG)
count_RPPOS <- subset(compounds, grepl("RPPOS", ID)) #3936
count_RPPOS <- count_RPPOS %>%
  mutate(across(everything(), ~na_if(., "")))
nrow(count_RPPOS)

HNEG_counts <- as.data.frame(table(count_HNEG$Class)) #832
RPPOS_counts <- as.data.frame(table(count_RPPOS$Class)) #2773

combined_counts <- merge(HNEG_counts, RPPOS_counts, by = "Var1", all = TRUE)

# Replace NAs with 0 to ensure correct summation
combined_counts[is.na(combined_counts)] <- 0

# Sum the 'Freq' columns
combined_counts$Total_Freq <- rowSums(combined_counts[, c("Freq.x", "Freq.y")])

# Optional: Clean up column names for better readability
colnames(combined_counts) <- c("Class", "Freq_HNEG", "Freq_RPPOS", "Freq")

# Convert counts to relative counts
combined_counts$Relative_Count <- combined_counts$Freq/(nrow(count_HNEG) + nrow(count_RPPOS))

combined_counts <- combined_counts%>%
  mutate(Class= reorder(Class, Relative_Count))

# Load ggplot2 for plotting
library(ggplot2)

# Create the bar plot
combined.all<-ggplot(combined_counts, aes(x = Relative_Count, y = Class)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Relative Count of Each Class", x = "Relative Count", y = "Class") +
  theme_minimal()

combined_counts$Potential_Count <- combined_counts$Relative_Count*249

important_counts <- as.data.frame(table(important_features$Class))
colnames(important_counts) <- c("Class", "Freq")

important_counts_all <- left_join(combined_counts,important_counts, by = "Class")


long_data <- important_counts_all %>%
  pivot_longer(cols = c(Potential_Count, Freq.y), 
               names_to = "Relative_Count_Type", 
               values_to = "Count")

plot.important<-ggplot(long_data[which(long_data$Freq.x > 10),], aes(x = Count, y = Class, fill = Relative_Count_Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "all important: Relative Count of Each Class", x = "Relative Count", y = "Class") +
  theme_minimal()


ggsave("plot.important.bar.pdf", plot.important, height = 10, width =10)


### Fisher exact test ----
perform_fisher_test <- function(data, total_dataset, total_subset, rare_num) {
    # Filter out rows where Freq.y is NA
  filtered_data <- data[data$Freq.x > rare_num & !is.na(data$Freq.y), ]
  
  # Initialize vectors to store results
  p_values <- numeric(nrow(filtered_data))
  odds_ratios <- numeric(nrow(filtered_data))
  
  # Loop through each row (each chemical class) to perform Fisher's Exact Test
  for (i in 1:nrow(filtered_data)) {
    # Construct the 2x2 table
    class_in_subset <- filtered_data$Freq.y[i]
    non_class_in_subset <- total_subset - class_in_subset
    class_in_dataset <- filtered_data$Freq.x[i] - class_in_subset
    non_class_in_dataset <- total_dataset - total_subset - class_in_dataset
    
    contingency_table <- matrix(c(class_in_subset, non_class_in_subset,
                                  class_in_dataset, non_class_in_dataset),
                                nrow = 2, byrow = TRUE)
    
    # Perform Fisher's Exact Test
    fisher_test_result <- fisher.test(contingency_table)
    
    # Store the p-value and odds ratio
    p_values[i] <- fisher_test_result$p.value
    odds_ratios[i] <- fisher_test_result$estimate
  }
  
  # Add the results back to the dataframe
  adjusted_p_values <- p.adjust(p_values, method = "BH")
  filtered_data$p_value <- p_values
  filtered_data$p_value_adj <- adjusted_p_values
  filtered_data$odds_ratio <- odds_ratios
  
  # Return the resulting dataframe with p-values and odds ratios
  return(filtered_data)
}

# Example usage of the function

# total_dataset = 4157
# total_subset = 100
important_all_result <- perform_fisher_test(important_counts_all, 4982, 198, 4)

write.csv(important_all_result, "important_all_result.csv")


important_features$Group <- factor(important_features$Group, 
                                   levels = c("Palsa_pos", "Bog_pos", "Fen_pos"))



count_vector <- c(48, 79, 71)

# try to compare with all features
result_list <- list()

group_levels <- c("Palsa_pos", "Bog_pos", "Fen_pos")
# Loop through all combinations of X.2 and Group levels
for (group in group_levels) {
  # Filter the data based on the current combination
  filtered_data <- important_features %>%
    filter(Group == group)
  
  # Count occurrences of each Class
  count_data <- as.data.frame(table(filtered_data$Class))
  
  # Rename columns
  colnames(count_data) <- c("Class", "Freq")
  
  join_table_name <- "combined_counts"
  join_table <- get(join_table_name, envir = .GlobalEnv)
  
  count_data <- left_join(join_table,count_data, by = "Class")
  
  # Store the result in the list
  result_list[[paste0(group)]] <- count_data
}


for (i in 1:3) {
  # Perform the Fisher test
  fisher_result <- perform_fisher_test(
    result_list[[i]],
    4982,
    count_vector[i],
    4
  )
  
  # Create a dynamic name for the result using the name of the list element
  result_name <- paste0(names(result_list)[i], "_result")
  
  # Assign the Fisher test result to a new object in the global environment
  assign(result_name, fisher_result, envir = .GlobalEnv)
}


# Filter the dataframes based on the condition p_value_adj < 0.05
filtered_bog <- Bog_pos_result %>% filter(p_value_adj < 0.05) %>% mutate(Source = "Bog")
filtered_fen <- Fen_pos_result %>% filter(p_value_adj < 0.05) %>% mutate(Source = "Fen")
filtered_palsa <- Palsa_pos_result %>% filter(p_value_adj < 0.05) %>% mutate(Source = "Palsa")

# Merge the dataframes
combined_data <- bind_rows(filtered_bog, filtered_fen, filtered_palsa)

# Create the point plot
class <- ggplot(combined_data, aes(x = log(odds_ratio), y = Class, color = Source)) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.7) +
  scale_color_manual(values = c("Palsa" = "#B8912F", "Bog" = "#5E813F", "Fen" = "#4273B0")) +
  labs(
    x = "Log(Odds Ratio)",
    y = "Class",
    size = "Feq.y",
    color = "Source"
  ) +
  theme_classic() +
  theme(legend.position = "none")




### test mz and RT -----
important_features$mz <- as.numeric(gsub("FT_([0-9.]+)_.+", "\\1", as.character(important_features$ID)))


important_features <- important_features %>%
  mutate(neutral_weight = case_when(
    grepl("RPPOS", ID) ~ mz - 1.00728,
    grepl("HNEG", ID) ~ mz + 1.00728
  ))

important_features$mz <- as.numeric(gsub("FT_([0-9.]+)_.+", "\\1", as.character(important_features$ID)))

compounds$mz<- as.numeric(gsub("FT_([0-9.]+)_.+", "\\1", as.character(compounds$ID)))
compounds <- compounds%>%
  mutate(neutral_weight = case_when(
    grepl("RPPOS", ID) ~ mz - 1.00728,
    grepl("HNEG", ID) ~ mz + 1.00728
  ))


important_features <- important_features %>%
  mutate(neutral_weight = case_when(
    grepl("RPPOS", ID) ~ mz - 1.00728,
    grepl("HNEG", ID) ~ mz + 1.00728
  ))

compounds$RT <- as.numeric(gsub("FT_[0-9.]+_([0-9.]+)_.+", "\\1", as.character(compounds$ID)))
important_features$RT<- as.numeric(gsub("FT_[0-9.]+_([0-9.]+)_.+", "\\1", as.character(important_features$ID)))


# density plot of neutral_weight and RT-----
important_features %>%
  filter(!((Bog < 0 & Palsa == 0 & Fen == 0) |  # Exclude Bog < 0 and Palsa = 0 & Fen = 0
             (Fen < 0 & Palsa == 0 & Bog == 0))) %>%  # Exclude Fen < 0 and Palsa = 0 & Bog = 0
  ggplot(., aes(x = neutral_weight)) +
  
  geom_density(aes(fill = case_when(
    Palsa > 0 ~ "Palsa_pos",
    Palsa < 0 ~ "Palsa_neg",
    Bog > 0 ~ "Bog_pos",
    Fen > 0 ~ "Fen_pos"
  )), alpha = 0.5, color = NA) +  # Remove outline for group densities
  
  # Add baseline density for compounds$neutral_weight
  geom_density(data = compounds, aes(x = neutral_weight), 
               color = "black", fill = NA) +  # Black outline, no fill for baseline
  # Customize fill and color scales
  scale_fill_manual(values = c("Palsa_pos" = "red", "Palsa_neg" = "darkred",
                               "Bog_pos" = "green", "Fen_pos" = "blue"))  +
  
  # Add titles and labels
  labs(title = "Density Plot of Neutral Weight across important feature Groups",
       x = "Neutral Weight",
       y = "Density") +
  
  # Apply minimal theme and adjust legend
  theme_minimal() +
  theme(legend.position = "bottom", legend.direction = "horizontal")

#RT plot.
important_RP_features_RT <- important_features %>%
  filter(grepl("RPPOS", ID))

important_RP_features_RT %>%
  filter(!((Bog < 0 & Palsa == 0 & Fen == 0) |  # Exclude Bog < 0 and Palsa = 0 & Fen = 0
             (Fen < 0 & Palsa == 0 & Bog == 0))) %>%  # Exclude Fen < 0 and Palsa = 0 & Bog = 0
  ggplot(., aes(x = RT)) +
  
  geom_density(aes(fill = case_when(
    Palsa > 0 ~ "Palsa_pos",
    Palsa < 0 ~ "Palsa_neg",
    Bog > 0 ~ "Bog_pos",
    Fen > 0 ~ "Fen_pos"
  )), alpha = 0.5, color = NA) +  # Remove outline for group densities +
  
  # Add density plot for important_features with no outline
  geom_density(alpha = 0.5, color = NA) +  # Remove outline for group densities
  
  # Add baseline density for compounds$RT
  geom_density(data = compounds %>% filter(grepl("RPPOS", ID)), aes(x = RT), 
               color = "black", fill = NA) +  # Black outline, no fill for baseline
  
  # Customize fill and color scales
  scale_fill_manual(values = c("Palsa_pos" = "red", "Palsa_neg" = "darkred",
                               "Bog_pos" = "green", "Fen_pos" = "blue")) +
  
  # Add titles and labels
  labs(title = "Density Plot Retention Time for Reverse Phase across Groups of important features",
       x = "Retention Time (RT)",
       y = "Density") +
  
  # Apply minimal theme and adjust legend
  theme_minimal() +
  theme(legend.position = "bottom", legend.direction = "horizontal")


#van Krevelen -----
extract_atom_count <- function(formula, element) {
  # Use regex to find the number following the element (e.g., C, H, O)
  match <- regmatches(formula, regexec(paste0(element, "([0-9]*)"), formula))
  if (length(match[[1]]) > 1 && match[[1]][2] != "") {
    as.numeric(match[[1]][2])
  } else {
    1  # If no number is found, assume it's 1 (e.g., C, H, O without explicit number)
  }
}

important_features <- important_features %>%
  mutate(
    # Extract number of atoms for C, H, and O
    C_count = sapply(Formula, extract_atom_count, element = "C"),
    H_count = sapply(Formula, extract_atom_count, element = "H"),
    O_count = sapply(Formula, extract_atom_count, element = "O"),
    
    # Calculate O:C and H:C ratios
    OC_ratio = O_count / C_count,
    HC_ratio = H_count / C_count
  )

classification <- tribble(
  ~Class, ~OC_low, ~OC_high, ~HC_low, ~HC_high,
  'Lipid', 0, 0.3, 1.5, 2.5,
  'Unsat. HC', 0, 0.125, 0.8, 1.5,
  'Cond. HC', 0, 0.95, 0.2, 0.8,
  'Protein', 0.3, 0.55, 1.5, 2.3,
  'Amino sugar', 0.55, 0.7, 1.5, 2.2,
  'Carbohydrate', 0.7, 1.5, 1.5, 2.5,
  'Lignin', 0.125, 0.65, 0.8, 1.5,
  'Tannin', 0.65, 1.1, 0.8, 1.5, 
) %>% 
  mutate(label_x = (OC_low + OC_high) / 2,
         label_y = HC_high - 0.1)

## Compound class rectangles (for plotting of Van Krevelen diagrams) 

class_rect <-  geom_rect(data = classification,
                         aes(xmin = OC_low,
                             xmax = OC_high,
                             ymin = HC_low,
                             ymax = HC_high),
                         color = 'black',
                         fill = NA,
                         linewidth = 1,
                         inherit.aes = FALSE, 
                         linetype = 'dashed')
rect_label <- geom_label(data = classification,
                         aes(x = label_x,
                             y = label_y,
                             label = Class),
                         inherit.aes = FALSE,
                         size = 3)


vk <- ggplot(important_features%>%
         filter(!grepl("Bad", Tags)), aes(x = OC_ratio, y = HC_ratio, alpha = 0.5, color = case_when(
  Palsa > 0 ~ "Palsa_pos",
  Bog > 0 ~ "Bog_pos",
  Fen > 0 ~ "Fen_pos",
))) +
  geom_point(size = 3, alpha = 0.7) +  # Add points with transparency
  labs(title = "Van Krevelen Diagram",
       x = "O:C Ratio",
       y = "H:C Ratio",
       color = "Group") +
  class_rect  +
  rect_label +
  theme_minimal() +  # Use a minimal theme for clarity
  scale_color_manual(values = c("Palsa_pos" = "#B8912F", "Palsa_neg" = "#EAD8A6",
                                "Bog_pos" = "#5E813F", "Fen_pos" = "#4273B0",
                                "Bog_neg" = "#B8D0AB", "Fen_neg" = "#AFC9E5")) +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  xlim(0, 1.5) +  # Limit x-axis from 0 to 2
  ylim(0, 3)    # Limit y-axis from 0 to 3




### nosc ----
nosc <- important_features%>%
  filter(!grepl("Bad", Tags)) %>%
  #Use OrgMassSpecR package to parse the chemical formula - Returns a named list
  mutate(breakdown = purrr::map(Formula, function(x) data.frame(OrgMassSpecR::ListFormula(x)))) %>%
  #Unnest that list into multiple columns, each named after the element
  unnest('breakdown') %>%
  #Use those columns to calculate the NOSC value
  mutate(NOSC = 4 - ((4*C + H - 2*O - 3*N - 2*S + 5*P)/C)) 

annotation_source <- read.csv("PERMA_annotationsources.csv")
non_full_match_ids <- annotation_source %>%
  filter(Annot..Source..Predicted.Compositions != "Full match") %>%
  pull(ID)  # Extract the IDs for non-full match formulas

nosc_all <- compounds%>%
  filter(!ID %in% c(HNEG_remove, RPPOS_remove))%>%
  filter(!ID %in% non_full_match_ids) %>%
  #Use OrgMassSpecR package to parse the chemical formula - Returns a named list
  mutate(breakdown = purrr::map(Formula, function(x) data.frame(OrgMassSpecR::ListFormula(x)))) %>%
  #Unnest that list into multiple columns, each named after the element
  unnest('breakdown') %>%
  #Use those columns to calculate the NOSC value 
  mutate(NOSC = 4 - ((4*C + H - 2*O - 3*N - 2*S + 5*P)/C)) %>%
  filter(NOSC > -1.6 & NOSC < 1.2 )




nosc_palsa_pos <- nosc %>%
  filter(NOSC > -1.6 & NOSC < 1.2) %>%
  filter(Palsa > 0)
# 
# nosc_palsa_neg <- nosc %>%
#   filter(NOSC > -1.6 & NOSC < 1.2) %>%
#   filter(Palsa < 0)

nosc_bog_pos <- nosc %>%
  filter(NOSC > -1.6 & NOSC < 1.2) %>%
  filter(Bog > 0)

nosc_fen_pos <- nosc %>%
  filter(NOSC > -1.6 & NOSC < 1.2) %>%
  filter(Fen > 0)

# Create empirical distribution function from full dataset
Fn <- ecdf(nosc_all$NOSC)


ks_palsa_pos<-ks.test(jitter(nosc_palsa_pos$NOSC), Fn)
# ks_palsa_neg<-ks.test(nosc_palsa_neg$NOSC, Fn)
ks_bog_pos<-ks.test(jitter(nosc_bog_pos$NOSC), Fn)
ks_fen_pos<-ks.test(jitter(nosc_fen_pos$NOSC), Fn)



library(ggpubr)



#NOSC FOR ALL
annotations <- c(
  paste0("D = ", round(ks_palsa_pos$statistic, 3), ", p = ", round(ks_palsa_pos$p.value, 3)),
  # paste0("D = ", round(ks_palsa_neg$statistic, 3), ", p = ", round(ks_palsa_neg$p.value, 3)),
  paste0("D = ", round(ks_bog_pos$statistic, 3), ", p = ", round(ks_bog_pos$p.value, 3)),
  paste0("D = ", round(ks_fen_pos$statistic, 3), ", p = ", round(ks_fen_pos$p.value, 3))
)

plot_nosc_all <- nosc %>%
  filter(NOSC > -1.6 & NOSC < 1.2) %>%
  filter(!((Bog < 0 & Palsa == 0 & Fen == 0) |  # Exclude Bog < 0 and Palsa = 0 & Fen = 0
             (Fen < 0 & Palsa == 0 & Bog == 0) |
             (Palsa < 0 & Bog == 0 & Fen == 0)
           )) %>%  # Exclude Fen < 0 and Palsa = 0 & Bog = 0
  ggplot(aes(x = factor(case_when(
    Palsa > 0 ~ "Palsa_pos",
    # Palsa < 0 ~ "Palsa_neg",
    Bog > 0 ~ "Bog_pos",
    Fen > 0 ~ "Fen_pos"
  ), levels = c("Palsa_pos","Bog_pos", "Fen_pos")),  # Reorder levels
  y = NOSC, fill = case_when(
    Palsa > 0 ~ "Palsa_pos",
    # Palsa < 0 ~ "Palsa_neg",
    Bog > 0 ~ "Bog_pos",
    Fen > 0 ~ "Fen_pos"
  ), color = case_when(
    Palsa > 0 ~ "Palsa_pos",
    # Palsa < 0 ~ "Palsa_neg",
    Bog > 0 ~ "Bog_pos",
    Fen > 0 ~ "Fen_pos"
  ))) +
  
  geom_violin(alpha = 0.5) +  # Violin plot for group-wise distribution
  
  geom_jitter(position = position_jitter(width = 0.2), size = 1.5, alpha = 0.7) +
  
  
  # Violin plot for nosc_all (overall distribution)
  geom_violin(data = nosc_all, aes(x = "Overall", y = NOSC), fill = "grey", color = "grey", alpha = 0.5) +
  
  # Add horizontal line at y = -0.3
  geom_hline(yintercept = 0, linetype = "solid", color = "darkgreen") +
  annotate("text", x = 3.5, y = 0.05, label = "Average for Plant Metabolites", color = "darkgreen", hjust = 0) +
  
  geom_hline(yintercept = -0.3, linetype = "solid", color = "black") +
  annotate("text", x = 3.5, y = -0.25, label = "Average for Soil Organic Matter", color = "black", hjust = 0) +
  
  
  # Customize fill and color scales
  scale_fill_manual(values = c("Palsa_pos" = "#B8912F",
                               "Bog_pos" = "#5E813F", "Fen_pos" = "#4273B0")) +
  
  # Customize fill and color scales
  scale_color_manual(values = c("Palsa_pos" = "#B8912F",
                                "Bog_pos" = "#5E813F", "Fen_pos" = "#4273B0")) +
  
  # Add titles and labels
  labs(title = "Violin Plot of NOSC across Groups of Predictive Features",
       x = "Group",
       y = "NOSC") +
  
  # Apply minimal theme and adjust legend
  theme_classic() +
  theme(legend.position = "none")  +
  geom_signif(comparisons = list(c("Palsa_pos", "Overall"), 
                                 c("Bog_pos", "Overall"),
                                 c("Fen_pos", "Overall")),
              annotations = annotations,  # Annotate with D and p-values
              y_position = c(1.6, 1.4, 1.2, 1),
              color = "black",# Adjust y positions for the lines
              tip_length = 0.1)  # Add some space for the lines

#NOSC FOR SULFER -----

nosc_s <- nosc %>%
  filter(NOSC > -1.6 & NOSC < 1.2) %>%
  filter(!((Bog < 0 & Palsa == 0 & Fen == 0) | 
             (Fen < 0 & Palsa == 0 & Bog == 0) |
             (Palsa < 0 & Bog == 0 & Fen == 0))) %>%
  filter(S >= 1) %>%
  mutate(group = case_when(
    Palsa > 0 ~ "Palsa_pos",
    Bog > 0 ~ "Bog_pos",
    Fen > 0 ~ "Fen_pos"
  ))

s_label <- nosc_s%>%
  group_by(group) %>%
  summarise(n = n()) %>%
  mutate(label = paste0(group, " (N = ", n, ")")) %>%
  dplyr::select(group, label)

plot_nosc_s <- nosc_s %>%
  ggplot(aes(x = factor(case_when(
    Palsa > 0 ~ "Palsa_pos",
    # Palsa < 0 ~ "Palsa_neg",
    Bog > 0 ~ "Bog_pos",
    Fen > 0 ~ "Fen_pos"
  ), levels = c("Palsa_pos","Bog_pos", "Fen_pos")),  # Reorder levels
  y = NOSC, fill = case_when(
    Palsa > 0 ~ "Palsa_pos",
    # Palsa < 0 ~ "Palsa_neg",
    Bog > 0 ~ "Bog_pos",
    Fen > 0 ~ "Fen_pos"
  ), color = case_when(
    Palsa > 0 ~ "Palsa_pos",
    # Palsa < 0 ~ "Palsa_neg",
    Bog > 0 ~ "Bog_pos",
    Fen > 0 ~ "Fen_pos"
  ))) +
  
  geom_boxplot(alpha = 0.5, outlier.shape = NA, color = "black") +  # Violin plot for group-wise distribution
  geom_jitter(position = position_jitter(width = 0.2), size = 1.5, alpha = 0.7) +
  # Customize fill and color scales
  scale_fill_manual(values = c("Palsa_pos" = "#B8912F","Bog_pos" = "#5E813F", "Fen_pos" = "#4273B0")) +
  # Customize fill and color scales
  scale_color_manual(values = c("Palsa_pos" = "#B8912F","Bog_pos" = "#5E813F", "Fen_pos" = "#4273B0")) +
  
  # Add titles and labels
  labs(title = "NOSC of S-containing Predictive Features",
       x = "Group",
       y = "NOSC") +
  
  # Apply minimal theme and adjust legend
  theme_classic() +
  theme(legend.position = "none")

nosc_s$group <- factor(nosc_s$group, levels = c("Palsa_pos", "Bog_pos", "Fen_pos"))

library(FSA)
library(rstatix)

# Perform Kruskal-Wallis test
kruskal_test <- kruskal_test(NOSC ~ group, data = nosc_s)

dunn_test <- dunnTest(NOSC ~ group, data = nosc_s, method = "bonferroni")
pvals <- dunn_test$res$P.adj
annotations <- ifelse(pvals < 0.05, paste0("P = ", round(pvals, 3)), "ns")

# Add significance annotations using geom_signif
plot_nosc_s <- plot_nosc_s +
  geom_signif(comparisons = list(c("Bog_pos", "Fen_pos"),
                                 c("Bog_pos","Palsa_pos"), 
                                 c("Fen_pos","Palsa_pos")),
              annotations = annotations,
              y_position = c(1.4, 1.55, 1.7),  # Adjust positions for clarity
              tip_length = 0.05, 
              color = "black")

#nosc-nitrogen
nosc_n <- nosc %>%
  filter(NOSC > -1.6 & NOSC < 1.2) %>%
  filter(!((Bog < 0 & Palsa == 0 & Fen == 0) | 
             (Fen < 0 & Palsa == 0 & Bog == 0) |
             (Palsa < 0 & Bog == 0 & Fen == 0))) %>%
  filter(N >= 1) %>%
  mutate(group = case_when(
    Palsa > 0 ~ "Palsa_pos",
    Bog > 0 ~ "Bog_pos",
    Fen > 0 ~ "Fen_pos"
  ))

n_label <- nosc_n%>%
  group_by(group) %>%
  summarise(n = n()) %>%
  mutate(label = paste0(group, " (N = ", n, ")")) %>%
  select(group, label)

plot_nosc_n <- nosc_n %>%
  ggplot(aes(x = factor(case_when(
    Palsa > 0 ~ "Palsa_pos",
    Bog > 0 ~ "Bog_pos",
    Fen > 0 ~ "Fen_pos"
  ), levels = c("Palsa_pos","Bog_pos", "Fen_pos")),  # Reorder levels
  y = NOSC, fill = case_when(
    Palsa > 0 ~ "Palsa_pos",
    Bog > 0 ~ "Bog_pos",
    Fen > 0 ~ "Fen_pos"
  ), color = case_when(
    Palsa > 0 ~ "Palsa_pos",
    Bog > 0 ~ "Bog_pos",
    Fen > 0 ~ "Fen_pos"
  ))) +
  
  geom_boxplot(alpha = 0.5, outlier.shape = NA, color = "black") +  # Violin plot for group-wise distribution
  geom_jitter(position = position_jitter(width = 0.2), size = 1.5, alpha = 0.7) +
  # Customize fill and color scales
  scale_fill_manual(values = c("Palsa_pos" = "#B8912F","Bog_pos" = "#5E813F", "Fen_pos" = "#4273B0")) +
  # Customize fill and color scales
  scale_color_manual(values = c("Palsa_pos" = "#B8912F","Bog_pos" = "#5E813F", "Fen_pos" = "#4273B0")) +
  
  # Add titles and labels
  labs(title = "NOSC of N-containing Predictive Features",
       x = "Group",
       y = "NOSC") +
  
  # Apply minimal theme and adjust legend
  theme_classic() +
  theme(legend.position = "none") 

nosc_n$group <- factor(nosc_n$group, levels = c("Palsa_pos", "Bog_pos", "Fen_pos"))

# Perform Kruskal-Wallis test
kruskal_test <- kruskal_test(NOSC ~ group, data = nosc_n)
dunn_test <- dunnTest(NOSC ~ group, data = nosc_n, method = "bonferroni")
pvals <- dunn_test$res$P.adj
annotations <- ifelse(pvals < 0.05, paste0("P = ", round(pvals, 3)), "ns")

# Add significance annotations using geom_signif
plot_nosc_n <- plot_nosc_n +
  geom_signif(comparisons = list(c("Bog_pos", "Fen_pos"),
                                 c("Bog_pos","Palsa_pos"), 
                                 c("Fen_pos","Palsa_pos")),
              annotations = annotations,
              y_position = c(1.55, 1.4, 1.7),  # Adjust positions for clarity
              tip_length = 0.05, 
              color = "black")




nosc_n_s <- ggarrange(plot_nosc_n, plot_nosc_s, ncol = 2,labels = c("3","4"))

class_nosc <- ggarrange(
  class, 
  plot_nosc_all, 
  nrow = 2,                # Arrange plots in 2 rows
  heights = c(1, 2),       # Custom heights for the rows
  labels = c("1", "2")     # Labels for the subplots
)

 plot_all <- ggarrange(
    class_nosc,               # First plot on the first row
    nosc_n_s,
  nrow = 2,heights = c(2, 1))     # Share a common legend across plots)



nosc_vk <- ggarrange(vk_class,plot_nosc_together, nrow = 2)
ggsave("nosc_vk.pdf", plot_all, width = 5, height = 8)





#nosc-phosphorus
nosc_p <- nosc %>%
  filter(NOSC > -1.6 & NOSC < 1.2) %>%
  filter(!((Bog < 0 & Palsa == 0 & Fen == 0) | 
             (Fen < 0 & Palsa == 0 & Bog == 0) |
             (Palsa < 0 & Bog == 0 & Fen == 0))) %>%
  filter(P >= 1) %>%
  mutate(group = case_when(
    Palsa > 0 ~ "Palsa_pos",
    Bog > 0 ~ "Bog_pos",
    Fen > 0 ~ "Fen_pos"
  ))

p_label <- nosc_p%>%
  group_by(group) %>%
  summarise(n = n()) %>%
  mutate(label = paste0(group, " (N = ", n, ")")) %>%
  select(group, label)

plot_nosc_p <- nosc_p %>%
  ggplot(aes(x = factor(case_when(
    Palsa > 0 ~ "Palsa_pos",
    Bog > 0 ~ "Bog_pos",
    Fen > 0 ~ "Fen_pos"
  ), levels = c("Palsa_pos","Bog_pos", "Fen_pos")),  # Reorder levels
  y = NOSC, fill = case_when(
    Palsa > 0 ~ "Palsa_pos",
    Bog > 0 ~ "Bog_pos",
    Fen > 0 ~ "Fen_pos"
  ), color = case_when(
    Palsa > 0 ~ "Palsa_pos",
    Bog > 0 ~ "Bog_pos",
    Fen > 0 ~ "Fen_pos"
  ))) +
  
  geom_boxplot(alpha = 0.5, outlier.shape = NA, color = "black") +  # Violin plot for group-wise distribution
  geom_jitter(position = position_jitter(width = 0.2), size = 1.5, alpha = 0.7) +
  # Customize fill and color scales
  scale_fill_manual(values = c("Palsa_pos" = "#B8912F","Bog_pos" = "#5E813F", "Fen_pos" = "#4273B0")) +
  # Customize fill and color scales
  scale_color_manual(values = c("Palsa_pos" = "#B8912F","Bog_pos" = "#5E813F", "Fen_pos" = "#4273B0")) +
  
  # Add titles and labels
  labs(title = "NOSC of P-containing Predictive Features",
       x = "Group",
       y = "NOSC") +
  
  # Apply minimal theme and adjust legend
  theme_linedraw() +
  theme(legend.position = "bottom", legend.direction = "horizontal") 

nosc_p$group <- factor(nosc_p$group, levels = c("Palsa_pos", "Bog_pos", "Fen_pos"))

# Perform Kruskal-Wallis test
kruskal_test <- kruskal_test(NOSC ~ group, data = nosc_p)
dunn_test <- dunnTest(NOSC ~ group, data = nosc_p, method = "bonferroni")
pvals <- dunn_test$res$P.adj
annotations <- ifelse(pvals < 0.05, paste0("P = ", round(pvals, 3)), "ns")

# Add significance annotations using geom_signif
plot_nosc_p <- plot_nosc_p +
  geom_signif(comparisons = list(c("Bog_pos", "Fen_pos"),
                                 c("Bog_pos","Palsa_pos"), 
                                 c("Fen_pos","Palsa_pos")),
              annotations = annotations,
              y_position = c(1.55, 1.4, 1.7),  # Adjust positions for clarity
              tip_length = 0.05, 
              color = "black")



#neutral weight
important_features %>%
  filter(!((Bog < 0 & Palsa == 0 & Fen == 0) |  # Exclude Bog < 0 and Palsa = 0 & Fen = 0
             (Fen < 0 & Palsa == 0 & Bog == 0))) %>%  # Exclude Fen < 0 and Palsa = 0 & Bog = 0
  ggplot(aes(x = factor(case_when(
    Palsa > 0 ~ "Palsa_pos",
    Palsa < 0 ~ "Palsa_neg",
    Bog > 0 ~ "Bog_pos",
    Fen > 0 ~ "Fen_pos"
  ), levels = c("Palsa_pos", "Palsa_neg", "Bog_pos", "Fen_pos")),  # Reorder levels
  y = neutral_weight, fill = case_when(
    Palsa > 0 ~ "Palsa_pos",
    Palsa < 0 ~ "Palsa_neg",
    Bog > 0 ~ "Bog_pos",
    Fen > 0 ~ "Fen_pos"
  ), color = case_when(
    Palsa > 0 ~ "Palsa_pos",
    Palsa < 0 ~ "Palsa_neg",
    Bog > 0 ~ "Bog_pos",
    Fen > 0 ~ "Fen_pos"
  ))) +
  
  geom_violin(alpha = 0.5, color = "black") +  # Violin plot for group-wise distribution
  
  geom_jitter(position = position_jitter(width = 0.2), size = 1.5, alpha = 0.7) +
  
  # Add baseline density for compounds$RT
  geom_violin(data = important_features, aes(x = "Overall", y =  neutral_weight), fill = "grey", color = "black", alpha = 0.5) +
  
  # Customize fill and color scales
  scale_fill_manual(values = c("Palsa_pos" = "red", "Palsa_neg" = "darkred",
                               "Bog_pos" = "green", "Fen_pos" = "blue")) +
  scale_color_manual(values = c("Palsa_pos" = "red", "Palsa_neg" = "darkred",
                               "Bog_pos" = "green", "Fen_pos" = "blue")) +
  
  # Add titles and labels
  labs(title = "Density Plot Retention Time for Reverse Phase across Groups of important features",
       x = "Group",
       y = "Neutral Weight") +
  
  # Apply minimal theme and adjust legend
  theme_minimal() +
  theme(legend.position = "bottom", legend.direction = "horizontal")

#oxidative states:




