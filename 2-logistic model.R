#logistical regresssion

library(dplyr)
library(lme4)
library(ggplot2)
library(caret)
library(tidyverse)
library(VennDiagram)
library(reprtree)
library(MLmetrics)
library(glmnet)
library(ggtern)
library(ggpubr)
library(reshape2)
###functions ----
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

view_result <- function(result){
  result %>%
    ggplot(aes(x = depth, y = percent_correct, color = habitat)) +
    geom_point(aes(size = total)) +  # Add points
    geom_line(aes(group = habitat), linewidth = 1) +  # Add lines
    scale_color_manual(values = habitat_colors) +  # Set manual colors
    labs(title = "Percentage of Correct Predictions by Habitat and Depth",
         x = "Depth",
         y = "Percent Match",
         color = "Habitat") +  # Label the legend
    theme_minimal()  # Use a minimal theme
}


#loading data ----
setwd("~/PERMA/R")
metadata <- read.csv("metadata.all.csv", row.names = 1)
metadata$habitat <- factor(metadata$habitat, levels = c("Palsa", "Bog", "Fen"))
depth_lookup <- data.frame(depth_category = c("1_5", "10_14", "20_24", "30_34", "40_44", "50_54", "60_64", "70_74", "80_84", "90_94"), 
                           depth_numeric = c(3, 12, 22, 32, 42, 52, 62, 72, 82, 92)
)
metadata$depth_n <- depth_lookup$depth_numeric[match(metadata$depth, depth_lookup$depth_category)]
metadata <- metadata[which(metadata$depth_n < 40),]

metadata$core <- paste0(metadata$Year,metadata$Habitat_plant,metadata$Sample_code)
habitat_colors <- c("Palsa" = "#B8912F", "Bog" = "#5E813F", "Fen" = "#4273B0", "Collapsed Palsa" = "orange" , "Poor Fen" = "blue")

metadata$pH <- as.numeric(metadata$pH)
metadata$Tsoil <- as.numeric(metadata$Tsoil)

compounds <- read.csv("compound_link.csv", header = T, row.names = 1) %>%
  rownames_to_column("ID")



# HNEG------
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


metadata.HNEG <- metadata[rownames(metadata) != "2014_S_2_20_24", ]
#use autochamber samples as training set
metadata.HNEG.ac <- subset(metadata.HNEG, grepl("Autochamber",metadata.HNEG$site))
#use avni and inoc as test set
metadata.HNEG.test <- subset(metadata.HNEG, !rownames(metadata.HNEG)%in% rownames(metadata.HNEG.ac))

HNEG.ac <- HNEG[rownames(metadata.HNEG.ac),]
metadata.HNEG.ac <- metadata.HNEG[rownames(HNEG.ac),]
HNEG.ac.pareto <- na.omit(HNEG.ac) %>%
  log_transform() %>%
  pareto_scale()


# RPPOS------
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

metadata.RPPOS <- metadata
#use autochamber samples as training set
metadata.RPPOS.ac <- subset(metadata.RPPOS, grepl("Autochamber",metadata.RPPOS$site))
#use avni and inoc as test set
metadata.RPPOS.test <- subset(metadata.RPPOS, !rownames(metadata.RPPOS)%in% rownames(metadata.RPPOS.ac))

RPPOS.ac <- RPPOS[rownames(metadata.RPPOS.ac),]
metadata.RPPOS.ac <- metadata.RPPOS[rownames(RPPOS.ac),]
RPPOS.ac.pareto <- na.omit(RPPOS.ac) %>%
  log_transform() %>%
  pareto_scale()




### include IncE, IncS, AvCj_P in training model -----


#HNEG
metadata.HNEG <- metadata[rownames(metadata) != "2014_S_2_20_24", ]
#use autochamber samples as training set
metadata.HNEG.ac <- subset(metadata.HNEG, grepl("Autochamber",metadata.HNEG$site))
#use avni and inoc as test set
metadata.HNEG.test <- subset(metadata.HNEG, !rownames(metadata.HNEG)%in% rownames(metadata.HNEG.ac))


rows_to_move <- metadata.HNEG.test %>%
  filter(Habitat_plant %in% c("AvCj_P", "IncE", "IncS"))

metadata.HNEG.ac <- bind_rows(metadata.HNEG.ac, rows_to_move)
#use avni and inoc as test set
metadata.HNEG.test <- metadata.HNEG.test %>%
  filter(!row.names(metadata.HNEG.test) %in% row.names(rows_to_move))

HNEG.ac <- HNEG[rownames(metadata.HNEG.ac),]
metadata.HNEG.ac <- metadata.HNEG[rownames(HNEG.ac),]
HNEG.ac.pareto <- na.omit(HNEG.ac) %>%
  log_transform() %>%
  pareto_scale()

set.seed(123)
train_control <- trainControl(method = "repeatedcv", number = 10, repeats=3, classProbs = TRUE, summaryFunction = multiClassSummary)

HNEG.all <- data.frame(habitat = metadata.HNEG.ac$habitat, depth = metadata.HNEG.ac$depth, HNEG.ac.pareto)

p = 0.7
HNEG.all$composite_var <- interaction(HNEG.all$habitat, HNEG.all$depth)
index <- createDataPartition(HNEG.all$composite_var, p = p, list = FALSE)
trainData <- HNEG.all[index, !names(HNEG.all) %in% c("depth", "composite_var")]
testData <- HNEG.all[-index, !names(HNEG.all) %in% c("depth", "composite_var")]

model_hilic <- train(habitat ~ ., 
               data = trainData, 
               metric = "Accuracy",
               method = "glmnet", 
               trControl = train_control,
               tuneGrid = expand.grid(alpha = 0.5,  
                                      lambda = 10^seq(-3, 1, length = 20)))
print(model_hilic)

predictions <- predict(model_hilic, newdata = testData)
#lambda        logLoss     AUC        prAUC      Accuracy   Kappa      Mean_F1    Mean_Sensitivity  Mean_Specificity  Mean_Pos_Pred_Value
#0.011288379  0.05306522  0.9995222  0.8709276  0.9914251  0.9870862  0.9911781  0.9908730         0.9956481         0.9925926            0.9960739

best_model <- model_hilic$finalModel
best_lambda <- model_hilic$bestTune$lambda 
best_alpha <- model_hilic$bestTune$alpha

# Confusion matrix for test set performance
confusionMatrix(predictions, testData$habitat)
confusionMatrix(predictions_adjusted, testData$habitat)

result_table <- model_hilic$results


metrics <- result_table[, c("lambda", "logLoss", "Accuracy")]

# Melt the data for ggplot2 (long format)
metrics_long <- melt(metrics, id.vars = "lambda", 
                     variable.name = "Metric", 
                     value.name = "Value")

# Plot metrics vs. lambda
hilic_tuning <- ggplot(metrics_long, aes(x = lambda, y = Value, color = Metric)) +
  geom_line() +
  geom_point() +
  scale_x_log10() +  # Log scale for lambda
  geom_vline(xintercept = best_lambda, color = "red", linetype = "dashed") +  # Add vertical line
  labs(title = paste("HILIC Tuning (Best Lambda =", signif(best_lambda, digits = 3), ")"),
       x = "Lambda (log scale)",
       y = "Metric Value") +
  theme_classic() +
  theme(legend.position = "bottom") +
  scale_color_brewer(palette = "Set1")


## model
HNEG.lr <- glmnet(HNEG.all[!names(HNEG.all) %in% c("habitat","depth", "composite_var")], HNEG.all$habitat, alpha = best_alpha, lambda = best_lambda, family = "multinomial")
coef_matrix_hn <- coef(HNEG.lr, s = best_lambda)

# Convert the coefficients to a data frame for easier interpretation
coef_df_hn <- do.call(cbind, lapply(coef_matrix_hn, as.matrix))
colnames(coef_df_hn) <- c( "Palsa","Bog", "Fen")

# run test
HNEG.test.pareto <- HNEG[rownames(metadata.HNEG.test),] %>%
  log_transform() %>%
  pareto_scale()%>%
  as.matrix(.)

predicted_p.HNEG <- predict(HNEG.lr, newx = HNEG.test.pareto, s = best_lambda, type = "response")
predicted_c.HNEG <- predict(HNEG.lr, newx = HNEG.test.pareto, s = best_lambda, type = "class")

df <- as.data.frame(predicted_p.HNEG)

metadata.HNEG.test <- subset(metadata.HNEG, !rownames(metadata.HNEG)%in% rownames(metadata.HNEG.ac))
metadata.HNEG.test$predict <- predicted_c.HNEG
metadata.HNEG.test <- cbind (metadata.HNEG.test,df)

#diagnose plot
label <- metadata.HNEG.test %>%
  filter(depth == "1_5")
label$name <- paste0(label$label, gsub("core", "", label$Sample_code))

site_color <- c("#5E813F","#4273B0","purple4","red4","red1", "#B8912F","#5E813F","#4273B0","#B8912F")
names(site_color) <-  c("IncS","IncE", "AvCe_F1", "AvCe_C", "AvC_P", "AvCj_P","Bog", "Fen","Palsa")

metadata.HNEG.test <- metadata.HNEG.test %>%
  mutate(Habitat_plant = case_when(
    Habitat_plant == "AvCe_F1" ~ "PT",
    Habitat_plant == "AvCe_C" ~ "ET",
    Habitat_plant == "AvC_P" ~ "RT",
    TRUE ~ Habitat_plant  # Keep other values unchanged
  ))

# Step 2: Pivot longer
metadata_long <- metadata.HNEG.test %>%
  pivot_longer(cols = c(Palsa.1, Bog.1, Fen.1),
               names_to = "category",
               values_to = "percentage")

metadata_long <- metadata_long %>%
  mutate(
    Habitat_plant = factor(Habitat_plant, levels = c("RT", "ET", "PT")),  # Set order for Habitat_plant
    Sample_code = factor(Sample_code, levels = c("core _1", "core _2", "core _3")),  # Set order for Sample_code
    combined = factor(
      paste(Habitat_plant, Sample_code, sep = "."),  # Combine Habitat_plant and Sample_code
      levels = c(
        "RT.core _1", "RT.core _2", "RT.core _3",
        "ET.core _1", "ET.core _2", "ET.core _3",
        "PT.core _1", "PT.core _2", "PT.core _3"
      )
    ),
    category = factor(category, levels = c("Fen.1", "Bog.1", "Palsa.1"))  # Order for stacking
  )

# Step 4: Create the stacking bar plot
result.HNEG.lr.sum <- metadata_long %>%
  ggplot(aes(x = combined, y = percentage, fill = category)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ depth, ncol = 1) +  # Facet by depth, align vertically
  scale_fill_manual(values = c("Fen.1" = "#4273B0", "Bog.1" = "#5E813F", "Palsa.1" = "#B8912F")) +  # Customize fill colors
  labs(
    title = "HNEG: %Correct Predictions",
    x = "Habitat (Sample Code)",
    y = "Percentage",
    fill = "Category"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

#RPPOS --------------------

metadata.RPPOS <- metadata
#use autochamber samples as training set
metadata.RPPOS.ac <- subset(metadata.RPPOS, grepl("Autochamber",metadata.RPPOS$site))
#use avni and inoc as test set
metadata.RPPOS.test <- subset(metadata.RPPOS, !rownames(metadata.RPPOS)%in% rownames(metadata.RPPOS.ac))


rows_to_move <- metadata.RPPOS.test %>%
  filter(Habitat_plant %in% c("AvCj_P", "IncE", "IncS"))

metadata.RPPOS.ac <- bind_rows(metadata.RPPOS.ac, rows_to_move)
#use avni and inoc as test set
metadata.RPPOS.test <- metadata.RPPOS.test %>%
  filter(!row.names(metadata.RPPOS.test) %in% row.names(rows_to_move))

RPPOS.ac <- RPPOS[rownames(metadata.RPPOS.ac),]
metadata.RPPOS.ac <- metadata.RPPOS[rownames(RPPOS.ac),]
RPPOS.ac.pareto <- na.omit(RPPOS.ac) %>%
  log_transform() %>%
  pareto_scale()

set.seed(123)
train_control <- trainControl(method = "repeatedcv", number = 10, repeats=3, classProbs = TRUE, summaryFunction = multiClassSummary)

RPPOS.all <- data.frame(habitat = metadata.RPPOS.ac$habitat, depth = metadata.RPPOS.ac$depth, RPPOS.ac.pareto)

p = 0.7
RPPOS.all$composite_var <- interaction(RPPOS.all$habitat, RPPOS.all$depth)
index <- createDataPartition(RPPOS.all$composite_var, p = p, list = FALSE)
trainData <- RPPOS.all[index, !names(RPPOS.all) %in% c("depth", "composite_var")]
testData <- RPPOS.all[-index, !names(RPPOS.all) %in% c("depth", "composite_var")]



rp.model <- train(habitat ~ ., 
                  data = trainData, 
                  metric = "Accuracy",
                  method = "glmnet", 
                  trControl = train_control,
                  tuneGrid = expand.grid(alpha = 0.5,  
                                         lambda = 10^seq(-3, 1, length = 20)))
print(rp.model)
#lambda        logLoss     AUC        prAUC      Accuracy   Kappa      Mean_F1    Mean_Sensitivity  Mean_Specificity  Mean_Pos_Pred_Value
#  0.029763514  0.06741333  0.9997082  0.8717654  0.9912330  0.9868356  0.9911446  0.9912698         0.9956481             0.9958633

predictions <- predict(rp.model, newdata = testData)

# Confusion matrix for test set performance
confusionMatrix(predictions, testData$habitat)

result_table <- rp.model$results

best_model <- rp.model$finalModel
best_lambda <- rp.model$bestTune$lambda 
best_alpha <- rp.model$bestTune$alpha


metrics <- result_table[, c("lambda", "logLoss", "Accuracy")]

# Melt the data for ggplot2 (long format)
metrics_long <- melt(metrics, id.vars = "lambda", 
                     variable.name = "Metric", 
                     value.name = "Value")

# Plot metrics vs. lambda
rp_tuning<-ggplot(metrics_long, aes(x = lambda, y = Value, color = Metric)) +
  geom_line() +
  geom_point() +
  scale_x_log10() +  # Log scale for lambda
  geom_vline(xintercept = best_lambda, color = "red", linetype = "dashed") +  # Add vertical line
  labs(title = paste("RP Tuning (Best Lambda =", signif(best_lambda, digits = 4), ")"),
       x = "Lambda (log scale)",
       y = "Metric Value") +
  theme_classic() +
  theme(legend.position = "bottom") +
  scale_color_brewer(palette = "Set1")





## model
RPPOS.lr <- glmnet(RPPOS.all[!names(RPPOS.all) %in% c("habitat","depth", "composite_var")], RPPOS.all$habitat, alpha = best_alpha, lambda = best_lambda, family = "multinomial")
coef_matrix_rp <- coef(RPPOS.lr, s = best_lambda)

# Convert the coefficients to a data frame for easier interpretation
coef_df_rp <- do.call(cbind, lapply(coef_matrix_rp, as.matrix))
colnames(coef_df_hn) <- c( "Palsa","Bog", "Fen")

intercepts <- sapply(c("Palsa", "Bog", "Fen"), function(habitat) {
  if (habitat %in% names(coef_matrix_hn)) {
    coef_matrix_hn[[habitat]]["(Intercept)", ]
  } else {
    NA  # Return NA if habitat is not found
  }
}, simplify = FALSE)

# Convert to a dataframe
intercepts_df <- do.call(rbind, intercepts)
rownames(intercepts_df) <- c("Palsa", "Bog", "Fen")



# run test
RPPOS.test.pareto <- RPPOS[rownames(metadata.RPPOS.test),] %>%
  log_transform() %>%
  pareto_scale()%>%
  as.matrix(.)


predicted_p.RPPOS <- predict(RPPOS.lr, newx = RPPOS.test.pareto, s = best_lambda, type = "response")
predicted_c.RPPOS <- predict(RPPOS.lr, newx = RPPOS.test.pareto, s = best_lambda, type = "class")

df <- as.data.frame(predicted_p.RPPOS)

metadata.RPPOS.test <- subset(metadata.RPPOS, !rownames(metadata.RPPOS)%in% rownames(metadata.RPPOS.ac))
metadata.RPPOS.test$predict <- predicted_c.RPPOS
metadata.RPPOS.test <- cbind (metadata.RPPOS.test,df)

metadata.RPPOS.test <- metadata.RPPOS.test %>%
  mutate(Habitat_plant = case_when(
    Habitat_plant == "AvCe_F1" ~ "PT",
    Habitat_plant == "AvCe_C" ~ "ET",
    Habitat_plant == "AvC_P" ~ "RT",
    TRUE ~ Habitat_plant  # Keep other values unchanged
  ))

# Step 2: Pivot longer
metadata_long <- metadata.RPPOS.test %>%
  pivot_longer(cols = c(Palsa.1, Bog.1, Fen.1),
               names_to = "category",
               values_to = "percentage")

metadata_long <- metadata_long %>%
  mutate(
    Habitat_plant = factor(Habitat_plant, levels = c("RT", "ET", "PT")),  # Set order for Habitat_plant
    Sample_code = factor(Sample_code, levels = c("core _1", "core _2", "core _3")),  # Set order for Sample_code
    combined = factor(
      paste(Habitat_plant, Sample_code, sep = "."),  # Combine Habitat_plant and Sample_code
      levels = c(
        "RT.core _1", "RT.core _2", "RT.core _3",
        "ET.core _1", "ET.core _2", "ET.core _3",
        "PT.core _1", "PT.core _2", "PT.core _3"
      )
    ),
    category = factor(category, levels = c("Fen.1", "Bog.1", "Palsa.1"))  # Order for stacking
  )

# Step 4: Create the stacking bar plot
result.RPPOS.lr.sum <- metadata_long %>%
  ggplot(aes(x = combined, y = percentage, fill = category)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ depth, ncol = 1) +  # Facet by depth, align vertically
  scale_fill_manual(values = c("Fen.1" = "#4273B0", "Bog.1" = "#5E813F", "Palsa.1" = "#B8912F")) +  # Customize fill colors
  labs(
    title = "RPPOS: %Correct Predictions",
    x = "Habitat (Sample Code)",
    y = "Percentage",
    fill = "Category"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

#binded
metadata.RPPOS.test <- metadata.RPPOS.test %>%
  mutate(source = "RPPOS")

metadata.HNEG.test <- metadata.HNEG.test %>%
  mutate(source = "HNEG")

# Combine the two datasets
metadata_combined <- rbind(metadata.RPPOS.test, metadata.HNEG.test)

# Pivot longer to prepare for plotting
metadata_long_combined <- metadata_combined %>%
  pivot_longer(cols = c(Palsa.1, Bog.1, Fen.1),
               names_to = "category",
               values_to = "percentage")

# Update factor levels and combined variable
metadata_long_combined <- metadata_long_combined %>%
  mutate(
    Habitat_plant = factor(Habitat_plant, levels = c("RT", "ET", "PT")),
    Sample_code = factor(Sample_code, levels = c("core _1", "core _2", "core _3")),
    combined = factor(
      paste(Habitat_plant, Sample_code, sep = "."),
      levels = c(
        "RT.core _1", "RT.core _2", "RT.core _3",
        "ET.core _1", "ET.core _2", "ET.core _3",
        "PT.core _1", "PT.core _2", "PT.core _3"
      )
    ),
    category = factor(category, levels = c("Palsa.1", "Bog.1","Fen.1"))
  )

metadata_long_combined <- metadata_long_combined %>%
  mutate(combined_all = paste(combined, source,sep = "-"))

# Create a stacked bar plot
combined_plot <- metadata_long_combined %>%
  ggplot(aes(x = combined_all, y = percentage, fill = category)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ depth, ncol = 1) +  # Facet by source and depth
  scale_fill_manual(values = c("Fen.1" = "#4273B0", "Bog.1" = "#5E813F", "Palsa.1" = "#B8912F")) +
  labs(
    title = "%Correct Predictions for RPPOS and HNEG",
    x = "Habitat (Sample Code)",
    y = "Percentage",
    fill = "Category"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


ggsave("prediction.pdf",combined_plot, width = 4, height = 10)


write.table(coef_df_rp, file = "Important_LR_RP_sig_addinc.csv", sep = ",")
write.table(coef_df_hn, file = "Important_LR_hn_sig_addinc.csv", sep = ",")

combined_plot <- ggarrange(rp_tuning, hilic_tuning, 
                           ncol = 2, nrow = 1, 
                           common.legend = TRUE, legend = "bottom")

# Save the combined plot
ggsave("tuning_plot.pdf", plot = combined_plot, width = 5, height = 3)