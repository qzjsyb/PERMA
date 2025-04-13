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
library(ggrepel)
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

#HN_sig <- read.csv("HN_shared_features_all.csv", stringsAsFactors = FALSE)

metadata.HNEG <- metadata[rownames(metadata) != "2014_S_2_20_24", ]
#use autochamber samples as training set
metadata.HNEG.ac <- subset(metadata.HNEG, grepl("Autochamber",metadata.HNEG$site))
#use avni and inoc as test set
metadata.HNEG.avni <- subset(metadata.HNEG, !rownames(metadata.HNEG)%in% rownames(metadata.HNEG.ac))

rows_to_move <- metadata.HNEG.avni %>%
  filter(Habitat_plant %in% c("AvCj_P", "IncE", "IncS"))

metadata.HNEG.ac <- bind_rows(metadata.HNEG.ac, rows_to_move)
#use avni and inoc as test set
metadata.HNEG.avni <- metadata.HNEG.avni %>%
  filter(!row.names(metadata.HNEG.avni) %in% row.names(rows_to_move))


HNEG.ac <- HNEG[rownames(metadata.HNEG.ac),]
metadata.HNEG.ac <- metadata.HNEG[rownames(HNEG.ac),]
HNEG.ac.pareto <- na.omit(HNEG.ac) %>%
  log_transform() %>%
  pareto_scale()

HNEG.ac.tic <- na.omit(HNEG.ac)  
HNEG.ac.tic <- HNEG.ac.tic / rowSums(HNEG.ac.tic)

HNEG.avni <- HNEG[rownames(metadata.HNEG.avni),]
HNEG.avni.pareto <- na.omit(HNEG.avni) %>%
  log_transform() %>%
  pareto_scale()

HNEG.avni.tic <- na.omit(HNEG.avni)  
HNEG.avni.tic <- HNEG.avni.tic / rowSums(HNEG.avni.tic)

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
metadata.RPPOS.avni <- subset(metadata.RPPOS, !rownames(metadata.RPPOS)%in% rownames(metadata.RPPOS.ac))


rows_to_move <- metadata.RPPOS.avni %>%
  filter(Habitat_plant %in% c("AvCj_P", "IncE", "IncS"))

metadata.RPPOS.ac <- bind_rows(metadata.RPPOS.ac, rows_to_move)
#use avni and inoc as test set
metadata.RPPOS.avni <- metadata.RPPOS.avni %>%
  filter(!row.names(metadata.RPPOS.avni) %in% row.names(rows_to_move))


RPPOS.ac <- RPPOS[rownames(metadata.RPPOS.ac),]
metadata.RPPOS.ac <- metadata.RPPOS[rownames(RPPOS.ac),]
RPPOS.ac.pareto <- na.omit(RPPOS.ac) %>%
  log_transform() %>%
  pareto_scale()

RPPOS.avni <- RPPOS[rownames(metadata.RPPOS.avni),]
RPPOS.avni.pareto <- na.omit(RPPOS.avni) %>%
  log_transform() %>%
  pareto_scale()


RPPOS.ac.tic <- na.omit(RPPOS.ac)  
RPPOS.ac.tic <- RPPOS.ac.tic / rowSums(RPPOS.ac.tic)
RPPOS.ac.log1ptic <- log1p(RPPOS.ac.tic)

RPPOS.avni <- RPPOS[rownames(metadata.RPPOS.avni.clean),]
RPPOS.avni.pareto <- na.omit(RPPOS.avni) %>%
  log_transform() %>%
  pareto_scale()

RPPOS.avni.tic <- na.omit(RPPOS.avni)  
RPPOS.avni.tic <- RPPOS.avni.tic / rowSums(RPPOS.avni.tic)
RPPOS.avni.log1ptic <- log1p(RPPOS.avni.tic)

#calculate similarity between all autochamber.
library(coop)

calculate_mu_plot <- function(data, color_line = "blue", subset_name) {
  # Step 1: Calculate similarity matrix
  # distance_matrix <- dist(data, method = "euclidean")
  # min_dist <- min(distance_matrix)
  # max_dist <- max(distance_matrix)
  # 
  # similarity_matrix <- 1 - ((distance_matrix - min_dist) / (max_dist - min_dist))
  # similarity_matrix <- as.matrix(similarity_matrix)
  # diag(similarity_matrix) <- 1
  #try cosine similarity
  #similarity_matrix <- tcosine(data)
  
  #try bray-curtis
  distance_matrix <- vegdist(data, method = "bray")  # Bray-Curtis dissimilarity
  
  
  # Step 2: Convert to similarity matrix
  similarity_matrix <- 1 - as.matrix(distance_matrix)  # Bray-Curtis similarity (1 - dissimilarity)
  diag(similarity_matrix) <- 1  # Set diagonal to 1
  
  n <- nrow(similarity_matrix)
  
  # Step 2: Calculate affinity matrix
  affinity_matrix <- matrix(0.5, n, n, dimnames = list(rownames(similarity_matrix), colnames(similarity_matrix)))
  
  calculate_affinity <- function(similarity_matrix) {
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        d_k <- similarity_matrix[i, ] - similarity_matrix[j, ]  # Compute differences
        d_k <- d_k[-c(i, j)]  # Remove self-comparisons
        
        rank_d <- rank(abs(d_k))  # Rank based on absolute differences
        rank_d[d_k < 0] <- -rank_d[d_k < 0]  # Assign original signs
        
        rank_sum <- sum(rank_d[rank_d > 0])  # Sum of positive ranks
        max_rank_sum <- (n-2)*(n-1)/2  # Maximum possible rank sum
        
        A_ij <- rank_sum / max_rank_sum  # Normalize rank sum
        affinity_matrix[i, j] <- A_ij
        affinity_matrix[j, i] <- 1 - A_ij  # Complementary affinity
      }
    }
    return(affinity_matrix)
  }
  
  affinity_matrix <- calculate_affinity(similarity_matrix)
  
  diag(similarity_matrix) <- NA
  diag(affinity_matrix) <- NA
  
  # Step 3: Compute mean similarity and mean affinity
  mean_similarity <- rowMeans(similarity_matrix, na.rm = TRUE)
  mean_affinity <- rowMeans(affinity_matrix, na.rm = TRUE)
  
  # Step 4: Calculate mu
  correlation <- cor(mean_similarity, mean_affinity, use = "complete.obs")
  std_dev_affinity <- sd(mean_affinity)
  std_dev_similarity <- sd(mean_similarity)
  
  mu_all <- correlation * (std_dev_affinity / std_dev_similarity)
  
  # Prepare data for plotting
  data_plot <- data.frame(
    mean_similarity = mean_similarity,
    mean_affinity = mean_affinity,
    sample = rownames(data)
  )
  
  # Add depth from metadata (assuming 'metadata' is a global variable with depth)
  matched_rows <- match(rownames(data_plot), rownames(metadata))
  data_plot$depth <- metadata$depth[matched_rows]
  data_plot$Habitat_plant <- metadata$Habitat_plant[matched_rows]
  
  # Step 5: Linear model and regression statistics
  model <- lm(mean_affinity ~ mean_similarity, data = data_plot)
  
  slope <- coef(model)[2]
  r_squared <- summary(model)$r.squared
  
  
  # Jackknife loop: drop each row one at a time and recalculate mu
  jackknife_mu_values <- numeric(nrow(data))
  for (i in 1:nrow(data)) {
    # Remove the i-th site (row) from the data
    reduced_matrix <- data[-i, ]
    # 
    # distance_matrix <- dist(reduced_matrix, method = "euclidean")
    # min_dist <- min(distance_matrix)
    # max_dist <- max(distance_matrix)
    # 
    # similarity_matrix <- 1 - ((distance_matrix - min_dist) / (max_dist - min_dist))
    # similarity_matrix <- as.matrix(similarity_matrix)
    # diag(similarity_matrix) <- 1
    #similarity_matrix <- tcosine(reduced_matrix)
    
    distance_matrix <- vegdist(reduced_matrix, method = "bray")  # Bray-Curtis dissimilarity
    
    # Step 2: Convert to similarity matrix
    similarity_matrix <- 1 - as.matrix(distance_matrix)  # Bray-Curtis similarity (1 - dissimilarity)
    diag(similarity_matrix) <- 1  # Set diagonal to 1
    
    n <- nrow(similarity_matrix)
    affinity_matrix <- matrix(0.5, n, n)
    affinity_matrix <- calculate_affinity(similarity_matrix)
    
    diag(similarity_matrix) <- NA
    diag(affinity_matrix) <- NA
    
    # Compute mean similarity and mean affinity for the reduced matrix
    mean_similarity <- rowMeans(similarity_matrix, na.rm = TRUE)
    mean_affinity <- rowMeans(affinity_matrix, na.rm = TRUE)
    
    # Calculate mu (correlation * ratio of standard deviations)
    correlation <- cor(mean_similarity, mean_affinity, use = "complete.obs")
    std_dev_affinity <- sd(mean_affinity)
    std_dev_similarity <- sd(mean_similarity)
    mu <- correlation * (std_dev_affinity / std_dev_similarity)
    
    # Store the mu value for the current jackknife iteration
    jackknife_mu_values[i] <- mu
  }
  
  jackknife_sd <- sd(jackknife_mu_values)
    
  habitat_colors <- c("#5E813F","#4273B0","purple4","red4","red1", "#B8912F")
  names(habitat_colors) <-  c("BOG","FEN", "AvCe_F1", "AvCe_C", "AvC_P", "PALSA")
  
  
  # Step 6: Create plot
  plot <- ggplot(data_plot, aes(x = mean_similarity, y = mean_affinity)) + 
    geom_point(aes(shape = as.factor(depth), color = Habitat_plant), size = 3) +
    geom_smooth(method = "lm", se = TRUE, color = color_line, linewidth = 1) +
    geom_hline(yintercept = mean(data_plot$mean_affinity, na.rm = TRUE) - 
                 sd(data_plot$mean_affinity, na.rm = TRUE), 
               linetype = "dashed", color = "grey", size = 1) +
    labs(
      title = paste(subset_name, 
                    "(μ=", round(mu_all, 3), "), sd=", round(jackknife_sd, 3)),
      x = "Mean Similarity",
      y = "Mean Affinity"
    ) +
    scale_color_manual(values = habitat_colors) + 
    theme_classic() +
    annotate("text", x = max(data_plot$mean_similarity) * 0.3, 
             y = max(data_plot$mean_affinity) * 0.9,
             label = paste("Slope =", round(slope, 3), 
                           "\nR² =", round(r_squared, 2)), 
             size = 5, color = color_line) +
    geom_text_repel(
      data = subset(data_plot, 
                    mean_affinity < (mean(mean_affinity, na.rm = TRUE) - 
                                       sd(mean_affinity, na.rm = TRUE))),
      aes(label = sample),  # Use the sample column directly
      size = 4, box.padding = 0.5, max.overlaps = 15, na.rm = TRUE
    )
  
  return(list(mu = mu_all, sd = jackknife_sd, plot = plot))
}

#use TIC scaled -----
palsa_data <- RPPOS.ac.tic[metadata.RPPOS.ac$habitat == "Palsa", ]
bog_data <- RPPOS.ac.tic[metadata.RPPOS.ac$habitat == "Bog", ]
fen_data <- RPPOS.ac.tic[metadata.RPPOS.ac$habitat == "Fen", ]

#all.mu.rp <- calculate_mu_plot(RPPOS.ac.pareto.shallow, color_line = "grey", subset_name = "all RP")
palsa.mu.rp.bray.alldepth <- calculate_mu_plot(palsa_data, color_line = "grey", subset_name = "palsa RP")
#4.425895 +- 0.02142638
bog.mu.rp.bray.alldepth <-calculate_mu_plot(bog_data, color_line = "grey", subset_name = "bog RP")
#3.850449 +- 0.02658061
fen.mu.rp.bray.alldepth <-calculate_mu_plot(fen_data, color_line = "grey", subset_name = "fen RP")
#2.899028 +- 0.02692167



mu_results_rp.habitat <- data.frame(year = integer(), habitat = character(), mu = numeric(), sd = numeric(), count = integer(), stringsAsFactors = FALSE)

# Define the years, excluding 2013, because no deep samples, 2022, because no palsa
years <- c(2013:2020, 2022:2023)

# Loop through each year and each habitat use raw
for (year in years) {
  # Subset data by habitat
  palsa_year <- palsa_data[grep(year, rownames(palsa_data)), ]
  bog_year <- bog_data[grep(year, rownames(bog_data)), ]
  fen_year <- fen_data[grep(year, rownames(fen_data)), ]
  
  # Function to calculate mu if data is available
  calculate_if_exists <- function(data, habitat) {
    if (nrow(data) > 0) {
      mu_result <- calculate_mu_plot(data, color_line = "grey", subset_name = paste0(year, " ", habitat, " RP"))
      return(data.frame(year = year, habitat = habitat, mu = mu_result$mu, sd = mu_result$sd, count = nrow(data)))
    } else {
      return(NULL)  # Return NULL if no data
    }
  }
  
  # Apply function to each habitat and store results
  results_list <- list(
    calculate_if_exists(palsa_year, "palsa"),
    calculate_if_exists(bog_year, "bog"),
    calculate_if_exists(fen_year, "fen")
  )
  
  # Bind non-null results to the dataframe
  mu_results_rp.habitat <- rbind(mu_results_rp.habitat, do.call(rbind, results_list[!sapply(results_list, is.null)]))
}

avcp <- RPPOS.avni.tic[grep("AvC_P", rownames(RPPOS.avni.tic)), ]
avcp_mu <- calculate_mu_plot(avcp, color_line = "grey", subset_name = "avcp")

avcec <- RPPOS.avni.tic[grep("AvCe_C", rownames(RPPOS.avni.tic)), ]
avcec_mu <- calculate_mu_plot(avcec, color_line = "grey", subset_name = "avcec")

avcef <- RPPOS.avni.tic[grep("AvCe_F", rownames(RPPOS.avni.tic)), ]
avcef_mu <- calculate_mu_plot(avcef, color_line = "grey", subset_name = "avcef")

subset_data.p <- subset(metadata.RPPOS.ac, Year == 2016 & Habitat_plant == "PALSA")
subset_data.p <- subset(metadata.RPPOS.ac, Year == 2017 & Habitat_plant == "PALSA")

subset_data.b <- subset(metadata.RPPOS.ac, Year == 2016 & Habitat_plant == "BOG")

subset_data.f <- subset(metadata.RPPOS.ac, Year == 2016 & Habitat_plant == "FEN")

mu_results_rp.habitat <- rbind(mu_results_rp.habitat,
                               data.frame(year = 2023, habitat = "AvC_P", mu = avcp_mu$mu, sd = avcp_mu$sd, count = nrow(avcp)),
                               data.frame(year = 2023, habitat = "AvCe_C", mu = avcec_mu$mu, sd = avcec_mu$sd, count = nrow(avcec)),
                               data.frame(year = 2023, habitat = "AvCe_F1", mu = avcef_mu$mu, sd = avcef_mu$sd, count = nrow(avcef)))

habitat_colors <- c("#5E813F", "#4273B0", "purple4", "red4", "red1", "#B8912F")
names(habitat_colors) <- c("bog", "fen", "AvCe_F1", "AvCe_C", "AvC_P", "palsa")

habitat_shades <- data.frame(
  habitat = c("palsa", "bog", "fen"),
  mu = c(palsa.mu.rp.bray.alldepth$mu, bog.mu.rp.bray.alldepth$mu, fen.mu.rp.bray.alldepth$mu),
  sd = c(palsa.mu.rp.bray.alldepth$sd, bog.mu.rp.bray.alldepth$sd, fen.mu.rp.bray.alldepth$sd)
)

habitat_shades <- habitat_shades %>%
  mutate(
    ymin = mu - sd,
    ymax = mu + sd
  )

rp_mu_plot_by_habitat <- ggplot(mu_results_rp.habitat, aes(x = year, y = mu, color = habitat, group = habitat)) +
  geom_rect(data = habitat_shades,
            aes(xmin = min(mu_results_rp.habitat$year),  # Adjust as per your data range
                xmax = max(mu_results_rp.habitat$year),  # Adjust as per your data range
                ymin = ymin,
                ymax = ymax,
                fill = habitat),
            alpha = 0.2) +  # Adjust transparency as needed+
  geom_line() +  # Lines for each habitat
  geom_point(aes(size = count / 3)) +  # Points for yearly mu
  scale_color_manual(values = habitat_colors) +  # Apply custom colors
  scale_fill_manual(values = habitat_colors)+  # Ensure 'habitat_colors' is defined
  scale_x_continuous(breaks = seq(min(mu_results_rp.habitat$year), max(mu_results_rp.habitat$year), by = 1)) + 
  scale_size_continuous(name = "#of Sample", labels = function(x) x * 3) +
  labs(color = "Habitat") +  # Add legend title for habitat
  theme_classic2()  # Optional, you can change the theme as desired


ggsave("rp_mu_plot_by_habitat.pdf", rp_mu_plot_by_habitat, width = 8, height = 4)


#HN, use TIC scaled, not included in main text ----

palsa_data <- HNEG.ac.tic[metadata.HNEG.ac$habitat == "Palsa", ]
bog_data <- HNEG.ac.tic[metadata.HNEG.ac$habitat == "Bog", ]
fen_data <- HNEG.ac.tic[metadata.HNEG.ac$habitat == "Fen", ]

#all.mu.hn <- calculate_mu_plot(HNEG.ac.pareto.shallow, color_line = "grey", subset_name = "all RP")
palsa.mu.hn.bray.alldepth <- calculate_mu_plot(palsa_data, color_line = "grey", subset_name = "palsa RP")
#3.339614 +- 0.03218098
bog.mu.hn.bray.alldepth <-calculate_mu_plot(bog_data, color_line = "grey", subset_name = "bog RP")
#3.599167 +- 0.02025399
fen.mu.hn.bray.alldepth <-calculate_mu_plot(fen_data, color_line = "grey", subset_name = "fen RP")
#2.529386 +- 0.02337224



mu_results_hn.habitat <- data.frame(year = integer(), habitat = character(), mu = numeric(), sd = numeric(), count = integer(), stringsAsFactors = FALSE)

# Define the years, excluding 2013, because no deep samples, 2022, because no palsa
years <- c(2013:2020, 2022:2023)

# Loop through each year and each habitat
for (year in years) {
  # Subset data by habitat
  palsa_year <- palsa_data[grep(year, rownames(palsa_data)), ]
  bog_year <- bog_data[grep(year, rownames(bog_data)), ]
  fen_year <- fen_data[grep(year, rownames(fen_data)), ]
  
  # Function to calculate mu if data is available
  calculate_if_exists <- function(data, habitat) {
    if (nrow(data) > 0) {
      mu_result <- calculate_mu_plot(data, color_line = "grey", subset_name = paste0(year, " ", habitat, " RP"))
      return(data.frame(year = year, habitat = habitat, mu = mu_result$mu, sd = mu_result$sd,count = nrow(data)))
    } else {
      return(NULL)  # Return NULL if no data
    }
  }
  
  # Apply function to each habitat and store results
  results_list <- list(
    calculate_if_exists(palsa_year, "palsa"),
    calculate_if_exists(bog_year, "bog"),
    calculate_if_exists(fen_year, "fen")
  )
  
  # Bind non-null results to the dataframe
  mu_results_hn.habitat <- rbind(mu_results_hn.habitat, do.call(rbind, results_list[!sapply(results_list, is.null)]))
}

avcp <- HNEG.avni.tic[grep("AvC_P", rownames(HNEG.avni.tic)), ]
avcp_mu <- calculate_mu_plot(avcp, color_line = "grey", subset_name = "avcp")

avcec <- HNEG.avni.tic[grep("AvCe_C", rownames(HNEG.avni.tic)), ]
avcec_mu <- calculate_mu_plot(avcec, color_line = "grey", subset_name = "avcec")

avcef <- HNEG.avni.tic[grep("AvCe_F", rownames(HNEG.avni.tic)), ]
avcef_mu <- calculate_mu_plot(avcef, color_line = "grey", subset_name = "avcef")

mu_results_hn.habitat <- rbind(mu_results_hn.habitat,
                               data.frame(year = 2023, habitat = "AvC_P", mu = avcp_mu$mu, sd = avcp_mu$sd, count = nrow(avcp)),
                               data.frame(year = 2023, habitat = "AvCe_C", mu = avcec_mu$mu, sd = avcec_mu$sd, count = nrow(avcec)),
                               data.frame(year = 2023, habitat = "AvCe_F1", mu = avcef_mu$mu, sd = avcef_mu$sd, count = nrow(avcef)))

habitat_colors <- c("#5E813F", "#4273B0", "purple4", "red4", "red1", "#B8912F")
names(habitat_colors) <- c("bog", "fen", "AvCe_F1", "AvCe_C", "AvC_P", "palsa")

habitat_shades <- data.frame(
  habitat = c("palsa", "bog", "fen"),
  mu = c(palsa.mu.hn.bray.alldepth$mu, bog.mu.hn.bray.alldepth$mu, fen.mu.hn.bray.alldepth$mu),
  sd = c(palsa.mu.hn.bray.alldepth$sd, bog.mu.hn.bray.alldepth$sd, fen.mu.hn.bray.alldepth$sd)
)

habitat_shades <- habitat_shades %>%
  mutate(
    ymin = mu - sd,
    ymax = mu + sd
  )

hn_mu_plot_by_habitat <- ggplot(mu_results_hn.habitat, aes(x = year, y = mu, color = habitat, group = habitat)) +
  geom_rect(data = habitat_shades,
            aes(xmin = min(mu_results_hn.habitat$year),  # Adjust as per your data range
                xmax = max(mu_results_hn.habitat$year),  # Adjust as per your data range
                ymin = ymin,
                ymax = ymax,
                fill = habitat),
            alpha = 0.2) +  # Adjust transparency as needed+
  geom_line() +  # Lines for each habitat
  geom_point(aes(size = count / 3)) +  # Points for yearly mu
  scale_color_manual(values = habitat_colors) +  # Apply custom colors
  scale_fill_manual(values = habitat_colors)+  # Ensure 'habitat_colors' is defined
  labs(color = "Habitat") +  # Add legend title for habitat
  theme_minimal()  # Optional, you can change the theme as desired


#used 2023 deep and shallow data.----
palsa_2023_deep <- RPPOS.ac.tic[metadata.RPPOS.ac$habitat == "Palsa" & metadata.RPPOS.ac$Year == 2023 & grepl("20_24|30_34", rownames(RPPOS.ac.tic)), ]
bog_2023_deep <- RPPOS.ac.tic[metadata.RPPOS.ac$habitat == "Bog"  & metadata.RPPOS.ac$Year == 2023 & grepl("20_24|30_34", rownames(RPPOS.ac.tic)), ]
fen_2023_deep <- RPPOS.ac.tic[metadata.RPPOS.ac$habitat == "Fen" & metadata.RPPOS.ac$Year == 2023 & grepl("20_24|30_34", rownames(RPPOS.ac.tic)), ]

palsa_2023_shallow <- RPPOS.ac.tic[metadata.RPPOS.ac$habitat == "Palsa" & metadata.RPPOS.ac$Year == 2023 & grepl("1_5|10_14", rownames(RPPOS.ac.tic)), ]
bog_2023_shallow <- RPPOS.ac.tic[metadata.RPPOS.ac$habitat == "Bog"  & metadata.RPPOS.ac$Year == 2023 & grepl("1_5|10_14", rownames(RPPOS.ac.tic)), ]
fen_2023_shallow <- RPPOS.ac.tic[metadata.RPPOS.ac$habitat == "Fen" & metadata.RPPOS.ac$Year == 2023 & grepl("1_5|10_14", rownames(RPPOS.ac.tic)), ]

palsa_2023 <- RPPOS.ac.tic[metadata.RPPOS.ac$habitat == "Palsa" & metadata.RPPOS.ac$Year == 2023, ]
bog_2023 <- RPPOS.ac.tic[metadata.RPPOS.ac$habitat == "Bog"  & metadata.RPPOS.ac$Year == 2023, ]
fen_2023 <- RPPOS.ac.tic[metadata.RPPOS.ac$habitat == "Fen" & metadata.RPPOS.ac$Year == 2023, ]



#all.mu.rp <- calculate_mu_plot(RPPOS.ac.pareto.deep, color_line = "grey", subset_name = "all RP")
palsa.mu.rp.bray.deep.2023 <- calculate_mu_plot(palsa_2023_deep, color_line = "grey", subset_name = "palsa RP")
bog.mu.rp.bray.deep.2023 <-calculate_mu_plot(bog_2023_deep, color_line = "grey", subset_name = "bog RP")
fen.mu.rp.bray.deep.2023 <-calculate_mu_plot(fen_2023_deep, color_line = "grey", subset_name = "fen RP")

palsa.mu.rp.bray.shallow.2023 <- calculate_mu_plot(palsa_2023_shallow, color_line = "grey", subset_name = "palsa RP")
bog.mu.rp.bray.shallow.2023 <-calculate_mu_plot(bog_2023_shallow, color_line = "grey", subset_name = "bog RP")
fen.mu.rp.bray.shallow.2023 <-calculate_mu_plot(fen_2023_shallow, color_line = "grey", subset_name = "fen RP")

palsa.mu.rp.bray.all.2023 <- calculate_mu_plot(palsa_2023, color_line = "grey", subset_name = "palsa RP")
bog.mu.rp.bray.all.2023 <-calculate_mu_plot(bog_2023, color_line = "grey", subset_name = "bog RP")
fen.mu.rp.bray.all.2023 <-calculate_mu_plot(fen_2023, color_line = "grey", subset_name = "fen RP")




avcp.deep <- RPPOS.avni.tic[grepl("AvC_P", rownames(RPPOS.avni.tic)) & grepl("20_24|30_34", rownames(RPPOS.avni.tic)), ]
avcp.deep_mu <- calculate_mu_plot(avcp.deep, color_line = "grey", subset_name = "avcp")

avcec.deep <- RPPOS.avni.tic[grepl("AvCe_C", rownames(RPPOS.avni.tic)) & grepl("20_24|30_34", rownames(RPPOS.avni.tic)), ]
avcec.deep_mu <- calculate_mu_plot(avcec.deep, color_line = "grey", subset_name = "avcec")

avcef.deep <- RPPOS.avni.tic[grepl("AvCe_F", rownames(RPPOS.avni.tic)) & grepl("20_24|30_34", rownames(RPPOS.avni.tic)), ]
avcef.deep_mu <- calculate_mu_plot(avcef.deep, color_line = "grey", subset_name = "avcef")

avcp.shallow <- RPPOS.avni.tic[grepl("AvC_P", rownames(RPPOS.avni.tic)) & grepl("1_5|10_14", rownames(RPPOS.avni.tic)), ]
avcp.shallow_mu <- calculate_mu_plot(avcp.shallow, color_line = "grey", subset_name = "avcp")

avcec.shallow <- RPPOS.avni.tic[grepl("AvCe_C", rownames(RPPOS.avni.tic)) & grepl("1_5|10_14", rownames(RPPOS.avni.tic)), ]
avcec.shallow_mu <- calculate_mu_plot(avcec.shallow, color_line = "grey", subset_name = "avcec")

avcef.shallow <- RPPOS.avni.tic[grepl("AvCe_F", rownames(RPPOS.avni.tic)) & grepl("1_5|10_14", rownames(RPPOS.avni.tic)), ]
avcef.shallow_mu <- calculate_mu_plot(avcef.shallow, color_line = "grey", subset_name = "avcef")


avcp <- RPPOS.avni.tic[grep("AvC_P", rownames(RPPOS.avni.tic)), ]
avcp_mu <- calculate_mu_plot(avcp, color_line = "grey", subset_name = "avcp")

avcec <- RPPOS.avni.tic[grep("AvCe_C", rownames(RPPOS.avni.tic)), ]
avcec_mu <- calculate_mu_plot(avcec, color_line = "grey", subset_name = "avcec")

avcef <- RPPOS.avni.tic[grep("AvCe_F", rownames(RPPOS.avni.tic)), ]
avcef_mu <- calculate_mu_plot(avcef, color_line = "grey", subset_name = "avcef")

mu_data_2023_shallow_deep <- data.frame(
  habitat = c("palsa", "bog", "fen", "palsa", "bog", "fen", "palsa", "bog", "fen", 
              "AvC_P", "AvCe_C", "AvCe_F1", "AvC_P", "AvCe_C", "AvCe_F1","AvC_P", "AvCe_C", "AvCe_F1"),
  mu = c(palsa.mu.rp.bray.deep.2023$mu, bog.mu.rp.bray.deep.2023$mu, fen.mu.rp.bray.deep.2023$mu, 
         palsa.mu.rp.bray.shallow.2023$mu, bog.mu.rp.bray.shallow.2023$mu, fen.mu.rp.bray.shallow.2023$mu,
         palsa.mu.rp.bray.all.2023$mu, bog.mu.rp.bray.all.2023$mu, fen.mu.rp.bray.all.2023$mu,
         avcp.deep_mu$mu, avcec.deep_mu$mu, avcef.deep_mu$mu, 
         avcp.shallow_mu$mu, avcec.shallow_mu$mu, avcef.shallow_mu$mu,
         avcp_mu$mu, avcec_mu$mu, avcef_mu$mu),
  sd = c(palsa.mu.rp.bray.deep.2023$sd, bog.mu.rp.bray.deep.2023$sd, fen.mu.rp.bray.deep.2023$sd, 
         palsa.mu.rp.bray.shallow.2023$sd, bog.mu.rp.bray.shallow.2023$sd, fen.mu.rp.bray.shallow.2023$sd,
         palsa.mu.rp.bray.all.2023$sd, bog.mu.rp.bray.all.2023$sd, fen.mu.rp.bray.all.2023$sd,
         avcp.deep_mu$sd, avcec.deep_mu$sd, avcef.deep_mu$sd, 
         avcp.shallow_mu$sd, avcec.shallow_mu$sd, avcef.shallow_mu$sd,
         avcp_mu$sd, avcec_mu$sd, avcef_mu$sd
         ),
  depth = rep(c("deep", "shallow", "all"), times = 2, each = 3)
)

ggplot(mu_data_2023_shallow_deep, aes(x = habitat, y = mu, color = habitat, shape = depth)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mu - sd, ymax = mu + sd), width = 0.2) +
  scale_color_manual(values = habitat_colors) +
  scale_shape_manual(values = c(1, 16, 17)) +  # 16 = filled circle (shallow), 17 = filled triangle (deep)
  theme_minimal() +
  labs(y = "Mu ± SD", x = "Habitat", title = "Comparison of Mu ± SD Across Groups")

shallowdeepmu <- ggplot(mu_data_2023_shallow_deep, aes(x = factor(habitat, levels = c("palsa", "bog", "AvC_P", "AvCe_C", "AvCe_F1", "fen")), 
                                      y = mu, fill = habitat)) +
  geom_line(aes(group = depth), size = 0.8, color = "black") +  # Black connecting line by depth
  geom_errorbar(aes(ymin = mu - sd, ymax = mu + sd), width = 0.1, color = "black") +  # Black error bars
  geom_point(size = 3, shape = 21, stroke = 1, color = "black") +  # Points with habitat fill and black outline
  scale_fill_manual(values = habitat_colors) +  # Use fill instead of color
  theme_classic() +
  facet_grid(rows = vars(factor(depth, levels = c("all", "shallow", "deep"))), scales = "free_y", switch = "x") +  # Correct facet order
  labs(y = "Mu ± SD", x = "Habitat", title = "Comparison of Mu ± SD Across Groups")


ggsave("shallowdeepmu.pdf", shallowdeepmu, width = 6, height = 3)
