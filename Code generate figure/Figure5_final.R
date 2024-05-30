library(ggplot2)
library(dplyr)
library(readr)

##f


# Read data
data <- read_csv("~/Documents/TACIT figure/dataset/Figure 6/PCF_Xenium_reg004.csv")

# Filter data based on X and Y coordinates
data <- data %>%
  filter(Y < 2800 & Y > 2000) %>%
  filter(X < 3700 & X > 3200)

#Louvain


# Define groups based on CD247 expression
Group1 <- ifelse(as.numeric(data$CD247) > 0, "CD247+", "CD247-")
Group2 <- ifelse(data$TACIT_Seurat_Xenium_PCF %in% c("CD8+ Effector T Cells", "CD4+ T Cells",
                                                     "CD8+ Exhausted T Cells", "Capillaries", "Venules", 
                                                     "M1 Macrophages", "M2 Macrophages", "B Cells", 
                                                     "Dendritic Cells", "Neutrophil"), 
                 data$TACIT_Seurat_Xenium_PCF, "Others")

group_ct_st_1 <- paste0(Group2, "::", Group1)
round(table(group_ct_st_1) / length(group_ct_st_1), 4)

# Define groups based on PDCD1 expression
Group1 <- ifelse(as.numeric(data$PDCD1) > 0, "CD274+", "CD274-")
Group2 <- ifelse(data$TACIT_Seurat_Xenium_PCF %in% c("CD8+ Effector T Cells", "CD4+ T Cells",
                                                     "CD8+ Exhausted T Cells", "Capillaries", "Venules", 
                                                     "M1 Macrophages", "M2 Macrophages", "B Cells", 
                                                     "Dendritic Cells", "Neutrophil"), 
                 data$TACIT_Seurat_Xenium_PCF, "Others")
Group2 <- ifelse(Group2 %in% c("CD8+ Effector T Cells", "CD4+ T Cells", "CD8+ Exhausted T Cells"), "T Cells", Group2)
Group2 <- ifelse(Group2 %in% c("M1 Macrophages", "M2 Macrophages", "Neutrophil"), "Mono/Mac", Group2)
Group2 <- ifelse(Group2 %in% c("Capillaries", "Venules"), "Capillaries & Venules", Group2)
Group2 <- ifelse(Group2 %in% c("T Cells", "Dendritic Cells", "Capillaries & Venules"), Group2, "Others")

group_ct_st_2 <- paste0(Group2, "::", Group1)

# Combine groups
Group_combine <- mapply(function(ct_st_1, ct_st_2) {
  if (ct_st_1 == "Others::PDCD1+" & ct_st_2 == "Others::CD274+" & ct_st_1 == "Others::PDCD1-" & ct_st_2 == "Others::CD274-") {
    return("Others")
  } else if (ct_st_1 != "Others::PDCD1+" & ct_st_2 == "Others::CD274+" & ct_st_1 != "Others::PDCD1-" & ct_st_2 == "Others::CD274-") {
    return(ct_st_1)
  } else if (ct_st_1 == "Others::PDCD1+" & ct_st_2 != "Others::CD274+" & ct_st_1 == "Others::PDCD1-" & ct_st_2 != "Others::CD274-") {
    return(ct_st_2)
  } else if (ct_st_1 != "Others::PDCD1+" & ct_st_2 != "Others::CD274+" & ct_st_1 != "Others::PDCD1-" & ct_st_2 != "Others::CD274-") {
    return(paste0(ct_st_1, "::", ct_st_2))
  } else {
    return("Others")
  }
}, group_ct_st_1, group_ct_st_2)

# Apply changes using gsub
Group_combine <- gsub("CD4\\+ T Cells::PDCD1- ::T Cells::CD274-", "CD4+ T Cells", Group_combine)
Group_combine <- gsub("CD4\\+ T Cells::PDCD1\\+ ::T Cells::CD274\\+", "CD4+ T Cells::PDCD1+::CD274+", Group_combine)
Group_combine <- gsub("Mono/Mac::PDCD1- ::Others::CD274-", "Mono/Mac", Group_combine)
Group_combine <- gsub("Others::PDCD1- ::Capillaries & Venules::CD274-", "Others", Group_combine)
Group_combine <- gsub("Others::PDCD1- ::Dendritic Cells::CD274-", "Dendritic Cells", Group_combine)
Group_combine <- gsub("Others::PDCD1- ::Others::CD274-", "Others", Group_combine)
Group_combine <- gsub("Others::PDCD1\\+ ::Others::CD274\\+", "Others", Group_combine)
Group_combine <- gsub("T Cells::PDCD1- ::T Cells::CD274-", "T Cells", Group_combine)
Group_combine <- gsub("T Cells::PDCD1\\+ ::T Cells::CD274\\+", "T Cells::PDCD1+::CD274+", Group_combine)

Group_combine <- ifelse(data$TACIT_Seurat_Xenium_PCF == "B Cells", "B Cells", Group_combine)

# Calculate the number of observations for each group
group_counts <- as.data.frame(table(Group_combine))
colnames(group_counts) <- c("group_ct_st", "count")

# Merge the counts back into the original data
data_plot_ct <- data.frame(X = data$X, Y = data$Y, group_ct_st = Group_combine)
data_plot_ct <- merge(data_plot_ct, group_counts, by = "group_ct_st")

# Create a new column with the group names and counts
data_plot_ct$group_with_count <- paste(data_plot_ct$group_ct_st, "(n =", data_plot_ct$count, ")")

round(table(data_plot_ct$group_with_count) / length(data_plot_ct$group_with_count), 4)

# Plot the data
ggplot(data_plot_ct, aes(x = X, y = Y)) +
  geom_point(aes(color = as.factor(group_with_count)), size = 1) +
  theme_classic(base_size = 15) +
  guides(color = guide_legend(override.aes = list(size = 3), title = "Group"))





###TACIT


# Define groups based on CD247 expression
Group1 <- ifelse(as.numeric(data$CD247) > 0, "CD247+", "CD247-")
Group2 <- ifelse(data$TACIT_PCF_Xenium %in% c("CD8+ Effector T Cells", "CD4+ T Cells",
                                                     "CD8+ Exhausted T Cells", "Capillaries", "Venules", 
                                                     "M1 Macrophages", "M2 Macrophages", "B Cells", 
                                                     "Dendritic Cells", "Neutrophil"), 
                 data$TACIT_PCF_Xenium, "Others")

group_ct_st_1 <- paste0(Group2, "::", Group1)
round(table(group_ct_st_1) / length(group_ct_st_1), 4)

# Define groups based on PDCD1 expression
Group1 <- ifelse(as.numeric(data$PDCD1) > 0, "CD274+", "CD274-")
Group2 <- ifelse(data$TACIT_PCF_Xenium %in% c("CD8+ Effector T Cells", "CD4+ T Cells",
                                                     "CD8+ Exhausted T Cells", "Capillaries", "Venules", 
                                                     "M1 Macrophages", "M2 Macrophages", "B Cells", 
                                                     "Dendritic Cells", "Neutrophil"), 
                 data$TACIT_PCF_Xenium, "Others")
Group2 <- ifelse(Group2 %in% c("CD8+ Effector T Cells", "CD4+ T Cells", "CD8+ Exhausted T Cells"), "T Cells", Group2)
Group2 <- ifelse(Group2 %in% c("M1 Macrophages", "M2 Macrophages", "Neutrophil"), "Mono/Mac", Group2)
Group2 <- ifelse(Group2 %in% c("Capillaries", "Venules"), "Capillaries & Venules", Group2)
Group2 <- ifelse(Group2 %in% c("T Cells", "Dendritic Cells", "Capillaries & Venules"), Group2, "Others")

group_ct_st_2 <- paste0(Group2, "::", Group1)

# Combine groups
Group_combine <- mapply(function(ct_st_1, ct_st_2) {
  if (ct_st_1 == "Others::PDCD1+" & ct_st_2 == "Others::CD274+" & ct_st_1 == "Others::PDCD1-" & ct_st_2 == "Others::CD274-") {
    return("Others")
  } else if (ct_st_1 != "Others::PDCD1+" & ct_st_2 == "Others::CD274+" & ct_st_1 != "Others::PDCD1-" & ct_st_2 == "Others::CD274-") {
    return(ct_st_1)
  } else if (ct_st_1 == "Others::PDCD1+" & ct_st_2 != "Others::CD274+" & ct_st_1 == "Others::PDCD1-" & ct_st_2 != "Others::CD274-") {
    return(ct_st_2)
  } else if (ct_st_1 != "Others::PDCD1+" & ct_st_2 != "Others::CD274+" & ct_st_1 != "Others::PDCD1-" & ct_st_2 != "Others::CD274-") {
    return(paste0(ct_st_1, "::", ct_st_2))
  } else {
    return("Others")
  }
}, group_ct_st_1, group_ct_st_2)

# Apply changes using gsub
Group_combine <- gsub("CD4\\+ T Cells::PDCD1- ::T Cells::CD274-", "CD4+ T Cells", Group_combine)
Group_combine <- gsub("CD4\\+ T Cells::PDCD1\\+ ::T Cells::CD274\\+", "CD4+ T Cells::PDCD1+::CD274+", Group_combine)
Group_combine <- gsub("Mono/Mac::PDCD1- ::Others::CD274-", "Mono/Mac", Group_combine)
Group_combine <- gsub("Others::PDCD1- ::Capillaries & Venules::CD274-", "Others", Group_combine)
Group_combine <- gsub("Others::PDCD1- ::Dendritic Cells::CD274-", "Dendritic Cells", Group_combine)
Group_combine <- gsub("Others::PDCD1- ::Others::CD274-", "Others", Group_combine)
Group_combine <- gsub("Others::PDCD1\\+ ::Others::CD274\\+", "Others", Group_combine)
Group_combine <- gsub("T Cells::PDCD1- ::T Cells::CD274-", "T Cells", Group_combine)
Group_combine <- gsub("T Cells::PDCD1\\+ ::T Cells::CD274\\+", "T Cells::PDCD1+::CD274+", Group_combine)

Group_combine <- ifelse(data$TACIT_PCF_Xenium == "B Cells", "B Cells", Group_combine)

# Calculate the number of observations for each group
group_counts <- as.data.frame(table(Group_combine))
colnames(group_counts) <- c("group_ct_st", "count")

# Merge the counts back into the original data
data_plot_ct <- data.frame(X = data$X, Y = data$Y, group_ct_st = Group_combine)
data_plot_ct <- merge(data_plot_ct, group_counts, by = "group_ct_st")

# Create a new column with the group names and counts
data_plot_ct$group_with_count <- paste(data_plot_ct$group_ct_st, "(n =", data_plot_ct$count, ")")

round(table(data_plot_ct$group_with_count) / length(data_plot_ct$group_with_count), 4)

# Plot the data
ggplot(data_plot_ct, aes(x = X, y = Y)) +
  geom_point(aes(color = as.factor(group_with_count)), size = 1) +
  theme_classic(base_size = 15) +
  guides(color = guide_legend(override.aes = list(size = 3), title = "Group"))




