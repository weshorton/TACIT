# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(ggforce)
library(tidyr)
library(pheatmap)
# Read data from CSV files
data=read.csv("~/Figure 3/result_TACIT_Xenium1_5_merge.csv")
Signature=readr::read_csv("~/Figure 3/Signature_Xenium_TACIT_Louvain.csv")
Signature=Signature[,-2]

#Color
user_colors <-  c(
  "Arterioles"="darkgrey",
  "B Cells"="#FFA500", 
  "Capillaries"="#aa6e28",
  "CD4+ T Cells"="#FF0000",
  "CD8+ Effector T cells"="#CC0000", 
  "CD8+ Exhausted T Cells"="#FF6347",
  "Dendritic Cells"="#FFD700",
  "Ductal Cells"="#00FF00",
  "Ductal Progenitors"="#008000",
  "Fibroblasts"="#f032e6",
  "IgA Plasma Cells"="#7B68EE",
  "IgG Plasma Cells"="purple",
  "Intermediate Epithelium"="powderblue",
  "Ionocytes"="lightgreen",
  "M1 Macrophages"="yellow3",
  "M2 Macrophages"="gold3",
  "Mucous Acinar Cells"="cyan",
  "Acinar Cells"="cyan",
  "Myoepithelium"="blue",
  "NK Cells"="darkred",
  "Pericytes"="#FFC0CB",
  "Regulatory T Cells"="brown2",
  "Seromucous Acinar Cells"="royalblue2",
  "Smooth Muscle"="azure3",
  "Others"="grey90",
  "T Cell Progenitors"="brown3",
  "Venules"="orange3","VEC"="orange3")

#-----------------------------
#f----------------------------
#-----------------------------


data_proportion=read.csv("~/Figure 3/cluster_proportions_scRNA.csv")

predicted=c(data$Seurat)
table(predicted)/length(predicted)
predicted_TACIT=c(data$TACIT_update)

# Example data setup - ensure this matches your actual data structure
predicted_Seurat <- data.frame(
  CellType = names(table(predicted)/length(predicted)),
  Seurat = as.numeric((table(predicted)/length(predicted)))
)

data_proportion <- data.frame(
  CellType =data_proportion$Cluster ,
  scRNA = data_proportion$Proportion.Freq
)

predicted_TACIT <- data.frame(
  CellType = names(table(predicted_TACIT)/length(predicted_TACIT)),
  TACIT = as.numeric((table(predicted_TACIT)/length(predicted_TACIT)))
)

predicted_Seurat_org <- data.frame(
  CellType = names(table(data$Seurat_orginial)/length(data$Seurat_orginial)),
  Seurat_org = as.numeric((table(data$Seurat_orginial)/length(data$Seurat_orginial)))
)



# Merge data frames based on cell type
merged_data <- merge(predicted_Seurat, data_proportion, "CellType")

merged_data <- merge(predicted_TACIT, merged_data, "CellType",all.y=T)

merged_data <- merge(predicted_Seurat_org, merged_data, "CellType",all.y=T)

merged_data$Seurat_org[is.na(merged_data$Seurat_org)==T]=0


#merged_data=merged_data[which(merged_data$scRNA<0.01),]

cor_TACIT=cor(merged_data$scRNA,merged_data$TACIT)
cor_Seurat=cor(merged_data$scRNA,merged_data$Seurat)
cor_Seurat_org=cor(merged_data$scRNA,merged_data$Seurat_org)

ggplot() +
  # Add points and linear model for TACIT vs. scRNA
  geom_point(data = merged_data, aes(x = scRNA, y = TACIT, color = "TACIT"), size = 4) +
  geom_smooth(data = merged_data, aes(x = scRNA, y = TACIT, color = "TACIT"), method = "lm", se = FALSE) +
  # Add points and linear model for Seurat vs. scRNA
  geom_point(data = merged_data, aes(x = scRNA, y = Seurat, color = "Seurat transfer"), size = 4) +
  geom_smooth(data = merged_data, aes(x = scRNA, y = Seurat, color = "Seurat transfer"), method = "lm", se = FALSE)+
  # Add points and linear model for Seurat vs. scRNA
  geom_point(data = merged_data, aes(x = scRNA, y = Seurat_org, color = "Seurat"), size = 4) +
  geom_smooth(data = merged_data, aes(x = scRNA, y = Seurat_org, color = "Seurat"), method = "lm", se = FALSE) +  labs(
    y = "",
    x = "",
    title = "",
    color = "Methods"  # Changing legend title for color
  ) +
  scale_color_manual(values = c("Seurat transfer" = "blue", "TACIT" = "red","Seurat"="purple")) +
  theme_classic(base_size = 25) +
  ylim(0, 0.25) +
  xlim(0, 0.45)



#-----------------------------
#g----------------------------
#-----------------------------

length(unique(data$Seurat))
length(unique(data$Seurat_orginial))
length(unique(data$TACIT_update))


# Create a data frame with the number of unique cell types identified by each method
data_unique_counts <- data.frame(
  Method = c("Seurat transfer", "Seurat", "TACIT"),
  Unique_Cell_Types = c(24, 16, 24)
)

# Plot the bar plot
ggplot(data_unique_counts, aes(x = Method, y = Unique_Cell_Types, fill = Method)) +
  geom_bar(stat = "identity") +
  labs(title = "Number of Unique Cell Types Identified by Each Method",
       x = "Method",
       y = "Number of Unique Cell Types") +
  theme_classic(base_size = 20) +
  scale_fill_manual(values = c("Seurat transfer" = "blue", "Seurat" = "purple", "TACIT" = "red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Adjust text angle for better readability
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")  # Remove legend as fill and x-axis are the same






#-----------------------------
#c----------------------------
#-----------------------------

# Group the data by the predicted label and calculate the mean for each column
data_plot <- data.frame(TACIT = c(data$TACIT_update), (data[, colnames(Signature)[-1]]))
colnames(data_plot)[-1] <- colnames(Signature)[-1]

# Calculate the 90th percentile for each cell type
mean_values_TACIT <- data_plot %>%
  group_by(TACIT) %>%
  summarise_all(~quantile(., 0.9))

mean_values_TACIT <- as.data.frame(mean_values_TACIT)
rownames(mean_values_TACIT) <- mean_values_TACIT$TACIT
mean_values_TACIT <- mean_values_TACIT[mean_values_TACIT$TACIT != "Others", ]

mean_values_TACIT <- as.data.frame(mean_values_TACIT[, -1])

# Reorder Signature based on custom order
Signature <- as.data.frame(Signature)
Signature_order <- reorder_columns_based_on_signature(Signature)

custom_order <- c(
  "T Cell Progenitors", "CD4+ T Cells", "CD8+ Effector T Cells", "CD8+ Exhausted T Cells", "Regulatory T Cells",  # T cells group
  "NK Cells", "Dendritic Cells", "B Cells", "IgG Plasma Cells", "IgA Plasma Cells",  # Other immune cells
  "M1 Macrophages", "M2 Macrophages",  # Macrophages
  "Intermediate Epithelium", "Ductal Progenitors", "Ductal Cells", "Myoepithelium", "Mucous Acinar Cells", "Seromucous Acinar Cells",  # Epithelial and related cells
  "Arterioles", "Venules", "Capillaries", "Pericytes", "Smooth Muscle",  # Vascular and muscle cells
  "Ionocytes", "Fibroblasts"  # Other cell types
)

# Use factor() to order the cell types vector according to the custom order
Signature$cell_type <- factor(Signature$cell_type, levels = custom_order)
Signature <- Signature[order(Signature$cell_type), ]

# Assuming mean_values_TACIT is a data frame or matrix and Signature_order is a data frame or list
ordered_index <- match(Signature_order$cell_type, rownames(mean_values_TACIT))

# Reorder mean_values_TACIT according to the ordered index
mean_values_TACIT <- mean_values_TACIT[ordered_index, ]

# Set breaks and colors for heatmap
my.breaks <- c(seq(-2, 0, by = 0.1), seq(0.1, 2, by = 0.1))
my.colors <- c(colorRampPalette(colors = c("blue", "white"))(length(my.breaks) / 2), colorRampPalette(colors = c("white", "red"))(length(my.breaks) / 2))

# Fill NA values in Signature_order with 0
Signature_order[is.na(Signature_order)] <- 0
Signature_order <- as.data.frame(Signature_order)
rownames(Signature_order) <- Signature_order$cell_type

# Plot the Signature heatmap
pheatmap(Signature_order[, -1],
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         show_colnames = TRUE,
         fontsize_col = 15,
         fontsize_row = 15)

# Scale and cap the values for mean_values_TACIT
aa <- scale(mean_values_TACIT[, colnames(Signature_order)[-1]])
aa[aa > 3] <- 3
rownames(aa) <- Signature_order$cell_type

# Plot the mean values heatmap
pheatmap(aa,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         show_colnames = TRUE,
         fontsize_col = 15,
         fontsize_row = 15)







#Seurat transfer

# Group the data by the predicted label and calculate the mean for each column
data_plot <- data.frame(TACIT = c(data$Seurat), (data[, colnames(Signature)[-1]]))
colnames(data_plot)[-1] <- colnames(Signature)[-1]

# Calculate the 90th percentile for each cell type
mean_values_TACIT <- data_plot %>%
  group_by(TACIT) %>%
  summarise_all(~quantile(., 0.9))

mean_values_TACIT <- as.data.frame(mean_values_TACIT)
rownames(mean_values_TACIT) <- mean_values_TACIT$TACIT
mean_values_TACIT <- mean_values_TACIT[mean_values_TACIT$TACIT != "Others", ]

mean_values_TACIT <- as.data.frame(mean_values_TACIT[, -1])

# Reorder Signature based on custom order
Signature <- as.data.frame(Signature)

# Use factor() to order the cell types vector according to the custom order
Signature$cell_type <- factor(Signature$cell_type, levels = custom_order)
Signature <- Signature[order(Signature$cell_type), ]

# Assuming mean_values_TACIT is a data frame or matrix and Signature_order is a data frame or list
ordered_index <- match(Signature_order$cell_type, rownames(mean_values_TACIT))

# Reorder mean_values_TACIT according to the ordered index
mean_values_TACIT <- mean_values_TACIT[ordered_index, ]

# Set breaks and colors for heatmap
my.breaks <- c(seq(-2, 0, by = 0.1), seq(0.1, 2, by = 0.1))
my.colors <- c(colorRampPalette(colors = c("blue", "white"))(length(my.breaks) / 2), colorRampPalette(colors = c("white", "red"))(length(my.breaks) / 2))

# Fill NA values in Signature_order with 0
Signature_order[is.na(Signature_order)] <- 0
Signature_order <- as.data.frame(Signature_order)
rownames(Signature_order) <- Signature_order$cell_type

# Plot the Signature heatmap
pheatmap(Signature_order[, -1],
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         show_colnames = TRUE,
         fontsize_col = 15,
         fontsize_row = 15)

# Scale and cap the values for mean_values_TACIT
aa <- scale(mean_values_TACIT[, colnames(Signature_order)[-1]])
aa[aa > 3] <- 3
rownames(aa) <- Signature_order$cell_type

# Plot the mean values heatmap
pheatmap(aa,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         show_colnames = TRUE,
         fontsize_col = 15,
         fontsize_row = 15)







#Seurat v5

# Group the data by the predicted label and calculate the mean for each column
data_plot <- data.frame(TACIT = c(data$Seurat_orginial), (data[, colnames(Signature)[-1]]))
colnames(data_plot)[-1] <- colnames(Signature)[-1]

# Calculate the 90th percentile for each cell type
mean_values_TACIT <- data_plot %>%
  group_by(TACIT) %>%
  summarise_all(~quantile(., 0.9))

mean_values_TACIT <- as.data.frame(mean_values_TACIT)
rownames(mean_values_TACIT) <- mean_values_TACIT$TACIT
mean_values_TACIT <- mean_values_TACIT[mean_values_TACIT$TACIT != "Others", ]

mean_values_TACIT <- as.data.frame(mean_values_TACIT[, -1])

# Reorder Signature based on custom order
Signature <- as.data.frame(Signature)

# Use factor() to order the cell types vector according to the custom order
Signature$cell_type <- factor(Signature$cell_type, levels = custom_order)
Signature <- Signature[order(Signature$cell_type), ]

# Assuming mean_values_TACIT is a data frame or matrix and Signature_order is a data frame or list
ordered_index <- match(Signature_order$cell_type, rownames(mean_values_TACIT))

# Reorder mean_values_TACIT according to the ordered index
mean_values_TACIT <- mean_values_TACIT[ordered_index, ]

# Set breaks and colors for heatmap
my.breaks <- c(seq(-2, 0, by = 0.1), seq(0.1, 2, by = 0.1))
my.colors <- c(colorRampPalette(colors = c("blue", "white"))(length(my.breaks) / 2), colorRampPalette(colors = c("white", "red"))(length(my.breaks) / 2))

# Fill NA values in Signature_order with 0
Signature_order[is.na(Signature_order)] <- 0
Signature_order <- as.data.frame(Signature_order)
rownames(Signature_order) <- Signature_order$cell_type

# Plot the Signature heatmap
pheatmap(Signature_order[, -1],
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         show_colnames = TRUE,
         fontsize_col = 15,
         fontsize_row = 15)

# Scale and cap the values for mean_values_TACIT
aa <- scale(mean_values_TACIT[, colnames(Signature_order)[-1]])
aa[aa > 3] <- 3
rownames(aa) <- Signature_order$cell_type

# Plot the mean values heatmap
pheatmap(aa,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         show_colnames = TRUE,
         fontsize_col = 15,
         fontsize_row = 15)




#-----------------------------
#a----------------------------
#-----------------------------


#TACIT
ggplot(data, aes(x = UMAP1, y = UMAP2, color = TACIT_update)) +
  geom_point(size = 0.1) + 
  scale_color_manual(values = user_colors) +
  theme_classic(base_size = 15)+ guides(color = guide_legend(override.aes = list(size = 3),title = "Cell type"))+theme(legend.position = "bottom")

#Seurat transfer
ggplot(data, aes(x = UMAP1, y = UMAP2, color = Seurat)) +
  geom_point(size = 0.1) + 
  scale_color_manual(values = user_colors) +
  theme_classic(base_size = 15)+ guides(color = guide_legend(override.aes = list(size = 3),title = "Cell type"))+theme(legend.position = "bottom")
#Seurat
ggplot(data, aes(x = UMAP1, y = UMAP2, color = Seurat_orginial)) +
  geom_point(size = 0.1) + 
  scale_color_manual(values = user_colors) +
  theme_classic(base_size = 15)+ guides(color = guide_legend(override.aes = list(size = 3),title = "Cell type"))+theme(legend.position = "bottom")



#-----------------------------
#e----------------------------
#-----------------------------

#Seurat

# Subset the data for Group 3 and within specified X range
data_plot_sub <- data[which(data$Group == "Group 3"), ]
data_plot_sub <- data_plot_sub[which(data_plot_sub$X < 5000 & data_plot_sub$X > 3000), ]

# Define group based on specific cell types
group <- ifelse(data_plot_sub$Seurat_orginial %in% c("T Cell Progenitors", "CD8+ Exhausted T Cells",
                                                     "CD4+ T Cells", "CD8+ Effector T Cells"), 
                data_plot_sub$Seurat_orginial, "Others")

# Plot the data with color differentiation for the defined groups
ggplot(data_plot_sub, aes(x = X, y = Y, color = group)) +
  geom_point(size = 0.1) + 
  scale_color_manual(values = c("Others" = "grey90", "T Cell Progenitors" = "orange", 
                                "CD4+ T Cells" = "darkgreen", "CD8+ Effector T Cells" = "#f032e6", 
                                "CD8+ Exhausted T Cells" = "royalblue2")) +
  theme_classic(base_size = 15) +
  guides(color = guide_legend(override.aes = list(size = 3), title = "Cell type")) +
  theme(legend.position = "none",
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
        axis.line.y = element_blank(), axis.line.x = element_blank()) +
  xlab("") + ylab("")

# Filter data for selected groups
data_plot_sub <- data_plot_sub[which(group %in% c("T Cell Progenitors", "CD8+ Exhausted T Cells",
                                                  "CD4+ T Cells", "CD8+ Effector T Cells")), ]

# Subset the Signature dataframe for selected cell types
Signature_sub <- Signature[which(Signature$cell_type %in% c("T Cell Progenitors", "CD8+ Exhausted T Cells",
                                                            "CD4+ T Cells", "CD8+ Effector T Cells")), -1]

# Select markers with non-zero sums
selected_marker <- names(colSums(Signature_sub))[which(as.numeric(colSums(Signature_sub)) > 0)]

# Re-define group based on specific cell types
group <- ifelse(data_plot_sub$Seurat_orginial %in% c("T Cell Progenitors", "CD8+ Exhausted T Cells",
                                                     "CD4+ T Cells", "CD8+ Effector T Cells"), 
                data_plot_sub$Seurat_orginial, "Others")

# Perform UMAP on the subset data
data_umap <- umap::umap(data_plot_sub[, colnames(Signature)[-1]], scale = TRUE, method = "umap-learn")
data_umap <- data_umap$layout
colnames(data_umap) <- c("UMAP1", "UMAP2")
data_umap <- as.data.frame(data_umap)

# Plot the UMAP result with color differentiation for the defined groups
ggplot(data_umap, aes(x = UMAP1, y = UMAP2, color = c(group))) +
  geom_point(size = 2) + 
  scale_color_manual(values = c("T Cell Progenitors" = "orange", "CD4+ T Cells" = "darkgreen", 
                                "CD8+ Effector T Cells" = "#f032e6", "CD8+ Exhausted T Cells" = "royalblue2")) +
  theme_classic(base_size = 15) +
  guides(color = guide_legend(override.aes = list(size = 3), title = "Cell type")) +
  theme(legend.position = "none",
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
        axis.line.y = element_blank(), axis.line.x = element_blank()) +
  xlab("") + ylab("")









#TACIT

# Subset the data for Group 3 and within specified X range
data_plot_sub <- data[which(data$Group == "Group 3"), ]
data_plot_sub <- data_plot_sub[which(data_plot_sub$X < 5000 & data_plot_sub$X > 3000), ]

# Define group based on specific cell types
group <- ifelse(data_plot_sub$TACIT_update %in% c("T Cell Progenitors", "CD8+ Exhausted T Cells",
                                                     "CD4+ T Cells", "CD8+ Effector T Cells"), 
                data_plot_sub$TACIT_update, "Others")

# Plot the data with color differentiation for the defined groups
ggplot(data_plot_sub, aes(x = X, y = Y, color = group)) +
  geom_point(size = 0.1) + 
  scale_color_manual(values = c("Others" = "grey90", "T Cell Progenitors" = "orange", 
                                "CD4+ T Cells" = "darkgreen", "CD8+ Effector T Cells" = "#f032e6", 
                                "CD8+ Exhausted T Cells" = "royalblue2")) +
  theme_classic(base_size = 15) +
  guides(color = guide_legend(override.aes = list(size = 3), title = "Cell type")) +
  theme(legend.position = "none",
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
        axis.line.y = element_blank(), axis.line.x = element_blank()) +
  xlab("") + ylab("")

# Filter data for selected groups
data_plot_sub <- data_plot_sub[which(group %in% c("T Cell Progenitors", "CD8+ Exhausted T Cells",
                                                  "CD4+ T Cells", "CD8+ Effector T Cells")), ]

# Subset the Signature dataframe for selected cell types
Signature_sub <- Signature[which(Signature$cell_type %in% c("T Cell Progenitors", "CD8+ Exhausted T Cells",
                                                            "CD4+ T Cells", "CD8+ Effector T Cells")), -1]

# Select markers with non-zero sums
selected_marker <- names(colSums(Signature_sub))[which(as.numeric(colSums(Signature_sub)) > 0)]

# Re-define group based on specific cell types
group <- ifelse(data_plot_sub$TACIT_update %in% c("T Cell Progenitors", "CD8+ Exhausted T Cells",
                                                     "CD4+ T Cells", "CD8+ Effector T Cells"), 
                data_plot_sub$TACIT_update, "Others")

# Perform UMAP on the subset data
data_umap <- umap::umap(data_plot_sub[, colnames(Signature)[-1]], scale = TRUE, method = "umap-learn")
data_umap <- data_umap$layout
colnames(data_umap) <- c("UMAP1", "UMAP2")
data_umap <- as.data.frame(data_umap)

# Plot the UMAP result with color differentiation for the defined groups
ggplot(data_umap, aes(x = UMAP1, y = UMAP2, color = c(group))) +
  geom_point(size = 2) + 
  scale_color_manual(values = c("T Cell Progenitors" = "orange", "CD4+ T Cells" = "darkgreen", 
                                "CD8+ Effector T Cells" = "#f032e6", "CD8+ Exhausted T Cells" = "royalblue2")) +
  theme_classic(base_size = 15) +
  guides(color = guide_legend(override.aes = list(size = 3), title = "Cell type")) +
  theme(legend.position = "none",
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
        axis.line.y = element_blank(), axis.line.x = element_blank()) +
  xlab("") + ylab("")









#Seurat transfer

# Subset the data for Group 3 and within specified X range
data_plot_sub <- data[which(data$Group == "Group 3"), ]
data_plot_sub <- data_plot_sub[which(data_plot_sub$X < 5000 & data_plot_sub$X > 3000), ]

# Define group based on specific cell types
group <- ifelse(data_plot_sub$Seurat %in% c("T Cell Progenitors", "CD8+ Exhausted T Cells",
                                                  "CD4+ T Cells", "CD8+ Effector T Cells"), 
                data_plot_sub$Seurat, "Others")

# Plot the data with color differentiation for the defined groups
ggplot(data_plot_sub, aes(x = X, y = Y, color = group)) +
  geom_point(size = 0.1) + 
  scale_color_manual(values = c("Others" = "grey90", "T Cell Progenitors" = "orange", 
                                "CD4+ T Cells" = "darkgreen", "CD8+ Effector T Cells" = "#f032e6", 
                                "CD8+ Exhausted T Cells" = "royalblue2")) +
  theme_classic(base_size = 15) +
  guides(color = guide_legend(override.aes = list(size = 3), title = "Cell type")) +
  theme(legend.position = "none",
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
        axis.line.y = element_blank(), axis.line.x = element_blank()) +
  xlab("") + ylab("")

# Filter data for selected groups
data_plot_sub <- data_plot_sub[which(group %in% c("T Cell Progenitors", "CD8+ Exhausted T Cells",
                                                  "CD4+ T Cells", "CD8+ Effector T Cells")), ]

# Subset the Signature dataframe for selected cell types
Signature_sub <- Signature[which(Signature$cell_type %in% c("T Cell Progenitors", "CD8+ Exhausted T Cells",
                                                            "CD4+ T Cells", "CD8+ Effector T Cells")), -1]

# Select markers with non-zero sums
selected_marker <- names(colSums(Signature_sub))[which(as.numeric(colSums(Signature_sub)) > 0)]

# Re-define group based on specific cell types
group <- ifelse(data_plot_sub$Seurat %in% c("T Cell Progenitors", "CD8+ Exhausted T Cells",
                                                  "CD4+ T Cells", "CD8+ Effector T Cells"), 
                data_plot_sub$Seurat, "Others")

# Perform UMAP on the subset data
data_umap <- umap::umap(data_plot_sub[, colnames(Signature)[-1]], scale = TRUE, method = "umap-learn")
data_umap <- data_umap$layout
colnames(data_umap) <- c("UMAP1", "UMAP2")
data_umap <- as.data.frame(data_umap)

# Plot the UMAP result with color differentiation for the defined groups
ggplot(data_umap, aes(x = UMAP1, y = UMAP2, color = c(group))) +
  geom_point(size = 2) + 
  scale_color_manual(values = c("T Cell Progenitors" = "orange", "CD4+ T Cells" = "darkgreen", 
                                "CD8+ Effector T Cells" = "#f032e6", "CD8+ Exhausted T Cells" = "royalblue2")) +
  theme_classic(base_size = 15) +
  guides(color = guide_legend(override.aes = list(size = 3), title = "Cell type")) +
  theme(legend.position = "none",
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
        axis.line.y = element_blank(), axis.line.x = element_blank()) +
  xlab("") + ylab("")


#-----------------------------
#e----------------------------
#-----------------------------

#TACIT

# Create a data frame combining TACIT outcomes and signature data
data_plot <- data.frame(TACIT = data$TACIT_update, (data[, colnames(Signature)[-1]]))
colnames(data_plot)[-1] <- colnames(Signature)[-1]

# Prepare the Signature data
colnames(Signature)[1] <- "cell_type"
Signature[is.na(Signature)] <- 0

# Calculate the median values for each cell type in TACIT
mean_values_TACIT <- data_plot %>%
  group_by(TACIT) %>%
  summarise_all(~quantile(., 0.9)) %>%
  as.data.frame()

rownames(mean_values_TACIT) <- mean_values_TACIT$TACIT
mean_values_TACIT <- mean_values_TACIT[which(mean_values_TACIT$TACIT != "Others"), ]
mean_values_TACIT <- as.data.frame(mean_values_TACIT[, -1])

# Reorder columns based on signature
Signature_order <- reorder_columns_based_on_signature(Signature = Signature)
ordered_index <- match(Signature_order$cell_type, rownames(mean_values_TACIT))
mean_values_TACIT <- mean_values_TACIT[ordered_index, ]

# Define color breaks and palette for heatmap
my.breaks <- c(seq(-2, 0, by = 0.1), seq(0.1, 2, by = 0.1))
my.colors <- c(colorRampPalette(colors = c("blue", "white"))(length(my.breaks) / 2), colorRampPalette(colors = c("white", "red"))(length(my.breaks) / 2))

# Prepare data for analysis
aa <- mean_values_TACIT[, colnames(Signature_order)[-1]]
rownames(aa) <- Signature_order$cell_type
aa <- aa - min(aa, na.rm = TRUE)
data_anb <- aa
ct <- Signature$cell_type
aa <- as.matrix(aa)

# Initialize matrix for storing p-values
mat <- matrix(NA, ncol = length(ct), nrow = ncol(data_anb))

# Loop over each marker and each cell type to perform t-tests
for (i in 1:ncol(data_anb)) {
  for (j in 1:length(ct)) {
    ct_choose <- ct[j]
    marker_choose <- colnames(data_anb)[i]
    
    select_ct_0 <- Signature$cell_type[Signature[, marker_choose] != 1]
    select_ct_0 <- select_ct_0[select_ct_0 != ct_choose]
    
    test_result <- tryCatch({
      t.test(mu = as.numeric(data_anb[rownames(data_anb) == ct_choose, marker_choose]),
             x = as.numeric(data_anb[rownames(data_anb) %in% select_ct_0, marker_choose]),
             alternative = "less")
    }, error = function(e) NA)
    
    mat[i, j] <- if (length(test_result) > 1) test_result$p.value else 1
  }
}

rownames(mat) <- colnames(data_anb)
colnames(mat) <- ct
aa_df <- as.data.frame(mat)
aa_df$Row_Name <- rownames(mat)

# Transform p-values to long format
p_val_long_format <- aa_df %>%
  pivot_longer(cols = -Row_Name, names_to = "Column_Name", values_to = "Value")

p_val_long_format$padjust <- p.adjust(p_val_long_format$Value, method = "BH")
p_val_long_format$ID <- paste0(p_val_long_format$Row_Name, "::", p_val_long_format$Column_Name)

# Prepare signature data in long format
S <- as.data.frame(Signature[, -1])
rownames(S) <- Signature$cell_type
S$Row_Name <- rownames(S)

long_S <- S %>%
  pivot_longer(cols = -Row_Name, names_to = "Column_Name", values_to = "Value")

long_S$ID <- paste0(long_S$Column_Name, "::", long_S$Row_Name)
colnames(long_S)[3] <- "Signature"

# Calculate log2 fold change for each marker and cell type
data_anb <- data[, colnames(Signature)[-1]]
data_anb <- as.data.frame(data_anb)
mat <- matrix(NA, ncol = length(ct), nrow = ncol(data_anb))

for (i in 1:ncol(data_anb)) {
  for (j in 1:length(ct)) {
    ct_choose <- ct[j]
    marker_choose <- colnames(data_anb)[i]
    
    select_ct_0 <- Signature$cell_type[which(Signature[, marker_choose] != 1)]
    select_ct_0 <- select_ct_0[which(select_ct_0 != ct_choose)]
    
    testt <- log2((mean(as.numeric(data_anb[which(data$TACIT_update == ct_choose), marker_choose]))+0.000001) /
                    (mean(as.numeric(data_anb[which(data$TACIT_update %in% select_ct_0), marker_choose]))+0.000001))
    
    mat[i, j] <- testt
  }
}

rownames(mat) <- colnames(data_anb)
colnames(mat) <- ct
mat[is.na(mat)] <- 0

aa_df <- as.data.frame(mat)
aa_df$Row_Name <- rownames(mat)

# Transform log2 fold change values to long format
logFC_long_format <- aa_df %>%
  pivot_longer(cols = -Row_Name, names_to = "Column_Name", values_to = "Value")

logFC_long_format$ID <- paste0(logFC_long_format$Row_Name, "::", logFC_long_format$Column_Name)
colnames(logFC_long_format)[3] <- "logFC"

# Merge logFC, signature, and adjusted p-values into a final dataset
data_final <- merge(logFC_long_format, long_S, by = "ID")
data_final2 <- merge(data_final, p_val_long_format, by = "ID")
data_final2$padjust[is.na(data_final2$padjust)==T]=1
# Classify data points based on significance and signature status
Group <- vector("character", nrow(data_final2))
for (i in 1:nrow(data_final2)) {
  if (data_final2$logFC[i] > 1 & data_final2$padjust[i] < 0.05 & data_final2$Signature[i] == 1) {
    Group[i] <- "Significant & Signature"
  } else if (data_final2$logFC[i] < 1 & data_final2$padjust[i] > 0.05 & data_final2$Signature[i] == 1) {
    Group[i] <- "Not Significant & Signature"
  } else if (data_final2$logFC[i] > 1 & data_final2$padjust[i] < 0.05 & data_final2$Signature[i] == 0) {
    Group[i] <- "Significant & Not Signature"
  } else if (data_final2$logFC[i] < 1 & data_final2$padjust[i] < 0.05 & data_final2$Signature[i] == 1) {
    Group[i] <- "Significant & Not Signature"
  } else {
    Group[i] <- "Others"
  }
}

# Plot the results
ggplot(data_final2, aes(x = logFC, y = -log10(padjust), color = Group)) +
  geom_point(alpha = 1) +
  scale_color_manual(values = c("Significant & Not Signature" = "black", "Others" = "gray90", "Significant & Signature" = "red", "Not Significant & Signature" = "blue")) +
  labs(x = "Average Log2 Fold Change", y = "-log10(p-value)", title = "TACIT") +
  theme_classic(base_size = 15)

data_final2$Group <- Group
data_final_TACIT <- data_final2







#Seurat V5

# Create a data frame combining TACIT outcomes and signature data
data_plot <- data.frame(TACIT = data$Seurat_orginial, (data[, colnames(Signature)[-1]]))
colnames(data_plot)[-1] <- colnames(Signature)[-1]

# Prepare the Signature data
colnames(Signature)[1] <- "cell_type"
Signature[is.na(Signature)] <- 0

# Calculate the median values for each cell type in TACIT
mean_values_TACIT <- data_plot %>%
  group_by(TACIT) %>%
  summarise_all(~quantile(., 0.9)) %>%
  as.data.frame()

rownames(mean_values_TACIT) <- mean_values_TACIT$TACIT
mean_values_TACIT <- mean_values_TACIT[which(mean_values_TACIT$TACIT != "Others"), ]
mean_values_TACIT <- as.data.frame(mean_values_TACIT[, -1])

# Reorder columns based on signature
Signature_order <- reorder_columns_based_on_signature(Signature = Signature)
ordered_index <- match(Signature_order$cell_type, rownames(mean_values_TACIT))
mean_values_TACIT <- mean_values_TACIT[ordered_index, ]

# Define color breaks and palette for heatmap
my.breaks <- c(seq(-2, 0, by = 0.1), seq(0.1, 2, by = 0.1))
my.colors <- c(colorRampPalette(colors = c("blue", "white"))(length(my.breaks) / 2), colorRampPalette(colors = c("white", "red"))(length(my.breaks) / 2))

# Prepare data for analysis
aa <- mean_values_TACIT[, colnames(Signature_order)[-1]]
rownames(aa) <- Signature_order$cell_type
aa <- aa - min(aa, na.rm = TRUE)
data_anb <- aa
ct <- Signature$cell_type
aa <- as.matrix(aa)

# Initialize matrix for storing p-values
mat <- matrix(NA, ncol = length(ct), nrow = ncol(data_anb))

# Loop over each marker and each cell type to perform t-tests
for (i in 1:ncol(data_anb)) {
  for (j in 1:length(ct)) {
    ct_choose <- ct[j]
    marker_choose <- colnames(data_anb)[i]
    
    select_ct_0 <- Signature$cell_type[Signature[, marker_choose] != 1]
    select_ct_0 <- select_ct_0[select_ct_0 != ct_choose]
    
    test_result <- tryCatch({
      t.test(mu = as.numeric(data_anb[rownames(data_anb) == ct_choose, marker_choose]),
             x = as.numeric(data_anb[rownames(data_anb) %in% select_ct_0, marker_choose]),
             alternative = "less")
    }, error = function(e) NA)
    
    mat[i, j] <- if (length(test_result) > 1) test_result$p.value else 1
  }
}

rownames(mat) <- colnames(data_anb)
colnames(mat) <- ct
aa_df <- as.data.frame(mat)
aa_df$Row_Name <- rownames(mat)

# Transform p-values to long format
p_val_long_format <- aa_df %>%
  pivot_longer(cols = -Row_Name, names_to = "Column_Name", values_to = "Value")

p_val_long_format$padjust <- p.adjust(p_val_long_format$Value, method = "BH")
p_val_long_format$ID <- paste0(p_val_long_format$Row_Name, "::", p_val_long_format$Column_Name)

# Prepare signature data in long format
S <- as.data.frame(Signature[, -1])
rownames(S) <- Signature$cell_type
S$Row_Name <- rownames(S)

long_S <- S %>%
  pivot_longer(cols = -Row_Name, names_to = "Column_Name", values_to = "Value")

long_S$ID <- paste0(long_S$Column_Name, "::", long_S$Row_Name)
colnames(long_S)[3] <- "Signature"

# Calculate log2 fold change for each marker and cell type
data_anb <- data[, colnames(Signature)[-1]]
data_anb <- as.data.frame(data_anb)
mat <- matrix(NA, ncol = length(ct), nrow = ncol(data_anb))

for (i in 1:ncol(data_anb)) {
  for (j in 1:length(ct)) {
    ct_choose <- ct[j]
    marker_choose <- colnames(data_anb)[i]
    
    select_ct_0 <- Signature$cell_type[which(Signature[, marker_choose] != 1)]
    select_ct_0 <- select_ct_0[which(select_ct_0 != ct_choose)]
    
    testt <- log2((mean(as.numeric(data_anb[which(data$Seurat_orginial == ct_choose), marker_choose]))+0.000001) /
                    (mean(as.numeric(data_anb[which(data$Seurat_orginial %in% select_ct_0), marker_choose]))+0.000001))
    
    mat[i, j] <- testt
  }
}

rownames(mat) <- colnames(data_anb)
colnames(mat) <- ct
mat[is.na(mat)] <- 0

aa_df <- as.data.frame(mat)
aa_df$Row_Name <- rownames(mat)

# Transform log2 fold change values to long format
logFC_long_format <- aa_df %>%
  pivot_longer(cols = -Row_Name, names_to = "Column_Name", values_to = "Value")

logFC_long_format$ID <- paste0(logFC_long_format$Row_Name, "::", logFC_long_format$Column_Name)
colnames(logFC_long_format)[3] <- "logFC"

# Merge logFC, signature, and adjusted p-values into a final dataset
data_final <- merge(logFC_long_format, long_S, by = "ID")
data_final2 <- merge(data_final, p_val_long_format, by = "ID")
data_final2$padjust[is.na(data_final2$padjust)==T]=1
# Classify data points based on significance and signature status
Group <- vector("character", nrow(data_final2))
for (i in 1:nrow(data_final2)) {
  if (data_final2$logFC[i] > 1 & data_final2$padjust[i] < 0.05 & data_final2$Signature[i] == 1) {
    Group[i] <- "Significant & Signature"
  } else if (data_final2$logFC[i] < 1 & data_final2$padjust[i] > 0.05 & data_final2$Signature[i] == 1) {
    Group[i] <- "Not Significant & Signature"
  } else if (data_final2$logFC[i] > 1 & data_final2$padjust[i] < 0.05 & data_final2$Signature[i] == 0) {
    Group[i] <- "Significant & Not Signature"
  } else if (data_final2$logFC[i] < 1 & data_final2$padjust[i] < 0.05 & data_final2$Signature[i] == 1) {
    Group[i] <- "Significant & Not Signature"
  } else {
    Group[i] <- "Others"
  }
}

# Plot the results
ggplot(data_final2, aes(x = logFC, y = -log10(padjust), color = Group)) +
  geom_point(alpha = 1) +
  scale_color_manual(values = c("Significant & Not Signature" = "black", "Others" = "gray90", "Significant & Signature" = "red", "Not Significant & Signature" = "blue")) +
  labs(x = "Average Log2 Fold Change", y = "-log10(p-value)", title = "TACIT") +
  theme_classic(base_size = 15)

data_final2$Group <- Group
data_final_Seurat <- data_final2





#Seurat transfer

# Create a data frame combining TACIT outcomes and signature data
data_plot <- data.frame(TACIT = data$Seurat, (data[, colnames(Signature)[-1]]))
colnames(data_plot)[-1] <- colnames(Signature)[-1]

# Prepare the Signature data
colnames(Signature)[1] <- "cell_type"
Signature[is.na(Signature)] <- 0

# Calculate the median values for each cell type in TACIT
mean_values_TACIT <- data_plot %>%
  group_by(TACIT) %>%
  summarise_all(~quantile(., 0.9)) %>%
  as.data.frame()

rownames(mean_values_TACIT) <- mean_values_TACIT$TACIT
mean_values_TACIT <- mean_values_TACIT[which(mean_values_TACIT$TACIT != "Others"), ]
mean_values_TACIT <- as.data.frame(mean_values_TACIT[, -1])

# Reorder columns based on signature
Signature_order <- reorder_columns_based_on_signature(Signature = Signature)
ordered_index <- match(Signature_order$cell_type, rownames(mean_values_TACIT))
mean_values_TACIT <- mean_values_TACIT[ordered_index, ]

# Define color breaks and palette for heatmap
my.breaks <- c(seq(-2, 0, by = 0.1), seq(0.1, 2, by = 0.1))
my.colors <- c(colorRampPalette(colors = c("blue", "white"))(length(my.breaks) / 2), colorRampPalette(colors = c("white", "red"))(length(my.breaks) / 2))

# Prepare data for analysis
aa <- mean_values_TACIT[, colnames(Signature_order)[-1]]
rownames(aa) <- Signature_order$cell_type
aa <- aa - min(aa, na.rm = TRUE)
data_anb <- aa
ct <- Signature$cell_type
aa <- as.matrix(aa)

# Initialize matrix for storing p-values
mat <- matrix(NA, ncol = length(ct), nrow = ncol(data_anb))

# Loop over each marker and each cell type to perform t-tests
for (i in 1:ncol(data_anb)) {
  for (j in 1:length(ct)) {
    ct_choose <- ct[j]
    marker_choose <- colnames(data_anb)[i]
    
    select_ct_0 <- Signature$cell_type[Signature[, marker_choose] != 1]
    select_ct_0 <- select_ct_0[select_ct_0 != ct_choose]
    
    test_result <- tryCatch({
      t.test(mu = as.numeric(data_anb[rownames(data_anb) == ct_choose, marker_choose]),
             x = as.numeric(data_anb[rownames(data_anb) %in% select_ct_0, marker_choose]),
             alternative = "less")
    }, error = function(e) NA)
    
    mat[i, j] <- if (length(test_result) > 1) test_result$p.value else 1
  }
}

rownames(mat) <- colnames(data_anb)
colnames(mat) <- ct
aa_df <- as.data.frame(mat)
aa_df$Row_Name <- rownames(mat)

# Transform p-values to long format
p_val_long_format <- aa_df %>%
  pivot_longer(cols = -Row_Name, names_to = "Column_Name", values_to = "Value")

p_val_long_format$padjust <- p.adjust(p_val_long_format$Value, method = "BH")
p_val_long_format$ID <- paste0(p_val_long_format$Row_Name, "::", p_val_long_format$Column_Name)

# Prepare signature data in long format
S <- as.data.frame(Signature[, -1])
rownames(S) <- Signature$cell_type
S$Row_Name <- rownames(S)

long_S <- S %>%
  pivot_longer(cols = -Row_Name, names_to = "Column_Name", values_to = "Value")

long_S$ID <- paste0(long_S$Column_Name, "::", long_S$Row_Name)
colnames(long_S)[3] <- "Signature"

# Calculate log2 fold change for each marker and cell type
data_anb <- data[, colnames(Signature)[-1]]
data_anb <- as.data.frame(data_anb)
mat <- matrix(NA, ncol = length(ct), nrow = ncol(data_anb))

for (i in 1:ncol(data_anb)) {
  for (j in 1:length(ct)) {
    ct_choose <- ct[j]
    marker_choose <- colnames(data_anb)[i]
    
    select_ct_0 <- Signature$cell_type[which(Signature[, marker_choose] != 1)]
    select_ct_0 <- select_ct_0[which(select_ct_0 != ct_choose)]
    
    testt <- log2((mean(as.numeric(data_anb[which(data$Seurat == ct_choose), marker_choose]))+0.000001) /
                    (mean(as.numeric(data_anb[which(data$Seurat %in% select_ct_0), marker_choose]))+0.000001))
    
    mat[i, j] <- testt
  }
}

rownames(mat) <- colnames(data_anb)
colnames(mat) <- ct
mat[is.na(mat)] <- 0

aa_df <- as.data.frame(mat)
aa_df$Row_Name <- rownames(mat)

# Transform log2 fold change values to long format
logFC_long_format <- aa_df %>%
  pivot_longer(cols = -Row_Name, names_to = "Column_Name", values_to = "Value")

logFC_long_format$ID <- paste0(logFC_long_format$Row_Name, "::", logFC_long_format$Column_Name)
colnames(logFC_long_format)[3] <- "logFC"

# Merge logFC, signature, and adjusted p-values into a final dataset
data_final <- merge(logFC_long_format, long_S, by = "ID")
data_final2 <- merge(data_final, p_val_long_format, by = "ID")
data_final2$padjust[is.na(data_final2$padjust)==T]=1
# Classify data points based on significance and signature status
Group <- vector("character", nrow(data_final2))
for (i in 1:nrow(data_final2)) {
  if (data_final2$logFC[i] > 1 & data_final2$padjust[i] < 0.05 & data_final2$Signature[i] == 1) {
    Group[i] <- "Significant & Signature"
  } else if (data_final2$logFC[i] < 1 & data_final2$padjust[i] > 0.05 & data_final2$Signature[i] == 1) {
    Group[i] <- "Not Significant & Signature"
  } else if (data_final2$logFC[i] > 1 & data_final2$padjust[i] < 0.05 & data_final2$Signature[i] == 0) {
    Group[i] <- "Significant & Not Signature"
  } else if (data_final2$logFC[i] < 1 & data_final2$padjust[i] < 0.05 & data_final2$Signature[i] == 1) {
    Group[i] <- "Significant & Not Signature"
  } else {
    Group[i] <- "Others"
  }
}

# Plot the results
ggplot(data_final2, aes(x = logFC, y = -log10(padjust), color = Group)) +
  geom_point(alpha = 1) +
  scale_color_manual(values = c("Significant & Not Signature" = "black", "Others" = "gray90", "Significant & Signature" = "red", "Not Significant & Signature" = "blue")) +
  labs(x = "Average Log2 Fold Change", y = "-log10(p-value)", title = "TACIT") +
  theme_classic(base_size = 15)

data_final2$Group <- Group
data_final_Seurat_transfer <- data_final2









data_final_Seurat$Methods="Seurat"
data_final_Seurat_transfer$Methods="Seurat transfer"
data_final_TACIT$Methods="TACIT"


data_final_plot=rbind(data_final_Seurat,data_final_Seurat_transfer,data_final_TACIT)

data_final_plot$Group_plot=ifelse(data_final_plot$Group=="Significant & Signature",data_final_plot$Methods,"Others")

data_final_plot_sub=data_final_plot[which(data_final_plot$Signature==1),]

# Filter the data to include only the specified methods
filtered_data <- data_final_plot_sub %>%
  filter(Methods %in% c("TACIT", "Seurat transfer", "Seurat", "Louvain","Reference"))

filtered_data$Methods=factor(filtered_data$Methods,level=c("Reference",
                                                           "TACIT","Seurat transfer","Seurat","Louvain"))

#logFC
ggplot(filtered_data, aes(x = Methods, y = logFC, fill = Methods)) +
  geom_boxplot(alpha = 1, outlier.size = 2, position = position_dodge(0.8)) +
  scale_fill_manual(values = c("TACIT" = "red", "CELESTA" = "purple", "Seurat transfer" = "green", "Seurat" = "blue", "Reference" = "yellow")) +
  theme_classic(base_size = 20) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  ylab("") +
  xlab("") +
  geom_signif(comparisons = list( c("TACIT", "Seurat transfer"), c("TACIT", "Seurat")),
              map_signif_level=TRUE,
              textsize = 3, vjust = -0.2)


#-log10(padjust)
ggplot(filtered_data, aes(x = Methods, y = -log10(padjust), fill = Methods)) +
  geom_boxplot(alpha = 1, outlier.size = 2, position = position_dodge(0.8)) +
  scale_fill_manual(values = c("TACIT" = "red", "CELESTA" = "purple", "Seurat transfer" = "green", "Seurat" = "blue", "Reference" = "yellow")) +
  theme_classic(base_size = 20) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  ylab("") +
  xlab("") +
  geom_signif(comparisons = list( c("TACIT", "Seurat transfer"), c("TACIT", "Seurat")),
              map_signif_level=TRUE,
              textsize = 3, vjust = -0.2)







