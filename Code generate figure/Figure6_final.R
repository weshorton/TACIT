# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(ggforce)
library(tidyr)
# Read data from CSV files
data1=read.csv("~/Documents/TACIT figure/dataset/Figure 6/PCF_Xenium_reg001.csv")
data2=read.csv("~/Documents/TACIT figure/dataset/Figure 6/PCF_Xenium_reg002.csv")
data3=read.csv("~/Documents/TACIT figure/dataset/Figure 6/PCF_Xenium_reg003.csv")
data4=read.csv("~/Documents/TACIT figure/dataset/Figure 6/PCF_Xenium_reg004.csv")
data5=read.csv("~/Documents/TACIT figure/dataset/Figure 6/PCF_Xenium_reg005.csv")
data6=read.csv("~/Documents/TACIT figure/dataset/Figure 6/PCF_Xenium_reg006.csv")





#-------------------
#-d-----------------
#-------------------

# Function to calculate correlations and prepare data
prepare_correlation_data <- function(data) {
  correlation_data <- data.frame(
    Gene_Pair = c("CD20-MS4A1", "CD4.x-CD4.y", "FOXP3.x-FOXP3.y", "CD3e-CD3E", 
                  "PD.1-PDCD1", "KRT14-Keratin_14", "CD31-PECAM1", "SMA-ACTA2", "CAV1-Caveolin"),
    Correlation = c(
      cor(data$CD20, data$MS4A1, method = "spearman"),
      cor(data$CD4.x, data$CD4.y, method = "spearman"),
      cor(data$FOXP3.x, data$FOXP3.y, method = "spearman"),
      cor(data$CD3e, data$CD3E, method = "spearman"),
      cor(data$PD.1, data$PDCD1, method = "spearman"),
      cor(data$KRT14, data$Keratin_14, method = "spearman"),
      cor(data$CD31, data$PECAM1, method = "spearman"),
      cor(data$SMA, data$ACTA2, method = "spearman"),
      cor(data$CAV1, data$Caveolin, method = "spearman")
    ),
    Category = c("Immune Marker", "Immune Marker", "Immune Marker", "Immune Marker", 
                 "Immune Marker", "Structural Marker", "Structural Marker", "Structural Marker", "Structural Marker")
  )
  return(correlation_data)
}



# List of datasets
data_list <- list(data1 = data1, data2 = data2, data3 = data3, data4 = data4, data5 = data5, data6 = data6)

# Prepare a combined data frame for all datasets
all_correlation_data <- do.call(rbind, lapply(data_list, prepare_correlation_data))

# Add dataset identifier
all_correlation_data$Dataset <- rep(names(data_list), each = 9)

# Perform the Wilcoxon rank-sum test
stat_test <- compare_means(Correlation ~ Category, data = all_correlation_data, method = "wilcox.test")

# Plot the combined boxplots with p-value annotation
ggplot(all_correlation_data, aes(x = Category, y = Correlation, fill = Category)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, size = 2) +
  labs(title = "",
       x = "Category",
       y = "Correlation") +
  theme_classic(base_size = 20) +
  scale_fill_manual(values = c("Immune Marker" = "blue", "Structural Marker" = "green")) +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 5)  +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())





#-------------------
#-e-----------------
#-------------------
# Define custom colors for different cell types
user_colors <- c(
  "Arterioles"="darkgrey",
  "B Cells"="#FFA500", 
  "Capillaries"="#aa6e28",
  "CD4+ T Cells"="#FF0000",
  "T Helper"="#FF0000",
  "CD8+ Effector T Cells"="#CC0000", 
  "CD8+ Exhausted T Cells"="#FF6347",
  "CD8+ T Cells exhausted"="#FF6347",
  "CD8+ T Cells"="#FF6347",
  "Dendritic Cells"="#FFD700",
  "Ductal Cells"="#00FF00",
  "Ductal Progenitors"="#008000",
  "Ductal Proliferating"="#008000",
  "Fibroblasts"="#f032e6",
  "Myofibroblast"="#f032e6",
  "IgA Plasma Cells"="#7B68EE",
  "IgG Plasma Cells"="purple",
  "Intermediate Epithelium"="powderblue",
  "Ionocytes"="lightgreen",
  "M1 Macrophages"="yellow3",
  "M2 Macrophages"="gold3",
  "Macrophage"="gold3",
  "Mucous Acinar Cells"="cyan",
  "Acinar Cells"="cyan",
  "Myoepithelium"="blue",
  "Myoepitelial Cells"="blue",
  "NK Cells"="darkred",
  "Pericytes"="#FFC0CB",
  "Memory T cell"="brown3",
  "Regulatory T Cells"="brown2",
  "Treg"="brown2",
  "Seromucous Acinar Cells"="royalblue2",
  "Smooth Muscle"="azure3",
  "T Cell Progenitors"="brown3",
  "Venules"="orange3",
  "VEC"="orange3",
  "VEC Progen"="orange1",
  "LECs"="orange",
  "Others"="grey90",
  "Neutrophil"="gold",
  "Lymphoid"="darkcyan"
)

# Filter the dataset for a specific region
data_sub <- data4[which(data4$Y < 2700 & data4$Y > 2000 & data4$X < 3700 & data4$X > 3300),]

# Plot Voronoi diagram for TACIT_PCF
ggplot(data_sub, aes(X, Y, group = -1L)) + 
  geom_voronoi_tile(aes(fill = TACIT_PCF), colour = 'black', max.radius = 40) +
  scale_fill_manual(values = user_colors) +
  theme_classic(base_size = 15) +
  labs(x = "X", y = "Y", title = "") + 
  guides(fill = guide_legend(override.aes = list(size = 3), title = "Cell type")) +
  theme(
    axis.title = element_blank(),        # Remove axis titles
    axis.text = element_blank(),         # Remove axis text
    axis.ticks = element_blank(),        # Remove axis ticks
    axis.line = element_blank(),         # Remove axis lines
    panel.background = element_blank(),  # Remove panel background
    panel.grid = element_blank()         # Remove grid lines
  ) + 
  ggtitle("PCF")

# Convert Xenium to match PCF cell types
data_sub$TACIT_Xenium_UPDATE <- ifelse(data_sub$TACIT_Xenium %in% c("Arterioles", "Venules"), "VEC", data_sub$TACIT_Xenium)
data_sub$TACIT_Xenium_UPDATE <- ifelse(data_sub$TACIT_Xenium_UPDATE %in% c("Ductal Cells", "Ductal Progenitors"), "Ductal Cells", data_sub$TACIT_Xenium_UPDATE)
data_sub$TACIT_Xenium_UPDATE <- ifelse(data_sub$TACIT_Xenium_UPDATE %in% c("Mucous Acinar Cells", "Seromucous Acinar Cells"), "Acinar Cells", data_sub$TACIT_Xenium_UPDATE)
data_sub$TACIT_Xenium_UPDATE <- ifelse(data_sub$TACIT_Xenium_UPDATE %in% c("CD8+ Exhausted T Cells"), "CD8+ T Cells", data_sub$TACIT_Xenium_UPDATE)
data_sub$TACIT_Xenium_UPDATE <- ifelse(data_sub$TACIT_Xenium_UPDATE %in% c("Myoepithelium", "Intermediate Epithelium"), "Myoepithelial", data_sub$TACIT_Xenium_UPDATE)
data_sub$TACIT_Xenium_UPDATE <- ifelse(data_sub$TACIT_Xenium_UPDATE %in% unique(data_sub$TACIT_PCF), data_sub$TACIT_Xenium_UPDATE, "Others")

# Plot Voronoi diagram for TACIT_Xenium_UPDATE
ggplot(data_sub, aes(X, Y, group = -1L)) + 
  geom_voronoi_tile(aes(fill = TACIT_Xenium_UPDATE), colour = 'black', max.radius = 40) +
  scale_fill_manual(values = user_colors) +
  theme_classic(base_size = 15) +
  labs(x = "X", y = "Y", title = "") + 
  guides(fill = guide_legend(override.aes = list(size = 3), title = "Cell type")) +
  theme(
    axis.title = element_blank(),        # Remove axis titles
    axis.text = element_blank(),         # Remove axis text
    axis.ticks = element_blank(),        # Remove axis ticks
    axis.line = element_blank(),         # Remove axis lines
    panel.background = element_blank(),  # Remove panel background
    panel.grid = element_blank()         # Remove grid lines
  ) + 
  ggtitle("Xenium")

# Identify agreement and disagreement between TACIT_PCF and TACIT_Xenium_UPDATE
data_sub$Group <- ifelse(data_sub$TACIT_PCF != data_sub$TACIT_Xenium_UPDATE, "Disagreement", "Agreement")

# Plot Voronoi diagram for agreement/disagreement
ggplot(data_sub, aes(X, Y, group = -1L)) + 
  geom_voronoi_tile(aes(fill = Group), colour = 'black', max.radius = 40) +
  scale_fill_manual(values = c("Disagreement" = "red", "Agreement" = "grey90")) +
  theme_classic(base_size = 15) +
  labs(x = "X", y = "Y", title = "") + 
  guides(fill = guide_legend(override.aes = list(size = 3), title = "Cell type")) +
  theme(
    axis.title = element_blank(),        # Remove axis titles
    axis.text = element_blank(),         # Remove axis text
    axis.ticks = element_blank(),        # Remove axis ticks
    axis.line = element_blank(),         # Remove axis lines
    panel.background = element_blank(),  # Remove panel background
    panel.grid = element_blank(),        # Remove grid lines
    legend.position = "none"             # Remove legend
  )





#-------------------
#-f-----------------
#-------------------


# Plot Voronoi diagram for TACIT_PCF
ggplot(data_sub, aes(X, Y, group = -1L)) + 
  geom_voronoi_tile(aes(fill = data_sub$TACIT_PCF_v2), colour = 'black', max.radius = 40) +
  scale_fill_manual(values = user_colors) +
  theme_classic(base_size = 15) +
  labs(x = "X", y = "Y", title = "") + 
  guides(fill = guide_legend(override.aes = list(size = 3), title = "Cell type")) +
  theme(
    axis.title = element_blank(),        # Remove axis titles
    axis.text = element_blank(),         # Remove axis text
    axis.ticks = element_blank(),        # Remove axis ticks
    axis.line = element_blank(),         # Remove axis lines
    panel.background = element_blank(),  # Remove panel background
    panel.grid = element_blank()         # Remove grid lines
  ) + 
  ggtitle("PCF")


# Plot Voronoi diagram for TACIT_Xenium_UPDATE
ggplot(data_sub, aes(X, Y, group = -1L)) + 
  geom_voronoi_tile(aes(fill = TACIT_Xenium_v2), colour = 'black', max.radius = 40) +
  scale_fill_manual(values = user_colors) +
  theme_classic(base_size = 15) +
  labs(x = "X", y = "Y", title = "") + 
  guides(fill = guide_legend(override.aes = list(size = 3), title = "Cell type")) +
  theme(
    axis.title = element_blank(),        # Remove axis titles
    axis.text = element_blank(),         # Remove axis text
    axis.ticks = element_blank(),        # Remove axis ticks
    axis.line = element_blank(),         # Remove axis lines
    panel.background = element_blank(),  # Remove panel background
    panel.grid = element_blank()         # Remove grid lines
  ) + 
  ggtitle("Xenium")

# Identify agreement and disagreement between TACIT_PCF and TACIT_Xenium_UPDATE
data_sub$Group <- ifelse(data_sub$TACIT_PCF_v2 != data_sub$TACIT_Xenium_v2, "Disagreement", "Agreement")

# Plot Voronoi diagram for agreement/disagreement
ggplot(data_sub, aes(X, Y, group = -1L)) + 
  geom_voronoi_tile(aes(fill = Group), colour = 'black', max.radius = 40) +
  scale_fill_manual(values = c("Disagreement" = "red", "Agreement" = "grey90")) +
  theme_classic(base_size = 15) +
  labs(x = "X", y = "Y", title = "") + 
  guides(fill = guide_legend(override.aes = list(size = 3), title = "Cell type")) +
  theme(
    axis.title = element_blank(),        # Remove axis titles
    axis.text = element_blank(),         # Remove axis text
    axis.ticks = element_blank(),        # Remove axis ticks
    axis.line = element_blank(),         # Remove axis lines
    panel.background = element_blank(),  # Remove panel background
    panel.grid = element_blank(),        # Remove grid lines
    legend.position = "none"             # Remove legend
  )



#-------------------
#-g-----------------
#-------------------



data_sub1=data4[which(data4$Y<2700&data4$Y>2000&data4$X<3700&data4$X>3300),]
data_sub2=data4[which(data4$Y<974&data4$Y>725&data4$X<1680&data4$X>1450),]
data_sub3=data4[which(data4$Y<3565&data4$Y>3180&data4$X<1150&data4$X>670),]
data_sub4=data4[which(data4$Y<4400&data4$Y>4075&data4$X<1430&data4$X>1160),]



data_sub_final=rbind(data.frame(data_sub1,TLS="TLS1"),
                     data.frame(data_sub2,TLS="TLS2"),
                     data.frame(data_sub3,TLS="TLS3"),
                     data.frame(data_sub4,TLS="TLS4"))



data_sub_final_plot=data.frame(Xenium=data_sub_final$TACIT_Xenium_v2,PCF=data_sub_final$TACIT_PCF_v2,TLS=data_sub_final$TLS)
data_sub_final_plot=data_sub_final_plot[which(data_sub_final_plot$Xenium!="Others"),]
data_sub_final_plot=data_sub_final_plot[which(data_sub_final_plot$PCF!="Others"),]


# Calculate proportions for each Xeniumwithin each TLS group
proportions_xenium <- data_sub_final_plot %>%
  group_by(TLS, Xenium) %>%
  summarise(count_xenium = n()) %>%
  ungroup() %>%
  group_by(TLS) %>%
  mutate(proportion_xenium = count_xenium / sum(count_xenium))

proportions_xenium$Group=paste0(proportions_xenium$TLS,"::",proportions_xenium$Xenium,sep="")

# Calculate proportions for each PCF within each TLS group
proportions_pcf <- data_sub_final_plot %>%
  group_by(TLS, PCF) %>%
  summarise(count_pcf = n()) %>%
  ungroup() %>%
  group_by(TLS) %>%
  mutate(proportion_pcf = count_pcf / sum(count_pcf))

proportions_pcf$Group=paste0(proportions_pcf$TLS,"::",proportions_pcf$PCF,sep="")


# Merge proportions by TLS
merged_proportions <- merge(proportions_xenium, proportions_pcf, by = "Group")


merged_proportions_sub=merged_proportions[which(merged_proportions$Xenium==unique(merged_proportions$Xenium)[3]),]
ggplot(merged_proportions_sub, aes(x = proportion_xenium, y = proportion_pcf, color = Xenium)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") + # Adjust label position
  labs(title = "",
       x = "Proportion Xenium",
       y = "Proportion PCF") +
  theme_classic(base_size = 20) +
  theme(legend.position = "bottom")+ylim(0,max(merged_proportions_sub$proportion_xenium,merged_proportions_sub$proportion_pcf))+xlim(0,max(merged_proportions_sub$proportion_xenium,merged_proportions_sub$proportion_pcf))







#-------------------
#-h-----------------
#-------------------



data_sub1=data4[which(data4$Y<2700&data4$Y>2000&data4$X<3700&data4$X>3300),]
data_sub2=data4[which(data4$Y<974&data4$Y>725&data4$X<1680&data4$X>1450),]
data_sub3=data4[which(data4$Y<3565&data4$Y>3180&data4$X<1150&data4$X>670),]
data_sub4=data4[which(data4$Y<4400&data4$Y>4075&data4$X<1430&data4$X>1160),]



data_sub_final=rbind(data.frame(data_sub1,TLS="TLS1"),
                     data.frame(data_sub2,TLS="TLS2"),
                     data.frame(data_sub3,TLS="TLS3"),
                     data.frame(data_sub4,TLS="TLS4"))

data_sub_final$TACIT_Xenium_UPDATE <- ifelse(data_sub_final$TACIT_Xenium %in% c("Arterioles", "Venules"), "VEC", data_sub_final$TACIT_Xenium)
data_sub_final$TACIT_Xenium_UPDATE <- ifelse(data_sub_final$TACIT_Xenium_UPDATE %in% c("Ductal Cells", "Ductal Progenitors"), "Ductal Cells", data_sub_final$TACIT_Xenium_UPDATE)
data_sub_final$TACIT_Xenium_UPDATE <- ifelse(data_sub_final$TACIT_Xenium_UPDATE %in% c("Mucous Acinar Cells", "Seromucous Acinar Cells"), "Acinar Cells", data_sub_final$TACIT_Xenium_UPDATE)
data_sub_final$TACIT_Xenium_UPDATE <- ifelse(data_sub_final$TACIT_Xenium_UPDATE %in% c("CD8+ Exhausted T Cells"), "CD8+ T Cells", data_sub_final$TACIT_Xenium_UPDATE)
data_sub_final$TACIT_Xenium_UPDATE <- ifelse(data_sub_final$TACIT_Xenium_UPDATE %in% c("Myoepithelium", "Intermediate Epithelium"), "Myoepithelial", data_sub_final$TACIT_Xenium_UPDATE)
data_sub_final$TACIT_Xenium_UPDATE <- ifelse(data_sub_final$TACIT_Xenium_UPDATE %in% unique(data_sub_final$TACIT_PCF), data_sub_final$TACIT_Xenium_UPDATE, "Others")





data_sub_final_plot=data.frame(Xenium=data_sub_final$TACIT_Xenium_UPDATE,PCF=data_sub_final$TACIT_PCF,TLS=data_sub_final$TLS)
data_sub_final_plot=data_sub_final_plot[which(data_sub_final_plot$Xenium%in%c("VEC")|
                                                data_sub_final_plot$PCF%in%c("VEC")),]
# data_sub_final_plot=data_sub_final_plot[which(data_sub_final_plot$Xenium!="Others"),]
# data_sub_final_plot=data_sub_final_plot[which(data_sub_final_plot$PCF!="Others"),]

data_sub_final_plot$score=ifelse(data_sub_final_plot$Xenium==data_sub_final_plot$PCF,"Agreement","Disagreement")

agreement_proportions_all <- data_sub_final_plot %>%
  group_by(TLS) %>%
  summarise(
    Total = n(),  # Total counts per TLS
    Agreement_Count = sum(score == "Agreement"),  # Count of "Agreement" per TLS
    Agreement_Prop = Agreement_Count / Total  # Proportion of "Agreement"
  )


data_sub_final_plot=data.frame(Xenium=data_sub_final$TACIT_Xenium_v2,PCF=data_sub_final$TACIT_PCF_v2,TLS=data_sub_final$TLS)
data_sub_final_plot=data_sub_final_plot[which(data_sub_final_plot$Xenium%in%c("VEC")|
                                                data_sub_final_plot$PCF%in%c("VEC")),]
# data_sub_final_plot=data_sub_final_plot[which(data_sub_final_plot$Xenium!="Others"),]
# data_sub_final_plot=data_sub_final_plot[which(data_sub_final_plot$PCF!="Others"),]


data_sub_final_plot$score=ifelse(data_sub_final_plot$Xenium==data_sub_final_plot$PCF,"Agreement","Disagreement")

agreement_proportions_sub <- data_sub_final_plot %>%
  group_by(TLS) %>%
  summarise(
    Total = n(),  # Total counts per TLS
    Agreement_Count = sum(score == "Agreement"),  # Count of "Agreement" per TLS
    Agreement_Prop = Agreement_Count / Total  # Proportion of "Agreement"
  )


agreement_proportions=rbind(data.frame(agreement_proportions_all,Group="All markers"),data.frame(agreement_proportions_sub,Group="Subset markers"))

# Create the boxplot
ggplot(agreement_proportions, aes(x = Group, y = Agreement_Prop, fill = Group)) +
  geom_boxplot() +
  labs(title = "",
       x = "",
       y = "Agreement Proportion",
       fill = "Group") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Adjust text angle for better readability
        plot.title = element_text(hjust = 0.5),legend.position = "none")+ylim(0,1)   +
  stat_compare_means(method = "t.test", label = "p.format", size = 5) 








data_sub_final_plot=data.frame(Xenium=data_sub_final$TACIT_Xenium_UPDATE,PCF=data_sub_final$TACIT_PCF,TLS=data_sub_final$TLS)
data_sub_final_plot=data_sub_final_plot[which(data_sub_final_plot$Xenium%in%c("B Cells")|
                                                data_sub_final_plot$PCF%in%c("B Cells")),]

data_sub_final_plot$score=ifelse(data_sub_final_plot$Xenium==data_sub_final_plot$PCF,"Agreement","Disagreement")

agreement_proportions_all <- data_sub_final_plot %>%
  group_by(TLS) %>%
  summarise(
    Total = n(),  # Total counts per TLS
    Agreement_Count = sum(score == "Agreement"),  # Count of "Agreement" per TLS
    Agreement_Prop = Agreement_Count / Total  # Proportion of "Agreement"
  )


data_sub_final_plot=data.frame(Xenium=data_sub_final$TACIT_Xenium_v2,PCF=data_sub_final$TACIT_PCF_v2,TLS=data_sub_final$TLS)
data_sub_final_plot=data_sub_final_plot[which(data_sub_final_plot$Xenium%in%c("B Cells")|
                                                data_sub_final_plot$PCF%in%c("B Cells")),]


data_sub_final_plot$score=ifelse(data_sub_final_plot$Xenium==data_sub_final_plot$PCF,"Agreement","Disagreement")

agreement_proportions_sub <- data_sub_final_plot %>%
  group_by(TLS) %>%
  summarise(
    Total = n(),  # Total counts per TLS
    Agreement_Count = sum(score == "Agreement"),  # Count of "Agreement" per TLS
    Agreement_Prop = Agreement_Count / Total  # Proportion of "Agreement"
  )


agreement_proportions=rbind(data.frame(agreement_proportions_all,Group="All markers"),data.frame(agreement_proportions_sub,Group="Subset markers"))

# Create the boxplot
ggplot(agreement_proportions, aes(x = Group, y = Agreement_Prop, fill = Group)) +
  geom_boxplot() +
  labs(title = "",
       x = "",
       y = "Agreement Proportion",
       fill = "Group") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Adjust text angle for better readability
        plot.title = element_text(hjust = 0.5),legend.position = "none")+ylim(0,1)  +
  stat_compare_means(method = "t.test", label = "p.format", size = 5) 









data_sub_final_plot=data.frame(Xenium=data_sub_final$TACIT_Xenium_UPDATE,PCF=data_sub_final$TACIT_PCF,TLS=data_sub_final$TLS)
data_sub_final_plot=data_sub_final_plot[which(data_sub_final_plot$Xenium%in%c("CD4+ T Cells")|
                                                data_sub_final_plot$PCF%in%c("CD4+ T Cells")),]

data_sub_final_plot$score=ifelse(data_sub_final_plot$Xenium==data_sub_final_plot$PCF,"Agreement","Disagreement")

agreement_proportions_all <- data_sub_final_plot %>%
  group_by(TLS) %>%
  summarise(
    Total = n(),  # Total counts per TLS
    Agreement_Count = sum(score == "Agreement"),  # Count of "Agreement" per TLS
    Agreement_Prop = Agreement_Count / Total  # Proportion of "Agreement"
  )


data_sub_final_plot=data.frame(Xenium=data_sub_final$TACIT_Xenium_v2,PCF=data_sub_final$TACIT_PCF_v2,TLS=data_sub_final$TLS)
data_sub_final_plot=data_sub_final_plot[which(data_sub_final_plot$Xenium%in%c("CD4+ T Cells")|
                                                data_sub_final_plot$PCF%in%c("CD4+ T Cells")),]


data_sub_final_plot$score=ifelse(data_sub_final_plot$Xenium==data_sub_final_plot$PCF,"Agreement","Disagreement")

agreement_proportions_sub <- data_sub_final_plot %>%
  group_by(TLS) %>%
  summarise(
    Total = n(),  # Total counts per TLS
    Agreement_Count = sum(score == "Agreement"),  # Count of "Agreement" per TLS
    Agreement_Prop = Agreement_Count / Total  # Proportion of "Agreement"
  )


agreement_proportions=rbind(data.frame(agreement_proportions_all,Group="All markers"),data.frame(agreement_proportions_sub,Group="Subset markers"))

# Create the boxplot
ggplot(agreement_proportions, aes(x = Group, y = Agreement_Prop, fill = Group)) +
  geom_boxplot() +
  labs(title = "",
       x = "",
       y = "Agreement Proportion",
       fill = "Group") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Adjust text angle for better readability
        plot.title = element_text(hjust = 0.5),legend.position = "none")+ylim(0,1)  +
  stat_compare_means(method = "t.test", label = "p.format", size = 5) 









data_sub_final_plot=data.frame(Xenium=data_sub_final$TACIT_Xenium_UPDATE,PCF=data_sub_final$TACIT_PCF,TLS=data_sub_final$TLS)
data_sub_final_plot=data_sub_final_plot[which(data_sub_final_plot$Xenium%in%c("Fibroblasts")|
                                                data_sub_final_plot$PCF%in%c("Fibroblasts")),]

data_sub_final_plot$score=ifelse(data_sub_final_plot$Xenium==data_sub_final_plot$PCF,"Agreement","Disagreement")

agreement_proportions_all <- data_sub_final_plot %>%
  group_by(TLS) %>%
  summarise(
    Total = n(),  # Total counts per TLS
    Agreement_Count = sum(score == "Agreement"),  # Count of "Agreement" per TLS
    Agreement_Prop = Agreement_Count / Total  # Proportion of "Agreement"
  )


data_sub_final_plot=data.frame(Xenium=data_sub_final$TACIT_Xenium_v2,PCF=data_sub_final$TACIT_PCF_v2,TLS=data_sub_final$TLS)
data_sub_final_plot=data_sub_final_plot[which(data_sub_final_plot$Xenium%in%c("Fibroblasts")|
                                                data_sub_final_plot$PCF%in%c("Fibroblasts")),]


data_sub_final_plot$score=ifelse(data_sub_final_plot$Xenium==data_sub_final_plot$PCF,"Agreement","Disagreement")

agreement_proportions_sub <- data_sub_final_plot %>%
  group_by(TLS) %>%
  summarise(
    Total = n(),  # Total counts per TLS
    Agreement_Count = sum(score == "Agreement"),  # Count of "Agreement" per TLS
    Agreement_Prop = Agreement_Count / Total  # Proportion of "Agreement"
  )


agreement_proportions=rbind(data.frame(agreement_proportions_all,Group="All markers"),data.frame(agreement_proportions_sub,Group="Subset markers"))

# Create the boxplot
ggplot(agreement_proportions, aes(x = Group, y = Agreement_Prop, fill = Group)) +
  geom_boxplot() +
  labs(title = "",
       x = "",
       y = "Agreement Proportion",
       fill = "Group") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Adjust text angle for better readability
        plot.title = element_text(hjust = 0.5),legend.position = "none")+ylim(0,1)  +
  stat_compare_means(method = "t.test", label = "p.format", size = 5) 





data_sub_final_plot=data.frame(Xenium=data_sub_final$TACIT_Xenium_UPDATE,PCF=data_sub_final$TACIT_PCF,TLS=data_sub_final$TLS)
data_sub_final_plot=data_sub_final_plot[which(data_sub_final_plot$Xenium%in%c("Acinar Cells")|
                                                data_sub_final_plot$PCF%in%c("Acinar Cells")),]

data_sub_final_plot$score=ifelse(data_sub_final_plot$Xenium==data_sub_final_plot$PCF,"Agreement","Disagreement")

agreement_proportions_all <- data_sub_final_plot %>%
  group_by(TLS) %>%
  summarise(
    Total = n(),  # Total counts per TLS
    Agreement_Count = sum(score == "Agreement"),  # Count of "Agreement" per TLS
    Agreement_Prop = Agreement_Count / Total  # Proportion of "Agreement"
  )


data_sub_final_plot=data.frame(Xenium=data_sub_final$TACIT_Xenium_v2,PCF=data_sub_final$TACIT_PCF_v2,TLS=data_sub_final$TLS)
data_sub_final_plot=data_sub_final_plot[which(data_sub_final_plot$Xenium%in%c("Acinar Cells")|
                                                data_sub_final_plot$PCF%in%c("Acinar Cells")),]


data_sub_final_plot$score=ifelse(data_sub_final_plot$Xenium==data_sub_final_plot$PCF,"Agreement","Disagreement")

agreement_proportions_sub <- data_sub_final_plot %>%
  group_by(TLS) %>%
  summarise(
    Total = n(),  # Total counts per TLS
    Agreement_Count = sum(score == "Agreement"),  # Count of "Agreement" per TLS
    Agreement_Prop = Agreement_Count / Total  # Proportion of "Agreement"
  )


agreement_proportions=rbind(data.frame(agreement_proportions_all,Group="All markers"),data.frame(agreement_proportions_sub,Group="Subset markers"))

# Create the boxplot
ggplot(agreement_proportions, aes(x = Group, y = Agreement_Prop, fill = Group)) +
  geom_boxplot() +
  labs(title = "",
       x = "",
       y = "Agreement Proportion",
       fill = "Group") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Adjust text angle for better readability
        plot.title = element_text(hjust = 0.5),legend.position = "none")+ylim(0,1)  +
  stat_compare_means(method = "t.test", label = "p.format", size = 5) 








data_sub_final_plot=data.frame(Xenium=data_sub_final$TACIT_Xenium_UPDATE,PCF=data_sub_final$TACIT_PCF,TLS=data_sub_final$TLS)
data_sub_final_plot=data_sub_final_plot[which(data_sub_final_plot$Xenium%in%c("Ductal Cells")|
                                                data_sub_final_plot$PCF%in%c("Ductal Cells")),]

data_sub_final_plot$score=ifelse(data_sub_final_plot$Xenium==data_sub_final_plot$PCF,"Agreement","Disagreement")

agreement_proportions_all <- data_sub_final_plot %>%
  group_by(TLS) %>%
  summarise(
    Total = n(),  # Total counts per TLS
    Agreement_Count = sum(score == "Agreement"),  # Count of "Agreement" per TLS
    Agreement_Prop = Agreement_Count / Total  # Proportion of "Agreement"
  )


data_sub_final_plot=data.frame(Xenium=data_sub_final$TACIT_Xenium_v2,PCF=data_sub_final$TACIT_PCF_v2,TLS=data_sub_final$TLS)
data_sub_final_plot=data_sub_final_plot[which(data_sub_final_plot$Xenium%in%c("Ductal Cells")|
                                                data_sub_final_plot$PCF%in%c("Ductal Cells")),]


data_sub_final_plot$score=ifelse(data_sub_final_plot$Xenium==data_sub_final_plot$PCF,"Agreement","Disagreement")

agreement_proportions_sub <- data_sub_final_plot %>%
  group_by(TLS) %>%
  summarise(
    Total = n(),  # Total counts per TLS
    Agreement_Count = sum(score == "Agreement"),  # Count of "Agreement" per TLS
    Agreement_Prop = Agreement_Count / Total  # Proportion of "Agreement"
  )


agreement_proportions=rbind(data.frame(agreement_proportions_all,Group="All markers"),data.frame(agreement_proportions_sub,Group="Subset markers"))
# Create the boxplot
ggplot(agreement_proportions, aes(x = Group, y = Agreement_Prop, fill = Group)) +
  geom_boxplot() +
  labs(title = "",
       x = "",
       y = "Agreement Proportion",
       fill = "Group") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Adjust text angle for better readability
        plot.title = element_text(hjust = 0.5),legend.position = "none")+ylim(0,1)  +
  stat_compare_means(method = "t.test", label = "p.format", size = 5) 











#-------------------
#-i-----------------
#-------------------




Group1=ifelse(as.numeric(data_sub4$PDCD1)>0,paste0("PDCD1","+",sep=" "),paste0("PDCD1","-",sep=" "))
Group2=ifelse(data_sub4$TACIT_Xenium_v2%in%c("CD4+ T Cells","B Cells","VEC"),data_sub4$TACIT_Xenium_v2,"Others")

group_ct_st_2=paste0(Group2,"::",Group1,sep="")

# Reclassify specific categories using levels
group_ct_st_2 <- ifelse(grepl("Others::PDCD1-", group_ct_st_2) | grepl("Others::PDCD1+", group_ct_st_2), "Others", group_ct_st_2)
group_ct_st_2 <- ifelse(grepl("VEC::PDCD1-", group_ct_st_2) | grepl("VEC::PDCD1+", group_ct_st_2), "VEC", group_ct_st_2)


# Check the transformed table
print(table(group_ct_st_2))


group_counts <- as.data.frame(table(group_ct_st_2))
colnames(group_counts) <- c("group_ct_st", "count")
# Merge the counts back into the original data
data_plot_ct <- data.frame(X=data_sub4$X,Y=data_sub4$Y,group_ct_st=group_ct_st_2)

data_plot_ct <- merge(data_plot_ct, group_counts, by = "group_ct_st")

# Create a new column with the group names and counts
data_plot_ct$group_with_count <- data_plot_ct$group_ct_st

colors <- c(
  "B Cells::PDCD1- " = "#1f77b4",  # Blue
  "B Cells::PDCD1+ " = "#ff7f0e",   # Orange
  "CD4+ T Cells::PDCD1- " = "#2ca02c",  # Green
  "CD4+ T Cells::PDCD1+ " = "#d62728",   # Red
  "Others" = "grey90",  # Purple
  "VEC" = "#8c564b"  # Brown
)


ggplot(data_plot_ct, aes(X, Y, group = -1L)) + 
  geom_voronoi_tile(aes(fill = group_with_count), colour = 'black', max.radius = 40) +
  theme_classic(base_size = 15) +
  labs(x = "X", y = "Y", title = "") + guides(fill = guide_legend(override.aes = list(size = 3),title = "Cell type"))+
  theme(
    axis.title = element_blank(),         # Remove axis titles
    axis.text = element_blank(),          # Remove axis text
    axis.ticks = element_blank(),         # Remove axis ticks
    axis.line = element_blank(),          # Remove axis lines
    panel.background = element_blank(),   # Remove panel background
    panel.grid = element_blank()          # Remove grid lines
  )+ggtitle("")+
  scale_fill_manual(values = colors)+theme(legend.position = "bottom")







Group1=ifelse(as.numeric(data_sub4$PDCD1)>0,paste0("PDCD1","+",sep=" "),paste0("PDCD1","-",sep=" "))
Group2=ifelse(data_sub4$TACIT_PCF_v2%in%c("CD4+ T Cells","B Cells","VEC"),data_sub4$TACIT_PCF_v2,"Others")

group_ct_st_2=paste0(Group2,"::",Group1,sep="")

# Reclassify specific categories using levels
group_ct_st_2 <- ifelse(grepl("Others::PDCD1-", group_ct_st_2) | grepl("Others::PDCD1+", group_ct_st_2), "Others", group_ct_st_2)
group_ct_st_2 <- ifelse(grepl("VEC::PDCD1-", group_ct_st_2) | grepl("VEC::PDCD1+", group_ct_st_2), "VEC", group_ct_st_2)


# Check the transformed table
print(table(group_ct_st_2))


group_counts <- as.data.frame(table(group_ct_st_2))
colnames(group_counts) <- c("group_ct_st", "count")
# Merge the counts back into the original data
data_plot_ct <- data.frame(X=data_sub4$X,Y=data_sub4$Y,group_ct_st=group_ct_st_2)

data_plot_ct <- merge(data_plot_ct, group_counts, by = "group_ct_st")

# Create a new column with the group names and counts
data_plot_ct$group_with_count <- data_plot_ct$group_ct_st

colors <- c(
  "B Cells::PDCD1- " = "#1f77b4",  # Blue
  "B Cells::PDCD1+ " = "#ff7f0e",   # Orange
  "CD4+ T Cells::PDCD1- " = "#2ca02c",  # Green
  "CD4+ T Cells::PDCD1+ " = "#d62728",   # Red
  "Others" = "grey90",  # Purple
  "VEC" = "#8c564b"  # Brown
)


ggplot(data_plot_ct, aes(X, Y, group = -1L)) + 
  geom_voronoi_tile(aes(fill = group_with_count), colour = 'black', max.radius = 40) +
  theme_classic(base_size = 15) +
  labs(x = "X", y = "Y", title = "") + guides(fill = guide_legend(override.aes = list(size = 3),title = "Cell type"))+
  theme(
    axis.title = element_blank(),         # Remove axis titles
    axis.text = element_blank(),          # Remove axis text
    axis.ticks = element_blank(),         # Remove axis ticks
    axis.line = element_blank(),          # Remove axis lines
    panel.background = element_blank(),   # Remove panel background
    panel.grid = element_blank()          # Remove grid lines
  )+ggtitle("")+
  scale_fill_manual(values = colors)+theme(legend.position = "bottom")






Group1=ifelse(as.numeric(data_sub4$PDCD1)>0,paste0("PDCD1","+",sep=" "),paste0("PDCD1","-",sep=" "))
Group2=ifelse(data_sub4$TACIT_PCF_Xenium_v2%in%c("CD4+ T Cells","B Cells","VEC"),data_sub4$TACIT_PCF_Xenium_v2,"Others")

group_ct_st_2=paste0(Group2,"::",Group1,sep="")

# Reclassify specific categories using levels
group_ct_st_2 <- ifelse(grepl("Others::PDCD1-", group_ct_st_2) | grepl("Others::PDCD1+", group_ct_st_2), "Others", group_ct_st_2)
group_ct_st_2 <- ifelse(grepl("VEC::PDCD1-", group_ct_st_2) | grepl("VEC::PDCD1+", group_ct_st_2), "VEC", group_ct_st_2)


# Check the transformed table
print(table(group_ct_st_2))


group_counts <- as.data.frame(table(group_ct_st_2))
colnames(group_counts) <- c("group_ct_st", "count")
# Merge the counts back into the original data
data_plot_ct <- data.frame(X=data_sub4$X,Y=data_sub4$Y,group_ct_st=group_ct_st_2)

data_plot_ct <- merge(data_plot_ct, group_counts, by = "group_ct_st")

# Create a new column with the group names and counts
data_plot_ct$group_with_count <- data_plot_ct$group_ct_st

colors <- c(
  "B Cells::PDCD1- " = "#1f77b4",  # Blue
  "B Cells::PDCD1+ " = "#ff7f0e",   # Orange
  "CD4+ T Cells::PDCD1- " = "#2ca02c",  # Green
  "CD4+ T Cells::PDCD1+ " = "#d62728",   # Red
  "Others" = "grey90",  # Purple
  "VEC" = "#8c564b"  # Brown
)


ggplot(data_plot_ct, aes(X, Y, group = -1L)) + 
  geom_voronoi_tile(aes(fill = group_with_count), colour = 'black', max.radius = 40) +
  theme_classic(base_size = 15) +
  labs(x = "X", y = "Y", title = "") + guides(fill = guide_legend(override.aes = list(size = 3),title = "Cell type"))+
  theme(
    axis.title = element_blank(),         # Remove axis titles
    axis.text = element_blank(),          # Remove axis text
    axis.ticks = element_blank(),         # Remove axis ticks
    axis.line = element_blank(),          # Remove axis lines
    panel.background = element_blank(),   # Remove panel background
    panel.grid = element_blank()          # Remove grid lines
  )+ggtitle("")+
  scale_fill_manual(values = colors)+theme(legend.position = "bottom")















Group1=ifelse(as.numeric(data_sub4$P_PD.1)>0,paste0("PD-1","+",sep=" "),paste0("PD-1","-",sep=" "))
Group2=ifelse(data_sub4$TACIT_Xenium_v2%in%c("CD4+ T Cells","B Cells","VEC"),data_sub4$TACIT_Xenium_v2,"Others")

group_ct_st_2=paste0(Group2,"::",Group1,sep="")

# Reclassify specific categories using levels
group_ct_st_2 <- ifelse(grepl("Others::PD-1-", group_ct_st_2) | grepl("Others::PD-1+", group_ct_st_2), "Others", group_ct_st_2)
group_ct_st_2 <- ifelse(grepl("VEC::PD-1-", group_ct_st_2) | grepl("VEC::PD-1+", group_ct_st_2), "VEC", group_ct_st_2)


# Check the transformed table
print(table(group_ct_st_2))


group_counts <- as.data.frame(table(group_ct_st_2))
colnames(group_counts) <- c("group_ct_st", "count")
# Merge the counts back into the original data
data_plot_ct <- data.frame(X=data_sub4$X,Y=data_sub4$Y,group_ct_st=group_ct_st_2)

data_plot_ct <- merge(data_plot_ct, group_counts, by = "group_ct_st")

# Create a new column with the group names and counts
data_plot_ct$group_with_count <- data_plot_ct$group_ct_st

colors <- c(
  "B Cells::PD-1- " = "#1f77b4",  # Blue
  "B Cells::PD-1+ " = "#ff7f0e",   # Orange
  "CD4+ T Cells::PD-1- " = "#2ca02c",  # Green
  "CD4+ T Cells::PD-1+ " = "#d62728",   # Red
  "Others" = "grey90",  # Purple
  "VEC" = "#8c564b"  # Brown
)


ggplot(data_plot_ct, aes(X, Y, group = -1L)) + 
  geom_voronoi_tile(aes(fill = group_with_count), colour = 'black', max.radius = 40) +
  theme_classic(base_size = 15) +
  labs(x = "X", y = "Y", title = "") + guides(fill = guide_legend(override.aes = list(size = 3),title = "Cell type"))+
  theme(
    axis.title = element_blank(),         # Remove axis titles
    axis.text = element_blank(),          # Remove axis text
    axis.ticks = element_blank(),         # Remove axis ticks
    axis.line = element_blank(),          # Remove axis lines
    panel.background = element_blank(),   # Remove panel background
    panel.grid = element_blank()          # Remove grid lines
  )+ggtitle("")+
  scale_fill_manual(values = colors)+theme(legend.position = "bottom")







Group1=ifelse(as.numeric(data_sub4$P_PD.1)>0,paste0("PD-1","+",sep=" "),paste0("PD-1","-",sep=" "))
Group2=ifelse(data_sub4$TACIT_PCF_v2%in%c("CD4+ T Cells","B Cells","VEC"),data_sub4$TACIT_PCF_v2,"Others")

group_ct_st_2=paste0(Group2,"::",Group1,sep="")

# Reclassify specific categories using levels
group_ct_st_2 <- ifelse(grepl("Others::PD-1-", group_ct_st_2) | grepl("Others::PD-1+", group_ct_st_2), "Others", group_ct_st_2)
group_ct_st_2 <- ifelse(grepl("VEC::PD-1-", group_ct_st_2) | grepl("VEC::PD-1+", group_ct_st_2), "VEC", group_ct_st_2)


# Check the transformed table
print(table(group_ct_st_2))


group_counts <- as.data.frame(table(group_ct_st_2))
colnames(group_counts) <- c("group_ct_st", "count")
# Merge the counts back into the original data
data_plot_ct <- data.frame(X=data_sub4$X,Y=data_sub4$Y,group_ct_st=group_ct_st_2)

data_plot_ct <- merge(data_plot_ct, group_counts, by = "group_ct_st")

# Create a new column with the group names and counts
data_plot_ct$group_with_count <- data_plot_ct$group_ct_st

colors <- c(
  "B Cells::PD-1- " = "#1f77b4",  # Blue
  "B Cells::PD-1+ " = "#ff7f0e",   # Orange
  "CD4+ T Cells::PD-1- " = "#2ca02c",  # Green
  "CD4+ T Cells::PD-1+ " = "#d62728",   # Red
  "Others" = "grey90",  # Purple
  "VEC" = "#8c564b"  # Brown
)


ggplot(data_plot_ct, aes(X, Y, group = -1L)) + 
  geom_voronoi_tile(aes(fill = group_with_count), colour = 'black', max.radius = 40) +
  theme_classic(base_size = 15) +
  labs(x = "X", y = "Y", title = "") + guides(fill = guide_legend(override.aes = list(size = 3),title = "Cell type"))+
  theme(
    axis.title = element_blank(),         # Remove axis titles
    axis.text = element_blank(),          # Remove axis text
    axis.ticks = element_blank(),         # Remove axis ticks
    axis.line = element_blank(),          # Remove axis lines
    panel.background = element_blank(),   # Remove panel background
    panel.grid = element_blank()          # Remove grid lines
  )+ggtitle("")+
  scale_fill_manual(values = colors)+theme(legend.position = "bottom")






Group1=ifelse(as.numeric(data_sub4$P_PD.1)>0,paste0("PD-1","+",sep=" "),paste0("PD-1","-",sep=" "))
Group2=ifelse(data_sub4$TACIT_PCF_Xenium_v2%in%c("CD4+ T Cells","B Cells","VEC"),data_sub4$TACIT_PCF_Xenium_v2,"Others")

group_ct_st_2=paste0(Group2,"::",Group1,sep="")

# Reclassify specific categories using levels
group_ct_st_2 <- ifelse(grepl("Others::PD-1-", group_ct_st_2) | grepl("Others::PD-1+", group_ct_st_2), "Others", group_ct_st_2)
group_ct_st_2 <- ifelse(grepl("VEC::PD-1-", group_ct_st_2) | grepl("VEC::PD-1+", group_ct_st_2), "VEC", group_ct_st_2)


# Check the transformed table
print(table(group_ct_st_2))


group_counts <- as.data.frame(table(group_ct_st_2))
colnames(group_counts) <- c("group_ct_st", "count")
# Merge the counts back into the original data
data_plot_ct <- data.frame(X=data_sub4$X,Y=data_sub4$Y,group_ct_st=group_ct_st_2)

data_plot_ct <- merge(data_plot_ct, group_counts, by = "group_ct_st")

# Create a new column with the group names and counts
data_plot_ct$group_with_count <- data_plot_ct$group_ct_st

colors <- c(
  "B Cells::PD-1- " = "#1f77b4",  # Blue
  "B Cells::PD-1+ " = "#ff7f0e",   # Orange
  "CD4+ T Cells::PD-1- " = "#2ca02c",  # Green
  "CD4+ T Cells::PD-1+ " = "#d62728",   # Red
  "Others" = "grey90",  # Purple
  "VEC" = "#8c564b"  # Brown
)


ggplot(data_plot_ct, aes(X, Y, group = -1L)) + 
  geom_voronoi_tile(aes(fill = group_with_count), colour = 'black', max.radius = 40) +
  theme_classic(base_size = 15) +
  labs(x = "X", y = "Y", title = "") + guides(fill = guide_legend(override.aes = list(size = 3),title = "Cell type"))+
  theme(
    axis.title = element_blank(),         # Remove axis titles
    axis.text = element_blank(),          # Remove axis text
    axis.ticks = element_blank(),         # Remove axis ticks
    axis.line = element_blank(),          # Remove axis lines
    panel.background = element_blank(),   # Remove panel background
    panel.grid = element_blank()          # Remove grid lines
  )+ggtitle("")+
  scale_fill_manual(values = colors)+theme(legend.position = "bottom")










