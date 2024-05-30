data=read.csv("~/Downloads/PCF_Xenium_reg004.csv")

data$TACIT_PCF_Xenium_v2=ifelse(data$TACIT_PCF_v2=="B Cells","B Cells",data$TACIT_PCF_Xenium_v2)

data_sub=data[which(data$Y<2700&data$Y>2000&data$X<3700&data$X>3300),]


user_colors <-  c(
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
  "Myfibroblast"="#f032e6",
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
  "Myoepithelial"="blue",
  "NK Cells"="darkred",
  "Pericytes"="#FFC0CB",
  "Memory T cell"="brown3",
  "Regulatory T Cells"="brown2",
  "Treg"="brown2",
  "Seromucous Acinar Cells"="royalblue2",
  "Smooth Muscle"="azure3",
  "T Cell Progenitors"="brown3",
  "Venules"="orange3","VEC"="orange3","VEC Progen"="orange1","LECs"="orange","Others"="grey90","Neutrophil"="gold","Lymphoid"="darkcyan")

ggplot(data_sub, aes(X, Y, group = -1L)) + 
  geom_voronoi_tile(aes(fill = Seurat_orginial), colour = 'black', max.radius = 40)+
  scale_fill_manual(values = user_colors) +
  theme_classic(base_size = 15) +
  labs(x = "X", y = "Y", title = "") + guides(fill = guide_legend(override.aes = list(size = 3),title = "Cell type"))+
  theme(
    axis.title = element_blank(),         # Remove axis titles
    axis.text = element_blank(),          # Remove axis text
    axis.ticks = element_blank(),         # Remove axis ticks
    axis.line = element_blank(),          # Remove axis lines
    panel.background = element_blank(),   # Remove panel background
    panel.grid = element_blank()          # Remove grid lines
  )+ggtitle("PCF")


ggplot(data_sub, aes(X, Y, group = -1L)) + 
  geom_voronoi_tile(aes(fill = data_sub$TACIT_Xenium_v2), colour = 'black', max.radius = 40)+
  scale_fill_manual(values = user_colors) +
  theme_classic(base_size = 15) +
  labs(x = "X", y = "Y", title = "") + guides(fill = guide_legend(override.aes = list(size = 3),title = "Cell type"))+
  theme(
    axis.title = element_blank(),         # Remove axis titles
    axis.text = element_blank(),          # Remove axis text
    axis.ticks = element_blank(),         # Remove axis ticks
    axis.line = element_blank(),          # Remove axis lines
    panel.background = element_blank(),   # Remove panel background
    panel.grid = element_blank()   ,legend.position = "none"       # Remove grid lines
  )+ggtitle("")

ggplot(data_sub, aes(X, Y, group = -1L)) + 
  geom_voronoi_tile(aes(fill = data_sub$TACIT_PCF_v2), colour = 'black', max.radius = 40)+
  scale_fill_manual(values = user_colors) +
  theme_classic(base_size = 15) +
  labs(x = "X", y = "Y", title = "") + guides(fill = guide_legend(override.aes = list(size = 3),title = "Cell type"))+
  theme(
    axis.title = element_blank(),         # Remove axis titles
    axis.text = element_blank(),          # Remove axis text
    axis.ticks = element_blank(),         # Remove axis ticks
    axis.line = element_blank(),          # Remove axis lines
    panel.background = element_blank(),   # Remove panel background
    panel.grid = element_blank()   ,legend.position = "none"       # Remove grid lines
  )+ggtitle("")


Group=ifelse(data_sub$TACIT_PCF!=data_sub$TACIT_Xenium_UPDATE,"Disagreement","Agreement")

ggplot(data_sub, aes(X, Y, group = -1L)) + 
  geom_voronoi_tile(aes(fill = Group), colour = 'black', max.radius = 40)+
  scale_fill_manual(values = c("Disagreement"="red","Agreement"="grey90")) +
  theme_classic(base_size = 15) +
  labs(x = "X", y = "Y", title = "") + guides(fill = guide_legend(override.aes = list(size = 3),title = "Cell type"))+
  theme(
    axis.title = element_blank(),         # Remove axis titles
    axis.text = element_blank(),          # Remove axis text
    axis.ticks = element_blank(),         # Remove axis ticks
    axis.line = element_blank(),          # Remove axis lines
    panel.background = element_blank(),   # Remove panel background
    panel.grid = element_blank(),legend.position = "none"          # Remove grid lines
  )





data_sub=data[which(data$Y<2700&data$Y>2000&data$X<3700&data$X>3300),]
data_sub$TACIT_Xenium_UPDATE=ifelse(data_sub$TACIT_Xenium%in%c("Arterioles","Venules"),"VEC",data_sub$TACIT_Xenium)
data_sub$TACIT_Xenium_UPDATE=ifelse(data_sub$TACIT_Xenium_UPDATE%in%c("Ductal Cells","Ductal Progenitors"),"Ductal Cells",data_sub$TACIT_Xenium_UPDATE)
data_sub$TACIT_Xenium_UPDATE=ifelse(data_sub$TACIT_Xenium_UPDATE%in%c("Mucous Acinar Cells","Seromucous Acinar Cells"),"Acinar Cells",data_sub$TACIT_Xenium_UPDATE)
data_sub$TACIT_Xenium_UPDATE=ifelse(data_sub$TACIT_Xenium_UPDATE%in%c("CD8+ Exhausted T Cells"),"CD8+ T Cells",data_sub$TACIT_Xenium_UPDATE)
data_sub$TACIT_Xenium_UPDATE=ifelse(data_sub$TACIT_Xenium_UPDATE%in%c("Myoepithelium","Intermediate Epithelium"),"Myoepithelial",data_sub$TACIT_Xenium_UPDATE)
data_sub$TACIT_Xenium_UPDATE=ifelse(data_sub$TACIT_Xenium_UPDATE%in%unique(data_sub$TACIT_PCF),data_sub$TACIT_Xenium_UPDATE,"Others")


ggplot(data_sub, aes(X, Y, group = -1L)) + 
  geom_voronoi_tile(aes(fill = TACIT_PCF), colour = 'black', max.radius = 40)+
  scale_fill_manual(values = user_colors) +
  theme_classic(base_size = 15) +
  labs(x = "X", y = "Y", title = "") + guides(fill = guide_legend(override.aes = list(size = 3),title = "Cell type"))+
  theme(
    axis.title = element_blank(),         # Remove axis titles
    axis.text = element_blank(),          # Remove axis text
    axis.ticks = element_blank(),         # Remove axis ticks
    axis.line = element_blank(),          # Remove axis lines
    panel.background = element_blank(),   # Remove panel background
    panel.grid = element_blank(),legend.position = "none"          # Remove grid lines
  )+ggtitle("")


ggplot(data_sub, aes(X, Y, group = -1L)) + 
  geom_voronoi_tile(aes(fill = data_sub$TACIT_PCF), colour = 'black', max.radius = 40)+
  scale_fill_manual(values = user_colors) +
  theme_classic(base_size = 15) +
  labs(x = "X", y = "Y", title = "") + guides(fill = guide_legend(override.aes = list(size = 3),title = "Cell type"))+
  theme(
    axis.title = element_blank(),         # Remove axis titles
    axis.text = element_blank(),          # Remove axis text
    axis.ticks = element_blank(),         # Remove axis ticks
    axis.line = element_blank(),          # Remove axis lines
    panel.background = element_blank(),   # Remove panel background
    panel.grid = element_blank(),legend.position = "none"          # Remove grid lines
  )+ggtitle("")




ggplot(data_sub, aes(X, Y, group = -1L)) + 
  geom_voronoi_tile(aes(fill = data_sub$TACIT_PCF_Xenium_v2), colour = 'black', max.radius = 40)+
  scale_fill_manual(values = user_colors) +
  theme_classic(base_size = 15) +
  labs(x = "X", y = "Y", title = "") + guides(fill = guide_legend(override.aes = list(size = 3),title = "Cell type"))+
  theme(
    axis.title = element_blank(),         # Remove axis titles
    axis.text = element_blank(),          # Remove axis text
    axis.ticks = element_blank(),         # Remove axis ticks
    axis.line = element_blank(),          # Remove axis lines
    panel.background = element_blank(),   # Remove panel background
    panel.grid = element_blank()          # Remove grid lines
  )+ggtitle("PCF & Xenium")




data_sub=data[which(data$Y<3565&data$Y>3180&data$X<1150&data$X>670),]


ggplot(data_sub, aes(X, Y, group = -1L)) + 
  geom_voronoi_tile(aes(fill = Seurat_orginial), colour = 'black', max.radius = 40)+
  scale_fill_manual(values = user_colors) +
  theme_classic(base_size = 15) +
  labs(x = "X", y = "Y", title = "") + guides(fill = guide_legend(override.aes = list(size = 3),title = "Cell type"))+
  theme(
    axis.title = element_blank(),         # Remove axis titles
    axis.text = element_blank(),          # Remove axis text
    axis.ticks = element_blank(),         # Remove axis ticks
    axis.line = element_blank(),          # Remove axis lines
    panel.background = element_blank(),   # Remove panel background
    panel.grid = element_blank()          # Remove grid lines
  )+ggtitle("PCF")


ggplot(data_sub, aes(X, Y, group = -1L)) + 
  geom_voronoi_tile(aes(fill = data_sub$TACIT_Xenium_v2), colour = 'black', max.radius = 40)+
  scale_fill_manual(values = user_colors) +
  theme_classic(base_size = 15) +
  labs(x = "X", y = "Y", title = "") + guides(fill = guide_legend(override.aes = list(size = 3),title = "Cell type"))+
  theme(
    axis.title = element_blank(),         # Remove axis titles
    axis.text = element_blank(),          # Remove axis text
    axis.ticks = element_blank(),         # Remove axis ticks
    axis.line = element_blank(),          # Remove axis lines
    panel.background = element_blank(),   # Remove panel background
    panel.grid = element_blank()          # Remove grid lines
  )+ggtitle("Xenium")




ggplot(data_sub, aes(X, Y, group = -1L)) + 
  geom_voronoi_tile(aes(fill = data_sub$TACIT_PCF_Xenium_v2), colour = 'black', max.radius = 40)+
  scale_fill_manual(values = user_colors) +
  theme_classic(base_size = 15) +
  labs(x = "X", y = "Y", title = "") + guides(fill = guide_legend(override.aes = list(size = 3),title = "Cell type"))+
  theme(
    axis.title = element_blank(),         # Remove axis titles
    axis.text = element_blank(),          # Remove axis text
    axis.ticks = element_blank(),         # Remove axis ticks
    axis.line = element_blank(),          # Remove axis lines
    panel.background = element_blank(),   # Remove panel background
    panel.grid = element_blank()          # Remove grid lines
  )+ggtitle("PCF & Xenium")







data_sub=data[which(data$Y<4400&data$Y>4075&data$X<1430&data$X>1160),]


ggplot(data_sub, aes(X, Y, group = -1L)) + 
  geom_voronoi_tile(aes(fill = TACIT_PCF_v2), colour = 'black', max.radius = 40)+
  scale_fill_manual(values = user_colors) +
  theme_classic(base_size = 15) +
  labs(x = "X", y = "Y", title = "") + guides(fill = guide_legend(override.aes = list(size = 3),title = "Cell type"))+
  theme(
    axis.title = element_blank(),         # Remove axis titles
    axis.text = element_blank(),          # Remove axis text
    axis.ticks = element_blank(),         # Remove axis ticks
    axis.line = element_blank(),          # Remove axis lines
    panel.background = element_blank(),   # Remove panel background
    panel.grid = element_blank()          # Remove grid lines
  )+ggtitle("")+theme(legend.position = "none")


ggplot(data_sub, aes(X, Y, group = -1L)) + 
  geom_voronoi_tile(aes(fill = data_sub$TACIT_PCF_Xenium_v2), colour = 'black', max.radius = 40)+
  scale_fill_manual(values = user_colors) +
  theme_classic(base_size = 15) +
  labs(x = "X", y = "Y", title = "") + guides(fill = guide_legend(override.aes = list(size = 3),title = "Cell type"))+
  theme(
    axis.title = element_blank(),         # Remove axis titles
    axis.text = element_blank(),          # Remove axis text
    axis.ticks = element_blank(),         # Remove axis ticks
    axis.line = element_blank(),          # Remove axis lines
    panel.background = element_blank(),   # Remove panel background
    panel.grid = element_blank()          # Remove grid lines
  )+ggtitle("")+theme(legend.position = "none")




ggplot(data_sub, aes(X, Y, group = -1L)) + 
  geom_voronoi_tile(aes(fill = data_sub$TACIT_PCF_Xenium_v2), colour = 'black', max.radius = 40)+
  scale_fill_manual(values = user_colors) +
  theme_classic(base_size = 15) +
  labs(x = "X", y = "Y", title = "") + guides(fill = guide_legend(override.aes = list(size = 3),title = "Cell type"))+
  theme(
    axis.title = element_blank(),         # Remove axis titles
    axis.text = element_blank(),          # Remove axis text
    axis.ticks = element_blank(),         # Remove axis ticks
    axis.line = element_blank(),          # Remove axis lines
    panel.background = element_blank(),   # Remove panel background
    panel.grid = element_blank()          # Remove grid lines
  )+ggtitle("PCF & Xenium")



















data1=read.csv("~/Downloads/PCF_Xenium_reg001.csv")
data2=read.csv("~/Downloads/PCF_Xenium_reg002.csv")
data3=read.csv("~/Downloads/PCF_Xenium_reg003.csv")
data4=read.csv("~/Downloads/PCF_Xenium_reg004.csv")
data5=read.csv("~/Downloads/PCF_Xenium_reg005.csv")
data6=read.csv("~/Downloads/PCF_Xenium_reg006.csv")



Signature=read.csv("~/Downloads/Signature_PCF_Xenium.csv")



# Scatter plot using ggplot2 for better visualization
data_plot <- data.frame(continuous_var=(data4$CD20), count_var=(data4$MS4A1))
# Calculate the correlation
correlation <- cor(data_plot$continuous_var, (data4$MS4A1), method = "pearson")

# Create the scatter plot with a regression line and correlation value
ggplot(data_plot, aes(x = continuous_var, y = count_var)) +
  geom_point(color = "blue") +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  labs(title = "",
       x = "PCF",
       y = "Xenium") +
  theme_minimal() +
  annotate("text", x = Inf, y = Inf, label = paste("Pearson's r =", round(correlation, 2)), 
           hjust = 1.1, vjust = 2, size = 5, color = "black")






Signature=read.csv("~/Downloads/Signature_PCF_Xenium.csv")
Signature=Signature[,-2]
library(dplyr)


# Define the new order of cell types
new_order <- c(
  "Seromucous Acinar Cells", "Mucous Acinar Cells","Myoepithelium","Intermediate Epithelium",
  "Ductal Cells", "Ductal Progenitors","Ionocytes",
  "Fibroblasts", "Myofibroblast",
  "Capillaries", "Venules", "Arterioles", "Pericytes",
   "Dendritic Cells", "T Cell Progenitors","CD4+ T Cells","Memory T cell",
  "CD8+ Effector T Cells","CD8+ Exhausted T Cells",  "M1 Macrophages","M2 Macrophages", "NK Cells", 
  "B Cells","IgA Plasma Cells","IgG Plasma Cells", 
  "Neutrophil",
  "Smooth Muscle", "Lymphoid"
)

# Reorder the Signature$cell_type vector
Signature$cell_type <- factor(Signature$cell_type, levels = new_order)
Signature=Signature[order(Signature$cell_type),]


library(dplyr)
# Group the data by the predicted label and calculate the mean for each column
data_plot=data.frame(TACIT=data4$TACIT_Xenium,(data4[,colnames(Signature)[-c(1)]]))
colnames(data_plot)[-1]=colnames(Signature)[-c(1)]





mean_values_TACIT <- data_plot %>%
  group_by(TACIT) %>%
  summarise_all(~quantile(., 0.5))
#    mean_values_TACIT[ , -1][mean_values_TACIT[ , -1] > 0] <- 1
mean_values_TACIT=as.data.frame(mean_values_TACIT)
rownames(mean_values_TACIT)=mean_values_TACIT$TACIT
mean_values_TACIT=mean_values_TACIT[which(mean_values_TACIT$TACIT!="Others"),]



Signature=reorder_columns_based_on_signature(Signature)

mean_values_TACIT <- as.data.frame((mean_values_TACIT[,-1]))

# Assuming mean_values_TACIT is a data frame or matrix and Signature_order is a data frame or list
ordered_index <- match(Signature$cell_type,rownames(mean_values_TACIT))

# Reorder mean_values_TACIT according to the ordered index
mean_values_TACIT <- mean_values_TACIT[ordered_index, ]

my.breaks <- c(seq(-2, 0, by=0.1),seq(0.1, 2, by=0.1))
my.colors <- c(colorRampPalette(colors = c("blue", "white"))(length(my.breaks)/2), colorRampPalette(colors = c("white", "red"))(length(my.breaks)/2))


aa=scale(mean_values_TACIT[, colnames(Signature)[-1]])
rownames(aa)=Signature$cell_type
pheatmap(
  aa,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  show_colnames = TRUE,
  fontsize_col = 15,
  fontsize_row = 15# Adjust this value as needed to increase x-axis text size
)

Signature_plot=as.matrix(Signature[,-1])
rownames(Signature_plot)=Signature$cell_type
pheatmap(
  Signature_plot,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  show_colnames = TRUE,
  fontsize_col = 15,
  fontsize_row = 15# Adjust this value as needed to increase x-axis text size
)







library(ggplot2)
library(ggpubr)

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



data_sub1=data[which(data$Y<2700&data$Y>2000&data$X<3700&data$X>3300),]
data_sub2=data[which(data$Y<974&data$Y>725&data$X<1680&data$X>1450),]
data_sub3=data[which(data$Y<3565&data$Y>3180&data$X<1150&data$X>670),]
data_sub4=data[which(data$Y<4400&data$Y>4075&data$X<1430&data$X>1160),]



# List of datasets
data_list <- list(data1 = data_sub1, data2 = data_sub2, data3 = data_sub3, data4 = data_sub3)
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


ggplot(all_correlation_data, aes(x = Gene_Pair, y = Correlation, fill = Gene_Pair)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, size = 2) +
  labs(title = "",
       x = "Category",
       y = "Correlation") +
  theme_classic(base_size = 20)    +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "none")







Group1=ifelse(as.numeric(data_sub4$PDCD1)>0,paste0("PDCD1","+",sep=" "),paste0("PDCD1","-",sep=" "))
Group2=ifelse(data_sub4$TACIT_PCF_v2%in%c("B Cells","CD4+ T Cells","VEC"),data_sub4$TACIT_PCF_v2,"Others")

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
  scale_fill_manual(values = colors)+theme(legend.position = "none")









Group1=ifelse(as.numeric(data_sub1$P_PD.1)>0,paste0("PD-1","+",sep=" "),paste0("PD-1","-",sep=" "))
Group2=ifelse(data_sub1$TACIT_Xenium_v2%in%c("CD4+ T Cells","B Cells",
                                          "VEC"),data_sub1$TACIT_Xenium_v2,"Others")

group_ct_st_2=paste0(Group2,"::",Group1,sep="")

# Reclassify specific categories using levels
group_ct_st_2 <- ifelse(grepl("Others::PD-1-", group_ct_st_2) | grepl("Others::PD-1+", group_ct_st_2), "Others", group_ct_st_2)
group_ct_st_2 <- ifelse(grepl("VEC::PD-1-", group_ct_st_2) | grepl("VEC::PD-1+", group_ct_st_2), "VEC", group_ct_st_2)


# Check the transformed table
print(table(group_ct_st_2))

group_counts <- as.data.frame(table(group_ct_st_2))
colnames(group_counts) <- c("group_ct_st", "count")
# Merge the counts back into the original data
data_plot_ct <- data.frame(X=data_sub1$X,Y=data_sub1$Y,group_ct_st=group_ct_st_2)

data_plot_ct <- merge(data_plot_ct, group_counts, by = "group_ct_st")

# Create a new column with the group names and counts
data_plot_ct$group_with_count <- paste(data_plot_ct$group_ct_st, "(n =", data_plot_ct$count, ")")



colors <- c(
  "B Cells::PD-1-  (n = 596 )" = "#1f77b4",  # Blue
  "B Cells::PD-1+  (n = 50 )" = "#ff7f0e",   # Orange
  "CD4+ T Cells::PD-1-  (n = 143 )" = "#2ca02c",  # Green
  "CD4+ T Cells::PD-1+  (n = 76 )" = "#d62728",   # Red
  "Others (n = 2011 )" = "grey90",  # Purple
  "VEC (n = 168 )" = "#8c564b"  # Brown
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
  )+ggtitle("Xenium")+
  scale_fill_manual(values = colors)








Group1=ifelse(as.numeric(data_sub1$P_PD.1)>0,paste0("PD-1","+",sep=" "),paste0("PD-1","-",sep=" "))
Group2=ifelse(data_sub1$TACIT_PCF_Xenium_v2%in%c("CD4+ T Cells","B Cells",
                                             "VEC"),data_sub1$TACIT_PCF_Xenium_v2,"Others")

group_ct_st_2=paste0(Group2,"::",Group1,sep="")

# Reclassify specific categories using levels
group_ct_st_2 <- ifelse(grepl("Others::PD-1-", group_ct_st_2) | grepl("Others::PD-1+", group_ct_st_2), "Others", group_ct_st_2)
group_ct_st_2 <- ifelse(grepl("VEC::PD-1-", group_ct_st_2) | grepl("VEC::PD-1+", group_ct_st_2), "VEC", group_ct_st_2)


# Check the transformed table
print(table(group_ct_st_2))

group_counts <- as.data.frame(table(group_ct_st_2))
colnames(group_counts) <- c("group_ct_st", "count")
# Merge the counts back into the original data
data_plot_ct <- data.frame(X=data_sub1$X,Y=data_sub1$Y,group_ct_st=group_ct_st_2)

data_plot_ct <- merge(data_plot_ct, group_counts, by = "group_ct_st")

# Create a new column with the group names and counts
data_plot_ct$group_with_count <- paste(data_plot_ct$group_ct_st, "(n =", data_plot_ct$count, ")")



colors <- c(
  "B Cells::PD-1-  (n = 754 )" = "#1f77b4",  # Blue
  "B Cells::PD-1+  (n = 98 )" = "#ff7f0e",   # Orange
  "CD4+ T Cells::PD-1-  (n = 165 )" = "#2ca02c",  # Green
  "CD4+ T Cells::PD-1+  (n = 135 )" = "#d62728",   # Red
  "Others (n = 1786 )" = "grey90",  # Purple
  "VEC (n = 106 )" = "#8c564b"  # Brown
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
  )+ggtitle("PCF & Xenium")+
  scale_fill_manual(values = colors)









Group1=ifelse(as.numeric(data_sub2$P_PD.1)>0,paste0("PD-1","+",sep=" "),paste0("PD-1","-",sep=" "))
Group2=ifelse(data_sub2$TACIT_PCF_Xenium_v2%in%c("CD4+ T Cells","B Cells",
                                                 "VEC"),data_sub2$TACIT_PCF_Xenium_v2,"Others")

group_ct_st_2=paste0(Group2,"::",Group1,sep="")

# Reclassify specific categories using levels
group_ct_st_2 <- ifelse(grepl("Others::PD-1-", group_ct_st_2) | grepl("Others::PD-1+", group_ct_st_2), "Others", group_ct_st_2)
group_ct_st_2 <- ifelse(grepl("VEC::PD-1-", group_ct_st_2) | grepl("VEC::PD-1+", group_ct_st_2), "VEC", group_ct_st_2)


# Check the transformed table
print(table(group_ct_st_2))

group_counts <- as.data.frame(table(group_ct_st_2))
colnames(group_counts) <- c("group_ct_st", "count")
# Merge the counts back into the original data
data_plot_ct <- data.frame(X=data_sub2$X,Y=data_sub2$Y,group_ct_st=group_ct_st_2)

data_plot_ct <- merge(data_plot_ct, group_counts, by = "group_ct_st")

# Create a new column with the group names and counts
data_plot_ct$group_with_count <- paste(data_plot_ct$group_ct_st, "(n =", data_plot_ct$count, ")")



colors <- c(
  "B Cells::PD-1-  (n = 215 )" = "#1f77b4",  # Blue
  "B Cells::PD-1+  (n = 131 )" = "#ff7f0e",   # Orange
  "CD4+ T Cells::PD-1-  (n = 66 )" = "#2ca02c",  # Green
  "CD4+ T Cells::PD-1+  (n = 75 )" = "#d62728",   # Red
  "Others (n = 260 )" = "grey90",  # Purple
  "VEC (n = 20 )" = "#8c564b"  # Brown
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
  )+ggtitle("PCF & Xenium")+
  scale_fill_manual(values = colors)












Group1=ifelse(as.numeric(data_sub2$P_PD.1)>0,paste0("PD-1","+",sep=" "),paste0("PD-1","-",sep=" "))
Group2=ifelse(data_sub2$TACIT_PCF_v2%in%c("CD4+ T Cells","B Cells",
                                                 "VEC"),data_sub2$TACIT_PCF_v2,"Others")

group_ct_st_2=paste0(Group2,"::",Group1,sep="")

# Reclassify specific categories using levels
group_ct_st_2 <- ifelse(grepl("Others::PD-1-", group_ct_st_2) | grepl("Others::PD-1+", group_ct_st_2), "Others", group_ct_st_2)
group_ct_st_2 <- ifelse(grepl("VEC::PD-1-", group_ct_st_2) | grepl("VEC::PD-1+", group_ct_st_2), "VEC", group_ct_st_2)


# Check the transformed table
print(table(group_ct_st_2))

group_counts <- as.data.frame(table(group_ct_st_2))
colnames(group_counts) <- c("group_ct_st", "count")
# Merge the counts back into the original data
data_plot_ct <- data.frame(X=data_sub2$X,Y=data_sub2$Y,group_ct_st=group_ct_st_2)

data_plot_ct <- merge(data_plot_ct, group_counts, by = "group_ct_st")

# Create a new column with the group names and counts
data_plot_ct$group_with_count <- paste(data_plot_ct$group_ct_st, "(n =", data_plot_ct$count, ")")

table(data_plot_ct$group_with_count)

colors <- c(
  "B Cells::PD-1-  (n = 200 )" = "#1f77b4",  # Blue
  "B Cells::PD-1+  (n = 130 )" = "#ff7f0e",   # Orange
  "CD4+ T Cells::PD-1-  (n = 59 )" = "#2ca02c",  # Green
  "CD4+ T Cells::PD-1+  (n = 49 )" = "#d62728",   # Red
  "Others (n = 312 )" = "grey90",  # Purple
  "VEC (n = 17 )" = "#8c564b"  # Brown
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
  )+ggtitle("PCF")+
  scale_fill_manual(values = colors)














Group1=ifelse(as.numeric(data_sub2$P_PD.1)>0,paste0("PD-1","+",sep=" "),paste0("PD-1","-",sep=" "))
Group2=ifelse(data_sub2$TACIT_Xenium_v2%in%c("CD4+ T Cells","B Cells",
                                          "VEC"),data_sub2$TACIT_Xenium_v2,"Others")

group_ct_st_2=paste0(Group2,"::",Group1,sep="")

# Reclassify specific categories using levels
group_ct_st_2 <- ifelse(grepl("Others::PD-1-", group_ct_st_2) | grepl("Others::PD-1+", group_ct_st_2), "Others", group_ct_st_2)
group_ct_st_2 <- ifelse(grepl("VEC::PD-1-", group_ct_st_2) | grepl("VEC::PD-1+", group_ct_st_2), "VEC", group_ct_st_2)


# Check the transformed table
print(table(group_ct_st_2))

group_counts <- as.data.frame(table(group_ct_st_2))
colnames(group_counts) <- c("group_ct_st", "count")
# Merge the counts back into the original data
data_plot_ct <- data.frame(X=data_sub2$X,Y=data_sub2$Y,group_ct_st=group_ct_st_2)

data_plot_ct <- merge(data_plot_ct, group_counts, by = "group_ct_st")

# Create a new column with the group names and counts
data_plot_ct$group_with_count <- paste(data_plot_ct$group_ct_st, "(n =", data_plot_ct$count, ")")

table(data_plot_ct$group_with_count)

colors <- c(
  "B Cells::PD-1-  (n = 208 )" = "#1f77b4",  # Blue
  "B Cells::PD-1+  (n = 117 )" = "#ff7f0e",   # Orange
  "CD4+ T Cells::PD-1-  (n = 40 )" = "#2ca02c",  # Green
  "CD4+ T Cells::PD-1+  (n = 27 )" = "#d62728",   # Red
  "Others (n = 319 )" = "grey90",  # Purple
  "VEC (n = 56 )" = "#8c564b"  # Brown
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
  )+ggtitle("Xenium")+
  scale_fill_manual(values = colors)







































Group1=ifelse(as.numeric(data_sub3$P_PD.1)>0,paste0("PD-1","+",sep=" "),paste0("PD-1","-",sep=" "))
Group2=ifelse(data_sub3$TACIT_PCF_Xenium_v2%in%c("CD4+ T Cells","B Cells",
                                                 "VEC"),data_sub3$TACIT_PCF_Xenium_v2,"Others")

group_ct_st_2=paste0(Group2,"::",Group1,sep="")

# Reclassify specific categories using levels
group_ct_st_2 <- ifelse(grepl("Others::PD-1-", group_ct_st_2) | grepl("Others::PD-1+", group_ct_st_2), "Others", group_ct_st_2)
group_ct_st_2 <- ifelse(grepl("VEC::PD-1-", group_ct_st_2) | grepl("VEC::PD-1+", group_ct_st_2), "VEC", group_ct_st_2)


# Check the transformed table
print(table(group_ct_st_2))

group_counts <- as.data.frame(table(group_ct_st_2))
colnames(group_counts) <- c("group_ct_st", "count")
# Merge the counts back into the original data
data_plot_ct <- data.frame(X=data_sub3$X,Y=data_sub3$Y,group_ct_st=group_ct_st_2)

data_plot_ct <- merge(data_plot_ct, group_counts, by = "group_ct_st")

# Create a new column with the group names and counts
data_plot_ct$group_with_count <- paste(data_plot_ct$group_ct_st, "(n =", data_plot_ct$count, ")")

table(data_plot_ct$group_with_count)

colors <- c(
  "B Cells::PD-1-  (n = 134 )" = "#1f77b4",  # Blue
  "B Cells::PD-1+  (n = 59 )" = "#ff7f0e",   # Orange
  "CD4+ T Cells::PD-1-  (n = 64 )" = "#2ca02c",  # Green
  "CD4+ T Cells::PD-1+  (n = 38 )" = "#d62728",   # Red
  "Others (n = 1164 )" = "grey90",  # Purple
  "VEC (n = 88 )" = "#8c564b"  # Brown
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
  )+ggtitle("PCF & Xenium")+
  scale_fill_manual(values = colors)












Group1=ifelse(as.numeric(data_sub3$P_PD.1)>0,paste0("PD-1","+",sep=" "),paste0("PD-1","-",sep=" "))
Group2=ifelse(data_sub3$TACIT_PCF_v2%in%c("CD4+ T Cells","B Cells",
                                          "VEC"),data_sub3$TACIT_PCF_v2,"Others")

group_ct_st_2=paste0(Group2,"::",Group1,sep="")

# Reclassify specific categories using levels
group_ct_st_2 <- ifelse(grepl("Others::PD-1-", group_ct_st_2) | grepl("Others::PD-1+", group_ct_st_2), "Others", group_ct_st_2)
group_ct_st_2 <- ifelse(grepl("VEC::PD-1-", group_ct_st_2) | grepl("VEC::PD-1+", group_ct_st_2), "VEC", group_ct_st_2)


# Check the transformed table
print(table(group_ct_st_2))

group_counts <- as.data.frame(table(group_ct_st_2))
colnames(group_counts) <- c("group_ct_st", "count")
# Merge the counts back into the data_sub3 data
data_plot_ct <- data.frame(X=data_sub3$X,Y=data_sub3$Y,group_ct_st=group_ct_st_2)

data_plot_ct <- merge(data_plot_ct, group_counts, by = "group_ct_st")

# Create a new column with the group names and counts
data_plot_ct$group_with_count <- paste(data_plot_ct$group_ct_st, "(n =", data_plot_ct$count, ")")

table(data_plot_ct$group_with_count)

colors <- c(
  "B Cells::PD-1-  (n = 58 )" = "#1f77b4",  # Blue
  "B Cells::PD-1+  (n = 24 )" = "#ff7f0e",   # Orange
  "CD4+ T Cells::PD-1-  (n = 39 )" = "#2ca02c",  # Green
  "CD4+ T Cells::PD-1+  (n = 16 )" = "#d62728",   # Red
  "Others (n = 1364 )" = "grey90",  # Purple
  "VEC (n = 46 )" = "#8c564b"  # Brown
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
  )+ggtitle("PCF")+
  scale_fill_manual(values = colors)














Group1=ifelse(as.numeric(data_sub3$P_PD.1)>0,paste0("PD-1","+",sep=" "),paste0("PD-1","-",sep=" "))
Group2=ifelse(data_sub3$TACIT_Xenium_v2%in%c("CD4+ T Cells","B Cells",
                                             "VEC"),data_sub3$TACIT_Xenium_v2,"Others")

group_ct_st_2=paste0(Group2,"::",Group1,sep="")

# Reclassify specific categories using levels
group_ct_st_2 <- ifelse(grepl("Others::PD-1-", group_ct_st_2) | grepl("Others::PD-1+", group_ct_st_2), "Others", group_ct_st_2)
group_ct_st_2 <- ifelse(grepl("VEC::PD-1-", group_ct_st_2) | grepl("VEC::PD-1+", group_ct_st_2), "VEC", group_ct_st_2)


# Check the transformed table
print(table(group_ct_st_2))

group_counts <- as.data.frame(table(group_ct_st_2))
colnames(group_counts) <- c("group_ct_st", "count")
# Merge the counts back into the original data
data_plot_ct <- data.frame(X=data_sub3$X,Y=data_sub3$Y,group_ct_st=group_ct_st_2)

data_plot_ct <- merge(data_plot_ct, group_counts, by = "group_ct_st")

# Create a new column with the group names and counts
data_plot_ct$group_with_count <- paste(data_plot_ct$group_ct_st, "(n =", data_plot_ct$count, ")")

table(data_plot_ct$group_with_count)

colors <- c(
  "B Cells::PD-1-  (n = 97 )" = "#1f77b4",  # Blue
  "B Cells::PD-1+  (n = 36 )" = "#ff7f0e",   # Orange
  "CD4+ T Cells::PD-1-  (n = 54 )" = "#2ca02c",  # Green
  "CD4+ T Cells::PD-1+  (n = 30 )" = "#d62728",   # Red
  "Others (n = 1225 )" = "grey90",  # Purple
  "VEC (n = 105 )" = "#8c564b"  # Brown
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
  )+ggtitle("Xenium")+
  scale_fill_manual(values = colors)

































Group1=ifelse(as.numeric(data_sub4$P_PD.1)>0,paste0("PD-1","+",sep=" "),paste0("PD-1","-",sep=" "))
Group2=ifelse(data_sub4$TACIT_PCF_Xenium_v2%in%c("CD4+ T Cells","B Cells",
                                                 "VEC"),data_sub4$TACIT_PCF_Xenium_v2,"Others")

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
data_plot_ct$group_with_count <- paste(data_plot_ct$group_ct_st, "(n =", data_plot_ct$count, ")")

table(data_plot_ct$group_with_count)

colors <- c(
  "B Cells::PD-1-  (n = 35 )" = "#1f77b4",  # Blue
  "B Cells::PD-1+  (n = 24 )" = "#ff7f0e",   # Orange
  "CD4+ T Cells::PD-1-  (n = 197 )" = "#2ca02c",  # Green
  "CD4+ T Cells::PD-1+  (n = 390 )" = "#d62728",   # Red
  "Others (n = 335 )" = "grey90",  # Purple
  "VEC (n = 83 )" = "#8c564b"  # Brown
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
  )+ggtitle("PCF & Xenium")+
  scale_fill_manual(values = colors)












Group1=ifelse(as.numeric(data_sub4$P_PD.1)>0,paste0("PD-1","+",sep=" "),paste0("PD-1","-",sep=" "))
Group2=ifelse(data_sub4$TACIT_PCF_v2%in%c("CD4+ T Cells","B Cells",
                                          "VEC"),data_sub4$TACIT_PCF_v2,"Others")

group_ct_st_2=paste0(Group2,"::",Group1,sep="")

# Reclassify specific categories using levels
group_ct_st_2 <- ifelse(grepl("Others::PD-1-", group_ct_st_2) | grepl("Others::PD-1+", group_ct_st_2), "Others", group_ct_st_2)
group_ct_st_2 <- ifelse(grepl("VEC::PD-1-", group_ct_st_2) | grepl("VEC::PD-1+", group_ct_st_2), "VEC", group_ct_st_2)


# Check the transformed table
print(table(group_ct_st_2))

group_counts <- as.data.frame(table(group_ct_st_2))
colnames(group_counts) <- c("group_ct_st", "count")
# Merge the counts back into the data_sub3 data
data_plot_ct <- data.frame(X=data_sub4$X,Y=data_sub4$Y,group_ct_st=group_ct_st_2)

data_plot_ct <- merge(data_plot_ct, group_counts, by = "group_ct_st")

# Create a new column with the group names and counts
data_plot_ct$group_with_count <- paste(data_plot_ct$group_ct_st, "(n =", data_plot_ct$count, ")")

table(data_plot_ct$group_with_count)

colors <- c(
  "B Cells::PD-1-  (n = 24 )" = "#1f77b4",  # Blue
  "B Cells::PD-1+  (n = 20 )" = "#ff7f0e",   # Orange
  "CD4+ T Cells::PD-1-  (n = 162 )" = "#2ca02c",  # Green
  "CD4+ T Cells::PD-1+  (n = 350 )" = "#d62728",   # Red
  "Others (n = 465 )" = "grey90",  # Purple
  "VEC (n = 43 )" = "#8c564b"  # Brown
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
  )+ggtitle("PCF")+
  scale_fill_manual(values = colors)














Group1=ifelse(as.numeric(data_sub4$P_PD.1)>0,paste0("PD-1","+",sep=" "),paste0("PD-1","-",sep=" "))
Group2=ifelse(data_sub4$TACIT_Xenium_v2%in%c("CD4+ T Cells","B Cells",
                                             "VEC"),data_sub4$TACIT_Xenium_v2,"Others")

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
data_plot_ct$group_with_count <- paste(data_plot_ct$group_ct_st, "(n =", data_plot_ct$count, ")")

table(data_plot_ct$group_with_count)

colors <- c(
  "B Cells::PD-1-  (n = 76 )" = "#1f77b4",  # Blue
  "B Cells::PD-1+  (n = 77 )" = "#ff7f0e",   # Orange
  "CD4+ T Cells::PD-1-  (n = 130 )" = "#2ca02c",  # Green
  "CD4+ T Cells::PD-1+  (n = 168 )" = "#d62728",   # Red
  "Others (n = 472 )" = "grey90",  # Purple
  "VEC (n = 141 )" = "#8c564b"  # Brown
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
  )+ggtitle("Xenium")+
  scale_fill_manual(values = colors)


























Group1=ifelse(as.numeric(data_sub1$PDCD1)>0,paste0("PDCD1","+",sep=" "),paste0("PDCD1","-",sep=" "))
Group2=ifelse(data_sub1$TACIT_PCF_v2%in%c("CD4+ T Cells","B Cells",
                                          "VEC"),data_sub1$TACIT_PCF_v2,"Others")

group_ct_st_2=paste0(Group2,"::",Group1,sep="")

# Reclassify specific categories using levels
group_ct_st_2 <- ifelse(grepl("Others::PDCD1-", group_ct_st_2) | grepl("Others::PDCD1+", group_ct_st_2), "Others", group_ct_st_2)
group_ct_st_2 <- ifelse(grepl("VEC::PDCD1-", group_ct_st_2) | grepl("VEC::PDCD1+", group_ct_st_2), "VEC", group_ct_st_2)


# Check the transformed table
print(table(group_ct_st_2))

group_counts <- as.data.frame(table(group_ct_st_2))
colnames(group_counts) <- c("group_ct_st", "count")
# Merge the counts back into the original data
data_plot_ct <- data.frame(X=data_sub1$X,Y=data_sub1$Y,group_ct_st=group_ct_st_2)

data_plot_ct <- merge(data_plot_ct, group_counts, by = "group_ct_st")

# Create a new column with the group names and counts
data_plot_ct$group_with_count <- paste(data_plot_ct$group_ct_st, "(n =", data_plot_ct$count, ")")

table(data_plot_ct$group_with_count)

colors <- c(
  "B Cells::PDCD1-  (n = 657 )" = "#1f77b4",  # Blue
  "B Cells::PDCD1+  (n = 1 )" = "#ff7f0e",   # Orange
  "CD4+ T Cells::PDCD1-  (n = 203 )" = "#2ca02c",  # Green
  "CD4+ T Cells::PDCD1+  (n = 3 )" = "#d62728",   # Red
  "Others (n = 2119 )" = "grey90",  # Purple
  "VEC (n = 61 )" = "#8c564b"  # Brown
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
  )+ggtitle("PCF")+
  scale_fill_manual(values = colors)









Group1=ifelse(as.numeric(data_sub1$PDCD1)>0,paste0("PDCD1","+",sep=" "),paste0("PDCD1","-",sep=" "))
Group2=ifelse(data_sub1$TACIT_Xenium_v2%in%c("CD4+ T Cells","B Cells",
                                             "VEC"),data_sub1$TACIT_Xenium_v2,"Others")

group_ct_st_2=paste0(Group2,"::",Group1,sep="")

# Reclassify specific categories using levels
group_ct_st_2 <- ifelse(grepl("Others::PDCD1-", group_ct_st_2) | grepl("Others::PDCD1+", group_ct_st_2), "Others", group_ct_st_2)
group_ct_st_2 <- ifelse(grepl("VEC::PDCD1-", group_ct_st_2) | grepl("VEC::PDCD1+", group_ct_st_2), "VEC", group_ct_st_2)


# Check the transformed table
print(table(group_ct_st_2))

group_counts <- as.data.frame(table(group_ct_st_2))
colnames(group_counts) <- c("group_ct_st", "count")
# Merge the counts back into the original data
data_plot_ct <- data.frame(X=data_sub1$X,Y=data_sub1$Y,group_ct_st=group_ct_st_2)

data_plot_ct <- merge(data_plot_ct, group_counts, by = "group_ct_st")

# Create a new column with the group names and counts
data_plot_ct$group_with_count <- paste(data_plot_ct$group_ct_st, "(n =", data_plot_ct$count, ")")

table(data_plot_ct$group_with_count)

colors <- c(
  "B Cells::PDCD1-  (n = 645 )" = "#1f77b4",  # Blue
  "B Cells::PDCD1+  (n = 1 )" = "#ff7f0e",   # Orange
  "CD4+ T Cells::PDCD1-  (n = 215 )" = "#2ca02c",  # Green
  "CD4+ T Cells::PDCD1+  (n = 4 )" = "#d62728",   # Red
  "Others (n = 2011 )" = "grey90",  # Purple
  "VEC (n = 168 )" = "#8c564b"  # Brown
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
  )+ggtitle("Xenium")+
  scale_fill_manual(values = colors)








Group1=ifelse(as.numeric(data_sub1$PDCD1)>0,paste0("PDCD1","+",sep=" "),paste0("PDCD1","-",sep=" "))
Group2=ifelse(data_sub1$TACIT_PCF_Xenium_v2%in%c("CD4+ T Cells","B Cells",
                                                 "VEC"),data_sub1$TACIT_PCF_Xenium_v2,"Others")

group_ct_st_2=paste0(Group2,"::",Group1,sep="")

# Reclassify specific categories using levels
group_ct_st_2 <- ifelse(grepl("Others::PDCD1-", group_ct_st_2) | grepl("Others::PDCD1+", group_ct_st_2), "Others", group_ct_st_2)
group_ct_st_2 <- ifelse(grepl("VEC::PDCD1-", group_ct_st_2) | grepl("VEC::PDCD1+", group_ct_st_2), "VEC", group_ct_st_2)


# Check the transformed table
print(table(group_ct_st_2))

group_counts <- as.data.frame(table(group_ct_st_2))
colnames(group_counts) <- c("group_ct_st", "count")
# Merge the counts back into the original data
data_plot_ct <- data.frame(X=data_sub1$X,Y=data_sub1$Y,group_ct_st=group_ct_st_2)

data_plot_ct <- merge(data_plot_ct, group_counts, by = "group_ct_st")

# Create a new column with the group names and counts
data_plot_ct$group_with_count <- paste(data_plot_ct$group_ct_st, "(n =", data_plot_ct$count, ")")

table(data_plot_ct$group_with_count)

colors <- c(
  "B Cells::PDCD1-  (n = 852 )" = "#1f77b4",  # Blue
  "B Cells::PDCD1+  (n = 98 )" = "#ff7f0e",   # Orange
  "CD4+ T Cells::PDCD1-  (n = 294 )" = "#2ca02c",  # Green
  "CD4+ T Cells::PDCD1+  (n = 6 )" = "#d62728",   # Red
  "Others (n = 1786 )" = "grey90",  # Purple
  "VEC (n = 106 )" = "#8c564b"  # Brown
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
  )+ggtitle("PCF & Xenium")+
  scale_fill_manual(values = colors)









Group1=ifelse(as.numeric(data_sub2$PDCD1)>0,paste0("PDCD1","+",sep=" "),paste0("PDCD1","-",sep=" "))
Group2=ifelse(data_sub2$TACIT_PCF_Xenium_v2%in%c("CD4+ T Cells","B Cells",
                                                 "VEC"),data_sub2$TACIT_PCF_Xenium_v2,"Others")

group_ct_st_2=paste0(Group2,"::",Group1,sep="")

# Reclassify specific categories using levels
group_ct_st_2 <- ifelse(grepl("Others::PDCD1-", group_ct_st_2) | grepl("Others::PDCD1+", group_ct_st_2), "Others", group_ct_st_2)
group_ct_st_2 <- ifelse(grepl("VEC::PDCD1-", group_ct_st_2) | grepl("VEC::PDCD1+", group_ct_st_2), "VEC", group_ct_st_2)


# Check the transformed table
print(table(group_ct_st_2))

group_counts <- as.data.frame(table(group_ct_st_2))
colnames(group_counts) <- c("group_ct_st", "count")
# Merge the counts back into the original data
data_plot_ct <- data.frame(X=data_sub2$X,Y=data_sub2$Y,group_ct_st=group_ct_st_2)

data_plot_ct <- merge(data_plot_ct, group_counts, by = "group_ct_st")

# Create a new column with the group names and counts
data_plot_ct$group_with_count <- paste(data_plot_ct$group_ct_st, "(n =", data_plot_ct$count, ")")

table(data_plot_ct$group_with_count)

colors <- c(
  "B Cells::PDCD1-  (n = 343 )" = "#1f77b4",  # Blue
  "B Cells::PDCD1+  (n = 3 )" = "#ff7f0e",   # Orange
  "CD4+ T Cells::PDCD1-  (n = 140 )" = "#2ca02c",  # Green
  "CD4+ T Cells::PDCD1+  (n = 1 )" = "#d62728",   # Red
  "Others (n = 260 )" = "grey90",  # Purple
  "VEC (n = 20 )" = "#8c564b"  # Brown
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
  )+ggtitle("PCF & Xenium")+
  scale_fill_manual(values = colors)












Group1=ifelse(as.numeric(data_sub2$PDCD1)>0,paste0("PDCD1","+",sep=" "),paste0("PDCD1","-",sep=" "))
Group2=ifelse(data_sub2$TACIT_PCF_v2%in%c("CD4+ T Cells","B Cells",
                                          "VEC"),data_sub2$TACIT_PCF_v2,"Others")

group_ct_st_2=paste0(Group2,"::",Group1,sep="")

# Reclassify specific categories using levels
group_ct_st_2 <- ifelse(grepl("Others::PDCD1-", group_ct_st_2) | grepl("Others::PDCD1+", group_ct_st_2), "Others", group_ct_st_2)
group_ct_st_2 <- ifelse(grepl("VEC::PDCD1-", group_ct_st_2) | grepl("VEC::PDCD1+", group_ct_st_2), "VEC", group_ct_st_2)


# Check the transformed table
print(table(group_ct_st_2))

group_counts <- as.data.frame(table(group_ct_st_2))
colnames(group_counts) <- c("group_ct_st", "count")
# Merge the counts back into the original data
data_plot_ct <- data.frame(X=data_sub2$X,Y=data_sub2$Y,group_ct_st=group_ct_st_2)

data_plot_ct <- merge(data_plot_ct, group_counts, by = "group_ct_st")

# Create a new column with the group names and counts
data_plot_ct$group_with_count <- paste(data_plot_ct$group_ct_st, "(n =", data_plot_ct$count, ")")

table(data_plot_ct$group_with_count)

colors <- c(
  "B Cells::PDCD1-  (n = 328 )" = "#1f77b4",  # Blue
  "B Cells::PDCD1+  (n = 2 )" = "#ff7f0e",   # Orange
  "CD4+ T Cells::PDCD1-  (n = 107 )" = "#2ca02c",  # Green
  "CD4+ T Cells::PDCD1+  (n = 1 )" = "#d62728",   # Red
  "Others (n = 312 )" = "grey90",  # Purple
  "VEC (n = 17 )" = "#8c564b"  # Brown
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
  )+ggtitle("PCF")+
  scale_fill_manual(values = colors)














Group1=ifelse(as.numeric(data_sub2$P_PD.1)>0,paste0("PDCD1","+",sep=" "),paste0("PDCD1","-",sep=" "))
Group2=ifelse(data_sub2$TACIT_Xenium_v2%in%c("CD4+ T Cells","B Cells",
                                             "VEC"),data_sub2$TACIT_Xenium_v2,"Others")

group_ct_st_2=paste0(Group2,"::",Group1,sep="")

# Reclassify specific categories using levels
group_ct_st_2 <- ifelse(grepl("Others::PDCD1-", group_ct_st_2) | grepl("Others::PDCD1+", group_ct_st_2), "Others", group_ct_st_2)
group_ct_st_2 <- ifelse(grepl("VEC::PDCD1-", group_ct_st_2) | grepl("VEC::PDCD1+", group_ct_st_2), "VEC", group_ct_st_2)


# Check the transformed table
print(table(group_ct_st_2))

group_counts <- as.data.frame(table(group_ct_st_2))
colnames(group_counts) <- c("group_ct_st", "count")
# Merge the counts back into the original data
data_plot_ct <- data.frame(X=data_sub2$X,Y=data_sub2$Y,group_ct_st=group_ct_st_2)

data_plot_ct <- merge(data_plot_ct, group_counts, by = "group_ct_st")

# Create a new column with the group names and counts
data_plot_ct$group_with_count <- paste(data_plot_ct$group_ct_st, "(n =", data_plot_ct$count, ")")

table(data_plot_ct$group_with_count)

colors <- c(
  "B Cells::PDCD1-  (n = 328 )" = "#1f77b4",  # Blue
  "B Cells::PDCD1+  (n = 2 )" = "#ff7f0e",   # Orange
  "CD4+ T Cells::PDCD1-  (n = 107 )" = "#2ca02c",  # Green
  "CD4+ T Cells::PDCD1+  (n = 1 )" = "#d62728",   # Red
  "Others (n = 312 )" = "grey90",  # Purple
  "VEC (n = 17 )" = "#8c564b"  # Brown
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
  )+ggtitle("Xenium")+
  scale_fill_manual(values = colors)







































Group1=ifelse(as.numeric(data_sub4$PDCD1)>0,paste0("PDCD1","+",sep=" "),paste0("PDCD1","-",sep=" "))
Group2=ifelse(data_sub4$TACIT_Xenium_v2%in%c("CD4+ T Cells"),data_sub4$TACIT_Xenium_v2,"Others")

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



















Group1=ifelse(as.numeric(data_sub3$PDCD1)>0,paste0("PDCD1","+",sep=" "),paste0("PDCD1","-",sep=" "))
Group2=ifelse(data_sub3$TACIT_PCF_v2%in%c("CD4+ T Cells","B Cells",
                                          "VEC"),data_sub3$TACIT_PCF_v2,"Others")

group_ct_st_2=paste0(Group2,"::",Group1,sep="")

# Reclassify specific categories using levels
group_ct_st_2 <- ifelse(grepl("Others::PDCD1-", group_ct_st_2) | grepl("Others::PDCD1+", group_ct_st_2), "Others", group_ct_st_2)
group_ct_st_2 <- ifelse(grepl("VEC::PDCD1-", group_ct_st_2) | grepl("VEC::PDCD1+", group_ct_st_2), "VEC", group_ct_st_2)


# Check the transformed table
print(table(group_ct_st_2))

group_counts <- as.data.frame(table(group_ct_st_2))
colnames(group_counts) <- c("group_ct_st", "count")
# Merge the counts back into the data_sub3 data
data_plot_ct <- data.frame(X=data_sub3$X,Y=data_sub3$Y,group_ct_st=group_ct_st_2)

data_plot_ct <- merge(data_plot_ct, group_counts, by = "group_ct_st")

# Create a new column with the group names and counts
data_plot_ct$group_with_count <- paste(data_plot_ct$group_ct_st, "(n =", data_plot_ct$count, ")")

table(data_plot_ct$group_with_count)

colors <- c(
  "B Cells::PDCD1-  (n = 82 )" = "#1f77b4",  # Blue
  "B Cells::PDCD1+  (n = 24 )" = "#ff7f0e",   # Orange
  "CD4+ T Cells::PDCD1-  (n = 54 )" = "#2ca02c",  # Green
  "CD4+ T Cells::PDCD1+  (n = 1 )" = "#d62728",   # Red
  "Others (n = 1364 )" = "grey90",  # Purple
  "VEC (n = 46 )" = "#8c564b"  # Brown
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
  )+ggtitle("PCF")+
  scale_fill_manual(values = colors)














Group1=ifelse(as.numeric(data_sub3$PDCD1)>0,paste0("PDCD1","+",sep=" "),paste0("PDCD1","-",sep=" "))
Group2=ifelse(data_sub3$TACIT_Xenium_v2%in%c("CD4+ T Cells","B Cells",
                                             "VEC"),data_sub3$TACIT_Xenium_v2,"Others")

group_ct_st_2=paste0(Group2,"::",Group1,sep="")

# Reclassify specific categories using levels
group_ct_st_2 <- ifelse(grepl("Others::PDCD1-", group_ct_st_2) | grepl("Others::PDCD1+", group_ct_st_2), "Others", group_ct_st_2)
group_ct_st_2 <- ifelse(grepl("VEC::PDCD1-", group_ct_st_2) | grepl("VEC::PDCD1+", group_ct_st_2), "VEC", group_ct_st_2)


# Check the transformed table
print(table(group_ct_st_2))

group_counts <- as.data.frame(table(group_ct_st_2))
colnames(group_counts) <- c("group_ct_st", "count")
# Merge the counts back into the original data
data_plot_ct <- data.frame(X=data_sub3$X,Y=data_sub3$Y,group_ct_st=group_ct_st_2)

data_plot_ct <- merge(data_plot_ct, group_counts, by = "group_ct_st")

# Create a new column with the group names and counts
data_plot_ct$group_with_count <- paste(data_plot_ct$group_ct_st, "(n =", data_plot_ct$count, ")")

table(data_plot_ct$group_with_count)

colors <- c(
  "B Cells::PDCD1-  (n = 133 )" = "#1f77b4",  # Blue
  "B Cells::PDCD1+  (n = 36 )" = "#ff7f0e",   # Orange
  "CD4+ T Cells::PDCD1-  (n = 82 )" = "#2ca02c",  # Green
  "CD4+ T Cells::PDCD1+  (n = 2 )" = "#d62728",   # Red
  "Others (n = 1225 )" = "grey90",  # Purple
  "VEC (n = 105 )" = "#8c564b"  # Brown
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
  )+ggtitle("Xenium")+
  scale_fill_manual(values = colors)

































Group1=ifelse(as.numeric(data_sub4$PDCD1)>0,paste0("PDCD1","+",sep=" "),paste0("PDCD1","-",sep=" "))
Group2=ifelse(data_sub4$TACIT_PCF_Xenium_v2%in%c("CD4+ T Cells","B Cells",
                                                 "VEC"),data_sub4$TACIT_PCF_Xenium_v2,"Others")

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
data_plot_ct$group_with_count <- paste(data_plot_ct$group_ct_st, "(n =", data_plot_ct$count, ")")

table(data_plot_ct$group_with_count)

colors <- c(
  "B Cells::PDCD1-  (n = 59 )" = "#1f77b4",  # Blue
  "B Cells::PDCD1+  (n = 24 )" = "#ff7f0e",   # Orange
  "CD4+ T Cells::PDCD1-  (n = 579 )" = "#2ca02c",  # Green
  "CD4+ T Cells::PDCD1+  (n = 8 )" = "#d62728",   # Red
  "Others (n = 335 )" = "grey90",  # Purple
  "VEC (n = 83 )" = "#8c564b"  # Brown
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
  )+ggtitle("PCF & Xenium")+
  scale_fill_manual(values = colors)












Group1=ifelse(as.numeric(data_sub4$PDCD1)>0,paste0("PDCD1","+",sep=" "),paste0("PDCD1","-",sep=" "))
Group2=ifelse(data_sub4$TACIT_PCF_v2%in%c("CD4+ T Cells","B Cells",
                                          "VEC"),data_sub4$TACIT_PCF_v2,"Others")

group_ct_st_2=paste0(Group2,"::",Group1,sep="")

# Reclassify specific categories using levels
group_ct_st_2 <- ifelse(grepl("Others::PDCD1-", group_ct_st_2) | grepl("Others::PDCD1+", group_ct_st_2), "Others", group_ct_st_2)
group_ct_st_2 <- ifelse(grepl("VEC::PDCD1-", group_ct_st_2) | grepl("VEC::PDCD1+", group_ct_st_2), "VEC", group_ct_st_2)


# Check the transformed table
print(table(group_ct_st_2))

group_counts <- as.data.frame(table(group_ct_st_2))
colnames(group_counts) <- c("group_ct_st", "count")
# Merge the counts back into the data_sub3 data
data_plot_ct <- data.frame(X=data_sub4$X,Y=data_sub4$Y,group_ct_st=group_ct_st_2)

data_plot_ct <- merge(data_plot_ct, group_counts, by = "group_ct_st")

# Create a new column with the group names and counts
data_plot_ct$group_with_count <- paste(data_plot_ct$group_ct_st, "(n =", data_plot_ct$count, ")")

table(data_plot_ct$group_with_count)

colors <- c(
  "B Cells::PDCD1-  (n = 44 )" = "#1f77b4",  # Blue
  "B Cells::PDCD1+  (n = 0 )" = "#ff7f0e",   # Orange
  "CD4+ T Cells::PDCD1-  (n = 505 )" = "#2ca02c",  # Green
  "CD4+ T Cells::PDCD1+  (n = 7 )" = "#d62728",   # Red
  "Others (n = 465 )" = "grey90",  # Purple
  "VEC (n = 43 )" = "#8c564b"  # Brown
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
  )+ggtitle("PCF")+
  scale_fill_manual(values = colors)














Group1=ifelse(as.numeric(data_sub4$PDCD1)>0,paste0("PDCD1","+",sep=" "),paste0("PDCD1","-",sep=" "))
Group2=ifelse(data_sub4$TACIT_Xenium_v2%in%c("CD4+ T Cells","B Cells",
                                             "VEC"),data_sub4$TACIT_Xenium_v2,"Others")

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
data_plot_ct$group_with_count <- paste(data_plot_ct$group_ct_st, "(n =", data_plot_ct$count, ")")

table(data_plot_ct$group_with_count)

colors <- c(
  "B Cells::PDCD1-  (n = 150 )" = "#1f77b4",  # Blue
  "B Cells::PDCD1+  (n = 3 )" = "#ff7f0e",   # Orange
  "CD4+ T Cells::PDCD1-  (n = 295 )" = "#2ca02c",  # Green
  "CD4+ T Cells::PDCD1+  (n = 3 )" = "#d62728",   # Red
  "Others (n = 472 )" = "grey90",  # Purple
  "VEC (n = 141 )" = "#8c564b"  # Brown
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
  )+ggtitle("Xenium")+
  scale_fill_manual(values = colors)




data_sub_final=rbind(data.frame(data_sub1,TLS="TLS1"),
                     data.frame(data_sub2,TLS="TLS2"),
                     data.frame(data_sub3,TLS="TLS3"),
                     data.frame(data_sub4,TLS="TLS4"))
data_sub_final=data_sub_final#[which(data_sub_final$TACIT_PCF_v2%in%c("CD4+ T Cells","B Cells",
                            #        "VEC")),]


data_sub_final_plot=data.frame(Xenium=data_sub_final$TACIT_Xenium_v2,PCF=data_sub_final$TACIT_PCF_v2,TLS=data_sub_final$TLS)
data_sub_final_plot=data_sub_final_plot[which(data_sub_final_plot$Xenium!="Others"),]
data_sub_final_plot=data_sub_final_plot[which(data_sub_final_plot$PCF!="Others"),]



library(dplyr)
library(tidyr)
library(ggplot2)

# Assume data_sub_final_plot is already loaded
# Calculate proportions
proportions <- data_sub_final_plot %>%
  pivot_longer(cols = c("Xenium", "PCF"), names_to = "Method", values_to = "Cell_Type") %>%
  count(TLS, Method, Cell_Type) %>%
  group_by(TLS, Method) %>%
  mutate(Total = sum(n)) %>%
  ungroup() %>%
  mutate(Proportion = n / Total) %>%
  select(TLS, Method, Cell_Type, Proportion)

# Reshape for plotting
proportions <- proportions %>%
  pivot_wider(names_from = Method, values_from = Proportion, names_prefix = "Proportion_")

proportions_PCF=proportions[which(proportions$Cell_Type=="Ductal Cells"),]

cor(proportions$Proportion_PCF,proportions$Proportion_Xenium,method = "spearman")

ggplot(proportions, aes(x = Proportion_Xenium, y = Proportion_PCF)) +
  geom_point(color = "black", size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "red")   +  # Display correlation
  theme_classic(base_size = 20) +
  labs(
    title = "",
    x = "Proportion by Xenium",
    y = "Proportion by PCF"
  )+  # Add a linear model fit line
  geom_text(aes(label = sprintf("R = %.2f", 0.99), x = Inf, y = Inf),
            hjust = 1.1, vjust = 1.1, color = "black", check_overlap = TRUE, size = 5) +ylim(0,0.22)+xlim(0,0.22)









data_sub_final_plot=data.frame(Xenium=data_sub_final$TACIT_1PCF_2Xenium,PCF=data_sub_final$TACIT_2PCF_1Xenium,TLS=data_sub_final$TLS)
data_sub_final_plot=data_sub_final_plot[which(data_sub_final_plot$Xenium!="Others"),]
data_sub_final_plot=data_sub_final_plot[which(data_sub_final_plot$PCF!="Others"),]
data_sub_final_plot$Xenium=ifelse(data_sub_final_plot$Xenium%in%c("Ductal Cells","Ductal Progenitors"),"Ductal Cells",data_sub_final_plot$Xenium)
data_sub_final_plot$PCF=ifelse(data_sub_final_plot$PCF%in%c("Ductal Cells","Ductal Progenitors"),"Ductal Cells",data_sub_final_plot$PCF)
data_sub_final_plot$Xenium=ifelse(data_sub_final_plot$Xenium%in%c("Mucous Acinar Cells","Seromucous Acinar Cells"),"Acinar Cells",data_sub_final_plot$Xenium)
data_sub_final_plot$PCF=ifelse(data_sub_final_plot$PCF%in%c("Mucous Acinar Cells","Seromucous Acinar Cells"),"Acinar Cells",data_sub_final_plot$PCF)


library(dplyr)
library(tidyr)
library(ggplot2)

# Assume data_sub_final_plot is already loaded
# Calculate proportions
proportions <- data_sub_final_plot %>%
  pivot_longer(cols = c("Xenium", "PCF"), names_to = "Method", values_to = "Cell_Type") %>%
  count(TLS, Method, Cell_Type) %>%
  group_by(TLS, Method) %>%
  mutate(Total = sum(n)) %>%
  ungroup() %>%
  mutate(Proportion = n / Total) %>%
  select(TLS, Method, Cell_Type, Proportion)

# Reshape for plotting
proportions <- proportions %>%
  pivot_wider(names_from = Method, values_from = Proportion, names_prefix = "Proportion_")

proportions_All=proportions[which(proportions$Cell_Type=="Venules"),]

cor(proportions_All$Proportion_PCF,proportions_All$Proportion_Xenium,method = "pearson")

ggplot(proportions_All, aes(x = Proportion_Xenium, y = Proportion_PCF)) +
  geom_point(color = "blue", size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "red")   +  # Display correlation
  theme_classic(base_size = 20) +
  labs(
    title = "",
    x = "Proportion by Xenium",
    y = "Proportion by PCF"
  )+  # Add a linear model fit line
  geom_text(aes(label = sprintf("R = %.2f", 0.99), x = Inf, y = Inf),
            hjust = 1.1, vjust = 1.1, color = "black", check_overlap = TRUE, size = 5) +ylim(0,0.22)+xlim(0,0.22)





data_sub_final_plot=data.frame(Xenium=data_sub_final$TACIT_Xenium,PCF=data_sub_final$TACIT_PCF,TLS=data_sub_final$TLS)
data_sub_final_plot=data_sub_final_plot[which(data_sub_final_plot$Xenium%in%c("Venules","Arterioles")|
                                                data_sub_final_plot$PCF%in%c("VEC")),]
data_sub_final_plot=data_sub_final_plot[which(data_sub_final_plot$Xenium!="Others"),]
data_sub_final_plot=data_sub_final_plot[which(data_sub_final_plot$PCF!="Others"),]

data_sub_final_plot$Xenium=ifelse(data_sub_final_plot$Xenium%in%c("Venules","Arterioles"),"VEC",data_sub_final_plot$Xenium)


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
data_sub_final_plot=data_sub_final_plot[which(data_sub_final_plot$Xenium!="Others"),]
data_sub_final_plot=data_sub_final_plot[which(data_sub_final_plot$PCF!="Others"),]


data_sub_final_plot$Xenium=ifelse(data_sub_final_plot$Xenium%in%c("Mucous Acinar Cells","Seromucous Acinar Cells"),"Acinar Cells",data_sub_final_plot$Xenium)
data_sub_final_plot$PCF=ifelse(data_sub_final_plot$PCF%in%c("Mucous Acinar Cells","Seromucous Acinar Cells"),"Acinar Cells",data_sub_final_plot$PCF)


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
        plot.title = element_text(hjust = 0.5),legend.position = "none")+ylim(0,1)  














data_sub_final$TACIT_Xenium_UPDATE=ifelse(data_sub_final$TACIT_Xenium%in%c("Arterioles","Venules"),"VEC",data_sub_final$TACIT_Xenium)
data_sub_final$TACIT_Xenium_UPDATE=ifelse(data_sub_final$TACIT_Xenium_UPDATE%in%c("Ductal Cells","Ductal Progenitors"),"Ductal Cells",data_sub_final$TACIT_Xenium_UPDATE)
data_sub_final$TACIT_Xenium_UPDATE=ifelse(data_sub_final$TACIT_Xenium_UPDATE%in%c("Mucous Acinar Cells","Seromucous Acinar Cells"),"Acinar Cells",data_sub_final$TACIT_Xenium_UPDATE)
data_sub_final$TACIT_Xenium_UPDATE=ifelse(data_sub_final$TACIT_Xenium_UPDATE%in%c("CD8+ Exhausted T Cells"),"CD8+ T Cells",data_sub_final$TACIT_Xenium_UPDATE)
data_sub_final$TACIT_Xenium_UPDATE=ifelse(data_sub_final$TACIT_Xenium_UPDATE%in%c("Myoepithelium","Intermediate Epithelium"),"Myoepithelial",data_sub_final$TACIT_Xenium_UPDATE)
data_sub_final$TACIT_Xenium_UPDATE=ifelse(data_sub_final$TACIT_Xenium_UPDATE%in%unique(data_sub_final$TACIT_PCF),data_sub_final$TACIT_Xenium_UPDATE,"Others")


data_sub_final_plot=data.frame(Xenium=data_sub_final$TACIT_1PCF_2Xenium,PCF=data_sub_final$TACIT_2PCF_1Xenium,TLS=data_sub_final$TLS)
data_sub_final_plot=data_sub_final_plot[which(data_sub_final_plot$Xenium!="Others"),]
data_sub_final_plot=data_sub_final_plot[which(data_sub_final_plot$PCF!="Others"),]
data_sub_final_plot$Xenium=ifelse(data_sub_final_plot$Xenium%in%c("Ductal Cells","Ductal Progenitors"),"Ductal Cells",data_sub_final_plot$Xenium)
data_sub_final_plot$PCF=ifelse(data_sub_final_plot$PCF%in%c("Ductal Cells","Ductal Progenitors"),"Ductal Cells",data_sub_final_plot$PCF)
data_sub_final_plot$Xenium=ifelse(data_sub_final_plot$Xenium%in%c("Mucous Acinar Cells","Seromucous Acinar Cells"),"Acinar Cells",data_sub_final_plot$Xenium)
data_sub_final_plot$PCF=ifelse(data_sub_final_plot$PCF%in%c("Mucous Acinar Cells","Seromucous Acinar Cells"),"Acinar Cells",data_sub_final_plot$PCF)





library(dplyr)
library(tidyr)
library(ggplot2)

# Assume data_sub_final_plot is already loaded
# Calculate proportions
proportions <- data_sub_final_plot %>%
  pivot_longer(cols = c("Xenium", "PCF"), names_to = "Method", values_to = "Cell_Type") %>%
  count(TLS, Method, Cell_Type) %>%
  group_by(TLS, Method) %>%
  mutate(Total = sum(n)) %>%
  ungroup() %>%
  mutate(Proportion = n / Total) %>%
  select(TLS, Method, Cell_Type, Proportion)

# Reshape for plotting
proportions <- proportions %>%
  pivot_wider(names_from = Method, values_from = Proportion, names_prefix = "Proportion_")

proportions_All=proportions[which(proportions$Cell_Type=="Venules"),]

cor(proportions_All$Proportion_PCF,proportions_All$Proportion_Xenium,method = "pearson")

ggplot(proportions_All, aes(x = Proportion_Xenium, y = Proportion_PCF)) +
  geom_point(color = "blue", size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "red")   +  # Display correlation
  theme_classic(base_size = 20) +
  labs(
    title = "",
    x = "Proportion by Xenium",
    y = "Proportion by PCF"
  )+  # Add a linear model fit line
  geom_text(aes(label = sprintf("R = %.2f", 0.99), x = Inf, y = Inf),
            hjust = 1.1, vjust = 1.1, color = "black", check_overlap = TRUE, size = 5) +ylim(0,0.22)+xlim(0,0.22)






















library(dplyr)

proportions <- data_sub_final_plot %>%
  group_by(TLS, ct) %>%
  summarise(
    total = n(),
    group_1 = sum(Group == 1),
    proportion = group_1 / total
  ) %>%
  ungroup()

print(proportions)







data_sub_final_plot2=data.frame(Group=ifelse(data_sub_final$PDCD1>0,1,0),ct=data_sub_final$TACIT_PCF_v2,TLS=data_sub_final$TLS)


library(dplyr)

proportions2 <- data_sub_final_plot2 %>%
  group_by(TLS, ct) %>%
  summarise(
    total = n(),
    group_1 = sum(Group ==1),
    proportion = group_1 / total
  ) %>%
  ungroup()

print(proportions2)



proportions_final=rbind(data.frame(proportions,marker="PD-1"),
                        data.frame(proportions2,marker="PDCD1"))


ggplot(proportions_final, aes(x = ct, y = proportion, fill = marker)) +
  geom_boxplot() +
  labs(title = "Proportion of Cell Types by Marker Positivity (PD-1 and PDCD1) in TLS",
       x = "",
       y = "Proportion") +
  theme_classic(base_size = 20)





data_sub_final=rbind(data.frame(data_sub1,TLS="TLS1"),
                     data.frame(data_sub2,TLS="TLS2"),
                     data.frame(data_sub3,TLS="TLS3"),
                     data.frame(data_sub4,TLS="TLS4"))

data_sub_final=data_sub_final


data_sub_final_plot=data.frame(Xenium=data_sub_final$TACIT_Xenium_v2,PCF=data_sub_final$TACIT_PCF_v2,TLS=data_sub_final$TLS)

confusionMatrix(factor(data_sub_final_plot$Xenium,levels = unique(data_sub_final_plot$Xenium)[-7]),factor(data_sub_final_plot$PCF,levels = unique(data_sub_final_plot$Xenium)[-7]))



# Load necessary libraries
library(dplyr)
library(ggplot2)

# Calculate proportions for each Xenium and PCF within each TLS group
proportions_xenium <- data_sub_final_plot %>%
  group_by(TLS, Xenium) %>%
  summarise(count_xenium = n()) %>%
  ungroup() %>%
  group_by(TLS) %>%
  mutate(proportion_xenium = count_xenium / sum(count_xenium))

proportions_xenium$Group=paste0(proportions_xenium$TLS,"::",proportions_xenium$Xenium,sep="")

proportions_xenium=proportions_xenium[which(proportions_xenium$Xenium%!="Others")),]


proportions_pcf <- data_sub_final_plot %>%
  group_by(TLS, PCF) %>%
  summarise(count_pcf = n()) %>%
  ungroup() %>%
  group_by(TLS) %>%
  mutate(proportion_pcf = count_pcf / sum(count_pcf))

proportions_pcf$Group=paste0(proportions_pcf$TLS,"::",proportions_pcf$PCF,sep="")

proportions_pcf=proportions_pcf[which(proportions_pcf$PCF!="Others"),]

# Merge proportions by TLS
merged_proportions <- merge(proportions_xenium, proportions_pcf, by = "Group")



# Draw the correlation plot
ggplot(merged_proportions, aes(x = proportion_xenium, y = proportion_pcf, color = Xenium)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Correlation Plot of Proportions",
       x = "Proportion Xenium",
       y = "Proportion PCF") +
  theme_classic(base_size = 20)



ggplot(merged_proportions, aes(x = proportion_xenium, y = proportion_pcf, color = Xenium)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  stat_cor(method = "pearson", label.x = 0.3, label.y = 0.3) + # Display correlation value on the plot
  labs(title = "Correlation Plot of Proportions",
       x = "Proportion Xenium",
       y = "Proportion PCF") +
  theme_minimal()


merged_proportions_sub=merged_proportions[which(merged_proportions$Xenium==unique(merged_proportions$Xenium)[8]),]
ggplot(merged_proportions_sub, aes(x = proportion_xenium, y = proportion_pcf, color = Xenium)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") + # Adjust label position
  labs(title = "",
       x = "Proportion Xenium",
       y = "Proportion PCF") +
  theme_classic(base_size = 20) +
  theme(legend.position = "bottom")+ylim(0,max(merged_proportions_sub$proportion_xenium,merged_proportions_sub$proportion_pcf))+xlim(0,max(merged_proportions_sub$proportion_xenium,merged_proportions_sub$proportion_pcf))












Group1=NULL
for (i in 1:nrow(data_sub1)) {
  if(data_sub1$TACIT_PCF_v2[i]%in%c("B Cells","CD4+ T Cells")){
    if(data_sub1$PDCD1[i]>0){
      Group1[i]=1
    }else{
      Group1[i]=0
    }
  }else{
    Group1[i]="Others"
  }
}


Group2=NULL
for (i in 1:nrow(data_sub1)) {
  if(data_sub1$TACIT_PCF_v2[i]%in%c("B Cells","CD4+ T Cells")){
    if(data_sub1$P_PD.1[i]>0){
      Group2[i]=1
    }else{
      Group2[i]=0
    }
  }else{
    Group2[i]="Others"
  }
}


group_ct_st_2=NULL
for (i in 1:length(Group2)) {
  if(Group2[i]!="Others"){
    if(Group1[i]==Group2[i]){
      group_ct_st_2[i]="Agreement"
    }else{
      group_ct_st_2[i]="Disagreement"
    }
  }else{
    group_ct_st_2[i]="Others"
  }
}

group_counts <- as.data.frame(table(group_ct_st_2))
colnames(group_counts) <- c("group_ct_st", "count")
# Merge the counts back into the original data
data_plot_ct <- data.frame(X=data_sub1$X,Y=data_sub1$Y,group_ct_st=group_ct_st_2)

data_plot_ct <- merge(data_plot_ct, group_counts, by = "group_ct_st")

# Create a new column with the group names and counts
data_plot_ct$group_with_count <- paste(data_plot_ct$group_ct_st, "(n =", data_plot_ct$count, ")")

table(data_plot_ct$group_with_count)

colors <- c(
  "Agreement (n = 745 )" = "blue",  # Blue
  "B Cells::PD-1+  (n = 44 )" = "#ff7f0e",   # Orange
  "CD4+ T Cells::PD-1-  (n = 133 )" = "#2ca02c",  # Green
  "Disagreement (n = 119 )" = "red",   # Red
  "Others (n = 2180 )" = "grey90",  # Purple
  "VEC (n = 61 )" = "#8c564b"  # Brown
)

ggplot(data_plot_ct, aes(X, Y, group = -1L)) + 
  geom_voronoi_tile(aes(fill = group_with_count), colour = 'black', max.radius = 40) +
  theme_classic(base_size = 15) +
  labs(x = "X", y = "Y", title = "") + guides(fill = guide_legend(override.aes = list(size = 3),title = "Group"))+
  theme(
    axis.title = element_blank(),         # Remove axis titles
    axis.text = element_blank(),          # Remove axis text
    axis.ticks = element_blank(),         # Remove axis ticks
    axis.line = element_blank(),          # Remove axis lines
    panel.background = element_blank(),   # Remove panel background
    panel.grid = element_blank()          # Remove grid lines
  )+ggtitle("")+
  scale_fill_manual(values = colors)










Group1=NULL
for (i in 1:nrow(data_sub2)) {
  if(data_sub2$TACIT_PCF_v2[i]%in%c("B Cells","CD4+ T Cells")){
    if(data_sub2$PDCD1[i]>0){
      Group1[i]=1
    }else{
      Group1[i]=0
    }
  }else{
    Group1[i]="Others"
  }
}


Group2=NULL
for (i in 1:nrow(data_sub2)) {
  if(data_sub2$TACIT_PCF_v2[i]%in%c("B Cells","CD4+ T Cells")){
    if(data_sub2$P_PD.1[i]>0){
      Group2[i]=1
    }else{
      Group2[i]=0
    }
  }else{
    Group2[i]="Others"
  }
}


group_ct_st_2=NULL
for (i in 1:length(Group2)) {
  if(Group2[i]!="Others"){
    if(Group1[i]==Group2[i]){
      group_ct_st_2[i]="Agreement"
    }else{
      group_ct_st_2[i]="Disagreement"
    }
  }else{
    group_ct_st_2[i]="Others"
  }
}

group_counts <- as.data.frame(table(group_ct_st_2))
colnames(group_counts) <- c("group_ct_st", "count")
# Merge the counts back into the original data
data_plot_ct <- data.frame(X=data_sub2$X,Y=data_sub2$Y,group_ct_st=group_ct_st_2)

data_plot_ct <- merge(data_plot_ct, group_counts, by = "group_ct_st")

# Create a new column with the group names and counts
data_plot_ct$group_with_count <- paste(data_plot_ct$group_ct_st, "(n =", data_plot_ct$count, ")")

table(data_plot_ct$group_with_count)

colors <- c(
  "Agreement (n = 258 )" = "blue",  # Blue
  "B Cells::PD-1+  (n = 180 )" = "#ff7f0e",   # Orange
  "CD4+ T Cells::PD-1-  (n = 133 )" = "#2ca02c",  # Green
  "Disagreement (n = 180 )" = "red",   # Red
  "Others (n = 329 )" = "grey90",  # Purple
  "VEC (n = 61 )" = "#8c564b"  # Brown
)

ggplot(data_plot_ct, aes(X, Y, group = -1L)) + 
  geom_voronoi_tile(aes(fill = group_with_count), colour = 'black', max.radius = 40) +
  theme_classic(base_size = 15) +
  labs(x = "X", y = "Y", title = "") + guides(fill = guide_legend(override.aes = list(size = 3),title = "Group"))+
  theme(
    axis.title = element_blank(),         # Remove axis titles
    axis.text = element_blank(),          # Remove axis text
    axis.ticks = element_blank(),         # Remove axis ticks
    axis.line = element_blank(),          # Remove axis lines
    panel.background = element_blank(),   # Remove panel background
    panel.grid = element_blank()          # Remove grid lines
  )+ggtitle("")+
  scale_fill_manual(values = colors)















Group1=NULL
for (i in 1:nrow(data_sub3)) {
  if(data_sub3$TACIT_PCF_v2[i]%in%c("B Cells","CD4+ T Cells")){
    if(data_sub3$PDCD1[i]>0){
      Group1[i]=1
    }else{
      Group1[i]=0
    }
  }else{
    Group1[i]="Others"
  }
}


Group2=NULL
for (i in 1:nrow(data_sub3)) {
  if(data_sub3$TACIT_PCF_v2[i]%in%c("B Cells","CD4+ T Cells")){
    if(data_sub3$P_PD.1[i]>0){
      Group2[i]=1
    }else{
      Group2[i]=0
    }
  }else{
    Group2[i]="Others"
  }
}


group_ct_st_2=NULL
for (i in 1:length(Group2)) {
  if(Group2[i]!="Others"){
    if(Group1[i]==Group2[i]){
      group_ct_st_2[i]="Agreement"
    }else{
      group_ct_st_2[i]="Disagreement"
    }
  }else{
    group_ct_st_2[i]="Others"
  }
}

group_counts <- as.data.frame(table(group_ct_st_2))
colnames(group_counts) <- c("group_ct_st", "count")
# Merge the counts back into the original data
data_plot_ct <- data.frame(X=data_sub3$X,Y=data_sub3$Y,group_ct_st=group_ct_st_2)

data_plot_ct <- merge(data_plot_ct, group_counts, by = "group_ct_st")

# Create a new column with the group names and counts
data_plot_ct$group_with_count <- paste(data_plot_ct$group_ct_st, "(n =", data_plot_ct$count, ")")

table(data_plot_ct$group_with_count)

colors <- c(
  "Agreement (n = 96 )" = "blue",  # Blue
  "B Cells::PD-1+  (n = 180 )" = "#ff7f0e",   # Orange
  "CD4+ T Cells::PD-1-  (n = 133 )" = "#2ca02c",  # Green
  "Disagreement (n = 41 )" = "red",   # Red
  "Others (n = 1410 )" = "grey90",  # Purple
  "VEC (n = 61 )" = "#8c564b"  # Brown
)

ggplot(data_plot_ct, aes(X, Y, group = -1L)) + 
  geom_voronoi_tile(aes(fill = group_with_count), colour = 'black', max.radius = 40) +
  theme_classic(base_size = 15) +
  labs(x = "X", y = "Y", title = "") + guides(fill = guide_legend(override.aes = list(size = 3),title = "Group"))+
  theme(
    axis.title = element_blank(),         # Remove axis titles
    axis.text = element_blank(),          # Remove axis text
    axis.ticks = element_blank(),         # Remove axis ticks
    axis.line = element_blank(),          # Remove axis lines
    panel.background = element_blank(),   # Remove panel background
    panel.grid = element_blank()          # Remove grid lines
  )+ggtitle("")+
  scale_fill_manual(values = colors)












Group1=NULL
for (i in 1:nrow(data_sub4)) {
  if(data_sub4$TACIT_PCF_v2[i]%in%c("B Cells","CD4+ T Cells")){
    if(data_sub4$PDCD1[i]>0){
      Group1[i]=1
    }else{
      Group1[i]=0
    }
  }else{
    Group1[i]="Others"
  }
}


Group2=NULL
for (i in 1:nrow(data_sub4)) {
  if(data_sub4$TACIT_PCF_v2[i]%in%c("B Cells","CD4+ T Cells")){
    if(data_sub4$P_PD.1[i]>0){
      Group2[i]=1
    }else{
      Group2[i]=0
    }
  }else{
    Group2[i]="Others"
  }
}


group_ct_st_2=NULL
for (i in 1:length(Group2)) {
  if(Group2[i]!="Others"){
    if(Group1[i]==Group2[i]){
      group_ct_st_2[i]="Agreement"
    }else{
      group_ct_st_2[i]="Disagreement"
    }
  }else{
    group_ct_st_2[i]="Others"
  }
}

group_counts <- as.data.frame(table(group_ct_st_2))
colnames(group_counts) <- c("group_ct_st", "count")
# Merge the counts back into the original data
data_plot_ct <- data.frame(X=data_sub4$X,Y=data_sub4$Y,group_ct_st=group_ct_st_2)

data_plot_ct <- merge(data_plot_ct, group_counts, by = "group_ct_st")

# Create a new column with the group names and counts
data_plot_ct$group_with_count <- paste(data_plot_ct$group_ct_st, "(n =", data_plot_ct$count, ")")

table(data_plot_ct$group_with_count)

colors <- c(
  "Agreement (n = 191 )" = "blue",  # Blue
  "B Cells::PD-1+  (n = 180 )" = "#ff7f0e",   # Orange
  "CD4+ T Cells::PD-1-  (n = 133 )" = "#2ca02c",  # Green
  "Disagreement (n = 365 )" = "red",   # Red
  "Others (n = 508 )" = "grey90",  # Purple
  "VEC (n = 61 )" = "#8c564b"  # Brown
)

ggplot(data_plot_ct, aes(X, Y, group = -1L)) + 
  geom_voronoi_tile(aes(fill = group_with_count), colour = 'black', max.radius = 40) +
  theme_classic(base_size = 15) +
  labs(x = "X", y = "Y", title = "") + guides(fill = guide_legend(override.aes = list(size = 3),title = "Group"))+
  theme(
    axis.title = element_blank(),         # Remove axis titles
    axis.text = element_blank(),          # Remove axis text
    axis.ticks = element_blank(),         # Remove axis ticks
    axis.line = element_blank(),          # Remove axis lines
    panel.background = element_blank(),   # Remove panel background
    panel.grid = element_blank()          # Remove grid lines
  )+ggtitle("")+
  scale_fill_manual(values = colors)

















Group1=NULL
for (i in 1:nrow(data_sub1)) {
  if(data_sub1$TACIT_Xenium_v2[i]%in%c("B Cells","CD4+ T Cells")){
    if(data_sub1$PDCD1[i]>0){
      Group1[i]=1
    }else{
      Group1[i]=0
    }
  }else{
    Group1[i]="Others"
  }
}


Group2=NULL
for (i in 1:nrow(data_sub1)) {
  if(data_sub1$TACIT_Xenium_v2[i]%in%c("B Cells","CD4+ T Cells")){
    if(data_sub1$P_PD.1[i]>0){
      Group2[i]=1
    }else{
      Group2[i]=0
    }
  }else{
    Group2[i]="Others"
  }
}


group_ct_st_2=NULL
for (i in 1:length(Group2)) {
  if(Group2[i]!="Others"){
    if(Group1[i]==Group2[i]){
      group_ct_st_2[i]="Agreement"
    }else{
      group_ct_st_2[i]="Disagreement"
    }
  }else{
    group_ct_st_2[i]="Others"
  }
}

group_counts <- as.data.frame(table(group_ct_st_2))
colnames(group_counts) <- c("group_ct_st", "count")
# Merge the counts back into the original data
data_plot_ct <- data.frame(X=data_sub1$X,Y=data_sub1$Y,group_ct_st=group_ct_st_2)

data_plot_ct <- merge(data_plot_ct, group_counts, by = "group_ct_st")

# Create a new column with the group names and counts
data_plot_ct$group_with_count <- paste(data_plot_ct$group_ct_st, "(n =", data_plot_ct$count, ")")

table(data_plot_ct$group_with_count)

colors <- c(
  "Agreement (n = 736 )" = "blue",  # Blue
  "B Cells::PD-1+  (n = 44 )" = "#ff7f0e",   # Orange
  "CD4+ T Cells::PD-1-  (n = 133 )" = "#2ca02c",  # Green
  "Disagreement (n = 129 )" = "red",   # Red
  "Others (n = 2179 )" = "grey90",  # Purple
  "VEC (n = 61 )" = "#8c564b"  # Brown
)

ggplot(data_plot_ct, aes(X, Y, group = -1L)) + 
  geom_voronoi_tile(aes(fill = group_with_count), colour = 'black', max.radius = 40) +
  theme_classic(base_size = 15) +
  labs(x = "X", y = "Y", title = "") + guides(fill = guide_legend(override.aes = list(size = 3),title = "Group"))+
  theme(
    axis.title = element_blank(),         # Remove axis titles
    axis.text = element_blank(),          # Remove axis text
    axis.ticks = element_blank(),         # Remove axis ticks
    axis.line = element_blank(),          # Remove axis lines
    panel.background = element_blank(),   # Remove panel background
    panel.grid = element_blank()          # Remove grid lines
  )+ggtitle("")+
  scale_fill_manual(values = colors)










Group1=NULL
for (i in 1:nrow(data_sub2)) {
  if(data_sub2$TACIT_Xenium_v2[i]%in%c("B Cells","CD4+ T Cells")){
    if(data_sub2$PDCD1[i]>0){
      Group1[i]=1
    }else{
      Group1[i]=0
    }
  }else{
    Group1[i]="Others"
  }
}


Group2=NULL
for (i in 1:nrow(data_sub2)) {
  if(data_sub2$TACIT_Xenium_v2[i]%in%c("B Cells","CD4+ T Cells")){
    if(data_sub2$P_PD.1[i]>0){
      Group2[i]=1
    }else{
      Group2[i]=0
    }
  }else{
    Group2[i]="Others"
  }
}


group_ct_st_2=NULL
for (i in 1:length(Group2)) {
  if(Group2[i]!="Others"){
    if(Group1[i]==Group2[i]){
      group_ct_st_2[i]="Agreement"
    }else{
      group_ct_st_2[i]="Disagreement"
    }
  }else{
    group_ct_st_2[i]="Others"
  }
}

group_counts <- as.data.frame(table(group_ct_st_2))
colnames(group_counts) <- c("group_ct_st", "count")
# Merge the counts back into the original data
data_plot_ct <- data.frame(X=data_sub2$X,Y=data_sub2$Y,group_ct_st=group_ct_st_2)

data_plot_ct <- merge(data_plot_ct, group_counts, by = "group_ct_st")

# Create a new column with the group names and counts
data_plot_ct$group_with_count <- paste(data_plot_ct$group_ct_st, "(n =", data_plot_ct$count, ")")

table(data_plot_ct$group_with_count)

colors <- c(
  "Agreement (n = 247 )" = "blue",  # Blue
  "B Cells::PD-1+  (n = 180 )" = "#ff7f0e",   # Orange
  "CD4+ T Cells::PD-1-  (n = 133 )" = "#2ca02c",  # Green
  "Disagreement (n = 145 )" = "red",   # Red
  "Others (n = 375 )" = "grey90",  # Purple
  "VEC (n = 61 )" = "#8c564b"  # Brown
)

ggplot(data_plot_ct, aes(X, Y, group = -1L)) + 
  geom_voronoi_tile(aes(fill = group_with_count), colour = 'black', max.radius = 40) +
  theme_classic(base_size = 15) +
  labs(x = "X", y = "Y", title = "") + guides(fill = guide_legend(override.aes = list(size = 3),title = "Group"))+
  theme(
    axis.title = element_blank(),         # Remove axis titles
    axis.text = element_blank(),          # Remove axis text
    axis.ticks = element_blank(),         # Remove axis ticks
    axis.line = element_blank(),          # Remove axis lines
    panel.background = element_blank(),   # Remove panel background
    panel.grid = element_blank()          # Remove grid lines
  )+ggtitle("")+
  scale_fill_manual(values = colors)















Group1=NULL
for (i in 1:nrow(data_sub3)) {
  if(data_sub3$TACIT_Xenium_v2[i]%in%c("B Cells","CD4+ T Cells")){
    if(data_sub3$PDCD1[i]>0){
      Group1[i]=1
    }else{
      Group1[i]=0
    }
  }else{
    Group1[i]="Others"
  }
}


Group2=NULL
for (i in 1:nrow(data_sub3)) {
  if(data_sub3$TACIT_Xenium_v2[i]%in%c("B Cells","CD4+ T Cells")){
    if(data_sub3$P_PD.1[i]>0){
      Group2[i]=1
    }else{
      Group2[i]=0
    }
  }else{
    Group2[i]="Others"
  }
}


group_ct_st_2=NULL
for (i in 1:length(Group2)) {
  if(Group2[i]!="Others"){
    if(Group1[i]==Group2[i]){
      group_ct_st_2[i]="Agreement"
    }else{
      group_ct_st_2[i]="Disagreement"
    }
  }else{
    group_ct_st_2[i]="Others"
  }
}

group_counts <- as.data.frame(table(group_ct_st_2))
colnames(group_counts) <- c("group_ct_st", "count")
# Merge the counts back into the original data
data_plot_ct <- data.frame(X=data_sub3$X,Y=data_sub3$Y,group_ct_st=group_ct_st_2)

data_plot_ct <- merge(data_plot_ct, group_counts, by = "group_ct_st")

# Create a new column with the group names and counts
data_plot_ct$group_with_count <- paste(data_plot_ct$group_ct_st, "(n =", data_plot_ct$count, ")")

table(data_plot_ct$group_with_count)

colors <- c(
  "Agreement (n = 149 )" = "blue",  # Blue
  "B Cells::PD-1+  (n = 180 )" = "#ff7f0e",   # Orange
  "CD4+ T Cells::PD-1-  (n = 133 )" = "#2ca02c",  # Green
  "Disagreement (n = 68 )" = "red",   # Red
  "Others (n = 1330 )" = "grey90",  # Purple
  "VEC (n = 61 )" = "#8c564b"  # Brown
)

ggplot(data_plot_ct, aes(X, Y, group = -1L)) + 
  geom_voronoi_tile(aes(fill = group_with_count), colour = 'black', max.radius = 40) +
  theme_classic(base_size = 15) +
  labs(x = "X", y = "Y", title = "") + guides(fill = guide_legend(override.aes = list(size = 3),title = "Group"))+
  theme(
    axis.title = element_blank(),         # Remove axis titles
    axis.text = element_blank(),          # Remove axis text
    axis.ticks = element_blank(),         # Remove axis ticks
    axis.line = element_blank(),          # Remove axis lines
    panel.background = element_blank(),   # Remove panel background
    panel.grid = element_blank()          # Remove grid lines
  )+ggtitle("")+
  scale_fill_manual(values = colors)












Group1=NULL
for (i in 1:nrow(data_sub4)) {
  if(data_sub4$TACIT_Xenium_v2[i]%in%c("B Cells","CD4+ T Cells")){
    if(data_sub4$PDCD1[i]>0){
      Group1[i]=1
    }else{
      Group1[i]=0
    }
  }else{
    Group1[i]="Others"
  }
}


Group2=NULL
for (i in 1:nrow(data_sub4)) {
  if(data_sub4$TACIT_Xenium_v2[i]%in%c("B Cells","CD4+ T Cells")){
    if(data_sub4$P_PD.1[i]>0){
      Group2[i]=1
    }else{
      Group2[i]=0
    }
  }else{
    Group2[i]="Others"
  }
}


group_ct_st_2=NULL
for (i in 1:length(Group2)) {
  if(Group2[i]!="Others"){
    if(Group1[i]==Group2[i]){
      group_ct_st_2[i]="Agreement"
    }else{
      group_ct_st_2[i]="Disagreement"
    }
  }else{
    group_ct_st_2[i]="Others"
  }
}

group_counts <- as.data.frame(table(group_ct_st_2))
colnames(group_counts) <- c("group_ct_st", "count")
# Merge the counts back into the original data
data_plot_ct <- data.frame(X=data_sub4$X,Y=data_sub4$Y,group_ct_st=group_ct_st_2)

data_plot_ct <- merge(data_plot_ct, group_counts, by = "group_ct_st")

# Create a new column with the group names and counts
data_plot_ct$group_with_count <- paste(data_plot_ct$group_ct_st, "(n =", data_plot_ct$count, ")")

table(data_plot_ct$group_with_count)

colors <- c(
  "Agreement (n = 208 )" = "blue",  # Blue
  "B Cells::PD-1+  (n = 180 )" = "#ff7f0e",   # Orange
  "CD4+ T Cells::PD-1-  (n = 133 )" = "#2ca02c",  # Green
  "Disagreement (n = 243 )" = "red",   # Red
  "Others (n = 613 )" = "grey90",  # Purple
  "VEC (n = 61 )" = "#8c564b"  # Brown
)

ggplot(data_plot_ct, aes(X, Y, group = -1L)) + 
  geom_voronoi_tile(aes(fill = group_with_count), colour = 'black', max.radius = 40) +
  theme_classic(base_size = 15) +
  labs(x = "X", y = "Y", title = "") + guides(fill = guide_legend(override.aes = list(size = 3),title = "Group"))+
  theme(
    axis.title = element_blank(),         # Remove axis titles
    axis.text = element_blank(),          # Remove axis text
    axis.ticks = element_blank(),         # Remove axis ticks
    axis.line = element_blank(),          # Remove axis lines
    panel.background = element_blank(),   # Remove panel background
    panel.grid = element_blank()          # Remove grid lines
  )+ggtitle("")+
  scale_fill_manual(values = colors)














Group1=NULL
for (i in 1:nrow(data_sub1)) {
  if(data_sub1$TACIT_PCF_Xenium_v2[i]%in%c("B Cells","CD4+ T Cells")){
    if(data_sub1$PDCD1[i]>0){
      Group1[i]=1
    }else{
      Group1[i]=0
    }
  }else{
    Group1[i]="Others"
  }
}


Group2=NULL
for (i in 1:nrow(data_sub1)) {
  if(data_sub1$TACIT_PCF_Xenium_v2[i]%in%c("B Cells","CD4+ T Cells")){
    if(data_sub1$P_PD.1[i]>0){
      Group2[i]=1
    }else{
      Group2[i]=0
    }
  }else{
    Group2[i]="Others"
  }
}


group_ct_st_2=NULL
for (i in 1:length(Group2)) {
  if(Group2[i]!="Others"){
    if(Group1[i]==Group2[i]){
      group_ct_st_2[i]="Agreement"
    }else{
      group_ct_st_2[i]="Disagreement"
    }
  }else{
    group_ct_st_2[i]="Others"
  }
}

group_counts <- as.data.frame(table(group_ct_st_2))
colnames(group_counts) <- c("group_ct_st", "count")
# Merge the counts back into the original data
data_plot_ct <- data.frame(X=data_sub1$X,Y=data_sub1$Y,group_ct_st=group_ct_st_2)

data_plot_ct <- merge(data_plot_ct, group_counts, by = "group_ct_st")

# Create a new column with the group names and counts
data_plot_ct$group_with_count <- paste(data_plot_ct$group_ct_st, "(n =", data_plot_ct$count, ")")

table(data_plot_ct$group_with_count)

colors <- c(
  "Agreement (n = 915 )" = "blue",  # Blue
  "B Cells::PD-1+  (n = 44 )" = "#ff7f0e",   # Orange
  "CD4+ T Cells::PD-1-  (n = 133 )" = "#2ca02c",  # Green
  "Disagreement (n = 237 )" = "red",   # Red
  "Others (n = 1892 )" = "grey90",  # Purple
  "VEC (n = 61 )" = "#8c564b"  # Brown
)

ggplot(data_plot_ct, aes(X, Y, group = -1L)) + 
  geom_voronoi_tile(aes(fill = group_with_count), colour = 'black', max.radius = 40) +
  theme_classic(base_size = 15) +
  labs(x = "X", y = "Y", title = "") + guides(fill = guide_legend(override.aes = list(size = 3),title = "Group"))+
  theme(
    axis.title = element_blank(),         # Remove axis titles
    axis.text = element_blank(),          # Remove axis text
    axis.ticks = element_blank(),         # Remove axis ticks
    axis.line = element_blank(),          # Remove axis lines
    panel.background = element_blank(),   # Remove panel background
    panel.grid = element_blank()          # Remove grid lines
  )+ggtitle("")+
  scale_fill_manual(values = colors)










Group1=NULL
for (i in 1:nrow(data_sub2)) {
  if(data_sub2$TACIT_PCF_Xenium_v2[i]%in%c("B Cells","CD4+ T Cells")){
    if(data_sub2$PDCD1[i]>0){
      Group1[i]=1
    }else{
      Group1[i]=0
    }
  }else{
    Group1[i]="Others"
  }
}


Group2=NULL
for (i in 1:nrow(data_sub2)) {
  if(data_sub2$TACIT_PCF_Xenium_v2[i]%in%c("B Cells","CD4+ T Cells")){
    if(data_sub2$P_PD.1[i]>0){
      Group2[i]=1
    }else{
      Group2[i]=0
    }
  }else{
    Group2[i]="Others"
  }
}


group_ct_st_2=NULL
for (i in 1:length(Group2)) {
  if(Group2[i]!="Others"){
    if(Group1[i]==Group2[i]){
      group_ct_st_2[i]="Agreement"
    }else{
      group_ct_st_2[i]="Disagreement"
    }
  }else{
    group_ct_st_2[i]="Others"
  }
}

group_counts <- as.data.frame(table(group_ct_st_2))
colnames(group_counts) <- c("group_ct_st", "count")
# Merge the counts back into the original data
data_plot_ct <- data.frame(X=data_sub2$X,Y=data_sub2$Y,group_ct_st=group_ct_st_2)

data_plot_ct <- merge(data_plot_ct, group_counts, by = "group_ct_st")

# Create a new column with the group names and counts
data_plot_ct$group_with_count <- paste(data_plot_ct$group_ct_st, "(n =", data_plot_ct$count, ")")

table(data_plot_ct$group_with_count)

colors <- c(
  "Agreement (n = 279 )" = "blue",  # Blue
  "B Cells::PD-1+  (n = 180 )" = "#ff7f0e",   # Orange
  "CD4+ T Cells::PD-1-  (n = 133 )" = "#2ca02c",  # Green
  "Disagreement (n = 208 )" = "red",   # Red
  "Others (n = 280 )" = "grey90",  # Purple
  "VEC (n = 61 )" = "#8c564b"  # Brown
)

ggplot(data_plot_ct, aes(X, Y, group = -1L)) + 
  geom_voronoi_tile(aes(fill = group_with_count), colour = 'black', max.radius = 40) +
  theme_classic(base_size = 15) +
  labs(x = "X", y = "Y", title = "") + guides(fill = guide_legend(override.aes = list(size = 3),title = "Group"))+
  theme(
    axis.title = element_blank(),         # Remove axis titles
    axis.text = element_blank(),          # Remove axis text
    axis.ticks = element_blank(),         # Remove axis ticks
    axis.line = element_blank(),          # Remove axis lines
    panel.background = element_blank(),   # Remove panel background
    panel.grid = element_blank()          # Remove grid lines
  )+ggtitle("")+
  scale_fill_manual(values = colors)















Group1=NULL
for (i in 1:nrow(data_sub3)) {
  if(data_sub3$TACIT_PCF_Xenium_v2[i]%in%c("B Cells","CD4+ T Cells")){
    if(data_sub3$PDCD1[i]>0){
      Group1[i]=1
    }else{
      Group1[i]=0
    }
  }else{
    Group1[i]="Others"
  }
}


Group2=NULL
for (i in 1:nrow(data_sub3)) {
  if(data_sub3$TACIT_PCF_Xenium_v2[i]%in%c("B Cells","CD4+ T Cells")){
    if(data_sub3$P_PD.1[i]>0){
      Group2[i]=1
    }else{
      Group2[i]=0
    }
  }else{
    Group2[i]="Others"
  }
}


group_ct_st_2=NULL
for (i in 1:length(Group2)) {
  if(Group2[i]!="Others"){
    if(Group1[i]==Group2[i]){
      group_ct_st_2[i]="Agreement"
    }else{
      group_ct_st_2[i]="Disagreement"
    }
  }else{
    group_ct_st_2[i]="Others"
  }
}

group_counts <- as.data.frame(table(group_ct_st_2))
colnames(group_counts) <- c("group_ct_st", "count")
# Merge the counts back into the original data
data_plot_ct <- data.frame(X=data_sub3$X,Y=data_sub3$Y,group_ct_st=group_ct_st_2)

data_plot_ct <- merge(data_plot_ct, group_counts, by = "group_ct_st")

# Create a new column with the group names and counts
data_plot_ct$group_with_count <- paste(data_plot_ct$group_ct_st, "(n =", data_plot_ct$count, ")")

table(data_plot_ct$group_with_count)

colors <- c(
  "Agreement (n = 197 )" = "blue",  # Blue
  "B Cells::PD-1+  (n = 180 )" = "#ff7f0e",   # Orange
  "CD4+ T Cells::PD-1-  (n = 133 )" = "#2ca02c",  # Green
  "Disagreement (n = 98 )" = "red",   # Red
  "Others (n = 1252 )" = "grey90",  # Purple
  "VEC (n = 61 )" = "#8c564b"  # Brown
)

ggplot(data_plot_ct, aes(X, Y, group = -1L)) + 
  geom_voronoi_tile(aes(fill = group_with_count), colour = 'black', max.radius = 40) +
  theme_classic(base_size = 15) +
  labs(x = "X", y = "Y", title = "") + guides(fill = guide_legend(override.aes = list(size = 3),title = "Group"))+
  theme(
    axis.title = element_blank(),         # Remove axis titles
    axis.text = element_blank(),          # Remove axis text
    axis.ticks = element_blank(),         # Remove axis ticks
    axis.line = element_blank(),          # Remove axis lines
    panel.background = element_blank(),   # Remove panel background
    panel.grid = element_blank()          # Remove grid lines
  )+ggtitle("")+
  scale_fill_manual(values = colors)












Group1=NULL
for (i in 1:nrow(data_sub4)) {
  if(data_sub4$TACIT_PCF_Xenium_v2[i]%in%c("B Cells","CD4+ T Cells")){
    if(data_sub4$PDCD1[i]>0){
      Group1[i]=1
    }else{
      Group1[i]=0
    }
  }else{
    Group1[i]="Others"
  }
}


Group2=NULL
for (i in 1:nrow(data_sub4)) {
  if(data_sub4$TACIT_PCF_Xenium_v2[i]%in%c("B Cells","CD4+ T Cells")){
    if(data_sub4$P_PD.1[i]>0){
      Group2[i]=1
    }else{
      Group2[i]=0
    }
  }else{
    Group2[i]="Others"
  }
}


group_ct_st_2=NULL
for (i in 1:length(Group2)) {
  if(Group2[i]!="Others"){
    if(Group1[i]==Group2[i]){
      group_ct_st_2[i]="Agreement"
    }else{
      group_ct_st_2[i]="Disagreement"
    }
  }else{
    group_ct_st_2[i]="Others"
  }
}

group_counts <- as.data.frame(table(group_ct_st_2))
colnames(group_counts) <- c("group_ct_st", "count")
# Merge the counts back into the original data
data_plot_ct <- data.frame(X=data_sub4$X,Y=data_sub4$Y,group_ct_st=group_ct_st_2)

data_plot_ct <- merge(data_plot_ct, group_counts, by = "group_ct_st")

# Create a new column with the group names and counts
data_plot_ct$group_with_count <- paste(data_plot_ct$group_ct_st, "(n =", data_plot_ct$count, ")")

table(data_plot_ct$group_with_count)

colors <- c(
  "Agreement (n = 236 )" = "blue",  # Blue
  "B Cells::PD-1+  (n = 180 )" = "#ff7f0e",   # Orange
  "CD4+ T Cells::PD-1-  (n = 133 )" = "#2ca02c",  # Green
  "Disagreement (n = 410 )" = "red",   # Red
  "Others (n = 418 )" = "grey90",  # Purple
  "VEC (n = 61 )" = "#8c564b"  # Brown
)

ggplot(data_plot_ct, aes(X, Y, group = -1L)) + 
  geom_voronoi_tile(aes(fill = group_with_count), colour = 'black', max.radius = 40) +
  theme_classic(base_size = 15) +
  labs(x = "X", y = "Y", title = "") + guides(fill = guide_legend(override.aes = list(size = 3),title = "Group"))+
  theme(
    axis.title = element_blank(),         # Remove axis titles
    axis.text = element_blank(),          # Remove axis text
    axis.ticks = element_blank(),         # Remove axis ticks
    axis.line = element_blank(),          # Remove axis lines
    panel.background = element_blank(),   # Remove panel background
    panel.grid = element_blank()          # Remove grid lines
  )+ggtitle("")+
  scale_fill_manual(values = colors)












