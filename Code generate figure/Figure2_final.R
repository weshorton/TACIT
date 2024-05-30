# Load necessary libraries
library(readr)
library(readxl)
library(uwot)
library(ggpubr)
library(ggplot2)
library(caret)
library(pROC)
library(pheatmap)
library(Seurat)
library(dplyr)
library(reticulate)
library(umap)
library(tidyr)

# Set working directory
setwd("~/Documents/Compare_TACIT/PCF/")

# Read data files
data <- read_csv("~/Documents/TACIT figure/dataset/Figure 2/CRC_DataSet.csv")
Signature <- readr::read_csv("~/Downloads/test_crc_signature.csv")
outcome <- read_csv("~/Documents/TACIT figure/dataset/Figure 2/CRC_outcome.csv")

# Prepare the Signature data
colnames(Signature)[1] = "cell_type"
Signature[is.na(Signature)] = 0  # Assign all NA as 0
Signature = as.data.frame(Signature)


#--------------------------
#a------------------------
#--------------------------

TACIT_val=calculate_metrics(outcome$TACIT[which(outcome$TACIT!="Others")],outcome$reference[which(outcome$TACIT!="Others")])
SCINA_val=calculate_metrics(outcome$SCINA[which(outcome$SCINA!="Others")],outcome$reference[which(outcome$SCINA!="Others")])
Celesta_val=calculate_metrics(outcome$CELESTA[which(outcome$CELESTA!="Others")],outcome$reference[which(outcome$CELESTA!="Others")])
Louvain_val002=calculate_metrics(outcome$Seurat_orginial[which(outcome$Seurat_orginial!="Others")],outcome$reference[which(outcome$Seurat_orginial!="Others")])

# Assuming you have the weighted recall, precision, and F1 scores for each method
data_plot <- data.frame(
  Method = c("TACIT","CELESTA","SCINA","Louvain"),
  WeightedRecall = c(TACIT_val$WeightedRecall,Celesta_val$WeightedRecall,SCINA_val$WeightedRecall,Louvain_val002$WeightedRecall),
  WeightedPrecision = c(TACIT_val$WeightedPrecision,Celesta_val$WeightedPrecision,SCINA_val$WeightedPrecision,Louvain_val002$WeightedPrecision),
  WeightedF1 = c(TACIT_val$WeightedF1,Celesta_val$WeightedF1,SCINA_val$WeightedF1,Louvain_val002$WeightedF1)
)


colnames(data_plot)=c("Method","Weighted Recall","Weighted Precision","Weighted F1")
# Transform data from wide to long format
long_data <- pivot_longer(data_plot, cols = -Method, names_to = "Metric", values_to = "Score")

long_data$Method=factor(long_data$Method,levels = c("TACIT","CELESTA","SCINA","Louvain"))
long_data$Metric=factor(long_data$Metric,levels = c("Weighted Recall","Weighted Precision","Weighted F1"))
# Plot
library(ggplot2)
ggplot(long_data, aes(x = Method, y = Score, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.7) +
  geom_text(aes(label = sprintf("%.2f", Score)), 
            position = position_dodge(width = 0.9), vjust = -0.25, color = "black", size = 6) +
  facet_wrap(~ Metric, scales = "free_y") +
  labs(x = "Cell Type", y = "Score", title = "") +
  scale_fill_manual(values = c("TACIT" = "red", "CELESTA" = "purple", "SCINA" = "green", "Louvain" = "blue")) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.6),
        legend.position = "right")+ylim(0,1)


#--------------------------
#b------------------------
#--------------------------

  
# Create a dataframe with the necessary columns
data_plot <- data.frame(
  ref = outcome$reference,
  CELESTA = outcome$CELESTA,
  SCINA = outcome$SCINA,
  Louvain = outcome$Seurat_orginial,
  TACIT = outcome$TACIT,
  TMA = data$`File Name`
)

# Convert the vector to a dataframe for easier manipulation
df <- data.frame(cell_type = data_plot$ref)

# Calculate the proportion of each cell type
proportions <- df %>%
  group_by(cell_type) %>%
  summarise(count = n()) %>%
  mutate(proportion = count / sum(count)) %>%
  filter(proportion < 0.01)


  
data_plot$ref=ifelse(data_plot$ref%in%c(proportions$cell_type),"Rare cell type","Not rare")
data_plot$CELESTA=ifelse(data_plot$CELESTA%in%c(proportions$cell_type),"Rare cell type","Not rare")
data_plot$TACIT=ifelse(data_plot$TACIT%in%c(proportions$cell_type),"Rare cell type","Not rare")
data_plot$Louvain=ifelse(data_plot$Louvain%in%c(proportions$cell_type),"Rare cell type","Not rare")
data_plot$SCINA=ifelse(data_plot$SCINA%in%c(proportions$cell_type),"Rare cell type","Not rare")




# Assuming your dataframe is named data_plot
# Convert factors to character to avoid drop of unused levels in subsequent filtering
data_plot <- data_plot %>% mutate(across(c(ref, TACIT,Louvain,SCINA,CELESTA), as.character))

# Summarizing counts for Immune cells (all)
summary_data <- data_plot %>% 
  group_by(TMA) %>% 
  summarise(Ref = sum(ref == "Rare cell type"),
            TACIT = sum(TACIT == "Rare cell type"),
            SCINA = sum(SCINA == "Rare cell type"),
            CELESTA = sum(CELESTA == "Rare cell type"),
            Louvain = sum(Louvain == "Rare cell type"))


method_name="TACIT"

# Calculating correlation
correlation <- round(cor(summary_data$Ref, summary_data[[method_name]], use = "complete.obs"), 2)

# Plot the correlation
ggplot(summary_data, aes(y = Ref, x = summary_data[[method_name]])) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red",size=2) + # x=y line
  theme_classic(base_size = 30) +
  labs(title = "",
       y = "Number of cells in Reference",
       x = paste("Number of cells in", method_name))  +
  theme(plot.title = element_text(size = 14)) +
  annotate("text", x = max(summary_data$Ref, na.rm = TRUE) * 0.2, y = max(summary_data[[method_name]], na.rm = TRUE), label = paste("R =", correlation), size = 1)+ geom_smooth(method = "lm", color = "blue", se = FALSE,size=2)





# Summarizing counts for Immune cells (all)
summary_data <- data_plot %>% 
  group_by(TMA) %>% 
  summarise(Ref = sum(ref != "Rare cell type"),
            TACIT = sum(TACIT != "Rare cell type"),
            SCINA = sum(SCINA != "Rare cell type"),
            CELESTA = sum(CELESTA != "Rare cell type"),
            Louvain = sum(Louvain != "Rare cell type"))


method_name="TACIT"

# Calculating correlation
correlation <- round(cor(summary_data$Ref, summary_data[[method_name]], use = "complete.obs"), 2)

# Plot the correlation
ggplot(summary_data, aes(y = Ref, x = summary_data[[method_name]])) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red",size=2) + # x=y line
  theme_classic(base_size = 30) +
  labs(title = "",
       y = "Number of cells in Reference",
       x = paste("Number of cells in", method_name))  +
  theme(plot.title = element_text(size = 14)) +
  annotate("text", x = max(summary_data$Ref, na.rm = TRUE) * 0.2, y = max(summary_data[[method_name]], na.rm = TRUE), label = paste("R =", correlation), size = 1)+ geom_smooth(method = "lm", color = "blue", se = FALSE,size=2)







#--------------------------
#c------------------------
#--------------------------

#-TACIT

# Create a data frame combining TACIT outcomes and signature data
data_plot <- data.frame(TACIT = outcome$TACIT, (data[, colnames(Signature)[-1]]))
colnames(data_plot)[-1] <- colnames(Signature)[-1]

# Prepare the Signature data
colnames(Signature)[1] <- "cell_type"
Signature[is.na(Signature)] <- 0

# Calculate the median values for each cell type in TACIT
mean_values_TACIT <- data_plot %>%
  group_by(TACIT) %>%
  summarise_all(~quantile(., 0.5)) %>%
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
    
    testt <- log2(mean(as.numeric(data_anb[which(outcome$TACIT == ct_choose), marker_choose])) /
                    mean(as.numeric(data_anb[which(outcome$TACIT %in% select_ct_0), marker_choose])))
    
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




#-SCINA

# Create a data frame combining TACIT outcomes and signature data
data_plot <- data.frame(TACIT = outcome$SCINA, (data[, colnames(Signature)[-1]]))
colnames(data_plot)[-1] <- colnames(Signature)[-1]

# Prepare the Signature data
colnames(Signature)[1] <- "cell_type"
Signature[is.na(Signature)] <- 0

# Calculate the median values for each cell type in TACIT
mean_values_TACIT <- data_plot %>%
  group_by(TACIT) %>%
  summarise_all(~quantile(., 0.5)) %>%
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
    
    testt <- log2(mean(as.numeric(data_anb[which(outcome$SCINA == ct_choose), marker_choose])) /
                    mean(as.numeric(data_anb[which(outcome$SCINA %in% select_ct_0), marker_choose])))
    
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
data_final_SCINA <- data_final2










#-CELESTA

# Create a data frame combining TACIT outcomes and signature data
data_plot <- data.frame(TACIT = outcome$CELESTA, (data[, colnames(Signature)[-1]]))
colnames(data_plot)[-1] <- colnames(Signature)[-1]

# Prepare the Signature data
colnames(Signature)[1] <- "cell_type"
Signature[is.na(Signature)] <- 0

# Calculate the median values for each cell type in TACIT
mean_values_TACIT <- data_plot %>%
  group_by(TACIT) %>%
  summarise_all(~quantile(., 0.5)) %>%
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
    
    testt <- log2(mean(as.numeric(data_anb[which(outcome$CELESTA == ct_choose), marker_choose])) /
                    mean(as.numeric(data_anb[which(outcome$CELESTA %in% select_ct_0), marker_choose])))
    
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
data_final_CELESTA <- data_final2






#-Louvain

# Create a data frame combining TACIT outcomes and signature data
data_plot <- data.frame(TACIT = outcome$Seurat_orginial, (data[, colnames(Signature)[-1]]))
colnames(data_plot)[-1] <- colnames(Signature)[-1]

# Prepare the Signature data
colnames(Signature)[1] <- "cell_type"
Signature[is.na(Signature)] <- 0

# Calculate the median values for each cell type in TACIT
mean_values_TACIT <- data_plot %>%
  group_by(TACIT) %>%
  summarise_all(~quantile(., 0.5)) %>%
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
    
    testt <- log2(mean(as.numeric(data_anb[which(outcome$Seurat_orginial == ct_choose), marker_choose])) /
                    mean(as.numeric(data_anb[which(outcome$Seurat_orginial %in% select_ct_0), marker_choose])))
    
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
data_final_Louvain <- data_final2










#-Reference

# Create a data frame combining TACIT outcomes and signature data
data_plot <- data.frame(TACIT = outcome$reference, (data[, colnames(Signature)[-1]]))
colnames(data_plot)[-1] <- colnames(Signature)[-1]

# Prepare the Signature data
colnames(Signature)[1] <- "cell_type"
Signature[is.na(Signature)] <- 0

# Calculate the median values for each cell type in TACIT
mean_values_TACIT <- data_plot %>%
  group_by(TACIT) %>%
  summarise_all(~quantile(., 0.5)) %>%
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
    
    testt <- log2(mean(as.numeric(data_anb[which(outcome$reference == ct_choose), marker_choose])) /
                    mean(as.numeric(data_anb[which(outcome$reference %in% select_ct_0), marker_choose])))
    
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
data_final_reference <- data_final2



data_final_CELESTA$Methods="CELESTA"
data_final_reference$Methods="Reference"
data_final_Louvain$Methods="Louvain"
data_final_SCINA$Methods="SCINA"
data_final_TACIT$Methods="TACIT"


data_final_plot=rbind(data_final_CELESTA,data_final_reference,data_final_Louvain,data_final_SCINA,data_final_TACIT)

data_final_plot$Group_plot=ifelse(data_final_plot$Group=="Significant & Signature",data_final_plot$Methods,"Others")

data_final_plot_sub=data_final_plot[which(data_final_plot$Signature==1),]

# Filter the data to include only the specified methods
filtered_data <- data_final_plot_sub %>%
  filter(Methods %in% c("TACIT", "CELESTA", "SCINA", "Louvain","Reference"))

filtered_data$Methods=factor(filtered_data$Methods,level=c("Reference",
                                                           "TACIT","CELESTA","SCINA","Louvain"))

#logFC
ggplot(filtered_data, aes(x = Methods, y = logFC, fill = Methods)) +
  geom_boxplot(alpha = 1, outlier.size = 2, position = position_dodge(0.8)) +
  scale_fill_manual(values = c("TACIT" = "red", "CELESTA" = "purple", "SCINA" = "green", "Louvain" = "blue", "Reference" = "yellow")) +
  theme_classic(base_size = 20) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  ylab("") +
  xlab("") +
  geom_signif(comparisons = list( c("TACIT", "SCINA"), c("TACIT", "Reference"), c("TACIT", "CELESTA"), c("TACIT", "Louvain")),
              map_signif_level=TRUE,
              textsize = 3, vjust = -0.2)


#-log10(padjust)
ggplot(filtered_data, aes(x = Methods, y = -log10(padjust), fill = Methods)) +
  geom_boxplot(alpha = 1, outlier.size = 2, position = position_dodge(0.8)) +
  scale_fill_manual(values = c("TACIT" = "red", "CELESTA" = "purple", "SCINA" = "green", "Louvain" = "blue", "Reference" = "yellow")) +
  theme_classic(base_size = 20) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  ylab("") +
  xlab("") +
  geom_signif(comparisons = list( c("TACIT", "SCINA"), c("TACIT", "Reference"), c("TACIT", "CELESTA"), c("TACIT", "Louvain")),
              map_signif_level=TRUE,
              textsize = 3, vjust = -0.2)






#--------------------------
#d------------------------
#--------------------------

data_sub=data[which(data$`File Name`=="reg008_B"),]
# Define a function to create a ggplot for given predicted outcome
plot_UMAP <- function(predicted_outcome, title) {
  df = data.frame(x = data_sub$X, y = data_sub$Y, predicted = predicted_outcome)
  ggplot(df, aes(x, y)) + 
    geom_point(aes(color = predicted), size = 3, alpha = 1) +
    theme_classic(base_size = 15) +
    theme(legend.position = "none") +
    labs(x = "X", y = "Y", title = title) +
    guides(fill = guide_legend(override.aes = list(size = 3), title = "Cell type")) +
    scale_color_manual(values = c("B cells" = "tan", "CD11c+ DCs" = "darkred", "CD3+ T cells" = "pink",
                                  "CD4+ T cells" = "violet", "CD4+ T cells CD45RO+" = "darkviolet",
                                  "macrophages" = "green", "CD8+ T cells" = "darkgreen",
                                  "granulocytes" = "blue", "lymphatics" = "skyblue", "NK cells" = "yellow",
                                  "Others" = "grey90", "immune cells" = "darkblue", "plasma cells" = "orange", 
                                  "stroma" = "gold", "Tregs" = "black", "tumor cells" = "red", 
                                  "vasculature" = "darkcyan", "smooth muscle" = "tan4"))
}

# Plot UMAP for each outcome
plot_UMAP(outcome$TACIT[which(data$`File Name` == "reg008_B")], "TACIT") + ggtitle("TACIT")
plot_UMAP(outcome$CELESTA[which(data$`File Name` == "reg008_B")], "CELESTA") + ggtitle("CELESTA")
plot_UMAP(outcome$Seurat_orginial[which(data$`File Name` == "reg008_B")], "Seurat Original") + ggtitle("Seurat Original")
plot_UMAP(outcome$reference[which(data$`File Name` == "reg008_B")], "Reference") + ggtitle("Reference")
plot_UMAP(outcome$SCINA[which(data$`File Name` == "reg008_B")], "SCINA") + ggtitle("SCINA")






#--------------------------
#e------------------------
#--------------------------

# Perform UMAP dimensionality reduction
data_umap = umap::umap(data[, colnames(Signature)[-1]], scale = TRUE, method = "umap-learn", metric = "correlation")
data_umap = data_umap$layout
colnames(data_umap) = c("UMAP1", "UMAP2")
data_umap = as.data.frame(data_umap)

# Plot UMAP with TACIT outcomes
ggplot(data_umap, aes(x = UMAP1, y = UMAP2, color = as.factor(outcome$TACIT))) +
  geom_point(size = 0.1) +
  scale_color_manual(values = c("B cells" = "tan", "CD11c+ DCs" = "darkred", "CD3+ T cells" = "pink",
                                "CD4+ T cells" = "violet", "CD4+ T cells CD45RO+" = "darkviolet",
                                "macrophages" = "green", "CD8+ T cells" = "darkgreen",
                                "granulocytes" = "blue", "lymphatics" = "skyblue", "NK cells" = "yellow",
                                "Others" = "grey90", "immune cells" = "darkblue", "plasma cells" = "orange", "stroma" = "gold",
                                "Tregs" = "black", "tumor cells" = "red", "vasculature" = "darkcyan", "smooth muscle" = "tan4")) +
  theme_classic(base_size = 25) +
  guides(color = guide_legend(override.aes = list(size = 3), title = "Cell type")) +
  theme(axis.line = element_blank(),  # Remove axis lines
        axis.ticks = element_blank(),  # Remove tick marks
        axis.text.x = element_blank(),  # Remove x axis text
        axis.text.y = element_blank()) +  # Remove y axis text
  ylab("") + xlab("") +  # Remove axis labels
  theme(legend.position = "none")  # Remove legend

# Plot UMAP with reference outcomes
ggplot(data_umap, aes(x = UMAP1, y = UMAP2, color = as.factor(outcome$reference))) +
  geom_point(size = 0.1) +
  scale_color_manual(values = c("B cells" = "tan", "CD11c+ DCs" = "darkred", "CD3+ T cells" = "pink",
                                "CD4+ T cells" = "violet", "CD4+ T cells CD45RO+" = "darkviolet",
                                "macrophages" = "green", "CD8+ T cells" = "darkgreen",
                                "granulocytes" = "blue", "lymphatics" = "skyblue", "NK cells" = "yellow",
                                "Others" = "grey90", "immune cells" = "darkblue", "plasma cells" = "orange", "stroma" = "gold",
                                "Tregs" = "black", "tumor cells" = "red", "vasculature" = "darkcyan", "smooth muscle" = "tan4")) +
  theme_classic(base_size = 25) +
  guides(color = guide_legend(override.aes = list(size = 3), title = "Cell type")) +
  theme(axis.line = element_blank(),  # Remove axis lines
        axis.ticks = element_blank(),  # Remove tick marks
        axis.text.x = element_blank(),  # Remove x axis text
        axis.text.y = element_blank()) +  # Remove y axis text
  ylab("") + xlab("") +  # Remove axis labels
  theme(legend.position = "none")  # Remove legend

# Plot UMAP with SCINA outcomes
ggplot(data_umap, aes(x = UMAP1, y = UMAP2, color = as.factor(outcome$SCINA))) +
  geom_point(size = 0.1) +
  scale_color_manual(values = c("B cells" = "tan", "CD11c+ DCs" = "darkred", "CD3+ T cells" = "pink",
                                "CD4+ T cells" = "violet", "CD4+ T cells CD45RO+" = "darkviolet",
                                "macrophages" = "green", "CD8+ T cells" = "darkgreen",
                                "granulocytes" = "blue", "lymphatics" = "skyblue", "NK cells" = "yellow",
                                "Others" = "grey90", "immune cells" = "darkblue", "plasma cells" = "orange", "stroma" = "gold",
                                "Tregs" = "black", "tumor cells" = "red", "vasculature" = "darkcyan", "smooth muscle" = "tan4")) +
  theme_classic(base_size = 25) +
  guides(color = guide_legend(override.aes = list(size = 3), title = "Cell type")) +
  theme(axis.line = element_blank(),  # Remove axis lines
        axis.ticks = element_blank(),  # Remove tick marks
        axis.text.x = element_blank(),  # Remove x axis text
        axis.text.y = element_blank()) +  # Remove y axis text
  ylab("") + xlab("") +  # Remove axis labels
  theme(legend.position = "none")  # Remove legend

# Plot UMAP with Seurat original outcomes
ggplot(data_umap, aes(x = UMAP1, y = UMAP2, color = as.factor(outcome$Seurat_orginial))) +
  geom_point(size = 0.1) +
  scale_color_manual(values = c("B cells" = "tan", "CD11c+ DCs" = "darkred", "CD3+ T cells" = "pink",
                                "CD4+ T cells" = "violet", "CD4+ T cells CD45RO+" = "darkviolet",
                                "macrophages" = "green", "CD8+ T cells" = "darkgreen",
                                "granulocytes" = "blue", "lymphatics" = "skyblue", "NK cells" = "yellow",
                                "Others" = "grey90", "immune cells" = "darkblue", "plasma cells" = "orange", "stroma" = "gold",
                                "Tregs" = "black", "tumor cells" = "red", "vasculature" = "darkcyan", "smooth muscle" = "tan4")) +
  theme_classic(base_size = 25) +
  guides(color = guide_legend(override.aes = list(size = 3), title = "Cell type")) +
  theme(axis.line = element_blank(),  # Remove axis lines
        axis.ticks = element_blank(),  # Remove tick marks
        axis.text.x = element_blank(),  # Remove x axis text
        axis.text.y = element_blank()) +  # Remove y axis text
  ylab("") + xlab("") +  # Remove axis labels
  theme(legend.position = "none")  # Remove legend






#--------------------------
#f------------------------
#--------------------------

#TACIT

# Combine TACIT outcomes with relevant data columns
data_plot <- data.frame(TACIT = outcome$TACIT, (data[, colnames(Signature)[-1]]))
colnames(data_plot)[-1] <- colnames(Signature)[-1]

# Calculate median values for each cell type in TACIT
mean_values_TACIT <- data_plot %>%
  group_by(TACIT) %>%
  summarise_all(~quantile(., 0.5)) %>%
  as.data.frame()

# Remove the 19th row and the "Others" cell type
mean_values_TACIT <- mean_values_TACIT[-19,]
mean_values_TACIT <- mean_values_TACIT[which(mean_values_TACIT$TACIT != "Others"),]
rownames(mean_values_TACIT) <- mean_values_TACIT$TACIT
mean_values_TACIT <- mean_values_TACIT[, -1]

# Reorder columns based on signature order
ordered_index <- match(Signature_order[, 1], rownames(mean_values_TACIT))
mean_values_TACIT <- mean_values_TACIT[ordered_index, ]

# Define color breaks and palette for heatmap
my.breaks <- c(seq(-2, 0, by = 0.1), seq(0.1, 2, by = 0.1))
my.colors <- c(colorRampPalette(colors = c("blue", "white"))(length(my.breaks) / 2), 
               colorRampPalette(colors = c("white", "red"))(length(my.breaks) / 2))

# Original to short name mapping for markers
short_names <- c(
  "CD31 - vasculature:Cyc_19_ch_3" = "CD31",
  "CD34 - vasculature:Cyc_20_ch_3" = "CD34",
  "Cytokeratin - epithelia:Cyc_10_ch_2" = "Cytokeratin",
  "Ki67 - proliferation:Cyc_5_ch_4" = "Ki67",
  "p53 - tumor suppressor:Cyc_3_ch_3" = "p53",
  "EGFR - signaling:Cyc_13_ch_3" = "EGFR",
  "Granzyme B - cytotoxicity:Cyc_13_ch_2" = "GranzymeB",
  "beta-catenin - Wnt signaling:Cyc_4_ch_4" = "BetaCatenin",
  "PD-L1 - checkpoint:Cyc_5_ch_3" = "PDL1",
  "aSMA - smooth muscle:Cyc_11_ch_2" = "aSMA",
  "CD44 - stroma:Cyc_2_ch_2" = "CD44",
  "Vimentin - cytoplasm:Cyc_8_ch_2" = "Vimentin",
  "Podoplanin - lymphatics:Cyc_19_ch_4" = "Podoplanin",
  "CD45 - hematopoietic cells:Cyc_4_ch_2" = "CD45",
  "CD15 - granulocytes:Cyc_14_ch_2" = "CD15",
  "CD11b - macrophages:Cyc_10_ch_3" = "CD11b",
  "HLA-DR - MHC-II:Cyc_5_ch_2" = "HLA-DR",
  "CD20 - B cells:Cyc_8_ch_3" = "CD20",
  "CD11c - DCs:Cyc_12_ch_3" = "CD11c",
  "CD163 - macrophages:Cyc_17_ch_3" = "CD163",
  "CD68 - macrophages:Cyc_18_ch_4" = "CD68",
  "CD38 - multifunctional:Cyc_20_ch_4" = "CD38",
  "CD138 - plasma cells:Cyc_21_ch_3" = "CD138",
  "CD56 - NK cells:Cyc_10_ch_4" = "CD56",
  "CD57 - NK cells:Cyc_17_ch_4" = "CD57",
  "CD3 - T cells:Cyc_16_ch_4" = "CD3",
  "CD7 - T cells:Cyc_16_ch_3" = "CD7",
  "CD8 - cytotoxic T cells:Cyc_3_ch_2" = "CD8",
  "FOXP3 - regulatory T cells:Cyc_2_ch_3" = "FOXP3",
  "CD4 - T helper cells:Cyc_6_ch_3" = "CD4",
  "CD45RO - memory cells:Cyc_18_ch_3" = "CD45RO"
)

# Rename columns in mean_values_TACIT and Signature_order based on short names
colnames(mean_values_TACIT) <- short_names[colnames(mean_values_TACIT)]
colnames(Signature_order) <- short_names[colnames(Signature_order)]

# Scale mean values for heatmap visualization
aa <- scale(mean_values_TACIT)
rownames(aa) <- Signature_order[, 1]

# Generate heatmap
pheatmap(aa, cluster_cols = FALSE,
         cluster_rows = FALSE,
         show_colnames = TRUE,
         fontsize_col = 12,
         fontsize_row = 12)  # Adjust fontsize_row as needed to increase x-axis text size






#SCINA

# Combine TACIT outcomes with relevant data columns
data_plot <- data.frame(TACIT = outcome$SCINA, (data[, colnames(Signature)[-1]]))
colnames(data_plot)[-1] <- colnames(Signature)[-1]

# Calculate median values for each cell type in TACIT
mean_values_TACIT <- data_plot %>%
  group_by(TACIT) %>%
  summarise_all(~quantile(., 0.5)) %>%
  as.data.frame()

# Remove the 19th row and the "Others" cell type
mean_values_TACIT <- mean_values_TACIT[-19,]
mean_values_TACIT <- mean_values_TACIT[which(mean_values_TACIT$TACIT != "Others"),]
rownames(mean_values_TACIT) <- mean_values_TACIT$TACIT
mean_values_TACIT <- mean_values_TACIT[, -1]

# Reorder columns based on signature order
ordered_index <- match(Signature_order[, 1], rownames(mean_values_TACIT))
mean_values_TACIT <- mean_values_TACIT[ordered_index, ]

# Define color breaks and palette for heatmap
my.breaks <- c(seq(-2, 0, by = 0.1), seq(0.1, 2, by = 0.1))
my.colors <- c(colorRampPalette(colors = c("blue", "white"))(length(my.breaks) / 2), 
               colorRampPalette(colors = c("white", "red"))(length(my.breaks) / 2))

# Original to short name mapping for markers
short_names <- c(
  "CD31 - vasculature:Cyc_19_ch_3" = "CD31",
  "CD34 - vasculature:Cyc_20_ch_3" = "CD34",
  "Cytokeratin - epithelia:Cyc_10_ch_2" = "Cytokeratin",
  "Ki67 - proliferation:Cyc_5_ch_4" = "Ki67",
  "p53 - tumor suppressor:Cyc_3_ch_3" = "p53",
  "EGFR - signaling:Cyc_13_ch_3" = "EGFR",
  "Granzyme B - cytotoxicity:Cyc_13_ch_2" = "GranzymeB",
  "beta-catenin - Wnt signaling:Cyc_4_ch_4" = "BetaCatenin",
  "PD-L1 - checkpoint:Cyc_5_ch_3" = "PDL1",
  "aSMA - smooth muscle:Cyc_11_ch_2" = "aSMA",
  "CD44 - stroma:Cyc_2_ch_2" = "CD44",
  "Vimentin - cytoplasm:Cyc_8_ch_2" = "Vimentin",
  "Podoplanin - lymphatics:Cyc_19_ch_4" = "Podoplanin",
  "CD45 - hematopoietic cells:Cyc_4_ch_2" = "CD45",
  "CD15 - granulocytes:Cyc_14_ch_2" = "CD15",
  "CD11b - macrophages:Cyc_10_ch_3" = "CD11b",
  "HLA-DR - MHC-II:Cyc_5_ch_2" = "HLA-DR",
  "CD20 - B cells:Cyc_8_ch_3" = "CD20",
  "CD11c - DCs:Cyc_12_ch_3" = "CD11c",
  "CD163 - macrophages:Cyc_17_ch_3" = "CD163",
  "CD68 - macrophages:Cyc_18_ch_4" = "CD68",
  "CD38 - multifunctional:Cyc_20_ch_4" = "CD38",
  "CD138 - plasma cells:Cyc_21_ch_3" = "CD138",
  "CD56 - NK cells:Cyc_10_ch_4" = "CD56",
  "CD57 - NK cells:Cyc_17_ch_4" = "CD57",
  "CD3 - T cells:Cyc_16_ch_4" = "CD3",
  "CD7 - T cells:Cyc_16_ch_3" = "CD7",
  "CD8 - cytotoxic T cells:Cyc_3_ch_2" = "CD8",
  "FOXP3 - regulatory T cells:Cyc_2_ch_3" = "FOXP3",
  "CD4 - T helper cells:Cyc_6_ch_3" = "CD4",
  "CD45RO - memory cells:Cyc_18_ch_3" = "CD45RO"
)

# Rename columns in mean_values_TACIT and Signature_order based on short names
colnames(mean_values_TACIT) <- short_names[colnames(mean_values_TACIT)]
colnames(Signature_order) <- short_names[colnames(Signature_order)]

# Scale mean values for heatmap visualization
aa <- scale(mean_values_TACIT)
rownames(aa) <- Signature_order[, 1]

# Generate heatmap
pheatmap(aa, cluster_cols = FALSE,
         cluster_rows = FALSE,
         show_colnames = TRUE,
         fontsize_col = 12,
         fontsize_row = 12)  # Adjust fontsize_row as needed to increase x-axis text size








#Louvain

# Combine TACIT outcomes with relevant data columns
data_plot <- data.frame(TACIT = outcome$Seurat_orginial, (data[, colnames(Signature)[-1]]))
colnames(data_plot)[-1] <- colnames(Signature)[-1]

# Calculate median values for each cell type in TACIT
mean_values_TACIT <- data_plot %>%
  group_by(TACIT) %>%
  summarise_all(~quantile(., 0.5)) %>%
  as.data.frame()

# Remove the 19th row and the "Others" cell type
mean_values_TACIT <- mean_values_TACIT[-19,]
mean_values_TACIT <- mean_values_TACIT[which(mean_values_TACIT$TACIT != "Others"),]
rownames(mean_values_TACIT) <- mean_values_TACIT$TACIT
mean_values_TACIT <- mean_values_TACIT[, -1]

# Reorder columns based on signature order
ordered_index <- match(Signature_order[, 1], rownames(mean_values_TACIT))
mean_values_TACIT <- mean_values_TACIT[ordered_index, ]

# Define color breaks and palette for heatmap
my.breaks <- c(seq(-2, 0, by = 0.1), seq(0.1, 2, by = 0.1))
my.colors <- c(colorRampPalette(colors = c("blue", "white"))(length(my.breaks) / 2), 
               colorRampPalette(colors = c("white", "red"))(length(my.breaks) / 2))

# Original to short name mapping for markers
short_names <- c(
  "CD31 - vasculature:Cyc_19_ch_3" = "CD31",
  "CD34 - vasculature:Cyc_20_ch_3" = "CD34",
  "Cytokeratin - epithelia:Cyc_10_ch_2" = "Cytokeratin",
  "Ki67 - proliferation:Cyc_5_ch_4" = "Ki67",
  "p53 - tumor suppressor:Cyc_3_ch_3" = "p53",
  "EGFR - signaling:Cyc_13_ch_3" = "EGFR",
  "Granzyme B - cytotoxicity:Cyc_13_ch_2" = "GranzymeB",
  "beta-catenin - Wnt signaling:Cyc_4_ch_4" = "BetaCatenin",
  "PD-L1 - checkpoint:Cyc_5_ch_3" = "PDL1",
  "aSMA - smooth muscle:Cyc_11_ch_2" = "aSMA",
  "CD44 - stroma:Cyc_2_ch_2" = "CD44",
  "Vimentin - cytoplasm:Cyc_8_ch_2" = "Vimentin",
  "Podoplanin - lymphatics:Cyc_19_ch_4" = "Podoplanin",
  "CD45 - hematopoietic cells:Cyc_4_ch_2" = "CD45",
  "CD15 - granulocytes:Cyc_14_ch_2" = "CD15",
  "CD11b - macrophages:Cyc_10_ch_3" = "CD11b",
  "HLA-DR - MHC-II:Cyc_5_ch_2" = "HLA-DR",
  "CD20 - B cells:Cyc_8_ch_3" = "CD20",
  "CD11c - DCs:Cyc_12_ch_3" = "CD11c",
  "CD163 - macrophages:Cyc_17_ch_3" = "CD163",
  "CD68 - macrophages:Cyc_18_ch_4" = "CD68",
  "CD38 - multifunctional:Cyc_20_ch_4" = "CD38",
  "CD138 - plasma cells:Cyc_21_ch_3" = "CD138",
  "CD56 - NK cells:Cyc_10_ch_4" = "CD56",
  "CD57 - NK cells:Cyc_17_ch_4" = "CD57",
  "CD3 - T cells:Cyc_16_ch_4" = "CD3",
  "CD7 - T cells:Cyc_16_ch_3" = "CD7",
  "CD8 - cytotoxic T cells:Cyc_3_ch_2" = "CD8",
  "FOXP3 - regulatory T cells:Cyc_2_ch_3" = "FOXP3",
  "CD4 - T helper cells:Cyc_6_ch_3" = "CD4",
  "CD45RO - memory cells:Cyc_18_ch_3" = "CD45RO"
)

# Rename columns in mean_values_TACIT and Signature_order based on short names
colnames(mean_values_TACIT) <- short_names[colnames(mean_values_TACIT)]
colnames(Signature_order) <- short_names[colnames(Signature_order)]

# Scale mean values for heatmap visualization
aa <- scale(mean_values_TACIT)
rownames(aa) <- Signature_order[, 1]

# Generate heatmap
pheatmap(aa, cluster_cols = FALSE,
         cluster_rows = FALSE,
         show_colnames = TRUE,
         fontsize_col = 12,
         fontsize_row = 12)  # Adjust fontsize_row as needed to increase x-axis text size






#reference

# Combine TACIT outcomes with relevant data columns
data_plot <- data.frame(TACIT = outcome$reference, (data[, colnames(Signature)[-1]]))
colnames(data_plot)[-1] <- colnames(Signature)[-1]

# Calculate median values for each cell type in TACIT
mean_values_TACIT <- data_plot %>%
  group_by(TACIT) %>%
  summarise_all(~quantile(., 0.5)) %>%
  as.data.frame()

# Remove the 19th row and the "Others" cell type
mean_values_TACIT <- mean_values_TACIT[-19,]
mean_values_TACIT <- mean_values_TACIT[which(mean_values_TACIT$TACIT != "Others"),]
rownames(mean_values_TACIT) <- mean_values_TACIT$TACIT
mean_values_TACIT <- mean_values_TACIT[, -1]

# Reorder columns based on signature order
ordered_index <- match(Signature_order[, 1], rownames(mean_values_TACIT))
mean_values_TACIT <- mean_values_TACIT[ordered_index, ]

# Define color breaks and palette for heatmap
my.breaks <- c(seq(-2, 0, by = 0.1), seq(0.1, 2, by = 0.1))
my.colors <- c(colorRampPalette(colors = c("blue", "white"))(length(my.breaks) / 2), 
               colorRampPalette(colors = c("white", "red"))(length(my.breaks) / 2))

# Original to short name mapping for markers
short_names <- c(
  "CD31 - vasculature:Cyc_19_ch_3" = "CD31",
  "CD34 - vasculature:Cyc_20_ch_3" = "CD34",
  "Cytokeratin - epithelia:Cyc_10_ch_2" = "Cytokeratin",
  "Ki67 - proliferation:Cyc_5_ch_4" = "Ki67",
  "p53 - tumor suppressor:Cyc_3_ch_3" = "p53",
  "EGFR - signaling:Cyc_13_ch_3" = "EGFR",
  "Granzyme B - cytotoxicity:Cyc_13_ch_2" = "GranzymeB",
  "beta-catenin - Wnt signaling:Cyc_4_ch_4" = "BetaCatenin",
  "PD-L1 - checkpoint:Cyc_5_ch_3" = "PDL1",
  "aSMA - smooth muscle:Cyc_11_ch_2" = "aSMA",
  "CD44 - stroma:Cyc_2_ch_2" = "CD44",
  "Vimentin - cytoplasm:Cyc_8_ch_2" = "Vimentin",
  "Podoplanin - lymphatics:Cyc_19_ch_4" = "Podoplanin",
  "CD45 - hematopoietic cells:Cyc_4_ch_2" = "CD45",
  "CD15 - granulocytes:Cyc_14_ch_2" = "CD15",
  "CD11b - macrophages:Cyc_10_ch_3" = "CD11b",
  "HLA-DR - MHC-II:Cyc_5_ch_2" = "HLA-DR",
  "CD20 - B cells:Cyc_8_ch_3" = "CD20",
  "CD11c - DCs:Cyc_12_ch_3" = "CD11c",
  "CD163 - macrophages:Cyc_17_ch_3" = "CD163",
  "CD68 - macrophages:Cyc_18_ch_4" = "CD68",
  "CD38 - multifunctional:Cyc_20_ch_4" = "CD38",
  "CD138 - plasma cells:Cyc_21_ch_3" = "CD138",
  "CD56 - NK cells:Cyc_10_ch_4" = "CD56",
  "CD57 - NK cells:Cyc_17_ch_4" = "CD57",
  "CD3 - T cells:Cyc_16_ch_4" = "CD3",
  "CD7 - T cells:Cyc_16_ch_3" = "CD7",
  "CD8 - cytotoxic T cells:Cyc_3_ch_2" = "CD8",
  "FOXP3 - regulatory T cells:Cyc_2_ch_3" = "FOXP3",
  "CD4 - T helper cells:Cyc_6_ch_3" = "CD4",
  "CD45RO - memory cells:Cyc_18_ch_3" = "CD45RO"
)

# Rename columns in mean_values_TACIT and Signature_order based on short names
colnames(mean_values_TACIT) <- short_names[colnames(mean_values_TACIT)]
colnames(Signature_order) <- short_names[colnames(Signature_order)]

# Scale mean values for heatmap visualization
aa <- scale(mean_values_TACIT)
rownames(aa) <- Signature_order[, 1]

# Generate heatmap
pheatmap(aa, cluster_cols = FALSE,
         cluster_rows = FALSE,
         show_colnames = TRUE,
         fontsize_col = 12,
         fontsize_row = 12)  # Adjust fontsize_row as needed to increase x-axis text size




#CELESTA

# Combine TACIT outcomes with relevant data columns
data_plot <- data.frame(TACIT = outcome$CELESTA, (data[, colnames(Signature)[-1]]))
colnames(data_plot)[-1] <- colnames(Signature)[-1]

# Calculate median values for each cell type in TACIT
mean_values_TACIT <- data_plot %>%
  group_by(TACIT) %>%
  summarise_all(~quantile(., 0.5)) %>%
  as.data.frame()

# Remove the 19th row and the "Others" cell type
mean_values_TACIT <- mean_values_TACIT[-19,]
mean_values_TACIT <- mean_values_TACIT[which(mean_values_TACIT$TACIT != "Others"),]
rownames(mean_values_TACIT) <- mean_values_TACIT$TACIT
mean_values_TACIT <- mean_values_TACIT[, -1]

# Reorder columns based on signature order
ordered_index <- match(Signature_order[, 1], rownames(mean_values_TACIT))
mean_values_TACIT <- mean_values_TACIT[ordered_index, ]

# Define color breaks and palette for heatmap
my.breaks <- c(seq(-2, 0, by = 0.1), seq(0.1, 2, by = 0.1))
my.colors <- c(colorRampPalette(colors = c("blue", "white"))(length(my.breaks) / 2), 
               colorRampPalette(colors = c("white", "red"))(length(my.breaks) / 2))

# Original to short name mapping for markers
short_names <- c(
  "CD31 - vasculature:Cyc_19_ch_3" = "CD31",
  "CD34 - vasculature:Cyc_20_ch_3" = "CD34",
  "Cytokeratin - epithelia:Cyc_10_ch_2" = "Cytokeratin",
  "Ki67 - proliferation:Cyc_5_ch_4" = "Ki67",
  "p53 - tumor suppressor:Cyc_3_ch_3" = "p53",
  "EGFR - signaling:Cyc_13_ch_3" = "EGFR",
  "Granzyme B - cytotoxicity:Cyc_13_ch_2" = "GranzymeB",
  "beta-catenin - Wnt signaling:Cyc_4_ch_4" = "BetaCatenin",
  "PD-L1 - checkpoint:Cyc_5_ch_3" = "PDL1",
  "aSMA - smooth muscle:Cyc_11_ch_2" = "aSMA",
  "CD44 - stroma:Cyc_2_ch_2" = "CD44",
  "Vimentin - cytoplasm:Cyc_8_ch_2" = "Vimentin",
  "Podoplanin - lymphatics:Cyc_19_ch_4" = "Podoplanin",
  "CD45 - hematopoietic cells:Cyc_4_ch_2" = "CD45",
  "CD15 - granulocytes:Cyc_14_ch_2" = "CD15",
  "CD11b - macrophages:Cyc_10_ch_3" = "CD11b",
  "HLA-DR - MHC-II:Cyc_5_ch_2" = "HLA-DR",
  "CD20 - B cells:Cyc_8_ch_3" = "CD20",
  "CD11c - DCs:Cyc_12_ch_3" = "CD11c",
  "CD163 - macrophages:Cyc_17_ch_3" = "CD163",
  "CD68 - macrophages:Cyc_18_ch_4" = "CD68",
  "CD38 - multifunctional:Cyc_20_ch_4" = "CD38",
  "CD138 - plasma cells:Cyc_21_ch_3" = "CD138",
  "CD56 - NK cells:Cyc_10_ch_4" = "CD56",
  "CD57 - NK cells:Cyc_17_ch_4" = "CD57",
  "CD3 - T cells:Cyc_16_ch_4" = "CD3",
  "CD7 - T cells:Cyc_16_ch_3" = "CD7",
  "CD8 - cytotoxic T cells:Cyc_3_ch_2" = "CD8",
  "FOXP3 - regulatory T cells:Cyc_2_ch_3" = "FOXP3",
  "CD4 - T helper cells:Cyc_6_ch_3" = "CD4",
  "CD45RO - memory cells:Cyc_18_ch_3" = "CD45RO"
)

# Rename columns in mean_values_TACIT and Signature_order based on short names
colnames(mean_values_TACIT) <- short_names[colnames(mean_values_TACIT)]
colnames(Signature_order) <- short_names[colnames(Signature_order)]

# Scale mean values for heatmap visualization
aa <- scale(mean_values_TACIT)
rownames(aa) <- Signature_order[, 1]

# Generate heatmap
pheatmap(aa, cluster_cols = FALSE,
         cluster_rows = FALSE,
         show_colnames = TRUE,
         fontsize_col = 12,
         fontsize_row = 12)  # Adjust fontsize_row as needed to increase x-axis text size







#--------------------------
#g------------------------
#--------------------------



data <- data.table::fread("~/Documents/TACIT figure/dataset/Figure 2/23_09_CODEX_HuBMAP_alldata_Dryad_merged.csv")
data=as.data.frame(data)
Signature <- read_excel("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Clean_up_disseration/data/BE data/test_BE_new.xlsx")
Signature=as.data.frame(Signature)
colnames(Signature)[1]="cell_type"
Signature[is.na(Signature)==T]=0 ##assign all NA as 0
Signature=as.data.frame(Signature)
Signature=Signature[,-2]

outcome <- readr::read_csv("~/Documents/TACIT figure/dataset/Figure 2/HI_outcome.csv")









#--------------------------
#g------------------------
#--------------------------

TACIT_val=calculate_metrics(outcome$TACIT[which(outcome$TACIT!="Others")],outcome$reference[which(outcome$TACIT!="Others")])
Louvain_val002=calculate_metrics(outcome$Louvain[which(outcome$Louvain!="Others")],outcome$reference[which(outcome$Louvain!="Others")])

# Assuming you have the weighted recall, precision, and F1 scores for each method
data_plot <- data.frame(
  Method = c("TACIT","Louvain"),
  WeightedRecall = c(TACIT_val$WeightedRecall,Louvain_val002$WeightedRecall),
  WeightedPrecision = c(TACIT_val$WeightedPrecision,Louvain_val002$WeightedPrecision),
  WeightedF1 = c(TACIT_val$WeightedF1,Louvain_val002$WeightedF1)
)


colnames(data_plot)=c("Method","Weighted Recall","Weighted Precision","Weighted F1")
# Transform data from wide to long format
long_data <- pivot_longer(data_plot, cols = -Method, names_to = "Metric", values_to = "Score")

long_data$Method=factor(long_data$Method,levels = c("TACIT","CELESTA","SCINA","Louvain"))
long_data$Metric=factor(long_data$Metric,levels = c("Weighted Recall","Weighted Precision","Weighted F1"))
# Plot
library(ggplot2)
ggplot(long_data, aes(x = Method, y = Score, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.7) +
  geom_text(aes(label = sprintf("%.2f", Score)), 
            position = position_dodge(width = 0.9), vjust = -0.25, color = "black", size = 6) +
  facet_wrap(~ Metric, scales = "free_y") +
  labs(x = "Cell Type", y = "Score", title = "") +
  scale_fill_manual(values = c("TACIT" = "red", "CELESTA" = "purple", "SCINA" = "green", "Louvain" = "blue")) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.6),
        legend.position = "right")+ylim(0,1)













#--------------------------
#h------------------------
#--------------------------


# Create a dataframe with the necessary columns
data_plot <- data.frame(
  ref = outcome$reference,
  Louvain = outcome$Louvain,
  TACIT = outcome$TACIT,
  TMA = data$unique_region
)

# Convert the vector to a dataframe for easier manipulation
df <- data.frame(cell_type = data_plot$ref)

# Calculate the proportion of each cell type
proportions <- df %>%
  group_by(cell_type) %>%
  summarise(count = n()) %>%
  mutate(proportion = count / sum(count)) %>%
  filter(proportion < 0.01)



data_plot$ref=ifelse(data_plot$ref%in%c(proportions$cell_type),"Rare cell type","Not rare")
data_plot$TACIT=ifelse(data_plot$TACIT%in%c(proportions$cell_type),"Rare cell type","Not rare")
data_plot$Louvain=ifelse(data_plot$Louvain%in%c(proportions$cell_type),"Rare cell type","Not rare")




# Assuming your dataframe is named data_plot
# Convert factors to character to avoid drop of unused levels in subsequent filtering
data_plot <- data_plot %>% mutate(across(c(ref, TACIT,Louvain), as.character))

# Summarizing counts for Immune cells (all)
summary_data <- data_plot %>% 
  group_by(TMA) %>% 
  summarise(Ref = sum(ref == "Rare cell type"),
            TACIT = sum(TACIT == "Rare cell type"),
            Louvain = sum(Louvain == "Rare cell type"))


method_name="TACIT"

# Calculating correlation
correlation <- round(cor(summary_data$Ref, summary_data[[method_name]], use = "complete.obs"), 2)

# Plot the correlation
ggplot(summary_data, aes(y = Ref, x = summary_data[[method_name]])) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red",size=2) + # x=y line
  theme_classic(base_size = 30) +
  labs(title = "",
       y = "Number of cells in Reference",
       x = paste("Number of cells in", method_name))  +
  theme(plot.title = element_text(size = 14)) +
  annotate("text", x = max(summary_data$Ref, na.rm = TRUE) * 0.2, y = max(summary_data[[method_name]], na.rm = TRUE), label = paste("R =", correlation), size = 1)+ geom_smooth(method = "lm", color = "blue", se = FALSE,size=2)





# Summarizing counts for Immune cells (all)
summary_data <- data_plot %>% 
  group_by(TMA) %>% 
  summarise(Ref = sum(ref != "Rare cell type"),
            TACIT = sum(TACIT != "Rare cell type"),
            Louvain = sum(Louvain != "Rare cell type"))


method_name="TACIT"

# Calculating correlation
correlation <- round(cor(summary_data$Ref, summary_data[[method_name]], use = "complete.obs"), 2)

# Plot the correlation
ggplot(summary_data, aes(y = Ref, x = summary_data[[method_name]])) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red",size=2) + # x=y line
  theme_classic(base_size = 30) +
  labs(title = "",
       y = "Number of cells in Reference",
       x = paste("Number of cells in", method_name))  +
  theme(plot.title = element_text(size = 14)) +
  annotate("text", x = max(summary_data$Ref, na.rm = TRUE) * 0.2, y = max(summary_data[[method_name]], na.rm = TRUE), label = paste("R =", correlation), size = 1)+ geom_smooth(method = "lm", color = "blue", se = FALSE,size=2)







#--------------------------
#i------------------------
#--------------------------

#-TACIT

# Create a data frame combining TACIT outcomes and signature data
data_plot <- data.frame(TACIT = outcome$TACIT, (data[, colnames(Signature)[-1]]))
colnames(data_plot)[-1] <- colnames(Signature)[-1]

# Prepare the Signature data
colnames(Signature)[1] <- "cell_type"
Signature[is.na(Signature)] <- 0

# Calculate the median values for each cell type in TACIT
mean_values_TACIT <- data_plot %>%
  group_by(TACIT) %>%
  summarise_all(~quantile(., 0.5)) %>%
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
    
    testt <- log2(mean(as.numeric(data_anb[which(outcome$TACIT == ct_choose), marker_choose])-min(as.numeric(data_anb[, marker_choose]))) /
                    mean(as.numeric(data_anb[which(outcome$TACIT %in% select_ct_0), marker_choose])-min(as.numeric(data_anb[, marker_choose]))))
    
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



#-Louvain

# Create a data frame combining TACIT outcomes and signature data
data_plot <- data.frame(TACIT = outcome$Louvain, (data[, colnames(Signature)[-1]]))
colnames(data_plot)[-1] <- colnames(Signature)[-1]

# Prepare the Signature data
colnames(Signature)[1] <- "cell_type"
Signature[is.na(Signature)] <- 0

# Calculate the median values for each cell type in TACIT
mean_values_TACIT <- data_plot %>%
  group_by(TACIT) %>%
  summarise_all(~quantile(., 0.5)) %>%
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
    
    testt <- log2(mean(as.numeric(data_anb[which(outcome$Louvain == ct_choose), marker_choose])-min(as.numeric(data_anb[, marker_choose]))) /
                    mean(as.numeric(data_anb[which(outcome$Louvain %in% select_ct_0), marker_choose])-min(as.numeric(data_anb[, marker_choose]))))
    
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
data_final_Louvain <- data_final2










#-Reference

# Create a data frame combining TACIT outcomes and signature data
data_plot <- data.frame(TACIT = outcome$reference, (data[, colnames(Signature)[-1]]))
colnames(data_plot)[-1] <- colnames(Signature)[-1]

# Prepare the Signature data
colnames(Signature)[1] <- "cell_type"
Signature[is.na(Signature)] <- 0

# Calculate the median values for each cell type in TACIT
mean_values_TACIT <- data_plot %>%
  group_by(TACIT) %>%
  summarise_all(~quantile(., 0.5)) %>%
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
    
    testt <- log2(mean(as.numeric(data_anb[which(outcome$reference == ct_choose), marker_choose])-min(as.numeric(data_anb[, marker_choose]))) /
                    mean(as.numeric(data_anb[which(outcome$reference %in% select_ct_0), marker_choose])-min(as.numeric(data_anb[, marker_choose]))))
    
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
data_final_reference <- data_final2




data_final_reference$Methods="Reference"
data_final_Louvain$Methods="Louvain"
data_final_TACIT$Methods="TACIT"


data_final_plot=rbind(data_final_reference,data_final_Louvain,data_final_TACIT)

data_final_plot$Group_plot=ifelse(data_final_plot$Group=="Significant & Signature",data_final_plot$Methods,"Others")

data_final_plot_sub=data_final_plot[which(data_final_plot$Signature==1),]

# Filter the data to include only the specified methods
filtered_data <- data_final_plot_sub %>%
  filter(Methods %in% c("TACIT", "Louvain","Reference"))

filtered_data$Methods=factor(filtered_data$Methods,level=c("Reference",
                                                           "TACIT","Louvain"))

#logFC
ggplot(filtered_data, aes(x = Methods, y = logFC, fill = Methods)) +
  geom_boxplot(alpha = 1, outlier.size = 2, position = position_dodge(0.8)) +
  scale_fill_manual(values = c("TACIT" = "red", "CELESTA" = "purple", "SCINA" = "green", "Louvain" = "blue", "Reference" = "yellow")) +
  theme_classic(base_size = 20) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  ylab("") +
  xlab("") +
  geom_signif(comparisons = list(  c("TACIT", "Reference"), c("TACIT", "Louvain")),
              map_signif_level=TRUE,
              textsize = 3, vjust = -0.2)


#-log10(padjust)
ggplot(filtered_data, aes(x = Methods, y = -log10(padjust), fill = Methods)) +
  geom_boxplot(alpha = 1, outlier.size = 2, position = position_dodge(0.8)) +
  scale_fill_manual(values = c("TACIT" = "red", "CELESTA" = "purple", "SCINA" = "green", "Louvain" = "blue", "Reference" = "yellow")) +
  theme_classic(base_size = 20) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  ylab("") +
  xlab("") +
  geom_signif(comparisons = list(  c("TACIT", "Reference"), c("TACIT", "Louvain")),
              map_signif_level=TRUE,
              textsize = 3, vjust = -0.2)












#--------------------------
#j------------------------
#--------------------------


# Subset data for the specific region "B009_Sigmoid"
data_sub <- data[which(data$unique_region == "B009_Sigmoid"), ]

# Plot TACIT predictions
df_tacit <- data.frame(x = data_sub$x, y = data_sub$y, predicted = outcome$TACIT[which(data$unique_region == "B009_Sigmoid")])
ggplot(df_tacit, aes(x = x, y = y, color = predicted)) +
  geom_point(size = 0.1) +
  theme_classic(base_size = 15) +
  theme(legend.position = "none") +
  labs(x = "X", y = "Y", title = "") +
  guides(color = guide_legend(override.aes = list(size = 3), title = "Cell type")) +
  scale_color_manual(values = c(
    "B" = "tan", "Endothelial" = "darkred", "CD7+ Immune" = "pink",
    "CD4+ T cell" = "violet", "DC" = "darkviolet",
    "M1 Macrophage" = "green", "CD8+ T" = "darkgreen",
    "Cycling TA" = "blue", "Lymphatic" = "skyblue", "NK" = "yellow",
    "Others" = "grey90", "Innate" = "darkblue", "Plasma" = "orange", "Stroma" = "gold",
    "Neutrophil" = "black", "Enterocyte" = "red", "vasculature" = "darkcyan", "Smooth muscle" = "tan4", "TA" = "red3", "Neuroendocrine" = "pink2", "Goblet" = "cyan", "M2 Macrophage" = "green3", "ICC" = "cyan3", "Nerve" = "yellow3",
    "Paneth" = "gold3"
  )) +
  ggtitle("") +
  theme(
    axis.line = element_blank(),  # Remove axis lines
    axis.ticks = element_blank(),  # Remove tick marks
    axis.text.x = element_blank(),  # Remove x axis text
    axis.text.y = element_blank()   # Remove y axis text
  ) + ylab("") + xlab("")

# Plot Reference outcomes
df_reference <- data.frame(x = data_sub$x, y = data_sub$y, predicted = outcome$reference[which(data$unique_region == "B009_Sigmoid")])
ggplot(df_reference, aes(x = x, y = y, color = predicted)) +
  geom_point(size = 0.1) +
  theme_classic(base_size = 15) +
  theme(legend.position = "none") +
  labs(x = "X", y = "Y", title = "") +
  guides(color = guide_legend(override.aes = list(size = 3), title = "Cell type")) +
  scale_color_manual(values = c(
    "B" = "tan", "Endothelial" = "darkred", "CD7+ Immune" = "pink",
    "CD4+ T cell" = "violet", "DC" = "darkviolet",
    "M1 Macrophage" = "green", "CD8+ T" = "darkgreen",
    "Cycling TA" = "blue", "Lymphatic" = "skyblue", "NK" = "yellow",
    "Others" = "grey90", "Innate" = "darkblue", "Plasma" = "orange", "Stroma" = "gold",
    "Neutrophil" = "black", "Enterocyte" = "red", "vasculature" = "darkcyan", "Smooth muscle" = "tan4", "TA" = "red3", "Neuroendocrine" = "pink2", "Goblet" = "cyan", "M2 Macrophage" = "green3", "ICC" = "cyan3", "Nerve" = "yellow3",
    "Paneth" = "gold3"
  )) +
  ggtitle("") +
  theme(
    axis.line = element_blank(),  # Remove axis lines
    axis.ticks = element_blank(),  # Remove tick marks
    axis.text.x = element_blank(),  # Remove x axis text
    axis.text.y = element_blank()   # Remove y axis text
  ) + ylab("") + xlab("")

# Plot Louvain predictions
df_louvain <- data.frame(x = data_sub$x, y = data_sub$y, predicted = outcome$Louvain[which(data$unique_region == "B009_Sigmoid")])
ggplot(df_louvain, aes(x = x, y = y, color = predicted)) +
  geom_point(size = 0.1) +
  theme_classic(base_size = 15) +
  theme(legend.position = "none") +
  labs(x = "X", y = "Y", title = "") +
  guides(color = guide_legend(override.aes = list(size = 3), title = "Cell type")) +
  scale_color_manual(values = c(
    "B" = "tan", "Endothelial" = "darkred", "CD7+ Immune" = "pink",
    "CD4+ T cell" = "violet", "DC" = "darkviolet",
    "M1 Macrophage" = "green", "CD8+ T" = "darkgreen",
    "Cycling TA" = "blue", "Lymphatic" = "skyblue", "NK" = "yellow",
    "Others" = "grey90", "Innate" = "darkblue", "Plasma" = "orange", "Stroma" = "gold",
    "Neutrophil" = "black", "Enterocyte" = "red", "vasculature" = "darkcyan", "Smooth muscle" = "tan4", "TA" = "red3", "Neuroendocrine" = "pink2", "Goblet" = "cyan", "M2 Macrophage" = "green3", "ICC" = "cyan3", "Nerve" = "yellow3",
    "Paneth" = "gold3"
  )) +
  ggtitle("") +
  theme(
    axis.line = element_blank(),  # Remove axis lines
    axis.ticks = element_blank(),  # Remove tick marks
    axis.text.x = element_blank(),  # Remove x axis text
    axis.text.y = element_blank()   # Remove y axis text
  ) + ylab("") + xlab("")


#--------------------------
#k------------------------
#--------------------------

# Perform UMAP dimensionality reduction on the data, excluding the first column (assumed to be non-numeric)
data_umap <- umap::umap(data[, colnames(Signature)[-1]], scale = TRUE, method = "umap-learn", metric = "correlation")
data_umap <- data_umap$layout

# Convert UMAP results to a data frame and set column names
colnames(data_umap) <- c("UMAP1", "UMAP2")
data_umap <- as.data.frame(data_umap)

# Define a common color palette for cell types
cell_type_colors <- c("B" = "tan", "Endothelial" = "darkred", "CD7+ Immune" = "pink",
                      "CD4+ T cell" = "violet", "DC" = "darkviolet",
                      "M1 Macrophage" = "green", "CD8+ T" = "darkgreen",
                      "Cycling TA" = "blue", "Lymphatic" = "skyblue", "NK" = "yellow",
                      "Others" = "grey90", "Innate" = "darkblue", "Plasma" = "orange", 
                      "Stroma" = "gold", "Neutrophil" = "black", "Enterocyte" = "red", 
                      "vasculature" = "darkcyan", "Smooth muscle" = "tan4", "TA" = "red3", 
                      "Neuroendocrine" = "pink2", "Goblet" = "cyan", "M2 Macrophage" = "green3", 
                      "ICC" = "cyan3", "Nerve" = "yellow3", "Paneth" = "gold3")

# Plot UMAP with TACIT outcomes
ggplot(data_umap, aes(x = UMAP1, y = UMAP2, color = outcome$TACIT)) +
  geom_point(size = 0.1) +
  scale_color_manual(values = cell_type_colors) +
  theme_classic(base_size = 15) +
  guides(color = guide_legend(override.aes = list(size = 3), title = "Cell type")) +
  ggtitle("") +
  theme(axis.line = element_blank(),         # Remove axis lines
        axis.ticks = element_blank(),        # Remove tick marks
        axis.text.x = element_blank(),       # Remove x axis text
        axis.text.y = element_blank(),       # Remove y axis text
        legend.position = "none") +          # Remove legend
  labs(x = "", y = "")

# Plot UMAP with reference outcomes
ggplot(data_umap, aes(x = UMAP1, y = UMAP2, color = outcome$reference)) +
  geom_point(size = 0.1) +
  scale_color_manual(values = cell_type_colors) +
  theme_classic(base_size = 15) +
  guides(color = guide_legend(override.aes = list(size = 3), title = "Cell type")) +
  ggtitle("") +
  theme(axis.line = element_blank(),         # Remove axis lines
        axis.ticks = element_blank(),        # Remove tick marks
        axis.text.x = element_blank(),       # Remove x axis text
        axis.text.y = element_blank(),       # Remove y axis text
        legend.position = "none") +          # Remove legend
  labs(x = "", y = "")

# Plot UMAP with Louvain outcomes
ggplot(data_umap, aes(x = UMAP1, y = UMAP2, color = outcome$Louvain)) +
  geom_point(size = 0.1) +
  scale_color_manual(values = cell_type_colors) +
  theme_classic(base_size = 15) +
  guides(color = guide_legend(override.aes = list(size = 3), title = "Cell type")) +
  ggtitle("") +
  theme(axis.line = element_blank(),         # Remove axis lines
        axis.ticks = element_blank(),        # Remove tick marks
        axis.text.x = element_blank(),       # Remove x axis text
        axis.text.y = element_blank(),       # Remove y axis text
        legend.position = "none") +          # Remove legend
  labs(x = "", y = "")





#--------------------------
#i------------------------
#--------------------------


# Reference Heatmap
# Group the data by the predicted label and calculate the median for each column
data_plot <- data.frame(TACIT = outcome$reference, (data[, colnames(Signature)[-1]]))
colnames(data_plot)[-1] <- colnames(Signature)[-1]

# Calculate median values for each cell type in TACIT
mean_values_TACIT <- data_plot %>%
  group_by(TACIT) %>%
  summarise_all(~quantile(., 0.5)) %>%
  as.data.frame()
rownames(mean_values_TACIT) <- mean_values_TACIT$TACIT
mean_values_TACIT <- mean_values_TACIT[which(mean_values_TACIT$TACIT != "Others"), ]
mean_values_TACIT <- mean_values_TACIT[, -1]

# Reorder columns based on signature order
Signature_order <- reorder_columns_based_on_signature(Signature)
ordered_index <- match(Signature_order$cell_type, rownames(mean_values_TACIT))
mean_values_TACIT <- mean_values_TACIT[ordered_index, ]

# Define color breaks and palette for heatmap
my.breaks <- c(seq(-2, 0, by = 0.1), seq(0.1, 2, by = 0.1))
my.colors <- c(colorRampPalette(colors = c("blue", "white"))(length(my.breaks) / 2), 
               colorRampPalette(colors = c("white", "red"))(length(my.breaks) / 2))

# Generate heatmap for Reference data
aa <- scale(mean_values_TACIT[, colnames(Signature_order)[-1]])
pheatmap(aa, cluster_cols = FALSE, cluster_rows = FALSE, show_colnames = FALSE, show_rownames = FALSE,
         fontsize_row = 15, legend = FALSE, fontsize_col = 15)

# TACIT Heatmap
# Group the data by the predicted label and calculate the median for each column
data_plot <- data.frame(TACIT = outcome$TACIT, (data[, colnames(Signature)[-1]]))
colnames(data_plot)[-1] <- colnames(Signature)[-1]

# Calculate median values for each cell type in TACIT
mean_values_TACIT <- data_plot %>%
  group_by(TACIT) %>%
  summarise_all(~quantile(., 0.5)) %>%
  as.data.frame()
rownames(mean_values_TACIT) <- mean_values_TACIT$TACIT
mean_values_TACIT <- mean_values_TACIT[which(mean_values_TACIT$TACIT != "Others"), ]
mean_values_TACIT <- mean_values_TACIT[, -1]

# Reorder columns based on signature order
Signature_order <- reorder_columns_based_on_signature(Signature)
ordered_index <- match(Signature_order$cell_type, rownames(mean_values_TACIT))
mean_values_TACIT <- mean_values_TACIT[ordered_index, ]

# Generate heatmap for TACIT data
aa <- scale(mean_values_TACIT[, colnames(Signature_order)[-1]])
pheatmap(aa, cluster_cols = FALSE, cluster_rows = FALSE, show_colnames = FALSE, show_rownames = FALSE,
         fontsize_row = 15, legend = FALSE, fontsize_col = 15)

# Louvain Heatmap
# Group the data by the predicted label and calculate the median for each column
data_plot <- data.frame(TACIT = outcome$Louvain, (data[, colnames(Signature)[-1]]))
colnames(data_plot)[-1] <- colnames(Signature)[-1]

# Calculate median values for each cell type in TACIT
mean_values_TACIT <- data_plot %>%
  group_by(TACIT) %>%
  summarise_all(~quantile(., 0.5)) %>%
  as.data.frame()
rownames(mean_values_TACIT) <- mean_values_TACIT$TACIT
mean_values_TACIT <- mean_values_TACIT[which(mean_values_TACIT$TACIT != "Others"), ]
mean_values_TACIT <- mean_values_TACIT[, -1]

# Reorder columns based on signature order
Signature_order <- reorder_columns_based_on_signature(Signature)
ordered_index <- match(Signature_order$cell_type, rownames(mean_values_TACIT))
mean_values_TACIT <- mean_values_TACIT[ordered_index, ]

# Generate heatmap for Louvain data
aa <- scale(mean_values_TACIT[, colnames(Signature_order)[-1]])
pheatmap(aa, cluster_cols = FALSE, cluster_rows = FALSE, show_colnames = FALSE, show_rownames = FALSE,
         fontsize_row = 15, legend = FALSE, fontsize_col = 15)


