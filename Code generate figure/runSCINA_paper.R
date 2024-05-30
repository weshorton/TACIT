#-----------------------------------------------------------------------------
###Goal: Obtain SCINA: Semi-supervised Category Identification and Assignment
#-----------------------------------------------------------------------------
###Cell type annotation Step: 
#### 1) Load the data
#### 2) Load signature matrix
#### 3) Quality Control
#### 4) Run SCINA
#-----------------------------------------------------------------------------

library(readr)
library(readxl)
library(dplyr)
library(tidyverse)
library('SCINA')

#-----------------------------------------------------------------------------

###Colon datasets (Figure 2 top)

#-----------------------------------------------------------------------------



# 1) Load the data
# Reading primary dataset for analysis
data <- readr::read_csv("CRC_DataSet.csv")

# 2) Load signature matrix
# Loading a CSV file that contains the signature matrix used for cell type identification
Signature <- readr::read_csv("Signature_CRC_TACIT_SCINA_Louvain.csv")

# 3) Quality Control
# Renaming certain cell types for consistency in nomenclature
data$ClusterName <- ifelse(data$ClusterName %in% c("CD11b+ CD68+ macrophages", "CD163+ macrophages",
                                                   "CD68+ macrophages", "CD68+ macrophages GzmB+",
                                                   "CD68+ CD163+ macrophages"), "macrophages", data$ClusterName)
data$ClusterName <- ifelse(data$ClusterName == "CD4+ T cells GATA3+", "CD4+ T cells", data$ClusterName)

# Filtering the data to exclude unwanted cell types based on their cluster names
data <- data[data$ClusterName %in% Signature$cell_type,]

# Print the dimensions of the filtered data to check how many rows and columns remain
print(dim(data))

# Replace all NA values in the signature matrix with 0 and convert to a dataframe
Signature[is.na(Signature)] <- 0 
Signature <- as.data.frame(Signature)

# Spatial validation: Visualizing the distribution of a specific marker in the dataset
marker <- "`CD31 - vasculature:Cyc_19_ch_3`"
ggplot(data, aes(x = X, y = Y)) +
  geom_point(aes_string(color = marker), size = 0.1) +
  theme_classic(base_size = 15) +
  ggtitle(marker) +
  scale_color_gradient(low = "grey90", high = "darkblue", limits = c(0, 100), na.value = "darkblue")

# 4) Run SCINA
# Setting up expressions by cell type for SCINA analysis
list_expression_by_cell_type <- list(
  "vasculature" = colnames(Signature[1, Signature[1, ] == 1 | is.na(Signature[1, ])]),
  "tumor cells" = colnames(Signature[2, Signature[2, ] == 1 | is.na(Signature[2, ])]),
  "stroma" = colnames(Signature[4, Signature[4, ] == 1 | is.na(Signature[4, ])]),
  "immune cells" = "CD45 - hematopoietic cells:Cyc_4_ch_2",
  "smooth muscle" = "aSMA - smooth muscle:Cyc_11_ch_2"
)

# Prepare expression data by removing the first column of Signature and apply log transformation
data_expression <- data[, colnames(Signature)[-1]]
exp_data_expression <- log(data_expression + 1)

# Execute SCINA to deduce cell types based on expression profiles
results <- SCINA(t(exp_data_expression), list_expression_by_cell_type, max_iter = 600, convergence_n = 30, 
                 convergence_rate = 0.999, sensitivity_cutoff = 0.9, rm_overlap = FALSE, allow_unknown = TRUE)

# Store the results in a dataframe along with the original IDs
data_result1 <- data.frame(cell_labels = results$cell_labels, ID = data$CellID)

# Identify immune cells using the initial SCINA results
# Create a new dataframe by filtering for immune cells or unknown cell types
exp_data_expression_input2 <- data.frame(exp_data_expression, label = data_result1$results.cell_labels)
exp_data_expression_input2 <- exp_data_expression_input2[which(exp_data_expression_input2$label %in% c("immune cells", "unknown")),]
colnames(exp_data_expression_input2)[1:31] <- colnames(exp_data_expression)

# Define expression signatures for different immune cell types excluding common hematopoietic marker
list_expression_by_cell_type <- list(
  "lymphatics" = setdiff(colnames(Signature[5, which(Signature[5,] == 1)]), "CD45 - hematopoietic cells:Cyc_4_ch_2"),
  "CD3+ T cells" = setdiff(colnames(Signature[13, which(Signature[13,] == 1)]), "CD45 - hematopoietic cells:Cyc_4_ch_2"),
  "granulocytes" = setdiff(colnames(Signature[7, which(Signature[7,] == 1)]), "CD45 - hematopoietic cells:Cyc_4_ch_2"),
  "B cells" = setdiff(colnames(Signature[8, which(Signature[8,] == 1)]), "CD45 - hematopoietic cells:Cyc_4_ch_2"),
  "CD11c+ DCs" = setdiff(colnames(Signature[9, which(Signature[9,] == 1)]), "CD45 - hematopoietic cells:Cyc_4_ch_2"),
  "macrophages" = setdiff(colnames(Signature[10, which(Signature[10,] == 1)]), "CD45 - hematopoietic cells:Cyc_4_ch_2"),
  "plasma cells" = setdiff(colnames(Signature[11, which(Signature[11,] == 1)]), "CD45 - hematopoietic cells:Cyc_4_ch_2"),
  "NK cells" = setdiff(colnames(Signature[12, which(Signature[12,] == 1)]), "CD45 - hematopoietic cells:Cyc_4_ch_2")
)

# Execute SCINA to further classify the filtered immune cells
results <- SCINA(t(exp_data_expression_input2[,-32]), list_expression_by_cell_type, max_iter = 100, convergence_n = 10, 
                 convergence_rate = 0.999, sensitivity_cutoff = 0.9, rm_overlap = FALSE, allow_unknown = TRUE)

# Display the table of SCINA classification results for immune cells
table(results$cell_labels)

# Prepare the results with IDs for merging later
data_result2 <- data.frame(results$cell_labels, ID = data_result1$ID[which(data_result1$results.cell_labels %in% c("immune cells", "unknown"))])

# Identify T cells within the previously classified immune cells
# Filter the expression data for T cells or unknown labels
exp_data_expression_input3 <- data.frame(exp_data_expression_input2[,-32], label = data_result2$results.cell_labels)
exp_data_expression_input3 <- exp_data_expression_input3[which(exp_data_expression_input3$label %in% c("CD3+ T cells", "unknown")),]
colnames(exp_data_expression_input3)[1:31] <- colnames(exp_data_expression)

# Define expression signatures for T cell subtypes
list_expression_by_cell_type <- list(
  "CD4+ T cells" = "CD4 - T helper cells:Cyc_6_ch_3",
  "CD8+ T cells" = "CD8 - cytotoxic T cells:Cyc_3_ch_2"
)

# Execute SCINA for T cell subtyping
results <- SCINA(t(exp_data_expression_input2[,-32]), list_expression_by_cell_type, max_iter = 100, convergence_n = 10, 
                 convergence_rate = 0.999, sensitivity_cutoff = 0.9, rm_overlap = FALSE, allow_unknown = TRUE)

# Display the table of SCINA classification results for T cells
table(results$cell_labels)

# Prepare the results with IDs for merging
data_result3 <- data.frame(results$cell_labels, ID = data_result2$ID[which(data_result1$results.cell_labels %in% c("CD3+ T cells", "unknown"))])

# Merge initial and refined classification results
final_data_result <- merge(data_result1, data_result2, by = "ID", all.x = TRUE)
head(final_data_result)

# Update final labels, resolving "unknown" labels and preferring refined labels over initial broad classifications
final_label <- ifelse(final_data_result$results.cell_labels.x %in% c("immune cells", "unknown"), final_data_result$results.cell_labels.y, final_data_result$results.cell_labels.x)
final_label <- ifelse(final_label == "unknown", "Others", final_label)

# Compile the final outcome with predicted and true labels
data_outcome <- data.frame(predicted = final_label, true = data$ClusterName)

# Write the final data outcome to a CSV file
write.csv(data_outcome, "result_CRC_SCINA.csv")








#-----------------------------------------------------------------------------

###Human intestine datasets (Figure 2 bottom)

#-----------------------------------------------------------------------------

#### 1) Load the data
data <-  readr::read_csv("23_09_CODEX_HuBMAP_alldata_Dryad_merged.csv")
#### 2) Load signature matrix
Signature <- readr::read_csv("Signature_HI_TACIT_SCINA_Louvain.csv")


#### 3) Quality Control
#Rename cell type
data$Cell.Type=ifelse(data$Cell.Type%in%c("CD57+ Enterocyte","CD66+ Enterocyte",
                                          "MUC1+ Enterocyte"),"Enterocyte",data$Cell.Type)

dim(data)

##assign all NA as 0
Signature[is.na(Signature)==T]=0 
Signature=as.data.frame(Signature)
##validate marker in spatial
data_sub=data[which(data$unique_region=="B011_Sigmoid"),]
marker=""
ggplot(data, aes(x = X, y = Y)) +
  geom_point(aes_string(color = marker), size = 0.1)  +
  theme_classic(base_size = 15) +
  theme() +ggtitle(marker)+
  scale_color_gradient(low = "grey90", high = "darkblue", limits = c(0, 10), na.value = "darkblue")

#### 4) Run SCINA

# Define expression signatures

list_expression_by_cell_type=list("B"=colnames(Signature[1,which(Signature[1,]==1|Signature[1,]=="NA")]),
                                  "CD4+ T cell"=colnames(Signature[2,which(Signature[2,]==1|Signature[2,]=="NA")]),
                                  "CD7+ Immune"=colnames(Signature[5,which(Signature[5,]==1|Signature[5,]=="NA")]),
                                  "Cycling TA"= colnames(Signature[7,which(Signature[7,]==1|Signature[7,]=="NA")]),
                                  "DC"=colnames(Signature[8,which(Signature[8,]==1|Signature[8,]=="NA")]),
                                  "Endothelial"=colnames(Signature[9,which(Signature[9,]==1|Signature[9,]=="NA")]),
                                  "Enterocyte"="Cytokeratin",
                                  "Goblet"=colnames(Signature[11,which(Signature[11,]==1|Signature[11,]=="NA")]),
                                  "ICC"= colnames(Signature[12,which(Signature[12,]==1|Signature[12,]=="NA")]),
                                  "Lymphatic"="Podoplanin",
                                  "M1 Macrophage"= colnames(Signature[14,which(Signature[14,]==1|Signature[14,]=="NA")]),
                                  "Nerve"="Synapto",
                                  "Neuroendocrine"=colnames(Signature[18,which(Signature[18,]==1|Signature[18,]=="NA")]),
                                  "Neutrophil"= colnames(Signature[19,which(Signature[19,]==1|Signature[19,]=="NA")]),
                                  "NK"="CD57",
                                  "Paneth"=colnames(Signature[21,which(Signature[21,]==1|Signature[21,]=="NA")]),
                                  "Plasma"= colnames(Signature[22,which(Signature[22,]==1|Signature[22,]=="NA")]),
                                  "Smooth muscle"="aSMA",
                                  "Stroma"="Vimentin",
                                  "CD57+ Enterocyte"=colnames(Signature[3,which(Signature[3,]==1|Signature[3,]=="NA")]),
                                  "CD66+ Enterocyte"=colnames(Signature[4,which(Signature[4,]==1|Signature[4,]=="NA")]),
                                  "CD8+ T"=colnames(Signature[6,which(Signature[6,]==1|Signature[6,]=="NA")]),
                                  "M2 Macrophage"=colnames(Signature[15,which(Signature[15,]==1|Signature[15,]=="NA")]),
                                  "MUC1+ Enterocyte"=colnames(Signature[16,which(Signature[16,]==1|Signature[16,]=="NA")])
)

# Prepare expression data
data_expression=data[,colnames(Signature)[-c(1)]]
# Run SCINA
results = SCINA(t(data_expression), list_expression_by_cell_type, max_iter = 600, convergence_n = 30, 
                convergence_rate = 0.999, sensitivity_cutoff = 0.9, rm_overlap=FALSE, allow_unknown=TRUE)

# Compile the final outcome with predicted and true labels
data_result1=data.frame(results$cell_labels,ID=data$V1)
data_outcome=data.frame(predicted=data_result1$results.cell_labels,true=data$Cell.Type)
#Output
write.csv(data_outcome,"~/Documents/Compare_TACIT/PCF/SCINA/result_HI.csv")





