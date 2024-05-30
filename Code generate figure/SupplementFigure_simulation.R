library(ggplot2)
library(tidyr)
library(dplyr)


#a

data=read_csv("~/Documents/TACIT figure/dataset/Figure 2/CRC_DataSet.csv")
data_predicted_return=read_csv("~/Documents/TACIT figure/dataset/Supplement_simulation/data_predicted_return_boostrap.csv")
data_predicted_return=data_predicted_return[,-1]
data_ID=read_csv("~/Documents/TACIT figure/dataset/Supplement_simulation/data_ID_boostrap.csv")
data_ID=data_ID[,-1]
threshold_norm=read_csv("~/Documents/TACIT figure/dataset/Supplement_simulation/threshold_boostrap.csv")
threshold_norm=threshold_norm[,-1]



# Calculate the standard deviation for each cell type
sd_values <- apply(threshold_norm[,-11], 1, sd)

# Combine the results into a new data frame
result <- data.frame(CellType = threshold_norm[,11], StandardDeviation = sd_values)

# 

# Assume 'threshold_norm' is loaded as a matrix
# Convert matrix to dataframe and tidy it
threshold_norm <- data.frame(threshold_norm)  # If it's not already a dataframe
threshold_norm$CellType <- rownames(threshold_norm)  # Ensure rownames are cell types


# Assuming 'threshold_norm' is your matrix and it has row names indicating cell types
threshold_norm <- data.frame(threshold_norm, row.names = NULL)  # Convert to dataframe without losing row names
threshold_norm$CellType <- rownames(threshold)  # Add cell types as a new column

threshold_long <- pivot_longer(threshold_norm,
                               cols = -CellType, 
                               names_to = "Measurement", 
                               values_to = "Value")

# Create box plots
ggplot(threshold_long, aes(x = CellType, y = Value, fill = CellType)) +
  geom_boxplot()  +  # Adjust text angle for better readability
  labs(title = "",
       x = "",
       y = "")+theme_classic(base_size = 20)+
  theme(axis.text.x = element_text(angle = 90),legend.position = "none") +
  scale_fill_manual(values = c("B cells"="tan","CD11c+ DCs"="darkred","CD3+ T cells"="pink",
                                "CD4+ T cells"="violet","CD4+ T cells CD45RO+"="darkviolet",
                                "macrophages"="green","CD8+ T cells"="darkgreen",
                                "granulocytes"="blue","lymphatics"="skyblue","NK cells"="yellow",
                                "Others"="grey90","immune cells"="darkblue","plasma cells"="orange","stroma"="gold",
                                "Tregs"="black","tumor cells"="red","vasculature"="darkcyan","smooth muscle"="tan4"))






#b


F1_value=NULL
Recall_value=NULL
Precision_value=NULL
mat_F1=matrix(NA,ncol = 10,nrow = 17)
mat_Recall=matrix(NA,ncol = 10,nrow = 17)
mat_Precision=matrix(NA,ncol = 10,nrow = 17)


for (i in 1:10) {
  pred_val=as.character(t(data_predicted_return[,i]))
  
  data_sub=data[which(data$CellID%in%as.numeric(t(data_ID[,i]))),]
  
  ref_val=data_sub$ClusterName

  val=calculate_metrics(predictions = pred_val[which(pred_val!="Others")],true_labels = ref_val[which(pred_val!="Others")])
  F1_value[i]=val$WeightedF1
  Precision_value[i]=val$WeightedPrecision
  Recall_value[i]=val$WeightedRecall
  mat_F1[,i]=val$F1Scores
  mat_Recall[,i]=val$Sensitivity
  mat_Precision[,i]=val$Precision
  
}
rownames(mat_F1)=names(val$Sensitivity)
rownames(mat_Recall)=names(val$Sensitivity)
rownames(mat_Precision)=names(val$Sensitivity)
rownames(mat_F1)[8]="Other immune"
rownames(mat_Recall)[8]="Other immune"
rownames(mat_Precision)[8]="Other immune"


# Function to calculate mean, standard deviation, and 95% confidence interval
calculate_stats <- function(values) {
  mean_value <- mean(values)
  sd_value <- sd(values)
  n <- length(values)
  error <- qt(0.975, df=n-1) * sd_value / sqrt(n)
  ci_lower <- mean_value - error
  ci_upper <- mean_value + error
  list(mean = mean_value, sd = sd_value, ci_lower = ci_lower, ci_upper = ci_upper)
}

# Calculate stats for F1_value, Recall_value, and Precision_value
F1_stats <- calculate_stats(F1_value)
Recall_stats <- calculate_stats(Recall_value)
Precision_stats <- calculate_stats(Precision_value)

# Print the results
print(F1_stats)
print(Recall_stats)
print(Precision_stats)



# Convert to data frame
mat_Recall <- as.data.frame(t(mat_Recall))

# Reshape data for ggplot2
data_long <- melt(mat_Recall, variable.name = "Cell_Type", value.name = "F1_Score")

# Plot using ggplot2
ggplot(data_long, aes(x = Cell_Type, y = F1_Score, fill = Cell_Type)) +
  geom_boxplot() +
  theme_classic(base_size = 20) +
  labs(title = "", x = "", y = "") +
  theme(axis.text.x = element_text(angle = 90), legend.position = "none")+ylim(0,1)




# Convert to data frame
mat_F1 <- as.data.frame(t(mat_F1))

# Reshape data for ggplot2
data_long <- melt(mat_F1, variable.name = "Cell_Type", value.name = "F1_Score")

# Plot using ggplot2
ggplot(data_long, aes(x = Cell_Type, y = F1_Score, fill = Cell_Type)) +
  geom_boxplot() +
  theme_classic(base_size = 20) +
  labs(title = "", x = "", y = "") +
  theme(axis.text.x = element_text(angle = 90), legend.position = "none")+ylim(0,1)


# Convert to data frame
mat_Precision <- as.data.frame(t(mat_Precision))

# Reshape data for ggplot2
data_long <- melt(mat_Precision, variable.name = "Cell_Type", value.name = "F1_Score")

# Plot using ggplot2
ggplot(data_long, aes(x = Cell_Type, y = F1_Score, fill = Cell_Type)) +
  geom_boxplot() +
  theme_classic(base_size = 20) +
  labs(title = "", x = "", y = "") +
  theme(axis.text.x = element_text(angle = 90), legend.position = "none")+ylim(0,1)



#b

# Combine data into a data frame
data_val <- data.frame(
  Index = 1:10,
  F1 = F1_value,
  Precision = Precision_value,
  Recall = Recall_value
)

# Reshape data for ggplot2
data_long <- reshape2::melt(data_val, id.vars = "Index", variable.name = "Metric", value.name = "Value")

# Plot using ggplot2
ggplot(data_long, aes(x = Index, y = Value, color = Metric)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(title = "F1, Precision, and Recall Values", x = "Index", y = "Value") +
  theme(legend.title = element_blank())



# Combine data into a data frame
data_val <- data.frame(
  F1 = F1_value,
  Precision = Precision_value,
  Recall = Recall_value
)

# Reshape data for ggplot2
data_long <- melt(data_val, variable.name = "Metric", value.name = "Value")

data_long$Metric=factor(data_long$Metric,levels = c("Recall","Precision","F1"))

# Plot using ggplot2
ggplot(data_long, aes(x = Metric, y = Value, fill = Metric)) +
  geom_boxplot() +
  theme_classic(base_size = 20) +
  labs(title = "", x = "", y = "") +
  theme(legend.position = "none")






data_resolution=read_csv("~/Documents/TACIT figure/dataset/Supplement_simulation/diff_res.csv")
data_resolution=data_resolution[,-1]
F1_value=NULL
Recall_value=NULL
Precision_value=NULL
for (i in 1:12) {
  pred_val=as.character(t(data_resolution[,i]))
  
  ref_val=data$ClusterName
  val=calculate_metrics(predictions = pred_val[which(pred_val!="Others")],true_labels = ref_val[which(pred_val!="Others")])
  F1_value[i]=val$WeightedF1
  Precision_value[i]=val$WeightedPrecision
  Recall_value[i]=val$WeightedRecall
  
}
colnames(data_resolution)
F1_value
Precision_value
Recall_value


# Given data
r_values <- colnames(data_resolution)
F1_value <- F1_value
Precision_value <- Precision_value
Recall_value <- Recall_value

# Convert r_values to numeric for ordering
r_numeric <- as.numeric(sub("r = ", "", r_values))

# Create a data frame
data_plot <- data.frame(
  r_values = factor(r_values, levels = r_values[order(r_numeric)]),
  F1 = F1_value,
  Precision = Precision_value,
  Recall = Recall_value
)

# Melt the data frame for ggplot
data_melt <- tidyr::pivot_longer(data_plot, cols = c("F1", "Precision", "Recall"), names_to = "Metric", values_to = "Value")

# Plot the data
ggplot(data_melt, aes(x = r_values, y = Value, color = Metric, group = Metric)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  labs(
    title = "",
    x = "",
    y = ""
  ) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values = c("F1" = "blue", "Precision" = "green", "Recall" = "red"))
















