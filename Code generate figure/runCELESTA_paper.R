#-----------------------------------------------------------------------------
###Goal: Obtain CELESTA: CELl typE identification with SpaTiAl information
#-----------------------------------------------------------------------------
###Cell type annotation Step: 
#### 1) Load the data
#### 2) Load signature matrix
#### 3) Quality Control
#### 4) Run CELESTA
#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------
library(CELESTA)
library(Rmixmod)
library(spdep)
library(ggplot2)
library(reshape2)
library(zeallot)
library(caret)
#-----------------------------------------------------------------------------

###Colon datasets (Figure 2 top)

#-----------------------------------------------------------------------------



#### 1) Load the data
data <-  readr::read_csv("CRC_DataSet.csv")
#### 2) Load signature matrix
prior_marker_info=read.csv("Signature_CRC_Celesta.csv")


#### 3) Quality Control
#Rename cell type
data$ClusterName=ifelse(data$ClusterName%in%c("CD11b+ CD68+ macrophages","CD163+ macrophages",
                                              "CD68+ macrophages","CD68+ macrophages GzmB+",
                                              "CD68+ CD163+ macrophages"),"macrophages",data$ClusterName)

data$ClusterName=ifelse(data$ClusterName%in%c("CD4+ T cells GATA3+"),"CD4+ T cells",data$ClusterName)


#Cells identified as immune/vasculature (n=2,153) and immune/tumor (n=1,797), 
#along with cells lacking a marker signature—including adipocytes (n=1,811), 
#nerves (n=659), undefined (n=6,524), monocytes (n=815), and cells categorized 
#as dirt (n=7,357)—were excluded from the analysis

data=data[which(data$ClusterName %in% Signature$cell_type),]
dim(data)

# Ensure the correct column names from prior_marker_info are used
marker_columns <- colnames(prior_marker_info)[-c(1, 2)]
data_sub <- data.frame(
  ClusterName = data$ClusterName,
  FileName = data$File.Name,
  CellID = data$CellID,
  X = data$X,
  Y = data$Y,
  TMA = data$TMA_AB,
  data[, marker_columns]
)

#### 4) run CELESTA

colnames(prior_marker_info)[1]=""
for (i in 1:length(unique(data$File.Name))) {
  imaging_data_sub_region=data_sub[which(data$File.Name==unique(data$File.Name)[i]),]
  
  CelestaObj <- CreateCelestaObject(project_title = "CRC",prior_marker_info,imaging_data_sub_region)
  
  CelestaObj <- FilterCells(CelestaObj,high_marker_threshold=0.9, low_marker_threshold=0.4)
  
  CelestaObj <-AssignCells(CelestaObj,max_iteration=10,cell_change_threshold=0.01,
                           high_expression_threshold_anchor=c(0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9),
                           low_expression_threshold_anchor=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
                           high_expression_threshold_index=c(0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3),
                           low_expression_threshold_index=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1))
  if(i==1){
    pred=CelestaObj@final_cell_type_assignment[,"Final cell type"]
    data_outcome_final=data.frame(pred,X=imaging_data_sub_region$X,Y=imaging_data_sub_region$Y,true=imaging_data_sub_region$ClusterName,ID=imaging_data_sub_region$CellID)
    
  }else{
    pred=CelestaObj@final_cell_type_assignment[,"Final cell type"]
    data_outcome=data.frame(pred,X=imaging_data_sub_region$X,Y=imaging_data_sub_region$Y,true=imaging_data_sub_region$ClusterName,ID=imaging_data_sub_region$CellID)
    data_outcome_final=rbind(data_outcome_final,data_outcome)
  }
  
  print(i)
}

# Calculating confusion matrix
confusionMatrix(data=factor(data_outcome_final$pred,levels = unique(data_outcome_final$true)),reference = factor(data_outcome_final$true,levels = unique(data_outcome_final$true)))

# Sort by ID after completing all operations
data_outcome_final=data_outcome_final[order(data_outcome_final$ID),]

write.csv(data_outcome_final,"result_CRC_CELESTA.csv")


#-----------------------------------------------------------------------------

###HI datasets (Figure 2 bottom)

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

#### 4) Run CELESTA

for (i in 1:length(unique(data$unique_region))) {
  imaging_data_sub_region=data[which(data$unique_region==unique(data$unique_region)[i]),]
  
  CelestaObj <- CreateCelestaObject(project_title = "HI",prior_marker_info,imaging_data_sub_region)
  
  CelestaObj <- FilterCells(CelestaObj,high_marker_threshold=0.9, low_marker_threshold=0.4)
  
  CelestaObj <-AssignCells(CelestaObj,max_iteration=10,cell_change_threshold=0.01,
                           high_expression_threshold_anchor=c(rep(0.7,22)),
                           low_expression_threshold_anchor=c(rep(1,22)),
                           high_expression_threshold_index=c(rep(0.3,22)),
                           low_expression_threshold_index=c(rep(1,22)))
  if(i==1){
    pred=CelestaObj@final_cell_type_assignment[,"Final cell type"]
    data_outcome_final=data.frame(pred,X=imaging_data_sub_region$X,Y=imaging_data_sub_region$Y,true=imaging_data_sub_region$ClusterName,ID=imaging_data_sub_region$CellID)
    
  }else{
    pred=CelestaObj@final_cell_type_assignment[,"Final cell type"]
    data_outcome=data.frame(pred,X=imaging_data_sub_region$X,Y=imaging_data_sub_region$Y,true=imaging_data_sub_region$ClusterName,ID=imaging_data_sub_region$CellID)
    data_outcome_final=rbind(data_outcome_final,data_outcome)
  }
  
  print(i)
}


data_outcome_final$pred[is.na(data_outcome_final$pred)==T]="Unknown"

write.csv(data_outcome_final,"result_HI_CELESTA.csv")










