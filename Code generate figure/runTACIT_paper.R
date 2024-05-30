#-----------------------------------------------------------------------------
###Goal: Obtain TACIT: Threshold-based Assignment for Cell Type Identifying in Spatial Omics
#-----------------------------------------------------------------------------
###Cell type annotation Step: 
#### 1) Load the data
#### 2) Load signature matrix
#### 3) Quality Control
#### 4) Run TACIT
#-----------------------------------------------------------------------------
library(ggplot2)
library(Seurat)
library(segmented)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------

###Colon datasets (Figure 2 top)

#-----------------------------------------------------------------------------


#### 1) Load the data
data <-  readr::read_csv("CRC_DataSet.csv")
#### 2) Load signature matrix
Signature <- readr::read_csv("Signature_CRC_TACIT_SCINA_Louvain.csv")


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

##assign all NA as 0
Signature[is.na(Signature)==T]=0 
Signature=as.data.frame(Signature)
##validate marker in spatial
marker="`CD31 - vasculature:Cyc_19_ch_3`"
ggplot(data, aes(x = X, y = Y)) +
  geom_point(aes_string(color = marker), size = 0.1)  +
  theme_classic(base_size = 15) +
  theme() +ggtitle(marker)+
  scale_color_gradient(low = "grey90", high = "darkblue", limits = c(0, 100), na.value = "darkblue")

#### 4) Run TACIT
CRC_TACIT=TACIT_update(data_anb = (data[,colnames(Signature)[-c(1)]]),r = 200,p=20,Signature = Signature)

#Get annotation
data_predicted=CRC_TACIT[3]
#Get annotation with cell ID
data_predicted=data.frame(cell_ID=data$CellID,TACIT=data_predicted$mem)
#Validate TACIT
confusionMatrix(data=factor(data_predicted$TACIT,levels = c(unique(data$ClusterName))),
                reference = factor(data$ClusterName,levels = c(unique(data$ClusterName))))

write.csv(data_predicted,"result_TACIT.csv")



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

#### 4) Run TACIT
CRC_TACIT=TACIT_update(data_anb = (data[,colnames(Signature)[-c(1)]]),r = 20,p=20,Signature = Signature)

#Get annotation
data_predicted=CRC_TACIT[3]
#Get annotation with cell ID
data_predicted=data.frame(TACIT=data_predicted$mem)
#Validate TACIT
confusionMatrix(data=factor(data_predicted$TACIT,levels = c(unique(data$Cell.Type))),
                reference = factor(data$Cell.Type,levels = c(unique(data$Cell.Type))))


write.csv(data_predicted,"result_HI.csv")



#-----------------------------------------------------------------------------

###MERFISH (Figure Supplement MERFISH)

#-----------------------------------------------------------------------------

#### 1) Load the data
data <-  readr::read_csv("Moffitt_and_Bambah-Mukku_et_al_merfish_all_cells.csv")
#### 2) Load signature matrix
Signature <- readr::read_csv("Signature_MERFISH_TACIT_Louvain.csv")


#### 3) Quality Control
#Remove ambiguous
data=data[which(data$Cell_class!="Ambiguous"),]

#Rename cell type
data$Cell_class=ifelse(data$Cell_class%in%c("Endothelial 1","Endothelial 2","Endothelial 3"),"Endothelial",data$Cell_class)
data$Cell_class=ifelse(data$Cell_class%in%c("OD Immature 1","OD Immature 2"),"OD Immature",data$Cell_class)
data$Cell_class=ifelse(data$Cell_class%in%c("OD Mature 1","OD Mature 2","OD Mature 3","OD Mature 4"),"OD Mature",data$Cell_class)

##assign all NA as 0
Signature[is.na(Signature)==T]=0 
Signature=as.data.frame(Signature)
##validate marker in spatial
marker=""
ggplot(data, aes(x = X, y = Y)) +
  geom_point(aes_string(color = marker), size = 0.1)  +
  theme_classic(base_size = 15) +
  theme() +ggtitle(marker)+
  scale_color_gradient(low = "grey90", high = "darkblue", limits = c(0, 10), na.value = "darkblue")


#### 4) Run TACIT
CRC_TACIT=TACIT_update(data_anb = (data[,colnames(Signature)[-c(1)]]),r = 20,p=50,Signature = Signature)

#Get annotation
data_predicted=CRC_TACIT[3]
#Get annotation with cell ID
data_predicted=data.frame(TACIT=data_predicted$mem)
#Validate TACIT
confusionMatrix(data=factor(data_predicted$TACIT,levels = c(unique(data$Cell_class))),
                reference = factor(data$Cell_class,levels = c(unique(data$Cell_class))))




write.csv(data_predicted,"result_MERFISH.csv")




#-----------------------------------------------------------------------------

###Xenium (Figure 3)

#-----------------------------------------------------------------------------

#### 1) Load the data
data <-  readr::read_csv("data_expression_Xenium.csv")
#### 2) Load signature matrix
Signature <- readr::read_csv("Signature_Xenium_TACIT_Louvain.csv")


#### 3) Quality Control
##assign all NA as 0
Signature[is.na(Signature)==T]=0 
Signature=as.data.frame(Signature)
##validate marker in spatial
marker=""
ggplot(data, aes(x = X, y = Y)) +
  geom_point(aes_string(color = marker), size = 0.1)  +
  theme_classic(base_size = 15) +
  theme() +ggtitle(marker)+
  scale_color_gradient(low = "grey90", high = "darkblue", limits = c(0, 10), na.value = "darkblue")


#### 4) Run TACIT
CRC_TACIT=TACIT_update(data_anb = (data[,colnames(Signature)[-c(1)]]),r = 20,p=50,Signature = Signature)

#Get annotation
data_predicted=CRC_TACIT[3]
#Get annotation with cell ID
data_predicted=data.frame(TACIT=data_predicted$mem)


write.csv(data_predicted,"result_Xenium.csv")








#-----------------------------------------------------------------------------

###PCF (Figure 4 & Figure 5)

#-----------------------------------------------------------------------------

#### 1) Load the data
data <-  readr::read_csv("region_1_PCF_TACIT.csv")
#### 2) Load signature matrix
Signature <- readr::read_csv("Signature_PCF_Figure_4.csv")


#### 3) Quality Control
##assign all NA as 0
Signature[is.na(Signature)==T]=0 
Signature=as.data.frame(Signature)
##validate marker in spatial
marker=""
ggplot(data, aes(x = X, y = Y)) +
  geom_point(aes_string(color = marker), size = 0.1)  +
  theme_classic(base_size = 15) +
  theme() +ggtitle(marker)+
  scale_color_gradient(low = "grey90", high = "darkblue", limits = c(0, 10), na.value = "darkblue")


#### 4) Run TACIT
CRC_TACIT=TACIT_update(data_anb = (data[,colnames(Signature)[-c(1)]]),r = 20,p=15,Signature = Signature)

#Get annotation
data_predicted=CRC_TACIT[3]
#Get annotation with cell ID
data_predicted=data.frame(TACIT=data_predicted$mem)


write.csv(data_predicted,"result_PCF_region.csv")

##Note
#Repeated this for all 6 region






#-----------------------------------------------------------------------------

###Xenium (Figure 4 & Figure 5)

#-----------------------------------------------------------------------------

#### 1) Load the data
data <-  readr::read_csv("region_1_Xenium_TACIT.csv")
#### 2) Load signature matrix
Signature <- readr::read_csv("Signature_Xenium_TACIT_Louvain.csv")


#### 3) Quality Control
##assign all NA as 0
Signature[is.na(Signature)==T]=0 
Signature=as.data.frame(Signature)
##validate marker in spatial
marker=""
ggplot(data, aes(x = X, y = Y)) +
  geom_point(aes_string(color = marker), size = 0.1)  +
  theme_classic(base_size = 15) +
  theme() +ggtitle(marker)+
  scale_color_gradient(low = "grey90", high = "darkblue", limits = c(0, 10), na.value = "darkblue")


#### 4) Run TACIT
CRC_TACIT=TACIT_update(data_anb = (data[,colnames(Signature)[-c(1)]]),r = 20,p=50,Signature = Signature)

#Get annotation
data_predicted=CRC_TACIT[3]
#Get annotation with cell ID
data_predicted=data.frame(TACIT=data_predicted$mem)


write.csv(data_predicted,"result_Xenium_region.csv")

##Note
#Repeated this for all 6 region













#-----------------------------------------------------------------------------

###PCF (Figure 6)

#-----------------------------------------------------------------------------

#### 1) Load the data
data <-  readr::read_csv("PCF_Xenium_region_001.csv")
#### 2) Load signature matrix
Signature <- readr::read_csv("Signature_PCF_figure6.csv")


#### 3) Quality Control
##assign all NA as 0
Signature[is.na(Signature)==T]=0 
Signature=as.data.frame(Signature)
##validate marker in spatial
marker=""
ggplot(data, aes(x = X, y = Y)) +
  geom_point(aes_string(color = marker), size = 0.1)  +
  theme_classic(base_size = 15) +
  theme() +ggtitle(marker)+
  scale_color_gradient(low = "grey90", high = "darkblue", limits = c(0, 10), na.value = "darkblue")


#### 4) Run TACIT
CRC_TACIT=TACIT_update(data_anb = (data[,colnames(Signature)[-c(1)]]),r = 20,p=10,Signature = Signature)

#Get annotation
data_predicted=CRC_TACIT[3]
#Get annotation with cell ID
data_predicted=data.frame(TACIT=data_predicted$mem)


write.csv(data_predicted,"result_PCF_Xenium_region_001_PCF.csv")

##Note
#Repeated this for all 6 region






#-----------------------------------------------------------------------------

###Xenium (Figure 6)

#-----------------------------------------------------------------------------

#### 1) Load the data
data <-  readr::read_csv("PCF_Xenium_region_001.csv")
#### 2) Load signature matrix
Signature <- readr::read_csv("Signature_Xenium_figure6.csv")


#### 3) Quality Control
##assign all NA as 0
Signature[is.na(Signature)==T]=0 
Signature=as.data.frame(Signature)
##validate marker in spatial
marker=""
ggplot(data, aes(x = X, y = Y)) +
  geom_point(aes_string(color = marker), size = 0.1)  +
  theme_classic(base_size = 15) +
  theme() +ggtitle(marker)+
  scale_color_gradient(low = "grey90", high = "darkblue", limits = c(0, 10), na.value = "darkblue")


#### 4) Run TACIT
CRC_TACIT=TACIT_update(data_anb = (data[,colnames(Signature)[-c(1)]]),r = 20,p=50,Signature = Signature)

#Get annotation
data_predicted=CRC_TACIT[3]
#Get annotation with cell ID
data_predicted=data.frame(TACIT=data_predicted$mem)


write.csv(data_predicted,"result_PCF_Xenium_region_001_Xenium.csv")

##Note
#Repeated this for all 6 region








#-----------------------------------------------------------------------------

###Xenium & PCF (Figure 6)

#-----------------------------------------------------------------------------

#### 1) Load the data
data <-  readr::read_csv("PCF_Xenium_region_001.csv")
#### 2) Load signature matrix
Signature <- readr::read_csv("Signature_PCF_Xenium_Figure6.csv")


#### 3) Quality Control
##assign all NA as 0
Signature[is.na(Signature)==T]=0 
Signature=as.data.frame(Signature)
##validate marker in spatial
marker=""
ggplot(data, aes(x = X, y = Y)) +
  geom_point(aes_string(color = marker), size = 0.1)  +
  theme_classic(base_size = 15) +
  theme() +ggtitle(marker)+
  scale_color_gradient(low = "grey90", high = "darkblue", limits = c(0, 10), na.value = "darkblue")


#### 4) Run TACIT
CRC_TACIT=TACIT_update(data_anb = (data[,colnames(Signature)[-c(1)]]),r = 20,p=50,Signature = Signature)

#Get annotation
data_predicted=CRC_TACIT[3]
#Get annotation with cell ID
data_predicted=data.frame(TACIT=data_predicted$mem)


write.csv(data_predicted,"result_PCF_Xenium_region_001_Xenium.csv")

##Note
#Repeated this for all 6 region











