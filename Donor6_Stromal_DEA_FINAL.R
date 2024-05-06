library(dplyr)
library(Seurat)
library(patchwork)
library(plotly)
library(data.table)
library(magrittr)
library(scProportionTest)
library(scCustomize)
library(qs)
library(ggplot2)

setwd("C:\\Users\\BASanto\\Desktop\\GUDMAP_Manuscripts\\scRNAseq_v_snRNAseq_MS\\Making_Final_UMAPs\\Data")

########
# Main #
########

# Load the main dataset
data_path <- paste(getwd(),"Donor_6_MS_Stromal_Compartment_FINAL.rds",sep="\\")
dataset <- "Donor_6_MS_Stromal_Compartment_FINAL.rds"
data <- readRDS(data_path)

# Set the active ident to be the cell type and reorder for figure generation
Idents(data) <- data@meta.data$cell_type
data$cell_type <- factor(data$cell_type,levels=c("Schwann","VSMC","SMC","Lymphatic","Arterial","Venous","EC","FB5","FB4","FB3","FB2","FB1"))

# Create ident for single nucleus and single cell subpopulations
Idents(data) <- data@meta.data$orig.ident
samples <- c("Run_Count_Donor-6_Dome","Run_Count_Donor-6_Dome-Nuc","Run_Count_Donor-6_Neck",
             "Run_Count_Donor-6_Neck-Nuc","Run_Count_Donor-6_UO","Run_Count_Donor-6_UO-Nuc",
             "Run_Count_Donor-6_UVJ","Run_Count_Donor-6_UVJ-Nuc")
new_labels <- c("SingleCell","SingleNuc","SingleCell","SingleNuc","SingleCell","SingleNuc",
                "SingleCell","SingleNuc")
for (i in 1:length(unique(Idents(data)))) {
  Idents(object=data,cells=WhichCells(object=data,ident=samples[i])) <- new_labels[i]
}
umap <- DimPlot(data,reduction="umap",label=FALSE,label.size=6)
umap
data@meta.data["sc_sn"] <- Idents(data)

# Plot sc and sn UMAPs
sc_select <- WhichCells(data,idents=c("SingleCell"))
sn_select <- WhichCells(data,idents=c("SingleNuc"))

sc <- DimPlot(data,label=FALSE,group.by="sc_sn",cells.highlight=list(sc_select),cols.highlight=c("#00BFC4"),cols="grey") + NoLegend() + NoAxes() + ggtitle(NULL)
sn <- DimPlot(data,label=FALSE,group.by="sc_sn",cells.highlight=list(sn_select),cols.highlight=c("#F87667"),cols="grey") + NoLegend() + NoAxes() + ggtitle(NULL)

#ggsave("donor6_stromal_umap_SC.jpg",plot=sc,width=20,height=20,units="cm",dpi=600)
#ggsave("donor6_stromal_umap_SN.jpg",plot=sn,width=20,height=20,units="cm",dpi=600)

# DEA
out_dir = ("C:\\Users\\BASanto\\Desktop\\GUDMAP_Manuscripts\\scRNAseq_v_snRNAseq_MS\\Making_Final_UMAPs\\Stromal_DEA")
if (file.exists(out_dir)) {} else {dir.create(out_dir)}

Idents(data) <- data@meta.data$sc_sn
cell_list <- unique(data@meta.data$cell_type)
min_pct <- 0.25
min_fc <- 1.5
max_q <- 0.05

has_enough <- vector()
sn_deg_list <- vector()
sc_deg_list <- vector()

cell_list <- list("FB1","FB2","FB3","FB4","FB5","EC","Arterial","Lymphatic","SMC","VSMC")

for (cell in cell_list) {
  
  # Run DEA for this cell type for singlenuc data
  sn_subset <- subset(data,idents="SingleNuc")
  Idents(sn_subset) <- sn_subset@meta.data$cell_type
  sn_subset@active.assay <- "RNA"
  sn_subset <- NormalizeData(sn_subset)
  sn_dea <- FindMarkers(sn_subset,ident.1=cell,min.pct=min_pct,logfc.threshold=min_fc)
  sn_dea['pct_diff']<- sn_dea["pct.1"]-sn_dea["pct.2"]
  sn_dea <- sn_dea[order(sn_dea$pct_diff,decreasing=TRUE),]
  sn_dea <- sn_dea[sn_dea$p_val_adj<max_q,]
  out_file <- paste(out_dir,"\\",cell,"_sn_degs.csv",sep="")
  #write.csv(sn_dea,out_file)
  
  # Run DEA for this cell type for singlecell data
  sc_subset <- subset(data,idents="SingleCell")
  Idents(sc_subset) <- sc_subset@meta.data$cell_type
  sc_subset@active.assay <- "RNA"
  sc_subset <- NormalizeData(sc_subset)
  sc_dea <- FindMarkers(sc_subset,ident.1=cell,min.pct=min_pct,logfc.threshold=min_fc)
  sc_dea['pct_diff']<- sc_dea["pct.1"]-sc_dea["pct.2"]
  sc_dea <- sc_dea[order(sc_dea$pct_diff,decreasing=TRUE),]
  sc_dea <- sc_dea[sc_dea$p_val_adj<max_q,]
  out_file <- paste(out_dir,"\\",cell,"_sc_degs.csv",sep="")
  #write.csv(sc_dea,out_file)
  
  if ((length(sn_dea$pct_diff)>=10)&(length(sc_dea$pct_diff)>=10)) {
    
    sn_ten <- sn_dea[1:10,]
    sc_ten <- sc_dea[1:10,]
    
    has_enough[length(has_enough)+1] <- cell
    sn_deg_list[(length(sn_deg_list)+1):(length(sn_deg_list)+10)] <- row.names(sn_ten)
    sc_deg_list[(length(sc_deg_list)+1):(length(sc_deg_list)+10)] <- row.names(sc_ten)
    
  }
}

Idents(data) <- data@meta.data$cell_type
data_subset <- subset(data,idents=c("FB1","FB2","FB3","FB4","FB5","EC","Arterial","SMC","VSMC","Lymphatic"))

data_subset@active.assay <- "RNA"
data_subset <- NormalizeData(data_subset)
all.genes <- rownames(data_subset)
data_subset <- ScaleData(data_subset,features=all.genes)

data_subset$cell_type <- factor(data_subset$cell_type,levels=c("FB1","FB2","FB3","FB4","FB5","EC","Arterial","Lymphatic","SMC","VSMC"))

cmap <- c("#CC3399","#FF00CC","#FF3399","#FF99CC","#660066","#FF9966","#FF3333","#FFCC00","#990066","#CC0033")

Idents(data_subset) <- data_subset@meta.data$sc_sn
data_subset_sn <- subset(data_subset,idents=c("SingleNuc"))
data_subset_sn@active.assay <- "RNA"
data_subset_sc <- subset(data_subset,idents=c("SingleCell"))
data_subset_sc@active.assay <- "RNA"

Idents(data_subset_sn) <- data_subset_sn@meta.data$cell_type
data_subset_sn <- NormalizeData(data_subset_sn)
all.genes <- rownames(data_subset_sn)
data_subset_sn <- ScaleData(data_subset_sn,features=all.genes)

Idents(data_subset_sc) <- data_subset_sc@meta.data$cell_type
data_subset_sc <- NormalizeData(data_subset_sc)
all.genes <- rownames(data_subset_sc)
data_subset_sc <- ScaleData(data_subset_sc,features=all.genes)

mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)

sn_master_map <- DoHeatmap(data_subset_sn,features=c(sn_deg_list),group.colors=cmap,label=FALSE) +  scale_fill_gradientn(colours = rev(mapal)) + NoLegend()
sn_master_map
#ggsave("donor6_stromal_DEA_SN_15Log2FC_PctDiff.jpg",plot=sn_master_map,width=25,height=25,units="cm",dpi=600)

sc_master_map <- DoHeatmap(data_subset_sc,features=c(sc_deg_list),group.colors=cmap,label=FALSE) +  scale_fill_gradientn(colours = rev(mapal)) + NoLegend()
sc_master_map
#ggsave("donor6_stromal_DEA_SC_15Log2FC_PctDiff.jpg",plot=sc_master_map,width=25,height=25,units="cm",dpi=600)

## END ##

