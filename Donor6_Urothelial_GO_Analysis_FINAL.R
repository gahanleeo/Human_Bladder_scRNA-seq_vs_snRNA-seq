library(dplyr)
library(Seurat)
library(patchwork)
library(plotly)
library(data.table)
library(magrittr)
library(stringr)
library(gprofiler2)
library(htmlwidgets)

setwd("C:\\Users\\BASanto\\Desktop\\GUDMAP_Manuscripts\\scRNAseq_v_snRNAseq_MS\\Making_Final_UMAPs")

########
# Data #
########

file_dir <- paste(getwd(),"\\Urothelial_DEA\\15_Log2FC",sep="")
file_dir <- dir(file_dir)
num_files <- length(file_dir)

output_dir <- paste(getwd(),"\\Urothelial_DEA_GSEA",sep="")
if (!dir.exists(output_dir)) {dir.create(output_dir)}

select <- c("GO:MF","GP:CC","GO:BP","WB")

data_labels <- vector()
for (i in 3:num_files) {
  
  # Load the DEGs
  this_file <- file_dir[i]
  this_subset <- str_split(this_file,".csv")
  this_subset <- this_subset[[1]][1]
  data_labels[length(data_labels)+1] <- this_subset
  
  this_type_degs <- paste(getwd(),"\\Urothelial_DEA\\15_Log2FC\\",this_file,sep="")
  this_type_degs <- read.csv(this_type_degs,header=TRUE)
  
  # Subset the DEGs
  degs <- this_type_degs$X
  
  # Run the GSEA on the DEGs
  gostres <- gost(query = c(degs), 
                  organism = "hsapiens", ordered_query = FALSE, 
                  multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                  measure_underrepresentation = FALSE, evcodes = FALSE, 
                  user_threshold = 0.05, correction_method = "g_SCS", 
                  domain_scope = "annotated", custom_bg = NULL, 
                  numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = TRUE)
  
  go_plot <- gostplot(gostres, capped = TRUE, interactive = TRUE)
  go_plot
  
  plot_path <- paste(output_dir,"\\",this_subset,".html",sep="")
  saveWidget(go_plot,file=plot_path)
  
  table_path <- paste(output_dir,"\\",this_subset,".pdf",sep="")
  res_df <- gostres$result
  res_df <- filter(res_df, source %in% select)
  res_df <- res_df[order(res_df$p_value,decreasing=FALSE),]
  gostres$result <- res_df
  go_table <- publish_gosttable(gostres, use_colors = TRUE, highlight_terms = gostres$result[c(1:100),],
                                show_columns = c("source", "term_name","intersection_size"),
                                filename = table_path, ggplot= TRUE)
}



