# This is fork from paper about bladder tissue analysis, a good reference!

### Codes for reproducing analysis in *Exploring the Utility of snRNA-seq in Profiling Human Bladder Tissue: A Comprehensive Comparison with scRNA-seq*
For scRNA-seq and snRNA-seq quality control filtering, data integration and batch correction, and clustering optimization, see the pipeline by Sona et al., forked under basanto/ scRNA_Analysis_Developing_for_Human_Bladder_Analysis. Raw and processed sequence and metadata are available on the Gene Expresson Omnibus (GEO) Database under accession GSE267964.

Paper link: 
https://www.sciencedirect.com/science/article/pii/S2589004224028554#fig2

### Key takeway:
# looking for GSEA analysis 

- Generated differential gene table
```
# parameter
min_pct <- 0.25
min_fc <- 1.5
max_q <- 0.05

has_enough <- vector()
sn_deg_list <- vector()
sc_deg_list <- vector()
  
  
  # Run DEA for this cell type for singlecell data
sc_subset <- subset(data,idents="SingleCell")
Idents(sc_subset) <- sc_subset@meta.data$cell_type
sc_subset@active.assay <- "RNA"
sc_subset <- NormalizeData(sc_subset)
sc_dea <- FindMarkers(sc_subset,ident.1=cell,min.pct=min_pct,logfc.threshold=min_fc)

# set pc_diff
sc_dea['pct_diff']<- sc_dea["pct.1"]-sc_dea["pct.2"]

# filter result
sc_dea <- sc_dea[order(sc_dea$pct_diff,decreasing=TRUE),]
sc_dea <- sc_dea[sc_dea$p_val_adj<max_q,]
  
```


