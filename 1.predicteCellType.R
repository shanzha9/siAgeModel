library(Seurat)
library(SeuratData)
library(ggplot2)
options(future.globals.maxSize = 8000 * 1024^2)


# read reference data
obj_ref <- readRDS('./scRNA-seqProcessedLabelledObject.rds')
obj_ref <- FindVariableFeatures(obj_ref, nfeatures = 2000)

# read query data
obj_query <- readRDS('./scRNA-seqProcessedObject.rds')
obj_query$orig.ident <- paste("sample_", obj_query$batch_key, sep = "")  # batch_key: column names of the batch: sample_1, sample_2, sample_3
obj_query <- NormalizeData(obj_query)

# transfer cell type labels
obj.anchors <- FindTransferAnchors(reference = obj_ref, query = obj_query, dims = 1:50,
                                        reference.reduction = "pca")
predictions <- TransferData(anchorset = obj.anchors, refdata = obj_ref$secondary_type, dims = 1:30)
obj_query <- AddMetaData(obj_query, metadata = predictions)
obj_query$secondary_type <- obj_query$predicted.id

# save the predicted cell type labels
saveRDS(obj_query, './predictedSeuratObject.rds', compress = F)

# plot the predicted cell type labels
features_secondary <- rev(c(
  "CD3G",
  "CCR7", "SELL", "LEF1", "IL7R", "TCF7",
  "CD4", "CD40LG", "GPR183",
  "AQP3", "CXCR3", "LIMS1",
  "ANXA1", "GATA3", "PLP2", "CRIP2",
  "FOXP3",
  "CD8A", "CD8B",
  "NKG7", "GNLY", "GZMK", "GZMB", "GZMA",
  "ITGAE", "HAVCR2", "TIGIT", "CTLA4",
  "CMC1", "DUSP2", "CCL5", # "HLA-DRB1",
  "ZNF683", "IKZF2",
  "SLC4A10", "TRAV1-2",
  "TRGV9", "TRDV2",
  "NCAM1", "KLRF1", "MKI67",
  "CD19", "MS4A1", "CD79A",
  "TCL1A", "IGHD", "IGHM",
  "BLK", "AIM2", "FCRL5",
  "JCHAIN", "TRAC", "IGKV4-1", "MZB1", "IGHA1",
  "CD14", "VCAN", "MX1",
  "FCGR3A", "S100A12", "MS4A7",
  "FCER1A", "CD1D",
  "SERPINF1", "LILRA4",
  "PPBP", "PF4"
))


p <- VlnPlot(obj_query, features = features_secondary, stack = T, pt.size = 0, flip = T, group.by = "secondary_type") + NoLegend()

ggsave('./featuresVlnplot.pdf', p, height = 16, width = 8.5)
