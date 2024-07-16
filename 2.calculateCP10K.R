library(Seurat)
library(stringr)

cp10k_path <- ""
if (!dir.exists(cp10k_path)) {
  dir.create(cp10k_path)
}

if (!dir.exists(paste0(cp10k_path, '/CPM_add_zero'))) {
  dir.create(paste0(cp10k_path, '/CPM_add_zero'))
}

if (!dir.exists(paste0(cp10k_path, '/CPM'))) {
  dir.create(paste0(cp10k_path, '/CPM'))
}

seurat_object <- readRDS("./predictedSeuratObject.rds")
table(seurat_object$secondary_type)
gene_list <- rownames(seurat_object@assays$RNA$counts)
sample_groups <- unique(seurat_object$orig.ident)

sample_counts_zero_expression_samplenocells <- as.data.frame(matrix(0.00,nrow = length(gene_list), ncol = 1))
rownames(sample_counts_zero_expression_samplenocells) <- as.character(gene_list)

celltypes_list <- unique(seurat_object$secondary_type)

for (i in celltypes_list) {
  cat(i, '\n')
  expression_from_cpm <- as.data.frame(matrix(,nrow = length(gene_list), ncol = 0))
  rownames(expression_from_cpm) <- as.character(gene_list)
  
  # zero expression cell data
  expression_from_cpm_add_zero <- as.data.frame(matrix(,nrow = length(gene_list), ncol = 0))
  rownames(expression_from_cpm_add_zero) <- as.character(gene_list)
  
  cache_obj <- subset(seurat_object, subset = secondary_type == i)
  cache_obj <- SetIdent(cache_obj, 
                        cells = rownames(cache_obj@meta.data),
                        value = cache_obj$orig.ident)
  
  for (j in sample_groups) {
    cat(j, '\n')
    if (j %in% unique(cache_obj$orig.ident)) {
      cache_obj_sample <- subset(cache_obj, subset = orig.ident == j)
      if (length(cache_obj_sample$orig.ident) > 1) {
        all_sample_counts <- as.data.frame(as.matrix(cache_obj_sample@assays$RNA$counts))
        sample_counts_normal <- apply(all_sample_counts, 2, function(x) (x/sum(x))*10000)
        sample_counts_normal_mean <- as.data.frame(apply(sample_counts_normal, 1, mean))
        colnames(sample_counts_normal_mean) <- j
        expression_from_cpm <- cbind(expression_from_cpm, sample_counts_normal_mean)
      } else {
        # there is only one cell is a certain sample
        # seurat4 or 5
        all_sample_counts <- as.data.frame(as.matrix(cache_obj_sample@assays[["RNA"]]@counts))
        # all_sample_counts <- as.data.frame(as.matrix(cache_obj_sample@assays[["RNA"]]@layers[["counts"]]))
        sample_counts_normal <- apply(all_sample_counts, 2, function(x) (x/sum(x))*10000)
        colnames(sample_counts_normal) <- j
        rownames(sample_counts_normal) <- gene_list
        expression_from_cpm <- cbind(expression_from_cpm, sample_counts_normal)
      }
    } else {
      # there is no cell is a certain sample
      cat('No cells in this sample, insert zero expression...\n')
      sample_counts_normal_mean <- sample_counts_zero_expression_samplenocells
    }
    colnames(sample_counts_normal_mean) <- j
    expression_from_cpm_add_zero <- cbind(expression_from_cpm_add_zero, sample_counts_normal_mean)
  }
  
  file_name <- paste0('./CPM_add_zero/', i, '_CPM_all_samples.csv')
  write.csv(expression_from_cpm_add_zero, file_name, row.names = T)
  file_name <- paste0('./CPM/', i, '_CPM_all_samples.csv')
  write.csv(expression_from_cpm, file_name, row.names = T)
}

cp10k_list <- list.files("~/data/age_review/7.ageModel/review_validation/cp10k_nm/CPM_add_zero", full.names = T,)
validation <- as.data.frame(matrix(, ncol = length(1), nrow = 0))

for (i in cp10k_list) {
  cat(i, '\n')
  celltype <- str_split(i, pattern = '/')[[1]][11] # 11 is the index of the celltype, need to change according to the path. Change it to ensure that that the celltype is in the format of '"B_Atypical_Memory_CPM_all_samples.csv"'
  celltype <- str_split(celltype, pattern = '_')[[1]]
  celltype <- celltype[1:(length(celltype) - 3)]
  celltype <- paste(celltype, collapse = '_')
  cat(celltype, '\n')
  data_celltype <- read.csv(i, row.names = 1, check.names = F)
  colnames(data_celltype) <- str_replace_all(colnames(data_celltype), '-', '_')
  colnames(data_celltype) <- str_replace_all(colnames(data_celltype), '[.]', '_')
  rownames(data_celltype) <- paste(celltype, rownames(data_celltype), sep = ".")
  validation <- rbind(validation, data_celltype)
}

write.csv(validation, './cp10k.csv', row.names = T)
