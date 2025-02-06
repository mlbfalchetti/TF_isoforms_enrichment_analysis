# TF_isoforms_enrichment_analysis

## Authors: Marcelo Falchetti, Stephanie Wang

- [ ] Create a Seurat object for every TF isoform (3,548 annotated splice isoforms)

- [ ] Get clusters (transcriptional states) for every TF isoform (check code block)

- [ ] For every cluster get a pseudobulk sample (https://satijalab.org/seurat/articles/de_vignette#perform-de-analysis-after-pseudobulking)

Output: For every TF isoform a matrix of:

- [ ] ncol = number of clusters
  
- [ ] nrow = number of genes
  
- [ ] matrix' cells = got from the AggregateExpression() function

```r
seurat <- CreateSeuratObject(
  counts = counts, 
  min.cells = ncol(counts) * 0.05, 
  min.features = 0
  )
seurat[["percent.mt"]] <- PercentageFeatureSet(
  seurat, 
  pattern = "^MT-"
  )
seurat <- subset(
  seurat, 
  subset = nFeature_RNA > 400 & nFeature_RNA < 7000 & percent.mt < 10
  )
seurat <- NormalizeData(
  seurat
  )
seurat <- FindVariableFeatures(
  seurat, 
  selection.method = "vst", 
  nfeatures = 2000
  )
seurat <- ScaleData(
  seurat, 
  features = rownames(seurat)
  )
seurat <- RunPCA(
  seurat, 
  features = VariableFeatures(object = seurat)
  )
seurat <- FindNeighbors(
  seurat, 
  dims = 1:10
  )
seurat <- FindClusters(
  seurat, 
  resolution = 0.5
  )
seurat <- RunUMAP(
  seurat, 
  dims = 1:10
  )
```

```r
seurat_pseudobulk <- AggregateExpression(
  seurat, 
  assays = "RNA", 
  return.seurat = T, 
  group.by = "seurat_clusters"
  )
```

```r
list_of_data <- list(tmp_01, tmp_02, tmp_03, tmp_04, tmp_05, tmp_06, tmp_07, tmp_08, tmp_09)
tmp <- reduce(list_of_data, full_join, by = "Gene")

tmp[is.na(tmp)] <- 0

sum_cols <- grep("^Sum\\.", names(tmp), value = TRUE)

corr_matrix <- matrix(NA, nrow = length(sum_cols), ncol = length(sum_cols), dimnames = list(sum_cols, sum_cols))
for (i in seq_along(sum_cols)) {
  for (j in seq_along(sum_cols)) {
    valid_genes <- !is.na(tmp[[sum_cols[i]]]) & !is.na(tmp[[sum_cols[j]]])
    if (sum(valid_genes) > 1) {
      corr_matrix[i, j] <- cor(tmp[[sum_cols[i]]][valid_genes], tmp[[sum_cols[j]]][valid_genes], method = "pearson")
    }
  }
}
colnames(corr_matrix) <- gsub("^Sum\\.", "", colnames(corr_matrix))
rownames(corr_matrix) <- gsub("^Sum\\.", "", rownames(corr_matrix))
corr_matrix <- melt(corr_matrix)
colnames(corr_matrix) <- c("Var1", "Var2", "Correlation")
ggplot(corr_matrix, aes(x = Var1, y = Var2, fill = Correlation)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", Correlation)), size = 4) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limit = c(-1, 1), space = "Lab", name = "Pearson\nCorrelation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.title.x = element_blank(), axis.title.y = element_blank()) +
  coord_fixed()
```
