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
