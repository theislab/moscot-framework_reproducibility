print("starting")

library(vctrs)
library(stringr)
library(anndata)
library(Signac)
library(Seurat)
library(anndata)
library(ggplot2)
library(tidyverse)
library(BSgenome.Mmusculus.UCSC.mm10)


seurat <- readRDS("/lustre/groups/ml01/workspace/moscot_paper/pancreas_revision/seurat_endocrine_annotated.rds")
BSgenome.Mmusculus.UCSC.mm10.renamed <- renameSeqlevels(BSgenome.Mmusculus.UCSC.mm10, value = str_replace(str_replace(seqnames(BSgenome.Mmusculus.UCSC.mm10), pattern = "chr", replacement = ""), pattern = "M", replacement = "MT"))
seurat <- RegionStats(seurat, genome = BSgenome.Mmusculus.UCSC.mm10.renamed, assay='atac')

seurat <- LinkPeaks(
  object = seurat,
  peak.assay = "atac",
  expression.assay = "scran",
  distance = 1e+07,
)

saveRDS(seurat, "/lustre/groups/ml01/workspace/moscot_paper/pancreas_revision/seurat_endocrine_with_links_long_range.rds")

