print("starting")

library(Seurat)
library(Signac)
library(ggplot2)
library(tidyverse)
#library(EnsDb.Mmusculus.v102)
#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

#BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")
#library(EnsDb.Mmusculus.v102)
#library(pastecs)
library(stringr)
library(BSgenome.Mmusculus.UCSC.mm10)
library(universalmotif)
library("chromVAR")


seurat <- readRDS("/lustre/groups/ml01/workspace/moscot_paper/pancreas_revision/seurat_endocrine_annotated_with_motifs.rds")

BSgenome.Mmusculus.UCSC.mm10.renamed <- renameSeqlevels(BSgenome.Mmusculus.UCSC.mm10, value = str_replace(str_replace(seqnames(BSgenome.Mmusculus.UCSC.mm10), pattern = "chr", replacement = ""), pattern = "M", replacement = "MT"))

seurat <- RunChromVAR(
  object = seurat,
  genome = BSgenome.Mmusculus.UCSC.mm10.renamed
)

write.csv(seurat[["chromvar"]]@data, "/lustre/groups/ml01/workspace/moscot_paper/pancreas_revision/cisBP_chromvar_annotations_reduced.csv")

saveRDS(seurat, "/lustre/groups/ml01/workspace/moscot_paper/pancreas_revision/cisBP_seurat_with_chromVAR_reduced.rds")