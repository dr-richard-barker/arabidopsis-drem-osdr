#!/usr/bin/env Rscript
# Widen the auto-decoder lever from the 4,000-HVG signature matrix to all genes.
#
# WHY: 06_weight_tf_prior.py can only re-weight an edge whose TF or target appears in
# the atlas signature matrix. That matrix holds 4,000 highly variable genes, so roughly
# a quarter of the prior's edges are atlas-informed. Genes outside the HVG set are
# cell-type-invariant by construction, so passing their edges through unchanged is
# defensible — but it is an argument, not a measurement. This script replaces the
# argument with the measurement by computing per-cluster means for EVERY gene.
#
# Output: cell_type_signatures_full.csv next to the atlas. 06_weight_tf_prior.py prefers
# that file automatically when it exists; nothing else needs changing.
#
#   Rscript scripts/build_full_signatures.R [/path/to/tropism_atlas]
#
# Requires Seurat and Matrix, and roughly 32 GB of RAM for the 3.2 GB atlas object.
# Derived from tropism_atlas/build_signature_matrix.R, with the HVG restriction removed.

suppressPackageStartupMessages({ library(Seurat); library(Matrix) })

args  <- commandArgs(trailingOnly = TRUE)
ATLAS_DIR <- if (length(args) >= 1) args[1] else file.path(Sys.getenv("HOME"), "Documents", "tropism_atlas")
ATLAS <- file.path(ATLAS_DIR, "GSE226097_global_integration_221009.rds")
OUT   <- file.path(ATLAS_DIR, "cell_type_signatures_full.csv")
ASSAY <- "RNA"; LAYER <- "data"; CLUSTER <- "orig.cluster"

if (!file.exists(ATLAS)) stop("atlas object not found: ", ATLAS)

cat("loading atlas (3.2 GB, this takes a few minutes) ...\n"); t0 <- Sys.time()
obj <- readRDS(ATLAS)
DefaultAssay(obj) <- ASSAY
cat(sprintf("loaded %d genes x %d cells in %.0fs\n", nrow(obj), ncol(obj),
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))

# Every gene, not just the variable ones — that is the entire point of this script.
expr <- GetAssayData(obj, assay = ASSAY, layer = LAYER)

clu <- as.character(obj@meta.data[[CLUSTER]])
clusters <- sort(unique(clu))
cat(sprintf("computing per-cluster means over %d clusters ...\n", length(clusters)))

ind <- sparse.model.matrix(~ 0 + factor(clu, levels = clusters))  # cells x clusters
colnames(ind) <- clusters
counts_per <- Matrix::colSums(ind)
sig <- as.matrix(expr %*% ind)
sig <- sweep(sig, 2, counts_per, "/")
rownames(sig) <- rownames(expr); colnames(sig) <- clusters

cat(sprintf("full signature matrix: %d genes x %d clusters\n", nrow(sig), ncol(sig)))
write.csv(sig, OUT)
cat("wrote", OUT, sprintf("(%.1f MB)\n", file.info(OUT)$size / 1e6))
cat("06_weight_tf_prior.py will now prefer this file over the HVG matrix.\n")
