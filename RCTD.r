library(Seurat)
library(ggplot2)
library(dplyr)
library(reticulate)
library(patchwork)
library(plyr)
library(tibble)
library(ggrepel)
library(reshape2)  
library(tidyr)
library(spacexr) 
set.seed(12)
celltypist = import("celltypist")
pandas <- import("pandas")
numpy = import("numpy")            
scanpy = import('scanpy')

main_dir <- 'HTAN/analysis/'
scRNA_data_dir <- file.path(main_dir, '/00.data/02.HTAPP_scRNA/CELLxGENE/')
slideseq_data_dir <- file.path(main_dir, '/00.data/03.HTAPP_slide_seq/CELLxGENE/CELLxGENE/')
slideseq_image_dir <- file.path(main_dir, '/00.data/03.HTAPP_slide_seq/CELLxGENE/HE/')
model <- 'HTAN/analysis/01.scRNA/05.1000_model.pkl'
sce_output_dir = file.path(main_dir, sprintf('02.HTAPP_scRNA/01.samples/%s',sample))
slideseq_output_dir = file.path(main_dir, sprintf('03.HTAPP_slide_seq/01.samples/%s',sample))
  
source('RenameGenesSeurat_v2.R')

sce_h5ad <- scanpy$read_h5ad(file.path(scRNA_data_dir, sprintf('%s_scRNA-seq.h5ad',sample)))
var = sce_h5ad$var  
sce_rds <- readRDS(file.path(scRNA_data_dir, sprintf('%s-scRNA-seq.rds',sample)))
sce <- RenameGenesSeurat_v2(sce_rds, 
                                newnames = var$feature_name,
                                gene.use = rownames(var),
                                de.assay = 'RNA')

sce$cell_type <- as.character(sce$cell_type) 
CAF_sce <- subset(sce, cell_type == 'fibroblast')

adata = CAF_sce
adata = scanpy$AnnData(X = numpy$array(t(as.matrix(adata[['RNA']]@counts))),
            obs = pandas$DataFrame(adata@meta.data),
            var = pandas$DataFrame(data.frame(gene = rownames(adata[['RNA']]@counts),
                                            row.names = rownames(adata[['RNA']]@counts))))                  
adata$raw = adata
adata = adata$raw$to_adata()  
library(dplyr)
scanpy$pp$normalize_total(adata, target_sum=1e4)
scanpy$pp$log1p(adata)
predictions = celltypist$annotate(adata, model = model)
CAF_sce = AddMetaData(CAF_sce, predictions$predicted_labels)
CAF_sce$predicted_labels = as.character(CAF_sce$predicted_labels)
sce$celltype <- sce$cell_type
sce$celltype <- as.character(sce$celltype)
sce@meta.data[match(colnames(CAF_sce),colnames(sce)), 'celltype'] <- as.character(CAF_sce$predicted_labels)
table(sce$celltype)

slide_h5ad <- scanpy$read_h5ad(file.path(slideseq_data_dir, sprintf('%s_Slide-seq.h5ad',sample)))
slide_rds <- readRDS(file.path(slideseq_data_dir, sprintf('%s.rds',sample)))

var = slide_h5ad$var  
slide_seq <- RenameGenesSeurat_v2(slide_rds, 
                                newnames = var$feature_name,
                                gene.use = rownames(var),
                                de.assay = 'RNA')
slide_seq$cell_type <- as.character(slide_seq$cell_type) 
CAF_slide_seq <- subset(slide_seq, cell_type=='fibroblast')
celltype_remove <-  table(CAF_sce$predicted_labels)[table(CAF_sce$predicted_labels) < 4]
CAF_sce <- subset(CAF_sce, predicted_labels %in% (celltype_remove |> names()), invert = TRUE)

counts <- CAF_sce[["RNA"]]@counts
counts <- counts[intersect(rownames(CAF_sce), rownames(CAF_slide_seq)),]


cluster <- as.factor(CAF_sce$predicted_labels)
names(cluster) <- colnames(CAF_sce)
nUMI <- CAF_sce$nCount_RNA
names(nUMI) <- colnames(CAF_sce)
reference <- Reference(counts, cluster, nUMI)

counts <- CAF_slide_seq@assays$RNA@counts

counts <- counts[intersect(rownames(CAF_sce), rownames(CAF_slide_seq)),]
coords <- CAF_slide_seq@reductions$spatial@cell.embeddings |> as.data.frame()

colnames(coords) <- c("x", "y")
coords[is.na(colnames(coords))] <- NULL
query <- SpatialRNA(coords, counts, colSums(counts))

all(rownames(query@counts) ==rownames(reference@counts) )
query@counts <- query@counts[rowSums(query@counts)>3,]
reference@counts <- reference@counts[rownames(query@counts),]

RCTD <- create.RCTD(query, reference, max_cores = 20, CELL_MIN_INSTANCE = 0, UMI_min = 0, counts_MIN = 0) 
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet") 

CAF_slide_seq <- AddMetaData(CAF_slide_seq, metadata = RCTD@results$results_df)
CAF_slide_seq@meta.data$first_type <- CAF_slide_seq@meta.data$first_type |> as.character()

slide_seq@meta.data$celltype <- slide_seq@meta.data$cell_type
slide_seq@meta.data[match(colnames(CAF_slide_seq),colnames(slide_seq)), 'celltype'] <- as.character(CAF_slide_seq$first_type)
metadata <- slide_seq@meta.data
RCTD_score <- RCTD@results$results_df |> data.frame()
RCTD_weights <- RCTD@results$weights |> data.frame()
RCTD_re <- cbind(RCTD_score, RCTD_weights)
