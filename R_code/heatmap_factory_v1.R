################################################################################
#
# heatmap factory
#
################################################################################

library(ComplexHeatmap)
library(stringr)
library(pals)

### data management ------------------------------------------------------------
# load data
## counts matrix; row names is genes, column names is samples
counts_matrix <- readRDS("~/Documents/BFX_proj/heatmap_factory/_input/Chow_PNAS_lognormcounts.rds")
## metadata for annotation; row names is samples
metadata <- readRDS("~/Documents/BFX_proj/heatmap_factory/_input/Chow_PNAS_meta_med.rds")
metadata_columns <- c("treatment", "sex")
## set gene lists; list of gene sets
gene_sets_ <- readRDS("~/Documents/BFX_proj/heatmap_factory/_input/PCA_results.RDS")[["PC1 Reactome_GSEA"]]
gene_sets <- list()
for(i in rownames(gene_sets_)){
  gene_sets[[gene_sets_[i, "Description"]]] <- unlist(str_split(gene_sets_[i, "core_enrichment_symbol"], "/"))
}
gene_sets <- gene_sets[1:5]
################################################################################

### heatmap factory function ---------------------------------------------------
heatmap_factory <- function(counts_matrix = counts_matrix,
                            metadata = metadata,
                            metadata_columns = metadata_columns,
                            gene_sets = gene_sets,
                            destination = NA){
  # check destination argument
  if(is.na(destination)){
    stop("no destination folder selected")
  }
  # check counts columns and metadata rows are in agreement
  if(!all(colnames(counts_matrix) == rownames(metadata))){
    stop("counts columns and metadata rows do not match")
  }
  # select metadata columns
  metadata_ <- metadata[, metadata_columns]
  # top annotation
  ## color legend 
  heat_cols_ <- list()
  ## legend text parameters
  legend_param_ <- list()
  for(i in colnames(metadata_)){
    heat_cols_[[i]] <- setNames(alphabet(length(unique(metadata_[, i]))),
                                unique(metadata_[, i]))
    legend_param_[[i]] <- list(title_gp = gpar(fontsize = 8),
                               labels_gp = gpar(fontsize = 8))
  }
  col_ann_ <- HeatmapAnnotation(df = metadata_,
                                col = heat_cols_,
                                annotation_legend_param = legend_param_)
  # for loop across all gene sets
  for(i in names(gene_sets)){
    genes_ <- gene_sets[[i]]
    counts_ <- counts_matrix[genes_[genes_ %in% rownames(counts_matrix)], ]
    heat_counts_ <- t(scale(t(counts_)))
    #heat_counts_[heat_counts_ < -2] <- -2
    #heat_counts_[heat_counts_ > 2] <- 2
    title_ <- str_replace_all(i, "/", "_")
    
    set.seed(415); hm_ <- Heatmap(heat_counts_,
                                  column_title = title_,
                                  col = colorRamp2(c(-2, 0, 2), c("purple", "black", "yellow")),
                                  top_annotation = col_ann_,
                                  name = "z-score",
                                  show_row_names = F,
                                  show_column_names = F,
                                  height = unit(14, "cm"),
                                  width = unit(12.5, "cm"),
                                  heatmap_legend_param = list(title_gp = gpar(fontsize = 8), 
                                                              labels_gp = gpar(fontsize = 8)))
    
    png(paste0(destination, title_, ".png"), height = 500, width = 500)
    draw(hm_)
    dev.off()
  }
}
################################################################################

# run heatmap factory
heatmap_factory(counts_matrix = counts_matrix,
                metadata = metadata,
                metadata_columns = metadata_columns,
                gene_sets = gene_sets,
                destination = "~/Documents/BFX_proj/heatmap_factory/_output/")
