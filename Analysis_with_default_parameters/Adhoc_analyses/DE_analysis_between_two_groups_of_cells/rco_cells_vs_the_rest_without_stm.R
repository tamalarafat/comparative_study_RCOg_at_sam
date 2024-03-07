# Specify the projects dir containing all the functions of scRNA-seq analysis pipeline
projects_dir <- "/home/ytamal2/Documents/2023/PhD_projects_Yasir/scExplorer/Functions/"

# Functions - Markers identification
list_1 <- list.files(paste0(projects_dir, "Library_handler"), pattern = "*.R$", full.names = TRUE) 
sapply(list_1, source, .GlobalEnv)

# Functions - Data manipulation
list_2 <- list.files(paste0(projects_dir, "Functions_matrix_manipulation"), pattern = "*.R$", full.names = TRUE) 
sapply(list_2, source, .GlobalEnv)

# Functions - Library and packages handler
list_3 <- list.files(paste0(projects_dir, "Functions_marker_identification"), pattern = "*.R$", full.names = TRUE) 
sapply(list_3, source, .GlobalEnv)

# Storing directory
res_dir = "/netscratch/dep_tsiantis/grp_laurent/tamal/2024/scRNA_projects/comparative_study_RCOg_at_sam/Analysis_with_default_parameters/Analysis_outputs"

load("/netscratch/dep_tsiantis/grp_laurent/tamal/2024/scRNA_projects/comparative_study_RCOg_at_sam/Analysis_with_default_parameters/Seurat_objects/seurat_object_of_K_50.RData")

Idents(integrated.data) <- "RNA_snn_res.0.2"

stm_cells = WhichCells(object = integrated.data, expression = AT1G62360 > 0)

rco_cells = WhichCells(object = integrated.data, expression = AT5G67651 > 0)

cells_to_keep = rownames(integrated.data@meta.data)[!rownames(integrated.data@meta.data) %in% stm_cells]

integrated.data$cell_ID = rownames(integrated.data@meta.data)


# Remove all the stm expressing cells before performing differential expression analysis between RCO-expressing cells and the rest

# Subset the seurat object
integrated.data = subset(integrated.data, subset = cell_ID %in% cells_to_keep)

cells_to_compare = cells_to_keep[!cells_to_keep %in% rco_cells]

rco_markers = FindMarkers(integrated.data, ident.1 = rco_cells, ident.2 = cells_to_compare, test.use = "wilcox", only.pos = FALSE, logfc.threshold = 0.001, min.pct = 0.001)

save(rco_cells_markers, file = "rco_cells_markers.RData")