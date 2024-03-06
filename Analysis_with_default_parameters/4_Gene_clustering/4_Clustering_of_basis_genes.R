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

###############################
###
# Convert the liger object to a seurat object to perform downstream analysis
###

storing_dir = "/netscratch/dep_tsiantis/grp_laurent/tamal/2024/scRNA_projects/comparative_study_RCOg_at_sam/Analysis_with_default_parameters"

cluster_basis_genes(seurat_object_dir = "/netscratch/dep_tsiantis/grp_laurent/tamal/2024/scRNA_projects/comparative_study_RCOg_at_sam/Analysis_with_default_parameters/Seurat_objects/seurat_object_with_all_factorized_K.RData", 
                    reduction_name_pattern = "^inmf",
                    store_data = TRUE,
                    store_dir = storing_dir,
                    store_folder = "Basis_objects")
