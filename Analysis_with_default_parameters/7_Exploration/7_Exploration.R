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

# Function to generate results for a seurat object with clustering results for multiple resolution parameter values
resolution_explorer <- function(seuratObject, 
                                store_dir = NULL, 
                                store_folder = "Your_choice_of_name"
){
  
  # Creating necessary storing space to store the results
  
  if (missing(store_dir)){
    store_dir = getwd()
  }
  
  if (!dir.exists(str_c(store_dir, "/", store_folder))){
    dir.create(str_c(store_dir, "/", store_folder), showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }
  
  temp_dir = str_c(store_dir, "/", store_folder, "/")
  
  # let's go over the resolution values
  temp = str_sort(names(seuratObject@meta.data)[str_detect(names(seuratObject@meta.data), pattern = "res")], numeric = TRUE)
  
  for (i in c(1:length(temp))){
    
    temp_str = substr(gsub("[^[:alnum:]]", "", temp[i]), start = nchar(gsub("[^[:alnum:]]", "", temp[i])) - 1, stop = nchar(gsub("[^[:alnum:]]", "", temp[i])))
    temp_str = str_c("RES_", gsub(pattern = "[A-Za-z]", replacement = "", temp_str))
    
    if (!dir.exists(str_c(temp_dir, "/", temp_str))){
      dir.create(str_c(temp_dir, "/", temp_str), showWarnings = TRUE, recursive = FALSE, mode = "0777")
    }
    
    temp_res_dir = str_c(temp_dir, "/", temp_str, "/")
    
    # Set the cluster identity to the seurat object
    Idents(seuratObject) <- temp[i]
    
    Idents(seuratObject) <- factor(Idents(seuratObject), levels = seq(0, length(levels(Idents(seuratObject))) - 1))
    
    replicate_proportion_per_cluster(seuratObject = seuratObject, store_dir = temp_res_dir, replicate_metadata_name = "Replicates", split_metadata_name = "Genotypes")
    
    group_proportion_per_cluster(seuratObject = seuratObject, store_dir = temp_res_dir, replicate_metadata_name = "Replicates", split_metadata_name = "Genotypes", rep_prop = FALSE)
    
    group_proportion_per_cluster(seuratObject = seuratObject, store_dir = temp_res_dir, replicate_metadata_name = "Replicates", split_metadata_name = "Genotypes", rep_prop = TRUE)
    
    feature_count_n_proportion(seuratObject = seuratObject, store_dir = temp_res_dir, gene_ID = "AT1G62360", gene_name = "STM")
    
    feature_count_n_proportion(seuratObject = seuratObject, store_dir = temp_res_dir, gene_ID = "AT5G67651", gene_name = "RCO")
    
    RDimension_plot(seuratObject = seuratObject, store_dir = temp_res_dir, store_folder = "Reduced_representation",dimension_reduction_name = "umap")
    
    # RDimension_plot(seuratObject = seuratObject, store_dir = temp_res_dir, store_folder = "Reduced_representation",dimension_reduction_name = "tsne")
    
    genes_feature_plot(seuratObject = seuratObject, store_dir = temp_res_dir, store_folder = "Genes_feature_plot", gene_IDs = c("AT1G62360", "AT5G67651"), gene_names = c("STM", "RCO"), reduction_name = "umap")
    
    genes_feature_plot(seuratObject = seuratObject, store_dir = temp_res_dir, store_folder = "Genes_feature_plot", gene_IDs = c("AT1G62360", "AT5G67651"), gene_names = c("STM", "RCO"), reduction_name = "tsne")
    
    markers_dot_plot(seuratObject = seuratObject, store_dir = temp_res_dir, marker_file = cell_type_markers, genes_ID_column = 1, genes_name_column = 2, group_clusters = TRUE)
    
    # side_by_side_proportion_comparison(seuratObject = seuratObject, store_dir = temp_res_dir, split_variable_name = "Genotypes")
    
    RDimension_split_group_visualization(seuratObject = seuratObject, store_dir = temp_res_dir, split_variable = "Genotypes")
    
    # cells_per_cluster_split_viz(seuratObject = seuratObject, store_dir = temp_res_dir, split_variable = "Genotypes")
    
    feature_expression_split_visualization(seuratObject = seuratObject, store_dir = temp_res_dir, split_variable = "Genotypes", gene_ID = "AT5G67651", gene_name = "RCO")
    
    feature_expression_split_visualization(seuratObject = seuratObject, store_dir = temp_res_dir, split_variable = "Genotypes", gene_ID = "AT1G62360", gene_name = "STM")
    
    # cluster_conserved_marker_finder(seuratObject = seuratObject, store_dir = temp_res_dir, grouping_variable = "Genotypes")
    #
    # between_groups_differentially_expressed_genes(seuratObject = seuratObject, store_dir = temp_res_dir, between_group_variable = "Genotypes")
    # 
    # group_n_cluster_DEG_finder(seuratObject = seuratObject, store_dir = temp_res_dir)
    # 
    # cluster_marker_finder(seuratObject = seuratObject, store_dir = temp_res_dir)
    
  }
  
}


# For a series of seurat objects
resolution_explorer_series_of_seurat_objects <- function(seurat_object_dir, 
                                                         file_name_pattern, 
                                                         store_dir = NULL,
                                                         store_folder = "Cells_genes_clusters",
                                                         ident_levels = NULL,
                                                         add_replicate_labels = NULL,
                                                         add_species_labels = NULL,
                                                         add_tissue_labels = NULL
){
  
  # Lets get the saved file names
  seuratObjects = str_sort(list.files(path = seurat_object_dir, pattern = file_name_pattern), numeric = TRUE)
  
  # Let's create a directory to store the GO annotation results for different factorization
  if (missing(store_dir)){
    store_dir = getwd()
  }
  
  # Let's create a directory to store the GO annotation results
  if (!dir.exists(str_c(store_dir, "/", store_folder))){
    dir.create(str_c(store_dir, "/", store_folder), showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }
  
  # storing the directory information in a temporary variable
  temp_main_dir = str_c(store_dir, "/", store_folder, "/")
  
  for (i in c(1:length(seuratObjects))){
    # Creating the file name with path location
    file_name = str_c(seurat_object_dir, "/", seuratObjects[i])
    
    # Loading the liger object and storing it to a temporary variable
    integrated.data = loadRData(file_name) # Load the RData file and assign it to a variable using the function loadRData
    
    # Dataset information
    integrated.data$Replicates <- integrated.data$orig.ident
    integrated.data$Replicates <- factor(integrated.data$Replicates, levels = c("S1", "S2", "R3", "R5", "R6", "R7", "R8", "R10"), labels = c("RCOg-SAM-1", "RCOg-SAM-2", "RCOg-Leaf-3", "RCOg-Leaf-5", "RCOg-Leaf-6", "RCOg-Leaf-7", "RCOg-Leaf-8", "RCOg-Leaf-10"))
    
    # Genotypes information
    integrated.data$Genotypes <- integrated.data$orig.ident
    integrated.data$Genotypes <- factor(integrated.data$Genotypes, levels = c("S1", "S2", "R3", "R5", "R6", "R7", "R8", "R10"), labels = c("SAM_RCOg", "SAM_RCOg", "Leaf_RCOg", "Leaf_RCOg", "Leaf_RCOg", "Leaf_RCOg", "Leaf_RCOg", "Leaf_RCOg"))
    
    # Tissue information
    integrated.data$Tissue <- integrated.data$orig.ident
    integrated.data$Tissue <- factor(integrated.data$Tissue, levels = c("S1", "S2", "R3", "R5", "R6", "R7", "R8", "R10"), labels = c("SAM", "SAM", "Leaf", "Leaf", "Leaf", "Leaf", "Leaf", "Leaf"))
    
    save(integrated.data, file = str_c(seurat_object_dir, "/", "seurat_object_of_K_", dim(integrated.data@reductions[["inmf"]])[2], ".RData"))
    
    # Create a directory to store the results of the particular factorization
    # Check if the folder exists (nested to previously created folder), if not create one to store the figures for the clustering using the coefficient matrix with K (any)
    # if (!dir.exists(str_c(temp_main_dir, "Factor_K_", ncol(integrated.data@reductions$inmf)))){
    #   dir.create(str_c(temp_main_dir, "Factor_K_", ncol(integrated.data@reductions$inmf)), showWarnings = TRUE, recursive = F, mode = "0777")
    # }
    
    # temp_factor_dir = str_c(temp_main_dir, "Factor_K_", ncol(integrated.data@reductions$inmf), "/")
    
    resolution_explorer(seuratObject = integrated.data, store_dir = temp_main_dir, store_folder = str_c("Factor_K_", ncol(integrated.data@reductions$inmf)))
    
  }
}

# Known cell type markers
cell_type_markers = read.csv("/netscratch/dep_tsiantis/grp_laurent/tamal/2022/Input_files/Additional_inputs/Known_markers_file/Curated_markers_thaliana.csv")

# Directory containing the seurat objects
seurat_files_dir = "/netscratch/dep_tsiantis/grp_laurent/tamal/2024/scRNA_projects/comparative_study_RCOg_at_sam/Analysis_with_default_parameters/Seurat_objects"

# Storing directory
storing_dir = "/netscratch/dep_tsiantis/grp_laurent/tamal/2024/scRNA_projects/comparative_study_RCOg_at_sam/Analysis_with_default_parameters/Analysis_outputs"

resolution_explorer_series_of_seurat_objects(seurat_object_dir = seurat_files_dir, 
                                             file_name_pattern = "seurat_object_of_", 
                                             store_dir = storing_dir, 
                                             store_folder = "Explore_resolutions")