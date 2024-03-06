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

###
# Optimize to a new K from previously factorized K (Higher than the new K) 
###

storing_dir = "/netscratch/dep_tsiantis/grp_laurent/tamal/2024/scRNA_projects/comparative_study_RCOg_at_sam/Analysis_with_default_parameters"

optimize_to_new_k(optimized_liger_object = "/netscratch/dep_tsiantis/grp_laurent/tamal/2024/scRNA_projects/comparative_study_RCOg_at_sam/Analysis_with_default_parameters/01_Data_preprocessing_and_integration/integrated_RCOg_at_sam_liger.RData", 
                  choice_of_new_K = c(30:49),
                  store_dir = storing_dir, 
                  store_folder = "Liger_objects")
