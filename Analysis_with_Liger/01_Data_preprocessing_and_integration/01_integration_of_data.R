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
# Load data tables
###

###
# Protoplasting induced genes
###

# Load the table containing the list of protoplasting-induced genes.
PP_genes_table = read.csv("/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/Protoplasting_genes/Col0.leaf.vs.protoplast_table_2pseudorep_final_August_2022.csv")

# Gene IDs - protoplasting-induced genes
PP_genes = PP_genes_table$GeneID

###
# WT A. thaliana
###

# Load data - WT OX 1st Experiment (leaf 5 and 6)
SAM_data_1E <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_RCOg_1_SAM_final/filtered_feature_bc_matrix/", gene.column = 1)

# Load data - WT OX 2nd Experiment (leaf 6 and 7)
SAM_data_2E <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_RCOg_2_SAM_final/filtered_feature_bc_matrix/", gene.column = 1)


###
# RCOg A. thaliana
###

# Data from 3rd experiment
RCO_data_3E <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_RCOg_RNA_3rd_ALL/filtered_feature_bc_matrix/")

# Data from 3rd experiment
RCO_data_5E <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_RCOg_RNA_5th_ALL_2/filtered_feature_bc_matrix/")

# Data from 3rd experiment
RCO_data_6E <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_RCOg_RNA_6th_ALL/filtered_feature_bc_matrix/")

# Data from 3rd experiment
RCO_data_7E <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_RCOg_RNA_7th_ALL_2/filtered_feature_bc_matrix/")

# Data from 3rd experiment
RCO_data_8E <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_RCOg_RNA_8th_New/filtered_feature_bc_matrix/", gene.column = 1)

# Data from 3rd experiment
RCO_data_10E <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_RCOg_RNA_10th_New/filtered_feature_bc_matrix/", gene.column = 1)


# All gene IDs - Arabidopsis Thaliana
thaliana_genes = rownames(SAM_data_1E)

# Remove protoplasting-induced genes from the total set of hirsuta genes
genes_to_keep = setdiff(thaliana_genes, PP_genes)

##### Remove the protoplasting induced genes
SAM_data_1E <- SAM_data_1E[genes_to_keep, ]
SAM_data_2E <- SAM_data_2E[genes_to_keep, ]

RCO_data_3E <- RCO_data_3E[genes_to_keep, ]
RCO_data_5E <- RCO_data_5E[genes_to_keep, ]
RCO_data_6E <- RCO_data_6E[genes_to_keep, ]
RCO_data_7E <- RCO_data_7E[genes_to_keep, ]
RCO_data_8E <- RCO_data_8E[genes_to_keep, ]
RCO_data_10E <- RCO_data_10E[genes_to_keep, ]


# Create seurat object and perform initial filtering - 
# 1. Remove genes if their expression was not detected in at least one cell out of every 500 ("min.cells"),
# 2. Remove cells if at least 200 genes were not detected to be expressed (min.features = 200),
# 3. Remove cells with a total count of more than 110000 (nCount_RNA > 110000).
# 4. Remove cells if 5% or more of the total count of a cell belongs to the mitochondiral genes.
# 5. Remove cells if 10% or more of the total count of a cell belongs to the chloroplast genes.

###
# SAM - 1 E
###

# First replicate - COL 1E - total cells 4850; filter out genes that are not detected in at least 13 cells
SAM_1E <- CreateSeuratObject(counts = SAM_data_1E, project = "SAM_1E", min.features = 200)

# Add metadata information to the seurat object
SAM_1E[[c("Species", "Replicates", "Genotype", "Tissue")]] <- c("Thaliana", "RCOg-SAM-1", "RCOg", "LS") # LS - Leaf and Apex

# Remove cells with a total count more than 110000
SAM_1E <- subset(SAM_1E, subset = nCount_RNA <= 110000)

# calculate the percentage of total counts belonging to the mitochondiral genes.
SAM_1E[["percent.mt"]] <- PercentageFeatureSet(SAM_1E, pattern = "^ATM")

# calculate the percentage of total counts belonging to the chloroplast genes.
SAM_1E[["percent.pt"]] <- PercentageFeatureSet(SAM_1E, pattern = "^ATC")

# Remove cells using the mitochondiral percentage and chloroplast percentage threshold
SAM_1E <- subset(SAM_1E, subset = percent.mt < 5 & percent.pt < 10)

# Normalize the data - log-normalization
SAM_1E <- NormalizeData(SAM_1E, verbose = FALSE)

# Find a set of highly avariable genes - 2000 HVGs
SAM_1E <- FindVariableFeatures(SAM_1E, selection.method = "vst", nfeatures = 2000)

# Extract the count table from the seurat object
SAM_W1 <- GetAssayData(SAM_1E, assay = "RNA", slot = "counts")

# Prepend a distinctive string to the cell labels 
colnames(SAM_W1) <- paste("S1", colnames(SAM_W1), sep = "_")

###
# SAM - 2 E
###

# First replicate - COL 2E - total cells 10760; filter out genes that are not detected in at least 21 cells
SAM_2E <- CreateSeuratObject(counts = SAM_data_2E, project = "SAM_2E", min.features = 200)

# Add metadata information to the seurat object
SAM_2E[[c("Species", "Replicates", "Genotype", "Tissue")]] <- c("Thaliana", "RCOg-SAM-2", "RCOg", "LS") # LS - Leaf and Apex

# Remove cells with a total count more than 110000
SAM_2E <- subset(SAM_2E, subset = nCount_RNA <= 110000)

# calculate the percentage of total counts belonging to the mitochondiral genes.
SAM_2E[["percent.mt"]] <- PercentageFeatureSet(SAM_2E, pattern = "^ATM")

# calculate the percentage of total counts belonging to the chloroplast genes.
SAM_2E[["percent.pt"]] <- PercentageFeatureSet(SAM_2E, pattern = "^ATC")

# Remove cells using the mitochondiral percentage and chloroplast percentage threshold
SAM_2E <- subset(SAM_2E, subset = percent.mt < 5 & percent.pt < 10)

# Normalize the data - log-normalization
SAM_2E <- NormalizeData(SAM_2E, verbose = FALSE)

# Find a set of highly avariable genes - 2000 HVGs
SAM_2E <- FindVariableFeatures(SAM_2E, selection.method = "vst", nfeatures = 2000)

# Extract the count table from the seurat object
SAM_W2 <- GetAssayData(SAM_2E, assay = "RNA", slot = "counts")

# Prepend a distinctive string to the cell labels 
colnames(SAM_W2) <- paste("S2", colnames(SAM_W2), sep = "_")


###
# RCOg - 3 E
###

# Total cells 4550; filter out genes that are not detected in at least 9 cells
RCO_3E <- CreateSeuratObject(counts = RCO_data_3E, project = "RCO_3E", min.features = 200)

# Add metadata information to the seurat object
RCO_3E[[c("Species", "Replicates", "Genotype", "Tissue")]] <- c("Thaliana", "RCOg-3", "RCO", "Leaf")

# Remove cells with a total count more than 110000
RCO_3E <- subset(RCO_3E, subset = nCount_RNA <= 110000)

# calculate the percentage of total counts belonging to the mitochondiral genes.
RCO_3E[["percent.mt"]] <- PercentageFeatureSet(RCO_3E, pattern = "^ATM")

# calculate the percentage of total counts belonging to the chloroplast genes.
RCO_3E[["percent.pt"]] <- PercentageFeatureSet(RCO_3E, pattern = "^ATC")

# Remove cells using the mitochondiral percentage and chloroplast percentage threshold
RCO_3E <- subset(RCO_3E, subset = percent.mt < 5 & percent.pt < 10)

# Normalize the data - log-normalization
RCO_3E <- NormalizeData(RCO_3E, verbose = FALSE)

# Find a set of highly avariable genes - 2000 HVGs
RCO_3E <- FindVariableFeatures(RCO_3E, selection.method = "vst", nfeatures = 2000)

# Extract the count matrix
RCO_G3 <- GetAssayData(RCO_3E, assay = "RNA", slot = "counts")

# Add a string (identifier) with the barcodes
colnames(RCO_G3) <- paste("R3", colnames(RCO_G3), sep = "_")


###
# RCOg - 5 E
###

# Total cells 7000; filter out genes that are not detected in at least 14 cells
RCO_5E <- CreateSeuratObject(counts = RCO_data_5E, project = "RCO_5E", min.features = 200)

# Add metadata information to the seurat object
RCO_5E[[c("Species", "Replicates", "Genotype", "Tissue")]] <- c("Thaliana", "RCOg-5", "RCO", "Leaf")

# Remove cells with a total count more than 110000
RCO_5E <- subset(RCO_5E, subset = nCount_RNA <= 110000)

# calculate the percentage of total counts belonging to the mitochondiral genes.
RCO_5E[["percent.mt"]] <- PercentageFeatureSet(RCO_5E, pattern = "^ATM")

# calculate the percentage of total counts belonging to the chloroplast genes.
RCO_5E[["percent.pt"]] <- PercentageFeatureSet(RCO_5E, pattern = "^ATC")

# Remove cells using the mitochondiral percentage and chloroplast percentage threshold
RCO_5E <- subset(RCO_5E, subset = percent.mt < 5 & percent.pt < 10)

# Normalize the data - log-normalization
RCO_5E <- NormalizeData(RCO_5E, verbose = FALSE)

# Find a set of highly avariable genes - 2000 HVGs
RCO_5E <- FindVariableFeatures(RCO_5E, selection.method = "vst", nfeatures = 2000)

# Extract the count matrix
RCO_G5 <- GetAssayData(RCO_5E, assay = "RNA", slot = "counts")

# Add a string (identifier) with the barcodes
colnames(RCO_G5) <- paste("R5", colnames(RCO_G5), sep = "_")


###
# RCOg - 6 E
###

# Total cells 4550; filter out genes that are not detected in at least 16 cells
RCO_6E <- CreateSeuratObject(counts = RCO_data_6E, project = "RCO_6E", min.features = 200)

# Add metadata information to the seurat object
RCO_6E[[c("Species", "Replicates", "Genotype", "Tissue")]] <- c("Thaliana", "RCOg-6", "RCO", "Leaf")

# Remove cells with a total count more than 110000
RCO_6E <- subset(RCO_6E, subset = nCount_RNA <= 110000)

# calculate the percentage of total counts belonging to the mitochondiral genes.
RCO_6E[["percent.mt"]] <- PercentageFeatureSet(RCO_6E, pattern = "^ATM")

# calculate the percentage of total counts belonging to the chloroplast genes.
RCO_6E[["percent.pt"]] <- PercentageFeatureSet(RCO_6E, pattern = "^ATC")

# Remove cells using the mitochondiral percentage and chloroplast percentage threshold
RCO_6E <- subset(RCO_6E, subset = percent.mt < 5 & percent.pt < 10)

# Normalize the data - log-normalization
RCO_6E <- NormalizeData(RCO_6E, verbose = FALSE)

# Find a set of highly avariable genes - 2000 HVGs
RCO_6E <- FindVariableFeatures(RCO_6E, selection.method = "vst", nfeatures = 2000)

# Extract the count matrix
RCO_G6 <- GetAssayData(RCO_6E, assay = "RNA", slot = "counts")

# Add a string (identifier) with the barcodes
colnames(RCO_G6) <- paste("R6", colnames(RCO_G6), sep = "_")


###
# RCOg - 7 E
###

# Total cells 4550; filter out genes that are not detected in at least 9 cells
RCO_7E <- CreateSeuratObject(counts = RCO_data_7E, project = "RCO_7E", min.features = 200)

# Add metadata information to the seurat object
RCO_7E[[c("Species", "Replicates", "Genotype", "Tissue")]] <- c("Thaliana", "RCOg-7", "RCO", "Leaf")

# Remove cells with a total count more than 110000
RCO_7E <- subset(RCO_7E, subset = nCount_RNA <= 110000)

# calculate the percentage of total counts belonging to the mitochondiral genes.
RCO_7E[["percent.mt"]] <- PercentageFeatureSet(RCO_7E, pattern = "^ATM")

# calculate the percentage of total counts belonging to the chloroplast genes.
RCO_7E[["percent.pt"]] <- PercentageFeatureSet(RCO_7E, pattern = "^ATC")

# Remove cells using the mitochondiral percentage and chloroplast percentage threshold
RCO_7E <- subset(RCO_7E, subset = percent.mt < 5 & percent.pt < 10)

# Normalize the data - log-normalization
RCO_7E <- NormalizeData(RCO_7E, verbose = FALSE)

# Find a set of highly avariable genes - 2000 HVGs
RCO_7E <- FindVariableFeatures(RCO_7E, selection.method = "vst", nfeatures = 2000)

# Extract the count matrix
RCO_G7 <- GetAssayData(RCO_7E, assay = "RNA", slot = "counts")

# Add a string (identifier) with the barcodes
colnames(RCO_G7) <- paste("R7", colnames(RCO_G7), sep = "_")


###
# RCOg - 8 E
###

# Total cells 4550; filter out genes that are not detected in at least 9 cells
RCO_8E <- CreateSeuratObject(counts = RCO_data_8E, project = "RCO_8E", min.features = 200)

# Add metadata information to the seurat object
RCO_8E[[c("Species", "Replicates", "Genotype", "Tissue")]] <- c("Thaliana", "RCOg-8", "RCO", "Leaf")

# Remove cells with a total count more than 110000
RCO_8E <- subset(RCO_8E, subset = nCount_RNA <= 110000)

# calculate the percentage of total counts belonging to the mitochondiral genes.
RCO_8E[["percent.mt"]] <- PercentageFeatureSet(RCO_8E, pattern = "^ATM")

# calculate the percentage of total counts belonging to the chloroplast genes.
RCO_8E[["percent.pt"]] <- PercentageFeatureSet(RCO_8E, pattern = "^ATC")

# Remove cells using the mitochondiral percentage and chloroplast percentage threshold
RCO_8E <- subset(RCO_8E, subset = percent.mt < 5 & percent.pt < 10)

# Normalize the data - log-normalization
RCO_8E <- NormalizeData(RCO_8E, verbose = FALSE)

# Find a set of highly avariable genes - 2000 HVGs
RCO_8E <- FindVariableFeatures(RCO_8E, selection.method = "vst", nfeatures = 2000)

# Extract the count matrix
RCO_G8 <- GetAssayData(RCO_8E, assay = "RNA", slot = "counts")

# Add a string (identifier) with the barcodes
colnames(RCO_G8) <- paste("R8", colnames(RCO_G8), sep = "_")


###
# RCOg - 10 E
###

# Total cells 4550; filter out genes that are not detected in at least 9 cells
RCO_10E <- CreateSeuratObject(counts = RCO_data_10E, project = "RCO_10E", min.features = 200)

# Add metadata information to the seurat object
RCO_10E[[c("Species", "Replicates", "Genotype", "Tissue")]] <- c("Thaliana", "RCOg-10", "RCO", "Leaf")

# Remove cells with a total count more than 110000
RCO_10E <- subset(RCO_10E, subset = nCount_RNA <= 110000)

# calculate the percentage of total counts belonging to the mitochondiral genes.
RCO_10E[["percent.mt"]] <- PercentageFeatureSet(RCO_10E, pattern = "^ATM")

# calculate the percentage of total counts belonging to the chloroplast genes.
RCO_10E[["percent.pt"]] <- PercentageFeatureSet(RCO_10E, pattern = "^ATC")

# Remove cells using the mitochondiral percentage and chloroplast percentage threshold
RCO_10E <- subset(RCO_10E, subset = percent.mt < 5 & percent.pt < 10)

# Normalize the data - log-normalization
RCO_10E <- NormalizeData(RCO_10E, verbose = FALSE)

# Find a set of highly avariable genes - 2000 HVGs
RCO_10E <- FindVariableFeatures(RCO_10E, selection.method = "vst", nfeatures = 2000)

# Extract the count matrix
RCO_G10 <- GetAssayData(RCO_10E, assay = "RNA", slot = "counts")

# Add a string (identifier) with the barcodes
colnames(RCO_G10) <- paste("R10", colnames(RCO_G10), sep = "_")

# Integration of the replicates - find anchors
anchFeatures <- SelectIntegrationFeatures(object.list = list(SAM_1E, SAM_2E, RCO_3E, RCO_5E, RCO_6E, RCO_7E, RCO_8E, RCO_10E))

fileGenerator(anchFeatures, "Seurat_HVG_standard.txt")


###
# Creating a Liger object, pre-processing, and performing integration
###

# Lets create the liger object
AT_genotypes <- createLiger(list(RS1 = SAM_W1, RS2 = SAM_W2, RG3 = RCO_G3, RG5 = RCO_G5, RG6 = RCO_G6, RG7 = RCO_G7, RG8 = RCO_G8, RG10 = RCO_G10), remove.missing = F)

# Add metadata information to the liger object for the replicates in the same way as it was added in the seurat object

# Add replicate information
AT_genotypes@cell.data$Replicates <- AT_genotypes@cell.data$dataset
AT_genotypes@cell.data$Replicates <- factor(AT_genotypes@cell.data$Replicates, levels = c("RS1", "RS2", "RG3", "RG5", "RG6", "RG7", "RG8", "RG10"), labels = c("RCOg-SAM-1", "RCOg-SAM-2", "RCOg-3", "RCOg-5", "RCOg-6", "RCOg-7", "RCOg-8", "RCOg-10"))

# Add species information
AT_genotypes@cell.data$Genotype <- str_sub(AT_genotypes@cell.data$dataset, 1, nchar(AT_genotypes@cell.data$dataset) - 1)
AT_genotypes@cell.data$Genotype <- factor(AT_genotypes@cell.data$Genotype, levels = c("RS", "RG"), labels = c("LS", "Leaf"))


AT_genotypes@norm.data <- list(RS1 = SAM_1E@assays$RNA@data, 
                               RS2 = SAM_2E@assays$RNA@data, 
                               RG3 = RCO_3E@assays$RNA@data,
                               RG5 = RCO_5E@assays$RNA@data,
                               RG6 = RCO_6E@assays$RNA@data,
                               RG7 = RCO_7E@assays$RNA@data,
                               RG8 = RCO_8E@assays$RNA@data,
                               RG10 = RCO_10E@assays$RNA@data)

# Assigning a set of highly variable genes
AT_genotypes@var.genes <- anchFeatures

# Scale the feature count
AT_genotypes <- scaleNotCenter(AT_genotypes)

# Check which datasets are we integrating
table(AT_genotypes@cell.data$dataset)

# Run liger integration - factorization of the matrices
AT_genotypes <- optimizeALS(AT_genotypes, k = 50, nrep = 10, lambda = 5)

# Quantile normalization of the data - integration in the shared space
AT_genotypes <- quantile_norm(AT_genotypes)

# Run liger implemented UMAP
AT_genotypes <- runUMAP(AT_genotypes)

Liger_object_K_50 <- AT_genotypes

# save the object
save(Liger_object_K_50, file = "integrated_RCOg_at_sam_liger.RData")

writeLines(capture.output(sessionInfo()), "Session_info_integrated_RCOg_at_sam_liger.txt")