# Note: scDesign3 requre R version V.4.2.0
# R lib path: /home/gayan/R/x86_64-pc-linux-gnu-library/4.3
# Independent Rscript and python script since loading R packages are different

# Launch Command:
# Rscript /home/mzr19001/data_simulation/generateSyntheticCounts.R \
#     "RNA_countmatrix_2500cell_4_True_top1000genes" \
#     "ATAC_countmatrix_2500cell_4_True_top5000peaks" \
#     "/home/mzr19001/data_simulation/Results" \
#     "/home/mzr19001/data_simulation/Results/synthetic" \
#     "default" \
#     "/home/mzr19001/data_simulation/Results/cell_types_2500_4_True.txt" \
#     "default" \
#     "1"

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("devtools", quietly = TRUE)) 
  install.packages("devtools")
if (!requireNamespace("scDesign3", quietly = TRUE)) {
devtools::install_github("SONGDONGYUAN1994/scDesign3")
}
if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
BiocManager::install("SingleCellExperiment")
}

# Load libraries
cat(sprintf("[scReadSim] Loading R packages...\n"))
suppressPackageStartupMessages(library(scDesign3))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(SingleCellExperiment))


scMultiOmics_runSyntheticCount <- function(samplename_RNA, samplename_ATAC, directory, out_directory, n_cell_new="default", celllabel_file="default", n_cluster="default", n_cores=1){
    cat(sprintf("[scReadSim] Reading RNA count matrix %s.txt...\n", samplename_RNA))
    # Load in scRNA-seq UMI count matrix
    RNA_mat <- read.table(sprintf("%s/%s.txt", directory, samplename_RNA), sep="\t",header = FALSE,  row.names = 1) %>% as.matrix()
    cat(sprintf("[scReadSim] Reading ATAC count matrix %s.txt...\n", samplename_ATAC))
    # Load in scATAC-seq count matrix 
    ATAC_mat <- read.table(sprintf("%s/%s.txt", directory, samplename_ATAC), sep="\t",header = FALSE,  row.names = 1) %>% as.matrix()
    
    # Extract dimensionality
    n_RNA <- nrow(RNA_mat)
    n_ATAC <- nrow(ATAC_mat)

    matrix_num <- RNA_mat
    # Extract cell labels
    cat(sprintf("[scReadSim] Loading cell label file %s...\n", celllabel_file))
        full_cell_label <- unlist(read.table(celllabel_file, header=FALSE))
        if (length(full_cell_label) == ncol(matrix_num)){
            colnames(matrix_num) <- full_cell_label
        } else {
            stop("[scReadSim] Number of cell labels differs from the cell number contained in the count matrix!\n ")
        }
    start_time <- proc.time()
    # Generate sce object
    mat <- rbind(RNA_mat, ATAC_mat)
    sce <- SingleCellExperiment(list(counts = mat))
    colData(sce)$cell_type <- as.factor(full_cell_label)
    n_cell_old <- ncol(mat)
    if (n_cell_new == "default"){
        n_cell_new <- n_cell_old
    }

    # Simulate
    set.seed(2022)
    new_sce <- scdesign3(sce = sce,
                        assay_use = "counts",
                        celltype = "cell_type",
                        pseudotime = NULL,
                        spatial = NULL, 
                        other_covariates = NULL, 
                        ncell = n_cell_new, 
                        mu_formula = "cell_type", 
                        sigma_formula = "cell_type", 
                        family_use = "nb", 
                        n_cores = n_cores, 
                        usebam = FALSE, 
                        corr_formula = "cell_type", 
                        copula = "gaussian", 
                        fastmvn = FALSE, 
                        DT = TRUE, 
                        pseudo_obs = FALSE, 
                        family_set = "gauss", 
                        nonnegative = TRUE, 
                        nonzerovar = FALSE, 
                        return_model = FALSE, 
                        parallelization = "mclapply", 
                        trace = FALSE
                        )
    end_time <- proc.time()

    cat(sprintf("[scReadSim - Benchmark] Simulation completed in %s seconds.\n", end_time - start_time))

    # Write out data files
    synthetic_RNA_mat <- new_sce$new_count[1:n_RNA, ]
    synthetic_ATAC_mat <- new_sce$new_count[(1+n_RNA):nrow(mat), ]
    if (n_cell_new == n_cell_old){
      synthetic_cell_label <- full_cell_label
    } else {
      synthetic_cell_label <- new_sce$new_covariate$cell_type
    }
    cat(sprintf("[scReadSim] Writing out synthetic cell labels to %s...\n", out_directory))
    write.table(synthetic_cell_label, sprintf("%s/scMultiOmics.scDesign3Simulated.CellTypeLabel.txt", out_directory), row.names = FALSE,col.names = FALSE, quote = FALSE)
    cat(sprintf("[scReadSim] Writing out synthetic count matrices to %s...\n", out_directory))
    write.table(synthetic_RNA_mat, sprintf("%s/%s.scMultiOmics.scDesign3Simulated.RNA.txt", out_directory, samplename_RNA), sep="\t", row.names = TRUE,col.names = FALSE)
    write.table(synthetic_ATAC_mat, sprintf("%s/%s.scMultiOmics.scDesign3Simulated.ATAC.txt", out_directory, samplename_ATAC), sep="\t", row.names = TRUE,col.names = FALSE)
    cat("[scReadSim] Done.\n")
}

args <- commandArgs(trailingOnly = TRUE)

# Assign parameters from command-line arguments
samplename_RNA <- args[1]
samplename_ATAC <- args[2]
directory <- args[3]
out_directory <- args[4]
n_cell_new <- args[5]
celllabel_file <- args[6]
n_cluster <- args[7]
n_cores <- as.numeric(args[8])

# Call the function
scMultiOmics_runSyntheticCount(
    samplename_RNA = samplename_RNA,
    samplename_ATAC = samplename_ATAC,
    directory = directory,
    out_directory = out_directory,
    n_cell_new = n_cell_new,
    celllabel_file = celllabel_file,
    n_cluster = n_cluster,
    n_cores = n_cores
)