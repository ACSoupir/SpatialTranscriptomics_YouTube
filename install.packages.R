renv::activate()
install.packages("Seurat") #had to run 2x
#restart R session
install.packages("igraph", type = "binary") #compilation fails
#retry Seurat
install.packages("Seurat")
#need to update rlang for purrr
install.packages("rlang")
#retry Seurat
install.packages("Seurat")
#devtools for github installs
install.packages("devtools")
#install biocmanager
install.packages("BiocManager")
#install dependencies for InSituType
BiocManager::install(c('SingleCellExperiment', 'SummarizedExperiment', 'sparseMatrixStats'))
install.packages("lsa")

#install insitutype for phenotyping
#reinstall R 4.1.2
#delete rlang
install.packages("rlang", dependencies = TRUE)
#cant get to work at the moment
devtools::install_github("https://github.com/Nanostring-Biostats/InSituType")
#worked?