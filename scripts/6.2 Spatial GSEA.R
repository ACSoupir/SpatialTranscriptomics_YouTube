
# Libraries ---------------------------------------------------------------

library(tidyverse)
library(pbmcapply)
#New package!!
#BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)


# Data --------------------------------------------------------------------

meta = readRDS("data/insitutyped_metadata.rds")
meta$unique_fov = paste0("slide", meta$slide_ID_numeric, "_fov", meta$fov)
all_fov = readRDS("data/all_fov_gsea.rds")


# Variables for enrichments -----------------------------------------------
#parameters for running the spatial enrichment
cores = 8
num_sds = 1
min_cells = 20
reps = 1000

all_cell_enrichments = lapply(all_fov, function(x){
  #extract cell level scores
  x$Cell_Sets %>%
    #transpose to make cells the rows and gene sets columns
    t() %>%
    #convert matrix back to data frame
    data.frame(check.names = FALSE)
}) %>%
  #make the list output a data frame
  do.call(bind_rows, .) %>%
  #make rownames column for cell IDs
  rownames_to_column("cell_ID")

#get unique FOVs
slide_fovs = meta$unique_fov %>% unique()



# Spatial Enrichment ------------------------------------------------------

results = pbmclapply(setNames(slide_fovs, slide_fovs), function(fov){
  cat(fov, "\n")
  #subset meta data and join only enrichment for those cells
  slide_meta = meta %>% 
    filter(unique_fov == !!fov) %>%
    select(x_FOV_px, y_FOV_px, cell_ID)
  slide_data = slide_meta %>%
    left_join(all_cell_enrichments, by = 'cell_ID')
  dist_mat = as.matrix(dist(slide_data[,1:2]))
  #iterate over gene sets
  pvals = lapply(setNames(gene_sets, gene_sets), function(gs){
    #select columns for 
    tmp_dat = slide_data[,c("x_FOV_px", 'y_FOV_px', gs)]
    #metrics
    m = mean(slide_data[[gs]], na.rm = TRUE)
    sd = sd(slide_data[[gs]], na.rm = TRUE)
    high = which(tmp_dat[[gs]] > (m + sd))
    rs = length(high)
    if(rs < min_cells)
      return(NA)
    #observed
    obs = sum(dist_mat[high, high])/2
    #permutations
    perms = sapply(1:reps, function(rep){
      cs = sample(1:nrow(tmp_dat), rs, replace = FALSE)
      sum(dist_mat[cs, cs])/2
    })
    #return p value
    return(sum(perms<obs)/1000)
  }) %>%
    #make vector
    do.call(c, .)
}, mc.cores = 8, mc.preschedule = FALSE)

#collapse vectors to a table per sample
results_simplified = lapply(names(results), function(samp){
  results[[samp]] %>% t() %>% data.frame(check.names = FALSE) %>% mutate(unique_fov = samp)
}) %>% 
  do.call(bind_rows, .)

#save output
#saveRDS(results_simplified, "data/spatialGSEA.rds")
results_simplified = readRDS("data/spatialGSEA.rds")


# Create Heatmap ----------------------------------------------------------

mat = results_simplified %>% 
  column_to_rownames("unique_fov") %>%
  as.matrix() %>%
  t()
mat[is.na(mat)] = 1
mat = -log2(mat + 0.0001)

Heatmap(mat, row_names_max_width = max_text_width(rownames(mat)), 
        column_title = "STenrich - Spatial Gene Set Enrichment Analysis",
        show_column_names = FALSE)





