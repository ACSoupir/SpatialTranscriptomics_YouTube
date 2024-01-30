
# Libraries ---------------------------------------------------------------

library(tidyverse)
library(msigdbr)
library(boot)
library(parallel)
library(ggpubr)
#New Libraries!!
#install as root
#apt-get update
#apt-get upgrade
#apt-get install -y libproj-dev
#apt-get install -y libgdal-dev
#apt-get install -y libudunits2-dev libgeos-dev
#devtools::install_github("r-spatial/spdep") must come from github - not cran
library(spdep)

# Data --------------------------------------------------------------------
#sample data
genes = readRDS("data/sct_gex.rds")
meta = readRDS("data/insitutyped_metadata.rds")
meta$unique_fov = paste0("slide", meta$slide_ID_numeric, "_fov", meta$fov)

#gene sets
all_gene_sets = msigdbr(species = "Homo sapiens", category = "H")
msigdbr_list = split(x = all_gene_sets$gene_symbol, f = all_gene_sets$gs_name)
ligand_pairs_available = readRDS(url("https://github.com/ZJUFanLab/CellTalkDB/raw/master/database/human_lr_pair.rds"))

#which ligand-receptors are in EMT
EMT_genes = intersect(msigdbr_list$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION, row.names(genes))
EMT_lr_pairs = ligand_pairs_available %>%
  filter(ligand_gene_symbol %in% EMT_genes, receptor_gene_symbol %in% EMT_genes)


# Bivariate Spatial Autocorrelation ---------------------------------------
#checking low number of nearby cells which makes most sense for communication
#doesn't make sense to check correlation for 100 nearest neighbors
nearest_neighbors = 3

#samples
fovs = unique(meta$unique_fov)
#lapply!!
morans_list = mclapply(setNames(fovs, fovs), function(fov){
  #subset to get the metadata and gene expression for the cells
  fov_meta = meta %>% dplyr::filter(unique_fov == !!fov)
  rm(meta)
  fov_expression = genes[,colnames(genes) %in% fov_meta$cell_id] %>%
    t() %>% data.frame(check.names = FALSE) %>%
    rownames_to_column("cell_id") %>%
    full_join(fov_meta %>%
                select(cell_id, x_FOV_px, y_FOV_px), 
              by = "cell_id") #silent the join
  rm(genes)
  #convert our gene expression to an `sf` object specifying the x and y locs
  spdf = st_as_sf(fov_expression,
                  coords = c("x_FOV_px", "y_FOV_px"))
  #using the new sf object, find the k-nearest neighbors
  knn = knearneigh(st_coordinates(spdf), k = nearest_neighbors)
  #convert the knn list to a neighbor list
  knn.nb = knn2nb(knn)
  #create a weight list from the neighbor list
  knn.nb.listw = nb2listw(knn.nb, style = "B")
  #style here is binary - either 1 or 0
  #can do row with "W" which will just be 0.33 [1/nearest_neighbor]
  
  lr_moranI_out = lapply(EMT_lr_pairs$lr_pair, function(pair){
    #grab the ligand and receptor genes making sure they exist
    ligand = EMT_lr_pairs %>% filter(lr_pair == !!pair) %>% pull(ligand_gene_symbol)
    receptor = EMT_lr_pairs %>% filter(lr_pair == !!pair) %>% pull(receptor_gene_symbol)
    #calculate bivariate morans I
    morans_I = moran_bv(spdf[[ligand]], 
                        spdf[[receptor]],
                        knn.nb.listw,
                        nsim = 1000)
    #calculate bootstrap confidence interval
    morans_ci = boot.ci(morans_I, conf = 0.95, type = 'basic')
    return(data.frame(original = morans_I$t0, #observed Morans I
                      lr_pair = pair, #ligand receptor pair
                      ci_low = morans_ci$basic[4], #low confidence interval
                      ci_high = morans_ci$basic[5])) #high confidence interval
  }) %>% 
    #bind to make dataframe wiht each row as l-r pair
    do.call(bind_rows, .) %>%
    #add column for the FOV
    mutate(unique_fov = fov)
}, mc.cores = 24)

saveRDS(morans_list, "data/ligand_receptor_CellTalkDB-knn3_EMT.rds")


# Find sample to plot -----------------------------------------------------

morans_df = do.call(bind_rows, morans_list) %>%
  arrange(desc(original))
head(morans_df)

meta_plot = meta %>%
  filter(unique_fov == "slide1_fov217")
pp = meta_plot %>%
  ggplot() + 
  geom_point(aes(x = x_FOV_px, y = y_FOV_px)) +
  coord_equal() + 
  theme_bw()
pp
#ligand expression SPP1
lpp = genes[,colnames(genes) %in% meta_plot$cell_id] %>%
  t() %>% data.frame(check.names = FALSE) %>%
  rownames_to_column("cell_id") %>%
  full_join(meta_plot %>%
              select(cell_id, x_FOV_px, y_FOV_px), 
            by = "cell_id") %>%
  ggplot() + 
  geom_point(aes(x = x_FOV_px, y = y_FOV_px, color = COL7A1)) +
  coord_equal() + 
  theme_bw()
lpp
#ligand expression ITGB3
rpp = genes[,colnames(genes) %in% meta_plot$cell_id] %>%
  t() %>% data.frame(check.names = FALSE) %>%
  rownames_to_column("cell_id") %>%
  full_join(meta_plot %>%
              select(cell_id, x_FOV_px, y_FOV_px), 
            by = "cell_id") %>%
  ggplot() + 
  geom_point(aes(x = x_FOV_px, y = y_FOV_px, color = ITGB1)) +
  coord_equal() + 
  theme_bw()
rpp
  