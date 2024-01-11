
# Libraries ---------------------------------------------------------------

library(tidyverse)
library(Seurat)

#new libraries
#BiocManager::install("GSVA")
library(GSVA)
#install.packages("msigdbr")
library(msigdbr)
#install.packages("pbmcapply")
library(pbmcapply)


# Read in Data ------------------------------------------------------------
#seurat cosmx data
# object = readRDS("data/HFC_insitutype.rds")
# genes = object@assays$SCT@scale.data %>%
#   as.matrix()
# #save gene expression
# saveRDS(genes, "data/sct_gex.rds")
#read in saved data
genes = readRDS("data/sct_gex.rds")
meta = readRDS("data/insitutyped_metadata.rds")

#gene sets
all_gene_sets = msigdbr(species = "Homo sapiens") %>%
  filter(gs_cat == "H") %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()
# gene_set_list = split(all_gene_sets$gene_symbol, all_gene_sets$gs_name)
# saveRDS(gene_set_list, "data/hallmark_geneset_list.rds")
gene_set_list = readRDS("data/hallmark_geneset_list.rds")

# Run Gene Set Enrichment -------------------------------------------------
#identify the unique samples
meta$unique_fov = paste0("slide", meta$slide_ID_numeric, "_fov", meta$fov)
fovs = unique(meta$unique_fov)
#conduct the gene set enrichment per fov for all hallmark gene sets
#NOTE: Must stay under memory limits or will return errors!
all_fov = pbmclapply(setNames(fovs, fovs), function(fov_id){
  #get fov cells and genes associated with those cells
  fov_meta = meta %>% filter(unique_fov == fov_id)
  fov_genes = genes[,colnames(genes) %in% fov_meta$cell_ID]
  #run the enrichment
  cell_gs = gsva(fov_genes, gene_set_list)
  #get the unique cell clusters from insutype
  phen = unique(fov_meta$unsup_insitu)
  phen = setNames(phen, phen)
  #average enrichment by cell type
  av_gs = lapply(phen, function(x){
    rowMeans(cell_gs[,fov_meta$cell_ID[fov_meta$unsup_insitu == x]] %>% data.frame())
  }) %>% do.call(bind_cols, .)
  row.names(av_gs) = row.names(cell_gs) %>%
    gsub("HALLMARK_", "", .) %>%
    gsub("_", " ", .)
  #return list of cell and phenotype level GSEA
  return(list(`Cell_Sets` = cell_gs,
              Phenotype_sets = av_gs))
}, mc.cores = 4)
#save gene set output
# saveRDS(all_fov, "data/all_fov_gsea.rds")
all_fov = readRDS("data/all_fov_gsea.rds")


# Plot --------------------------------------------------------------------
#gsea for first FOV
fov_gsea = all_fov$slide1_fov100$Cell_Sets %>%
  data.frame(check.names = FALSE) %>%
  rownames_to_column("Hallmark Gene Set") %>%
  gather("cell_ID", "score", -`Hallmark Gene Set`)
fov_meta = meta %>%
  filter(cell_ID %in% fov_gsea$cell_ID)

fov_plotting_dat = fov_meta %>% 
  select(unique_fov, #sample ID
         x_FOV_px, y_FOV_px, #cell location
         unsup_insitu, cell_ID) %>%
  full_join(fov_gsea)

#pick a couple gene sets
unique(fov_plotting_dat$`Hallmark Gene Set`)
fov_plotting_dat %>%
  filter(`Hallmark Gene Set` %in% c("HALLMARK_APOPTOSIS")) %>%
  ggplot() +
  geom_point(aes(x = x_FOV_px, y = y_FOV_px, color = score), size = 3) +
  scale_colour_gradient(low = "blue", high = "red") +
  facet_wrap(~`Hallmark Gene Set`) +
  coord_equal() +
  theme_bw()
