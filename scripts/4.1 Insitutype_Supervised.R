
# Libraries ---------------------------------------------------------------

library(InSituType)
library(Seurat)
library(tidyverse)

# load data ---------------------------------------------------------------

object = readRDS("data/HFC_sctransform.rds")

#load in profiles
#https://github.com/Nanostring-Biostats/CellProfileLibrary/tree/master
load(url("https://github.com/Nanostring-Biostats/CellProfileLibrary/raw/master/Human/Adult/Brain_AllenBrainAtlas.RData"))

#prep sample data
#must be count data!
counts = object@assays$RNA@counts %>% 
  as.matrix() %>% t()
counts[1:10, 1:10]

negmean = object@assays$negprobes@counts %>%
  as.matrix() %>% t() %>%
  Matrix::rowMeans()
head(negmean)


# Get cohorts -------------------------------------------------------------
#can create cell cohorts many different ways but will just use immunofluorescence
cohort = fastCohorting(object@meta.data[,colnames(object@meta.data) %>% grep("Max|Mean", ., value= TRUE)],
                       gaussian_transform = TRUE)
table(cohort)


# Perform Supervised Clustering -------------------------------------------
#using counts, false positives, our cell groups, and profiles, annotate cells
supervised_clusters = insitutypeML(x = counts,
                                   neg = negmean,
                                   cohort = cohort,
                                   reference_profiles = as.matrix(profile_matrix))
#some genes aren't in reference
#example is ITGA5 in expression but called ITGAV in reference
#ITGA10 in expression and ITGAX in reference
#I think this comes down to this particular _run_ and not reference
#if your data has many missing, might require a gene symbol review

#add phenotypes to the metadata for plotting
object$InSituTypeIDs = supervised_clusters$clust

# Visualize ---------------------------------------------------------------

#heatmap
heatmap(sweep(supervised_clusters$profiles, 1, 
              pmax(apply(supervised_clusters$profiles, 1, max), 0.2), "/"),
        scale = "none", 
        main = "Mean cell expression by phenotype",
        margins = c(10, 5))

#UMAP
DimPlot(object, group.by = "InSituTypeIDs", reduction = "UMAP")


#flight plot to view posterio$ probabilities
cols <- c('#8DD3C7','#BEBADA','#FB8072','#80B1D3','#FDB462',
          '#B3DE69','#FCCDE5','#D9D9D9','#BC80BD','#CCEBC5',
          '#FFED6F','#E41A1C','#377EB8','#4DAF4A','#984EA3',
          '#FF7F00','#FFFF33','#A65628','#F781BF','#999999')
#get same number of colors as cell types
cols <- cols[seq_along(unique(supervised_clusters$clust))]
#give colors the cell type names
names(cols) <- unique(supervised_clusters$clust)
#plot flightpath plot
fp <- flightpath_plot(flightpath_result = NULL,
                      insitutype_result = supervised_clusters,
                      col = cols[supervised_clusters$clust])
print(fp)
#we want higher posterior probabilities which indicates greater confidence in asignment

#plot 2d of first FOV
object@meta.data %>%
  filter(slide_ID_numeric == 1,
         fov == 1) %>% 
  ggplot() +
  geom_point(aes(x = x_FOV_px, y = y_FOV_px, color = InSituTypeIDs)) +
  theme_bw() +
  coord_fixed()
