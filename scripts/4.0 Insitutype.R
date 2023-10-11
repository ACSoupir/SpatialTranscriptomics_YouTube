
# Libraries ---------------------------------------------------------------

library(InSituType)
library(Seurat)
library(tidyverse)


# Get data ----------------------------------------------------------------

object = readRDS("data/HFC_sctransform.rds")

#working data
counts = object@assays$RNA@counts %>%
  as.matrix() %>% t()

neg = object@assays$negprobes@counts %>% 
  as.matrix() %>% t()

#calculate average negative probe per cell
negmean = rowMeans(neg)
head(negmean)

#prep immunofluorecence data
immunofluorescence = object@meta.data %>%
  select(contains(c("Mean", "Max"))) %>%
  select(-contains("DAPI"))

#get per cell background
negmean.per.totcount <- mean(rowMeans(neg)) / mean(rowSums(counts))
per.cell.bg <- rowSums(counts) * negmean.per.totcount
head(per.cell.bg)


# Calculate cohorts -------------------------------------------------------
#will leave n_cohorts for fully unsupervised
cohort <- fastCohorting(immunofluorescence,
                        gaussian_transform = TRUE)
#check clusters and cohort numbers
table(cohort)


# Unsupervised Clustering -------------------------------------------------

unsup <- insitutype(
  x = counts,
  neg = negmean,
  cohort = cohort,
  #can use the per cell background from above
  bg = per.cell.bg,
  #we see 14 from louvain clustering, set from 10 to 20
  n_clusts = 10:20,
  #NULL value runs unsupervised clustering; entering a matrix here would run semi-supervised clustering.
  reference_profiles = NULL,
  n_phase1 = 200,
  n_phase2 = 500,
  n_phase3 = 2000,
  n_starts = 1,
  max_iters = 5
)

str(unsup)

#look at expression profiles
heatmap(sweep(unsup$profiles, 1, pmax(apply(unsup$profiles, 1, max), .2), "/"), scale = "none",
        main = "Cluster mean expression profiles")


# plot assignments --------------------------------------------------------

cols <-
  c(
    '#8DD3C7',
    '#BEBADA',
    '#FB8072',
    '#80B1D3',
    '#FDB462',
    '#B3DE69',
    '#FCCDE5',
    '#D9D9D9',
    '#BC80BD',
    '#CCEBC5',
    '#FFED6F',
    '#E41A1C',
    '#377EB8',
    '#4DAF4A',
    '#984EA3',
    '#FF7F00',
    '#FFFF33',
    '#A65628',
    '#F781BF',
    '#999999'
  )

cols <- cols[seq_along(unique(unsup$clust))]
names(cols) <- unique(unsup$clust)

par(mfrow = c(1, 1))
par(mar = c(0, 0, 3, 0))

# plot(mini_nsclc$x, mini_nsclc$y, pch = 16, cex = .75, asp = 1, cex.main = 0.75,
#      main = "cells in physical space",
#      col = cols[unsup$clust], xlab = "", ylab = "", xaxt = "n", yaxt = "n")

plot(object@reductions$UMAP@cell.embeddings, pch = 16, cex = .75, asp = 1, cex.main = 0.75,
     main = "cells in UMAP space",     
     col = cols[unsup$clust], xlab = "", ylab = "", xaxt = "n", yaxt = "n")
legend("bottomleft", pch = 16, col = cols, legend = names(cols), cex = 0.7)



# flight plot! ------------------------------------------------------------

# make the flightpath plot
fp <- flightpath_plot(flightpath_result = NULL, insitutype_result = unsup, col = cols[unsup$clust])
class(fp)
#> [1] "gg"     "ggplot"
print(fp)



# Plot FOV from before ----------------------------------------------------

object$unsup_insitu = unsup$clust
metadata = object@meta.data

slide1_fov1 = metadata %>%
  filter(fov == 1,
         slide_ID_numeric == 1)

alex_theme = theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

#plot slide 1 FOV 1
slide1_fov1 %>%
  ggplot() +
  geom_point(aes(x = x_FOV_px, y = y_FOV_px, color = unsup_insitu)) +
  coord_fixed() +
  labs(title = "New InSituType Clusters") + 
  alex_theme

slide1_fov1 %>%
  ggplot() +
  geom_point(aes(x = x_FOV_px, y = y_FOV_px, color = RNA_nbclust_clusters_refined)) +
  coord_fixed() +
  labs(title = "Nanostring Provided") +
  alex_theme


# save object -------------------------------------------------------------

saveRDS(object, "data/HFC_insitutype")
