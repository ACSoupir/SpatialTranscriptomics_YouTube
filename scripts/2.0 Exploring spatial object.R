
# Libraries ---------------------------------------------------------------

library(Seurat)
library(tidyverse)


# Load Data ---------------------------------------------------------------

object = readRDS("data/SeuratObj_withTranscripts.RDS")


# Metadata ----------------------------------------------------------------

metadata = object@meta.data


# Explore data ------------------------------------------------------------

slide_fov = metadata %>%
  select(fov, slide_ID_numeric) %>%
  distinct()

slide1_fov1 = metadata %>%
  filter(fov == 1,
         slide_ID_numeric == 1)

#each pixel = 0.180 um

#plot slide 1 FOV 1
slide1_fov1 %>%
  ggplot() +
  geom_point(aes(x = x_FOV_px, y = y_FOV_px, color = RNA_nbclust_clusters_refined)) +
  coord_fixed()

slide1_fov1 %>%
  ggplot() +
  geom_point(aes(x = x_FOV_px, y = y_FOV_px, color = RNA_nbclust_clusters_refined)) +
  coord_fixed() +
  labs(title = "Slide 1 - FOV 1: Human Frontal Cortex (Nanostring CosMx SMI)") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))


# Exploring Cell 2 of FOV 1 on Slide 1 ------------------------------------

slide1_fov1_cell2 = object@misc$transcriptCoords %>%
  filter(slideID == 1,
         fov == 1,
         CellId == 2) 

#plottign cell 2 transcript locations
slide1_fov1_cell2 %>%
  ggplot() +
  geom_point(aes(x = x_FOV_px, y = y_FOV_px, color = target)) +
  coord_fixed() +
  theme(legend.position = "none")
