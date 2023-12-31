
# Libraries ---------------------------------------------------------------

library(InSituType)
library(Seurat)
library(tidyverse)

# Data --------------------------------------------------------------------

object = readRDS("data/HFC_insitutype.rds")

#extract metadata and then save, remove object
metadata = object@meta.data
#saveRDS(metadata, "data/insitutyped_metadata.rds")
rm(object)


# Prep Data ---------------------------------------------------------------

#create a unique id for the slide - fov combinations
metadata$unique_fov = paste0("slide", metadata$slide_ID_numeric, "_fov", metadata$fov)

#extract and make long data frame of relevant spatial information
spatial_data = metadata %>%
  select(unique_fov, #sample ID
         x_FOV_px, y_FOV_px, #cell location
         unsup_insitu) %>% #cell cluster
  mutate(positive = 1) %>% #add column for positivity of the cells
  arrange(unsup_insitu) %>% #arrange cells clusters in alphabetical order
  pivot_wider(names_from = "unsup_insitu",
              values_from = "positive",
              values_fill = 0)

#convert to spatial list
spatial_list = split(spatial_data, spatial_data$unique_fov)

#create sample summary level data
cell_types = letters[1:14]
summary_df = lapply(spatial_list, function(spat){
  #need to collapse to the nubmer of positive of each cell type and 
  #percent of cells postive for each cell type
  spat %>% 
    group_by(unique_fov) %>% #can remove cell locations here for sample level information
    summarise(total_cells = n(),
              across(!!cell_types, sum)) %>% #get cell numbers
    mutate(across(!!cell_types, ~ .x /total_cells * 100, .names = "percent_{col}")) #percent of total
}) %>%
  do.call(bind_rows, .) #bind all samples to one data frame

#add column to summary_df for 'patient' level
#can be the same as the unique_fov column
#is linker with spatial time between summary_df and clinical
summary_df$patient_id = summary_df$unique_fov

#create the clinical data
#this is also sample level data

clinical = metadata %>% 
  select(unique_fov, fov, slide_ID_numeric, assay_type, Run_name, tissue,
         Panel, version, contains("qc"), median_negprobes:percOfDataFromErrorPerCell) %>% 
  distinct() %>%
  mutate(patient_id = unique_fov, .before = 1)


# SpatialTIME - Build Object ----------------------------------------------

if(!require(spatialTIME)){
  print("Missing spatialTIME")
  install.packages("spatialTIME")
  library(spatialTIME)
}

#create mif
spat_obj = create_mif(clinical_data = clinical,
                      sample_data = summary_df,
                      spatial_list = spatial_list,
                      patient_id = "patient_id",
                      sample_id = "unique_fov")
#view contents
spat_obj

# Run Nearest Neighbor G(r) Function --------------------------------------
#find range of x and y for appropriate function range
ranges = spatial_data %>%
  group_by(unique_fov) %>%
  summarise(xrange = max(x_FOV_px) - min(x_FOV_px),
            yrange = max(y_FOV_px) - min(y_FOV_px))
view(ranges)

#since minimum size is ~1400, a range of 1-200 will be more than fine
#NOTE: for permutations, can take a VERY long time to run them all
# the more permutations, the more accurate the CSR estimate
# if calculating Ripley's K(r), permutations aren't needed and can use exact measure
#   - see Ripley's K(r) documentation for more information
#NOTE: for number of workers, the environment is duplicated (parallel package)
# need to keep in mind the size of the environment and available ram
# can remove extra objects like metadata and spatial_data to cut down on size
spat_obj2 = NN_G(mif = spat_obj, #object with all 
                 mnames = cell_types, #the columns of the cell types we want to measure
                 r_range = 0:200, #radius to calculate summary function
                 num_permutations = 100, #number of permutations to estimate CSR
                 edge_correction = "rs", #edge correction - reduced sampling
                 keep_perm_dis = FALSE, #whether to keep iterations fo CSR estimates
                 workers = 16, #number of CPU cores to use for calculating G(r)
                 overwrite = TRUE, #if G(r) exists in object, replace
                 xloc = "x_FOV_px", yloc = "y_FOV_px") #cell location columns
#took abotu 17 minutes to run with ~12GB of RAM

#saveRDS(spat_obj2, "data/nng_mif.rds")

# Visualize Spatial Summary NN G(r) ---------------------------------------

#extract the derived table from object
nng_derived = spat_obj2$derived$univariate_NN
head(nng_derived)

#plot for first FOV
sum_plot = nng_derived %>%
  tibble() %>%
  filter(unique_fov == "slide1_fov1") %>%
  ggplot() +
  geom_line(aes(x = r, y = `Degree of Clustering Permutation`,
                color = Marker)) +
  labs(title = "Slide 1 - FOV 1 Nearest Neighbor G(r)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
point_plot = spat_obj2$spatial$slide1_fov1 %>%
  gather("Marker", "Positive", -unique_fov, -x_FOV_px, -y_FOV_px) %>%
  filter(Positive == 1) %>%
  ggplot() +
  geom_point(aes(x = x_FOV_px, y = y_FOV_px,
                color = Marker), size = 3) +
  labs(title = "Slide 1 - FOV 1 Point Plot") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_equal() #show correct scale of sample

#show together
ggpubr::ggarrange(sum_plot, point_plot,
                  nrow = 2, ncol = 1)








