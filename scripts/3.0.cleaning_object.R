
# Libraries ---------------------------------------------------------------

library(Seurat)
library(tidyverse)

# Load Data ---------------------------------------------------------------

object = readRDS("data/SeuratObj_withTranscripts.RDS")

#check size of the transcript coordinates
print(object.size(object@misc$transcriptCoords), units = "MB")

#save compressed
transcripts = object@misc$transcriptCoords
data.table::fwrite(transcripts, "data/trancript_locations.csv.gz")

#remove transcript locations from object
object@misc = list(NULL)
rm(transcripts)

#clean environment
gc(full = TRUE)

#save cleaned object
saveRDS(object, "data/HFC_no_transcript.rds")

#read in our no transcript object
object = readRDS("data/HFC_no_transcripts.rds")


# Assays in Object --------------------------------------------------------

Assays(object)
object@assays$RNA_normalized = NULL
gc(full = TRUE)

#what reductions are in object
Reductions(object)
object@reductions = list(NULL)

#save our fresh object
saveRDS(object, "data/HFC_fresh.rds")
