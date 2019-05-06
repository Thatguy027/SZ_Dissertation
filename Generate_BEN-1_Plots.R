try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)))

# data folders
arsenic_data <- "Arsenic/Data/"
etoposide_data <- "Etoposide/Data/"
ben1_data <- "201809_rappaport/Data/ben1_data/Final_Tables/"

# source base theme and colors for plots
source("Scripts/Figure_Themes.R")

# source scripts 
source("Scripts/PD_Talk_Functions.R")


########################################################################################################################
# BEN-1
########################################################################################################################

# # # LOAD AND PROCESS INDELS
ben1_variants <- data.table::fread(glue::glue("{ben1_data}TS9_ben-1_variants.tsv"))%>%
  na.omit()

gwa_mappings <-  data.table::fread(file = glue::glue("{ben1_data}TS5_GWA_processed_marker_mapping.tsv"))%>%
  na.omit()%>%
  dplyr::mutate(snpGT = ifelse(allele==-1,"REF", "ALT"))%>%
  dplyr::select(snpMarker = marker, strain, trait, value, snpGT)%>%
  dplyr::filter(trait == "q90.tof")%>%
  dplyr::left_join(.,ben1_variants, by = "strain")