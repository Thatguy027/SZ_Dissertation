try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)))

# data folders
arsenic_data <- "Arsenic/Data/"
etoposide_data <- "Etoposide/Data/"
ben1_data <- "201809_rappaport/Data/ben1_data/Final_Tables/"
outlier_data <- "Data/Strain_Outliers/Data/"

# source base theme and colors for plots
source("Scripts/Figure_Themes.R")

# source scripts 
source("Scripts/PD_Talk_Functions.R")

library("ape")

ancestry.colours <- c("gold2", "plum4","darkorange1", 
                      "lightskyblue2", "firebrick","burlywood3", "gray51", 
                      "springgreen4", "lightpink2", "deepskyblue4", "black", 
                      "mediumpurple4","orange", "maroon", "yellow3", "brown4", 
                      "yellow4", "sienna4","chocolate", "gray19")
########################################################################################################################
# outlier regions
########################################################################################################################

outlier_regions <- data.table::fread(glue::glue("{outlier_data}/Population_Outlier_Bins.tsv"), 
                                   col.names = c("CHROM",	"START_BIN",	"END_BIN",	"COUNT",	"MID_BIN",	"index",	"rank",	"outlier",	"direction",	"GENOMIC_REGION",	"STRAIN")) %>%
  tidyr::unite(MARKER_BIN, CHROM, MID_BIN, sep = ":", remove = F) %>%
  # tidyr::unite(STRAIN_RANK, STRAIN, rank, sep = "_", remove = F) %>%
  dplyr::select(CHROM, MID_BIN, MARKER_BIN, rank, STRAIN) %>%
  dplyr::arrange(STRAIN, CHROM, MID_BIN, rank) %>%
  dplyr::distinct(CHROM, MID_BIN, MARKER_BIN, STRAIN, .keep_all = T) %>%
  # dplyr::filter(rank < 5) %>%
  tidyr::spread(STRAIN, rank)

outlier_regions[is.na(outlier_regions)] <- 0

outlier_regions_toclust <- outlier_regions %>%
  dplyr::select(-CHROM, -MARKER_BIN, -MID_BIN) %>%
  t()

dd <- dist(scale(outlier_regions_toclust), method = "euclidean")
hc <- hclust(dd, method = "ward.D2")

clus6 = cutree(hc, 20)
# plot(as.phylo(hc), hang = -1, cex = 0.6, type = "unrooted", no.margin = TRUE)
plot(as.phylo(hc), type = "fan", tip.color = ancestry.colours[clus6],
     label.offset = 1, cex = 0.7)


outlier_pca <- prcomp(scale(outlier_regions_toclust))
outlier_pca_df <- data.frame(strain = colnames(outlier_regions[,4:ncol(outlier_regions)]),
                             outlier_pca$x)

ggplot(outlier_pca_df)+
  aes(x = PC4, y = PC7)+
  geom_point()

View(outlier_pca_df[,1:4])
