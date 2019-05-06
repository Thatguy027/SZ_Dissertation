try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)))

# data folders
etoposide_data <- "Etoposide/Data/"
etoposide_cegwas2_data <- "Etoposide/cegwas2/data/"

# source base theme and colors for plots
source("Scripts/Figure_Themes.R")

# source scripts 
source("Scripts/PD_Talk_Functions.R")


########################################################################################################################
# ETOPOSIDE
########################################################################################################################

######################################################################################################################## Linkage LOD PLOT

etoposide_linkage <- data.table::fread(glue::glue("{etoposide_data}TableS3_LM-Etoposide-medianTOF.csv"))

etoposide_linkage %>%
  maxlodplot_edit() +
  base_theme +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.border = element_rect(fill = NA))

ggsave(filename = "Plots/etoposide_medianTOF_Linkage.png", height = 4, width = 12, dpi = 400)
ggsave(filename = "Plots/etoposide_medianTOF_Linkage.pdf", height = 4, width = 12, dpi = 400)

######################################################################################################################## Linkage LOD PLOT

etoposide_linkage_pheno <- data.table::fread(glue::glue("{etoposide_data}TableS2_LM-phenotypes.csv"))

etoposide_linkage_pheno <- etoposide_linkage_pheno %>%
  dplyr::filter(trait == "median.TOF") %>%
  dplyr::select(-condition)

data("N2xCB4856cross")

blankcross <- N2xCB4856cross

# Plot Linkage Mapping - PxG
etoposide_cross <- linkagemapping::mergepheno(blankcross, etoposide_linkage_pheno, set = 2)

corrected_map <- dplyr::filter(etoposide_linkage, trait == "etoposide.median.TOF") 
corrected_map$trait <- ".median.TOF"

pxgplot_edit(cross = etoposide_cross, map = corrected_map)+ 
  scale_x_discrete(breaks=c("N2", "CB4856"),labels=c("N2", "CB4856"))+
  labs(y = "Etoposide Response", x = "RIAIL Genotype at Peak Marker") +
  base_theme +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = axis_color)) +
  labs(title = "")
 
ggsave(filename = "Plots/etoposide_medianTOF_Linkage_PxG.png", height = 6, width = 8, dpi = 400)
ggsave(filename = "Plots/etoposide_medianTOF_Linkage_PxG.pdf", height = 6, width = 8, dpi = 400)

######################################################################################################################## NILs

etoposide_nil_pheno <- data.table::fread(glue::glue("{etoposide_data}TableS8_NIL-phenotypes.csv"))

etoposide_nil_pheno %>%
  dplyr::filter(trait == "median.TOF",
                strain %in% c("N2","CB4856","ECA220","ECA219","ECA216"))%>%
  dplyr::mutate(strain1 = factor(strain, levels = c("N2","CB4856","ECA216","ECA219","ECA220"),
                                 labels =c("N2",
                                           "CB4856", 
                                           "N2>CB4856\n II:11.43 - 12.11 Mb ", 
                                           "CB4856>N2\n II:11.64 - 11.9 Mb",
                                           "N2>CB4856\n II:12.01 - 12.1 Mb ")))%>%
  ggplot(.) +
  aes(x = factor(strain1), 
      y = phenotype, 
      fill=strain) +
  geom_beeswarm(alpha = point_alpha,
                priority = "density",
                cex = 1.2)+
  geom_boxplot(outlier.colour = NA, alpha = boxplot_alpha)+
  scale_fill_manual(values = c("N2" = "#F9A227","CB4856" = "#2790F9",
                               "ECA220" = "#F9A227","ECA219" = "#2790F9","ECA216" = "#ECA216"))+
  base_theme +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = axis_color),
        axis.title.x = element_blank(),
        plot.title = element_blank()) +
  labs( y = paste0("Animal Length"))

ggsave(filename = "Plots/etoposide_medianTOF_NIL.png", height = 6, width = 12, dpi = 400)
ggsave(filename = "Plots/etoposide_medianTOF_NIL.pdf", height = 6, width = 12, dpi = 400)

######################################################################################################################## GWA MANHATTAN PLOT

etoposide_gwa <- data.table::fread(glue::glue("{etoposide_cegwas2_data}etoposide_median.TOF_processed_mapping.tsv"))
independent_tests <- 827
gwas_fine_mapping <- readr::read_tsv(glue::glue("{etoposide_cegwas2_data}etoposide_median.TOF_processed_mapping.tsv")) 
# geno_matrix <- readr::read_tsv(glue::glue("{etoposide_data}Figure 2-source data 8.zip"))%>%
#   na.omit()

cegwas2_manplot(etoposide_gwa, eigen_cutoff = -log10(0.05/independent_tests))[[1]] + 
  base_theme +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.border = element_rect(fill = NA), 
        plot.title = element_blank())

ggsave(filename = "Plots/etoposide_medianTOF_GWA.png", height = 4, width = 12, dpi = 400)
ggsave(filename = "Plots/etoposide_medianTOF_GWA.pdf", height = 4, width = 12, dpi = 400)

######################################################################################################################## GWA PxG

peak_pos <- na.omit(etoposide_gwa) %>%
  dplyr::filter(CHROM == "II", POS > 10e6) %>%
  dplyr::mutate(facet_marker = paste0(CHROM, ":", peakPOS)) %>%
  dplyr::pull(facet_marker) %>%
  unique()

pxg_df <- na.omit(etoposide_gwa) %>%
  dplyr::filter(CHROM == "II", POS > 10e6) %>%
  dplyr::mutate(facet_marker = paste0(CHROM, ":", peakPOS)) %>%
  dplyr::group_by(allele, facet_marker)%>%
  dplyr::mutate(mean_pheno = mean(as.numeric(value), na.rm = T))%>%
  dplyr::mutate(n2_cb = case_when(
    strain == "N2" ~ "N2",
    strain == "CB4856" ~ "CB4856", 
    TRUE ~ "Other"
  )) 

pxg_df %>%
  ggplot()+
  aes(x = factor(allele, levels = c(-1,1), labels = c("REF","ALT")))+
  geom_beeswarm(cex=1.2,
                priority='density',
                aes(y = as.numeric(value),
                    fill = n2_cb,
                    size = n2_cb),
                shape = 21, 
                alpha = point_alpha,
                data = dplyr::filter(pxg_df, !strain %in% c("CB4856", "N2")))+
  geom_boxplot(aes(y=as.numeric(value)), alpha = boxplot_alpha, fill = "gray70", outlier.colour = NA) +
  geom_beeswarm(cex=1.2,
                priority='density',
                aes(y = as.numeric(value),
                    fill = n2_cb,
                    size = n2_cb),
                shape = 21, 
                data = dplyr::filter(pxg_df, strain %in% c("CB4856", "N2")))+
  scale_fill_manual(values=strain_colors)+
  scale_size_manual(values=c(point_highlight_size,point_highlight_size,point_size))+
  labs(y = "Animal Length",
       x = glue::glue("Genotype at {peak_pos}")) +
  base_theme +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = axis_color))

ggsave(filename = "Plots/etoposide_medianTOF_GWA_PxG.png", height = 6, width = 8, dpi = 400)
ggsave(filename = "Plots/etoposide_medianTOF_GWA_PxG.pdf", height = 6, width = 8, dpi = 400)

######################################################################################################################## GWA PEAK LD
gm <- readr::read_tsv(glue::glue("{etoposide_cegwas2_data}Genotype_Matrix.tsv"))
etoposide_gwa <- data.table::fread(glue::glue("{etoposide_cegwas2_data}etoposide_median.TOF_processed_mapping.tsv"))

LD_output <- Plot_Peak_LD(etoposide_gwa, gm)

LD_output[[1]] + 
  base_theme + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())

ggsave(filename = "Plots/etoposide_medianTOF_GWA_PeakLD.png", height = 10, width = 14, dpi = 400)
ggsave(filename = "Plots/etoposide_medianTOF_GWA_PeakLD.pdf", height = 10, width = 14, dpi = 400)

######################################################################################################################## GWA Fine mapping
# something is wrong with the compression, need to compress to push to github an uncompress to load and 
etoposide_fine_mapping <-  data.table::fread(glue::glue("unzip -p {etoposide_cegwas2_data}etoposide_median.TOF_snpeff_genes.tsv.zip")) 

snpeff_fine <- etoposide_fine_mapping %>%
  dplyr::filter(CHROM == "II", POS > 10e6) %>%
  dplyr::select(MARKER, POS, STRAIN, REF,ALT, TGT = STRAIN_GENOTYPE, VARIANT_IMPACT,
                VARIANT_LD_WITH_PEAK_MARKER, PEAK_MARKER, QTL_INTERVAL_START,
                QTL_INTERVAL_END, VARIANT_LOG10p)

snpeff_fine$VARIANT_IMPACT[is.na(snpeff_fine$VARIANT_IMPACT)] <- "INTERGENIC"

LD_genotypes <- snpeff_fine %>%
  dplyr::filter(STRAIN == "CB4856") %>%
  dplyr::mutate(cb_alt = ifelse(REF == TGT, "CB4856 REF", "CB4856 ALT")) %>%
  dplyr::mutate(tidy_marker = gsub("_",":",MARKER))

peak_roi_marker <- LD_genotypes %>%
  dplyr::filter(tidy_marker == PEAK_MARKER)

LD_genotypes%>%
  na.omit() %>%
  ggplot() +
  aes(x = POS/1e6) +
  geom_vline(aes(xintercept = 11.872440), 
             color = "red",
             linetype = 2) +
  geom_vline(aes(xintercept = 11.64), color = "gray60")+
  geom_vline(aes(xintercept = 12.11), color = "gray60")+
  geom_point(aes(fill = factor(VARIANT_IMPACT,
                               levels = rev(c("INTERGENIC", "MODIFIER", "LOW", "MODERATE", "HIGH")), ordered = T), 
                 y = VARIANT_LOG10p), 
             size = point_size,
             shape = 21)+
  facet_grid(.~cb_alt) + 
  scale_fill_viridis_d(name = "Variant\nImpact", direction = -1) +
  theme_bw(15)+
  base_theme + 
  labs(x = "Genomic Position (Mb)",
       y = expression(-log[10](italic(p))))+
  theme(panel.grid.major = element_blank(),
        axis.line = element_line(colour = axis_color))

ggsave(filename = "Plots/Etoposide_medianTOF_GWA_FineMap.png", height = 6, width = 16, dpi = 400)
ggsave(filename = "Plots/Etoposide_medianTOF_GWA_FineMap.pdf", height = 6, width = 16, dpi = 400)

######################################################################################################################## SWAP - Etoposide only

etoposide_swap <- data.table::fread(glue::glue("{etoposide_data}TableS12_Allele-swap-all-drugs.csv")) 

etoposide_swap%>%
  dplyr::filter(trait == "mean.TOF",
                condition == "Etoposide",
                strain!="Bristol(M)ECA400",
                strain!="Bristol(M)ECA401",
                strain!="Hawaii(Q)ECA547",
                strain!="Hawaii(Q)ECA549")%>%
  dplyr::mutate(strain1 = factor(strain, levels = c("Bristol", "Bristol(M)ECA402" ,"Hawaii", "Hawaii(Q)ECA548" ),
                                 labels =c("N2\nTOP-2(Q762)", "N2\nTOP-2(M762)", "CB4856\nTOP-2(M762)","CB4856\nTOP-2(Q762)")))%>%
  ggplot(.) +
  aes(x = factor(strain1), 
      y = phenotype, 
      fill= strain) +
  geom_beeswarm(cex=1.2,priority='density', alpha = point_alpha, size = point_size)+
  geom_boxplot(outlier.colour = NA, alpha = boxplot_alpha)+
  scale_fill_manual(values = c("Bristol" = "#F9A227", "Hawaii" = "#2790F9",
                               "Bristol(M)ECA402" = "gray50","Hawaii(Q)ECA548" = "gray50"))+
  base_theme +  theme(legend.position = "none",
                      panel.grid.major = element_blank(),
                      axis.line = element_line(colour = axis_color),
                      axis.title.x = element_blank()) +
  labs( y = paste0("Animal Length"))

ggsave(filename = "Plots/etoposide_medianTOF_SWAP.png", height = 6, width = 10, dpi = 400)
ggsave(filename = "Plots/etoposide_medianTOF_SWAP.pdf", height = 6, width = 10, dpi = 400)

######################################################################################################################## SWAP - All poisons

poison_swap <- data.table::fread(glue::glue("{etoposide_data}CRISPR_swap_regress_phenotype.csv")) 

######################################################################################################################## amsacrine
poison_swap%>%
  dplyr::filter(trait == "median.TOF",
                condition == "Amsacrine",
                strain!="ECA400", 
                strain != "ECA549", 
                strain != "ECA401", 
                strain != "ECA548")%>%
  dplyr::mutate(strain1 = factor(strain, levels = c("N2", "ECA400", "ECA401", "ECA402" ,"CB4856","ECA547", "ECA548" ,"ECA549"),
                                 labels =c("N2\nTOP-2(Q762)", "N2\nTOP-2(M762)", "N2\nTOP-2(M762)","N2\nTOP-2(M762)",
                                           "CB4856\nTOP-2(M762)","CB4856\nTOP-2(Q762)","CB4856\nTOP-2(Q762)","CB4856\nTOP-2(Q762)"))) %>%
  ggplot(.) +
  aes(x = factor(strain1), 
      y = phenotype, 
      fill= strain) +
  geom_beeswarm(cex=1.2,priority='density', alpha = point_alpha, size = point_size)+
  geom_boxplot(outlier.colour = NA, alpha = boxplot_alpha)+
  scale_fill_manual(values = c("N2" = "#F9A227", "CB4856" = "#2790F9",
                               "ECA400" = "gray50","ECA401" = "gray50","ECA402" = "gray50",
                               "ECA547" = "gray50","ECA548"= "gray50" ,
                               "ECA549"= "gray50"))+
  facet_grid(.~condition, scales = "free") +
  base_theme +  theme(legend.position = "none",
                      panel.grid.major = element_blank(),
                      axis.line = element_line(colour = axis_color),
                      axis.title.x = element_blank()) +
  labs( y = paste0("Animal Length"))

ggsave(filename = "Plots/amsacrine_medianTOF_SWAP.png", height = 6, width = 10, dpi = 400)
ggsave(filename = "Plots/amsacrine_medianTOF_SWAP.pdf", height = 6, width = 10, dpi = 400)

######################################################################################################################## xk469
poison_swap%>%
  dplyr::filter(trait == "median.TOF",
                condition == "XK469",
                strain!="ECA400", 
                strain != "ECA549", 
                strain != "ECA401", 
                strain != "ECA548")%>%
  dplyr::mutate(strain1 = factor(strain, levels = c("N2", "ECA400", "ECA401", "ECA402" ,"CB4856","ECA547", "ECA548" ,"ECA549"),
                                 labels =c("N2\nTOP-2(Q762)", "N2\nTOP-2(M762)", "N2\nTOP-2(M762)","N2\nTOP-2(M762)",
                                           "CB4856\nTOP-2(M762)","CB4856\nTOP-2(Q762)","CB4856\nTOP-2(Q762)","CB4856\nTOP-2(Q762)"))) %>%
  ggplot(.) +
  aes(x = factor(strain1), 
      y = phenotype, 
      fill= strain) +
  geom_beeswarm(cex=1.2,priority='density', alpha = point_alpha, size = point_size)+
  geom_boxplot(outlier.colour = NA, alpha = boxplot_alpha)+
  scale_fill_manual(values = c("N2" = "#F9A227", "CB4856" = "#2790F9",
                               "ECA400" = "gray50","ECA401" = "gray50","ECA402" = "gray50",
                               "ECA547" = "gray50","ECA548"= "gray50" ,
                               "ECA549"= "gray50"))+
  facet_grid(.~condition, scales = "free") +
  base_theme +  theme(legend.position = "none",
                      panel.grid.major = element_blank(),
                      axis.line = element_line(colour = axis_color),
                      axis.title.x = element_blank()) +
  labs( y = paste0("Animal Length"))

ggsave(filename = "Plots/XK469_medianTOF_SWAP.png", height = 6, width = 10, dpi = 400)
ggsave(filename = "Plots/XK469_medianTOF_SWAP.pdf", height = 6, width = 10, dpi = 400)

######################################################################################################################## teniposide
poison_swap%>%
  dplyr::filter(trait == "q90.TOF",
                condition == "Teniposide",
                strain!="ECA400", 
                strain != "ECA549", 
                strain != "ECA401", 
                strain != "ECA548")%>%
  dplyr::mutate(strain1 = factor(strain, levels = c("N2", "ECA400", "ECA401", "ECA402" ,"CB4856","ECA547", "ECA548" ,"ECA549"),
                                 labels =c("N2\nTOP-2(Q762)", "N2\nTOP-2(M762)", "N2\nTOP-2(M762)","N2\nTOP-2(M762)",
                                           "CB4856\nTOP-2(M762)","CB4856\nTOP-2(Q762)","CB4856\nTOP-2(Q762)","CB4856\nTOP-2(Q762)"))) %>%
  ggplot(.) +
  aes(x = factor(strain1), 
      y = phenotype, 
      fill= strain) +
  geom_beeswarm(cex=1.2,priority='density', alpha = point_alpha, size = point_size)+
  geom_boxplot(outlier.colour = NA, alpha = boxplot_alpha)+
  scale_fill_manual(values = c("N2" = "#F9A227", "CB4856" = "#2790F9",
                               "ECA400" = "gray50","ECA401" = "gray50","ECA402" = "gray50",
                               "ECA547" = "gray50","ECA548"= "gray50" ,
                               "ECA549"= "gray50"))+
  facet_grid(.~condition, scales = "free") +
  base_theme +  theme(legend.position = "none",
                      panel.grid.major = element_blank(),
                      axis.line = element_line(colour = axis_color),
                      axis.title.x = element_blank()) +
  labs( y = paste0("Animal Length"))

ggsave(filename = "Plots/Teniposide_medianTOF_SWAP.png", height = 6, width = 10, dpi = 400)
ggsave(filename = "Plots/Teniposide_medianTOF_SWAP.pdf", height = 6, width = 10, dpi = 400)

######################################################################################################################## dactinomycin

dact_swap <- data.table::fread(glue::glue("{etoposide_data}top2_crisprswap_dactinomycin_assayregressed.csv")) 

dact_swap%>%
  dplyr::filter(trait == "norm.n",
                strain!= "ECA400", 
                strain != "ECA401", 
                strain != "ECA548")%>%
  dplyr::mutate(strain1 = factor(strain, levels = c("N2", "ECA400", "ECA401", "ECA402" ,"CB4856","ECA547", "ECA548" ,"ECA549"),
                                 labels =c("N2\nTOP-2(Q762)", "N2\nTOP-2(M762)", "N2\nTOP-2(M762)","N2\nTOP-2(M762)",
                                           "CB4856\nTOP-2(M762)","CB4856\nTOP-2(Q762)","CB4856\nTOP-2(Q762)","CB4856\nTOP-2(Q762)"))) %>%
  ggplot(.) +
  aes(x = factor(strain1), 
      y = phenotype, 
      fill= strain) +
  geom_beeswarm(cex=1.2,priority='density', alpha = point_alpha, size = point_size)+
  geom_boxplot(outlier.colour = NA, alpha = boxplot_alpha)+
  scale_fill_manual(values = c("N2" = "#F9A227", "CB4856" = "#2790F9",
                               "ECA400" = "gray50","ECA401" = "gray50","ECA402" = "gray50",
                               "ECA547" = "gray50","ECA548"= "gray50" ,
                               "ECA549"= "gray50"))+
  facet_grid(.~condition, scales = "free") +
  base_theme +  theme(legend.position = "none",
                      panel.grid.major = element_blank(),
                      axis.line = element_line(colour = axis_color),
                      axis.title.x = element_blank()) +
  labs( y = paste0("Brood Size"))

ggsave(filename = "Plots/Dactionmycin_Brood_SWAP.png", height = 6, width = 10, dpi = 400)
ggsave(filename = "Plots/Dactionmycin_Brood_SWAP.pdf", height = 6, width = 10, dpi = 400)

######################################################################################################################## human data

human_swap <- data.table::fread(glue::glue("{etoposide_data}20170130_human-cell_reads.csv")) 

reads <- data.frame(t(human_swap))%>%
  dplyr::select(-X4)
colnames(reads) <- c("guide", "condition","amplicon","top2a_wt","top2a_edit","top2b_wt","top2b_edit")
reads <- reads[2:nrow(reads),]
reads$top2a_wt <- as.numeric(as.character(reads$top2a_wt))
reads$top2a_edit <- as.numeric(as.character(reads$top2a_edit))
reads$top2b_wt <- as.numeric(as.character(reads$top2b_wt))
reads$top2b_edit <- as.numeric(as.character(reads$top2b_edit))


# # # #  TOP2A
reads <- dplyr::mutate(reads, fraction_edited = ifelse(amplicon == "TOP2A", top2a_edit/ top2a_wt, top2b_edit/ top2b_wt))

top2a_o <- dplyr::filter(reads, guide == "No guide", 
                         # guide == "TOP2A", 
                         condition %in%c("Cellline Control D1"),
                         amplicon == "TOP2A")%>%  
  dplyr::select( condition, amplicon, fraction_edited)%>%
  tidyr::spread(condition, fraction_edited)

colnames(top2a_o) <- c("amplicon","0")

top2a_ctrl <- dplyr::filter(reads, guide != "No guide", 
                            guide == "TOP2A", 
                            condition %in%c("Early Time Point", "No Drug"),
                            amplicon == "TOP2A")%>%
  dplyr::select(guide, condition, amplicon, fraction_edited)%>%
  tidyr::spread(condition, fraction_edited)

colnames(top2a_ctrl) <- c("guide","amplicon","1","2")

top2a_ctrl <- top2a_ctrl%>%
  dplyr::left_join(.,top2a_o,by="amplicon")

top2a_expt <- reads%>%
  dplyr::filter( guide == "TOP2A", 
                 !condition %in%c("Early Time Point"),
                 amplicon == "TOP2A")%>%
  dplyr::select(guide, condition, amplicon, fraction_edited)%>%
  dplyr::left_join(.,top2a_ctrl, by = c("amplicon","guide"))%>%
  tidyr::gather(timepoint, base_fraction, -condition, -amplicon, -guide, -fraction_edited)%>%
  dplyr::mutate(fraction = ifelse(timepoint == "0", base_fraction,
                                  ifelse(timepoint == "1", base_fraction,
                                         ifelse(timepoint == "2", fraction_edited, NA))))%>%
  dplyr::mutate(fold = ifelse(timepoint == "0", base_fraction,
                              ifelse(timepoint == "1", base_fraction,
                                     ifelse(timepoint == "2", fraction_edited/base_fraction, NA))))%>%
  dplyr::mutate(l_frac = log10(fraction))

top2a_expt%>%
  dplyr::filter(timepoint != "0", grepl("Mid|No|Actin", condition))%>%
  dplyr::mutate(xlabels = ifelse(timepoint == 1, "Pre-treatment", "Post-treatment"))%>%
  ggplot(.)+
  aes(x = factor(xlabels, levels = c("Pre-treatment", "Post-treatment"), 
                 labels = c("Pre-treatment", "Post-treatment"), 
                 ordered = T), 
      y = (fraction), group = condition, color = condition)+
  geom_line(size = 2)+
  scale_color_manual(values=c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#990000", "#D55E00", "#CC79A7"),
                     name="Experimental\nCondition",
                     breaks=c("No Drug","Amsacrine Hydrochloride Mid", "Actinomycin D Low", "Etoposide Mid","Teniposide Mid","XK469 Mid"),
                     labels=c("No Drug","Amsacrine", "Dactinomycin", "Etoposide",  "Teniposide", "XK469"))+
  scale_y_log10(limits=c(0.0001,0.1))+
  base_theme+
  scale_x_discrete(expand = c(0,.1))+  
  theme(panel.grid.major = element_blank(),
        axis.line = element_line(colour = axis_color),
        axis.title.x = element_blank()) +
  labs(y = "Fraction Edited Reads")

ggsave(filename = "Plots/Human_SWAP.png", height = 6, width = 8, dpi = 400)
ggsave(filename = "Plots/Human_SWAP.pdf", height = 6, width = 8, dpi = 400)

######################################################################################################################## PopGene

# args <- c("II","10050022","12062611","Ce330_annotated.vcf.gz","WS245_exons.gff", "Etoposide")

load("Data/Etoposide_II_10050022-12062611_Diversity_Statistics.Rda")

qtl_start <- 10050022
qtl_end <- 12062611
top2_start <- 11875145
top2_end <- 11880588

td_df <- neutrality_df %>%
  dplyr::filter(statistic %in% c("Fay.Wu.H","Zeng.E","Tajima.D")) %>%
  dplyr::group_by(statistic) %>%
  dplyr::mutate(scaled_value = scale(value))

td_df %>%
  ggplot()+
  aes(x = WindowPosition/1e6, y = value, color = statistic)+
  geom_point(size = point_size, alpha = point_alpha)+
  scale_color_manual(values = c("#BE0032","#0067A5","#222222"))+
  facet_grid(statistic~., scales = "free")+
  base_theme +
  geom_vline(aes(xintercept = top2_start/1e6), color = "#E68FAC")+
  geom_vline(aes(xintercept = top2_end/1e6), color = "#E68FAC")+
  theme(panel.grid.major = element_blank(),
        legend.position = "none",
        axis.line = element_line(colour = axis_color),
        axis.title.y = element_blank()) +
  labs(x = "Genomic Position (Mb)")

ggsave(filename = "Plots/etoposide_popgen.png", height = 8, width = 12, dpi = 400)
ggsave(filename = "Plots/Human_SWAP.pdf", height = 6, width = 8, dpi = 400)

load("Data/etoposide_test_interval_10050022-12062611_Diversity_Statistics.Rda")

td_df <- neutrality_df %>%
  dplyr::filter(statistic %in% c("Fay.Wu.H","Zeng.E","Tajima.D")) %>%
  dplyr::group_by(statistic) %>%
  dplyr::mutate(scaled_value = scale(value))

td_df %>%
  dplyr::filter(WindowPosition > qtl_start) %>%
  dplyr::filter(WindowPosition < qtl_end) %>%
  ggplot()+
  aes(x = WindowPosition/1e6, y = value, color = statistic)+
  annotate("rect",
           xmin=dbt_start/1e6, 
           xmax=dbt_end/1e6, 
           ymin=-Inf, 
           ymax=Inf, 
           alpha=0.5, 
           fill="cyan") +
  geom_point(size = point_size, alpha = point_alpha)+
  scale_color_manual(values = c("#BE0032","#0067A5","#222222"))+
  facet_grid(statistic~., scales = "free")+
  base_theme +
  theme(panel.grid.major = element_blank(),
        legend.position = "none",
        axis.line = element_line(colour = axis_color),
        axis.title.x = element_blank()) +
  labs(x = "Genomic Position (Mb)")

td_df %>%
  dplyr::filter(WindowPosition > qtl_start) %>%
  dplyr::filter(WindowPosition < qtl_end) %>%
  ggplot()+
  aes(x = WindowPosition/1e6, y = value, color = statistic)+
  annotate("rect",
           xmin=dbt_start/1e6, 
           xmax=dbt_end/1e6, 
           ymin=-Inf, 
           ymax=Inf, 
           alpha=0.5, 
           fill="cyan") +
  geom_point(size = point_size, alpha = point_alpha)+
  scale_color_manual(values = col_blind_colors)+
  facet_grid(statistic~., scales = "free")+
  base_theme +
  theme(panel.grid.major = element_blank(),
        legend.position = "none",
        axis.line = element_line(colour = axis_color),
        axis.title.x = element_blank()) +
  labs(x = "Genomic Position (Mb)")+
  xlim(c((dbt_start-50000)/1e6, (dbt_end+50000)/1e6))




