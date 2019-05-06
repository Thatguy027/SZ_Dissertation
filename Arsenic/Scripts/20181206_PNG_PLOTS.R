# set to location of files
main.dir <- "~/Dropbox/AndersenLab/LabFolders/Stefan/Manuscripts/Arsenic/"
data.dir <- paste0(main.dir,"Data/")
final.dir <- paste0(main.dir,"Final_tables/")
plot.dir <- paste0(main.dir,"Plots/")


resub.dir <- paste0(main.dir,"Resubmission/")

resub.data <- paste0(resub.dir,"Data/")
resub.plots <- paste0(resub.dir,"Plots/")
resub.gwas <- paste0(resub.dir,"cegwas2_GWA_Results/")
resub.plots.png <- paste0(resub.dir,"Plots/PNG/")

dir.create(resub.data)
dir.create(resub.plots)
dir.create(resub.plots.png)

script.dir <- paste0(resub.dir,"Scripts/")

source(paste0(script.dir,"20181009_processing_functions.R"))

########################################################################################################################
# DOSE RESPONSE EXPERIMENT - PLOTS - START
########################################################################################################################
dr_traits <- data.table::fread(glue::glue("{resub.data}Supplemental_Data1_DR-traits_PCA_per_Condition.tsv"))
dr_trait_cor <- data.table::fread(glue::glue("{resub.data}Supplemental_Data2_DR-trait-correlation.tsv"))
dr_loadings <- data.table::fread(glue::glue("{resub.data}Supplemental_Data3_DR-PC-Loadings.tsv"))

dr_complete_pcs <- data.table::fread(glue::glue("{resub.data}Supplemental_Data4_DR-PCA_Complete_Experiment.tsv"))
dr_complete_loadings <- data.table::fread(glue::glue("{resub.data}Supplemental_Data5_DR-PCA_Complete_Experiment.tsv"))


# plot DR - correlations and loadings - per condition

for(cond in 1:length(unique(dr_traits$Condition))){
  
  analysis_condition <- unique(dr_traits$Condition)[cond]
  
  if(analysis_condition == "arsenic1"){
    tidy_condition <- "250 µM Arsenic"
    tidy_subplot <- "B"
    plot_space <- 15
  } else if (analysis_condition == "arsenic2") {
    tidy_condition <- "500 µM Arsenic"
    tidy_subplot <- "C"
    plot_space <- 15
  } else if (analysis_condition == "arsenic3") {
    tidy_condition <- "1000 µM Arsenic"
    tidy_subplot <- "D"
    plot_space <- 15
  } else if (analysis_condition == "arsenic4") {
    tidy_condition <- "2000 µM Arsenic"
    tidy_subplot <- "E"
    plot_space <- 15
  } else {
    tidy_condition <- "Water"
    tidy_subplot <- "A"
    plot_space <- 40
  }
  
  
  cond_trait_df <- dr_traits %>%
    dplyr::filter(Condition == analysis_condition)
  
  cond_cor_df <- dr_trait_cor %>%
    dplyr::filter(Condition == analysis_condition)
  
  cond_loadings_df <- dr_loadings %>%
    dplyr::filter(Condition == analysis_condition)
  
  
  trait_order <- dplyr::filter(cond_cor_df, 
                               trait_a == "PC1",
                               !grepl("PC", trait_b)) %>%
    dplyr::arrange(desc(abs(trait_cor)))
  
  dr_cor_plot <- cond_cor_df %>%
    dplyr::filter(!grepl("PC", trait_b),
                  !grepl("PC", trait_a)) %>%
    ggplot() +
    aes(x = factor(trait_a,
                   levels = trait_order$trait_b), 
        y = factor(trait_b,
                   levels = trait_order$trait_b), 
        fill = trait_cor)+
    geom_tile()+
    scale_fill_gradient2(low = "blue", 
                         high = "yellow", 
                         mid = "purple", 
                         midpoint = 0, 
                         limit = c(-1,1), 
                         space = "Lab", 
                         name="Pearson\nCorrelation") +
    arsenic.theme +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    labs(title = tidy_condition)
  
  dr_pc_loading_plot <- ggplot(cond_loadings_df) +
    aes(x = factor(Trait,
                   levels = trait_order$trait_b), 
        y = factor(PC, levels = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")), 
        fill = loading)+
    geom_tile()+
    scale_fill_gradient2(low = "blue",
                         high = "yellow",
                         mid = "purple", 
                         space = "Lab", 
                         name="Loading") +
    arsenic.theme +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),plot.margin = margin(10, 10, 10, plot_space))
  
  
  cowplot::plot_grid(dr_cor_plot, dr_pc_loading_plot, 
                     labels = tidy_subplot, 
                     ncol = 1, 
                     align = 'v',
                     vjust = 1, 
                     label_size = 18)
  
  ggsave(glue::glue("{resub.plots.png}FS1{tidy_subplot}_Trait-correlation_loadings.png"), 
         height = 12, 
         width = 20, dpi = 300)
  
}

# plot dose response
dr_traits%>%
  dplyr::filter(Trait %in% c("median.TOF","norm.n")) %>%
  dplyr::bind_rows(.,dr_complete_pcs) %>%
  dplyr::filter(Trait %in% c("median.TOF","norm.n", "PC1")) %>%
  dplyr::mutate(flip_pc = ifelse(Trait == "PC1", -Value, Value)) %>%
  dplyr::mutate(Tidy_trait = ifelse(Trait == "median.TOF", "Animal Size", 
                                    ifelse(Trait == "norm.n","Brood Size", "PC1"))) %>%
  ggplot()+
  aes(x = factor(Condition, 
                 levels = c("Water", "arsenic1", "arsenic2", "arsenic3","arsenic4"),
                 labels = c("0", "250", "500", "1000","2000")), 
      y = flip_pc, 
      fill = Strain)+
  geom_boxplot(outlier.colour = NA)+
  facet_grid(Tidy_trait ~ ., scales = "free_y") +
  scale_fill_manual(values = c("blue", "cadetblue3", "hotpink3", "orange"), name = "Strain")+
  theme_bw()+
  labs(x = "Arsenic (µM)") + 
  arsenic.theme +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 16),
        axis.title.y = element_blank())

ggsave(glue::glue("{resub.plots.png}FS2_DR.png"), 
       height = 10, 
       width = 10, dpi = 300)


# plot effect sizes for DR
effects_df <- data.table::fread(glue::glue("{resub.data}Supplemental_Data6_DR-PCA_Effect-Sizes.tsv"))
h2_df <- data.table::fread(glue::glue("{resub.data}Supplemental_Data7_DR-Heritability_estimates.tsv"))

plt_df <- h2_df%>%
  dplyr::filter(condition != "Water") %>%
  dplyr::left_join(.,effects_df, by = c("condition", "trait")) %>%
  dplyr::filter(strain1 == "ALL") %>%
  dplyr::mutate(cpoint = ifelse(trait == "PC1", 1,
                                ifelse(trait == "norm.n",3,
                                       ifelse(trait == "median.TOF", 2, 4))),
                tidy_cond = factor(condition, 
                                   levels = c( "arsenic1", "arsenic2", "arsenic3","arsenic4"),
                                   labels = c( "250 µM", "500 µM", "1000 µM","2000 µM"))) 

ggplot()+
  aes(x = partial.omegasq, y = DR_H2)+
  geom_point(color = "gray50", 
             data = dplyr::filter(plt_df, !trait %in% c("PC1","norm.n","median.TOF")))+
  geom_point(aes(color = factor(cpoint, labels=c("PC1","Animal Length","Brood Size")),
                 size = factor(cpoint, labels=c("PC1","Animal Length","Brood Size")),
                 fill = factor(cpoint, labels=c("PC1","Animal Length","Brood Size")),
                 shape = factor(cpoint, labels=c("PC1","Animal Length","Brood Size"))),
             data = dplyr::filter(plt_df, trait %in% c("PC1","norm.n","median.TOF"))) +
  facet_grid(~tidy_cond) +
  theme_bw(18)+
  scale_shape_manual(values=c(21,21,21,1))+
  scale_color_manual(values=c("black","black","black","gray50"))+
  scale_fill_manual(values=c("red","hotpink3","cadetblue3","gray50"))+
  scale_size_manual(values=c(3,3,3,1))+
  labs(x = expression(italic(omega)[p]^2), y = "Broad-Sense Heritability",
       color="Trait",
       size = "Trait",
       fill="Trait",
       shape="Trait")+
  theme(strip.background = element_blank(),
        legend.position = "top",
        panel.grid.minor = element_blank())

ggsave(glue::glue("{resub.plots.png}FS3_DR_H2_v_Effect_Sizes.png"), 
       height = 8, 
       width = 16, dpi = 300)

########################################################################################################################
# LINKAGE MAPPING - PLOTS - START
########################################################################################################################

linkage_pheno <- data.table::fread(glue::glue("{resub.data}Supplemental_Data8_Linkage_Phenotypes.tsv"))
riail_trait_cor <- data.table::fread(glue::glue("{resub.data}Supplemental_Data9_Linkage_trait_correlations.tsv"))
riail_pc_loadings <- data.table::fread(glue::glue("{resub.data}Supplemental_Data10_Linkage_PC_Loadings.tsv"))

trait_order <- riail_pc_loadings %>%
  dplyr::filter(PC == "PC1") %>%
  dplyr::arrange(desc(abs(loading)))

RIAIL_cor_plot <- riail_trait_cor %>%
  dplyr::filter(!grepl("PC", trait_b),
                !grepl("PC", trait_a)) %>%
  ggplot() +
  aes(x = factor(trait_a,
                 levels = trait_order$Trait), 
      y = factor(trait_b,
                 levels = trait_order$Trait), 
      fill = trait_cor)+
  geom_tile()+
  scale_fill_gradient2(low = "blue", 
                       high = "yellow", 
                       mid = "purple", 
                       midpoint = 0, 
                       limit = c(-1,1), 
                       space = "Lab", 
                       name="Pearson\nCorrelation") +
  arsenic.theme +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())


RIAIL_pc_loading_plot <- riail_pc_loadings%>%
  dplyr::filter(PC %in% c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")) %>%
  ggplot() + 
  aes(x = factor(Trait,
                 levels = trait_order$Trait), 
      y = factor(PC, levels = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")), 
      fill = loading)+
  geom_tile()+
  scale_fill_gradient2(low = "blue",
                       high = "yellow",
                       mid = "purple",
                       space = "Lab",
                       name="Loading") +
  arsenic.theme +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())


cowplot::plot_grid(RIAIL_cor_plot, 
                   RIAIL_pc_loading_plot, 
                   labels = "AUTO", 
                   ncol = 1, 
                   align = 'v',
                   vjust = 1, 
                   label_size = 18)

ggsave(glue::glue("{resub.plots.png}FS4_RIAIL_TRAIT_COR_LOADINGS.png"), 
       height = 12, 
       width = 20, dpi = 300)

# Plot Linkage Mapping
linkage_pheno <- data.table::fread(glue::glue("{resub.data}Supplemental_Data8_Linkage_Phenotypes.tsv"))
lods <- data.table::fread(glue::glue("{resub.data}Supplemental_Data11_Linkage_Annotated_LODs.tsv"))

arsenic_linkage <- linkage_pheno %>%
  dplyr::rename(strain = Strain, phenotype = Value, trait = Trait) %>%
  dplyr::select(-Condition)

data("N2xCB4856cross")

blankcross <- N2xCB4856cross

# Plot Linkage Mapping - PxG
arsenic_cross <- linkagemapping::mergepheno(blankcross, arsenic_linkage, set = 2)

pxg_linkage_PC <- pxgplot_edit(arsenic_cross, dplyr::filter(lods, trait == ".PC1"))+ 
  scale_x_discrete(breaks=c("N2", "CB4856"),labels=c("N2", "CB4856"))+
  labs(y = paste0("PC1"))

pxg_linkage_length <- pxgplot_edit(arsenic_cross, dplyr::filter(lods,  trait == ".median.TOF"))+ 
  scale_x_discrete(breaks=c("N2", "CB4856"),labels=c("N2", "CB4856"))+
  labs(y = paste0("Animal Length"))

pxg_linkage_brood <- pxgplot_edit(arsenic_cross, dplyr::filter(lods,  trait == ".norm.n"))+ 
  scale_x_discrete(breaks=c("N2", "CB4856"),labels=c("N2", "CB4856"))+
  labs(y = paste0("Brood Size"))

plot_grid(pxg_linkage_PC,
          pxg_linkage_length, 
          pxg_linkage_brood, labels = "AUTO", ncol = 1, align = 'v', label_size = 18)

ggsave(glue::glue("{resub.plots.png}FS5_Linkage_PxG.png"), 
       height = 8, 
       width = 14, dpi = 300)

lods %>%
  dplyr::filter(trait %in% c(".median.TOF",".norm.n", ".PC1")) %>%
  maxlodplot_edit2() +
  arsenic.theme +
  theme(axis.title.x = element_text(size=16),axis.title.y = element_text(size=16),
        strip.background.x = element_blank())

ggsave(glue::glue("{resub.plots.png}FS6_Linkage_PC1-brood-length.png"), 
       height = 4, 
       width = 12, dpi = 300)
########################################################################################################################
# PLOT LINKAGE MAPPING - END
########################################################################################################################

riail_h2_df <- data.table::fread(glue::glue("{resub.data}Supplemental_Data12_Linkage_H2_ESTIMATES.tsv"))

broad_estimates <- ggplot(riail_h2_df)+
  geom_point(aes(x = H2_realized, y = H2_expectation, color = QTL_II))+
  scale_color_manual(values = c("red","black")) +
  theme_classic(16) +
  geom_abline(slope=1,intercept = 0, linetype = 2, alpha = 0.5)+
  ggrepel::geom_label_repel(aes(x = H2_realized, y = H2_expectation, label = tidy_trait),
                            min.segment.length = 0.01,
                            nudge_y=0.2,
                            segment.alpha = 0.5,
                            data = dplyr::filter(riail_h2_df, trait %in% c("PC1", "norm.n","median.TOF")) %>%
                              dplyr::mutate(tidy_ve = round(QTL_II_VE/H2_expectation, digits = 3),
                                            tidy_trait = case_when(
                                              trait == "PC1" ~ glue::glue("PC 1"),
                                              trait == "norm.n" ~ glue::glue("Brood Size"),
                                              trait == "median.TOF" ~ glue::glue("Animal Length"),
                                              TRUE ~ "NA"
                                            ))) +
  labs(x = "H2: A", y = "H2: E(A)", color = "II:7296342 QTL")+
  xlim(c(0,1))+
  ylim(c(0,1))

narrow_estimates <- ggplot(riail_h2_df)+
  geom_point(aes(x = h2_realized, y = h2_expectation, color = QTL_II))+
  scale_color_manual(values = c("red","black")) +
  theme_classic(16) +
  geom_abline(slope=1,intercept = 0, linetype = 2, alpha = 0.5)+
  ggrepel::geom_label_repel(aes(x = h2_realized, y = h2_expectation, label = tidy_trait),
                            min.segment.length = 0.01,
                            nudge_y=0.2,
                            segment.alpha = 0.5,
                            data = dplyr::filter(riail_h2_df, trait %in% c("PC1", "norm.n","median.TOF")) %>%
                              dplyr::mutate(tidy_ve = round(QTL_II_VE/H2_expectation, digits = 3),
                                            tidy_trait = case_when(
                                              trait == "PC1" ~ glue::glue("PC 1"),
                                              trait == "norm.n" ~ glue::glue("Brood Size"),
                                              trait == "median.TOF" ~ glue::glue("Animal Length"),
                                              TRUE ~ "NA"
                                            ))) +
  labs(x = "H2: A", y = "H2: E(A)", color = "II:7296342 QTL")+
  xlim(c(0,1))+
  ylim(c(0,1))

broad_narrow_realized <- ggplot(riail_h2_df)+
  geom_point(aes(x = H2_realized, y = h2_realized, color = QTL_II))+
  scale_color_manual(values = c("red","black")) +
  theme_classic(16) +
  geom_abline(slope=1,intercept = 0, linetype = 2, alpha = 0.5)+
  ggrepel::geom_label_repel(aes(x = H2_realized, y = h2_realized, label = tidy_trait),
                            min.segment.length = 0.01,
                            nudge_y=0.2,
                            segment.alpha = 0.5,
                            data = dplyr::filter(riail_h2_df, trait %in% c("PC1", "norm.n","median.TOF")) %>%
                              dplyr::mutate(tidy_ve = round(QTL_II_VE/H2_expectation, digits = 3),
                                            tidy_trait = case_when(
                                              trait == "PC1" ~ glue::glue("PC 1"),
                                              trait == "norm.n" ~ glue::glue("Brood Size"),
                                              trait == "median.TOF" ~ glue::glue("Animal Length"),
                                              TRUE ~ "NA"
                                            ))) +
  labs(x = "H2: A", y = "h2: A", color = "II:7296342 QTL") +
  xlim(c(0,1))+
  ylim(c(0,1))

broad_narrow_expected <- ggplot(riail_h2_df)+
  geom_point(aes(x = H2_expectation, y = h2_expectation, color = QTL_II))+
  scale_color_manual(values = c("red","black")) +
  theme_classic(16) +
  ggrepel::geom_label_repel(aes(x = H2_expectation, y = h2_expectation, label = tidy_trait),
                            min.segment.length = 0.01,
                            nudge_y=0.2,
                            segment.alpha = 0.5,
                            data = dplyr::filter(riail_h2_df, trait %in% c("PC1", "norm.n","median.TOF")) %>%
                              dplyr::mutate(tidy_ve = round(QTL_II_VE/H2_expectation, digits = 3),
                                            tidy_trait = case_when(
                                              trait == "PC1" ~ glue::glue("PC 1"),
                                              trait == "norm.n" ~ glue::glue("Brood Size"),
                                              trait == "median.TOF" ~ glue::glue("Animal Length"),
                                              TRUE ~ "NA"
                                            ))) +
  geom_abline(slope=1,intercept = 0, linetype = 2, alpha = 0.5)+
  labs(x = "H2: E(A)", y = "h2: E(A)", color = "II:7296342 QTL") +
  xlim(c(0,1))+
  ylim(c(0,1))

legend <- cowplot::get_legend(broad_estimates)

plots <- cowplot::plot_grid(broad_estimates + theme(legend.position = "none"),
                            narrow_estimates + theme(legend.position = "none"),
                            broad_narrow_realized + theme(legend.position = "none"),
                            broad_narrow_expected + theme(legend.position = "none"),
                            labels = "AUTO", 
                            align = 'v', 
                            label_size = 18, 
                            ncol = 2)

linkage_h2_figure <- plot_grid( plots, legend, rel_widths = c(3, .3))

ggsave(plot = linkage_h2_figure,
       glue::glue("{resub.plots.png}FS7_Linkage_heritability_estimates.png"), 
       height = 12, 
       width = 16, dpi = 300)

########################################################################################################################
# LINKAGE MAPPING HERITABILITY ESTIMATES - END
########################################################################################################################

# LINKAGE QTL SUMMARY
lods <- data.table::fread(glue::glue("{resub.data}Supplemental_Data10_Linkage_Annotated_LODs.tsv"))

linkage_peaks <- na.omit(lods)
linkage_peaks$trait <- gsub("^.","",linkage_peaks$trait)

goodtraits <- linkage_peaks %>%
  ungroup() %>%
  dplyr::group_by(chr) %>%
  dplyr::arrange(chr, ci_l_pos)

goodtraits$pos <- as.numeric(goodtraits$pos)
goodtraits$ci_l_pos <- as.numeric(goodtraits$ci_l_pos)
goodtraits$ci_r_pos <- as.numeric(goodtraits$ci_r_pos)
goodtraits$lod <- as.numeric(goodtraits$lod)

#Set chromosome boundaries
newrows <- goodtraits[1,]
newrows$pos <- as.numeric(newrows$pos)
newrows$ci_l_pos <- as.numeric(newrows$ci_l_pos)
newrows$ci_r_pos <- as.numeric(newrows$ci_r_pos)
newrows$lod <- as.numeric(newrows$lod)

newrows[1,] = c(NA,"I",1,"pc1",1,NA,NA,NA,NA,NA,0,NA,14972282)
newrows[2,] = c(NA,"II",1,"pc1",NA,NA,NA,NA,NA,NA,0,NA,15173999)
newrows[3,] = c(NA,"III",1,"pc1",NA,NA,NA,NA,NA,NA,0,NA,13829314)
newrows[4,] = c(NA,"IV",1,"pc1",NA,NA,NA,NA,NA,NA,0,NA,17450860)
newrows[5,] = c(NA,"V",1,"pc1",NA,NA,NA,NA,NA,NA,0,NA,20914693)
newrows[6,] = c(NA,"X",1,"pc1",NA,NA,NA,NA,NA,NA,0,NA,17748731)
newrows$pos <- as.numeric(newrows$pos)
newrows$ci_l_pos <- as.numeric(newrows$ci_l_pos)
newrows$ci_r_pos <- as.numeric(newrows$ci_r_pos)
newrows$lod <- as.numeric(newrows$lod)
newrows$trait <- ""
newrows <- dplyr::ungroup(newrows)

goodtraits$pos <- as.numeric(goodtraits$pos)
goodtraits$ci_l_pos <- as.numeric(goodtraits$ci_l_pos)
goodtraits$ci_r_pos <- as.numeric(goodtraits$ci_r_pos)
goodtraits$lod <- as.numeric(goodtraits$lod)

linkage_summary_plot<-goodtraits%>%
  dplyr::ungroup()%>%
  ggplot()+
  aes(x=pos/1E6, y=trait)+
  theme_bw() +
  scale_fill_gradient(high = "#D7263D", low = "#0072B2",
                      name = "LOD")+
  scale_color_gradient(high = "#D7263D", low = "#0072B2",
                       name ="LOD")+
  geom_segment(aes(x = ci_l_pos/1e6, y = trait, xend = ci_r_pos/1e6, yend = trait, color = lod), size = 2, alpha = 1) +
  geom_segment(data=newrows,aes(x = 0, y = trait, xend = ci_r_pos/1e6, yend = trait), size = 2.5, alpha = 0) +
  geom_point(aes(fill=lod),colour = "black",size = 2, alpha = 1, shape = 25)+
  xlab("Genomic Position (Mb)") + ylab("") +
  theme(strip.background = element_rect(colour = "black", fill = "white",
                                        size = 0.75, linetype = "solid")) +
  theme_bw(18) + 
  theme(strip.background = element_blank(),
        panel.background = element_blank(),
        strip.text.y = element_blank(),
        panel.spacing = unit(0.5, "lines"))+
  facet_grid(. ~ chr, scales = "free", space = "free")+
  ggplot2::labs(x = "Genomic Position (Mb)")

ggsave(plot = linkage_summary_plot,
       glue::glue("{resub.plots.png}FS8_Linkage_QTL_Summary.png"), 
       height = 14, 
       width = 16, dpi = 300)

########################################################################################################################
# LINKAGE MAPPING - PLOTS - START
########################################################################################################################

lods <- data.table::fread(glue::glue("{resub.data}Supplemental_Data11_Linkage_Annotated_LODs.tsv"))
nil_pheno <- data.table::fread(glue::glue("{resub.data}Supplemental_Data13_NIL_control_regressed.tsv"))
linkage_pheno <- data.table::fread(glue::glue("{resub.data}Supplemental_Data8_Linkage_Phenotypes.tsv"))

lod_plot <- lods %>%
  dplyr::filter(trait == ".PC1") %>%
  maxlodplot_edit() +
  theme_bw(18) +
  arsenic.theme +
  theme(plot.title = element_blank(),
        strip.background.x = element_blank(),
        axis.text.x = ggplot2::element_text(size = 12),
        axis.text.y = ggplot2::element_text(size = 12),
        panel.grid = element_blank())

arsenic_linkage <- linkage_pheno %>%
  dplyr::rename(strain = Strain, phenotype = Value, trait = Trait) %>%
  dplyr::select(-Condition)

data("N2xCB4856cross")

blankcross <- N2xCB4856cross

# PC trait correlation plot for RIAILs
riail_to_lm <- linkage_pheno %>%
  dplyr::select(-Condition) %>%
  dplyr::ungroup() %>%
  dplyr::filter(Trait %in% c("norm.n", "median.TOF", "PC1")) %>%
  tidyr::spread(Trait, Value) %>%
  tidyr::gather(CorTrait, value, -Strain, -PC1) 


lm_ls <- list()
for(tra in 1:length(unique(riail_to_lm$CorTrait))) {
  riail_trait_to_lm <- riail_to_lm %>%
    dplyr::filter(CorTrait == unique(riail_to_lm$CorTrait)[tra]) 
  
  lm_ls[[tra]] <- riail_trait_to_lm %>%
    tidyr::nest(CorTrait) %>%
    mutate(model = purrr::map(data, ~ glm(PC1 ~ value, data = .x)),
           adj.r.squared = purrr::map_dbl(model, ~ signif(1 - (.x$deviance/.$null.deviance), 3)),
           intercept = purrr::map_dbl(model, ~ signif(.x$coef[[1]],3)),
           slope = purrr::map_dbl(model, ~ signif(.x$coef[[2]], 3)),
           pvalue = purrr::map_dbl(model, ~ signif(summary(.x)$coef[2,4], 3)) 
    )    %>%
    select(-data, -model) %>% 
    left_join(riail_to_lm)
}

riail_lm <- dplyr::bind_rows(lm_ls)

riail_p_cor_plot <- riail_lm %>%
  ggplot()+
  aes(x = value, y = PC1,color = CorTrait)+
  geom_hline(yintercept = 0, color = "gray50", alpha = 0.5) +
  geom_vline(xintercept = 0, color = "gray50", alpha = 0.5) +
  geom_point(aes( )) + 
  geom_smooth(method='glm',formula=y~x, se = F) +
  scale_color_manual(values = c("hotpink3", "cadetblue3"), 
                     breaks = c("median.TOF","norm.n"),
                     labels = c("Length", "Brood Size"))  +
  theme_bw(18) +
  labs(x = "Raw Phenotype Value", color = NULL, y = "PC 1") +
  arsenic.theme +
  geom_text(aes(50, -15, label = paste0("y = ", unique(slope), 
                                        "x + ", unique(intercept), 
                                        "\n R^2 = ", unique(adj.r.squared), "; p = ", unique(pvalue))),
            data = dplyr::filter(riail_lm, CorTrait == "norm.n")) +
  geom_text(aes(-55, 15, label = paste0("y = ", unique(slope), 
                                        "x + ", unique(intercept), 
                                        "\n R^2 = ", unique(adj.r.squared), "; p = ", unique(pvalue))),
            data = dplyr::filter(riail_lm, CorTrait == "median.TOF")) +
  theme(panel.grid = element_blank(),
        legend.position = "top",
        legend.background = element_rect(size = 0.5, color = "black"))


NIL_pc <- nil_pheno%>%
  dplyr::filter(Trait == "PC1",
                Strain%in%c("N2","CB4856","ECA414","ECA434"))%>%
  dplyr::mutate(strain1 = factor(Strain, levels = c("N2","CB4856","ECA414","ECA434","ECA581",
                                                    "ECA582","ECA589","ECA590","ECA591"),
                                 labels =c("N2",
                                           "CB4856", 
                                           "CB4856>N2\n II:5.75 - 8.02Mb", 
                                           "CB4856>N2\n II:7.83 - 9.66Mb",
                                           "N2(C78S)\nA", "N2(C78S)\nB",
                                           "CB4856(S78C)\nA","CB4856(S78C)\nB",
                                           "CB4856(S78C)\nC")))%>%
  ggplot(.) +
  aes(x = factor(strain1), 
      y = -Value, 
      fill=Strain) +
  geom_beeswarm(alpha = .4,priority = "density",cex = 1.2)+
  geom_boxplot(outlier.colour = NA, alpha = 0.7)+
  scale_fill_manual(values = c("N2" = "orange","CB4856" = "blue",
                               "ECA414" = "blue","ECA434" = "blue",
                               "ECA581" = "gray50","ECA582" = "gray50",
                               "ECA589" = "gray50","ECA590" = "gray50",
                               "ECA591" = "gray50"))+
  theme_bw()+
  arsenic.theme +
  labs( y = paste0("PC 1"))+
  theme(legend.position="none",
        strip.background = element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        panel.grid.minor = element_blank())


bottom_row <- plot_grid(lemon::reposition_legend(riail_p_cor_plot, 'bottom right'), 
                        NIL_pc, 
                        labels = c("B", "C"), 
                        rel_widths = c(0.75,1),
                        align = 'v', 
                        label_size = 18, 
                        ncol = 2)

figure1 <- plot_grid(lod_plot, 
                     bottom_row, 
                     labels = c("A", ""), 
                     rel_heights = c(1,1),
                     align = 'v', 
                     label_size = 18, 
                     ncol = 1)


ggsave(plot = figure1,
       glue::glue("{resub.plots.png}F1_Linkage_NILs.png"), 
       height = 8, 
       width = 16, dpi = 300)

cor.test(x = dplyr::filter(nil_pheno, Strain %in% c("N2", "CB4856"), Trait == "PC1") %>% dplyr::pull(Strain),
         y = dplyr::filter(nil_pheno, Strain %in% c("N2", "CB4856"), Trait == "PC1") %>% dplyr::pull(Value),
         method = "spearman")

aov_res <- aov(dplyr::filter(nil_pheno, Strain %in% c("N2", "CB4856", "ECA414", "ECA434"), Trait == "PC1") %>% dplyr::pull(Value) ~ 
                 dplyr::filter(nil_pheno, Strain %in% c("N2", "CB4856", "ECA414", "ECA434"), Trait == "PC1") %>% dplyr::pull(Strain))
summary(aov_res)
tuk <- TukeyHSD(aov_res)

sjstats::anova_stats(car::Anova(aov(Value ~ Strain, data = dplyr::filter(nil_pheno, Strain %in% c("N2", "CB4856"), Trait == "PC1"))))
sjstats::anova_stats(car::Anova(aov(Value ~ Strain, data = dplyr::filter(nil_pheno, Strain %in% c("N2", "ECA414"), Trait == "PC1"))))
sjstats::anova_stats(car::Anova(aov(Value ~ Strain, data = dplyr::filter(nil_pheno, Strain %in% c("N2", "ECA434"), Trait == "PC1"))))
# N2 - CB
# 3.66
# N2 - left NIL
# 2.97 - 81% recapitulation
# N2 - right NIL
# 3.081 - 84% recapitulation

# Brood size and animal length NILs

nil_pheno%>%
  dplyr::filter(Trait %in% c("norm.n", "median.TOF"),
                Strain%in%c("N2","CB4856","ECA414","ECA434"))%>%
  dplyr::mutate(strain1 = factor(Strain, levels = c("N2","CB4856","ECA414","ECA434","ECA581",
                                                    "ECA582","ECA589","ECA590","ECA591"),
                                 labels =c("N2",
                                           "CB4856", 
                                           "CB4856>N2\n II:5.75 - 8.02Mb", 
                                           "CB4856>N2\n II:7.83 - 9.66Mb",
                                           "N2(C78S)\nA", "N2(C78S)\nB",
                                           "CB4856(S78C)\nA","CB4856(S78C)\nB",
                                           "CB4856(S78C)\nC")),
                tidy_trait = factor(Trait, levels = c("norm.n", "median.TOF"), 
                                    labels = c("Brood Size", "Animal Length")))%>%
  ggplot(.) +
  aes(x = factor(strain1), 
      y = Value, 
      fill=Strain) +
  geom_beeswarm(alpha = .4,priority = "density",cex = 1.2)+
  geom_boxplot(outlier.colour = NA, alpha = 0.7)+
  scale_fill_manual(values = c("N2" = "orange","CB4856" = "blue",
                               "ECA414" = "blue","ECA434" = "blue",
                               "ECA581" = "gray50","ECA582" = "gray50",
                               "ECA589" = "gray50","ECA590" = "gray50",
                               "ECA591" = "gray50"))+
  theme_bw()+
  facet_grid(tidy_trait~.)+
  arsenic.theme +
  labs( y = paste0("PC 1"))+
  theme(legend.position="none",
        strip.background = element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        panel.grid.minor = element_blank())

ggsave(glue::glue("{resub.plots.png}FS9_NIL_length_brood_recapitulation.png"), 
       height = 6, 
       width = 10, dpi = 300)


# NIL TRAIT CORRELATION AND LOADINGS

nil_pheno <- data.table::fread(glue::glue("{resub.data}Supplemental_Data13_NIL_control_regressed.tsv"))
nil_trait_cor <- data.table::fread(glue::glue("{resub.data}Supplemental_Data14_NIL_Trait_Correlations.tsv"))
nil_pc_loading <- data.table::fread(glue::glue("{resub.data}Supplemental_Data15_NIL_PC_Loadings.tsv"))

trait_order <- nil_pc_loading %>%
  dplyr::filter(PC == "PC1") %>%
  dplyr::arrange(desc(abs(loading)))

NIL_cor_plot <- nil_trait_cor %>%
  dplyr::filter(!grepl("PC", trait_b),
                !grepl("PC", trait_a)) %>%
  ggplot() +
  aes(x = factor(trait_a,
                 levels = trait_order$Trait), 
      y = factor(trait_b,
                 levels = trait_order$Trait), 
      fill = trait_cor)+
  geom_tile()+
  scale_fill_gradient2(low = "blue", 
                       high = "yellow", 
                       mid = "purple", 
                       midpoint = 0, 
                       limit = c(-1,1), 
                       space = "Lab", 
                       name="Pearson\nCorrelation") +
  arsenic.theme +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())


NIL_pc_loading_plot <- nil_pc_loading%>%
  dplyr::filter(PC %in% c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")) %>%
  ggplot() + 
  aes(x = factor(Trait,
                 levels = trait_order$Trait), 
      y = factor(PC, levels = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")), 
      fill = loading)+
  geom_tile()+
  scale_fill_gradient2(low = "blue",
                       high = "yellow",
                       mid = "purple",
                       space = "Lab",
                       name="Loading") +
  arsenic.theme +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())


cowplot::plot_grid(NIL_cor_plot, 
                   NIL_pc_loading_plot, 
                   labels = "AUTO", 
                   ncol = 1, 
                   align = 'v',
                   vjust = 1, 
                   label_size = 18)

ggsave(glue::glue("{resub.plots.png}FS10_NIL_SWAP_TRAIT_COR_LOADINGS.png"), 
       height = 12, 
       width = 20, dpi = 300)


# PC trait correlation plot

nil_to_lm <- nil_pheno %>%
  dplyr::filter(Strain %in% c("N2", "CB4856", "ECA434", "ECA414"))  %>%
  dplyr::filter(Trait %in% c("norm.n", "median.TOF", "PC1")) %>%
  tidyr::spread(Trait, Value) %>%
  tidyr::gather(CorTrait, value, -(Strain:u_id), -PC1) 

lm_ls <- list()
for(tra in 1:length(unique(nil_to_lm$CorTrait))) {
  nil_trait_to_lm <- nil_to_lm %>%
    dplyr::filter(CorTrait == unique(nil_to_lm$CorTrait)[tra]) 
  
  lm_ls[[tra]] <- nil_trait_to_lm %>%
    tidyr::nest(CorTrait) %>%
    mutate(model = purrr::map(data, ~ glm(PC1 ~ value, data = .x)),
           adj.r.squared = purrr::map_dbl(model, ~ signif(1 - (.x$deviance/.$null.deviance), 3)),
           intercept = purrr::map_dbl(model, ~ signif(.x$coef[[1]],3)),
           slope = purrr::map_dbl(model, ~ signif(.x$coef[[2]], 3)),
           pvalue = purrr::map_dbl(model, ~ signif(summary(.x)$coef[2,4], 3)) 
    )    %>%
    select(-data, -model) %>% 
    left_join(nil_to_lm)
}

nil_lm <- dplyr::bind_rows(lm_ls)

nil_plot <- nil_lm %>%
  ggplot()+
  aes(x = value, y = PC1,color = CorTrait)+
  geom_hline(yintercept = 0, color = "gray50", alpha = 0.5) +
  geom_vline(xintercept = 0, color = "gray50", alpha = 0.5) +
  geom_point(aes( )) + 
  geom_smooth(method='glm',formula=y~x, se = F) +
  scale_color_manual(values = c("hotpink3", "cadetblue3"), 
                     breaks = c("median.TOF","norm.n"),
                     labels = c("Length", "Brood\nSize"))  +
  theme_bw(18) +
  labs(x = "Phenotype Value", color = "Trait", y = "PC 1") +
  arsenic.theme +
  geom_text(aes(45, 10, 
                label = paste("Adj R2 = ", adj.r.squared, "\n",
                              "Intercept =",intercept, "\n",
                              "Slope =", slope, "\n",
                              "P =", pvalue)), data = dplyr::filter(nil_lm, CorTrait == "norm.n")) +
  geom_text(aes(-35, -10, 
                label = paste("Adj R2 = ", adj.r.squared, "\n",
                              "Intercept =",intercept, "\n",
                              "Slope =", slope, "\n",
                              "P =", pvalue)), data = dplyr::filter(nil_lm, CorTrait == "median.TOF"))+
  theme(panel.grid = element_blank())


ggsave(plot = nil_plot,
       glue::glue("{resub.plots.png}FS11_NIL_PC_length_brood_correlation.png"), 
       height = 6, 
       width = 10, dpi = 300)

########################################################################################################################
# LINKAGE MAPPING - PLOTS - END
########################################################################################################################

########################################################################################################################
# GWA MAPPING - PLOT TRAIT CORRELATIONS - START
########################################################################################################################

gwas_traits <- data.table::fread(glue::glue("{resub.data}Supplemental_Data15_GWAS_ALL_TRAITS.tsv"))
gwas_loadings <- data.table::fread(glue::glue("{resub.data}Supplemental_Data16_GWAS_PC_Loadings.tsv"))
gwas_trait_cor <- data.table::fread(glue::glue("{resub.data}Supplemental_Data18_GWAS_Trait_Correlations.tsv"))

trait_order <- gwas_loadings %>%
  dplyr::filter(PC == "PC1") %>%
  dplyr::arrange(desc(abs(Loading)))

gwas_cor_plot <- gwas_trait_cor %>%
  dplyr::filter(!grepl("PC", trait_b),
                !grepl("PC", trait_a)) %>%
  ggplot() +
  aes(x = factor(trait_a,
                 levels = trait_order$Trait), 
      y = factor(trait_b,
                 levels = trait_order$Trait), 
      fill = trait_cor)+
  geom_tile()+
  scale_fill_gradient2(low = "blue", 
                       high = "yellow", 
                       mid = "purple", 
                       midpoint = 0, 
                       limit = c(-1,1), 
                       space = "Lab", 
                       name="Pearson\nCorrelation") +
  arsenic.theme +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())


gwas_pc_loading_plot <- gwas_loadings%>%
  dplyr::filter(PC %in% c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")) %>%
  ggplot() + 
  aes(x = factor(Trait,
                 levels = trait_order$Trait), 
      y = factor(PC, levels = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")), 
      fill = Loading)+
  geom_tile()+
  scale_fill_gradient2(low = "blue",
                       high = "yellow",
                       mid = "purple",
                       space = "Lab",
                       name="Loading") +
  arsenic.theme +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())


cowplot::plot_grid(gwas_cor_plot, 
                   gwas_pc_loading_plot, 
                   labels = "AUTO", 
                   ncol = 1, 
                   align = 'v',
                   vjust = 1, 
                   label_size = 18)

ggsave(glue::glue("{resub.plots.png}FS12_GWAS_Trait-correlation_loadings.png"), 
       height = 12, 
       width = 20, dpi = 300)

########################################################################################################################
# GWA MAPPING - PLOT TRAIT CORRELATIONS - END
########################################################################################################################


########################################################################################################################
# PROCESS GWA MAPPING PHENOTYPES - GENOMIC HERITABILITY - DR H2 - EFFECT SIZE - PLOTS - START
########################################################################################################################

gwas_h2 <- data.table::fread(glue::glue("{resub.data}Supplemental_Data20_GWAS_Genomic_H2_Estimates.tsv"))
h2_df <- data.table::fread(glue::glue("{resub.data}Supplemental_Data7_DR-Heritability_estimates.tsv"))

pt_df <- h2_df %>%
  left_join(., gwas_h2, by = "trait") %>%
  dplyr::filter(condition!="Water")%>%
  dplyr::mutate(cpoint = ifelse(trait == "PC1", 1,
                                ifelse(trait == "norm.n",3,
                                       ifelse(trait == "median.TOF", 2, 4))),
                tidy_cond = factor(condition, 
                                   levels = c( "arsenic1", "arsenic2", "arsenic3","arsenic4"),
                                   labels = c( "250 µM", "500 µM", "1000 µM","2000 µM")))

ggplot(pt_df)+
  aes(y = narrow, x = Broad_H2_add_epp_model)+
  geom_point(color = "gray50", 
             data = dplyr::filter(pt_df, !trait %in% c("PC1","norm.n","median.TOF"),
                                  condition == "arsenic3")) +
  geom_point(aes(color = factor(cpoint, labels=c("PC1","Animal Length","Brood Size")),
                 size = factor(cpoint, labels=c("PC1","Animal Length","Brood Size")),
                 fill = factor(cpoint, labels=c("PC1","Animal Length","Brood Size")),
                 shape = factor(cpoint, labels=c("PC1","Animal Length","Brood Size"))),
             data = dplyr::filter(pt_df, trait %in% c("PC1","norm.n","median.TOF"),
                                  condition == "arsenic3"))+
  scale_shape_manual(values=c(21,21,21,1))+
  scale_color_manual(values=c("black","black","black","gray50"))+
  scale_fill_manual(values=c("red","hotpink3","cadetblue3","gray50"))+
  scale_size_manual(values=c(3,3,3,1))+
  labs(y = expression(italic(h)^2), x = expression(italic(H)^2),
       color="Trait",
       size = "Trait",
       fill="Trait",
       shape="Trait")+
  theme_bw(18)+
  ylim(c(0,1))+
  scale_x_continuous(breaks=c(0,0.5,1), limits = c(0,1))+
  geom_abline(slope=1,intercept = 0, linetype = 2, alpha = 0.5)+
  theme(strip.background = element_blank(),
        legend.position = "top",
        panel.grid.minor = element_blank())

ggsave(glue::glue("{resub.plots.png}FS13_GWAS-H2_v_h2.png"), 
       height = 8, 
       width = 10, dpi = 300)

########################################################################################################################
# PROCESS GWA MAPPING PHENOTYPES - GENOMIC HERITABILITY - DR H2 - EFFECT SIZE - PLOTS - END
########################################################################################################################

########################################################################################################################
# GWA MAPPING - PLOT FIGURE 2 - START
########################################################################################################################

gwas_mapping <- data.table::fread(glue::glue("{resub.data}Supplemental_Data21_PC1_processed_mapping.tsv"))
independent_tests <- data.table::fread(glue::glue("{resub.data}Supplemental_Data23_total_independent_tests.txt"))
gwas_fine_mapping <- data.table::fread(glue::glue("{resub.data}Supplemental_Data22_PC1_snpeff_genes.tsv")) 
geno_matrix <- readr::read_tsv(glue::glue("{resub.data}Supplemental_Data19_Genotype_Matrix.tsv"))%>%
  na.omit()

qtl_efffect <- na.omit(gwas_mapping)

Za <- diag(nrow(qtl_efffect))
Ze <- diag(nrow(qtl_efffect))

A <- sommer::A.mat(t(geno_matrix[,qtl_efffect$strain]))
E <- sommer::E.mat(t(geno_matrix[,qtl_efffect$strain]))

# Find variance explained by QTL
gwas_mod <- sommer::mmer(Y=qtl_efffect$value,
                         Z=list(add = list(Z=Za,K=A), 
                                epp = list(Z=Ze,K=E)),
                         IMP = T, init.equal = F)

suma <- summary(gwas_mod)$var.comp.table

H2_no_qtl <- sum(suma[1:2,1])/sum(suma[,1])
h2_no_qtl <- sum(suma[1,1])/sum(suma[,1])

gwas_mod <- sommer::mmer(Y=qtl_efffect$value, 
                         X=cbind(1,qtl_efffect$allele),
                         Z=list(add = list(Z=Za,K=A), 
                                epp = list(Z=Ze,K=E)),
                         IMP = T, init.equal = F)

suma <- summary(gwas_mod)$var.comp.table

H2_qtl <- sum(suma[1:2,1])/sum(suma[,1])
h2_qtl <- sum(suma[1,1])/sum(suma[,1])

H2_no_qtl-H2_qtl
# 0.22

base_manplot <- cegwas2_manplot(gwas_mapping, eigen_cutoff = -log10(0.05/independent_tests$V1))

pc1_manplot <- base_manplot[[1]]+
  theme(plot.title = element_blank())

peak_pos <- na.omit(gwas_mapping) %>%
  dplyr::mutate(facet_marker = paste0(CHROM, ":", peakPOS)) %>%
  dplyr::pull(facet_marker) %>%
  unique()

pxg_split_plot <- na.omit(gwas_mapping) %>%
  dplyr::mutate(facet_marker = paste0(CHROM, ":", peakPOS)) %>%
  dplyr::group_by(allele, facet_marker)%>%
  dplyr::mutate(mean_pheno = mean(as.numeric(value), na.rm = T))%>%
  dplyr::mutate(n2_cb = case_when(
    strain == "N2" ~ "1",
    strain == "CB4856" ~ "2", 
    TRUE ~ "3"
  )) %>%
  ggplot()+
  aes(x = factor(allele, levels = c(-1,1), labels = c("REF","ALT")))+
  geom_beeswarm(cex=2,priority='density',
                aes(y = as.numeric(value),
                    fill = n2_cb,
                    size = n2_cb),
                shape = 21)+
  geom_boxplot(aes(y=as.numeric(value)), alpha = 0.5, fill = "gray70") +
  scale_fill_manual(values=c("orange","blue", "gray50"))+
  scale_size_manual(values=c(4,4,2))+
  theme_bw(15)+
  labs(y = "PC1",
       x = glue::glue("Genotype at {peak_pos}")) +
  arsenic.theme +
  theme(strip.background.x = element_blank())+
  theme(legend.position = "none")


snpeff_fine <- gwas_fine_mapping %>%
  dplyr::select(MARKER, POS, STRAIN, REF,ALT, TGT = STRAIN_GENOTYPE, VARIANT_IMPACT,
                VARIANT_LD_WITH_PEAK_MARKER, PEAK_MARKER, QTL_INTERVAL_START,
                QTL_INTERVAL_END, VARIANT_LOG10p)

snpeff_fine$VARIANT_IMPACT[is.na(snpeff_fine$VARIANT_IMPACT)] <- "INTERGENIC"

LD_genotypes <- snpeff_fine %>%
  dplyr::filter(STRAIN == "CB4856", TGT == ALT)  %>%
  dplyr::mutate(tidy_marker = gsub("_",":",MARKER))

peak_roi_marker <- LD_genotypes %>%
  dplyr::filter(tidy_marker == PEAK_MARKER)

finemap_plot <- LD_genotypes%>%
  na.omit() %>%
  ggplot() +
  aes(x = POS/1e6) +
  geom_vline(aes(xintercept = 7931252/1e6), 
             color = "red",
             linetype = 2) +
  geom_vline(aes(xintercept = 7.83), color = "gray60")+
  geom_vline(aes(xintercept = 8.02), color = "gray60")+
  geom_point(aes(fill = factor(VARIANT_IMPACT,
                               levels = rev(c("INTERGENIC", "MODIFIER", "LOW", "MODERATE", "HIGH"))), 
                 y = VARIANT_LOG10p), 
             size = 2,
             shape = 21)+
  scale_fill_viridis_d(name = "Variant\nImpact", direction = -1) +
  theme_bw(15)+
  arsenic.theme + 
  labs(x = "Genomic Position (Mb)",
       y = expression(-log[10](italic(p))))+
  theme(strip.background.x  = element_blank())

bottom <- plot_grid(pxg_split_plot, 
                    finemap_plot, 
                    labels = c("B", "C"), 
                    rel_widths = c(0.5,1),
                    align = 'v', 
                    label_size = 18, 
                    ncol = 2)

figure2 <- plot_grid(pc1_manplot, 
                     bottom, 
                     labels = c("A", ""), 
                     rel_heights = c(1,1),
                     align = 'v', 
                     label_size = 18, 
                     ncol = 1)

ggsave(plot = figure2,
       glue::glue("{resub.plots.png}F2_GWAS.png"), 
       height = 8, 
       width = 12, dpi = 300)

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
  geom_vline(aes(xintercept = 7931252/1e6), 
             color = "red",
             linetype = 2) +
  geom_vline(aes(xintercept = 7.83), color = "gray60")+
  geom_vline(aes(xintercept = 8.02), color = "gray60")+
  geom_point(aes(fill = factor(VARIANT_IMPACT,
                               levels = rev(c("INTERGENIC", "MODIFIER", "LOW", "MODERATE", "HIGH"))), 
                 y = VARIANT_LOG10p), 
             size = 2,
             shape = 21)+
  facet_grid(.~cb_alt) + 
  scale_fill_viridis_d(name = "Variant\nImpact", direction = -1) +
  theme_bw(15)+
  arsenic.theme + 
  labs(x = "Genomic Position (Mb)",
       y = expression(-log[10](italic(p))))+
  theme(strip.background.x  = element_blank())

ggsave(glue::glue("{resub.plots.png}FS15_GWAS_FINE_MAP.png"), 
       height = 8, 
       width = 12, dpi = 300)

########################################################################################################################
# GWA MAPPING - PLOT FIGURE 2 - END
########################################################################################################################

########################################################################################################################
# DBT-1 C78S SWAP - PLOTS and EFFECT SIZES- START
########################################################################################################################

swap_pheno <- data.table::fread(glue::glue("{resub.data}Supplemental_Data13_NIL_control_regressed.tsv"))

# calculate % effect size of swaps
test_swap <- swap_pheno%>%
  dplyr::filter(Trait == "PC1",
                Strain != "ECA591", 
                Strain != "ECA414", 
                Strain != "ECA434", 
                Strain != "ECA582", 
                Strain != "ECA589")%>%
  dplyr::mutate(strain1 = factor(Strain, 
                                 levels = c("N2","ECA581","CB4856","ECA590"),
                                 labels =c("N2\nDBT-1(C78)", 
                                           "N2\nDBT-1 (S78)",
                                           "CB4856\nDBT-1 (S78)",
                                           "CB4856\nDBT-1(C78)" ),
                                 ordered = T))


sjstats::anova_stats(car::Anova(aov(Value ~ Strain, data = dplyr::filter(test_swap, Strain %in% c("N2", "CB4856")))))
sjstats::anova_stats(car::Anova(aov(Value ~ Strain, data = dplyr::filter(test_swap, Strain %in% c("N2", "ECA581")))))
sjstats::anova_stats(car::Anova(aov(Value ~ Strain, data = dplyr::filter(test_swap, Strain %in% c("CB4856", "ECA590")))))
# n2 - 581
# 0.69
# # cb - 590
# 0.65

swap_pheno%>%
  dplyr::filter(Trait == "PC1",
                Strain != "ECA591", 
                Strain != "ECA414", 
                Strain != "ECA434", 
                Strain != "ECA582", 
                Strain != "ECA589")%>%
  dplyr::mutate(strain1 = factor(Strain, 
                                 levels = c("N2","ECA581","CB4856","ECA590"),
                                 labels =c("N2\nDBT-1(C78)", 
                                           "N2\nDBT-1 (S78)",
                                           "CB4856\nDBT-1 (S78)",
                                           "CB4856\nDBT-1(C78)" ),
                                 ordered = T))%>%
  ggplot(.) +
  aes(x = factor(strain1), 
      y = Value, 
      fill=Strain) +
  geom_beeswarm(cex=1.2,priority='density', size = 2, alpha = 0.5)+
  geom_boxplot(outlier.colour = NA, alpha = 0.75)+
  scale_fill_manual(values = c("N2" = "orange","CB4856" = "blue",
                               "ECA414" = "blue","ECA434" = "blue",
                               "ECA581" = "gray50","ECA582" = "gray50",
                               "ECA589" = "gray50","ECA590" = "gray50",
                               "ECA591" = "gray50"))+
  theme_bw()+
  arsenic.theme+
  theme(legend.position="none",
        axis.title.x = ggplot2::element_text(size = 0, face = "bold", color = "black", vjust = -0.3),
        axis.title.y = ggplot2::element_text(size = 16, face = "bold", color = "black", vjust = -0.3))+
  labs( y = paste0("PC 1"))

ggsave(glue::glue("{resub.plots.png}F3_Allele_Swap_PC1.png"), 
       height = 4, 
       width = 8, dpi = 300)


# PC trait correlation plot
nil_to_lm <- swap_pheno %>%
  dplyr::filter(Strain != "ECA591", 
                Strain != "ECA414", 
                Strain != "ECA434", 
                Strain != "ECA582", 
                Strain != "ECA589") %>%
  dplyr::filter(Trait %in% c("norm.n", "median.TOF", "PC1")) %>%
  tidyr::spread(Trait, Value) %>%
  tidyr::gather(CorTrait, value, -(Strain:u_id), -PC1) 


lm_ls <- list()
for(tra in 1:length(unique(nil_to_lm$CorTrait))) {
  nil_trait_to_lm <- nil_to_lm %>%
    dplyr::filter(CorTrait == unique(nil_to_lm$CorTrait)[tra]) 
  
  lm_ls[[tra]] <- nil_trait_to_lm %>%
    tidyr::nest(CorTrait) %>%
    mutate(model = purrr::map(data, ~ lm(PC1 ~ value, data = .x)),
           adj.r.squared = purrr::map_dbl(model, ~ signif(summary(.x)$adj.r.squared, 5)),
           intercept = purrr::map_dbl(model, ~ signif(.x$coef[[1]],5)),
           slope = purrr::map_dbl(model, ~ signif(.x$coef[[2]], 5)),
           pvalue = purrr::map_dbl(model, ~ signif(summary(.x)$coef[2,4], 5)) 
    )    %>%
    dplyr::select(-data, -model) %>% 
    dplyr::left_join(nil_to_lm)
}

nil_lm <- dplyr::bind_rows(lm_ls)

nil_lm %>%
  ggplot()+
  aes(x = value, y = PC1,color = CorTrait)+
  geom_point(aes( )) + 
  geom_smooth(method='lm',formula=y~x, se = F) +
  scale_color_manual(values = c("hotpink3", "cadetblue3"), 
                     breaks = c("median.TOF","norm.n"),
                     labels = c("Length", "Brood\nSize"))  +
  theme_bw(18) +
  labs(x = "Phenotype Value", color = "Trait", y = "PC 1") +
  arsenic.theme +
  geom_text(aes(45, 10, 
                label = paste("Adj R2 = ", adj.r.squared, "\n",
                              "Intercept =",intercept, "\n",
                              "Slope =", slope, "\n",
                              "P =", pvalue)), data = dplyr::filter(nil_lm, CorTrait == "norm.n")) +
  geom_text(aes(-35, -10, 
                label = paste("Adj R2 = ", adj.r.squared, "\n",
                              "Intercept =",intercept, "\n",
                              "Slope =", slope, "\n",
                              "P =", pvalue)), data = dplyr::filter(nil_lm, CorTrait == "median.TOF"))+
  theme(panel.grid.minor = element_blank())


ggsave(glue::glue("{resub.plots.png}FS17_Allele_Swap_PC1_Brood_Size_Correlation.png"), 
       height = 8, 
       width = 12, dpi = 300)

swap_pheno%>%
  dplyr::filter(Trait %in% c("median.TOF", "norm.n"),
                Strain != "ECA591", 
                Strain != "ECA414", 
                Strain != "ECA434", 
                Strain != "ECA582", 
                Strain != "ECA589")%>%
  dplyr::mutate(strain1 = factor(Strain, 
                                 levels = c("N2","ECA581","CB4856","ECA590"),
                                 labels =c("N2\nDBT-1(C78)", 
                                           "N2\nDBT-1 (S78)",
                                           "CB4856\nDBT-1 (S78)",
                                           "CB4856\nDBT-1(C78)" ),
                                 ordered = T),
                tidy_trait = ifelse(Trait == "norm.n", "Brood Size", "Animal Length"))%>%
  ggplot(.) +
  aes(x = factor(strain1), 
      y = Value, 
      fill=Strain) +
  geom_beeswarm(cex=1.2,priority='density', size = 2, alpha = 0.5)+
  geom_boxplot(outlier.colour = NA, alpha = 0.75)+
  scale_fill_manual(values = c("N2" = "orange","CB4856" = "blue",
                               "ECA414" = "blue","ECA434" = "blue",
                               "ECA581" = "gray50","ECA582" = "gray50",
                               "ECA589" = "gray50","ECA590" = "gray50",
                               "ECA591" = "gray50"))+
  theme_bw()+
  facet_grid(tidy_trait ~.) +
  arsenic.theme+
  theme(legend.position="none",
        axis.title.x = ggplot2::element_text(size = 0, face = "bold", color = "black", vjust = -0.3),
        axis.title.y = ggplot2::element_text(size = 16, face = "bold", color = "black", vjust = -0.3),
        strip.background = element_blank(),
        strip.text.y = element_text(size = 16))+
  labs( y = paste0("Phenotype Value"))


ggsave(glue::glue("{resub.plots.png}FS18_Allele_Swap_Brood_Length.png"), 
       height = 6, 
       width = 10, dpi = 300)

########################################################################################################################
# DBT-1 C78S SWAP - PLOTS and EFFECT SIZES- END
########################################################################################################################

########################################################################################################################
# DBT-1 C78S C15ISO RESCUE PHENOTYPES - PLOT PC1 v Brood and Length - START
########################################################################################################################

rescue_pheno <- data.table::fread( glue::glue("{resub.data}Supplemental_Data23_Rescue_control_regressed.tsv"))

# PC trait correlation plot
rescue_to_lm <- rescue_pheno %>%
  dplyr::filter(Trait %in% c("norm.n", "median.TOF", "PC1"),
                Condition %in% c("Arsenic",
                                 "C15ISO64", 
                                 "ArsenicC15ISO64")) %>%
  tidyr::spread(Trait, Value) %>%
  tidyr::gather(CorTrait, value, -(Strain:conc), -PC1) 


lm_ls <- list()
for(tra in 1:length(unique(rescue_to_lm$CorTrait))) {
  rescue_trait_to_lm <- rescue_to_lm %>%
    dplyr::filter(CorTrait == unique(rescue_to_lm$CorTrait)[tra]) %>%
    dplyr::distinct(Strain, Condition, value, .keep_all = T)
  
  lm_ls[[tra]] <- rescue_trait_to_lm %>%
    tidyr::nest(CorTrait) %>%
    mutate(model = purrr::map(data, ~ lm(PC1 ~ value, data = .x)),
           adj.r.squared = purrr::map_dbl(model, ~ signif(summary(.x)$adj.r.squared, 5)),
           intercept = purrr::map_dbl(model, ~ signif(.x$coef[[1]],5)),
           slope = purrr::map_dbl(model, ~ signif(.x$coef[[2]], 5)),
           pvalue = purrr::map_dbl(model, ~ signif(summary(.x)$coef[2,4], 5)) 
    )    %>%
    dplyr::select(-data, -model) %>% 
    dplyr::left_join(rescue_trait_to_lm)
}

rescue_lm <- dplyr::bind_rows(lm_ls)

rescue_lm %>%
  ggplot()+
  aes(x = value, 
      y = PC1,
      color = CorTrait)+
  geom_point(aes( )) + 
  geom_smooth(method='lm',formula=y~x, se = F) +
  scale_color_manual(values = c("hotpink3", "cadetblue3"), 
                     breaks = c("median.TOF","norm.n"),
                     labels = c("Length", "Brood\nSize"))  +
  theme_bw(18) +
  labs(x = "Phenotype Value", color = "Trait", y = "PC 1") +
  arsenic.theme +
  geom_text(aes(5, 10, 
                label = paste("Adj R2 = ", adj.r.squared, "\n",
                              "Intercept =",intercept, "\n",
                              "Slope =", slope, "\n",
                              "P =", pvalue)),
            data = dplyr::filter(rescue_lm, CorTrait == "norm.n")) +
  geom_text(aes(200, -10, 
                label = paste("Adj R2 = ", adj.r.squared, "\n",
                              "Intercept =",intercept, "\n",
                              "Slope =", slope, "\n",
                              "P =", pvalue)), 
            data = dplyr::filter(rescue_lm, CorTrait == "median.TOF"))+
  theme(panel.grid.minor = element_blank())


ggsave(glue::glue("{resub.plots.png}FS19_Rescue_PC1_Brood_Size_Correlation.png"), 
       height = 10, 
       width = 16, dpi = 300)

########################################################################################################################
# DBT-1 C78S C15ISO RESCUE PHENOTYPES - PLOT PC1 v Brood and Length - END
########################################################################################################################

########################################################################################################################
# DBT-1 C78S C15ISO RESCUE PHENOTYPES - PLOT TRAIT CORRELATIONS AND LOADINGS - START
########################################################################################################################

rescue_trait_cor <- data.table::fread( glue::glue("{resub.data}Supplemental_Data24_Rescue_Trait_Correlations.tsv"))
rescue_pc_loading <- data.table::fread( glue::glue("{resub.data}Supplemental_Data25_Rescue_PC_Loadings.tsv"))

trait_order <- rescue_pc_loading %>%
  dplyr::filter(PC == "PC1") %>%
  dplyr::arrange(desc(abs(loading)))

rescue_cor_plot <- rescue_trait_cor %>%
  dplyr::filter(!grepl("PC", trait_b),
                !grepl("PC", trait_a)) %>%
  ggplot() +
  aes(x = factor(trait_a,
                 levels = trait_order$Trait), 
      y = factor(trait_b,
                 levels = trait_order$Trait), 
      fill = trait_cor)+
  geom_tile()+
  scale_fill_gradient2(low = "blue", 
                       high = "yellow", 
                       mid = "purple", 
                       midpoint = 0, 
                       limit = c(-1,1), 
                       space = "Lab", 
                       name="Pearson\nCorrelation") +
  arsenic.theme +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())


rescue_pc_loading_plot <- rescue_pc_loading%>%
  dplyr::filter(PC %in% c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")) %>%
  ggplot() + 
  aes(x = factor(Trait,
                 levels = trait_order$Trait), 
      y = factor(PC, levels = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")), 
      fill = loading)+
  geom_tile()+
  scale_fill_gradient2(low = "blue",
                       high = "yellow",
                       mid = "purple",
                       space = "Lab",
                       name="Loading") +
  arsenic.theme +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())


cowplot::plot_grid(rescue_cor_plot, 
                   rescue_pc_loading_plot, 
                   labels = "AUTO", 
                   ncol = 1, 
                   align = 'v',
                   vjust = 1, 
                   label_size = 18)

ggsave(glue::glue("{resub.plots.png}FS20_Rescue_Trait-correlation_loadings.png"), 
       height = 12, 
       width = 20, dpi = 300)

########################################################################################################################
# DBT-1 C78S C15ISO RESCUE PHENOTYPES - PLOT TRAIT CORRELATIONS AND LOADINGS - END
########################################################################################################################

########################################################################################################################
# DBT-1 C78S C15ISO RESCUE PHENOTYPES - BOXPLOTS - START
########################################################################################################################

rescue_pheno_pr <- rescue_pheno %>%
  dplyr::group_by(Strain, Condition, Trait)%>%
  dplyr::mutate(mph = median(Value),
                sph = sd(Value))%>%
  dplyr::mutate(flag_h = 2*sph+mph,
                flag_l = mph-2*sph)%>%
  dplyr::mutate(cut_h =ifelse(Value >= 2*sph+mph, "YES", "NO"),
                cut_l =ifelse(Value <= mph-2*sph, "YES", "NO"))%>%
  dplyr::filter(cut_h != "YES" , cut_l !="YES")%>%
  dplyr::ungroup()%>%
  dplyr::select(-cut_l, -cut_h,-flag_l,-flag_h,-sph,-mph)

mc <- "64"
cc <- c("Arsenic",
        "C15ISO64", 
        "ArsenicC15ISO64")

boxplot_plt(df = rescue_pheno_pr,
            trt = "PC1",
            cond = cc,
            fancy_name = paste0("PC 1"),
            strains = c("N2",
                        "CB4856",
                        "ECA581",
                        "ECA590"),
            fancy_strains = c("Bristol",
                              "Hawaii",
                              "Bristol\n(C78S)",
                              "Hawaii\n(S78C)"),
            ordered_conditions = c("EtOH", 
                                   "C15ISO12",
                                   "C15ISO24", 
                                   "C15ISO48",
                                   "C15ISO64", 
                                   "C15ISO100", 
                                   'Arsenic', 
                                   "ArsenicC15ISO12",
                                   "ArsenicC15ISO24",
                                   "ArsenicC15ISO48",
                                   "ArsenicC15ISO64",
                                   "ArsenicC15ISO100"),
            fancy_ordered_conditions = c("EtOH", 
                                         "C15ISO\n12",
                                         "C15ISO\n24", 
                                         "C15ISO\n48",
                                         "C15ISO\n64", 
                                         "C15ISO\n100", 
                                         'Arsenic', 
                                         "Arsenic\nC15ISO\n12",
                                         "Arsenic\nC15ISO\n24",
                                         "Arsenic\nC15ISO\n48",
                                         "Arsenic\nC15ISO\n64",
                                         "Arsenic\nC15ISO\n100"),
            r_conc = mc)+ theme(axis.text.x = ggplot2::element_text(size = 12),
                                axis.title.y = element_text(size =16))

ggsave(glue::glue("{resub.plots.png}FS21_Rescue_PC1.png"), 
       height = 6, 
       width = 18, dpi = 300)

mc <- "64"
cc <- c("Arsenic",
        "ArsenicC15ISO64")

boxplot_plt(df = rescue_pheno_pr,
            trt = "PC1",
            cond = cc,
            fancy_name = paste0("PC 1"),
            strains = c("N2",
                        "CB4856",
                        "ECA581",
                        "ECA590"),
            fancy_strains = c("Bristol",
                              "Hawaii",
                              "Bristol\n(C78S)",
                              "Hawaii\n(S78C)"),
            ordered_conditions = c("EtOH", 
                                   "C15ISO12",
                                   "C15ISO24", 
                                   "C15ISO48",
                                   "C15ISO64", 
                                   "C15ISO100", 
                                   'Arsenic', 
                                   "ArsenicC15ISO12",
                                   "ArsenicC15ISO24",
                                   "ArsenicC15ISO48",
                                   "ArsenicC15ISO64",
                                   "ArsenicC15ISO100"),
            fancy_ordered_conditions = c("EtOH", 
                                         "C15ISO\n12",
                                         "C15ISO\n24", 
                                         "C15ISO\n48",
                                         "C15ISO\n64", 
                                         "C15ISO\n100", 
                                         'Arsenic', 
                                         "Arsenic\nC15ISO\n12",
                                         "Arsenic\nC15ISO\n24",
                                         "Arsenic\nC15ISO\n48",
                                         "Arsenic\nC15ISO\n64",
                                         "Arsenic\nC15ISO\n100"),
            r_conc = mc)+ theme(axis.text.x = ggplot2::element_text(size = 12),
                                axis.title.y = element_text(size =16))

ggsave(glue::glue("{resub.plots.png}F4C_Rescue_PC1.png"), 
       height = 6, 
       width = 12, dpi = 300)

########################################################################################################################
# DBT-1 C78S C15ISO RESCUE PHENOTYPES - BOXPLOTS - END
########################################################################################################################

########################################################################################################################
# PROCESS AND PLOT HUMAN COMPETITION EXPERIMENT - START
########################################################################################################################

human_reads <- readr::read_csv(glue::glue("{resub.dir}Human_Experiment/Human_Reads.csv"))

pr_reads <- human_reads %>%
  dplyr::filter(Arsenic_Concentration %in% c("0uM","5uM")) %>%
  na.omit() %>%
  dplyr::mutate(perc_edit = Guide1_edit/Guide1_total) %>%
  dplyr::select(Edit, Replicate, Arsenic_Concentration, perc_edit) 

cntrl <- pr_reads %>%
  dplyr::filter(Arsenic_Concentration == "0uM") %>%
  dplyr::rename(control_edit = perc_edit) %>%
  dplyr::select(-Arsenic_Concentration) %>%
  dplyr::group_by(Edit, Replicate) %>%
  dplyr::mutate(control_means = mean(control_edit))

arsenic_5 <- pr_reads %>%
  dplyr::filter(Arsenic_Concentration == "5uM") %>%
  dplyr::rename(arsenic_edit = perc_edit) %>%
  dplyr::left_join(., cntrl, by = c("Edit", "Replicate")) %>%
  dplyr::mutate(delta_cntrl = (arsenic_edit - control_means)*100) %>%
  dplyr::distinct(Edit, delta_cntrl, .keep_all=T) %>%
  dplyr::group_by(Edit) %>%
  dplyr::summarise(mean_delta = mean(delta_cntrl),
                   sd_delta = sd(delta_cntrl))

ggplot(arsenic_5)+
  aes(x = factor(Edit, levels = c("W84C","S112C","R113C")), y = mean_delta, fill = factor(Edit, levels = c("W84C","S112C","R113C")))+
  geom_bar(stat="identity", color = "black", size = 1) +
  geom_errorbar(aes(ymin=mean_delta, ymax=mean_delta+sd_delta), width=.2, size =1) +
  theme_classic(18) +
  scale_fill_manual(values = c("hotpink3", "cadetblue3","black")) +
  labs(y = "Percent Edit Enrichment in Arsenic") +
  theme(axis.title.x = element_blank(),
        legend.position = "none")

ggsave(glue::glue("{resub.plots.png}FXX_Human_Cell_Barplot.png"), 
       height = 6, 
       width = 8, dpi = 300)

########################################################################################################################
# PROCESS HUMAN COMPETITION EXPERIMENT - END
########################################################################################################################
