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
# INTRO 
########################################################################################################################

######################################################################################################################## Phenotypic distribution logo

pheno <- data.table::fread(glue::glue("{arsenic_data}Figure 1-source data 8.tsv")) %>%
  dplyr::filter(Trait == "mean.TOF") %>%
  na.omit() %>%
  dplyr::select(strain = Strain, value = Value) %>%
  dplyr::distinct()

pheno%>%
  na.omit() %>%
  dplyr::arrange(value)%>%
  dplyr::mutate(strain2 = factor(strain, levels = unique(strain), labels = unique(strain), ordered = T))%>%
  dplyr::mutate(norm_pheno_temp = ifelse(value == min(value), 0, 1))%>%
  dplyr::mutate(delta_pheno = ifelse(norm_pheno_temp == 0, 0, abs(dplyr::lag(value) - value)))%>%
  dplyr::mutate(norm_pheno = cumsum(delta_pheno)) %>%
  dplyr::mutate(final_pheno = norm_pheno/max(norm_pheno)) %>% 
  ggplot()+
  aes(x = final_pheno,stat(count))+
  geom_density(fill = "#C6C6C6") +
  base_theme +
  labs(y = "Trait")+
  scale_x_continuous(breaks = c(0.1,0.9), 
                     labels = c("Sensitive", "Resistant"),limits = c(0,1))+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 24),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(filename = "Plots/Arsenic/Distribution.pdf", height = 6, width = 10, dpi = 400)

######################################################################################################################## RIAIL GENOTYPES

df <- readr::read_tsv("Data/gt_hmm_fill.tsv") %>%
  dplyr::group_by(sample) %>%
  dplyr::filter(chrom != "MtDNA") %>%
  dplyr::mutate(gt = ifelse(gt == 1, "N2", "CB4856")) %>%
  dplyr::mutate(low_sites = ifelse(sites < 100, TRUE, FALSE)) %>%
  dplyr::filter(grepl("QX", sample)) %>%
  dplyr::mutate(num_riail = as.numeric(gsub("QX", "",sample))) %>%
  dplyr::filter(num_riail > 249)

plot_riail_geno(df) +
  theme(axis.ticks.y = element_blank(),
        plot.background = element_rect(fill = background_color),
        panel.background = element_rect(fill = background_color, colour = NA),
        text = element_text(family = axes_title_font, size = axes_text_size),
        axis.text = element_text(family = number_font,size = rel(0.8), colour = "grey30", margin = unit(0.1, "cm")))

ggsave(filename = "Plots/Arsenic/RIAIL_Genotypes.png", height = 8, width = 12, dpi = 400)


######################################################################################################################## Example Dose Response

arsenic_DR <- data.table::fread(glue::glue("{arsenic_data}Figure 1-source data 1.tsv"))
arsenic_DR <- data.table::fread(glue::glue("{arsenic_data}Figure 1-source data 1.tsv"))

dr_plt_df <- arsenic_DR%>%
  dplyr::filter(Trait %in% "mean.TOF", Strain %in% c("N2", "CB4856"), Condition %in% c("Water", "arsenic1")) %>%
  dplyr::filter(Trait %in% c("mean.TOF", "PC1")) %>%
  dplyr::mutate(flip_pc = ifelse(Trait == "PC1", -Value, Value)) %>%
  dplyr::mutate(Tidy_trait = ifelse(Trait == "mean.TOF", "Animal Size", 
                                    ifelse(Trait == "norm.n"," Brood Size",
                                           ifelse(Trait == "mean.norm.EXT", "Optical Density",
                                                  ifelse(Trait =="mean.norm.yellow", "Fluorescence","PC 1")))),
                tidy_cond = factor(Condition, 
                                   levels = c("Water", "arsenic1", "arsenic2", "arsenic3","arsenic4"),
                                   labels = c("Control", "Toxin", "500", "1000","2000"))) 


dr_plt_df%>%
  ggplot()+
  aes(x = tidy_cond, 
      y = flip_pc, 
      fill = Strain)+
  geom_beeswarm(alpha = point_alpha,
                priority = "density",
                cex = 1.2,dodge.width=.8)+
  geom_boxplot(outlier.colour = NA, alpha = boxplot_alpha)+
  scale_fill_manual(values = strain_colors, name = "Strain")+
  scale_color_manual(values = strain_colors, name = "Strain")+
  theme_bw()+
  labs(x = "Arsenic (µM)", y = "Animal Length") + 
  base_theme +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = axis_color),
        axis.title.x = element_blank())

ggsave(filename = "Plots/Arsenic/Drug_Effect_Example.png", height = 6, width = 8, dpi = 400)

########################################################################################################################
# ARSENIC
########################################################################################################################

######################################################################################################################## Dose Response

arsenic_DR <- data.table::fread(glue::glue("{arsenic_data}Figure 1-source data 1.tsv"))
arsenic_DR <- data.table::fread(glue::glue("{arsenic_data}Figure 1-source data 1.tsv"))

dr_plt_df <- arsenic_DR%>%
  dplyr::filter(Trait %in% "mean.TOF") %>%
  dplyr::filter(Trait %in% c("mean.TOF", "PC1")) %>%
  dplyr::mutate(flip_pc = ifelse(Trait == "PC1", -Value, Value)) %>%
  dplyr::mutate(Tidy_trait = ifelse(Trait == "mean.TOF", "Animal Size", 
                                    ifelse(Trait == "norm.n"," Brood Size",
                                           ifelse(Trait == "mean.norm.EXT", "Optical Density",
                                                  ifelse(Trait =="mean.norm.yellow", "Fluorescence","PC 1")))),
                tidy_cond = factor(Condition, 
                                   levels = c("Water", "arsenic1", "arsenic2", "arsenic3","arsenic4"),
                                   labels = c("0", "250", "500", "1000","2000"))) 


dr_plt_df%>%
  ggplot()+
  aes(x = tidy_cond, 
      y = flip_pc, 
      fill = Strain)+
  geom_boxplot(outlier.colour = NA)+
  scale_fill_manual(values = strain_colors, name = "Strain")+
  theme_bw()+
  labs(x = "Arsenic (µM)", y = "Animal Length") + 
  base_theme +
  theme(panel.grid.major = element_blank(),
        axis.line = element_line(colour = axis_color),
        axis.title.x = element_blank())

ggsave(filename = "Plots/Arsenic/Arsenic_Animal_Length_DR.png", height = 6, width = 12, dpi = 400)
ggsave(filename = "Plots/Arsenic/Arsenic_Animal_Length_DR.pdf", height = 6, width = 12, dpi = 400)

######################################################################################################################## Linkage LOD PLOT

arsenic_linkage <- data.table::fread(glue::glue("{arsenic_data}Figure 1-source data 11.tsv"))

arsenic_linkage %>%
  dplyr::filter(trait == ".PC1") %>%
  maxlodplot_edit() +
  base_theme +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.border = element_rect(fill = NA))

ggsave(filename = "Plots/Arsenic/Arsenic_PC1_Linkage.png", height = 4, width = 12, dpi = 400)
ggsave(filename = "Plots/Arsenic/Arsenic_PC1_Linkage.pdf", height = 4, width = 12, dpi = 400)

######################################################################################################################## Linkage LOD PLOT

arsenic_linkage_pheno <- data.table::fread(glue::glue("{arsenic_data}Figure 1-source data 8.tsv"))

arsenic_linkage_pheno <- arsenic_linkage_pheno %>%
  dplyr::rename(strain = Strain, phenotype = Value, trait = Trait) %>%
  dplyr::select(-Condition)

data("N2xCB4856cross")

blankcross <- N2xCB4856cross

# Plot Linkage Mapping - PxG
arsenic_cross <- linkagemapping::mergepheno(blankcross, arsenic_linkage_pheno, set = 2)

pxgplot_edit(arsenic_cross, dplyr::filter(arsenic_linkage, trait == ".PC1"))+ 
  scale_x_discrete(breaks=c("N2", "CB4856"),labels=c("N2", "CB4856"))+
  labs(y = "Animal Length", x = "RIAIL Genotype at Peak Marker") +
  base_theme +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = axis_color)) +
  labs(title = "")

ggsave(filename = "Plots/Arsenic/Arsenic_PC1_Linkage_PxG.png", height = 6, width = 12, dpi = 400)
ggsave(filename = "Plots/Arsenic/Arsenic_PC1_Linkage_PxG.pdf", height = 6, width = 12, dpi = 400)

riail_bar_pheno <- rial_bar_plot(arsenic_cross,
                                 dplyr::filter(arsenic_linkage, trait == ".PC1"),
                                 color_by_genotype = F )

riail_bar_pheno[[2]] +
  base_theme +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = axis_color),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

ggsave(filename = "Plots/Arsenic/Arsenic_PC1_Linkage_BAR_gray.png", height = 6, width = 12, dpi = 400)
ggsave(filename = "Plots/Arsenic/Arsenic_PC1_Linkage_BAR_gray.pdf", height = 6, width = 12, dpi = 400)

riail_bar_pheno <- rial_bar_plot(arsenic_cross,
                                 dplyr::filter(arsenic_linkage, trait == ".PC1"),
                                 color_by_genotype = T )

riail_bar_pheno[[1]] +
  base_theme +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = axis_color),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

ggsave(filename = "Plots/Arsenic/Arsenic_PC1_Linkage_BAR_GENO_Color.png", height = 6, width = 12, dpi = 400)
ggsave(filename = "Plots/Arsenic/Arsenic_PC1_Linkage_BAR_GENO_Color.pdf", height = 6, width = 12, dpi = 400)

######################################################################################################################## PC correlation

riail_to_lm <- arsenic_linkage_pheno %>%
  dplyr::ungroup() %>%
  dplyr::filter(trait %in% c("norm.n", "mean.TOF", "PC1")) %>%
  tidyr::spread(trait, phenotype) %>%
  tidyr::gather(CorTrait, value, -strain, -PC1) 


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
    dplyr::select(-data, -model) %>% 
    left_join(riail_to_lm)
}

riail_lm <- dplyr::bind_rows(lm_ls)

riail_lm %>%
  ggplot()+
  aes(x = value, y = PC1,color = CorTrait)+
  geom_hline(yintercept = 0, color = "gray50", alpha = 0.5) +
  geom_vline(xintercept = 0, color = "gray50", alpha = 0.5) +
  geom_point(aes( )) + 
  geom_smooth(method='glm',formula=y~x, se = F) +
  scale_color_manual(values = c("hotpink3", "cadetblue3"), 
                     breaks = c("mean.TOF","norm.n"),
                     labels = c("Length", "Brood Size"))  +
  theme_bw(18) +
  labs(x = "Raw Phenotype Value", color = NULL, y = "PC 1") +
  base_theme +
  theme(panel.grid = element_blank(),
        legend.position = "top",
        legend.background = element_rect(size = 0.5, color = "black"))

ggsave(filename = "Plots/Arsenic/Arsenic_PC1_TOF-BROOD_Cor.png", height = 6, width = 8, dpi = 400)
ggsave(filename = "Plots/Arsenic/Arsenic_PC1_TOF-BROOD_Cor.pdf", height = 6, width = 8, dpi = 400)

######################################################################################################################## NILs

arsenic_nil_pheno <- data.table::fread(glue::glue("{arsenic_data}Figure 1-source data 13.tsv"))

arsenic_nil_pheno %>%
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
  geom_beeswarm(alpha = point_alpha,
                priority = "density",
                cex = 1.2)+
  geom_boxplot(outlier.colour = NA, alpha = boxplot_alpha)+
  scale_fill_manual(values = c("N2" = "#F9A227","CB4856" = "#2790F9",
                               "ECA414" = "#2790F9","ECA434" = "#2790F9"))+
  base_theme +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = axis_color),
        axis.title.x = element_blank(),
        plot.title = element_blank()) +
  labs( y = paste0("PC 1"))

ggsave(filename = "Plots/Arsenic/Arsenic_PC1_NIL.png", height = 6, width = 12, dpi = 400)
ggsave(filename = "Plots/Arsenic/Arsenic_PC1_NIL.pdf", height = 6, width = 12, dpi = 400)

######################################################################################################################## GWA MANHATTAN PLOT

arsenic_gwa <- data.table::fread(glue::glue("{arsenic_data}Figure 2-source data 4.tsv"))
independent_tests <- 500
arsenic_fine_mapping <- readr::read_tsv(glue::glue("{arsenic_data}Supplemental_Data23_PC1_snpeff_genes.tsv")) 
# geno_matrix <- readr::read_tsv(glue::glue("{arsenic_data}Figure 2-source data 8.zip"))%>%
#   na.omit()

cegwas2_manplot(plot_df = arsenic_gwa, eigen_cutoff = -log10(0.05/independent_tests),mapped_cutoff = "BF")[[1]] + 
  base_theme +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.border = element_rect(fill = NA)) +
  labs(title = "")

ggsave(filename = "Plots/Arsenic/Arsenic_PC1_GWA.png", height = 4, width = 12, dpi = 400)
ggsave(filename = "Plots/Arsenic/Arsenic_PC1_GWA.pdf", height = 4, width = 12, dpi = 400)

######################################################################################################################## GWA PxG

peak_pos <- na.omit(arsenic_gwa) %>%
  dplyr::filter(CHROM == "II") %>%
  dplyr::mutate(facet_marker = paste0(CHROM, ":", peakPOS)) %>%
  dplyr::pull(facet_marker) %>%
  unique()

pxg_df <- na.omit(arsenic_gwa) %>%
  dplyr::filter(CHROM == "II") %>%
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
  labs(y = "PC 1",
       x = glue::glue("Genotype at {peak_pos}")) +
  base_theme +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = axis_color))

ggsave(filename = "Plots/Arsenic/Arsenic_PC1_GWA_PxG.png", height = 6, width = 8, dpi = 400)
ggsave(filename = "Plots/Arsenic/Arsenic_PC1_GWA_PxG.pdf", height = 6, width = 8, dpi = 400)

arsenic_phenos <- na.omit(arsenic_gwa) %>%
  dplyr::select(strain, phenotype = value, allele) %>%
  dplyr::distinct(strain,phenotype,allele)%>%
  dplyr::mutate(clean_geno = ifelse(strain == "N2", "N2",
                                    ifelse(strain == "CB4856","CB4856",
                                           ifelse(allele == -1, "REF","ALT")))) %>%
  dplyr::arrange(phenotype) %>%
  dplyr::mutate(norm_pheno_temp = ifelse(phenotype == min(phenotype), 0, 1))%>%
  dplyr::mutate(delta_pheno = ifelse(norm_pheno_temp == 0, 0, abs(dplyr::lag(phenotype) - phenotype)))%>%
  dplyr::mutate(norm_pheno = cumsum(delta_pheno)) %>%
  dplyr::mutate(final_pheno = norm_pheno/max(norm_pheno)) %>%
  dplyr::select(phenotype = final_pheno, clean_geno) %>%
  dplyr::mutate(strain = factor(1:n()))

arsenic_phenos%>%
  ggplot()+
  aes(x=strain, y= phenotype, fill = clean_geno)+
  geom_bar(stat="identity", color = "black", size = .1) +
  labs(x = "Wild isolate", y = paste0("Arsenic sensititivity")) +
  scale_fill_manual(values=strain_colors)+
  base_theme +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line = element_line(colour = axis_color))

ggsave(filename = "Plots/Arsenic/Arsenic_TOF_GWA_PHENO_BAR.png", height = 6, width = 12, dpi = 400)
ggsave(filename = "Plots/Arsenic/Arsenic_TOF_GWA_PHENO_BAR.pdf", height = 6, width = 12, dpi = 400)

######################################################################################################################## GWA PEAK LD
gm <- readr::read_tsv(glue::glue("{arsenic_data}Figure 2-source data 5.tsv"))
arsenic_gwa <- data.table::fread(glue::glue("{arsenic_data}Figure 2-source data 4.tsv"))

LD_output <- Plot_Peak_LD(arsenic_gwa, gm)

LD_output[[1]] + 
  base_theme + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())

ggsave(filename = "Plots/Arsenic/Arsenic_PC1_GWA_PeakLD.png", height = 8, width = 12, dpi = 400)
ggsave(filename = "Plots/Arsenic/Arsenic_PC1_GWA_PeakLD.pdf", height = 8, width = 12, dpi = 400)

######################################################################################################################## GWA Fine mapping

arsenic_fine_mapping <-  readr::read_tsv(glue::glue("{arsenic_data}Supplemental_Data23_PC1_snpeff_genes.tsv")) 

snpeff_fine <- arsenic_fine_mapping %>%
  dplyr::filter(CHROM == "II") %>%
  dplyr::select(MARKER, POS, STRAIN, REF,ALT, TGT = STRAIN_GENOTYPE, VARIANT_IMPACT,
                VARIANT_LD_WITH_PEAK_MARKER, PEAK_MARKER, QTL_INTERVAL_START,QTL_INTERVAL_END, VARIANT_LOG10p)

snpeff_fine$VARIANT_IMPACT[is.na(snpeff_fine$VARIANT_IMPACT)] <- "INTERGENIC"

LD_genotypes <- snpeff_fine %>%
  dplyr::filter(STRAIN == "CB4856") %>%
  dplyr::mutate(cb_alt = ifelse(REF == TGT, "N2", "CB4856")) %>%
  dplyr::mutate(tidy_marker = gsub("_",":",MARKER))

peak_roi_marker <- LD_genotypes %>%
  dplyr::filter(tidy_marker == PEAK_MARKER)

LD_genotypes%>%
  # dplyr::filter(cb_alt == "CB4856") %>%
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
             size = 4,
             shape = 21)+
  scale_fill_viridis_d(name = "Variant\nImpact", direction = -1) +
  theme_bw(15)+
  base_theme + 
  labs(x = "Genomic Position (Mb)",
       y = expression(-log[10](italic(p))))+
  theme(panel.grid.major = element_blank(),
        axis.line = element_line(colour = axis_color))+
  xlim(c(7.775,8.175))

ggsave(filename = "Plots/Arsenic/Arsenic_PC1_GWA_FineMap.png", height = 6, width = 16, dpi = 400)
ggsave(filename = "Plots/Arsenic/Arsenic_PC1_GWA_FineMap.pdf", height = 6, width = 16, dpi = 400)

LD_genotypes%>%
  dplyr::filter(cb_alt == "CB4856") %>%
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
             size = 4,
             shape = 21)+
  scale_fill_viridis_d(name = "Variant\nImpact", direction = -1) +
  theme_bw(15)+
  base_theme + 
  labs(x = "Genomic Position (Mb)",
       y = expression(-log[10](italic(p))))+
  theme(panel.grid.major = element_blank(),
        axis.line = element_line(colour = axis_color))

ggsave(filename = "Plots/Arsenic/Arsenic_PC1_CB_GWA_FineMap.png", height = 6, width = 16, dpi = 400)
ggsave(filename = "Plots/Arsenic/Arsenic_PC1_CB_GWA_FineMap.pdf", height = 6, width = 16, dpi = 400)


LD_genotypes%>%
  na.omit() %>%
  ggplot() +
  aes(x = POS/1e6) +
  geom_vline(aes(xintercept = 7931252/1e6), 
             color = "red",
             linetype = 2) +
  geom_vline(aes(xintercept = 7.83), color = "gray60")+
  geom_vline(aes(xintercept = 8.02), color = "gray60")+
  geom_point(aes(fill = cb_alt, 
                 y = VARIANT_LOG10p), 
             size = point_size,
             shape = 21)+
  scale_fill_manual(values = strain_colors, name = "CB4856 Allele Status", breaks = c("ALT", "REF")) +
  base_theme + 
  labs(x = "Genomic Position (Mb)",
       y = expression(-log[10](italic(p))))+
  theme(panel.grid.major = element_blank(),
        axis.line = element_line(colour = axis_color))+
  xlim(c(7.75,8.175))

ggsave(filename = "Plots/Arsenic/Arsenic_TOF_GWA_FineMap_N2_CB_COLORS.png", height = 6, width = 12, dpi = 400)

######################################################################################################################## SWAP


arsenic_swap <- data.table::fread(glue::glue("{arsenic_data}Figure 1-source data 13.tsv"))

arsenic_swap%>%
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
      y = -Value, 
      fill= Strain) +
  geom_beeswarm(cex=1.2,priority='density', alpha = point_alpha, size = point_size)+
  geom_boxplot(outlier.colour = NA, alpha = boxplot_alpha)+
  scale_fill_manual(values = c("N2" = "#F9A227", "CB4856" = "#2790F9",
                               "ECA581" = "gray50","ECA590" = "gray50"))+
  base_theme +  theme(legend.position = "none",
                      panel.grid.major = element_blank(),
                      axis.line = element_line(colour = axis_color),
                      axis.title.x = element_blank()) +
  labs( y = paste0("PC 1"))

ggsave(filename = "Plots/Arsenic/Arsenic_PC1_SWAP.png", height = 6, width = 10, dpi = 400)
ggsave(filename = "Plots/Arsenic/Arsenic_PC1_SWAP.pdf", height = 6, width = 10, dpi = 400)

######################################################################################################################## Metabolites ratio in arsenic


metabolites <- data.table::fread(glue::glue("{arsenic_data}Supplemental_data26_Processed_Metabolite_Measurements.tsv"))

# second experiment to get more replicates of ECA581
L1_fa <- readr::read_csv(glue::glue("{arsenic_data}Supplemental_Data27_L1_FA_N2_ECA581.csv")) 

controls_new <- L1_fa %>% 
  tidyr::gather(FA,value,-Strain,-Condition,-Replicate) %>%
  dplyr::filter(Condition == "Water")%>%
  dplyr::select(-Condition) %>%
  dplyr::rename(control_value = value)

arsenic100_new <- L1_fa %>% 
  tidyr::gather(FA,value,-Strain,-Condition,-Replicate) %>%
  dplyr::filter(Condition == "Arsenic")%>%
  dplyr::left_join(.,controls_new, by = c("Strain", "FA", "Replicate")) %>%
  dplyr::rowwise()%>%
  dplyr::mutate(delta_control = value - control_value) %>%
  dplyr::filter(FA %in% c("17_ratio", "15_ratio"))%>%
  dplyr::select(Strain, FA, control_value, delta_control)

arsenic100_new$FA <- gsub("17_ratio", "C17iso/C17n", arsenic100_new$FA)
arsenic100_new$FA <- gsub("15_ratio", "C15iso/C15n", arsenic100_new$FA)

rC15 <- metabolites%>%
  dplyr::filter(compound %in% c("rC15", "rC17"), concentration != "200") %>%
  dplyr::filter(strain %in% c("590", "CB"))

rC15$compound <- gsub("rC15", "C15iso/C15n", rC15$compound)
rC15$compound <- gsub("rC17", "C17iso/C17n", rC15$compound)
rC15$concentration <- gsub("Mock", "Water", rC15$concentration)
rC15$concentration <- gsub("100", "Arsenic", rC15$concentration)

controls <- dplyr::filter(rC15, concentration == "Water")%>%
  dplyr::group_by(strain, compound, replicate)%>%
  dplyr::select(-concentration) %>%
  dplyr::rename(control_value = value)

arsenic100 <- dplyr::filter(rC15, concentration == "Arsenic")%>%
  dplyr::left_join(.,controls, by = c("strain", "compound","replicate"))%>%
  dplyr::rowwise()%>%
  dplyr::mutate(delta_control = value - control_value) %>%
  dplyr::select(Strain = strain, FA = compound, control_value, delta_control) 

complete_arsenic <- dplyr::bind_rows(arsenic100,arsenic100_new) %>%
  dplyr::group_by(FA) %>%
  dplyr::mutate(scale_delta = scale(delta_control)) %>%
  dplyr::mutate(expt = ifelse(Strain %in% c("N2","ECA581"), "N2 Background", "CB4856 Background")) %>%
  dplyr::group_by(Strain, FA,expt) %>%
  dplyr::mutate(mean_scale_value = mean(scale_delta),
                sd_scale_value = sd(scale_delta),
                mean_value = mean(delta_control),
                sd_value = sd(delta_control),
                tidy_strain = factor(Strain, 
                                     levels = c("N2","ECA581","CB", "590"), 
                                     labels = c("N2\nDBT-1(C78)", "N2\nDBT-1(S78)","CB4856\nDBT-1(S78)", "CB4856\nDBT-1(C78)")))

complete_arsenic %>%
  ggplot()+
  aes(x = tidy_strain)+
  geom_bar(aes(y = mean_value, fill = tidy_strain), 
           color = "black", 
           stat = "identity", 
           position = position_dodge())+
  geom_errorbar(aes(ymin=mean_value-abs(sd_value), 
                    ymax=mean_value+abs(sd_value)),
                width=.2)+
  facet_wrap(expt~FA, scales = "free")+
  scale_fill_manual(values=c("#F9A227","gray50","#2790F9","gray50"))+
  theme_bw(16)+
  base_theme +  theme(legend.position = "none",
                      panel.grid.major = element_blank(),
                      axis.line = element_line(colour = axis_color),
                      axis.title.x = element_blank()) +
  labs(y = "Arsenic - Control")

ggsave(filename = "Plots/Arsenic/Arsenic_PC1_ISOratio_Arsenic.png", height = 10, width = 10, dpi = 400)
ggsave(filename = "Plots/Arsenic/Arsenic_PC1_ISOratio_Arsenic.pdf", height = 10, width = 10, dpi = 400)

######################################################################################################################## Metabolites straight in arsenic
controls_new <- L1_fa %>% 
  tidyr::gather(FA,value,-Strain,-Condition,-Replicate) %>%
  dplyr::filter(Condition == "Water")%>%
  dplyr::select(-Condition) %>%
  dplyr::rename(control_value = value)

arsenic100_new <- L1_fa %>% 
  tidyr::gather(FA,value,-Strain,-Condition,-Replicate) %>%
  dplyr::filter(Condition == "Arsenic")%>%
  dplyr::left_join(.,controls_new, by = c("Strain", "FA", "Replicate")) %>%
  dplyr::rowwise()%>%
  dplyr::mutate(delta_control = value - control_value) %>%
  dplyr::filter(FA %in% c("C17n", "C15n"))%>%
  dplyr::select(Strain, FA, control_value, delta_control)

rC15 <- metabolites%>%
  dplyr::filter(compound %in% c("C17n", "C15n"), concentration != "200") %>%
  dplyr::filter(strain %in% c("590", "CB"))

rC15$concentration <- gsub("Mock", "Water", rC15$concentration)
rC15$concentration <- gsub("100", "Arsenic", rC15$concentration)

controls <- dplyr::filter(rC15, concentration == "Water")%>%
  dplyr::group_by(strain, compound, replicate)%>%
  dplyr::select(-concentration) %>%
  dplyr::rename(control_value = value)

arsenic100 <- dplyr::filter(rC15, concentration == "Arsenic")%>%
  dplyr::left_join(.,controls, by = c("strain", "compound","replicate"))%>%
  dplyr::rowwise()%>%
  dplyr::mutate(delta_control = value - control_value) %>%
  dplyr::select(Strain = strain, FA = compound, control_value, delta_control) 

complete_arsenic <- dplyr::bind_rows(arsenic100,arsenic100_new) %>%
  dplyr::group_by(FA) %>%
  dplyr::mutate(scale_delta = scale(delta_control)) %>%
  dplyr::mutate(expt = ifelse(Strain %in% c("N2","ECA581"), "N2 Background", "CB4856 Background")) %>%
  dplyr::group_by(Strain, FA,expt) %>%
  dplyr::mutate(mean_scale_value = mean(scale_delta),
                sd_scale_value = sd(scale_delta),
                mean_value = mean(delta_control),
                sd_value = sd(delta_control),
                tidy_strain = factor(Strain, 
                                     levels = c("N2","ECA581","CB", "590"), 
                                     labels = c("N2\nDBT-1(C78)", "N2\nDBT-1(S78)","CB4856\nDBT-1(S78)", "CB4856\nDBT-1(C78)")))

complete_arsenic %>%
  ggplot()+
  aes(x = tidy_strain)+
  geom_bar(aes(y = mean_value, fill = tidy_strain), 
           color = "black", 
           stat = "identity", 
           position = position_dodge())+
  geom_errorbar(aes(ymin=mean_value-abs(sd_value), ymax=mean_value+abs(sd_value)),
                width=.2)+
  facet_wrap(expt~FA, scales = "free")+
  scale_fill_manual(values=c("#F9A227","gray50","#2790F9","gray50"))+
  theme_bw(16)+
  base_theme +  theme(legend.position = "none",
                      panel.grid.major = element_blank(),
                      axis.line = element_line(colour = axis_color),
                      axis.title.x = element_blank()) +
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        strip.text.y = element_text(face ="bold"))+
  labs(y = "Arsenic - Control")

ggsave(filename = "Plots/Arsenic/Arsenic_PC1_StraightChain_Arsenic.png", height = 10, width = 10, dpi = 400)
ggsave(filename = "Plots/Arsenic/Arsenic_PC1_StraightChain_Arsenic.pdf", height = 10, width = 10, dpi = 400)

######################################################################################################################## Metabolites ratios in control

controls_new <- L1_fa %>% 
  tidyr::gather(FA,value,-Strain,-Condition,-Replicate) %>%
  dplyr::filter(Condition == "Water")%>%
  dplyr::select(-Condition) %>%
  dplyr::rename(control_value = value) %>%
  dplyr::filter(FA %in% c("15_ratio", "17_ratio"))%>%
  dplyr::select(Strain, FA, control_value)

controls_new$FA <- gsub("17_ratio", "C17iso/C17n", controls_new$FA)
controls_new$FA <- gsub("15_ratio", "C15iso/C15n", controls_new$FA)

rC15 <- metabolites%>%
  dplyr::filter(compound %in% c("rC17", "rC15"), concentration == "Mock") %>%
  dplyr::filter(strain %in% c("590", "CB"))

rC15$compound <- gsub("rC15", "C15iso/C15n", rC15$compound)
rC15$compound <- gsub("rC17", "C17iso/C17n", rC15$compound)
rC15$concentration <- gsub("Mock", "Water", rC15$concentration)

controls <- dplyr::filter(rC15, concentration == "Water")%>%
  dplyr::group_by(strain, compound, replicate)%>%
  dplyr::select(-concentration) %>%
  dplyr::rename(control_value = value)%>%
  dplyr::ungroup()%>%
  dplyr::select(Strain = strain, FA = compound, control_value)


complete_arsenic <- dplyr::bind_rows(controls,controls_new) %>%
  dplyr::mutate(expt = ifelse(Strain %in% c("N2","ECA581"), "N2 Background", "CB4856 Background")) %>%
  dplyr::group_by(Strain, FA,expt) %>%
  dplyr::mutate(mean_value = mean(control_value),
                sd_value = sd(control_value),
                tidy_strain = factor(Strain, 
                                     levels = c("N2","ECA581","CB", "590"), 
                                     labels = c("N2\nDBT-1(C78)", "N2\nDBT-1(S78)","CB4856\nDBT-1(S78)", "CB4856\nDBT-1(C78)")))

complete_arsenic %>%
  ggplot()+
  aes(x = tidy_strain)+
  geom_bar(aes(y = mean_value, fill = tidy_strain), 
           color = "black", 
           stat = "identity", 
           position = position_dodge())+
  geom_errorbar(aes(ymin=mean_value-abs(sd_value), ymax=mean_value+abs(sd_value)),
                width=.2)+
  facet_wrap(expt~FA, scales = "free")+
  scale_fill_manual(values=c("orange","gray50","blue","gray50"))+
  scale_fill_manual(values=c("#F9A227","gray50","#2790F9","gray50"))+
  theme_bw(16)+
  base_theme +  theme(legend.position = "none",
                      panel.grid.major = element_blank(),
                      axis.line = element_line(colour = axis_color),
                      axis.title.x = element_blank()) +
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        strip.text.y = element_text(face ="bold"))+
  labs(y = "Metabolite Levels")

ggsave(filename = "Plots/Arsenic/Arsenic_PC1_ISO_Control.png", height = 10, width = 10, dpi = 400)
ggsave(filename = "Plots/Arsenic/Arsenic_PC1_ISO_Control.pdf", height = 10, width = 10, dpi = 400)

######################################################################################################################## RESCUE
arsenic_rescue <- data.table::fread(glue::glue("{arsenic_data}Figure 4-source data 5.tsv"))

rescue_pheno_pr <- arsenic_rescue %>%
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
        "ArsenicC15ISO64")

boxplot_plt(df = rescue_pheno_pr,
            trt = "PC1",
            cond = cc,
            fancy_name = "PC 1",
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
            r_conc = mc)+ 
  base_theme +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = axis_color),
        axis.title.x = element_blank()) 

ggsave(filename = "Plots/Arsenic/Arsenic_PC1_Rescue.png", height = 8, width = 14, dpi = 400)
ggsave(filename = "Plots/Arsenic/Arsenic_PC1_Rescue.pdf", height = 8, width = 14, dpi = 400)

######################################################################################################################## Human data

human_reads <- readr::read_csv(glue::glue("{arsenic_data}Figure 5-source data 1.csv"))

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
  base_theme+
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = axis_color),
        axis.title.x = element_blank()) +
  scale_fill_manual(values = c("hotpink3", "cadetblue3","black")) +
  labs(y = "Percent Edit Enrichment in Arsenic") +
  theme(axis.title.x = element_blank(),
        legend.position = "none")


ggsave(glue::glue("Plots/Arsenic/Arsenic_Human_Cell.png"), height = 6, width = 8)
ggsave(glue::glue("Plots/Arsenic/Arsenic_Human_Cell.pdf"), height = 6, width = 8)

######################################################################################################################## PopGene

# Rscript --vanilla Interval_Popgen.R II 7598325 8210489 Ce330_annotated.vcf.gz WS245_exons.gff Arsenic 249_samples.txt

load("Run_Popgen/Arsenic_II_7598325-8210489_Diversity_Statistics.Rda")

qtl_start <- 7598325
qtl_end <- 8210489
dbt_start <- 7942463
dbt_end <- 7945206

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
  geom_vline(aes(xintercept = dbt_start/1e6), color = "#E68FAC")+
  geom_vline(aes(xintercept = dbt_end/1e6), color = "#E68FAC")+
  theme(panel.grid.major = element_blank(),
        legend.position = "none",
        axis.line = element_line(colour = axis_color),
        axis.title.y = element_blank()) +
  labs(x = "Genomic Position (Mb)") +
  xlim(qtl_start/1e6,qtl_end/1e6)

ggsave(filename = "Plots/Arsenic/Arsenic_popgen.png", height = 8, width = 12, dpi = 400)
ggsave(filename = "Plots/Arsenic/Arsenic_popgen.pdf", height = 8, width = 12, dpi = 400)

######################################################################################################################## SUBSTRATE PLOTS

strains_330 %>%
  dplyr::group_by(substrate,TGT)%>%
  dplyr::mutate(gt_ct = n())%>%
  dplyr::filter(substrate !="None")%>%
  dplyr::group_by(substrate) %>%
  dplyr::mutate(substrate_ct = n()) %>%
  dplyr::distinct(substrate, TGT, .keep_all =T) %>%
  dplyr::mutate(perc_gt = gt_ct/substrate_ct*100) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(TGT),desc(perc_gt)) %>%
  dplyr::mutate(fac_sub = factor(substrate, levels = substrate, labels = substrate)) %>%
  ggplot() +
  aes(y = perc_gt, x = fac_sub, fill = TGT)+ 
  geom_bar(stat="identity") +
  geom_text(aes(label=gt_ct, color = TGT), 
            position=position_stack(), 
            size = 12, family = "Itim", hjust = 1) +
  scale_fill_manual(values = strain_colors) +
  scale_color_manual(values = c(background_color, axis_color)) +
  coord_flip()+
  base_theme+
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = axis_color),
        axis.title.x = element_blank())+
  scale_y_continuous(limits = c(0,101), expand = c(0, 0))+
  labs(x = "Substrate")

strains_330 %>%
  dplyr::group_by(long, lat) %>%
  dplyr::mutate(same_loc = n()) %>%
  dplyr::group_by(substrate,TGT, same_loc)%>%
  dplyr::mutate(gt_ct = n())%>%
  dplyr::filter(substrate !="None")%>%
  dplyr::group_by(substrate) %>%
  dplyr::mutate(substrate_ct = n()) %>%
  dplyr::distinct(substrate, TGT, .keep_all =T) %>%
  dplyr::mutate(perc_gt = gt_ct/substrate_ct*100) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(TGT),desc(perc_gt)) %>%
  dplyr::mutate(fac_sub = factor(substrate, levels = substrate, labels = substrate)) %>%
  dplyr::mutate(sameL = ifelse(same_loc > 1, "Same Location", "Different Location")) %>%
  ggplot() +
  aes(y = perc_gt, x = fac_sub, fill = TGT)+ 
  geom_bar(stat="identity") +
  geom_text(aes(label=gt_ct),position=position_stack(), size = 12, family = "Itim", hjust = 1) +
  scale_fill_manual(values = strain_colors) +
  facet_grid(.~sameL) +
  coord_flip()+
  base_theme+
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = axis_color),
        axis.title.x = element_blank())+
  scale_y_continuous(limits = c(0,101), expand = c(0, 0))+
  labs(x = "Substrate")

