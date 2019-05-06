# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ##
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ##
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ # Process Data and Generate Figures for Arsenic Manuscript # ~ ~ ~ ~ #
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## 
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ##

# set to location of files
main.dir <- "~/Dropbox/AndersenLab/LabFolders/Stefan/Manuscripts/Arsenic/"
data.dir <- paste0(main.dir,"Data/")
final.dir <- paste0(main.dir,"Final_tables/")
plot.dir <- paste0(main.dir,"Plots/")
script.dir <- paste0(main.dir,"Scripts/")

# load functions used in this script
source(paste0(script.dir,"processing_functions.R"))

trait_of_interest <- "median.TOF"
comparison_trait <- "q90.TOF"
condition_of_interest <- "arsenictrioxide"
control_of_interest <- "1percwater"

# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ##
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ##
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ # DOSE RESPINSE EXPERIMENT # ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ #
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ##
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ##

biopruned <- data.table::fread(paste0(final.dir, "Supplemental_data1_DR_Summarized.tsv"))

# plot Dose response experiement by subtracting mean control value
# brood
pr_df <- biopruned%>%
  dplyr::ungroup()%>%
  dplyr::select(-date, -round,-experiment,-assay)%>%
  tidyr::gather(trait, value, -condition,-plate,-control,-strain,-row,-col)%>%
  dplyr::filter(trait == "norm.n", grepl("arsenic|wat",condition,ignore.case = T))%>%
  dplyr::filter(!is.na(condition))

control_means <- pr_df%>%
  dplyr::filter(condition == "Water")%>%
  dplyr::select(strain, value)%>%
  dplyr::group_by(strain)%>%
  dplyr::summarise(m_ctrl = mean(value, na.rm = T))

subtract_controls <- pr_df%>%
  dplyr::left_join(.,control_means, by = "strain")%>%
  dplyr::mutate(subtracted_pheno = value - m_ctrl)


brood_DR <- subtract_controls%>%
  ggplot()+
  aes(x = factor(condition, 
                 levels = c("Water", "arsenic1", "arsenic2", "arsenic3","arsenic4"),
                 labels = c("0", "250", "500", "1000","2000")), 
      y = subtracted_pheno, 
      fill = strain)+
  geom_boxplot(outlier.colour = NA)+
  scale_fill_manual(values = c("blue", "cadetblue3", "hotpink3", "orange"), name = "Strain")+
  theme_bw()+
  labs(x = "Arsenic (µM)",
       y = "Brood Size") + arsenic.theme + theme(axis.title.x = element_blank(),
                                                 axis.text.x = element_blank(),
                                                 legend.position = "top",
                                                 legend.text=element_text(size=14),
                                                 legend.title = element_text(size = 16))

# length
pr_df <- biopruned%>%
  dplyr::ungroup()%>%
  dplyr::select(-date, -round,-experiment,-assay)%>%
  tidyr::gather(trait, value, -condition,-plate,-control,-strain,-row,-col)%>%
  dplyr::filter(trait == "median.TOF", grepl("arsenic|wat",condition,ignore.case = T))%>%
  dplyr::filter(!is.na(condition))

control_means <- pr_df%>%
  dplyr::filter(condition == "Water")%>%
  dplyr::select(strain, value)%>%
  dplyr::group_by(strain)%>%
  dplyr::summarise(m_ctrl = mean(value, na.rm = T))

subtract_controls <- pr_df%>%
  dplyr::left_join(.,control_means, by = "strain")%>%
  dplyr::mutate(subtracted_pheno = value - m_ctrl)

TOF_DR <- subtract_controls%>%
  ggplot()+
  aes(x = factor(condition, 
                 levels = c("Water", "arsenic1", "arsenic2", "arsenic3","arsenic4"),
                 labels = c("0", "250", "500", "1000","2000")), 
      y = subtracted_pheno, 
      fill = strain)+
  geom_boxplot(outlier.colour = NA)+
  # geom_jitter()+
  scale_fill_manual(values = c("blue", "cadetblue3", "hotpink3", "orange"), name = "Strain")+
  theme_bw()+
  labs(x = "Arsenic (µM)",
       y = "Animal Length") + arsenic.theme + theme(legend.position = "none")


cowplot::plot_grid(brood_DR, TOF_DR, labels = "AUTO", ncol = 1, align = 'v', label_size = 18)

ggsave(paste0(plot.dir,paste0("FS1_DR.pdf")), 
       height = 8, 
       width = 10)


# ~ Conclusion --> 1000µM is the first concentration with a drop in all strain brood size phenotypes. 
# Use this concentration for everything.

# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ##
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## 
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ # PLOT Linkage mapping and NIL figure # ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ # 
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ##
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ##

# regenerate cross object for pxg plot
assayregressed <- data.table::fread(paste0(final.dir,"Supplemental_data3_LINKAGE_assayRegressed.tsv"))
#Prune based on bins
RIAILs2BAMFpruned <- bamf_prune(assayregressed, drop=T)
#Regress out control strains
RIAILs2regressed <- regress(RIAILs2BAMFpruned)
arsenic_linkage <- RIAILs2regressed %>%
  dplyr::filter(condition == condition_of_interest, trait == trait_of_interest | trait == "norm.n")
data("N2xCB4856cross")
blankcross <- N2xCB4856cross
arsenic_cross <- linkagemapping::mergepheno(blankcross, arsenic_linkage, set = 2)

annotatedmap <- data.table::fread(file = paste0(final.dir,"Supplemental_data5_LINKAGE_annotatedLODs_with_brood.tsv"))


pxg_linkage_TOF <- pxgplot_edit(arsenic_cross, dplyr::filter(annotatedmap, 
                                                             trait == paste0("arsenictrioxide.",trait_of_interest)))+ 
  scale_x_discrete(breaks=c("N2", "CB4856"),labels=c("N2", "CB4856"))+
  labs(y = paste0("Animal Length"))

pxg_linkage_brood <- pxgplot_edit(arsenic_cross, dplyr::filter(annotatedmap, 
                                                               trait == paste0("arsenictrioxide.","norm.n")))+ 
  scale_x_discrete(breaks=c("N2", "CB4856"),labels=c("N2", "CB4856"))+
  labs(y = paste0("Brood Size"))

plot_grid(pxg_linkage_brood, pxg_linkage_TOF, labels = "AUTO", ncol = 2, align = 'v', label_size = 18)

ggsave(paste0(plot.dir,paste0("FS2_LINKAGE_PxG.pdf")), 
       height = 4, 
       width = 14)

LOD_plot <- maxlodplot_edit2(annotatedmap)+theme(axis.title.x = element_text(size=16),
                                                 axis.title.y = element_text(size=16))

# load nil data for plot
controlregressed <- data.table::fread(file = paste0(final.dir,"TS8_NILs_SWAPs_controlRegressed.tsv"))

NIL_TOF<-controlregressed%>%
  dplyr::filter(trait == trait_of_interest,
                strain%in%c("N2","CB4856","ECA414","ECA434"))%>%
  dplyr::mutate(strain1 = factor(strain, levels = c("N2","CB4856","ECA414","ECA434","ECA581",
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
      y = phenotype, 
      fill=strain) +
  geom_jitter(width = 0.25, alpha = .4)+
  geom_boxplot(outlier.colour = NA, alpha = 0.7)+
  scale_fill_manual(values = c("N2" = "orange","CB4856" = "blue",
                               "ECA414" = "blue","ECA434" = "blue",
                               "ECA581" = "gray50","ECA582" = "gray50",
                               "ECA589" = "gray50","ECA590" = "gray50",
                               "ECA591" = "gray50"))+
  theme_bw()+
  labs( y = paste0("Animal Length"))+
  theme(axis.text.x = ggplot2::element_text(size = 12),
        axis.text.y = ggplot2::element_text(size = 12),
        axis.title.y = ggplot2::element_text(size = 16, face = "bold", color = "black", vjust = -0.3),
        strip.text.x = element_text(size = 16, face = "bold"),
        legend.position="none",
        axis.title.x = ggplot2::element_text(size = 0, face = "bold", color = "black", vjust = -0.3))


NIL_brood <- controlregressed%>%
  dplyr::filter(trait == "norm.n",
                strain%in%c("N2","CB4856","ECA414","ECA434"))%>%
  dplyr::mutate(strain1 = factor(strain, levels = c("N2","CB4856","ECA414","ECA434","ECA581",
                                                    "ECA582","ECA589","ECA590","ECA591"),
                                 labels =c("N2",
                                           "CB4856", 
                                           "CB4856>N2\n II:5.75 - 8.02Mb", 
                                           "CB4856>N2\n II:7.83 - 9.66Mb" ,
                                           "N2(C78S)\nA", 
                                           "N2(C78S)\nB" ,
                                           "CB4856(S78C)\nA",
                                           "CB4856(S78C)\nB",
                                           "CB4856(S78C)\nC")))%>%
  ggplot(.) +
  aes(x = factor(strain1), 
      y = phenotype, 
      fill=strain) +
  geom_jitter(width = 0.25, alpha = .4)+
  geom_boxplot(outlier.colour = NA, alpha = 0.7)+
  scale_fill_manual(values = c("N2" = "orange","CB4856" = "blue",
                               "ECA414" = "blue","ECA434" = "blue",
                               "ECA581" = "gray50","ECA582" = "gray50",
                               "ECA589" = "gray50","ECA590" = "gray50",
                               "ECA591" = "gray50"))+
  theme_bw()+
  theme(axis.text.x = ggplot2::element_text(size = 12),
        axis.text.y = ggplot2::element_text(size = 12),
        axis.title.y = ggplot2::element_text(size = 16, face = "bold", color = "black", vjust = -0.3),
        strip.text.x = element_text(size = 16, face = "bold"),
        legend.position="none",
        axis.title.x = ggplot2::element_text(size = 0, face = "bold", color = "black", vjust = -0.3))+
  labs( y = paste0("Brood Size"))

bottom_row <- plot_grid(NIL_brood, NIL_TOF, labels = c('B', 'C'), align = 'h', rel_widths = c(1, 1), label_size = 18)
plot_grid(LOD_plot, bottom_row, labels = c('A', ''), ncol = 1, rel_heights = c(1, 1), label_size = 18)

ggsave(paste0(plot.dir,paste0("F1_LINKAGE.pdf")), 
       height = 8, 
       width = 14)

# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## 
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## 
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ # PCA GWAS # ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ #  
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ##
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ##

pca_var_df <- data.table::fread(file = paste0(final.dir,"TS16_GWAS_PCA_VE.tsv"))

pca_plot <- ggplot(pca_var_df)+
  aes(x = PC, y = Variance)+
  geom_point()+
  theme_bw()+
  theme(legend.position = 'none')+
  labs(x="Principal Component",y = "Cumulative Variance Explained")+
  theme(legend.position = "none")+
  ylim(c(0.5,1))+
  xlim(c(0,20))+
  geom_hline(aes(yintercept = 0.95), color = "red", linetype =4, alpha = 0.7)+
  arsenic.theme

plot_grid(pca_plot, labels = c(''), ncol = 1)

ggsave(paste0(plot.dir, "FS3_PCA_variance_explained.pdf"), 
       height = 6, 
       width = 6)

# --- > 10 PCs capture >95% of the variation < --- #

# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ #
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ #
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ # PLOT PCA MAPPINGS # ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ # 
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ #
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ #

pr_pca_maps <- data.table::fread(file = paste0(final.dir,"Supplemental_data17_GWAS_PCA_processed_mappings.tsv"))

# plot PC1 mapping
pca_gwa_plot<- pr_pca_maps%>%
  dplyr::filter(CHROM!="MtDNA")%>%
  manplot_edit()

# phenotype x genotype plot for PC1
pca_pxg <- pr_pca_maps%>%
  na.omit()%>%
  dplyr::mutate(col_st = ifelse(strain == "N2", "orange",
                                ifelse(strain == "CB4856", "blue","red")))%>%
  ggplot(.)+
  aes(x=factor(as.character(allele), labels = c("REF","ALT")), y = value)+
  geom_boxplot(fill = "gray90", outlier.colour = NA)+
  geom_jitter(width = .25,  fill = "red",alpha = 0.5, shape = 21, color = "black")+
  scale_fill_manual(values = c("blue","orange","red"))+
  scale_size_manual(values = c(3,3,2))+
  scale_alpha_manual(values = c(1,1,.5))+
  facet_grid(~marker)+
  theme_bw()+
  theme(legend.position = 'none')+
  labs(x="",y = "Principal Component 1")+
  theme(legend.position = "none") + arsenic.theme

# PC1 fine-mapping

pca_genes <- data.table::fread(file = paste0(final.dir,"Supplemental_data18_GWAS_PCA_fine_mapping.tsv"))

fine_map_plot <- pca_genes%>%
  dplyr::filter(strain=="CB4856", trait == "pc1")%>%
  ggplot(.)+
  aes(x=POS/1e6, y = -log10(corrected_spearman_cor_p))+
  geom_vline(aes(xintercept = 7.83), color = "gray60")+
  geom_vline(aes(xintercept = 8.02), color = "gray60")+
  geom_point(aes(fill = GT),alpha = 0.7, shape = 21, color = "black", size =2)+
  scale_fill_manual(values=c("blue","orange"))+
  theme_bw()+
  theme(legend.position = 'none')+
  labs(x="Genomic Position (Mb)",y = expression(-log[10](italic(p))))+
  theme(legend.position = "none") +
  arsenic.theme

bottom_row <- plot_grid(pca_pxg, fine_map_plot, 
                        labels = c('B', 'C'), align = 'h', 
                        rel_widths = c(.75, 1), 
                        label_size = 18)

plot_grid(pca_gwa_plot[[1]] + 
          theme(axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16)), 
          bottom_row, labels = c('A', ''), ncol = 1, rel_heights = c(1, 1), label_size = 18)

ggsave(paste0(plot.dir,paste0("F2_GWA.pdf")), 
       height = 8, 
       width = 14)

# visualize PC1 correlation with animal length trait median.TOF

g8_regressed <- data.table::fread(file = paste0(final.dir,"TS13_GWAS_controlRegressed.tsv"))
arsenic_pca_df <- data.table::fread(file = paste0(final.dir,"TS15_GWAS_PCA_trait_DF.tsv"))

pc_size <- g8_regressed%>%
  dplyr::filter(trait == trait_of_interest)%>%
  dplyr::select(strain,trait,phenotype)%>%
  dplyr::left_join(.,arsenic_pca_df)

mod <- lm(PC1~phenotype,data=pc_size)
arsenic_eq <- transform(pc_size, Fitted = fitted(mod))
eq <- substitute(italic(r)~"="~rvalue*","~italic(p)~"="~pvalue, 
                 list(rvalue = sprintf("%.2f",sign(coef(mod)[2])*sqrt(summary(mod)$r.squared)), 
                      pvalue = format(summary(mod)$coefficients[2,4], digits = 3)))

dftext <- data.frame(PC1 = 18, phenotype = -50, eq = as.character(as.expression(eq)))

arsenic_eq%>%
  ggplot()+
  aes(x = phenotype, y = PC1)+
  geom_point()+
  theme_bw()+ 
  geom_text(aes(label = eq), data = dftext, parse = TRUE)+
  geom_smooth(se=FALSE, method = "lm")+
  theme(legend.position = "none")+
  labs(x=paste0("Animal Length (", trait_of_interest,")"), y = paste0("PC1"))+
  arsenic.theme

ggsave(paste0(plot.dir,"FS4.pdf"), 
       height = 6, 
       width = 8)

# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ##
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ##
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ # ALLELE SWAP EXPERIMENT # ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ #
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## 
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## 

controlregressed <- data.table::fread(file = paste0(final.dir,"TS8_NILs_SWAPs_controlRegressed.tsv"))

swap_size <- controlregressed%>%
  dplyr::filter(trait == trait_of_interest,
                strain != "ECA591", 
                strain != "ECA414", 
                strain != "ECA434", 
                strain != "ECA582", 
                strain != "ECA589")%>%
  dplyr::mutate(strain1 = factor(strain, 
                                 levels = c("N2","ECA581","CB4856","ECA590"),
                                 labels =c("N2\nDBT-1(C78)", 
                                           "N2\nDBT-1 (S78)",
                                           "CB4856\nDBT-1 (S78)",
                                           "CB4856\nDBT-1(C78)" ),
                                 ordered = T))%>%
  ggplot(.) +
  aes(x = factor(strain1), 
      y = phenotype, 
      fill=strain) +
  geom_jitter(width = 0.25, alpha = .4)+
  geom_boxplot(outlier.colour = NA, alpha = 0.7)+
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
  labs( y = paste0("Animal Length"))

swap_brood <- controlregressed%>%
  dplyr::filter(trait == "norm.n",
                strain != "ECA591", 
                strain != "ECA414", 
                strain != "ECA434", 
                strain != "ECA582", 
                strain != "ECA589")%>%
  dplyr::mutate(strain1 = factor(strain, 
                                 levels = c("N2","ECA581","CB4856","ECA590"),
                                 labels =c("N2\nDBT-1(C78)", 
                                           "N2\nDBT-1 (S78)",
                                           "CB4856\nDBT-1 (S78)",
                                           "CB4856\nDBT-1(C78)" ),
                                 ordered = T))%>%
  ggplot(.) +
  aes(x = factor(strain1), 
      y = phenotype, 
      fill=strain) +
  geom_jitter(width = 0.25, alpha = .4)+
  geom_boxplot(outlier.colour = NA, alpha = 0.7)+
  scale_fill_manual(values = c("N2" = "orange","CB4856" = "blue",
                               "ECA414" = "blue","ECA434" = "blue",
                               "ECA581" = "gray50","ECA582" = "gray50",
                               "ECA589" = "gray50","ECA590" = "gray50",
                               "ECA591" = "gray50"))+
  theme_bw()+
  arsenic.theme+
  theme(legend.position="none",
        axis.title.x = ggplot2::element_text(size = 0, face = "bold", color = "black", vjust = -0.3),
        axis.text.x = ggplot2::element_text(size = 0),
        axis.title.y = ggplot2::element_text(size = 14, face = "bold", color = "black", vjust = -0.3))+
  labs( y = paste0("Brood Size"))

plot_grid(swap_brood, swap_size, labels = "AUTO", ncol = 1, align = 'v', label_size = 18)

ggsave(paste0(plot.dir,paste0("F3_SWAP.pdf")), 
       height = 8, 
       width = 10)

# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ##
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ##
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ # ALLELE SWAP EXPERIMENT - PC1 as a trait # ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ #
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ##
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ##

controlregressed <- data.table::fread(file = paste0(final.dir,"TS8_NILs_SWAPs_controlRegressed.tsv"))
swap_pc_df <- data.table::fread(paste0(final.dir,"TS19_SWAP_PCA_Phenotypes.tsv"))

swap_pc1 <- swap_pc_df%>%
  tidyr::gather(trait, phenotype,-strain, -plate,-row,-col)%>%
  dplyr::filter(trait == "PC1",
                strain != "ECA591", 
                strain != "ECA414", 
                strain != "ECA434", 
                strain != "ECA582", 
                strain != "ECA589")%>%
  dplyr::mutate(strain1 = factor(strain, 
                                 levels = c("N2","ECA581","CB4856","ECA590"),
                                 labels =c("N2", 
                                           "N2\nDBT-1 (C78S)",
                                           "CB4856",
                                           "CB4856\nDBT-1(S78C)" ),
                                 ordered = T))%>%
  ggplot(.) +
  aes(x = factor(strain1), 
      y = phenotype, 
      fill=strain) +
  geom_jitter(width = 0.25, alpha = .4)+
  geom_boxplot(outlier.colour = NA, alpha = 0.7)+
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
  labs( y = paste0("PC1"))

# Show correlation between PC1 and median.TOF

swap_pc_median <- controlregressed%>%
  dplyr::filter(trait == trait_of_interest)%>%
  dplyr::select(strain, plate, row, col, phenotype)%>%
  dplyr::left_join(.,swap_pc_df, by = c("strain","plate","row","col"))

mod <- lm(PC1~phenotype,data=swap_pc_median)
# arsenic_eq <- transform(pc_size, Fitted = fitted(mod))
eq <- substitute(italic(r)~"="~rvalue*","~italic(p)~"="~pvalue, 
                 list(rvalue = sprintf("%.2f",sign(coef(mod)[2])*sqrt(summary(mod)$r.squared)), 
                      pvalue = format(summary(mod)$coefficients[2,4], digits = 3)))

dftext <- data.frame(PC1 = 18, phenotype = 25, eq = as.character(as.expression(eq)))

swap_pc_v_tof <- swap_pc_median%>%
  ggplot()+
  aes(x = phenotype, y = PC1)+
  geom_point(aes(color = strain))+
  theme_bw()+ 
  geom_text(aes(label = eq), data = dftext, parse = TRUE)+
  geom_smooth(se=FALSE, method = "lm", color = "red", alpha = 0.25)+
  scale_color_manual(values = c("N2" = "orange","CB4856" = "blue",
                                "ECA414" = "blue","ECA434" = "blue",
                                "ECA581" = "gray50","ECA582" = "gray50",
                                "ECA589" = "gray50","ECA590" = "gray50",
                                "ECA591" = "gray50"))+
  theme_bw()+
  arsenic.theme+
  theme(legend.position="none",
        axis.title.x = ggplot2::element_text(size = 16, face = "bold", color = "black", vjust = -0.3),
        axis.title.y = ggplot2::element_text(size = 16, face = "bold", color = "black", vjust = -0.3))+
  labs(x=paste0("Animal Length"), y = paste0("PC1"))

swap_pc_brood <- controlregressed%>%
  dplyr::filter(trait == "norm.n")%>%
  dplyr::select(strain, plate, row, col, phenotype)%>%
  dplyr::left_join(.,swap_pc_df, by = c("strain","plate","row","col"))

mod <- lm(PC1~phenotype,data=swap_pc_brood)
# arsenic_eq <- transform(pc_size, Fitted = fitted(mod))
eq <- substitute(italic(r)~"="~rvalue*","~italic(p)~"="~pvalue, 
                 list(rvalue = sprintf("%.2f",sign(coef(mod)[2])*sqrt(summary(mod)$r.squared)), 
                      pvalue = format(summary(mod)$coefficients[2,4], digits = 3)))

dftext <- data.frame(PC1 = 18, phenotype = 25, eq = as.character(as.expression(eq)))

swap_pc_v_brood <- swap_pc_brood%>%
  ggplot()+
  aes(x = phenotype, y = PC1)+
  geom_point(aes(color = strain))+
  theme_bw()+ 
  geom_text(aes(label = eq), data = dftext, parse = TRUE)+
  geom_smooth(se=FALSE, method = "lm", color = "red", alpha = 0.25)+
  scale_color_manual(values = c("N2" = "orange","CB4856" = "blue",
                                "ECA414" = "blue","ECA434" = "blue",
                                "ECA581" = "gray50","ECA582" = "gray50",
                                "ECA589" = "gray50","ECA590" = "gray50",
                                "ECA591" = "gray50"))+
  theme_bw()+
  arsenic.theme+
  theme(legend.position="none",
        axis.title.x = ggplot2::element_text(size = 16, face = "bold", color = "black", vjust = -0.3),
        axis.title.y = ggplot2::element_text(size = 16, face = "bold", color = "black", vjust = -0.3))+
  labs(x=paste0("Brood Size"), y = paste0("PC1"))

bottom_row <- plot_grid(swap_pc_v_brood, 
                        swap_pc_v_tof, 
                        labels = c('B', 'C'), 
                        align = 'h', 
                        rel_widths = c(1, 1), 
                        label_size = 18)

plot_grid(swap_pc1, 
          bottom_row, 
          labels = c('A', ''), 
          ncol = 1, 
          rel_heights = c(1, 1), 
          label_size = 18)

ggsave(paste0(plot.dir,paste0("FS5.pdf")), 
       height = 8, 
       width = 10)


# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ##
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ##
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ # Metabolite measurements in C. elegans   # ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ #
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ##
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ##

# measure metabolites

metabolites <- data.table::fread(paste0(final.dir, "Supplemental_data22_Metabolite_Measurements.tsv"))

rC15 <- metabolites%>%
  dplyr::filter(compound %in% c("rC15"), concentration != "200") 

met_ratio <- rC15%>%
  ggplot()+
  aes(x = factor(strain, 
                 levels = c("N2","CB","581","590"), 
                 labels = c("N2", "CB4856","N2\nDBT-1(S78)", "CB4856\nDBT-1(C78)")), 
      y = value, 
      fill = concentration)+
  geom_point(shape =21, size = 3)+
  scale_fill_manual(values = c("hotpink3","cadetblue3"))+
  theme_classic()+
  arsenic.theme+ theme(legend.position = "None",
                       axis.title.x = ggplot2::element_text(size = 0),
                       axis.title.y = ggplot2::element_text(size = 16))+
  labs(y = "C15iso/C15straight")

#c15iso
for( s in unique(metabolites$strain)){
  rC15 <- metabolites%>%
    dplyr::filter(compound %in% c("C15iso"), concentration != "200") %>%
    dplyr::filter(strain==s)
  
  print(s)
  fit <- aov(value ~ concentration , data=rC15)
  print(TukeyHSD(fit))
}
#c17iso
for( s in unique(metabolites$strain)){
  rC15 <- metabolites%>%
    dplyr::filter(compound %in% c("C17iso"), concentration != "200") %>%
    dplyr::filter(strain==s)
  
  print(s)
  fit <- aov(value ~ concentration , data=rC15)
  print(TukeyHSD(fit))
}
#rC17
for( s in unique(metabolites$strain)){
  rC15 <- metabolites%>%
    dplyr::filter(compound %in% c("rC17"), concentration != "200") %>%
    dplyr::filter(strain==s)
  
  print(s)
  fit <- aov(value ~ concentration , data=rC15)
  print(TukeyHSD(fit))
}
#rc15
for( s in unique(metabolites$strain)){
  rC15 <- metabolites%>%
    dplyr::filter(compound %in% c("rC15"), concentration != "200") %>%
    dplyr::filter(strain==s)
  
  print(s)
  fit <- aov(value ~ concentration , data=rC15)
  print(TukeyHSD(fit))
}
#c15n
for( s in unique(metabolites$strain)){
  rC15 <- metabolites%>%
    dplyr::filter(compound %in% c("C15n"), concentration != "200") %>%
    dplyr::filter(strain==s)
  
  print(s)
  fit <- aov(value ~ concentration , data=rC15)
  print(TukeyHSD(fit))
}

top_row <- plot_grid(met_ratio, 
                     met_ratio, 
                     labels = c('A', 'B'), 
                     align = 'h', 
                     rel_widths = c(1, 1), 
                     label_size = 18)


# c17

rC17 <- metabolites%>%
  dplyr::filter(compound %in% c("rC17")) 

rC17%>%
  ggplot()+
  aes(x = factor(strain, 
                 levels = c("N2","CB","581","590"), 
                 labels = c("N2", "CB4856","N2\nDBT-1(S78)", 
                            "CB4856\nDBT-1(C78)")), 
      y = value, 
      fill = concentration)+
  geom_point(shape =21, size = 3)+
  scale_fill_manual(values = c("hotpink3","cadetblue3"))+
  theme_classic()+
  arsenic.theme+ theme(legend.position = "None",
                       axis.title.x = ggplot2::element_text(size = 0),
                       axis.title.y = ggplot2::element_text(size = 14))+
  labs(y = "C17iso/C17straight")

ggsave(paste0(plot.dir,paste0("FS6.pdf")), 
       height = 4, 
       width = 6)



#  #  # # ## # # iso metabolites

c17 <- metabolites%>%
  dplyr::filter(compound %in% c("C17iso"))

c17plot <- c17%>%
  ggplot()+
  aes(x = factor(strain, levels = c("N2","CB","581","590"), 
                 labels = c("N2", "CB4856","N2\nDBT-1(S78)", "CB4856\nDBT-1(C78)")), 
      y = value, 
      fill = concentration)+
  geom_point(shape =21, size = 3)+
  scale_fill_manual(values = c("hotpink3","cadetblue3"))+
  theme_classic()+
  arsenic.theme+ theme(legend.position = "None",
                       axis.title.x = ggplot2::element_text(size = 0),
                       axis.title.y = ggplot2::element_text(size = 14))+
  labs(y = "C17ISO")

c15 <- metabolites%>%
  dplyr::filter(compound %in% c("C15iso"))

c15plot <- c15%>%
  ggplot()+
  aes(x = factor(strain, levels = c("N2","CB","581","590"), 
                 labels = c("N2", "CB4856", "N2\nDBT-1(S78)", "CB4856\nDBT-1(C78)")), 
      y = value, 
      fill = concentration)+
  geom_point(shape =21, size = 3)+
  scale_fill_manual(values = c("hotpink3","cadetblue3"))+
  theme_classic()+
  arsenic.theme+ theme(legend.position = "None",
                       axis.title.x = ggplot2::element_text(size = 0),
                       axis.title.y = ggplot2::element_text(size = 14))+
  labs(y = "C15ISO")

plot_grid(c15plot, 
          c17plot, 
          labels = "AUTO", 
          ncol = 1,
          align = 'v',
          label_size = 16)

ggsave(paste0(plot.dir,paste0("FS7.pdf")), 
       height = 8, 
       width = 10)


# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ #
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ #
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ # C15ISO RESCUE EXPERIMENT # ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ # 
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ #
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ #


processed_biopruned <- data.table::fread(file = paste0(final.dir,"TS21_RESCUE_processed.tsv"))

boxplot_plt(df = processed_biopruned,
            trt = trait_of_interest,
            cond = c("Arsenic",
                     "C15ISO64", 
                     "ArsenicC15ISO64"),
            fancy_name = paste0("Animal Length"),
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
            r_conc = "64")+ theme(axis.text.x = ggplot2::element_text(size = 12),
                                  axis.title.y = element_text(size =16))

ggsave(paste0(plot.dir,paste0("FS9.pdf")), 
       height = 4, 
       width = 14)



main_display_rescue <- boxplot_plt(df = processed_biopruned,
                                   trt = trait_of_interest,
                                   cond = c("Arsenic",
                                            "ArsenicC15ISO64"),
                                   fancy_name = paste0("Animal Length"),
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
                                   r_conc = "64")+
  arsenic.theme+ theme(legend.position = "None",
                       axis.title.x = ggplot2::element_text(size = 0),
                       axis.title.y = ggplot2::element_text(size = 16))

plot_grid(top_row, 
          main_display_rescue, 
          labels = c('', 'C'), 
          ncol = 1, 
          rel_heights = c(0.7, 1), 
          label_size = 18)

ggsave(paste0(plot.dir,paste0("F4.pdf")), 
       height = 8, 
       width = 12)

# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## 
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## 
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ #  Human results  # ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ # 
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ##
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## 

celldata <- data.table::fread(paste0(final.dir,"TS23_Human_cell_allele_swap.csv"))

zeros <- celldata%>%
  dplyr::filter(Concentration == 0)%>%
  dplyr::group_by(Edit, Concentration)%>%
  dplyr::mutate(m_f_zero = mean(Fraction_edit))%>%
  dplyr::ungroup()%>%
  dplyr::select(Edit, Fraction_edit_zero = Fraction_edit, m_f_zero)

subtract_zero <- celldata%>%
  dplyr::filter(Concentration == 5)%>%
  dplyr::group_by(Edit, Concentration)%>%
  dplyr::mutate(m_f_five = mean(Fraction_edit))%>%
  dplyr::ungroup()%>%
  dplyr::select(Edit, Fraction_edit_five = Fraction_edit, m_f_five) %>%
  dplyr::left_join(., zeros, by = "Edit") %>%
  dplyr::mutate(c_zero = Fraction_edit_zero - m_f_zero,
                c_five = Fraction_edit_five - m_f_zero)%>%
  dplyr::select(Edit, c_zero, c_five)%>%
  dplyr::distinct(Edit, c_zero, c_five)%>%
  dplyr::group_by(Edit)%>%
  dplyr::mutate(m_c_five = mean(c_five),
                m_c_zero = mean(c_zero))

subtract_zero%>% 
  tidyr::gather(Conc, Fraction_edit, -Edit, -m_c_zero, -m_c_five)%>%
  dplyr::mutate(Concentration = ifelse(Conc == "c_five", 5, 0))%>%
  dplyr::distinct(Edit, Fraction_edit, Concentration, .keep_all =T)%>%
  dplyr::select(-Conc) %>%
  tidyr::gather(Means, mean_value, -Edit, -Fraction_edit, -Concentration) %>%
  dplyr::mutate(Mean_Concentration = ifelse(Means == "m_c_five", 5, 0))%>%
  dplyr::filter(Mean_Concentration == Concentration)%>%
  ggplot()+
  aes(x = factor(Concentration), y = Fraction_edit, 
      color = factor(Edit,
                     levels = c("W84C","S112C","R113C")), 
      fill = factor(Edit,
                    levels = c("W84C","S112C","R113C")))+
  geom_point(size =3, shape = 22, aes(x = factor(Concentration), y = mean_value))+
  geom_line( aes(x = factor(Concentration), y = mean_value, group = Edit))+
  scale_color_manual(values = c("hotpink3", "cadetblue3","black"), name = "Edit:")+
  scale_fill_manual(values = c("hotpink3", "cadetblue3","black"), name = "Edit:")+
  labs(x = "Arsenic Concentration (µM)", y = "Percent Edited Cells\nRelative to Control")+
  theme_classic(18)+
  arsenic.theme+
  theme(legend.position = "top",
        legend.text = element_text(size = 16))+
  ylim(c(0,5))+
  scale_x_discrete(expand=c(.1,0)) 


ggsave(paste0(plot.dir,"F5.pdf"), 
       height = 5, 
       width = 6)


# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ##
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## 
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ # Tajima's D across GWAS interval # ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ #
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## 
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ##

dbt <- tajimas_d_temp(vcf_path = paste0(final.dir),
                      vcf_name = "Supplemental_data25_WI.20160408.impute.vcf.gz", 
                      chromosome = "II", 
                      interval_start = 7.43e6, 
                      interval_end = 8.33e6,
                      site_of_interest = 7942463,
                      slide_distance = 10,
                      window_size = 100)

dbt

ggsave(paste0(plot.dir,paste0("FS11.pdf")), 
       height = 4, 
       width = 12)

# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## 
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ##
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ #  C78S world-wide allele distribution # ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~  # 
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## 
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## 

isolation_info <- data.table::fread(paste0(final.dir,"TS26_Strain_Isolation_Info.tsv"))

world <- map_data(map = "world")
world <- world[world$region != "Antarctica",] # intercourse antarctica

# world pop
ggplot()+ geom_map(data=world, map=world,
                   aes(x=long, y=lat, map_id=region),
                   color="white", fill="#7f7f7f", size=0.05, alpha=1/4)+
  geom_point(data = isolation_info, aes(x=long, y=lat, fill=GT), shape =21, alpha = 0.7) + 
  scale_fill_manual(values = c("blue", "orange"),name = "Cys78Ser")+
  theme_map()

ggsave(paste0(plot.dir,"FS12.pdf"), 
       height = 5, 
       width = 10)
