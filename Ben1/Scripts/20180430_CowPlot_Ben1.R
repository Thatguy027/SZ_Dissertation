# Finalized figures for BEN-1 paper

# calculate ka/ks with http://services.cbu.uib.no/tools/kaks
# use "seqdump.txt" in data/custom_Popgen folder
# ka/ks Ce - Cbr = 

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # 
# set to location of files
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # 
main.dir <- "~/Dropbox/HTA/Results/2017_anthelmintic_gwa/analysis/ben1_paper/"
data.dir <- paste0(main.dir,"Data/")
plot.dir <- paste0(main.dir,"Plots/")
script.dir <- paste0(main.dir,"Scripts/")

source(paste0(script.dir,"ben1_processing_functions.R"))

trait_of_interest <- "q90.TOF"
PCA_of_interest <- "PC1"
comparison_trait <- "q90.EXT"
condition_of_interest <- "albendazole"
control_oof_interest <- "DMSO"

plot.theme <-  theme(axis.text.x = ggplot2::element_text(size = 14),
                     axis.text.y = ggplot2::element_text(size = 14),
                     axis.title.x = ggplot2::element_text(size = 18, face = "bold", color = "black", vjust = -0.3),
                     axis.title.y = ggplot2::element_text(size = 18, face = "bold", color = "black", vjust = -0.3),
                     strip.text.x = element_text(size = 14, face = "bold"))

library(cowplot)

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # 
# PLOT DOSE RESPONSES 
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

dr_data <- data.table::fread(paste0(data.dir,"DR_Processed.tsv"))

dr_data$strain <- gsub("Ju", "JU", dr_data$strain)
pr_df <- dr_data%>%
  dplyr::ungroup()%>%
  dplyr::filter(trait == "norm.n", grepl("alb|DMSO",condition,ignore.case = T))%>%
  dplyr::filter(!is.na(condition))

control_means <- pr_df%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::select(strain, value)%>%
  dplyr::group_by(strain)%>%
  dplyr::summarise(m_ctrl = mean(value, na.rm = T))

subtract_controls <- pr_df%>%
  dplyr::left_join(.,control_means, by = "strain")%>%
  dplyr::mutate(subtracted_pheno = value - m_ctrl)


brood_DR <- subtract_controls%>%
  ggplot()+
  aes(x = factor(condition, 
                 levels = c("DMSO", "albendazole3125", "albendazole625", 
                            "albendazole125",  "albendazole25"),
                 labels = c("0", "3.125", "6.25", "12.5", "25")), 
      y = subtracted_pheno, 
      fill = strain)+
  geom_boxplot(outlier.colour = NA)+
  scale_fill_manual(values = c("blue", "cadetblue3", "hotpink3", "orange"), name = "Strain")+
  theme_bw()+
  labs(x = "Albendazole Concentration (µM)",
       y = "Brood Size") + plot.theme + theme(legend.position = "none")

# length
pr_df <- dr_data%>%
  dplyr::ungroup()%>%
  dplyr::filter(trait == trait_of_interest, grepl("alb|DMSO",condition,ignore.case = T))%>%
  dplyr::filter(!is.na(condition))

control_means <- pr_df%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::select(strain, value)%>%
  dplyr::group_by(strain)%>%
  dplyr::summarise(m_ctrl = mean(value, na.rm = T))

subtract_controls <- pr_df%>%
  dplyr::left_join(.,control_means, by = "strain")%>%
  dplyr::mutate(subtracted_pheno = value - m_ctrl)


size_DR <- subtract_controls%>%
  ggplot()+
  aes(x = factor(condition, 
                 levels = c("DMSO", "albendazole3125", "albendazole625", 
                            "albendazole125",  "albendazole25"),
                 labels = c("0", "3.125", "6.25", "12.5", "25")), 
      y = subtracted_pheno, 
      fill = strain)+
  geom_boxplot(outlier.colour = NA)+
  scale_fill_manual(values = c("blue", "cadetblue3", "hotpink3", "orange"), name = "Strain")+
  theme_bw()+
  labs(x = "Albendazole (µM)",
       y = "Animal Length") + plot.theme + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  theme(legend.position = "top", 
        legend.background = element_rect(color = "black", size = 1, linetype = "solid"),
        legend.direction = "horizontal",
        legend.justification = "right")


plot_grid( size_DR,brood_DR, labels = "AUTO", ncol = 1, align = 'v')


ggsave(paste0(plot.dir,"FS1_DR_",condition_of_interest,"_",trait_of_interest,".png"), 
       dpi = 300,
       height = 6, 
       width = 10)

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # 
# Correlation between brood size and animal length for dose response experiments.
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

dr_data$strain <- gsub("Ju", "JU", dr_data$strain)
pr_df <- dr_data%>%
  dplyr::ungroup()%>%
  dplyr::filter(trait == "norm.n" | trait == trait_of_interest, grepl("alb|DMSO",condition,ignore.case = T))%>%
  dplyr::filter(!is.na(condition))

control_means <- pr_df%>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::select(strain, trait, value)%>%
  dplyr::group_by(strain, trait)%>%
  dplyr::summarise(m_ctrl = mean(value, na.rm = T))

subtract_controls <- pr_df%>%
  dplyr::left_join(.,control_means, by = c("strain","trait"))%>%
  dplyr::mutate(subtracted_pheno = value - m_ctrl)%>%
  dplyr::select(condition, strain, plate, row,col, trait, subtracted_pheno) %>%
  dplyr::ungroup()%>%
  # dplyr::mutate(i = row_number())%>%
  tidyr::spread(trait, subtracted_pheno)


subtract_controls%>%
  ggplot()+
  aes(x = q90.TOF, 
      y = norm.n, 
      fill = strain)+
  geom_point( shape = 21, size = 2)+
  scale_fill_manual(values = c("blue", "cadetblue3", "hotpink3", "orange"), name = "Strain:")+
  theme_bw()+
  labs(x = "Animal Length",
       y = "Brood Size") + plot.theme +
  theme(legend.position = "top", 
        legend.background = element_rect(color = "black", size = 1, linetype = "solid"),
        legend.direction = "horizontal",
        legend.justification = "right")

ggsave(paste0(plot.dir,"FS2_DR_correlation_",trait_of_interest,"_v_broodsize",".png"), 
       dpi = 300,
       height = 6, 
       width = 8)

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # 
# PLOT GWAS and BURDEN
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

pr_maps <- data.table::fread(paste0(data.dir,"GWAS_processed_mapping.tsv"))

snv_manplot <- manplot_edit(pr_maps)[[2]] +
  theme(axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank())

load(paste0(data.dir, "q90TOF.VariableThresholdPrice_processed.Rda"))

burden_manplot <- q90burden_pr%>%
  ggplot()+
  aes(x = POS/1e6, y = Stat, size = NumVar, alpha = 0.5, color = significant)+
  geom_point()+
  scale_color_manual(values=c("black","red"))+
  facet_grid(.~CHROM, scales = "free", space = "free")+
  theme_bw()+
  plot.theme+
  ggplot2::theme(legend.position = "none")+
  labs(x = "Genomic Position (Mb)", y = "Test Statisitic")

plot_grid(snv_manplot, burden_manplot, labels = "AUTO", ncol = 1, align = 'v')

ggsave(paste0(plot.dir,"F1_GWA-BURDEN_",condition_of_interest,"_",trait_of_interest,".png"), 
       dpi = 300,
       height = 6, 
       width = 12)

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # 
# PLOT GWAS and BURDEN
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

unique(pr_maps$trait)

animal_length_pxg <- pr_maps %>%
  dplyr::filter(trait == "q90.tof")%>%
  na.omit()%>%
  dplyr::arrange(CHROM,POS)%>%
  dplyr::mutate(marker2 = factor(marker, levels = unique(marker), labels = unique(marker)))%>%
  ggplot(.)+
  aes(x=factor(as.character(allele), labels = c("REF","ALT")), y = value)+
  geom_boxplot(fill = "gray90", outlier.colour = NA)+
  geom_jitter(width = .25,  fill = "red",alpha = 0.5, shape = 21, color = "black")+
  scale_fill_manual(values = c("blue","orange","red"))+
  scale_size_manual(values = c(3,3,2))+
  scale_alpha_manual(values = c(1,1,.5))+
  facet_grid(~marker2)+
  theme_bw()+
  theme(legend.position = 'none')+
  labs(x="",y = "Animal Length")+
  theme(legend.position = "none") + plot.theme

animal_length_qtl_ld <- modified_ld_plot(plot_df = pr_maps, trait = "q90.TOF")

bottom_row <- plot_grid(animal_length_qtl_ld, scale = 1)

plot_grid(animal_length_pxg, bottom_row, labels = "AUTO", ncol = 1, align = 'v')

ggsave(paste0(plot.dir,"FS3_GWA-SNV_PXG_LD_",condition_of_interest,"_",trait_of_interest,".png"), 
       dpi = 300,
       height = 6, 
       width = 8)

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # 
# BAR PLOT OF BEN-1 VARIATION
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

# # # LOAD AND PROCESS INDELS
ben1_indels <- data.table::fread(paste0(data.dir, "ben1_Indels/20171021_ben1_indels.csv"))%>%
  na.omit()

pr_indels <- ben1_indels %>%
  tidyr::unite(marker, Type, Start, End, sep ="_",remove=F)%>%
  dplyr::filter(comments == "frame shift" | grepl("Exon", location) | Type == "trans")%>%
  dplyr::filter(!grepl("Intron",location,ignore.case = T))%>%
  dplyr::mutate(GT = ifelse(marker == "_NA_NA", "REF", "ALT"))%>%
  dplyr::select(marker, strain = Strain, GT,start = Start,end =End)%>%
  dplyr::distinct(strain, marker, GT,.keep_all=T)

# # # LOAD AND PROCESS SNPS
ben1_snps <- cegwas::snpeff("ben-1",severity = "ALL",elements = "ALL")

ben1_snps_pr <- ben1_snps%>%
  dplyr::select(POS, strain, GT, nt_change, effect, impact)%>%
  tidyr::unite(marker, effect,impact,nt_change,POS,sep="_",remove=F)%>%
  dplyr::select(marker, strain, GT,start = POS)%>%
  dplyr::mutate(end = start)%>%
  dplyr::filter(!grepl("MODIF|syn",marker))

# # # COMBINE SNPS AND INDELS
ben1_variants <- dplyr::bind_rows(ben1_snps_pr,pr_indels)%>%
  dplyr::filter(GT=="ALT")

gwa_mappings <- data.table::fread(file = paste0(data.dir,"GWAS_processed_mapping.tsv"))%>%
  na.omit()%>%
  dplyr::mutate(snpGT = ifelse(allele==-1,"REF", "ALT"))%>%
  dplyr::select(snpMarker = marker, strain, trait, value, snpGT)%>%
  dplyr::filter(trait == "q90.tof")%>%
  dplyr::left_join(.,ben1_variants, by = "strain")

gwa_mappings$snpMarker <- gsub("_",":", gwa_mappings$snpMarker)
colors <- c("#2121D9","#DF0101","#04B404","#FF9326","#A945FF","#0089B2","#B26314","#610B5E","#FE2E9A","#BFF217")
gwa_mappings <- gwa_mappings%>%
  dplyr::mutate(length = end-start)%>%
  dplyr::mutate(ben1_prediction = ifelse(is.na(GT), "A", 
                                         ifelse(grepl("del", marker),"Deletion",
                                                ifelse(grepl("ins", marker,ignore.case = T),"Insertion",
                                                       ifelse(grepl("inv", marker),"Inversion",
                                                              ifelse(grepl("stop", marker),"Stop Gained",
                                                                     ifelse(grepl("trans", marker),"Transposon Insertion",
                                                                            ifelse(grepl("missense", marker),"Missense","Splice Donor"))))))))




gwa_mappings_bar <- gwa_mappings %>%
  dplyr::distinct(strain, value, marker, ben1_prediction)%>%
  dplyr::arrange(value)%>%
  dplyr::mutate(strain2 = factor(strain, levels = unique(strain), labels = unique(strain), ordered = T))%>%
  dplyr::distinct(strain, value, .keep_all=T)%>%
  dplyr::mutate(norm_pheno_temp = ifelse(value == min(value), 0, 1))%>%
  dplyr::mutate(delta_pheno = ifelse(norm_pheno_temp == 0, 0, 
                                     ifelse(norm_pheno_temp == 1 & value < 0 , abs(lag(value) - value) , 
                                            ifelse(norm_pheno_temp == 1 & value > 0 ,  value - lag(value) , NA ))))%>%
  dplyr::mutate(norm_pheno = cumsum(delta_pheno)) %>%
  dplyr::mutate(final_pheno = norm_pheno/max(norm_pheno))

ben1_all <- gwa_mappings_bar%>%
  dplyr::mutate(build_GT = ifelse(ben1_prediction %in% c("Missense","Stop Gained", "Splice Donor","Deletion", "Insertion", "Inversion", "Transposon Insertion"), ben1_prediction, "A"))

ben1_all$build_GT[is.na(ben1_all$build_GT)] <- "A"

library(lemon)

ben1_variation_bar <- plot_bar(ben1_all) + plot.theme + 
  theme(axis.text.x = ggplot2::element_blank(),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = )+
  theme(legend.position =c(0,1), 
        legend.background = element_rect(color = "black", size = 1, linetype = "solid"),
        legend.direction = "horizontal",
        legend.justification = "left")

reposition_legend(ben1_variation_bar, 'top left')

ggsave(plot = reposition_legend(ben1_variation_bar, 'top left'),
       paste0(plot.dir,"F2_ben1_resistance_barplot",condition_of_interest,"_",trait_of_interest,".png"), 
       dpi = 300,
       height = 4, 
       width = 10)
ggsave(plot = reposition_legend(ben1_variation_bar, 'top left'),
       paste0(plot.dir,"F2_ben1_resistance_barplot",condition_of_interest,"_",trait_of_interest,".pdf"), 
       dpi = 300,
       height = 4, 
       width = 10)

pr_resid_maps <- data.table::fread( paste0(data.dir,"GWAS_ben1_variants_regressed_from_",trait_of_interest,"_processed_mapping.tsv"))
manhattan_plots <- manplot_edit(pr_resid_maps)

ben1_resid_manplot <- manhattan_plots[[2]]+theme_bw(15) +plot.theme + theme(legend.position = "none")

ben1_resid_genes <- data.table::fread( paste0(data.dir,"GWAS_finemapping_genes_ben1_regressed.tsv"))


fine_map <- ben1_resid_genes %>%
  dplyr::filter(strain == "EG4725",
                effect %in% c("missense_variant",
                              "splice_donor_variant&intron_variant",
                              "stop_gained",
                              "splice_acceptor_variant&intron_variant",
                              "splice_region_variant&non_coding_exon_variant",
                              "missense_variant&splice_region_variant",
                              "splice_region_variant"))%>%
  dplyr::mutate(tidy_effect = factor(effect, levels = c("missense_variant",
                                                        "splice_donor_variant&intron_variant",
                                                        "stop_gained",
                                                        "splice_acceptor_variant&intron_variant",
                                                        "splice_region_variant&non_coding_exon_variant",
                                                        "missense_variant&splice_region_variant",
                                                        "splice_region_variant"),
                                     labels = c("Missense",
                                                "Splice Donor\nIntron",
                                                "Stop Gained",
                                                "Splice Acceptor\nIntron",
                                                "Splice Region\nExon",
                                                "Splice Region\nMissense",
                                                "Splice Region")))%>%
  ggplot()+
  aes(x = POS/1e6, y = -log10(corrected_spearman_cor_p), 
      fill=tidy_effect,color = tidy_effect, shape = tidy_effect)+
  geom_point(size =2)+
  ggplot2::labs(x = "Genomic Position (Mb)",
                y = expression(bold(-log[10](bolditalic(p)))))+
  scale_color_brewer(palette="Set1", name = "Effect:")+
  scale_fill_brewer(palette="Set1", name = "Effect:")+
  scale_shape_manual(values = c(16:25), name = "Effect:")+
  theme_bw(15)+
  plot.theme+
  theme(axis.title.y = element_blank())


plot_grid(ben1_resid_manplot, fine_map, labels = c("A","B"), 
          label_size = 18, ncol = 2, align = 'v', rel_widths = c(2,1))

ggsave(paste0(plot.dir,"F5_Ben-1_variation_regressed_manplot_finemap_",condition_of_interest,"_",trait_of_interest,".png"), 
       dpi = 300,
       height = 6, 
       width = 18)
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # 
# CHRX PxG
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

regress_ben1 <- data.table::fread(paste0(data.dir,"GWAS_",trait_of_interest,"_with_ben1_covariate.tsv"))

x_pg_df <- pr_resid_maps %>%
       dplyr::filter(trait == "regressed_trait")%>%
       na.omit()%>%
       dplyr::select(qtl_marker=marker, strain, reg_val = value, qtl_allele = allele)%>%
       dplyr::left_join(regress_ben1,.,by ="strain")

x_pg_df%>%
  ggplot()+
  aes(x = factor(qtl_allele, levels = c(-1,1),
                 labels = c("REF", "ALT")), y = reg_val)+
  geom_boxplot(outlier.colour = NA)+
  geom_jitter(width = 0.25, aes(fill = factor(covariate)), shape = 21, size =2)+
  scale_fill_manual(values = c("cadetblue3","hotpink3"), labels = c("REF", "ALT"), name=expression(italic(ben-1))) +
  theme_bw(15)+
  labs(x="QTL allele", y = "Regressed Animal Length")+
  theme(axis.title.x = ggplot2::element_text(size = 18, face = "bold", color = "black", vjust = -0.3), 
        axis.title.y = ggplot2::element_text(size = 18, face = "bold", color = "black"))

ggsave(paste0(plot.dir,"GWAS_Xqtl_split_colored_by_ben1_variation.png"), 
       dpi = 300,
       height = 6, 
       width = 8)

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # 
# BEN-1 HAPLOTYPE GENE MODEL PLOT
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

# # # LOAD AND PROCESS INDELS
ben1_indels <- data.table::fread(paste0(data.dir, "ben1_Indels/20171021_ben1_indels.csv"))%>%
  na.omit()

pr_indels <- ben1_indels %>%
  tidyr::unite(marker, Type, Start, End, sep ="_",remove=F)%>%
  dplyr::filter(comments == "frame shift" | grepl("Exon", location) | Type == "trans")%>%
  dplyr::filter(!grepl("Intron",location,ignore.case = T))%>%
  dplyr::mutate(GT = ifelse(marker == "_NA_NA", "REF", "ALT"))%>%
  dplyr::select(marker, strain = Strain, GT,start = Start,end =End)%>%
  dplyr::distinct(strain, marker, GT,.keep_all=T)

# # # LOAD AND PROCESS SNPS
ben1_snps <- cegwas::snpeff("ben-1",severity = "ALL",elements = "ALL")

ben1_snps_pr <- ben1_snps%>%
  dplyr::select(POS, strain, GT, nt_change, effect, impact)%>%
  tidyr::unite(marker, effect,impact,nt_change,POS,sep="_",remove=F)%>%
  dplyr::select(marker, strain, GT,start = POS)%>%
  dplyr::mutate(end = start)%>%
  dplyr::filter(!grepl("MODIF|syn",marker))

# # # COMBINE SNPS AND INDELS
ben1_variants <- dplyr::bind_rows(ben1_snps_pr,pr_indels)%>%
  dplyr::filter(GT=="ALT")


for_plot <- data.table::fread(file = paste0(data.dir,"GWAS_processed_mapping.tsv"))%>%
  na.omit()%>%
  dplyr::mutate(snpGT = ifelse(allele==-1,"REF", "ALT"))%>%
  # dplyr::select(snpMarker = marker, strain, trait, value, snpGT)%>%
  dplyr::filter(trait == "q90.tof")%>%
  na.omit()%>%
  dplyr::filter(peak_id==1)%>%
  dplyr::select(strain, phenotype = value) %>%
  dplyr::left_join(.,ben1_variants,by="strain")


# construct ben-1 model
ben1_gene_info <- data.frame(feature = c("5'UTR", "Exon", "Intron", "Exon", "Intron", "Exon", "Intron", "Exon", "Intron", "Exon", "3'UTR"),
                             start = c(3541630,3541595,3541432,3540411,3540065,3540004,3539841,3539784,3539310,3538501,3538325),
                             stop = c(3541596,3541433,3540412,3540064,3540005,3539842,3539785,3539311,3538502,3538326,3538293),
                             color = c("blue","orange","gray","orange","gray","orange","gray","orange","gray","orange","blue"),
                             min_y = min(for_plot$phenotype)-10,
                             max_y = max(for_plot$phenotype)+10,
                             gene = "ben-1",
                             strand = "Watson")%>%
  dplyr::mutate(size = start-stop)


plot_df<- for_plot%>%
  dplyr::filter(GT != "HET", !is.na(GT), GT == "ALT")%>%
  dplyr::filter(!grepl("intron_variant_MODIFIER|modif|synon|del_3542405_3542407",marker,ignore.case = T))%>%
  dplyr::arrange(start)%>%
  dplyr::select(strain, phenotype, marker)%>%
  dplyr::mutate(GT = 1)%>%
  tidyr::spread(marker, GT)

plot_df[is.na(plot_df)] <- 0

ben1Cluster <- kmeans(plot_df[, 4:ncol(plot_df)], 26, nstart = 20)
max(ben1Cluster$cluster)

ben1haps <- data.frame(cluster= ben1Cluster$cluster, plot_df)%>%
  dplyr::group_by(cluster)%>%
  dplyr::mutate(mean_pheno = mean(phenotype))%>%
  dplyr::ungroup()%>%
  dplyr::arrange((mean_pheno))%>%
  dplyr::mutate(f_clust = factor(cluster, levels = unique(cluster), ordered = T))%>%
  tidyr::gather(marker, GT, -strain, -phenotype, -cluster, -f_clust,-mean_pheno)

marker_starts <- dplyr::select(for_plot,marker, start,end)%>%
  dplyr::distinct()

marker_starts$marker <- gsub(">","\\.",marker_starts$marker)
marker_starts$marker <- gsub("-","\\.",marker_starts$marker)
marker_starts$marker <- gsub("\\+","\\.",marker_starts$marker)
marker_starts$marker <- gsub("&","\\.",marker_starts$marker)

ben1haps1 <- dplyr::left_join(ben1haps, marker_starts, by = "marker")%>%
  dplyr::mutate(plot_shape = ifelse(grepl("ins",marker,ignore.case = T), "insertion",
                                    ifelse(grepl("del",marker), "deletion",
                                           ifelse(grepl("missense",marker), "missense",
                                                  ifelse(grepl("stop",marker), "stop_gained",
                                                         ifelse(grepl("spli",marker),"spli",
                                                                ifelse(grepl("MODIF",marker),"mod",
                                                                       ifelse(grepl("inv",marker),"inv", 
                                                                              ifelse(grepl("trans",marker),"trans", "other")))))))))%>%
  dplyr::filter(GT!=0)%>%
  dplyr::group_by(cluster)%>%
  dplyr::mutate(hap_strains= paste(unique(strain), collapse = " "))%>%
  dplyr::distinct(cluster,marker,.keep_all=T)%>%
  dplyr::ungroup()%>%
  dplyr::arrange(mean_pheno)

pt_size <- 5
ben1_high_var<-ggplot(ben1haps1)+
  aes(x=start, 
      y = factor(f_clust,labels=unique(hap_strains),levels=unique(cluster), ordered=T), 
      fill = factor(GT), 
      shape = plot_shape,
      size = plot_shape)+
  scale_y_discrete(labels = scales::wrap_format(25),position = "right")+
  geom_segment(aes(x = as.numeric(start), 
                   xend = as.numeric(end), 
                   yend = factor(f_clust,labels=unique(hap_strains),levels=unique(cluster), ordered=T)), 
               size =5, 
               color = "red")+
  geom_point( color ="black")+
  scale_shape_manual(values = c("insertion" = 73, 
                                "deletion" = 2, 
                                "missense" = 77, 
                                "stop_gained" = 13, 
                                "mod"=96,
                                "spli" = 83,
                                "trans" = 84,
                                "inv" = 37,
                                other= 16))+
  scale_size_manual(values = c("insertion" = pt_size,
                               "deletion" = pt_size,
                               "missense" = pt_size,
                               "stop_gained" = pt_size,
                               "mod"=pt_size,
                               "spli" = pt_size,
                               "trans"=pt_size,
                               "inv" = pt_size,
                               "other"= pt_size))+
  theme_bw()+
  labs(y="",x = "Genomic Position (Mb)")+
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=14, face="bold", color="black"),
        axis.title.y = element_text(size=14, face="bold", color="black",vjust=1),
        plot.title = element_text(size=10, face="bold",vjust=1),
        legend.title=element_blank(),
        legend.position = "none")+
  geom_vline(aes(xintercept=c(3541595)), linetype="dotdash", alpha = .25, color = "blue")+
  geom_vline(aes(xintercept=c(3541432)), linetype="dotdash", alpha = .25, color = "blue")+
  geom_vline(aes(xintercept=c(3540411)), linetype="dotdash", alpha = .25, color = "blue")+
  geom_vline(aes(xintercept=c(3540065)), linetype="dotdash", alpha = .25, color = "blue")+
  geom_vline(aes(xintercept=c(3540004)), linetype="dotdash", alpha = .25, color = "blue")+
  geom_vline(aes(xintercept=c(3539841)), linetype="dotdash", alpha = .25, color = "blue")+
  geom_vline(aes(xintercept=c(3539784)), linetype="dotdash", alpha = .25, color = "blue")+
  geom_vline(aes(xintercept=c(3539310)), linetype="dotdash", alpha = .25, color = "blue")+
  geom_vline(aes(xintercept=c(3538501)), linetype="dotdash", alpha = .25, color = "blue") +
  geom_vline(aes(xintercept=c(3538325)), linetype="dotdash", alpha = .25, color = "blue")

missense_df <- ben1_snps%>%
  dplyr::select(CHROM,POS, aa_change,effect)%>%
  dplyr::filter(effect == "missense_variant")%>%
  dplyr::distinct(aa_change, .keep_all =T)%>%
  dplyr::mutate(plot_aa = gsub(pattern = "p\\.","",aa_change))

missense_df$letter <- c("D404N", "M257I", "F200Y", "A185P","S145F","Q131L","E69G","K58E")

missense_df_low <- dplyr::filter(missense_df, letter %in% c("M257I","A185P","Q131L","K58E"))
missense_df_high <- dplyr::filter(missense_df, !(letter %in% c("M257I","A185P","Q131L","K58E")))

ben1_model <- gene_model(ben1_gene_info)+
  ylim(-4,4)+
  geom_segment(aes(x = as.numeric(POS), y = -1, xend = as.numeric(POS), yend = -2),data = missense_df_low, size =2)+
  geom_segment(aes(x = as.numeric(POS), y = 1, xend = as.numeric(POS), yend = 2),data = missense_df_high, size =2)+
  geom_label(aes(x = POS, y = -3, label = letter), data = missense_df_low,size = 2)+
  geom_label(aes(x = POS, y = 3, label = letter), data = missense_df_high,size = 2)


cowplot::ggdraw() +
  cowplot::draw_plot(ben1_high_var, 0, 0, 1, .8)+
  cowplot::draw_label(label = expression(bolditalic('ben-1')), 0.5, y = .92)+
  cowplot::draw_plot(ben1_model, .13, .81, .52, .1)+
  draw_plot_label(c("A", "B"), c(0, 0), c(1, 0.8), size = 14)

ggsave(paste0(plot.dir,"F5_Ben1_haps_by_phenotype.png"),
       dpi = 300,
       height = 12, 
       width = 18)

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # 
# BEN-1 phylogeny + BAR
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

library(ape)
library(phangorn)
library(seqinr)
library(ggtree)
library(distory)
library(GGally)
library(data.table)
library(pophelper)
library(RColorBrewer)
# source("https://bioconductor.org/biocLite.R")
# biocLite("ggtree")

# # # LOAD AND PROCESS INDELS
gwa_mappings <- data.table::fread(paste0(data.dir, "Processed_Ben1_variants_with_pheno.tsv"))

base.dir <- "~/Dropbox/AndersenLab/LabFolders/Stefan/indels/Post_calling_analysis/"

tree <-  ape::read.tree( file = paste0(data.dir, "genome.tree"))

tree <- ggtree(tree, 
               branch.length = "rate",layout = "circular") 

ben1_variation_matrix <- gwa_mappings%>%
  dplyr::ungroup()%>%
  dplyr::select(strain,marker)%>%
  dplyr::distinct()

missing_strains <- data.frame(strain = row.names(cegwas::kinship)[!row.names(cegwas::kinship)%in%ben1_variation_matrix$strain],
                              marker = NA, stringsAsFactors = F)
colors <- c("#2121D9","#DF0101","#04B404","#FF9326","#A945FF","#0089B2","#B26314","#610B5E","#FE2E9A","#BFF217")

ben1_strain_markers <- ben1_variation_matrix %>%
  dplyr::bind_rows(.,missing_strains)%>%
  dplyr::distinct()%>%
  dplyr::filter(!grepl("missense_variant_MODERATE_c.599T>A_353", marker))%>%
  dplyr::mutate(ben1_prediction = ifelse(is.na(marker), NA, 
                                         ifelse(grepl("del", marker),"Deletion",
                                                ifelse(grepl("ins", marker,ignore.case = T),"Insertion",
                                                       ifelse(grepl("inv", marker),"Inversion",
                                                              ifelse(grepl("stop", marker),"Stop Gained",
                                                                     ifelse(grepl("trans", marker),"Transposon Insertion",
                                                                            ifelse(grepl("missense", marker),"Missense","Splice Donor"))))))))

ben1_tree <- tree %<+% ben1_strain_markers+ 
  geom_tiplab(aes(color=ben1_prediction))+
  scale_color_manual(values=colors,name = expression(paste("Variation at ", italic("ben-1"))))


ben1_variation_bar <- plot_bar_phylo(ben1_all) + 
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        axis.ticks.y = element_blank(),
        legend.background = element_rect(color = "black", size = 0.5, linetype = "solid"))+
  theme(legend.position = "none")

# ben1_variation_bar <- reposition_legend(ben1_variation_bar, 'bottom right')

plot_grid(ben1_variation_bar,
          ben1_tree, labels = "AUTO", ncol = 2, align = 'v', rel_widths = c(1,2))

ggsave(paste0(plot.dir,"F-ALT_Ben1_bar_w_phylo.png"),
       dpi = 300,
       height = 26, 
       width = 10)


# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # 
# BEN-1 phylogeny + WORLD MAP
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
library(ggplot2)  # FYI you need v2.0
library(dplyr)    # yes, i could have not done this and just used 'subset' instead of 'filter'
library(ggalt)    # devtools::install_github("hrbrmstr/ggalt")
library(ggthemes) # theme_map and tableau colors

world <- map_data("world")
world <- world[world$region != "Antarctica",] # intercourse antarctica

phenos <- gwa_mappings %>%
  dplyr::select(strain,value)

isolation_info <-readr::read_tsv("https://elegansvariation.org/strain/strain_data.tsv")%>%
  dplyr::filter(reference_strain == "True")%>%
  dplyr::select(strain = isotype, long = longitude, lat = latitude, landscape, substrate)%>%
  dplyr::filter(lat != "None")%>%
  dplyr::left_join(ben1_strain_markers,.,by="strain")%>%
  dplyr::filter(!is.na(lat)) %>%
  dplyr::left_join(.,phenos,by="strain") %>%
  dplyr::distinct(strain, value, .keep_all = T)%>%
  dplyr::filter(!is.na(value))

isolation_info$lat <- as.numeric(isolation_info$lat)
isolation_info$long <- as.numeric(isolation_info$long)

write.table(isolation_info, 
            file = paste0(data.dir,"isolation_infor_with_variants.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F)

map <- ggplot()+ geom_map(data=world, map=world,
                          aes(x=long, y=lat, map_id=region),
                          color="white", fill="#7f7f7f", size=0.05, alpha=1)+
  geom_point(data=dplyr::filter(isolation_info, !is.na(marker)),  
             aes(x=long, y=lat, fill=ben1_prediction), size = 2, shape = 21)+
  scale_fill_manual(values=colors,name = expression(paste("Variation at ", italic("ben-1"))))+
  theme_map()+ 
  theme(legend.background = element_rect(fill= 'white',
                                         size=0.5, 
                                         linetype="solid", 
                                         colour ="black"),
        legend.text = element_text(colour="black", size=14),
        legend.title = element_text(colour="black", size=16, face = "bold"),
        legend.position = "bottom")

map <- reposition_legend(map,position = 'bottom left')

ben1_tree <- tree %<+% ben1_strain_markers+ 
  geom_tippoint(aes( color=ben1_prediction), size = 2)+
  scale_color_manual(values=colors,name = expression(paste("Variation at ", italic("ben-1")))) + coord_flip() + scale_x_reverse()
ben1_tree

cowplot::plot_grid(map,
                   ben1_tree, labels = "AUTO", ncol = 1, align = 'v', rel_heights = c(1,0.5), label_size = 18)

ggsave(paste0(plot.dir,"F6_map_phylo.png"),
       dpi = 300,
       height = 12, 
       width = 20)

ancestry.colours <- c("gold2", "plum4", "darkorange1", 
                      "lightskyblue2", "firebrick", "burlywood3","gray51", 
                      "springgreen4", "lightpink2","deepskyblue4", "black", 
                      "mediumpurple4", "orange", "maroon","yellow3", "brown4", 
                      "yellow4", "sienna4", "chocolate", "gray19")

ben_by_subs <-readr::read_tsv("https://elegansvariation.org/strain/strain_data.tsv")%>%
  dplyr::filter(reference_strain == "True")%>%
  dplyr::select(strain = isotype, long = longitude, lat = latitude, landscape, substrate)%>%
  dplyr::filter(lat != "None")%>%
  dplyr::left_join(ben1_strain_markers,.,by="strain")%>%
  dplyr::filter(!is.na(lat)) %>%
  dplyr::left_join(.,phenos,by="strain") %>%
  dplyr::distinct(strain, value, .keep_all = T) %>%
  dplyr::mutate(LoF = ifelse(is.na(marker) | strain == "CB4856", "REF","ALT")) %>%
  dplyr::filter(landscape != "None")


ben_by_subs$landscape <- gsub(" zoo","",ben_by_subs$landscape)
ben_by_subs$landscape <- gsub(" ","\n",ben_by_subs$landscape)

sample_loc_plot <- ggplot(ben_by_subs)+
  aes(x = LoF, fill = factor(landscape))+
  geom_bar(position="fill",color = "black")+
  theme_classic(20)+
  labs(x = expression(paste("Variation at ", italic("ben-1"))),
       y = "Frequency")+
  scale_fill_brewer(palette = "Set1",name = "Sampling\nLocation")+
  plot.theme

top_row <- cowplot::plot_grid(map,
                              sample_loc_plot, labels = "AUTO", ncol = 2, align = 'h', rel_widths = c(1,0.5), label_size = 18)

cowplot::plot_grid(top_row,
                   ben1_tree, labels = c("","C"), ncol = 1, align = 'v', rel_heights = c(1,0.5), label_size = 18)


ggsave(paste0(plot.dir,"F6_map_phylo_ALT.png"),
       dpi = 300,
       height = 12, 
       width = 20)


ben_by_subs <-readr::read_tsv("https://elegansvariation.org/strain/strain_data.tsv")%>%
  dplyr::filter(reference_strain == "True")%>%
  dplyr::select(strain = isotype, long = longitude, lat = latitude, landscape, substrate)%>%
  dplyr::filter(lat != "None")%>%
  dplyr::left_join(ben1_strain_markers,.,by="strain")%>%
  dplyr::filter(!is.na(lat)) %>%
  dplyr::left_join(.,phenos,by="strain") %>%
  dplyr::distinct(strain, value, .keep_all = T) %>%
  dplyr::mutate(LoF = ifelse(is.na(marker) | strain == "CB4856", "REF","ALT")) %>%
  dplyr::filter(landscape != "None")


ben_by_subs$landscape <- gsub(" zoo","",ben_by_subs$landscape)
ben_by_subs$landscape <- gsub(" ","\n",ben_by_subs$landscape)

sample_loc_plot <- ggplot(ben_by_subs)+
  aes(x = LoF, fill = factor(landscape))+
  geom_bar(position="fill",color = "black")+
  theme_classic(20)+
  labs(x = expression(paste("Variation at ", italic("ben-1"))),
       y = "Frequency")+
  scale_fill_brewer(palette = "Set1",name = "Sampling\nLocation")+
  plot.theme

substrate_type <- ggplot(ben_by_subs)+
  aes(x = LoF, fill = factor(substrate))+
  geom_bar(position="fill",color = "black")+
  scale_fill_manual(values=ancestry.colours,name = "Sampling\nSubstrate")+
  theme_classic(20)+
  labs(x = expression(paste("Variation at ", italic("ben-1"))),
       y = "Frequency")+
  plot.theme+
  theme(axis.title.y = element_blank())

cowplot::plot_grid(sample_loc_plot,
                                substrate_type, labels = c("A","B"), ncol = 2, align = 'h', rel_widths = c(1,1), label_size = 18)


ggsave(paste0(plot.dir,"F4_sampling_location_substrate.png"),
       dpi = 300,
       height = 6, 
       width = 10)


top_right <- cowplot::plot_grid(sample_loc_plot,
                                substrate_type, labels = c("B","C"), ncol = 2, align = 'h', rel_widths = c(1,1), label_size = 18)

top_row <- cowplot::plot_grid(map,
                              top_right, labels = c("A",""), ncol = 2, align = 'h', rel_widths = c(1,0.75), label_size = 18)

cowplot::plot_grid(top_row,
                   ben1_tree, labels = c("","D"), ncol = 1, align = 'v', rel_heights = c(1,0.5), label_size = 18)

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # 
# F200Y vs ben-1 deletion crispr
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

pr_df <- data.table::fread(file = paste0(data.dir,"F200Y_del_sorter_assay_regressed.tsv"))

N2_means <- pr_df%>%
  dplyr::filter(strain == "N2")%>%
  dplyr::select(genotype, phenotype)%>%
  dplyr::group_by(genotype)%>%
  dplyr::summarise(del_ctrl = mean(phenotype, na.rm = T))

PAM_means <- pr_df%>%
  dplyr::filter(strain %in% c("ECA880", "ECA881"))%>%
  dplyr::select(genotype, phenotype)%>%
  dplyr::group_by(genotype)%>%
  dplyr::summarise(swap_ctrl = mean(phenotype, na.rm = T))


subtract_controls <- pr_df%>%
  dplyr::left_join(.,N2_means, by = "genotype")%>%
  dplyr::left_join(.,PAM_means, by = "genotype")%>%
  dplyr::mutate(subtracted_pheno = ifelse(genotype %in% c("del"), abs(phenotype - del_ctrl),
                                          abs(phenotype - swap_ctrl)))

HTA_plot <- subtract_controls%>%
  dplyr::filter(trait == trait_of_interest, condition %in% condition_of_interest)%>%
  dplyr::mutate(strain1 = factor(group, levels = c("N2","del","PAM","F200Y"),
                                 labels =c("N2","del","PAM","F200Y")))%>%
  ggplot(.) +
  aes(x = factor(strain1), 
      y = subtracted_pheno, 
      fill=group) +
  geom_boxplot(outlier.colour = NA, alpha = 0.7)+
  scale_fill_manual(values = c("blue", "cadetblue3", "hotpink3", "orange"))+
  theme_bw()+
  labs( y = "Animal length")+
  theme(axis.text.x = element_text(size=10, face="bold", color="black"),
        axis.text.y = element_text(size=12, face="bold", color="black"),
        axis.title.x = element_text(size=0, face="bold", color="black", vjust=-.3),
        axis.title.y = element_text(size=16, face="bold", color="black"),
        strip.text.x = element_text(size=24, face="bold", color="black"),
        strip.text.y = element_text(size=16, face="bold", color="black"),
        plot.title = element_text(size=16, face="bold", vjust = 1),
        legend.position="none",
        panel.background = element_rect( color="black",size=1.2),
        strip.background = element_rect(color = "black", size = 1.2),
        panel.border = element_rect( colour = "black"))

##### TukeyHSD ##### 
run_TukeyHSD(subtract_controls,trait_of_interest)

############ Figure 6 B: competition assay plot ########

CA_v2 <- data.table::fread(file = paste0(data.dir,"F200Y_del_competition.csv"))%>%
  dplyr::filter(!strain == "N2")


CA_plot <-ggplot(CA_v2, aes(x=week, y=frac_mean, fill=condition, color=strain))+
  scale_y_continuous(limits = c(30, 100))+
  geom_errorbar(aes(ymin=frac_mean - standdev, ymax=frac_mean + standdev), width=.1, 
                position=position_dodge(0.05))+
  geom_line(aes(linetype=condition))+
  geom_point()+
  scale_color_manual(values = c("blue", "orange"))+
  theme_bw()+
  labs( x = "Generation")+
  labs( y = "Allele Frequency (%)")+
  theme(axis.text.x = element_text(size=10, face="bold", color="black"),
        axis.text.y = element_text(size=12, face="bold", color="black"),
        axis.title.x = element_text(size=16, face="bold", color="black", vjust=-.3),
        axis.title.y = element_text(size=16, face="bold", color="black"),
        strip.text.x = element_text(size=24, face="bold", color="black"),
        strip.text.y = element_text(size=16, face="bold", color="black"),
        plot.title = element_text(size=16, face="bold", vjust = 1),
        legend.position="top",
        panel.background = element_rect( color="black",size=1.2),
        strip.background = element_rect(color = "black", size = 1.2),
        panel.border = element_rect( colour = "black"))

plot_grid( HTA_plot,CA_plot, labels = "AUTO", ncol = 2, align = 'v')

ggsave(paste0(plot.dir,"figure6.png"), 
       dpi = 300,
       height = 6, 
       width = 10)


# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # 
# tubulin tajima d calculations
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

ys <- c(-2.6, 0.5)

# III:3537688..3541628
st = 3537688
en = 3541628
ben1_d <- tajimas_d_temp(vcf_path = paste0(data.dir,"custom_Popgen") ,vcf_name = "WI.20170531.impute.vcf.gz", 
                         chromosome = "III", interval_start = st-1e5, interval_end = en+1e5, site_of_interest = st,
                         slide_distance = 1,window_size = 100, outgroup = "XZ1516")
ben1td <- ben1_d[[2]]+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())+
  ggplot2::geom_vline(ggplot2::aes(xintercept = 3541628/1e+06), color = "red", alpha = 0.7, size = 1)+
  ylim(ys)

# III:10740868..10742461
st = 10740868 
en = 10742461
tbb1_d <- tajimas_d_temp(vcf_path = paste0(data.dir,"custom_Popgen") ,vcf_name = "WI.20170531.impute.vcf.gz", 
                         chromosome = "III", interval_start = st-1e5, interval_end = en+1e5 ,site_of_interest = st,
                         slide_distance = 1,window_size = 100, outgroup = "XZ1516")
tbb1td <- tbb1_d[[2]]+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())+
  ggplot2::geom_vline(ggplot2::aes(xintercept = 10742461/1e+06), color = "red", alpha = 0.7, size = 1)+
  ylim(ys)

# III:4015769..4017643
st = 4015769
en = 4017643

tbb2_d <- tajimas_d_temp(vcf_path = paste0(data.dir,"custom_Popgen") ,vcf_name = "WI.20170531.impute.vcf.gz", 
                         chromosome = "III", interval_start = st-1e5, interval_end = en+1e5 ,site_of_interest = st,
                         slide_distance = 1,window_size = 100, outgroup = "XZ1516") 
tbb2td <- tbb2_d[[2]]+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())+
  ggplot2::geom_vline(ggplot2::aes(xintercept = 4017643/1e+06), color = "red", alpha = 0.7, size = 1)+
  ylim(ys)

# X:9434452..9436893
st = 9434452
en = 9436893
tbb4 <- tajimas_d_temp(vcf_path = paste0(data.dir,"custom_Popgen") ,vcf_name = "WI.20170531.impute.vcf.gz",
                       chromosome = "X", interval_start = st-1e5, interval_end = en+1e5 ,site_of_interest = st,
                       slide_distance = 1,window_size = 100, outgroup = "XZ1516") 
tbb4td <- tbb4[[2]]+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())+
  ggplot2::geom_vline(ggplot2::aes(xintercept = 9436893/1e+06), color = "red", alpha = 0.7, size = 1)+
  ylim(ys)


#V:12261804..12263368
st = 12261804
en = 12263368

tbb6 <- tajimas_d_temp(vcf_path = paste0(data.dir,"custom_Popgen") ,vcf_name = "WI.20170531.impute.vcf.gz",
                       chromosome = "V", interval_start = st-1e5, interval_end = en+1e5 ,site_of_interest = st,
                       slide_distance = 1,window_size = 100, outgroup = "XZ1516") 
tbb6td <- tbb6[[2]]+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())+
  ggplot2::geom_vline(ggplot2::aes(xintercept = 12263368/1e+06), color = "red", alpha = 0.7, size = 1)+
  ylim(ys)


# X:7774859..7776689
st = 7774859
en = 7776689
mec7 <- tajimas_d_temp(vcf_path = paste0(data.dir,"custom_Popgen") ,vcf_name = "WI.20170531.impute.vcf.gz",
                       chromosome = "X", interval_start = st-1e5, interval_end = en+1e5 ,site_of_interest = st,
                       slide_distance = 1,window_size = 100, outgroup = "XZ1516") 
mec7td <- mec7[[2]]+
  theme(axis.title.y = element_blank())+
  ggplot2::geom_vline(ggplot2::aes(xintercept = 7776689/1e+06), color = "red", alpha = 0.7, size = 1)+
  ylim(ys)


# combine all plots
plots <- cowplot::plot_grid(ben1td,tbb1td,tbb2td,tbb4td,tbb6td,mec7td, ncol = 1, labels = "AUTO")
plots

# now add the title
title <- ggdraw() + draw_label("Tajima's D for Beta Tubulins", fontface='bold')
plot_w_title <- plot_grid(title, plots, ncol=1, rel_heights=c(0.05, 1))

ytitle <- ggdraw() + draw_label("Tajima's D", fontface='bold', angle = 90)

plot_grid(ytitle, plot_w_title, ncol=2, rel_widths=c(0.025, 1))

ggsave(paste0(plot.dir,"TajimaDfigure.png"), 
       dpi = 300,
       height = 10, 
       width = 10)

ggsave(paste0(plot.dir,"TajimaDfigure.pdf"), 
       dpi = 300,
       height = 10, 
       width = 10)
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # 
# ben-1 gene tajima d with manually curated variants
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

ben1_snps_pr <- ben1_snps%>%
  dplyr::select(POS, strain, GT, nt_change, effect, impact)%>%
  tidyr::unite(marker, effect,impact,nt_change,POS,sep="_",remove=F)%>%
  dplyr::select(marker, strain, GT,start = POS)%>%
  dplyr::mutate(end = start)%>%
  dplyr::filter(!grepl("MODIF|syn",marker))

# # # COMBINE SNPS AND INDELS
ben1_variants <- dplyr::bind_rows(ben1_snps_pr,pr_indels)

ben1_variants <- dplyr::bind_rows(ben1_snps_pr,pr_indels)%>%
  # dplyr::filter(GT=="ALT")%>%
  dplyr::mutate(GTn = ifelse(GT=="ALT",1,0)) %>%
  dplyr::select(-start,-end, -GT)%>%
  tidyr::spread(marker, GTn)%>%
  dplyr::select(-strain)

ben1_variants[is.na(ben1_variants)] <- 0

non_syn <- gene_level_TajimasD(ben1_variants)
non_syn

ben1_snps_pr <- ben1_snps%>%
  dplyr::select(POS, strain, GT, nt_change, effect, impact)%>%
  tidyr::unite(marker, effect,impact,nt_change,POS,sep="_",remove=F)%>%
  dplyr::select(marker, strain, GT,start = POS)%>%
  dplyr::mutate(end = start)

ben1_variants <- dplyr::bind_rows(ben1_snps_pr,pr_indels)%>%
  # dplyr::filter(GT=="ALT")%>%
  dplyr::mutate(GTn = ifelse(GT=="ALT",1,0)) %>%
  dplyr::select(-start,-end, -GT)%>%
  tidyr::spread(marker, GTn)%>%
  dplyr::select(-strain)

ben1_variants[is.na(ben1_variants)] <- 0

syn <- gene_level_TajimasD(ben1_variants)
syn


ben1_variants <- ben1_snps%>%
  dplyr::select(POS, strain, GT, nt_change, effect, impact)%>%
  tidyr::unite(marker, effect,impact,nt_change,POS,sep="_",remove=F)%>%
  dplyr::select(marker, strain, GT,start = POS)%>%
  dplyr::mutate(end = start)%>%
  dplyr::filter(!grepl("MODIF",marker))%>%
  dplyr::mutate(GTn = ifelse(GT=="ALT",1,0)) %>%
  dplyr::select(-start,-end, -GT)%>%
  tidyr::spread(marker, GTn)%>%
  dplyr::select(-strain)


# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # 
# ben-1 neutrality statistics with phylo and map
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
### TREE
tree <-  ape::read.tree( file = paste0(data.dir, "genome.tree"))

tree <- ggtree(tree, 
               branch.length = "rate",
               layout = "circular") 

ben1_variation_matrix <- gwa_mappings%>%
  dplyr::ungroup()%>%
  dplyr::select(strain,marker)%>%
  # na.omit()%>%
  dplyr::distinct()

missing_strains <- data.frame(strain = row.names(cegwas::kinship)[!row.names(cegwas::kinship)%in%ben1_variation_matrix$strain],
                              marker = NA, stringsAsFactors = F)
colors <- c("#2121D9","#DF0101","#04B404","#FF9326","#A945FF","#0089B2","#B26314","#610B5E","#FE2E9A","#BFF217")

ben1_strain_markers <- ben1_variation_matrix %>%
  dplyr::bind_rows(.,missing_strains)%>%
  dplyr::distinct()%>%
  dplyr::filter(!grepl("missense_variant_MODERATE_c.599T>A_353", marker))%>%
  dplyr::mutate(ben1_prediction = ifelse(is.na(marker), NA, 
                                         ifelse(grepl("del", marker),"Deletion",
                                                ifelse(grepl("ins", marker,ignore.case = T),"Insertion",
                                                       ifelse(grepl("inv", marker),"Inversion",
                                                              ifelse(grepl("stop", marker),"Stop Gained",
                                                                     ifelse(grepl("trans", marker),"Transposon Insertion",
                                                                            ifelse(grepl("missense", marker),"Missense","Splice Donor"))))))))

ben1_tree <- tree %<+% ben1_strain_markers+ 
  geom_tippoint(aes( color=ben1_prediction), size = 2)+
  scale_color_manual(values=colors,name = expression(paste("Variation at ", italic("ben-1")))) 

ben1_tree <- tree %<+% ben1_strain_markers+ 
  geom_tippoint(aes( color=ben1_prediction), size = 2)+
  scale_color_manual(values=colors,name = expression(paste("Variation at ", italic("ben-1")))) + coord_flip() + scale_x_reverse()

# # # MAP 
world <- map_data("world")
world <- world[world$region != "Antarctica",] # intercourse antarctica

isolation_info <- data.table::fread(paste0(data.dir,"isolation_infor_with_variants.tsv")) 

map <- ggplot()+ geom_map(data=world, map=world,
                          aes(x=long, y=lat, map_id=region),
                          color="white", fill="#7f7f7f", size=0.05, alpha=1)+
  geom_point(data=dplyr::filter(isolation_info, !is.na(marker)),  
             aes(x=long, y=lat, fill=ben1_prediction), size = 2, shape = 21)+
  scale_fill_manual(values=colors,name = expression(paste("Variation at ", italic("ben-1"))))+
  theme_map()+ 
  theme(legend.background = element_rect(fill= 'white',
                                         size=0.5, 
                                         linetype="solid", 
                                         colour ="black"),
        legend.text = element_text(colour="black", size=14),
        legend.title = element_text(colour="black", size=16, face = "bold"),
        legend.position = "bottom")


### NEUTRALITY STATS

library(PopGenome)
vcf_name = "WI.20170531.impute.vcf.gz"
slide_distance = 1
window_size = 100

chr1 <- c(1,15072434)
chr2 <- c(1,15279421)
chr3 <- c(1,13783801)
chr4 <- c(1,17493829)
chr5 <- c(1,20924180)
chr6 <- c(1,17718942)
chr7 <- c(1,13794)

chr.lengths <- list(chr1,chr2,chr3,chr4,chr5,chr6,chr7)
chroms <- c("I","II","III","IV","V","X")

w.CHROM <- 3

setwd("~/Dropbox/HTA/Results/2017_anthelmintic_gwa/analysis/ben1_paper/Data/custom_Popgen/")

st = 3537688
en = 3541628

gen <- PopGenome::readVCF(vcf_name, 
                          numcols = 1000, 
                          tid = chroms[w.CHROM], 
                          frompos =  st-1e5, 
                          topos =  en+1e5, 
                          approx = F,
                          gffpath = "protein_coding.gff",
                          include.unknown=TRUE)


# Set most diverged strain as the outgroup
gen2 <- PopGenome::set.outgroup(gen, c("XZ1516"),  diploid = FALSE)

s_gen <- PopGenome::sliding.window.transform(gen2, 
                                             width = window_size, 
                                             jump = slide_distance, 
                                             whole.data = FALSE)

s_gen_stats <- PopGenome::detail.stats(s_gen)

s_gen_stats <- PopGenome::neutrality.stats(s_gen_stats,detail = T)


gene_names <- get.feature.names(object = s_gen_stats,chr = chroms[w.CHROM],gff.file = "protein_coding.gff")

s_gen_stats <- PopGenome::diversity.stats(s_gen_stats,pi = T)

get.sum.data(s_gen_stats)

# . . . . 
# . . . . SAVE GENOME OBJECT
# . . . . 

save(s_gen_stats,file = paste0(data.dir,"20180602ben1_popgenome.Rda"))

# . . . . 
# . . . . DEFINE SLIDE WINDOWS
# . . . . 

load(paste0(data.dir,"20180602ben1_popgenome.Rda"))

windowStarts <- data.frame(snp_index = 1:length(colnames(s_gen_stats@BIG.BIAL[[1]])),
                           position = colnames(s_gen_stats@BIG.BIAL[[1]]))

slide_index <- cbind(data.frame(lapply(s_gen_stats@SLIDE.POS,  function(x) as.numeric(floor(mean(x))))))%>%
  tidyr::gather(temp, snp_index)%>%
  dplyr::select(-temp)%>%
  dplyr::left_join(., windowStarts, by = "snp_index")

# . . . . 
# . . . . ANALYZYE NEUTRALITY/DIVERSITY STATS
# . . . .


neutrality_ls <- list()
for(popnumber in 1){
  popname <- LETTERS[[popnumber]]
  neutrality_ls[[popnumber]] <- data.frame(PopGenome::get.neutrality(s_gen_stats, theta = T, stats = T)[[popnumber]],
                                           PopGenome::get.diversity(s_gen_stats)[[popnumber]],
                                           PI = s_gen_stats@Pi[,popnumber])%>%
    dplyr::mutate(Population = popname,
                  WindowPosition = slide_index$position)%>%
    dplyr::select(-Pi)
}

neutrality_df <- dplyr::bind_rows(neutrality_ls)
neutrality_df <- neutrality_df[,colSums(is.na(neutrality_df))<nrow(neutrality_df)]
neutrality_df <- tidyr::gather(neutrality_df, statistic, value, -Population, -WindowPosition)


plt_df <- neutrality_df%>%
  dplyr::filter(statistic %in%c("Tajima.D","Fay.Wu.H","Zeng.E"))%>%
  dplyr::mutate(f_stats = ifelse(statistic == "Tajima.D", "Tajima's D",
                                 ifelse(statistic == "Fay.Wu.H", "Fay and Wu's H","Zeng's E")))
statspt <- plt_df%>%
  ggplot()+
  aes(x = WindowPosition/1e6, y = as.numeric(value), fill = f_stats)+
  scale_color_manual(values = "blue")+
  geom_point(shape = 21, size =2)+
  scale_fill_manual(values=c("hotpink3","black","cadetblue3"), name = "")+
  theme_bw(18)+
  theme(axis.text.x = ggplot2::element_text(size = 16),
        axis.text.y = ggplot2::element_text(size = 16),
        axis.title.x = ggplot2::element_text(size = 18, face = "bold", color = "black", vjust = -0.3),
        axis.title.y = ggplot2::element_blank(), strip.text.y = element_text(size = 14, face = "bold"))+
  labs(x = "Genomic Position (Mb)")+
  geom_vline(aes(xintercept=c(3.541595)), linetype="dotdash", alpha = 0.5, color = "black")+
  geom_vline(aes(xintercept=c(3.538293)), linetype="dotdash", alpha = 0.5, color = "black") +
  theme(strip.background = element_blank(),legend.position = "top")


topplot <- cowplot::plot_grid(statspt,map,labels = c("A","B"), rel_widths = c(0.5,1),label_size = 18)

bottomplot <- cowplot::plot_grid(ben1_tree, labels = c("C"), label_size = 18)

cowplot::plot_grid(topplot,bottomplot, nrow = 2)


ggsave(paste0(plot.dir,"neut_map_phylo.png"),
       height = 10,
       width  = 14,
       dpi = 300)



