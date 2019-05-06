# Finalized figures for BEN-1 paper

# calculate ka/ks with http://services.cbu.uib.no/tools/kaks
# use "seqdump.txt" in data/custom_Popgen folder


# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # 
# set to location of files
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # 
main.dir <- "~/Dropbox/HTA/Results/2017_anthelmintic_gwa/analysis/ben1_paper/"
data.dir <- paste0(main.dir,"Raw_data/")
plot.dir <- paste0(main.dir,"Final_Figures/TIFFs/")
script.dir <- paste0(main.dir,"Scripts/")
final.dir <- paste0(main.dir,"Final_Tables/")


source(paste0(script.dir,"ben1_processing_functions.R"))

trait_of_interest <- "q90.TOF"
condition_of_interest <- "albendazole"
control_of_interest <- "DMSO"

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # 
# FIGURE S1 - PLOT DOSE RESPONSES 
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

dr_data <- data.table::fread(paste0(final.dir,"TS2_DR_Processed.tsv"))

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
       y = "Brood Size") + 
  plot.theme + 
  theme(legend.position = "none")

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
        axis.text.x = element_blank(),
        legend.position = "top", 
        legend.background = element_rect(color = "black", size = 0.5, linetype = "solid"),
        legend.text = element_text(size =10),
        legend.direction = "horizontal",
        legend.justification = "right")

plot_grid( size_DR,brood_DR, labels = "AUTO", ncol = 1, align = 'v', label_size = 12)

ggsave(paste0(plot.dir,"S1_fig.pdf"), 
       dpi = 300,
       height = 4.5, 
       width = 7.5)

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # 
#  FIGURE 1 -PLOT GWAS and BURDEN
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

pr_maps <- data.table::fread(paste0(final.dir,"TS5_GWA_processed_marker_mapping.tsv"))%>%
  dplyr::filter(trait != "norm.n")

snv_manplot <- manplot_edit(pr_maps)[[1]]+
  plot.theme + 
  theme(axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank())

q90burden_pr <- data.table::fread(paste0(final.dir, "TS8_burden_test_mappings.tsv"))

burden_manplot <- q90burden_pr%>%
  ggplot()+
  aes(x = POS/1e6, y = Stat, size = NumVar, alpha = 0.5, color = significant)+
  geom_point(size = 0.25)+
  scale_color_manual(values=c("black","red"))+
  facet_grid(.~CHROM, scales = "free", space = "free")+
  theme_bw()+
  plot.theme+
  ggplot2::theme(legend.position = "none",
                 strip.background = element_blank(),
                 strip.text = element_blank())+
  labs(x = "Genomic Position (Mb)", y = "Test Statisitic")

plot_grid(snv_manplot, burden_manplot, labels = "AUTO", ncol = 1, align = 'v', label_size = 11, rel_heights = c(1,1))

ggsave(paste0(plot.dir,"Fig1.tiff"), 
       dpi = 300,
       height = 3.5, 
       width = 7, compression = "lzw")

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # 
#  FIGURE S2 - PLOT GWAS and BURDEN
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

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

plot_grid(animal_length_pxg, bottom_row, labels = "AUTO", ncol = 1, align = 'v', label_size = 18)

ggsave(paste0(plot.dir,"S2_fig.pdf"), 
       dpi = 300,
       height = 6, 
       width = 8)

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # 
#  FIGURE 2A - BAR PLOT OF BEN-1 VARIATION
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

gwa_mappings <- data.table::fread(file = paste0(final.dir,"TS5_GWA_processed_marker_mapping.tsv"))%>%
  dplyr::filter(trait!="norm.n")%>%
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
  dplyr::arrange(marker)%>%
  dplyr::distinct(strain, value, .keep_all=T)%>%
  dplyr::arrange(value)%>%
  dplyr::mutate(strain2 = factor(strain, levels = unique(strain), labels = unique(strain), ordered = T))%>%
  dplyr::mutate(norm_pheno_temp = ifelse(value == min(value), 0, 1))%>%
  dplyr::mutate(delta_pheno = ifelse(norm_pheno_temp == 0, 0, abs(dplyr::lag(value) - value)))%>%
  dplyr::mutate(norm_pheno = cumsum(delta_pheno)) %>%
  dplyr::mutate(final_pheno = norm_pheno/max(norm_pheno))

ben1_all <- gwa_mappings_bar%>%
  dplyr::mutate(build_GT = ifelse(ben1_prediction %in% c("Missense","Stop Gained", "Splice Donor","Deletion", "Insertion", "Inversion", "Transposon Insertion"), ben1_prediction, "A"))

ben1_all$build_GT[is.na(ben1_all$build_GT)] <- "A"

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


ggsave(plot = reposition_legend(ben1_variation_bar, 'top left'),
       paste0(plot.dir,"F2A_ben1_resistance_barplot",condition_of_interest,"_",trait_of_interest,".png"), 
       dpi = 300,
       height = 4, 
       width = 10)
ggsave(plot = reposition_legend(ben1_variation_bar, 'top left'),
       paste0(plot.dir,"F2A_ben1_resistance_barplot",condition_of_interest,"_",trait_of_interest,".pdf"), 
       dpi = 300,
       height = 4, 
       width = 10)


# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # 
# FIGURE 2B - BEN-1 HAPLOTYPE GENE MODEL PLOT
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

# THINGS TO MANUALLY ADJUST IN AN SVG SOFTWARE
# COLORS OF VARIANTS
# GROUP THESE STRAINS TOGETHER - ED3005 JU1808 MY18 MY795 - BC MY795 HAS A VARIANT UNRELIABLE CALLED

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


for_plot <- data.table::fread(file = paste0(final.dir,"TS5_GWA_processed_marker_mapping.tsv"))%>%
  dplyr::filter(trait!="norm.n")%>%
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

ggsave(paste0(plot.dir,"F2B_Ben1_haps_by_phenotype.png"),
       dpi = 300,
       height = 12, 
       width = 18)

ggsave(paste0(plot.dir,"F2B_Ben1_haps_by_phenotype.pdf"),
       dpi = 300,
       height = 12, 
       width = 18)

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # 
# FIGURE S3 - Generated CHIMERA
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # 
# TEXT INFORMATION -  ben-1 gene tajima d with manually curated variants
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

# # # COMBINE SNPS AND INDELS - GET THIS INFORMATION FROM PREVIOUS SECTION

ben1_snps_pr <- ben1_snps%>%
  dplyr::select(POS, strain, GT, nt_change, effect, impact)%>%
  tidyr::unite(marker, effect,impact,nt_change,POS,sep="_",remove=F)%>%
  dplyr::select(marker, strain, GT,start = POS)%>%
  dplyr::mutate(end = start)%>%
  dplyr::filter(!grepl("MODIF|syn",marker)) #  REMOVE MODIFIER (INTRON) and SYNONYMOUS VARIANTS

ben1_variants <- dplyr::bind_rows(ben1_snps_pr,pr_indels)%>%
  dplyr::mutate(GTn = ifelse(GT=="ALT",1,0)) %>%
  dplyr::select(-start,-end, -GT)%>%
  tidyr::spread(marker, GTn)%>%
  dplyr::select(-strain)

ben1_variants[is.na(ben1_variants)] <- 0

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # For coding variation
non_syn <- gene_level_TajimasD(ben1_variants)
non_syn
# Wattersons_Theta Average_PWD TajimasD Singletons FuLiF
# 1             5.25        0.39    -2.59         26 -8.09

ben1_snps_pr <- ben1_snps%>%
  dplyr::select(POS, strain, GT, nt_change, effect, impact)%>%
  tidyr::unite(marker, effect,impact,nt_change,POS,sep="_",remove=F)%>%
  dplyr::select(marker, strain, GT,start = POS)%>%
  dplyr::mutate(end = start)

ben1_variants <- dplyr::bind_rows(ben1_snps_pr,pr_indels)%>%
  dplyr::mutate(GTn = ifelse(GT=="ALT",1,0)) %>%
  dplyr::select(-start,-end, -GT)%>%
  tidyr::spread(marker, GTn)%>%
  dplyr::select(-strain)

ben1_variants[is.na(ben1_variants)] <- 0

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # For coding + MODIFIER + SYNONYMOUS variation
syn <- gene_level_TajimasD(ben1_variants)
syn
# Wattersons_Theta Average_PWD TajimasD Singletons FuLiF
# 1            12.64        4.15    -2.02         42 -6.07


# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
# FIGURE 3 - ben-1 neutrality statistics with phylo and map
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
### TREE
tree <-  ape::read.tree( file = paste0(final.dir, "genome.tree"))

tree <- ggtree(tree, 
               branch.length = "rate", size = .2) 


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
                                                                     ifelse(grepl("trans", marker),"Transposon\nInsertion",
                                                                            ifelse(grepl("missense", marker),"Missense","Splice Donor"))))))))

ben1_tree <- tree %<+% na.omit(ben1_strain_markers)+ 
  geom_tippoint(aes( color=ben1_prediction), size = 0.3)+
  scale_color_manual(values=colors,
                     name = expression(paste(italic("ben-1"), " variation")),
                     labels = c("Deletion", "Insertion", "Inversion", "Missense",
                                "Splice Donor", "Stop Gained", "Transposon\nInsertion", ""),
                     guide = guide_legend(title.position = "top", 
                                          ncol=2))+
  coord_flip() + 
  scale_x_reverse()+
  theme(legend.background = element_rect(fill= 'white',
                                         size=0.25, 
                                         linetype="solid", 
                                         colour ="black"),
        legend.text = element_text(colour="black", size=8),
        legend.title = element_text(colour="black", size=10, face = "bold"),
        legend.justification = c(1, 0), 
        legend.position = c(0.3, .1),
        legend.box = "horizontal",
        legend.key.size = unit(.1, "cm"))

# # # MAP 
world <- map_data("world")
world <- world[world$region != "Antarctica",] # intercourse antarctica

isolation_info <- readr::read_tsv("https://elegansvariation.org/strain/strain_data.tsv")%>%
  dplyr::left_join(ben1_strain_markers,.,by="strain")%>%
  dplyr::select(strain, marker, ben1_prediction, long=longitude, lat=latitude, landscape, substrate)

isolation_info$long <- as.numeric(isolation_info$long)
isolation_info$lat <- as.numeric(isolation_info$lat)

map <-ggplot()+ geom_map(data=world, map=world,
                         aes(x=long, y=lat, map_id=region),
                         color="white", fill="#7f7f7f", 
                         size=0.05, alpha=1)+
  geom_point(data=dplyr::filter(isolation_info, !is.na(marker)),  
             aes(x=long, y=lat, color = ben1_prediction), 
             size = 0.5)+
  scale_color_manual(values=colors,
                     name = expression(paste(italic("ben-1"), " variation")),
                     guide = guide_legend(title.position = "top", 
                                          ncol=2))+
  theme_map()+ 
  theme(legend.position = "none")


### NEUTRALITY STATS

# . . . . 
# . . . . DEFINE SLIDE WINDOWS
# . . . . 

load(paste0(final.dir,"TS10_ben1_popgen.Rda"))

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
                                 ifelse(statistic == "Fay.Wu.H", "Fay & Wu's H","Zeng's E")))
statspt <- plt_df%>%
  ggplot()+
  aes(x = WindowPosition/1e6, 
      y = as.numeric(value), 
      color = f_stats)+
  geom_line(size = 0.25)+
  scale_color_manual(values=c("hotpink3","black","cadetblue3"), name = "",
                     guide = guide_legend(label.hjust = 0,
                                          label.position = "left",
                                          reverse = T,
                                          override.aes = list(shape = 15,
                                                              size = 1)))+
  theme_bw()+
  plot.theme +
  theme(axis.title.y = ggplot2::element_blank(),
        legend.text = element_text(colour="black", size=8),
        legend.background = element_rect(color = "black", size = 0.25))+
  labs(x = "Genomic Position (Mb)")+
  geom_vline(aes(xintercept=c(3.541595)),
             linetype="dotdash", alpha = 0.5, color = "black")+
  geom_vline(aes(xintercept=c(3.538293)), 
             linetype="dotdash", alpha = 0.5, color = "black") +
  theme(legend.justification = c(0, 0), 
        legend.position = c(0.01, 0.01),
        legend.title = element_blank(),
        legend.key.height = unit(0.01,"cm"),
        legend.key.width = unit(0.01,"cm"))

topplot <- cowplot::plot_grid(statspt,
                              map,
                              labels = c("A","B"), 
                              rel_widths = c(0.6,1),label_size = 11)

bottomplot <- cowplot::plot_grid(ben1_tree, labels = c("C"), label_size = 11)

cowplot::plot_grid(topplot,bottomplot, 
                   rel_heights = c(1,.7),
                   nrow = 2)



ggsave(paste0(plot.dir,"Fig3.pdf"),
       height = 4.645,
       width  = 6.5,
       dpi = 300)

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
# FIGURE S4 - Sampling location and statistics
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

phenos <- gwa_mappings %>%
  dplyr::select(strain,value)

ben_by_subs <-readr::read_tsv("https://elegansvariation.org/strain/strain_data.tsv")%>%
  dplyr::filter(reference_strain == "True",
                release != "20180413")%>%
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
ben_by_subs$substrate <- gsub(" ","\n",ben_by_subs$substrate)

sample_loc_plot <- ggplot(ben_by_subs)+
  aes(x = LoF, fill = factor(landscape))+
  geom_bar(position="fill",color = "black")+
  theme_classic()+
  labs(x = expression(paste("Variation at ", italic("ben-1"))),
       y = "Frequency")+
  scale_fill_brewer(palette = "Set1",name = "Sampling\nLocation")+
  plot.theme

substrate_type <- ggplot(ben_by_subs)+
  aes(x = LoF, fill = factor(substrate))+
  geom_bar(position="fill",color = "black")+
  scale_fill_manual(values=ancestry.colours,
                    name = "Sampling\nSubstrate",
                    guide = guide_legend(ncol=2))+
  theme_classic()+
  labs(x = expression(paste("Variation at ", italic("ben-1"))),
       y = "Frequency")+
  plot.theme+
  theme(axis.title.y = element_blank())

cowplot::plot_grid(sample_loc_plot,
                   substrate_type, labels = c("A","B"), 
                   ncol = 2, 
                   align = 'h', 
                   rel_widths = c(1,1), 
                   label_size = 11)

ggsave(paste0(plot.dir,"S4_fig.pdf"),
       dpi = 300,
       height = 4.5, 
       width = 7.5)

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
# FIGURE 4 - HTA and competition assay of F200Y and DEL alleles
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

pr_df <- data.table::fread(file = paste0(final.dir,"TS13_HTA_ben-1_regressed.tsv"))

main_figure <- pr_df%>%
  dplyr::filter(strain %in% c("N2", "ECA883", "ECA919"))
# 
# N2_means <- pr_df%>%
#   dplyr::filter(strain == "N2")%>%
#   dplyr::select(genotype, phenotype)%>%
#   dplyr::group_by(genotype)%>%
#   dplyr::summarise(del_ctrl = mean(phenotype, na.rm = T))
# 
# PAM_means <- pr_df%>%
#   dplyr::filter(strain %in% c("ECA880", "ECA881"))%>%
#   dplyr::select(genotype, phenotype)%>%
#   dplyr::group_by(genotype)%>%
#   dplyr::summarise(swap_ctrl = mean(phenotype, na.rm = T))
# 
# subtract_controls <- pr_df%>%
#   dplyr::left_join(.,N2_means, by = "genotype")%>%
#   dplyr::left_join(.,PAM_means, by = "genotype")%>%
#   dplyr::mutate(subtracted_pheno = ifelse(genotype %in% c("del"), (phenotype - del_ctrl),
#                                           (phenotype - swap_ctrl)))

HTA_plot <- main_figure%>%
  dplyr::filter(trait == trait_of_interest, condition == "Albendazole")%>%
  dplyr::mutate(strain1 = factor(group, levels = c("N2","del","F200Y"),
                                 labels =c("N2","Del","F200Y")))%>%
  ggplot(.) +
  aes(x = factor(strain1), 
      y = phenotype, 
      fill=group) +
  geom_boxplot(outlier.colour = NA, alpha = 0.7)+
  scale_fill_manual(values = c("cadetblue2","hotpink2","grey75"))+
  theme_bw()+
  labs( y = "Animal length")+
  plot.theme +
  theme(axis.title.x = element_blank(),
        legend.position="none")

##### TukeyHSD ##### 
run_TukeyHSD(main_figure,trait_of_interest)

############ Figure 4 B: competition assay plot ########

competition_assay <- data.table::fread(file = paste0(final.dir,"TS14_competition_assay.csv"))%>%
  dplyr::filter(TargetType == "Ch1Unknown")%>%
  dplyr::select(Condition, Replicate, Generation, FractionalAbundance)%>%
  dplyr::mutate(id = paste(Condition,Generation, sep = "-")) %>%
  na.omit()%>% 
  dplyr::group_by(Condition, Generation)%>% 
  dplyr::mutate(Mean =  mean(100-FractionalAbundance)) %>% 
  dplyr::mutate(SND =  sd(100-FractionalAbundance))%>%
  dplyr::ungroup()%>%
  dplyr::select(id, FractionalAbundance)

outliers  <-  competition_assay%>%
  dplyr::group_by(id)%>%
  dplyr::transmute_if(is.numeric, isnt_out_funs)

competition_assay_outliers <- dplyr::bind_cols(competition_assay,outliers)%>%
  tidyr::separate(id, c("Condition", "Generation"), sep = "-")%>%
  # dplyr::mutate(outlier = ifelse(!z | !mad | !tukey, "OUTLIER", "OK"))%>%
  dplyr::ungroup() %>%
  dplyr::filter(!(!z | !mad | !tukey))%>%
  dplyr::group_by(Condition, Generation)%>%
  dplyr::mutate(Mean =  mean(100-FractionalAbundance)) %>% 
  dplyr::mutate(SND =  sd(100-FractionalAbundance))%>%
  tidyr::separate(Condition, into = c("Strain", "Condition"))%>% 
  dplyr::filter(!Strain == "N2")

CA_plot <-competition_assay_outliers%>%
  dplyr::ungroup()%>%
  dplyr::mutate(id = paste0(Strain,Condition))%>%
  ggplot(aes(x=Generation, y=Mean, fill=Condition, color=Strain))+
  scale_y_continuous(limits = c(40, 100))+
  geom_errorbar(aes(ymin=Mean - SND, ymax=Mean + SND), 
                width=.1, 
                position=position_dodge(0.05))+
  geom_line(aes(linetype=Condition, group=id))+
  geom_point()+
  scale_color_manual(values = c("cadetblue3", "hotpink3"),
                     guide= guide_legend(nrow=2,
                                         keyheight = 1))+
  scale_linetype_manual(values = c(1, 2),
                        guide= guide_legend(nrow=2, 
                                            keywidth = 2,
                                            keyheight = 1))+
  theme_bw()+
  theme(legend.title = element_text(colour="black", size = 10, face = "bold"))+
  theme(legend.text = element_text(colour="black", size = 8))+
  labs( x = "Generation")+
  labs( y = expression(bold(paste("Allele Frequency of ", bolditalic("ben-1"), " Edits (%)"))))+
  plot.theme+
  theme(legend.position="top")


plot_grid( HTA_plot,CA_plot, labels = "AUTO", ncol = 2, align = 'v', label_size = 11)

ggsave(paste0(plot.dir,"Fig4.tiff"), 
       dpi = 300,
       height = 3.75, 
       width = 7.5)

ggsave(paste0(plot.dir,"Fig4.pdf"), 
       dpi = 300,
       height = 3.75, 
       width = 7.5)

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
# FIGURE S5 - HTA and competition assay of F200Y and DEL alleles - all strains
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

pr_df <- data.table::fread(file = paste0(final.dir,"TS13_HTA_ben-1_regressed.tsv"))


pr_df%>%
  dplyr::filter(trait == trait_of_interest, condition == "Albendazole")%>%
  ggplot(.) +
  aes(x = factor(strain, levels = c("N2","ECA882","ECA883","ECA884","ECA917","ECA918","ECA919","ECA920","ECA921" )), 
      y = phenotype, 
      fill=group) +
  geom_boxplot(outlier.colour = NA, alpha = 0.7)+
  scale_fill_manual(values = c("cadetblue2","hotpink2","grey75","grey75"))+
  theme_bw()+
  labs( y = "Animal length")+
  plot.theme +
  theme(axis.title.x = element_blank(),
        legend.position="none")

ggsave(paste0(plot.dir,"S5_fig.pdf"), 
       dpi = 300,
       height = 3.75, 
       width = 7.5)

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # 
# FIGURE 5 - CHRX manhattan plot and fine mapping
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

pr_resid_maps <- data.table::fread( paste0(final.dir,"TS17_GWA_marker_mappings_ben1_regressed.tsv"))
manhattan_plots <- manplot_edit(pr_resid_maps)

ben1_resid_manplot <- manhattan_plots[[1]] +
  plot.theme + 
  theme(legend.position = "none")

ben1_resid_genes <- data.table::fread( paste0(final.dir,"TS18_GWA_marker_fine_mappings_ben1_regressed.tsv"))

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
  geom_point(size =1)+
  ggplot2::labs(x = "Genomic Position (Mb)",
                y = expression(bold(-log[10](bolditalic(p)))))+
  scale_color_brewer(palette="Set1", name = "Effect:")+
  scale_fill_brewer(palette="Set1", name = "Effect:")+
  scale_shape_manual(values = c(16:25), name = "Effect:")+
  theme_bw()+
  plot.theme+
  theme(axis.title.y = element_blank(),
  legend.text = element_text(colour="black", size=8),
  legend.title = element_text(colour="black", size=10, face = "bold"),
  legend.key.size = unit(.1, "cm"))

plot_grid(ben1_resid_manplot, fine_map, labels = c("A","B"), 
          label_size = 11, ncol = 2, align = 'v', rel_widths = c(1.5,1))

ggsave(paste0(plot.dir,"Fig5.tiff"), 
       dpi = 300,
       height = 2.5, 
       width = 7.5)
ggsave(paste0(plot.dir,"Fig5.pdf"), 
       dpi = 300,
       height = 2.5, 
       width = 7.5)

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # 
# FIGURE S6 - CHRX burden
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

benreg_burden <- read.table(paste0(final.dir, "TS21_GWA_ben1_regressed.VariableThresholdPrice_processed.tsv"),header = T)

benreg_burden%>%
  ggplot()+
  aes(x = POS/1e6, y = Stat, size = NumVar, alpha = 0.5, color = significant)+
  geom_point(size = 0.3)+
  scale_color_manual(values=c("black","red"))+
  facet_grid(.~CHROM, scales = "free", space = "free")+
  theme_bw()+
  plot.theme +
  ggplot2::theme(legend.position = "none",
                 strip.background = element_blank())+
  labs(x = "Genomic Position (Mb)", y = "Test Statisitic")

ggsave(paste0(plot.dir,"S6_fig.pdf"), 
       height = 2.5, 
       width = 7.5)


# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # 
# FIGURE S7 - single marker PxG, colored by ben-1 variants
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

gwa_mappings <- data.table::fread(paste0(final.dir,"TS5_GWA_processed_marker_mapping.tsv"))%>%
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

gwa_mappings <- gwa_mappings %>%
  tidyr::separate(snpMarker, into = c("CHROM","POS"),sep = ":",remove=F,convert=T)%>%
  dplyr::arrange(CHROM,POS)%>%
  dplyr::mutate(marker2 = factor(snpMarker, levels = unique(snpMarker), labels = unique(snpMarker)))%>%
  dplyr::distinct(marker2, value, .keep_all=T)

ggplot(gwa_mappings)+
  aes(x = snpGT, y = value)+
  geom_boxplot(outlier.colour = NA)+
  facet_grid(.~marker2)+
  geom_jitter(shape = 21, color= "black", fill = "gray90", 
              alpha = 0.5, data = dplyr::filter(gwa_mappings, is.na(marker)),
              size = 0.3)+
  geom_jitter(aes(fill=ben1_prediction), data = na.omit(gwa_mappings), 
              shape = 21, color= "black", size = 1)+
  scale_fill_manual(values=colors,name = expression(paste("Variation at ", italic("ben-1"))))+
  theme_bw()+
  ggplot2::theme(panel.background = ggplot2::element_rect(color = "black",size = 1.2),
                 strip.background = element_blank(),
                 panel.border = element_rect(size = 0.1))+
  labs(x = "SNV Genotype at QTL", y = paste0("Animal Length (", trait_of_interest,")"))

ggsave(paste0(plot.dir,"S7_fig.pdf"), 
       height = 2.5, 
       width = 7.5)


# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # 
# FIGURE S8 - CHRX PxG - LD
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
pr_resid_maps <- data.table::fread( paste0(final.dir,"TS17_GWA_marker_mappings_ben1_regressed.tsv"))
regress_ben1 <- data.table::fread(paste0(final.dir,"TS16_ben1_regressed_with_ben1_covariate.tsv"))

x_pg_df <- pr_resid_maps %>%
  dplyr::filter(trait == "regressed_trait")%>%
  na.omit()%>%
  dplyr::select(qtl_marker=marker, strain, reg_val = value, qtl_allele = allele )%>%
  dplyr::left_join(regress_ben1,.,by ="strain")%>%
  dplyr::mutate(qtl_marker1 = factor(qtl_marker, levels = unique(qtl_marker)))

x_pg_df%>%
  ggplot()+
  aes(x = factor(qtl_allele, levels = c(-1,1),
                 labels = c("REF", "ALT")), y = reg_val)+
  geom_boxplot(outlier.colour = NA)+
  geom_jitter(width = 0.25, aes(fill = factor(ben1_variant_trait)), shape = 21, size =2)+
  scale_fill_manual(values = c("cadetblue3","hotpink3"), labels = c("REF", "ALT"), name=expression(italic(ben-1))) +
  theme_bw()+
  plot.theme +
  theme(strip.background = element_blank(),
        panel.border = element_rect(size = 0.3))+
  facet_grid(.~qtl_marker1)+
  labs(x="QTL allele", y = "Regressed Animal Length")+
  theme(axis.title.x = ggplot2::element_text(size = 12, face = "bold", color = "black", vjust = -0.3), 
        axis.title.y = ggplot2::element_text(size = 12, face = "bold", color = "black"))


ggsave(paste0(plot.dir,"S8_fig.pdf"), 
       dpi = 300,
       height = 7, 
       width = 7)


# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # 
# FIGURE S9 - phenotyping WN2002 
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

# read wild strain HTA data 


assayregressed_WI <- data.table::fread(paste0(final.dir,"TS23_WI_rephenotype_HTA_processed.tsv"))

assayregressed_WI%>%
  dplyr::filter(trait == trait_of_interest, condition == "Albendazole")%>%
  ggplot(.) +
  aes(x = factor(strain, levels = c("N2", "JU2141", "WN2002","JU2581")), 
      y = phenotype, 
      fill=strain) +
  geom_boxplot(outlier.colour = NA, alpha = 0.7)+
  theme_bw()+
  labs( y = "Animal length")+
  scale_fill_manual(values=c("hotpink3","cadetblue3","orange","gray75"),name="Strain")+
  plot.theme +
  theme(axis.title.x = element_blank(),
        legend.position="none")

ggsave(paste0(plot.dir,"S9_fig.pdf"), 
       dpi = 300,
       height = 4.5, 
       width = 7.5)


# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # 
# FIGURE S10 - Generated in Jalview
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #


# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # 
# FIGURE S11 - Tajima's D for tubulins
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
ys <- c(-2.6, 0.5)
# III:3537688..3541628
st = 3537688
en = 3541628
ben1_d <- tajimas_d_temp(vcf_path = paste0(final.dir) ,vcf_name = "WI.20170531.impute.vcf.gz", 
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
tbb1_d <- tajimas_d_temp(vcf_path =paste0(final.dir) ,vcf_name = "WI.20170531.impute.vcf.gz", 
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

tbb2_d <- tajimas_d_temp(vcf_path = paste0(final.dir),vcf_name = "WI.20170531.impute.vcf.gz", 
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
tbb4 <- tajimas_d_temp(vcf_path = paste0(final.dir) ,vcf_name = "WI.20170531.impute.vcf.gz",
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

tbb6 <- tajimas_d_temp(vcf_path = paste0(final.dir),vcf_name = "WI.20170531.impute.vcf.gz",
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
mec7 <- tajimas_d_temp(vcf_path = paste0(final.dir) ,vcf_name = "WI.20170531.impute.vcf.gz",
                       chromosome = "X", interval_start = st-1e5, interval_end = en+1e5 ,site_of_interest = st,
                       slide_distance = 1,window_size = 100, outgroup = "XZ1516") 
mec7td <- mec7[[2]]+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size =16))+
  ggplot2::geom_vline(ggplot2::aes(xintercept = 7776689/1e+06), color = "red", alpha = 0.7, size = 1)+
  ylim(ys)

# combine all plots
plots <- cowplot::plot_grid(ben1td,tbb1td,tbb2td,tbb4td,tbb6td,mec7td, ncol = 1, labels = "AUTO")
plots

# now add the title
title <- ggdraw() + draw_label("Tajima's D for beta-tubulins", fontface='bold', size = 18)
plot_w_title <- plot_grid(title, plots, ncol=1, rel_heights=c(0.05, 1))

ytitle <- ggdraw() + draw_label("Tajima's D", fontface='bold', angle = 90, size = 16)

plot_grid(ytitle, plot_w_title, ncol=2, rel_widths=c(0.025, 1))

ggsave(paste0(plot.dir,"S11_fig.pdf"), 
       dpi = 300,
       height = 7, 
       width = 7)

ggsave(paste0(plot.dir,"TajimaDfigure.pdf"), 
       dpi = 300,
       height = 10, 
       width = 10)

