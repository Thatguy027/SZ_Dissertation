try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)))

# data folders
albendazole_data <- "Ben1/"

# source base theme and colors for plots
source("Scripts/Figure_Themes.R")

# source scripts 
source("Scripts/PD_Talk_Functions.R")

trait_of_interest <- "q90.TOF"
condition_of_interest <- "albendazole"
control_of_interest <- "DMSO"


######################################################################################################################## Dose Response

dr_data <- data.table::fread(glue::glue("{albendazole_data}Final_Tables/TS2_DR_Processed.tsv"))

dr_data$strain <- gsub("Ju", "JU", dr_data$strain)

# length
pr_df <- dr_data%>%
  dplyr::ungroup()%>%
  dplyr::filter(trait == trait_of_interest, grepl("alb|DMSO",condition,ignore.case = T))%>%
  dplyr::filter(!is.na(condition)) %>%
  dplyr::arrange(value) %>%
  dplyr::mutate(strain2 = factor(strain, levels = unique(strain), labels = unique(strain), ordered = T))%>%
  dplyr::mutate(norm_pheno_temp = ifelse(value == min(value), 0, 1))%>%
  dplyr::mutate(delta_pheno = ifelse(norm_pheno_temp == 0, 0, abs(dplyr::lag(value) - value)))%>%
  dplyr::mutate(norm_pheno = cumsum(delta_pheno)) %>%
  dplyr::mutate(final_pheno = norm_pheno/max(norm_pheno))

pr_df%>%
  ggplot()+
  aes(x = factor(condition, 
                 levels = c("DMSO", "albendazole3125", "albendazole625", 
                            "albendazole125",  "albendazole25"),
                 labels = c("0", "3.125", "6.25", "12.5", "25")), 
      y = final_pheno, 
      fill = strain)+
  geom_boxplot(outlier.colour = NA)+
  scale_fill_manual(values = c("blue", "cadetblue3", "hotpink3", "orange"), name = "Strain")+
  theme_bw()+
  labs(x = "Albendazole (µM)",
       y = "Relative Animal Length") + 
  base_theme+
  theme(panel.grid.major = element_blank(),
        axis.line = element_line(colour = axis_color))


ggsave(filename = "Plots/Albendazole_DoseResponse.png", height = 4.5, width = 9, dpi = 400)
ggsave(filename = "Plots/Albendazole_DoseResponse.pdf", height = 4.5, width = 9, dpi = 400)


######################################################################################################################## GWAS

pr_maps <- data.table::fread(glue::glue("{albendazole_data}Final_Tables/Albendazole_q90.TOF_processed_mapping.tsv"))

independent_tests <- 1159.75

cegwas2_manplot(plot_df = pr_maps, eigen_cutoff = -log10(0.05/independent_tests), mapped_cutoff = "BF")[[1]] +
  base_theme +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.border = element_rect(fill = NA),
        plot.title = element_blank()) +
  ggplot2::labs(x = "Genomic Position (Mb)",
                y = expression(-log[10](italic(p))))

ggsave(filename = "Plots/Albendazole/ABZ_q90TOF_GWA.png", height = 4, width = 12, dpi = 400)
ggsave(filename = "Plots/Albendazole/ABZ_q90TOF_GWA.pdf", height = 4, width = 12, dpi = 400)

q90burden_pr <- data.table::fread(glue::glue( "{albendazole_data}Final_Tables/Albendazole_q90.TOF.Skat.assoc"))%>%
  dplyr::rowwise()%>%
  dplyr::mutate(CHROM = strsplit(strsplit(RANGE,split = ",")[[1]][1],split = ":")[[1]][1],
                POS = as.numeric(strsplit(strsplit(strsplit(RANGE,split = ",")[[1]][1],split = ":")[[1]][2],split = "-")[[1]][1]),
                endPOS = as.numeric(strsplit(strsplit(strsplit(RANGE,split = ",")[[1]][1],split = ":")[[1]][2],split = "-")[[1]][2]))%>%
  dplyr::mutate(size = abs(POS-endPOS))%>%
  dplyr::filter(CHROM!="MtDNA")%>%
  dplyr::ungroup()%>%
  dplyr::filter(NumVar > 1)%>%
  dplyr::mutate(significant = ifelse(PermPvalue < .05/n(), TRUE,FALSE ))%>%
  dplyr::arrange(PermPvalue)

q90burden_pr%>%
  ggplot()+
  aes(x = POS/1e6, y = -log10(as.numeric(Pvalue)), alpha = significant, color = significant)+
  geom_point()+
  scale_color_manual(values=c("black","#D7263D"))+
  scale_alpha_manual(values = c(0.5,1))+
  facet_grid(.~CHROM, scales = "free", space = "free")+
  ggplot2::theme_bw() +
  base_theme +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.border = element_rect(fill = NA)) +
  labs(x = "Genomic Position (Mb)", 
       y = expression(-log[10](italic(p))))

ggsave("Plots/Albendazole/ABZ_q90TOF_SKAT.png",  dpi = 600, height = 4, width = 12)
ggsave("Plots/Albendazole/ABZ_q90TOF_SKAT.pdf",  dpi = 600, height = 4, width = 12)

######################################################################################################################## abz_colored by ben1 variants
ben1_indels <- data.table::fread(glue::glue("{albendazole_data}Final_Tables/ben1_variants.csv"))%>%
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

gwa_mappings <- pr_maps %>%
  na.omit()%>%
  dplyr::mutate(snpGT = ifelse(allele==-1,"REF", "ALT"))%>%
  dplyr::select(snpMarker = marker, strain, trait, value, snpGT)%>%
  dplyr::left_join(.,ben1_variants, by = "strain")

gwa_mappings$snpMarker <- gsub("_",":", gwa_mappings$snpMarker)

colors <- c("Deletion" = "#E69F00","Insertion" = "#56B4E9","Inversion" = "#D55E00","Stop Gained" = "#F0E442",
            "Transposon\nInsertion" = "#0072B2","Missense" = "#009E73","Splice Donor" = "#CC79A7","None" = "#999999")

gwa_mappings <- gwa_mappings%>%
  dplyr::mutate(length = end-start)%>%
  dplyr::mutate(ben1_prediction = ifelse(is.na(GT), "A", 
                                         ifelse(grepl("del", marker),"Deletion",
                                                ifelse(grepl("ins", marker,ignore.case = T),"Insertion",
                                                       ifelse(grepl("inv", marker),"Inversion",
                                                              ifelse(grepl("stop", marker),"Stop Gained",
                                                                     ifelse(grepl("trans", marker),"Transposon\nInsertion",
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
  dplyr::mutate(build_GT = ifelse(ben1_prediction %in% c("Missense","Stop Gained", "Splice Donor","Deletion", "Insertion", "Inversion", "Transposon\nInsertion"), ben1_prediction, "A"))


plot_bar(ben1_all, pt_colors = colors) + 
  base_theme +
  theme(axis.text.x = ggplot2::element_blank(),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        axis.line = element_line(colour = "black"))

ggsave("Plots/Albendazole/ABZ_pheno_ben1_colors.png", dpi = 300, height = 6, width = 10)
ggsave("Plots/Albendazole/ABZ_pheno_ben1_colors.pdf", dpi = 300, height = 6, width = 10)

######################################################################################################################## Ben-1 world distribution

ben1_variation_matrix <- gwa_mappings%>%
  dplyr::ungroup()%>%
  dplyr::select(strain,marker)%>%
  dplyr::distinct()

missing_strains <- data.frame(strain = row.names(cegwas::kinship)[!row.names(cegwas::kinship)%in%ben1_variation_matrix$strain],
                              marker = NA, stringsAsFactors = F)

colors <- c("Deletion" = "#E69F00","Insertion" = "#56B4E9","Inversion" = "#D55E00","Stop Gained" = "#F0E442",
            "Transposon\nInsertion" = "#0072B2","Missense" = "#009E73","Splice Donor" = "#CC79A7","None" = "#999999")

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


world <- map_data("world")
world <- world[world$region != "Antarctica",] # intercourse antarctica

isolation_info <- readr::read_tsv("https://elegansvariation.org/strain/strain_data.tsv")%>%
  dplyr::left_join(ben1_strain_markers,.,by="strain")%>%
  dplyr::select(strain, marker, ben1_prediction, long=longitude, lat=latitude, landscape, substrate)

isolation_info$long <- as.numeric(isolation_info$long)
isolation_info$lat <- as.numeric(isolation_info$lat)

map <- ggplot()+ geom_map(data=world, map=world,
                          aes(x=long, y=lat, map_id=region),
                          color=axis_color, fill=background_color, size= 0.5, alpha=1)+
  geom_point(data=dplyr::filter(isolation_info, !is.na(marker)),  
             aes(x=long, y=lat, fill=ben1_prediction), size = 4, shape = 21)+
  scale_fill_manual(values=colors,name = expression(paste("Variation at ", italic("ben-1"))))+
  theme_map()+
  theme(panel.background = element_rect(fill = background_color, colour = NA),
        text = element_text(family = axes_title_font, size = axes_text_size),
        legend.text = element_text(size = 22, family = number_font),
        legend.title = element_text(size = 22),
        legend.key.size = unit(1, "cm"),
        legend.background = element_rect(fill = background_color)) 

ggsave("Plots/Albendazole/ben-1_world_dist.png", height = 10,width  = 20, dpi = 300)
ggsave("Plots/Albendazole/ben-1_world_dist.pdf", height = 10,width  = 20, dpi = 300)



######################################################################################################################## Ben-1 popgen
# Rscript --vanilla Interval_Popgen.R II 2000000 5000000 Ce330_annotated.vcf.gz WS245_exons.gff BEN1 249_samples.txt

tree <-  ape::read.tree( file = glue::glue("{albendazole_data}Final_Tables/genome.tree"))

tree <- ggtree(tree, branch.length = "rate") +xlim(-.075,.475)

tree

ben1_variation_matrix <- gwa_mappings%>%
  dplyr::ungroup()%>%
  dplyr::select(strain,marker)%>%
  dplyr::distinct()

missing_strains <- data.frame(strain = row.names(cegwas::kinship)[!row.names(cegwas::kinship)%in%ben1_variation_matrix$strain],
                              marker = NA, stringsAsFactors = F)

colors <- c("Deletion" = "#E69F00","Insertion" = "#56B4E9","Inversion" = "#D55E00","Stop Gained" = "#F0E442",
            "Transposon\nInsertion" = "#0072B2","Missense" = "#009E73","Splice Donor" = "#CC79A7","None" = "#999999")

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

tree %<+% na.omit(ben1_strain_markers)+ 
  geom_tippoint(aes( color=ben1_prediction), size = 2.5)+
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

ggsave("Plots/Albendazole/ben-1_tree.png", height = 4 ,width  = 20, dpi = 300)
ggsave("Plots/Albendazole/ben-1_tree.pdf", height = 4,width  = 20, dpi = 300)

########################################################################################################################


load("Run_Popgen/BEN1_II_Statistics.Rda")

ben1_start <- 3537688
ben1_end <- 3541628

td_df <- neutrality_df %>%
  dplyr::filter(statistic %in% c("Fay.Wu.H","Zeng.E","Tajima.D")) %>%
  dplyr::group_by(statistic) %>%
  dplyr::mutate(scaled_value = scale(value))

td_df %>%
  ggplot()+
  aes(x = WindowPosition/1e6, y = value, color = statistic)+
  geom_point(size = point_size, alpha = point_alpha)+
  scale_color_manual(values = c("#BE0032","#0067A5","#222222"))+
  # facet_grid(statistic~., scales = "free")+
  base_theme +
  geom_vline(aes(xintercept = ben1_start/1e6), color = "gray60", linetype='dashed')+
  geom_vline(aes(xintercept = ben1_end/1e6), color = "gray60", linetype='dashed')+
  theme(panel.grid.major = element_blank(),
        axis.line = element_line(colour = axis_color),
        axis.title.y = element_blank()) +
  labs(x = "Genomic Position (Mb)", color = "Statistic")+xlim(c(3.427688, 3.651628))

ggsave(filename = "Plots/Albendazole/ben1_popgen.png", height = 8, width = 12, dpi = 400)
ggsave(filename = "Plots/Albendazole/ben1_popgen.pdf", height = 8, width = 12, dpi = 400)


######################################################################################################################## HTA assay of F200Y and DEL alleles

pr_df <- data.table::fread(glue::glue("{albendazole_data}Final_Tables/TS13_HTA_ben-1_regressed.tsv"))

main_figure <- pr_df%>%
  dplyr::filter(strain %in% c("N2", "ECA883", "ECA919"))%>%
  dplyr::arrange(phenotype)%>%
  dplyr::mutate(strain2 = factor(strain, levels = unique(strain), labels = unique(strain), ordered = T))%>%
  dplyr::mutate(norm_pheno_temp = ifelse(phenotype == min(phenotype), 0, 1))%>%
  dplyr::mutate(delta_pheno = ifelse(norm_pheno_temp == 0, 0, abs(dplyr::lag(phenotype) - phenotype)))%>%
  dplyr::mutate(norm_pheno = cumsum(delta_pheno)) %>%
  dplyr::mutate(final_pheno = norm_pheno/max(norm_pheno))

main_figure%>%
  dplyr::filter(trait == trait_of_interest, condition == "Albendazole")%>%
  dplyr::mutate(strain1 = factor(group, levels = c("N2","del","F200Y"),
                                 labels =c("N2","Del","F200Y")))%>%
  ggplot(.) +
  aes(x = factor(strain1), 
      y = final_pheno, 
      fill=group) +
  geom_boxplot(outlier.colour = NA, alpha = 1)+
  scale_fill_manual(values = c("#E69F00", "#009E73","#999999"))+
  theme_bw()+
  labs( y = "Relative ABZ resistance")+
  base_theme+
  theme(axis.title.x = element_blank(),
        legend.position="none",
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.line = element_line(colour = axis_color))

ggsave(filename = "Plots/Albendazole/ABZ_allele_HTA.png", height = 6, width = 8, dpi = 400)
ggsave(filename = "Plots/Albendazole/ABZ_allele_HTA.pdf", height = 6, width = 8, dpi = 400)

######################################################################################################################## competition assay of F200Y and DEL alleles

competition_assay <- data.table::fread(glue::glue("{albendazole_data}Final_Tables/TS14_competition_assay.csv"))%>%
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
  dplyr::filter(Strain == "N2")

competition_assay_outliers%>%
  dplyr::ungroup()%>%
  dplyr::mutate(id = paste0(Strain,Condition))%>%
  ggplot(aes(x=Generation, y=Mean, fill=Condition, color=Strain))+
  geom_hline(aes(yintercept = 50), color = "gray70", linetype = "dashed")+
  geom_errorbar(aes(ymin=Mean - SND, ymax=Mean + SND, alpha=Condition), width=.1,size =1)+
  geom_line(aes(alpha=Condition, group=id), size =2)+
  scale_color_manual(values = c("#E69F00", "#009E73","#999999"))+
  scale_alpha_manual(values = c(1,0.2))+
  theme_bw()+
  theme(legend.title = element_text(colour="black", size = 15, face = "bold"))+
  theme(legend.text = element_text(colour="black", size = 14))+
  labs( x = "Generation")+
  labs( y = expression(bold(paste("Allele Frequency of ", bolditalic("ben-1"), " Edits (%)"))))+
  base_theme+
  theme(legend.position="none",
        panel.grid.major.x = element_blank())+
  theme(legend.position="none",
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text = element_text(size = 24),
        axis.line = element_line(colour = axis_color))+
  scale_y_continuous(breaks = c(0,50,100), labels = c(0, 0.5, 1), limits = c(0,100))

ggsave(filename = "Plots/Albendazole/ABZ_competition.png", height = 6, width = 10, dpi = 400)
ggsave(filename = "Plots/Albendazole/ABZ_competition.pdf", height = 6, width = 10, dpi = 400)

######################################################################################################################## Ben-1 removed


pr_maps <- data.table::fread(glue::glue("{albendazole_data}Final_Tables/Albendazole_mean.TOF_processed_mapping.tsv"))

independent_tests <- 1159.75

cegwas2_manplot(plot_df = pr_maps, eigen_cutoff = -log10(0.05/independent_tests), mapped_cutoff = "EIGEN")[[1]] +
  base_theme +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.border = element_rect(fill = NA),
        plot.title = element_blank()) +
  ggplot2::labs(x = "Genomic Position (Mb)",
                y = expression(-log[10](italic(p))))

ggsave(filename = "Plots/Albendazole/ABZ_meanTOF_ben1v_Removed_GWA.png", height = 4, width = 12, dpi = 400)
ggsave(filename = "Plots/Albendazole/ABZ_meanTOF_ben1v_Removed_GWA.pdf", height = 4, width = 12, dpi = 400)





# random
isolation_info <- readr::read_tsv("https://elegansvariation.org/strain/strain_data.tsv")

strains_249 <- isolation_info%>%
  dplyr::filter(release %in% c("20160408", "20170531"), reference_strain == "True")%>%
  dplyr::select(strain = isotype, long = longitude, lat = latitude) %>%
  dplyr::pull(strain)


ben1 <- cegwas2::query_vcf("ben-1", vcf = "Data/Ce330_annotated.vcf.gz", impact = "ALL") %>%
  dplyr::filter(nchar(REF) == nchar(ALT)) 

ben1_syn <- dplyr::filter(ben1, effect == "synonymous_variant") %>%
  dplyr::distinct(CHROM,POS, .keep_all = T)

ben1_nonsyn <- dplyr::filter(ben1, effect %in% c("missense_variant","stop_gained","splice_donor_variant&intron_variant")) %>%
  dplyr::distinct(CHROM,POS, .keep_all = T)
