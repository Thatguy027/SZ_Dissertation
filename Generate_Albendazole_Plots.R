try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)))

# data folders
albendazole_data <- "201809_rappaport/Data/ben1_data/"

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

cegwas2_manplot(plot_df = pr_maps, eigen_cutoff = -log10(0.05/independent_tests), mapped_cutoff = "EIGEN")[[1]] +
  base_theme +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.border = element_rect(fill = NA),
        plot.title = element_blank()) +
  ggplot2::labs(x = "Genomic Position (Mb)",
                y = expression(-log[10](italic(p))))

ggsave(filename = "Plots/ABZ_q90TOF_GWA.png", height = 4, width = 12, dpi = 400)

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

ggsave("Plots/ABZ_q90TOF_SKAT.png", 
       dpi = 600,
       height = 4, 
       width = 12)


######################################################################################################################## Ben-1 regressed




pr_maps <- data.table::fread(glue::glue("{albendazole_data}Final_Tables/TS17_GWA_marker_mappings_ben1_regressed.tsv"))

independent_tests <- 8000

cegwas2_manplot(plot_df = pr_maps, eigen_cutoff = -log10(0.05/independent_tests), mapped_cutoff = "BF")[[1]] +
  base_theme +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.border = element_rect(fill = NA),
        plot.title = element_blank()) +
  ggplot2::labs(x = "Genomic Position (Mb)",
                y = expression(-log[10](italic(p))))

ggsave(filename = "Plots/ABZ_q90TOF_ben1v_Removed_GWA.png", height = 4, width = 12, dpi = 400)

######################################################################################################################## Ben-1 removed




pr_maps <- data.table::fread(glue::glue("{albendazole_data}Final_Tables/ben1v_rem_Albendazole_mean.TOF_processed_mapping.tsv"))

independent_tests <- 1159.75

cegwas2_manplot(plot_df = pr_maps, eigen_cutoff = -log10(0.05/independent_tests), mapped_cutoff = "EIGEN")[[1]] +
  base_theme +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.border = element_rect(fill = NA),
        plot.title = element_blank()) +
  ggplot2::labs(x = "Genomic Position (Mb)",
                y = expression(-log[10](italic(p))))

ggsave(filename = "Plots/ABZ_q90TOF_ben1v_Removed_GWA.png", height = 4, width = 12, dpi = 400)


########################################################################################################################


load("Data/ben1_III_3000000-4000000_Diversity_Statistics.Rda")

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
  facet_grid(statistic~., scales = "free")+
  base_theme +
  geom_vline(aes(xintercept = ben1_start/1e6), color = "#E68FAC")+
  geom_vline(aes(xintercept = ben1_end/1e6), color = "#E68FAC")+
  theme(panel.grid.major = element_blank(),
        legend.position = "none",
        axis.line = element_line(colour = axis_color),
        axis.title.y = element_blank()) +
  labs(x = "Genomic Position (Mb)")



# 249 strains

load("Data/BEN1_249_III_3000000-4000000_Diversity_Statistics.Rda")

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
  facet_grid(statistic~., scales = "free")+
  base_theme +
  geom_vline(aes(xintercept = ben1_start/1e6), color = "#E68FAC")+
  geom_vline(aes(xintercept = ben1_end/1e6), color = "#E68FAC")+
  theme(panel.grid.major = element_blank(),
        legend.position = "none",
        axis.line = element_line(colour = axis_color),
        axis.title.y = element_blank()) +
  labs(x = "Genomic Position (Mb)") 

# 249 strains - zoom

load("Data/BEN1_249_III_3000000-4000000_Diversity_Statistics.Rda")

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
  facet_grid(statistic~., scales = "free")+
  base_theme +
  geom_vline(aes(xintercept = ben1_start/1e6), color = "#E68FAC")+
  geom_vline(aes(xintercept = ben1_end/1e6), color = "#E68FAC")+
  theme(panel.grid.major = element_blank(),
        legend.position = "none",
        axis.line = element_line(colour = axis_color),
        axis.title.y = element_blank()) +
  labs(x = "Genomic Position (Mb)") +xlim(c(3.427688, 3.651628))

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
