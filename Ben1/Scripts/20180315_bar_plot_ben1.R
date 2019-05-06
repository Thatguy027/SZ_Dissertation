

# set to location of files
main.dir <- "~/Dropbox/HTA/Results/2017_anthelmintic_gwa/analysis/ben1_paper/"
data.dir <- paste0(main.dir,"Data/")
plot.dir <- paste0(main.dir,"Plots/")
script.dir <- paste0(main.dir,"Scripts/")

source(paste0(script.dir,"ben1_processing_functions.R"))

trait_of_interest <- "q90.TOF"
PCA_of_interest <- "PC1"
comparison_trait <- "q90.EXT"
condition_of_interest <- "albendazole"
control_of_interest <- "DMSO"



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



plot_bar <- function(df){
ggplot(df)+
  aes(x=strain2, y=final_pheno, fill = factor(build_GT,
                                              levels = sort(unique(build_GT)), 
                                              labels= c("None", sort(unique(build_GT)[2:8]) )))+
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("None" = "gray75","Deletion" = colors[1],
                             "Insertion" = colors[2],
                             "Inversion" = colors[3],
                             "Missense" = colors[4],
                             "Splice Donor" = colors[5],
                             "Stop Gained" = colors[6],
                             "Transposon Insertion" = colors[7]),name = expression(paste("Variation at ", italic("ben-1"))))+
  theme_bw()+
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 0),
                 axis.text.y = ggplot2::element_text(size = 16),
                 axis.title.x = ggplot2::element_text(size = 20, face = "bold", color = "black", vjust = -0.3), 
                 axis.title.y = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                 strip.text.x = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                 strip.text.y = ggplot2::element_text(size = 20, face = "bold", color = "black"),
                                                      plot.title = ggplot2::element_text(size = 24, face = "bold", vjust = 1), 
                                                      # panel.background = ggplot2::element_rect(color = "black",size = 1.2),
                                                      panel.grid.major = element_blank(), 
                                                      panel.grid.minor = element_blank(),
                                                      panel.background = element_blank())+
                                                      labs(x = "Strain", y = paste0("Resistance to Albendazole"))
                                                      
} 

#### generate different dataframes/plots #### 

#none #

ben1_none <- gwa_mappings_bar%>%
  dplyr::mutate(build_GT = ifelse(ben1_prediction == "Missense", "A", "A"))

ben1_missense$build_GT[is.na(ben1_missense$build_GT)] <- "A"

plot_bar(ben1_none)

ggsave(paste0(plot.dir,"GWA_",trait_of_interest,"_barplot_non_colored.png"), 
       height = 4, 
       width = 12)


# missense #
ben1_missense <- gwa_mappings_bar%>%
  dplyr::mutate(build_GT = ifelse(ben1_prediction == "Missense", "Missense", "A"))

ben1_missense$build_GT[is.na(ben1_missense$build_GT)] <- "A"

plot_bar(ben1_missense)

ggsave(paste0(plot.dir,"GWA_",trait_of_interest,"_barplot_colored_by_missense.png"), 
height = 4, 
width = 12)


# SNVs #
ben1_SNV<- gwa_mappings_bar%>%
  dplyr::mutate(build_GT = ifelse(ben1_prediction %in% c("Missense", "Stop Gained", "Splice Donor"), ben1_prediction, "A"))

ben1_SNV$build_GT[is.na(ben1_SNV$build_GT)] <- "A"

plot_bar(ben1_SNV)

ggsave(paste0(plot.dir,"GWA_",trait_of_interest,"_barplot_colored_by_SNV.png"), 
       height = 4, 
       width = 12)

# + deletions # 

ben1_SNV_del <- gwa_mappings_bar%>%
  dplyr::mutate(build_GT = ifelse(ben1_prediction %in% c("Missense","Stop Gained", "Splice Donor","Deletion"), ben1_prediction, "A"))

ben1_SNV_del$build_GT[is.na(ben1_SNV_del$build_GT)] <- "A"

plot_bar(ben1_SNV_del)

ggsave(paste0(plot.dir,"GWA_",trait_of_interest,"_barplot_colored_by_SNV_del.png"), 
       height = 4, 
       width = 12)

# + insertions # 

ben1_SNV_del_ins <- gwa_mappings_bar%>%
  dplyr::mutate(build_GT = ifelse(ben1_prediction %in% c("Missense","Stop Gained", "Splice Donor","Deletion", "Insertion"), ben1_prediction, "A"))

ben1_SNV_del_ins$build_GT[is.na(ben1_SNV_del_ins$build_GT)] <- "A"

plot_bar(ben1_SNV_del_ins)

ggsave(paste0(plot.dir,"GWA_",trait_of_interest,"_barplot_colored_by_SNV_del_ins.png"), 
       height = 4, 
       width = 12)


# all #

ben1_all <- gwa_mappings_bar%>%
  dplyr::mutate(build_GT = ifelse(ben1_prediction %in% c("Missense","Stop Gained", "Splice Donor","Deletion", "Insertion", "Inversion", "Transposon Insertion"), ben1_prediction, "A"))

ben1_all$build_GT[is.na(ben1_all$build_GT)] <- "A"

plot_bar(ben1_all)

ggsave(paste0(plot.dir,"GWA_",trait_of_interest,"_barplot_colored_by_all.png"), 
       height = 4, 
       width = 12)







#ggsave(paste0(plot.dir,"GWA_",trait_of_interest,"_barplot_colored_by_ben1.pdf"), 
                                                      height = 4, 
                                                      width = 12)















