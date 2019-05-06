# Figure script for ben-1 manuscript



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
  control_oof_interest <- "DMSO"

# # # # # # ## # # # # # ## # # # # # ## # # # # # ## # # # # # ## # # # # # ## # # # # # ## # # # # # ## # # # # # ## # # # # # ## # # # # # ## # # # # # ## # # # # # # DOSE RESPONSE EXPERIMENT
dirs <- paste0(data.dir,"DoseResponse/")

raw <- read_data(dirs)
raw_nocontam <- remove_contamination(raw)
summedraw <- sumplate(raw_nocontam, directories = FALSE, quantiles = TRUE)
biopruned <- bioprune(summedraw)

pr_df <- biopruned%>%
  dplyr::ungroup()%>%
  dplyr::select(-date, -round,-experiment,-assay)%>%
  tidyr::gather(trait, value, -condition,-control,-strain, -plate, -row, -col)%>%
  dplyr::filter(!is.na(condition))

write.table(pr_df , 
            file = paste0(data.dir,"DR_Processed.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F)

pr_df%>%
  dplyr::filter(trait == trait_of_interest)%>%
  # dplyr::
  ggplot()+
  aes(x = factor(condition, 
                 levels = c("DMSO", "albendazole3125", "albendazole625", "albendazole125",  "albendazole25"),
                 labels = c("0", "3.125", "6.25", "12.5", "25")), 
      y = value, 
      fill = strain)+
  geom_boxplot()+
  labs(y = paste0("Animal Length (", trait_of_interest,")"),
       x = paste0("Albendazole Concentration (µM)"))+
  theme_bw()+
  scale_fill_manual(values = c("blue","cadetblue3","hotpink3","orange"), 
                    labels = c("CB4856", "DL238", "JU775", "N2"),
                    name = "Strain")+
  theme(axis.text.x = ggplot2::element_text(size = 14),
        axis.text.y = ggplot2::element_text(size = 14),
        # legend.position = "none",
        axis.title.x = ggplot2::element_text(size = 18, face = "bold", color = "black", vjust = -0.3),
        axis.title.y = ggplot2::element_text(size = 18, face = "bold", color = "black", vjust = -0.3))

ggsave(paste0(plot.dir,"DR_",condition_of_interest,"_",trait_of_interest,".pdf"), 
       height = 4, 
       width = 10)

# # # # # # ## # # # # # ## # # # # # ## # # # # # ## # # # # # ## # # # # # ## # # # # # ## # # # # # ## # # # # # ## # # # # # ## # # # # # ## # # # # # ## # # # # # # PCA DOSE RESPONSE

ben1_pca <- biopruned%>%
  dplyr::ungroup()%>%
  dplyr::select(-(date:assay))%>%
  dplyr::filter(!is.na(strain))

ben1_pca_strains <- ben1_pca$strain
ben1_pca_condition <- ben1_pca$condition

ben1_pca1 <- ben1_pca%>%
  dplyr::select(-(plate:col))

row.names(ben1_pca1) <- ben1_pca_strains

# NEED to scale phenotypes to set mean = 0 and standard deviation = 1 for a all traits so PC loadings are not driven by traits with largest variances
scaled_ben1_pca <- data.frame(scale(ben1_pca1))

sapply(scaled_ben1_pca,mean)
sapply(scaled_ben1_pca,sd)

ben1_pca_analysis <- prcomp(scaled_ben1_pca)
# LOOK AT LOADINGS
sort(ben1_pca_analysis$rotation[,1])
sort(ben1_pca_analysis$rotation[,2])
sort(ben1_pca_analysis$rotation[,3])

# PLOT VE BY PC
vars <- apply(ben1_pca_analysis$x, 2, var)  
props <- vars / sum(vars)
pca_var_df <- data.frame(PC = as.numeric(gsub("PC","",names(props))), Variance = cumsum(props))

ggplot(pca_var_df)+
  aes(x = PC, y = Variance)+
  geom_point()+
  theme_bw()+
  theme(legend.position = 'none')+
  labs(x="Principal Component",y = "Cummulative Variance Explained")+
  theme(axis.text.x = ggplot2::element_text(size = 14),
        axis.text.y = ggplot2::element_text(size = 14),
        legend.position = "none",
        axis.title.x = ggplot2::element_text(size = 18, face = "bold", color = "black", vjust = -0.3),
        axis.title.y = ggplot2::element_text(size = 18, face = "bold", color = "black", vjust = -0.3),
        strip.text.x = element_text(size = 14, face = "bold"))+
  ylim(c(0,1))

pc_df <- data.frame(strain = ben1_pca_strains, condition = ben1_pca_condition, PC1 = ben1_pca_analysis$x[,1])

pc_df%>%
  ggplot()+
  aes(x = factor(condition, 
                 levels = c("DMSO", "albendazole3125", "albendazole625", "albendazole125",  "albendazole25"),
                 labels = c("0", "3.125", "6.25", "12.5", "25")), 
      y = PC1, 
      fill = strain)+
  geom_boxplot(outlier.colour = NA)+
  labs(y = PCA_of_interest,
       x = paste0("Albendazole Concentration \u03BCM"))+
  theme_bw()+
  scale_fill_manual(values = c("blue","cadetblue3","hotpink3","orange"), 
                    labels = c("CB4856", "DL238", "JU775", "N2"),
                    name = "Strain")+
  theme(axis.text.x = ggplot2::element_text(size = 14),
        axis.text.y = ggplot2::element_text(size = 14),
        # legend.position = "none",
        axis.title.x = ggplot2::element_text(size = 18, face = "bold", color = "black", vjust = -0.3),
        axis.title.y = ggplot2::element_text(size = 18, face = "bold", color = "black", vjust = -0.3))

ggsave(paste0(plot.dir,"DR_",condition_of_interest,"_",PCA_of_interest,".pdf"), device=cairo_pdf, 
       height = 4, 
       width = 10)

# # # # # # ## # # # # # ## # # # # # ## # # # # # ## # # # # # ## # # # # # ## # # # # # ## # # # # # ## # # # # # ## # # # # # ## # # # # # ## # # # # # ## # # # # # # GWAS DATA

# # # # # # # EXPERIMENT 
# # replicates
# A1A 
# A1B
# B1A
# B1B
# # # # For every strain four independent bleaches were performed over two days 
# # # # 
# # # # The first letter corresponds to the day. 
# # # # Only one drug prep was made per day... So two days have independent drug preps. 
# # # # Within each day there are independent bleaches. 
# # # # Number corresponds to strain set
# # # # The second letter corresponds to bleach on given drug prep

# Define a vector of your experiement directories
dirs <- c(paste0(data.dir,"GWAS/20170605_GWAA1B/"),
          paste0(data.dir,"GWAS/20170605_GWAA1A/"),
          paste0(data.dir,"GWAS/20170612_GWAA2A/"),
          paste0(data.dir,"GWAS/20170612_GWAA2B/"),
          paste0(data.dir,"GWAS/20170619_GWAA3A/"),
          paste0(data.dir,"GWAS/20170619_GWAA3B/"),
          paste0(data.dir,"GWAS/20170626_GWAA4A/"),
          paste0(data.dir,"GWAS/20170626_GWAA4B/"),
          paste0(data.dir,"GWAS/20170717_GWAA5A/"),
          paste0(data.dir,"GWAS/20170717_GWAA5B/"))


# # Read in the data
raw_a <- easysorter::read_data(dirs)

dirs <- c(paste0(data.dir,"GWAS/20170606_GWAB1A/"),
          paste0(data.dir,"GWAS/20170606_GWAB1B/"),
          paste0(data.dir,"GWAS/20170613_GWAB2A/"),
          paste0(data.dir,"GWAS/20170613_GWAB2B/"),
          paste0(data.dir,"GWAS/20170620_GWAB3A/"),
          paste0(data.dir,"GWAS/20170620_GWAB3B/"),
          paste0(data.dir,"GWAS/20170627_GWAB4A/"),
          paste0(data.dir,"GWAS/20170627_GWAB4B/"),
          paste0(data.dir,"GWAS/20170718_GWAB5A/"),
          paste0(data.dir,"GWAS/20170718_GWAB5B/"))

# Read in the data
raw_b <- easysorter::read_data(dirs)

# Remove all data from the contaminated wells
raw_nocontam_a <- easysorter::remove_contamination(raw_a)
raw_nocontam_b <- easysorter::remove_contamination(raw_b)

# Summarize the data
summedraw_a <- easysorter::sumplate(raw_nocontam_a, directories = TRUE, quantiles = TRUE)
summedraw_b <- easysorter::sumplate(raw_nocontam_b, directories = TRUE, quantiles = TRUE)

#Prune based on biological impossibilities
biopruned_a <- easysorter::bioprune(summedraw_a)
biopruned_b <- easysorter::bioprune(summedraw_b)

# removes any remaining wash wells. 
all_data <- dplyr::bind_rows(list(biopruned_a,biopruned_b))%>%
  dplyr::filter(!is.na(strain))

write.table(all_data , 
            file = paste0(data.dir,"GWAS_summarized.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F)

# # # # IDENTICAL UP TO HERE
# all_data <- data.table::fread(paste0(data.dir,"GWAS_summarized.tsv"))

# FILTER OUTLIERS
long_data <- all_data %>% 
  tidyr::gather(., trait, phenotype, -date, -experiment,-round,-assay,-plate,-condition,-control,-strain, -row, -col)%>%
  data.frame()%>%
  dplyr::filter(!grepl("red|yellow|green", trait))


outlier_replicates <- long_data %>%
  dplyr::group_by(condition, round, trait, strain)%>%
  dplyr::mutate(n_reps = n())%>%
  dplyr::filter(n_reps > 3)%>%
  # dplyr::select(-n_reps)%>%
  dplyr::mutate(med_ph = median(phenotype, na.rm = T),
                sd_ph = sd(phenotype, na.rm= T))%>%
  dplyr::rowwise()%>%
  dplyr::mutate(outlier = ifelse((phenotype > med_ph + 1.8*sd_ph) | (phenotype < med_ph - 1.8*sd_ph), "OUTLIER", "OK"))%>%
  ungroup()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # VISUALIZE OUTLIERS

outlier_replicates%>%
  dplyr::filter(trait == trait_of_interest)%>%
  ggplot(.)+
  aes(x = strain, y = phenotype)+
  geom_boxplot(outlier.colour = NA)+
  geom_jitter(aes(color = outlier), width = 0)+
  labs(y = paste0("Animal Length (", trait_of_interest,")"),
       x = paste0("Albendazole Concentration \u03BCM"))+
  theme_bw()+
  facet_grid(condition~.)+
  scale_fill_manual(values = c("blue","cadetblue3","hotpink3","orange"), 
                    labels = c("CB4856", "DL238", "JU775", "N2"),
                    name = "Strain")+
  theme(axis.text.x = ggplot2::element_text(size = 14),
        axis.text.y = ggplot2::element_text(size = 14),
        # legend.position = "none",
        axis.title.x = ggplot2::element_text(size = 18, face = "bold", color = "black", vjust = -0.3),
        axis.title.y = ggplot2::element_text(size = 18, face = "bold", color = "black", vjust = -0.3))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # VISUALIZE OUTLIERS


outlier_replicates1 <- dplyr::ungroup(outlier_replicates)%>%
  tidyr::unite(outlier_id, experiment, round,assay,condition, trait, plate,strain, row, col)%>%
  dplyr::select(outlier_id, outlier)

long_data2 <- dplyr::ungroup(long_data)%>%
  tidyr::unite(outlier_id, experiment, round,assay,condition, trait, plate,strain, row, col)%>%
  dplyr::left_join(.,outlier_replicates1, by = "outlier_id")%>%
  dplyr::ungroup()%>%
  dplyr::filter(outlier == "OK")%>%
  tidyr::separate(outlier_id, into = c("experiment", "round","assay","condition", "trait", "plate","strain", "row", "col"),sep = "_")%>%
  dplyr::group_by(condition, round, trait)%>%
  tidyr::unite(id , date,condition,round,strain,phenotype,plate,row,col,trait,condition,assay,remove = F, sep = "_")

row.names(long_data2) <- long_data2$id

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  REGRESSION
resid_data <- long_data2%>%
  do(augment(lm(phenotype ~ experiment + assay, .)))%>%
  dplyr::select(condition,round,trait,experiment,assay,phenotype, resid_phenotype = .resid)

pr_residual_df <- long_data2%>%
  dplyr::ungroup()%>%
  dplyr::left_join(.,resid_data,by=c("condition","round","trait","experiment","assay","phenotype"))%>%
  dplyr::distinct(date,condition,round,strain,phenotype,plate,row,col,trait,condition,assay,.keep_all=T)

DMSO <-pr_residual_df%>%
  dplyr::filter(condition=="DMSO")%>%
  dplyr::group_by(strain,trait)%>%
  dplyr::summarise(ctrl_phenotype = mean(phenotype),
                   ctrl_resid_phenotype = mean(resid_phenotype))

drug <-pr_residual_df%>%
  dplyr::filter(condition!="DMSO")%>%
  dplyr::left_join(.,DMSO,by=c("trait",'strain'))

ctrl_regressed <- drug %>%
  dplyr::select(condition,trait,strain,resid_phenotype,ctrl_resid_phenotype)%>%
  group_by(condition,trait)%>%
  do(augment(lm(resid_phenotype ~ ctrl_resid_phenotype, .)))%>%
  dplyr::left_join(.,drug,by=c('condition','trait','resid_phenotype','ctrl_resid_phenotype'))%>%
  dplyr::select(condition, trait,strain, phenotype = .resid)%>%
  dplyr::group_by(condition,trait,strain)%>%
  dplyr::summarise(phenotype = mean(phenotype, na.rm = T)) %>%
  tidyr::unite(cond_trait, condition, trait) %>%
  dplyr::rename(trait = cond_trait)%>%
  tidyr::spread(trait, phenotype)

write.table(ctrl_regressed, 
            file = paste0(data.dir,"GWAS_processed_phenotype.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # PCA 

gwa_pca_df <- na.omit(ctrl_regressed)

gwa_pca_strains <- gwa_pca_df$strain

gwa_pca <- gwa_pca_df%>%
  dplyr::ungroup()%>%
  dplyr::select(-strain)

row.names(gwa_pca) <- gwa_pca_strains

# NEED to scale phenotypes to set mean = 0 and standard deviation = 1 for a all traits so PC loadings are not driven by traits with largest variances
scaled_gwa_pca <- data.frame(scale(gwa_pca))

sapply(scaled_gwa_pca,mean)
sapply(scaled_gwa_pca,sd)

ben1_pca_analysis <- prcomp(scaled_gwa_pca)
# LOOK AT LOADINGS
sort(ben1_pca_analysis$rotation[,1])
sort(ben1_pca_analysis$rotation[,2])
sort(ben1_pca_analysis$rotation[,3])

# PLOT VE BY PC
vars <- apply(ben1_pca_analysis$x, 2, var)  
props <- vars / sum(vars)
pca_var_df <- data.frame(PC = as.numeric(gsub("PC","",names(props))), Variance = cumsum(props))

gwa_pca_traits <- data.frame(strain = gwa_pca_strains, ben1_pca_analysis$x[,1:3])

gwa_raw_traits <- ctrl_regressed%>%
  tidyr::gather(cond_trait, phenotype, -strain)%>%
  tidyr::separate(cond_trait, into = c("condition", "trait"), sep = "_")%>%
  dplyr::select(-condition)%>%
  dplyr::filter(trait == trait_of_interest)%>%
  tidyr::spread(trait, phenotype)%>%
  dplyr::left_join(.,gwa_pca_traits, by = "strain")

write.table(gwa_raw_traits, 
            file = paste0(data.dir,"GWAS_traits_for_mapping.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F)

pr_pheno <- cegwas::process_pheno(gwa_raw_traits)
maps <- cegwas::gwas_mappings(pr_pheno)
pr_maps <- cegwas::process_mappings(maps, phenotype_df = pr_pheno)

manhattan_plots <- manplot_edit(pr_maps)

manhattan_plots[[1]]

ggsave(paste0(plot.dir,"GWA_PC1","_manplot.pdf"), 
       height = 4, 
       width = 12)

manhattan_plots[[2]]

ggsave(paste0(plot.dir,"GWA_",trait_of_interest,"_manplot.pdf"), 
       height = 4, 
       width = 12)

write.table(pr_maps, 
            file = paste0(data.dir,"GWAS_processed_mapping.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F)

genes <- process_correlations(variant_correlation(pr_maps, condition_trait = F))

write.table(genes, 
            file = paste0(data.dir,"GWAS_finemapping_genes_MOD_HIGH.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F)

# # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # # make phenotype for  testing

traits <- gwa_raw_traits%>%
  tidyr::gather(trait, value, -strain)%>%
  dplyr::mutate(Fam = "elegans", Sample = strain, Paternal = 0, Maternal = 0, Sex = 2)%>%
  tidyr::spread(trait,value)%>%
  dplyr::select(-strain)

traits[is.na(traits)] <- -9

burden_raw <- dplyr::select(traits, Fam:Sex, q90.TOF)

write.table(burden_raw, 
            file = paste0(data.dir,"GWAS_",trait_of_interest,".ped"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F)

burden_pc1 <- dplyr::select(traits, Fam:Sex, PC1)

write.table(burden_pc1, 
            file = paste0(data.dir,"GWAS_","PC1",".ped"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F)

burden_pc2 <- dplyr::select(traits, Fam:Sex, PC2)

write.table(burden_pc2, 
            file = paste0(data.dir,"GWAS_","PC2",".ped"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F)

burden_pc3 <- dplyr::select(traits, Fam:Sex, PC3)

write.table(burden_pc3, 
            file = paste0(data.dir,"GWAS_","PC3",".ped"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F)


# # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # # BURDEN TESTING IN TERMINAL

#vcfdir=/Users/stefan/Dropbox/AndersenLab/LabFolders/Stefan/burden/snp_indel_snpEff.vcf.gz
#refflatdir=/Users/stefan/Dropbox/AndersenLab/LabFolders/Stefan/burden/refFlat.ws245.txt

# docker run -i -t -v $(pwd):/home -w /home zhanxw/rvtests /rvtests/executable/rvtest \
# --pheno ben1_paper/GWAS_q90.TOF.ped --inVcf snp_indel_snpEff.vcf.gz --freqUpper .05 \
# --freqLower 0.003 --out q90TOF --geneFile refFlat.ws245.txt  --vt price
# 
# docker run -i -t -v $(pwd):/home -w /home zhanxw/rvtests /rvtests/executable/rvtest --pheno ben1_paper/GWAS_q90.TOF.ped --inVcf WI.20170531.impute.vcf.gz --freqUpper .05 --freqLower 0.003 --out q90TOF --geneFile refFlat.ws245.txt  --vt price --kernel skat
# docker run -i -t -v $(pwd):/home -w /home zhanxw/rvtests /rvtests/executable/rvtest --pheno ben1_paper/GWAS_PC1.ped --inVcf WI.20170531.impute.vcf.gz --freqUpper .05 --freqLower 0.003 --out PC1 --geneFile refFlat.ws245.txt  --vt price --kernel skat
# docker run -i -t -v $(pwd):/home -w /home zhanxw/rvtests /rvtests/executable/rvtest --pheno ben1_paper/GWAS_PC2.ped --inVcf WI.20170531.impute.vcf.gz --freqUpper .05 --freqLower 0.003 --out PC2 --geneFile refFlat.ws245.txt  --vt price --kernel skat
# docker run -i -t -v $(pwd):/home -w /home zhanxw/rvtests /rvtests/executable/rvtest --pheno ben1_paper/GWAS_PC3.ped --inVcf WI.20170531.impute.vcf.gz --freqUpper .05 --freqLower 0.003 --out PC3 --geneFile refFlat.ws245.txt  --vt price --kernel skat
# # # #VT PRICE
q90burden <- read.table(paste0(data.dir, "q90TOF.VariableThresholdPrice.assoc"),header = T)
q90burden$Gene <- as.character(q90burden$Gene)
q90burden$RANGE <- as.character(q90burden$RANGE)

q90burden<-q90burden%>%
  dplyr::rowwise()%>%
  dplyr::mutate(CHROM = strsplit(strsplit(RANGE,split = ",")[[1]][1],split = ":")[[1]][1],
                POS = as.numeric(strsplit(strsplit(strsplit(RANGE,split = ",")[[1]][1],split = ":")[[1]][2],split = "-")[[1]][1]),
                endPOS = as.numeric(strsplit(strsplit(strsplit(RANGE,split = ",")[[1]][1],split = ":")[[1]][2],split = "-")[[1]][2]))%>%
  dplyr::mutate(size = abs(POS-endPOS))

q90burden_pr<- q90burden%>%
  ungroup()%>%
  dplyr::mutate(significant = ifelse(PermPvalue < .05/n(), TRUE,FALSE ))%>%
  dplyr::filter(CHROM!="MtDNA",NumVar>1,size >500)

save(q90burden_pr, file = (paste0(data.dir, "q90TOF.VariableThresholdPrice_processed.Rda")))

q90burden_pr%>%
  ggplot()+
  aes(x = POS/1e6, y = Stat, size = NumVar, alpha = 0.5, color = significant)+
  geom_point()+
  scale_color_manual(values=c("black","red"))+
  facet_grid(.~CHROM, scales = "free", space = "free")+
  theme_bw()+
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 16),
                 axis.text.y = ggplot2::element_text(size = 16),
                 axis.title.x = ggplot2::element_text(size = 20, face = "bold", color = "black", vjust = -0.3), 
                 axis.title.y = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                 strip.text.x = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                 strip.text.y = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                 plot.title = ggplot2::element_text(size = 24, face = "bold", vjust = 1), 
                 panel.background = ggplot2::element_rect(color = "black",size = 1.2),
                 legend.position = "none")+
  labs(x = "Genomic Position (Mb)", y = "Test Statisitic")

ggsave(paste0(plot.dir,"GWA_",trait_of_interest,"_BURDEN_VTprice_manplot.pdf"), 
       height = 4, 
       width = 12)

pc1burden <- read.table(paste0(data.dir, "PC1.VariableThresholdPrice.assoc"),header = T)
pc1burden$Gene <- as.character(pc1burden$Gene)
pc1burden$RANGE <- as.character(pc1burden$RANGE)

pc1burden<-pc1burden%>%
  dplyr::rowwise()%>%
  dplyr::mutate(CHROM = strsplit(strsplit(RANGE,split = ",")[[1]][1],split = ":")[[1]][1],
                POS = as.numeric(strsplit(strsplit(strsplit(RANGE,split = ",")[[1]][1],split = ":")[[1]][2],split = "-")[[1]][1]),
                endPOS = as.numeric(strsplit(strsplit(strsplit(RANGE,split = ",")[[1]][1],split = ":")[[1]][2],split = "-")[[1]][2]))%>%
  dplyr::mutate(size = abs(POS-endPOS))

pc1burden_pr<- pc1burden%>%
  ungroup()%>%
  dplyr::mutate(significant = ifelse(PermPvalue < .05/n(), TRUE,FALSE ))%>%
  dplyr::filter(CHROM!="MtDNA",NumVar>1,size >500)

pc1burden_pr%>%
  ggplot()+
  aes(x = POS/1e6, y = Stat, size = NumVar, alpha = 0.5, color = significant)+
  geom_point()+
  scale_color_manual(values=c("black","red"))+
  facet_grid(.~CHROM, scales = "free", space = "free")+
  theme_bw()+
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 16),
                 axis.text.y = ggplot2::element_text(size = 16),
                 axis.title.x = ggplot2::element_text(size = 20, face = "bold", color = "black", vjust = -0.3), 
                 axis.title.y = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                 strip.text.x = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                 strip.text.y = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                 plot.title = ggplot2::element_text(size = 24, face = "bold", vjust = 1), 
                 panel.background = ggplot2::element_rect(color = "black",size = 1.2),
                 legend.position = "none")+
  labs(x = "Genomic Position (Mb)", y = "Test Statisitic")

ggsave(paste0(plot.dir,"GWA_PC1","_BURDEN_VTprice_manplot.pdf"), 
       height = 4, 
       width = 12)


# # # # SKAT
# process burden mappings
q90burdenSkat <- read.table(paste0(data.dir, "q90TOF.Skat.assoc"),header = T)
q90burdenSkat$Gene <- as.character(q90burdenSkat$Gene)
q90burdenSkat$RANGE <- as.character(q90burdenSkat$RANGE)

q90burdenSkat<-q90burdenSkat%>%
  dplyr::rowwise()%>%
  dplyr::mutate(CHROM = strsplit(strsplit(RANGE,split = ",")[[1]][1],split = ":")[[1]][1],
                POS = as.numeric(strsplit(strsplit(strsplit(RANGE,split = ",")[[1]][1],split = ":")[[1]][2],split = "-")[[1]][1]),
                endPOS = as.numeric(strsplit(strsplit(strsplit(RANGE,split = ",")[[1]][1],split = ":")[[1]][2],split = "-")[[1]][2]))%>%
  dplyr::mutate(size = abs(POS-endPOS))

skat <- dplyr::filter(q90burdenSkat,CHROM!="MtDNA")%>%
  dplyr::ungroup()%>%
  dplyr::filter(NumVar > 1)%>%
  dplyr::mutate(significant = ifelse(Pvalue < .05/n(), TRUE,FALSE ))

skat%>%
  # dplyr::filter(significant == TRUE)%>%
  ggplot()+
  aes(x = POS/1e6, y = -log10(Pvalue), alpha = 0.5, color = significant)+
  geom_point()+
  scale_color_manual(values=c("black","red"))+
  facet_grid(.~CHROM, scales = "free", space = "free")+
  theme_bw()+
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 16),
                 axis.text.y = ggplot2::element_text(size = 16),
                 axis.title.x = ggplot2::element_text(size = 20, face = "bold", color = "black", vjust = -0.3), 
                 axis.title.y = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                 strip.text.x = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                 strip.text.y = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                 plot.title = ggplot2::element_text(size = 24, face = "bold", vjust = 1), 
                 panel.background = ggplot2::element_rect(color = "black",size = 1.2),
                 legend.position = "none")+
  labs(x = "Genomic Position (Mb)", y = expression(-log[10](italic(p))))

ggsave(paste0(plot.dir,"GWA_",trait_of_interest,"_BURDEN_SKAT_manplot.pdf"), 
       height = 4, 
       width = 12)


# # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # # MANUAL CURATION OF ben-1 VARIANTS

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
  dplyr::mutate(ben1_prediction = ifelse(is.na(GT), NA, 
                                         ifelse(grepl("del", marker),"Deletion",
                                                ifelse(grepl("ins", marker,ignore.case = T),"Insertion",
                                                       ifelse(grepl("inv", marker),"Inversion",
                                                              ifelse(grepl("stop", marker),"Stop Gained",
                                                                     ifelse(grepl("trans", marker),"Transposon Insertion",
                                                                            ifelse(grepl("missense", marker),"Missense","Splice Donor"))))))))
ggplot(gwa_mappings)+
  aes(x = snpGT, y = value)+
  geom_boxplot(outlier.colour = NA)+
  facet_grid(.~snpMarker)+
  geom_jitter(shape = 21, color= "black", fill = "gray90", alpha = 0.5, data = dplyr::filter(gwa_mappings, is.na(marker)))+
  geom_jitter(aes(fill=ben1_prediction), data = na.omit(gwa_mappings), shape = 21, color= "black", size = 2)+
  scale_fill_manual(values=colors,name = expression(paste("Variation at ", italic("ben-1"))))+
  theme_bw()+
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 16),
                 axis.text.y = ggplot2::element_text(size = 16),
                 axis.title.x = ggplot2::element_text(size = 20, face = "bold", color = "black", vjust = -0.3), 
                 axis.title.y = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                 strip.text.x = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                 strip.text.y = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                 plot.title = ggplot2::element_text(size = 24, face = "bold", vjust = 1), 
                 panel.background = ggplot2::element_rect(color = "black",size = 1.2))+
  labs(x = "SNV Genotype at QTL", y = paste0("Animal Length (", trait_of_interest,")"))
  

ggsave(paste0(plot.dir,"GWA_",trait_of_interest,"_boxplot_colored_by_ben1.pdf"), 
       height = 4, 
       width = 12)

# # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # # Arranged bar plot
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


ggplot(gwa_mappings_bar)+
  aes(x=strain2, y=final_pheno, fill = factor(ben1_prediction,
                                              levels = sort(unique(ben1_prediction)), 
                                              labels= c("None", sort(unique(ben1_prediction)[2:8]) )))+
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("gray50",colors),name = expression(paste("Variation at ", italic("ben-1"))))+
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

ggsave(paste0(plot.dir,"GWA_",trait_of_interest,"_barplot_colored_by_ben1.pdf"), 
       height = 4, 
       width = 12)

ggplot(gwa_mappings_bar)+
  aes(x=strain2, y=final_pheno, fill = factor(ben1_prediction,
                                              levels = sort(unique(ben1_prediction)), 
                                              labels= c("None", sort(unique(ben1_prediction)[2:8]) )))+
  geom_bar(stat="identity" ) +
  scale_fill_manual(values=c("gray50","gray50","gray50","gray50","gray50","gray50","gray50","gray50","gray50"),name = expression(paste("Variation at ", italic("ben-1"))))+
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


ggsave(paste0(plot.dir,"GWA_",trait_of_interest,"_barplot_nocolors.pdf"), 
       height = 4, 
       width = 12)
# # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # # REGRESSION OF ben-1 variation and remapping

regress_ben1 <- gwa_mappings %>%
  dplyr::mutate(covariate = ifelse(is.na(GT) | strain == "CB4856", 0, 1))%>%
  dplyr::distinct(strain, value, covariate)

cor(regress_ben1$value, regress_ben1$covariate)^2
# [1] 0.7382114 variance explained by putative ben-1 LOF


residual_df <- data.frame(strain = regress_ben1$strain, regressed_trait = residuals(lm(value ~ covariate, data = regress_ben1)), ben1_variant_trait = regress_ben1$covariate)

pr_resid_pheno <- cegwas::process_pheno(residual_df)
resid_maps <- cegwas::gwas_mappings(pr_resid_pheno)
pr_resid_maps <- cegwas::process_mappings(resid_maps, phenotype_df = pr_resid_pheno, BF=5)

manhattan_plots <- manplot_edit(pr_resid_maps)

manhattan_plots[[1]]

ggsave(paste0(plot.dir,"GWA_ben1_variant_as_trait","_manplot.pdf"), 
       height = 4, 
       width = 12)

manhattan_plots[[2]]

ggsave(paste0(plot.dir,"GWA_ALL_ben1_variants_regressed_from_",trait_of_interest,"_manplot.pdf"), 
       height = 4, 
       width = 12)

write.table(pr_resid_maps, 
            file = paste0(data.dir,"GWAS_ben1_variants_regressed_from_",trait_of_interest,"_processed_mapping.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F)

genes <- process_correlations(variant_correlation(dplyr::filter(pr_resid_maps,trait =="regressed_trait"), condition_trait = F,variant_severity = "ALL"))

write.table(genes, 
            file = paste0(data.dir,"GWAS_finemapping_genes_ben1_regressed.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F)


# # # ## ben-1 regression burden test

traits <- residual_df%>%
  dplyr::select(-ben1_variant_trait)%>%
  dplyr::mutate(Fam = "elegans", Sample = strain, Paternal = 0, Maternal = 0, Sex = 2)%>%
  dplyr::select(-strain)

traits[is.na(traits)] <- -9

burden_raw <- dplyr::select(traits, Fam:Sex, regressed_trait)

write.table(burden_raw, 
            file = paste0(data.dir,"GWAS_Ben-1_Regressed.ped"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F)

benreg_burden <- read.table(paste0(data.dir, "GWAS_Ben-1_Regressed.VariableThresholdPrice.assoc"),header = T)
benreg_burden$Gene <- as.character(benreg_burden$Gene)
benreg_burden$RANGE <- as.character(benreg_burden$RANGE)

benreg_burden<-benreg_burden%>%
  dplyr::rowwise()%>%
  dplyr::mutate(CHROM = strsplit(strsplit(RANGE,split = ",")[[1]][1],split = ":")[[1]][1],
                POS = as.numeric(strsplit(strsplit(strsplit(RANGE,split = ",")[[1]][1],split = ":")[[1]][2],split = "-")[[1]][1]),
                endPOS = as.numeric(strsplit(strsplit(strsplit(RANGE,split = ",")[[1]][1],split = ":")[[1]][2],split = "-")[[1]][2]))%>%
  dplyr::mutate(size = abs(POS-endPOS))

benreg_burden_pr<- benreg_burden%>%
  ungroup()%>%
  dplyr::mutate(significant = ifelse(PermPvalue < .05/n(), TRUE,FALSE ))%>%
  dplyr::filter(CHROM!="MtDNA",NumVar>1,size >500)

save(benreg_burden_pr, file = (paste0(data.dir, "GWAS_Ben-1_Regressed.VariableThresholdPrice_processed.Rda")))

benreg_burden_pr%>%
  ggplot()+
  aes(x = POS/1e6, y = Stat, size = NumVar, alpha = 0.5, color = significant)+
  geom_point()+
  scale_color_manual(values=c("black","red"))+
  facet_grid(.~CHROM, scales = "free", space = "free")+
  theme_bw()+
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 16),
                 axis.text.y = ggplot2::element_text(size = 16),
                 axis.title.x = ggplot2::element_text(size = 20, face = "bold", color = "black", vjust = -0.3), 
                 axis.title.y = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                 strip.text.x = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                 strip.text.y = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                 plot.title = ggplot2::element_text(size = 24, face = "bold", vjust = 1), 
                 panel.background = ggplot2::element_rect(color = "black",size = 1.2),
                 legend.position = "none")+
  labs(x = "Genomic Position (Mb)", y = "Test Statisitic")

ggsave(paste0(plot.dir,"GWA_Ben1regressed_BURDEN_VTprice_manplot.png"), 
       height = 4, 
       width = 12)

# process burden mappings
q90regSkat <- read.table(paste0(data.dir, "GWAS_Ben-1_Regressed.Skat.assoc"),header = T)
q90regSkat$Gene <- as.character(q90regSkat$Gene)
q90regSkat$RANGE <- as.character(q90regSkat$RANGE)

q90regSkat<-q90regSkat%>%
  dplyr::rowwise()%>%
  dplyr::mutate(CHROM = strsplit(strsplit(RANGE,split = ",")[[1]][1],split = ":")[[1]][1],
                POS = as.numeric(strsplit(strsplit(strsplit(RANGE,split = ",")[[1]][1],split = ":")[[1]][2],split = "-")[[1]][1]),
                endPOS = as.numeric(strsplit(strsplit(strsplit(RANGE,split = ",")[[1]][1],split = ":")[[1]][2],split = "-")[[1]][2]))%>%
  dplyr::mutate(size = abs(POS-endPOS))

skat <- dplyr::filter(q90regSkat,CHROM!="MtDNA")%>%
  dplyr::ungroup()%>%
  dplyr::filter(NumVar > 1)%>%
  dplyr::mutate(significant = ifelse(Pvalue < .05/n(), TRUE,FALSE ))

skat%>%
  ggplot()+
  aes(x = POS/1e6, y = -log10(Pvalue), alpha = 0.5, color = significant)+
  geom_point()+
  scale_color_manual(values=c("black","red"))+
  facet_grid(.~CHROM, scales = "free", space = "free")+
  theme_bw()+
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 16),
                 axis.text.y = ggplot2::element_text(size = 16),
                 axis.title.x = ggplot2::element_text(size = 20, face = "bold", color = "black", vjust = -0.3), 
                 axis.title.y = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                 strip.text.x = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                 strip.text.y = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                 plot.title = ggplot2::element_text(size = 24, face = "bold", vjust = 1), 
                 panel.background = ggplot2::element_rect(color = "black",size = 1.2),
                 legend.position = "none")+
  labs(x = "Genomic Position (Mb)", y = expression(-log[10](italic(p))))

ggsave(paste0(plot.dir,"GWA_Ben-1_Regressed_BURDEN_SKAT_manplot.pdf"), 
       height = 4, 
       width = 12)

# # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # # # # # # ben-1 haplotype plot


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

# pr_plot_df <- group_by(plot_df, .groups = colnames(plot_df[3:ncol(plot_df)]))

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
  scale_y_discrete(labels = scales::wrap_format(25))+
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
  # geom_vline(aes(xintercept=c(3541630)), linetype="dotdash", alpha = .25, color = "blue")+
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
  # dplyr::filter(grepl("missense",marker), start < 3542000)%>%
  dplyr::distinct(aa_change, .keep_all =T)%>%
  dplyr::mutate(plot_aa = gsub(pattern = "p\\.","",aa_change))

missense_df$letter <- c("D404N", "M257I", "F200Y", "A185P","S145F","Q131L","E69G","K58E")

missense_df_low <- dplyr::filter(missense_df, letter %in% c("M257I","A185P","Q131L","K58E"))
missense_df_high <- dplyr::filter(missense_df, !(letter %in% c("M257I","A185P","Q131L","K58E")))

ben1_model <- gene_model(ben1_gene_info)+
  ylim(-4,4)+
  geom_segment(aes(x = as.numeric(POS), y = -1, xend = as.numeric(POS), yend = -2),data = missense_df_low, size =2)+
  geom_segment(aes(x = as.numeric(POS), y = 1, xend = as.numeric(POS), yend = 2),data = missense_df_high, size =2)+
  # ggrepel::geom_label_repel(aes(x = POS, y = -1.1, label = letter), data = missense_df_low,size = 2)+
  # ggrepel::geom_label_repel(aes(x = POS, y = 1.1, label = letter), data = missense_df_high,size = 2)
  geom_label(aes(x = POS, y = -3, label = letter), data = missense_df_low,size = 2)+
  geom_label(aes(x = POS, y = 3, label = letter), data = missense_df_high,size = 2)




cowplot::ggdraw() +
  # cowplot::draw_plot(pxg_pl, .75, 0.069, .25, .73)+
  cowplot::draw_plot(ben1_high_var, 0, 0, 1, .8)+
  cowplot::draw_label(label = expression(bolditalic('ben-1')), 0.5, y = .92)+
  cowplot::draw_plot(ben1_model, .13, .81, .52, .1)

ggsave(paste0(plot.dir,"Ben1_haps_by_phenotype.pdf"), 
       height = 12, 
       width = 18)

# # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # # ben-1 phylo

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
control_oof_interest <- "DMSO"

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


write.table(gwa_mappings, paste0(data.dir, "Processed_Ben1_variants_with_pheno.tsv"), quote = F, col.names = T, row.names = F)

# Load admixture populations
base.dir <- "~/Dropbox/AndersenLab/LabFolders/Stefan/indels/Post_calling_analysis/"

load( file = paste0(base.dir, "Data/snp_distance_matrix.Rda"))

snps_NJ <- NJ(dm_snps)

tree <- ggtree(snps_NJ, 
               layout = "circular", branch.length = "none") 

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



tree %<+% ben1_strain_markers+ 
  geom_tiplab2(aes(color=ben1_prediction))+
  scale_color_manual(values=colors,name = expression(paste("Variation at ", italic("ben-1"))))+
  # theme(plot.margin=unit(c(1.5,1.5,1.5,1.5),"cm"))+
  theme(legend.position=c(0.9, 0.8))

ggsave(paste0(plot.dir,"Ben1_phylo.pdf"), 
       height = 18, 
       width = 18)


# # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # # World map
library(ggplot2)  # FYI you need v2.0
library(dplyr)    # yes, i could have not done this and just used 'subset' instead of 'filter'
library(ggalt)    # devtools::install_github("hrbrmstr/ggalt")
library(ggthemes) # theme_map and tableau colors

world <- map_data("world")
world <- world[world$region != "Antarctica",] # intercourse antarctica


  
isolation_info <-readr::read_tsv("https://elegansvariation.org/strain/strain_data.tsv")%>%
  dplyr::filter(reference_strain == "True")%>%
  dplyr::select(strain = isotype, long = longitude, lat = latitude)%>%
  dplyr::filter(lat != "None")%>%
  dplyr::left_join(ben1_strain_markers,.,by="strain")%>%
  dplyr::filter(!is.na(lat))

isolation_info$lat <- as.numeric(isolation_info$lat)
isolation_info$long <- as.numeric(isolation_info$long)

ggplot()+ geom_map(data=world, map=world,
                   aes(x=long, y=lat, map_id=region),
                   color="white", fill="#7f7f7f", size=0.05, alpha=1)+
  geom_point(data=dplyr::filter(isolation_info, !is.na(marker)),  
             aes(x=long, y=lat, color=ben1_prediction)) + 
  scale_color_manual(values=colors,name = expression(paste("Variation at ", italic("ben-1"))))+
  ggrepel::geom_text_repel(data =dplyr::filter(isolation_info, !is.na(marker)),
                              aes(x=long, y=lat, label = strain),
                              point.padding = 0.001) + 
  theme_map()


ggsave(paste0(plot.dir,paste0("World_MAP_ben-1_variants_labels.pdf")), 
       height = 10, 
       width = 20)

ggplot()+ geom_map(data=world, map=world,
                   aes(x=long, y=lat, map_id=region),
                   color="white", fill="#7f7f7f", size=0.05, alpha=1)+
  geom_point(data=dplyr::filter(isolation_info, !is.na(marker)),  
             aes(x=long, y=lat, fill=ben1_prediction), size = 5, alpha = 0.75, shape = 21)+
  scale_fill_manual(values=colors,name = expression(paste("Variation at ", italic("ben-1"))))+
  theme_map()


ggsave(paste0(plot.dir,paste0("World_MAP_ben-1_variants_no-labels.pdf")), 
       height = 10, 
       width = 20)

ggsave(paste0(plot.dir,paste0("World_MAP_ben-1_variants_no-labels.png")), 
       height = 10, 
       width = 20,
       dpi = 300)


ggplot() +
  geom_map(data=world, map=world,
           aes(x=long, y=lat, map_id=region),
           color="white", fill="#7f7f7f", size=0.05, alpha=1)+
    geom_point(data=dplyr::filter(isolation_info, !is.na(marker)),  
             aes(x=long, y=lat, fill=ben1_prediction), size = 5, alpha = 0.75, shape = 21) +
  coord_cartesian(xlim = c(-9,20), ylim = c(32,60))+
  scale_fill_manual(values=colors,name = expression(paste("Variation at ", italic("ben-1"))))+
  theme_map()

ggsave(paste0(plot.dir,paste0("World_MAP_ben-1_variants_Europe_no-labels.png")), 
       height = 10, 
       width = 10,
       dpi = 300)

# # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # # beta tubulin tajimas D

tajimas_d_temp <- function (vcf_path = paste0(path.package("cegwas"), "/"), vcf_name = paste0("WI.", 
                                                                                              vcf_version, ".impute.vcf.gz"), chromosome = "II", interval_start = 11021073, 
                            interval_end = 12008179, window_size = 300, slide_distance = 100, 
                            samples = sort(colnames(snps[, 5:ncol(snps)])), outgroup = "XZ1516", 
                            site_of_interest = 11875145) 
{
  setwd(vcf_path)
  gen <- PopGenome::readVCF(vcf_name, numcols = 10000, tid = chromosome, 
                            frompos = interval_start, topos = interval_end, samplenames = samples, 
                            approx = T)
  gen1 <- PopGenome::set.populations(gen, list(samples), diploid = TRUE)
  gen2 <- PopGenome::set.outgroup(gen1, outgroup, diploid = TRUE)
  s_gen <- PopGenome::sliding.window.transform(gen2, width = window_size, 
                                               jump = slide_distance, whole.data = FALSE)
  test <- data.frame(snps = 1:length(as.numeric(colnames(as.data.frame(s_gen@BIG.BIAL[[1]][, 
                                                                                           ])))), position = as.numeric(colnames(as.data.frame(s_gen@BIG.BIAL[[1]][, 
                                                                                                                                                                   ]))))
  slides <- list()
  for (i in 1:length(s_gen@SLIDE.POS)) {
    slides[[i]] <- data.frame(window = i, snps = s_gen@SLIDE.POS[[i]])
  }
  slides <- dplyr::rbind_all(slides) %>% dplyr::left_join(data.frame(test), 
                                                          ., by = "snps")
  s_gen_stats <- PopGenome::neutrality.stats(s_gen)
  td <- data.frame(window = 1:length(s_gen@SLIDE.POS), Td = s_gen_stats@Tajima.D) %>% 
    dplyr::left_join(slides, ., by = "window") %>% dplyr::group_by(window) %>% 
    dplyr::mutate(swindow = min(as.numeric(position)), ewindow = max(as.numeric(position))) %>% 
    dplyr::mutate(midwindow = (swindow + ewindow)/2) %>% 
    dplyr::rename(Td = pop.1) %>% dplyr::distinct(Td, window, 
                                                  .keep_all = T)
  tajimas_d_plot <- ggplot2::ggplot(td) + ggplot2::aes(x = position/1e+06, 
                                                       y = Td) + ggplot2::geom_point(size = 2) + ggplot2::theme_bw() + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 24, 
                                                       face = "bold", color = "black"), axis.text.y = ggplot2::element_text(size = 24, 
                                                                                                                            face = "bold", color = "black"), axis.title.x = ggplot2::element_text(size = 24, 
                                                                                                                                                                                                  face = "bold", color = "black", vjust = -0.3), axis.title.y = ggplot2::element_text(size = 24, 
                                                                                                                                                                                                                                                                                      face = "bold", color = "black"), strip.text.x = ggplot2::element_text(size = 24, 
                                                                                                                                                                                                                                                                                                                                                            face = "bold", color = "black"), strip.text.y = ggplot2::element_text(size = 16, 
                                                                                                                                                                                                                                                                                                                                                                                                                                  face = "bold", color = "black"), plot.title = ggplot2::element_text(size = 24, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      face = "bold", vjust = 1), legend.position = "none", 
                   panel.background = ggplot2::element_rect(color = "black", 
                                                            size = 1.2), strip.background = ggplot2::element_rect(color = "black", 
                                                                                                                  size = 1.2)) + ggplot2::labs(x = "Genomic Position (Mb)", 
                                                                                                                                               y = "Tajima's D")
  if (!is.null(site_of_interest)) {
    tajimas_d_plot <- tajimas_d_plot + ggplot2::geom_vline(ggplot2::aes(xintercept = site_of_interest/1e+06), 
                                                           color = "red", alpha = 0.7, size = 2)
  }
  return(list(td, tajimas_d_plot))
}


library(cegwas)
library(dplyr)
library(ggplot2)
ben1_d <- tajimas_d_temp(vcf_path = "~/Dropbox/AndersenLab/LabFolders/Stefan/vcfs/",vcf_name = "gatk_snpEff_diploid.vcf.gz", 
                         chromosome = "III", interval_start = 3437688, interval_end = 3641628,site_of_interest = 3539628,slide_distance = 1,window_size = 100)
ben1_d

ggsave("~/Dropbox/HTA/Results/2017_anthelmintic_gwa/analysis/ben1_paper/Plots/TajimaD_ben1_100window_1slide.pdf",
       height = 5,
       width = 12)


# tbb1
tbb1_d <- tajimas_d_temp(vcf_path = "~/Dropbox/AndersenLab/LabFolders/Stefan/vcfs/",vcf_name = "gatk_snpEff_diploid.vcf.gz", 
                         chromosome = "III", interval_start = 10730868, interval_end = 10782461,site_of_interest = 10741461,slide_distance = 1,window_size = 100)
tbb1_d
ggsave("~/Dropbox/HTA/Results/2017_anthelmintic_gwa/analysis/ben1_paper/Plots/TajimaD_tbb1_100window_1slide.pdf",
       height = 5,
       width = 12)
# tbb2
tbb2 <- tajimas_d_temp(vcf_path = "~/Dropbox/AndersenLab/LabFolders/Stefan/vcfs/",vcf_name = "gatk_snpEff_diploid.vcf.gz",
                       chromosome = "III", interval_start = 3975769, interval_end = 4059643,site_of_interest = 4016643,slide_distance = 1,window_size = 100) 
tbb2
ggsave("~/Dropbox/HTA/Results/2017_anthelmintic_gwa/analysis/ben1_paper/Plots/TajimaD_tbb2_100window_1slide.pdf",
       height = 5,
       width = 12)

# tbb4 X:9434452..9436893
tbb4 <- tajimas_d_temp(vcf_path = "~/Dropbox/AndersenLab/LabFolders/Stefan/vcfs/",vcf_name = "gatk_snpEff_diploid.vcf.gz",
                       chromosome = "X", interval_start = 9334452, interval_end = 9536893,site_of_interest = 9436893,slide_distance = 1,window_size = 100) 
tbb4
ggsave("~/Dropbox/HTA/Results/2017_anthelmintic_gwa/analysis/ben1_paper/Plots/TajimaD_tbb4_100window_1slide.pdf",
       height = 5,
       width = 12)
# tbb6 V: 12261804 12263368
tbb6 <- tajimas_d_temp(vcf_path = "~/Dropbox/AndersenLab/LabFolders/Stefan/vcfs/",vcf_name = "gatk_snpEff_diploid.vcf.gz",
                       chromosome = "V", interval_start = 12161804, interval_end = 12363368,site_of_interest = 12263368,slide_distance = 1,window_size = 100) 
tbb6
ggsave("~/Dropbox/HTA/Results/2017_anthelmintic_gwa/analysis/ben1_paper/Plots/TajimaD_tbb6_100window_1slide.pdf",
       height = 5,
       width = 12)
# mec7  X:7774859 7776689
mec7 <- tajimas_d_temp(vcf_path = "~/Dropbox/AndersenLab/LabFolders/Stefan/vcfs/",vcf_name = "gatk_snpEff_diploid.vcf.gz",
                       chromosome = "X", interval_start = 7674859, interval_end = 7876689,site_of_interest = 7776689,slide_distance = 1,window_size = 100) 
mec7
ggsave("~/Dropbox/HTA/Results/2017_anthelmintic_gwa/analysis/ben1_paper/Plots/TajimaD_mec7_100window_1slide.pdf",
       height = 5,
       width = 12)

# Ben-1 locus popgenome

library(PopGenome)
vcf_name = "ben1_delly_strelka.vcf.gz"
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

gen <- PopGenome::readVCF(vcf_name, 
                          numcols = 1000, 
                          tid = chroms[w.CHROM], 
                          frompos = chr.lengths[[w.CHROM]][1], 
                          topos = 3550315, 
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

save(s_gen_stats,file = paste0(data.dir,"ben1_popgenome.Rda"))


# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # 
# F200Y vs ben-1 deletion crispr
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
# ASSAY A
raw <- read_data(paste0(data.dir,"20180312_ben1_edits_A/"))
raw_nocontam <- remove_contamination(raw)
summedraw <- sumplate(raw_nocontam, quantiles = TRUE)

write.table(summedraw, 
            file = paste0(data.dir,"F200Y_del_sorter_summarized_A.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F)

summedraw_A <- data.table::fread(paste0(data.dir,"F200Y_del_sorter_summarized_A.tsv"))

biopruned_A <- bioprune(summedraw_A)%>%
  tidyr::gather(trait, value, -date, -experiment, -round, -assay, -plate, -condition, -control, -strain, -row, -col )%>%
  dplyr::filter(!is.na(condition))%>%
  dplyr::group_by(strain, assay, condition, trait)%>%
  dplyr::mutate(mph = median(value),
                sph = sd(value))%>%
  dplyr::mutate(flag_h = 2*sph+mph,
                flag_l = mph-2*sph)%>%
  dplyr::mutate(cut_h =ifelse(value >= 2*sph+mph, "YES", "NO"),
                cut_l =ifelse(value <= mph-2*sph, "YES", "NO"))%>%
  dplyr::filter(cut_h != "YES" , cut_l !="YES")%>%
  dplyr::ungroup()%>%
  dplyr::select(-cut_l, -cut_h,-flag_l,-flag_h,-sph,-mph,-date,-experiment,-round)%>%
  dplyr::rename(phenotype = value)

biopruned_A <- biopruned_A%>%
  dplyr::mutate(assay = "A")

# ASSAY B
raw <- read_data(paste0(data.dir,"20180313_ben1_edits_B/"))
raw_nocontam <- remove_contamination(raw)
summedraw <- sumplate(raw_nocontam, quantiles = TRUE)

write.table(summedraw, 
            file = paste0(data.dir,"F200Y_del_sorter_summarized_B.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F)

summedraw_B <- data.table::fread(paste0(data.dir,"F200Y_del_sorter_summarized_B.tsv"))

biopruned_B <- bioprune(summedraw_B)%>%
  tidyr::gather(trait, value, -date, -experiment, -round, -assay, -plate, -condition, -control, -strain, -row, -col )%>%
  dplyr::filter(!is.na(condition))%>%
  dplyr::group_by(strain, assay, condition, trait)%>%
  dplyr::mutate(mph = median(value),
                sph = sd(value))%>%
  dplyr::mutate(flag_h = 2*sph+mph,
                flag_l = mph-2*sph)%>%
  dplyr::mutate(cut_h =ifelse(value >= 2*sph+mph, "YES", "NO"),
                cut_l =ifelse(value <= mph-2*sph, "YES", "NO"))%>%
  dplyr::filter(cut_h != "YES" , cut_l !="YES")%>%
  dplyr::ungroup()%>%
  dplyr::select(-cut_l, -cut_h,-flag_l,-flag_h,-sph,-mph,-date,-experiment,-round)%>%
  dplyr::rename(phenotype = value)

biopruned_B <- biopruned_B%>%
  dplyr::mutate(assay = "B")

# combine dataframes and remove any remaining wash wells 
all_data <- dplyr::bind_rows(list(biopruned_A,biopruned_B))%>%
  dplyr::filter(!is.na(strain))

#### control regression ######

assayregressed <- easysorter::regress(all_data, assay = TRUE)%>%
  dplyr::group_by(assay,strain, condition, trait)%>%
  dplyr::mutate(mph = median(phenotype),
                sph = sd(phenotype))%>%
  dplyr::mutate(flag_h = 2*sph+mph,
                flag_l = mph-2*sph)%>%
  dplyr::mutate(cut_h =ifelse(phenotype >= 2*sph+mph, "YES", "NO"),
                cut_l =ifelse(phenotype <= mph-2*sph, "YES", "NO"))%>%
  dplyr::filter(cut_h != "YES" , cut_l !="YES")

pr_df <- assayregressed%>%
  dplyr::ungroup()%>%
  dplyr::filter(trait == "q90.TOF", grepl("alb",condition,ignore.case = T))%>%
  dplyr::filter(!strain =="PTM")%>%
  dplyr::filter(!is.na(condition))%>%
  dplyr::mutate(genotype = ifelse(strain %in% c("ECA880","ECA881","ECA917","ECA918","ECA919","ECA920","ECA921"), "F200Y",
                                  "del"))%>%
  dplyr::mutate(group = ifelse(strain %in% c("N2"), "N2",
                               ifelse(strain %in% c("ECA880","ECA881"), "PAM",
                                      ifelse(strain %in% c("ECA917","ECA918","ECA919","ECA920","ECA921"),"F200Y",
                                             "del"))))

write.table(pr_df, 
            file = paste0(data.dir,"F200Y_del_sorter_assay_regressed.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F)



