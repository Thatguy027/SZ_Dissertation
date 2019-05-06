# set to location of files
main.dir <- "~/Dropbox/AndersenLab/LabFolders/Stefan/Manuscripts/Arsenic/"
data.dir <- paste0(main.dir,"Data/")
final.dir <- paste0(main.dir,"Final_tables/")
plot.dir <- paste0(main.dir,"Plots/")


resub.dir <- paste0(main.dir,"Resubmission/")

resub.data <- paste0(resub.dir,"Data/")
resub.plots <- paste0(resub.dir,"Plots/")
resub.gwas <- paste0(resub.dir,"cegwas2_GWA_Results/")


dir.create(resub.data)
dir.create(resub.plots)

script.dir <- paste0(resub.dir,"Scripts/")

source(paste0(script.dir,"20181009_processing_functions.R"))

########################################################################################################################
# PROCESS DOSE RESPONSE EXPERIMENT - START
########################################################################################################################

dirs <- paste0(data.dir,"DoseResponse_new/")

# easysorter analysis
raw <- read_data(dirs)
raw_nocontam <- remove_contamination(raw)
summedraw <- sumplate(raw_nocontam, directories = FALSE, quantiles = TRUE)
biopruned <- bioprune(summedraw)

# subtract phenotypes from mean of controls
biopruned_pr <- biopruned %>%
  dplyr::ungroup() %>%
  dplyr::select(-(date:assay)) %>%
  dplyr::filter(!is.na(condition)) %>%
  tidyr::gather(trait, value, -(plate:col))

control_dr <- biopruned_pr %>%
  dplyr::filter(condition == "Water") %>%
  dplyr::group_by(strain, trait) %>%
  dplyr::summarise(control_value = mean(value, na.rm = T))

# RUN PCA ON ENTIRE EXPERIMENT
subtract_dr_to_pc <- biopruned_pr %>%
  dplyr::left_join(., control_dr, by = c("strain", "trait")) %>%
  dplyr::mutate(delta_control = value - control_value) %>%
  dplyr::select(strain, condition, trait, delta_control, plate, row, col)%>%
  dplyr::mutate(u_strain = paste(strain, condition, plate,row,col,sep="_"))%>%
  tidyr::spread(trait, delta_control) %>%
  dplyr::select(-contains("green"),-contains("yellow"), -contains("red"),-contains("cv"),-contains("iqr"), -contains("var"))

dr_arsenic_pca <- data.frame(scale(subtract_dr_to_pc[,7:ncol(subtract_dr_to_pc)]))
dr_arsenic_pca_analysis <- prcomp(dr_arsenic_pca)

pc_ve_90 <- which(cumsum(dr_arsenic_pca_analysis$sdev^2/sum(dr_arsenic_pca_analysis$sdev^2)) > 0.9)[1]

dr_complete_pc_df <- data.frame(strain = subtract_dr_to_pc$strain,
                                condition = subtract_dr_to_pc$condition,
                                plate = subtract_dr_to_pc$plate,
                                row = subtract_dr_to_pc$row,
                                col = subtract_dr_to_pc$col,
                                dr_arsenic_pca_analysis$x[,1:pc_ve_90]) %>%
  dplyr::rename(Strain = strain, Condition = condition) %>%
  dplyr::select(-(plate:col)) %>%
  tidyr::gather(Trait, Value, -Strain, -Condition)

# extract PC loadings
dr_complete_rotation_df <- data.frame(trait = row.names(dr_arsenic_pca_analysis$rotation),
                                      loadings = dr_arsenic_pca_analysis$rotation[,1:pc_ve_90])

colnames(dr_complete_rotation_df) <- gsub("loadings.", "", colnames(dr_complete_rotation_df))

dr_complete_rotation_df <- dr_complete_rotation_df %>%
  tidyr::gather(PC, loading, -trait) 

# RUN PCA FOR EACH CONDITION

subtract_dr <- biopruned_pr %>%
  dplyr::left_join(., control_dr, by = c("strain", "trait")) %>%
  dplyr::mutate(delta_control = value - control_value) %>%
  dplyr::select(strain, condition, trait, delta_control, plate, row, col)%>%
  dplyr::mutate(u_strain = paste(strain,condition,plate,row,col,sep="_"))%>%
  tidyr::spread(trait,delta_control)

trait_list <- list()
cor_list <- list()
loading_list <- list()

for(cond in 1:length(unique(subtract_dr$condition))){
  
  analysis_condition <- unique(subtract_dr$condition)[cond]
  
  # run PCA on DR traits
  cond_df <- subtract_dr %>%
    dplyr::filter(condition == analysis_condition) %>%
    dplyr::select(-contains("cv"),-contains("iqr"), -contains("var"))
  
  dr_arsenic_pca <- data.frame(scale(cond_df[,7:ncol(cond_df)]))
  dr_arsenic_pca_analysis <- prcomp(dr_arsenic_pca)
  
  pc_ve_90 <- which(cumsum(dr_arsenic_pca_analysis$sdev^2/sum(dr_arsenic_pca_analysis$sdev^2)) > 0.9)[1]
  
  dr_pc_df <- data.frame(strain = cond_df$strain,
                         condition = cond_df$condition,
                         plate = cond_df$plate,
                         row = cond_df$row,
                         col = cond_df$col,
                         dr_arsenic_pca_analysis$x[,1:pc_ve_90])
  
  # combine PC with raw traits
  dr_all_traits <- dplyr::left_join(cond_df, dr_pc_df, by = c("strain", "plate", "condition","row", "col")) %>%
    dplyr::select(-u_strain, -plate, -row, -col)
  
  # extract PC loadings
  dr_rotation_df <- data.frame(trait = row.names(dr_arsenic_pca_analysis$rotation),
                               loadings = dr_arsenic_pca_analysis$rotation[,1:pc_ve_90])
  
  colnames(dr_rotation_df) <- gsub("loadings.", "", colnames(dr_rotation_df))
  
  dr_rotation_df <- dr_rotation_df %>%
    tidyr::gather(PC, loading, -trait) %>%
    dplyr::mutate(Condition = analysis_condition) %>%
    dplyr::rename(Trait = trait)
  
  # Find trait correlations
  trait_correlations <- cor(dr_all_traits[,3:ncol(dr_all_traits)])%>%
    data.frame() 
  
  trait_correlations$trait_b <- row.names(trait_correlations)
  
  dr_trait_correlations <- trait_correlations%>%
    dplyr::mutate(Condition = analysis_condition) %>%
    tidyr::gather(trait_a, trait_cor, -trait_b, -Condition) 
  
  # make DR trait data long
  dr_all_traits_long <- dr_all_traits %>%
    tidyr::gather(Trait, Value, -condition, -strain) %>%
    dplyr::rename(Strain = strain, Condition = condition)
  
  trait_list[[cond]] <- dr_all_traits_long
  cor_list[[cond]] <- dr_trait_correlations
  loading_list[[cond]] <- dr_rotation_df
}

dr_all_traits_long <- dplyr::bind_rows(trait_list)
dr_trait_correlations <- dplyr::bind_rows(cor_list)
dr_rotation_df <- dplyr::bind_rows(loading_list)

########################################################################################################################
# PROCESS DOSE RESPONSE EXPERIMENT - END
########################################################################################################################

########################################################################################################################
# DOSE RESPONSE EXPERIMENT - EFFECT SIZES - START
########################################################################################################################

# Calculate effect sizes of traits

dr_var_part <- dr_all_traits_long%>%
  dplyr::filter(!grepl("PC",Trait)) %>%
  dplyr::bind_rows(.,dr_complete_pc_df) %>%
  dplyr::mutate(cond_trait = paste0(Condition, "_",Trait))

effect_sizes <- list()
for(cond_tr in 1:length(unique(dr_var_part$cond_trait))){
  
  test_cond_trait <- unique(dr_var_part$cond_trait)[cond_tr]
  test_trait <- dr_var_part %>%
    dplyr::filter(cond_trait == test_cond_trait)
  
  lm_dr <- lm(Value ~ Strain, data = test_trait)
  complete_condtr <- data.frame(anova_stats(lm_dr),
                                condition = strsplit(test_cond_trait,split = "_")[[1]][1],
                                trait = strsplit(test_cond_trait,split = "_")[[1]][2],
                                strain1 = "ALL",
                                strain2 = "ALL")
  
  combos <- combn(unique(test_trait$Strain),m = 2)
  
  c_ls <- list()
  for(s in 1:ncol(combos)){
    two_compare <- test_trait %>%
      dplyr::filter(Strain %in% combos[,s])
    
    lm_dr <- lm(Value ~ Strain, data = two_compare)
    c_ls[[s]] <- data.frame(anova_stats(lm_dr),
                            condition = strsplit(test_cond_trait,split = "_")[[1]][1],
                            trait = strsplit(test_cond_trait,split = "_")[[1]][2],
                            strain1 = combos[1,s],
                            strain2 = combos[2,s])
  }
  
  effect_sizes[[cond_tr]] <- dplyr::bind_rows(c_ls) %>%
    dplyr::bind_rows(.,complete_condtr) %>%
    dplyr::filter(term!="Residuals", term!="(Intercept)", condition != "Water")
}

effects_df <- dplyr::bind_rows(effect_sizes) %>%
  dplyr::filter(term!="Residuals", term!="(Intercept)", condition != "Water")

########################################################################################################################
# DOSE RESPONSE EXPERIMENT - EFFECT SIZES - END
########################################################################################################################

########################################################################################################################
# DOSE RESPONSE EXPERIMENT - HERITABILITY ESTIMATES - START
########################################################################################################################

# Calculate heritability of DR data
pr_dr <- dr_all_traits_long %>%
  dplyr::group_by(Strain, Condition, Trait) %>%
  dplyr::filter(Value > (mean(Value, na.rm = T)) - 2*sd(Value, na.rm = T),
                Value < (mean(Value, na.rm = T)) + 2*sd(Value, na.rm = T))%>%
  dplyr::ungroup() %>%
  dplyr::mutate(cond_trait = paste0(Condition, "_",Trait))

herits <- list()

for(cond_tr in 1:length(unique(pr_dr$cond_trait))){
  
  test_cond_trait <- unique(pr_dr$cond_trait)[cond_tr]
  
  for_calc <- pr_dr %>%
    dplyr::filter(cond_trait == test_cond_trait) %>%
    dplyr::select(strain = Strain, value = Value)
  
  reffMod = lme4::lmer(value ~ (1|strain), data=for_calc)
  
  Var_Random_effect <- as.numeric(lme4::VarCorr(reffMod))
  Var_Residual <- attr(lme4::VarCorr(reffMod), "sc")^2
  H2 <- Var_Random_effect/(Var_Random_effect+Var_Residual)
  cVg <- (100*sqrt(Var_Random_effect))/mean(reffMod@resp$mu)
  vVe <- (100*sqrt(Var_Residual))/mean(reffMod@resp$mu)
  
  herits[[cond_tr]] <- data.frame(DR_H2 = H2,
                                  Strain_Effect = Var_Random_effect,
                                  Residual_Var = Var_Residual,
                                  Coef_Gv = cVg,
                                  Coef_Rv = vVe,
                                  condition = strsplit(test_cond_trait,split = "_")[[1]][1],
                                  trait = strsplit(test_cond_trait,split = "_")[[1]][2])
  
}

broad_h2_df <- dplyr::bind_rows(herits)

########################################################################################################################
# DOSE RESPONSE EXPERIMENT - HERITABILITY ESTIMATES - END
########################################################################################################################

########################################################################################################################
# SAVE DOSE RESPONSE EXPERIMENT DATA
########################################################################################################################
# save DR trait data
write.table(dr_all_traits_long,
            file = glue::glue("{resub.data}Supplemental_Data1_DR-TRAITS_per_Condition.tsv"),
            sep = "\t",
            quote = F,
            col.names = T,
            row.names = F)

# Save DR H2 estimates
write.table(broad_h2_df,
            file = glue::glue("{resub.data}Supplemental_Data2_DR-Heritability_estimates.tsv"),
            sep = "\t",
            quote = F,
            col.names = T,
            row.names = F)

# save DR Effect size data
write.table(effects_df,
            file = glue::glue("{resub.data}Supplemental_Data3_DR-PCA_Effect-Sizes.tsv"),
            sep = "\t",
            quote = F,
            col.names = T,
            row.names = F)

# save DR PC loadings
write.table(dr_rotation_df,
            file = glue::glue("{resub.data}Supplemental_Data4_DR-PCA-LOADINGS_Per_Condition.tsv"),
            sep = "\t",
            quote = F,
            col.names = T,
            row.names = F)

# save DR trait data
write.table(dr_complete_pc_df,
            file = glue::glue("{resub.data}Supplemental_Data5_DR-PCA-TRAITS_Complete_Experiment.tsv"),
            sep = "\t",
            quote = F,
            col.names = T,
            row.names = F)

# save DR trait data
write.table(dr_complete_rotation_df,
            file = glue::glue("{resub.data}Supplemental_Data6_DR-PCA-LOADINGS_Complete_Experiment.tsv"),
            sep = "\t",
            quote = F,
            col.names = T,
            row.names = F)

# save DR correlation data
write.table(dr_trait_correlations,
            file = glue::glue("{resub.data}Supplemental_Data7_DR-TRAIT_COR.tsv"),
            sep = "\t",
            quote = F,
            col.names = T,
            row.names = F)
########################################################################################################################
# SAVE DOSE RESPONSE EXPERIMENT DATA
########################################################################################################################

dr_traits <- data.table::fread(glue::glue("{resub.data}Supplemental_Data1_DR-TRAITS_per_Condition.tsv"))
dr_trait_cor <- data.table::fread(glue::glue("{resub.data}Supplemental_Data7_DR-TRAIT_COR.tsv"))
dr_loadings <- data.table::fread(glue::glue("{resub.data}Supplemental_Data4_DR-PCA-LOADINGS_Per_Condition.tsv"))

dr_complete_pcs <- data.table::fread(glue::glue("{resub.data}Supplemental_Data5_DR-PCA-TRAITS_Complete_Experiment.tsv"))
dr_complete_loadings <- data.table::fread(glue::glue("{resub.data}Supplemental_Data6_DR-PCA-LOADINGS_Complete_Experiment.tsv"))

effects_df <- data.table::fread(glue::glue("{resub.data}Supplemental_Data3_DR-PCA_Effect-Sizes.tsv"))
h2_df <- data.table::fread(glue::glue("{resub.data}Supplemental_Data2_DR-Heritability_estimates.tsv"))
########################################################################################################################
# DOSE RESPONSE EXPERIMENT - PLOTS - START
########################################################################################################################
# plot dose response
dr_plt_df <- dr_traits%>%
  dplyr::filter(Trait %in% c("median.TOF","norm.n")) %>%
  dplyr::bind_rows(.,dr_complete_pcs) %>%
  dplyr::filter(Trait %in% c("median.TOF","norm.n", "PC1")) %>%
  dplyr::mutate(flip_pc = ifelse(Trait == "PC1", -Value, Value)) %>%
  dplyr::mutate(Tidy_trait = ifelse(Trait == "median.TOF", "Animal Size", 
                                    ifelse(Trait == "norm.n"," Brood Size", "PC 1")),
                tidy_cond = factor(Condition, 
                                   levels = c("Water", "arsenic1", "arsenic2", "arsenic3","arsenic4"),
                                   labels = c("0", "250", "500", "1000","2000"))) 

dr_plots <- list()
for(condit in 1:length(unique(dr_plt_df$Tidy_trait))){
  
  
  dr_plots[[condit]] <- dr_plt_df%>%
    dplyr::filter(Tidy_trait == unique(dr_plt_df$Tidy_trait)[condit]) %>%
    ggplot()+
    aes(x = tidy_cond, 
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
}

legend <- cowplot::get_legend(dr_plots[[3]])
dr_plots <- cowplot::plot_grid(dr_plots[[2]]+theme(axis.title.x = element_blank(),
                                                   axis.text.x = element_blank(),
                                                   legend.position = "none"), 
                               dr_plots[[1]]+theme(axis.title.x = element_blank(),
                                                   axis.text.x = element_blank(),
                                                   legend.position = "none"), 
                               dr_plots[[3]]+theme(legend.position = "none"),
                               labels = "AUTO", 
                               ncol = 1, 
                               align = 'v',
                               vjust = 1, 
                               label_size = 18)

plot_grid( dr_plots, legend, rel_widths = c(3, .3))

ggsave(glue::glue("{resub.plots}FS1_DR.pdf"), 
       height = 10, 
       width = 10)

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
  
  ggsave(glue::glue("{resub.plots}FS2{tidy_subplot}_Trait-correlation_loadings.pdf"), 
         height = 12, 
         width = 20)
  
}

# plot effect sizes for DR
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

ggsave(glue::glue("{resub.plots}FS3_DR_H2_v_Effect_Sizes.pdf"), 
       height = 8, 
       width = 16)

########################################################################################################################
# PROCESS LINKAGE MAPPING - START
########################################################################################################################

#Load directories
dirs <- c(paste0(data.dir,"LINKAGE/20140630_RIAILs2a/"), 
          paste0(data.dir,"LINKAGE/20140701_RIAILs2b/"), 
          paste0(data.dir,"LINKAGE/20140707_RIAILs2c/"), 
          paste0(data.dir,"LINKAGE/20140708_RIAILs2d/"), 
          paste0(data.dir,"LINKAGE/20140714_RIAILs2e/"), 
          paste0(data.dir,"LINKAGE/20140715_RIAILs2f/"), 
          paste0(data.dir,"LINKAGE/20140721_RIAILs2g/"), 
          paste0(data.dir,"LINKAGE/20140722_RIAILs2h/"), 
          paste0(data.dir,"LINKAGE/20140728_RIAILs2i/"), 
          paste0(data.dir,"LINKAGE/20140729_RIAILs2j/"))

#Read data from all directories
raw <- read_data(dirs)

arsenic_linkage_raw <- list(data.frame(raw[[1]][[1]], day = "score"),
                            data.frame(raw[[1]][[2]], day = "sort"),
                            data.frame(raw[[2]][[1]], day = "score"),
                            data.frame(raw[[2]][[2]], day = "sort"),
                            data.frame(raw[[3]][[1]], day = "score"),
                            data.frame(raw[[3]][[2]], day = "sort"),
                            data.frame(raw[[4]][[1]], day = "score"),
                            data.frame(raw[[4]][[2]], day = "sort"),
                            data.frame(raw[[5]][[1]], day = "score"),
                            data.frame(raw[[5]][[2]], day = "sort"),
                            data.frame(raw[[6]][[1]], day = "score"),
                            data.frame(raw[[6]][[2]], day = "sort"),
                            data.frame(raw[[7]][[1]], day = "score"),
                            data.frame(raw[[7]][[2]], day = "sort"),
                            data.frame(raw[[8]][[1]], day = "score"),
                            data.frame(raw[[8]][[2]], day = "sort"),
                            data.frame(raw[[9]][[1]], day = "score"),
                            data.frame(raw[[9]][[2]], day = "sort"),
                            data.frame(raw[[10]][[1]], day = "score"),
                            data.frame(raw[[10]][[2]], day = "sort"))%>%
  dplyr::bind_rows()

#Remove contamination
raw_noncontam <- remove_contamination(raw)
#Summarize data
summary <- sumplate(raw_noncontam, directories=T, quantiles=T)
#Prune based on biological impossibilities
biopruned <- bioprune(summary)
#Regress out the effect of the assay
assayregressed <- regress(biopruned, assay=TRUE)

#Prune based on bins
RIAILs2BAMFpruned <- bamf_prune(assayregressed, drop=T)
#Regress out control strains
RIAILs2regressed <- regress(RIAILs2BAMFpruned)

# remove rockman riails
arsenic_linkage_processed <- RIAILs2regressed %>%
  dplyr::mutate(n_strain = as.numeric(gsub("QX","",strain)))%>%
  dplyr::filter(n_strain>239)%>%
  dplyr::select(-n_strain) %>%
  dplyr::select(Condition = condition, 
                Strain = strain, 
                Trait = trait, 
                Value = phenotype)

# RUN PCA ON LINKAGE EXPERIMENT
linkage_to_pc <- arsenic_linkage_processed %>%
  tidyr::spread(Trait, Value) %>%
  na.omit() %>% # remove NAs for PC
  dplyr::select(-contains("cv"),-contains("iqr"), -contains("var"))

linkage_arsenic_pca <- data.frame(scale(linkage_to_pc[,3:ncol(linkage_to_pc)]))
linkage_arsenic_pca_analysis <- prcomp(linkage_arsenic_pca)

pc_ve_90 <- which(cumsum(linkage_arsenic_pca_analysis$sdev^2/sum(linkage_arsenic_pca_analysis$sdev^2)) > 0.9)[1]
# VE by PC
#  [1] 0.7087478 0.7843706 0.8314552 0.8672650 0.8942172 0.9123562 0.9267090 0.9363419 0.9456650 0.9539984 0.9601834 0.9654195

linkage_pc_df <- data.frame(Strain = linkage_to_pc$Strain,
                            Condition = linkage_to_pc$Condition,
                            linkage_arsenic_pca_analysis$x[,1:pc_ve_90]) %>%
  tidyr::gather(Trait, Value, -Strain, -Condition)

# combine raw and PC phenotypes
all_linkage_traits <- dplyr::bind_rows(arsenic_linkage_processed, linkage_pc_df)

# extract PC loadings
linkage_rotation_df <- data.frame(trait = row.names(linkage_arsenic_pca_analysis$rotation),
                                  loadings = linkage_arsenic_pca_analysis$rotation[,1:pc_ve_90])

colnames(linkage_rotation_df) <- gsub("loadings.", "", colnames(linkage_rotation_df))

linkage_rotation_df <- linkage_rotation_df %>%
  tidyr::gather(PC, loading, -trait) %>%
  dplyr::rename(Trait = trait)

# Find trait correlations
RIAIL_trait_correlations <- cor(linkage_arsenic_pca)%>%
  data.frame() 

RIAIL_trait_correlations$trait_b <- row.names(RIAIL_trait_correlations)

RIAIL_trait_correlations <- RIAIL_trait_correlations%>%
  tidyr::gather(trait_a, trait_cor, -trait_b) 


write.table(all_linkage_traits, 
            file = glue::glue("{resub.data}Supplemental_Data8_Linkage_Phenotypes.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F)

write.table(RIAIL_trait_correlations, 
            file = glue::glue("{resub.data}Supplemental_Data9_Linkage_trait_correlations.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F)

write.table(linkage_rotation_df, 
            file = glue::glue("{resub.data}Supplemental_Data10_Linkage_PC_Loadings.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F)

# PLOT RIAIL trait correlations and loadings

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

ggsave(glue::glue("{resub.plots}FS4_RIAIL_TRAIT_COR_LOADINGS.pdf"), 
       height = 12, 
       width = 20)

# Perform Linkage Mapping
arsenic_linkage <- all_linkage_traits %>%
  dplyr::rename(strain = Strain, phenotype = Value, trait = Trait) %>%
  dplyr::select(-Condition)

data("N2xCB4856cross")
blankcross <- N2xCB4856cross

arsenic_cross <- linkagemapping::mergepheno(blankcross, arsenic_linkage, set = 2)
map <- linkagemapping::fsearch(arsenic_cross, permutations = 1000, thresh = "GWER")
annotatedmap <- linkagemapping::annotate_lods(map, arsenic_cross)

write.table(annotatedmap, 
            file = glue::glue("{resub.data}Supplemental_Data11_Linkage_Annotated_LODs.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F)

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

ggsave(glue::glue("{resub.plots}FS5_Linkage_PxG.pdf"), 
       height = 8, 
       width = 14)

lods %>%
  dplyr::filter(trait %in% c(".median.TOF",".norm.n", ".PC1")) %>%
  maxlodplot_edit2() +
  arsenic.theme +
  theme(axis.title.x = element_text(size=16),axis.title.y = element_text(size=16),
        strip.background.x = element_blank(),
        panel.grid.major = element_blank())

ggsave(glue::glue("{resub.plots}FS6_Linkage_PC1-brood-length.pdf"), 
       height = 4, 
       width = 12)
########################################################################################################################
# PROCESS LINKAGE MAPPING - END
########################################################################################################################

########################################################################################################################
# LINKAGE MAPPING HERITABILITY ESTIMATES - START
########################################################################################################################
linkage_pheno <- data.table::fread(glue::glue("{resub.data}Supplemental_Data8_Linkage_Phenotypes.tsv"))
lods <- data.table::fread(glue::glue("{resub.data}Supplemental_Data11_Linkage_Annotated_LODs.tsv"))

arsenic_linkage <- linkage_pheno %>%
  dplyr::rename(strain = Strain, phenotype = Value, trait = Trait) %>%
  dplyr::select(-Condition)

data("N2xCB4856cross")

blankcross <- N2xCB4856cross

########################################################################################################################
# relatedness matrix - bloom method
g.chr=lapply(N2xCB4856cross$geno, function(x) x$argmax)
#no point having duplicated markers
g.chr=lapply(g.chr, function(x) x[,-(duplicated(x, MARGIN=2))])
# scale markers
g.chr.s=lapply(g.chr, scale)
g.s=do.call('cbind', g.chr.s)
rownames(g.s)=as.character(N2xCB4856cross$pheno$strain)

# calculated additive relatedness matrix
A=tcrossprod(g.s)/ncol(g.s)
#sometimes I see people divide by the mean of the trace of the crossproduct, makes little difference here
colnames(A)=as.character(N2xCB4856cross$pheno$strain)
rownames(A)=as.character(N2xCB4856cross$pheno$strain)
AA=A*A
########################################################################################################################

# Plot Linkage Mapping - PxG
arsenic_cross <- linkagemapping::mergepheno(blankcross, arsenic_linkage, set = 2)

riail_geno <- linkagemapping::extract_genotype(arsenic_cross)

riail_h2_ls <- list()
for(trt in 1:length(unique(colnames(arsenic_cross$pheno[3:ncol(arsenic_cross$pheno)])))){
  
  anal_trt <- unique(colnames(arsenic_cross$pheno[3:ncol(arsenic_cross$pheno)]))[trt]
  # anal_trt <- gsub("^.","",anal_trt)
  
  trait_peaks <- dplyr::filter(lods, trait == anal_trt) %>%
    na.omit()
  
  if(ncol(trait_peaks) > 0) {
    qtl_ve <- sum(trait_peaks$var_exp, na.rm = T)
    
    if("II" %in% trait_peaks$chr & 7296342 %in% trait_peaks$pos) {
      qtl_II <- TRUE
      qtl_II_ve <-  dplyr::filter(trait_peaks, chr == "II", pos == 7296342) %>%dplyr::pull(var_exp)
    } else {
      qtl_II <- FALSE
      qtl_II_ve <-  NA
    }
  } else {
    qtl_ve <- NA
    qtl_II <- FALSE
    qtl_II_ve <-  NA
  }
  
  riail_p <- arsenic_cross$pheno[240:nrow(arsenic_cross$pheno),trt+2]
  
  strains.not.na <- arsenic_cross$pheno$strain[240:nrow(arsenic_cross$pheno)][which(!is.na(riail_p))]
  
  riail_g <- riail_geno[240:nrow(riail_geno), ]
  
  phenod_geno <- riail_g[which(!is.na(riail_p)), ]
  riail_p <- riail_p[which(!is.na(riail_p))]
  
  y <- riail_p 
  ######################################################################################################################
  # Bloom method
  Z=matrix(0,length(y), ncol(A))
  Z[cbind(1:length(y), match(strains.not.na, colnames(A)))]=1
  ZA=Z%*%A%*%t(Z)
  ZAA=Z%*%AA%*%t(Z)
  
  rr_bloom=regress::regress(y~1, ~ZA+ZAA, verbose=T)
  rs_bloom=rr_bloom$sigma
  
  H2_bloom <- sum(rs_bloom[1:2])/sum(rs_bloom)
  h2_bloom <- sum(rs_bloom[1])/sum(rs_bloom)
  ######################################################################################################################
  
  ######################################################################################################################
  # Sommer method
  Za <- diag(length(y))
  Ze <- diag(length(y))
  
  A_sommer <- sommer::A.mat(phenod_geno)
  E_sommer <- sommer::E.mat(phenod_geno)
  
  # with epistatic term
  complete_add <- sommer::mmer(Y=y, Z=list(add = list(Z=Za,K=A_sommer), 
                                           epp = list(Z=Ze,K=E_sommer)))
  
  suma <- summary(complete_add)$var.comp.table
  
  H2_sommer <- sum(suma[1:2,1])/sum(suma[,1])
  h2_sommer <- sum(suma[1,1])/sum(suma[,1])
  ######################################################################################################################
  
  riail_h2_ls[[trt]] <- data.frame(trait = anal_trt,
                                   H2_realized = H2_sommer,
                                   h2_realized = h2_sommer,
                                   H2_expectation = H2_bloom,
                                   h2_expectation = h2_bloom,
                                   QTL_VE = qtl_ve,
                                   QTL_II = qtl_II,
                                   QTL_II_VE = qtl_II_ve)
}

riail_h2_df <- dplyr::bind_rows(riail_h2_ls)

riail_h2_df$trait <- gsub("^.","",riail_h2_df$trait)

write.table(riail_h2_df, 
            file = glue::glue("{resub.data}Supplemental_Data12_Linkage_H2_ESTIMATES.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F)

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
       glue::glue("{resub.plots}FS7_Linkage_heritability_estimates.pdf"), 
       height = 12, 
       width = 16)

########################################################################################################################
# LINKAGE MAPPING HERITABILITY ESTIMATES - END
########################################################################################################################

# LINKAGE QTL SUMMARY
lods <- data.table::fread(glue::glue("{resub.data}Supplemental_Data11_Linkage_Annotated_LODs.tsv"))

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
       glue::glue("{resub.plots}FS8_Linkage_QTL_Summary.pdf"), 
       height = 14, 
       width = 16)


########################################################################################################################
# PROCESS NIL PHENOTYPES - START
########################################################################################################################
# 1000µM arscenic trioxide

raw <- easysorter::read_data(paste0(data.dir,"NILs_n_SWAPs/20170403_dbtswap/"))

raw_nocontam <- remove_contamination(raw)

summedraw <- sumplate(raw_nocontam, quantiles = TRUE)%>%
  dplyr::mutate(f.L1L2L3 = f.L1+f.L2L3)

summedraw1 <- dplyr::mutate(summedraw, n_sorted = n/norm.n)%>%
  tidyr::gather(trait, value, -date, -experiment, -round, -assay, -plate, 
                -condition, -control, -strain, -row, -col, -n_sorted )%>%
  dplyr::select(condition, strain,plate, row, col, control,trait, n_sorted)

biopruned <- bioprune(summedraw)

controlregressed <- easysorter::regress(biopruned, assay = FALSE)
controlregressed <- dplyr::left_join(controlregressed,summedraw1,
                                     by = c("condition", "strain","plate", "row", "col", 
                                            "control","trait","date","experiment","round","assay"))


nil_to_pc <- controlregressed %>%
  dplyr::ungroup() %>%
  dplyr::select(strain, condition, plate, row, col, trait, phenotype) %>%
  dplyr::mutate(u_strain = paste(strain,condition,plate,row,col,sep="_"))%>%
  tidyr::spread(trait, phenotype) %>%
  na.omit() %>%
  dplyr::select(-contains("green"),-contains("red"),-contains("yellow"), -contains("var"), -contains("iqr"), -contains("cv"))

nil_arsenic_pca <- data.frame(scale(nil_to_pc[,7:ncol(nil_to_pc)]))
nil_arsenic_pca_analysis <- prcomp(nil_arsenic_pca)

pc_ve_90 <- which(cumsum(nil_arsenic_pca_analysis$sdev^2/sum(nil_arsenic_pca_analysis$sdev^2)) > 0.9)[1]

nil_pc_df <- data.frame(strain = nil_to_pc$strain,
                        condition = nil_to_pc$condition,
                        u_id = nil_to_pc$u_strain,
                        nil_arsenic_pca_analysis$x[,1:pc_ve_90]) %>%
  tidyr::gather(trait, value, -strain, -condition, -u_id)


nil_all_traits <- controlregressed %>%
  dplyr::ungroup()%>%
  dplyr::mutate(u_id = paste(strain,condition,plate,row,col,sep="_"))%>%
  dplyr::select(strain, condition, u_id, trait, value = phenotype) %>%
  dplyr::bind_rows(., nil_pc_df) %>%
  dplyr::rename(Strain = strain, Condition = condition, Trait = trait, Value = value)

# extract PC loadings
NIL_rotation_df <- data.frame(trait = row.names(nil_arsenic_pca_analysis$rotation),
                              nil_arsenic_pca_analysis$rotation[,1:pc_ve_90])

NIL_rotation_df <- NIL_rotation_df %>%
  tidyr::gather(PC, loading, -trait) %>%
  dplyr::rename(Trait = trait)

# Find trait correlations
NIL_trait_correlations <- cor(nil_to_pc[,7:ncol(nil_to_pc)])%>%
  data.frame() 

NIL_trait_correlations$trait_b <- row.names(NIL_trait_correlations)

NIL_trait_correlations <- NIL_trait_correlations%>%
  tidyr::gather(trait_a, trait_cor, -trait_b) 


write.table(nil_all_traits, 
            file = glue::glue("{resub.data}Supplemental_Data13_NIL_control_regressed.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F )

write.table(NIL_trait_correlations, 
            file = glue::glue("{resub.data}Supplemental_Data14_NIL_Trait_Correlations.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F )

write.table(NIL_rotation_df, 
            file = glue::glue("{resub.data}Supplemental_Data15_NIL_PC_Loadings.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F )

# nil genotypes are supplementa data 16

########################################################################################################################
# PROCESS NIL PHENOTYPES - END
########################################################################################################################




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

# Plot Linkage Mapping - PxG
# arsenic_cross <- linkagemapping::mergepheno(blankcross, arsenic_linkage, set = 2)
# 
# pxg_linkage_PC_II <- pxgplot_edit_1peak(arsenic_cross, dplyr::filter(lods, trait == ".PC1", chr == "II"))+ 
#   scale_x_discrete(breaks=c("N2", "CB4856"),labels=c("N2", "CB4856"))+
#   labs(y = paste0("PC 1"), x = "RIAIL Genotype at II:7296342") +
#   arsenic.theme +
#   theme(strip.background = element_blank(),
#         strip.text = element_blank(),
#         panel.grid.minor = element_blank())

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
  # geom_text(aes(50, -15, label = paste0("y = ", unique(slope), 
  #                                       "x + ", unique(intercept), 
  #                                       "\n R^2 = ", unique(adj.r.squared), "; p = ", unique(pvalue))),
  #           data = dplyr::filter(riail_lm, CorTrait == "norm.n")) +
  # geom_text(aes(-55, 15, label = paste0("y = ", unique(slope), 
  #                                       "x + ", unique(intercept), 
  #                                       "\n R^2 = ", unique(adj.r.squared), "; p = ", unique(pvalue))),
  #           data = dplyr::filter(riail_lm, CorTrait == "median.TOF")) +
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
       glue::glue("{resub.plots}F1_Linkage_NILs.pdf"), 
       height = 8, 
       width = 16)

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
# 4.007
# N2 - left NIL
# 3.609 - 90% recapitulation
# N2 - right NIL
# 3.792 - 94.6% recapitulation

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

ggsave(glue::glue("{resub.plots}FS9_NIL_length_brood_recapitulation.pdf"), 
       height = 6, 
       width = 10)


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

ggsave(glue::glue("{resub.plots}FS10_NIL_SWAP_TRAIT_COR_LOADINGS.pdf"), 
       height = 12, 
       width = 20)


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
  geom_text(aes(45, -10, 
                label = paste("Adj R2 = ", adj.r.squared, "\n",
                              "Intercept =",intercept, "\n",
                              "Slope =", slope, "\n",
                              "P =", pvalue)), data = dplyr::filter(nil_lm, CorTrait == "norm.n")) +
  geom_text(aes(-35, 10, 
                label = paste("Adj R2 = ", adj.r.squared, "\n",
                              "Intercept =",intercept, "\n",
                              "Slope =", slope, "\n",
                              "P =", pvalue)), data = dplyr::filter(nil_lm, CorTrait == "median.TOF"))+
  theme(panel.grid = element_blank())


ggsave(plot = nil_plot,
       glue::glue("{resub.plots}FS11_NIL_PC_length_brood_correlation.pdf"), 
       height = 6, 
       width = 10)

########################################################################################################################
# LINKAGE MAPPING - PLOTS - END
########################################################################################################################


########################################################################################################################
# PROCESS GWA MAPPING PHENOTYPES - START
########################################################################################################################

g8 <- easysorter::read_data(c(paste0(data.dir,"GWAS/20140617_arsenicGWAS8a/"),
                              paste0(data.dir,"GWAS/20140617_arsenicGWAS8b/")))

raw_nocontam <- easysorter::remove_contamination(g8)
summedraw_a <- easysorter::sumplate(raw_nocontam, quantiles = TRUE, directories = T)
summedraw_a <- easysorter::bioprune(summedraw_a)
summedraw_a <- easysorter::bamf_prune(summedraw_a)

g8_prune <- summedraw_a%>%
  dplyr::filter(bamfoutlier1 != TRUE,bamfoutlier2 != TRUE,bamfoutlier2 != TRUE)%>%
  dplyr::select(-bamfoutlier1,-bamfoutlier2,-bamfoutlier3)

# # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # assay regression 

g8_prune_assay <- easysorter::regress(g8_prune, assay = TRUE)

# # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # control regression

g8_regressed <- easysorter::regress(g8_prune_assay, assay = F)

# # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # ## ## ## ## # PCA MAPPINGS

arsenic_pca <- g8_regressed%>%
  dplyr::ungroup()%>%
  dplyr::select(Strain = strain, trait, phenotype)%>%
  tidyr::spread(trait, phenotype) %>%
  dplyr::select(-contains("cv"),-contains("iqr"), -contains("var"))
  

row.names(arsenic_pca) <- arsenic_pca$Strain
strain_order <- arsenic_pca$Strain

# GWAS PC analysis
arsenic_pca_ls <- run_pca(arsenic_pca, format = "w", drop_method = "t", n_pc = NA)

arsenic_pca_df <- arsenic_pca_ls[[1]]
arsenic_pca_loadings <- arsenic_pca_ls[[2]]
arsenic_pca_ve <- arsenic_pca_ls[[3]]

gwas_all_traits <- g8_regressed %>%
  dplyr::ungroup() %>%
  dplyr::select(Strain = strain, Trait = trait, Value = phenotype) %>%
  dplyr::bind_rows(., arsenic_pca_df) %>%
  tidyr::spread(Trait, Value) %>%
  dplyr::rename(strain = Strain) %>%
  dplyr::select(-contains("var"),-contains("cv"),-contains("iqr"))

gwas_all_traits_fixnames <- cegwas2::process_phenotypes(gwas_all_traits)

write.table(gwas_all_traits_fixnames, 
            file = glue::glue("{resub.data}Supplemental_Data17_GWAS_ALL_TRAITS.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F )

write.table(arsenic_pca_loadings, 
            file = glue::glue("{resub.data}Supplemental_Data18_GWAS_PC_Loadings.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F )

# write.table(arsenic_pca_ve, 
#             file = glue::glue("{resub.data}Supplemental_Data19_GWAS_PC_VE.tsv"), 
#             sep = "\t", 
#             quote = F, 
#             col.names = T, 
#             row.names = F )

# Find trait correlations
gwas_trait_correlations <- cor(arsenic_pca[,2:ncol(arsenic_pca)],use = "complete.obs")%>%
  data.frame() 

gwas_trait_correlations$trait_b <- row.names(gwas_trait_correlations)

gwas_trait_correlations <- gwas_trait_correlations%>%
  tidyr::gather(trait_a, trait_cor, -trait_b)

write.table(gwas_trait_correlations, 
            file = glue::glue("{resub.data}Supplemental_Data19_GWAS_Trait_Correlations.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F )

dir.create(glue::glue("{resub.data}Individual_Arsenic_Traits"))

for(trait in colnames(gwas_all_traits_fixnames[2:ncol(gwas_all_traits_fixnames)])){
  
  col_name <- quo_name(enquo(trait))
  
  temp_pheno <- gwas_all_traits %>%
    dplyr::select(strain, trait) 
  
  write.table(temp_pheno, 
              file = glue::glue("{resub.data}Individual_Arsenic_Traits/{trait}.tsv"), 
              sep = "\t", 
              quote = F, 
              col.names = T, 
              row.names = F )
  
}

########################################################################################################################
# PROCESS GWA MAPPING PHENOTYPES - END
########################################################################################################################

########################################################################################################################
# mappings are run using cegwas2-nf, which can be found on github: AndersenLab/cegwas2-nf
########################################################################################################################

########################################################################################################################
# PROCESS GWA MAPPING PHENOTYPES - GENOMIC HERITABILITY - START
########################################################################################################################
geno_matrix <- readr::read_tsv(glue::glue("{resub.data}Supplemental_Data21_Genotype_Matrix.tsv"))%>%
  na.omit()

new_pheno <- readr::read_tsv(file = glue::glue("{resub.data}Supplemental_Data17_GWAS_ALL_TRAITS.tsv"))%>%
  as.data.frame() %>%
  tidyr::gather(trait,value, -strain)

narrow_h2 <- list()
for(i in 1:length(unique(new_pheno$trait))){
  
  h2_trait <- unique(new_pheno$trait)[i]
  
  strains <- new_pheno %>%
    dplyr::filter(trait == h2_trait ) %>%
    dplyr::arrange(strain)%>%
    dplyr::pull(strain)
  
  y <- new_pheno %>%
    dplyr::filter(trait == h2_trait ) %>%
    dplyr::arrange(strain) %>%
    dplyr::pull(value)
  
  Za <- diag(length(y))
  Ze <- diag(length(y))
  
  A <- sommer::A.mat(t(geno_matrix[,strains]))
  E <- sommer::E.mat(t(geno_matrix[,strains]))
  
  # additive only
  complete_broad <- sommer::mmer(Y=y, Z=list(add = list(Z=Za,K=A)),
                                 IMP = T, init.equal = F)
  
  suma <- summary(complete_broad)$var.comp.table
  
  H2 <- sum(suma[1,1])/sum(suma[,1])
  
  # with epistatic term
  complete_add <- sommer::mmer(Y=y, Z=list(add = list(Z=Za,K=A), 
                                           epp = list(Z=Ze,K=E)),
                               IMP = T, init.equal = F)
  
  suma <- summary(complete_add)$var.comp.table
  
  H2_2 <- sum(suma[1:2,1])/sum(suma[,1])
  h2_2 <- sum(suma[1,1])/sum(suma[,1])
  
  narrow_h2[[i]] <- data.frame(Broad_H2_add_model = H2,
                               narrow = h2_2,
                               Broad_H2_add_epp_model = H2_2,
                               trait = as.character(h2_trait))
}

gwas_h2 <- data.frame(dplyr::bind_rows(narrow_h2))

write.table(gwas_h2, 
            file = glue::glue("{resub.data}Supplemental_Data22_GWAS_Genomic_H2_Estimates.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F )

########################################################################################################################
# PROCESS GWA MAPPING PHENOTYPES - GENOMIC HERITABILITY - END
########################################################################################################################

########################################################################################################################
# GWA MAPPING - PLOT TRAIT CORRELATIONS - START
########################################################################################################################

gwas_traits <- data.table::fread(glue::glue("{resub.data}Supplemental_Data17_GWAS_ALL_TRAITS.tsv"))
gwas_loadings <- data.table::fread(glue::glue("{resub.data}Supplemental_Data18_GWAS_PC_Loadings.tsv"))
gwas_trait_cor <- data.table::fread(glue::glue("{resub.data}Supplemental_Data19_GWAS_Trait_Correlations.tsv"))

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

ggsave(glue::glue("{resub.plots}FS12_GWAS_Trait-correlation_loadings.pdf"), 
       height = 12, 
       width = 20)

########################################################################################################################
# GWA MAPPING - PLOT TRAIT CORRELATIONS - END
########################################################################################################################


########################################################################################################################
# PROCESS GWA MAPPING PHENOTYPES - GENOMIC HERITABILITY - DR H2 - EFFECT SIZE - PLOTS - START
########################################################################################################################

gwas_h2 <- data.table::fread(glue::glue("{resub.data}Supplemental_Data22_GWAS_Genomic_H2_Estimates.tsv"))
h2_df <- data.table::fread(glue::glue("{resub.data}Supplemental_Data2_DR-Heritability_estimates.tsv"))

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

ggsave(glue::glue("{resub.plots}FS13_GWAS-H2_v_h2.pdf"), 
       height = 8, 
       width = 10)

########################################################################################################################
# PROCESS GWA MAPPING PHENOTYPES - GENOMIC HERITABILITY - DR H2 - EFFECT SIZE - PLOTS - END
########################################################################################################################



########################################################################################################################
# GWA MAPPING - PLOT FIGURE 2 - START
########################################################################################################################

gwas_mapping <- data.table::fread(glue::glue("{resub.data}Supplemental_Data20_PC1_processed_mapping.tsv"))
independent_tests <- data.table::fread(glue::glue("{resub.data}Supplemental_Data24_total_independent_tests.txt"))
gwas_fine_mapping <- data.table::fread(glue::glue("{resub.data}Supplemental_Data23_PC1_snpeff_genes.tsv")) 
geno_matrix <- readr::read_tsv(glue::glue("{resub.data}Supplemental_Data21_Genotype_Matrix.tsv"))%>%
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

(H2_no_qtl-H2_qtl)/H2_no_qtl
# 0.211

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
       glue::glue("{resub.plots}F2_GWAS.pdf"), 
       height = 8, 
       width = 12)

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

ggsave(glue::glue("{resub.plots}FS15_GWAS_FINE_MAP.pdf"), 
       height = 8, 
       width = 12)

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
# 3.02/3.999
# 0.7551888
# # cb - 590
# 2.557/3.999
# 0.6394099

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
        axis.title.y = ggplot2::element_text(size = 16, face = "bold", color = "black", vjust = -0.3),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())+
  labs( y = paste0("PC 1"))

ggsave(glue::glue("{resub.plots}F3_Allele_Swap_PC1.pdf"), 
       height = 4, 
       width = 6)

sig_swap <- swap_pheno%>%
  dplyr::filter(Trait == "PC1",
                Strain != "ECA591", 
                Strain != "ECA414", 
                Strain != "ECA434", 
                Strain != "ECA582", 
                Strain != "ECA589") %>%
  dplyr::select(strain = Strain, trait = Trait, phenotype = Value)

run_TukeyHSD(sig_swap,trt = "PC1")
# $`stat_df$strain`
# diff        lwr        upr    p adj
# ECA581-CB4856  -1.701308  -2.604441 -0.7981744 1.42e-05
# ECA590-CB4856  -8.444066  -9.341742 -7.5463896 0.00e+00
# N2-CB4856     -10.608607 -11.506283 -9.7109313 0.00e+00
# ECA590-ECA581  -6.742758  -7.645891 -5.8396250 0.00e+00
# N2-ECA581      -8.907300  -9.810433 -8.0041668 0.00e+00
# N2-ECA590      -2.164542  -3.062218 -1.2668657 0.00e+00

# brood and size
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


ggsave(glue::glue("{resub.plots}FS16_Allele_Swap_Brood_Length.pdf"), 
       height = 6, 
       width = 10)

sig_swap <- swap_pheno%>%
  dplyr::filter(Trait %in% c("median.TOF", "norm.n"),
                Strain != "ECA591", 
                Strain != "ECA414", 
                Strain != "ECA434", 
                Strain != "ECA582", 
                Strain != "ECA589") %>%
  dplyr::select(strain = Strain, trait = Trait, phenotype = Value)

run_TukeyHSD(sig_swap,trt = "median.TOF")
# $`stat_df$strain`
# diff        lwr        upr     p adj
# ECA581-CB4856   2.95946  -3.356933   9.275854 0.6174821
# ECA590-CB4856 -33.69528 -39.973503 -27.417048 0.0000000
# N2-CB4856     -54.35146 -60.629687 -48.073232 0.0000000
# ECA590-ECA581 -36.65474 -42.971129 -30.338342 0.0000000
# N2-ECA581     -57.31092 -63.627314 -50.994527 0.0000000
# N2-ECA590     -20.65618 -26.934412 -14.377957 0.0000000
run_TukeyHSD(sig_swap,trt = "norm.n")
# $`stat_df$strain`
# diff       lwr        upr    p adj
# ECA581-CB4856 -13.42928 -20.58601  -6.272544 1.55e-05
# ECA590-CB4856 -48.44585 -55.55934 -41.332359 0.00e+00
# N2-CB4856     -27.37680 -34.49029 -20.263315 0.00e+00
# ECA590-ECA581 -35.01657 -42.17330 -27.859838 0.00e+00
# N2-ECA581     -13.94753 -21.10426  -6.790795 6.70e-06
# N2-ECA590      21.06904  13.95555  28.182533 0.00e+00


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
  geom_text(aes(45, -5, 
                label = paste("Adj R2 = ", adj.r.squared, "\n",
                              "Intercept =",intercept, "\n",
                              "Slope =", slope, "\n",
                              "P =", pvalue)), data = dplyr::filter(nil_lm, CorTrait == "norm.n")) +
  geom_text(aes(-35, 7, 
                label = paste("Adj R2 = ", adj.r.squared, "\n",
                              "Intercept =",intercept, "\n",
                              "Slope =", slope, "\n",
                              "P =", pvalue)), data = dplyr::filter(nil_lm, CorTrait == "median.TOF"))+
  theme(panel.grid.minor = element_blank())


ggsave(glue::glue("{resub.plots}FS17_Allele_Swap_PC1_Brood_Size_Correlation.pdf"), 
       height = 8, 
       width = 12)

########################################################################################################################
# DBT-1 C78S SWAP - PLOTS and EFFECT SIZES- END
########################################################################################################################

########################################################################################################################
# Fatty acid measurements in L1 animals with arsenic
########################################################################################################################

metabolites <- data.table::fread(glue::glue("{resub.data}Supplemental_Data25_Metabolite_measurements.csv"))%>%
  tidyr::separate(Sample, into = c("strain", "replicate", "concentration"), sep = "-")%>%
  dplyr::select(-TIC)%>%
  dplyr::filter(!(strain == "581" & replicate == "C" & concentration == "Mock"))%>%
  dplyr::mutate(rC15 = C15iso/C15n,
                rC17 = C17iso/C17n)%>%
  tidyr::gather(compound, value, -strain, -replicate, -concentration) %>%
  dplyr::filter(concentration != "200")

write.table(metabolites, 
            file = glue::glue("{resub.data}Supplemental_data26_Processed_Metabolite_Measurements.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F)

metabolites <- data.table::fread(glue::glue("{resub.data}Supplemental_data26_Processed_Metabolite_Measurements.tsv"))

# second experiment to get more replicates of ECA581
L1_fa <- readr::read_csv(glue::glue("{resub.data}Supplemental_Data27_L1_FA_N2_ECA581.csv")) 

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
  geom_errorbar(aes(ymin=mean_value-abs(sd_value), ymax=mean_value+abs(sd_value)),
                width=.2)+
  facet_wrap(expt~FA, scales = "free")+
  scale_fill_manual(values=c("orange","gray50","blue","gray50"))+
  theme_bw(16)+
  arsenic.theme+
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        strip.text.y = element_text(face ="bold"))+
  labs(y = "Arsenic - Control")


ggsave(glue::glue("{resub.plots}F4B_ISO_FA_measurements.pdf"), 
       height = 8, 
       width = 8)

# cb and cb swap experiment - C15ISO
stat_df <- complete_arsenic %>%
  dplyr::ungroup() %>%
  dplyr::filter(grepl("15",FA),
                Strain %in% c("CB", "590"))%>%
  dplyr::select(strain = Strain, phenotype= delta_control)

TukeyHSD(aov(stat_df$phenotype ~ stat_df$strain))
#0.0427733
stat_df%>%
  group_by(strain)%>%
  dplyr::summarise(mfa=mean(phenotype))
# fold change
(0.172-0.0457)/0.0457
# 2.763676

# cb and cb swap experiment - C17ISO
stat_df <- complete_arsenic %>%
  dplyr::ungroup() %>%
  dplyr::filter(grepl("17",FA),
                Strain %in% c("CB", "590"))%>%
  dplyr::select(strain = Strain, phenotype= delta_control)

TukeyHSD(aov(stat_df$phenotype ~ stat_df$strain))
#0.164721
stat_df%>%
  group_by(strain)%>%
  dplyr::summarise(mfa=mean(phenotype))
# fold change
(1.92-0.756)/0.756
# 1.539683

# n2 and n2 swap experiment - C17ISO
stat_df <- complete_arsenic %>%
  dplyr::ungroup() %>%
  dplyr::filter(grepl("15",FA),
                Strain %in% c("N2", "ECA581"))%>%
  dplyr::select(strain = Strain, phenotype= delta_control)

TukeyHSD(aov(stat_df$phenotype ~ stat_df$strain))
#0.0357889
stat_df%>%
  group_by(strain)%>%
  dplyr::summarise(mfa=mean(phenotype))
# fold change
(0.0410-0.00477)/0.00477
# 7.595388

# n2 and n2 swap experiment - C17ISO
stat_df <- complete_arsenic %>%
  dplyr::ungroup() %>%
  dplyr::filter(grepl("17",FA),
                Strain %in% c("N2", "ECA581"))%>%
  dplyr::select(strain = Strain, phenotype= delta_control)

TukeyHSD(aov(stat_df$phenotype ~ stat_df$strain))
#0.003747
stat_df%>%
  group_by(strain)%>%
  dplyr::summarise(mfa=mean(phenotype))
# fold change
abs(-0.221-0.0392)/0.0392
# 6.637755

### PLOT C17 for CB and SWAP

c17 <- metabolites%>%
  dplyr::filter(compound %in% c("C17iso"),
                strain %in% c("CB","590"),
                concentration != "200") 

c17%>%
  ggplot()+
  aes(x = factor(strain, 
                 levels = c("N2","CB","581","590"), 
                 labels = c("N2", "CB4856","N2\nDBT-1(S78)", "CB4856\nDBT-1(C78)")), 
      y = (value), 
      fill = concentration)+
  geom_point(shape =21, size = 3)+
  scale_fill_manual(values = c("hotpink3","cadetblue3"))+
  theme_classic()+
  arsenic.theme+ theme(legend.position = "None",
                       axis.title.x = ggplot2::element_text(size = 0),
                       axis.title.y = ggplot2::element_text(size = 16))+
  labs(y = "C17iso")

ggsave(glue::glue("{resub.plots}FS18_C17ISO_CB_SBswap.pdf"), 
       height = 4, 
       width = 6)

for( s in unique(c17$strain)){
  rC15 <- metabolites%>%
    dplyr::filter(compound %in% c("C17iso"), concentration != "200") %>%
    dplyr::filter(strain==s)
  
  print(s)
  fit <- aov(value ~ concentration , data=rC15)
  print(TukeyHSD(fit))
}

# cb iso
# $concentration
# diff      lwr      upr     p adj
# Mock-100 13357867 -4190114 30905849 0.1020881
#swap iso
# $concentration
# diff     lwr      upr     p adj
# Mock-100 54496753 9154091 99839415 0.0289175
##### PLOT STRAIGHT CHAIN

metabolites <- data.table::fread(glue::glue("{resub.data}Supplemental_data26_Processed_Metabolite_Measurements.tsv"))

# second experiment to get more replicates of ECA581
L1_fa <- readr::read_csv(glue::glue("{resub.data}Supplemental_Data27_L1_FA_N2_ECA581.csv")) 

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
  scale_fill_manual(values=c("orange","gray50","blue","gray50"))+
  theme_bw(16)+
  arsenic.theme+
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        strip.text.y = element_text(face ="bold"))+
  labs(y = "Arsenic - Control")

ggsave(glue::glue("{resub.plots}FS19_StraightChain_FA_measurements.pdf"), 
       height = 8, 
       width = 8)


# cb and cb swap experiment - C15ISO
stat_df <- complete_arsenic %>%
  dplyr::ungroup() %>%
  dplyr::filter(grepl("15",FA),
                Strain %in% c("CB", "590"))%>%
  dplyr::select(strain = Strain, phenotype= delta_control)

TukeyHSD(aov(stat_df$phenotype ~ stat_df$strain))
#0.7670935

# cb and cb swap experiment - C17ISO
stat_df <- complete_arsenic %>%
  dplyr::ungroup() %>%
  dplyr::filter(grepl("17",FA),
                Strain %in% c("CB", "590"))%>%
  dplyr::select(strain = Strain, phenotype= delta_control)

TukeyHSD(aov(stat_df$phenotype ~ stat_df$strain))
#0.1464522

# n2 and n2 swap experiment - C17ISO
stat_df <- complete_arsenic %>%
  dplyr::ungroup() %>%
  dplyr::filter(grepl("15",FA),
                Strain %in% c("N2", "ECA581"))%>%
  dplyr::select(strain = Strain, phenotype= delta_control)

TukeyHSD(aov(stat_df$phenotype ~ stat_df$strain))
#0.1042516

# n2 and n2 swap experiment - C17ISO
stat_df <- complete_arsenic %>%
  dplyr::ungroup() %>%
  dplyr::filter(grepl("17",FA),
                Strain %in% c("N2", "ECA581"))%>%
  dplyr::select(strain = Strain, phenotype= delta_control)

TukeyHSD(aov(stat_df$phenotype ~ stat_df$strain))
#0.1805046

#### PLOT RATIOS in CONTROL CONDITION MOTHER FUCKER

metabolites <- data.table::fread(glue::glue("{resub.data}Supplemental_data26_Processed_Metabolite_Measurements.tsv"))

# second experiment to get more replicates of ECA581
L1_fa <- readr::read_csv(glue::glue("{resub.data}Supplemental_Data27_L1_FA_N2_ECA581.csv")) 

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
  theme_bw(16)+
  arsenic.theme+
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        strip.text.y = element_text(face ="bold"))+
  labs(y = "Metabolite Levels")

ggsave(glue::glue("{resub.plots}FS20_iso_ratio_control.pdf"), 
       height = 8, 
       width = 8)

stat_df <- complete_arsenic %>%
  dplyr::ungroup() %>%
  dplyr::filter(grepl("15",FA),
                Strain %in% c("N2", "ECA581"))%>%
  dplyr::select(strain = Strain, phenotype= control_value)

TukeyHSD(aov(stat_df$phenotype ~ stat_df$strain))
#0.1239674

stat_df <- complete_arsenic %>%
  dplyr::ungroup() %>%
  dplyr::filter(grepl("17",FA),
                Strain %in% c("N2", "ECA581"))%>%
  dplyr::select(strain = Strain, phenotype= control_value)

TukeyHSD(aov(stat_df$phenotype ~ stat_df$strain))
#0.0044667

stat_df <- complete_arsenic %>%
  dplyr::ungroup() %>%
  dplyr::filter(grepl("15",FA),
                Strain %in% c("CB", "590"))%>%
  dplyr::select(strain = Strain, phenotype= control_value)

TukeyHSD(aov(stat_df$phenotype ~ stat_df$strain))
#0.0168749

stat_df <- complete_arsenic %>%
  dplyr::ungroup() %>%
  dplyr::filter(grepl("17",FA),
                Strain %in% c("CB", "590"))%>%
  dplyr::select(strain = Strain, phenotype= control_value)

TukeyHSD(aov(stat_df$phenotype ~ stat_df$strain))
#0.0342525

# # # # # # # ## # ISO AND STRAIGHT IN CONTROL CONDITIONS

metabolites <- data.table::fread(glue::glue("{resub.data}Supplemental_data26_Processed_Metabolite_Measurements.tsv"))

# second experiment to get more replicates of ECA581
L1_fa <- readr::read_csv(glue::glue("{resub.data}Supplemental_Data27_L1_FA_N2_ECA581.csv")) 

controls_new <- L1_fa %>% 
  tidyr::gather(FA,value,-Strain,-Condition,-Replicate) %>%
  dplyr::filter(Condition == "Water")%>%
  dplyr::select(-Condition) %>%
  dplyr::rename(control_value = value) %>%
  dplyr::filter(FA %in% c("C15_branched", "C17iso", "C17n", "C15n"))%>%
  dplyr::select(Strain, FA, control_value)

controls_new$FA <- gsub("C15_branched", "C15iso", controls_new$FA)
# controls_new$FA <- gsub("15_ratio", "C15iso/C15n", controls_new$FA)

rC15 <- metabolites%>%
  dplyr::filter(compound %in% c("C15iso", "C17iso", "C17n", "C15n"), concentration == "Mock") %>%
  dplyr::filter(strain %in% c("590", "CB"))

# rC15$compound <- gsub("rC15", "C15iso/C15n", rC15$compound)
# rC15$compound <- gsub("rC17", "C17iso/C17n", rC15$compound)
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
  facet_wrap(expt~FA, scales = "free",ncol = 2)+
  scale_fill_manual(values=c("orange","gray50","blue","gray50"))+
  theme_bw(16)+
  arsenic.theme+
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        strip.text.y = element_text(face ="bold"))+
  labs(y = "Metabolite Level")


ggsave(glue::glue("{resub.plots}FS21_iso_and_straight_control.pdf"), 
       height = 14, 
       width = 8)

stat_df <- complete_arsenic %>%
  dplyr::ungroup() %>%
  dplyr::filter(grepl("15n",FA),
                Strain %in% c("N2", "ECA581"))%>%
  dplyr::select(strain = Strain, phenotype= control_value)

TukeyHSD(aov(stat_df$phenotype ~ stat_df$strain))
#0.5817993

stat_df <- complete_arsenic %>%
  dplyr::ungroup() %>%
  dplyr::filter(grepl("17n",FA),
                Strain %in% c("N2", "ECA581"))%>%
  dplyr::select(strain = Strain, phenotype= control_value)

TukeyHSD(aov(stat_df$phenotype ~ stat_df$strain))
#0.35827

stat_df <- complete_arsenic %>%
  dplyr::ungroup() %>%
  dplyr::filter(grepl("15n",FA),
                Strain %in% c("CB", "590"))%>%
  dplyr::select(strain = Strain, phenotype= control_value)

TukeyHSD(aov(stat_df$phenotype ~ stat_df$strain))
#0.0787388

stat_df <- complete_arsenic %>%
  dplyr::ungroup() %>%
  dplyr::filter(grepl("17n",FA),
                Strain %in% c("CB", "590"))%>%
  dplyr::select(strain = Strain, phenotype= control_value)

TukeyHSD(aov(stat_df$phenotype ~ stat_df$strain))
#0.0086572

stat_df <- complete_arsenic %>%
  dplyr::ungroup() %>%
  dplyr::filter(grepl("15iso",FA),
                Strain %in% c("N2", "ECA581"))%>%
  dplyr::select(strain = Strain, phenotype= control_value)

TukeyHSD(aov(stat_df$phenotype ~ stat_df$strain))
#0.0265059

stat_df <- complete_arsenic %>%
  dplyr::ungroup() %>%
  dplyr::filter(grepl("17iso",FA),
                Strain %in% c("N2", "ECA581"))%>%
  dplyr::select(strain = Strain, phenotype= control_value)

TukeyHSD(aov(stat_df$phenotype ~ stat_df$strain))
#0.0022501

stat_df <- complete_arsenic %>%
  dplyr::ungroup() %>%
  dplyr::filter(grepl("15iso",FA),
                Strain %in% c("CB", "590"))%>%
  dplyr::select(strain = Strain, phenotype= control_value)

TukeyHSD(aov(stat_df$phenotype ~ stat_df$strain))
#0.0036201

stat_df <- complete_arsenic %>%
  dplyr::ungroup() %>%
  dplyr::filter(grepl("17iso",FA),
                Strain %in% c("CB", "590"))%>%
  dplyr::select(strain = Strain, phenotype= control_value)

TukeyHSD(aov(stat_df$phenotype ~ stat_df$strain))
#0.0086572

########################################################################################################################
# Fatty acid measurements in L4 animals
########################################################################################################################

L4_fa <- readr::read_csv(glue::glue("{resub.data}Supplemental_Data28_YA_Fatty_Acid_Measurements.csv")) %>%
  dplyr::mutate(`C17iso/C18straight` = C17iso/C18,
                `C15iso/C18straight` = C15iso/C18) %>%
  dplyr::select(Strain, `C17iso/C18straight`, `C15iso/C18straight`)%>%
  tidyr::gather(ratio, value, -Strain)

ggplot(L4_fa)+
  aes(x = factor(Strain, 
                 levels = c("N2","ECA581","CB4856","ECA590"),
                 labels =c("N2\nDBT-1(C78)", 
                           "N2\nDBT-1 (S78)",
                           "CB4856\nDBT-1 (S78)",
                           "CB4856\nDBT-1(C78)" ),
                 ordered = T), 
      y = value,
      fill= factor(Strain, 
                   levels = c("N2","ECA581","CB4856","ECA590"),
                   labels =c("N2\nDBT-1(C78)", 
                             "N2\nDBT-1 (S78)",
                             "CB4856\nDBT-1 (S78)",
                             "CB4856\nDBT-1(C78)" ),
                   ordered = T))+
  geom_boxplot()+
  geom_point()+
  scale_fill_manual(values=c("orange","gray50","blue","gray50"))+
  arsenic.theme + 
  facet_grid(ratio~., scales = "free")+
  labs(y = "Ratio Value")+
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold",vjust = 1))

ggsave(glue::glue("{resub.plots}FS22_L4_C15iso_to_C18_ratio.pdf"), 
       height = 6, 
       width = 8)

########################################################################################################################
# PROCESS DBT-1 C78S C15ISO RESCUE PHENOTYPES - START
########################################################################################################################

raw <- easysorter::read_data(paste0(data.dir,"RESCUE/20171024_arsenicsupp/"))

raw_nocontam <- remove_contamination(raw)
summedraw <- sumplate(raw_nocontam, quantiles = TRUE)

biopruned <- bioprune(summedraw)%>%
  dplyr::select(-contains("green"),-contains("red"),-contains("yellow"), -contains("var"), -contains("iqr"), -contains("cv")) %>%
  tidyr::gather(trait, value, -date, -experiment, -round, -assay, -plate, -condition, -control, -strain, -row, -col )%>%
  dplyr::rename(phenotype = value)%>%
  dplyr::mutate(conc = ifelse(grepl("24",condition), ".24 µmol",
                              ifelse(grepl("48",condition), ".48 µmol",
                                     ifelse(grepl("64",condition), ".64 µmol",
                                            ifelse(grepl("100",condition), "1 µmol",
                                                   ifelse(grepl("12",condition), ".12 µmol",NA))))))

repeated_arsenic <- list()
supp_concentrations <- c(".24 µmol",".48 µmol",".64 µmol","1 µmol",".12 µmol")

for(i in 1:length(supp_concentrations)){
  repeated_arsenic[[i]] <- dplyr::filter(biopruned, condition %in% c("EtOH","Arsenic"))%>%
    dplyr::mutate(conc = supp_concentrations[i])
}

repeated_arsenic_df <- dplyr::bind_rows(repeated_arsenic)

rescue_processed <- dplyr::filter(biopruned, !(condition %in% c("EtOH","Arsenic")))%>%
  dplyr::ungroup() %>%
  dplyr::bind_rows(.,repeated_arsenic_df) %>%
  dplyr::select(-assay) %>%
  dplyr::mutate(u_strain = paste(strain, condition, plate, row, col, conc, sep = "_")) %>%
  dplyr::select(Strain = strain, Condition = condition, Trait = trait, Value = phenotype, row, col, control, u_strain, conc) %>%
  tidyr::spread(Trait, Value) %>%
  na.omit()

rescue_arsenic_pca <- data.frame(scale(rescue_processed[,8:ncol(rescue_processed)]))
rescue_arsenic_pca_analysis <- prcomp(rescue_arsenic_pca)

pc_ve_90 <- which(cumsum(rescue_arsenic_pca_analysis$sdev^2/sum(rescue_arsenic_pca_analysis$sdev^2)) > 0.9)[1]

rescue_pc_df <- data.frame(Strain = rescue_processed$Strain,
                           Condition = rescue_processed$Condition,
                           u_strain = rescue_processed$u_strain,
                           row = rescue_processed$row,
                           col = rescue_processed$col,
                           control = rescue_processed$control,
                           conc = rescue_processed$conc,
                           rescue_arsenic_pca_analysis$x[,1:pc_ve_90]) %>%
  tidyr::gather(Trait, Value, -Strain, -Condition, -u_strain, -row, -col, -control, -conc)

rescue_all_traits <- dplyr::filter(biopruned, !(condition %in% c("EtOH","Arsenic")))%>%
  dplyr::ungroup() %>%
  dplyr::bind_rows(.,repeated_arsenic_df) %>%
  dplyr::select(-assay) %>%
  dplyr::mutate(u_strain = paste(strain, condition, plate, row, col, conc, sep = "_")) %>%
  dplyr::select(Strain = strain, Condition = condition, u_strain,row, col,control, Trait = trait, Value = phenotype, conc) %>%
  na.omit()  %>%
  dplyr::bind_rows(., rescue_pc_df)

# extract PC loadings
rescue_rotation_df <- data.frame(trait = row.names(rescue_arsenic_pca_analysis$rotation),
                                 rescue_arsenic_pca_analysis$rotation[,1:pc_ve_90])

rescue_rotation_df <- rescue_rotation_df %>%
  tidyr::gather(PC, loading, -trait) %>%
  dplyr::rename(Trait = trait)

# Find trait correlations
rescue_trait_correlations <- cor(rescue_processed[,8:ncol(rescue_processed)])%>%
  data.frame() 

rescue_trait_correlations$trait_b <- row.names(rescue_trait_correlations)

rescue_trait_correlations <- rescue_trait_correlations%>%
  tidyr::gather(trait_a, trait_cor, -trait_b) 


write.table(rescue_all_traits, 
            file = glue::glue("{resub.data}Supplemental_Data29_Rescue_control_regressed.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F )

write.table(rescue_trait_correlations, 
            file = glue::glue("{resub.data}Supplemental_Data30_Rescue_Trait_Correlations.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F )

write.table(rescue_rotation_df, 
            file = glue::glue("{resub.data}Supplemental_Data31_Rescue_PC_Loadings.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F )

########################################################################################################################
# PROCESS DBT-1 C78S C15ISO RESCUE PHENOTYPES - END
########################################################################################################################


########################################################################################################################
# DBT-1 C78S C15ISO RESCUE PHENOTYPES - PLOT PC1 v Brood and Length - START
########################################################################################################################

rescue_pheno <- data.table::fread( glue::glue("{resub.data}Supplemental_Data29_Rescue_control_regressed.tsv"))

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
  geom_text(aes(15, 9, 
                label = paste("Adj R2 = ", adj.r.squared, "\n",
                              "Intercept =",intercept, "\n",
                              "Slope =", slope, "\n",
                              "P =", pvalue)),
            data = dplyr::filter(rescue_lm, CorTrait == "norm.n")) +
  geom_text(aes(200, -5, 
                label = paste("Adj R2 = ", adj.r.squared, "\n",
                              "Intercept =",intercept, "\n",
                              "Slope =", slope, "\n",
                              "P =", pvalue)), 
            data = dplyr::filter(rescue_lm, CorTrait == "median.TOF"))+
  theme(panel.grid.minor = element_blank())


ggsave(glue::glue("{resub.plots}FS24_Rescue_PC1_Brood_Size_Correlation.pdf"), 
       height = 10, 
       width = 16)

########################################################################################################################
# DBT-1 C78S C15ISO RESCUE PHENOTYPES - PLOT PC1 v Brood and Length - END
########################################################################################################################

########################################################################################################################
# DBT-1 C78S C15ISO RESCUE PHENOTYPES - PLOT TRAIT CORRELATIONS AND LOADINGS - START
########################################################################################################################

rescue_trait_cor <- data.table::fread( glue::glue("{resub.data}Supplemental_Data30_Rescue_Trait_Correlations.tsv"))
rescue_pc_loading <- data.table::fread( glue::glue("{resub.data}Supplemental_Data31_Rescue_PC_Loadings.tsv"))

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

ggsave(glue::glue("{resub.plots}FS23_Rescue_Trait-correlation_loadings.pdf"), 
       height = 12, 
       width = 20)

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

ggsave(glue::glue("{resub.plots}FS25_Rescue_PC1_full_experiment.pdf"), 
       height = 6, 
       width = 18)

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

ggsave(glue::glue("{resub.plots}F4C_Rescue_PC1.pdf"), 
       height = 6, 
       width = 12)

rescue_pheno_pr %>%
  dplyr::filter(Trait == "PC1") %>%
  dplyr::filter(Condition %in% c("ArsenicC15ISO64","Arsenic")) %>%
  dplyr::select(Strain, Condition, Value) %>%
  dplyr::distinct() %>%
  dplyr::filter(Strain %in% c("N2","ECA581")) %>%
  dplyr::group_by(Strain, Condition) %>%
  dplyr::summarise(mv = mean(Value))

# c15iso rescues 56.3% of the phenotypic difference, between n2 and n2 swap 
(5.82 - 3.91)/(5.82 - 2.43)

rescue_pheno_pr %>%
  dplyr::filter(Trait == "PC1") %>%
  dplyr::filter(Condition %in% c("ArsenicC15ISO64","Arsenic")) %>%
  dplyr::select(Strain, Condition, Value) %>%
  dplyr::distinct() %>%
  dplyr::filter(Strain %in% c("CB4856","ECA590")) %>%
  dplyr::group_by(Strain, Condition) %>%
  dplyr::summarise(mv = mean(Value))

# c15iso rescues 35.9% of the phenotypic difference , between cb and cb swap
abs(6.69 - 5.21)/abs(6.69 - 2.57)

########################################################################################################################
# DBT-1 C78S C15ISO RESCUE PHENOTYPES - BOXPLOTS - END
########################################################################################################################



########################################################################################################################
# PROCESS AND PLOT HUMAN COMPETITION EXPERIMENT - START
########################################################################################################################

human_reads <- readr::read_csv(glue::glue("{resub.data}Supplemental_Data32_Cell_Line_Swap.csv"))

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

ggsave(glue::glue("{resub.plots}F5B_Human_Cell_Barplot.pdf"), 
       height = 6, 
       width = 8)

human_reads <- human_reads %>%
  dplyr::mutate(id = paste(Edit, Replicate, Arsenic_Concentration,c(1:2), sep = "_")) %>%
  na.omit()

fet_ls <- list()
for(ft in 1:length(unique(human_reads$id))) {
  
  expid <- unique(human_reads$id)[ft]
  
  aa_edit <- strsplit(expid, split = "_")[[1]][1]
  replicate_n <- strsplit(expid, split = "_")[[1]][2]
  arsenic_c <- strsplit(expid, split = "_")[[1]][3]
  well_rep <- as.numeric(strsplit(expid, split = "_")[[1]][4])
  
  c_wt <- dplyr::filter(human_reads, 
                        Edit == aa_edit,
                        Replicate == replicate_n,
                        Arsenic_Concentration == "0uM") %>% dplyr::pull(Guide1_wt)
  
  c_ed <- dplyr::filter(human_reads, 
                        Edit == aa_edit,
                        Replicate == replicate_n,
                        Arsenic_Concentration == "0uM") %>% dplyr::pull(Guide1_edit)
  
  a_wt <- dplyr::filter(human_reads, 
                        Edit == aa_edit,
                        Replicate == replicate_n,
                        Arsenic_Concentration == arsenic_c) %>% dplyr::pull(Guide1_wt)
  
  a_ed <- dplyr::filter(human_reads, 
                        Edit == aa_edit,
                        Replicate == replicate_n,
                        Arsenic_Concentration == arsenic_c) %>% dplyr::pull(Guide1_edit)
  
  test_table <- matrix(c(c_wt[well_rep],c_ed[well_rep],
                         a_wt[well_rep],a_ed[well_rep]),
                       nrow = 2,
                       dimnames = list(c("Control", "Arsenic"),
                                       c("Unedited", "Edited")))
  
  fet_ls[[ft]] <- data.frame(p_val = fisher.test(test_table, alternative = "greater")$p.value,
                             Edit = aa_edit,
                             Conc = arsenic_c,
                             Rep = well_rep)
  
}

fet_df <- dplyr::bind_rows(fet_ls)

write.table(fet_df, 
            file = glue::glue("{resub.data}Supplemental_Data33_Human_Arsenic_Exact_test.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F )

########################################################################################################################
# PROCESS HUMAN COMPETITION EXPERIMENT - END
########################################################################################################################


# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ##
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## 
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ # Tajima's D across GWAS interval # ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ #
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## 
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ##

# downloaded CeNDR release 20160408 and renamed TS23_WI.20160408.impute.vcf.gz
# 
# dbt <- tajimas_d_temp(vcf_path = paste0(final.dir),
#                       vcf_name = "TS24_WI.20160408.impute.vcf.gz", 
#                       chromosome = "II", 
#                       interval_start = 7.43e6, 
#                       interval_end = 8.33e6,
#                       site_of_interest = 7942463,
#                       slide_distance = 10,
#                       window_size = 500)
# 
# write.table(dbt[[1]],
#             file = paste0(final.dir,"TS25_TajimasD.tsv"),
#             sep = "\t",
#             quote = F,
#             col.names = T,
#             row.names = F)

# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## 
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ##
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ #  C78S world-wide allele distribution # ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~  # 
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## 
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## 

df <- gs_key("1V6YHzblaDph01sFDI8YK_fP0H7sVebHQTXypGdiQIjI") %>%
  gs_read()

new_pheno <- readr::read_tsv(file = glue::glue("{resub.data}Supplemental_Data17_GWAS_ALL_TRAITS.tsv"))%>%
  as.data.frame() %>%
  tidyr::gather(trait,value, -strain)

dbt_c78 <- cegwas2::query_vcf("dbt-1") %>%
  dplyr::filter(aa_change == "p.Cys78Ser") %>%
  dplyr::filter(SAMPLE %in% new_pheno$strain) %>%
  dplyr::mutate(GT = ifelse(a1 == REF, "REF", "ALT"))%>%
  dplyr::select(strain = SAMPLE, GT )

isolation_info <-df%>%
  dplyr::filter(reference_strain == 1)%>%
  dplyr::select(strain = isotype, long = longitude, lat = latitude)%>%
  dplyr::filter(strain %in% new_pheno$strain)%>%
  dplyr::filter(lat != "None")%>%
  dplyr::left_join(dbt_c78,.,by="strain")%>%
  dplyr::filter(!is.na(lat))

isolation_info$lat <- as.numeric(isolation_info$lat)
isolation_info$long <- as.numeric(isolation_info$long)

write.table(fet_df, 
            file = glue::glue("{resub.data}Supplemental_Data35_Strain_Isolation_Info.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F )


isolation_info <- data.table::fread(glue::glue("{resub.data}Supplemental_Data35_Strain_Isolation_Info.tsv"))

world <- map_data(map = "world")
world <- world[world$region != "Antarctica",] # intercourse antarctica

# world pop
ggplot()+ geom_map(data=world, map=world,
                   aes(x=long, y=lat, map_id=region),
                   color="white", fill="#7f7f7f", size=0.05, alpha=1/4)+
  geom_point(data = isolation_info, aes(x=long, y=lat, fill=GT), shape =21, alpha = 0.7) + 
  scale_fill_manual(values = c("blue", "orange"),name = "Cys78Ser")+
  theme_map()

ggsave(glue::glue("{resub.plots}FS28_World_Distribution.pdf"), 
       height = 5, 
       width = 10)