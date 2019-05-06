# Figure script for arsenic manuscript

# set to location of files
main.dir <- "~/Dropbox/AndersenLab/LabFolders/Stefan/Manuscripts/Arsenic/"
data.dir <- paste0(main.dir,"Data/")
final.dir <- paste0(main.dir,"Final_tables/")
plot.dir <- paste0(main.dir,"Plots/")
script.dir <- paste0(main.dir,"Scripts/")

source(paste0(script.dir,"processing_functions.R"))

trait_of_interest <- "median.TOF"
comparison_trait <- "q90.TOF"
condition_of_interest <- "arsenictrioxide"
control_oof_interest <- "1percwater"

# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ #
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ #
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ # DOSE RESPINSE EXPERIMENT # ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ #
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ #
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ #

dirs <- paste0(data.dir,"DoseResponse_new/")

raw <- read_data(dirs)
raw_nocontam <- remove_contamination(raw)
summedraw <- sumplate(raw_nocontam, directories = FALSE, quantiles = TRUE)
biopruned <- bioprune(summedraw)

write.table(biopruned, 
            file = paste0(final.dir,"TS1_DR_Summarized.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F)

# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ #
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ #
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ # Linkage mapping # ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ # 
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ # PROCESS RIAIL DATA # ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ # 
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ #
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ #

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

write.table(arsenic_linkage_raw, 
            file = paste0(final.dir,"TS2_LINKAGE_rawSorter.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F)

#Remove contamination
raw_noncontam <- remove_contamination(raw)
#Summarize data
summary <- sumplate(raw_noncontam, directories=T, quantiles=T)
#Prune based on biological impossibilities
biopruned <- bioprune(summary)
#Regress out the effect of the assay
assayregressed <- regress(biopruned, assay=TRUE)

write.table(assayregressed, 
            file = paste0(final.dir,"TS3_LINKAGE_assayRegressed.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F)

#Prune based on bins
RIAILs2BAMFpruned <- bamf_prune(assayregressed, drop=T)
#Regress out control strains
RIAILs2regressed <- regress(RIAILs2BAMFpruned)

arsenic_linkage_processed <- RIAILs2regressed %>%
  dplyr::filter(condition == condition_of_interest)%>%
  dplyr::mutate(n_strain = as.numeric(gsub("QX","",strain)))%>%
  dplyr::filter(n_strain>239)%>%
  dplyr::select(-n_strain)

#Save final regressed df and be able to load it in the future
write.table(arsenic_linkage_processed, 
            file = paste0(final.dir,"TS4_LINKAGE_controlRegressed.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F)

# # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## #  PERFORM LINKAGE MAPPING
arsenic_linkage <- RIAILs2regressed %>%
  dplyr::filter(condition == condition_of_interest, trait == trait_of_interest | trait == "norm.n")
data("N2xCB4856cross")
blankcross <- N2xCB4856cross
arsenic_cross <- linkagemapping::mergepheno(blankcross, arsenic_linkage, set = 2)
map <- linkagemapping::fsearch(arsenic_cross, permutations = 1000, thresh = "GWER")
annotatedmap <- linkagemapping::annotate_lods(map, arsenic_cross)

write.table(annotatedmap, 
            file = paste0(final.dir,"TS5_LINKAGE_annotatedLODs_with_brood.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F)

# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ #
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ #
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ # NIL PHENOTYPES # ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ # 
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ #
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ #
# 1000µM arscenic trioxide

raw <- easysorter::read_data(paste0(data.dir,"NILs_n_SWAPs/20170403_dbtswap/"))
raw_save <- list(data.frame(raw[[1]], day = "score"),
                 data.frame(raw[[2]], day = "sort"))%>%
  dplyr::bind_rows()

write.table(raw_save, 
            file = paste0(final.dir,"TS6_NILs_SWAPs_rawSorter.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F)

raw_nocontam <- remove_contamination(raw)

summedraw <- sumplate(raw_nocontam, quantiles = TRUE)%>%
  dplyr::mutate(f.L1L2L3 = f.L1+f.L2L3)

summedraw1 <- dplyr::mutate(summedraw, n_sorted = n/norm.n)%>%
  tidyr::gather(trait, value, -date, -experiment, -round, -assay, -plate, 
                -condition, -control, -strain, -row, -col, -n_sorted )%>%
  dplyr::select(condition, strain,plate, row, col, control,trait, n_sorted)

biopruned <- bioprune(summedraw)

write.table(biopruned, 
            file = paste0(final.dir,"TS7_NILs_SWAPs_Summarized.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F)

controlregressed <- easysorter::regress(biopruned, assay = FALSE)
controlregressed <- dplyr::left_join(controlregressed,summedraw1,
                                     by = c("condition", "strain","plate", "row", "col", 
                                            "control","trait","date","experiment","round","assay"))

write.table(controlregressed, 
            file = paste0(final.dir,"TS8_NILs_SWAPs_controlRegressed.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F)

# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ #
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ #
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ # LINKAGE - NIL QTL SUMMARY # ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ # 
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ #
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ #

annotatedmap <- data.table::fread(paste0(final.dir,"TS5_LINKAGE_annotatedLODs_with_brood.tsv"))

# define brood size QTL region
region <- paste0(na.omit(dplyr::filter(annotatedmap, trait == "arsenictrioxide.norm.n"))$chr[1],":",
                 na.omit(dplyr::filter(annotatedmap, trait == "arsenictrioxide.norm.n"))$ci_l_pos[1],"-",
                 na.omit(dplyr::filter(annotatedmap, trait == "arsenictrioxide.norm.n"))$ci_r_pos[1])

# 1243 in brood size qtl
# 683 with mod-high
# 91 with high
# 15769 SNPs
linkage_summary <- interval_summary(region)

# define left-right bounds of median.TOF QTL
left <- na.omit(annotatedmap)$ci_l_pos[1]
right <- na.omit(annotatedmap)$ci_r_pos[1]

cb_variants <- snpeff(region, severity = "ALL")

cb_variants<-cb_variants %>%
  dplyr::filter(strain=="CB4856", GT =="ALT")

# 567 genes with variation in CB for larger brood size QTL
length(unique(cb_variants$gene_id))

size_cb_variants <- cb_variants %>%
  dplyr::filter(strain=="CB4856", GT =="ALT", POS > left, POS < right)

# 373 genes with variants in CB for animal length
length(unique(size_cb_variants$gene_id))

cb_variants_mod<-cb_variants %>%
  dplyr::filter(strain=="CB4856", GT =="ALT", impact %in% c("MODERATE", "HIGH"))

# 187 genes with high-moderate effect in larger interval
length(unique(cb_variants_mod$gene_id))

size_cb_variants_mod<-cb_variants_mod %>%
  dplyr::filter(strain=="CB4856", GT =="ALT", impact %in% c("MODERATE", "HIGH"), POS > left, POS < right)

# 187 genes with high-moderate effect in larger interval
length(unique(size_cb_variants_mod$gene_id))

# 127 genes with high-moderate effect in larger interval
length(unique(cb_variants_mod$gene_id))

size_nil<-cb_variants %>%
  dplyr::filter(strain=="CB4856", GT =="ALT", POS > 7.83e6, POS < 8.02e6)

# 30 genes with high-moderate effect in larger interval
length(unique(size_nil$gene_id))

size_nil<-cb_variants %>%
  dplyr::filter(strain=="CB4856", GT =="ALT", impact %in% c("MODERATE", "HIGH"), POS > 7.83e6, POS < 8.02e6)

# 30 genes with high-moderate effect in larger interval
length(unique(size_nil$gene_id))

# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ #
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ #
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ # GWAS PHENOTYPES # ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ # 
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ #
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ #

g8 <- easysorter::read_data(c(paste0(data.dir,"GWAS/20140617_arsenicGWAS8a/"),
                              paste0(data.dir,"GWAS/20140617_arsenicGWAS8b/")))

arsenic_gwas_raw <- list(data.frame(g8[[1]][[1]], day = "score"),
                         data.frame(g8[[1]][[2]], day = "sort"),
                         data.frame(g8[[2]][[1]], day = "score"),
                         data.frame(g8[[2]][[2]], day = "sort"))%>%
  dplyr::bind_rows()

write.table(arsenic_gwas_raw , 
            file = paste0(final.dir,"TS10_GWAS_rawSorter.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F)

raw_nocontam <- easysorter::remove_contamination(g8)
summedraw_a <- easysorter::sumplate(raw_nocontam, quantiles = TRUE, directories = T)
summedraw_a <- easysorter::bioprune(summedraw_a)
summedraw_a <- easysorter::bamf_prune(summedraw_a)

g8_prune <- summedraw_a%>%
  dplyr::filter(bamfoutlier1 != TRUE,bamfoutlier2 != TRUE,bamfoutlier2 != TRUE)%>%
  dplyr::select(-bamfoutlier1,-bamfoutlier2,-bamfoutlier3)

write.table(g8_prune, 
            file = paste0(final.dir,"TS11_GWAS_Summarized_Pruned.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F)

# # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # assay regression 

g8_prune_assay <- easysorter::regress(g8_prune, assay = TRUE)

write.table(g8_prune_assay, 
            file = paste0(final.dir,"TS12_GWAS_assayRegressed.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F)

# # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # control regression

g8_regressed <- easysorter::regress(g8_prune_assay, assay = F)

write.table(g8_regressed, 
            file = paste0(final.dir,"TS13_GWAS_controlRegressed.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F)

# # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # ## ## ## ## # PERFORM MAPPINGS

arsenic_raw <- g8_regressed%>%
  dplyr::filter(trait == trait_of_interest | trait == "norm.n")%>%
  dplyr::ungroup()%>%
  dplyr::select(strain,trait, phenotype)%>%
  tidyr::spread(trait, phenotype)

# colnames(arsenic_raw) <- c("strain","arsenic")

pr_pheno <- cegwas::process_pheno(arsenic_raw)
maps <- cegwas::gwas_mappings(pr_pheno, mapping_snp_set = F)
pr_maps <- cegwas::process_mappings(maps, pr_pheno)

write.table(pr_maps,
            file = paste0(final.dir,"TS14_GWAS_Brood_Size_processed_mappings.tsv"),
            sep = "\t",
            quote = F,
            col.names = T,
            row.names = F)

# # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # ## ## ## ## # PCA MAPPINGS

arsenic_pca <- g8_regressed%>%
  dplyr::ungroup()%>%
  dplyr::select(strain, trait, phenotype)%>%
  tidyr::spread(trait, phenotype)

row.names(arsenic_pca) <- arsenic_pca$strain
strain_order <- arsenic_pca$strain

arsenic_pca <- arsenic_pca%>%
  dplyr::select(-strain, -q90.norm.green, -q90.norm.yellow, 
                -q25.norm.red, -q10.yellow, -q10.norm.yellow, -q10.norm.red, -mean.norm.green, -mean.norm.yellow)


scaled_arsenic_pca <- data.frame(scale(arsenic_pca))
arsenic_pca_analysis <- prcomp(scaled_arsenic_pca)

arsenic_pca_df <- data.frame(arsenic_pca_analysis$x[,1:10])
arsenic_pca_df$strain <- strain_order
arsenic_pca_df <- dplyr::select(arsenic_pca_df, strain, PC1:PC10)

write.table(arsenic_pca_df,
            file = paste0(final.dir,"TS15_GWAS_PCA_trait_DF.tsv"),
            sep = "\t",
            quote = F,
            col.names = T,
            row.names = F)

vars <- apply(arsenic_pca_analysis$x, 2, var)  
props <- vars / sum(vars)
pca_var_df <- data.frame(PC = as.numeric(gsub("PC","",names(props))), Variance = cumsum(props))

write.table(pca_var_df,
            file = paste0(final.dir,"TS16_GWAS_PCA_VE.tsv"),
            sep = "\t",
            quote = F,
            col.names = T,
            row.names = F)

pca_pheno <- cegwas::process_pheno(arsenic_pca_df)
pca_maps <- cegwas::gwas_mappings(pca_pheno, mapping_snp_set = F)
pr_pca_maps <- cegwas::process_mappings(pca_maps, pca_pheno)

pr_pca_maps_pc1 <- dplyr::filter(pr_pca_maps, trait == "pc1")

write.table(pr_pca_maps_pc1,
            file = paste0(final.dir,"TS17_GWAS_PCA_processed_mappings.tsv"),
            sep = "\t",
            quote = F,
            col.names = T,
            row.names = F)

# PC1 fine-mapping
pca_genes <- cegwas::process_correlations(cegwas::variant_correlation(pr_pca_maps,
                                                                      variant_severity = "ALL",
                                                                      condition_trait = F))

pca_genes_pc1 <- dplyr::filter(pca_genes, trait == "pc1")

write.table(pca_genes_pc1,
            file = paste0(final.dir,"TS18_GWAS_PCA_fine_mapping.tsv"),
            sep = "\t",
            quote = F,
            col.names = T,
            row.names = F)


# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ #
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ #
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ # SWAP PCA ANALYSIS # ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ # 
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ #
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ #

controlregressed <- data.table::fread(file = paste0(final.dir,"TS8_NILs_SWAPs_controlRegressed.tsv"))

swap_pc <- controlregressed%>%
  dplyr::select(strain, trait, phenotype,plate,row,col)%>%
  dplyr::mutate(u_strain = paste(strain,plate,row,col,sep="_"))%>%
  tidyr::spread(trait,phenotype)

scaled_arsenic_pca <- data.frame(scale(swap_pc[,6:ncol(swap_pc)]))
swap_arsenic_pca_analysis <- prcomp(scaled_arsenic_pca)
swap_pc_df <- data.frame(strain = swap_pc$strain,
                         plate = swap_pc$plate,
                         row = swap_pc$row,
                         col = swap_pc$col,
                         swap_arsenic_pca_analysis$x[,1:10])

write.table(swap_pc_df,
            file = paste0(final.dir,"TS19_SWAP_PCA_Phenotypes.tsv"),
            sep = "\t",
            quote = F,
            col.names = T,
            row.names = F)

# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ #
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ #
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ # C15ISO RESCUE EXPERIMENT  # ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ # 
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ #
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ #

raw <- easysorter::read_data(paste0(data.dir,"RESCUE/20171024_arsenicsupp/"))
raw_save <- list(data.frame(raw[[1]], day = "score"),
                 data.frame(raw[[2]], day = "sort"))%>%
  dplyr::bind_rows()

write.table(raw_save, 
            file = paste0(final.dir,"TS20_RESCUE_rawSorter.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F)

raw_nocontam <- remove_contamination(raw)
summedraw <- sumplate(raw_nocontam, quantiles = TRUE)

biopruned <- bioprune(summedraw)%>%
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

processed_biopruned <- dplyr::filter(biopruned, !(condition %in% c("EtOH","Arsenic")))%>%
  dplyr::bind_rows(.,repeated_arsenic_df)

write.table(processed_biopruned, 
            file = paste0(final.dir,"TS21_RESCUE_processed.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F)

# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ #
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ #
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ # Metabolite Measurements  # ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ # 
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ #
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ #

metabolites <- data.table::fread(paste0(data.dir, "Metabolite_measurements.tsv"))%>%
  tidyr::separate(Sample, into = c("strain", "replicate", "concentration"), sep = "-")%>%
  dplyr::select(-TIC)%>%
  dplyr::filter(!(strain == "581" & replicate == "C" & concentration == "Mock"))%>%
  dplyr::mutate(rC15 = C15iso/C15n,
                rC17 = C17iso/C17n)%>%
  tidyr::gather(compound, value, -strain, -replicate, -concentration) %>%
  dplyr::filter(concentration != "200")

write.table(metabolites, 
            file = paste0(final.dir,"TS22_Metabolite_Measurements.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F)

# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ##
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## 
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ # Tajima's D across GWAS interval # ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ #
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## 
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ##

# downloaded CeNDR release 20160408 and renamed TS23_WI.20160408.impute.vcf.gz

dbt <- tajimas_d_temp(vcf_path = paste0(final.dir),
                      vcf_name = "Supplemental_data25_WI.20160408.impute.vcf.gz", 
                      chromosome = "II", 
                      interval_start = 7.43e6, 
                      interval_end = 8.33e6,
                      site_of_interest = 7942463,
                      slide_distance = 10,
                      window_size = 50)

write.table(dbt[[1]],
            file = paste0(final.dir,"TS25_TajimasD.tsv"),
            sep = "\t",
            quote = F,
            col.names = T,
            row.names = F)

# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## 
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ##
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ #  C78S world-wide allele distribution # ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~  # 
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## 
# ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## ~ ~ ~ ~ ~ ~ ~ ~~ ~ ~ ~ ~ ~ ~ ~ ## 

g8_regressed <- data.table::fread(paste0(final.dir,"TS13_GWAS_controlRegressed.tsv"))

world <- map_data(map = "world")
world <- world[world$region != "Antarctica",] # intercourse antarctica

dbt_c78 <- cegwas2::query_vcf("dbt-1") %>%
  dplyr::filter(aa_change == "p.Cys78Ser") %>%
  dplyr::filter(SAMPLE %in% g8_regressed$strain) %>%
  dplyr::mutate(GT = ifelse(a1 == REF, "REF", "ALT"))%>%
  dplyr::select(strain = SAMPLE, GT )

df <- gs_key("1V6YHzblaDph01sFDI8YK_fP0H7sVebHQTXypGdiQIjI") %>%
  gs_read()

isolation_info <-df%>%
  dplyr::filter(reference_strain == 1)%>%
  dplyr::select(strain = isotype, long = longitude, lat = latitude)%>%
  dplyr::filter(strain %in% g8_regressed$strain)%>%
  dplyr::filter(lat != "None")%>%
  dplyr::left_join(dbt_c78,.,by="strain")%>%
  dplyr::filter(!is.na(lat))

isolation_info$lat <- as.numeric(isolation_info$lat)
isolation_info$long <- as.numeric(isolation_info$long)

write.table(isolation_info,
            file = paste0(final.dir,"TS26_Strain_Isolation_Info.tsv"),
            sep = "\t",
            quote = F,
            col.names = T,
            row.names = F)
