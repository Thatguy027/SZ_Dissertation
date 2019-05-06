# Table script for ben-1 manuscript


# set to location of files
main.dir <- "~/Dropbox/HTA/Results/2017_anthelmintic_gwa/analysis/ben1_paper/"
data.dir <- paste0(main.dir,"Raw_data/")
plot.dir <- paste0(main.dir,"Plots/")
script.dir <- paste0(main.dir,"Scripts/")
final.dir <- paste0(main.dir,"Final_Tables/")

source(paste0(script.dir,"ben1_processing_functions.R"))

trait_of_interest <- "q90.TOF"
condition_of_interest <- "albendazole"
control_of_interest <- "DMSO"

# # # # # ##  DOSE RESPONSE EXPERIMENT (TS1_DR_raw.Rda, TS2_DR_Processed.tsv) ## # # # # #

dirs <- paste0(data.dir,"DoseResponse/")

raw <- read_data(dirs)

save(raw, file = paste0(final.dir,"TS1_DR_raw.Rda"))

raw_nocontam <- remove_contamination(raw)
summedraw <- sumplate(raw_nocontam, directories = FALSE, quantiles = TRUE)
biopruned <- bioprune(summedraw)

pr_df <- biopruned%>%
  dplyr::ungroup()%>%
  dplyr::select(-date, -round,-experiment,-assay)%>%
  tidyr::gather(trait, value, -condition,-control,-strain, -plate, -row, -col)%>%
  dplyr::filter(!is.na(condition))

write.table(pr_df , 
            file = paste0(final.dir,"TS2_DR_Processed.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F)

# # # # # ## HTA DATA - WILD STRAINS (TS3_HTA_raw.Rda, TS4_HTA_processed.tsv) ## # # # # # 

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

GWA_raw <- c(raw_a,raw_b)

save(GWA_raw, file = paste0(final.dir,"TS3_HTA_raw.Rda"))

raw_ls <-list()
for(blk in 1:length(GWA_raw)){
  raw_tof1 <- dplyr::select(GWA_raw[[blk]][[1]], date:TOF)
  raw_tof2 <- dplyr::select(GWA_raw[[blk]][[2]], date:TOF)
  
  raw_ls[[blk]] <- list(raw_tof1, raw_tof2)
}
save(raw_ls, file = paste0(final.dir,"TS3_HTA_raw_tof.Rda"))

# load(paste0(final.dir,"TS3_HTA_raw.Rda"))

# Remove all data from the contaminated wells
raw_nocontam_a <- easysorter::remove_contamination(raw_a)
raw_nocontam_b <- easysorter::remove_contamination(raw_b)

# Summarize the data
summedraw_a <- easysorter::sumplate(raw_nocontam_a, directories = TRUE, quantiles = TRUE)
summedraw_b <- easysorter::sumplate(raw_nocontam_b, directories = TRUE, quantiles = TRUE)

# Prune based on biological impossibilities
biopruned_a <- easysorter::bioprune(summedraw_a)
biopruned_b <- easysorter::bioprune(summedraw_b)

# removes any remaining wash wells. 
all_data <- dplyr::bind_rows(list(biopruned_a,biopruned_b))%>%
  dplyr::filter(!is.na(strain))

# FILTER OUTLIERS
long_data <- all_data %>% 
  tidyr::gather(., trait, phenotype, -date, -experiment,-round,-assay,-plate,-condition,-control,-strain, -row, -col)%>%
  data.frame()%>%
  dplyr::filter(!grepl("red|yellow|green", trait))

outlier_replicates <- long_data %>%
  dplyr::group_by(condition, round, trait, strain)%>%
  dplyr::mutate(n_reps = n())%>%
  dplyr::filter(n_reps > 3)%>%
  dplyr::mutate(med_ph = median(phenotype, na.rm = T),
                sd_ph = sd(phenotype, na.rm= T))%>%
  dplyr::rowwise()%>%
  dplyr::mutate(outlier = ifelse((phenotype > med_ph + 1.8*sd_ph) | (phenotype < med_ph - 1.8*sd_ph), "OUTLIER", "OK"))%>%
  ungroup()

# REMOVE OUTLIERS

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

# REGRESSION
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
            file = paste0(final.dir,"TS4_HTA_processed.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F)


# # # # # ## GWA single marker mapping (TS5_GWA_processed_marker_mapping.tsv, TS6_GWA_finemapping_genes_MOD_HIGH.tsv) ## # # # # #

ctrl_regressed <- data.table::fread(paste0(final.dir,"TS4_HTA_processed.tsv"))

gwa_raw_traits <- ctrl_regressed%>%
  tidyr::gather(cond_trait, phenotype, -strain)%>%
  tidyr::separate(cond_trait, into = c("condition", "trait"), sep = "_")%>%
  dplyr::select(-condition)%>%
  dplyr::filter(trait == trait_of_interest | trait == "norm.n")%>%
  tidyr::spread(trait, phenotype)

pr_pheno <- cegwas::process_pheno(gwa_raw_traits)
maps <- cegwas::gwas_mappings(pr_pheno)
pr_maps <- cegwas::process_mappings(maps, phenotype_df = pr_pheno)

brood_prmaps <- maps %>%
  dplyr::filter(trait == "norm.n")%>%
  cegwas::process_mappings(., phenotype_df = pr_pheno)

comb_maps <- dplyr::rbind_all(list(pr_maps, brood_prmaps))

write.table(comb_maps, 
            file = paste0(final.dir,"TS5_GWA_processed_marker_mapping.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F)

genes <- process_correlations(variant_correlation(pr_maps, condition_trait = F))

write.table(genes, 
            file = paste0(final.dir,"TS6_GWA_finemapping_genes_MOD_HIGH.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F)

fine_map_small <- dplyr::distinct(genes, CHROM, POS, REF, ALT, nt_change, aa_change, gene_name, effect, corrected_spearman_cor_p)

write.table(fine_map_small, 
            file = paste0(final.dir,"TS6_GWA_finemapping_genes_MOD_HIGH_small.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F)

# # # # # ## BURDEN TESTING (TS7_burden_test_mappings.Rda) ## # # # # #

phenos <- data.table::fread(paste0(final.dir, "TS4_HTA_processed.tsv"))
traits <- phenos%>%
  dplyr::select(strain,Albendazole_q90.TOF)%>%
  tidyr::gather(trait, value, -strain)%>%
  dplyr::mutate(Fam = "elegans", Sample = strain, Paternal = 0, Maternal = 0, Sex = 2)%>%
  tidyr::spread(trait,value)%>%
  dplyr::select(-strain)

traits[is.na(traits)] <- -9

burden_raw <- dplyr::select(traits, Fam:Sex, Albendazole_q90.TOF)

write.table(burden_raw, 
            file = paste0(final.dir,"TS7_GWAS_",trait_of_interest,".ped"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F)

#vcfdir=/Users/stefan/Dropbox/AndersenLab/LabFolders/Stefan/burden/snp_indel_snpEff.vcf.gz
#refflatdir=/Users/stefan/Dropbox/AndersenLab/LabFolders/Stefan/burden/refFlat.ws245.txt

# docker run -i -t -v $(pwd):/home -w /home zhanxw/rvtests /rvtests/executable/rvtest \
# --pheno ben1_paper/GWAS_q90.TOF.ped --inVcf snp_indel_snpEff.vcf.gz --freqUpper .05 \
# --freqLower 0.003 --out q90TOF --geneFile refFlat.ws245.txt  --vt price
# 
# docker run -i -t -v $(pwd):/home -w /home zhanxw/rvtests /rvtests/executable/rvtest --pheno ben1_paper/GWAS_q90.TOF.ped --inVcf WI.20170531.impute.vcf.gz --freqUpper .05 --freqLower 0.003 --out q90TOF --geneFile refFlat.ws245.txt  --vt price --kernel skat
# # # #VT PRICE
q90burden <- read.table(paste0(final.dir, "q90TOF.VariableThresholdPrice.assoc"),header = T)
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

write.table(q90burden_pr, 
            file = paste0(final.dir,"TS8_burden_test_mappings.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F)

# # # # # ## MANUAL CURATION OF ben-1 VARIANTS (TS8_ben-1_variants.tsv) ## # # # # # 

# LOAD AND PROCESS INDELS
ben1_indels <- data.table::fread(paste0(data.dir, "ben1_Indels/20171021_ben1_indels.csv"))%>%
  na.omit()

pr_indels <- ben1_indels %>%
  tidyr::unite(marker, Type, Start, End, sep ="_",remove=F)%>%
  dplyr::filter(comments == "frame shift" | grepl("Exon", location) | Type == "trans")%>%
  dplyr::filter(!grepl("Intron",location,ignore.case = T))%>%
  dplyr::mutate(GT = ifelse(marker == "_NA_NA", "REF", "ALT"))%>%
  dplyr::select(marker, strain = Strain, GT,start = Start,end =End)%>%
  dplyr::distinct(strain, marker, GT,.keep_all=T)

# LOAD AND PROCESS SNPS
ben1_snps <- cegwas::snpeff("ben-1",severity = "ALL",elements = "ALL")

ben1_snps_pr <- ben1_snps%>%
  dplyr::select(POS, strain, GT, nt_change, effect, impact)%>%
  tidyr::unite(marker, effect,impact,nt_change,POS,sep="_",remove=F)%>%
  dplyr::select(marker, strain, GT,start = POS)%>%
  dplyr::mutate(end = start)%>%
  dplyr::filter(!grepl("MODIF|syn",marker))

# COMBINE SNPS AND INDELS
ben1_variants <- dplyr::bind_rows(ben1_snps_pr,pr_indels)%>%
  dplyr::filter(GT=="ALT")

write.table(ben1_variants, 
            file = paste0(final.dir,"TS9_ben-1_variants.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F)

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # 
# ben-1 gene tajima d with manually curated variants
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

ben1_indels <- data.table::fread(paste0(data.dir, "ben1_Indels/20171021_ben1_indels.csv"))%>%
  na.omit()
pr_indels <- ben1_indels %>%
  tidyr::unite(marker, Type, Start, End, sep ="_",remove=F)%>%
  dplyr::filter(comments == "frame shift" | grepl("Exon", location) | Type == "trans")%>%
  dplyr::filter(!grepl("Intron",location,ignore.case = T))%>%
  dplyr::mutate(GT = ifelse(marker == "_NA_NA", "REF", "ALT"))%>%
  dplyr::select(marker, strain = Strain, GT,start = Start,end =End)%>%
  dplyr::distinct(strain, marker, GT,.keep_all=T)


ben1_snps <- cegwas::snpeff("ben-1",severity = "ALL",elements = "ALL")
ben1_snps_pr <- ben1_snps%>%
  dplyr::select(POS, strain, GT, nt_change, effect, impact)%>%
  tidyr::unite(marker, effect,impact,nt_change,POS,sep="_",remove=F)%>%
  dplyr::select(marker, strain, GT,start = POS)%>%
  dplyr::mutate(end = start)%>%
  dplyr::filter(!grepl("MODIF|syn",marker))

# # # COMBINE SNPS AND INDELS
ben1_variants <- dplyr::bind_rows(ben1_snps_pr,pr_indels)

ben1_variants <- dplyr::bind_rows(ben1_snps_pr,pr_indels)%>%
  dplyr::mutate(GTn = ifelse(GT=="ALT",1,0)) %>%
  dplyr::select(-start,-end, -GT)%>%
  tidyr::spread(marker, GTn)%>%
  dplyr::select(-strain)

ben1_variants[is.na(ben1_variants)] <- 0

non_syn <- gene_level_TajimasD(ben1_variants)
non_syn
# > non_syn
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

syn <- gene_level_TajimasD(ben1_variants)
syn
# > syn
# Wattersons_Theta Average_PWD TajimasD Singletons FuLiF
# 1            12.64        4.15    -2.02         42 -6.07


### NEUTRALITY STATS

library(PopGenome)
vcf_name = "WI.20170531.impute.vcf.gz"
slide_distance = 50
window_size = 2000

chr3 <- c(1,13783801)

chroms <- c("I","II","III","IV","V","X")

w.CHROM <- 3

setwd(final.dir)

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
                                             type = 2)

s_gen_stats <- PopGenome::detail.stats(s_gen)
s_gen_stats <- PopGenome::neutrality.stats(s_gen_stats,detail = T)
gene_names <- get.feature.names(object = s_gen_stats,chr = chroms[w.CHROM],gff.file = "protein_coding.gff")
s_gen_stats <- PopGenome::diversity.stats(s_gen_stats,pi = T)
get.sum.data(s_gen_stats)

# . . . . 
# . . . . SAVE GENOME OBJECT
# . . . . 

save(s_gen_stats,file = paste0(final.dir,"TS10_ben1_popgen.Rda"))

# . . . . 
# . . . . DEFINE SLIDE WINDOWS
# . . . . 

load(paste0(final.dir,"20180602ben1_popgenome.Rda"))

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

write.table(plt_df, 
            file = paste0(final.dir,"TS11_ben1_fayH_tajimaD_zengE.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F)

# # # # # ## CRISPR F200Y vs ben-1 deletion (TSxx_HTA_ben-1_raw.Rda, TSxx_HTA_ben-1_regressed.tsv, (TSxx_competition_assay_raw.csv), TSxx_competition_assay_processed.tsv) ## # # # # # 

# read ASSAY 
HTA_raw<- read_data(paste0(data.dir,"20180312_ben1_edits/"))

# Remove all data from the contaminated wells
raw_nocontam <- remove_contamination(HTA_raw)

# Summarize the data
summedraw <- sumplate(raw_nocontam, quantiles = TRUE)

# too big for plos pathogens
# save(HTA_raw, file = paste0(final.dir,"TS12_HTA_ben-1_raw.Rda"))

write.table(summedraw, 
            file = paste0(final.dir,"TS12_HTA_ben-1_summarized.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F)

# Prune based on biological impossibilities
biopruned <- bioprune(summedraw)%>%
  dplyr::filter(strain %in% c("N2","ECA882","ECA883","ECA884","ECA917","ECA918","ECA919","ECA920","ECA921" )) %>%
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

# control regression

assayregressed <- easysorter::regress(biopruned, assay = FALSE)%>%
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
            file = paste0(final.dir,"TS13_HTA_ben-1_regressed.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F)

# competition assay 

competition_assay <- data.table::fread(file = paste0(data.dir,"competition_assay.csv"))%>% ### TSxx_competition_assay_raw.csv
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


write.table(competition_assay_outliers, 
            file = paste0(final.dir,"TS14_competition_assay_processed.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F)

# # # # # ## REGRESSION OF ben-1 variation and remapping  (TSxx_GWA_marker_mappings_ben1_regressed.tsv, TSxx_GWA_marker_fine_mappings_ben1_regressed.tsv, TSxx_burden_mapping_ben-1_regressed.Rda)## # # # # # 

# LOAD AND PROCESS INDELS
ben1_indels <- data.table::fread(paste0(data.dir, "ben1_Indels/20171021_ben1_indels.csv"))%>%
  na.omit()

pr_indels <- ben1_indels %>%
  tidyr::unite(marker, Type, Start, End, sep ="_",remove=F)%>%
  dplyr::filter(comments == "frame shift" | grepl("Exon", location) | Type == "trans")%>%
  dplyr::filter(!grepl("Intron",location,ignore.case = T))%>%
  dplyr::mutate(GT = ifelse(marker == "_NA_NA", "REF", "ALT"))%>%
  dplyr::select(marker, strain = Strain, GT,start = Start,end =End)%>%
  dplyr::distinct(strain, marker, GT,.keep_all=T)

# LOAD AND PROCESS SNPS
ben1_snps <- cegwas::snpeff("ben-1",severity = "ALL",elements = "ALL")

ben1_snps_pr <- ben1_snps%>%
  dplyr::select(POS, strain, GT, nt_change, effect, impact)%>%
  tidyr::unite(marker, effect,impact,nt_change,POS,sep="_",remove=F)%>%
  dplyr::select(marker, strain, GT,start = POS)%>%
  dplyr::mutate(end = start)%>%
  dplyr::filter(!grepl("MODIF|syn",marker))

# COMBINE SNPS AND INDELS
ben1_variants <- dplyr::bind_rows(ben1_snps_pr,pr_indels)%>%
  dplyr::filter(GT=="ALT")

gwa_mappings <- data.table::fread(file = paste0(final.dir,"TS5_GWA_processed_marker_mapping.tsv"))%>%
  na.omit()%>%
  dplyr::mutate(snpGT = ifelse(allele==-1,"REF", "ALT"))%>%
  dplyr::select(snpMarker = marker, strain, trait, value, snpGT)%>%
  dplyr::filter(trait == "q90.tof")%>%
  dplyr::left_join(.,ben1_variants, by = "strain")


regress_ben1 <- gwa_mappings %>%
  dplyr::mutate(covariate = ifelse(is.na(GT) | strain == "CB4856", 0, 1))%>%
  dplyr::distinct(strain, value, covariate)

residual_df <- data.frame(strain = regress_ben1$strain, regressed_trait = residuals(lm(value ~ covariate, data = regress_ben1)), ben1_variant_trait = regress_ben1$covariate)

write.table(residual_df, 
            file = paste0(final.dir,"TS16_ben1_regressed_with_ben1_covariate.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F)

pr_resid_pheno <- cegwas::process_pheno(residual_df)
resid_maps <- cegwas::gwas_mappings(pr_resid_pheno)
pr_resid_maps <- cegwas::process_mappings(resid_maps, phenotype_df = pr_resid_pheno)

pr_resid_maps <- dplyr::filter(pr_resid_maps, trait == "regressed_trait")

write.table(pr_resid_maps, 
            file = paste0(final.dir,"TS17_GWA_marker_mappings_ben1_regressed.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F)

genes <- process_correlations(variant_correlation(dplyr::filter(pr_resid_maps,trait =="regressed_trait"), condition_trait = F,variant_severity = "ALL"))

write.table(genes, 
            file = paste0(final.dir,"TS18_GWA_marker_fine_mappings_ben1_regressed.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F)

# ben-1 regression burden test  

traits <- residual_df%>%
  dplyr::select(-ben1_variant_trait)%>%
  dplyr::mutate(Fam = "elegans", Sample = strain, Paternal = 0, Maternal = 0, Sex = 2)%>%
  dplyr::select(-strain)

traits[is.na(traits)] <- -9

burden_raw <- dplyr::select(traits, Fam:Sex, regressed_trait)

write.table(burden_raw, 
            file = paste0(final.dir,"TS19_GWA_ben1_regressed.ped"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F)

benreg_burden <- read.table(paste0(final.dir, "TS20_GWA_ben1_regressed.VariableThresholdPrice.assoc"),header = T)
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


write.table(benreg_burden_pr, 
            file = paste0(final.dir,"TS21_GWA_ben1_regressed.VariableThresholdPrice_processed.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F)

# beta tubulin tajima's d

# III:3537688..3541628
st = 3537688
en = 3541628
ben1_d <- tajimas_d_temp(vcf_path = paste0(final.dir) ,vcf_name = "WI.20170531.impute.vcf.gz", 
                         chromosome = "III", interval_start = st-1e5, interval_end = en+1e5, site_of_interest = st,
                         slide_distance = 1,window_size = 100, outgroup = "XZ1516")

# III:10740868..10742461
st = 10740868 
en = 10742461
tbb1_d <- tajimas_d_temp(vcf_path = paste0(final.dir) ,vcf_name = "WI.20170531.impute.vcf.gz", 
                         chromosome = "III", interval_start = st-1e5, interval_end = en+1e5 ,site_of_interest = st,
                         slide_distance = 1,window_size = 100, outgroup = "XZ1516")

# III:4015769..4017643
st = 4015769
en = 4017643

tbb2_d <- tajimas_d_temp(vcf_path = paste0(final.dir) ,vcf_name = "WI.20170531.impute.vcf.gz", 
                         chromosome = "III", interval_start = st-1e5, interval_end = en+1e5 ,site_of_interest = st,
                         slide_distance = 1,window_size = 100, outgroup = "XZ1516") 

# X:9434452..9436893
st = 9434452
en = 9436893
tbb4 <- tajimas_d_temp(vcf_path = paste0(final.dir) ,vcf_name = "WI.20170531.impute.vcf.gz",
                       chromosome = "X", interval_start = st-1e5, interval_end = en+1e5 ,site_of_interest = st,
                       slide_distance = 1,window_size = 100, outgroup = "XZ1516") 

#V:12261804..12263368
st = 12261804
en = 12263368

tbb6 <- tajimas_d_temp(vcf_path = paste0(final.dir) ,vcf_name = "WI.20170531.impute.vcf.gz",
                       chromosome = "V", interval_start = st-1e5, interval_end = en+1e5 ,site_of_interest = st,
                       slide_distance = 1,window_size = 100, outgroup = "XZ1516") 

# X:7774859..7776689
st = 7774859
en = 7776689
mec7 <- tajimas_d_temp(vcf_path = paste0(final.dir) ,vcf_name = "WI.20170531.impute.vcf.gz",
                       chromosome = "X", interval_start = st-1e5, interval_end = en+1e5 ,site_of_interest = st,
                       slide_distance = 1,window_size = 100, outgroup = "XZ1516") 

# combine all plots
tubulin_td <- dplyr::rbind_all(list(ben1_d[[1]],tbb1_d[[1]],tbb2_d[[1]],tbb4[[1]],tbb6[[1]],mec7[[1]]))

write.table(tubulin_td, 
            file = paste0(final.dir,"TS22_tubulin_tajima_d.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F)


# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # 
# FIGURE S11 - phenotyping WN2002 
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

raw <- read_data(paste0(raw_data.dir,"20180313_WI_ABZ/"))
raw_nocontam <- remove_contamination(raw)
summedraw <- sumplate(raw_nocontam, quantiles = TRUE)

summedraw_WI <- summedraw%>%
  dplyr::filter(strain %in% c("N2", "JU2141", "WN2002", "JU2581"))

biopruned_WI <- bioprune(summedraw_WI)%>%
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


assayregressed_WI <- easysorter::regress(biopruned_WI, assay = FALSE)%>%
  dplyr::group_by(strain, condition, trait)%>%
  dplyr::mutate(mph = median(phenotype),
                sph = sd(phenotype))%>%
  dplyr::mutate(flag_h = 2*sph+mph,
                flag_l = mph-2*sph)%>%
  dplyr::mutate(cut_h =ifelse(phenotype >= 2*sph+mph, "YES", "NO"),
                cut_l =ifelse(phenotype <= mph-2*sph, "YES", "NO"))%>%
  dplyr::filter(cut_h != "YES" , cut_l !="YES")

write.table(assayregressed_WI, 
            file = paste0(final.dir,"TS23_WI_rephenotype_HTA_processed.tsv"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = F)
