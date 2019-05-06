# Load easysorter
# devtools::install_github("AndersenLab/easysorter")
# devtools::install_github("AndersenLab/COPASutils")
# devtools::install_github("AndersenLab/easysorter")
library(easysorter)
library(broom)
library(ggplot2)
library(tidyr)
library(dplyr)
library(directlabels)

readSorter <- function(file, tofmin=0, tofmax=10000, extmin=0, extmax=10000, reflx=TRUE)  {
  data <- read.delim(file=file, header=T, na.strings=c("n/a"), as.is=T, stringsAsFactors=F)
  if(is.null(data$EXT) & reflx){
    stop("This file appears to have come from a machine with an LP Sampler, not a ReFLx module. Please set `reflx` = FALSE and try again.")
  }
  data <- data[!is.na(data$TOF),]
  data <- data[,!is.na(data[1,])]
  data <- data[(data$TOF>=tofmin & data$TOF<=tofmax) | data$TOF == -1,]
  if(reflx){
    data <- data[(data$EXT>=extmin & data$EXT<=extmax) | data$EXT == -1,]
    data$Column <- as.factor(data$Column)
    data$Row <- as.factor(data$Row)
  } else {
    data <- data[(data$Extinction>=extmin & data$Extinction<=extmax) | data$Extinction == -1,]
    data$Column <- as.factor(data$Column)
    data$Row <- as.factor(data$Row)
  }
  #removed the following two lines because in cases where a plate was interrupted during sort/score and two files exist for a certain plate the Row and Column data were defaulting to start at row A and col 1. This messes with plate stitching.
  #levels(data$Row) <- LETTERS[1:8]
  #levels(data$Column) <- 1:12
  return(data)
}

raw <- read_data("~/Dropbox/HTA/Results/20180312_ben1_edits_A/")

raw_nocontam <- remove_contamination(raw)

summedraw <- sumplate(raw_nocontam, quantiles = TRUE)

save(summedraw, file = "~/Dropbox/HTA/Results/20180312_ben1_edits_A/analysis/data/20180315_easysorter_processed_ben1.Rda")



load("~/Dropbox/HTA/Results/20180312_ben1_edits_A/analysis/data/20180315_easysorter_processed_ben1.Rda")


biopruned_A <- bioprune(summedraw)%>%
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


load("~/Dropbox/HTA/Results/20180313_ben1_edits_B/analysis/data/20180315_easysorter_processed_ben1.Rda")

biopruned_B <- bioprune(summedraw)%>%
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

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # 
# subtract controls 
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

# dr_data <- assayregressed
# 
# 
# del <- c("ECA882","ECA883","ECA884")
# 
# F200Y <- c("ECA917","ECA918","ECA919","ECA920","ECA921")
# 
# N2 <- c("N2")
# 
# PAM <- c("ECA880","ECA881")


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
  

plot_results <- function(df, cond,trt){
  df%>%
    dplyr::filter(trait == trt, condition %in% cond)%>%
    dplyr::mutate(strain1 = factor(group, levels = c("N2","del","PAM","F200Y"),
                                   labels =c("N2","del","PAM","F200Y")))%>%
    ggplot(.) +
    aes(x = factor(strain1), 
        y = subtracted_pheno, 
        fill=group) +
    # geom_jitter(width = 0.25, alpha = .4)+
    geom_boxplot(outlier.colour = NA, alpha = 0.7)+
    # facet_grid(.~condition)+
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
}

plt_trt <- "q90.TOF"
HTA_plot <- plot_results(subtract_controls, c("Albendazole",plt_trt), plt_trt)




##### TukeyHSD ##### 


run_TukeyHSD <- function(df, trt){
  stat_df <- df %>%
    dplyr::filter(trait == trt)%>%
    dplyr::select(group, subtracted_pheno)
  
  aov_res <- aov(stat_df$subtracted_pheno ~ stat_df$group)
  summary(aov_res)
  tuk <- TukeyHSD(aov_res)
  
  psig=as.numeric(apply(tuk$`stat_df$group`[,2:3],1,prod)>=0)+1
  op=par(mar=c(4.2,9,3.8,2))
  plot(tuk,col=psig,yaxt="n")
  for (j in 1:length(psig)){
    axis(2,at=j,labels=rownames(tuk$`stat_df$group`)[length(psig)-j+1],
         las=1,cex.axis=.8,col.axis=psig[length(psig)-j+1])
  }
  par(op)
  
  pwtuk <- TukeyHSD(aov_res)
  
  return(pwtuk)
}

plt_trt <- "q90.TOF"
run_TukeyHSD(subtract_controls,plt_trt)


############ Figure 6 B: competition assay plot ########


df <-CA_v2%>%
  dplyr::filter(!strain == "N2")


CA_plot <-ggplot(df, aes(x=week, y=frac_mean, fill=condition, color=strain))+
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


# ggsave(paste("~/Dropbox/HTA/Results/20180312_ben1_edits_A/analysis/plots/",plt_trt,"_pxg.pdf",sep=""),
#        height = 6,
#        width = 8)
