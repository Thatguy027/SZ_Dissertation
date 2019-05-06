library(easysorter)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cegwas)
library(PopGenome)
library(linkagemapping)
library(directlabels)
library(broom)
library(ggalt)    # devtools::install_github("hrbrmstr/ggalt")
library(ggthemes) # theme_map and tableau colors
library(cowplot)
library(googlesheets)
library(ggbeeswarm)
library(sjstats)
library(car)
library(tidyverse)
library(ggpmisc)
library(lemon)

arsenic.theme <-  theme(axis.text.x = ggplot2::element_text(size = 14),
                        axis.text.y = ggplot2::element_text(size = 14),
                        axis.title.x = ggplot2::element_text(size = 16, face = "bold", color = "black", vjust = -0.3),
                        axis.title.y = ggplot2::element_text(size = 16, face = "bold", color = "black", vjust = -0.3),
                        strip.text.x = element_text(size = 14, face = "bold"))


get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}

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

manplot_edit <- function(plot_df, bf_line_color = "gray") {
  plot_traits <- unique(plot_df$trait)
  plots <- lapply(plot_traits, function(i) {
    plot_df %>%
      dplyr::filter(trait == i) %>%
      dplyr::distinct(marker, .keep_all = T) %>%
      ggplot2::ggplot(.) +
      ggplot2::aes(x = POS/1e6, y = log10p) +
      ggplot2::scale_color_manual(values = c("black","red","blue")) +
      ggplot2::geom_rect(ggplot2::aes(xmin = startPOS/1e6, 
                                      xmax = endPOS/1e6, 
                                      ymin = 0, 
                                      ymax = Inf, 
                                      fill = "blue"), 
                         color = "blue",fill = "cyan",linetype = 2, 
                         alpha=.3)+
      ggplot2::geom_hline(ggplot2::aes(yintercept = BF),
                          color = bf_line_color, 
                          alpha = .75,  
                          size = 1) +
      ggplot2::geom_point( ggplot2::aes(color= factor(aboveBF)) ) +
      ggplot2::facet_grid( . ~ CHROM, scales = "free_x" , space = "free_x") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(size = 16),
                     axis.text.y = ggplot2::element_text(size = 16),
                     axis.title.x = ggplot2::element_text(size = 20, face = "bold", color = "black", vjust = -0.3), 
                     axis.title.y = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                     strip.text.x = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                     strip.text.y = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                     plot.title = ggplot2::element_text(size = 24, face = "bold", vjust = 1), 
                     panel.background = ggplot2::element_rect(color = "black",size = 1.2),
                     legend.position = "none") +
      ggplot2::labs(x = "Genomic Position (Mb)",
                    y = expression(-log[10](italic(p))))
  })
  plots
}

maxlodplot_edit <- function(map){
  map1 <- map %>%
    dplyr::group_by(marker) %>%
    dplyr::filter(lod == max(lod))
  
  cis <- map %>%
    dplyr::group_by(marker) %>%
    dplyr::mutate(maxlod=max(lod))%>%
    dplyr::group_by(iteration) %>%
    dplyr::filter(!is.na(var_exp)) %>%
    dplyr::do(head(., n=1))
  
  
  if(nrow(cis) == 0) {
    plot <- ggplot2::ggplot(map1, ggplot2::aes(x = genotype, y = pheno)) + ggplot2::geom_blank()
    return(plot)
  }
  
  map1 <- cidefiner(cis, map1)
  
  plot <- ggplot2::ggplot(map1) +
    ggplot2::aes(x = pos/1e6, y = lod)
  
  if(nrow(cis) != 0) {
    plot <- plot + 
      ggplot2::geom_ribbon(ggplot2::aes(x = pos/1e6, ymin = 0, ymax = ci_lod),
                           fill = "blue", alpha = 0.5) +
      ggplot2::geom_point(data = cis, ggplot2::aes(x=pos/1e6, y=(1.05*maxlod)),
                          fill ="red", shape=25, size=3.2, show_guide = FALSE) +
      ggplot2::geom_text(data = cis,
                         ggplot2::aes(x=pos/1e6,
                                      y=(ifelse(lod < 10, 1.5*lod, 1.1*lod)),
                                      label = paste0(100*round(var_exp, digits = 4),"%")),
                         colour = "black", size=3)
  }
  
  plot <- plot + ggplot2::geom_line(size = 1, alpha = 0.85) +
    ggplot2::facet_grid(.~chr, scales ="free", space = "free") +
    ggplot2::labs(x = "Genomic Position (Mb)", y = "LOD") +
    ggplot2::scale_colour_discrete(name="Mapping\nIteration") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 16),
                   axis.text.y = ggplot2::element_text(size = 16),
                   axis.title.x = ggplot2::element_text(size = 20, face = "bold", color = "black", vjust = -0.3), 
                   axis.title.y = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                   strip.text.x = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                   strip.text.y = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                   plot.title = ggplot2::element_text(size = 24, face = "bold", vjust = 1), 
                   panel.background = ggplot2::element_rect(color = "black",size = 1.2))
  return(plot)
}

maxlodplot_edit2 <- function(map){
  
  traits <- unique(map$trait)
  
  map_trait_1 <- map %>%
    dplyr::group_by(marker)%>%
    dplyr::filter(trait == traits[1])%>%
    dplyr::filter(lod == max(lod))
  
  cis_trait_1 <- map %>%
    dplyr::group_by(marker) %>%
    dplyr::filter(trait == traits[1])%>%
    dplyr::mutate(maxlod=max(lod))%>%
    dplyr::group_by(iteration) %>%
    dplyr::filter(!is.na(var_exp)) %>%
    dplyr::do(head(., n=1))
  
  map_trait_2 <- map %>%
    dplyr::group_by(marker)%>%
    dplyr::filter(trait == traits[2])%>%
    dplyr::filter(lod == max(lod))

  cis_trait_2 <- map %>%
    dplyr::group_by(marker) %>%
    dplyr::filter(trait == traits[2])%>%
    dplyr::mutate(maxlod=max(lod))%>%
    dplyr::group_by(iteration) %>%
    dplyr::filter(!is.na(var_exp)) %>%
    dplyr::do(head(., n=1))
  
  map_trait_3 <- map %>%
    dplyr::group_by(marker)%>%
    dplyr::filter(trait == traits[3])%>%
    dplyr::filter(lod == max(lod))
  
  cis_trait_3 <- map %>%
    dplyr::group_by(marker) %>%
    dplyr::filter(trait == traits[3])%>%
    dplyr::mutate(maxlod=max(lod))%>%
    dplyr::group_by(iteration) %>%
    dplyr::filter(!is.na(var_exp)) %>%
    dplyr::do(head(., n=1))
  
  
  if(nrow(cis_trait_1) == 0 | nrow(cis_trait_2) == 0) {
    plot <- ggplot2::ggplot(map1, ggplot2::aes(x = genotype, y = pheno)) + ggplot2::geom_blank()
    return(plot)
  }
  
  map_trait_1 <- cidefiner(cis_trait_1, map_trait_1)
  map_trait_2 <- cidefiner(cis_trait_2, map_trait_2)
  map_trait_3 <- cidefiner(cis_trait_3, map_trait_3)
  
  map1 <- dplyr::bind_rows(map_trait_1,map_trait_2, map_trait_3)
  
  plot <- ggplot2::ggplot(map1) +
    ggplot2::aes(x = pos/1e6, y = lod)
    
  
  if(nrow(cis_trait_1) != 0 | nrow(cis_trait_2) != 0 | nrow(cis_trait_3) != 0) {
    plot <- plot + 
      ggplot2::geom_ribbon(ggplot2::aes(x = pos/1e6, ymin = 0, ymax = ci_lod),
                           fill = "hotpink3", alpha = 0.5, data = map_trait_1)  +
      ggplot2::geom_ribbon(ggplot2::aes(x = pos/1e6, ymin = 0, ymax = ci_lod),
                           fill = "cadetblue3", alpha = 0.5, data = map_trait_2) +
      ggplot2::geom_ribbon(ggplot2::aes(x = pos/1e6, ymin = 0, ymax = ci_lod),
                           fill = "black", alpha = 0.5, data = map_trait_3)
  }
  
  plot <- plot + 
    ggplot2::geom_line(size = 1, alpha = 0.85, color = "hotpink3",data = map_trait_1) +
    ggplot2::geom_line(size = 1, alpha = 0.85, color = "cadetblue3",data = map_trait_2) +
    ggplot2::geom_line(size = 1, alpha = 0.85, color = "black",data = map_trait_3) +
    ggplot2::facet_grid(.~chr, scales ="free", space = "free") +
    ggplot2::labs(x = "Genomic Position (Mb)", y = "LOD") +
    ggplot2::scale_colour_manual(values = c("hotpink3","black","cadetblue3", "black"),
                                 labels = c("Animal Length", "Brood Size"),
                                 name="Phenotype") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 16),
                   axis.text.y = ggplot2::element_text(size = 16),
                   axis.title.x = ggplot2::element_text(size = 20, face = "bold", color = "black", vjust = -0.3), 
                   axis.title.y = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                   strip.text.x = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                   strip.text.y = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                   plot.title = ggplot2::element_text(size = 24, face = "bold", vjust = 1), 
                   panel.background = ggplot2::element_rect(color = "black",size = 1.2))
  return(plot)
}


cidefiner <- function(cis, map) {
  ci_lod <- vapply(1:nrow(map), function(marker) {
    pos <- map$pos[marker]
    chr <- map$chr[marker]
    lod <- map$lod[marker]
    s <- sum(chr == cis$chr & (pos >= cis$ci_l_pos & pos <= cis$ci_r_pos))
    inci <- as.logical(s)
    cilodscore <- ifelse(inci, lod, 0)
    return(cilodscore)
  }, numeric(1))
  map$ci_lod <- ci_lod
  return(map)
}

pxgplot_edit <- function(cross, map, parent="N2xCB4856") {
  peaks <- map %>% 
    dplyr::group_by(iteration) %>%
    dplyr::filter(!is.na(var_exp)) %>%
    dplyr::do(head(., n=1))
  
  if(nrow(peaks) == 0) {
    stop("No QTL identified")
  }
  
  uniquemarkers <- gsub("-", "\\.", unique(peaks$marker))
  
  colnames(cross$pheno) <- gsub("-", "\\.", colnames(cross$pheno))
  
  pheno <- cross$pheno %>%
    dplyr::select_(map$trait[1])
  geno <- data.frame(extract_genotype(cross)) %>%
    dplyr::select(which(colnames(.) %in% uniquemarkers)) %>%
    data.frame(., pheno)
  
  colnames(geno)[1:(ncol(geno)-1)] <- sapply(colnames(geno)[1:(ncol(geno)-1)],
                                             function(marker) {
                                               paste(
                                                 unlist(
                                                   peaks[
                                                     peaks$marker == 
                                                       gsub("\\.",
                                                            "-",
                                                            marker),
                                                     c("chr", "pos")]),
                                                 collapse = ":")
                                             })
  colnames(geno)[ncol(geno)] <- "pheno"
  
  split <- tidyr::gather(geno, marker, genotype, -pheno)
  
  split$genotype <- sapply(split$genotype, function(x){
    if(is.na(x)) {
      return(NA)
    }
    if(parent=="N2xCB4856") {
      if(x == -1) {
        "N2"
      } else {
        "CB4856"
      }
    } else if(parent=="N2xLSJ2") {
      if(x == -1) {
        "N2"
      } else {
        "LSJ2"
      }
    } else if(parent=="AF16xHK104") {
      if(x==-1) {
        "AF16"
      } else {
        "HK104"
      }
    }
  })
  
  split$genotype <- factor(split$genotype, levels = c("N2","CB4856","LSJ2","AF16","HK104"))
  
  ggplot2::ggplot(split) +
    ggbeeswarm::geom_beeswarm(ggplot2::aes(x = genotype, y = pheno), alpha = .4,priority = "density",cex = 1.2) +
    ggplot2::geom_boxplot(ggplot2::aes(x = genotype, y = pheno, fill = genotype), outlier.shape = NA,alpha=0.7) +
    ggplot2::scale_fill_manual(values = c("N2" = "orange", "CB4856" = "blue", "LSJ2" = "green", "AF16" = "indianred", "HK104"= "gold")) +
    ggplot2::facet_wrap(~ marker, ncol = 5) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 16), 
                   axis.text.y = ggplot2::element_text(size = 16), 
                   axis.title.x = ggplot2::element_text(size = 0,  face = "bold", color = "black", vjust = -0.3), 
                   axis.title.y = ggplot2::element_text(size = 16, face = "bold", color = "black"), 
                   strip.text.x = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                   strip.text.y = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                   plot.title = ggplot2::element_text(size = 0,  face = "bold", vjust = 1), 
                   panel.background = ggplot2::element_rect(color = "black",  size = 1.2),
                   legend.position = "none") +
    ggplot2::ggtitle(peaks$trait[1]) +
    ggplot2::labs(x = "Genotype", y = "Phenotype")
}

tajimas_d_temp <- function (vcf_path = paste0(path.package("cegwas"), "/"), 
                            vcf_name = paste0("WI.", vcf_version, ".impute.vcf.gz"), 
                            chromosome = "II", interval_start = 11021073, 
                            interval_end = 12008179, window_size = 300, slide_distance = 100, 
                            samples = colnames(snps[, 5:ncol(snps)]), outgroup = "N2", 
                            site_of_interest = 11875145) 
{
  setwd(vcf_path)
  gen <- PopGenome::readVCF(vcf_name, numcols = 10000, tid = chromosome, 
                            frompos = interval_start, topos = interval_end, samplenames = samples, 
                            approx = F)
  gen1 <- PopGenome::set.populations(gen, list(samples), diploid = FALSE)
  gen2 <- PopGenome::set.outgroup(gen1, outgroup, diploid = FALSE)
  s_gen <- PopGenome::sliding.window.transform(gen2, width = window_size, 
                                               jump = slide_distance, 
                                               whole.data = FALSE, type = 2)
  test <- data.frame(snps = 1:length(as.numeric(colnames(as.data.frame(s_gen@BIG.BIAL[[1]][,])))), 
                     position = as.numeric(colnames(as.data.frame(s_gen@BIG.BIAL[[1]][, ]))))
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
    dplyr::rename(Td = pop.1) %>% 
    dplyr::distinct(Td, window,  .keep_all = T)
  
  tajimas_d_plot <- ggplot2::ggplot(td) + 
    ggplot2::aes(x = position/1e+06, y = Td) + 
    ggplot2::geom_point(size = 1, alpha = 0.75) + 
    ggplot2::theme_bw() + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 16),
                   axis.text.y = ggplot2::element_text(size = 16),
                   axis.title.x = ggplot2::element_text(size = 20, face = "bold", color = "black", vjust = -0.3), 
                   axis.title.y = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                   strip.text.x = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                   strip.text.y = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                   plot.title = ggplot2::element_text(size = 24,face = "bold", vjust = 1), 
                   legend.position = "none", 
                   panel.background = ggplot2::element_rect(color = "black", size = 1.2), 
                   strip.background = ggplot2::element_rect(color = "black", size = 1.2)) + 
    ggplot2::labs(x = "Genomic Position (Mb)",  y = "Tajima's D")
  
  if (!is.null(site_of_interest)) {
    tajimas_d_plot <- tajimas_d_plot + 
      ggplot2::geom_vline(ggplot2::aes(xintercept = site_of_interest/1e+06), color = "red", alpha = 0.7, size = 1)
  }
  return(list(td, tajimas_d_plot))
}

rxn_norm_plt <- function(df, trt, cond, fancy_name, 
                         strains, fancy_strains, 
                         ordered_conditions, fancy_ordered_conditions){
  
  
  plot_concentrations <- df%>%
    dplyr::filter(trait == trt, condition %in% cond, grepl("C15",condition))
  
  plot_concentrations <- unique(plot_concentrations$conc)
  
  df%>%
    dplyr::filter(trait == trt, condition %in% cond, conc %in% plot_concentrations)%>%
    dplyr::mutate(strain1 = factor(strain, levels = strains,
                                   labels = fancy_strains),
                  condition1 = factor(condition, 
                                      levels = ordered_conditions,
                                      labels = fancy_ordered_conditions))%>%
    dplyr::mutate(strain_cond = paste(strain, condition, sep = "\n"))%>%
    dplyr::group_by(strain_cond)%>%
    dplyr::mutate(m_ph = mean(phenotype),
                  sd_ph = sd(phenotype))%>%
    ggplot(.) +
    aes(x = factor(condition1), 
        y = m_ph, 
        color = strain1,
        group = strain1,
        alpha = strain1) +
    geom_line(aes(linetype = strain1))+
    geom_point()+
    geom_errorbar(aes(ymin=m_ph-sd_ph, ymax=m_ph+sd_ph), width=.1)+
    scale_color_manual(values = c("Bristol" = "orange","Hawaii" = "blue",
                                  "Bristol\n(C78S)" = "orange","Bristol\n(C78S)" = "orange",
                                  "Hawaii\n(S78C)" = "blue","Hawaii\n(S78C)" = "blue"),
                       name="Strain")+
    scale_linetype_manual(values = c("Bristol" = 1,"Hawaii" = 1,
                                     "Bristol\n(C78S)" = 4,"Bristol\n(C78S)" = 4,
                                     "Hawaii\n(S78C)" = 4,"Hawaii\n(S78C)" = 4),
                          name="Strain")+
    scale_alpha_manual(values = c("Bristol" = 1,"Hawaii" = 1,
                                  "Bristol\n(C78S)" = .5,"Bristol\n(C78S)" = .5,
                                  "Hawaii\n(S78C)" = .5,"Hawaii\n(S78C)" = .5),
                       name="Strain")+
    theme_bw()+
    facet_grid(.~conc, scales = "free")+
    labs( y = trt)+
    theme(axis.text.x = element_text(size=10, face="bold", color="black"),
          axis.text.y = element_text(size=12, face="bold", color="black"),
          axis.title.x = element_text(size=0, face="bold", color="black", vjust=-.3),
          axis.title.y = element_text(size=16, face="bold", color="black"),
          strip.text.x = element_text(size=24, face="bold", color="black"),
          strip.text.y = element_text(size=16, face="bold", color="black"),
          plot.title = element_text(size=16, face="bold", vjust = 1),
          # legend.position="none",
          panel.background = element_rect( color="black",size=1.2),
          strip.background = element_rect(color = "black", size = 1.2),
          panel.border = element_rect( colour = "black"))+
    labs(y = fancy_name)
}

boxplot_plt <- function(df, 
                        trt, 
                        cond,
                        fancy_name, 
                        strains,
                        fancy_strains, 
                        ordered_conditions, 
                        fancy_ordered_conditions,
                        r_conc){
  
  plot_concentrations <- df%>%
    dplyr::filter(Trait == trt, 
                  Condition %in% cond, 
                  grepl("C15", Condition))
  
  plot_concentrations <- unique(plot_concentrations$conc)
  
  df%>%
    dplyr::filter(Trait == trt, 
                  Condition %in% cond, 
                  conc %in% plot_concentrations)%>%
    dplyr::mutate(strain1 = factor(Strain, 
                                   levels = strains,
                                   labels = fancy_strains),
                  condition1 = factor(Condition, 
                                      levels = ordered_conditions,
                                      labels = fancy_ordered_conditions))%>%
    dplyr::mutate(strain_cond = paste(Strain, Condition, sep = "\n"))%>%
    dplyr::mutate(strain_cond_reorder = factor(strain_cond, 
                                               levels = c(paste0("N2\nC15ISO",r_conc),
                                                          "N2\nArsenic",
                                                          paste0("N2\nArsenicC15ISO",r_conc),
                                                          paste0("ECA581\nC15ISO",r_conc),
                                                          "ECA581\nArsenic",
                                                          paste0("ECA581\nArsenicC15ISO",r_conc),
                                                          paste0("CB4856\nC15ISO",r_conc),
                                                          "CB4856\nArsenic",
                                                          paste0("CB4856\nArsenicC15ISO",r_conc),
                                                          paste0("ECA590\nC15ISO",r_conc),
                                                          "ECA590\nArsenic",
                                                          paste0("ECA590\nArsenicC15ISO",r_conc)),
                                               labels = c(paste("N2\nC15ISO"),
                                                          "N2\nArsenic",
                                                          paste0("N2\nArsenic\nC15ISO"),
                                                          paste0("N2\nDBT-1(C78S)\nC15ISO"),
                                                          "N2\nDBT-1(C78S)\nArsenic",
                                                          paste0("N2\nDBT-1(C78S)\nArsenic\nC15ISO"),
                                                          paste0("CB4856\nC15ISO"),
                                                          "CB4856\nArsenic",
                                                          paste0("CB4856\nArsenic\nC15ISO"),
                                                          paste0("CB4856\nDBT-1(S78C)\nC15ISO"),
                                                          "CB4856\nDBT-1(S78C)\nArsenic",
                                                          paste0("CB4856\nDBT-1(S78C)\nArsenic\nC15ISO"))))%>%
    ggplot(.) +
    aes(x = factor(strain_cond_reorder), 
        y = Value,
        fill = strain1) +
    geom_beeswarm(priority = "density", cex = 0.6, size = 2, alpha = 0.75) +
    geom_boxplot(outlier.colour = NA, alpha = .7)+
    scale_fill_manual(values = c("Bristol" = "orange","Hawaii" = "blue",
                                 "Bristol\n(C78S)" = "gray","Bristol\n(C78S)" = "gray",
                                 "Hawaii\n(S78C)" = "gray","Hawaii\n(S78C)" = "gray"),
                      name="Strain")+
    theme_bw()+
    # facet_grid(.~conc, scales = "free")+
    labs( y = trt)+
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 16), 
                   axis.text.y = ggplot2::element_text(size = 16), 
                   axis.title.x = ggplot2::element_text(size = 0,  face = "bold", color = "black", vjust = -0.3), 
                   axis.title.y = ggplot2::element_text(size = 16, face = "bold", color = "black"), 
                   strip.text.x = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                   strip.text.y = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                   plot.title = ggplot2::element_text(size = 0,  face = "bold", vjust = 1), 
                   panel.background = ggplot2::element_rect(color = "black",  size = 1.2),
                   legend.position = "none")+
    labs(y = fancy_name)
}


run_TukeyHSD <- function(df, trt){
  stat_df <- df %>%
    dplyr::filter(trait == trt)%>%
    dplyr::select(strain, phenotype)
  
  aov_res <- aov(stat_df$phenotype ~ stat_df$strain)
  summary(aov_res)
  tuk <- TukeyHSD(aov_res)
  
  psig=as.numeric(apply(tuk$`stat_df$strain`[,2:3],1,prod)>=0)+1
  op=par(mar=c(4.2,9,3.8,2))
  plot(tuk,col=psig,yaxt="n")
  for (j in 1:length(psig)){
    axis(2,at=j,labels=rownames(tuk$`stat_df$strain`)[length(psig)-j+1],
         las=1,cex.axis=.8,col.axis=psig[length(psig)-j+1])
  }
  par(op)
  
  pwtuk <- TukeyHSD(aov_res)
  
  return(pwtuk)
}

run_pca <- function(data, format = "w", drop_method = "t", n_pc = 10){
  
  if( format == "l" ){
    print("Data is in long format")
    
    trait_col <- quo(colnames(data)[2])
    value_col <- quo(colnames(data)[3])
    
    n_traits <- data %>%
      dplyr::pull(!!trait_col) %>%
      unique() %>%
      length()
    
    traits <- data %>%
      dplyr::pull(!!trait_col) %>%
      unique() 
    
    wd <- data %>%
      tidyr::spread(!!trait_col, !!value_col)
    
  } else {
    print("Data is in wide format")
    
    wd <- data
    
    n_traits <- ncol(data) - 1
    
    traits <- colnames(data)[2:ncol(data)]
  }
  
  # deal with NAs
  n_obs <- nrow(wd)
  
  if( n_obs != nrow(na.omit(wd)) ) {
    
    if(drop_method == "t"){
      
      na_cols <- colnames(wd)[colSums(is.na(wd)) > 0]
      
      wd <- wd %>%
        dplyr::select(-na_cols)
      
      print(glue::glue("Traits removed because they have missing data: {na_cols}"))
      
      cols_removed <- length(na_cols)
      
    } else if ( drop_method == "s" ){
      
      dropped_strains <- row.names(wd)[rowSums(is.na(wd)) > 0]
      
      wd <- na.omit(wd) 
      
      print(glue::glue("Strains removed because they have missing data: {dropped_strains}"))
      
      cols_removed <- 0
    }
  }
  
  # 
  # # if we wanted to deal with meta data
  # if( (ncol(wd)-1) == n_traits - cols_removed ) {
  #   
  #   print("No meta data")
  #   
  meta_cols <- 1
  #   
  # } else {
  #   
  #   n_meta_cols <- (ncol(wd)-1) - (n_traits - cols_removed)
  #   
  #   meta_cols <- colnames(wd)[!colnames(wd) %in% colnames(data)[1] | colnames(wd) %in% ]
  #   
  # }
  
  to_pca <- data.frame(scale(wd[,(meta_cols+1):ncol(wd)]))
  
  pca_object <- prcomp(to_pca)
  
  if(is.numeric(n_pc)){
    
    n_pc <- n_pc 
    
  } else {
    
    t_pc_ve <- sum((pca_object$sdev)^2)
    
    pc_ve <- cumsum((pca_object$sdev)^2)/t_pc_ve
    
    n_pc <- which(pc_ve > 0.9)[1]
    
  }
  
  PC_df <- data.frame(Strain = wd[,1],
                      pca_object$x[,1:n_pc]) %>%
    tidyr::gather(Trait, Value, -Strain)
  
  # extract PC loadings
  loading_df <- data.frame(Trait = row.names(pca_object$rotation),
                           pca_object$rotation[,1:n_pc]) 
  
  loading_df <- loading_df %>%
    tidyr::gather(PC, Loading, -Trait)
  
  vars <- apply(pca_object$x, 2, var)  
  props <- vars / sum(vars)
  pca_var_df <- data.frame(PC = as.numeric(gsub("PC","",names(props))), Variance = cumsum(props))
  
  list(PC_df, loading_df, pca_var_df)
  
}

cegwas2_manplot <- function(plot_df, 
                            bf_line_color = "gray",
                            eigen_line_color = "gray",
                            eigen_cutoff = independent_test_cutoff) {
  plot_traits <- unique(plot_df$trait)
  plots <- lapply(plot_traits, function(i) {
    
    bf_cut <- -log10(0.05/length(unique(plot_df$marker)))
    
    plot_df_pr <- plot_df %>%
      dplyr::filter(trait == i,
                    CHROM != "MtDNA") %>%
      dplyr::distinct(marker, .keep_all = T) %>%
      dplyr::mutate(EIGEN_CUTOFF = eigen_cutoff,
                    BF = bf_cut) %>%
      dplyr::mutate(EIGEN_SIG = ifelse(log10p > bf_cut, "1", 
                                       ifelse(log10p > EIGEN_CUTOFF, "2", "0")) )
    
    plot_df_pr  %>%
      ggplot2::ggplot(.) +
      ggplot2::aes(x = POS/1e6, y = log10p) +
      ggplot2::scale_color_manual(values = c("0" = "black", 
                                             "1" = "red",
                                             "2" = "hotpink3")) +
      ggplot2::geom_rect(ggplot2::aes(xmin = startPOS/1e6, 
                                      xmax = endPOS/1e6, 
                                      ymin = 0, 
                                      ymax = Inf, 
                                      fill = "hotpink3"), 
                         color = "hotpink",linetype = 2, 
                         alpha=.3, data = dplyr::filter(plot_df_pr, EIGEN_SIG!="1") %>% na.omit())+
      ggplot2::geom_rect(ggplot2::aes(xmin = startPOS/1e6, 
                                      xmax = endPOS/1e6, 
                                      ymin = 0, 
                                      ymax = Inf, 
                                      fill = "blue"), 
                         color = "blue",fill = "cyan",linetype = 2, 
                         alpha=.3, data = dplyr::filter(plot_df_pr, EIGEN_SIG=="1") %>% na.omit())+
      ggplot2::geom_hline(ggplot2::aes(yintercept = bf_cut),
                          color = bf_line_color, 
                          alpha = .75,  
                          size = 1) +
      ggplot2::geom_hline(ggplot2::aes(yintercept = EIGEN_CUTOFF),
                          color = eigen_line_color, 
                          alpha = .75,  
                          size = 1,
                          linetype = 2) +
      ggplot2::geom_point( ggplot2::aes(color= factor(EIGEN_SIG)) ) +
      ggplot2::facet_grid( . ~ CHROM, scales = "free_x" , space = "free_x") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(size = 16),
                     axis.text.y = ggplot2::element_text(size = 16),
                     axis.title.x = ggplot2::element_text(size = 20, face = "bold", color = "black", vjust = -0.3), 
                     axis.title.y = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                     strip.text.x = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                     strip.text.y = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                     plot.title = ggplot2::element_text(size = 24, face = "bold", vjust = 1), 
                     panel.background = ggplot2::element_rect(color = "black",size = 1.2),
                     legend.position = "none",
                     strip.background = element_blank()) +
      ggplot2::labs(x = "Genomic Position (Mb)",
                    y = expression(-log[10](italic(p))),
                    title = plot_traits)
  })
  plots
}




pxgplot_edit_1peak <- function(cross, map, parent="N2xCB4856") {
  peaks <- map %>% 
    dplyr::group_by(iteration) %>%
    dplyr::filter(!is.na(var_exp)) %>%
    dplyr::do(head(., n=1))
  
  if(nrow(peaks) == 0) {
    stop("No QTL identified")
  }
  
  uniquemarkers <- gsub("-", "\\.", unique(peaks$marker))
  
  colnames(cross$pheno) <- gsub("-", "\\.", colnames(cross$pheno))
  
  pheno <- cross$pheno %>%
    dplyr::select_(map$trait[1])
  geno <- data.frame(extract_genotype(cross)) %>%
    dplyr::select(which(colnames(.) %in% uniquemarkers)) %>%
    data.frame(., pheno)
  
  colnames(geno)[1:(ncol(geno)-1)] <- sapply(colnames(geno)[1:(ncol(geno)-1)],
                                             function(marker) {
                                               paste(
                                                 unlist(
                                                   peaks[
                                                     peaks$marker == 
                                                       gsub("\\.",
                                                            "-",
                                                            marker),
                                                     c("chr", "pos")]),
                                                 collapse = ":")
                                             })
  colnames(geno)[ncol(geno)] <- "pheno"
  
  split <- tidyr::gather(geno, marker, genotype, -pheno)
  
  split$genotype <- sapply(split$genotype, function(x){
    if(is.na(x)) {
      return(NA)
    }
    if(parent=="N2xCB4856") {
      if(x == -1) {
        "N2"
      } else {
        "CB4856"
      }
    } else if(parent=="N2xLSJ2") {
      if(x == -1) {
        "N2"
      } else {
        "LSJ2"
      }
    } else if(parent=="AF16xHK104") {
      if(x==-1) {
        "AF16"
      } else {
        "HK104"
      }
    }
  })
  
  split$genotype <- factor(split$genotype, levels = c("N2","CB4856","LSJ2","AF16","HK104"))
  
  ggplot2::ggplot(split) +
    ggbeeswarm::geom_beeswarm(ggplot2::aes(x = genotype, y = pheno), alpha = .4,priority = "density",cex = 1.2) +
    ggplot2::geom_boxplot(ggplot2::aes(x = genotype, y = pheno, fill = genotype), outlier.shape = NA,alpha=0.7) +
    ggplot2::scale_fill_manual(values = c("N2" = "orange", "CB4856" = "blue", "LSJ2" = "green", "AF16" = "indianred", "HK104"= "gold")) +
    # ggplot2::facet_wrap(~ marker, ncol = 5) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 16), 
                   axis.text.y = ggplot2::element_text(size = 16), 
                   axis.title.x = ggplot2::element_text(size = 0,  face = "bold", color = "black", vjust = -0.3), 
                   axis.title.y = ggplot2::element_text(size = 16, face = "bold", color = "black"), 
                   strip.text.x = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                   strip.text.y = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                   plot.title = ggplot2::element_text(size = 0,  face = "bold", vjust = 1), 
                   # panel.background = ggplot2::element_rect(color = "black",  size = 1.2),
                   legend.position = "none") +
    ggplot2::ggtitle(peaks$trait[1]) +
    ggplot2::labs(x = "Genotype", y = "Phenotype")
}


maxlodplot_edit3 <- function(map){
  
  traits <- unique(map$trait)
  
  map_trait_1 <- map %>%
    dplyr::group_by(marker)%>%
    dplyr::filter(trait == traits[1])%>%
    dplyr::filter(lod == max(lod))
  
  cis_trait_1 <- map %>%
    dplyr::group_by(marker) %>%
    dplyr::filter(trait == traits[1])%>%
    dplyr::mutate(maxlod=max(lod))%>%
    dplyr::group_by(iteration) %>%
    dplyr::filter(!is.na(var_exp)) %>%
    dplyr::do(head(., n=1))
  
  map_trait_2 <- map %>%
    dplyr::group_by(marker)%>%
    dplyr::filter(trait == traits[2])%>%
    dplyr::filter(lod == max(lod))
  
  cis_trait_2 <- map %>%
    dplyr::group_by(marker) %>%
    dplyr::filter(trait == traits[2])%>%
    dplyr::mutate(maxlod=max(lod))%>%
    dplyr::group_by(iteration) %>%
    dplyr::filter(!is.na(var_exp)) %>%
    dplyr::do(head(., n=1))
  
  map_trait_3 <- map %>%
    dplyr::group_by(marker)%>%
    dplyr::filter(trait == traits[3])%>%
    dplyr::filter(lod == max(lod))
  
  cis_trait_3 <- map %>%
    dplyr::group_by(marker) %>%
    dplyr::filter(trait == traits[3])%>%
    dplyr::mutate(maxlod=max(lod))%>%
    dplyr::group_by(iteration) %>%
    dplyr::filter(!is.na(var_exp)) %>%
    dplyr::do(head(., n=1))
  
  map_trait_4 <- map %>%
    dplyr::group_by(marker)%>%
    dplyr::filter(trait == traits[4])%>%
    dplyr::filter(lod == max(lod))
  
  cis_trait_4 <- map %>%
    dplyr::group_by(marker) %>%
    dplyr::filter(trait == traits[4])%>%
    dplyr::mutate(maxlod=max(lod))%>%
    dplyr::group_by(iteration) %>%
    dplyr::filter(!is.na(var_exp)) %>%
    dplyr::do(head(., n=1))
  
  # map_trait_5 <- map %>%
  #   dplyr::group_by(marker)%>%
  #   dplyr::filter(trait == traits[5])%>%
  #   dplyr::filter(lod == max(lod))
  # 
  # cis_trait_5 <- map %>%
  #   dplyr::group_by(marker) %>%
  #   dplyr::filter(trait == traits[5])%>%
  #   dplyr::mutate(maxlod=max(lod))%>%
  #   dplyr::group_by(iteration) %>%
  #   dplyr::filter(!is.na(var_exp)) %>%
  #   dplyr::do(head(., n=1))
  
  
  if(nrow(cis_trait_1) == 0 | nrow(cis_trait_2) == 0) {
    plot <- ggplot2::ggplot(map1, ggplot2::aes(x = genotype, y = pheno)) + ggplot2::geom_blank()
    return(plot)
  }
  
  map_trait_1 <- cidefiner(cis_trait_1, map_trait_1)
  map_trait_2 <- cidefiner(cis_trait_2, map_trait_2)
  map_trait_3 <- cidefiner(cis_trait_3, map_trait_3)
  map_trait_4 <- cidefiner(cis_trait_4, map_trait_4)
  # map_trait_5 <- cidefiner(cis_trait_5, map_trait_5)
  
  map1 <- dplyr::bind_rows(map_trait_1,map_trait_2, map_trait_3,map_trait_4)
  
  plot <- ggplot2::ggplot(map1) +
    ggplot2::aes(x = pos/1e6, y = lod)
  
  
  if(nrow(cis_trait_1) != 0 | nrow(cis_trait_2) != 0 | nrow(cis_trait_3) != 0) {
    plot <- plot + 
      ggplot2::geom_ribbon(ggplot2::aes(x = pos/1e6, ymin = 0, ymax = ci_lod),
                           fill = "hotpink3", alpha = 0.5, data = map_trait_2)  +
      ggplot2::geom_ribbon(ggplot2::aes(x = pos/1e6, ymin = 0, ymax = ci_lod),
                           fill = "cadetblue3", alpha = 0.5, data = map_trait_3) +
      ggplot2::geom_ribbon(ggplot2::aes(x = pos/1e6, ymin = 0, ymax = ci_lod),
                           fill = "red", alpha = 0.5, data = map_trait_1)
      ggplot2::geom_ribbon(ggplot2::aes(x = pos/1e6, ymin = 0, ymax = ci_lod),
                         fill = "orange", alpha = 0.5, data = map_trait_4)
      # ggplot2::geom_ribbon(ggplot2::aes(x = pos/1e6, ymin = 0, ymax = ci_lod),
      #                      fill = "black", alpha = 0.5, data = map_trait_5)
  }
  
  plot <- plot + 
    ggplot2::geom_line(size = 1, alpha = 0.85, color = "hotpink3",data = map_trait_2) + #size
    ggplot2::geom_line(size = 1, alpha = 0.85, color = "cadetblue3",data = map_trait_3) +#brood
    ggplot2::geom_line(size = 1, alpha = 0.85, color = "red",data = map_trait_1) +#ext
    ggplot2::geom_line(size = 1, alpha = 0.85, color = "orange",data = map_trait_4) +#fluor
    # ggplot2::geom_line(size = 1, alpha = 0.85, color = "black",data = map_trait_5) +
    ggplot2::facet_grid(.~chr, scales ="free", space = "free") +
    ggplot2::labs(x = "Genomic Position (Mb)", y = "LOD") +
    # ggplot2::scale_colour_manual(values = c("hotpink3","black","cadetblue3", "red"),
    #                              labels = c("Animal Length", "Brood Size", "Optical Density"),
    #                              name="Phenotype") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 16),
                   axis.text.y = ggplot2::element_text(size = 16),
                   axis.title.x = ggplot2::element_text(size = 20, face = "bold", color = "black", vjust = -0.3), 
                   axis.title.y = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                   strip.text.x = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                   strip.text.y = ggplot2::element_text(size = 20, face = "bold", color = "black"), 
                   plot.title = ggplot2::element_text(size = 24, face = "bold", vjust = 1), 
                   panel.background = ggplot2::element_rect(color = "black",size = 1.2))
  return(plot)
}

