library(easysorter)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cegwas)
library(PopGenome)
library(directlabels)
library(broom)
library(grDevices)
library(cowplot)
library(ggtree)
library(ggalt)    # devtools::install_github("hrbrmstr/ggalt")
library(ggthemes) # theme_map and tableau colors
library(lemon)


ancestry.colours <- c("gold2", "plum4", "darkorange1", 
                      "lightskyblue2", "firebrick", "burlywood3","gray51", 
                      "springgreen4", "lightpink2","deepskyblue4", "black", 
                      "mediumpurple4", "orange", "maroon","yellow3", "brown4", 
                      "yellow4", "sienna4", "chocolate", "gray19")

plot.theme <-  theme(axis.text.x = ggplot2::element_text(size = 8),
                     axis.text.y = ggplot2::element_text(size = 8),
                     axis.title.x = ggplot2::element_text(size = 10, face = "bold", color = "black", vjust = -0.3),
                     axis.title.y = ggplot2::element_text(size = 10, face = "bold", color = "black", vjust = -0.3),
                     strip.text.x = element_text(size = 10, face = "bold"))

read_file <- function(file, tofmin = 60, tofmax = 2000, extmin = 0,
                      extmax = 10000, SVM = TRUE, levels = 2){
  
  # Read the raw sorter files and make the row names
  
  plate <- readSorter(file, tofmin, tofmax, extmin, extmax)
  modplate <- with(plate,
                   data.frame(row=Row,
                              col=as.factor(Column),
                              sort=Status.sort,
                              TOF=TOF,
                              EXT=EXT,
                              time=Time.Stamp,
                              green=Green,
                              yellow=Yellow,
                              red=Red))
  
  # Extract the time so that it is realtive to the first worm sorted
  
  modplate <- modplate %>%
    dplyr::group_by(row, col) %>%
    dplyr::do(COPASutils::extractTime(.))
  modplate <- data.frame(modplate)
  
  # Normalize the optical values by time of flight
  
  modplate[, 10:13] <- apply(modplate[, c(5, 7:9)], 2,
                             function(x) x / modplate$TOF)
  colnames(modplate)[10:13] <- c("norm.EXT", "norm.green", "norm.yellow",
                                 "norm.red")
  
  # Handle the SVM predictions if requested
  
  if(SVM){
    plateprediction <- kernlab::predict(
      COPASutils::bubbleSVMmodel_noProfiler,
      modplate[,3:length(modplate)],
      type="probabilities")
    modplate$object <- plateprediction[, "1"]
    modplate$call50 <- factor(as.numeric(modplate$object > 0.5),
                              levels=c(1, 0), labels=c("object", "bubble"))
  }
  
  # Calculate the life stage values based on the size of the worms
  
  modplate$stage <- ifelse(modplate$TOF >= 60 & modplate$TOF < 90, "L1",
                           ifelse(modplate$TOF >= 90 & modplate$TOF < 200,
                                  "L2/L3",
                                  ifelse(modplate$TOF >= 200
                                         & modplate$TOF < 300, "L4",
                                         ifelse(modplate$TOF >= 300,
                                                "adult", NA))))
  
  # Convert integer values to numerics
  
  modplate[, as.vector(which(lapply(modplate, class) == "integer"))] <- lapply(
    modplate[, as.vector(which(lapply(modplate, class) == "integer"))],
    as.numeric)
  
  # Get info about the plate using the new_info function 
  
  plateinfo <- new_info(file, levels)
  
  # Get the template base directory
  
  templatedir <- strsplit(file, "/")[[1]]
  templatedir <- templatedir[-c(length(templatedir), length(templatedir) - 1)]
  templatedir <- paste0(templatedir, collapse = "/")
  templatedir <- paste0(templatedir, "/")
  
  # Get the template file paths
  
  strainsfile <- paste0(templatedir, "strains/",
                        plateinfo$straintemplate[1], ".csv")
  conditionsfile <- paste0(templatedir, "conditions/",
                           plateinfo$conditiontemplate[1],
                           ".csv")
  controlsfile <- paste0(templatedir, "controls/",
                         plateinfo$controltemplate[1], ".csv")
  contamfile <- paste0(templatedir,
                       "contamination/",
                       sprintf("p%02d", plateinfo$plate[1]),
                       "_contamination.csv")
  
  # Read all of the templates
  
  strains <- read_template(strainsfile, type="strains")
  conditions  <- read_template(conditionsfile, type="conditions")
  controls <- read_template(controlsfile, type="controls")
  contam <- read_template(contamfile, type="contam")
  
  # Join all of the metadata and template info to the plate data
  
  modplate <- cbind(plateinfo[,1:5], modplate)
  modplate <- dplyr::left_join(modplate, strains, by = c("row", "col"))
  modplate <- dplyr::left_join(modplate, conditions, by = c("row", "col"))
  modplate <- dplyr::left_join(modplate, controls, by = c("row", "col"))
  modplate <- dplyr::left_join(modplate, contam, by = c("row", "col"))
  
  return(modplate)
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
                         alpha=.3, size = 0.3)+
      ggplot2::geom_hline(ggplot2::aes(yintercept = BF),
                          color = bf_line_color, 
                          alpha = .75,  
                          size = 0.5) +
      ggplot2::geom_point( ggplot2::aes(color= factor(aboveBF)), size = 0.25 ) +
      ggplot2::facet_grid( . ~ CHROM, scales = "free_x" , space = "free_x") +
      ggplot2::theme_bw() +
      plot.theme +
      theme(strip.background = element_blank(),
            legend.position = "none") +
      ggplot2::labs(x = "Genomic Position (Mb)",
                    y = expression(bold(-log[10](bolditalic(p)))))
  })
  plots
}

plot_bar <- function(df){
  ggplot(df)+
    aes(x=strain2, y=final_pheno, fill = factor(build_GT,
                                                levels = sort(unique(build_GT)), 
                                                labels= c("None", sort(unique(build_GT)[2:8]) )))+
    geom_bar(stat="identity", color = "black", size = .2) +
    scale_fill_manual(values=c("None" = "gray75","Deletion" = colors[1],
                               "Insertion" = colors[2],
                               "Inversion" = colors[3],
                               "Missense" = colors[4],
                               "Splice Donor" = colors[5],
                               "Stop Gained" = colors[6],
                               "Transposon Insertion" = colors[7]),name = expression(paste("Variation at ", italic("ben-1"))))+
    theme_bw()+
    labs(x = "Strain", y = paste0("Albendazole Resistance"))
  
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
                                      y=(1.2*lod),
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
    ggplot2::geom_jitter(ggplot2::aes(x = genotype, y = pheno), alpha = .4, width = .25) +
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
                                               whole.data = FALSE)
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
    ggplot2::geom_line(size = 1) + 
    ggplot2::theme_bw() + 
    plot.theme +
    ggplot2::theme(legend.position = "none", 
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

boxplot_plt <- function(df, trt, cond, fancy_name, 
                        strains, fancy_strains, 
                        ordered_conditions, fancy_ordered_conditions,
                        r_conc){
  
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
        y = phenotype,
        fill = strain1) +
    geom_jitter(width = .25,alpha =.4)+
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


gene_model <- function(df, 
                       gene_color = "blue",
                       intron_color = "black",
                       utr3_color = "red",
                       utr5_color = "gray60",
                       gene_alpha = 0.5){
  
  strand <- unique(df$strand)
  
  if(strand == "Watson"){
    
    exons <- df%>%
      dplyr::filter(feature == "Exon")
    
    introns <- df%>%
      dplyr::filter(feature == "Intron")%>%
      dplyr::mutate(halfpt = (start+stop )/2)
    
    startUTR <- df %>%
      dplyr::filter(feature == "5'UTR")
    
    endUTR <- df %>%
      dplyr::filter(feature == "3'UTR")%>%
      dplyr::select(start, stop)%>%
      tidyr::gather(loc, x)%>%
      dplyr::mutate(y = ifelse(loc == "start", 1, 0))
    
    endUTR_a <- dplyr::filter(endUTR, loc == "start")
    endUTR_a$y <- -1
    
    endUTR_pt <- bind_rows(endUTR, endUTR_a)
    
    ggplot(exons)+
      geom_rect( aes(xmin =  stop, xmax = start, ymin = -1 , ymax = 1), fill = gene_color, color = "black",alpha = gene_alpha)+
      geom_segment(aes(x = stop, y = 1, xend = halfpt, yend = 2), data = introns, color = intron_color)+
      geom_segment(aes(x = start, y = 1, xend = halfpt, yend = 2), data = introns, color = intron_color)+
      geom_rect(aes(xmin =  stop, xmax = start, ymin = -1 , ymax = 1), data = startUTR, fill = utr5_color, color = "black",alpha = gene_alpha)+
      geom_polygon(aes(x = x, y = y), fill = utr3_color,color = "black", data = endUTR_pt,alpha = gene_alpha)+
      theme_void()
    
  } else {
    
    print("CRICK STRAND STILL NEEDS TO BE SORTED OUT")
    # exons <- df%>%
    #   dplyr::filter(feature == "Exon")
    # 
    # introns <- df%>%
    #   dplyr::filter(feature == "Intron")%>%
    #   dplyr::mutate(halfpt = (start+stop )/2)
    # 
    # startUTR <- df %>%
    #   dplyr::filter(feature == "5'UTR")
    # 
    # endUTR <- df %>%
    #   dplyr::filter(feature == "3'UTR")%>%
    #   dplyr::select(start, stop)%>%
    #   tidyr::gather(loc, x)%>%
    #   dplyr::mutate(y = ifelse(loc == "stop", 1, 0))
    # 
    # endUTR_a <- dplyr::filter(endUTR, loc == "start")
    # endUTR_a$y <- -1
    # 
    # endUTR_pt <- bind_rows(endUTR, endUTR_a)
    # 
    # ggplot(exons)+
    #   geom_rect( aes(xmin =  start, xmax = stop, ymin = -1 , ymax = 1), fill = "blue", color = "black")+
    #   geom_segment(aes(x = stop, y = 1, xend = halfpt, yend = 2), data = introns, color = "black")+
    #   geom_segment(aes(x = start, y = 1, xend = halfpt, yend = 2), data = introns, color = "black")+
    #   geom_rect(aes(xmin =  stop, xmax = start, ymin = -1 , ymax = 1), data = startUTR, fill = "gray60", color = "black")+
    #   geom_polygon(aes(x = x, y = y), fill = "red",color = "black", data = endUTR_pt)+
    #   theme_void()
    
  }
}

modified_ld_plot <- function (plot_df, trait = NULL) 
{
  if (is.null(trait)) {
    snp_df <- plot_df %>% na.omit()
  }
  else {
    snp_df <- dplyr::filter(plot_df, trait == trait) %>% 
      na.omit()
  }
  ld_snps <- dplyr::filter(snps, CHROM %in% snp_df$CHROM, POS %in% 
                             snp_df$POS)
  ld_snps <- data.frame(snp_id = paste(ld_snps$CHROM, ld_snps$POS, 
                                       sep = "_"), data.frame(ld_snps)[, 5:ncol(ld_snps)])
  sn <- list()
  for (i in 1:nrow(ld_snps)) {
    sn[[i]] <- genetics::genotype(as.character(gsub(1, "T/T", 
                                                    gsub(-1, "A/A", ld_snps[i, 4:ncol(ld_snps)]))))
  }
  test <- data.frame(sn)
  colnames(test) <- (ld_snps$snp_id)
  if (ncol(test) == 1) {
    print("Only one significant SNP, not calculating LD")
  }
  else {
    ldcalc <- t(genetics::LD(test)[[4]])
    diag(ldcalc) <- 1
    LDs <- tbl_df(data.frame(ldcalc) %>% dplyr::add_rownames(var = "SNP1")) %>% 
      tidyr::gather(SNP2, corr, -SNP1) %>% dplyr::arrange(SNP1) %>% 
      tidyr::separate(SNP1, sep = "_", into = c("CHROM1", 
                                                "POS1"), remove = F) %>% dplyr::arrange(CHROM1, 
                                                                                        as.numeric(POS1))
    ldplot <- ggplot2::ggplot(LDs) + ggplot2::aes(x = factor(SNP1, levels = unique(SNP1), ordered = T), 
                                                  y = factor(SNP2, levels = unique(SNP1), ordered = T)) + 
      ggplot2::geom_tile(ggplot2::aes(fill = corr)) + 
      ggplot2::geom_text(ggplot2::aes(label = signif(corr, 3)), fontface = "bold", size = 8) + 
      ggplot2::theme(axis.text.x = ggplot2::element_text(size = 10, face = "bold", color = "black"), 
                     axis.text.y = ggplot2::element_text(size = 10,face = "bold", color = "black"), 
                     axis.title.x = ggplot2::element_text(size = 0, face = "bold", color = "black", vjust = -0.3), 
                     axis.title.y = ggplot2::element_text(size = 0, face = "bold", color = "black"), legend.position = "none") + 
      scale_x_discrete(labels = function(x) {
        gsub("_", ":", x)
      }, expand = c(0, 0)) + scale_y_discrete(position = "right", 
                                              labels = function(x) {
                                                gsub("_", ":", x)
                                              }, expand = c(0, 0)) + scale_fill_continuous(high = "#FF0000", 
                                                                                           low = "white", na.value = "white")
    return(ldplot)
  }
}


plot_bar_phylo <- function(df){
  
  df <- df%>%
    dplyr::mutate(strain3=ifelse(build_GT=="A","",strain))
  
  lb <- df$strain3
  
  df%>%
    ggplot()+
    aes(x=strain2, y=final_pheno, fill = factor(build_GT,
                                                levels = sort(unique(build_GT)), 
                                                labels= c("None", sort(unique(build_GT)[2:8]) )))+
    geom_bar(stat="identity", color = "black", size = .2) +
    scale_fill_manual(values=c("None" = "gray75","Deletion" = colors[1],
                               "Insertion" = colors[2],
                               "Inversion" = colors[3],
                               "Missense" = colors[4],
                               "Splice Donor" = colors[5],
                               "Stop Gained" = colors[6],
                               "Transposon Insertion" = colors[7]),name = expression(paste("Variation at ", italic("ben-1"))))+
    theme_bw()+
    labs(x = "Strain", y = paste0("Albendazole Resistance"))+
    coord_flip()+
    theme(axis.text.x = ggplot2::element_text(size = 14),
          axis.title.x = ggplot2::element_text(size = 18, face = "bold", color = "black", vjust = -0.3),
          axis.title.y = ggplot2::element_text(size = 18, face = "bold", color = "black", vjust = -0.3),
          strip.text.x = element_text(size = 14, face = "bold"))+
    theme(axis.text.y = element_text(size=ifelse(lb=="",0,12)))
  
} 

run_TukeyHSD <- function(df, trt){
  stat_df <- df %>%
    dplyr::filter(trait == trt)%>%
    dplyr::select(group, phenotype)
  
  aov_res <- aov(stat_df$phenotype ~ stat_df$group)
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

# gt_matrix is a genotype matrix data frame that has markers as columns and samples as rows
# do not include column of strain names
gene_level_TajimasD <- function(gt_matrix){
  
  # mean pairwise differences
  d <- dist(gt_matrix, method="manhattan")
  pwd = mean(d)
  
  # number of sequences
  ns <- nrow(gt_matrix)
  # number of segregating sites
  s <- ncol(gt_matrix)
  
  # initialize vector of values
  con=vector(length=10)
  
  # from - https://ocw.mit.edu/courses/health-sciences-and-technology/hst-508-quantitative-genomics-fall-2005/study-materials/tajimad1.pdf
  # a1
  con[1]=sum(1/c(1:(ns-1)))
  # a2
  con[2]=sum((1/c(1:(ns-1)))^2)
  # b1
  con[3]=(ns+1)/(3*ns-3)
  # b2
  con[4]=2*(ns*ns+ns+3)/(9*ns*ns-9*ns)
  # c1
  con[5]=con[3]-1/con[1]
  # c2
  con[6]=con[4]-(ns+2)/(ns*con[1])+con[2]/(con[1]*con[1])
  # e1 
  con[7]=con[5]/con[1]
  # e2
  con[8]=con[6]/(con[1]*con[1]+con[2])
  # 
  con[9]=2*(ns*con[1]-2*(ns-1))/((ns-1)*(ns-2))
  # 
  con[10]=con[9]+(ns-2)/((ns-1)*(ns-1))+((2/(ns-1))*(1.5-(2*(con[1]+1/ns)-3)/(ns-2)-1/ns))
  # 
  con[11]=(ns*ns*con[2]/((ns-1)*(ns-1))+con[1]*con[1]*con[10]-2*ns*con[1]*(con[1]+1)/((ns-1)*(ns-1)))/(con[1]*con[1]+con[2])
  #
  con[12]=ns/(ns-1)*(con[1]-ns/(ns-1))-con[11]
  
  # Wattersons_Theta <- round(s/con[1],2) # Watterson's Theta
  # Average_PWD <- round(pwd,2) # average pairwise differences
  # TajimasD <- round((pwd-s/con[1])/sqrt(con[7]*s+con[8]*s*(s-1)),2) # Tajima's D
  # 
  fq<-apply(gt_matrix, 2, sum, na.rm=T)
  n.na<-apply(is.na(gt_matrix), 2, sum)
  n.sing <- sum(fq==1)+sum(fq==(gt_matrix-n.na-1))
  
  # FuLiF <- round((ns/(ns-1)*s-con[1]*n.sing)/sqrt(con[12]*s+con[11]*s*s),2) # Fu and Li's D
  
  gene_summary <- data.frame(Wattersons_Theta = round(s/con[1],2),
                             Average_PWD = round(pwd,2) ,
                             TajimasD = round((pwd-s/con[1])/sqrt(con[7]*s+con[8]*s*(s-1)),2),
                             Singletons = n.sing,
                             FuLiF = round((ns/(ns-1)*s-con[1]*n.sing)/sqrt(con[12]*s+con[11]*s*s),2))
  return(gene_summary)
}




isnt_out_tukey <- function(x, k = 1.5, na.rm = TRUE) {
  quar <- quantile(x, probs = c(0.25, 0.75), na.rm = na.rm)
  iqr <- diff(quar)
  
  (quar[1] - k * iqr <= x) & (x <= quar[2] + k * iqr)
}
isnt_out_z <- function(x, thres = 3, na.rm = TRUE) {
  abs(x - mean(x, na.rm = na.rm)) <= thres * sd(x, na.rm = na.rm)
}
isnt_out_mad <- function(x, thres = 3, na.rm = TRUE) {
  abs(x - median(x, na.rm = na.rm)) <= thres * mad(x, na.rm = na.rm)
}

isnt_out_funs <- funs(
  z = isnt_out_z,
  mad = isnt_out_mad,
  tukey = isnt_out_tukey
)

