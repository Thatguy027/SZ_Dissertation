plot_riail_geno <- function(riail_gt) {
  
  df$index <- dplyr::group_indices(riail_gt)
  
  strain_index <- df$sample
  names(strain_index) <- df$index + 0.5
  
  r_gt <- ggplot(df,  aes(xmin = start, xmax = end, ymin = index, ymax = index + 1, fill = gt)) +
    geom_rect(aes(alpha = low_sites)) +
    scale_alpha_discrete(range = c(1.0, 0.65)) +
    scale_fill_manual(values = strain_colors) +
    facet_grid(.~chrom, scales="free", space="free") +
    scale_x_continuous(labels = function(x) { x/1e6 }, expand = c(0,0)) +
    scale_y_continuous(breaks = unique(df$index) + 0.5, labels = function(x) { strain_index[as.character(x)] }, expand = c(0,0)) + 
    theme(strip.background = element_blank(),
          axis.text.y = element_blank(),
          legend.position = "None")
  return(r_gt)
}



maxlodplot_edit <- function(map){
  
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
      ggplot2::geom_point(data = cis, ggplot2::aes(x=pos/1e6, y=(1+maxlod)),
                          fill ="red", shape=25, size=2, show_guide = FALSE) +
      ggplot2::geom_text(data = cis,
                         ggplot2::aes(x=pos/1e6,
                                      y=(ifelse(lod < 10, 2+lod, 2+lod)),
                                      label = paste0(100*round(var_exp, digits = 4),"%")),
                         colour = "black", size=3)
  }
  
  plot <- plot + ggplot2::geom_line(size = 0.75, alpha = 1, color = axis_color) +
    ggplot2::facet_grid(.~chr, scales ="free", space = "free") +
    ggplot2::labs(x = "Genomic Position (Mb)", y = "LOD") +
    ggplot2::scale_colour_discrete(name="Mapping\nIteration")
  return(plot)
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
    ggbeeswarm::geom_beeswarm(ggplot2::aes(x = genotype, y = pheno), alpha = point_alpha, priority = "density",cex = 1.2) +
    ggplot2::geom_boxplot(ggplot2::aes(x = genotype, y = pheno, fill = genotype), outlier.shape = NA,alpha = boxplot_alpha) +
    ggplot2::scale_fill_manual(values = strain_colors) +
    ggplot2::facet_wrap(~ marker, ncol = 5) +
    ggplot2::ggtitle(peaks$trait[1]) +
    ggplot2::labs(x = "Genotype", y = "Phenotype")
}


cegwas2_manplot <- function(plot_df, 
                            bf_line_color = "gray",
                            eigen_line_color = "gray",
                            eigen_cutoff = independent_test_cutoff,
                            mapped_cutoff = "BF") {
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
    
    pt <- plot_df_pr  %>%
      ggplot2::ggplot(.) +
      ggplot2::aes(x = POS/1e6, y = log10p) +
      ggplot2::scale_color_manual(values = c("0" = "black", 
                                             "1" = "red",
                                             "2" = "hotpink3")) +
      ggplot2::scale_alpha_manual(values = c("0" = 0.5, 
                                             "1" = 1,
                                             "2" = 1)) +
      ggplot2::geom_hline(ggplot2::aes(yintercept = bf_cut),
                          color = bf_line_color, 
                          alpha = .75,  
                          size = 1) +
      ggplot2::geom_hline(ggplot2::aes(yintercept = EIGEN_CUTOFF),
                          color = eigen_line_color, 
                          alpha = .75,  
                          size = 1,
                          linetype = 2) +
      ggplot2::facet_grid( . ~ CHROM, scales = "free_x" , space = "free_x") +
      ggplot2::theme_bw() +
      ggplot2::labs(x = "Genomic Position (Mb)",
                    y = expression(-log[10](italic(p))),
                    title = plot_traits)
    if(mapped_cutoff == "BF") {
      pt <- pt +
        ggplot2::geom_rect(ggplot2::aes(xmin = startPOS/1e6, 
                                        xmax = endPOS/1e6, 
                                        ymin = 0, 
                                        ymax = Inf, 
                                        fill = "blue"), 
                           color = "blue",fill = "cyan",linetype = 2, 
                           alpha=.3, data = dplyr::filter(plot_df_pr, EIGEN_SIG=="1") %>% na.omit()) 
      if(grepl("2", unique(plot_df_pr$EIGEN_SIG) )){
       pt <- pt + ggplot2::geom_rect(ggplot2::aes(xmin = startPOS/1e6, 
                                             xmax = endPOS/1e6, 
                                             ymin = 0, 
                                             ymax = Inf, 
                                             fill = "hotpink3"), 
                                color = "hotpink",linetype = 2, 
                                alpha=.3, data = dplyr::filter(plot_df_pr, EIGEN_SIG!="1") %>% na.omit())
      }
      pt <- pt + ggplot2::geom_point( ggplot2::aes(color= factor(EIGEN_SIG), alpha = factor(EIGEN_SIG)), size = 1 ) 
    } else {
      pt +  ggplot2::geom_rect(ggplot2::aes(xmin = startPOS/1e6, 
                                            xmax = endPOS/1e6, 
                                            ymin = 0, 
                                            ymax = Inf, 
                                            fill = "hotpink3"), 
                               color = "hotpink",linetype = 2, 
                               alpha=.3, data = dplyr::filter(plot_df_pr, EIGEN_SIG!="1") %>% na.omit()) +
        ggplot2::geom_point( ggplot2::aes(color= factor(EIGEN_SIG), alpha = factor(EIGEN_SIG)), size = 1 )
    }
  })
  plots
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
    geom_beeswarm(priority = "density", cex = 0.6, size = point_size, alpha = point_alpha) +
    geom_boxplot(outlier.colour = NA, alpha = boxplot_alpha)+
    scale_fill_manual(values = c("Bristol" = "#F9A227","Hawaii" = "#2790F9",
                                 "Bristol\n(C78S)" = "gray","Bristol\n(C78S)" = "gray",
                                 "Hawaii\n(S78C)" = "gray","Hawaii\n(S78C)" = "gray"),
                      name="Strain")+
    labs( y = trt)+
    labs(y = fancy_name)
}

Plot_Peak_LD <- function(processed_mapping, genotype_matrix){
  snp_df <- processed_mapping %>% na.omit()
  ld_snps <- dplyr::filter(genotype_matrix, CHROM %in% snp_df$CHROM, POS %in% snp_df$POS)
  ld_snps <- data.frame(snp_id = paste(ld_snps$CHROM, ld_snps$POS, sep = "_"), data.frame(ld_snps)[, 5:ncol(ld_snps)])
  sn <- list()
  for (i in 1:nrow(ld_snps)) {
    sn[[i]] <- genetics::genotype(as.character(gsub(1, "T/T", gsub(-1, "A/A", ld_snps[i, 4:ncol(ld_snps)]))))
  }
  
  test <- data.frame(sn)
  colnames(test) <- (ld_snps$snp_id)
  ldcalc <- t(genetics::LD(test)[[4]])^2
  diag(ldcalc) <- 1
  LDs <- tbl_df(data.frame(ldcalc) %>% dplyr::add_rownames(var = "SNP1")) %>% 
    tidyr::gather(SNP2, r2, -SNP1) %>% dplyr::arrange(SNP1) %>% 
    tidyr::separate(SNP1, sep = "_", into = c("CHROM1", "POS1"), remove = F) %>% 
    dplyr::arrange(CHROM1,  as.numeric(POS1))
  
  ldplot <- ggplot2::ggplot(LDs) + 
    ggplot2::aes(x = factor(SNP1, levels = unique(SNP1), ordered = T), 
                 y = factor(SNP2, levels = unique(SNP1), ordered = T)) + 
    ggplot2::geom_tile(ggplot2::aes(fill = r2), color = background_color) + 
    ggplot2::geom_text(ggplot2::aes(label = signif(r2, 3)), fontface = "bold", size = 4.5) + 
    scale_x_discrete(labels = function(x) {
      gsub("_", ":", x)
    }, expand = c(0, 0)) + 
    scale_y_discrete(position = "right",  
                     labels = function(x) {gsub("_", ":", x)}, 
                     expand = c(0, 0)) + 
    scale_fill_continuous(high = "#FF0000", low = "white", na.value = "white")
  
  return(list(ldplot, LDs))
}

rial_bar_plot <- function(cross, map, parent="N2xCB4856", color_by_genotype = TRUE) {
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
    dplyr::select_( phenotype = map$trait[1])
  
  geno <- data.frame(extract_genotype(cross)) %>%
    dplyr::select(which(colnames(.) %in% uniquemarkers)) %>%
    data.frame(., pheno) %>%
    na.omit() %>%
    dplyr::arrange(phenotype) %>%
    dplyr::mutate(norm_pheno_temp = ifelse(phenotype == min(phenotype), 0, 1))%>%
    dplyr::mutate(delta_pheno = ifelse(norm_pheno_temp == 0, 0, abs(dplyr::lag(phenotype) - phenotype)))%>%
    dplyr::mutate(norm_pheno = cumsum(delta_pheno)) %>%
    dplyr::mutate(final_pheno = norm_pheno/max(norm_pheno)) %>%
    dplyr::select(-(phenotype:norm_pheno))
  
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
  
  split <- geno%>%
    dplyr::mutate(strain = factor(1:n()))%>%
    tidyr::gather(marker, genotype, -pheno, -strain)
  
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
  
  plots <- list()
  for(qtl in 1:length(unique(split$marker))){
    plots[[qtl]] <- split%>%
      dplyr::filter(marker == unique(split$marker)[qtl]) %>%
      ggplot()+
      aes(x=strain, y= pheno, fill = genotype)+
      geom_bar(stat="identity", color = "black", size = .1) +
      facet_grid(marker~.)+
      labs(x = "Recombinant inbred line", y = paste0("Arsenic sensititivity")) 
    
    if(color_by_genotype == TRUE){
      plots[[qtl]] <- plots[[qtl]] + scale_fill_manual(values= strain_colors,
                        name = expression(paste("Genotype at QTL")))
    } else {
      plots[[qtl]] <- plots[[qtl]] + scale_fill_manual(values= c("#999999","#999999"),
                        name = expression(paste("Genotype at QTL")))
    }
  }
  return(plots)
}

