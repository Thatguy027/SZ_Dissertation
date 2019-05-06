library(cegwas)
library(tidyverse)
library(rio)

tajimas_d_temp <- function (vcf_path = paste0(path.package("cegwas"), "/"), vcf_name = paste0("WI.", 
                                                                                              vcf_version, ".impute.vcf.gz"), chromosome = "II", interval_start = 11021073, 
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



ben1_d <- tajimas_d_temp(vcf_path = "~/Dropbox/AndersenLab/Reagents/WormReagents/_SEQ/WI/WI-20170531/vcf/",vcf_name = "WI.20170531.impute.vcf.gz", 
                         chromosome = "III", interval_start = 3437688, interval_end = 3641628,site_of_interest = 3539628,slide_distance = 10,window_size = 500)
ben1_d
ggsave("~/Dropbox/HTA/Results/2017_anthelmintic_gwa/analysis/ben1_paper/Plots/ben1_TajimasD.pdf",height= 4,width = 24)

# tbb1
tbb1_d <- tajimas_d_temp(vcf_path = "~/Dropbox/AndersenLab/Reagents/WormReagents/_SEQ/WI/WI-20170531/vcf/",vcf_name = "WI.20170531.impute.vcf.gz", 
                         chromosome = "III", interval_start = 10730868, interval_end = 10752461,site_of_interest = 10741461,slide_distance = 5,window_size = 10)
tbb1_d
ggsave("~/Dropbox/HTA/Results/2017_anthelmintic_gwa/analysis/ben1_paper/Plots/tbb1_TajimasD.pdf",height= 4,width = 24)

# tbb2
tbb2 <- tajimas_d_temp(vcf_name = "WI.20160408.impute.vcf.gz", chromosome = "III", interval_start = 3905769, interval_end = 4027643,site_of_interest = 4016643,slide_distance = 5,window_size = 10) 
tbb2
ggsave("~/Dropbox/HTA/Results/2017_anthelmintic_gwa/analysis/ben1_paper/Plots/ben1_TajimasD.pdf",height= 4,width = 24)

# tbb4 X:9434452..9436893
tbb4 <- tajimas_d_temp(vcf_name = "WI.20160408.impute.vcf.gz", chromosome = "X", interval_start = 9334452, interval_end = 9536893,site_of_interest = 9436893,slide_distance = 5,window_size = 10) 
tbb4
ggsave("~/Dropbox/HTA/Results/2017_anthelmintic_gwa/analysis/ben1_paper/Plots/tbb4_TajimasD.pdf",height= 4,width = 24)

# tbb6 V: 12261804 12263368
tbb6 <- tajimas_d_temp(vcf_name = "WI.20160408.impute.vcf.gz", chromosome = "V", interval_start = 12161804, interval_end = 12363368,site_of_interest = 12263368,slide_distance = 5,window_size = 10) 
tbb6
ggsave("~/Dropbox/HTA/Results/2017_anthelmintic_gwa/analysis/ben1_paper/Plots/tbb6_TajimasD.pdf",height= 4,width = 24)

# mec7  X:7774859 7776689
mec7 <- tajimas_d_temp(vcf_name = "WI.20160408.impute.vcf.gz", chromosome = "X", interval_start = 7674859, interval_end = 7876689,site_of_interest = 7776689,slide_distance = 5,window_size = 10) 
mec7
ggsave("~/Dropbox/HTA/Results/2017_anthelmintic_gwa/analysis/ben1_paper/Plots/mec7_TajimasD.pdf",height= 4,width = 24)


