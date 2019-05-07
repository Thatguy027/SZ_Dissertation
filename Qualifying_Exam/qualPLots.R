library(dplyr)
library(ggplot2)
library(data.table)
library(tidyr)
library(stringr)
library(broom)
library(compute.es)

presentation <- theme(axis.text.x = element_text(size=16, face="bold", color="black"),
                      axis.text.y = element_text(size=16, face="bold", color="black"),
                      axis.title.x = element_text(size=20, face="bold", color="black"),
                      axis.title.y = element_text(size=20, face="bold", color="black"),
                      strip.text.x = element_text(size=20,face="bold", color="black"),
                      strip.text.y = element_text(size=20,angle=0, face="bold", color="black"),
                      plot.title = element_text(size=24, face="bold"))

# Figure 1 A & B

source("~/Dropbox/AndersenLab/RCode/Stefan/LinkagePlotFunction.R")
load("~/Dropbox/AndersenLab/RCode/GWAS/Data/processed_3it_linkage.Rda")

linkagePlot("etoposide_resid.n")
ggsave(filename = "~/Dropbox/AndersenLab/LabFolders/Stefan/Qualifying_Exam/figures/LOD_etop_resid.n.pdf",
       width=14,
       height=6)

linkagePlot("amsacrine_resid.median.TOF")
ggsave(filename = "~/Dropbox/AndersenLab/LabFolders/Stefan/Qualifying_Exam/figures/LOD_amsa_resid_med-length.pdf",
       width=14,
       height=6)

linkagePlot("amsacrine_resid.q90.EXT")
ggsave(filename = "~/Dropbox/AndersenLab/LabMeetings/Stefan/20150413_OncDevBio/LOD_amsa_resid_q90ext.pdf",
       width=14,
       height=6)

# Table 1 - Number of genes in interval. Find intervals in Mapping reports

source("~/Dropbox/AndersenLab/RCode/Stefan/Interval_Genes.R")
# amsacrine
ama <- browsePXGvariantSPLITS("V",10124453 , 10839540)
# etoposide
etop <- browsePXGvariantSPLITS("II",11492171 , 12008179)
# etoposide
etopV <- browsePXGvariantSPLITS("V",9903872 , 12125475)

etopNamsa <- data.frame(etopV$NAME[!(etopV$NAME%in%ama$NAME)])

# Figure 2 - Phenotype by Genotype Plots
source("~/Dropbox/AndersenLab/RCode/Stefan/linkage_PxG_Function.R")
load("~/Dropbox/AndersenLab/LabFolders/Stefan/Mapping_Reports/Data/Linkage_Genotypes.Rda")
load("~/Dropbox/AndersenLab/LabFolders/Stefan/Mapping_Reports/Data/Linkage_Phenotypes.Rda")


PxGlinkage("UCE2-2047","etoposide_resid.n")
ggsave(filename = "~/Dropbox/AndersenLab/LabFolders/Stefan/Qualifying_Exam/figures/pxg_linkage_etop_chr2.pdf",
       width=10,
       height=6)

PxGlinkage("UCE5-1861","etoposide_resid.n")
ggsave(filename = "~/Dropbox/AndersenLab/LabFolders/Stefan/Qualifying_Exam/figures/pxg_linkage_etop_chrv.pdf",
       width=10,
       height=6)

PxGlinkage("UCE5-1839","amsacrine_resid.median.TOF")
ggsave(filename = "~/Dropbox/AndersenLab/LabFolders/Stefan/Qualifying_Exam/figures/pxg_linkage_amsacrine_chrv.pdf",
       width=10,
       height=6)

####### plot that splits base on two variants
condtrt <- "etoposide_resid.n"
peakMarker <- c("UCE2-2047","UCE5-1861")

temp <- str_split_fixed(condtrt,pattern="_",n=2)
trt <- temp[,2]
condition <- temp[,1]


gen1 <- Lgeno %>%
  filter(snp%in%peakMarker[1])%>%
  arrange(strain)%>%
  select(strain,geno1=geno)

gen2 <- Lgeno %>%
  filter(snp%in%peakMarker[2])%>%
  arrange(strain)%>%
  select(strain,geno2=geno)%>%
  left_join(.,gen1,by="strain")


Lphen%>%
  filter(grepl(condition,trait))%>%
  filter(grepl(trt,trait))%>%
  left_join(.,gen2,by="strain")%>%
  mutate(color = ifelse(strain=="N2", "N2",
                        ifelse(strain == "CB4856", "CB4856",
                               ifelse(geno1 == 1&geno2==1, "RIAILs-N2-N2",
                                      ifelse(geno1 == 2&geno2==2, "RIAILs-CB-CB",
                                             ifelse(geno1 == 1&geno2==2, "RIAILs-N2-CB",
                                                    ifelse(geno1 == 2&geno2==1, "RIAILs-CB-N2",0)))))),
         label = ifelse(strain=="N2", "N2",
                        ifelse(strain == "CB4856", "CB4856",
                               ifelse(geno1 == 1&geno2==1, "RIAILs-N2-N2",
                                      ifelse(geno1 == 2&geno2==2, "RIAILs-CB-CB",
                                             ifelse(geno1 == 1&geno2==2, "RIAILs-N2-CB",
                                                    ifelse(geno1 == 2&geno2==1, "RIAILs-CB-N2",0)))))))%>%
  filter(!is.na(color))%>%
  ggplot(.)+
  aes(x = color, y = value,fill = ifelse(strain=="N2", "N2",
                                         ifelse(strain == "CB4856", "CB4856",
                                                ifelse(geno1 == 1, "RIAILs-N2",
                                                       ifelse(geno1 == 2, "RIAILs-CB",0)))),
      color = ifelse(strain=="N2", "N2",
                     ifelse(strain == "CB4856", "CB4856",
                            ifelse(geno1 == 1, "RIAILs-N2",
                                   ifelse(geno1 == 2, "RIAILs-CB",0)))))+
  scale_fill_manual(values = c("N2"="orange","CB4856"="blue","RIAILs-CB"="blue","RIAILs-N2"="orange"), name = "Genotype")+
  scale_color_manual(values = c("N2"="black","CB4856"="black","RIAILs-CB"="#666666","RIAILs-N2"="#666666"),name="Genotype")+
  geom_boxplot(outlier.shape=NA)+
  theme_bw()+
  theme(axis.text.x = element_text(size=16, face="bold", color="black"),
        axis.text.y = element_text(size=16, face="bold", color="black"),
        axis.title.x = element_text(size=20, face="bold", color="black"),
        axis.title.y = element_text(size=20, face="bold", color="black",vjust=1),
        strip.text.x = element_text(size=20, face="bold", color="black"),
        strip.text.y = element_text(size=20, face="bold", color="black"),
        plot.title = element_text(size=24, face="bold",vjust=1),
        legend.title = element_text(size=14),
        panel.border = element_rect(size=2))+
  geom_jitter(alpha=.75, size=2)


# Figure 3 A - Dominance Results
load("~/Dropbox/AndersenLab/LabFolders/Stefan/Qualifying_Exam/figures/processed_resid_Dominance_Data.Rda")

regressed %>%
  dplyr::filter(trait == "n")%>%
  dplyr::mutate(strain1 = factor(strain, levels = c("CB4856" ,"N2"   , "N2GFP",  "CB4856GFP","GFP"   ),
                          labels = c("CB4856","N2", "N2/N2-GFP", "CB4856/N2-GFP","N2-GFP/N2-GFP")))%>%
  # filter(strain1 != "N2")%>%
  ggplot(.)+
  aes(x = strain1, y = .resid, color = strain, fill = strain)+
  scale_color_manual(values = c("black", "black","black", "black","black"))+
  scale_fill_manual(values = c("blue", "blue", "orange","orange","orange"))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(size =4)+
  theme_bw()+
  presentation+ 
  theme(legend.position="none",
        panel.border = element_rect(size=1, colour = "black"))+
  labs(x = "Strain", y = "Animal Size")

ggsave(filename = "~/Dropbox/AndersenLab/LabFolders/Stefan/Qualifying_Exam/figures/etoposide_top2_dom_n.pdf",
       width=10,
       height=6)

test <- regressed%>%filter(trait == "n")

am = lm(test$Evalue ~ test$strain-1)
aov.a=anova(am)
tssq = sum(aov.a[,2])
a.var.exp = aov.a[1:(nrow(aov.a)-1),2]/tssq  
a.eff.size= as.vector(coefficients(am))

test <- regressed%>%filter(trait == "n" & strain %in% c("CB4856","N2"))

am = lm(test$Evalue ~ test$strain-1)
aov.a=anova(am)
tssq = sum(aov.a[,2])
a.var.exp = aov.a[1:(nrow(aov.a)-1),2]/tssq  
a.eff.size= as.vector(coefficients(am))


test <- regressed%>%filter(trait == "n" & strain %in% c("CB4856","GFP"))

am = lm(test$.resid ~ test$strain-1)
aov.a=anova(am)
tssq = sum(aov.a[,2])
a.var.exp = aov.a[1:(nrow(aov.a)-1),2]/tssq  
a.eff.size= as.vector(coefficients(am))

test <- regressed%>%filter(trait == "n" & strain %in% c("CB4856","CB4856GFP"))

am = lm(test$.resid ~ test$strain-1)
aov.a=anova(am)
tssq = sum(aov.a[,2])
a.var.exp = aov.a[1:(nrow(aov.a)-1),2]/tssq  
a.eff.size= as.vector(coefficients(am))

test <- regressed%>%filter(trait == "n" & strain %in% c("CB4856","N2GFP"))

am = lm(test$.resid ~ test$strain-1)
aov.a=anova(am)
tssq = sum(aov.a[,2])
a.var.exp = aov.a[1:(nrow(aov.a)-1),2]/tssq  
a.eff.size= as.vector(coefficients(am))


# Figure 3 B Complementation results

load("~/Dropbox/AndersenLab/LabFolders/Stefan/Qualifying_Exam/figures/Etoposide_complementation_top2_processed_Data.Rda")

regressed %>%
  filter(trait == "n")%>%
  # filter(!(strain %in% c("CB4856min1","N2min1")))%>%
  mutate(strain1 = factor(strain, levels = c("CB4856" ,"N2"  ,  "CB4856top2" , "N2top2","CB4856min1","N2min1"   ),
                          labels = c("CB4856","N2",  "CB4856/∆top-2","N2/∆top-2","CB4856min1","N2min1")))%>%
  ggplot(.)+
  aes(x = strain1, y = .resid, color = strain, fill = strain)+
  scale_color_manual(values = c("black", "black","black", "black","black", "black"))+
  scale_fill_manual(values = c("blue", "blue", "gray","orange","orange","gray"))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter()+
  theme_bw()+
  presentation+ 
  theme(legend.position="none",
        panel.border = element_rect(size=1, colour = "black"))+
  labs(x = "Strain", y = "Brood Size")


regressed %>%
  filter(trait == "n")%>%
  filter(!(strain %in% c("CB4856min1","N2min1")))%>%
  mutate(strain1 = factor(strain, levels = c("CB4856" ,"N2"  ,  "CB4856top2" , "N2top2"  ),
                          labels = c("CB4856","N2",  "CB4856/∆top-2","N2/∆top-2")))%>%
  ggplot(.)+
  aes(x = strain1, y = .resid, color = strain, fill = strain)+
  scale_color_manual(values = c("black", "black","black", "black"))+
  scale_fill_manual(values = c("blue", "gray", "orange","gray"))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter()+
  theme_bw()+
  presentation+ 
  theme(legend.position="none",
        panel.border = element_rect(size=1, colour = "black"))+
  labs(x = "Strain", y = "Brood Size")
ggsave(filename = "~/Dropbox/AndersenLab/LabFolders/Stefan/Qualifying_Exam/figures/etoposide_top2_comp_nv1.pdf",
       width=12,
       height=6)

# # # # npp-3

load("~/Dropbox/AndersenLab/LabFolders/Stefan/WetLab/Experiments/20150413_etoposide_npp3_Comp/regressed.Rda")

regressed %>%
  filter(trait == "n")%>%
  # filter(!(strain %in% c("CB4856min1","N2min1")))%>%
  mutate(strain1 = factor(strain, levels = c("CB4856" ,"N2"  ,  "CB4856npp3" , "N2npp3" ,"CB4856min1","N2min1"  ),
                          labels = c("CB4856","N2",  "CB4856/∆npp-3","N2/∆npp-3","CB4856min1","N2min1")))%>%
  ggplot(.)+
  aes(x = strain1, y = .resid, color = strain, fill = strain)+
  scale_color_manual(values = c("black", "black","black", "black","black", "black"))+
  scale_fill_manual(values = c("blue", "blue", "gray","orange","orange","gray"))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter()+
  theme_bw()+
  presentation+ 
  theme(legend.position="none",
        panel.border = element_rect(size=1, colour = "black"))+
  labs(x = "Strain", y = "Brood Size")

regressed %>%
  filter(trait == "n")%>%
  filter(!(strain %in% c("CB4856min1","N2min1")))%>%
  mutate(strain1 = factor(strain, levels = c("CB4856" ,"N2"  ,  "CB4856npp3" , "N2npp3"),
                          labels = c("CB4856","N2",  "CB4856/∆npp-3","N2/∆npp-3")))%>%
  ggplot(.)+
  aes(x = strain1, y = .resid, color = strain, fill = strain)+
  scale_color_manual(values = c("black", "black","black", "black"))+
  scale_fill_manual(values = c("blue", "gray", "orange","gray"))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter()+
  theme_bw()+
  presentation+ 
  theme(legend.position="none",
        panel.border = element_rect(size=1, colour = "black"))+
  labs(x = "Strain", y = "Brood Size")
ggsave(filename = "~/Dropbox/AndersenLab/LabFolders/Stefan/Qualifying_Exam/figures/etoposide_npp3_comp_n.pdf",
       width=12,
       height=6)

# p values
prune <- regressed%>%
  filter(trait == "n")%>%
  filter(!(strain %in% c("CB4856min1","N2min1")))%>%
  select(trait,.resid,strain)

pairwise.t.test(prune$.resid,prune$strain,pool.sd = F,p.adj = "bonf")


# Appendix Figure - GWAS

load("~/Dropbox/AndersenLab/LabFolders/Stefan/Mapping_Reports/Data/processed_resid_GWAS_mappings.Rda")

Gmaps%>%
    filter(grepl("etoposide",pheno))%>%
    filter(grepl("resid.q90.TOF", pheno))%>%
    ggplot(.)+
    aes(x=pos/1e6,y=log10p,color=ifelse(log10p>BF,"1","0"))+
    geom_point()+
    scale_color_manual(values=c("1"="red","0"="black"))+
    facet_grid(.~chr, scales="free")+
    theme_bw()+
    geom_hline(aes(yintercept=BF), color = "gray",linetype="dashed")+
    theme(axis.text.x = element_text(size=16, face="bold", color="black"),
          axis.text.y = element_text(size=16, face="bold", color="black"),
          axis.title.x = element_text(size=20, face="bold", color="black", vjust=-.3),
          axis.title.y = element_text(size=20, face="bold", color="black"),
          strip.text.x = element_text(size=20, face="bold", color="black"),
          strip.text.y = element_text(size=20, face="bold", color="black"),
          plot.title = element_text(size=24, face="bold", vjust = 1),
          legend.position="none",
          panel.border = element_rect(size=1, colour = "black"),
          strip.background = element_rect(size = 1))+
    labs(x = "Genomic Position (Mb)", y = expression(-log[10](p)))
ggsave(filename = "~/Dropbox/AndersenLab/LabFolders/Stefan/Qualifying_Exam/figures/q90TOF_etop_manplot.pdf",
       width=14,
       height=6)

Gmaps%>%
  filter(grepl("amsacrine",pheno))%>%
  filter(grepl("resid.q90.TOF", pheno))%>%
  ggplot(.)+
  aes(x=pos/1e6,y=log10p,color=ifelse(log10p>BF,"1","0"))+
  geom_point()+
  scale_color_manual(values=c("1"="red","0"="black"))+
  facet_grid(.~chr, scales="free")+
  theme_bw()+
  geom_hline(aes(yintercept=BF), color = "gray",linetype="dashed")+
  theme(axis.text.x = element_text(size=16, face="bold", color="black"),
        axis.text.y = element_text(size=16, face="bold", color="black"),
        axis.title.x = element_text(size=20, face="bold", color="black", vjust=-.3),
        axis.title.y = element_text(size=20, face="bold", color="black"),
        strip.text.x = element_text(size=20, face="bold", color="black"),
        strip.text.y = element_text(size=20, face="bold", color="black"),
        plot.title = element_text(size=24, face="bold", vjust = 1),
        legend.position="none",
        panel.border = element_rect(size=1, colour = "black"),
        strip.background = element_rect(size = 1))+
  labs(x = "Genomic Position (Mb)", y = expression(-log[10](p)))

# histogram of GWAS data

load("/Users/stefan/Dropbox/AndersenLab/RCode/Stefan/good_gwasMappingsINlinkage_phenotypes.Rda")

hm%>%
  filter(condition == "etoposide" & trait == "q90.TOF")%>%
  ggplot(.)+
  aes(x=pheno)+
  geom_histogram(color = "black",fill="black")+
  theme_bw()+
  presentation+
  labs(x="Body Size (µM)", y = "Count")
ggsave("~/Dropbox/AndersenLab/LabFolders/Stefan/Qualifying_Exam/figures/etop_q90hist.pdf",
       width = 10,
       height= 6)

## map variants

library(devtools)
source_gist("da9a54435b5d7ee4909d")
map_alleles(CHROMOSOME = "II", 
            POSITION = 14524396,
            ref_color = "#B20000",
            alt_color = "#B20000")

## Andersen RIAIL genotype

library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(grid)

genos <- data.frame(fread("~/Dropbox/AndersenLab/RCode/GWAS/Ancillary/N2xCB4856_598RIAILs_gen.csv",
                          header = T))

load("~/Dropbox/AndersenLab/RCode/GWAS/Ancillary/marker_pos_conversion.Rda")

test <- do.call(cbind, lapply(genos[,4:ncol(genos)], function(x){
  temp <- data.frame(gen = x)
  temp <- mutate(temp, num = ifelse(gen =="AA",1,
                                    ifelse(gen=="AB",0,NA)))
  temp <- select(temp, -gen)
}))

test1 <- data.frame(genos[,1:3], test)
colnames(test1) <- c("id","chr","CM",seq(1,ncol(test1)-3))
test1 <- left_join(test1, conv, by= "id")

riails <-c(240:598)
ch <- c("I","II","III","IV","V","X")

test1 %>%
  gather(riail,geno, -id, -chr,-CM,-pos)%>%
  filter(riail %in% riails)%>%
  mutate(cbgen = ifelse(geno==1,0,
                        ifelse(geno==0,1,NA)))%>%
  filter(chr %in% ch)%>%
  group_by(riail)%>%
  mutate(marker=seq(1:n()))%>%
  ungroup()%>%
  group_by(riail,chr)%>%
  mutate(cst = seq(1:n()))%>%
  mutate(li = ifelse(cst==1,marker,0))%>%
  ggplot(.)+
  aes(x = marker, y = geno)+
  facet_grid(riail~., scales = "free")+
  geom_ribbon(aes(ymin = 0, ymax = geno), fill ="orange")+
  geom_ribbon(aes(ymin = 0, ymax = cbgen), fill ="blue")+
  geom_vline(aes(xintercept=li))+
  theme_bw()+
  theme(axis.text.x = element_text(size=16, face="bold", color="black"),
        axis.text.y = element_text(size=0, face="bold", color="black"),
        axis.title.x = element_text(size=20, face="bold", color="black"),
        axis.title.y = element_text(size=20, face="bold", color="black"),
        strip.text.x = element_text(size=20,face="bold", color="black"),
        strip.text.y = element_text(size=0, angle =0, face="bold", color="black"),
        plot.title = element_text(size=24, face="bold"),
        panel.margin = unit(0, "lines"),
        strip.background = element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_blank(),
        legend.position = "none",
        axis.ticks = element_blank())+
  labs(x = "Marker", y = "RIAIL")

ggsave("~/Dropbox/AndersenLab/LabFolders/Stefan/Qualifying_Exam/figures/ECA_RIAIL_GENO1.pdf",
       width = 12,
       height= 6)


# # # # ## #normal distribution


x <- seq(-4, 4, length=100)
hx <- dnorm(x)

degf <- c(1, 3, 8, 30)
colors <- c("red", "blue", "darkgreen", "gold", "black")
labels <- c("df=1", "df=3", "df=8", "df=30", "normal")


df <- data.frame(x = x, dens = hx)

ggplot(df)+
  aes(x = x, y = dens)+
  geom_line()+
  theme_bw()+
  presentation


ggsave("~/Dropbox/AndersenLab/LabFolders/Stefan/Qualifying_Exam/figures/NormalDist.pdf",
       width = 12,
       height= 8)

