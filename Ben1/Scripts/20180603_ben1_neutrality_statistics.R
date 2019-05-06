library(PopGenome)
vcf_name = "WI.20170531.impute.vcf.gz"
slide_distance = 10
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

st = 3537688
en = 3541628

gen <- PopGenome::readVCF(vcf_name, 
                          numcols = 1000, 
                          tid = chroms[w.CHROM], 
                          frompos =  st-1e5, 
                          topos =  en+1e5, 
                          approx = F,
                          gffpath = "wormbase.gff",
                          include.unknown=TRUE)

# Set most diverged strain as the outgroup

gen2 <- PopGenome::set.outgroup(gen, c("XZ1516"),  diploid = FALSE)

genes <- set.synnonsyn(gen2, ref.chr="III.fa")
genes <- splitting.data(genes, subsites="gene")
genes <- neutrality.stats(genes, subsites="nonsyn",detail = T)
nonsynTaj <- genes@Tajima.D
nonsynH <- genes@Fay.Wu.H
nonsynE<- genes@Zeng.E
nonsynFuD<- genes@Fu.Li.D

genes <- PopGenome::neutrality.stats(genes, subsites="syn",detail = T)
synTaj <- genes@Tajima.D
synH <- genes@Fay.Wu.H
synE<- genes@Zeng.E
synFuD<- genes@Fu.Li.D

genes <- PopGenome::neutrality.stats(genes, subsites="intron",detail = T)
intronTaj <- genes@Tajima.D
intronH <- genes@Fay.Wu.H
intronE<- genes@Zeng.E
intronFuD<- genes@Fu.Li.D

genes <- PopGenome::neutrality.stats(genes, subsites="coding",detail = T)
codingTaj <- genes@Tajima.D
codingH <- genes@Fay.Wu.H
codingE<- genes@Zeng.E
codingFuD<- genes@Fu.Li.D

gene_names <- get.feature.names(object = genes,chr = chroms[w.CHROM],gff.file = "wormbase.gff")
Td <- data.frame(TajimaD_non = nonsynTaj, 
                 FayWuH_non = nonsynH,
                 ZengE_non = nonsynE,
                 FuLiD_non = nonsynFuD,
                 TajimaD_syn = synTaj, 
                 FayWuH_syn = synH,
                 ZengE_syn = synE,
                 FuLiD_syn = synFuD,
                 TajimaD_intron = intronTaj, 
                 FayWuH_intron = intronH,
                 ZengE_intron = intronE,
                 FuLiD_intron = intronFuD,
                 TajimaD_coding = codingTaj, 
                 FayWuH_coding = codingH,
                 ZengE_coding = codingE,
                 FuLiD_coding = codingFuD,
                 gene_id = gene_names)%>%
  na.omit()

colnames(Td) <- c("TajimaD_non","FayWuH_non","ZengE_non","FuLiD_non",
                  "TajimaD_syn","FayWuH_syn","ZengE_syn","FuLiD_syn",
                  "TajimaD_intron","FayWuH_intron","ZengE_intron","FuLiD_intron",
                  "TajimaD_coding","FayWuH_coding","ZengE_coding","FuLiD_coding","gene_id")

