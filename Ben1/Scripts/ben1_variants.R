library(cegwas)
library(readr)

ben1 <- snpeff("ben-1",severity = "ALL",elements = "ALL")

save(ben1, file = "~/Dropbox/HTA/Results/2017_anthelmintic_gwa/analysis/ben1_paper/Data/ben1variants.Rda")

readr::write_tsv(ben1, path = "~/Dropbox/HTA/Results/2017_anthelmintic_gwa/analysis/ben1_paper/Data/ben1variants.tsv",
          na = "NA", append = FALSE, col_names = T)


ben1 <- snpeff("III:3536000-3546000",severity = "ALL",elements = "ALL")
