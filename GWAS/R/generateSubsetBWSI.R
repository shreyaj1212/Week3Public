## Generate a smaller dataset for testing

## Read all of the data and sample 1000 SNPs for each
## Do the controls once only

library(tidyverse)
outdir <- "./Data/BWSI_100set/"
basedir <- "./Data/"

diseases <- c ("T1D", "T2D", "HT", "CAD", "RA", "CD", "BD", "58C", "NBS")
#diseases <- c ("58C", "NBS")

#diseases <- c ("T1D")
#controls <- c ("58C", "NBS")

chromosomes <- c(1:22)

walk (diseases, function (d) {
  walk (chromosomes, function (c){
    c = ifelse (as.numeric(c) < 10, str_glue("0{c}"), c)
    df = loadFile(d,c,basedir)
    df = df %>% slice(1:100)
    write_delim(df, str_glue("{outdir}/{d}/Affx_gt_{d}_Chiamo_{c}.tped.gz"), 
                delim = "\t", col_names = F)
  })
})


a <- loadFile("T1D", "01", outdir)
