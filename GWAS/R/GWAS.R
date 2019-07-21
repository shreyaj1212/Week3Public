source("./R/functions.R")
library(HardyWeinberg)
library(tidyverse)
library(foreach)
library(parallel)

GWASD <- function (disease, basedir, outdir) {
  #start time
  strt<-Sys.time()
  
  # getDoParWorkers()
  # 
  # cl<-makeCluster(2)
  # registerDoParallel(cl)
  
  #start time
  strt<-Sys.time()
  
  results <- foreach(i=1:22, 
                     .combine='rbind', 
                     .export=c('GWAS','gwas', 'homozygousOrMono', 
                               'computeAlleleFrequency','hweChisq', 
                               'chitest', 'loadFile','loadSnpInfoFile'),
                     .packages=c('tidyverse', 'HardyWeinberg','data.table', 
                                 'foreach')) %dopar% 
    GWAS(disease,i, basedir, outdir)
  
  
  # write out the result file
  dim(results)
  print(head(results))
  colnames(results) <-  c("SNP","RSID","CHR", "MajA", "MinA", "MAFD","MAFC","OR", 
                          "P", "CHI2", "HWE" )
  
  resultFile <-  str_glue("{outdir}/{disease}.tsv")
  write.table(results,file = resultFile,quote = FALSE,dec = ".",sep = "\t",
              col.names = NA, row.names = T )
 
  return(results)
  
}

# GWAS Function ######
GWAS <- function (disease, chrom, basedir,outdir) {
  print (str_glue("Disease: {disease}, chrom: {chrom}"))
  # for chromosome num < 10 need to append 0 for file format
  if (as.numeric(chrom) < 10) {
    chrom <- sprintf('0%s',chrom)
  }
  
  # create output directory if it doesnt exist
  if (!(file.exists(outdir))) {
    dir.create( outdir, showWarnings = FALSE, recursive = FALSE,mode = "0777" )
    Sys.chmod(outdir, mode = "0777", use_umask = TRUE)
  }
  
  # load the data
  controlData <- loadFile("58C", chrom, basedir)
  controlData2 <- loadFile("NBS", chrom, basedir)
  dData <- loadFile(disease, chrom, basedir)
  snpfile <- paste(basedir, "/58C/snps", sep = "")
  snpInfo <- loadSnpInfoFile(chrom, snpfile)
  
  # control and t2d format: chNum, snpId, start, finish, <genotype>
  cntrl <- inner_join(controlData, controlData2, by=c("X1", "X2", "X3", "X4"))
  controls <- left_join(cntrl, snpInfo, by=c("X2" = "snpId")) 
  dData <- left_join(dData, snpInfo, by=c("X2" = "snpId"))
  
  alC <- controls %>% select (-c(1:4), -rsid)
  alD <- dData %>% select (-c(1:4), -rsid)
  
  headC <- controls %>% select (c(1:4), rsid)
  headD <- dData %>% select (c(1:4), rsid)
  colnames(headD) <- c('chr_num', 'snpId', 'start', 'finish', "rsid")
  colnames(headC) <- c('chr_num', 'snpId', 'start', 'finish', "rsid")
  
  ## For testing only
  ## Comment when done
   # alC <- alC %>% slice (1:100)
   # alD <- alD %>% slice (1:100)
   # headC <- headC %>% slice (1:100)
   # headD <- headD %>% slice (1:100)

  library(foreach)
  results <- foreach (i=1:nrow(headC), 
                      .combine = 'rbind', 
                      .export=c('gwas', 'homozygousOrMono', 
                                'computeAlleleFrequency','hweChisq', 'chitest')) %dopar%
             gwas(alC[i,], alD[i,], headC[i,], headD[i,])
  

  dim(results)
  print(head(results))
  colnames(results) <-  c("SNP","RSID","CHR", "MajA", "MinA", "MAFD","MAFC","OR", 
                            "P", "HWE" )
  
  resultFile <-  str_glue("{outdir}/{disease}_{chrom}.tsv")
  
  write.table(results, file = resultFile, quote = FALSE, 
              dec = ".",  sep = "\t", col.names = NA,  row.names = T )
  return(results)
}




