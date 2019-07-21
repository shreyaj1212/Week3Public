## Function: load the data

loadFile <- function(disease, chNum, filepath="") {
  # build filename
  filename <- str_glue("Affx_gt_{disease}_Chiamo_{chNum}.tped.gz")
  
  # load file
  if (filepath=="") {
    filepath <- getwd()
  }
  fData <- read_delim(str_glue("{filepath}/{disease}/{filename}"), delim="\t",
                    col_names=F);
#  fData <- read_tsv (str_glue("{filepath}/{disease}/{filename}"), col_names = F)
  
  # get the snpinfo
  # assumes that the snpinfo in a sub-directory snps from the curdir
  return (fData)
}

# load snpinfo file
# filename:snps_chNum
loadSnpInfoFile <- function(chNum, filepath="") {
  snpInfo <- read_tsv(str_glue("{filepath}/snps_{chNum}"), col_names=F)
  
  # only interested in the last 2 columns
  snpInfo <- snpInfo %>% select (snpId=X4, rsid=X5)
  return (snpInfo)
}

#' Title computeAlleleFrequency
#'
#' @param one SNP for all of the individuals
#'
#' @return allele frequency
#' @export
#'
#' @examples 
computeAlleleFrequency <- function(snp) {
  alleles <- unlist(map(snp,function(genotype) {strsplit(genotype, " ")} ))
  t <- table(alleles)
  
  ## get the minor and major alleles
  ## minor and major alleles 
  minorAllele = which.min(t)
  majorAllele <- which.max(t)
  
  maf <- 0
  
  maf <- ifelse (unname(minorAllele) == 1, 
                 unname(t[1]/(sum(t[1], t[2]))),
                 unname(t[2]/(sum(t[1], t[2])))
          )
  
  return(list(t, maf))
}

# HWE
hweChisq <- function (p, o) {
  # o = observed genotype
  # p = the maf for the control group
  t <- sum(o)
  e <- cbind( round(((p^2) * t)), round((2*p*(1-p)*t)), round((((1-p)^2)*t)))
  z <- map2(o,e, function (obs,expd) { return(((obs-expd)^2)/expd)})
  return(pchisq(sum(unlist(z)), df = 1, lower.tail = F))
}

homozygousOrMono <- function (geno) {
  return (ifelse (length(geno) != 3, TRUE, FALSE))
}


gwas <- function (alC, alD, hdC, hdD) {
  if (!hdC["snpId"] == hdD["snpId"])
    return()  

  genotype.C <- table(as.character(alC))
  genotype.D <- table(as.character(alD))

  results <-  data.frame(SNP = as.character(),
                         rsid = as.character (),
                         chrom_num = as.integer(),
                         min_al = as.character(),
                         maj_al = as.character(),
                         MAFD = as.numeric(),
                         MAFC = as.numeric(),
                         or = as.numeric(),
                         or_pval = as.numeric(),
                         or_chi2 = as.numeric(),
                         hw_pval = as.numeric()) 
  
 
  # filtering out data that does not "conform"
  # homozygous, mono-allelic
  # if chr num is missing
  # or if ID is missing
  if ((!is.na(hdC$chr_num)) & (!is.na(hdD$chr_num))) {
    if ((hdD$chr_num != 0) &  (hdC$chr_num != 0) &  
        (!homozygousOrMono(genotype.C)) &  (!homozygousOrMono (genotype.D))) {
      
      # frequencies etc
      cf <- computeAlleleFrequency(alC)
      df <- computeAlleleFrequency(alD)
      tab <- bind_rows(df[[1]], cf[[1]])
      row.names(tab) <- c("Disease", "Control")
     
      ## get the minor allele for disease
      minorAllele <- which.min(df[[1]])
      
      # The first allele is the minor allele
      if (unname(minorAllele) == 1){ 
        or <- unname((tab[1,1]*tab[2,2]) / (tab[1,2] * tab[2,1]))
      } else{ # The second allele is the minor allele.
        or <- unname((tab[1,2]*tab[2,1]) / (tab[1,1]*tab[2,2]))
      }
     
      chiTest <- chisq.test(tab, correct=F)
      or_pval <-chiTest$p.value 
      chi2 <- chiTest$statistic
      
      # calculate the pvalue of HW equilibrium
      hwe <- as.numeric(hweChisq(cf[[2]], genotype.C))
     
      results<- results %>% add_row(SNP = as.character(hdC["snpId"]),
                            rsid = as.character(hdC["rsid"]),
                            chrom_num = as.integer(hdC["chr_num"]),
                            min_al = as.character(names(minorAllele)),
                            maj_al = as.character(names(which.max(cf[[1]]))),
                            MAFD = round(as.numeric(df[[2]]), 2),
                            MAFC = round(as.numeric(cf[[2]]), 2),
                            or = as.numeric(or),
                            or_pval = as.numeric(or_pval),
                            or_chi2 = as.numeric (chi2),
                            hw_pval = as.numeric(hwe)) 
      
      
      return(results)
    }
  }
}

gif <- function (fData) {
  lambda = median(fData$CHI2)/qchisq(0.5, 1)
  return(lambda)
}


filterData <- function (file) {
  df <- read_tsv(file)
  ## filter those that deviate from HWE
  ## filter those that have MAFC < 1%
  ## filter those that are not significant
  
  # Identify significant snps
  snps_significant = df %>% 
                     filter(P <= 1e-7 & MAFC >= 0.01 & HWE >= 0.05) 
  return (snps_significant)
}

myqqplot <- function (p) {
 
  p = case_when(p > 30 ~ 30, p<=30 ~ p)
  pe = sort(-log10(seq(1, 1/length(p), length.out = length(p))))
  plot(pe, p, xlim=c(0,5), ylim=c(0,30), 
       ylab = "Observed -log10(p-value) (values above 30 set at 30)",
       xlab="Expected -log10(p-value)", 
       main = "QQplot T2D")
}

manhattanPlot <- function (fData, x, y, p_threshold = 0.05) {
  p_bonferroni <- p_threshold/nrow(fData)

  fData <- fData %>%
           filter(MAFC > 0.01 & HWE >= 0.05)
  
  f1 = fData %>%
       arrange(CHR, SNP) 
  
  axisdf = f1 %>% 
           group_by(CHR) %>%
           summarize(center = (max(ID) + min(ID))/2)
  
  
  ggplot(f1, aes_string(x=x, y=str_glue("-log10({y})"))) +
    # Show all points
    geom_point( aes(color=as.factor(CHR)), alpha=0.7, size=0.5) +
    scale_color_manual(values = rep(c("lightblue", "darkblue"), 22 )) +
    
    # custom X axis:
    scale_x_continuous(label = axisdf$CHR, breaks = axisdf$center ) +
    scale_y_continuous(expand = c(0, 0) ) +
    
    # Y limit
    coord_cartesian(ylim = c(0, 20))+
    
    geom_hline(yintercept = -log10(p_bonferroni), size=0.5) +
    
    theme_bw() +
    theme( 
      legend.position="none",
      #  panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ) +
    labs(title="Manhattan plot WTCCC T2D", x="Chromosomes", y="-log10(P)")
}

#merge files
smooshFiles <- function (filepath, disease) {
  filenames=list.files(path=filepath, pattern = disease, full.names=TRUE)
  datalist = lapply(filenames, function(x) {read_tsv (file=x)})
  Reduce(function(x,y) {rbind(x,y)}, datalist)
}