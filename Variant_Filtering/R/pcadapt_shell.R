#This script was written by Tom Oosting and is intended to identify outlier loci using pcadapt which can be run op a HPC and incorporated into BGS pipelines

#Import Packages
library(tidyverse)
library(readr)
library(optparse)
library(pcadapt)
library(qvalue)
library(glue)
library(runner)
library(rlist)

#import parameters
option_list <- list(
  make_option(c('--plink')  ,action='store',type='character',default=NULL  ,help='path to plink bed file without the .bed extension, also expects fam and bim file with the same filename'),
  make_option(c('--out')    ,action='store',type='character',default=NULL  ,help='path and file extension used for all output'),
  make_option(c('--K0')     ,action='store',type='numeric'  ,default=10    ,help='number K PCs to determine the number of informative PCs'),
  make_option(c('--K')      ,action='store',type='numeric'  ,default=1     ,help='Number of PCs that are used to identify outliers'),
  make_option(c('--maf')    ,action='store',type='numeric'  ,default=0.05  ,help='minimaf allele frequency applied'),
  make_option(c('--q')      ,action='store',type='numeric'  ,default=0.05  ,help='q-value cutoff used to identify outliers '),
  make_option(c('--slw')    ,action='store',type='numeric'  ,default=50000 ,help='q-value cutoff used to identify outliers '),
  make_option(c('--minNsnp'),action='store',type='numeric'  ,default=5     ,help='q-value cutoff used to identify outliers '),
  make_option(c('--mode')   ,action='store',type='character',default="full",help='which part of the analyses should be run, options are full (default), scree, or outlier'))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser) 

bed_ext <- opt$plink
out_ext <- opt$out
K0      <- opt$K0
K       <- opt$K
maf     <- opt$maf
q       <- opt$q
mode    <- opt$mode
slw     <- opt$slw
minNsnp <- opt$minNsnp

### functions
chr_info <- function(SNP_df = NULL){chr_info <- SNP_df                   %>% 
                     group_by(CHR)                                       %>% 
                     summarise(Length = max(POS))                        %>% 
                     arrange(factor(CHR, levels = unique(SNP_df$CHR)))   %>% 
                     mutate(tot = cumsum(Length)-Length, 
                            LG  = c(1:length(unique(SNP_df$CHR))))
return(chr_info)}

#load plink data
bed <- read.pcadapt(paste0(bed_ext,".bed"), type = "bed")
fam <- read_tsv(paste0(bed_ext,".fam"),col_names = FALSE)
bim <- read_tsv(paste0(bed_ext,".bim"),col_names = FALSE)
colnames(bim) <- c("CHR","LOC","Cmn","POS","REF","ALT")

#get chromosome information 
#requires columns CHR and POS
chr_info <- chr_info(bim)    

#add chr info to SNPinfo
#cumulative base pair (BPcum) is added for manhattan plotting
bim <- left_join(bim,chr_info[,c("CHR","LG","tot")]) %>% 
  arrange(LG,POS)                               %>% 
  mutate(BPcum = tot + POS)

#pcadapt - stage 1
#run pcadapt with high K (e.g. number of sample locations +5)
#determine right number of PC to retain by applying Cattell’s rule:
#retain number of PCs that correspond to eigenvalues to the left of the straight line.
if(mode == "scree" | mode == "full"){
  pcadapt <- pcadapt(input = bed, K = K0)
  
  png(glue("{out_ext}_screeplot.png"),res = 300, units = "in", width = 10, height = 6)
  plot(pcadapt, option = "screeplot")
  dev.off()
  png(glue("{out_ext}_scores.png"),res = 300, units = "in", width = 10, height = 6)
  plot(pcadapt, option = "scores", i = 1, j =2, pop = fam$X1)
  dev.off()
}


#pcadapt - stage 2
if(mode == "outlier" | mode == "full"){
  #run pcadapt for K PCs
  pcadapt_K <- pcadapt(input = bed, K = K, min.maf = maf)
  
  #summary plots
  png(glue("{out_ext}_K{K}_scores.png"),res = 300, units = "in", width = 7, height = 6)
  plot(pcadapt_K, option = "scores", pop = fam$X1)
  dev.off()
  
  png(glue("{out_ext}_K{K}_QQplot.png"),res = 300, units = "in", width = 7, height = 6)
  plot(pcadapt_K, option = "qqplot")
  dev.off()
  
  png(glue("{out_ext}_K{K}_pval_hist.png"),res = 300, units = "in", width = 10, height = 6)
  hist(pcadapt_K$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
  dev.off()
  
  png(glue("{out_ext}_K{K}_stat.distribution.png"),res = 300, units = "in", width = 7, height = 6)
  plot(pcadapt_K, option = "stat.distribution")
  dev.off()
  
  #obtain q-values
  qval  <- qvalue(pcadapt_K$pvalues)
  bim$p <- pcadapt_K$pvalues
  bim$q <- qval$qvalues
  
  #print summary
  summary(qval)
  capture.output(summary(qval), file = glue("{out_ext}_qvalue_sum.txt"))
  
  png(glue("{out_ext}_K{K}_q{q}_qvalue_sum.png"),res = 300, units = "in", width = 10, height = 6)
  plot(qval)
  dev.off()
  
  png(glue("{out_ext}_K{K}_q{q}_qvalue_hist.png"),res = 300, units = "in", width = 10, height = 6)
  hist(qval)
  dev.off()
  
  #select all significant outliers
  bim_outliers <- bim[which(qval$qvalues < q),]
  if(nrow(bim_outliers) != 0){
    write_tsv(bim_outliers, file = glue("{out_ext}_K{K}_q{q}_all_outliers.tsv"), col_names = TRUE)
  # filter for independent loci
  # LGs containing outliers
  LGs <- unique(bim_outliers$CHR)
  #create emtpy list
  runner_list <- list()
  for(i in LGs){
    LG_bim <- bim_outliers %>% filter(CHR == i)  
    
    
    runner_LOC <- runner(x   = LG_bim
                        ,idx = "BPcum"  
                        ,f   = function(x){x[["LOC"]][[which.min(x[["q"]])]] }     
                        ,k   = slw
                        ,at  = seq(from =1, to = max(LG_bim$BPcum+slw), by = slw)
                        ,na_pad = FALSE)
    
    runner_N   <-   runner(x   = LG_bim
                          ,idx = "BPcum"  
                          ,f   = function(x){nrow(x)}     
                          ,k   = slw
                          ,at  = seq(from =1, to = max(LG_bim$BPcum+slw), by = slw)
                          ,na_pad = FALSE)
    
    runner_df <- data.frame(LOC = runner_LOC,
                            N   = runner_N) %>% 
      filter(!is.na(LOC))        %>%  #filter out sliding windows without outliers
      filter(N >= minNsnp) #filter out outliers from sliding window that do not have the minimum number of SNPs
    
    if(nrow(runner_df != 0)){
      LG_bim_selected <- inner_join(LG_bim,runner_df, by = c("LOC"="LOC"))
      runner_list[[paste0(i,"_runner")]] <- LG_bim_selected
    } 
  }
  
  bim_independent  <- as.data.frame(rlist::list.rbind(runner_list))
  write_tsv(bim_independent, file = glue("{out_ext}_K{K}_q{q}_independent_outliers.tsv"), col_names = TRUE)
  
  #Manhattan plot
  png(glue("{out_ext}_K{K}_q{q}_Manhattan.png"),res = 300, units = "in", width = 14, height = 6)
  plot(x    = bim$BPcum, y= -log10(bim$q),
       xlab = "base pairs (cumulative)",
       ylab = "log10(q-value)",
       main = glue("Manhattan plot showing outlier loci, {nrow(bim_outliers)} significant outliers (blue), {nrow(bim_independent)} independent significant outliers (green)"))
  points(x = bim_outliers$BPcum,
         y = -log10(bim_outliers$q), 
         col = "blue")
  points(x = bim_independent$BPcum, 
         y = -log10(bim_independent$q), 
         col = "green")
  dev.off()
  } else {
    print("no outliers were found")
  }
}



