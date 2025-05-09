library(data.table)
library(dplyr)
library(stats)
library(optparse)
library(ggplot2)
library(patchwork)
library(rlist)
library(readr)


option_list <- list(
  make_option(c('--GT_file')     ,action='store',type='character',default=NULL,help='Genotypes obtained from vcftools'),
  make_option(c('--AD_file')     ,action='store',type='character',default=NULL,help='Allelic depths obtained from vcftools'),
  make_option(c('--remove')      ,action='store',type='character',default="./no.txt",help='...'),
  make_option(c('--conf.level')  ,action='store',type='numeric'  ,default=0.95,help='...'),
  make_option(c('--probability') ,action='store',type='numeric'  ,default=0.50,help='...'),
  make_option(c('--plots')       ,action='store',type='logical'  ,default=FALSE,help='set to TRUE for creating summary plots'),
  make_option(c('--by_chrom')    ,action='store',type='logical'  ,default=FALSE,help='if TRUE, calculation will be done per chromosome to prevent memory issues'),
  make_option(c('--out_file')    ,action='store',type='character',default=NULL,help='name output file, without extension'))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser) 

gt_file   <- opt$GT_file
ad_file   <- opt$AD_file
out       <- opt$out_file
rm        <- opt$remove
conf      <- opt$conf.level
prb       <- opt$probability
by_chrom  <- opt$by_chrom

###
#gt_file   <- "C:/Users/oostinto/Desktop/raapoi/LG2_snapper_382_tmp.GT.FORMAT"
#ad_file   <- "C:/Users/oostinto/Desktop/raapoi/LG2_snapper_382_tmp.AD.FORMAT"
#out       <- "C:/Users/oostinto/Desktop/raapoi/test_AB"
#rm        <- NULL
#conf      <- 0.99
#prb       <- 0.5
#by_chrom  <- TRUE
###

#read in data
f <- function(x, pos) subset(x, CHROM != "XYZ")
gt <- read_tsv_chunked(gt_file, chunk_size = 100000, callback =  DataFrameCallback$new(f),col_names = TRUE, progress = TRUE )
ad <- read_tsv_chunked(ad_file, chunk_size = 100000, callback =  DataFrameCallback$new(f),col_names = TRUE, progress = TRUE )

plots_dir <- paste0(dirname(out),"/plots")
dir.create(plots_dir)
plots_out <- paste0(plots_dir,"/",basename(out))

#remove individuals from analyses
#only done when file is provided
if(file.exists(rm)){
  print("removing individuals")
  rm_df  <- fread(rm, header = FALSE)
  rm_IND <- rm_df$V1
  rm_IND <- intersect(rm_IND,colnames(gt))
  gt <- gt %>% dplyr::select(-c(all_of(rm_IND)))
  ad <- ad %>% dplyr::select(-c(all_of(rm_IND)))
}

IND   <- colnames(gt)[-1:-2]
N_IND <- length(IND)
N     <- nrow(gt)
signf <- 1-conf
paste0("performing analyses on ",N_IND," individuals")


if(isTRUE(by_chrom)){
    print("performing AB analysis per LG")
    AB_df_ls <- list()
    AB_exclude_ls <- list()
    
    LGs <- unique(gt$CHROM)
    for(LG in LGs){
      
    gt_LG <- gt %>% filter(CHROM == LG) %>% as.matrix()
    ad_LG <- ad %>% filter(CHROM == LG) %>% as.matrix()
    #replace possible phased genotype notations
    gt_LG <-  gsub("\\|","\\/",gt_LG)
    
    G1_LG   <- which(gt_LG == "0/0")
    G2_LG   <- which(gt_LG == "1/1")
    G3_LG   <- which(gt_LG == "./.")
    Gtot_LG <- sort(c(G1_LG,G2_LG,G3_LG))
    #set allelic depths of non-heterozygote sites to NA
    ad_LG[Gtot_LG] <- NA
    
    #create matrixis for AD for ref and alt & remove CHR and POS column
    ad_ref_LG <- matrix(as.numeric(gsub(",\\d*","",ad_LG)[,-1:-2]),ncol = N_IND)
    ad_alt_LG <- matrix(as.numeric(gsub("\\d*,","",ad_LG)[,-1:-2]),ncol = N_IND)
     
    # create data frame with allelic depth of the ref and alt allele over all samples with heterozygote genotypes
    AB_df_LG <-  data.frame(CHROM        = gt_LG[,1],
                            POS          = gt_LG[,2],
                            sumRef       = rowSums(ad_ref_LG , na.rm = TRUE),
                            sumAlt       = rowSums(ad_alt_LG , na.rm = TRUE),
                            meanRef      = rowMeans(ad_ref_LG, na.rm = TRUE),
                            meanAlt      = rowMeans(ad_alt_LG, na.rm = TRUE)) %>%
                     mutate(meanDif      = (sumRef-sumAlt)/2,
                            meanDif_abs  = abs(meanDif),
                            meanFrac_ref = sumRef/(sumRef+sumAlt),
                            N            = sumRef/meanRef)
    
    #perform binomial test to test for a significant deviation from equal proportions
    binom   <- function(x,y) stats::binom.test(c(x,y), p = prb, alternative='two.sided', conf.level = conf )$p.value  
    p.value_LG <- apply(AB_df_LG[,c('sumRef','sumAlt')], 1, function(z) binom(z[1],z[2]) )
    AB_df_LG$binomp <- p.value_LG
    #estimate p.value using a false discovery rate
    AB_df_LG$binompFDR <- p.adjust(AB_df_LG$binomp,method = "fdr")
    AB_df_LG <- AB_df_LG %>% mutate(signf = if_else(binompFDR < signf, TRUE, FALSE))
    AB_exclude_LG <- AB_df_LG %>% dplyr::filter(binompFDR < signf) %>% dplyr::select(CHROM,POS)
    #add results from LG to list
    AB_df_ls[[LG]] <- AB_df_LG 
    AB_exclude_ls[[LG]] <- AB_exclude_LG 
    }
    #bind dfs in the list
    AB_df <- rlist::list.rbind(AB_df_ls)
    AB_exclude <- rlist::list.rbind(AB_exclude_ls)
  
} else {
    print("performing AB analysis over entire dataset")
    ad <- as.matrix(ad)
    gt <- as.matrix(gt)
    #replace possible phased genotype notations
    gt <-  gsub("\\|","\\/",gt)
    
    #extract locations non-heterozygote sites
    G1   <- which(gt == "0/0")
    G2   <- which(gt == "1/1")
    G3   <- which(gt == "./.")
    Gtot <- sort(c(G1,G2,G3))
    #set allelic depths of non-heterozygote sites to NA
    ad[Gtot] <- NA
    
    #create matrixis for AD for ref and alt & remove CHR and POS column
    ad_ref <- matrix(as.numeric(gsub(",\\d*","",ad)[,-1:-2]),ncol = N_IND)
    ad_alt <- matrix(as.numeric(gsub("\\d*,","",ad)[,-1:-2]),ncol = N_IND)
    
    # create data frame with allelic depth of the ref and alt allele over all samples with heterozygote genotypes
    AB_df <-  data.frame(CHROM        = gt[,1],
                         POS          = gt[,2],
                         sumRef       = rowSums(ad_ref , na.rm = TRUE),
                         sumAlt       = rowSums(ad_alt , na.rm = TRUE),
                         meanRef      = rowMeans(ad_ref, na.rm = TRUE),
                         meanAlt      = rowMeans(ad_alt, na.rm = TRUE)) %>%
                  mutate(meanDif      = (sumRef-sumAlt)/2,
                         meanDif_abs  = abs(meanDif),
                         meanFrac_ref = sumRef/(sumRef+sumAlt),
                         N            = sumRef/meanRef)
    
    #perform binomial test to test for a significant deviation from equal proportions
    binom   <- function(x,y) stats::binom.test(c(x,y), p = prb, alternative='two.sided', conf.level = conf )$p.value  
    p.value <- apply(AB_df[,c('sumRef','sumAlt')], 1, function(z) binom(z[1],z[2]) )
    AB_df$binomp <- p.value
    #estimate p.value using a false discovery rate
    AB_df$binompFDR <- p.adjust(AB_df$binomp,method = "fdr")
    AB_df <- AB_df %>% mutate(signf = if_else(binompFDR < signf, TRUE, FALSE))
    AB_exclude <- AB_df %>% dplyr::filter(binompFDR < signf) %>% dplyr::select(CHROM,POS)
}

#write output to table
write.table(AB_df     , file = paste0(plots_out,"_SNP_info_pval",signf,".tsv"), row.names=FALSE,  quote=FALSE , sep = "\t")
write.table(AB_exclude, file = paste0(out,".exclude_pval",signf,".list"), row.names=FALSE, col.names = FALSE, quote=FALSE , sep = "\t")

# if set to true
#create plots
if(opt$plots){
  #stats by SNP
  AB_df_signf  <- AB_df %>% filter(signf)
  #stats by sample
  AD_sample_mat <-  t(apply(ad_ref,1, function(x) x/mean(x, na.rm = TRUE)))
  AD_sample_df  <- data.frame(IND  = IND,
                              diff = colMeans(AD_sample_mat, na.rm = TRUE))
  AD_sample_df_bias <- AD_sample_df %>% filter(diff >= 1.5)
  N_bias <- nrow(AD_sample_df_bias)
  IND_bias <- AD_sample_df_bias$IND
  
  
  p1 <-ggplot(AB_df,aes(x=meanFrac_ref,y=binompFDR, color = N))+
    geom_point()+
    geom_hline(yintercept= signf, color= "red")+
    ylab("p-value")+
    xlab("Proportion ref allele")+
    theme_bw()
  
  p2 <-ggplot(AB_df,aes(x=meanFrac_ref,y=binompFDR, color = N))+
    geom_point()+
    geom_hline(yintercept= signf, color= "red")+
    ylab("p-value")+
    xlab("Proportion ref allele")+
    coord_cartesian(xlim = c(0.3,0.7), ylim = c(0,0.06))+
    theme_bw()
  
  p3 <- ggplot(AB_df,aes(x=meanDif,y=binompFDR, color = N))+
    geom_point()+
    geom_hline(yintercept= signf, color= "red")+
    ylab("p-value")+
    xlab("N AD deviation")+
    theme_bw()
  
  p4 <- ggplot(AB_df,aes(x=meanDif,y=binompFDR, color = N))+
    geom_point()+
    geom_hline(yintercept= signf, color= "red")+
    ylab("p-value")+
    xlab("N AD deviation")+
    coord_cartesian(xlim = c(-100,100), ylim = c(0,0.06))+
    theme_bw()
	
 print("making plot 1") 
 t1 <- (p1+p2)/(p3+p4) + 
        plot_layout(guides = "collect")+
        plot_annotation(title = paste0("Allelic imbalance detected for ",nrow(AB_df_signf)," SNPs out of ",N),
                        subtitle = paste0("SNPs identified using using a significance value of <",signf))
        ggsave(plot = t1, filename = paste0(plots_out,"_stats_SNPs_pval",signf,".png"), dpi = 300, width = 6, height = 6, units = "in")
  
  print("making plot 2") 
  p5 <- ggplot(AD_sample_df,aes(x=IND,y=diff))+
    geom_point(color = "red")+
    geom_text(data=AD_sample_df_bias,aes(x=IND,y=diff, label = IND))+
    ggtitle(label = paste0(N_bias ," sample(s) detected that could cause a bias"), subtitle = "Consider removing for AB analyses ")+
    xlab("Individuals")+
    ylab("relative AD from mean")+
    theme(axis.text.x = element_blank())
  ggsave(plot = p5, filename = paste0(plots_out,"_stats_samples_pval",signf,".png"), dpi = 300, width = 8, height = 4, units = "in")
  
  p6 <- ggplot(AB_df,aes(x=meanFrac_ref,y=N, color = signf))+
    geom_point()+
    theme(legend.position = c(0.15, 0.9))
  
  d1 <- ggplot(AB_df,aes(x=meanFrac_ref, fill = signf))+
    geom_density(alpha = 0.4) + 
    theme_void() + 
    theme(legend.position = "none")
  
  d2 <- ggplot(AB_df,aes(x=N, fill = signf))+
    geom_density(alpha = 0.4) + 
    theme_void() + 
    theme(legend.position = "none")+ 
    coord_flip()
  
  print("making plot 3") 
  t2 <- d1 + plot_spacer() + p6 + d2 + 
        plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))
        ggsave(plot = t2, filename = paste0(plots_out,"_stats_2d_density_pval",signf,".png"), dpi = 300, width = 6, height = 6,  units = "in")
}


