#functions are written or modified by T. Oosting
#If you use these functions please cite Oosting et al., 2023, and if you use the PCA, PCoA, and DAPC functions also cite Therkildsen et al., 2019.
#packages
library(dplyr)
library(stringr)
library(ggplot2)
library(raster)
library(gdistance)
library(reshape2)
library(boot)
library(SNPRelate)
library(vcfR)

#obtain vector for ordering individuals and group by factor
#e.g. sample_info$ORD <- order_INDS(IND = sample_info$IND, FACTOR = sample_info$POP)
order_INDS <- function(IND = NULL, FACTOR = NULL ){
  tmp_df <- data.frame(IND = IND,
                       FAC = FACTOR)
  
  ORD <- tmp_df %>% arrange(FAC, .by_group = TRUE) %>% dplyr::select(IND)
  ORD2 <- order(tmp_df$FAC, ORD$IND)
  
  tmp_df2 <- data.frame(IND = 1:nrow(ORD),
                        ORD = ORD2) %>% 
                        arrange(ORD)
  return(tmp_df2$IND)
}


#Rfunctions
#####reading in vcf or gds file into R (SNPrelate)#####
snpgdsReadGDS <- function(vcf_file = NULL, gds_file = NULL){
  if(is.null(gds_file)){
    print("no gds file suplied, looking for gds file with similar name as vcf file")
    gds_file <- str_replace(vcf_file,".vcf.gz",".gds")
    if(file.exists(gds_file)){
      print("gds file found, loading now")
      gds <- SNPRelate::snpgdsOpen(gds_file,readonly = FALSE)
      return(gds)
    } else {
      print("no gds file, converting vcf to gds in same directory")
      SNPRelate::snpgdsVCF2GDS(vcf_file, gds_file, method="biallelic.only")
      gds <- SNPRelate::snpgdsOpen(gds_file,readonly = FALSE)
      gdsfmt::add.gdsn(gds,"snp.id", paste0(read.gdsn(index.gdsn(gds, "snp.chromosome")),
                                            ":",
                                            read.gdsn(index.gdsn(gds, "snp.position"))) ,
                       replace = TRUE)
      return(gds)
    } 
  } else {
    print("gds file supplied")
    if(file.exists(gds_file)){
      print("gds file found, loading now")
      gds <- SNPRelate::snpgdsOpen(gds_file,readonly = FALSE)
      return(gds)
    } else {
      print("no gds file, converting vcf to gds if vcf is suplied")
      SNPRelate::snpgdsVCF2GDS(vcf_file, gds_file, method="biallelic.only")
      gds <- SNPRelate::snpgdsOpen(gds_file,readonly = FALSE)
      gdsfmt::add.gdsn(gds,"snp.id", paste0(read.gdsn(index.gdsn(gds, "snp.chromosome")),
                                            ":",
                                            read.gdsn(index.gdsn(gds, "snp.position"))) ,
                       replace = TRUE)
      return(gds)
    }
  }
}

#####getting SNP information from gds file (SNPrelate)#####

#Obtain level of heterozygosity per SNP
snpgdsSNPHet <- function(gds       = NULL,
                         sample.id = NULL,
                         snp.id    = NULL,
						             margin    = 2 ){apply(SNPRelate::snpgdsGetGeno(gdsobj = gds,
                                                                        sample.id = sample.id,
                                                                        snp.id = snp.id),
                                                                        margin,
                                                                        #verbose = FALSE,
                                                                        function(y){length(which(y==1))/length(which(!is.na(y)))})}

#Get reference allele count (AC) per SNP
snpgdsSNPACRef <- function(gds       = NULL,
                           sample.id = NULL,
                           snp.id    = NULL,
	   					             margin    = 2 ){apply(SNPRelate::snpgdsGetGeno(gdsobj = gds,
                                                                          sample.id = sample.id,
                                                                          snp.id = snp.id),
                                                                          margin,
                                                                          #verbose = FALSE,
                                                                          function(y){(length(which(y==2))*2)+length(which(y==1))})}

#Get frequency alternative allele
snpgdsSNPACAlt <- function(gds       = NULL,
                           sample.id = NULL,
                           snp.id    = NULL,
						               margin    = 2 ){apply(SNPRelate::snpgdsGetGeno(gdsobj = gds,
                                                                          sample.id = sample.id,
                                                                          snp.id = snp.id),
                                                                          margin,
                                                                          #verbose = FALSE,
                                                                          function(y){(length(which(y==0))*2)+length(which(y==1))})}

#Get alternative allele count (AC) per SNP
snpgdsSNPRef_freq <- function(gds    = NULL,
                           sample.id = NULL,
                           snp.id    = NULL,
						               margin    = 2 ){apply(SNPRelate::snpgdsGetGeno(gdsobj = gds,
                                                                          sample.id = sample.id,
                                                                          snp.id = snp.id),
                                                                          margin,
                                                                          #verbose = FALSE,
                                                                          function(y){((length(which(y==2))*2)+length(which(y==1)))/(length(which(!is.na(y)))*2)})}

#Get frequency reference allele
snpgdsSNPAlt_freq <- function(gds    = NULL,
                           sample.id = NULL,
                           snp.id    = NULL,
						               margin    = 2 ){apply(SNPRelate::snpgdsGetGeno(gdsobj = gds,
                                                                          sample.id = sample.id,
                                                                          snp.id = snp.id),
                                                                          margin,
                                                                          #verbose = FALSE,
                                                                          function(y){((length(which(y==0))*2)+length(which(y==1)))/(length(which(!is.na(y)))*2)})}

#Get frequency missing data
snpgdsSNPMissing <- function(gds       = NULL,
                             sample.id = NULL,
                             snp.id    = NULL,
						                 margin    = 2 ){apply(SNPRelate::snpgdsGetGeno(gdsobj = gds,
                                                                            sample.id = sample.id,
                                                                            snp.id = snp.id),
                                                                            margin,
                                                                            #verbose = FALSE,
                                                                            function(y){length(which(is.na(y)))/length(y)})}

#Obtain fraction of homozygote genotypes of reference allele per SNP
snpgdsSNPhom_ref <- function(gds       = NULL,
                             sample.id = NULL,
                             snp.id    = NULL,
                             margin    = 2 ){apply(SNPRelate::snpgdsGetGeno(gdsobj = gds,
                                                                            sample.id = sample.id,
                                                                            snp.id = snp.id),
                                                                            margin,
                                                                            #verbose = FALSE,
                                                                            function(y){length(which(y==0))/length(which(!is.na(y)))})}

#Obtain fraction of homozygote genotypes of alternative allele per SNP
snpgdsSNPhom_alt <- function(gds       = NULL,
                             sample.id = NULL,
                             snp.id    = NULL,
                             margin    = 2 ){apply(SNPRelate::snpgdsGetGeno(gdsobj = gds,
                                                                            sample.id = sample.id,
                                                                            snp.id = snp.id),
                                                                            margin,
                                                                            #verbose = FALSE,
                                                                            function(y){length(which(y==2))/length(which(!is.na(y)))})}



#SNP summary table from gds SNPRelate object
#function returns data frame with following statistics:
# Chromosome                             (CHR)
# Linkage group                          (LG) i.e. numerical variant to CHR
# snp.id                                 (LOC)
# snp.position                           (POS)
# Reference allele                       (REF)
# Alternative allele                     (ALT)
# fraction of missing data               (missing)
# minimum allele frequency               (MAF)
# heterozygosity                         (HET)
# Allele count reference allele          (AC_ref)
# Allele count alternative allele        (AC_alt)
# Frequency reference allele             (REF_freq)
# Frequency alternative allele           (ALT_freq)
# p-value for Hardy-Weinberg equilibrium (HWE_pval)
snpgdsSNPsum <- function(gds_snprelate = NULL,
                         sample.id     = NULL,
                         snp.id        = NULL,
                         extended      = FALSE){SNP_info <- data.frame(CHR = read.gdsn(index.gdsn(gds_snprelate, "snp.chromosome")),
                                                                       LG  = as.numeric( str_remove_all(read.gdsn(index.gdsn(gds_snprelate, "snp.chromosome")), "[:alpha:]|[:punct:]" )),
                                                                       LOC = paste0(read.gdsn(index.gdsn(gds_snprelate, "snp.chromosome")),
                                                                                    ":",
                                                                                    read.gdsn(index.gdsn(gds_snprelate, "snp.position"))),
                                                                       POS = read.gdsn(index.gdsn(gds_snprelate, "snp.position")),
                                                                       REF = stringr::str_extract(read.gdsn(index.gdsn(gds_snprelate, "snp.allele")),"^\\w"),
                                                                       ALT = stringr::str_extract(read.gdsn(index.gdsn(gds_snprelate, "snp.allele")),"\\w$"))

                              if(!is.null(snp.id)){SNP_info <- SNP_info[which(SNP_info$LOC %in% snp.id),]}
                              
                              if(extended){
                              SNP_info$MISSING  = SNPRelate::snpgdsSNPRateFreq(gds_snprelate, sample.id = sample.id, snp.id = snp.id)$MissingRate
                              SNP_info$MAF      = SNPRelate::snpgdsSNPRateFreq(gds_snprelate, sample.id = sample.id, snp.id = snp.id)$MinorFreq
                              SNP_info$HET      = snpgdsSNPHet(gds = gds_snprelate, sample.id = sample.id, snp.id = snp.id, margin = 2)
                              SNP_info$HOM_ref  = snpgdsSNPhom_ref(gds = gds_snprelate, sample.id = sample.id, snp.id = snp.id, margin = 2)
                              SNP_info$HOM_alt  = snpgdsSNPhom_alt(gds = gds_snprelate, sample.id = sample.id, snp.id = snp.id, margin = 2)
                              #SNP_info$AC_ref   = snpgdsSNPACRef(gds = gds_snprelate, sample.id = sample.id, snp.id = snp.id, margin = 2)
                              #SNP_info$AC_alt   = snpgdsSNPACAlt(gds = gds_snprelate, sample.id = sample.id, snp.id = snp.id, margin = 2)
                              SNP_info$REF_freq = snpgdsSNPRef_freq(gds_snprelate, sample.id = sample.id, snp.id = snp.id, margin = 2)
                              SNP_info$ALT_freq = snpgdsSNPAlt_freq(gds_snprelate, sample.id = sample.id, snp.id = snp.id, margin = 2)
                              SNP_info$MISS     = snpgdsSNPMissing(gds_snprelate, sample.id = sample.id, snp.id = snp.id, margin = 2)
                              SNP_info$HWE_pval = SNPRelate::snpgdsHWE(gds_snprelate, sample.id = sample.id, snp.id = snp.id)
                              }
                              return(SNP_info)
}

#same as function above but it requires a vector indicating populations
#this function will provide information about all SNPs per population
snpgdsSNPsum_byPOP <- function(gds       = NULL,
                               sample.id = NULL,
                               snp.id    = NULL,
                               pop.id    = NULL){ if(is.null(sample.id)){sample.id <- read.gdsn(index.gdsn(gds, "sample.id"))}
  
                                                  tmp_df1 <- data.frame(matrix(ncol = 16, nrow = 0))
                                                  colnames(tmp_df1) <- c("CHR","LG","LOC","POS","REF","ALT","MISSING","MAF","HET","Hom_ref","Hom_alt","REF_freq","ALT_freq","MISS","HWE_pval","POP")
                                                  
                                                  pops <- unique(pop.id)
                                                  populations <- split(sample.id, pop.id)
                                                  
                                                  for (i in pops) {
                                                    tmp_df2 <- snpgdsSNPsum(gds, sample.id = populations[[i]], snp.id = snp.id, extended = TRUE)
                                                    tmp_df2$POP <- i
                                                    tmp_df1 <- rbind(tmp_df1,tmp_df2)}
                                                  return(tmp_df1)
}

#get either maf or mac per populations, the last column in the returned df indicates the number of populatios that have the alternative allele
#this analyses assumes that the alternative allele is the minar allele, i.e. always uses the ALT column for the vcf/gds file
#this this way maf/mac is always giving for the same allele per POP.
#function needs:
# gds    = gds file
# pop.id = population vector
# mac    = logical operator the return allele counts when true

snpgdsSNPmaf_byPOP <- function(gds       = NULL,
                               mac       = FALSE,
                               pop.id    = NULL,
                               sample.id = NULL){
                               
                               if(is.null(sample.id)){
                                 sample.id <- read.gdsn(index.gdsn(gds, "sample.id"))
                               }
  
                               populations <- split(sample.id, pop.id)
                               
                               N <- lengths(populations)
                               P <- names(populations)
                               Npop <- length(N)
                               
                               tmp_df <- data.frame(LOC = read.gdsn(index.gdsn(gds, "snp.id")))
                               #get frequencies of the reference alleles
                               freq_REF <- snpgdsSNPList(gds)
                               #determine whether the REF or ALT allele has the highest freq
                               # 0 = REF allele has highest freq
                               # 1 = ALT allele has highest freq
                               # this sets the same reference allele for all populations regardless of the frequency
                               major_allele <- rep(0,nrow(freq_REF))
                               major_allele[which(freq_REF$afreq < 0.5)] <- 1
                               
                               major_freq <- abs( major_allele - freq_REF$afreq )
                               
                               for (i in P) {
                                 #if mac is TRUE alternative allele counts are given instead of frequency   
                                 if(!mac){
                                   X1 <- snpgdsSNPList(gds, sample.id = populations[[i]])
                                   X2 <- abs( major_allele - X1$afreq )
                                   X3 <- 1-X2
                                 } else {
                                   X1 <- length(populations[[i]])*2*major_allele
                                   X2 <- snpgdsSNPACAlt(gds = gds, sample.id = populations[[i]])
                                   X3 <- abs(X2-X1)
                                 }
                                 tmp_df <- cbind(tmp_df,X3)
                               }
                               
                               X <- apply(tmp_df[,-1],1,function(z){length( which(z > 0))})
                               tmp_df <- cbind(tmp_df,X)
                               colnames(tmp_df) <- c("LOC",paste0(P,"(",N,")"),"N")
                               
                               return(tmp_df)
}


#####getting sample information from gds file (SNPrelate)#####

#Obtain level of heterozygosity per sample
snpgdsINDHet <- function(gds       = NULL,
                         sample.id = NULL,
                         snp.id    = NULL){apply(SNPRelate::snpgdsGetGeno(gdsobj = gds,
                                                                          sample.id = sample.id,
                                                                          snp.id = snp.id),
                                                                          1,
                                                                          #verbose = FALSE,
                                                                          function(y){length(which(y==1))/length(which(!is.na(y)))})}



snpgdsINDMiss <- function(gds       = NULL,
                          sample.id = NULL,
                          snp.id    = NULL){apply(SNPRelate::snpgdsGetGeno(gdsobj = gds,
                                                                           sample.id = sample.id,
                                                                           snp.id = snp.id),
                                                                           1,
                                                                           #verbose = FALSE,
                                                                           function(y){length(which(is.na(y)))/length(y)})}






#####estimating pairwise FST between sample locations#####

#estimate pwFSTs using SNPrelate
# 2 parametes required: gds & pop.id
# requires gds object generated via SNPrelate (not SeqArray!)
# vector for populations should reflect the order of samples in the gds object!!!
snpgdsPWFST <- function(gds = NULL, pop.id = NULL, snp.id = NULL, sample.id = NULL, out = NULL, reps = 1000, signf = 0.05,
                        exclude.monomorphic = TRUE, autosome.only = FALSE, FST.method = "W&C84", maf = 0.05, missing.rate = 0.95){
  #set vector of samples included FST estimates
  if(is.null(sample.id)){sample.id <- read.gdsn(index.gdsn(gds, "sample.id"))}
  #extract population information 
  pops <- as.character(unique(pop.id))
  N_pops <- length(pops)
  populations <- split(sample.id,  pop.id)
  pop_combs <- as.data.frame(t(utils::combn(pops,2)))
  pop_combs_df <- as.data.frame(t(utils::combn(pops,2)))
  pop_combs_df <- pop_combs %>% mutate(merge = paste(V1,V2,sep = "."))
  N_comb <- nrow(pop_combs)
  pop_df <- data.frame(IND = sample.id,
                       POP = pop.id)
  #print population information
  print("starting pairwise FST estimation", quote = F)
  print("analyses includes:", quote = F)
  print(paste(length(sample.id),"individuals"), quote = F)
  print(paste(N_pops, "populations"), quote = F)
  print(paste(N_comb,"population combinations"), quote = F)
  #create empty data frame to print output
  pwFST_df <- data.frame(matrix(ncol = 16,nrow = 0))
  #loop over population combinations
  for(i in 1:N_comb){
    #get populations
    pop1       <- pop_combs_df$V1[i] 
    pop2       <- pop_combs_df$V2[i]
    indo_2pops <- pop_df %>% dplyr::filter(POP %in% c(pop1, pop2))
    #print information
    print("#######################################################", quote = F)
    print(paste("estimating",i,"out of",N_comb))
    print(paste("population1:",pop1))
    print(paste("population2:",pop2))
    #estimate PWFst 
    fst_2pop <- snpgdsFst(gdsobj         = gds,
                          sample.id      = indo_2pops$IND,
                          population     = factor(indo_2pops$POP,levels = c(pop1,pop2)),
                          method         = FST.method,
                          autosome.only  = autosome.only,
                          remove.monosnp = exclude.monomorphic,
                          maf            = maf,
                          missing.rate   = missing.rate)
    N_SNPs <- length(fst_2pop$FstSNP)
    # bootstrapping for CI over non-weighted FST estimates form all SNPs
    confid <- 1-signf
    # internal function to obtain the mean
    Bmean <- function(data,indices){
      d <- data[indices] # allows boot to select sample 
      return(mean(d))}
    results <- boot(data= fst_2pop$FstSNP, statistic=Bmean, R=reps)
    #create QQplots (only when out is specified)
    if(!is.null(out)){
      dir.create(paste0(dirname(out),"/QQplots"), recursive = TRUE, showWarnings = FALSE)
      png(paste0(dirname(out),"/QQplots/",basename(out),"_qqplot_",pop1,"_",pop2,"_",reps,"bootstrap.png"), width = 14, height = 7, units = "in", res = 300)
      plot(results)   
      dev.off()}
    
    #obtain pval
    Ttest <- t.test(results$t)
    pval  <- Ttest$p.value
    
    #obtain confidence intervals
    CIs <- boot::boot.ci(results, type=c("norm", "basic", "perc"), conf = confid)
    CI_vec <- c(CIs$normal[2],CIs$normal[3],CIs$basic[4],CIs$basic[5],CIs$percent[4],CIs$percent[5],pval)
    
    data <- as.character(c(pop1,pop2,fst_2pop$Fst,summary(fst_2pop$FstSNP),CI_vec,N_SNPs))
    pwFST_df <- pwFST_df %>% mutate_all(as.character)
    pwFST_df <- rbind(pwFST_df,data)
  }
  #assign column names
  colnames(pwFST_df) <- c("pop1", "pop2","FST_weighted", "Min", "1st_Qu", "Median", "Mean", "3rd_Qu", "Max","CI95%_min_norm","CI95%_max_norm","CI95%_min_basic","CI95%_max_basic","CI95%_min_perc","CI95%_max_perc", "pval","N_snps")
  pwFST_df[,3:15] <- sapply(pwFST_df[,3:15],as.numeric)
  pwFST_df[,3:15] <- pwFST_df[,3:15] %>% rowwise() %>% round(digits = 7)
  pwFST_df <- pwFST_df %>% mutate(CI_range = paste(`CI95%_min_norm`,`CI95%_max_norm`, sep = " - "))
  #print pwFST table
  
  if(!is.null(out)){
    write_tsv(pwFST_df,paste0(out,"_pwFST.tsv"))
  }
  return(pwFST_df)
}

#convert output to matrix for plotting with pheatmap
snpgdsFSTdf2mat <- function(pwFST_df=NULL, order=NULL){
                            ### convert pwFST data frame to square matrix via distance matrix
                            ## set up storage matrix
                            # get names for row and columns
                            nameVals <- sort(unique(unlist(pwFST_df[1:2])))
                            # set factor for pops
                            if(!is.null(order)){
                              print("setting matrix to custom order")
                              nameVals <- nameVals[order]
                            }
                            # construct 0 matrix of correct dimensions with row and column names
                            myMat <- matrix(0, length(nameVals), length(nameVals), dimnames = list(nameVals, nameVals))
                            # fill in the matrix with matrix indexing on row and column names
                            myMat[as.matrix(pwFST_df[c("pop1", "pop2")])] <- pwFST_df[["FST_weighted"]]
                            myMat <- myMat + t(myMat)
                            return(myMat)
}


#####converting gds to bed file#####
snpgdsSNPbim <- function(gds_snprelate = NULL,
                         sample.id     = NULL,
                         snp.id        = NULL){SNP_info <- data.frame(chr = read.gdsn(index.gdsn(gds_snprelate, "snp.chromosome")),
                                                                      id  = paste0(read.gdsn(index.gdsn(gds_snprelate, "snp.chromosome")),
                                                                                   ":",
                                                                                   read.gdsn(index.gdsn(gds_snprelate, "snp.position"))),
                                                                      posg  = 0,
                                                                      pos = read.gdsn(index.gdsn(gds_snprelate, "snp.position")),
                                                                      ref = stringr::str_extract(read.gdsn(index.gdsn(gds_snprelate, "snp.allele")),"^\\w"),
                                                                      alt = stringr::str_extract(read.gdsn(index.gdsn(gds_snprelate, "snp.allele")),"\\w$"))
                         return(SNP_info)
}


snpgdsgds2bed <- function(gds = NULL, pop_df = NULL, out = NULL, exclude = NULL ){
  library(genio)
  
  ### get bim file ###
  bim <- snpgdsSNPbim(gds = gds)
  
  #check if CHR is numeric - change if not
  if(is.numeric(class(bim$chr))){
    print("CHR is in numeric format, no modifications made to CHR info field")
  } else {
    print("CHR is non-numeric, will replace CHR names with numeric sequence, original CHR name is still indicated in the variant identifier")
    chr_old <- unique(bim$chr)
    chr_new <- seq(length(chr_old))
    bim$chr <- chr_new[match(bim$chr, chr_old)]
  }
  
  ### get fam file ###
  if(!is.null(pop_df)){
    print("pop_file provided, generating fam file")
    print("looking for IND and POP fields")
    
    #vector if samples in the order they occr in the vcf/gds
    INDS_ordered <- read.gdsn(index.gdsn(gds, "sample.id"))
    
    if(!is.null(exclude)){
      print("sample to exclude were provided, removing now...")
      INDS_ordered <- setdiff(INDS_ordered,exclude)
      length(INDS_ordered)
    }
    
    #read pop file and filter the file for names present in vcf/gds
    pop_df <- pop_df %>% filter(IND %in% INDS_ordered) %>% dplyr::arrange(IND)
    #order samples in the order they occur in the vcf/gds
    pop_df[match(INDS_ordered, pop_df$IND),]
    #create df
    fam <- data.frame(fam   =  pop_df$POP,
                      id    =  pop_df$IND,
                      pat   =  0,
                      mat   =  0,
                      sex   =  0,
                      pheno = -9)
    
    #write_tsv(fam,file = glue("{out_ext}.fam"), col_names = F)  
  } else {
    INDS_ordered <- read.gdsn(index.gdsn(gds, "sample.id"))
    fam <- data.frame(fam   =  INDS_ordered,
                      id    =  INDS_ordered,
                      pat   =  0,
                      mat   =  0,
                      sex   =  0,
                      pheno = -9)
  }
  
  ### get bed file ###
  GT <- snpgdsGetGeno(gds)
  bed  <- t(GT)
  GT_cols <- which(read.gdsn(index.gdsn(gds, "sample.id")) %in% INDS_ordered)
  bed  <- bed[,GT_cols]
  
  ### write files ###
  #write new plink
  genio::write_plink(file = out, X = bed, fam = fam, bim = bim)
  
  
}

#####creating objects for Manhattan plots#####
#collect chromosome information
#date frame requires following columns:
# Chromosome   as "CHR"
# SNP position as "POS"
chr_info <- function(SNP_df = NULL){chr_info <- SNP_df %>% group_by(CHR) %>% summarise(Length = max(POS)) %>% 
  arrange(factor(CHR, levels = unique(SNP_df$CHR)))         %>% 
  mutate(tot = cumsum(Length)-Length, 
         LG  = c(1:length(unique(SNP_df$CHR))))
return(chr_info)
}

#create data frame that contains information for where chromosome start and stop with in plot
#requires: 
# CHR
# LG
# BPcum
axisdf  <- function(SNP_df){ axisdf <- SNP_df %>% group_by(CHR) %>% summarize(LG     =  first(LG),
                                                                              center = (max(BPcum) + min(BPcum) ) / 2, 
                                                                              start  =  min(BPcum),
                                                                              stop   =  max(BPcum),
                                                                              ymin   =  -Inf,
                                                                              ymax   =  Inf) %>% arrange(LG)
return(axisdf)
}
#create x-scale
x_scale <- function(axisdf){ x_scale <- scale_x_continuous(label=axisdf$CHR,breaks = axisdf$center,expand=c(0,0))
return(x_scale)
}


#create rectangles for background in ggplot
bg_rect <- function(axisdf){ bg_rect <- annotate("rect", xmin = axisdf$start,xmax = axisdf$stop, 
                                                 ymin = axisdf$ymin, ymax = axisdf$ymax, 
                                                 alpha = .5, 
                                                 fill = rep(c("darkgrey","white"),(nrow(axisdf)+1)/2)[1:nrow(axisdf)])
return(bg_rect)
}


#####
## This script contains some essential functions for the individual level PCA and DAPC analysis and visualization
# These functions use a covariance matrix as the input
# Individual ID and population labels should also be supplied
# It can be imported to any R script using: source("individual_pca_functions.R")
# Please read each function for further details
######geom_enterotype
geom_enterotype <- function(mapping = NULL, data = NULL, stat = "identity",  position = "identity",
                            alpha = 0.3, prop = 0.5, ..., lineend = "butt", linejoin = "round", 
                            linemitre = 1, arrow = NULL, na.rm = FALSE, parse = FALSE, 
                            nudge_x = 0, nudge_y = 0, label.padding = unit(0.15, "lines"), 
                            label.r = unit(0.15, "lines"), label.size = 0.1, 
                            show.legend = TRUE, inherit.aes = TRUE, show.point=T, show.label=T, show.ellipse=T, show.line=T) {
  ## This function is used to create an enterotype plot from two dimemsional data. Features can be turned on and off using show.point, show.label, show.ellipse, and show.line options. 
  library(ggplot2)
  # create new stat and geom for PCA scatterplot with ellipses
  StatEllipse <- ggproto("StatEllipse", Stat, 
                         required_aes = c("x", "y"), 
                         compute_group = function(., data, scales, level = 0.75, segments = 51, ...) {
                           library(MASS)
                           dfn <- 2
                           dfd <- length(data$x) - 1
                           if (dfd < 3) {
                             ellipse <- rbind(c(NA, NA))
                           } else {
                             v <- cov.trob(cbind(data$x, data$y))
                             shape <- v$cov
                             center <- v$center
                             radius <- sqrt(dfn * qf(level, dfn, dfd))
                             angles <- (0:segments) * 2 * base::pi/segments
                             unit.circle <- cbind(cos(angles), sin(angles))
                             ellipse <- t(center + radius * t(unit.circle %*% chol(shape)))
                           }
                           ellipse <- as.data.frame(ellipse)
                           colnames(ellipse) <- c("x", "y")
                           return(ellipse)
                         })
  
  # write new ggproto 
  GeomEllipse <- ggproto("GeomEllipse", Geom, 
                         draw_group = function(data, panel_scales, coord) {
                           n <- nrow(data)
                           if (n == 1) 
                             return(zeroGrob())
                           munched <- coord_munch(coord, data, panel_scales)
                           munched <- munched[order(munched$group), ]
                           first_idx <- !duplicated(munched$group)
                           first_rows <- munched[first_idx, ]
                           grid::pathGrob(munched$x, munched$y, default.units = "native", 
                                          id = munched$group, 
                                          gp = grid::gpar(col = first_rows$colour, 
                                                          fill = alpha(first_rows$fill, first_rows$alpha), lwd = first_rows$size * .pt, lty = first_rows$linetype))
                         }, 
                         default_aes = aes(colour = "NA", fill = "grey20", size = 0.5, linetype = 1, alpha = NA, prop = 0.5), 
                         handle_na = function(data, params) {
                           data
                         }, 
                         required_aes = c("x", "y"), 
                         draw_key = draw_key_path
  )
  
  # create a new stat for PCA scatterplot with lines which totally directs to the center
  StatConline <- ggproto("StatConline", Stat, 
                         compute_group = function(data, scales) {
                           library(miscTools)
                           library(MASS)
                           df <- data.frame(data$x,data$y)
                           mat <- as.matrix(df)
                           center <- cov.trob(df)$center
                           names(center)<- NULL 
                           mat_insert <- insertRow(mat, 2, center )
                           for(i in 1:nrow(mat)) {
                             mat_insert <- insertRow( mat_insert, 2*i, center )
                             next
                           }
                           mat_insert <- mat_insert[-c(2:3),]
                           rownames(mat_insert) <- NULL
                           mat_insert <- as.data.frame(mat_insert,center)
                           colnames(mat_insert) =c("x","y")
                           return(mat_insert)
                         },
                         required_aes = c("x", "y")
                         
  )
  
  # create a new stat for PCA scatterplot with center labels
  StatLabel <- ggproto("StatLabel" ,Stat,
                       compute_group = function(data, scales) {
                         library(MASS)
                         df <- data.frame(data$x,data$y)
                         center <- cov.trob(df)$center
                         names(center)<- NULL 
                         center <- t(as.data.frame(center))
                         center <- as.data.frame(cbind(center))
                         colnames(center) <- c("x","y")
                         rownames(center) <- NULL
                         return(center)
                       },
                       required_aes = c("x", "y")
  )
  
  
  layer1 <- layer(data = data, mapping = mapping, stat = stat, geom = GeomPoint, 
                  position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
                  params = list(na.rm = na.rm, ...))
  layer2 <- layer(stat = StatEllipse, data = data, mapping = mapping, geom = GeomEllipse, position = position, show.legend = FALSE, 
                  inherit.aes = inherit.aes, params = list(na.rm = na.rm, prop = prop, alpha = alpha, ...))
  layer3 <- layer(data = data, mapping = mapping, stat =  StatConline, geom = GeomPath, 
                  position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
                  params = list(lineend = lineend, linejoin = linejoin, 
                                linemitre = linemitre, arrow = arrow, alpha=0.3, na.rm = na.rm, ...))
  if (!missing(nudge_x) || !missing(nudge_y)) {
    if (!missing(position)) {
      stop("Specify either `position` or `nudge_x`/`nudge_y`", 
           call. = FALSE)
    }
    position <- position_nudge(nudge_x, nudge_y)
  }
  layer4 <- layer(data = data, mapping = mapping, stat = StatLabel, geom = GeomLabel, 
                  position = position_jitter(width = 0.015, height = 0.015, seed = 21), show.legend = FALSE, inherit.aes = inherit.aes, 
                  params = list(parse = parse, label.padding = label.padding, 
                                label.r = label.r, label.size = label.size, na.rm = na.rm, ...))
  return(list(layer1,layer2,layer3,layer4)[c(show.point, show.ellipse, show.line, show.label)])
}

## PCA #############

PCA <- function(cov_matrix, ind_label, pop_label, x_axis, y_axis, show.point=T, show.label=T, show.ellipse=T, show.line=T, alpha=0, index_exclude=vector())
{
  ## This function takes a covariance matrix and performs PCA. 
  # cov_matrix: a square covariance matrix generated by most pca softwares
  # ind_label: a vector in the same order and length as cov_matrix; it contains the individual labels of the individuals represented in the covariance matrix
  # pop_label: a vector in the same order and length as cov_matrix; it contains the population labels of the individuals represented in the covariance matrix
  # x_axis: an integer that determines which principal component to plot on the x axis
  # y_axis: an integer that determines which principal component to plot on the y axis
  # show.point: whether to show individual points
  # show.label: whether to show population labels
  # show.ellipse: whether to show population-specific ellipses
  # show.line: whether to show lines connecting population means with each individual point
  # alpha: the transparency of ellipses
  # index_exclude: the indices of individuals to exclude from the analysis
  index_include <- setdiff(seq_along(ind_label), index_exclude)
  m <- as.matrix(cov_matrix)
  m[is.na(m)]<- median(m, na.rm = T)
  m<-m[index_include, index_include] ## Remove 4SJH, four 3Ps individuals, and contaminated ones
  e <- eigen(m)
  e_value<-e$values
  x_variance<-e_value[x_axis]/sum(e_value)*100
  y_variance<-e_value[y_axis]/sum(e_value)*100
  e <- as.data.frame(e$vectors)
  e <- cbind(ind_label[index_include], pop_label[index_include], e) ## with the above individuals removed
  #colnames(e)[3:331]<-paste0("PC",1:329)
  colnames(e)[3:(dim(e)[1])]<-paste0("PC",1:(dim(e)[1]-2)) ## with the above individuals removed
  colnames(e)[1:2]<-c("individual", "population")
  assign("pca_table", e, .GlobalEnv)
  
  PCA_plot<-ggplot(data=e[,],aes(x=e[,x_axis+2], y=e[,y_axis+2], color=population,label=population, shape=population)) + 
    geom_enterotype(alpha=alpha, show.point=show.point, show.label=show.label, show.ellipse=show.ellipse, show.line=show.line) +
    scale_shape_manual(values = c(rep(c(15,16,17,18),7), 15, 16)) +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank()
    ) +
    xlab(paste0("PC", x_axis, "(",round(x_variance,2),"%)")) +
    ylab(paste0("PC", y_axis ,"(",round(y_variance,2),"%)")) 
  print(PCA_plot)
}

#PCA but set colors manually
PCA_manual_colors <- function(cov_matrix, ind_label, pop_label, x_axis, y_axis, manual_colors, manual_shapes , show.point=T, show.label=T, show.ellipse=T, show.line=T, alpha=0,index_exclude=vector(), legend.pos = "right", legend.name = NULL)
{
  ## This function takes a covariance matrix and performs PCA. 
  # cov_matrix: a square covariance matrix generated by most pca softwares
  # ind_label: a vector in the same order and length as cov_matrix; it contains the individual labels of the individuals represented in the covariance matrix
  # pop_label: a vector in the same order and length as cov_matrix; it contains the population labels of the individuals represented in the covariance matrix
  # x_axis: an integer that determines which principal component to plot on the x axis
  # y_axis: an integer that determines which principal component to plot on the y axis
  # show.point: whether to show individual points
  # show.label: whether to show population labels
  # show.ellipse: whether to show population-specific ellipses
  # show.line: whether to show lines connecting population means with each individual point
  # alpha: the transparency of ellipses
  # index_exclude: the indices of individuals to exclude from the analysis
  index_include <- setdiff(seq_along(ind_label), index_exclude)
  m <- as.matrix(cov_matrix)
  m[is.na(m)]<- median(m, na.rm = T)
  m<-m[index_include, index_include] ## Remove 4SJH, four 3Ps individuals, and contaminated ones
  e <- eigen(m)
  e_value<-e$values
  x_variance<-e_value[x_axis]/sum(e_value)*100
  y_variance<-e_value[y_axis]/sum(e_value)*100
  e <- as.data.frame(e$vectors)
  e <- cbind(ind_label[index_include], pop_label[index_include], e) ## with the above individuals removed
  #colnames(e)[3:331]<-paste0("PC",1:329)
  colnames(e)[3:(dim(e)[1])]<-paste0("PC",1:(dim(e)[1]-2)) ## with the above individuals removed
  colnames(e)[1:2]<-c("individual", "population")
  assign("pca_table", e, .GlobalEnv)
  
  PCA_plot<-ggplot(data=e[,],aes(x=e[,x_axis+2], y=e[,y_axis+2], color=population,label=population, shape=population)) + 
    geom_enterotype(alpha=alpha, show.point=show.point, show.label=show.label, show.ellipse=show.ellipse, show.line=show.line) +
    labs(col = legend.name, shape = legend.name) +
    scale_shape_manual(values = manual_shapes) +
    scale_color_manual(values = manual_colors) +
    theme_bw() +
    theme(legend.position = legend.pos,
          axis.text.x = element_blank(),
          axis.text.y = element_blank()) +
    xlab(paste0("PC", x_axis, "(",round(x_variance,2),"%)")) +
    ylab(paste0("PC", y_axis ,"(",round(y_variance,2),"%)")) 
  print(PCA_plot)
}

## DAPC ########

DAPC <- function (n=50, x_axis, y_axis, show.point=T, show.label=T, show.ellipse=T, show.line=T, alpha=0) {
  ## This function should follow immediately after a PCA function. It takes the PCA output and performs a linear discriminant analysis
  # n: the number of principal components to use for discriminant analysis
  # x_axis: an integer that determines which linear discriminant to plot on the x axis
  # y_axis: an integer that determines which linear discriminant to plot on the y axis
  # show.point: whether to show individual points
  # show.label: whether to show population labels
  # show.ellipse: whether to show population-specific ellipses
  # show.line: whether to show lines connecting population means with each individual point
  # alpha: the transparency of ellipses
  fit <- lda(population ~ ., data=pca_table[,2:(n+1)], na.action="na.omit", CV=F, output = "Scatterplot")
  plda <- predict(object = fit,
                  newdata = pca_table[,2:(n+1)])
  prop.lda = fit$svd^2/sum(fit$svd^2)
  dataset = data.frame(group = pca_table[,2], lda = plda$x)
  DAPC_plot<- ggplot(dataset, aes(dataset[,1+x_axis], dataset[,1+y_axis], color= group, label=group, shape=group)) + 
    geom_enterotype(show.point=show.point, show.label=show.label, show.ellipse=show.ellipse, show.line=show.line, alpha=alpha) +
    scale_shape_manual(values = c(rep(c(15,16,17,18),7), 15, 16)) +
    labs(x = paste0("LD", x_axis," (", memisc::percent(prop.lda[x_axis]), ")", sep=""),
         y = paste0("LD", y_axis, " (", memisc::percent(prop.lda[y_axis]), ")", sep=""))
  print(DAPC_plot)
}

## PCoA ########

PCoA <- function(dist_matrix, ind_label, pop_label, k, x_axis, y_axis, show.point=T, show.label=T, show.ellipse=T, show.line=T, alpha=0, index_exclude=vector())
{
  ## This function takes a pairwise distance matrix and performs PCoA
  # dist_matrix: a square distance matrix generated by most pca softwares
  # ind_label: a vector in the same order and length as dist_matrix; it contains the individual labels of the individuals represented in the covariance matrix
  # pop_label: a vector in the same order and length as dist_matrix; it contains the population labels of the individuals represented in the covariance matrix
  # x_axis: an integer that determines which principal component to plot on the x axis
  # y_axis: an integer that determines which principal component to plot on the y axis
  # show.point: whether to show individual points
  # show.label: whether to show population labels
  # show.ellipse: whether to show population-specific ellipses
  # show.line: whether to show lines connecting population means with each individual point
  # alpha: the transparency of ellipses
  # index_exclude: the indices of individuals to exclude from the analysis
  
  index_include <- setdiff(seq_along(ind_label), index_exclude)
  m <- as.matrix(dist_matrix)
  m[is.na(m)]<- median(m, na.rm = T)
  m <- m[index_include, index_include] ## Remove 4SJH, four 3Ps individuals, and contaminated ones
  mds <- cmdscale(as.dist(m), k=k)
  mds <- as.data.frame(mds)
  colnames(mds) <- paste0("dist_", 1:k)
  mds <- cbind(ind_label[index_include], pop_label[index_include], mds)
  colnames(mds)[1:2]<-c("individual", "population")
  eigen_value <- cmdscale(as.dist(m), k=k, eig = T)$eig
  var_explained <- round(eigen_value/sum(eigen_value)*100, 2)
  assign("pcoa_table", mds, .GlobalEnv)
  assign("var_explained", var_explained, .GlobalEnv)
  
  PCoA_plot<-ggplot(data=mds[,], aes(x=mds[,x_axis+2], y=mds[,y_axis+2], color=population,label=population, shape=population)) + 
    geom_enterotype(alpha=alpha, show.point=show.point, show.label=show.label, show.ellipse=show.ellipse, show.line=show.line) +
    scale_shape_manual(values = c(rep(c(15,16,17,18),7), 15, 16)) +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank()
    ) +
    xlab(paste0("PCo", x_axis, "(",var_explained[x_axis],"%)")) +
    ylab(paste0("PCo", y_axis ,"(",var_explained[y_axis],"%)")) 
  print(PCoA_plot)
}


#functions to help plotting with ggplot
#primarily contains functions that return objects for making manhattam plots

#############################################################################################################################################################################
#############################################################################################################################################################################
#####                                                                                                                                                                   #####
#####                                                     functions for efficient plotting                                                             #####
#####                                                                                                                                                                   #####
#############################################################################################################################################################################
#############################################################################################################################################################################

saveToPDF <- function(...) {
    d = dev.copy(pdf,...)
    dev.off(d)
}

saveToPNG <- function(...) {
    d = dev.copy(png,...)
    dev.off(d)
}

#############################################################################################################################################################################
#############################################################################################################################################################################
#####                                                                                                                                                                   #####
#####                                                     functions of creating objects for Manhattan lots                                                             #####
#####                                                                                                                                                                   #####
#############################################################################################################################################################################
#############################################################################################################################################################################



#collect chromosome information
#date frame requires following columns:
# Chromosome   as "CHR"
# SNP position as "POS"
chr_info <- function(SNP_df = NULL){chr_info <- SNP_df %>% group_by(CHR) %>% summarise(Length = max(POS)) %>% 
                                                arrange(factor(CHR, levels = unique(SNP_df$CHR)))         %>% 
                                                mutate(tot = cumsum(Length)-Length, 
                                                       LG  = c(1:length(unique(SNP_df$CHR))))
                                                return(chr_info)
}

#create data frame that contains information for where chromosome start and stop with in plot
#requires: 
# CHR
# LG
# BPcum
axisdf  <- function(SNP_df){ axisdf <- SNP_df %>% group_by(CHR) %>% summarize(LG     =  first(LG),
                                                                              center = (max(BPcum) + min(BPcum) ) / 2, 
                                                                              start  =  min(BPcum),
                                                                              stop   =  max(BPcum),
                                                                              ymin   =  -Inf,
                                                                              ymax   =  Inf) %>% arrange(LG)
                                                                              return(axisdf)
}
#create x-scale
x_scale <- function(axisdf){ x_scale <- scale_x_continuous(label=axisdf$CHR,breaks = axisdf$center,expand=c(0,0))
                            return(x_scale)
}


#create rectangles for background in ggplot
bg_rect <- function(axisdf){ bg_rect <- annotate("rect", xmin = axisdf$start,xmax = axisdf$stop, 
                                                         ymin = axisdf$ymin, ymax = axisdf$ymax, 
                                                         alpha = .5, 
                                                         fill = rep(c("darkgrey","white"),(nrow(axisdf)+1)/2)[1:nrow(axisdf)])
                                                         return(bg_rect)
}

#make manhattan plot
manhattan_outliers <- function(snp_df=NULL,outlier_df=NULL,out="manhattan_plot.png"){
  #internal functions:
  chr_info <- function(SNP_df = NULL){chr_info <- SNP_df %>% group_by(CHR) %>% summarise(Length = max(POS)) %>% 
    arrange(factor(CHR, levels = unique(SNP_df$CHR)))         %>% 
    mutate(tot = cumsum(Length)-Length, 
           LG  = c(1:length(unique(SNP_df$CHR))))
  return(chr_info)}
  axisdf  <- function(SNP_df){ axisdf <- SNP_df %>% group_by(CHR) %>% summarize(LG     =  first(LG),
                                                                                center = (max(BPcum) + min(BPcum) ) / 2, 
                                                                                start  =  min(BPcum),
                                                                                stop   =  max(BPcum),
                                                                                ymin   =  -Inf,
                                                                                ymax   =  Inf) %>% arrange(LG)
  return(axisdf)}
  x_scale <- function(axisdf){ x_scale <- scale_x_continuous(label=axisdf$CHR,breaks = axisdf$center,expand=c(0,0))
  return(x_scale)}
  bg_rect <- function(axisdf){ bg_rect <- annotate("rect", xmin = axisdf$start,xmax = axisdf$stop, 
                                                   ymin = axisdf$ymin, ymax = axisdf$ymax, 
                                                   alpha = .5, 
                                                   fill = rep(c("darkgrey","white"),(nrow(axisdf)+1)/2)[1:nrow(axisdf)])
  return(bg_rect)}
  #get chromosome information 
  chr_info <- chr_info(snp_df)    
  #add chr info to SNPinfo
  snp_df <- left_join(snp_df,chr_info[,c("CHR","LG","tot")]) %>% 
    arrange(LG,POS)                               %>% 
    mutate(BPcum = tot + POS)
  #ggplot objects
  axisdf <- axisdf(snp_df)
  bg_rect <- bg_rect(axisdf)
  x_scale <- x_scale(axisdf)
  
  #merge dataframes
  outlier_df$snp <- "nonsignf"
  outlier_df$snp[which(outlier_df$OutlierFlag)] <- "signf"
  outlier_df$snp[which(outlier_df$selected)] <- "select"
  outlier_df$snp <- factor(outlier_df$snp, levels = c("nonsignf","signf","select"))
  snp_df      <- left_join(snp_df,outlier_df) %>% filter(!is.na(FST))
  N_nignf     <- length(which(snp_df$OutlierFlag))
  N_selected  <- length(which(snp_df$selected))
  #make plot
  plot <- ggplot(snp_df,aes(x=BPcum,y=FST,color=snp)) + 
    bg_rect +
    geom_point() + 
    labs(title = "Manhattan plot",
         subtitle = paste(N_nignf,"significant SNPs,",N_selected,"selected"),
         x = "Chromosome",
         y = expression(italic(F)[ST])) +
    xlab("Chromosome")+
    coord_cartesian(ylim = c(min(snp_df$FST,na.rm =T)
                             ,max(snp_df$FST,na.rm =T)))+
    scale_color_manual(values = c("darkgrey","yellow","red")) +
    theme_bw() + 
    x_scale + 
    theme(panel.grid  = element_blank(),
          plot.margin = unit(c(0,0,0.0,0), units = "in"),
          axis.text.x = element_text(angle = 0, size = 12))
  ggsave(plot = plot, filename = out, width = 16, height = 4, dpi = 300, units = "in")
  #return(plot)
}

### distance
#function that estimates the ditance betrween locations taking into account boundaries between terrestrial and marine habitat
#required are:
# a shapefile of the terrestrial habitat of the area of interest
# extent of the area of interest in a vector as: c(min_long, max_long, min_lath, max_lath)
# GPS coordinates of samples or sample locations
# if interested in marine habitats set marine = TRUE! this flips the matrix
# if plots = TRUE, plot showing extent of area is shown
distance <- function(name       = NULL,
                     long       = NULL,
                     lath       = NULL,
                     shapefile  = NULL,
                     extent     = NULL,
                     gridsize   = 0.1,
                     directions = 8,
                     marine     = FALSE,
                     plot       = FALSE){#collect GPS information samples
                                         GPS_info <- data.frame(NAME = name,
                                                                LONG = long,
                                                                LATH = lath)  
                                          
                                         GPS_info <- GPS_info %>% group_by(NAME) %>% summarise(LONG = mean(LONG), LATH = mean(LATH))
                                          
                                         #create raster
                                         r <- raster(extent(extent), res=gridsize)
                                         rr <- rasterize(shapefile, r)
                                         
                                         #if marine is TRUE flip raster  
                                         if(marine){
                                            roce <- rr
                                            roce[is.na(roce[])] <- -99
                                            roce[roce[]>=1] <- NA
                                            roce[roce[]==-99] <- 1
                                          } else { roce <- rr }
                                         
                                         if(plot){
                                            plot(roce,ext=c(extent))
                                         }
                                         
                                         #create transition object
                                         troce <- transition(roce, mean, directions = directions) #less than 1 min
                                         #corrext for local distances
                                         troceC <- geoCorrection(troce, "c", scl = F)
                                         #put GPS coordinates in matrix
                                         pC_region <- as.matrix(GPS_info[c("LONG","LATH")])
                                         
                                         pC_region <- as.matrix(GPS_info[,-1])
                                         #calculate least-cost distance
                                         cosDist_reg <- costDistance(troceC, pC_region)
                                         #convert to Km and round to the nearest Km
                                         cosDist_reg <- cosDist_reg/1000
                                         cosDist_reg <- round(cosDist_reg,0)
                                         #put distances in square-matrix
                                         cosDist_reg_mat <- as.data.frame(as.matrix(cosDist_reg))
                                         colnames(cosDist_reg_mat) <- GPS_info$NAME
                                         rownames(cosDist_reg_mat) <- GPS_info$NAME
                                         cosDist_reg_mat <- as.matrix(cosDist_reg_mat)
                                         
                                         N_inf <-  length(which(is.infinite(cosDist_reg_mat[1,])))
                                         if(N_inf != 0){
                                           print(paste("warning,",N_inf,"individual(s) have GPS coordinates located on excluded grids"))
                                           print(paste("please look at individual(s):", row.names(cosDist_reg_mat)[  which(is.infinite(cosDist_reg_mat[1,]))]))
                                           print("this can be caused by wrong GPS coordinates, or coordinates too close to borders of of the grid edges")
                                           print("for now, slightly move GPS coorddinates so individuals are analysed correctly")
                                         }
                                         
                                         
                                         cosDist_reg_df <- reshape2::melt(as.matrix(cosDist_reg_mat),varnames = c("pop1","pop2"))
                                         colnames(cosDist_reg_df) <- c("pop1","pop2","distance")
                                         
                                         
                                         
                                         distance_list <- list()
                                         distance_list$matrix <- cosDist_reg_mat
                                         distance_list$dataframe <- cosDist_reg_df
                                         distance_list$raster <- roce
                                         return(distance_list)
}

### convert data frame to square matrix via distance matrix
df2squarematrix <- function(pop1 = NULL, pop2 = NULL, value = NULL){
  df <- data.frame(pop1 = pop1, pop2 = pop2, value = value)
  ## set up storage matrix
  # get names for row and columns
  nameVals <- sort(unique(unlist(df[1:2])))
  # construct 0 matrix of correct dimensions with row and column names
  myMat <- matrix(0, length(nameVals), length(nameVals), dimnames = list(nameVals, nameVals))
  # fill in the matrix with matrix indexing on row and column names
  myMat[as.matrix(df[c("pop1", "pop2")])] <- df[["value"]]
  #myMat <- myMat + t(myMat) #doubled all my values should remove after double checking
  return(myMat)
}

###
plotCorRes <- function(cor_mat, pop=NULL, ord=NULL, superpop=NULL,
                       title="Correlation of residuals", min_z=NA,max_z=NA, 
                       cex.main=1.5, cex.lab=1.5, cex.legend=1.5, color_palette=c("#001260", "#EAEDE9", "#601200"),
                       pop_labels = c(T,T), plot_legend = T, adjlab = 0.1, rotatelabpop=0, rotatelabsuperpop=0,lineswidth=1, lineswidthsuperpop=2,
                       adjlabsuperpop=0.16,cex.lab.2 = 1.5){
  
  op <- par(mfrow=c(1,1) ,mar=c(5,4,4,2) +0.1,xpd=F, oma=c(0,0,0,0))
  on.exit(par(op))
  
  N <- dim(cor_mat)[1]
  

  if(is.null(ord)&!is.null(pop)) ord <- order(pop)
  if(is.null(ord)&is.null(pop)) ord <- 1:nrow(cor_mat)

  if(is.null(pop)){
      pop <- rep(" ", nrow(cor_mat))
      lineswidth <- 0
  }
    
  pop<-pop[ord]
  
  N_pop <- vapply(unique(pop[ord]), function(x) sum(pop==x),1)
  
  cor_mat <- cor_mat[ord,ord]
  
  ## Set lower part of matrix as population mean correlation
  mean_cors <- matrix(ncol=length(unique(pop)), nrow=length(unique(pop)))
  colnames(mean_cors) <- unique(pop)
  rownames(mean_cors) <- unique(pop)
  
  for(i1 in 1:(length(unique(pop)))){
    for(i2 in 1:(length(unique(pop)))){
      p1 <- unique(pop)[i1]
      p2 <- unique(pop)[i2]
      mean_cors[i1,i2]<- mean(cor_mat[which(pop==p1),
                                      which(pop==p2)][!is.na(cor_mat[which(pop==p1),
                                                                     which(pop==p2)])])
      
    }
  }
  
  for(i1 in 1:(N-1)){
    for(i2 in (i1+1):N){
      cor_mat[i1, i2] <- mean_cors[pop[i2], pop[i1]]
      
    }
  }
  
  z_lims <- c(min_z, max_z)
  
    if(all(is.na(z_lims))) z_lims <- c(-max(abs(cor_mat[!is.na(cor_mat)])),
                                         max(abs(cor_mat[!is.na(cor_mat)])))
  #if(all(is.null(z_lims))) max_z <- max(abs(cor_mat[!is.na(cor_mat)]))

    
  if(any(is.na(z_lims))) z_lims <- c(-z_lims[!is.na(z_lims)], z_lims[!is.na(z_lims)])
  #if(any(is.null(z_lims))) max_z <- z_lims[!is.null(z_lims)]

    min_z <- z_lims[1]
    max_z <- z_lims[2]
    
  diag(cor_mat) <- 10
  nHalf <- 10
  
  # make sure col palette is centered on 0
  Min <- min_z
  Max <- max_z
  Thresh <- 0
  
  ## Make vector of colors for values below threshold
  rc1 <- colorRampPalette(colors = color_palette[1:2], space="Lab")(nHalf)    
  ## Make vector of colors for values above threshold
  rc2 <- colorRampPalette(colors = color_palette[2:3], space="Lab")(nHalf)
  rampcols <- c(rc1, rc2)
  
  rampcols[c(nHalf, nHalf+1)] <- rgb(t(col2rgb(color_palette[2])), maxColorValue=256) 
  
  rb1 <- seq(Min, Thresh, length.out=nHalf+1)
  rb2 <- seq(Thresh, Max, length.out=nHalf+1)[-1]
  rampbreaks <- c(rb1, rb2)
  
  rlegend <- as.raster(matrix(rampcols, ncol=1)[length(rampcols):1,])
  if(plot_legend){
    layout(matrix(1:2,ncol=2), width = c(4,1),height = c(1,1))
    par(mar=c(5,4,4,0),oma=c(1,4.5,2,0))
  }else
    par(mar=c(5,4,4,5),oma=c(1,4.5,2,0))
  image(t(cor_mat), col=rampcols, breaks=rampbreaks,
        yaxt="n",xaxt="n", zlim=c(min_z,max_z),useRaster=T,
        main=title, 
        oldstyle=T,cex.main=cex.main,xpd=NA)
  image(ifelse(t(cor_mat>max_z),1,NA),col="darkred",add=T)
  if(min(cor_mat)<min_z) image(ifelse(t(cor_mat<min_z),1,NA),col="darkslateblue",add=T)
  image(ifelse(t(cor_mat==10),1,NA),col="black",add=T)
  
  # put pop info
  if(pop_labels[2])
    text(sort(tapply(1:length(pop),pop,mean)/length(pop)),-adjlab,unique(pop),xpd=NA,cex=cex.lab, srt=rotatelabpop)
  if(pop_labels[1])
    text(-adjlab,sort(tapply(1:length(pop),pop,mean)/length(pop)),unique(pop),xpd=NA, cex=cex.lab,srt=90-rotatelabpop)
  abline(v=grconvertX(cumsum(sapply(unique(pop),function(x){sum(pop==x)}))/N,"npc","user"),
         col=1,lwd=lineswidth,xpd=F)
  abline(h=grconvertY(cumsum(sapply(unique(pop),function(x){sum(pop==x)}))/N, "npc", "user"),
         col=1,lwd=lineswidth,xpd=F)
  
  # put superpop if not null
    if(!is.null(superpop)){
        superpop <- superpop[ord]
    if(pop_labels[2])
      text(sort(tapply(1:length(superpop),superpop,mean)/length(superpop)),-adjlabsuperpop,unique(superpop),xpd=NA,cex=cex.lab.2, srt=rotatelabsuperpop, font=2)
    if(pop_labels[1])
      text(-adjlabsuperpop,sort(tapply(1:length(superpop),superpop,mean)/length(superpop)),unique(superpop),xpd=NA, cex=cex.lab.2,srt=90-rotatelabsuperpop,font=2)
    abline(v=grconvertX(cumsum(sapply(unique(superpop),function(x){sum(superpop==x)}))/N,"npc","user"),
           col=1,lwd=lineswidthsuperpop,xpd=F)
    abline(h=grconvertY(cumsum(sapply(unique(superpop),function(x){sum(superpop==x)}))/N, "npc", "user"),
           col=1,lwd=lineswidthsuperpop,xpd=F)
  }
  
  if(plot_legend){
    par(mar=c(5,0.5,4,2))
    plot(c(0,1),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')    
    
    rasterImage(rlegend, 0, 0.25, 0.4,0.75)
    text(x=0.8, y = c(0.25,0.5, 0.75),
         labels = c(-max(abs(min_z),abs(max_z)), 0, max(abs(min_z),abs(max_z))),
         cex=cex.legend,xpd=NA)
  }
}


orderInds <- function(q=NULL, pop=NULL, popord=NULL){
  # Function to order individuals for admixture and evalAdmix plots. 
  # recommended is to use pop, then if q is given it will order within pop by admixture proporiton. poporder allows to pre-specify order of populations
  # if only q is given will group individuals by main cluster they are assigned
  
  ordpop <- function(x, pop, q){
    idx <- which(pop==x)
    main_k <- which.max(apply(as.matrix(q[idx,]),2,mean))
    ord <- order(q[idx,main_k])
    idx[ord]
  } 
  
  if(!is.null(pop)){
    
    if(is.null(popord)) popord <- unique(pop)
    
    if(!is.null(q)){ 
      
      ord <- unlist(sapply(popord, ordpop, pop=pop, q=q))
      
    } else if (is.null(q)) {
      
      ord <- unlist(sapply(popord, function(x) which(pop==x)))
      
    }
  } else if (is.null(pop)&!is.null(q)) {
    
    # get index of k with max value per individual
    main_k <- apply(q,1, which.max)
    
    # get max q per indivdiual
    main_q <- q[cbind(1:nrow(q),main_k)]
    
    ord <- order(main_k, main_q)
    
  } else {stop("Need at least an argument to order.")}

  return(ord)
  
}


orderK <- function(q, refinds= NULL,refpops = NULL, pop=NULL){
  # Function to order ancestral populations, useful to keep cluster colors in admix plot the same when comparing results across different k values
  # if you give refinds will use maximum Q value of each individual to define clusters
  # if you give refpops (must also give pops) will use maximum mean admixture proportions within inds from pop to define clusters
  # if any refpops or refinds have same cluster as maximum, the admixture plot will look really bad (you will lose a cluster and another will be twice)
  
  k <- ncol(q)
  kord <- integer(0)
  
  if(is.null(refinds)){
  refpops <- refpops[1:k]
  
  for(p in refpops){
    
    kord <- c(kord, which.max(apply(q[pop==p,],2,mean)))
    
  }
  } else {
    
    refinds <- refinds[1:k]
    
    for(i in refinds){
      
      kord <- c(kord, which.max(q[i,]))
    }
  }
  
    # if(any(rowSums(q[,kord]!=1))) warning("reordered admixture proportions don't sum to 1, make sure every refind or refpop defines a unique cluster.")

    return(kord)
}


plotAdmix <- function(q, pop=NULL, ord=NULL, inds=NULL,
                      colorpal= c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999"),
                      main=paste("Admixture proportions assuming K =",k),
                      cex.main=1.5, cex.lab=1, rotatelab=0,padj=0, cex.inds=1){
  # simple function to plot admixture proprotions, just to make sure the ordering of individuals is handled as in plotCorRes.
  
  k <- ncol(q)
  
  if(k>length(colorpal))
    warning("not enought colors for all Ks in palette.")
  
  # if(!is.null(ord)) if(!ord) ord <- 1:nrow(q)
  
  if(is.null(ord)&!is.null(pop)) ord <- order(pop)
  if(is.null(ord)&is.null(pop)) ord <- 1:nrow(q)
  
  barplot(t(q)[,ord], col=colorpal, space=0, border=NA, cex.axis=1.2,cex.lab=1.8,
          ylab="Admixture proportions", xlab="", main=main, cex.main=cex.main,xpd=NA)
  
  if(!is.null(inds)){
    text(x = 1:nrow(q) - 0.5,-0.1, inds[ord],xpd=NA,srt=90, cex=cex.inds)
  }
  
  if(!is.null(pop)){
    
    text(sort(tapply(1:length(pop),pop[ord],mean)),-0.05-padj,unique(pop[ord]),xpd=NA, srt=rotatelab, cex=cex.lab)
    abline(v=1:nrow(q), col="white", lwd=0.2)
    abline(v=cumsum(sapply(unique(pop[ord]),function(x){sum(pop[ord]==x)})),col=1,lwd=1.2)
    
  }
  
}


