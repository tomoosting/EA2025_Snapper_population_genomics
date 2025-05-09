#vcf2bim
library(SNPRelate)
library(glue)
library(tidyverse)
library(gdsfmt)
library(optparse)
library(genio)

#arguments
option_list <- list(
  make_option(c('--gds_file') ,action='store',type='character',default=NULL,help='gds file'),
  make_option(c('--pop_file') ,action='store',type='character',default=NULL,help='tsv poplation file with atleast IND, POP, and name data set culumns'),
  make_option(c('--out')      ,action='store',type='character',default=NULL,help='file name without .bed, .bim or .fam'),
  make_option(c('--exclude')  ,action='store',type='character',default=NULL,help='vector of sample names to exlude: c("IND1","IND2","IND#")')
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser) 

gds_file <- c(opt$gds_file)
pop_file <- opt$pop_file
out      <- opt$out
exclude  <- opt$exclude

###############################################################################################
#functions
#open gds file
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


#get bim file from gds
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

###############################################################################################

#open gds file
gds <- snpgdsReadGDS(gds_file = gds_file)

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
if(!is.null(pop_file)){
  print("pop_file provided, generating fam file")
  print("looking for IND and POP fields")
  
  #read pop_file, needs IND and POP!
  pop_df <- read_tsv(pop_file)
  
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

#close gds file
gdsfmt::closefn.gds(gds)