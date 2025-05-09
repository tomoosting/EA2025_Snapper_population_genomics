# load packages
library(dplyr)
library(optparse)
#library(SeqArray)
library(SNPRelate)
library(stringr)
library(adegenet)
library(vcfR)

#arguments
option_list <- list(
  make_option(c('--gzvcf')         ,action='store',type='character',default=NULL,help='file containing all SNP from vcf file'),
  make_option(c('--snprelate_out') ,action='store',type='character',default=NULL,help='...'),
  make_option(c('--seqarray_out')  ,action='store',type='character',default=NULL,help='...'),
  make_option(c('--vcfR_out')      ,action='store',type='character',default=NULL,help='...'),
  make_option(c('--genlight_out')  ,action='store',type='character',default=NULL,help='...')
  )

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser) 

###gds for snprelate
if(!is.null(opt$snprelate_out)){
	SNPRelate::snpgdsVCF2GDS(c(opt$gzvcf), paste0(opt$snprelate_out,".gds") , method="biallelic.only")  
	gds_snprelate <- SNPRelate::snpgdsOpen(paste0(opt$snprelate_out,".gds"),readonly = FALSE)
	gdsfmt::add.gdsn(gds_snprelate,"snp.id", paste0(read.gdsn(index.gdsn(gds_snprelate, "snp.chromosome")),
													":",
                                                    read.gdsn(index.gdsn(gds_snprelate, "snp.position"))) ,
													replace = TRUE)
	gdsfmt::closefn.gds(gds_snprelate)
} 

###gds for SeqArray
if(!is.null(opt$seqarray_out)){
	SeqArray::seqVCF2GDS(opt$gzvcf, paste0(seqarray_out,".gds"))
	gds_seqarray <- SeqArray::seqOpen(paste0(seqarray_out,".gds"), readonly = FALSE)
	SeqArray::seqAddValue(gds_seqarray,
						varnm = "variant.id", 
						paste0(seqGetData(gds_seqarray, "chromosome"),":",seqGetData(gds_seqarray, "position")),
						replace = T)
	gdsfmt::closefn.gds(gds_seqarray)
}	

###vcfR or genlight object
if(!is.null(opt$vcfR_out) | !is.null(opt$genlight_out)){
	vcfR <- vcfR::read.vcfR(opt$gzvcf)
	
	if(!is.null(opt$vcfR_out)){ 
		save(vcfR,file = paste0(opt$vcfR_out,".vcfR.R"))
	}
	
	if(!is.null(opt$genlight_out)){
		genlight <- vcfR::vcfR2genlight(vcfR)
		save(genlight,file = paste0(opt$genlight_out,".genlight.R"))
	}
}

