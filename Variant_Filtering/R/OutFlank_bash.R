library(tidyverse, quietly = TRUE)
library(vcfR, quietly = TRUE)
library(qvalue, quietly = TRUE)
library(OutFLANK, quietly = TRUE)
library(glue, quietly = TRUE)
library(runner, quietly = TRUE)
library(SNPRelate, quietly = TRUE)
library(data.table, quietly = TRUE)
library(optparse, quietly = TRUE)
library(R.utils, quietly = TRUE)
library(scales, quietly = TRUE)
library(ggsci, quietly = TRUE)

#####################################################################################################################################
# import parameters
option_list <- list(
  #input
  make_option(c('--vcf_file')      ,action='store',type='character',default=NULL   ,help='gzipped VCF file '),
  make_option(c('--gds_file')      ,action='store',type='character',default=NULL   ,help='gds file, does not need to be supplied when vcf is given as input'),
  make_option(c('--pop_file')      ,action='store',type='character',default=NULL   ,help='population file with 2 nessesary columns, IND and POP'),
  #outout  
  make_option(c('--fst_file')      ,action='store',type='character',default=NULL   ,help='path to fst estimates'),
  make_option(c('--compress')      ,action='store',type='logical'  ,default=FALSE  ,help='if true, gzip fst estimates'),
  make_option(c('--out')           ,action='store',type='character',default=NULL   ,help='path to output directiry'),
  #plots
  make_option(c('--plots_out')     ,action='store',type='character',default=NULL   ,help='path to output directiry'),
  make_option(c('--fst_plot')      ,action='store',type='logical'  ,default=TRUE   ,help='if TRUE plot summary information on FST estimates'),
  make_option(c('--qval_plot')     ,action='store',type='logical'  ,default=TRUE   ,help='if TRUE plot summary information on qvalues'),
  make_option(c('--PCA_plot')      ,action='store',type='logical'  ,default=TRUE   ,help='if TRUE plot PCA'),
  make_option(c('--Manhattan_plot'),action='store',type='logical'  ,default=TRUE   ,help='if TRUE plot Manhattan plot'),
  #parameters
  make_option(c('--by_chr')        ,action='store',type='logical'  ,default=FALSE   ,help='if TRUE estimate fst estimates per chromosome, useful for large datasets'),
  make_option(c('--minNPOP')       ,action='store',type='numeric'  ,default=NULL   ,help='minimum number of populations that has are polymorphic for a given snp to be included in the analyses'),
  make_option(c('--qval')          ,action='store',type='numeric'  ,default=0.05   ,help='qval; see OutFlank manual for details'),
  make_option(c('--LTrim')         ,action='store',type='numeric'  ,default=0.05   ,help='left trim; ; see OutFlank manual for details'),
  make_option(c('--TRim')          ,action='store',type='numeric'  ,default=0.05   ,help='right triml ; see OutFlank manual for details'),
  make_option(c('--Hmin')          ,action='store',type='numeric'  ,default=0.10   ,help='minimal heterozygosity; ; see OutFlank manual for details'),
  make_option(c('--slw')           ,action='store',type='numeric'  ,default=100000 ,help='size non-overlapping sliding window to select independently segregating sites'),
  make_option(c('--filter_method') ,action='store',type='character',default="FST"  ,help='either FST or LD; use FST to select SNP with higest FST value per sliding window, or LD to use linkage disequilibrium filtering'),
  make_option(c('--min_snps')      ,action='store',type='numeric'  ,default=1      ,help='when filter_method is FST, set minimum number of outliers required before outliers are excepted'),
  make_option(c('--ld.threshold')  ,action='store',type='numeric'  ,default=0.2    ,help='when filter_method is LD , set rsquared threshold for LD based filtering of outliers')
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

#input
vcf_file       = c(opt$vcf_file)
gds_file       = c(opt$gds_file)
pop_file       = opt$pop_file
#output
fst_file	   = opt$fst_file
out            = opt$out
compress	   = opt$compress
#plots
plots_out      = opt$plots_out
fst_plot	   = opt$fst_plot
qval_plot      = opt$qval_plot
PCA_plot       = opt$PCA_plot
Manhattan_plot = opt$Manhattan_plot
#parameters
by_chr		   = opt$by_chr
minNPOP		   = opt$minNPOP  
qval           = opt$qval
Hmin           = opt$Hmin    
LTrim          = opt$LTrim   
RTrim          = opt$TRim    
sliding_window = opt$slw
filter_method  = opt$filter_method      
min_snps       = opt$min_snps
ld.threshold   = opt$ld.threshold       

source("./R/OutFlank_Rfunctions.R") #needs to be in same folder as this script
source("./R/general_Rfunctions.R")  #needs to be in same folder as this script

#####################################################################################################################################
#extension to genotype files
#out_dir <- dirname(out)
#out_ext <- basename(out)
#dir.create(glue("{out_dir}/fst"), recursive = TRUE)
#dir.create(glue("{out_dir}/plots"), recursive = TRUE)
#dir.create(glue("{out_dir}/outliers"), recursive = TRUE)
#fst_ext       <- glue::glue("{out_dir}/fst/{out_ext}")
#plots_ext     <- glue::glue("{out_dir}/plots/{out_ext}")
#outliers_ext  <- glue::glue("{out_dir}/outliers/{out_ext}")

#load pop file
sample_info <- read_tsv(glue::glue(pop_file))
#make set color df
Npops      <- length(unique(sample_info$POP))
hex_codes2 <- hue_pal()(Npops)
shapes     <- c(rep(c(15,16,17,18),10))[1:Npops]

#set colors and node shpes for pops
pop_cols <- sample_info   %>% 
            group_by(POP) %>% 
            summarise(comparison = first(POP)) %>% #had GENETIC_CLUSTER here but needs to work without that variable 
            mutate(color = hex_codes2,
                   shape = shapes)
				   
#####################################################################################################################################
#load data
print("loading genotype data")
gds <- snpgdsReadGDS(vcf_file = vcf_file, gds_file = gds_file)

#get SNP info
print("get SNP info")
bim <- snpgdsSNPsum(gds = gds)

#get chromosome information 
#requires columns CHR and POS
chr_info <- chr_info(bim)    

#add chr info to SNPinfo
#cumilative base pair (BPcum) is added for manhattan plotting
bim <- left_join(bim,chr_info[,c("CHR","LG","tot")]) %>% 
       arrange(LG,POS)                               %>% 
       mutate(BPcum = tot + POS)

#select sample_info
sample_info <- sample_info %>% filter(IND %in% read.gdsn(index.gdsn(gds, "sample.id")))
pops   <- as.character(unique(sample_info$POP))
N_pops <- length(pops)

pop_cols_select <- pop_cols %>% filter(POP %in% sample_info$POP)
colors <- pop_cols_select$color
shapes <- pop_cols_select$shape

#####################################################################################################################################
#estimate fst values
#saves output to file and can reload files created
print("estimating FST values")

if(isTRUE(compress))
 {print("compress is TRUE, gzip fst file")
  fst_file <- paste0(fst_file,"_fst.tsv.gz")} else {paste0(fst_file,"_fst.tsv")}

fst <- OutFlankgetFST(gds       = gds,
                      pop.id    = sample_info$POP,
                      by_chr    = by_chr,
                      fst_file  = fst_file,
                      plots_out = plots_out,
					  fst_plot  = fst_plot)

#####################################################################################################################################
#get maf per POP
#last column shows the number of populations contain the rare allele
#use this to filter out SNPs that are polymorphic in only 1 or 2 populations and might bias your analyses					  
if(!is.null(minNPOP)){
print("minNPOP was set, extracting ALT allele counts per POP and filetering sites")
mac_POP        <- snpgdsSNPmaf_byPOP(gds = gds, pop.id = sample_info$POP, mac = TRUE)
mac_POP$LOC    <- bim$LOC #waarom is dit ineens nodig??
#get SNP.ids from sites that are found in 3 or more sample locations
mac_POP_select <- mac_POP$LOC[which(mac_POP$N >= minNPOP)]

#add line that filters out unwanted sites based on filtering criteria above 
fst_selected <- fst %>% filter(LOC %in% mac_POP_select)
} else { fst_selected <- fst  } 


#####################################################################################################################################
#run outflank function
print("running OutFlank")
OutFlank_list <- OutFLANK(FstDataFrame      = fst_selected,
                          NumberOfSamples   = N_pops,
                          qthreshold        = qval,
                          Hmin              = Hmin,
                          LeftTrimFraction  = LTrim,
                          RightTrimFraction = RTrim)

print("selecting outliers")
outlier_plots <- if(!is.null(plots_out)){glue::glue("{plots_out}_qval{qval}_Hmin{Hmin}")} else {NULL}
outliers <- OutFlankgetoutliers(gds            = gds,
                                OutFlank_list  = OutFlank_list,
								plots_out      = outlier_plots,
								qval_plot      = qval_plot,
								out            = glue::glue("{out}_qval{qval}_Hmin{Hmin}"),
                                sliding_window = sliding_window,
                                filter_method  = filter_method,
                                min_snps       = min_snps,     #when method is "FST"
                                ld.threshold   = ld.threshold) #when method is "LD"

if(!is.null(plots_out) & isTRUE(PCA_plot)){
print("estimate Genetic Relationship Matrix")
GRM <- snpgdsGRM(gdsobj        = gds,
                 sample.id     = sample_info$IND,
                 snp.id        = outliers$LOC,
                 method        = "Eigenstrat",
                 autosome.only = FALSE)

print("plotting PCA")
png(filename = glue::glue("{plots_out}_qval{qval}_Hmin{Hmin}_GRM.png"), width = 8, height = 4, units = "in", res = 300)
PCA_manual_colors(GRM$grm, sample_info$IND ,sample_info$POP,
            1, 2, show.point=T, show.label=F, manual_colors = colors, manual_shapes = shapes,
            show.ellipse=T, show.line=T, alpha=0)         
dev.off()
}
#####################################################################################################################################
#create manhattan plot 
if(!is.null(plots_out) & isTRUE(Manhattan_plot)){
print("creating Manhattan plot")
comp_cols <- scale_color_manual(values = ggsci::pal_npg("nrc")(5)[c(3)]) 

axisdf <- axisdf(bim)
bg_rect <- bg_rect(axisdf)
x_scale <- x_scale(axisdf)

theme_last <- theme(panel.grid  = element_blank(),
                    plot.margin = unit(c(0,0,0.0,0), units = "in"),
                    axis.text.x = element_text(angle = 0, size = 12))

test <- left_join(fst_selected,bim)

ggplot(test,aes(x=BPcum)) + bg_rect + 
  geom_point(aes(y=FST)) +
  geom_point(data = outliers, aes(x=BPcum,y=FST),color="green")+
  xlab("Chromosome")+
  coord_cartesian(ylim = c(min(test$FST,na.rm =T),max(test$FST,na.rm =T)))+
  ylab(expression(italic(F)[ST]))+ comp_cols + theme_bw() + x_scale + theme_last
ggsave(filename = glue::glue("{plots_out}_qval{qval}_Hmin{Hmin}_Manhattan.png"), width = 12, height = 6, units = "in", dpi = 300  )
}
#####################################################################################################################################
print("closing gds file")
SNPRelate::snpgdsClose(gds)
