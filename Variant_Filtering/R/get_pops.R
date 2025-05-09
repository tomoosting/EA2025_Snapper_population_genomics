# load packages
library(dplyr)
library(optparse)
library(stringr)
library(glue)
library(data.table)
library(readr)

#arguments
option_list <- list(
  make_option(c('--pop_file')     ,action='store',type='character',default=NULL ,help='poplation file with atleast IND, POP, and name data set culumns'),
  make_option(c('--set_name')     ,action='store',type='character',default=NULL ,help='column name this contains TRUE and FALSE vales for which individuals to include'),
  make_option(c('--out')          ,action='store',type='character',default=NULL ,help='output path'),
  make_option(c('--split_by')     ,action='store',type='character',default=NULL ,help='seperate ouputfiles per population, indicate column what group names'),
  make_option(c('--colnames')     ,action='store',type='logical'  ,default=TRUE ,help='include colnames'),
  make_option(c('--samples_only') ,action='store',type='logical'  ,default=FALSE,help='only output sample names')
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser) 

infile <- fread(opt$pop_file)
name  <- as.character(opt$set_name)

if(is.null(opt$split_by)){
  
  print("writing single output file")
  
  if(opt$samples_only){
	  subset <- infile[,c("IND")][which(infile[[name]])]
  } else {
	  subset <- infile[,c("IND","POP")][which(infile[[name]])]
  }
  write_tsv(subset, file = glue("{opt$out}.tsv"), col_names = opt$colnames)
  
} else {
  
  print("writing one output file per population")
  group <- as.character(opt$split_by)
  subset <- infile[which(infile[[name]]),]
  subset[[group]] <- as.factor(subset[[group]])
  
  if(opt$samples_only){
	  groups_list <- split(subset[,c("IND")],subset[[group]])
  } else {
	  groups_list <- split(subset[,c("IND","POP")],subset[[group]])
  }
  groups <- as.character(unique(subset[[group]]))
  for(i in groups){write_tsv(groups_list[[i]], file = glue("{opt$out}_{i}.tsv"),col_names = opt$colnames)}
}
