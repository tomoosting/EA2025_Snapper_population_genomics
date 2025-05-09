#custom fucntion for running OutFlank
#required packages
library(OutFLANK, quietly = TRUE)
library(rlist, quietly = TRUE)
library(readr, quietly = TRUE)
library(glue, quietly = TRUE)
library(runner, quietly = TRUE)

#############################################################################################################################################################################
#############################################################################################################################################################################
#####                                                                                                                                                                   #####
#####                                                               functions for running OutFlank                                                                      #####
#####                                                                                                                                                                   #####
#############################################################################################################################################################################
#############################################################################################################################################################################

#estimate FST values per CHR using gds file
OutFlankgetFST <- function(gds       = NULL,
                           pop.id    = NULL,
                           fst_file  = NULL,
                           by_chr    = FALSE,
                           plots_out = NULL,
						   fst_plot  = TRUE){if(file.exists(fst_file)){
                                                  print("FST file present - loading file and FST estimates will not be recalculated")
                                                  fst_estimates <- read_tsv(fst_file, progress = TRUE)
                                                  return(fst_estimates)
                                            } else {
                                                  #obtain SNP info
                                                  CHR <- unique(read.gdsn(index.gdsn(gds, "snp.chromosome")))
                                                  N_CHR <- length(CHR)
                                                  SNP_info <- data.frame(CHR = read.gdsn(index.gdsn(gds, "snp.chromosome")),
                                                                         LG  = as.numeric( str_remove_all(read.gdsn(index.gdsn(gds, "snp.chromosome")), "[:alpha:]|[:punct:]" )),
                                                                         LOC = paste0(read.gdsn(index.gdsn(gds, "snp.chromosome")),
                                                                                      ":",
                                                                                      read.gdsn(index.gdsn(gds, "snp.position"))),
                                                                         POS = read.gdsn(index.gdsn(gds, "snp.position")))
                                                  #extract genotype matrix
                                                  GT <- snpgdsGetGeno(gds)
                                                  if(!by_chr){
                                                  print("estimating fst values for total dataset, if this crashes try setting by by_chr to TRUE")
												  #estimate FST for entire data set  
                                                  G <- GT
                                                  G[is.na(G)] <-  9
                                                  G[G== 2]    <- -2
                                                  G[G== 0]    <-  2
                                                  G[G==-2]    <-  0
                                                  snp <- G
                                                  #estimate fst
                                                  fst <- MakeDiploidFSTMat(snp, locusNames = SNP_info$LOC, popNames = pop.id)
                                                  fst_estimates <- fst %>% mutate(POS       = SNP_info$POS
                                                                                 ,FST_ratio = FST/FSTNoCorr
                                                                                 ,LOC       = SNP_info$LOC
                                                                                 ,CHR       = SNP_info$CHR
                                                                                 ,LG        = SNP_info$LG)
                                                  } else {
                                                  print("performing fst estimates per linkage group")
                                                  #make empty list for FST values to be saved
                                                  fst_list <- list()
                                                  #loop over CHRs for estimating FST
                                                     for(i in CHR){
                                                         print(glue("calculating fst for {i}"))
                                                         #extract genotype info from GT matrix  
                                                         LG_rows <- which(SNP_info$CHR == i)
                                                         LG_GT  <- GT[,LG_rows]
                                                         LG_info <- SNP_info %>% filter(CHR == i)
                                                         #create matrix per CHR
                                                         G <- LG_GT
                                                         G[is.na(G)] <-  9
                                                         G[G== 2]    <- -2
                                                         G[G== 0]    <-  2
                                                         G[G==-2]    <-  0
                                                         snp <- G
                                                         #estimate fst
                                                         fst <- MakeDiploidFSTMat(snp, locusNames = LG_info$LOC, popNames = pop.id)
                                                         fst <- fst %>% mutate(POS       = LG_info$POS
                                                                              ,FST_ratio = FST/FSTNoCorr
                                                                              ,LOC       = LG_info$LOC
                                                                              ,CHR       = LG_info$CHR
                                                                              ,LG        = LG_info$LG)
                                                         fst_list[[paste0("fst_",i)]] <- fst
                                                         }     
                                                     #bind all fst estimates in one table
                                                     fst_estimates <- rlist::list.rbind(fst_list)
                                                  }
												  print("writing fst estimates to file")
                                                  write_tsv(fst_estimates,file = fst_file)
                                                  ### summary plot
                                                  if(!is.null(plots_out) & isTRUE(fst_plot)){
                                                       png(glue::glue("{plots_out}_FST_correction.png"),res = 300, units = "in", width = 12, height = 5)
                                                       par(mfrow=c(1,2))
                                                       plot(fst_estimates$FST, fst_estimates$FSTNoCorr, 
                                                            xlab="FST", ylab="uncorrected FST",
                                                            pch=20, col="blue") ; abline(0,1)
                                                       hist(fst_estimates$FSTNoCorr,xlim=c(0,round(max(fst_estimates$FSTNoCorr,na.rm = TRUE),3))  ,breaks=seq(0,round(max(fst_estimates$FSTNoCorr,na.rm = TRUE),3)+0.01, by=0.001))  
                                                       dev.off()
                                                  }
                                                  return(fst_estimates)
                                          }                                           
                             }

#select outliers using the output from the output obtained from the OutFlank script
OutFlankgetoutliers <- function(gds            = NULL,
                                OutFlank_list  = NULL,
                                out            = NULL,
                                sliding_window = 20000,
                                filter_method  = "FST",
                                min_snps       = 1,
                                ld.threshold   = 0.2,
                                plots_out      = NULL,
								qval_plot      = TRUE){bim      <- snpgdsSNPsum(gds = gds)                   # get SNP info from gds file 
                                                       # get CHR info
                                                       chr_info <- bim %>% group_by(CHR) %>% summarise(Length = max(POS)) %>% 
                                                                   arrange(factor(CHR, levels = unique(bim$CHR) ))        %>% 
                                                                   mutate(tot = cumsum(Length)-Length, 
                                                                          LG  = c(1:length(unique(bim$CHR))))
											   
                                                       bim      <- left_join(bim,chr_info[,c("CHR","LG","tot")]) %>% 
                                                                   arrange(LG,POS)                               %>% 
                                                                   mutate(BPcum = tot + POS)                 # get cumulative base pair count for Manhattan plots
                                                       OutFlank <- left_join(OutFlank_list$results,bim)      # join OutFlank results with SNP info
                                                       OutFlank <- dplyr::filter(OutFlank,!base::is.na(FST)) # remove SNP with NA as FST value
                                                       
                                                       outliers <- ungroup(OutFlank[which(OutFlank$OutlierFlag==TRUE),]) # get df with SNP identified as outliers
                                                       N_outliers <- nrow(outliers)                                      # get number of outliers
                                                       
                                                       if(N_outliers == 0){
                                                         print(glue("0 outliers were found"))
                                                       } else {
                                                         print(glue("{N_outliers} outlier(s) were found"))
                                                         
                                                         write_tsv(outliers, file = glue::glue("{out}_all_outliers_N{N_outliers}.tsv"))
                                                         # select outlier SNPs per sliding window
                                                         if(filter_method == "FST"){
                                                           
                                                           print("applying FST method for selecting outliers per sliding window")
                                                           #create list to save outliers per CHR to
                                                           runner_list <- list()
                                                           for(i in chr_info$CHR){
                                                             LG_outliers <- outliers %>% filter(CHR == i)  
                                                             #select highest FST outlier per sliding window
                                                             runner_LOC             <-   runner(x   = LG_outliers
                                                                                               ,idx = "BPcum"  
                                                                                               ,f   = function(x){x[["LOC"]][[which.max(x[["FST"]])]] }     
                                                                                               ,k   = sliding_window
                                                                                               ,at  = seq(from =1, to = max(outliers$BPcum,na.rm = TRUE)+sliding_window, by = sliding_window)
                                                                                               ,na_pad = FALSE)
                                                             #get number of outliers per sliding window
                                                             runner_N               <-   runner(x   = LG_outliers
                                                                                               ,idx = "BPcum"  
                                                                                               ,f   = function(x){nrow(x)}     
                                                                                               ,k   = sliding_window
                                                                                               ,at  = seq(from =1, to = max(outliers$BPcum,na.rm = TRUE)+sliding_window, by = sliding_window)
                                                                                               ,na_pad = FALSE)
                                                             
                                                             runner_df <- data.frame(LOC = runner_LOC,
                                                                                     N   = runner_N) %>% 
                                                               filter(!is.na(LOC))        %>%  #filter out sliding windows without outliers
                                                               filter(N >= min_snps) #filter out outliers from sliding window that do not have the minimum number of SNPs
                                                             
                                                             if(nrow(runner_df != 0)){
                                                               LG_outliers_selected <- inner_join(LG_outliers,runner_df, by = c("LOC"="LOC"))
                                                               runner_list[[paste0(i,"_runner")]] <- LG_outliers_selected} 
                                                           }
                                                           #bind dfs from list into one df
                                                           selected_outliers  <- as.data.frame(rlist::list.rbind(runner_list))
                                                           N_selected <- nrow(selected_outliers)
                                                           write_tsv(selected_outliers, file = glue::glue("{out}_selected_outliers_N{N_selected}.tsv"))
                                                           #return(selected_outliers)
                                                           
                                                         } else if(filter_method == "LD"){
                                                           print("applying LD method for selecting outliers per sliding window")
                                                           
                                                           LD_vec <- unlist(SNPRelate::snpgdsLDpruning(gdsobj =  gds,
                                                                                                       sample.id = NULL,
                                                                                                       snp.id = outliers$LocusName,
                                                                                                       autosome.only = FALSE,
                                                                                                       slide.max.bp = sliding_window,
                                                                                                       ld.threshold = ld.threshold))
                                                           
                                                           selected_outliers <- outliers %>% filter(LocusName %in% unname(LD_vec))
                                                           N_selected <- nrow(selected_outliers)
                                                           write_tsv(selected_outliers, file = glue::glue("{out}_selected_outliers_N{N_selected}.tsv"))
                                                           #return(selected_outliers)
                                                           
                                                         } else {
                                                           print("error: filter_method has to be either 'FST' or 'LD' ")
                                                         }
                                                         
                                                         #create plots if plots_out is given
                                                         if(!is.null(plots_out) & isTRUE(qval_plot)){
                                                           print("created summary plots")
                                                           #plot qvalue information
                                                           qvals <- qvalue(OutFlank$pvaluesRightTail)
                                                           
                                                           png(glue::glue("{plots_out}_qvalue_summary_plot.png"),res = 300, units = "in", width = 16, height = 10)
                                                              print(plot(qvals))
                                                           dev.off()  
                                                           
                                                           png(glue::glue("{plots_out}_qvalue_summary_histogram.png"),res = 300, units = "in", width = 16, height = 10)
                                                             print(hist(qvals))  
                                                           dev.off()   
                                                           
                                                           png(glue::glue("{plots_out}_{N_selected}SNPs_outflank_summary.png"),res = 300, units = "in", width = 16, height = 10)
                                                             par(mfrow=c(2,2))
                                                             # outlier by He
                                                             plot(OutFlank$He, OutFlank$FST, pch=20, col="grey")
                                                             points(outliers$He, outliers$FST, pch=21, col="blue")
                                                             points(selected_outliers$He, selected_outliers$FST, pch=21, col="green")
                                                             # outlier along genome
                                                             plot(OutFlank$BPcum, OutFlank$FST, pch=20, col="grey")
                                                             points(outliers$BPcum, outliers$FST, pch=21, col="blue")
                                                             points(selected_outliers$BPcum, selected_outliers$FST, pch=21, col="green")
                                                             # output fit check
                                                             OutFLANKResultsPlotter(OutFlank_list, withOutliers = TRUE,
                                                                                    NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001,
                                                                                    Zoom = FALSE, RightZoomFraction = 0.05, titletext = NULL)
                                                             # pvalue histogram
                                                             hist(OutFlank$pvaluesRightTail)
                                                           dev.off()  
                                                         }
                                                           return(selected_outliers)
                                                       }
                                                     } 
