# Author
# Tom Oosting
# 2025-05-10

#load packages
library(optparse)
library(glue)
library(ggplot2)
library(tidyverse)

option_list <- list(make_option(c("--ld_file")   , type = "character", default = NULL  , metavar = "character", help = ""), 
                    make_option(c("--bin_size")  , type = "integer"  , default = 100   , metavar = "integer"  , help = ""), 
                    make_option(c("--r2_thresh") , type = "double"   , default = 0.2   , metavar = "double"   , help = ""), 
                    make_option(c("--thin_size") , type = "integer"  , default = 10000 , metavar = "integer"  , help = "")) 

# Parse the command-line options to Global environment
list2env(parse_args(OptionParser(option_list = option_list)), envir = .GlobalEnv)

# insert own path the file
# ld_file <- "C:/Users/oostinto/Desktop/to_repo/snapper_norm_neutral_plink.ld"

#set output file extension
ld_out <- str_remove(ld_file,".ld$")

# load dataframe
decay_df <- read_tsv(file = ld_file) %>% mutate(dis = abs(BP_A-BP_B))

#QQplot
png(filename = glue("{ld_out}_qqnorm.png"))
  qqnorm(decay_df$R2)  
dev.off()

#create dataframe
decay_binned <- decay_df %>%  group_by(group = cut(dis, breaks = seq(0, max(dis), bin_size))) %>%  summarise(R2 = mean(R2))
sd           <- decay_df %>%  group_by(group = cut(dis, breaks = seq(0, max(dis), bin_size))) %>%  summarise(sd = sd(R2))
quantiles    <- decay_df %>%  group_by(group = cut(dis, breaks = seq(0, max(dis), bin_size))) %>%  summarise(quantiles = quantile(R2, probs = c(0.95)))
decay_binned$sd <- sd$sd
decay_binned$bin <- seq(bin_size/2,10^6,bin_size)[1:nrow(decay_binned)]
decay_binned <- decay_binned %>% mutate(min = R2-sd , max = R2+sd) %>% filter(!is.na(sd) & bin < 20000)

#plot
ggplot(decay_binned)+
  geom_ribbon(aes(x=bin,ymin=min,ymax=max),fill = "lightblue", alpha=0.5)+
  geom_point(aes(x=bin,y=R2))+
  geom_hline(yintercept = r2_thresh, linetype = "dashed", color = "red", size = 1)+
  geom_vline(xintercept = thin_size, linetype = "dashed", color = "black")+
  xlab("Distance (bp)") + 
  ylab(expression(italic(r)^2))+
  theme_bw()
ggsave(filename = glue("{ld_out}.png"))