# Rscript visualize-plotable.R input.plotable output.png
library(tidyverse)  

args <- commandArgs(trailingOnly = TRUE)
if (length(args) <2) {
 cat("Please provide an input file and an output file\n"+
       "Usage Rscript visualize-plotable.R input.plotable output.png")
  quit("no", 1)
}


#
# some parameters; feel free to modify
#
mindeletion=10 # minimum length of the internal deletions
width=8       # plot width 
height=5      # plot height
dpi=300       # plot dpi
#
# end parameters
#

# get parameters
file<-args[1]
outfile<-args[2]

# debug 
file<-"/Users/robertkofler/gh/teplotter/test2/mdg1#LTR_Gypsy_te"




data <- read_tsv(file,col_names = FALSE,cols(.default = col_character()))


# split of coverage
cov <- data |> filter(str_detect(X3, "cov"))
cov <- cov |> rename(seqid=X1,sampleid=X2,feature=X3,pos=X4,cov=X5) 
cov <- cov |> mutate(pos = as.double(pos),cov= as.double(cov))



# split of ambcoverge
ambcov <- data |> filter(str_detect(X3, "ambcov"))
ambcov <- ambcov |> rename(seqid=X1,sampleid=X2,feature=X3,pos=X4,cov=X5) 
ambcov <- ambcov |>mutate(pos = as.double(pos),cov= as.double(cov))

# split of snps
snp <- data |> filter(str_detect(X3, "snp"))
snp <- snp |> rename(seqid=X1,sampleid=X2,feature=X3,pos=X4,refc=X5,base=X6,count=X7) 
snp <- snp |>  mutate(pos = as.double(pos),count= as.double(count))


# split of deletion
deletion <- data |> filter(str_detect(X3, "del"))
deletion <- deletion |> rename(seqid=X1,sampleid=X2,feature=X3,start=X4,end=X5,startcov=X6,endcov=X7,count=X8) 
deletion <- deletion |> mutate(start = as.double(start),end= as.double(end),startcov = as.double(startcov),endcov= as.double(endcov),count= as.double(count))

# split of insertion
insertion <- data |> filter(str_detect(X3, "ins"))
insertion <- insertion |> rename(seqid=X1,sampleid=X2,feature=X3,pos=X4,length=X5,count=X6) 
insertion <- insertion |> mutate(pos = as.double(pos), length= as.double(length), count= as.double(count))

# prepare insertions
# filter min size of insertion
deletion<- deletion |> filter(end-start>mindeletion)
# size of scaling
deletion$scale=log(deletion$count)

theme_set(theme_bw())
plo<-ggplot()+
  geom_polygon(data = coverage, mapping = aes(x = pos, y = cov), fill = 'grey', color = 'grey') +
  geom_polygon(data = ambcoverage, aes(x = pos, y = ambcov), fill = 'lightgrey', color = 'lightgrey')+
  geom_curve(data = deletion, mapping = aes(x = start, y = startcov, xend = end, yend = endcov, linewidth = scale),  curvature = -0.15, ncp=5,show.legend = FALSE)+
  scale_linewidth(range = c(0.3, 2))+xlab("position") + ylab("coverage")+
  geom_bar(data=snp,aes(x=pos,y=count,fill=base),stat="identity",width=2)+
  geom_bar(data=insertion,aes(x=pos,y=count),stat="identity",color="grey50",width=4)

# faceting
nseq<-n_distinct(cov$seqid)
nsample<-n_distinct(cov$sampleid)
if (nseq > 1 & nsample>1) {
  plo<-plo+facet_grid(seqid~sampleid)
} else if (nseq>1){
  plo<-plo+facet_grid(seqid~.)
}else if (nsample>1){
  plo<-plo+facet_grid(.~sampleid)
}

ggsave(outfile, plot = plo, width = width, height = height, dpi = dpi)
plot(plo)

