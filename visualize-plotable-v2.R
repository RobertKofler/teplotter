# Rscript ultra_simple.R Sophie
library(tidyverse)  

#args <- commandArgs(trailingOnly = TRUE)
#if (length(args) == 0) {
 # cat("Please provide a file\n")
  #quit("no", 1)
#}
#file<-args[1]
file<-"/Users/robertkofler/gh/teplotter/test2/mdg1#LTR_Gypsy_te"
mininsertion=5

data <- read_tsv(file,col_names = FALSE,cols(.default = col_character()))

# split of coverage
cov <- data |> filter(str_detect(X3, "cov"))
cov <- cov |> rename(seqid=X1,sampleid=X2,feature=X3,pos=X4,cov=X5) 
cov <- cov |>mutate(pos = as.double(pos),cov= as.double(cov))


# split of ambcoverge
ambcov <- data |> filter(str_detect(X3, "ambcov"))
ambcov <- ambcov |> rename(seqid=X1,sampleid=X2,feature=X3,pos=X4,cov=X5) 
ambcov <- ambcov |>mutate(pos = as.double(pos),cov= as.double(cov))

# split of snps
snp <- data |> filter(str_detect(X3, "snp"))
snp <- snp |> rename(seqid=X1,sampleid=X2,feature=X3,pos=X4,refc=X5,ac=X6,tc=X7,cc=X8,gc=X9) 
snp <- snp |>  mutate(pos = as.double(pos),ac= as.double(ac),tc= as.double(tc),cc= as.double(cc),gc= as.double(gc))

# split of insertions
insertion <- data |> filter(str_detect(X3, "ins"))
insertion <- insertion |> rename(seqid=X1,sampleid=X2,feature=X3,start=X4,end=X5,startcov=X6,endcov=X7,count=X8) 
insertion <- insertion |> mutate(start = as.double(start),end= as.double(end),startcov = as.double(startcov),endcov= as.double(endcov),count= as.double(count))

# split of deletion
deletion <- data |> filter(str_detect(X3, "del"))
deletion <- deletion |> rename(seqid=X1,sampleid=X2,feature=X3,pos=X4,length=X5,count=X6) 
deletion <- deletion |> mutate(pos = as.double(pos), length= as.double(length), count= as.double(count))

# prepare insertions
# filter min size of insertion
insertion<- insertion |> filter(end-start>mininsertion)
# size of scaling
insertion$scale=log(insertion$count)

theme_set(theme_bw())
plo<-ggplot()+
  geom_polygon(data = coverage, mapping = aes(x = pos, y = cov), fill = 'grey', color = 'grey') +
  geom_polygon(data = ambcoverage, aes(x = pos, y = ambcov), fill = 'lightgrey', color = 'lightgrey')+
  geom_curve(data = insertion, mapping = aes(x = start, y = startcov, xend = end, yend = endcov, linewidth = scale),  curvature = -0.15, ncp=5,show.legend = FALSE)+
  scale_linewidth(range = c(0.3, 2))+xlab("position") + ylab("coverage")
plot(plo)

