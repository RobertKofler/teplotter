# Rscript ultra_simple.R Sophie
library(tidyverse)  

#args <- commandArgs(trailingOnly = TRUE)
#if (length(args) == 0) {
 # cat("Please provide a file\n")
  #quit("no", 1)
#}
#file<-args[1]
file<-"/Users/robertkofler/gh/teplotter/test/mdg1#LTR_Gypsy_te"

lines <- read_lines(file)|>
  str_trim() |>
  discard(~ .x == "" | str_starts(.x, fixed("#"))) |>   # skip empty & comments
  str_split("\\s+")  |> keep(~ length(.x) >= 3)

# mdg1#LTR/Gypsy_te	test	ambcov	7478	0.0
# mdg1#LTR/Gypsy_te	test	ambcov	7479	0.0
coverage <- lines |>
  keep(~ .x[3] == "cov") |>
  map_dfr(~ tibble(seqid = .x[1], sampleid = .x[2], feature = .x[3], pos=as.double(.x[4]), cov =  as.double(.x[5])))

ambcoverage <- lines |>
  keep(~ .x[3] == "ambcov") |>
  map_dfr(~ tibble(seqid = .x[1], sampleid = .x[2], feature = .x[3],pos=as.double(.x[4]), ambcov =  as.double(.x[5])))

#mdg1#LTR/Gypsy_te	test	snp	7102	T	1.22	0	0.07	0.0
#mdg1#LTR/Gypsy_te	test	snp	7128	A	0	0.0	0.0	4.01
snp <- lines |>
  keep(~ .x[3] == "snp") |>
  map_dfr(~ tibble(seqid = .x[1], sampleid = .x[2], feature = .x[3], pos =as.double(.x[4]), refc=.x[5], ac= as.double(.x[6]), tc= as.double(.x[7]),cc= as.double(.x[8]),gc= as.double(.x[9])))

# mdg1#LTR/Gypsy_te	test	ins	7415	7417	33.78	32.49	1.26
# mdg1#LTR/Gypsy_te	test	ins	7458	7460	15.04	13.71	0.99
insertion <- lines |>
  keep(~ .x[3] == "ins") |>
  map_dfr(~ tibble(seqid = .x[1], sampleid = .x[2], feature = .x[3], start =  as.double(.x[4]),end =  as.double(.x[5]),startcov =  as.double(.x[6]),encov =  as.double(.x[7]),count=  as.double(.x[8])))

# mdg1#LTR/Gypsy_te	test	del	7274	3	1.16
deletion <- lines |>
  keep(~ .x[3] == "del") |>
  map_dfr(~ tibble(seqid = .x[1], sampleid = .x[2], feature = .x[3], pos =  as.double(.x[4]),length =  as.double(.x[5]),count =  as.double(.x[6])))

theme_set(theme_bw())
plo<-ggplot()+
  geom_polygon(data = coverage, mapping = aes(x = pos, y = cov), fill = 'grey', color = 'grey') +
  geom_polygon(data = ambcoverage, aes(x = pos, y = ambcov), fill = 'lightgrey', color = 'lightgrey')
plot(plo)
  
expr_spli_plot <- ggplot() +
  geom_polygon(data = s, mapping = aes(x = pos, y = cov), fill = 'grey', color = 'grey') +
  geom_polygon(data = as, aes(x = pos, y = -cov), fill = 'lightgrey', color = 'lightgrey') +
  geom_curve(data = a_s, mapping = aes(x = start, y = cov.x, xend = end, yend = cov.y, linewidth = size_scaled), curvature = -0.15, ncp = 10, show.legend = FALSE) +
  geom_curve(data = a_as, mapping = aes(x = start, y = -cov.x, xend = end, yend = -cov.y, linewidth = size_scaled), curvature = 0.15, ncp = 10, show.legend = FALSE) +
  geom_hline(yintercept = 0, color = "black") + # Add black line at Y = 0
  facet_grid(time ~ rep) +
  scale_size(range = c(0.1, 1)) + # Adjust the range for line thickness
  coord_cartesian(ylim = c(min(a$cov) -1, max(a$cov) + 1)) + # Adjust y-axis limits
  xlab("position") +
  ylab("expression [rpm]") +
  theme(
    panel.grid.major = element_blank(), # Remove major gridlines
    panel.grid.minor = element_blank(), # Remove minor gridlines
    panel.border = element_rect(color = "black", fill = NA, size = 1), # Add border around each plot
    strip.background = element_rect(fill = "grey80", color = "black"), # Shaded grey box for facet titles
    strip.text = element_text(size = 10),
    panel.spacing = unit(0.4, "lines"), # Reduce space between plots
    legend.position = "bottom", # Move legend to the bottom
    legend.direction = "horizontal" # Make the legend horizontal
  )
