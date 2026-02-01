# Rscript ultra_simple.R Sophie
library(tidyverse)   # or just stringr + purrr

#args <- commandArgs(trailingOnly = TRUE)
#if (length(args) == 0) {
 # cat("Please provide a file\n")
  #quit("no", 1)
#}
#file<-args[1]
file<-"/Users/robertkofler/dev/data/Dmel01.normalized"


lines <- readLines(file, warn = FALSE)

data_raw <- lines |>
  str_split("[\t|]") 

snp_data <- data_raw |>
  keep(~ length(.x) >= 2) |>
  keep(~ str_detect(.x[[2]], "snp")) |>
  map_dfr(~ as_tibble_row(set_names(.x, paste0("V", seq_along(.x))))) |>
  mutate(across(everything(), str_trim))          # optional: clean whitespace

