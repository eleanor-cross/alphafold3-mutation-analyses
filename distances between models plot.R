library(dplyr)
library(readxl)
library(stringr)
library(tidyr)
library(ggplot2)
library(pals)
library(ggnewscale)

rm(list = ls())
setwd("C:/users/ezraa/onedrive/desktop/PURA")
directory <- "C:/users/ezraa/onedrive/desktop/PURA/PURA folded mutations/"
csv_files <- list.files(path = directory, pattern = "*distances.csv", full.names = TRUE, recursive = TRUE)

read_function = function(name){
  d = read.csv(name) %>% mutate(original = str_extract(name, "[^/]+$") %>% gsub(pattern = '.csv', replacement = ''))
  return(d)
}

aas <- c("ala", "arg", "asn", "asp", "cys", "glu", "gln", "gly", "his", "ile", "leu", "lys", "met", "phe", "pro", "ser", "thr", "trp", "tyr", "val")

distances <- csv_files %>% lapply(read_function) %>% bind_rows() %>%
  mutate(name = sapply(FUN = gsub, X = original, pattern = "_between.*|pura_|fold_", replacement = ""))


ggplot() + 
  geom_tile(data = distances, aes(y = name, x = resi, fill = average_distances)) + 
  scale_fill_gradientn(colors = pals::coolwarm(100)) +  
  theme_minimal() + 
  theme(legend.position = 'bottom', legend.title.position = 'top', panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
