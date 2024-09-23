library(dplyr)
library(readxl)
library(stringr)
library(tidyr)
library(ggplot2)
library(pals)
library(ggnewscale)

install.packages("vsr
")

rm(list = ls())
setwd("C:/users/ezraa/onedrive/desktop/PURA")
directory <- "C:/users/ezraa/onedrive/desktop/PURA/PURA folded mutations/"
csv_files <- list.files(path = directory, pattern = "*.csv", full.names = TRUE, recursive = TRUE)

csv_files
