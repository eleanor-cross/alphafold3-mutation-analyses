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
csv_files <- list.files(path = directory, pattern = "*.csv", full.names = TRUE, recursive = TRUE)
csv_files = csv_files[grepl("between_models", csv_files)==FALSE]

read_function = function(name){
  d = read.csv(name) %>% mutate(original = str_extract(name, "[^/]+$") %>% gsub(pattern = '.csv', replacement = ''))
  return(d)
}
distances <- csv_files %>% lapply(read_function) %>% bind_rows() %>%
  mutate(name = sapply(FUN = gsub, X = original, pattern = "_model.*|pura_|fold_", replacement = "")) 

distances = distances %>% select(-c(x,y,z))

list = unique(distances$original)

fix = function(x){
  return(gsub(x, pattern = "\\[|\\]", replacement = "") %>% as.numeric() %>% abs())
}

domains = data.frame(matrix(nrow = 322, ncol = 0))
domains$resi = c(1:322)
domains$Domain = NA
domains$Domain[c(60:125)] = "PUR I"
domains$Domain[c(142:213)] = "PUR II"
domains$Domain[c(215:281)] = "PUR III"
domains = domains[!is.na(domains$Domain),]

I_color = "#FF0000"
II_color = "#046C9A"
III_color = "#00A08A"


make_plot = function(i){
  one = distances[distances$original==list[i],]
  one = suppressWarnings(one %>%
  mutate(across(c('X_8chw_distances','X_8cht_distances','fold_wt_pura_model_0_distances'),
                ~fix(`.`),
                .names = 'd_{.col}')) %>%
    mutate(ind = c(1:nrow(one))))
  missing = one[is.na(one$d_fold_wt_pura_model_0_distances),] %>%
    mutate(Termination = "")
  
  p = suppressWarnings(ggplot() + 
  geom_tile(data = one,
            aes(x = resi, 
                y = 3,
                fill = d_fold_wt_pura_model_0_distances)) +
  geom_tile(data = one,
            aes(x = resi, 
                y = 5, 
                fill = d_X_8chw_distances)) +
  geom_tile(data = one, 
            aes(x = resi, 
                y = 7, 
                fill = d_X_8cht_distances)) +
  scale_fill_gradientn(name = "Distance", colours = coolwarm(100), na.value = 'transparent')+
  theme_minimal() +
  coord_equal(ratio = 40) + 
  ggnewscale::new_scale_fill() + 
  geom_tile(data = domains,
            aes(
            x = resi, 
            y = 1, 
            fill = Domain)
  ) +
    scale_fill_manual(values = c(I_color,II_color,III_color)) +
  xlim(0,325) + 
  ylim(0,8) +
  theme(legend.position = 'bottom',
        legend.title.position = 'top',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_y_continuous(breaks = c(1,3,5,7),
                     labels = c("Domains", "WT model 0", "8chw", "8cht")) +
  xlab("Residue") + 
  ylab("") +
  coord_fixed(ratio = 25))
  
  if(nrow(missing)>0){
    p = p + 
      ggnewscale::new_scale_fill() + 
      geom_tile(data = missing,
                aes(x = ind, 
                    y = 3,
                    fill = Termination))+
      scale_fill_manual(values = "black") 
  }
  name = list[i]
  folder = gsub(pattern = "_model.*", replacement = "", name)
  filename = paste0(directory,folder,'/',name,'/',name,'_distance_plots.png')
  
  suppressWarnings(ggsave(filename, plot = p, width = 8, height = 8, dpi = 300))
  print(paste('Wrote ', filename))
}

iterate = c(1:length(list))
sapply(FUN = make_plot, iterate)


