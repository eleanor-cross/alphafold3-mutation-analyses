library(dplyr)
library(readxl)
library(stringr)
library(tidyr)
library(ggplot2)
library(pals)
library(ggnewscale)
library(data.table)

rm(list = ls())
setwd("C:/users/ezraa/onedrive/desktop/PURA")
directory <- "C:/users/ezraa/onedrive/desktop/PURA/PURA folded mutations/"
csv_files <- list.files(path = directory, pattern = "*.csv", full.names = TRUE, recursive = TRUE)
csv_files = csv_files[grepl("between_models", csv_files)==FALSE]

read_function = function(name){
  d = fread(name) %>% mutate(original = str_extract(name, "[^/]+$") %>% gsub(pattern = '.csv', replacement = ''))
  return(d)
}
d <- csv_files %>% lapply(read_function) %>% bind_rows() %>%
  mutate(name = sapply(FUN = gsub, X = original, pattern = "_model.*|pura_|fold_", replacement = "")) 

d = d %>% select(-c(x,y,z))

fix = function(x){
  return(gsub(x, pattern = "\\[|\\]", replacement = "") %>% as.numeric() %>% abs())
}

d = d %>% rename(`8chw` = `_8chw_distances`, `8cht` = `_8cht_distances`, wt = `fold_wt_pura_model_0_distances`) %>% 
  mutate(suppressWarnings(mutate(across(c('8chw','8cht','wt'), ~fix(`.`),.names = 'd_{.col}'))))

d$type = NA
d = d %>% mutate(type = ifelse(grepl(pattern = "fs", original)==TRUE, 
                               "fs",
                               ifelse(grepl("ter|Ter|x|X", original)==TRUE,
                                      "ns",
                                      "ms"
                                      )
                               )
                 )


d = d %>% filter(grepl("fold_leu54cysfster24_", original)==FALSE & 
                   grepl("fold_lys97ter_", original)==FALSE & 
                   grepl("fold_phe73ser", original)==FALSE
)

types = d %>% select(c(original, name, type)) %>% filter(grepl("leu114del", name)==FALSE) %>% unique() %>%  mutate(mutation_location = gsub("[^0-9]", "", name)) %>%  group_by(type) %>% arrange(mutation_location) %>% mutate(index = row_number()) %>% ungroup()


d2 = left_join(d[grepl("leu114del",d$name)==FALSE,],types)
# now have an index column 
  
range(d2$d_8chw, na.rm = T)
range(d2$d_8cht, na.rm = T)
range(d2$d_wt, na.rm = T)

d2[grepl('leu114del',d2$name)==TRUE,]

testing_plot = ggplot() + geom_tile(data = d2, 
                     aes(x = resi, 
                         y = index, 
                         fill = d_wt)) + 
  scale_fill_gradientn(colors = coolwarm(100),
                       values = scales::rescale(seq(0,max(d2$d_8cht, na.rm = T), length.out = 100))) + 
  facet_grid(~type)
testing_plot

theme = theme(legend.position = 'bottom',
              legend.title.position = 'top',
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_rect(fill = 'white'),
              axis.title.y = element_blank(), 
              axis.ticks.y = element_blank(),
              axis.text.y = element_blank(),
              axis.title.x = element_blank())

tw = types[types$index %in% seq(3,max(types$index), by = 5),] %>% pivot_wider(names_from = type, values_from = index) 

d2 = d2 %>% group_by(name) %>% 
  mutate(ymin = min(index,na.rm = TRUE)-.5,
         ymax = max(index, na.rm = TRUE)+.5,
         xmin = min(resi, na.rm = TRUE)-.5,
         xmax = max(resi, na.rm = TRUE)+.5)

rect_df = d2 %>% select(c(name,type, xmin,xmax,ymin,ymax)) %>% unique()

ms = ggplot() + geom_tile(data = d2[d2$type=="ms",],
                          aes(x = resi, 
                              y = index,
                              fill = d_wt)) + 
  scale_fill_gradientn(name = "WT Distance", 
                       colors = coolwarm(100), 
                       values = scales::rescale(seq(0,max(d2$d_8cht, na.rm = T), 
                                                    length.out = 100))) +
  xlim(-50,330) +
  geom_text(data = tw, size = 3,
            aes(y = ms, x=0, label = name), hjust = 1) + geom_rect(data = rect_df[rect_df$type=="ms",],
            color = "black", 
            fill = "transparent",
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) +
  theme

ns = ggplot() + geom_tile(data = d2[d2$type=="ns",],
                     aes(x = resi, 
                         y = index,
                         fill = d_wt)) + 
  scale_fill_gradientn(name = "WT Distance", 
                       colors = coolwarm(100), 
                       values = scales::rescale(seq(0,max(d2$d_8cht, na.rm = T), 
                                                    length.out = 100))) +
  xlim(-50,330) +
  geom_text(data = tw, size = 3,
            aes(y = ns, x=0, label = name), hjust = 1) + geom_rect(data = rect_df[rect_df$type=="ns",],
                                                                   color = "black", 
                                                                   fill = "transparent",
                                                                   aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) +
  theme

fs = ggplot() + geom_tile(data = d2[d2$type=="fs",],
                          aes(x = resi, 
                              y = index,
                              fill = d_wt)) + 
  scale_fill_gradientn(name = "WT Distance", 
                       colors = coolwarm(100), 
                       values = scales::rescale(seq(0,max(d2$d_8cht, na.rm = T), 
                                                    length.out = 100))) +
  xlim(-50,330) +
  geom_text(data = tw, size = 3,
            aes(y = fs, x=0, label = name), hjust = 1) + geom_rect(data = rect_df[rect_df$type=="fs",],
                                                                   color = "black", 
                                                                   fill = "transparent",
                                                                   aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) +
  theme



library(ggpubr)
ggarrange(ms,ns,fs, nrow = 1, common.legend = TRUE)


# View(d[grepl("73ser",name)==TRUE],)
rect_df[rect_df$type=="ms",]
# dm = d %>% select(c(resi, d_wt, original)) %>% pivot_wider(names_from = resi, values_from = d_wt) %>% select(-original) %>% as.matrix()

unique(d$name[d$type=="ms"])


