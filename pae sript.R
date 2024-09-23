library(jsonlite)
library(dplyr)
rm(list = ls())
d = read_json("C:\\Users\\ezraa\\OneDrive\\Desktop\\PURA\\PURA folded mutations\\fold_pura_lys160glu\\fold_pura_lys160glu_full_data_0.json")

d = d$..JSON 

pae = d[[1]][['pae']]

paedf = data.frame(matrix(nrow = 322,ncol = 0))

i = 1
for(i in 1:length(pae)){
  paedf[,paste0('resi',i)] = unlist(pae[i])
}
rm(i)

rm(pae,d,p)

# paem = paedf %>% as.matrix()
# heatmap(paem)

library(tidyr)

pael = paedf %>% 
  mutate(`Scored Residue` = 1:nrow(paedf)) %>%
  pivot_longer(cols = -`Scored Residue`) %>% 
  mutate(`Aligned Residue` = gsub(pattern = "[A-Za-z]", replacement = "", name)) 

pael = pael %>% rename(`Expected Position Error (Ångströms)`=value)

rm(paedf)

order = unique(pael$`Aligned Residue`) %>% rev()

pael$`Aligned Residue` = factor(pael$`Aligned Residue`, levels = order, ordered = TRUE)


p = ggplot() + 
  geom_tile(data = pael, 
            aes(x = `Scored Residue`, 
                y = `Aligned Residue`, 
                fill = `Expected Position Error (Ångströms)`)) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(),
        legend.position = 'bottom'# Remove minor gridlines
        ) +
  scale_y_discrete(breaks = function(y) y[c(1,seq(23, 300, by = 25),322)]) +
  scale_x_continuous(breaks = c(seq(0,300,by = 25),322)) 
  

windows()
p


