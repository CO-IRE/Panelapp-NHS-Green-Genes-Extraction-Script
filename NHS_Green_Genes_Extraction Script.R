library(rvest)
library(purrr)
library(tidyverse)
library(dplyr)

update.packages()
website <- "https://panelapp.genomicsengland.co.uk/panels/"
page <- read_html(website)

c_ref <- page %>%
  html_nodes("a") %>% # find all links
  html_attr("href")

df_ref <- tibble(ref = c_ref, id = NA) %>%
  filter(str_detect(ref, 'download')) %>%
  mutate(ref = str_remove(ref, '/panels/')) %>%
  mutate(id = ref) %>%
  mutate(id = str_remove(id, '/download/01234/'))

# Linux command - if you are using Windows, please make sure that you create a new folder with the name 'gene_panel'
# and remove the "system('mkdir gene_panel')" line
#system('mkdir gene_panel')
#setwd('gene_panel')

walk2(df_ref$ref, df_ref$id, function(a, b)
  download.file(url = paste0(website, a), destfile = paste0('gene_panel_', b))
)

#Here
files_panel <- list.files()

As <- as.data.frame(files_panel)

for (i in files_panel)
  
  
{
  panel_total_i <- files_panel[i:i+1] %>% map_dfr(~ read_tsv(.x) %>% 
                                                select(`Entity Name`, `Entity type`, `Gene Symbol`, `Sources(; separated)`, 
                                                       Level4, Phenotypes) %>% 
                                                mutate(source = .x) 
  
}

panel_total <- files_panel %>% map_dfr(~ read_tsv(.x)) %>% 
                                         select(`Entity Name`, `Entity type`, `Gene Symbol`, `Sources(; separated)`, 
                                                Level4, Phenotypes) %>% 
                                         mutate(source = .x) 

panel_total2 <- files_panel[1:5] %>% map_dfr(~ read_tsv(.x) %>% 
                                              select(`Entity Name`, `Entity type`, `Gene Symbol`, `Sources(; separated)`, 
                                                     Level4, Phenotypes) %>% 
                                              mutate(source = .x) )

panel_total3 <- files_panel[6:8] %>% map_dfr(~ read_tsv(.x) %>% 
                                               select(`Entity Name`, `Entity type`, `Gene Symbol`, `Sources(; separated)`, 
                                                      Level4, Phenotypes) %>% 
                                               mutate(source = .x) )


length(files_panel)

install.packages("installr"); library(installr)
updateR()

  
  for (k in length(files_panel)) { 
    
  
k %>% map_dfr(~ read_tsv(.x)) %>% 
    select(`Entity Name`, `Entity type`, `Gene Symbol`, `Sources(; separated)`, 
           Level4, Phenotypes)# %>% 
  #  mutate_(source = .x) 

  
}


y <- NULL

for (k in length(files_panel)) {
  
tmp <-  k %>% map_dfr(~ read_tsv(.x)) %>% 
    select(`Entity Name`, `Entity type`, `Gene Symbol`, `Sources(; separated)`, 
           Level4, Phenotypes)# %>% 
  #  mutate_(source = .x) 


  
}



Attempt <- as.data.frame(files_panel)

is.vector(files_panel)

result<-rbind(panel_total, panel_total_2)

nrow(Attempt)

for(row in 1:nrow(Attempt)) {
   
  map_dfr(~ read_tsv(.x)) %>% 
    select(`Entity Name`, `Entity type`, `Gene Symbol`, `Sources(; separated)`, 
           Level4, Phenotypes) %>% 
    mutate(source = .x) 
  }


#Do five at a time and then add their rows to panel_total
#Could [] files_panel and add + to panel_total

# Filtering out genes with a evidence level (red - amber)
panel_total <- panel_total %>%
  rename(entity_name = `Entity Name`,
         entity_type = `Entity type`,
         gene = `Gene Symbol`,
         sources = `Sources(; separated)`) %>%
  filter(entity_type == 'gene') %>%  # optional - we can include regions in our analysis
  filter(str_detect(sources, 'Expert Review Green')) %>%
  select(gene, Level4, -sources, source, Phenotypes)

setwd('..')

write_tsv(panel_total, 'panel_genes.tsv')
