# for github.com/davemcg/OGVFB_RNAseq

library(tidyverse)
si <- read_tsv('data/sample_info2.tsv')

si %>% filter(!is.na(Run), Organism == 'Human') %>% select(Sample, Run) %>%
  mutate(R1 = glue::glue('{Run}_1.fastq.gz'),
         R2 = glue::glue('{Run}_2.fastq.gz')) %>%
  select(-Run) %>%
  pivot_longer(cols = c(R1,R2)) %>%
  select(-name) %>%
  write_tsv('~/Desktop/x.tsv')


si %>% filter(!is.na(Run), Organism == 'Mouse') %>%
  dplyr::select(Sample, Run) %>%
  mutate(R1 = glue::glue('{Run}_1.fastq.gz'),
         R2 = glue::glue('{Run}_2.fastq.gz')) %>%
  dplyr::select(-Run) %>%
  pivot_longer(cols = c(R1,R2)) %>%
  dplyr::select(-name) %>%
  write_tsv('~/Desktop/x2.tsv')
