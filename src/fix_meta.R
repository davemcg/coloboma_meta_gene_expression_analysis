meta <- read_tsv('~/git/coloboma_meta_gene_expression_analysis/data/t_Info.tsv')
ERP123356 <- read_tsv('data/ERP123356.txt')
ERP123357 <- read_tsv('data/ERP123357.txt')
erp <- bind_rows(ERP123356, ERP123357)

meta <- meta %>% mutate(RawCore = str_extract(RawFile, '^[a-zA-Z]+\\d+') %>% gsub('_','',.))
meta_erp <- meta %>%
  left_join(erp %>% select(RawCore = run_accession, sample_accession)) %>%
  filter(grepl('ERR',RawFile)) %>%
  filter(!grepl('_2',RawFile)) %>%
  mutate(Sample = sample_accession)

meta_corrected <- bind_rows(meta %>%
                              filter(!grepl('ERR',RawFile)),
                            meta_erp) %>%
  mutate(Sample = case_when(is.na(Sample) ~ RawCore,
                            TRUE ~ Sample)) %>%
  select(-sample_accession)

# fix the chick data
meta_chick <- meta_corrected %>%
  filter(Organism == 'Chick') %>%
  left_join(SRP079981 %>%
              select(-Organism), by = c('Sample' = 'GEO_Accession (exp)')) %>%
  select(RawFile:Run) %>%
  mutate(RawFile = Run) %>%
  select(-RawFile, -RawCore)

meta3 <- bind_rows(meta_corrected %>%
                     filter(Organism != 'Chick') %>%
                     mutate(Run = RawCore),
                   meta_chick) %>%
  select(Run, Sample:Other) %>%
  mutate(Run = case_when(!grepl('GSM', Run) ~ Run)) %>%
  mutate(Layout = case_when(Accession == 'GSE84916' ~ 'Single',
                            Accession == 'GSE159822' ~ 'Paired',
                            grepl('MTAB', Accession) ~ 'Paired'))
write_tsv(meta3, file = 'data/sample_info2.tsv')
