library(tidyverse)

GenStats="/bigdata/littlab/arajewski/FULTranscriptomes/multiqc_data/multiqc_general_stats.txt"
RSeQC="/bigdata/littlab/arajewski/FULTranscriptomes/multiqc_data/multiqc_rseqc_read_distribution.txt"
QCNames="/bigdata/littlab/arajewski/FULTranscriptomes/QC_Names.csv"

read_tsv(RSeQC,
                  col_types = cols()) %>% 
  select(1, 28:31,) %>% 
    mutate(Sample = gsub("RSeQC_", "", Sample)) %>%
  mutate(PctExon = rowSums(across(2:4), na.rm=TRUE),
         PctIntron = introns_tag_pct,
         PctOther = (100 - (PctExon + PctIntron))) %>%
  inner_join(read_tsv(GenStats),
             by = c("Sample")) %>%
  inner_join(read_csv(QCNames),
             by = c("Sample")) %>% 
  select(c(16:19,15,14,9,6:8)) %>% 
  rename(MeanCov = `QualiMap_mqc-generalstats-qualimap-mean_coverage`,
         PctMap = `STAR_mqc-generalstats-star-uniquely_mapped_percent`,
         NumMap = `STAR_mqc-generalstats-star-uniquely_mapped`) %>% 
  mutate(Stage = factor(Stage, levels = c("1", "2", "3", "Tr","4"))) %>% 
  arrange(Species, Stage, Replicate) %>% 
  write_csv("Tables/MappingStatistics.csv")
