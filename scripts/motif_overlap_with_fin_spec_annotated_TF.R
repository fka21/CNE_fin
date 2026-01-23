library(biomaRt)

tbl1 <- read.csv("~/Downloads/zfin_search_results.csv")
tbl2 <- read_tsv("~/Downloads/Danio_rerio_TF.txt")
tbl3 <- read.csv("~/Downloads/DatabaseExtract_v_1.01.csv")

ensembl <- useMart("ensembl", dataset = "drerio_gene_ensembl")

# Query the genes using biomaRt
genes_info <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description", "hsapiens_homolog_ensembl_gene"),
                    filters = "external_gene_name",
                    values = tbl1$name,
                    mart = ensembl)

fin_tf <- filter(genes_info, genes_info$ensembl_gene_id %in% tbl2$Ensembl)
fin_tf$HGNC <- tbl3$HGNC.symbol[match(fin_tf$hsapiens_homolog_ensembl_gene, tbl3$Ensembl.ID)]

motif <- read_tsv("~/Documents/Projects/CNE_fin/downstream_analysis/enriched_motifs.tsv")

motif %>% 
  mutate(
    fin_spec = case_when(
      motif_alt_id %in% fin_tf$HGNC ~ "Fin specific TF",
      TRUE ~ "Non-specific"
    )
  ) %>% 
  group_by(fin_spec, id) %>% 
  summarise(count = n()) %>%
  mutate(percentage = count / sum(count)) %>%
  ggplot(aes(x = fin_spec, y = percentage, fill = id)) +
  geom_col(position = 'fill', color = "black", alpha = 0.9) +
  geom_text(aes(label = scales::percent(percentage), group = interaction(fin_spec, id)),
            position = position_fill(vjust = 0.5), color = "black") +
  scale_y_continuous(labels = scales::percent_format(scale = 100)) +  # Display percentages on the y-axis
  theme_bw() +
  scale_fill_manual(values = c("#1f78b4", "#a6cee3", "#ff7f00", "#fdbf6f")) +
  labs(x = NULL, y = "Percentage", fill = "ID")

ggsave("~/Documents/Projects/CNE_fin/downstream_analysis/fin_specific_enriched_motifs.png",
       units = 'in', width = 6, height = 5, dpi = 720)


# Print the result
print(genes_info)
