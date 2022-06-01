library(GO.db)
library(org.Hs.eg.db)
library(tidyverse)
library(scales)
library(readxl)
library(here)
library(cowplot)
library(ggrepel)
library(limma)
library(edgeR)

source(here("lib/contrast_scatter_plot.R"))
source(here("lib/limma_test.R"))

data_file <- here("data/all_crapome_data.xlsx")

out_dir <- here("results/prober-MS", Sys.Date())
dir.create(out_dir, recursive = T)


contrast_dir <- file.path(out_dir, "prober-contrasts")
dir.create(contrast_dir, recursive = T)

sheets <- excel_sheets(data_file)
data <- map(set_names(sheets), ~ read_excel(data_file, .))

go_terms <- unlist(Term(GOTERM))
chromatin_terms <- map(c("chromatin", "DNA binding", "DNA-templated",
                         "DNA methylation", "DNA demethylation",
                         "DNA helicase", "DNA repair", "DNA replication",
                         "histone", "nucleosome",
                         "transcription"),
                       ~ str_subset(go_terms, .)) %>% unlist()

ms_protein_list <- map(data, ~ .$PROTID) %>% unlist() %>% unique() %>%
  stringi::stri_replace_all(regex = "(C.+)ORF(\\d+)", replacement = "$1orf$2")

ms_proteins_symbols <-
  AnnotationDbi::select(org.Hs.eg.db, keys = ms_protein_list, keytype = "SYMBOL", columns = c("SYMBOL", "ENTREZID")) %>%
  group_by(SYMBOL) %>%
  slice(1) %>%
  ungroup()

unmapped_symbols <- ms_proteins_symbols %>% filter(is.na(ENTREZID)) %>% pull(SYMBOL)

ms_proteins_aliases <-
  AnnotationDbi::select(org.Hs.eg.db, keys = unmapped_symbols, keytype = "ALIAS", columns = c("ALIAS", "SYMBOL", "ENTREZID")) %>%
  group_by(SYMBOL) %>%
  slice(1) %>%
  ungroup()


ms_protein_entrez <- bind_rows(
  ms_proteins_symbols %>% mutate(PROTID = SYMBOL),
  ms_proteins_aliases %>% mutate(PROTID = ALIAS)
) %>% filter(!is.na(ENTREZID))


protein_go_terms <-
  AnnotationDbi::select(org.Hs.eg.db, keytypes = "ENTREZID", keys = ms_protein_entrez$ENTREZID,
                        columns = c("ENTREZID", "GO", "ONTOLOGY")) %>%
  mutate(TERM = go_terms[GO]) %>%
  inner_join(ms_protein_entrez, by = "ENTREZID")


protein_chromatin_go_terms <- protein_go_terms %>%
  filter(TERM %in% chromatin_terms)

ms_protein_entrez %>%
  mutate(CHROMATIN_GO = PROTID %in% protein_chromatin_go_terms$PROTID) %>%
  left_join(protein_chromatin_go_terms %>%
              group_by(PROTID) %>% summarise(CHROMATIN_GO_TERMS = paste(TERM, collapse = "|"))) %>%
  write_csv(file.path(out_dir, "prober_chromatin_go_terms.csv"))

saint_fc_dir <- file.path(out_dir, "prober-saint-vs-fc")
dir.create(saint_fc_dir, recursive = T)

fc_axis_limits <-
  c("NFkB(-TNF)" = 50,
    "NFkB (+TNF)" = 50,
    "NFkB(+TNF) RelA-KD" = 50,
    "STAT1(+IFN)" = 30,
    "STAT1(-IFN)" = 30,
    "YY1" = 31,
    "YY1-KD" = 31,
    "rS7296179 C-allele" = 30,
    "rS7296179 G-allele" = 30,
    "hTERT  228MT" = 30,
    "hTERT  228WT" = 30,
    "hTERT  250WT" = 40,
    "hTERT C250MT" = 40,
    "hTERT 161WT" = 42,
    "hTERT A161MT" = 42)

imap(data,
     ~ saint_vs_foldchange_scatterplot(
       .x, protein_labels = protein_chromatin_go_terms$PROTID,
       fc_axis_limit = fc_axis_limits[.y], fc_threshold = 1,
       plot_path = file.path(saint_fc_dir, paste0(str_replace_all(.y, "\\s+", "_"), ".pdf"))))


chromatin_prot_enrichment <-
  tibble(NAME = names(data), data = data) %>%
  unnest(data) %>%
  select(NAME, PROTID, BAIT_SP) %>%
  mutate(chromatin_term = PROTID %in% protein_chromatin_go_terms$PROTID) %>%
  group_by(NAME) %>%
  summarise(A = sum(BAIT_SP >= 0.9 & chromatin_term),
            B = sum(BAIT_SP >= 0.9 & !chromatin_term),
            C = sum(BAIT_SP < 0.9 & chromatin_term),
            D = sum(BAIT_SP < 0.9 & !chromatin_term)) %>%
  mutate(fisher = pmap(., function(A, B, C, D, ...)
    broom::tidy(fisher.test(matrix(c(A, B, C, D), nrow = 2),
                            alternative = "greater")))) %>%
  unnest(fisher)

write_csv(chromatin_prot_enrichment,
          file.path(out_dir, "prober_chromatin_go_term_enrichment_saint0.9.csv"))





contrast_scatter_plot(data = data, contrast = c("YY1", "YY1-KD"),
                      xlab = "YY1 PROBER-MS", ylab = "YY1-KD vs Control (log2FC)",
                      protein_labels = protein_chromatin_go_terms$PROTID, fc_threshold = 1,
                      plot_path = file.path(contrast_dir, paste0("YY1-KD_vs_Ctrl", ".pdf")))

contrast_scatter_plot(data = data, contrast = c("STAT1(-IFN)", "STAT1(+IFN)"), x_data = 2,
                      xlab = "STAT1 +IFN PROBER-MS", ylab = "STAT1 +IFN vs -IFN (log2FC)",
                      protein_labels = protein_chromatin_go_terms$PROTID, fc_threshold = 1,
                      plot_path = file.path(contrast_dir, paste0("STAT1-IFN_vs_Ctrl", ".pdf")))


contrast_scatter_plot(data = data, contrast = c("NFkB(-TNF)", "NFkB (+TNF)"), x_data = 2,
                      xlab = "NFkB +TNF PROBER-MS", ylab = "NFkB +TNF vs -TNF (log2FC)",
                      protein_labels = protein_chromatin_go_terms$PROTID, fc_threshold = 1,
                      plot_path = file.path(contrast_dir, paste0("NFkB-TNF_vs_Ctrl", ".pdf")))

contrast_scatter_plot(data = data, contrast = c("NFkB (+TNF)", "NFkB(+TNF) RelA-KD"), x_data = 1,
                      xlab = "NFkB +TNF PROBER-MS", ylab = "NFkB +TNF RelA-KD vs Control (log2FC)",
                      protein_labels = protein_chromatin_go_terms$PROTID, fc_threshold = 1,
                      plot_path = file.path(contrast_dir, paste0("NFkB-TNF_RelA-KD_vs_Ctrl", ".pdf")))


contrast_scatter_plot(data = data, contrast = c("rS7296179 C-allele", "rS7296179 G-allele"), x_data = 1,
                      xlab = "rS7296179 C-allele PROBER-MS", ylab = "rS7296179 G-allele vs C-allele (log2FC)",
                      protein_labels = protein_chromatin_go_terms$PROTID, fc_threshold = 1,
                      plot_path = file.path(contrast_dir, paste0("rS7296179_G-allele_vs_C-allele", ".pdf")))



contrast_scatter_plot_MA(data = data, contrast = c("YY1", "YY1-KD"),
                         xlab = "Average FC", ylab = "YY1-KD vs Control (log2FC)",
                         protein_labels = protein_chromatin_go_terms$PROTID, fc_threshold = 1,
                         plot_path = file.path(contrast_dir, paste0("YY1-KD_vs_Ctrl", "_MA.pdf")))

contrast_scatter_plot_MA(data = data, contrast = c("STAT1(-IFN)", "STAT1(+IFN)"),
                         xlab = "Average FC", ylab = "STAT1 +IFN vs -IFN (log2FC)",
                         protein_labels = protein_chromatin_go_terms$PROTID, fc_threshold = 1,
                         plot_path = file.path(contrast_dir, paste0("STAT1-IFN_vs_Ctrl", "_MA.pdf")))


contrast_scatter_plot_MA(data = data, contrast = c("NFkB(-TNF)", "NFkB (+TNF)"),
                         xlab = "Average FC", ylab = "NFkB +TNF vs -TNF (log2FC)",
                         protein_labels = protein_chromatin_go_terms$PROTID, fc_threshold = 1,
                         plot_path = file.path(contrast_dir, paste0("NFkB-TNF_vs_Ctrl", "_MA.pdf")))

contrast_scatter_plot_MA(data = data, contrast = c("NFkB (+TNF)", "NFkB(+TNF) RelA-KD"),
                         xlab = "Average FC", ylab = "NFkB +TNF RelA-KD vs Control (log2FC)",
                         protein_labels = protein_chromatin_go_terms$PROTID, fc_threshold = 1,
                         plot_path = file.path(contrast_dir, paste0("NFkB-TNF_RelA-KD_vs_Ctrl", "_MA.pdf")))


contrast_scatter_plot_MA(data = data, contrast = c("rS7296179 C-allele", "rS7296179 G-allele"),
                         xlab = "Average FC", ylab = "rS7296179 G-allele vs C-allele (log2FC)",
                         protein_labels = protein_chromatin_go_terms$PROTID, fc_threshold = 1,
                         plot_path = file.path(contrast_dir, paste0("rS7296179_G-allele_vs_C-allele", "_MA.pdf")))




# remove background points

contrast_scatter_plot_MA(data = data, contrast = c("YY1", "YY1-KD"),
                         xlab = "Average FC", ylab = "YY1-KD vs Control (log2FC)",
                         protein_labels = protein_chromatin_go_terms$PROTID, fc_threshold = 1,
                         show_background = F,
                         plot_path = file.path(contrast_dir, paste0("YY1-KD_vs_Ctrl", "_MA-no-bkg.pdf")))

contrast_scatter_plot_MA(data = data, contrast = c("STAT1(-IFN)", "STAT1(+IFN)"),
                         xlab = "Average FC", ylab = "STAT1 +IFN vs -IFN (log2FC)",
                         protein_labels = protein_chromatin_go_terms$PROTID, fc_threshold = 1,
                         show_background = F,
                         plot_path = file.path(contrast_dir, paste0("STAT1-IFN_vs_Ctrl", "_MA-no-bkg.pdf")))


contrast_scatter_plot_MA(data = data, contrast = c("NFkB(-TNF)", "NFkB (+TNF)"),
                         xlab = "Average FC", ylab = "NFkB +TNF vs -TNF (log2FC)",
                         protein_labels = protein_chromatin_go_terms$PROTID, fc_threshold = 1,
                         show_background = F,
                         plot_path = file.path(contrast_dir, paste0("NFkB-TNF_vs_Ctrl", "_MA-no-bkg.pdf")))

contrast_scatter_plot_MA(data = data, contrast = c("NFkB (+TNF)", "NFkB(+TNF) RelA-KD"),
                         xlab = "Average FC", ylab = "NFkB +TNF RelA-KD vs Control (log2FC)",
                         protein_labels = protein_chromatin_go_terms$PROTID, fc_threshold = 1,
                         show_background = F,
                         plot_path = file.path(contrast_dir, paste0("NFkB-TNF_RelA-KD_vs_Ctrl", "_MA-no-bkg.pdf")))


contrast_scatter_plot_MA(data = data, contrast = c("rS7296179 C-allele", "rS7296179 G-allele"),
                         xlab = "Average FC", ylab = "rS7296179 G-allele vs C-allele (log2FC)",
                         protein_labels = protein_chromatin_go_terms$PROTID, fc_threshold = 1,
                         show_background = F,
                         plot_path = file.path(contrast_dir, paste0("rS7296179_G-allele_vs_C-allele", "_MA-no-bkg.pdf")))







yy1_kd_highlights <-
  c("YY1",
    "YY2",
    "NFRKB",
    "RUVBL1",
    "KANSL1",
    "KANSL3",
    "CHD8",
    "TAF4",
    "TAF9B",
    "ZHX1",
    "ZHX2",
    "ZHX3",
    "SALL3")

nfkb_stim_highlights <-
  c("RELA",
    "REL",
    "BICRA",
    "CHD8",
    "NCOA2",
    "NCOR1",
    "BCOR",
    "KMT2A",
    "SMARCC1",
    "SUPT16H",
    "MYBBP1A",
    "GLO1")

stat_stim_highights <-
  c("STAT1",
    "CREBBP",
    "KMT2D",
    "ERF",
    "NCOR2",
    "ELF1",
    "ELF2",
    "GABPA",
    "FOXP1",
    "FOXP4")

rela_kd_highlights <-
  c("RELA",
    "REL",
    "BICRA",
    "CHD8",
    "NCOA2",
    "NCOR1",
    "BCOR",
    "KMT2A",
    "SMARCC1",
    "SUPT16H",
    "MYBBP1A",
    "GLO1")

snp_highlights <-
  c("RELA",
    "REL",
    "NFKB1",
    "EP300",
    "NCOA2",
    "EP400")



nfkb_tnf_saint_genes <- unique(c(data$`NFkB (+TNF)`  %>% filter(BAIT_SP >= 0.9) %>% pull(PROTID),
                                 data$`NFkB(-TNF)`  %>% filter(BAIT_SP >= 0.9) %>% pull(PROTID)))

contrast_limma_analysis(data, quant_norm = T,
                         contrast_labels = c("NFkB(-TNF)", "NFkB (+TNF)"),
                         bait1 = "S\\d",
                         bait2 = "NFKB",
                         condition1 = "NOTNF",
                         condition2 = "WITHTNF",
                         contrast_name = "baitBAIT2.conditionCOND2",
                         fdr_cutoff = 0.25,
                         # p_cutoff = 0.05,
                         saint_cutoff = 0,
                         pseudo_count = 2,
                         min_protein_count = 10,
                         # protein_labels = protein_chromatin_go_terms$PROTID,
                         protein_highlights = protein_chromatin_go_terms$PROTID,
                         # protein_highlights = nfkb_tnf_saint_genes,
                         logFC_label = "NFkB +TNFa vs vehicle (log2FC)",
                         plot_path = file.path(contrast_dir, "NFkB-TNF_vs_vehicle_full"))

# contrast_limma_analysis2(data, quant_norm = T,
#                          contrast_labels = c("NFkB(-TNF)", "NFkB (+TNF)"),
#                          controls = "S\\d",
#                          condition1 = "NFKB_NOTNF",
#                          condition2 = "NFKB_WITHTNF",
#                          fdr_cutoff = 0.1,
#                          saint_cutoff = 0,
#                          pseudo_count = 2,
#                          min_protein_count = 10,
#                          protein_labels = protein_chromatin_go_terms$PROTID,
#                          protein_highlights = nfkb_stim_highlights,
#                          logFC_label = "NFkB +TNFa vs vehicle (log2FC)",
#                          plot_path = file.path(contrast_dir, "NFkB-TNF_vs_vehicle_groups"))
#
#
# contrast_limma_analysis2(data, quant_norm = T,
#                          contrast_labels = c("NFkB(-TNF)", "NFkB (+TNF)"),
#                          condition1 = "NFKB_NOTNF",
#                          condition2 = "NFKB_WITHTNF",
#                          fdr_cutoff = 0.1,
#                          saint_cutoff = 0,
#                          pseudo_count = 2,
#                          min_protein_count = 10,
#                          protein_labels = protein_chromatin_go_terms$PROTID,
#                          protein_highlights = nfkb_stim_highlights,
#                          logFC_label = "NFkB +TNFa vs vehicle (log2FC)",
#                          plot_path = file.path(contrast_dir, "NFkB-TNF_vs_vehicle_noctrls"))


stat_ifn_saint_genes <- unique(c(data$`STAT1(+IFN)`  %>% filter(BAIT_SP >= 0.9) %>% pull(PROTID),
                                 data$`STAT1(-IFN)`  %>% filter(BAIT_SP >= 0.9) %>% pull(PROTID)))

contrast_limma_analysis(data, quant_norm = T,
                         contrast_labels = c("STAT1(-IFN)", "STAT1(+IFN)"),
                         bait1 = "S\\d",
                         bait2 = "STAT1",
                         condition1 = "NOIFN",
                         condition2 = "WITHIFN",
                         contrast_name = "baitBAIT2.conditionCOND2",
                         # protein_labels = protein_chromatin_go_terms$PROTID,
                         protein_highlights = protein_chromatin_go_terms$PROTID,
                         # protein_highlights = stat_ifn_saint_genes,

                         min_protein_count = 10,
                         pseudo_count = 2,
                         fdr_cutoff = 0.25,
                         # p_cutoff = 0.05,
                         saint_cutoff = 0,
                         logFC_label = "STAT1 +IFNg vs vehicle (log2FC)",
                         plot_path = file.path(contrast_dir, "STAT1-IFN_vs_vehicle_full"))



yy1_kd_saint_genes <- unique(c(data$YY1  %>% filter(BAIT_SP >= 0.9) %>% pull(PROTID),
                               data$`YY1-KD`  %>% filter(BAIT_SP >= 0.9) %>% pull(PROTID)))

contrast_limma_analysis(data, quant_norm = T,
                         contrast_labels = c("YY1", "YY1-KD"),
                         controls = "SCR\\d",
                         condition1 = "YY1",
                         condition2 = "YY1-KD",
                         min_protein_count = 5,
                         pseudo_count = 2,
                         fdr_cutoff = 0.25,
                         # p_cutoff = 0.05,
                         saint_cutoff = 0,
                         # protein_labels = protein_chromatin_go_terms$PROTID,
                         protein_highlights = protein_chromatin_go_terms$PROTID,
                         # protein_highlights = yy1_kd_saint_genes,
                         logFC_label = "YY1-KD vs mock (log2FC)",
                         plot_path = file.path(contrast_dir, "YY1-KD_vs_mock_groups"))



rela_kd_saint_genes <- unique(c(data$`NFkB (+TNF)`  %>% filter(BAIT_SP >= 0.9) %>% pull(PROTID),
                                data$`NFkB(+TNF) RelA-KD`  %>% filter(BAIT_SP >= 0.9) %>% pull(PROTID)))

contrast_limma_analysis(data, quant_norm = T,
                         contrast_labels = c("NFkB (+TNF)", "NFkB(+TNF) RelA-KD"),
                         controls = "S\\d_(WITH|NO)TNF",
                         condition1 = "NFKB_WITHTNF",
                         condition2 = "SHRNA",
                         min_protein_count = 5,
                         pseudo_count = 2,
                         fdr_cutoff = 0.25,
                         # p_cutoff = 0.05,
                         saint_cutoff = 0,
                         # protein_labels = protein_chromatin_go_terms$PROTID,
                         protein_highlights = protein_chromatin_go_terms$PROTID,
                         # protein_highlights = rela_kd_saint_genes,
                         logFC_label = "RELA-KD vs mock (log2FC)",
                         plot_path = file.path(contrast_dir, "NFkB-TNF_RelA-KD_vs_mock_groups"))

rs7296179_saint_genes <- unique(c(data$`rS7296179 G-allele`  %>% filter(BAIT_SP >= 0.9) %>% pull(PROTID),
                                  data$`rS7296179 C-allele`  %>% filter(BAIT_SP >= 0.9) %>% pull(PROTID)))

contrast_limma_analysis(data = data, quant_norm = T,
                         contrast_labels = c("rS7296179 G-allele", "rS7296179 C-allele"),
                         controls = "S\\d",
                         condition1 = "RS7296179-LOW",
                         condition2 = "RS7296179_HIGH",
                         logFC_label = "rS7296179 C-allele vs G-allele (log2FC)",
                         min_protein_count = 10,
                         pseudo_count = 2,
                         fdr_cutoff = 0.25,
                         # p_cutoff = 0.05,
                         saint_cutoff = 0,
                         # protein_labels = protein_chromatin_go_terms$PROTID,
                         protein_highlights = protein_chromatin_go_terms$PROTID,
                         # protein_highlights = rs7296179_saint_genes,
                         plot_path = file.path(contrast_dir, "rS7296179_C-allele_vs_G-allele_groups"))



tert161_saint_genes <- unique(c(data$`hTERT 161WT`  %>% filter(BAIT_SP >= 0.9) %>% pull(PROTID),
                                data$`hTERT A161MT`  %>% filter(BAIT_SP >= 0.9) %>% pull(PROTID)))

contrast_limma_analysis(data = data, quant_norm = T,
                         contrast_labels = c("hTERT 161WT", "hTERT A161MT"),
                         controls = "SCR\\d",
                         condition1 = "TRIMER-WT-A161C",
                         condition2 = "TRIMER-MUT-A161C",
                         logFC_label = "hTERT A161C-MUT vs WT (log2FC)",
                         min_protein_count = 5,
                         pseudo_count = 2,
                         fdr_cutoff = 0.25,
                         # p_cutoff = 0.05,
                         p_limit = 2e-4,
                         saint_cutoff = 0,
                         # protein_labels = protein_chromatin_go_terms$PROTID,
                         protein_highlights = protein_chromatin_go_terms$PROTID,
                         # protein_highlights = c("ELF1", "ELF2", "ETV3", "ELK1", "ERF",
                         #                        tert161_saint_genes),
                         plot_path = file.path(contrast_dir, "hTERT_A161C-MUT_vs_WT_groups"))


tert228_saint_genes <- unique(c(data$`hTERT  228WT`  %>% filter(BAIT_SP >= 0.9) %>% pull(PROTID),
                                data$`hTERT 228MT`  %>% filter(BAIT_SP >= 0.9) %>% pull(PROTID)))

contrast_limma_analysis(data = data, quant_norm = T,
                         contrast_labels = c("hTERT  228WT", "hTERT 228MT"),
                         controls = "SCR\\d",
                         condition1 = "TRIMER-WT-C228",
                         condition2 = "TRIMER-MUT-C228T",
                         logFC_label = "hTERT C228T-MUT vs WT (log2FC)",
                         min_protein_count = 5,
                         pseudo_count = 2,
                         fdr_cutoff = 0.25,
                         # p_cutoff = 0.05,
                         p_limit = 3e-4,
                         saint_cutoff = 0,
                         # protein_labels = protein_chromatin_go_terms$PROTID,
                         protein_highlights = protein_chromatin_go_terms$PROTID,
                         # protein_highlights = c("ELF1", "ELF2", "ETV3", "ELK1", "ERF",
                         #                        tert228_saint_genes),
                         plot_path = file.path(contrast_dir, "hTERT_C228T-MUT_vs_WT_groups"))


tert250_saint_genes <- unique(c(data$`hTERT  250WT`  %>% filter(BAIT_SP >= 0.9) %>% pull(PROTID),
                                data$`hTERT C250MT`  %>% filter(BAIT_SP >= 0.9) %>% pull(PROTID)))

contrast_limma_analysis(data = data, quant_norm = T,
                         contrast_labels = c("hTERT  250WT", "hTERT C250MT"),
                         controls = "SCR\\d",
                         condition1 = "TRIMER-WT-C250T",
                         condition2 = "TRIMER-MUT-C250T",
                         logFC_label = "hTERT C250T-MUT vs WT (log2FC)",
                         min_protein_count = 5,
                         pseudo_count = 2,
                         fdr_cutoff = 0.25,
                         # p_cutoff = 0.01,
                         p_limit = 3e-5,
                         saint_cutoff = 0,
                         # protein_labels = protein_chromatin_go_terms$PROTID,
                         protein_highlights = protein_chromatin_go_terms$PROTID,
                         # protein_highlights = c("ELF1", "ELF2", "ETV3", "ELK1", "ERF",
                         #                        tert250_saint_genes),
                         plot_path = file.path(contrast_dir, "hTERT_C250T-MUT_vs_WT_groups"))

#####
# additional MS data

ms_data <- list.files(here("data/MS"), pattern = "^[A-Za-z0-9].+\\.xlsx", full.names = T) %>%
  set_names(., str_replace(basename(.), "\\.xlsx$", "")) %>% map(read_excel)

ms_protein_list <- map(ms_data, ~ .$PROTID) %>% unlist() %>% unique() %>%
  stringi::stri_replace_all(regex = "(C.+)ORF(\\d+)", replacement = "$1orf$2") %>%
  c(ms_protein_list)

ms_proteins_symbols <-
  AnnotationDbi::select(org.Hs.eg.db, keys = ms_protein_list, keytype = "SYMBOL", columns = c("SYMBOL", "ENTREZID")) %>%
  group_by(SYMBOL) %>%
  slice(1) %>%
  ungroup()

unmapped_symbols <- ms_proteins_symbols %>% filter(is.na(ENTREZID)) %>% pull(SYMBOL)

ms_proteins_aliases <-
  AnnotationDbi::select(org.Hs.eg.db, keys = unmapped_symbols, keytype = "ALIAS", columns = c("ALIAS", "SYMBOL", "ENTREZID")) %>%
  group_by(SYMBOL) %>%
  slice(1) %>%
  ungroup()


ms_protein_entrez <- bind_rows(
  ms_proteins_symbols %>% mutate(PROTID = SYMBOL),
  ms_proteins_aliases %>% mutate(PROTID = ALIAS)
) %>% filter(!is.na(ENTREZID))


protein_go_terms <-
  AnnotationDbi::select(org.Hs.eg.db, keytypes = "ENTREZID", keys = ms_protein_entrez$ENTREZID,
                        columns = c("ENTREZID", "GO", "ONTOLOGY")) %>%
  mutate(TERM = go_terms[GO]) %>%
  inner_join(ms_protein_entrez, by = "ENTREZID")


protein_chromatin_go_terms <- protein_go_terms %>%
  filter(TERM %in% chromatin_terms)

ms_protein_entrez %>%
  mutate(CHROMATIN_GO = PROTID %in% protein_chromatin_go_terms$PROTID) %>%
  left_join(protein_chromatin_go_terms %>%
              group_by(PROTID) %>% summarise(CHROMATIN_GO_TERMS = paste(TERM, collapse = "|"))) %>%
  write_csv(file.path(out_dir, "prober_chromatin_go_terms.csv"))


## YY1 U2OS vs 293T

yy1_celltype_genes <- unique(c(ms_data$`U2OS-YY1`  %>% filter(BAIT_SP >= 0.9) %>% pull(PROTID),
                                  data$YY1  %>% filter(BAIT_SP >= 0.9) %>% pull(PROTID)))

list(YY1_293T = data$YY1 %>% rename_with(~ str_replace(., "-R?(\\d)$", "_293T-R\\1"), cols = starts_with("IP_")),
     YY1_U2OS = ms_data$`U2OS-YY1` %>% rename_with(~ str_replace(., "-R?(\\d)$", "_U2OS-R\\1"), cols = starts_with("IP_"))) %>%
  contrast_limma_analysis(quant_norm = T,
                           contrast_labels = c("YY1_U2OS", "YY1_293T"),
                           bait1 = "SCR\\d?",
                           bait2 = "YY1",
                           condition1 = "293T",
                           condition2 = "U2OS",
                           min_protein_count = 5,
                           pseudo_count = 2,
                           fdr_cutoff = 0.25,
                           fdr_line = 0.25,
                           # p_cutoff = 0.05,
                           saint_cutoff = 0,
                           protein_highlights = protein_chromatin_go_terms$PROTID,
                           label_size = 1.75,
                           point_sizes = c(0.75, 1.5),
                           logFC_label = "YY1 U2OS vs 293T (log2FC)",
                           plot_path = file.path(contrast_dir, "YY1-U2OS_vs_293T_full"))


bs1_celltype_genes <- unique(c(ms_data$`U2OS-26-9`  %>% filter(BAIT_SP >= 0.9) %>% pull(PROTID),
                               ms_data$`26-9_SAINT`  %>% filter(BAIT_SP >= 0.9) %>% pull(PROTID)))

list(BS1_293T = ms_data$`26-9_SAINT` %>% rename_with(~ str_replace(., "-R?(\\d)$", "_293T-R\\1"), cols = starts_with("IP_")),
     BS1_U2OS = ms_data$`U2OS-26-9` %>% rename_with(~ str_replace(., "-R?(\\d)$", "_U2OS-R\\1"), cols = starts_with("IP_"))) %>%
  contrast_limma_analysis(quant_norm = T,
                           contrast_labels = c("BS1_U2OS", "BS1_293T"),
                           bait1 = "SCR\\d?",
                           bait2 = "26-9",
                           condition1 = "293T",
                           condition2 = "U2OS",
                           min_protein_count = 5,
                           pseudo_count = 2,
                           fdr_cutoff = 0.25,
                           fdr_line = 0.25,
                           # p_cutoff = 0.05,
                           saint_cutoff = 0,
                           protein_highlights = protein_chromatin_go_terms$PROTID,
                           label_size = 1.75,
                           point_sizes = c(0.75, 1.5),
                           logFC_label = "BS1 U2OS vs 293T (log2FC)",
                           plot_path = file.path(contrast_dir, "BS1-U2OS_vs_293T_full"))


### DNA-Bead Pull Downs

yy1_dna_pd_genes <- unique(c(ms_data$`YY1PD_vs_beads+SCR`%>% filter(BAIT_SP >= 0.9) %>% pull(PROTID),
                             data$YY1  %>% filter(BAIT_SP >= 0.9) %>% pull(PROTID)))

list(YY1_prober = data$YY1 %>% rename_with(~ str_replace(., "-R?(\\d)$", "_prober-R\\1"), cols = starts_with("IP_")),
     YY1_dna_pd = ms_data$`YY1PD_vs_beads+SCR` %>% rename_with(~ str_replace(., "-R?(\\d)$", "_PD-R\\1"), cols = starts_with("IP_"))) %>%
  contrast_limma_analysis(quant_norm = T,
                           contrast_labels = c("YY1_dna_pd", "YY1_prober"),
                           bait1 = "SCR\\d?",
                           bait2 = "YY1",
                           condition1 = "PD",
                           condition2 = "prober",
                           min_protein_count = 5,
                           pseudo_count = 2,
                           fdr_cutoff = 0.25,
                           fdr_line = 0.25,
                           # p_cutoff = 0.05,
                           saint_cutoff = 0,
                           protein_highlights = protein_chromatin_go_terms$PROTID,
                           label_size = 1.75,
                           point_sizes = c(0.75, 1.5),
                           logFC_label = "YY1 PROBER vs DNA PD (log2FC)",
                           plot_path = file.path(contrast_dir, "YY1-PROBER_vs_DNAPD_full"))


nfkbnotnf_dna_pd_genes <- unique(c(ms_data$`NFkBPD–TNF` %>% filter(BAIT_SP >= 0.9) %>% pull(PROTID),
                             data$`NFkB(-TNF)`  %>% filter(BAIT_SP >= 0.9) %>% pull(PROTID)))

list(NFKB_prober = data$`NFkB(-TNF)` %>% rename_with(~ str_replace(., "[_-]R?(\\d)$", "_prober-R\\1"), cols = starts_with("IP_")),
     NFKB_dna_pd = ms_data$`NFkBPD–TNF` %>%
       rename_with(~ str_replace(., "RELA", "NFKB_NOTNF_R"), cols = starts_with("IP_")) %>%
       rename_with(~ str_replace(., "SCR", "S1_NOTNF_R"), cols = starts_with("IP_")) %>%
                     rename_with(~ str_replace(., "[_-]R?(\\d)$", "_PD-R\\1"), cols = starts_with("IP_"))) %>%
  contrast_limma_analysis(quant_norm = T,
                           contrast_labels = c("NFKB_dna_pd", "NFKB_prober"),
                           bait1 = "S\\d?_NOTNF",
                           bait2 = "NFKB_NOTNF",
                           condition1 = "PD",
                           condition2 = "prober",
                           min_protein_count = 5,
                           pseudo_count = 2,
                           fdr_cutoff = 0.25,
                           fdr_line = 0.25,
                           # p_cutoff = 0.05,
                           saint_cutoff = 0,
                           protein_highlights = protein_chromatin_go_terms$PROTID,
                           label_size = 1.75,
                           point_sizes = c(0.75, 1.5),
                           logFC_label = "NFKB-TNF PROBER vs DNA PD (log2FC)",
                           plot_path = file.path(contrast_dir, "NFKB-TNF-PROBER_vs_DNAPD_full"))


nfkbtnf_dna_pd_genes <- unique(c(ms_data$`NFkBPD+TNF` %>% filter(BAIT_SP >= 0.9) %>% pull(PROTID),
                              data$`NFkB (+TNF)`  %>% filter(BAIT_SP >= 0.9) %>% pull(PROTID)))

list(NFKB_prober = data$`NFkB (+TNF)` %>% rename_with(~ str_replace(., "[_-]R?(\\d)$", "_prober-R\\1"), cols = starts_with("IP_")),
     NFKB_dna_pd = ms_data$`NFkBPD+TNF` %>%
       rename_with(~ str_replace(., "RELA\\+TNF", "NFKB_WITHTNF"), cols = starts_with("IP_")) %>%
       rename_with(~ str_replace(., "SCR\\+TNF", "S1_WITHTNF"), cols = starts_with("IP_")) %>%
       rename_with(~ str_replace(., "[_-]R?(\\d)$", "_PD-R\\1"), cols = starts_with("IP_"))) %>%
  contrast_limma_analysis(quant_norm = T,
                           contrast_labels = c("NFKB_dna_pd", "NFKB_prober"),
                           bait1 = "S\\d?_WITHTNF",
                           bait2 = "NFKB_WITHTNF",
                           condition1 = "PD",
                           condition2 = "prober",
                           min_protein_count = 5,
                           pseudo_count = 2,
                           fdr_cutoff = 0.25,
                           fdr_line = 0.25,
                           # p_cutoff = 0.05,
                           saint_cutoff = 0,
                           protein_highlights = protein_chromatin_go_terms$PROTID,
                           label_size = 1.75,
                           point_sizes = c(0.75, 1.5),
                           logFC_label = "NFKB+TNF PROBER vs DNA PD (log2FC)",
                           plot_path = file.path(contrast_dir, "NFKB+TNF-PROBER_vs_DNAPD_full"))



nfkb_dna_pd_genes <- unique(c(ms_data$`NFkBPD+TNF` %>% filter(BAIT_SP >= 0.9) %>% pull(PROTID),
                              ms_data$`NFkBPD–TNF` %>% filter(BAIT_SP >= 0.9) %>% pull(PROTID)))

list(NFKB_notnf = ms_data$`NFkBPD–TNF` %>%
       rename_with(~ str_replace(., "RELA", "NFKB_NOTNF_R"), cols = starts_with("IP_")) %>%
       rename_with(~ str_replace(., "SCR", "S1_NOTNF_R"), cols = starts_with("IP_")),
     NFKB_tnf = ms_data$`NFkBPD+TNF` %>%
       rename_with(~ str_replace(., "RELA\\+TNF", "NFKB_WITHTNF"), cols = starts_with("IP_")) %>%
       rename_with(~ str_replace(., "SCR\\+TNF", "S1_WITHTNF"), cols = starts_with("IP_"))) %>%
  contrast_limma_analysis(quant_norm = T,
                           contrast_labels = c("NFKB_notnf", "NFKB_tnf"),
                           bait1 = "S\\d",
                           bait2 = "NFKB",
                           condition1 = "NOTNF",
                           condition2 = "WITHTNF",
                           min_protein_count = 5,
                           pseudo_count = 2,
                           fdr_cutoff = 0.25,
                           fdr_line = 0.25,
                           # p_cutoff = 0.05,
                           saint_cutoff = 0,
                           protein_highlights = protein_chromatin_go_terms$PROTID,
                           label_size = 1.75,
                           point_sizes = c(0.75, 1.5),
                           logFC_label = "NFKB DNA-PD +TNF vs -TNF (log2FC)",
                           plot_path = file.path(contrast_dir, "NFKB-DNAPD_+TNF_vs_-TNF_full"))



### BioID

yy1_bioid_genes <- unique(c(ms_data$`bioid-YY1` %>% filter(BAIT_SP >= 0.9) %>% pull(PROTID),
                             data$YY1  %>% filter(BAIT_SP >= 0.9) %>% pull(PROTID)))

list(YY1_prober = data$YY1 %>%
       rename_with(~ str_replace(., "SCR", "CTRL"), cols = starts_with("IP_")) %>%
       rename_with(~ str_replace(., "[-_]R?(\\d)$", "_prober-R\\1"), cols = starts_with("IP_")),
     YY1_bioid = ms_data$`bioid-YY1` %>%
       rename_with(~ str_replace(., "NLS_BASU_GFP", "CTRL1"), cols = starts_with("IP_")) %>%
       # rename_with(~ str_replace(., "EGFP_BASU", "CTRL2"), cols = starts_with("IP_")) %>%
       rename_with(~ str_replace(., "[-_]R?(\\d)$", "_bioid-R\\1"), cols = starts_with("IP_"))) %>%
  contrast_limma_analysis(quant_norm = T,
                           contrast_labels = c("YY1_prober", "YY1_bioid"),
                           bait1 = "CTRL\\d?",
                           bait2 = "YY1",
                           condition1 = "bioid",
                           condition2 = "prober",
                           min_protein_count = 5,
                           pseudo_count = 2,
                           fdr_cutoff = 0.25,
                           fdr_line = 0.25,
                           # p_cutoff = 0.05,
                           saint_cutoff = 0,
                           protein_highlights = protein_chromatin_go_terms$PROTID,
                           label_size = 1.75,
                           point_sizes = c(0.75, 1.5),
                           logFC_label = "YY1 PROBER vs BioID (log2FC)",
                           plot_path = file.path(contrast_dir, "YY1-PROBER_vs_BioID_full"))


nfkbnotnf_bioid_genes <- unique(c(ms_data$`bioid-RelA-TNF` %>% filter(BAIT_SP >= 0.9) %>% pull(PROTID),
                                   data$`NFkB(-TNF)`  %>% filter(BAIT_SP >= 0.9) %>% pull(PROTID)))

list(NFKB_prober = data$`NFkB(-TNF)` %>%
       rename_with(~ str_replace(., "S(\\d)", "CTRL\\1"), cols = starts_with("IP_")) %>%
       rename_with(~ str_replace(., "[_-]R?(\\d)$", "_prober-R\\1"), cols = starts_with("IP_")),
     NFKB_bioid = ms_data$`bioid-RelA-TNF` %>%
       rename_with(~ str_replace(., "RELA", "NFKB_NOTNF"), cols = starts_with("IP_")) %>%
       # rename_with(~ str_replace(., "NLS_BASU_GFP", "CTRL1_NOTNF"), cols = starts_with("IP_")) %>%
       rename_with(~ str_replace(., "EGFP_BASU", "CTRL2_NOTNF"), cols = starts_with("IP_")) %>%
       rename_with(~ str_replace(., "[-_]R?(\\d)$", "_bioid-R\\1"), cols = starts_with("IP_"))) %>%
  contrast_limma_analysis(quant_norm = T,
                           contrast_labels = c("NFKB_bioid", "NFKB_prober"),
                           bait1 = "CTRL\\d?_NOTNF",
                           bait2 = "NFKB_NOTNF",
                           condition1 = "bioid",
                           condition2 = "prober",
                           min_protein_count = 5,
                           pseudo_count = 2,
                           fdr_cutoff = 0.25,
                           fdr_line = 0.25,
                           # p_cutoff = 0.05,
                           saint_cutoff = 0,
                           protein_highlights = protein_chromatin_go_terms$PROTID,
                           label_size = 1.75,
                           point_sizes = c(0.75, 1.5),
                           logFC_label = "NFKB-TNF PROBER vs BioID (log2FC)",
                           plot_path = file.path(contrast_dir, "NFKB-TNF-PROBER_vs_BioID_full"))




nfkbtnf_bioid_genes <- unique(c(ms_data$`bioid-RelA+TNF` %>% filter(BAIT_SP >= 0.9) %>% pull(PROTID),
                                data$`NFkB (+TNF)`  %>% filter(BAIT_SP >= 0.9) %>% pull(PROTID)))

list(NFKB_prober = data$`NFkB (+TNF)` %>%
       rename_with(~ str_replace(., "S(\\d)", "CTRL\\1"), cols = starts_with("IP_")) %>%
       rename_with(~ str_replace(., "[_-]R?(\\d)$", "_prober-R\\1"), cols = starts_with("IP_")),
     NFKB_bioid = ms_data$`bioid-RelA+TNF` %>%
       rename_with(~ str_replace(., "RELA_TNF", "NFKB_WITHTNF"), cols = starts_with("IP_")) %>%
       rename_with(~ str_replace(., "NLS_BASU_GFP(_TNF)?", "CTRL1_WITHTNF"), cols = starts_with("IP_")) %>%
       # rename_with(~ str_replace(., "EGFP_BASU", "CTRL2_WITHTNF"), cols = starts_with("IP_")) %>%
       rename_with(~ str_replace(., "[-_]R?(\\d)$", "_bioid-R\\1"), cols = starts_with("IP_"))) %>%
  contrast_limma_analysis(quant_norm = T,
                           contrast_labels = c("NFKB_bioid", "NFKB_prober"),
                           bait1 = "CTRL\\d?_WITHTNF",
                           bait2 = "NFKB_WITHTNF",
                           condition1 = "bioid",
                           condition2 = "prober",
                           min_protein_count = 5,
                           pseudo_count = 2,
                           fdr_cutoff = 0.25,
                           fdr_line = 0.25,
                           # p_cutoff = 0.05,
                           saint_cutoff = 0,
                           protein_highlights = protein_chromatin_go_terms$PROTID,
                           label_size = 1.75,
                           point_sizes = c(0.75, 1.5),
                           logFC_label = "NFKB+TNF PROBER vs BioID (log2FC)",
                           plot_path = file.path(contrast_dir, "NFKB+TNF-PROBER_vs_BioID_full"))




nfkb_bioid_genes <- unique(c(ms_data$`NFkBPD+TNF` %>% filter(BAIT_SP >= 0.9) %>% pull(PROTID),
                              ms_data$`NFkBPD–TNF` %>% filter(BAIT_SP >= 0.9) %>% pull(PROTID)))

list(NFKB_notnf = ms_data$`bioid-RelA-TNF` %>%
       rename_with(~ str_replace(., "RELA", "NFKB_NOTNF"), cols = starts_with("IP_")) %>%
       rename_with(~ str_replace(., "NLS_BASU_GFP", "CTRL1_NOTNF"), cols = starts_with("IP_")) %>%
       rename_with(~ str_replace(., "EGFP_BASU", "CTRL2_NOTNF"), cols = starts_with("IP_")),
     NFKB_tnf = ms_data$`bioid-RelA+TNF` %>%
       rename_with(~ str_replace(., "RELA_TNF", "NFKB_WITHTNF"), cols = starts_with("IP_")) %>%
       rename_with(~ str_replace(., "NLS_BASU_GFP(_TNF)?", "CTRL1_WITHTNF"), cols = starts_with("IP_")) %>%
       rename_with(~ str_replace(., "EGFP_BASU(_TNF)?", "CTRL2_WITHTNF"), cols = starts_with("IP_"))) %>%
  contrast_limma_analysis(quant_norm = T,
                           contrast_labels = c("NFKB_notnf", "NFKB_tnf"),
                           bait1 = "CTRL\\d",
                           bait2 = "NFKB",
                           condition1 = "NOTNF",
                           condition2 = "WITHTNF",
                           min_protein_count = 5,
                           pseudo_count = 2,
                           fdr_cutoff = 0.25,
                           fdr_line = 0.25,
                           # p_cutoff = 0.05,
                           saint_cutoff = 0,
                           protein_highlights = protein_chromatin_go_terms$PROTID,
                           protein_labels = c("NFKBIB", "NFKBIA", "EP300", "RELB",
                                              "CREBBP", "NFKBIE", "NFKB2", "NFKB1", "RELA"),
                           label_size = 1.75,
                           point_sizes = c(0.75, 1.5),
                           logFC_label = "NFKB BioID +TNF vs -TNF (log2FC)",
                           plot_path = file.path(contrast_dir, "NFKB-BioID_+TNF_vs_-TNF_full"))



### scatters
scatter_dir <- file.path(out_dir, "saint_scatter")
dir.create(scatter_dir, recursive = T)

saint_scatter <- function(saint_x, saint_y, plot_dir, plot_name = "SAINT_scatter",
                          label_x = "SAINT FC x", label_y = "SAINT FC y",
                          label_n = 15,
                          label_both = Inf,
                          additional_labels = c()) {

  df <- full_join(saint_x %>% select(PROTID, BAIT_FC_A, BAIT_SP),
                  saint_y %>% select(PROTID, BAIT_FC_A, BAIT_SP),
                  by = "PROTID") %>%
    mutate(across(starts_with("BAIT"), ~ ifelse(is.na(.), 0, .))) %>%
    mutate(saint = case_when(BAIT_SP.x >= 0.9 & BAIT_SP.y >= 0.9 ~ "both",
                             BAIT_SP.x >= 0.9 ~ "x",
                             BAIT_SP.y >= 0.9 ~ "y",
                             TRUE ~ NA_character_)) %>%
    group_by(saint) %>%
    mutate(rank.x = min_rank(-BAIT_FC_A.x), rank.y = min_rank(-BAIT_FC_A.y),
           rank.both = min_rank(-sqrt(BAIT_FC_A.x*BAIT_FC_A.y))) %>%
    ungroup() %>%
    mutate(label = case_when( PROTID %in% additional_labels ~ PROTID,
                              saint == "both" & rank.both <= label_both ~ PROTID,
                              saint == "x" & rank.x <= label_n ~ PROTID,
                              saint == "y" & rank.y <= label_n ~ PROTID,
                              TRUE ~ "")) %>%
    arrange(label, !is.na(saint))

  ggplot(df, aes(x = BAIT_FC_A.x, y = BAIT_FC_A.y)) +
    # geom_point(color = "grey50", size = 0.75,
    #            data = df %>% filter(is.na(saint))) +
    geom_point(aes(color = saint, size = ifelse(is.na(saint), 0.5, 1.25)), show.legend = F) +
    geom_text_repel(aes(label = label), point.size = 0.5,
                    data = df %>% filter(label != "" | !is.na(saint)),
                    # nudge_x = 0.1, nudge_y = 0.1,
                    # point.padding = unit(0.05, "lines"),
                    # box.padding = unit(0.2, "lines"),
                    force = 2, force_pull = 2,
                    max.time = 5,
                    size = 2.5, fontface = "bold",
                    show.legend = F, max.overlaps = Inf,
                    min.segment.length = 0, segment.size = 0.2) +
    scale_size_identity() +
    scale_size_identity(aesthetics = "point.size") +
    scale_x_continuous(trans = pseudo_log_trans(base = 2, sigma = 4),
                       expand = expansion(mult = c(0.1, 0.05)),
                       breaks = function(x) unique(c(0, 1, 2, 5, 10, extended_breaks(n = 5)(x)))) +
    scale_y_continuous(trans = pseudo_log_trans(base = 2, sigma = 4),
                       expand = expansion(mult = c(0.1, 0.05)),
                       breaks = function(x) unique(c(0, 1, 2, 5, 10, extended_breaks(n = 5)(x)))) +
    expand_limits(x = max(df$BAIT_FC_A.x, df$BAIT_FC_A.y), y = max(df$BAIT_FC_A.x, df$BAIT_FC_A.y)) +
    labs(x = label_x, y = label_y) +
    coord_equal() +
    theme_cowplot()
  ggsave(file.path(plot_dir, paste0(plot_name, ".pdf")), width = 5, height = 5)
}

saint_scatter(ms_data$`bioid-YY1`, data$YY1, scatter_dir,
              label_x = "BioID YY1", label_y = "PROBER YY1",
              plot_name = "YY1_PROBER_vs_BioID")



saint_scatter(ms_data$`YY1PD_vs_beads+SCR`, data$YY1, scatter_dir,
              label_x = "DNA PD YY1 ", label_y = "PROBER YY1",
              plot_name = "YY1_PROBER_vs_DNA-PD")



saint_scatter(ms_data$`26-11-SAINT`, ms_data$`26-9_SAINT`, scatter_dir,
              label_x = "YY1 BS2", label_y = "YY1 BS1",
              plot_name = "PROBER_YY1_BS1_vs_BS2")
saint_scatter(ms_data$`26-14-SAINT`, ms_data$`26-9_SAINT`, scatter_dir,
              label_x = "YY1 BS3", label_y = "YY1 BS1",
              plot_name = "PROBER_YY1_BS1_vs_BS3")
saint_scatter(ms_data$`26-14-SAINT`, ms_data$`26-11-SAINT`, scatter_dir,
              label_x = "YY1 BS3", label_y = "YY1 BS2",
              plot_name = "PROBER_YY1_BS2_vs_BS3")


saint_scatter(ms_data$`26-9_SAINT`, ms_data$`U2OS-26-9`, scatter_dir,
              label_x = "293T YY1 BS1", label_y = "U2OS YY1 BS1",
              plot_name = "YY1_BS1_U2OS_vs_293T")

saint_scatter(data$YY1, ms_data$`U2OS-YY1`, scatter_dir,
              label_x = "293T YY1", label_y = "U2OS YY1",
              plot_name = "YY1_U2OS_vs_293T")


saint_scatter(ms_data$`bioid-RelA+TNF`, data$`NFkB (+TNF)`, scatter_dir,
              label_x = "BioID RelA +TNF", label_y = "PROBER NFkB +TNF",
              plot_name = "NFKB+TNF_PROBER_vs_BioID",
              additional_labels = c("NFKBIB", "NFKBIA", "EP300", "RELB",
                                    "CREBBP", "NFKBIE", "NFKB2", "NFKB1", "RELA"))

saint_scatter(ms_data$`bioid-RelA-TNF`, data$`NFkB(-TNF)`, scatter_dir,
              label_x = "BioID RelA -TNF", label_y = "PROBER NFkB -TNF",
              plot_name = "NFKB-TNF_PROBER_vs_BioID",
              additional_labels = c("NFKBIB", "NFKBIA", "EP300", "RELB",
                                    "CREBBP", "NFKBIE", "NFKB2", "NFKB1", "RELA"))

saint_scatter(ms_data$`NFkBPD+TNF`, data$`NFkB (+TNF)`, scatter_dir,
              label_x = "DNA PD NFkB+TNF", label_y = "PROBER NFkB+TNF",
              plot_name = "NFKB+TNF_PROBER_vs_DNAPD",
              additional_labels = c("NFKBIB", "NFKBIA", "EP300", "RELB",
                                    "CREBBP", "NFKBIE", "NFKB2", "NFKB1", "RELA"))

saint_scatter(ms_data$`NFkBPD–TNF`, data$`NFkB(-TNF)`, scatter_dir,
              label_x = "DNA PD NFkB-TNF", label_y = "PROBER NFkB-TNF",
              plot_name = "NFKB-TNF_PROBER_vs_DNAPD",
              additional_labels = c("NFKBIB", "NFKBIA", "EP300", "RELB",
                                    "CREBBP", "NFKBIE", "NFKB2", "NFKB1", "RELA"))



saint_scatter(ms_data$`NFkBPD–TNF`, ms_data$`NFkBPD+TNF`, scatter_dir,
              label_x = "DNA-PD NFkB -TNF", label_y = "DNA-PD NFkB +TNF",
              plot_name = "NFKB_DNAPD_+TNF_vs_-TNF",
              additional_labels = c("NFKBIB", "NFKBIA", "EP300", "RELB",
                                    "CREBBP", "NFKBIE", "NFKB2", "NFKB1", "RELA"))


saint_scatter(ms_data$`bioid-RelA-TNF`, ms_data$`bioid-RelA+TNF`, scatter_dir,
              label_both = 25,
              label_x = "BioID RelA -TNF", label_y = "BioID RelA +TNF",
              plot_name = "RelA_BioID_+TNF_vs_-TNF",
              additional_labels = c("NFKBIB", "NFKBIA", "EP300", "RELB",
                                    "CREBBP", "NFKBIE", "NFKB2", "NFKB1", "RELA"))


saint_scatter(data$`NFkB(-TNF)`, data$`NFkB (+TNF)`, scatter_dir,
              label_x = "PROBER NFkB -TNF", label_y = "PROBER NFkB +TNF",
              plot_name = "NFKB_PROBER_+TNF_vs_-TNF",
              additional_labels = c("NFKBIB", "NFKBIA", "EP300", "RELB",
                                    "CREBBP", "NFKBIE", "NFKB2", "NFKB1", "RELA"))




