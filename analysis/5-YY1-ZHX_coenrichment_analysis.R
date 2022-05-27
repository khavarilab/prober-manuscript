
library(GenomicRanges)
library(plyranges)
library(ComplexHeatmap)
library(tidyverse)
library(here)
library(scales)
library(cowplot)
library(eulerr)

out_dir <- file.path("results/ZHX_overlaps", Sys.Date())
dir.create(out_dir, recursive = T, showWarnings = F)


## Data download and wrangling

# DNAse-seq and ChIP-seq peak files were downloaded from the ENCODE portal (Sloan et al., 2016; Davis et al., 2018)
# with the following accession identifiers:
# ENCFF433TIR, ENCFF821KDJ, ENCFF711IED, ENCFF209DJG, ENCFF389ELU, ENCFF671FGG, ENCFF938BND, ENCFF632ZHY, ENCFF371NVU, ENCFF042DUJ, ENCFF933XNU, ENCFF662LUH, ENCFF787LPL, ENCFF101MTI, ENCFF376CAG, ENCFF158NBU, ENCFF041ODC, ENCFF252PLM, and ENCFF286EMA.
# For each cell line, a union peak set was created by merging peaks from all DNAse-seq experiments and ChIP-seq experiments with the selected targets. A union ChIP-seq peak set was created if there were multiple experiments present for each target in each cell line.

encode_metadata <- read_tsv(here("data/encode_yy1_zhx/metadata.tsv"))
encode_data_dir <- here("data/encode_yy1_zhx/")

zhx_cell_lines <- encode_metadata %>%
  filter(str_detect(Assay, "ChIP-seq"), `Output type` == "conservative IDR thresholded peaks") %>%
  group_by(`Biosample term name`) %>%
  filter(any(str_detect(`Experiment target`, "ZHX")),
         any(str_detect(`Experiment target`, "YY1"))) %>%
  ungroup


zhx_chip_peaks <- zhx_cell_lines %>%
  transmute(`File accession`, `Experiment accession`,
            `Experiment target`, `Biosample term name`,
            peaks = map(
              file.path(encode_data_dir, paste0(`File accession`, ".bed.gz")),
              read_narrowpeaks))

dnase_peaks <- encode_metadata %>%
  filter(str_detect(Assay, "DNase-seq"),
         `Biosample term name` %in% zhx_cell_lines$`Biosample term name`) %>%
  mutate(peaks = map(
    file.path(encode_data_dir, paste0(`File accession`, ".bed.gz")),
    read_narrowpeaks))

merged_dnase_peaks <-
  dnase_peaks %>% group_by(`Biosample term name`) %>%
  summarise(dnase_peaks = list(reduce(peaks, function(x, y) GenomicRanges::union(x[x %over% y], y[y %over% x])))) %>%
  ungroup


merged_dnase_peaks <-
  dnase_peaks %>% group_by(`Biosample term name`) %>%
  summarise(dnase_peaks = list(purrr::reduce(peaks, GenomicRanges::union))) %>%
  ungroup

merged_chip_peaks <-
  zhx_chip_peaks %>%
  group_by(`Biosample term name`, `Experiment target`) %>%
  summarise(peaks = list(purrr::reduce(peaks, GenomicRanges::union))) %>%
  ungroup

chip_overlap_with_dnase_peaks <- merged_chip_peaks %>%
  left_join(merged_dnase_peaks) %>%
  transmute(`Biosample term name`, `Experiment target`,
            n_chip_peaks = map_int(peaks, length),
            n_dnase_peaks = map_int(dnase_peaks, length),
            n_overlap_peaks = map2_int(peaks, dnase_peaks, ~ sum(.x %over% .y)),
            frac_overlap = n_overlap_peaks / n_chip_peaks)

### Calculate number for ChIP-seq peaks overlapping DNAse accessible sites

chip_overlap_with_dnase_peaks %>%
  mutate(target = str_replace(`Experiment target`, "-human", "")) %>%
  ggplot(aes(x = target, y = frac_overlap)) +
  facet_grid(~ `Biosample term name`, scales = "free_x", space = "free_x") +
  geom_col() +
  geom_text(aes(label = n_chip_peaks), vjust = - 0.5) +
  scale_y_continuous(limits = c(0, 1.05), expand = c(0, 0)) +
  labs(y = "Fraction of ChIP-seq peaks\noverlaping DNAse-seq peaks") +
  theme_cowplot() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave(file.path(out_dir, "ChIP-seq_peaks_in_accessible_sites.pdf"),
       width = 6, height = 4)

## Calculate pairwise statistics for ChIP-seq peaks in each cell line


biosample_merged_peaks <-
  bind_rows(merged_dnase_peaks %>% rename(peaks = dnase_peaks),
            merged_chip_peaks) %>%
  group_by(`Biosample term name`) %>%
  summarise(merged_peaks = list(purrr::reduce(peaks, GenomicRanges::union))) %>%
  ungroup


pairwise_chip_peaks_overlaps <- merged_chip_peaks %>%
  full_join(., ., by = "Biosample term name") %>%
  left_join(merged_dnase_peaks) %>%
  left_join(biosample_merged_peaks)

pairwise_chip_peaks_summary <- pairwise_chip_peaks_overlaps %>%
  filter(`Experiment target.x` > `Experiment target.y`) %>%
  transmute(`Biosample term name`, `Experiment target.x`, `Experiment target.y`,
            n_dnase_peaks = map_dbl(dnase_peaks, length),
            n_peaks.x = map_dbl(peaks.x, length),
            n_peaks.y = map_dbl(peaks.y, length),
            n_merged_peaks = map_dbl(merged_peaks, length),
            n_overlaps.x = map2_dbl(merged_peaks, peaks.x, ~ sum(.x %over% .y)),
            n_overlaps.y = map2_dbl(merged_peaks, peaks.y, ~ sum(.x %over% .y)),
            n_overlaps.xy = pmap_dbl(list(merged_peaks, peaks.x, peaks.y),
                                     ~ sum(..1 %over% ..2 & ..1 %over% ..3)))

pairwise_chip_peaks_test <- pairwise_chip_peaks_summary %>%
  mutate(fisher_test = pmap(list(n_merged_peaks, n_overlaps.x, n_overlaps.y, n_overlaps.xy),
                            ~ fisher.test(matrix(c(..1 - ..2 - ..3 + ..4, ..2 - ..4,
                                                   ..3 - ..4, ..4), nrow = 2), alternative = "greater") %>%
                              broom::tidy()),
         lambda_pois = n_overlaps.x * n_overlaps.y / n_merged_peaks,
         # pois_enrichment = n_overlaps.xy / pois_lambda,
         log_pvalue_pois = map2_dbl(n_overlaps.xy, lambda_pois,
                                    ~ ppois(.x, .y, lower.tail = FALSE, log.p = T )),
         log10_pvalue_pois = log_pvalue_pois / log(10),
         pvalue_pois = exp(log_pvalue_pois)) %>%
  unnest(fisher_test)

write_csv(pairwise_chip_peaks_test,
          file.path(out_dir, "pairwise_overlap_statistics.csv"))


## Calculate overlap statistics


batch_overlaps <- function(query, subjects, subject_labels) {
  map2_dfc(subjects, subject_labels, ~ list(query %over% .x) %>% set_names(.y))
}


biosample_labeled_peaks <- biosample_merged_peaks %>%
  left_join(merged_chip_peaks %>% group_by(`Biosample term name`) %>%
              summarise(target = list(str_replace(`Experiment target`, "-human", "")),
                        peaks = list(peaks))) %>%
  mutate(overlaps = pmap(list(query = merged_peaks, subjects = peaks, subject_labels = target),
                         batch_overlaps)) %>%
  mutate(euler = map(overlaps, ~ euler(filter_all(., any_vars(.)))),
         euler_labels = map(overlaps, ~ set_names(paste0(colnames(.), "\n(", colSums(.), ")"), colnames(.))))

target_colors <- c(brewer_pal(palette = "Set1")(2)[c(1,2,2)] %>%
                     colorspace::lighten(0.4) %>% set_names(c("YY1", "ZHX1", "ZHX2")),
                   brewer_pal(palette = "Pastel2")(3) %>%
                     colorspace::lighten(0.8) %>% set_names(c("CTCF", "SMC3", "RAD21")))

venn_plots <- biosample_labeled_peaks %>%
  dplyr::rename(Biosample = `Biosample term name`) %>%
  transmute(Biosample,
            plots = pmap(.,
                         function(euler, euler_labels, Biosample, ...) {
                           plot(euler, main = Biosample,
                                labels = euler_labels,
                                # quantities = T,
                                fills = target_colors[names(euler_labels)])
                         }))

pwalk(venn_plots, function(plots, Biosample) {
  pdf(file.path(out_dir, paste0(Biosample, "_euler_venn.pdf")),
      width = 5, height = 5)
  print(plots)
  dev.off()
})


calc_deviation_comb_mat <- function(comb_mat) {

  S <- set_size(comb_mat)
  n <- sum(comb_size(comb_mat))

  df <- tibble(combination_id = comb_name(comb_mat),
               comb_index = str_split(combination_id, "") %>%
                 map(~ as.logical(as.numeric(.))),
               combination = map_chr(comb_index, ~ paste(names(S)[.], collapse = "|")),
               observed = comb_size(comb_mat),
               included_prob = map_dbl(comb_index, ~ prod(S[.]/n)),
               excluded_prob = map_dbl(comb_index, ~ prod(1 - S[!.]/n)),
               expected_prob = included_prob * excluded_prob,
               expected = expected_prob * n,
               lfc = log2(observed / expected),
               deviation = (observed - expected) / expected,
               log_pvalue = ppois(observed, expected, lower.tail = F, log.p = T),
               log10_pvalue = log_pvalue / log(10),
               pvalue = exp(log_pvalue))

  return(df)
}

make_upset_plot <- function(comb_mat, deviation, deviation_calc="lfc", max_comb_size=3) {

  deviation_subset <- deviation[comb_degree(comb_mat) > 0 &
                                  comb_degree(comb_mat) <= max_comb_size]

  deviation_rescaled <- rescale_mid(deviation_subset, mid = 0)
  deviation_colors <-
    gradient_n_pal(rev(brewer_pal(palette = "RdBu")(9)),
                   values = c(-Inf, seq(min(deviation_rescaled), 0.5, length.out = 4)[-4],
                              seq(0.5, max(deviation_rescaled), length.out = 4), Inf))(deviation_rescaled)

  comb_mat_subset <- comb_mat[comb_degree(comb_mat) > 0 & comb_degree(comb_mat) <= max_comb_size]


  if (deviation_calc == "lfc") {
    anno_name = "log2(Obs/Exp)"
  } else if (deviation_calc == "deviation") {
    anno_name = "(Obs-Exp)/Exp"
  }


  do.call("columnAnnotation",
          args = c(set_names(list(anno_barplot(deviation_subset,
                                               gp = gpar(fill = deviation_colors))),
                             anno_name),
                   list(annotation_name_side = "left",
                        annotation_name_rot = 90,
                        annotation_height = unit(40, "mm")))) %v%
    UpSet(comb_mat_subset,
          set_order = order(set_size(comb_mat_subset), decreasing = F),
          comb_order = order(deviation_subset, decreasing = T),
          top_annotation = upset_top_annotation(comb_mat_subset, names = "Intersection size",
                                                annotation_name_rot = 90,
                                                height = unit(30, "mm")))
}

biosample_peaks_stats <- biosample_labeled_peaks %>%
  dplyr::rename(Biosample = `Biosample term name`) %>%
  mutate(comb_mat = map(overlaps, make_comb_mat),
         stats = map(comb_mat, calc_deviation_comb_mat))

biosample_peaks_stats %>%
  pwalk(function(Biosample, stats, ...) {
    stats %>% select(-combination_id, - comb_index) %>%
      write_csv(file.path(out_dir, paste0(Biosample, "_upset_stats.csv")))
  })

biosample_peaks_stats %>%
  mutate(deviation = map(stats, ~ .[["lfc"]])) %>%
  pwalk(function(overlaps, comb_mat, deviation, Biosample, ...) {
    pdf(file.path(out_dir, paste0(Biosample, "_upset_plot_lfc.pdf")),
        width = 5, height = 5)
    print(make_upset_plot(comb_mat, deviation,
                          deviation_calc = "lfc"))
    dev.off()
  })


biosample_peaks_stats %>%
  mutate(deviation = map(stats, ~ .[["deviation"]])) %>%
  pwalk(function(overlaps, comb_mat, deviation, Biosample, ...) {
    pdf(file.path(out_dir, paste0(Biosample, "_upset_plot_deviation.pdf")),
        width = 5, height = 5)
    print(make_upset_plot(comb_mat, deviation,
                          deviation_calc = "deviation"))
    dev.off()
  })



