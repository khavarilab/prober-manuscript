library(BSgenome.Hsapiens.UCSC.hg38)
library(TFBSTools)
library(plyranges)
library(here)
library(magrittr)
library(tidyverse)


out_dir <- here(file.path("results/endogenous_YY1_sites_293T", Sys.Date()))
dir.create(out_dir, showWarnings = F, recursive = T)

cores <- 12

YY1_binding_seqs <- DNAStringSet(
  c("seqA" = "GCCGCCATCTTG",
    "seqB" = "GCCGCCATTTTG",
    "seqC" = "CCGGCCATCTTG"
  ))

YY1_binding_seqs_revCom <- reverseComplement(YY1_binding_seqs)

YY1_chip_peaks <- read_narrowpeaks(here("data/ENCODE_YY1_293T/ENCFF395WOT.bed.gz"),
                                   genome_info = "hg38")

YY1_chip_signal <- read_bigwig(here("data/ENCODE_YY1_293T/ENCFF258PGM.bigWig"),
                               genome_info = "hg38")

YY1_chip_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, YY1_chip_peaks)

YY1_sequence_matches_for <- vcountPDict(YY1_binding_seqs, YY1_chip_seqs) %>%
  set_rownames(names(YY1_binding_seqs))

YY1_sequence_matches_rev <- vcountPDict(YY1_binding_seqs_revCom, YY1_chip_seqs) %>%
  set_rownames(names(YY1_binding_seqs_revCom))

YY1_sequence_matches <- YY1_sequence_matches_for + YY1_sequence_matches_rev

YY1_peaks_matching_seq <- YY1_chip_peaks[colSums(YY1_sequence_matches) > 0]
YY1_seqs_matching <- YY1_chip_seqs[colSums(YY1_sequence_matches) > 0]


YY1_peaks_table <- YY1_peaks_matching_seq %>% as_tibble() %>%
  mutate(peakSequence = as.character(YY1_seqs_matching),
         YY1_peak = as.character(row_number()))

YY1_position_matches_for <- map_dfr(YY1_peaks_table$peakSequence, ~ matchPDict(YY1_binding_seqs, DNAString(.)) %>%
                                  as.list() %>%
                                  map_dfr(as.data.frame, .id = "YY1_seq_name"),
                                .id = "YY1_peak") %>%
  select(YY1_peak, YY1_seq_name, YY1_seq_start = start, YY1_seq_end = end) %>%
  mutate(YY1_seq_strand = "+")

YY1_position_matches_rev <- map_dfr(YY1_peaks_table$peakSequence, ~ matchPDict(YY1_binding_seqs_revCom, DNAString(.)) %>%
                                      as.list() %>%
                                      map_dfr(as.data.frame, .id = "YY1_seq_name"),
                                    .id = "YY1_peak") %>%
  select(YY1_peak, YY1_seq_name, YY1_seq_start = start, YY1_seq_end = end) %>%
  mutate(YY1_seq_strand = "-")

flank_len <- 12
YY1_results_table <- left_join(YY1_peaks_table,
                               bind_rows(YY1_position_matches_for, YY1_position_matches_rev)) %>%
  mutate(YY1_binding_seq = as.character(getSeq(BSgenome.Hsapiens.UCSC.hg38, names = seqnames,
                                  start = start + YY1_seq_start - 1 - flank_len,
                                  end = start + YY1_seq_end -1 + flank_len)) %>% set_names(NULL)) %>%
  mutate(YY1_seq_id = as.character(row_number()))


yy1_gr <- YY1_results_table %>% transmute(YY1_seq_id,
                                          seqnames,
                                          peak_start = start,
                                          peak_end = end,
                                          start = peak_start + YY1_seq_start - 1 - flank_len,
                                          end = peak_start + YY1_seq_end - 1 + flank_len) %>%
  as_granges()


yy1_signal <- rtracklayer::import(here("data/ENCODE_YY1_293T/ENCFF258PGM.bigWig"), format = "bigWig",
       selection = BigWigSelection(yy1_gr))

yy1_binding_signal <- yy1_gr %>%
  join_overlap_intersect(yy1_signal) %>%
  plyranges::group_by(YY1_seq_id) %>%
  plyranges::summarise(
    meanSignal = weighted.mean(score, width),
    maxSignal = max(score)) %>%
  as_tibble()


homer_pwms <- motifbreakR::homer %>% convert_motifs("TFBSTools-PWMatrix") %>%
  set_names(map_chr(., ~ name(.))) %>% do.call(PWMatrixList, .)

hocomoco_pwms <- motifbreakR::hocomoco %>% convert_motifs("TFBSTools-PWMatrix") %>%
  set_names(map_chr(., ~ name(.))) %>% do.call(PWMatrixList, .)

hocomoco_search <- searchSeq(hocomoco_pwms, DNAStringSet(YY1_results_table$YY1_binding_seq %>% set_names(YY1_results_table$YY1_seq_id)),
                             min.score = 0.9, mc.cores = cores) %>% as.list()

hocomoco_results <- hocomoco_search[map_lgl(hocomoco_search, ~ length(.) > 0)] %>%
  map_dfr(as.data.frame) %>%
  select(YY1_seq_id = seqnames, motif_start = start, motif_end = end, motif_strand = strand,
         absScore, relScore, motif_ID = ID, motif_TF= TF, motif_seq = siteSeqs)


YY1_results_table_motifs <- left_join(YY1_results_table, yy1_binding_signal) %>%
  left_join(hocomoco_results) %>%
  select(YY1_peak, YY1_seq_id, YY1_binding_seq, meanSignal, maxSignal,
         chr = seqnames, start, end, width, signalValue, qValue,
         YY1_seq_start, YY1_seq_end, YY1_seq_strand,
         motif_ID, motif_TF, absScore, relScore, motif_start, motif_end, motif_strand, motif_seq, peakSequence)

write_csv(YY1_results_table_motifs,
          file.path(out_dir, "YY1_binding_sites_293T_with_motifs.csv"))

