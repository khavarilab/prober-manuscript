saint_vs_foldchange_scatterplot <- function(data, protein_labels, plot_path,
                                            fc_axis_limit = NA,
                                            saint_threshold = 0.9, fc_threshold = 3) {
  ref_column <- which(colnames(data) == "BAIT_iREF")

  data <- data %>%
    mutate(COUNT_1 = data[[ref_column+1]],
           COUNT_2 = data[[ref_column+2]],
           above_threshold = BAIT_SP >= saint_threshold &
             BAIT_FC_A >= fc_threshold,
           group = case_when(above_threshold & PROTID %in% protein_labels ~ 3,
                             above_threshold ~ 2,
                             TRUE ~ 1),
           label = case_when(group == 3 ~ PROTID,
                             group == 2 ~ "",
                             TRUE ~ NA_character_)) %>%
    filter(COUNT_1 > 0 & COUNT_2 > 0) %>%
    arrange(group)

  fc_axis_limit <- ifelse(is.na(fc_axis_limit),
                          max(data$BAIT_FC_A), fc_axis_limit)

  ggplot(data, aes(y = BAIT_SP, x = BAIT_FC_A, label = label)) +
    geom_hline(yintercept = saint_threshold, color = "grey50", lty = 2) +
    geom_point(aes(color = as.factor(group), size = !is.na(label)),
               show.legend = F) +
    geom_text_repel(size = 3, fontface = "bold", show.legend = F,
                    nudge_y = -0.01, max.overlaps = 1000, min.segment.length = 0) +
    labs(x = "Fold Change",
         y = "SAINT") +
    scale_color_manual(values = c("grey70", ggsci::pal_npg()(2)[2:1])) +
    scale_size_manual(values = c(1.5, 2.5)) +
    scale_x_continuous(limits = c(0, fc_axis_limit)) +
    theme_cowplot()
  ggsave(plot_path, width = 6, height = 5)

}




contrast_scatter_plot <- function(data, contrast, xlab, ylab, protein_labels, plot_path,
                                  x_data = 1, saint_threshold = 0.9, fc_threshold = 3) {

  ref_column_1 <- which(colnames(data[[1]]) == "BAIT_iREF")
  ref_column_2 <- which(colnames(data[[2]]) == "BAIT_iREF")


  contrast_data <- full_join(
    data[[contrast[1]]] %>%
      mutate(MIN_COUNT_1 = pmin(data[[contrast[1]]][[ref_column_1+1]],
                                data[[contrast[1]]][[ref_column_1+2]], na.rm = T)) %>%
      select(PROTID = PROTID,
             BAIT_FC_1 = BAIT_FC_A,
             BAIT_SP_1 = BAIT_SP,
             MIN_COUNT_1),
    data[[contrast[2]]] %>%
      mutate(MIN_COUNT_2 = pmin(data[[contrast[2]]][[ref_column_2+1]],
                                data[[contrast[2]]][[ref_column_2+2]], na.rm = T)) %>%
      select(PROTID = PROTID,
             BAIT_FC_2 = BAIT_FC_A,
             BAIT_SP_2 = BAIT_SP,
             MIN_COUNT_2),
    by = "PROTID") %>%
    mutate(across(starts_with("BAIT_FC_"), ~ ifelse(is.na(.), 1, .)),
           across(starts_with("BAIT_SP_"), ~ ifelse(is.na(.), 0, .)),
           across(starts_with("MIN_COUNT_"), ~ ifelse(is.na(.), 0, .)),
           max_SAINT = pmax(BAIT_SP_1, BAIT_SP_2, na.rm = T),
           max_FC = pmax(BAIT_FC_1, BAIT_FC_2, na.rm = T),
           above_threshold = max_FC >= fc_threshold &
             max_SAINT >= saint_threshold,
           group = case_when(above_threshold & PROTID %in% protein_labels ~ 3,
                             above_threshold ~ 2,
                             TRUE ~ 1),
           label = case_when(group == 3 ~ PROTID,
                             group == 2 ~ "",
                             TRUE ~ NA_character_)) %>%
    filter(MIN_COUNT_1 > 0 | MIN_COUNT_2 > 0) %>%
    arrange(group)


  x_data <- sym(paste0("BAIT_FC_", x_data))

  ggplot(contrast_data, aes(x = !!x_data, y = log2(BAIT_FC_2 / BAIT_FC_1),
                            label = label)) +
    geom_hline(yintercept = 0, color = "grey50", lty = 2) +
    geom_point(aes(color = as.factor(group), size = !is.na(label)),
               show.legend = F) +
    geom_text_repel(size = 3, fontface = "bold", show.legend = F, max.overlaps = 1000, min.segment.length = 0) +
    labs(x = xlab,
         y = ylab) +
    scale_color_manual(values = c("grey70", ggsci::pal_npg()(2)[2:1])) +
    scale_size_manual(values = c(1.5, 2.5)) +
    theme_cowplot()
  ggsave(plot_path, width = 6, height = 5)

}


contrast_scatter_plot_MA <- function(data, contrast, xlab, ylab, protein_labels, plot_path,
                                  x_data = 1, saint_threshold = 0.9, fc_threshold = 3,
                                  show_background = T) {

  ref_column_1 <- which(colnames(data[[1]]) == "BAIT_iREF")
  ref_column_2 <- which(colnames(data[[2]]) == "BAIT_iREF")


  contrast_data <- full_join(
    data[[contrast[1]]] %>%
      mutate(MIN_COUNT_1 = pmin(data[[contrast[1]]][[ref_column_1+1]],
                                data[[contrast[1]]][[ref_column_1+2]], na.rm = T)) %>%
      select(PROTID = PROTID,
             BAIT_FC_1 = BAIT_FC_A,
             BAIT_SP_1 = BAIT_SP,
             MIN_COUNT_1),
    data[[contrast[2]]] %>%
      mutate(MIN_COUNT_2 = pmin(data[[contrast[2]]][[ref_column_2+1]],
                                data[[contrast[2]]][[ref_column_2+2]], na.rm = T)) %>%
      select(PROTID = PROTID,
             BAIT_FC_2 = BAIT_FC_A,
             BAIT_SP_2 = BAIT_SP,
             MIN_COUNT_2),
    by = "PROTID") %>%
    mutate(across(starts_with("BAIT_FC_"), ~ ifelse(is.na(.), 1, .)),
           across(starts_with("BAIT_SP_"), ~ ifelse(is.na(.), 0, .)),
           across(starts_with("MIN_COUNT_"), ~ ifelse(is.na(.), 0, .)),
           max_SAINT = pmax(BAIT_SP_1, BAIT_SP_2, na.rm = T),
           max_FC = pmax(BAIT_FC_1, BAIT_FC_2, na.rm = T),
           above_threshold = max_FC >= fc_threshold &
             max_SAINT >= saint_threshold,
           group = case_when(above_threshold & PROTID %in% protein_labels ~ 3,
                             above_threshold ~ 2,
                             TRUE ~ 1),
           label = case_when(group == 3 ~ PROTID,
                             group == 2 ~ "",
                             TRUE ~ NA_character_)) %>%
    filter(MIN_COUNT_1 > 0 | MIN_COUNT_2 > 0) %>%
    arrange(group)

  if (!show_background) {
    contrast_data <- contrast_data %>%
      filter(above_threshold)
  }

  ggplot(contrast_data, aes(x = sqrt(BAIT_FC_1 * BAIT_FC_2), y = log2(BAIT_FC_2 / BAIT_FC_1),
                            label = label)) +
    geom_hline(yintercept = 0, color = "grey50", lty = 2) +
    geom_point(aes(color = as.factor(group), size = !is.na(label)),
               show.legend = F) +
    geom_text_repel(size = 3, fontface = "bold", show.legend = F, max.overlaps = 1000, min.segment.length = 0) +
    labs(x = xlab,
         y = ylab) +
    scale_color_manual(values = c("grey70", ggsci::pal_npg()(2)[2:1]) %>% set_names(1:3)) +
    scale_size_manual(values = c(1.5, 2.5) %>% set_names(c(F, T))) +
    expand_limits(x = 0) +
    theme_cowplot()
  ggsave(plot_path, width = 6, height = 5)

}




contrast_scatter_plot_saint <- function(data, contrast, xlab, ylab, protein_labels, plot_path,
                                  x_data = 1, saint_threshold = 0.9, fc_threshold = 3) {

  ref_column_1 <- which(colnames(data[[1]]) == "BAIT_iREF")
  ref_column_2 <- which(colnames(data[[2]]) == "BAIT_iREF")


  contrast_data <- full_join(
    data[[contrast[1]]] %>%
      mutate(MIN_COUNT_1 = pmin(data[[contrast[1]]][[ref_column_1+1]],
                                data[[contrast[1]]][[ref_column_1+2]], na.rm = T)) %>%
      select(PROTID = PROTID,
             BAIT_FC_1 = BAIT_FC_A,
             BAIT_SP_1 = BAIT_SP,
             MIN_COUNT_1),
    data[[contrast[2]]] %>%
      mutate(MIN_COUNT_2 = pmin(data[[contrast[2]]][[ref_column_2+1]],
                                data[[contrast[2]]][[ref_column_2+2]], na.rm = T)) %>%
      select(PROTID = PROTID,
             BAIT_FC_2 = BAIT_FC_A,
             BAIT_SP_2 = BAIT_SP,
             MIN_COUNT_2),
    by = "PROTID") %>%
    mutate(across(starts_with("BAIT_FC_"), ~ ifelse(is.na(.), 1, .)),
           across(starts_with("BAIT_SP_"), ~ ifelse(is.na(.), 0, .)),
           across(starts_with("MIN_COUNT_"), ~ ifelse(is.na(.), 0, .)),
           max_SAINT = pmax(BAIT_SP_1, BAIT_SP_2, na.rm = T),
           max_FC = pmax(BAIT_FC_1, BAIT_FC_2, na.rm = T),
           above_threshold = max_FC >= fc_threshold &
             max_SAINT >= saint_threshold,
           group = case_when(above_threshold & PROTID %in% protein_labels ~ 3,
                             above_threshold ~ 2,
                             TRUE ~ 1),
           label = case_when(group == 3 ~ PROTID,
                             group == 2 ~ "",
                             TRUE ~ NA_character_)) %>%
    filter(MIN_COUNT_1 > 0 | MIN_COUNT_2 > 0) %>%
    arrange(group)


  x_data <- sym(paste0("BAIT_SP_", x_data))

  ggplot(contrast_data, aes(x = !!x_data, y = log2(BAIT_FC_2 / BAIT_FC_1),
                            label = label)) +
    geom_hline(yintercept = 0, color = "grey50", lty = 2) +
    geom_point(aes(color = as.factor(group), size = !is.na(label)),
               show.legend = F) +
    geom_text_repel(size = 3, fontface = "bold", show.legend = F, max.overlaps = 1000, min.segment.length = 0) +
    labs(x = xlab,
         y = ylab) +
    scale_color_manual(values = c("grey70", ggsci::pal_npg()(2)[2:1])) +
    scale_size_manual(values = c(1.5, 2.5)) +
    theme_cowplot()
  ggsave(plot_path, width = 6, height = 5)

}

#
# contrast_limma_analysis <-
#   function(data, contrast_labels,
#            condition1 = NA, condition2 = NA,
#            bait1, bait2, plot_path,
#            quant_norm = T,
#            contrast_name = "baitBAIT2.conditionCOND2",
#            protein_labels, protein_highlights = c(),
#            fdr_cutoff = 0.2,
#            p_cutoff = NA,
#            saint_cutoff = 0.9,
#            logFC_label = "log2FC", pval_label = "-log10 p-value",
#            meanFC_label = "Average FC") {
#
#     if (!is.na(condition1)) {
#       condition1_bait1_expr <- paste0("IP_", bait1, "_", condition1,  "[_-]?R\\d")
#       condition1_bait2_expr <- paste0("IP_", bait2, "_", condition1,  "[_-]?R\\d")
#     } else {
#       condition1_bait1_expr <- paste0("IP_", bait1, "[_-]?R\\d")
#       condition1_bait2_expr <- paste0("IP_", bait2, "[_-]?R\\d")
#     }
#
#     if (!is.na(condition2)) {
#       condition2_bait1_expr <- paste0("IP_", bait1, "_", condition2,  "[_-]?R\\d")
#       condition2_bait2_expr <- paste0("IP_", bait2, "_", condition2,  "[_-]?R\\d")
#     } else {
#       condition2_bait1_expr <- paste0("IP_", bait1, "[_-]?R\\d")
#       condition2_bait2_expr <- paste0("IP_", bait2, "[_-]?R\\d")
#     }
#
#
#
#
#
#     #
#     # condition1_bait_expr <- "IP_NFKB_WITHTNF_R\\d"
#     # condition2_bait_expr <- "IP_NFKB_NOTNF_R\\d"
#     # condition1_scr_expr <- "IP_S\\d_WITHTNF_R\\d"
#     # condition2_scr_expr <- "IP_S\\d_NOTNF_R\\d"
#
#
#     contrast_data <- full_join(
#       data[[contrast_labels[1]]] %>%
#         dplyr::rename_with(function(x) paste0("COUNTS_BAIT1_COND1_", seq_along(x)),
#                            matches(condition1_bait1_expr)) %>%
#         dplyr::rename_with(function(x) paste0("COUNTS_BAIT2_COND1_", seq_along(x)),
#                            matches(condition1_bait2_expr)) %>%
#         dplyr::select(PROTID, starts_with("COUNTS"),
#                       BAIT_FC_1 = BAIT_FC_A, BAIT_SP_1 = BAIT_SP),
#       data[[contrast_labels[2]]] %>%
#         dplyr::rename_with(function(x) paste0("COUNTS_BAIT1_COND2_", seq_along(x)),
#                            matches(condition2_bait1_expr)) %>%
#         dplyr::rename_with(function(x) paste0("COUNTS_BAIT2_COND2_", seq_along(x)),
#                            matches(condition2_bait2_expr)) %>%
#         dplyr::select(PROTID, starts_with("COUNTS_"),
#                       BAIT_FC_2 = BAIT_FC_A, BAIT_SP_2 = BAIT_SP)) %>%
#       mutate(across(starts_with("COUNTS_"), ~ ifelse(is.na(.), 0, .)),
#              across(starts_with("BAIT_FC_"), ~ ifelse(is.na(.), 1, .)))
#
#     if (is.na(condition2)) {
#       contrast_data <- contrast_data %>%
#         select(!starts_with("COUNTS_BAIT1_COND2_"))
#     }
#
#     count_data <- contrast_data %>%
#       select(PROTID, starts_with("COUNTS_")) %>%
#       as.data.frame() %>%
#       magrittr::set_rownames(.$PROTID) %>%
#       select(-PROTID) %>% as.matrix()
#
#
#     if (! is.na(condition2)) {
#       sample_data <- tibble(sample = colnames(count_data)) %>%
#         extract(sample, c("bait", "condition"),
#                 "COUNTS_(.+)_(.+)_\\d", remove = F) %>%
#         mutate(bait = factor(bait, levels = c("BAIT1", "BAIT2")),
#                condition = factor(condition, levels = c("COND1", "COND2"))) %>%
#         as.data.frame() %>%
#         magrittr::set_rownames(.$sample) %>%
#         select(-sample)
#
#       design <- model.matrix(~ bait * condition + 0, sample_data)
#       colnames(design) <- str_replace(colnames(design), ":", ".")
#
#     } else {
#       sample_data <- tibble(sample = colnames(count_data)) %>%
#         mutate(condition = case_when(str_detect(sample, "BAIT1_COND1") ~ "CTRL",
#                                      str_detect(sample, "BAIT2_COND1") ~ "COND1",
#                                      str_detect(sample, "BAIT2_COND2") ~ "COND2")) %>%
#         as.data.frame() %>%
#         magrittr::set_rownames(.$sample) %>%
#         select(-sample)
#
#       design <- model.matrix(~ condition + 0, sample_data)
#
#     }
#
#
#     data_obj <- DGEList(count_data)
#
#     keep <- filterByExpr(data_obj, design = design,
#                  min.count = 1, min.total.count = 10)
#     data_obj <- data_obj[keep,]
#
#     # count_data <- count_data[rowSums(count_data > 0) > 1,]
#     # count_data <- count_data[rowSums(count_data) > 5, ]
#
#
#     data_obj <- calcNormFactors(data_obj)
#
#     data_obj <- cpm(data_obj, log = T, prior.count = 10)
#
#
#     if (quant_norm) {
#       # log_transformed_data <- log(data_obj + 1)
#
#       data_obj <- preprocessCore::normalize.quantiles(data_obj)
#
#       # data_obj <- DGEList(normalized_log_transformed_data)
#     }
#
#
#
#     contrast <- makeContrasts(contrasts = contrast_name,
#                               levels = design)
#
#
#     fit1 <- lmFit(data_obj, design = design)
#
#     fit2 <- contrasts.fit(fit1, contrasts = contrast)
#
#     fit3 <- eBayes(fit2, trend = T, proportion = 0.1)
#
#     limma_results <- topTable(fit3, number = Inf) %>%
#       as_tibble(rownames = "PROTID")
#
#     if( length(setdiff(protein_highlights, limma_results$PROTID)) > 0) {
#       print(setdiff(protein_highlights, limma_results$PROTID))
#     }
#
#     results <-  limma_results %>%
#       left_join(contrast_data) %>%
#       mutate(above_threshold = (!is.na(fdr_cutoff) & adj.P.Val <= fdr_cutoff) |
#                (!is.na(p_cutoff) & P.Value < p_cutoff),
#              group = case_when(above_threshold & PROTID %in% protein_labels ~ 3,
#                                above_threshold ~ 2,
#                                TRUE ~ 1) %>% factor(levels = 1:3),
#              label = case_when(group == 3  &
#                                  (BAIT_SP_1 >saint_cutoff | BAIT_SP_2 > saint_cutoff) ~ PROTID,
#                                group == 2 | group == 3 ~ "",
#                                TRUE ~ NA_character_)) %>%
#       arrange(group)
#
#     write_csv(results, paste0(plot_path, "_limma-results.csv"))
#
#     volcano_gg <- ggplot(results, aes(x = logFC, y = -log10(P.Value),
#                         label = label)) +
#       geom_vline(xintercept = 0, color = "grey50", lty = 2) +
#       geom_point(aes(color = group, size = !is.na(label)),
#                  show.legend = F) +
#       geom_text_repel(size = 3, fontface = "bold", show.legend = F, max.overlaps = 1000, min.segment.length = 0) +
#       scale_color_manual(values = c("grey70", ggsci::pal_npg()(2)[c(2,1)]) %>% set_names(1:3)) +
#       scale_size_manual(values = c(1.5, 2.5) %>% set_names(c(F, T))) +
#       labs(x = logFC_label,
#            y = pval_label) +
#       theme_cowplot() +
#       ggsave(paste0(plot_path, "_limma-volcano.pdf"), width = 6, height = 5)
#
#     ma_gg <- ggplot(results,
#                     aes(x = sqrt(BAIT_FC_1 * BAIT_FC_2),
#                         y = logFC,
#                                         label = label)) +
#                       geom_hline(yintercept = 0, color = "grey50", lty = 2) +
#                       geom_point(aes(color = group, size = !is.na(label)),
#                                  show.legend = F) +
#                       geom_text_repel(size = 3, fontface = "bold", show.legend = F, max.overlaps = 1000, min.segment.length = 0) +
#       scale_color_manual(values = c("grey70", ggsci::pal_npg()(2)[c(2,1)]) %>% set_names(1:3)) +
#       scale_size_manual(values = c(1.5, 2.5) %>% set_names(c(F, T))) +
#                       labs(x = meanFC_label,
#                            y = logFC_label) +
#                       theme_cowplot() +
#                       ggsave(paste0(plot_path, "_limma-MA.pdf"), width = 6, height = 5)
#
#     return(list(volcano = volcano_gg,
#                 ma = ma_gg))
#
#   }
#



