
contrast_limma_analysis <-
  function(data, contrast_labels,
           condition1, condition2,
           bait1 = NA, bait2 = NA,
           controls = NA,
           plot_path,
           quant_norm = F,
           contrast_name = ifelse(is.na(bait1) || is.na(bait2),
                                  "conditionCOND2 - conditionCOND1",
                                  "baitBAIT2.conditionCOND2"),
           protein_labels = c(), protein_highlights = c(),
           only_label_highlights = T,
           min_protein_count  = 1,
           pseudo_count = 1,
           fdr_cutoff = 0.25,
           fdr_line = 0.25,
           p_cutoff = NA,
           p_limit = NA,
           saint_cutoff = 0.9,
           label_size = 3,
           point_sizes = c(1.5, 2.5),
           logFC_label = "log2FC", pval_label = "-log10 p-value",
           meanFC_label = "Average FC over scramble controls") {


    saint_results <- full_join(
        data[[contrast_labels[1]]] %>%
          dplyr::select(PROTID,
                        BAIT_FC_1 = BAIT_FC_A, BAIT_SP_1 = BAIT_SP),
        data[[contrast_labels[2]]] %>%
          dplyr::select(PROTID,
                        BAIT_FC_2 = BAIT_FC_A, BAIT_SP_2 = BAIT_SP)) %>%
      mutate(across(starts_with("BAIT_FC_"), ~ ifelse(is.na(.) | . < 1, 1, .)),
             across(starts_with("BAIT_SP_"), ~ ifelse(is.na(.), 0, .)))

    proteins_fc1 <- saint_results %>%
      filter(BAIT_FC_1 >= 1.5 | BAIT_FC_2 >= 1.5) %>% pull(PROTID)

    contrast_data <- full_join(
      data[[contrast_labels[1]]] %>%
        dplyr::select(PROTID, starts_with("IP_")),
      data[[contrast_labels[2]]] %>%
        dplyr::select(PROTID, starts_with("IP_"))) %>%
      mutate(across(starts_with("IP_"), ~ ifelse(is.na(.), 0, .)))


    if (any(duplicated(contrast_data$PROTID))) {
      stop("Non-identical count columns with the same name, or multiple rows with the same PROTIDs")
    }






    if (!is.na(bait1) && !is.na(bait2)) {
      condition1_bait1_expr <- paste0("IP_", bait1, "_", condition1,  "[_-]?R\\d")
      condition1_bait2_expr <- paste0("IP_", bait2, "_", condition1,  "[_-]?R\\d")
      condition2_bait1_expr <- paste0("IP_", bait1, "_", condition2,  "[_-]?R\\d")
      condition2_bait2_expr <- paste0("IP_", bait2, "_", condition2,  "[_-]?R\\d")

      contrast_data <- contrast_data %>%
        dplyr::rename_with(function(x) paste0("COUNTS_BAIT1_COND1_", seq_along(x)),
                         matches(condition1_bait1_expr)) %>%
        dplyr::rename_with(function(x) paste0("COUNTS_BAIT2_COND1_", seq_along(x)),
                           matches(condition1_bait2_expr)) %>%
        dplyr::rename_with(function(x) paste0("COUNTS_BAIT1_COND2_", seq_along(x)),
                           matches(condition2_bait1_expr)) %>%
        dplyr::rename_with(function(x) paste0("COUNTS_BAIT2_COND2_", seq_along(x)),
                           matches(condition2_bait2_expr)) %>%
        dplyr::select(PROTID, starts_with("COUNTS_"))


      sample_data <- tibble(sample = colnames(contrast_data)) %>%
        extract(sample, c("bait", "condition"),
                "COUNTS_(.+)_(.+)_\\d", remove = F) %>%
        mutate(bait = factor(bait, levels = c("BAIT1", "BAIT2")),
               condition = factor(condition, levels = c("COND1", "COND2"))) %>%
        as.data.frame() %>%
        magrittr::set_rownames(.$sample) %>%
        select(-sample)

      design <- model.matrix(~ bait * condition + 0, sample_data)
      colnames(design) <- str_replace(colnames(design), ":", ".")

    } else {
      condition1_expr <- paste0("IP_", condition1, "[_-]?R\\d")
      condition2_expr <- paste0("IP_", condition2, "[_-]?R\\d")

      contrast_data <- contrast_data %>%
        dplyr::rename_with(function(x) paste0("COUNTS_COND1_", seq_along(x)),
                           matches(condition1_expr)) %>%
        dplyr::rename_with(function(x) paste0("COUNTS_COND2_", seq_along(x)),
                           matches(condition2_expr))

      if (!is.na(controls)) {
        controls_expr <- paste0("^IP_", controls)
        contrast_data <- contrast_data %>%
          dplyr::rename_with(function(x) paste0("COUNTS_CTRLS_", seq_along(x)),
                             matches(controls_expr))
      }

      contrast_data <- contrast_data %>%
        dplyr::select(PROTID, starts_with("COUNTS_"))

      sample_data <- tibble(sample = colnames(contrast_data)) %>%
        mutate(condition = case_when(str_detect(sample, "CTRLS") ~ "CTRLS",
                                     str_detect(sample, "COND1") ~ "COND1",
                                     str_detect(sample, "COND2") ~ "COND2")) %>%
        as.data.frame() %>%
        magrittr::set_rownames(.$sample) %>%
        select(-sample)

      design <- model.matrix(~ condition + 0, sample_data)
    }

    count_matrix <- contrast_data %>%
      select(PROTID, starts_with("COUNTS_")) %>%
      as.data.frame() %>%
      magrittr::set_rownames(.$PROTID) %>%
      select(-PROTID) %>% as.matrix()

    data_obj <- DGEList(count_matrix)

    keep_expressed <- filterByExpr(data_obj, design = design,
                         min.count = 0, min.total.count = min_protein_count)
    keep_fc <- rownames(data_obj) %in% proteins_fc1

    data_obj <- calcNormFactors(data_obj)

    data_obj <- data_obj[keep_expressed & keep_fc, , keep.lib.sizes = T]

    # count_data <- count_data[rowSums(count_data > 0) > 1,]
    # count_data <- count_data[rowSums(count_data) > 5, ]



    data_obj <- cpm(data_obj, log = T, prior.count = pseudo_count)


    if (quant_norm) {
      # log_transformed_data <- log(data_obj + 1)

      data_obj <- preprocessCore::normalize.quantiles(data_obj) %>%
        magrittr::set_colnames(colnames(data_obj)) %>%
        magrittr::set_rownames(rownames(data_obj))

      # data_obj <- DGEList(normalized_log_transformed_data)
    }



    contrast <- makeContrasts(contrasts = contrast_name,
                              levels = design)


    fit1 <- lmFit(data_obj, design = design)

    fit2 <- contrasts.fit(fit1, contrasts = contrast)

    fit3 <- eBayes(fit2, trend = T, robust = T)

    pdf(paste0(plot_path, "_limma-SAplot.pdf"), width = 6, height = 6)
    plotSA(fit3)
    dev.off()

    if (!is.na(p_cutoff)) {
      MA_highlights <- fit3$p.value <= p_cutoff
    } else {
      MA_highlights <- p.adjust(fit3$p.value, method = 'fdr') <= fdr_cutoff
    }


    pdf(paste0(plot_path, "_limma-MAplot.pdf"), width = 6, height = 6)
    plotMD(fit3, status = MA_highlights, values = T)
    dev.off()

    limma_results <- topTable(fit3, number = Inf) %>%
      as_tibble(rownames = "PROTID")

    if( length(setdiff(protein_labels, limma_results$PROTID)) > 0) {
      print(setdiff(protein_labels, limma_results$PROTID))
    }

    results <-  limma_results %>%
      left_join(contrast_data) %>%
      left_join(saint_results) %>%
      mutate(above_threshold = (!is.na(fdr_cutoff) & adj.P.Val <= fdr_cutoff) |
               (!is.na(p_cutoff) & P.Value <= p_cutoff),
             group = case_when(above_threshold & PROTID %in% protein_highlights ~ 3,
                               above_threshold ~ 2,
                               PROTID %in% protein_labels & PROTID %in% protein_highlights ~ 3,
                               PROTID %in% protein_labels & ! PROTID %in% protein_highlights ~ 2,
                               TRUE ~ 1) %>% factor(levels = 1:4),
             label = case_when(PROTID %in% protein_labels ~ PROTID,
                               group == 3  &
                                 (BAIT_SP_1 >= saint_cutoff | BAIT_SP_2 >= saint_cutoff) ~ PROTID,
                               !only_label_highlights & above_threshold & group == 2 ~ PROTID,
                               group == 4 ~ PROTID,
                               group == 2 | group == 3 ~ "",
                               TRUE ~ NA_character_)) %>%
      arrange(group)

    write_csv(results, paste0(plot_path, "_limma-results.csv"))



    if (!is.na(fdr_line)) {
      below_fdr_line <- results %>%
        distinct(P.Value, adj.P.Val) %>%
        filter(adj.P.Val < fdr_line) %>%
        slice_max(P.Value, n = 1)
      above_fdr_line <- results %>%
        distinct(P.Value, adj.P.Val) %>%
        filter(adj.P.Val >= fdr_line) %>%
        slice_min(P.Value, n = 1)

      fdr_line_p <- (below_fdr_line$P.Value + above_fdr_line$P.Value) / 2
    } else {
      fdr_line_p <- NA
    }

    volcano_gg <- results %>%
      mutate(label = ifelse(!is.na(p_limit) & P.Value < p_limit,
                            paste0(label, "\n(p = ",
                                   scientific(P.Value, digits = 2), ")"),
                            label),
             P.Value = ifelse(!is.na(p_limit) & P.Value < p_limit,
                              p_limit, P.Value)) %>%
      ggplot(aes(x = logFC, y = -log10(P.Value),
                                      label = label)) +
      geom_vline(xintercept = 0, color = "grey50", lty = 2) +
      geom_hline(yintercept = -log10(fdr_line_p), color = "grey50", lty = 2) +
      geom_point(aes(color = group, size = !is.na(label)),
                 show.legend = F) +
      geom_text_repel(size = label_size, fontface = "bold", show.legend = F, max.overlaps = 1000,
                      min.segment.length = 0, segment.size = 0.2) +
      scale_color_manual(values = c("grey70", ggsci::pal_npg()(5)[c(2,1,5)]) %>% set_names(1:4)) +
      scale_size_manual(values = point_sizes %>% set_names(c(F, T))) +
      labs(x = logFC_label,
           y = pval_label) +
      theme_cowplot()
    ggsave(paste0(plot_path, "_limma-volcano.pdf"), width = 6, height = 5)

    ma_gg <- ggplot(results,
                    aes(x = sqrt(BAIT_FC_1 * BAIT_FC_2),
                        y = logFC,
                        label = label)) +
      geom_hline(yintercept = 0, color = "grey50", lty = 2) +
      geom_point(aes(color = group, size = !is.na(label)),
                 show.legend = F) +
      geom_text_repel(size = label_size, fontface = "bold", show.legend = F, max.overlaps = 1000,
                      segment.size = 0.2, min.segment.length = 0) +
      scale_color_manual(values = c("grey70", ggsci::pal_npg()(5)[c(2,1,5)]) %>% set_names(1:4)) +
      scale_size_manual(values = point_sizes %>% set_names(c(F, T))) +
      labs(x = meanFC_label,
           y = logFC_label) +
      expand_limits(x = 1) +
      scale_x_continuous(trans = "log2") +
      theme_cowplot()
    ggsave(paste0(plot_path, "_limma-MA-FC.pdf"), width = 6, height = 5)

    return(list(volcano = volcano_gg,
                ma = ma_gg))

  }




