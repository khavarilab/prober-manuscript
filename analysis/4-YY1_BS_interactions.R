
library(tidyverse)
library(readxl)
library(tidygraph)
library(here)
library(scales)
library(cowplot)
library(ggrepel)
library(ggraph)
library(ggforce)

out_dir <- file.path("results/string_interactions", Sys.Date())
dir.create(out_dir, recursive = T)

rev_data <- list.files(here("data/rev_MS"), pattern = "^[A-Za-z0-9].+\\.xlsx", full.names = T) %>%
  set_names(., str_replace(basename(.), "\\.xlsx$", "")) %>% map(read_excel)


bs_data <-
  bind_rows(rev_data$`26-9_SAINT` %>%
              filter(BAIT_SP >= 0.9) %>%
              select(PROTID) %>% mutate(site = "BS1"),
            rev_data$`26-11-SAINT` %>%
              filter(BAIT_SP >= 0.9) %>%
              select(PROTID) %>% mutate(site = "BS2"),
            rev_data$`26-14-SAINT` %>%
              filter(BAIT_SP >= 0.9) %>%
              select(PROTID) %>% mutate(site = "BS3"))


bs_node_df <- bs_data %>% mutate(dummy = T) %>%
  pivot_wider(names_from = "site", values_from = "dummy", values_fill = F)

remap_names <- c("GLTSCR1" = "BICRA",
                 "QARS" = "QARS1",
                 "WHSC1L1" = "NSD3")

string_interactions <- read_tsv(here("data/string/YY1_BS_string_interactions.tsv")) %>%
  rename(PROTID.x = `#node1`, PROTID.y = node2) %>%
  mutate(PROTID.x = ifelse(PROTID.x %in% names(remap_names), remap_names[PROTID.x], PROTID.x),
         PROTID.y = ifelse(PROTID.y %in% names(remap_names), remap_names[PROTID.y], PROTID.y)) %>%
  mutate(tmp.x = PROTID.x,
         tmp.y = PROTID.y,
         PROTID.x = ifelse(tmp.x < tmp.y, tmp.x, tmp.y),
         PROTID.y = ifelse(tmp.x < tmp.y, tmp.y, tmp.x)) %>%
  filter(experimentally_determined_interaction > 0)


interaction_sums <- bind_rows(string_interactions %>% select(PROTID_1 = PROTID.x, PROTID_2 = PROTID.y, experimentally_determined_interaction),
                              string_interactions %>% select(PROTID_1 = PROTID.y, PROTID_2 = PROTID.x, experimentally_determined_interaction)) %>%
  group_by(PROTID_1) %>%
  summarise(sum = sum(experimentally_determined_interaction))

string_interactions_scaled <- string_interactions %>%
  left_join(interaction_sums %>% rename(sum.x = sum), by = c("PROTID.x" = "PROTID_1")) %>%
  left_join(interaction_sums %>% rename(sum.y = sum), by = c("PROTID.y" = "PROTID_1")) %>%
  mutate(scaled_interaction = experimentally_determined_interaction / max(c(sum.x, sum.y)))

# View(string_interactions_scaled %>% select(PROTID.x, PROTID.y, experimentally_determined_interaction, sum.x, sum.y, scaled_interaction))


bs_interaction_data <- bs_data %>%
  group_by(site) %>% mutate(n = n()) %>% ungroup %>%
  inner_join(. , ., by = c("site", "n")) %>%
  complete(PROTID.x, PROTID.y) %>%
  # mutate(shared_binding_site = ifelse(site.x == site.y, 1, 0)) %>%
  filter(PROTID.x < PROTID.y) %>%
  # select(PROTID.x, PROTID.y, site) %>%
  # distinct() %>%
  full_join(string_interactions_scaled) %>%
  mutate(across(c(experimentally_determined_interaction, scaled_interaction), .fns = ~ ifelse(is.na(.), 0, .))) %>%
  distinct(PROTID.x, PROTID.y, site, experimentally_determined_interaction, scaled_interaction)


ks.test(bs_interaction_data %>% filter(is.na(site)) %>% pull(experimentally_determined_interaction),
        bs_interaction_data %>% filter(site == "BS1") %>% pull(experimentally_determined_interaction),
        alternative = "greater") %>%
  print()
ks.test(bs_interaction_data %>% filter(is.na(site)) %>% pull(experimentally_determined_interaction),
        bs_interaction_data %>% filter(site == "BS2") %>% pull(experimentally_determined_interaction),
        alternative = "greater") %>%
  print()
ks.test(bs_interaction_data %>% filter(is.na(site)) %>% pull(experimentally_determined_interaction),
        bs_interaction_data %>% filter(site == "BS3") %>% pull(experimentally_determined_interaction),
        alternative = "greater") %>%
  print()

colorpal <- brewer_pal(palette = "Set1")(3)

bs_interaction_data %>%
  ggplot(aes(x = experimentally_determined_interaction, color = site)) +
  stat_ecdf() +
  coord_cartesian(ylim = c(0.8, 1)) +
  scale_color_manual(values = colorpal, na.value = "grey50") +
  labs(x = "STRING interaction score", y = "Cumulative fraction") +
  theme_cowplot() +
  theme(legend.position = c(1, 0), legend.justification = c(1, 0))
ggsave(file.path(out_dir, "string_ecdf.pdf"), width = 4, height = 3)


bs_edge_df <- bs_interaction_data %>%
  group_by(PROTID.x, PROTID.y, experimentally_determined_interaction, scaled_interaction) %>%
  summarise(shared_binding_site = n_distinct(site, na.rm = T)) %>%
  ungroup() %>%
  mutate(edge_weight = shared_binding_site + 10*scaled_interaction) %>%
  select(PROTID.x, PROTID.y, experimentally_determined_interaction, scaled_interaction, edge_weight) %>%
  filter(edge_weight > 0, PROTID.x %in% bs_node_df$PROTID, PROTID.y %in% bs_node_df$PROTID)




bs_graph <- tbl_graph(nodes = bs_node_df, edges = bs_edge_df, node_key = "PROTID", directed = F)

initial_layout <- create_layout(bs_graph %>%
                                  activate(edges) %>%
                                  filter(scaled_interaction > 0),
                                layout = "fr")#, weights = .E()$edge_weight)


ggraph(bs_graph, layout = "fr", coords = as.matrix(initial_layout[,c("x", "y")]),
       weights = .E()$edge_weight, start.temp = 0.5, niter = 1000) +
  geom_edge_link(aes(edge_width = .E()$experimentally_determined_interaction,
                     filter = .E()$experimentally_determined_interaction > 0),
                 color = "grey70") +
  geom_node_point(size = 1) +
  geom_mark_circle(aes(x, y, filter = BS1), fill = colorpal[1], color = colorpal[1], alpha = 0.15) +
  geom_mark_circle(aes(x, y, filter = BS2), fill = colorpal[2], color = colorpal[2], alpha = 0.15) +
  geom_mark_circle(aes(x, y, filter = BS3), fill = colorpal[3], color = colorpal[3], alpha = 0.15) +
  geom_text_repel(aes(x, y, label = .N()$PROTID), segment.colour = NA, size = 4, fontface = "bold",
                  box.padding = unit(0.15, "lines"), max.overlaps = Inf) +
  scale_edge_width_continuous(range = c(0.1, 1), name = "STRING") +
  scale_x_continuous(expand = expansion(mult = .15)) +
  scale_y_continuous(expand = expansion(mult = .15)) +
  coord_equal() +
  theme_nothing() +
  theme(legend.position = "bottom")
ggsave(file.path(out_dir, "YY1_binding_site_network.pdf"), width = 8, height = 8)

