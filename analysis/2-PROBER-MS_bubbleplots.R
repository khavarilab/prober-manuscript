library(tidyverse)
library(scales)
library(readxl)
library(here)
library(cowplot)


out_dir <- here(file.path("results/bubbleplots", Sys.Date()))
dir.create(out_dir, recursive = T)


bubbleplot_data <- read_excel(here("data/bubbleplot_dec20.xlsx"))


bubbleplot_data %>%
  mutate(Protein = factor(Protein, levels = Protein)) %>%
  pivot_longer(-Protein, names_sep = "-", names_to = c("stat", "bait"), values_to = "value") %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  filter(FC > 1) %>%
  mutate(bait = factor(bait, levels = c("YY1", "NFkB", "STAT1"),
                       labels = c("YY1", "NF-kB (+TNFa)", "STAT1 (+IFNy)"))) %>%
  ggplot(aes(x = fct_rev(Protein), y = fct_rev(bait))) +
  geom_point(aes(color = FC, size = SAINT)) +
  scale_size_continuous(range = c(0, 4), guide = guide_legend(reverse = T),
                        trans = exp_trans(base = 10)) +
  expand_limits(color = 1, size = 0) +
  scale_x_discrete(position = "top") +
  scale_color_distiller(palette = "YlOrRd", direction = 1, trans = "log2",
                        breaks = c(2, 4, 8, 16, 32)) +
  theme(panel.background = element_rect(color = "black", fill = "white"),
        panel.grid = element_line(size = 0.2, color = "grey80"),
        axis.text.x = element_text(angle = 45, hjust = 0),
        axis.title = element_blank(), legend.box = "horizontal")

ggsave(file.path(out_dir, "bubbleplot.pdf"), width = 10, height = 2.5)
