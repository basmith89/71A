library(tidyverse)
library(ggplot2)
library(ggh4x)

df <- read.csv("GSA_scores__for_histo.csv", check.names = FALSE)

long_df <- df %>%
  pivot_longer(
    cols = -Pathway,
    names_to = "Condition",
    values_to = "Score"
  ) %>%
  separate(Condition, into = c("DPC", "Group"), sep = "-") %>%
  filter(
    Pathway %in% c(
      "Astrocyte Function",
      "Oligodendrocyte Function",
      "Microglia Function",
      "Neurons and Neurotransmission"
    )
  ) %>%
  mutate(
    DPC = factor(DPC, levels = c("1DPC", "2DPC", "3DPC")),
    Group = factor(Group, levels = c("Bp", "Yp"))
  )

gsa_comp_plot <- ggplot(long_df, aes(x = DPC, y = Score, color = Group)) +
  geom_point(
    position = position_dodge(width = 0.4),
    size = 3
  ) +
  facet_wrap(~ Pathway, scales = "free_y") +
  labs(
    title = "Neural and Glial Pathway Scores by DPC and Group",
    y = "Pathway Score"
  ) +
  theme_grey()


#new
#
#
# Build per-pathway color scales
pathway_scales <- long_df %>%
  group_by(Pathway) %>%
  summarise(
    breaks = list(sort(pretty(Score))),
    .groups = "drop"
  ) %>%
  mutate(
    colors = map(breaks, ~ ifelse(.x > 0, "red", ifelse(.x < 0, "blue", "black")))
  )

scale_list <- map2(
  pathway_scales$breaks,
  pathway_scales$colors,
  ~ scale_y_continuous(
    breaks = .x,
    limits = range(.x),  # force limits to match only actual breaks
    guide = guide_axis(theme = theme(axis.text.y = element_text(color = .y)))
  )
)
names(scale_list) <- pathway_scales$Pathway

gsa_comp_plot <- ggplot(long_df, aes(x = DPC, y = Score, color = Group)) +
  geom_point(
    position = position_dodge(width = 0.4),
    size = 3
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  facet_wrap(~ Pathway, scales = "free_y") +
  facetted_pos_scales(y = scale_list) +  # ggh4x: per-facet y scales
  labs(
    title = "",
    y = "Directed Global Significance Score"
  ) +
  theme_grey()



#plot to specific mm for pub
ggsave(
  filename = "GSA_comp_plot.png",
  plot = gsa_comp_plot,
  width = 85 * (14/5),
  height = 53 * (14/5),
  units = "mm",
  dpi = 300
)
