library(tidyverse)
library(ggplot2)
library(ggsignif)

# Read the CSV
df <- read.csv("your_raw_data_file.csv",
               stringsAsFactors = FALSE,
               check.names = FALSE)

annotation_df <- read.csv("stats_plot_df.csv",
                          stringsAsFactors = FALSE,
                          check.names = FALSE)

# Rename first column to something explicit
colnames(df)[1] <- "Cytokine"

df[df == ""] <- NA

long_df <- df %>%
  pivot_longer(
    cols = -Cytokine,
    names_to = "Day",
    values_to = "Value"
  ) %>%
  mutate(
    Day = str_extract(Day, "^[^\\.]+")  # removes replicate suffixes
  )


####Plot a single cytokine like IL-6
il6_df <- long_df %>%
  filter(Cytokine == "IL-6") %>%
  mutate(
    Day = factor(
      Day,
      levels = c("Control", setdiff(unique(Day), "Control"))
    )
  ) 

ggplot(il6_df, aes(x = Day, y = Value)) +
  geom_boxplot(outlier.shape = 16) +
  labs(
    title = "IL-6 Levels Across Days",
    x = "Day Post Challenge",
    y = "IL-6 Expression"
  ) +
  geom_jitter(
    aes(color = Day),           # Color points
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), 
    size = 2,                     # Size of the points for visibility
    alpha = 0.6                   # Slight transparency to avoid clutter
  ) +
  theme_minimal() +
  geom_signif(comparisons = list(c("Control", "1DPC"),
                                 c("Control", "2DPC"),
                                  c("Control", "3DPC")),
              map_signif_level = TRUE,
              step_increase = .08)
              #y_position = c(2,2.1,2.2))


#filter cytokines of interest
filtered_df <- long_df %>%
  filter(
    Cytokine %in% c(
      "IL-6", "IL-1 alpha", "IL-27", "MIP-2 alpha (CXCL2)",
      "G-CSF", "TNF alpha", "MCP-1 (CCL2)",
      "GRO-alpha (CXCL1)", "RANTES (CCL5)"
    )
  ) %>%
  mutate(
    Day = factor(
      Day,
      levels = c("Control", setdiff(unique(Day), "Control"))
    )
  )

#plot the filtered cytokines of interest in a facet wrapped group
cyto_fil_plot <- ggplot(filtered_df, aes(x = Day, y = Value)) +
  geom_boxplot(outlier.shape = 16) +
  #geom_jitter(
  #  aes(color = Day),
  #  position = position_jitterdodge(
  #    jitter.width = 0.2,
  #    dodge.width = 0.8
  #  ),
  #  size = 2,
  #  alpha = 0.6
  #) +
  geom_signif(
    data = annotation_df,
    aes(xmin = start, xmax = end, annotations = label, y_position = y),
    textsize = 4, vjust = -0.1,
    manual = TRUE) +
  facet_wrap(~ Cytokine, scales = "free_y") +
  scale_y_continuous(
    expand = expansion(mult = c(0.00, 0.20))
  ) +
  labs(
    x = "Day Post Challenge",
    y = "Log10[pg/µL]"
  ) +
  theme_bw() 

#plot to specific mm for pub
ggsave(
  filename = "outfile.png",
  plot = cyto_fil_plot,
  width = 85 * (14/5),
  height = 64 * (14/5),
  units = "mm",
  dpi = 300
)

### plot all cytokines
ggplot(long_df, aes(x = Day, y = Value)) +
  geom_boxplot(outlier.shape = 16) +
  facet_wrap(~ Cytokine, scales = "free_y") +
  theme_minimal()
