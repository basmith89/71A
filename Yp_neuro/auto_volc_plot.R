library(tidyverse)
library(scales)
library(ggrepel)
library(patchwork)

data_dir <- "Your_data_dir"
setwd("wd")

#csv file intake
# note rev() I use this to output in the correct order when display multiple figs in one.
# this would be dependent on what filenames are, but was a quick hack in my specific case.
csv_files <- rev(list.files(
  path = data_dir,
  pattern = "D3\\.csv$", #### edit this line to filter files
  full.names = TRUE
))

cols <- c("Up" = "red", "Down" = "#26b3ff", "Ns" = "grey")
sizes <- c("Up" = 2, "Down" = 2, "ns" = 1) 
alphas <- c("Up" = 1, "Down" = 1, "ns" = 0.5)

#store filtering as a function
process_df <- function(df) {
  df %>%
    mutate(
      gene_type = case_when(
        Log2.fold.change >= 1 & BY.p.value <= 0.05  ~ "Up",
        Log2.fold.change <= -1 & BY.p.value <= 0.05 ~ "Down",
        TRUE                                       ~ "Ns"
      )
    )
}

#store plotting as a function
plot_volcano <- function(df, sig_genes, up_genes, down_genes, title = NULL) {
  
  sig_interest_genes <- df %>%
    filter(
      Gene_ID %in% sig_genes,
      abs(Log2.fold.change) >= 1,
      BY.p.value <= 0.05
    )
  
  up_int_genes <- df %>%
    filter(
      Gene_ID %in% up_genes,
      abs(Log2.fold.change) >= 1,
      BY.p.value <= 0.05
    )
  
  down_int_genes <- df %>%
    filter(
      Gene_ID %in% down_genes,
      abs(Log2.fold.change) >= 1,
      BY.p.value <= 0.05
    )
  
  ggplot(df, aes(Log2.fold.change, -log10(P.value))) +
    geom_point(aes(colour = gene_type),
               alpha = 0.7,
               shape = 16,
               size = 1) +
    geom_point(data = up_int_genes,
               shape = 21,
               size = 2,
               fill = "firebrick",
               colour = "black") +
    geom_point(data = down_int_genes,
               shape = 21,
               size = 2,
               fill = "steelblue",
               colour = "black") +
    geom_hline(yintercept = -log10(0.05),
               linetype = "dashed") +
    geom_vline(xintercept = c(log2(0.5), log2(2)),
               linetype = "dashed") +
    geom_label_repel(
      data = sig_interest_genes,
      aes(label = Gene_ID),
      force = 2,
      nudge_y = 1
    ) +
    scale_colour_manual(values = cols) +
    scale_x_continuous(breaks = seq(-11, 11.1, 2),
                       limits = c(-11, 11.1)) +
    labs(
      color = "Diff Exp",
      title = title
    )
}

#lists of interested genes for plotting
astro_genes <- c("Lcn2", "Agt", "Fbln5", "Ptgs2", "Gbp2",
                 "Hspb1", "Serpina3n", "Serping1")

micro_sig <- c("Gpr34", "Irf1", "Srgn", "P2ry12",
               "Cx3cr1", "C1qc", "C3", "C1qa", "C1qb")

micro_up   <- c("Irf1", "Srgn", "C1qc", "C3", "C1qa", "C1qb")
micro_down <- c("Gpr34", "P2ry12", "Cx3cr1")

oligo_genes <- c("Opalin", "Sox10", "Ugt8a", "Pllp")


#run the loop on all files with .csv in directory
#this will run and save all files individually
for (file in csv_files) {
  
  df <- read.csv(file, header = TRUE) %>%
    process_df()
  
  base_name <- tools::file_path_sans_ext(basename(file))
  
  p1 <- plot_volcano(
    df,
    sig_genes = astro_genes,
    up_genes  = astro_genes,
    down_genes = character(0),
    title = paste(base_name, "- Astrocytes")
  )
  
  p2 <- plot_volcano(
    df,
    sig_genes = micro_sig,
    up_genes  = micro_up,
    down_genes = micro_down,
    title = paste(base_name, "- Microglia")
  )
  
  p3 <- plot_volcano(
    df,
    sig_genes = oligo_genes,
    up_genes  = character(0),
    down_genes = oligo_genes,
    title = paste(base_name, "- Oligodendrocytes")
  )
#}  
  #print(p1)
  #print(p2)
  #print(p3)
  
  ## Optional: save
    ggsave(paste0(base_name, "_astro.png"), p1, width = 6, height = 5)
    ggsave(paste0(base_name, "_micro.png"), p2, width = 6, height = 5)
    ggsave(paste0(base_name, "_oligo.png"), p3, width = 6, height = 5)
}


#################################################
###Same thing but now combine the into one fig###
#################################################

#counter necessary for labeling with the apropriate letter labels in automation
counter <- 0

for (file in csv_files) {
  
  df <- read.csv(file, header = TRUE) %>%
    process_df()
  
  base_name <- tools::file_path_sans_ext(basename(file))
  
  p1 <- plot_volcano(
    df,
    sig_genes  = astro_genes,
    up_genes   = astro_genes,
    down_genes = character(0),
    title = paste(base_name, "- Astrocytes")
  )
  
  p2 <- plot_volcano(
    df,
    sig_genes  = micro_sig,
    up_genes   = micro_up,
    down_genes = micro_down,
    title = paste(base_name, "- Microglia")
  )
  
  p3 <- plot_volcano(
    df,
    sig_genes  = oligo_genes,
    up_genes   = character(0),
    down_genes = oligo_genes,
    title = paste(base_name, "- Oligodendrocytes")
  )
  
  p1 <- p1 + labs(title = NULL)
  p2 <- p2 + labs(title = NULL)
  p3 <- p3 + labs(title = NULL)
  
  # Remove all group legends (color key)
  p1 <- p1 + theme(legend.position = "none")
  p2 <- p2 + theme(legend.position = "none")
  p3 <- p3 + theme(legend.position = "none")
  
  # Keep only y-axis title on p1, remove from p2 and p3
  p1 <- p1 + labs(y = expression(-Log[10](italic(P))))
  p2 <- p2 + labs(y = NULL)
  p3 <- p3 + labs(y = NULL)
  
  # Keep only x-axis title on p2, remove from p1 and p3
  p1 <- p1 + labs(x = NULL)
  p2 <- p2 + labs(x = "Log2 Fold Change")
  p3 <- p3 + labs(x = NULL)
  
  #counter code, 3 plots per file
  tags <- LETTERS[counter * 3 + 1:3]
  counter <- counter + 1
  
  # Combine
  combined <- p1 + p2 + p3 +
    plot_layout(nrow = 1) +
    plot_annotation(tag_levels = list(tags)) #labels each figure plot in alph order starting with "X" letter.
  
  
  # Combine into a single grouped figure and save
 # combined <- p1 + p2 + p3 +
  #  plot_layout(nrow = 1, guides = "collect") +
   # plot_annotation(tag_levels = "A")
  
  ggsave(
    paste0(base_name, "_combined.png"),
    combined,
    width = 18,
    height = 5
  )
}
