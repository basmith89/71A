library(tidyverse)
library(scales)
library(ggrepel)

df <- read.csv("PATH", header =TRUE, sep = ",")

df <- df %>%
  mutate(gene_type = case_when(Log2.fold.change >= 1 & BY.p.value <= 0.05 ~ "up",
                               Log2.fold.change <= -1 & BY.p.value <= 0.05 ~ "down",
                               TRUE ~ "ns"))   


###Label vol plot w/ genes of interest
sig_interest_genes <- df %>%
  filter(Gene_ID %in% c("Opalin", "Sox10", "Ugt8a", "Pllp", "C1qc", "C3", "Lcn2", "Agt"),
         (Log2.fold.change >= 1 | Log2.fold.change <= -1),
         BY.p.value <= 0.05
  )

up_int_genes <- df %>%
  filter(Gene_ID %in% c("C1qc", "C3", "Lcn2", "Agt"),
         (Log2.fold.change >= 1 | Log2.fold.change <= -1),
         BY.p.value <= 0.05
  )

down_int_genes <- df %>%
  filter(Gene_ID %in% c("Opalin", "Sox10", "Ugt8a", "Pllp"),
         (Log2.fold.change >= 1 | Log2.fold.change <= -1),
         BY.p.value <= 0.05
  )

ggplot(data = df,
       aes(x = Log2.fold.change,
           y = -log10(P.value))) + 
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
  labs(color = "Diff Exp") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed") +
  geom_label_repel(data = sig_interest_genes, # Add labels last to appear as the top layer  
                   aes(label = Gene_ID),
                   force = 2,
                   nudge_y = 1) +
  scale_colour_manual(values = cols) + 
  scale_x_continuous(breaks = c(seq(-11, 11.1, 2)),     
                     limits = c(-11, 11.1))  




#####Copied above for microglia genes########
sig_interest_genes <- df %>%
  filter(Gene_ID %in% c("Gpr34", "Irf1", "Srgn", "P2ry12", "Cx3cr1"),
         (Log2.fold.change >= 1 | Log2.fold.change <= -1),
         BY.p.value <= 0.05
  )

up_int_genes <- df %>%
  filter(Gene_ID %in% c("Irf1", "Srgn"),
         (Log2.fold.change >= 1 | Log2.fold.change <= -1),
         BY.p.value <= 0.05
  )

down_int_genes <- df %>%
  filter(Gene_ID %in% c("Gpr34", "P2ry12", "Cx3cr1"),
         (Log2.fold.change >= 1 | Log2.fold.change <= -1),
         BY.p.value <= 0.05
  )

ggplot(data = df,
       aes(x = Log2.fold.change,
           y = -log10(P.value))) + 
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
  labs(color = "Diff Exp") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed") +
  geom_label_repel(data = sig_interest_genes, # Add labels last to appear as the top layer  
                   aes(label = Gene_ID),
                   force = 2,
                   nudge_y = 1) +
  scale_colour_manual(values = cols) + 
  scale_x_continuous(breaks = c(seq(-11, 11, 2)),     
                     limits = c(-11, 11))  







####DEGs color coded by p adj value
cols <- c("up" = "red", "down" = "#26b3ff", "ns" = "grey") 
sizes <- c("up" = 2, "down" = 2, "ns" = 1) 
alphas <- c("up" = 1, "down" = 1, "ns" = 0.5)

df %>%
  ggplot(aes(x = Log2.fold.change,
             y = -log10(P.value),
             fill = gene_type,    
             size = gene_type,
             alpha = gene_type)) + 
  geom_point(shape = 21, # Specify shape and color as fixed local parameters    
             colour = "black") + 
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed") +
  scale_fill_manual(values = cols) + # Modify point color
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  scale_x_continuous(breaks = c(seq(-11, 11, 2)),       
                     limits = c(-11, 11))  

#simple vol plot
vol_plot <- df %>%
  ggplot(aes(x = Log2.fold.change,
         y = -log10(BY.p.value))) +
  geom_point()

vol_plot +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1),
             linetype = "dashed")

  df_filtered <- df %>% filter(BY.p.value != 1) 