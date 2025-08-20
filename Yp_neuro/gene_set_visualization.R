library(tidyr)
library(dplyr)
library(ggplot2)


df <- read.csv("~/Documents/USAMRIID/code/Yp_neuro_code/Day_1/gene_set_summary_counts.txt", header =TRUE, sep = "\t")

#pivot data
df_long <- df %>% 
  pivot_longer(
    cols = !Gene.set, 
    names_to = "DEG_Group", 
    values_to = "Count"
  )

#lollipop
ggplot(df_long)+
  geom_linerange(aes(x = reorder(Gene.set, Count), ymin = 0, ymax = Count, colour = DEG_Group), 
                 size = 1, position = position_dodge(width = 0.5)) +
  geom_point(aes(x = reorder(Gene.set, Count), y = Count, colour = DEG_Group),
             size = 2, position = position_dodge(width = 0.5)) +
  coord_flip() +
  xlab("Gene Set Pathways")+
  ylab("Differentially Expressed Pathway Counts Day 1")+
  theme(axis.title = element_text(face = "bold"), legend.position = c(0.85,0.15))

#stacked barplot
ggplot(df_long, aes(fill= DEG_Group, y=Count, x=Gene.set)) +
         geom_bar(position = "stack", stat = "identity")