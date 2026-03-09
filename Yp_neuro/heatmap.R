library(ComplexHeatmap)
library(dplyr)
library(circlize)

heat_df_yp <- read.csv("heatmap_scores_file_bug1.csv", header = TRUE, sep = ",")
heat_df_bp <- read.csv("heatmap_scores_filebug2.csv", header = TRUE, sep = ",")

#filter the undirected data
heat_df_yp <- heat_df_yp %>% select(-c(Undirected.DPI..differential.expression.in.1.vs..baseline.of.CTRL, Undirected.DPI..differential.expression.in.2.vs..baseline.of.CTRL, Undirected.DPI..differential.expression.in.3.vs..baseline.of.CTRL))

rownames(heat_df_yp) <- heat_df_yp[, 1] #make first col rownames
heat_df_yp <- heat_df_yp[, -1] #delete the old first column

heat_df_yp <- heat_df_yp %>%
  rename("1DPC vs. Control" = Directed.DPI..differential.expression.in.1.vs..baseline.of.CTRL, "2DPC vs. Control" = Directed.DPI..differential.expression.in.2.vs..baseline.of.CTRL, "3DPC vs. Control" = Directed.DPI..differential.expression.in.3.vs..baseline.of.CTRL)

col_range <- colorRamp2(c(-5, 0, 5), c("blue", "white", "red"))

heatmap_yp <- Heatmap(heat_df_yp, name = "Directed Global Significance Score",
        col = col_range,
        column_order = order(as.numeric(gsub("column", "", colnames(heat_df_yp)))),
        column_title = "Y.pestis CO92 - BALB/c")



#####Repeat for Bp
#filter the undirected data
heat_df_bp <- heat_df_bp %>% select(-c(Undirected.DPI..differential.expression.in.1.vs..baseline.of.CTRL, Undirected.DPI..differential.expression.in.2.vs..baseline.of.CTRL, Undirected.DPI..differential.expression.in.3.vs..baseline.of.CTRL))

rownames(heat_df_bp) <- heat_df_bp[, 1] #make first col rownames
heat_df_bp <- heat_df_bp[, -1] #delete the old first column

heat_df_bp <- heat_df_bp %>%
  rename("1DPC vs. Control" = Directed.DPI..differential.expression.in.1.vs..baseline.of.CTRL, "2DPC vs. Control" = Directed.DPI..differential.expression.in.2.vs..baseline.of.CTRL, "3DPC vs. Control" = Directed.DPI..differential.expression.in.3.vs..baseline.of.CTRL)

col_range <- colorRamp2(c(-10, 0, 10), c("blue", "white", "red"))
heatmap_bp <- Heatmap(heat_df_bp, name = "Directed Global Significance Score",
        col = col_range,
        column_order = order(as.numeric(gsub("column", "", colnames(heat_df_bp)))),
        column_title = "B. pseudomallei ATS2021 C57BL/6")


png(file="heatmap_yp.png", width = 85 * (14/5),height=53 * (14/5),units="mm",res=300)
draw(heatmap_yp)
png(file="heatmap_bp.png", width = 85 * (14/5),height=53 * (14/5),units="mm",res=300)
draw(heatmap_bp)
dev.off
