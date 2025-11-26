setwd("C:/Users/AayudhDas/OneDrive - StratusTX/SINGLE CELL/ST-101_Res_1")
list.files()
library(ggplot2)

gene_data <- read.csv("cluster9_14.csv")
head(gene_data)

gene_data$Gene <- ifelse(
  is.na(gene_data$Function) | trimws(gene_data$Function) == "",
  gene_data$Gene,
  paste0(gene_data$Gene, " | ", trimws(gene_data$Function))
)

# Order genes from high to low (GJA5 at top)
gene_data$Gene <- factor(gene_data$Gene, 
                         levels = rev(gene_data$Gene[order(gene_data$avg_log2FC, decreasing = TRUE)]))

# Count negative genes to add dashed line
negative_genes <- sum(gene_data$avg_log2FC < 0)
dotted_line_position <- negative_genes + 0.5

# Save PNG
png("gene_dotplot_single_line_right_y_axis.png", width = 13, height = 10, units = 'in', res = 300)

ggplot(gene_data, aes(x = 0, y = as.numeric(Gene))) +
  geom_hline(yintercept = dotted_line_position, linetype = "dashed", color = "black", size = 0.8) +
  geom_point(aes(size = abs(avg_log2FC), color = avg_log2FC > 0), alpha = 0.8) +
  scale_y_continuous(
    breaks = 1:length(levels(gene_data$Gene)),
    labels = levels(gene_data$Gene),
    sec.axis = dup_axis() # <- Move labels to right-hand side
  ) +
  scale_size_continuous(range = c(1, 8), name = "log2FC", 
                        breaks = c(2, 4, 6, 8), 
                        labels = c("2", "4", "6", "8")) +
  scale_color_manual(values = c("FALSE" = "#717C7D", "TRUE" = "#D50032"), guide = "none") +
  theme_classic() +
  theme(
    text = element_text(family = "Arial Black", size = 14),
    axis.text.y.right = element_text(color = "black", face = "italic", size = 11), # Right-side labels
    axis.text.y.left = element_blank(),  # Hide left labels
    axis.ticks.y.left = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(color = "black", face = "bold", size = 14),
    axis.line.x = element_line(color = "black", size = 0.8),
    axis.line.y.right = element_line(color = "black", size = 0.8),  # Right y-axis line
    axis.line.y.left = element_blank(),  # Hide left y-axis line
    legend.position = "left",
    legend.title = element_text(size = 12, face = "italic"),
    legend.text = element_text(size = 12),
    plot.margin = margin(20, 40, 20, 20)
  ) +
  labs(x = "log2FC", y = NULL, title = NULL)

dev.off()
