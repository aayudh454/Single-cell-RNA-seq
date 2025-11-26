setwd("C:/Users/AayudhDas/OneDrive - StratusTX/SINGLE CELL/ST-101_Res_1")
library(alphahull)
library(dplyr)
library(ggplot2)
library(Seurat)
library(cowplot)
library(tictoc)
library(svglite)
library(dplyr)
library(reshape2)

load("step3_cluster_ST-101_res1.RData")

theme_custom <- theme_classic() +
  theme(
    text = element_text(family = "Arial", face = "bold", colour = "black", size = 10),
    plot.title = element_text(hjust = 0.5, face = "italic", family = "Arial", colour = "black", size = 14),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank()
  )

umap_res1 <- DimPlot(seuset2, reduction = "umap", label = TRUE,label.size = 6,repel = TRUE) + theme_custom + NoLegend() + ggtitle("Res: 1.0; Clusters = 19")
umap_res1

umap_df <- Embeddings(seuset2, "umap") %>% as.data.frame()
colnames(umap_df)[1:2] <- c("x","y")
df_umap <- umap_df %>% mutate(cluster = as.character(Idents(seuset2)))

# --- target points (16,1,13,18) ---
targets <- c("16","1","13","18")
target  <- df_umap %>% filter(cluster %in% targets)

# --- concave hull via alpha-shape ---
ah <- ashape(target$x, target$y, alpha = 2)   # smaller alpha = tighter; try 1.5–3
edges <- as.matrix(ah$edges[, c("ind1","ind2")])  # indices into ah$x
coords <- ah$x                                  # matrix [x,y] of the source points

# --- order edges into a single loop (no igraph needed) ---
order_cycle <- function(ed){
  ed <- ed[!duplicated(t(apply(ed, 1, sort))), , drop = FALSE]  # dedupe undirected edges
  path <- c(ed[1,1], ed[1,2])                                   # start with first edge
  ed   <- ed[-1, , drop = FALSE]
  # walk until we return to start or run out
  repeat {
    last <- tail(path, 1)
    i <- which(ed[,1] == last | ed[,2] == last)
    if (length(i) == 0) break
    i <- i[1]
    nextv <- if (ed[i,1] == last) ed[i,2] else ed[i,1]
    # if we’ve looped back to the start, close and stop
    if (nextv == path[1]) { path <- c(path, nextv); break }
    path <- c(path, nextv)
    ed <- ed[-i, , drop = FALSE]
    if (nrow(ed) == 0) break
  }
  unique(path)
}

idx <- order_cycle(edges)

# Fallback: if chaining failed, use convex hull to avoid messy lines
if (length(idx) < 3) idx <- chull(coords[,1], coords[,2])

hull_df <- as.data.frame(coords[idx, , drop = FALSE])
colnames(hull_df) <- c("x","y")

# --- plot clean polygon around 16/1/13/18 only ---
umap_res1 +
  geom_polygon(data = hull_df, aes(x = x, y = y),
               fill = NA, colour = "#404444", linewidth = 0.75,linetype = "dotted")

#####----------two-------------------------------------------------------------------------
# --- helper to compute polygon from given cluster ids ---
get_hull_df <- function(df_umap, clusters, alpha = 2) {
  target <- df_umap %>% filter(cluster %in% clusters)
  ah <- ashape(target$x, target$y, alpha = alpha)
  edges <- as.matrix(ah$edges[, c("ind1","ind2")])
  coords <- ah$x
  
  order_cycle <- function(ed){
    ed <- ed[!duplicated(t(apply(ed, 1, sort))), , drop = FALSE]
    path <- c(ed[1,1], ed[1,2])
    ed   <- ed[-1, , drop = FALSE]
    repeat {
      last <- tail(path, 1)
      i <- which(ed[,1] == last | ed[,2] == last)
      if (length(i) == 0) break
      i <- i[1]
      nextv <- if (ed[i,1] == last) ed[i,2] else ed[i,1]
      if (nextv == path[1]) { path <- c(path, nextv); break }
      path <- c(path, nextv)
      ed <- ed[-i, , drop = FALSE]
      if (nrow(ed) == 0) break
    }
    unique(path)
  }
  
  idx <- order_cycle(edges)
  if (length(idx) < 3) idx <- chull(coords[,1], coords[,2])
  hull_df <- as.data.frame(coords[idx, , drop = FALSE])
  colnames(hull_df) <- c("x","y")
  return(hull_df)
}

# --- polygon for clusters 16,1,13,18 ---
hull_df1 <- get_hull_df(df_umap, c("16","1","13","18"), alpha = 2)

# --- polygon for clusters 9,14 ---
hull_df2 <- get_hull_df(df_umap, c("9","14"), alpha = 2)

hull_df3 <- get_hull_df(df_umap, c("10","7","6"), alpha = 2)

# --- polygon for clusters CD34- ---
hull_df4 <- get_hull_df(df_umap, c("0","5","12","4"), alpha = 0.75)


# --- plot both polygons ---
png("Res1_polygons.png", width = 7, height = 6.2, units = "in", res = 300)
umap_res1 +
  geom_polygon(data = hull_df1, aes(x = x, y = y),
               fill = NA, colour = "#404444", linewidth = 1, linetype = "dashed") +
  geom_polygon(data = hull_df2, aes(x = x, y = y),
               fill = NA, colour = "#D50032", linewidth = 1, linetype = "solid")+
  geom_polygon(data = hull_df3, aes(x = x, y = y),
               fill = NA, colour = "black", linewidth = 1, linetype = "solid")+
  geom_polygon(data = hull_df4, aes(x = x, y = y),
              fill = NA, colour = "black", linewidth = 1, linetype = "dotdash")
dev.off()
