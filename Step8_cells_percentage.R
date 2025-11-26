# --- Paths & Input ---
setwd("C:/Users/AayudhDas/OneDrive - StratusTX/SINGLE CELL/ST-101_Res_1")
list.files()
meta_matrix_all <- read.csv("raw_counts_expr_df_res1_res2.csv", stringsAsFactors = FALSE)
head(meta_matrix_all)

# Keep only required columns
meta_matrix_all<- meta_matrix_all[, c("cell_id","cluster_res1","cluster_res2", "CD34", "THY1", "PTPRC", "PECAM1")]

# Rename columns: THY1 → CD90, PTPRC → CD45, PECAM1 → CD31
colnames(meta_matrix_all) <- c("cell_id", "cluster_res1","cluster_res2","CD34", "CD90", "CD45", "CD31")

# Check first few rowss
head(meta_matrix_all)

cd34_zero_rows <- sum(meta_matrix_all$CD34 == 0, na.rm = TRUE)

cat("Rows with CD34 = 0:", cd34_zero_rows, "\n")


# Filter rows where CD34 > 0
meta_matrix_all_cd34 <- meta_matrix_all[
  meta_matrix_all$CD34 > 0,
]

# Filter rows where both CD34 > 0 and CD90 > 0
meta_matrix_all_cd34_cd90 <- meta_matrix_all[
  meta_matrix_all$CD34 > 0 & meta_matrix_all$CD90 > 0,
]

meta_matrix_all_CD90 <- meta_matrix_all[
  meta_matrix_all$CD90 > 0,
]

meta_matrix_all_CD45 <- meta_matrix_all[
  meta_matrix_all$CD45 > 0,
]


meta_matrix_all_cd34_cd45 <- meta_matrix_all[
  meta_matrix_all$CD34 > 0 & meta_matrix_all$CD45 > 0,
]

total_rows <- nrow(meta_matrix_all)
cd34_rows <- nrow(meta_matrix_all_cd34)
cd90_rows <- nrow(meta_matrix_all_CD90)
cd45_rows <- nrow(meta_matrix_all_CD45)
cd34_cd90_rows <- nrow(meta_matrix_all_cd34_cd90)
cd34_cd45_rows <- nrow(meta_matrix_all_cd34_cd45)

cat("Percentage of CD34+ cells:", round((cd34_rows / total_rows) * 100, 2), "%\n")
cat("Percentage of CD34+CD90+ cells:", round((cd34_cd90_rows / total_rows) * 100, 2), "%\n")
cat("Percentage of CD90+ cells:", round((cd90_rows / total_rows) * 100, 2), "%\n")
cat("Percentage of CD45+ cells:", round((cd45_rows / total_rows) * 100, 2), "%\n")
cat("Percentage of CD34+CD45+ cells:", round((cd34_cd45_rows / total_rows) * 100, 2), "%\n")

#------------------ENG1--------------------------------------------------------------

meta_matrix_eng1 <- meta_matrix_all[grepl("ENG1_ENG1", meta_matrix_all$cell_id), ]

# Check the first few rows
head(meta_matrix_eng1)

meta_matrix_eng1_cd34 <- meta_matrix_eng1[
  meta_matrix_eng1$CD34 > 0,
]

# --- CD34+CD90+ filtering on ENG1 dataset ---
meta_matrix_eng1_cd34_cd90 <- meta_matrix_eng1[
  meta_matrix_eng1$CD34 > 0 & meta_matrix_eng1$CD90 > 0,
]

# --- Row counts ---
total_rows_eng1 <- nrow(meta_matrix_eng1)
cd34_rows_eng1 <- nrow(meta_matrix_eng1_cd34)
cd34_cd90_rows_eng1 <- nrow(meta_matrix_eng1_cd34_cd90)

# --- Percentages ---
cat("ENG1 dataset total rows:", total_rows_eng1, "\n")
cat("ENG1 CD34+ cells:", cd34_rows_eng1, 
    "(", round((cd34_rows_eng1 / total_rows_eng1) * 100, 2), "% )\n")
cat("ENG1 CD34+CD90+ cells:", cd34_cd90_rows_eng1, 
    "(", round((cd34_cd90_rows_eng1 / total_rows_eng1) * 100, 2), "% )\n")

####---------------gmp1------------------
meta_matrix_gmp1 <- meta_matrix_all[grepl("ENG1_GMP1", meta_matrix_all$cell_id), ]

# Check the first few rows
head(meta_matrix_gmp1)

meta_matrix_gmp1_cd34 <- meta_matrix_gmp1[
  meta_matrix_gmp1$CD34 > 0,
]

# --- CD34+CD90+ filtering on ENG1 dataset ---
meta_matrix_gmp1_cd34_cd90 <- meta_matrix_gmp1[
  meta_matrix_gmp1$CD34 > 0 & meta_matrix_gmp1$CD90 > 0,
]

# --- Row counts ---
total_rows_eng1 <- nrow(meta_matrix_gmp1)
cd34_rows_eng1 <- nrow(meta_matrix_gmp1_cd34)
cd34_cd90_rows_eng1 <- nrow(meta_matrix_gmp1_cd34_cd90)

# --- Percentages ---
cat("ENG1 dataset total rows:", total_rows_eng1, "\n")
cat("ENG1 CD34+ cells:", cd34_rows_eng1, 
    "(", round((cd34_rows_eng1 / total_rows_eng1) * 100, 2), "% )\n")
cat("ENG1 CD34+CD90+ cells:", cd34_cd90_rows_eng1, 
    "(", round((cd34_cd90_rows_eng1 / total_rows_eng1) * 100, 2), "% )\n")

############-----------PD160--------------------------------
meta_matrix_pd160 <- meta_matrix_all[grepl("PD160", meta_matrix_all$cell_id), ]

# Check the first few rows
head(meta_matrix_pd160)

meta_matrix_pd160_cd34 <- meta_matrix_pd160[
  meta_matrix_pd160$CD34 > 0,
]

# --- CD34+CD90+ filtering on ENG1 dataset ---
meta_matrix_pd160_cd34_cd90 <- meta_matrix_pd160[
  meta_matrix_pd160$CD34 > 0 & meta_matrix_pd160$CD90 > 0,
]

# --- Row counts ---
total_rows_eng1 <- nrow(meta_matrix_pd160)
cd34_rows_eng1 <- nrow(meta_matrix_pd160_cd34)
cd34_cd90_rows_eng1 <- nrow(meta_matrix_pd160_cd34_cd90)

# --- Percentages ---
cat("ENG1 dataset total rows:", total_rows_eng1, "\n")
cat("ENG1 CD34+ cells:", cd34_rows_eng1, 
    "(", round((cd34_rows_eng1 / total_rows_eng1) * 100, 2), "% )\n")
cat("ENG1 CD34+CD90+ cells:", cd34_cd90_rows_eng1, 
    "(", round((cd34_cd90_rows_eng1 / total_rows_eng1) * 100, 2), "% )\n")

######---------------table-------------
# --- Helper function ---
get_stats <- function(df) {
  total <- nrow(df)
  cd34_n <- sum(df$CD34 > 0, na.rm = TRUE)
  cd34_cd90_n <- sum(df$CD34 > 0 & df$CD90 > 0, na.rm = TRUE)
  data.frame(
    `Total Cells` = total,
    `CD34+` = sprintf("%d (%.2f%%)", cd34_n, 100 * cd34_n / total),
    `CD34+CD90+` = sprintf("%d (%.2f%%)", cd34_cd90_n, 100 * cd34_cd90_n / total),
    check.names = FALSE
  )
}

# Subsets
meta_matrix_eng1  <- meta_matrix_all[grepl("ENG1_ENG1",  meta_matrix_all$cell_id), ]
meta_matrix_gmp1  <- meta_matrix_all[grepl("ENG1_GMP1",  meta_matrix_all$cell_id), ]
meta_matrix_pd160 <- meta_matrix_all[grepl("PD160", meta_matrix_all$cell_id), ]

# Combined (ENG1 + GMP1 + PD160)
meta_matrix_combined <- rbind(meta_matrix_eng1, meta_matrix_gmp1, meta_matrix_pd160)

# Build table
tbl <- rbind(
  `PD160+ENG1+GMP1` = get_stats(meta_matrix_combined),
  `ENG1`            = get_stats(meta_matrix_eng1),
  `GMP1`            = get_stats(meta_matrix_gmp1),
  `PD-160`          = get_stats(meta_matrix_pd160)
)

# Print
print(tbl)
write.csv(tbl,"percentage_cells.csv")

