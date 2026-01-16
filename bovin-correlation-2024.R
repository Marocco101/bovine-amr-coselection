# ==============================================================================
# Project: Bovine AMR Co-selection Analysis 
# Date: Jan 2026
# Description: 
#   Ecological analysis of potential co-selection between Colistin resistance 
#   and alternative antibiotics (TMP-SMX, Tetracycline) in bovine E. coli 
#   using RESAPATH surveillance data.
# ==============================================================================

# 1. Load Libraries
library(dplyr)
library(ggplot2)
library(ggrepel)

# 2. Define Data Loading Function ("Surgical Extraction")
#    Note: Handles specific CSV format issues (semicolon delimiter, quotes) 
#    encountered in RESAPATH data exports.
surgical_load <- function(filepath) {
  
  message(paste("Processing:", filepath))
  
  # Read file as raw text lines to avoid parsing errors
  raw_lines <- readLines(filepath, warn = FALSE, skipNul = TRUE)
  
  if(length(raw_lines) == 0) stop("File is empty!")
  
  # Remove quotes and split by semicolon
  clean_lines <- gsub("\"", "", raw_lines)
  split_data <- strsplit(clean_lines, ";")
  
  # Reconstruct data frame
  header <- split_data[[1]]
  body   <- split_data[-1]
  
  n_cols <- length(header)
  body_fixed <- lapply(body, function(x) {
    length(x) <- n_cols
    return(x)
  })
  
  df <- data.frame(do.call(rbind, body_fixed), stringsAsFactors = FALSE)
  colnames(df) <- header
  
  # Extract specific columns: 2 (Department) and 4 (Resistance Rate)
  # RESAPATH format assumption: Col 2 = 'dpt', Col 4 = 'pourR'
  # Check if columns exist
  if(ncol(df) < 4) stop("Error: Data structure mismatch (less than 4 columns).")
  
  df_clean <- df[, c(2, 4)]
  colnames(df_clean) <- c("departement", "rate")
  
  # Clean and convert data types
  df_final <- df_clean %>%
    mutate(rate = as.character(rate)) %>%
    mutate(rate = na_if(rate, "-")) %>%       # Handle missing data
    mutate(rate = gsub(",", ".", rate)) %>%   # Handle French decimals
    mutate(rate = as.numeric(rate)) %>%
    na.omit()
  
  message(paste(" -> Successfully extracted:", nrow(df_final), "rows"))
  return(df_final)
}

# ==============================================================================
# Analysis 1: Colistin vs TMP-SMX (Trimethoprim-Sulfonamides)
# ==============================================================================

# File paths
file_col <- "data/bovin_colistin_2024_fixed.csv"
file_tmp <- "data/bovin_trimsalfa_2024_fixed.csv"

# Load and prepare data
df_col <- surgical_load(file_col) %>% rename(rate_colistin = rate)
df_tmp <- surgical_load(file_tmp) %>% rename(rate_tmp = rate)

# Merge datasets by Department
df_merged_tmp <- inner_join(df_col, df_tmp, by = "departement")
print(paste("Analysis 1 (TMP-SMX) dataset size:", nrow(df_merged_tmp)))

# Statistical Test (Pearson Correlation)
cor_res_tmp <- cor.test(df_merged_tmp$rate_tmp, df_merged_tmp$rate_colistin)
print(cor_res_tmp)

# Visualization
p1 <- ggplot(df_merged_tmp, aes(x = rate_tmp, y = rate_colistin)) +
  geom_point(color = "#2c3e50", alpha = 0.7, size = 3) +
  geom_smooth(method = "lm", color = "#e74c3c", fill = "#fadbd8") +
  geom_text_repel(aes(label = departement), size = 3, max.overlaps = 15) +
  labs(
    title = "Correlation: TMP-SMX vs Colistin Resistance",
    subtitle = paste("Bovine E. coli (2024) | Pearson r =", round(cor_res_tmp$estimate, 2), 
                     ", p =", format.pval(cor_res_tmp$p.value, digits=3)),
    x = "TMP-SMX Resistance (%)", y = "Colistin Resistance (%)",
    caption = "Data Source: RESAPATH"
  ) +
  theme_minimal()

print(p1)
ggsave("plots/correlation_plot_tmp.png", plot = p1, width = 8, height = 6, bg = "white")


# ==============================================================================
# Analysis 2: Colistin vs Tetracycline
# ==============================================================================

# File path for Tetracycline
file_tet <- "data/bovin_tetra_2024_fixed.csv"

# Load data
df_tet <- surgical_load(file_tet) %>% rename(rate_tetra = rate)

# Merge datasets
df_merged_tet <- inner_join(df_col, df_tet, by = "departement")
print(paste("Analysis 2 (Tetracycline) dataset size:", nrow(df_merged_tet)))

# Statistical Test
cor_res_tet <- cor.test(df_merged_tet$rate_tetra, df_merged_tet$rate_colistin)
print(cor_res_tet)

# Visualization
p2 <- ggplot(df_merged_tet, aes(x = rate_tetra, y = rate_colistin)) +
  geom_point(color = "#27ae60", alpha = 0.7, size = 3) + # Green points
  geom_smooth(method = "lm", color = "#2c3e50", fill = "#d5d8dc") +
  geom_text_repel(aes(label = departement), size = 3, max.overlaps = 15) +
  labs(
    title = "Correlation: Tetracycline vs Colistin Resistance",
    subtitle = paste("Bovine E. coli (2024) | Pearson r =", round(cor_res_tet$estimate, 2), 
                     ", p =", format.pval(cor_res_tet$p.value, digits=3)),
    x = "Tetracycline Resistance (%)", 
    y = "Colistin Resistance (%)",
    caption = "Data Source: RESAPATH"
  ) +
  theme_minimal()

print(p2)
ggsave("plots/correlation_plot_tetra.png", plot = p2, width = 8, height = 6, bg = "white") 


