# libraries
suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# arguments
parser$add_argument("-f", "--genotype_file", type="character")
parser$add_argument("-n", "--number_of_clusters", type="integer")
parser$add_argument("-od", "--output_dir", type="character")

# get command line options
args <- parser$parse_args()

# for visualization of dendrograms
library(factoextra)
# for getting clusters n
library(NbClust)


print("Starting Workflow")

# main function
main <- function() {
print("Importing Dataframe")
  # args <- commandArgs(trailingOnly = TRUE)
  filename <- args$genotype_file
  sample_table <-
    read.delim(
      filename,
      header = T,
      stringsAsFactors = F,
      sep = "\t",
      row.names = "v")
      
  print("Manipulating Data")
  sample_table <- sample_table[, colSums(is.na(sample_table)) == 0]
  # transform to numeric
  sample_df <-
    as.data.frame(sapply(sample_table, as.numeric), 
                  row.names = rownames(sample_table))
  # extract controls
  sample_df_2 = sample_df[grep("CONTROL",
                               colnames(sample_df),
                               invert = TRUE ,
                               value = T)]
  # jaccard index -----------------------------------------------------------
  
  # do all possible pair combinations
  combinations = combn(colnames(sample_df_2), 2)
  (combinations)[0:6]
  
  # generate similarity matrix (samples vs samples)
  sim_mat = matrix(
    NA,
    ncol = ncol(sample_df_2),
    nrow = ncol(sample_df_2),
    dimnames = list(colnames(sample_df_2), colnames(sample_df_2))
  )
  
  for (i in 1:ncol(combinations)) {
    sim_mat[combinations[2, i], combinations[1, i]] =
      (sum(sample_df_2[, combinations[1, i]] != 0  &
             sample_df_2[, combinations[1, i]] ==
             sample_df_2[, combinations[2, i]])) /
      sum(sample_df_2[, combinations[2, i]] != 0 |
            sample_df_2[, combinations[1, i]] != 0)
    sim_mat[combinations[1, i], combinations[2, i]] =
      (sum(sample_df_2[, combinations[1, i]] != 0  &
             sample_df_2[, combinations[1, i]] ==
             sample_df_2[, combinations[2, i]])) /
      sum(sample_df_2[, combinations[2, i]] != 0 |
            sample_df_2[, combinations[1, i]] != 0)
  }
  diag(sim_mat) = 1
  sim_mat[is.nan(sim_mat) | is.na(sim_mat)] <- 0
  dis_mat = 1 - sim_mat
  
  # hierarchical clustering -------------------------------------------------
  
  # ball or hartigan
  index = "ball"
  
  res <-
    NbClust(
      data = dis_mat,
      diss = as.dist(dis_mat),
      distance = NULL,
      method = "ward.D2",
      index = index,
      min.nc = args$number_of_clusters,
    )
  fit = hclust(d = as.dist(dis_mat), method = "ward.D2")
  
  print("Creating Plots")
  # using fviz
  suppressWarnings(fviz_dend(
    x = as.dendrogram(fit),
    k = res$Best.nc[1][[1]],
    k_colors =rainbow(length(1:(res$Best.nc[1][[1]]))),
    rect = TRUE,
    rect_fill = TRUE,
    cex = 0.6,
    color_labels_by_k = F,
    labels_track_height = 0.2
  ))
  output_dir <- args$output_dir
  file_type <- 'tiff'
  ggsave(
    #filename = args$output_file,
    filename = paste0(output_dir, "dendrogram_output.", file_type),
    device = "tiff",
    plot = last_plot(),
    path = NULL,
    scale = 1,
    width = 20,
    height = 12,
    units = "in",
    dpi = 300,
    limitsize = TRUE,
    bg = NULL,
  )

# barplot accumulation of variants ----------------------------------------

  # convert str to int
  matrix_counts <- as.data.frame(t(sample_df), stringsAsFactors = F)
  
  # transform all homozygous genotypes to perform rowsums
  matrix_counts[matrix_counts == 2] <- 1
  matrix_counts <-
    as.data.frame(sapply(matrix_counts, as.numeric), 
                  row.names = rownames(matrix_counts))
  
  matrix_counts$total <- rowSums(matrix_counts)
  matrix_counts_ <- matrix_counts[!grepl(rownames(matrix_counts), 
                                         pattern = "CONTROL"), ]
  
  ggplot(matrix_counts_,
         (aes(x = (total), ))) +
    geom_histogram(colour = "cadetblue", binwidth = 1) +
    geom_text(
      aes(label = after_stat(count)),
      stat = "count",
      position = position_dodge(width = 1),
      vjust = -0.5,
      size = 6.5
    ) +
    xlab("Nº variants accumulated") + ylab("Nº of patients") +
    theme_light(base_size = 24) +
    theme(axis.text.x = element_text(hjust = 0.5),
          axis.title.x = element_text(vjust = -1.0)) +
    scale_fill_manual(values = "paleturquoise") +
    scale_x_continuous(breaks = c(unique(matrix_counts_$total)))
  
  file_type <- 'tiff'
  ggsave(
    filename = paste0(output_dir, "barplot_output.", file_type),
    device = "tiff",
    plot = last_plot(),
    path = NULL,
    scale = 1,
    width = 20,
    height = 12,
    units = "in",
    dpi = 300,
    limitsize = TRUE,
    bg = NULL,
  )

  print("Worflow Finished")
}

main()
