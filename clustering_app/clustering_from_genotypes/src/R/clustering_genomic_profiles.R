########### Description #######################################################
# jan 17 2023 IHC

# Apply Hierarchical clustering non supervised for variants

# clean global enviroment 
rm(list = ls())

# rda
load("clustering_genomic_profiles.Rda")
save.image("clustering_genomic_profiles.Rda")

########## Libraries ##########################################################

# for visualization of dendrograms
library(factoextra)
library(dendextend)
# for getting clusters n
library(NbClust)
# for heatmap.2 function
library(gplots)

# loading data ------------------------------------------------------------

# loading sample table with variant profiles:
sample_table <-
  read.delim(
    # "genotypes_sigexp_variants.tsv",
    "genotypes_exp_val_genes.tsv",
    header = T,
    stringsAsFactors = F,
    sep = "\t",
    row.names = "v"
  )
head(sample_table)

# Delete columns not required
sample_table <- sample_table[,-151]
head(sample_table)

# using dplyr
library(dplyr)
sample_table %>% 
  select_if(~ !any(is.na(.)))

# base R
sample_table <- sample_table[, colSums(is.na(sample_table)) == 0]

# transform to numeric
sample_df <-
  as.data.frame(sapply(sample_table, as.numeric), 
                row.names = rownames(sample_table))
# check:
sample_df[sample_df != 2 & sample_df != 1 & sample_df != 0]
sample_df[sample_df != 1 & sample_df != 0] # is there homozygous?
head(sample_df)

# transpose if required (variants should go in rows, samples in cols)
mat_df = as.data.frame(t(mat_df))
mat_df

# if you don't want to extract controls from the study
sample_df_2 <- sample_df

# if you want to extract controls
sample_df_2 = sample_df[grep("CONTROL",
                       colnames(sample_df),
                       invert = TRUE ,
                       value = T)]
colnames(sample_df_2)


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
head(sim_mat)

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
fit = hclust(d = as.dist(dis_mat), method = "ward.D2")
plot(fit)

# check distribution of values without controls
dis_mat_2 <-
  dis_mat[!grepl(colnames(dis_mat), pattern = "CONTROL"), !grepl(colnames(dis_mat), pattern = "CONTROL")]
hist(dis_mat_2)
boxplot(dis_mat_2[!grepl(colnames(dis_mat_2), pattern = "CONTROL")])

# jaccard genotypes with controls

### add partition to mat original

# with jaccard 66 (0/1/2) with controls
dis_mat_2 <- dis_mat

# hierarchical clustering -------------------------------------------------

# ball or hartigan
index = "ball"

res <-
  NbClust(
    data = dis_mat_2,
    diss = as.dist(dis_mat_2),
    distance = NULL,
    method = "ward.D2",
    index = index,
    min.nc = 3
  )
res$Best.nc
# 3 clusters

res$All.index
res$Best.partition

fit = hclust(d = as.dist(dis_mat_2), method = "ward.D2")
cutree(fit, k = 4)

all(res$Best.partition == cutree(fit, k = 4))

plot(fit, cex = 0.6) # plot tree
rect.hclust(fit, k = 4, border = 2:5) # add rectangle


# clustering customization -------------------------------------------------

# non supervised herarquical clustering
fit = hclust(as.dist(dis_mat), method = "ward.D2")
# save order of hclust
patients_order = colnames(sample_df_2)[fit$order]

# att of nodes for visualization
nodePar <-
  list(
    lab.cex = 0.6,
    cex = 0.7,
    col = "brown1",
    pch = c(NA, 19)
  )

plot(
  as.dendrogram(fit),
  nodePar = nodePar,
  horiz = FALSE,
  xlab = "Height",
  # hang = -1,
  edgePar = list(
    col = c("deepskyblue4", "deepskyblue2"),
    lwd = 2:1
  )
)
# cex = 0.1) # jerarquico, top-bottom

# save dendogram
# save(fit, file = "dendrogram.Rda")

# use dendextend package to colour better the dendrogram
dend1 <- color_branches(fit, k = 3)
dend1 <- color_labels(fit, k = 3)
plot(dend1)
set.seed(5665)

dev.off()

# simple separation by class
# controls
color_labels(as.dendrogram(fit),
             labels =
               fit$labels[grep(fit$labels, pattern = "CONTROL")],
             col = c("yellow")) %>% plot
color_labels(as.dendrogram(fit),
             labels = fit$labels[grep(fit$labels, pattern = "FOO")],
             col = c("red")) %>% plot
color_labels(as.dendrogram(fit),
             labels = fit$labels[grep(fit$labels, pattern = "FOP")],
             col = c("black")) %>% plot

rainbow(length(1:(res$Best.nc[1][[1]])))

# using fviz
fviz_dend(
  x = as.dendrogram(fit),
  k = res$Best.nc[1][[1]],
  k_colors =rainbow(length(1:(res$Best.nc[1][[1]]))),
  # k_colors = c("#2E9FDF", "#4c4c4c", "#E7B800"),
  # k_colors = c("#4c4c4c", "#4c4c4c", "#4c4c4c"), 
  # k_colors = c("#4c4c4c", "#00AFBB"), 2 clusters
  rect = TRUE,
  # rect_border = c("#2E9FDF", "#00AFBB", "#E7B800"),
  # rect_border = c("#FF7F50", "#00AFBB", "#E7B800"),
  rect_fill = TRUE,
  cex = 0.6,
  # main = "Dendrogram - ward.D2",
  # sub = "Genomic stratification",
  xlab = "samples",
  ylab = "height",
  color_labels_by_k = F,
  labels_track_height = 0.2
)


# deeper clustering study -------------------------------------------------

mycl <-
  cutree(fit, h = max(0.55)) # h = max(fit$height/0.5)) # 3 clusters

mycl <- cutree(fit, k = 4)
# cutree returns a vector of cluster membership
# in the order of the original data rows
# examine it
mycl

# examine the cluster membership by it's order
# in the dendrogram
mycl[fit$order]

# get clusters
cluster_A <- names(mycl[mycl == 4])
cluster_B <- names(mycl[mycl == 2])
cluster_C <- names(mycl[mycl == 1])
cluster_D <- names(mycl[mycl == 3])

cluster_num <- c(4, 2, 1, 3)
lapply(cluster_num, function(i){
  print(length(mycl[mycl == i]))
})

# save Rda
# save(mycl, file = "clusters_dendrogram_top6.Rda")

clusterCols <- rainbow(length(unique(mycl)))
myClusterSideBar <- clusterCols[mycl]
# myheatcol <- rev(redgreen(75))

# draw the heat map + dendrogram
heatmap.2(
  dis_mat,
  main = "Hierarchical Cluster",
  Rowv = as.dendrogram(fit),
  Colv = NA,
  dendrogram = "row",
  scale = "row",
  col = "heat.colors",
  density.info = "none",
  trace = "none",
  RowSideColors
)

# grab a cluster
dist_mat <- dis_mat[fit$order,]
cluster1 <- dis_mat[mycl == 2,]
View(cluster1)
# or simply add the cluster ID to your data
foo <- cbind(dist_mat, clusterID = mycl)

# examine the data with cluster ids attached, and ordered like the heat map
View(foo[fit$order,])

# change order of clusters ----------------------------------------------------------------

fit = hclust(d = as.dist(dis_mat), method = "ward.D2")

# change order
dd <- fit$order

a <- fit$order[1:17]
b <- fit$order[18:49]
c <- fit$order[50:150]

new_order <- vector()
new_order <- append(a, c)
new_order <- append(new_order, b)
fit$order <- new_order

plot(
  fit,
  hang = -1,
  ylab = "Height",
  xlab = "Samples",
  cex = 0.75
)

dend <- as.dendrogram(fit)
order.dendrogram(dend)

fviz_dend(
  x = fit,
  k = 3,
  # k_colors = c("#2E9FDF", "#00AFBB", "#E7B800"),
  k_colors = c("#4c4c4c", "#4c4c4c", "#4c4c4c"),
  rect = TRUE,
  # rect_border = c("#2E9FDF", "#00AFBB", "#E7B800"),
  rect_border = c("#FF7F50", "#00AFBB", "#E7B800"),
  rect_fill = TRUE,
  cex = 0.7,
  #main = "Dendrograma - ward.D2",
  xlab = "samples",
  ylab = "height",
  # sub = "Genomic stratification",
  color_labels_by_k = F,
  labels_track_height = 0.2
)

ggdendrogram(fit, rotate = F, size = 2)

dend <- as.dendrogram(fit)


# barplot accumulation of variants -------------------------------------------------------------------

# create a dataframe samples\group
temp <- list()
temp <- lapply(1:(res$Best.nc[1][[1]]), function(i){
  fit$labels[res$Best.partition == i ]
})
temp

samples_groups <- do.call('rbind', temp) 
samples_groups

samples_groups <- data.frame(res$Best.partition)
colnames(samples_groups) <- "group"
head(samples_groups)
write.table(samples_groups, "./clusters_val_exp.tsv", sep = "\t",
            row.names = T, col.names = T, quote = F)

# extract order based on dendrogram
a <- fit$order[1:17]
a
A <- fit$labels[a]
A

b <- fit$order[18:118]
b
B <- fit$labels[b]
B

View(dis_mat_2[, A])
V <- colSums(dis_mat_2[, A])
V
sum(V)

View(dis_mat_2[, B])
V2 <- colSums(dis_mat_2[, B])
V2
sum(V2)

# convert str to int
matrix_counts <- as.data.frame(t(sample_df), stringsAsFactors = F)

# transform all homozygous genotypes to perform rowsums
matrix_counts[matrix_counts == 2] <- 1

matrix_counts <-
  as.data.frame(sapply(matrix_counts, as.numeric), row.names = rownames(matrix_counts))

matrix_counts$total <- rowSums(matrix_counts)
matrix_counts$total

matrix_counts_ <- matrix_counts[!grepl(rownames(matrix_counts), pattern = "CONTROL"), ]
matrix_counts_

# order by subtypes dendrogram
all <- append(a, b)
all
# tosave <- t(matrix_counts_[all, c("chr22_17450952_A_G", "total")])
tosave <- t(matrix_counts_[all, c("total")])
tosave

# write.table(tosave, "../pof/number_of_variants_by_patients.tsv", sep = "\t",
#             row.names = F, col.names = T, quote = F)

# mat_df4 <- read.table(
#   "workspace/projects/pof/results/20_for_paper/number_of_variants_by_patients.tsv",
#   header = T)

dev.off()

ggplot(matrix_counts_,
       (aes(x = (total), ))) +
  geom_histogram(colour = "black", binwidth = 1) +
  geom_text(
    aes(label = after_stat(count)),
    stat = "count",
    position = position_dodge(width = 1),
    vjust = -0.5,
    size = 6.5
  ) +
  xlab("Number of variants") + ylab("Number of patients") +
  theme_light(base_size = 24) +
  theme(axis.text.x = element_text(hjust = 0.5),
        axis.title.x = element_text(vjust = -1.0)) +
  scale_fill_manual(values = "paleturquoise") +
  # scale_x_continuous(breaks = c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32))
  scale_x_continuous(breaks = c(unique(matrix_counts_$total)))
# +   scale_x_continuous(breaks = c(seq(0,30)))

