library(reshape2)
library(dendsort)
library(dendextend)
library(dendextendRcpp)
library(dplyr)


#cluster
distance_CDRc <- acast(Combined_CDR_AbTU_dis10_2, V1 ~ V2, mean, value.var = "V3")
dist_CDRc <- as.dist(distance_CDRc)
fit_CDRc <- hclust(dist_CDRc, method = "complete")

#make type dend
dend_CDRc <- as.dendrogram(fit_CDRc)

#color branches by label
Ab_CDRc_list <- pull(Ab_s_combined_CDRs, V1)

dend_CDRc %>% sort_smallest() %>% set("by_labels_branches_col", value = Ab_CDRc_list, TF_values = c(6,Inf)) -> painted_Ab_CDRc

#hang dendogram for better viewing
hang.dendrogram(painted_Ab_CDRc, hang = 0) -> painted_Ab_CDRc
plot(painted_Ab_CDRc, leaflab = "none", ylab ="Edit Distance Between CDRcs", main = "Clustering of 8,119 Full Length sdAbs")
abline(h=14.5)


Combined_CDR_AbTU_dis10_2 <- filter(Combined_CDR_AbTU_dis10_2, V3 > 0)
Combined_CDR_AbTU_dis10_2 <- filter(Combined_CDR_AbTU_dis10_2, V3 < 14)
