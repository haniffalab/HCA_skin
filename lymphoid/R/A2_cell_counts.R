# Reference: Montoro, D.T., Haber, A.L., Biton, M. et al. A revised airway epithelial hierarchy includes CFTR-expressing ionocytes. Nature 560, 319–324 (2018) doi:10.1038/s41586-018-0393-7
# Method: Testing for difference in labelled fraction for pulse-seq and conventional lineage tracing
# Link: https://www.nature.com/articles/s41586-018-0393-7
# Adapted from: https://github.com/adamh-broad/single_cell_airway/blob/d7741835738ed58a9f93f2febb4dd302d9a9455b/PulseSeq.md
###############################################################################
library(reshape2)
library(MASS)
###############################################################################
# x is a dataframe with the following columns:
# x$variable # celltype
# x$value # cell count
# x$batch # condition_donor
# x$timepoint # condition (e.g. adult, fetal)
x <- read.csv("../t_cell_counts_df_site.csv", row.names = 1)
x$variable <- factor(x$variable)
# Lines marked with '!' need to be adapted

# Define all your conditions that group your replicates (donors)
cond_list <- sort(as.vector(unique(x$timepoint)))
#  [1] "E_L_Tc"            "E_L_Tc IL13 IL22"  "E_L_Th"
#  [4] "E_L_Treg"          "E_NL_Tc"           "E_NL_Tc IL13 IL22"
#  [7] "E_NL_Th"           "E_NL_Treg"         "P_L_Tc"
# [10] "P_L_Tc17_Th17"     "P_L_Th"            "P_L_Treg"
# [13] "P_NL_Tc"           "P_NL_Tc17_Th17"    "P_NL_Th"
# [16] "P_NL_Treg"         "S_NL_Tc"           "S_NL_Th"
# [19] "S_NL_Treg"

# Define conditions you want to test:
cond_comparisons <- list(
  c(9, 13), # test, ref
  c(10, 14),
  c(11, 15),
  c(12, 16),
  c(1, 5),
  c(2, 6),
  c(3, 7),
  c(4, 8),
  c(9, 17), # healthy
  c(11, 18),
  c(12, 19),
  c(1, 17),
  c(3, 18),
  c(4, 19)
  )

pv_cond <- list()
# pv_cond_9vs13 <- vector() # P Tc
# pv_cond_10vs14 <- vector() # P Tc17_Th17
# pv_cond_11vs15 <- vector() # P Th
# pv_cond_12vs16 <- vector() # P Treg
#
# pv_cond_1vs5 <- vector() # E Tc
# pv_cond_2vs6 <- vector() # E Tc17_Th17
# pv_cond_3vs7 <- vector() # E Th
# pv_cond_4vs8 <- vector() # E Treg
###############################################################################
for(k in 1:length(cond_comparisons)) {
  print(k)
  print(cond_comparisons[[k]][1])
  pv_cond[[k]] <- vector()
}
# Compute the p-values for increasing labeled fraction (using negative binomial regression)
for(i in 1:length(unique(x$variable))) {
    cell = levels(x$variable)[i]
    totals = aggregate(x$value, by = list(x$batch), FUN = sum)
    colnames(totals) <- c("batch", "total")
    d = data.frame(subset(x, variable == cell))
    d <- merge(d, totals, by = "batch")

    if(any(d$total == 0)) {
      warning("Removing sample that didn't detect this celltype")
      d = subset(d, total > 0)
      }

    for(j in 1:length(cond_comparisons)) {
      print(j)
      reference <- cond_list[cond_comparisons[[j]][2]] # sorry
      test <- cond_list[cond_comparisons[[j]][1]]

      x1 = d[d$timepoint %in% c(reference, test), ]

      if(dim(x1)[1] == 0 || length(unique(x1$timepoint)) != 2) { # if there are no both groups
        pv_cond[[j]][i] <- "No data"
        print(paste0("No counts for ", cell))
        print(x1)
      } else {
        nb = MASS::glm.nb(formula = value ~ timepoint + offset(log(as.numeric(x1$total))), data = x1, maxit = 1000)#, control=glm.control(trace = 3))
        pv_cond[[j]][i] = anova(nb, test = "LRT")$`Pr(>Chi)`[2] # !
      }
    }
}


name_vector <- vector()
for(j in 1:length(cond_comparisons)) {
  print(j)
  reference <- cond_list[cond_comparisons[[j]][2]] # sorry
  test <- cond_list[cond_comparisons[[j]][1]]
  name_vector[j] <- paste0(test, "__vs__", reference)
}

pv <- do.call(cbind.data.frame, pv_cond)
rownames(pv) = levels(x$variable)
rownames(pv) <- c("Not expanded", "Expanded")
colnames(pv) <- name_vector
write.table(pv, file = "NB_pvals_results_site_all.txt", sep = "\t", quote = F)
write.csv(t(pv), file = "NB_pvals_results_site_all_transposed.csv")


# write.table(pv, file = "NB_pvals_results_site.txt", sep = "\t", quote = F)
# write.csv(t(pv), file = "NB_pvals_results_site_transposed.csv")
