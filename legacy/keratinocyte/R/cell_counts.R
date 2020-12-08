# Reference: Montoro, D.T., Haber, A.L., Biton, M. et al. A revised airway epithelial hierarchy includes CFTR-expressing ionocytes. Nature 560, 319–324 (2018) doi:10.1038/s41586-018-0393-7
# Method: Testing for difference in labelled fraction for pulse-seq and conventional lineage tracing
# Link: https://www.nature.com/articles/s41586-018-0393-7
# Source: https://github.com/adamh-broad/single_cell_airway/blob/d7741835738ed58a9f93f2febb4dd302d9a9455b/PulseSeq.md
###############################################################################
library(reshape2)
library(MASS)
###############################################################################
# x is a dataframe with the following columns:
# x$variable # celltype
# x$value # cell count
# x$batch # condition_donor
# x$timepoint # condition (e.g. adult, fetal)
x <- read.csv("../cell_counts_df.csv", row.names = 1)

# Define all your conditions that group your replicates (donors)
cond_1 <- "P_NL"
cond_2 <- "P_L"
cond_3 <- "E_NL"
cond_4 <- "E_L"

# Define as many empty variables as many conditions you want to test:
pv_cond_2vs1 <- vector()
pv_cond_4vs3 <- vector()
###############################################################################
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

    { # Copy this block as many times as many condition-pairs you test
    reference <- cond_1 # !
    test <- cond_2 # !

    x1 = d[d$timepoint %in% c(reference, test), ]
    nb = MASS::glm.nb(formula = value ~ timepoint + offset(log(as.numeric(x1$total))), data = x1, maxit = 1000)#, control=glm.control(trace = 3))
    pv_cond_2vs1[i] = anova(nb, test = "LRT")$`Pr(>Chi)`[2] # !
    }

    { # Copy this block as many times as many condition-pairs you test
    reference <- cond_3 # !
    test <- cond_4 # !

    x1 = d[d$timepoint %in% c(reference, test), ]
    nb = MASS::glm.nb(formula = value ~ timepoint + offset(log(as.numeric(x1$total))), data = x1, maxit = 1000)#, control=glm.control(trace = 3))
    pv_cond_4vs3[i] = anova(nb, test = "LRT")$`Pr(>Chi)`[2] # !
    }

}

pv = cbind.data.frame(pv_cond_2vs1, pv_cond_4vs3) # !
rownames(pv) = levels(x$variable)
write.table(pv, file = "NB_pvals_results.txt", sep = "\t", quote = F)
