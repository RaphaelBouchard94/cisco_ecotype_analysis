rm(list = ls())

###############################################################################
# 03.2 - Sex ratio between ecotypes (chi-square test)
#
# CHANGES vs the previous version
#  1. Uses true_eco (genetic ecotype), like every other analysis. The old
#     script used the `ecotype` column. Both are printed so you can compare.
#  2. Keeps only F and M, so the table is 2 x 2 (1 df).
#  3. Prints a ready-to-paste result with df and an unrounded p-value.
#     For 2 x 2 tables chisq.test() applies Yates' continuity correction by
#     default; mention it in the Methods.
###############################################################################

library(tidyverse)

data <- read.csv("data/data_with_gillrakers.csv", header = TRUE, sep = ",")

data_fm <- data %>% filter(sex %in% c("F", "M"))

sex_eco <- table(ecotype = data_fm$true_eco, sex = data_fm$sex)
print(sex_eco)

csq <- chisq.test(sex_eco)
csq
csq$expected

cat(sprintf("Sex ratio vs genetic ecotype: chi2(%d) = %.3f, p = %.3f (n = %d)\n",
            csq$parameter, csq$statistic, csq$p.value, sum(sex_eco)))

# For comparison: the column the old script used
if ("ecotype" %in% names(data_fm)) {
  csq_old <- chisq.test(table(data_fm$ecotype, data_fm$sex))
  cat(sprintf("Same test with the `ecotype` column (old script): chi2(%d) = %.3f, p = %.3f\n",
              csq_old$parameter, csq_old$statistic, csq_old$p.value))
}
