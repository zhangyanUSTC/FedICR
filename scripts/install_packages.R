packages <- c("MASS", "clue", "doParallel", "ggplot2", "gridExtra", "foreach",
              "cowplot", "combinat", "matrixStats", "mclust", "Matrix", "fossil", "haven")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
cat("Required packages are available.\n")
