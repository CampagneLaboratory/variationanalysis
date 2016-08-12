# Use this script to look at features that are very different between two groups of predictions (i.e., correct, wrong).
require("data.table")
require("dtplyr")
require("tidyr")
require("dplyr")

wrong <- read.delim("/Users/rct66/data/VariationAnalysisfeatures/rightg95features", header=TRUE)
right <- read.delim("/Users/rct66/data/VariationAnalysisfeatures/wrongl05features", header=TRUE)

rightG = right
wrongG = wrong

rightG$id <- TRUE
wrongG$id <- FALSE



both <- union(rightG,wrongG)



both %>%
  summarise_each(funs(t.test(.[id == TRUE],.[id == FALSE])$p.value), vars = normGermlineReverseCount0:invSomaticForwardCount4)