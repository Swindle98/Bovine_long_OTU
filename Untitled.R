library(MetaLonDA)
library(readr)

#Load data
Bracken.table <- readr::read_tsv("./bracken_db1.tsv")
meta.table <- read.csv("./data/3in1FromMetadata.csv")
long.meta <- subset(meta.table, Study == "Longitudinal")


#Adding metadata
n.group <- 3 # Number of groups represented in the table
n.samples <- 1 # Number of samples in each group. 1 sample was taken from each pen per timepoint.
n.time <- 14
count <- max(Bracken.table)
group <- 
# running MetaLonDa

output.metalonda <- metalondaAll(
  Count = count,
  group = 
  
)