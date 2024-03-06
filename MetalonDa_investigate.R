library(MetaLonDA)
library(readr)
library(dplyr)

#Load data
Bracken.table <- as.data.frame(readr::read_tsv("./bracken_db1.tsv"))
rownames(Bracken.table) <- Bracken.table$taxonomy_id
Bracken.table$"taxonomy_id" <- NULL
meta.table <- read.csv("./data/3in1FromMetadata.csv")
long.meta <- subset(meta.table, Study == "Longitudinal ")


# get vectors functions

count.vector <- function(sample.data){
  return(length(sample.data))
}



get.group <- function(meta.table, LabID){
  group <- meta.table[LabID, "Corral.ID"]
  print(class(group))
  if  (group %in% c("20", "21")){
    return("study")
  }
  else{
    return("control")
  }
}




get.vectors <- function(count.table, meta.table, count.index.columns=1){
  
  #print(count.table)
  rownames(meta.table) <- meta.table$LabID # sets the index to == LabID value for that row
  meta.table$"LabID" <- NULL # removes the LabID column
  
  samples <- colnames(count.table) # gets the column names of the count table
  #samples <- samples[-count.index.columns] #drops the header column for the taxonomy ID default is 1
  print(count.vector(samples))
  
  ID <- group <- timepoints <- c() # creates empty vectors for group and timepoints
  
  #group <- rep("Long", times = count.vector(samples))
               
  for (sample in samples){
    #for each sample in the samples vector, it gets the LabID, group, and timepoint from the metadata table
    # and appends them to the LabID, group, and timepoints vectors
    
    print(sample)
    LabID <- sub(".*-(L\\d+)_.*", "\\1", sample)
    group <- c(group, get.group(meta.table, LabID))
    timepoints <- c(timepoints, meta.table[LabID, "Day"])
    ID <- c(ID, paste0("c", meta.table[LabID, "Corral.ID"]))
  }
  return(list(group.vector = group, timepoint.vector = timepoints, ID.vector = ID))
}


# executing the get vectors functions 

vectors <- get.vectors(Bracken.table, long.meta)
print(vectors$ID.vector)
print(count.vector(vectors$ID.vector))
print(vectors$group.vector)
print(count.vector(vectors$group.vector))
print(vectors$timepoint.vector)
print(count.vector(vectors$timepoint.vector))


# running MetaLonDA

count = Bracken.table
group = vectors$group.vector
time = vectors$timepoint.vector
ID = vectors$ID.vector


result <- metalondaAll(
  Count = count, 
  Group = group, 
  Time = time, 
  ID = ID, 
  n.perm = 1000,
  fit.method = "nbinomial",
  pvalue.threshold = 0.05,
  adjust.method = "BH",
  time.unit = "days",
  norm.method = "none",
  prefix = "Output",
  ylabel = "Normalized Count",
  col = c("blue", "firebrick")
)
