library(MetaLonDA)
library(readr)
library(dplyr)



# Data import functions

taxprofiler.output.import <- function (tsv){
  #Imports files that have come out of taxprofiler and sets the row names to the species name or taxonomy_id
  #and removes the name and taxonomy_id columns. (Prefering species name over taxonomy_id)
  #returns dataframe
  
  df.table <- as.data.frame(readr::read_tsv(tsv))
  
  if ('name' %in% colnames(df.table)){
    # If the file contains a column with the name 'name' it sets the row names to the species names
    # and removes the name and taxonomy_id columns.
    
    print("[info] setting row names to species names")
    df.table <- df.table[!is.na(df.table$name),]
    rownames(df.table) <- df.table$name
    df.table$name <- NULL
    df.table$taxonomy_id <- NULL
  }
  else if ('taxonomy_id' %in% colnames(df.table)){
    print("[info] setting rownames to taxonomy_id")
    rownames(df.table) <- df.table$taxonomy_id
    df.table$taxonomy_id <- NULL
  }
  else{
    stop("The file does not contain a column with the name 'name' or 'taxonomy_id'")
  }
  return(df.table)
}
  

filter.num.features <- function(count.table, top.x.features){
  #Filters the count table to only include the top x features
  #returns dataframe
  
  print("[info] filtering count table to only include the top x features")
  return(count.table[order(rowSums(count.table), decreasing = TRUE)[1:top.x.features],])
}


# get vectors functions

data.frame.filter <- function(count.table, meta.table, groups){
  
  group.1 = groups[1]
  group.2 = groups[2]
  
  #Get the labIDs for the groups of interest
  meta.table.for.groups.of.interest <- meta.table[meta.table$Corral.ID %in% c(group.1, group.2),]
  LabIDs <- meta.table.for.groups.of.interest$LabID
  desired.col.names <- c()
  
  #Filter the count table to only include the samples from the groups of interest
  for (LabID in LabIDs){
    print(LabID)
    column <- grep(LabID, colnames(count.table))
    print(column)
    desired.col.names <- c(desired.col.names, column)
  }
  return(count.table[,desired.col.names])
  
}


count.vector <- function(sample.data){
  #Helper function that counts the number of items in a vector.
  #Used to keep code readable. 
  return(length(sample.data))
}


get.group <- function(meta.table, LabID, groups){
  group.1 = groups[1]
  group.2 = groups[2]
  group <- meta.table[LabID, "Corral.ID"]
  print(class(group))
  if  (group %in% c(group.1)){
    return(group.1)
  }
  
  else if (group %in% c(group.2)){
      return(group.2)
  }

  else{
    stop("More than the two desired groups are present in the metadata table")
  }
}


get.vectors <- function(count.table, meta.table, groups){
  # !!! Currently only works for taxprofiler outputs !!!
  
  #print(count.table)
  rownames(meta.table) <- meta.table$LabID # sets the index to == LabID value for that row
  meta.table$"LabID" <- NULL # removes the LabID column
  
  samples <- colnames(count.table) # gets the column names of the count table
  
  print(count.vector(samples))
  
  ID <- group <- timepoints <- c() # creates empty vectors for group and timepoints
  
  #group <- rep("Long", times = count.vector(samples))
               
  for (sample in samples){
    #for each sample in the samples vector, it gets the LabID, group, and timepoint from the metadata table
    # and appends them to the LabID, group, and timepoints vectors
    
    print(sample)
    LabID <- sub(".*-(L\\d+)_.*", "\\1", sample)
    group <- c(group, get.group(meta.table, LabID, groups))
    timepoints <- c(timepoints, meta.table[LabID, "Day"])
    ID <- c(ID, paste0("c", meta.table[LabID, "Corral.ID"]))
  }
  return(list(group.vector = group, timepoint.vector = timepoints, ID.vector = ID))
}


get.vector.stats <- function(count, group, time, ID){
  print(paste0("Count samples: ", ncol(count)))
  print(paste0("Group: ", count.vector(group)))
  print(paste0("Time: ", count.vector(time)))
  print(paste0("ID: ", count.vector(ID)))
}

#Load data
Bracken.table <- taxprofiler.output.import("./data/bracken_db1.tsv")
mOTU.table <- taxprofiler.output.import("./data/motus_db_mOTU.tsv")
meta.table <- read.csv("./data/3in1FromMetadata.csv")
long.meta <- subset(meta.table, Study == "Longitudinal ")


# running MetaLonDA


group.1 <- "20"
group.2 <- "21"
groups <- c(group.1, group.2)

count = mOTU.table
top.x.features <- 250 # Number of top features to be selected taxprofiler output tables are ranked by counts.
count <- filter.num.features(count, top.x.features) 

group.filtered.count <- data.frame.filter(count, long.meta, groups)


vectors = get.vectors(group.filtered.count, long.meta, groups)
group = vectors$group.vector
time = vectors$timepoint.vector
ID = vectors$ID.vector

get.vector.stats(group.filtered.count, group, time, ID)

result <- metalondaAll(
  Count = count, 
  Group = group, 
  Time = time, 
  ID = ID, 
  parall = TRUE,
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
