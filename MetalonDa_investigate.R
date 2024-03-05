library(MetaLonDA)
library(readr)
library(dplyr)

#Load data
Bracken.table <- readr::read_tsv("./bracken_db1.tsv")
rownames(Bracken.table) <- Bracken.table$taxonomy_id
Bracken.table$"taxonomy_id" <- NULL
meta.table <- read.csv("./data/3in1FromMetadata.csv")
long.meta <- subset(meta.table, Study == "Longitudinal ")

#Sorting for each group (corrral)


get.samples.for.group <- function(group, meta.table){
  group.table <- subset(meta.table, Corral.ID == group)
  group.table <- select(group.table, "LabID", "CGR.sample.ID")
  rownames(group.table) <- group.table$LabID
  group.table$"LabID" <- NULL
  return(group.table)
}

rename.samples <- function(sample.data){
  #print(sample.data)
  renamed <- c()
  for (index in rownames(sample.data)){
    # For loop takes the index of the sample.data and gets the LabID and CGR.sample.ID,
    # then renames the file to the format "CGR.sample.ID-LabID_pe_CGR.sample.ID-LabID_db1.bracken"
    # fitting the outputs for taxprofiler. 
    #print(paste("index value:", index))
    sample.ID <- index
    CGR.ID <- sample.data[index,"CGR.sample.ID"]
    #print(paste("name", index, "ID:", CGR.ID))
    taxprofiler.output.name <- paste0(CGR.ID,"-",sample.ID,"_pe_",CGR.ID,"-", sample.ID, "_db1.bracken")
    renamed <- c(renamed, taxprofiler.output.name)
  }
  return(renamed)
}

count.vector <- function(sample.data){
  return(length(sample.data))
}

get.timepoints <- function(corral, meta.table){
  dataframe <- subset(meta.table, Corral.ID == corral)
  timepoints <- dataframe[["Day"]]
  return(timepoints)
}


# running functions

c20.samples <- rename.samples(get.samples.for.group(20, long.meta))
c20.count <- count.vector(c20.samples)
c20.timepoints <- get.timepoints(20, long.meta)
print(c20.timepoints)
print(count.vector(c20.timepoints))

c21.samples <- rename.samples(get.samples.for.group(21, long.meta))
c23.samples <- rename.samples(get.samples.for.group(23, long.meta))
print(c20.samples)
print(c20.count)

print(c21.samples)
print(c23.samples)

## Attempt 2 at getting vectors

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