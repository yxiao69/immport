
get_biosample <- function(studies="SDY80") {
#
# Function to get biosample.txt files for an arbitrary ImmPort study 
# Example: 
#
# > studies <- c("SDY80","SDY112")
# > sapply(studies,get_biosample(x))
#

  library(RMySQL)
  library(dplyr)
  library(readr)
  
  # The IP number of our MySQL databaase on the GCP 
  
  host <- "35.229.69.7"
  
  # Connect to our Google study base
  
  con <- dbConnect(RMySQL::MySQL(),user='study', 
                   password='allstudies', 
                   dbname='shared_data', 
                   host=host)
  
  for (ii in 1:length(studies)) { 
    
    # dbplyr does lazy connections so it doesn't actually fetch the data at first
    
    biosample <- tbl(con,"biosample")
    
    output <- biosample %>% filter(study_accession == studies[ii])
    #output_folder="~/Downloads"
   
    fname <- paste(studies[ii],"biosample.txt",sep="_")
    #out<- paste(output_folder,studies[ii],sep="/")
    #dir.create(out)
    #fname <- paste(output_folder,fname,sep="/")
    print(fname)
    # Because dbplyr does lazy evaluation we have to use the collection 
    # function to force it to pull the data out of the data base
    
    write_tsv(collect(output),fname)
    msg <- paste0("Processed ",studies[ii]," to filename: ",fname)
    print(msg)
  }
  dbDisconnect(con)
}

#sumdat=read_tsv("/Users/xiaoyi/Documents/2019spring/data mining/studies_by_data_types.tsv")
#studies <- as.vector(sumdat$study)


