

##fcs,elisa,elispot
fcs.convertor<-function (names,col1,col2,type) {
  flist <- list()
  library(readr)
  for (ii in 1:length(names)) {
    tmp <- CreateConnection(names[ii])
    flist[[ii]] <- tmp$getDataset(type)
  }
  outname=list()
  # Write out files 
  for (ii in 1:length(flist)) {
    # Write out the files
    outname[[ii]] <- paste(names[ii],type,".tsv",sep="_")
    #  fullpath <- paste0(path,outname,sep="")
    if (nrow(flist[[ii]]) == 0) {
      str <- paste(outname[[ii]],"does not have any rows to write",sep=" ")
     # print(str)
    } else {
      write_tsv(flist[[ii]],path=outname[[ii]])
    }
    
  } #end write files
  
 # files <- list.files(pattern=paste0("*",type,"*"))
 
files=as.character(outname)
print(files)
  for (ii in 1:length(flist)) {
    fcs <- data.frame(flist[[ii]])
    ids <- sapply(strsplit(fcs$`participant_id`,"\\."),`[`,1)
    observation_ID <- gsub(" ","",paste(ids,format(fcs$`study_time_collected`,nsmall=1),sep="_"))
    
    # All we need is the first 10 columns from the data frame
    
    fcs.mat <- cbind(observation_ID,fcs)
    fcs.mat$observation_ID <- as.character(fcs.mat$observation_ID)
    fcs.mat$`Participant ID` <- ids 
    data.matrix <- fcs.mat %>% 
      group_by(UQ(as.name(col1)),observation_ID) %>%#`population_definition_reported`
      summarize(med=median(UQ(as.name(col2)))) %>% #`population_cell_number`
      spread(observation_ID,med) 
    
    fname <- paste(strsplit(files[ii],"_")[[1]][1],type,"data_matrix.tsv",sep="_")
    write_tsv(data.matrix,fname)
  }
  
  for (ii in 1:length(files)) {
    fcs <-data.frame(flist[[ii]])
    ids <- sapply(strsplit(fcs$`participant_id`,"\\."),`[`,1)
    observation_ID <- gsub(" ","",paste(ids,format(fcs$`study_time_collected`,nsmall=1),sep="_"))
    
    # All we need is the first 10 columns from the data frame
    
    fcs.tmp <- cbind(observation_ID,fcs)
    fcs.tmp$observation_ID <- as.character(fcs.tmp$observation_ID)
    fcs.tmp$`Participant ID` <- ids 
    
    fcs.2 <- fcs.tmp[,1:8]
    
    # Are there duplicates ?
    sum(duplicated(fcs.2))
    
    # Get rid of duplicated rows (except the first occurrence thereof)
    fcs.3 <- fcs.2[!duplicated(fcs.2),]
    
    fname <- paste(strsplit(files[ii],"_")[[1]][1],type,"annotation.tsv",sep="_")
    write_tsv(fcs.3,fname)
  }
}

## biosample

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
    

    
    biosample <- tbl(con,"biosample")
    
    output <- biosample %>% filter(study_accession == studies[ii])

    
    fname <- paste(studies[ii],"biosample.txt",sep="_")

    print(fname)
    
    
    write_tsv(collect(output),fname)
    msg <- paste0("Processed ",studies[ii]," to filename: ",fname)
    print(msg)
  }
  dbDisconnect(con)
}

## gene expression matrix
library(Biobase)
library(readr)

expression_fetch <- function(study="SDY80") {
  tmp <- CreateConnection(study)
  fmat=tmp$listGEMatrices()
  if (!is.null(fmat)){
    flist=as.vector(fmat[,name])
    

    
    for( i in 1:length(flist)){
      outmat=paste0(flist[i],'_Expression_Matrices.tsv')
      obj<-tmp$getGEMatrix(flist[i])
      mat<-data.frame(exprs(obj))
      mat_rowname=data.frame(cbind(gene=row.names(mat),mat),stringsAsFactors = FALSE)
      #row.names(mat_rowname)=NULL
      write_tsv(mat_rowname,path=outmat)
      
    }
  }else{print(paste0("No gene expression data in ", study))
  }
}

## antibody

antibody_fetch <- function(study="SDY80",type="neut_ab_titer") {
  #
  # Take a study name and a data type and trues to fetch it 
  #
  # INPUT: study - a name of a valid study
  #        type - "neut_ab_titer" or "hai"
  #
  # OUTPUT: two output files:
  #         1) an annotation file prefixed with the study name
  #         2) a matrix file prefixed with the study name
  #
  # Example:
  # setwd("~/Downloads") # Change this to a working directory
  #
  # sdy135.hai  <- antibody_fetch(study="SDY315",type="hai")
  # sdy135.neut <- antibody_fetch(study="SDY315",type="neut_ab_titer")
  # 
  # OR could use in a call to sapply as a loop
  #
  # studies <- c("SDY63", "SDY80", "SDY112", "SDY144", "SDY180")
  # neut <- sapply(studies,function(x) antibody_fetch(x,type="neut_ab_titer"))
  # 
  
  library(forcats)
  library(ImmuneSpaceR)
  library(tidyr)
  
  # This uses the Puthon code sample Shuzhao placed at 
  #  https://storage.cloud.google.com/immport/Shuzhao_example_SDY80/README.txt
  # 
  
  tmp <- CreateConnection(study)
  df  <- tmp$getDataset(type)
  
  # Check to see if anything was returned for the desired data type
  
  if (nrow(df) > 0) {
    virus <- df$virus
    #    unique(virus)
    
    # So as with the fcs data we want to create a variable called observation_id that is based on a combination of
    # the participant id and time point. 
    
    ids <- sapply(strsplit(df$`participant_id`,"\\."),`[`,1)
    observation_id <- gsub(" ","",
                           paste(ids,format(df$`study_time_collected`,nsmall=1),
                                 sep="_"))
    
    # Create a new data frame that has the observation_id as the first column. 
    
    df.2 <- data.frame(cbind(observation_id=observation_id,df),stringsAsFactors = FALSE)
    
    
    # Pick out the desired row names as specified by the code sample
    # https://storage.cloud.google.com/immport/Shuzhao_example_SDY80/README.txt
    
    df.3 <- df.2[,c("observation_id","participant_id","age_reported","gender","race","cohort",
                    "study_time_collected","study_time_collected_unit")]
    
    # Get just the the non duplicated rows
    
    df.4 <- df.3[!duplicated(df.3),]
    
    # Write this out as the annotation file
    
    anno.out <- paste(study,type,sep="_")
    
    anno.out <- paste(anno.out,"annotation.txt",sep="_")
    
    write_tsv(df.4,path=anno.out)
    
    # Now let's write out the matrix portion. We will use df.2 since it has all the data
    # Now, the code at https://storage.cloud.google.com/immport/Shuzhao_example_SDY80/README.txt 
    # wants to collape the virus names into a normalized set of four names. We can use the forcats
    # package to simplify this. 
    
    # Make virus a factor which is required by the forcats package
    
    df.2$virus <- as.factor(df.2$virus)
    
    # Normalize Virus strains by collapsing the existing virus names into a smaller set.
    # Notice that this specific to the study and data type and would need to be adjusted so
    # this code will fail for example if you tried to access SDY180 hai data
    
    #  virus.tmp <- fct_collapse(df.2$virus,
    #                            Brisbane_A=c("A/Brisbane/59/2007"),
    #                            California=c("A/California/7/2009","A_Ca_07_swine"),
    #                            Uruguay=c("A/Uruguay/716/2007","A_Uruguay_716_2007"),
    #                            Brisbane_B=c("B/Brisbane/60/2008","B_Brisbane_60_2001"))
    #
    
    virus.tmp <- fct_collapse(df.2$virus,
                              Brisbane_A=c("A/Brisbane/59/2007","A/Brisbane/10/2007"),
                              California=c("A/California/7/2009","A_Ca_07_swine"),
                              Uruguay=c("A/Uruguay/716/2007","A_Uruguay_716_2007"),
                              Brisbane_B=c("B/Brisbane/60/2008","B_Brisbane_60_2001","B/Brisbane/03/2007"),
                              Perth=c("A/Perth/16/2009","A/Perth/19/09 H3N2"),
                              Wisconsin=c("B/Wisconsin/01/2010"),
                              Victoria=c("A/Victoria/361/2011","A/Victoria/3/1975"),
                              Massachusetts="B/Massachusetts/02/2012",
                              Florida="B/Florida/4/2006",
                              South_Dakota="A/South Dakota/06/2007",
                              Lee="B/Lee/1940",
                              Influenza_A_HIN1=c("A/PR/8/34 (H1N1)","influenza A H1N1"),
                              Influenza_B="influenza B")
    
    # Create a new column at the end of the data frame with the newly collpsed normalized virus names
    # Now we want to create a table of the virus names by the observation ideas and then sum the 
    # value_preferred column
    
    if (type == "neut_ab_titer") {
      df.2$neut_ab_titer <- virus.tmp
      df.2 %>% group_by(neut_ab_titer,observation_id) %>% summarize(sum=sum(value_preferred))
      
      out <- df.2 %>% group_by(neut_ab_titer,observation_id) %>% 
        summarize(sum=sum(value_preferred)) %>% 
        spread(observation_id,sum)
      
    } else {
      df.2$hai <- virus.tmp
      df.2 %>% group_by(hai,observation_id) %>% summarize(sum=sum(value_preferred))
      
      out <- df.2 %>% group_by(hai,observation_id) %>% 
        summarize(sum=sum(value_preferred)) %>% 
        spread(observation_id,sum)
      
    }
    
    matrix.out <- paste(study,type,sep="_")
    matrix.out <- paste(matrix.out,"data_matrix.txt",sep="_")
    write_tsv(out,matrix.out)
    #    return(out)
  } else {
    # If there is nothing to process keep going
    print(paste("Study",study,"has no data of type",type,sep=" "))
  }
}
 
data.explorer <- function(studies=study.names,data_types=data.types,write=TRUE) {
  #
  # Example: 
  # > outmat <- data.explorer(studies=study.names,data_types=data.types,write=TRUE)


  
  library(ImmuneSpaceR)
  library(dplyr)
  library(readr)
  
  # Create a matrix whose row names are the studies and the column names
  # are the data_types
  
  outmat <- matrix("Y",nrow=length(studies),ncol=length(data_types))
  
  row.names(outmat) <- studies
  colnames(outmat) <- data_types
  
  # For each study (row) - connect to ImmuneSpace
  for (ii in 1:length(studies)) {
    tmp <- CreateConnection(studies[ii])
    
    # Now for each study type try to see if that data type exists
    # If it does then put a "Y" in the correspinding element. If not
    # then put "N"
    
    for (jj in 1:length(data_types)) {
      print(paste("Processing",data_types[jj],sep=" "))
      ds <- tmp$getDataset(data_types[jj])
      print(nrow(ds))
      if (nrow(ds) == 0) {
        outmat[ii,jj] = "N"
      } else {
        outmat[ii,jj] = "Y"
      }
    }
  }
  
  # Okay, let's write out the file which most users will want to read using
  # Excel so let's make it tab delimited. Also, we have to do some conversion
  # on the matrix into a data frame.
  
  om <- data.frame(outmat,stringsAsFactors = FALSE)
  om <- cbind(study=rownames(om),om)
  rownames(om) <- NULL
  
  if (write) {
    write_tsv(om,path="studies_by_data_types.tsv")
  }
  return(om)
}

## function to extract for given study 
## notice one type each time 
convertor <- function(types='elispot',study=c('SDY80','SDY296'),path="~/Downloads") {
  library(ImmuneSpaceR)
  library(tidyr)
  library(dplyr)
  #find exist data type for study
  outmat <- data.explorer(studies=study,data_types=types,write=FALSE)
  mainDir=path
  
for (i in 1:length(study)){
  subDir=study[i]
  if (file.exists(file.path(mainDir, subDir))) {
    print(paste0(subDir," exists in ",mainDir," and is a directory"))
  } else{
    dir.create(file.path(mainDir, subDir))
  }
  setwd(file.path(mainDir, subDir))
  #print the exit data type 
  print(outmat)
for ( type in types) {
  # Determine file type
  print(outmat[type])
  if (outmat[type]=='Y'){
  if(type=='gene_expression_files') {
    
      
      expre<-sapply(study,function(x) expression_fetch(x))
      sapply(study,function(x)get_biosample(x))
      
  } else if (type=='neut_ab_titer'|type=='hai') {
    
     
      sapply(study,function(x) antibody_fetch(x,type=type))
      
  } else { 
    
      if(type=='fcs_analyzed_result') {
      
           col1='population_name_reported' #group by col1
           col2='population_cell_number'
       
       } else if(type=='elispot') {
         
           col1='analyte' #group by col1
           col2='spot_number_reported' #the value take median 
          
       }  else if(type=='elisa') {
           col1='analyte' #group by col1
           col2='value_reported' #the value take median
       }
  
    fcs.convertor(study,col1,col2,type)
  
    }
  }
}

  file=list.files()
  for (i in 1:length(file)){
    print(file[i])
    f=read.table(file =file[i] , sep = '\t', header = TRUE,fill=TRUE)
    if("participant_id"%in%colnames(f)){
      f[["participant_id"]]=gsub('\\.[0-9]+','',f[["participant_id"]])
      write_tsv(f,path=file[i])
    }
  }
  
}
}

#create output folder
dir.create("~/Downloads/test.com")


#generate file 
convertor(type=c("elisa","elispot","fcs_analyzed_result", "hai","neut_ab_titer", "gene_expression_files"),study=c('SDY269'),path="~/Downloads/test.com")
