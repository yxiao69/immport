#organize  folder into study 
library(ImmuneSpaceR)
mainDir <- "~/Downloads"


sumdat=read_tsv("/Users/xiaoyi/Documents/2019spring/data mining/studies_by_data_types.tsv")
studys <- as.vector(sumdat$study)
names<-colnames(sumdat)[-c(1,7)]
for ( i in 1:length(studys)){
  for(j in 1:length(names)){
subDir <- paste0(names[j],'.out')
newDir<-studys[i]
#dir.create(file.path(mainDir, newDir))
setwd(file.path(mainDir, subDir))

f1=paste(studys[i],names[j],'.tsv',sep = '_')
  f2=paste(studys[i],names[j],'annotation.tsv',sep = '_')
  f3=paste(studys[i],names[j],'data_matrix.tsv',sep = '_')
 file.rename(paste0(studys[i],'_fcs_data_matrix.tsv'),to=f3)
file.copy(from=c(f1,f2,f3), to=file.path(mainDir, newDir))
  }
}

