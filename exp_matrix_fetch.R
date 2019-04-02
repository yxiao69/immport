## extract expression matrix
#organize  folder into study 
#library(ImmuneSpaceR)
#mainDir <- "~/Downloads"
#sumdat=read_tsv("/Users/xiaoyi/Documents/2019spring/data mining/studies_by_data_types.tsv")
#studys <- as.vector(sumdat$study)

library(Biobase)
library(readr)
#studys <- as.vector(sumdat$study)[12:20]
expression_fetch <- function(study="SDY80") {
  tmp <- CreateConnection(study)
  fmat=tmp$listGEMatrices()
  if (!is.null(fmat)){
    flist=as.vector(fmat[,name])
    
    #mainDir <- "~/Downloads/immport_result/"
    #subDir <- study
    #dir.create(file.path(mainDir, subDir))
    #setwd(file.path(mainDir, subDir))

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

#expre<-sapply(studys,function(x) expression_fetch(x))