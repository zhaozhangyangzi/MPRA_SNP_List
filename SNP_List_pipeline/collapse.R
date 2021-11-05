# this script is to collapse the significant file
sigFile <- list.files(pattern="*sig_LD_comb.txt")

sigData <- read.delim(sigFile,stringsAsFactors=F,header=F)

causalSNP <- sigData[,1]
tab.causalSNP <- table(causalSNP)

if (max(tab.causalSNP)==1){
     data.tmp <- data.frame(sigData[,1],sigData[,2:4],sigData[,2:4],stringsAsFactors=F)
     write.table(sigData,file=paste0('collapse-',sigFile),sep='\t',quote = F,col.names = F, row.names = F)
} else{
     collapseIndex <- which(tab.causalSNP>1)
	 sigData1 <- sigData[-which(sigData[,1] %in% names(tab.causalSNP)[collapseIndex]),]
     sigData1 <- data.frame(sigData1[,1],sigData1[,2:4],sigData1[,2:4],stringsAsFactors=F)
     for (i in 1:length(collapseIndex)){
		  causalSNP.tmp <- names(tab.causalSNP)[collapseIndex[i]]
		  data.tmp <- sigData[grep(causalSNP.tmp,sigData[,1]),]
		  
		  r2.list <- data.tmp[,2]
		  if (sum(!is.na(r2.list))!=0){		  
		    top.r2 <- max(data.tmp[,2],na.rm=T)
		    top.leadSNP <- unique(data.tmp[which.max(data.tmp[,2]),3])
		    top.cohort <- unique(data.tmp[which.max(data.tmp[,2]),4])
		  } else{
		    top.r2 <- unique(data.tmp[,2])
			top.leadSNP <- unique(data.tmp[,3])
			top.cohort <- data.tmp[1,4]
		  }
		  collapse.r2 <- paste(unique(data.tmp[,2]),collapse=',')
		  collapse.leadSNP <- paste(unique(data.tmp[,3]),collapse=',')
		  collapse.cohort <- paste(unique(data.tmp[,4]),collapse=',')
		  
		  collapse.tmp <- c(causalSNP.tmp,top.r2,top.leadSNP,top.cohort,collapse.r2,collapse.leadSNP,collapse.cohort)
		  sigData1 <- rbind(sigData1,collapse.tmp)
}
     write.table(sigData1,file=paste0('collapse-',sigFile),sep='\t',quote = F,col.names = F, row.names = F)
}

