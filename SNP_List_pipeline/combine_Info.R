### combine the GTEx with linked SNP file

##### read in GTEx data
GTEx.skin.nonExp <- read.delim("../data/GTEx/Skin_Not_Sun_Exposed_Suprapubic.v8.signif_variant_gene_pairs_SNPID.txt",
                               header=F,stringsAsFactors=F,sep = " ")
GTEx.skin.Exp <- read.delim("../data/GTEx/Skin_Sun_Exposed_Lower_leg.v8.signif_variant_gene_pairs_SNPID.txt",
                            header=F,stringsAsFactors=F,sep = " ")
GTEx.all <- unique(rbind(GTEx.skin.nonExp,GTEx.skin.Exp))
dim(GTEx.all)
length(unique(GTEx.all[,7]))
GTEx.eGene <- unique(GTEx.all[,c(2:5,7,9)])

## get the symbol
library("org.Hs.eg.db") # remember to install it if you don't have it already
geneID <- unlist(sapply(GTEx.eGene[,6],function(x)strsplit(x,"\\.")[[1]][1]))
symbols <- mapIds(org.Hs.eg.db, keys = geneID, keytype = "ENSEMBL", column="SYMBOL")
GTEx.eGene <- data.frame(GTEx.eGene,symbols)
write.table(GTEx.eGene,file="GTEx_skin.txt",row.names = F,col.names = F,quote=F,sep='\t')

SNPID.filter <- read.delim("SNPID_filtered_all.txt",header=F,stringsAsFactors = F)
GTEx.eGene.filter <- GTEx.eGene[which(as.character(GTEx.eGene[,5]) %in% SNPID.filter[,1]),]

#library("biomaRt")
#ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
#geneSymbol <- getBM(attributes='hgnc_symbol', 
#      filters = 'ensembl_gene_id', 
#      values = as.character(GTEx.eGene[,5]), 
#      mart = ensembl)
#GTEx.eGene.filter[which(is.na(GTEx.eGene.filter[,7])),7] <- ""
GTEx.eGene.uniq <- aggregate(GTEx.eGene.filter[,5:7], list(GTEx.eGene.filter[,5]), paste, collapse=",")
GTEx.eGene.uniq <- GTEx.eGene.uniq[,c(1,3,4)]

##### read in the lindedin SNP file
fileName <- list.files("./GGR_filtered/",full.names = T)

for (i in 1:length(fileName)){
  linkedSNP <- read.delim(fileName[i],header=F,stringsAsFactors = F,sep='\t')
  linkedSNP_GTEx <- data.frame(linkedSNP,GTEx.eGene.uniq[match(as.character(linkedSNP[,3]),as.character(GTEx.eGene.uniq[,1])),])
  write.table(linkedSNP_GTEx,file=paste0(fileName[i],".txt"),col.names = F,row.names = F,quote=F,sep='\t')
}

###############################
## combine multiple files
## strict
file.strict <- list.files("./GGR_filtered/",full.names = T,pattern = "*strict*")[c(2,4,6)]
tmp <- read.delim(file.strict[1],header=F,stringsAsFactors = F)[,c(1,2,3,4,5,7,8)]
tmp <- data.frame(tmp,"ATAC")
colnames(tmp)[8] <- ""

tmp1 <- read.delim(file.strict[2],header=F,stringsAsFactors = F)[,c(1,2,3,4,5,7,8)]
tmp1 <- data.frame(tmp1,"H3K27ac")
colnames(tmp1)[8] <- ""

tmp2 <- read.delim(file.strict[3],header=F,stringsAsFactors = F)[,c(1,2,3,4,5,7,8)]
tmp2 <- data.frame(tmp2,"H3K4me1")
colnames(tmp2)[8] <- ""

tmp.all <- rbind(tmp,tmp1,tmp2)
tmp1 <- unique(tmp.all[,1:7])
tmp2 <- unique(tmp.all[,c(3,8)])
tmp2.aggregate <- aggregate(tmp2, list(tmp2[,1]), paste, collapse=",")
tmp2.aggregate <- tmp2.aggregate[,c(1,3)]
tmp.all1 <- data.frame(tmp1,tmp2.aggregate[match(as.character(tmp1[,3]),as.character(tmp2.aggregate[,1])),2])
index_linkedSNP <- read.delim("sig_LD_comb_uniq.txt",header=F,stringsAsFactors = F)
index_linkedSNP_uniq <- aggregate(index_linkedSNP, list(index_linkedSNP[,1,3]), max, na.rm=T)[,2:4] # linkedSNP, LD, indexSNP

tmp.all.index <- data.frame(tmp.all1,index_linkedSNP_uniq[match(as.character(tmp.all1[,3]),as.character(index_linkedSNP_uniq[,1])),])
indexSNP_info <- read.delim("../data/GWAS-skin-12112019.txt",stringsAsFactors = F)

i=1
tmp.all.index.info <- data.frame(tmp.all.index[i,],indexSNP_info[grep(as.character(tmp.all.index[i,11]),as.character(indexSNP_info[,1])),])

for (i in 2:nrow(tmp.all.index)){
  
  tmp <- data.frame(tmp.all.index[i,],indexSNP_info[grep(as.character(tmp.all.index[i,11]),as.character(indexSNP_info[,1])),])
  tmp.all.index.info <- rbind(tmp.all.index.info,tmp)
}

final.tab1 <- data.frame(LinkedSNP_IDs=as.character(tmp.all.index.info[,3]),
                         LinkedSNP_Position=paste0("chr",tmp.all.index.info[,1],":",tmp.all.index.info[,2]),
                         LinkedSNP_Ref_allele=tmp.all.index.info[,4],
                         LinkedSNP_Alt_allele=tmp.all.index.info[,5],
                         IndexSNP_ID=as.character(tmp.all.index.info[,11]),
                         IndexSNP_Position=paste0("chr",tmp.all.index.info[,13],":",tmp.all.index.info[,14]),
                         IndexSNP_Report_Gene=tmp.all.index.info[,15],
                         LD_score=tmp.all.index.info[,10],
                         Histone_Type=tmp.all.index.info[,8],
                         Disease=tmp.all.index.info[,19],
                         IndexSNP_Strongest_risk_allele=tmp.all.index.info[,16],
                         p_value=tmp.all.index.info[,17],
                         OR_or_Beta=tmp.all.index.info[,18],
                         eGene_symbol=tmp.all.index.info[,7],
                         eGene_ID=tmp.all.index.info[,6],check.names=F,stringsAsFactors = F)

final.tab1$LD_score[which(final.tab1$LinkedSNP_IDs==final.tab1$IndexSNP_ID)] <- NA

sum(final.tab1[,12]>5e-8,na.rm=T)                        
final.tab1.sig <- final.tab1[which(final.tab1[,12]<5e-8),]                         
write.table(final.tab1,file="Skin_diseases_strict_all.txt",quote=F,row.names = F,sep='\t',na = "")      
write.table(final.tab1.sig,file="Skin_diseases_strict_sig.txt",quote=F,row.names = F,sep='\t',na = "")      

#########################################
## 50bp
file.50bp <- list.files("./GGR_filtered/",full.names = T,pattern = "*50bp*")[c(2,4,6)]
tmp <- read.delim(file.50bp[1],header=F,stringsAsFactors = F)[,c(1,2,3,4,5,7,8)]
tmp <- data.frame(tmp,"ATAC")
colnames(tmp)[8] <- ""

tmp1 <- read.delim(file.50bp[2],header=F,stringsAsFactors = F)[,c(1,2,3,4,5,7,8)]
tmp1 <- data.frame(tmp1,"H3K27ac")
colnames(tmp1)[8] <- ""

tmp2 <- read.delim(file.50bp[3],header=F,stringsAsFactors = F)[,c(1,2,3,4,5,7,8)]
tmp2 <- data.frame(tmp2,"H3K4me1")
colnames(tmp2)[8] <- ""

tmp.all <- rbind(tmp,tmp1,tmp2)
tmp1 <- unique(tmp.all[,1:7])
tmp2 <- unique(tmp.all[,c(3,8)])
tmp2.aggregate <- aggregate(tmp2, list(tmp2[,1]), paste, collapse=",")
tmp2.aggregate <- tmp2.aggregate[,c(1,3)]
tmp.all1 <- data.frame(tmp1,tmp2.aggregate[match(as.character(tmp1[,3]),as.character(tmp2.aggregate[,1])),2])
index_linkedSNP <- read.delim("sig_LD_comb_uniq.txt",header=F,stringsAsFactors = F)
index_linkedSNP_uniq <- aggregate(index_linkedSNP, list(index_linkedSNP[,1,3]), max, na.rm=T)[,2:4] # linkedSNP, LD, indexSNP

tmp.all.index <- data.frame(tmp.all1,index_linkedSNP_uniq[match(as.character(tmp.all1[,3]),as.character(index_linkedSNP_uniq[,1])),])
indexSNP_info <- read.delim("../data/GWAS-skin-12112019.txt",stringsAsFactors = F)

i=1
tmp.all.index.info <- data.frame(tmp.all.index[i,],indexSNP_info[grep(as.character(tmp.all.index[i,11]),as.character(indexSNP_info[,1])),])

for (i in 2:nrow(tmp.all.index)){
  
  tmp <- data.frame(tmp.all.index[i,],indexSNP_info[grep(as.character(tmp.all.index[i,11]),as.character(indexSNP_info[,1])),])
  tmp.all.index.info <- rbind(tmp.all.index.info,tmp)
}

final.tab1 <- data.frame(LinkedSNP_IDs=as.character(tmp.all.index.info[,3]),
                         LinkedSNP_Position=paste0("chr",tmp.all.index.info[,1],":",tmp.all.index.info[,2]),
                         LinkedSNP_Ref_allele=tmp.all.index.info[,4],
                         LinkedSNP_Alt_allele=tmp.all.index.info[,5],
                         IndexSNP_ID=as.character(tmp.all.index.info[,11]),
                         IndexSNP_Position=paste0("chr",tmp.all.index.info[,13],":",tmp.all.index.info[,14]),
                         IndexSNP_Report_Gene=tmp.all.index.info[,15],
                         LD_score=tmp.all.index.info[,10],
                         Histone_Type=tmp.all.index.info[,8],
                         Disease=tmp.all.index.info[,19],
                         IndexSNP_Strongest_risk_allele=tmp.all.index.info[,16],
                         p_value=tmp.all.index.info[,17],
                         OR_or_Beta=tmp.all.index.info[,18],
                         eGene_symbol=tmp.all.index.info[,7],
                         eGene_ID=tmp.all.index.info[,6],check.names=F,stringsAsFactors = F)

final.tab1$LD_score[which(final.tab1$LinkedSNP_IDs==final.tab1$IndexSNP_ID)] <- NA

sum(final.tab1[,12]>5e-8,na.rm=T)                        
final.tab1.sig <- final.tab1[which(final.tab1[,12]<5e-8),]                         
write.table(final.tab1,file="Skin_diseases_50bp_all.txt",quote=F,row.names = F,sep='\t',na = "")      
write.table(final.tab1.sig,file="Skin_diseases_50bp_sig.txt",quote=F,row.names = F,sep='\t',na = "")      

########################################
## 100bp
file.100bp <- list.files("./GGR_filtered/",full.names = T,pattern = "*100bp*")[c(2,4,6)]
tmp <- read.delim(file.100bp[1],header=F,stringsAsFactors = F)[,c(1,2,3,4,5,7,8)]
tmp <- data.frame(tmp,"ATAC")
colnames(tmp)[8] <- ""

tmp1 <- read.delim(file.100bp[2],header=F,stringsAsFactors = F)[,c(1,2,3,4,5,7,8)]
tmp1 <- data.frame(tmp1,"H3K27ac")
colnames(tmp1)[8] <- ""

tmp2 <- read.delim(file.100bp[3],header=F,stringsAsFactors = F)[,c(1,2,3,4,5,7,8)]
tmp2 <- data.frame(tmp2,"H3K4me1")
colnames(tmp2)[8] <- ""

tmp.all <- rbind(tmp,tmp1,tmp2)
tmp1 <- unique(tmp.all[,1:7])
tmp2 <- unique(tmp.all[,c(3,8)])
tmp2.aggregate <- aggregate(tmp2, list(tmp2[,1]), paste, collapse=",")
tmp2.aggregate <- tmp2.aggregate[,c(1,3)]
tmp.all1 <- data.frame(tmp1,tmp2.aggregate[match(as.character(tmp1[,3]),as.character(tmp2.aggregate[,1])),2])
index_linkedSNP <- read.delim("sig_LD_comb_uniq.txt",header=F,stringsAsFactors = F)
index_linkedSNP_uniq <- aggregate(index_linkedSNP, list(index_linkedSNP[,1,3]), max, na.rm=T)[,2:4] # linkedSNP, LD, indexSNP

tmp.all.index <- data.frame(tmp.all1,index_linkedSNP_uniq[match(as.character(tmp.all1[,3]),as.character(index_linkedSNP_uniq[,1])),])
indexSNP_info <- read.delim("../data/GWAS-skin-12112019.txt",stringsAsFactors = F)

i=1
tmp.all.index.info <- data.frame(tmp.all.index[i,],indexSNP_info[grep(as.character(tmp.all.index[i,11]),as.character(indexSNP_info[,1])),])

for (i in 2:nrow(tmp.all.index)){
  
  tmp <- data.frame(tmp.all.index[i,],indexSNP_info[grep(as.character(tmp.all.index[i,11]),as.character(indexSNP_info[,1])),])
  tmp.all.index.info <- rbind(tmp.all.index.info,tmp)
}

final.tab1 <- data.frame(LinkedSNP_IDs=as.character(tmp.all.index.info[,3]),
                         LinkedSNP_Position=paste0("chr",tmp.all.index.info[,1],":",tmp.all.index.info[,2]),
                         LinkedSNP_Ref_allele=tmp.all.index.info[,4],
                         LinkedSNP_Alt_allele=tmp.all.index.info[,5],
                         IndexSNP_ID=as.character(tmp.all.index.info[,11]),
                         IndexSNP_Position=paste0("chr",tmp.all.index.info[,13],":",tmp.all.index.info[,14]),
                         IndexSNP_Report_Gene=tmp.all.index.info[,15],
                         LD_score=tmp.all.index.info[,10],
                         Histone_Type=tmp.all.index.info[,8],
                         Disease=tmp.all.index.info[,19],
                         IndexSNP_Strongest_risk_allele=tmp.all.index.info[,16],
                         p_value=tmp.all.index.info[,17],
                         OR_or_Beta=tmp.all.index.info[,18],
                         eGene_symbol=tmp.all.index.info[,7],
                         eGene_ID=tmp.all.index.info[,6],check.names=F,stringsAsFactors = F)

final.tab1$LD_score[which(final.tab1$LinkedSNP_IDs==final.tab1$IndexSNP_ID)] <- NA

sum(final.tab1[,12]>5e-8,na.rm=T)                        
final.tab1.sig <- final.tab1[which(final.tab1[,12]<5e-8),]                         
write.table(final.tab1,file="Skin_diseases_100bp_all.txt",quote=F,row.names = F,sep='\t',na = "")      
write.table(final.tab1.sig,file="Skin_diseases_100bp_sig.txt",quote=F,row.names = F,sep='\t',na = "")      

#######################################################
## 200bp
file.200bp <- list.files("./GGR_filtered/",full.names = T,pattern = "*200bp*")[c(2,4,6)]
tmp <- read.delim(file.200bp[1],header=F,stringsAsFactors = F)[,c(1,2,3,4,5,7,8)]
tmp <- data.frame(tmp,"ATAC")
colnames(tmp)[8] <- ""

tmp1 <- read.delim(file.200bp[2],header=F,stringsAsFactors = F)[,c(1,2,3,4,5,7,8)]
tmp1 <- data.frame(tmp1,"H3K27ac")
colnames(tmp1)[8] <- ""

tmp2 <- read.delim(file.200bp[3],header=F,stringsAsFactors = F)[,c(1,2,3,4,5,7,8)]
tmp2 <- data.frame(tmp2,"H3K4me1")
colnames(tmp2)[8] <- ""

tmp.all <- rbind(tmp,tmp1,tmp2)
tmp1 <- unique(tmp.all[,1:7])
tmp2 <- unique(tmp.all[,c(3,8)])
tmp2.aggregate <- aggregate(tmp2, list(tmp2[,1]), paste, collapse=",")
tmp2.aggregate <- tmp2.aggregate[,c(1,3)]
tmp.all1 <- data.frame(tmp1,tmp2.aggregate[match(as.character(tmp1[,3]),as.character(tmp2.aggregate[,1])),2])
index_linkedSNP <- read.delim("sig_LD_comb_uniq.txt",header=F,stringsAsFactors = F)
index_linkedSNP_uniq <- aggregate(index_linkedSNP, list(index_linkedSNP[,1,3]), max, na.rm=T)[,2:4] # linkedSNP, LD, indexSNP

tmp.all.index <- data.frame(tmp.all1,index_linkedSNP_uniq[match(as.character(tmp.all1[,3]),as.character(index_linkedSNP_uniq[,1])),])
indexSNP_info <- read.delim("../data/GWAS-skin-12112019.txt",stringsAsFactors = F)

i=1
tmp.all.index.info <- data.frame(tmp.all.index[i,],indexSNP_info[grep(as.character(tmp.all.index[i,11]),as.character(indexSNP_info[,1])),])

for (i in 2:nrow(tmp.all.index)){
  
  tmp <- data.frame(tmp.all.index[i,],indexSNP_info[grep(as.character(tmp.all.index[i,11]),as.character(indexSNP_info[,1])),])
  tmp.all.index.info <- rbind(tmp.all.index.info,tmp)
}

final.tab1 <- data.frame(LinkedSNP_IDs=as.character(tmp.all.index.info[,3]),
                         LinkedSNP_Position=paste0("chr",tmp.all.index.info[,1],":",tmp.all.index.info[,2]),
                         LinkedSNP_Ref_allele=tmp.all.index.info[,4],
                         LinkedSNP_Alt_allele=tmp.all.index.info[,5],
                         IndexSNP_ID=as.character(tmp.all.index.info[,11]),
                         IndexSNP_Position=paste0("chr",tmp.all.index.info[,13],":",tmp.all.index.info[,14]),
                         IndexSNP_Report_Gene=tmp.all.index.info[,15],
                         LD_score=tmp.all.index.info[,10],
                         Histone_Type=tmp.all.index.info[,8],
                         Disease=tmp.all.index.info[,19],
                         IndexSNP_Strongest_risk_allele=tmp.all.index.info[,16],
                         p_value=tmp.all.index.info[,17],
                         OR_or_Beta=tmp.all.index.info[,18],
                         eGene_symbol=tmp.all.index.info[,7],
                         eGene_ID=tmp.all.index.info[,6],check.names=F,stringsAsFactors = F)

final.tab1$LD_score[which(final.tab1$LinkedSNP_IDs==final.tab1$IndexSNP_ID)] <- NA

sum(final.tab1[,12]>5e-8,na.rm=T)                        
final.tab1.sig <- final.tab1[which(final.tab1[,12]<5e-8),]                         
write.table(final.tab1,file="Skin_diseases_200bp_all.txt",quote=F,row.names = F,sep='\t',na = "")      
write.table(final.tab1.sig,file="Skin_diseases_200bp_sig.txt",quote=F,row.names = F,sep='\t',na = "")      

#############################################
## 500bp
file.500bp <- list.files("./GGR_filtered/",full.names = T,pattern = "*500bp*")[c(2,4,6)]
tmp <- read.delim(file.500bp[1],header=F,stringsAsFactors = F)[,c(1,2,3,4,5,7,8)]
tmp <- data.frame(tmp,"ATAC")
colnames(tmp)[8] <- ""

tmp1 <- read.delim(file.500bp[2],header=F,stringsAsFactors = F)[,c(1,2,3,4,5,7,8)]
tmp1 <- data.frame(tmp1,"H3K27ac")
colnames(tmp1)[8] <- ""

tmp2 <- read.delim(file.500bp[3],header=F,stringsAsFactors = F)[,c(1,2,3,4,5,7,8)]
tmp2 <- data.frame(tmp2,"H3K4me1")
colnames(tmp2)[8] <- ""

tmp.all <- rbind(tmp,tmp1,tmp2)
tmp1 <- unique(tmp.all[,1:7])
tmp2 <- unique(tmp.all[,c(3,8)])
tmp2.aggregate <- aggregate(tmp2, list(tmp2[,1]), paste, collapse=",")
tmp2.aggregate <- tmp2.aggregate[,c(1,3)]
tmp.all1 <- data.frame(tmp1,tmp2.aggregate[match(as.character(tmp1[,3]),as.character(tmp2.aggregate[,1])),2])
index_linkedSNP <- read.delim("sig_LD_comb_uniq.txt",header=F,stringsAsFactors = F)
index_linkedSNP_uniq <- aggregate(index_linkedSNP, list(index_linkedSNP[,1,3]), max, na.rm=T)[,2:4] # linkedSNP, LD, indexSNP

tmp.all.index <- data.frame(tmp.all1,index_linkedSNP_uniq[match(as.character(tmp.all1[,3]),as.character(index_linkedSNP_uniq[,1])),])
indexSNP_info <- read.delim("../data/GWAS-skin-12112019.txt",stringsAsFactors = F)

i=1
tmp.all.index.info <- data.frame(tmp.all.index[i,],indexSNP_info[grep(as.character(tmp.all.index[i,11]),as.character(indexSNP_info[,1])),])

for (i in 2:nrow(tmp.all.index)){
  
  tmp <- data.frame(tmp.all.index[i,],indexSNP_info[grep(as.character(tmp.all.index[i,11]),as.character(indexSNP_info[,1])),])
  tmp.all.index.info <- rbind(tmp.all.index.info,tmp)
}

final.tab1 <- data.frame(LinkedSNP_IDs=as.character(tmp.all.index.info[,3]),
                         LinkedSNP_Position=paste0("chr",tmp.all.index.info[,1],":",tmp.all.index.info[,2]),
                         LinkedSNP_Ref_allele=tmp.all.index.info[,4],
                         LinkedSNP_Alt_allele=tmp.all.index.info[,5],
                         IndexSNP_ID=as.character(tmp.all.index.info[,11]),
                         IndexSNP_Position=paste0("chr",tmp.all.index.info[,13],":",tmp.all.index.info[,14]),
                         IndexSNP_Report_Gene=tmp.all.index.info[,15],
                         LD_score=tmp.all.index.info[,10],
                         Histone_Type=tmp.all.index.info[,8],
                         Disease=tmp.all.index.info[,19],
                         IndexSNP_Strongest_risk_allele=tmp.all.index.info[,16],
                         p_value=tmp.all.index.info[,17],
                         OR_or_Beta=tmp.all.index.info[,18],
                         eGene_symbol=tmp.all.index.info[,7],
                         eGene_ID=tmp.all.index.info[,6],check.names=F,stringsAsFactors = F)

final.tab1$LD_score[which(final.tab1$LinkedSNP_IDs==final.tab1$IndexSNP_ID)] <- NA

sum(final.tab1[,12]>5e-8,na.rm=T)                        
final.tab1.sig <- final.tab1[which(final.tab1[,12]<5e-8),]                         
write.table(final.tab1,file="Skin_diseases_500bp_all.txt",quote=F,row.names = F,sep='\t',na = "")      
write.table(final.tab1.sig,file="Skin_diseases_500bp_sig.txt",quote=F,row.names = F,sep='\t',na = "")      

##############################################
## 1kb
file.1000bp <- list.files("./GGR_filtered/",full.names = T,pattern = "*1000bp*")[c(2,4,6)]
tmp <- read.delim(file.1000bp[1],header=F,stringsAsFactors = F)[,c(1,2,3,4,5,7,8)]
tmp <- data.frame(tmp,"ATAC")
colnames(tmp)[8] <- ""

tmp1 <- read.delim(file.1000bp[2],header=F,stringsAsFactors = F)[,c(1,2,3,4,5,7,8)]
tmp1 <- data.frame(tmp1,"H3K27ac")
colnames(tmp1)[8] <- ""

tmp2 <- read.delim(file.1000bp[3],header=F,stringsAsFactors = F)[,c(1,2,3,4,5,7,8)]
tmp2 <- data.frame(tmp2,"H3K4me1")
colnames(tmp2)[8] <- ""

tmp.all <- rbind(tmp,tmp1,tmp2)
tmp1 <- unique(tmp.all[,1:7])
tmp2 <- unique(tmp.all[,c(3,8)])
tmp2.aggregate <- aggregate(tmp2, list(tmp2[,1]), paste, collapse=",")
tmp2.aggregate <- tmp2.aggregate[,c(1,3)]
tmp.all1 <- data.frame(tmp1,tmp2.aggregate[match(as.character(tmp1[,3]),as.character(tmp2.aggregate[,1])),2])
index_linkedSNP <- read.delim("sig_LD_comb_uniq.txt",header=F,stringsAsFactors = F)
index_linkedSNP_uniq <- aggregate(index_linkedSNP, list(index_linkedSNP[,1,3]), max, na.rm=T)[,2:4] # linkedSNP, LD, indexSNP

tmp.all.index <- data.frame(tmp.all1,index_linkedSNP_uniq[match(as.character(tmp.all1[,3]),as.character(index_linkedSNP_uniq[,1])),])
indexSNP_info <- read.delim("../data/GWAS-skin-12112019.txt",stringsAsFactors = F)

i=1
tmp.all.index.info <- data.frame(tmp.all.index[i,],indexSNP_info[grep(as.character(tmp.all.index[i,11]),as.character(indexSNP_info[,1])),])

for (i in 2:nrow(tmp.all.index)){
  
  tmp <- data.frame(tmp.all.index[i,],indexSNP_info[grep(as.character(tmp.all.index[i,11]),as.character(indexSNP_info[,1])),])
  tmp.all.index.info <- rbind(tmp.all.index.info,tmp)
}

final.tab1 <- data.frame(LinkedSNP_IDs=as.character(tmp.all.index.info[,3]),
                         LinkedSNP_Position=paste0("chr",tmp.all.index.info[,1],":",tmp.all.index.info[,2]),
                         LinkedSNP_Ref_allele=tmp.all.index.info[,4],
                         LinkedSNP_Alt_allele=tmp.all.index.info[,5],
                         IndexSNP_ID=as.character(tmp.all.index.info[,11]),
                         IndexSNP_Position=paste0("chr",tmp.all.index.info[,13],":",tmp.all.index.info[,14]),
                         IndexSNP_Report_Gene=tmp.all.index.info[,15],
                         LD_score=tmp.all.index.info[,10],
                         Histone_Type=tmp.all.index.info[,8],
                         Disease=tmp.all.index.info[,19],
                         IndexSNP_Strongest_risk_allele=tmp.all.index.info[,16],
                         p_value=tmp.all.index.info[,17],
                         OR_or_Beta=tmp.all.index.info[,18],
                         eGene_symbol=tmp.all.index.info[,7],
                         eGene_ID=tmp.all.index.info[,6],check.names=F,stringsAsFactors = F)

final.tab1$LD_score[which(final.tab1$LinkedSNP_IDs==final.tab1$IndexSNP_ID)] <- NA

sum(final.tab1[,12]>5e-8,na.rm=T)                        
final.tab1.sig <- final.tab1[which(final.tab1[,12]<5e-8),]                         
write.table(final.tab1,file="Skin_diseases_1000bp_all.txt",quote=F,row.names = F,sep='\t',na = "")      
write.table(final.tab1.sig,file="Skin_diseases_1000bp_sig.txt",quote=F,row.names = F,sep='\t',na = "")      






