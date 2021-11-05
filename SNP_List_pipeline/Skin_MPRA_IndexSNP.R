## this script is for skin MPRA proj

# get the SNP from GWAS catelog
setwd("/Users/12107/Dropbox/David/skin_MPRA/")
GWAScatelog <- read.delim("./Pub_Data/gwas_catalog_v1.0-associations_e96_r2019-11-21.tsv",check.names = F,stringsAsFactors = F)
head(GWAScatelog)

## disease
diseaseList <- c("Psoriasis","Psoriatic arthritis","Rosacea","Acne","Atopic Dermatitis",
                 "Lupus Erythematosus","Cutaneous Squamous Cell Carcinoma","Leprosy",
                 "Sarcoidosis","Keloid","Male pattern baldness","Androgenetic Alopecia",
                 "Alopecia Areata","Cutaneous Basal Cell Carcinoma",
                 "Cutaneous Squamous Cell Carcinoma","Melanoma","Vitiligo",
                 "Systemic Sclerosis","Stevens-Johnson and Toxic Epidermolysis Necrolysis",
                 "Cutaneous Inflammation","Behçet's Disease","Skin youthfulness ",
                 "skin pigmentation","Skin color")

GWAScatelog.Psoriasis <- GWAScatelog[grep("Psoriasis",GWAScatelog$`DISEASE/TRAIT`,ignore.case = TRUE),]
GWAScatelog.Psoriatic.arthritis <- GWAScatelog[grep("Psoriatic arthritis",GWAScatelog$`DISEASE/TRAIT`,ignore.case = TRUE),]
GWAScatelog.Rosacea <- GWAScatelog[grep("Rosacea",GWAScatelog$`DISEASE/TRAIT`,ignore.case = TRUE),]
GWAScatelog.Acne <- GWAScatelog[grep("Acne",GWAScatelog$`DISEASE/TRAIT`,ignore.case = TRUE),]
GWAScatelog.Acne <- GWAScatelog.Acne[-grep("Menstruation",GWAScatelog.Acne$`DISEASE/TRAIT`),]
GWAScatelog.Atopic.Dermatitis <- GWAScatelog[grep("Atopic dermatitis",GWAScatelog$`DISEASE/TRAIT`,ignore.case = TRUE),]
GWAScatelog.lupus.Erythematosus <- GWAScatelog[grep("Cutaneous lupus erythematosus",GWAScatelog$`DISEASE/TRAIT`,ignore.case = TRUE),]
GWAScatelog.Cutaneous.Squamous.Cell.Carcinoma <- GWAScatelog[grep("Cutaneous squamous cell carcinoma",GWAScatelog$`DISEASE/TRAIT`,ignore.case = TRUE),] 
#GWAScatelog.psoriasis <- GWAScatelog[grep("Cutaneous psoriasis",GWAScatelog$`DISEASE/TRAIT`,ignore.case = TRUE),] 
GWAScatelog.Leprosy <- GWAScatelog[grep("Leprosy",GWAScatelog$`DISEASE/TRAIT`,ignore.case = TRUE),] 
GWAScatelog.Sarcoidosis <- GWAScatelog[grep("Sarcoidosis",GWAScatelog$`DISEASE/TRAIT`,ignore.case = TRUE),] 
GWAScatelog.Keloid <- GWAScatelog[grep("Keloid",GWAScatelog$`DISEASE/TRAIT`,ignore.case = TRUE),] 
GWAScatelog.Keloid2 <- data.frame(matrix(NA,ncol=ncol(GWAScatelog.Keloid),nrow=2))
colnames(GWAScatelog.Keloid2) <- colnames(GWAScatelog.Keloid)
GWAScatelog.Keloid2$SNPS <- c("rs940187","rs1511412")
GWAScatelog.Keloid2$`DISEASE/TRAIT` <- "Keloid"
GWAScatelog.bald <- GWAScatelog[grep("Male-pattern baldness",GWAScatelog$`DISEASE/TRAIT`,ignore.case = TRUE),]
GWAScatelog.Androgenetic.Alopecia <- GWAScatelog[grep("Androgenetic Alopecia",GWAScatelog$`DISEASE/TRAIT`,ignore.case = TRUE),] 
GWAScatelog.Alopecia.Areata <- GWAScatelog[grep("Alopecia Areata",GWAScatelog$`DISEASE/TRAIT`,ignore.case = TRUE),]
GWAScatelog.Basal <- GWAScatelog[grep("Basal Cell Carcinoma",GWAScatelog$`DISEASE/TRAIT`,ignore.case = TRUE),] 
GWAScatelog.Melanoma <- GWAScatelog[grep("melanoma",GWAScatelog$`DISEASE/TRAIT`,ignore.case = TRUE),] 
GWAScatelog.Melanoma <- GWAScatelog.Melanoma[-which(GWAScatelog.Melanoma[,"DISEASE/TRAIT"]=="Non-melanoma skin cancer"),]
GWAScatelog.Vitiligo <- GWAScatelog[grep("Vitiligo",GWAScatelog$`DISEASE/TRAIT`,ignore.case = TRUE),] 
GWAScatelog.Systemic.Sclerosis <- GWAScatelog[grep("Systemic Sclerosis",GWAScatelog$`DISEASE/TRAIT`,ignore.case = TRUE),] 
GWAScatelog.Stevens <- GWAScatelog[grep("Stevens|Johnson|Toxic Epidermolysis Necrolysis",GWAScatelog$`DISEASE/TRAIT`,ignore.case = TRUE),] 
GWAScatelog.inflamation <- GWAScatelog[grep("Inflammatory skin disease",GWAScatelog$`DISEASE/TRAIT`,ignore.case = TRUE),]
GWAScatelog.Behcet <- GWAScatelog[grep("Behcet",GWAScatelog$`DISEASE/TRAIT`,ignore.case = TRUE),] 
GWAScatelog.youthfulness <-  data.frame(matrix(NA,ncol=ncol(GWAScatelog.Keloid),nrow=3))
colnames(GWAScatelog.youthfulness) <- colnames(GWAScatelog.Keloid)
GWAScatelog.youthfulness$SNPS <- c("rs6975107","rs318125","rs7616661")
GWAScatelog.youthfulness$`DISEASE/TRAIT` <- "Skin youthfulness"
#GWAScatelog.youthfulness <- GWAScatelog[grep("youthful",GWAScatelog$`DISEASE/TRAIT`,ignore.case = TRUE),] 
GWAScatelog.color <- GWAScatelog[grep("pigment",GWAScatelog$`DISEASE/TRAIT`,ignore.case = TRUE),] 

## SNP from paper
grep("28537254",GWAScatelog$PUBMEDID) # already in GWAS catelog
grep("30315195",GWAScatelog$PUBMEDID) # already in GWAS catelog
grep("30915103",GWAScatelog$PUBMEDID) # need to be included
grep("30908599",GWAScatelog$PUBMEDID) # need to be included

## 30915103
paper.30915103 <- read.csv("./Pub_Data/30915103/Table_2.csv",skip=5,stringsAsFactors = F,
                           header=F)
data.30915103 <- as.data.frame(matrix(NA,ncol=ncol(GWAScatelog.color),nrow=nrow(paper.30915103)))
colnames(data.30915103) <- colnames(GWAScatelog.color)
data.30915103$PUBMEDID <- "30915103"
data.30915103$`DISEASE/TRAIT` <- "Atopic Dermatitis"
data.30915103$`INITIAL SAMPLE SIZE` <- "1012cases and 1362controls"
data.30915103$`REPLICATION SAMPLE SIZE` <- "1634cases and 1263controls"
data.30915103$CHR_ID <- paper.30915103$V2
data.30915103$`REPORTED GENE(S)` <- paper.30915103$V3
data.30915103$SNPS <- paper.30915103$V1
data.30915103$`P-VALUE`<- paper.30915103$V20
data.30915103$`OR or BETA` <- paper.30915103$V21

## 30908599
## none


GWAS.SNP.all <- data.frame(rbind(GWAScatelog.Psoriasis[,c("SNPS","CHR_ID","CHR_POS","REPORTED GENE(S)","STRONGEST SNP-RISK ALLELE","P-VALUE","OR or BETA")],
                                 GWAScatelog.Psoriatic.arthritis[,c("SNPS","CHR_ID","CHR_POS","REPORTED GENE(S)","STRONGEST SNP-RISK ALLELE","P-VALUE","OR or BETA")],
                                 GWAScatelog.Rosacea[,c("SNPS","CHR_ID","CHR_POS","REPORTED GENE(S)","STRONGEST SNP-RISK ALLELE","P-VALUE","OR or BETA")],
                                 GWAScatelog.Acne[,c("SNPS","CHR_ID","CHR_POS","REPORTED GENE(S)","STRONGEST SNP-RISK ALLELE","P-VALUE","OR or BETA")],
                                 GWAScatelog.Atopic.Dermatitis[,c("SNPS","CHR_ID","CHR_POS","REPORTED GENE(S)","STRONGEST SNP-RISK ALLELE","P-VALUE","OR or BETA")],
                                 GWAScatelog.lupus.Erythematosus[,c("SNPS","CHR_ID","CHR_POS","REPORTED GENE(S)","STRONGEST SNP-RISK ALLELE","P-VALUE","OR or BETA")],
                                 GWAScatelog.Cutaneous.Squamous.Cell.Carcinoma[,c("SNPS","CHR_ID","CHR_POS","REPORTED GENE(S)","STRONGEST SNP-RISK ALLELE","P-VALUE","OR or BETA")],
                                 GWAScatelog.Leprosy[,c("SNPS","CHR_ID","CHR_POS","REPORTED GENE(S)","STRONGEST SNP-RISK ALLELE","P-VALUE","OR or BETA")],
                                 GWAScatelog.Sarcoidosis[,c("SNPS","CHR_ID","CHR_POS","REPORTED GENE(S)","STRONGEST SNP-RISK ALLELE","P-VALUE","OR or BETA")],
                                 GWAScatelog.Keloid[,c("SNPS","CHR_ID","CHR_POS","REPORTED GENE(S)","STRONGEST SNP-RISK ALLELE","P-VALUE","OR or BETA")],
                                 GWAScatelog.Keloid2[,c("SNPS","CHR_ID","CHR_POS","REPORTED GENE(S)","STRONGEST SNP-RISK ALLELE","P-VALUE","OR or BETA")],
                                 GWAScatelog.bald[,c("SNPS","CHR_ID","CHR_POS","REPORTED GENE(S)","STRONGEST SNP-RISK ALLELE","P-VALUE","OR or BETA")],
                                 GWAScatelog.Alopecia.Areata[,c("SNPS","CHR_ID","CHR_POS","REPORTED GENE(S)","STRONGEST SNP-RISK ALLELE","P-VALUE","OR or BETA")],
                                 GWAScatelog.Basal[,c("SNPS","CHR_ID","CHR_POS","REPORTED GENE(S)","STRONGEST SNP-RISK ALLELE","P-VALUE","OR or BETA")],
                                 GWAScatelog.Melanoma[,c("SNPS","CHR_ID","CHR_POS","REPORTED GENE(S)","STRONGEST SNP-RISK ALLELE","P-VALUE","OR or BETA")],
                                 GWAScatelog.Vitiligo[,c("SNPS","CHR_ID","CHR_POS","REPORTED GENE(S)","STRONGEST SNP-RISK ALLELE","P-VALUE","OR or BETA")],
                                 GWAScatelog.Systemic.Sclerosis[,c("SNPS","CHR_ID","CHR_POS","REPORTED GENE(S)","STRONGEST SNP-RISK ALLELE","P-VALUE","OR or BETA")],
                                 GWAScatelog.Stevens[,c("SNPS","CHR_ID","CHR_POS","REPORTED GENE(S)","STRONGEST SNP-RISK ALLELE","P-VALUE","OR or BETA")],
                                 GWAScatelog.inflamation[,c("SNPS","CHR_ID","CHR_POS","REPORTED GENE(S)","STRONGEST SNP-RISK ALLELE","P-VALUE","OR or BETA")],
                                 GWAScatelog.Behcet[,c("SNPS","CHR_ID","CHR_POS","REPORTED GENE(S)","STRONGEST SNP-RISK ALLELE","P-VALUE","OR or BETA")],
                                 GWAScatelog.youthfulness[,c("SNPS","CHR_ID","CHR_POS","REPORTED GENE(S)","STRONGEST SNP-RISK ALLELE","P-VALUE","OR or BETA")],
                                 GWAScatelog.color[,c("SNPS","CHR_ID","CHR_POS","REPORTED GENE(S)","STRONGEST SNP-RISK ALLELE","P-VALUE","OR or BETA")],
                                 data.30915103[,c("SNPS","CHR_ID","CHR_POS","REPORTED GENE(S)","STRONGEST SNP-RISK ALLELE","P-VALUE","OR or BETA")]),
                           Disease=c(rep("Psoriasis",nrow(GWAScatelog.Psoriasis)),
                                     rep("Psoriatic arthritis",nrow(GWAScatelog.Psoriatic.arthritis)),
                                     rep("Rosacea",nrow(GWAScatelog.Rosacea)),
                                     rep("Acne",nrow(GWAScatelog.Acne)),
                                     rep("Atopic Dermatitis",nrow(GWAScatelog.Atopic.Dermatitis)),
                                     rep("Lupus Erythematosus",nrow(GWAScatelog.lupus.Erythematosus)),
                                     rep("Cutaneous Squamous Cell Carcinoma",nrow(GWAScatelog.Cutaneous.Squamous.Cell.Carcinoma)),
                                     rep("Leprosy",nrow(GWAScatelog.Leprosy)),
                                     rep("Sarcoidosis",nrow(GWAScatelog.Sarcoidosis)),
                                     rep("Keloid",(nrow(GWAScatelog.Keloid)+2)),
                                     rep("Male pattern baldness",nrow(GWAScatelog.bald)),
                                     rep("Alopecia Areata",nrow(GWAScatelog.Alopecia.Areata)),
                                     rep("Cutaneous Basal Cell Carcinoma",nrow(GWAScatelog.Basal)),
                                     rep("Melanoma",nrow(GWAScatelog.Melanoma)),
                                     rep("Vitiligo",nrow(GWAScatelog.Vitiligo)),
                                     rep("Systemic Sclerosis",nrow(GWAScatelog.Systemic.Sclerosis)),
                                     rep("Stevens-Johnson and Toxic Epidermolysis Necrolysis",nrow(GWAScatelog.Stevens)),
                                     rep("Cutaneous Inflammation",nrow(GWAScatelog.inflamation)),
                                     rep("Behceta Disease",nrow(GWAScatelog.Behcet)),
                                     rep("Skin youthfulness",nrow(GWAScatelog.youthfulness)),
                                     rep("skin pigmentation",nrow(GWAScatelog.color)),
                                     rep("Atopic Dermatitis",nrow(data.30915103))))
length(unique(GWAS.SNP.all$SNP)) #2700
GWAS.SNP.all.rs <- GWAS.SNP.all[grep("rs",GWAS.SNP.all$SNPS),]
GWAS.SNP.all.rs <- GWAS.SNP.all.rs[-grep("x",GWAS.SNP.all$SNPS),]
length(unique(GWAS.SNP.all.rs$SNPS)) #2376
## from linux script: 2430 unique index SNPs
## fliter based on p value, filter based on GTEx data, filter based on ATAC in GGR data.
## start from eQTL then get SNP list 
## merge those two lists

GWAS.SNP.all.sig <- GWAS.SNP.all[which(GWAS.SNP.all$P.VALUE < (5e-8)),]
length(unique(GWAS.SNP.all.sig$SNPS)) # 2082
GWAS.SNP.all.rs.sig <- GWAS.SNP.all.rs[which(GWAS.SNP.all.rs$P.VALUE < (5e-8)),]
length(unique(GWAS.SNP.all.rs.sig$SNPS)) # 1895

#save.image(file="./RData/leadSNP.RData")

write.table(GWAS.SNP.all,file="./Results/GWAS-skin-12112019.txt",sep='\t',quote = F, row.names=F)
write.table(GWAS.SNP.all.sig,file="./Results/GWAS-skin-sig-12112019.txt",sep='\t',quote = F, row.names=F)

write.table(GWAS.SNP.all.rs,file="./Results/GWAS-skin-with_rsID-12112019.txt",sep='\t',quote = F, row.names=F)
write.table(GWAS.SNP.all.rs.sig,file="./Results/GWAS-skin-with_rsID-sig-12112019.txt",sep='\t',quote = F, row.names=F)

leadGWAS.all <- read.delim("GWAS-skin-with_rsID-12112019.txt",check.names=F,stringsAsFactors = F)
leadGWAS_SNPID <- c(leadGWAS.all$SNPS,"rs11923593","rs763035")
write.table(leadGWAS_SNPID,file="leadSNP-all-SNPID.txt",quote = F, row.names=F,col.names=F)



