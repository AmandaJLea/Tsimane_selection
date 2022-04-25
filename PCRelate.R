
# https://www.bioconductor.org/packages/devel/bioc/vignettes/GENESIS/inst/doc/pcair.html
# GENESIS uses PC-AiR for population structure inference that is robust to known or cryptic relatedness, and it uses PC-Relate for accurate relatedness estimation in the presence of population structure, admixutre, and departures from Hardy-Weinberg equilibrium.
# The PC-Relate method is used to accurately estimate measures of recent genetic relatedness in samples with unknown or unspecified population structure without requiring reference population panels. PC-Relate uses ancestry representative principal components to account for sample ancestry differences and provide estimates that are robust to population structure, ancestry admixture, and departures from Hardy-Weinberg equilibirum.

########################
# run KING
########################

cd /Genomics/ayroleslab2/alea/tsimane/transcriptome/genotyped_reads

~/programs/plink2 --bfile /Genomics/ayroleslab2/alea/tsimane/int_data_files/TSI1_MEG_filt2 --maf 0.05 --max-alleles 2 --snps-only --out TSI1_MEG_filt2 --make-king-table --hwe 0.000001 

~/programs/plink2 --bfile /Genomics/ayroleslab2/alea/tsimane/int_data_files/TSI1_MEG_filt2 --maf 0.05 --max-alleles 2 --snps-only --out TSI1_MEG_filt2 --make-king square0 --hwe 0.000001 

########################
# PC-AiR - old sample names
########################

setwd('~/Desktop')
library(SNPRelate)
library(GENESIS)
library(GWASTools)

snpgdsBED2GDS(bed.fn = "TSI1_MEG_filt2.bed", 
              bim.fn = "TSI1_MEG_filt2.bim", 
              fam.fn = "TSI1_MEG_filt2.fam", 
              out.gdsfn = "TSI1_MEG_filt2.gds")

# LD pruning
gds <- snpgdsOpen('TSI1_MEG_filt2.gds')
snpset <- snpgdsLDpruning(gds, method="corr", slide.max.bp=10e6, 
                          ld.threshold=sqrt(0.1), verbose=FALSE)
pruned <- unlist(snpset, use.names=FALSE)
length(pruned)
snpgdsClose(gds)

geno <- GdsGenotypeReader(filename = "TSI1_MEG_filt2.gds")
genoData <- GenotypeData(geno)

# KING robust estimates
# king <- snpgdsIBDKING(gds)
# KINGmat <- kingToMatrix(king)
library(data.table)
KING<-fread('TSI1_MEG_filt2.kin0')     
id<-unique((c(KING$IID1,KING$IID2)))
KINGmat<-diag(ncol=length(id),nrow=length(id))
rownames(KINGmat)<-id;colnames(KINGmat)<-id
 
for (i in 1:dim(KING)[1]){
k<-which(colnames(KINGmat)==KING$IID1[i] )
j<-which(rownames(KINGmat)==KING$IID2[i] )
KINGmat[k,j ]<-KING$KINSHIP[i] 
KINGmat[j,k]<-KING$KINSHIP[i] }
write.table(KINGmat,'TSI1_MEG_filt2_KING_matrix.txt',sep='\t')

# related ids
KINGmat=read.delim('TSI1_MEG_filt2_KING_matrix.txt')
KINGmat2<-as.matrix(KINGmat)
rownames(KINGmat2)<-id;colnames(KINGmat2)<-id

part <- pcairPartition(kinobj = KINGmat2, divobj = KINGmat2)
mypcair <- pcair(genoData, kinobj = KINGmat2, divobj = KINGmat2, snp.include = pruned, unrel.set = part$unrels)
plot(mypcair)

########################
# PC-Relate - old sample names
########################

genoData2 <- GenotypeBlockIterator(genoData, snpInclude=pruned)

mypcrelate <- pcrelate(genoData2, pcs = mypcair$vectors[,1:2], 
                       training.set = mypcair$unrels)

smoothScatter(mypcrelate$kinBtwn$k0, mypcrelate$kinBtwn$kin, xlab="k0", ylab="kinship")
iids <- as.character(getScanID(genoData2))
matrix<-pcrelateToMatrix(mypcrelate, sample.include = iids, thresh = 2^(-11/2), scaleKin = 2)

write.table(as.matrix(matrix),'TSI1_MEG_filt2_PCRelate_matrix.txt',sep='\t')

matrix<-read.delim('TSI1_MEG_filt2_PCRelate_matrix.txt')
part <- pcairPartition(kinobj = as.matrix(matrix), divobj = as.matrix(matrix),kin.thresh=0.125)
write.table(as.matrix(part$unrels),'TSI1_MEG_filt2_PCRelate_unrelated.txt',sep='\t')

########################
# compare
########################

res1<-as.data.frame(mypcrelate$kinBtwn)
res1$pair<-paste(res1$ID1,res1$ID2,sep='_')

KING$pair<-paste(KING$IID2,KING$IID1,sep='_')
res2<-merge(res1,KING,by='pair')
smoothScatter(res2$KINSHIP,res2$kin, ylab="PC-Relate", xlab="KING")
x=c(0,1);y=c(0,1);abline(lm(y~x),lty=2)

res3<-subset(res2,KINSHIP>0.0625 | kin>0.0625)
plot(res3$KINSHIP,res3$kin, ylab="PC-Relate", xlab="KING")
x=c(0,1);y=c(0,1);abline(lm(y~x),lty=2)

