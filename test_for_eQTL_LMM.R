#########
# genotype data
#########

path_out=/Genomics/ayroleslab2/alea/tsimane/transcriptome

for k in {2..22}; do
~/programs/plink2 --bfile $path_out/TOPMedImputed_clean_ALLchr_LDfilt_MOS --chr $k --recode A-transpose --out $path_out/TOPMedImputed_clean_chr${k}_LDfilt_MOS;
~/programs/plink2 --bfile $path_out/TOPMedImputed_clean_ALLchr_LDfilt_TSI --chr $k --recode A-transpose --out $path_out/TOPMedImputed_clean_chr${k}_LDfilt_TSI;
done

#########
# for cluster
#########

# for f in `cat chr_names.txt` ; do cat eQTL.R | sed -e s/CHROMNAME/$f/g > eQTL.$f.R; done
# for f in `cat chr_names.txt` ; do cat eQTL.sh | sed -e s/CHROMNAME/$f/g > eQTL.$f.sh; done
# rm commands1.sh; touch commands1.sh
# for f in `cat chr_names.txt` ; do echo 'sh eQTL.'$f'.sh' >> commands1.sh; done

#########
# load data
#########

library(data.table)
library(EMMREML)

setwd('/Genomics/ayroleslab2/alea/tsimane/transcriptome')
tot=read.delim('3Dec20_all_total_counts.txt',header=F)
align=read.delim('27Oct20_CT17_alignment_stats.txt')
meta=read.csv('27Oct20_CT17_rna_techdata_v2.csv',stringsAsFactors=F)
meta=merge(meta,align,by.x='Cole.ID',by.y='Sample')

# reorder and remove US
meta=merge(meta,tot,by.x='Cole.ID',by.y='V2')
meta=subset(meta,popn!='USA' & popn!='')
meta=meta[order(meta$Cole.ID),]
meta$popn2<-0
meta$popn2[which(meta$popn=='Tsimane')]<-1

meta2=subset(meta,sequencing_batch==1 & popn=='Tsimane')
meta1=subset(meta,sequencing_batch==2 & popn=='Moseten')

# read in more data
exp1A=read.delim('3Dec20_corrected_exp_Moseten_noSVA.txt')
exp2A=read.delim('3Dec20_corrected_exp_Tsimane_noSVA.txt')

kmatrix2=read.delim('output/TOPMedImputed_clean_ALLchr_TSI_rel.cXX.txt',header=F)
kmatrix1=read.delim('output/TOPMedImputed_clean_ALLchr_MOS_rel.cXX.txt',header=F)

geno2=fread('TOPMedImputed_clean_CHROMNAME_LDfilt_TSI.traw')
geno1=fread('TOPMedImputed_clean_CHROMNAME_LDfilt_MOS.traw')

genesnp1=fread('TOPMedImputed_clean_ALLchr_LDfilt_MOS.SNPs_genes',header=F)
genesnp2=fread('TOPMedImputed_clean_ALLchr_LDfilt_TSI.SNPs_genes',header=F)

genesnp1=subset(genesnp1,V1=='CHROMNAME')
genesnp2=subset(genesnp2,V1=='CHROMNAME')
genesnp1$V11<-paste(genesnp1$V1,genesnp1$V2,genesnp1$V4,genesnp1$V5,sep=':')
genesnp2$V11<-paste(genesnp2$V1,genesnp2$V2,genesnp2$V4,genesnp2$V5,sep=':')

# match files
# geno and kmatrix are matched
meta2$pid3<-paste(meta2$pid2,meta2$pid2,sep='_')
x<-match(names(geno2[,-c(1:6)]),meta2$pid3)
meta2_resort<-meta2[x,]
exp2A_resort<-exp2A[,x]
meta2_resort$age_visit[which(is.na(meta2_resort$age_visit))]<-median(meta2_resort$age_visit,na.rm=T)
Z_TSI<-diag(dim(meta2_resort)[1])
identical(meta2_resort$pid3,names(geno2[,-c(1:6)]))
identical(meta2_resort$Cole.ID,names(exp2A_resort))

meta1$pid3<-paste(meta1$pid2,meta1$pid2,sep='_')
x<-match(names(geno1[,-c(1:6)]),meta1$pid3)
meta1_resort<-meta1[x,]
exp1A_resort<-exp1A[,x]
meta1_resort$age_visit[which(is.na(meta1_resort$age_visit))]<-median(meta1_resort$age_visit,na.rm=T)
Z_MOS<-diag(dim(meta1_resort)[1])
identical(meta1_resort$pid3,names(geno1[,-c(1:6)]))
identical(meta1_resort$Cole.ID,names(exp1A_resort))

#################
# run eQTL pipeline
#################

res_tsimane=matrix(nrow=nrow(genesnp2),ncol=9)

for (i in 1:nrow(genesnp2)){
x<-which(geno2$SNP %in% genesnp2$V11[i])
y<-which(rownames(exp2A_resort) %in% genesnp2$V9[i])

meta2_resort$snp<-t(geno2[x,-c(1:6)])
design<-model.matrix(~meta2_resort$male + meta2_resort$age_visit + meta2_resort$snp)

emma=emmreml(y=t(exp2A_resort[y,]),X=design,Z=as.matrix(Z_TSI),K=as.matrix(kmatrix2),varbetahat=T,varuhat=T,PEVuhat=T,test=T)

res_tsimane[i,]=t(c(emma$beta[2:4],emma$pvalbeta[,"none"][2:4],emma$varbetahat[2:4]))
print(i) }

write.table(res_tsimane,'3Dec20_eQTL_TSI_CHROMNAME.txt',row.names=F,sep='\t')

##

res_moseten=matrix(nrow=nrow(genesnp1),ncol=9)

for (i in 1:nrow(genesnp1)){
x<-which(geno1$SNP %in% genesnp1$V11[i])
y<-which(rownames(exp1A_resort) %in% genesnp1$V9[i])

meta1_resort$snp<-t(geno1[x,-c(1:6)])
design<-model.matrix(~meta1_resort$male + meta1_resort$age_visit + meta1_resort$snp)

emma=emmreml(y=t(exp1A_resort[y,]),X=design,Z=as.matrix(Z_MOS),K=as.matrix(kmatrix1),varbetahat=T,varuhat=T,PEVuhat=T,test=T)

res_moseten[i,]=t(c(emma$beta[2:4],emma$pvalbeta[,"none"][2:4],emma$varbetahat[2:4]))
print(i) }

write.table(res_moseten,'3Dec20_eQTL_MOS_CHROMNAME.txt',row.names=F,sep='\t')

#################
# run eQTL pipeline - permute
#################

res_tsimane1=matrix(nrow=nrow(genesnp2),ncol=10)
res_tsimane2=matrix(nrow=nrow(genesnp2),ncol=10)
res_tsimane3=matrix(nrow=nrow(genesnp2),ncol=10)

for (k in 1:10){

geno2_perm<-geno2[,c(1:6,sample(7:dim(geno2)[2])),with=F]
for (i in 1:nrow(genesnp2)){
x<-which(geno2_perm$SNP %in% genesnp2$V11[i])
y<-which(rownames(exp2A_resort) %in% genesnp2$V9[i])

meta2_resort$snp<-t(geno2_perm[x,-c(1:6)])
design<-model.matrix(~meta2_resort$male + meta2_resort$age_visit + meta2_resort$snp)

emma=emmreml(y=t(exp2A_resort[y,]),X=design,Z=as.matrix(Z_TSI),K=as.matrix(kmatrix2),varbetahat=T,varuhat=T,PEVuhat=T,test=T)

res_tsimane1[i,k]=emma$beta[4]
res_tsimane2[i,k]=emma$pvalbeta[,"none"][4]
res_tsimane3[i,k]=emma$varbetahat[4] }
print(k) }
write.table(res_tsimane1,'3Dec20_eQTL_TSI_permbeta_CHROMNAME.txt',row.names=F,sep='\t')
write.table(res_tsimane2,'3Dec20_eQTL_TSI_permpval_CHROMNAME.txt',row.names=F,sep='\t')
write.table(res_tsimane3,'3Dec20_eQTL_TSI_permse_CHROMNAME.txt',row.names=F,sep='\t')

##

res_moseten1=matrix(nrow=nrow(genesnp1),ncol=10)
res_moseten2=matrix(nrow=nrow(genesnp1),ncol=10)
res_moseten3=matrix(nrow=nrow(genesnp1),ncol=10)

for (k in 1:10){

geno1_perm<-geno1[,c(1:6,sample(7:dim(geno1)[2])),with=F]
for (i in 1:nrow(genesnp1)){
x<-which(geno1_perm$SNP %in% genesnp1$V11[i])
y<-which(rownames(exp1A_resort) %in% genesnp1$V9[i])

meta1_resort$snp<-t(geno1_perm[x,-c(1:6)])
design<-model.matrix(~meta1_resort$male + meta1_resort$age_visit + meta1_resort$snp)

emma=emmreml(y=t(exp1A_resort[y,]),X=design,Z=as.matrix(Z_MOS),K=as.matrix(kmatrix1),varbetahat=T,varuhat=T,PEVuhat=T,test=T)

res_moseten1[i,k]=emma$beta[4]
res_moseten2[i,k]=emma$pvalbeta[,"none"][4]
res_moseten3[i,k]=emma$varbetahat[4] }
print(k) }
write.table(res_moseten1,'3Dec20_eQTL_MOS_permbeta_CHROMNAME.txt',row.names=F,sep='\t')
write.table(res_moseten2,'3Dec20_eQTL_MOS_permpval_CHROMNAME.txt',row.names=F,sep='\t')
write.table(res_moseten3,'3Dec20_eQTL_MOS_permse_CHROMNAME.txt',row.names=F,sep='\t')
