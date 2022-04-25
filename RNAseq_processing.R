setwd('/Genomics/ayroleslab2/alea/tsimane/transcriptome')
library(edgeR)
library(limma)

lengths=read.delim('/Genomics/ayroleslab2/alea/ref_genomes/hg38/hg38_PCgene_start_end.txt',header=F,stringsAsFactors=F)
data=read.delim('/Genomics/ayroleslab2/alea/tsimane/transcriptome/counts/28Oct20_all_ReadsPerGene.txt',header=F,stringsAsFactors=F)
names=read.delim('/Genomics/ayroleslab2/alea/tsimane/transcriptome/counts/28Oct20_sample_names.txt',header=F)
lengths$length<-lengths$V3-lengths$V2

data2=data[which(data$V1 %in% lengths$V4),]
remove<-c('ENSG00000244734',  'ENSG00000188536', 'ENSG00000206172' )

globin_counts<-data2[which(data2$V1 %in% remove),]
data2=data2[-which(data2$V1 %in% remove),]

tot_counts<-apply(data2[,-1],2,sum)
tot_counts2=matrix(nrow=dim(data2)[1],ncol=length(tot_counts))
for (i in 1:length(tot_counts)){tot_counts2[,i]<-tot_counts[i] }

data3=10^6 *(data2[,-1]/tot_counts2)
med_TPM<-apply(data3,1,median)

data4=data3[which(med_TPM>3),]
names(data4)<-names$V1
row.names(data4)<-data2$V1[which(med_TPM>3)]
write.table(data4,'3Dec20_all_ReadsPerGene_norm.txt',sep='\t',quote=F)

data4=data2[which(med_TPM>3),-1]
names(data4)<-names$V1
row.names(data4)<-data2$V1[which(med_TPM>3)]
write.table(data4,'3Dec20_all_ReadsPerGene_filt.txt',sep='\t',quote=F)

write.table(cbind(tot_counts,names(data4)),'3Dec20_all_total_counts.txt',sep='\t',col.names=F,row.names=F)

# globin prop
data2=data[which(data$V1 %in% lengths$V4),]
globin_counts2<-as.data.frame(t(globin_counts[,-1]))
globin_counts2$tot_counts<-colSums(data2[,-1])
globin_counts2$globin_counts<-apply(globin_counts2[,1:3],1,sum)
globin_counts2$globin_prop<-globin_counts2$globin_counts/(globin_counts2$tot_counts)
globin_counts2$id<-names$V1
write.table(globin_counts2,'3Dec20_globin_counts.txt',row.names=F,sep='\t')

#################
# PCA on Tsimane and Moseten together
#################

library(edgeR)
library(limma)
library(corrplot)
library(ggplot2)
library(ggpubr)

setwd('/Users/amandalea/Dropbox/Amanda_files/future_lab/Tsimane/current_pipeline_RNAseq/data')
data=read.delim('3Dec20_all_ReadsPerGene_filt.txt')
tot=read.delim('3Dec20_all_total_counts.txt',header=F)
meta=read.csv('27Oct20_CT17_rna_techdata_v2.csv',stringsAsFactors=F)
align=read.delim('27Oct20_CT17_alignment_stats.txt')
meta=merge(meta,align,by.x='Cole.ID',by.y='Sample')

# reorder and remove US
meta=merge(meta,tot,by.x='Cole.ID',by.y='V2')
meta=subset(meta,popn!='USA' & popn!='')
meta=meta[order(meta$Cole.ID),]
meta$popn2<-0
meta$popn2[which(meta$popn=='Tsimane')]<-1

data2=data[,meta$Cole.ID]
identical(as.character(meta$Cole.ID),as.character(names(data2)))

## Normalize using limma
voom_RNA <- voom(calcNormFactors(DGEList(counts=data2,lib.size = meta$V1)),plot=FALSE)

## PCA
pca<-prcomp(t(voom_RNA$E), scale = TRUE, center = TRUE)
# Variability of each principal component: pr.var
pr.var <- pca$sdev^2
# Variance explained by each principal component: pve
pve <- pr.var / sum(pr.var)
par(mfrow=c(1,1))
plot(1:20,pve[1:20],xlab='PC',ylab='PVE',bty='n')

meta$PC1<-pca$x[,1]
meta$PC2<-pca$x[,2]

p1<-ggplot(meta, aes(PC1, PC2, colour = factor(sequencing_batch))) + geom_point()+theme_bw(13)+ theme(legend.position="top")
p2<-ggplot(meta, aes(PC1, PC2, colour = popn)) + geom_point()+theme_bw(13)+ theme(legend.position="top")
ggarrange(p1, p2)

#################
# PCA on Tsimane and Moseten seperately
#################

meta1=subset(meta,sequencing_batch==1 & popn=='Tsimane')
meta2=subset(meta,sequencing_batch==2 & popn=='Moseten')

counts1=data[,meta1$Cole.ID]
identical(as.character(meta1$Cole.ID),as.character(names(counts1)))

counts2=data[,meta2$Cole.ID]
identical(as.character(meta2$Cole.ID),as.character(names(counts2)))

voom_RNA1 <- voom(calcNormFactors(DGEList(counts=counts1,lib.size = meta1$V1)),plot=FALSE)
voom_RNA2 <- voom(calcNormFactors(DGEList(counts=counts2,lib.size = meta2$V1)),plot=FALSE)

## PCA
pca1<-prcomp(t(voom_RNA1$E), scale = TRUE, center = TRUE)
pca2<-prcomp(t(voom_RNA2$E), scale = TRUE, center = TRUE)

## Use limma + SVA to remove technical effects 
library(sva)
meta1$age_visit[which(is.na(meta1$age_visit))]<-median(meta1$age_visit,na.rm=T)
meta2$age_visit[which(is.na(meta2$age_visit))]<-median(meta2$age_visit,na.rm=T)

mod = model.matrix(~male + age_visit ,data=meta1)
mod0 = model.matrix(~1, data=meta1)
sva1 = sva(voom_RNA1$E,mod,mod0,n.sv=10)

mod = model.matrix(~male + age_visit ,data=meta2)
mod0 = model.matrix(~1, data=meta2)
sva2 = sva(voom_RNA2$E,mod,mod0,n.sv=10)

# control for surrgogate variables plus major technical effects
mod1<-model.matrix(~sva1$sv )
fit <-lmFit(voom_RNA1,mod1)
fit <- eBayes(fit)
resid<-as.data.frame(residuals.MArrayLM(object=fit, voom_RNA1))
names(resid)<-names(counts1)
write.table(resid,'3Dec20_corrected_exp_Tsimane.txt',sep='\t')

mod1<-model.matrix(~sva2$sv  )
fit <-lmFit(voom_RNA2,mod1)
fit <- eBayes(fit)
resid<-as.data.frame(residuals.MArrayLM(object=fit, voom_RNA2))
names(resid)<-names(counts2)
write.table(resid,'3Dec20_corrected_exp_Moseten.txt',sep='\t')

## NO SVA

mod1<-model.matrix(~ meta1$Uniquely.mapped.. + meta1$RNA.Conc_ng.ul + meta1$V1 +meta1$X260.280.Ratio + meta1$Avg.aligned.length + meta1$Total.mapped..)
fit <-lmFit(voom_RNA1,mod1)
fit <- eBayes(fit)
resid<-as.data.frame(residuals.MArrayLM(object=fit, voom_RNA1))
names(resid)<-names(counts1)
write.table(resid,'3Dec20_corrected_exp_Tsimane_noSVA.txt',sep='\t')

mod1<-model.matrix(~ meta2$Uniquely.mapped.. + meta2$RNA.Conc_ng.ul + meta2$V1 +meta2$X260.280.Ratio + meta2$Avg.aligned.length + meta2$Total.mapped..)
fit <-lmFit(voom_RNA2,mod1)
fit <- eBayes(fit)
resid<-as.data.frame(residuals.MArrayLM(object=fit, voom_RNA2))
names(resid)<-names(counts2)
write.table(resid,'3Dec20_corrected_exp_Moseten_noSVA.txt',sep='\t')
