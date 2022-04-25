setwd("/Users/alea/Dropbox/Amanda_files/lea_lab/Tsimane/current_pipeline_and_data_array/Tsimane_only/demography")

########
# observed data
########

library(data.table)
library(qqman)
library(ggplot2)
library(tidyr)

genes=read.delim('~/Dropbox/Amanda_files/lea_lab/Tsimane/data/14Sep20_genes_near_50kb_windows.bed',header=F)
keep=read.delim('~/Dropbox/Amanda_files/lea_lab/Tsimane/current_pipeline_and_data_array/26Jan21_merged_sites_keep.bed',header=F)
keep$SNP<-paste(keep$V1,keep$V2,sep=':')

# results
pbs=read.delim('~/Dropbox/Amanda_files/lea_lab/Tsimane/current_pipeline_and_data_array/26Jan21_PBS_results_post_RFmix.txt')
pbs$SNP2<-paste("chr",paste(pbs$CHR,pbs$LOC,sep=':'),sep='')

iHS=read.delim('~/Dropbox/Amanda_files/lea_lab/Tsimane/current_pipeline_and_data_array/26Jan21_iHS_XP_post_RFmix.txt',header=F)
iHS$SNP<-paste("chr",paste(iHS$V2,iHS$V1,sep=':'),sep='')

iHS_filt<-subset(iHS,SNP %in% keep$SNP)
PBS_filt<-subset(pbs,SNP2 %in% keep$SNP)
both=merge(iHS_filt,PBS_filt,by.x='SNP',by.y='SNP2',all.x=T,all.y=T)

# FCS
both=both[complete.cases(both),]
both<-both[order(both$PBS,decreasing = TRUE),]
both$PBS_rank<-1:dim(both)[1]

both<-both[order(abs(both$V5),decreasing = TRUE),]
both$iHS_rank<-1:dim(both)[1]

both<-both[order(abs(both$V8),decreasing = TRUE),]
both$XP_rank<-1:dim(both)[1]

both$FCS<- (-log10(both$iHS_rank/dim(both)[1]) + -log10(both$PBS_rank/dim(both)[1]) + -log10(both$XP_rank/dim(both)[1]))

candidates=read.delim('11Nov21_candidate_regions_hg38.txt',header=F)
i=1
candidate_snps<-subset(both,CHR==candidates$V1[i] & LOC>candidates$V2[i] & LOC<candidates$V3[i])
candidate_snps$candidate<-i

for (i in 2:dim(candidates)[1]){
tmp<-subset(both,CHR==candidates$V1[i] & LOC>candidates$V2[i] & LOC<candidates$V3[i])
tmp$candidate<-i
candidate_snps=rbind(candidate_snps,tmp)
}

########
# simulated data
########

# PBS

i=1
pbs_all=read.delim(paste("8Nov21_simulation_PBS_results",i,".txt",sep=''))
pbs_all$iter<-i
	
	for (i in 2:100) {
	pbs=read.delim(paste("8Nov21_simulation_PBS_results",i,".txt",sep=''))
	pbs$iter<-i
	pbs_all<-rbind(pbs_all,pbs) }

# iHS and XP

i=1
ihs_all=read.delim(paste("8Nov21_simulation_iHS_XP_results",i,".txt",sep=''),header=F)
ihs_all$iter<-i
	
	for (i in 2:100) {
	ihs=read.delim(paste("8Nov21_simulation_iHS_XP_results",i,".txt",sep=''),header=F)
	ihs$iter<-i
	ihs_all<-rbind(ihs_all,ihs) }
	
ihs_all$SNP2<-paste(1,ihs_all$V1,ihs_all$iter,sep='_')
pbs_all$SNP2<-paste(pbs_all$SNP,pbs_all$iter,sep='_')

both=merge(ihs_all,pbs_all,by='SNP2')
#both=subset(both,V5>-10 & V8>-10 & PBS>-10)

# FCS
both=both[complete.cases(both),]
both<-both[order(both$PBS,decreasing = TRUE),]
both$PBS_rank<-1:dim(both)[1]

both<-both[order(abs(both$V5),decreasing = TRUE),]
both$iHS_rank<-1:dim(both)[1]

both<-both[order(abs(both$V8),decreasing = TRUE),]
both$XP_rank<-1:dim(both)[1]

both$FCS<- (-log10(both$iHS_rank/dim(both)[1]) + -log10(both$PBS_rank/dim(both)[1]) + -log10(both$XP_rank/dim(both)[1]))
	
plot(density(both$FCS),lwd=2,xlab='FCS: simulated data')
tmp<-as.data.frame(aggregate(candidate_snps$FCS~candidate_snps$candidate,FUN=mean))	
makeTransparent<-function(someColor, alpha=100)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
    blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}
points(tmp[,2],jitter(rep(.1,21),factor=8),pch=20,cex=2,col=makeTransparent('steelblue'))	

quantile<-c()
for (i in 1:21){
tmp2<-subset(both,FCS>tmp[i,2])
quantile<-c(quantile,dim(tmp2)[1]/dim(both)[1])
}