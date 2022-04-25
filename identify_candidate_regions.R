###########
# sites to keep
########### 

mask=/Genomics/ayroleslab2/alea/ref_genomes/public_datasets/20141020.strict_mask.whole_genome.bed
zcat /Genomics/ayroleslab2/alea/tsimane/selection_files/phased_haps_Jan2021/TSI1_MEG_filt3_w1000G_ALLchr.masked.merged.vcf.gz | awk '{OFS="\t"; print "chr"$1,$2,$2+1}' | grep -v "#" > /Genomics/ayroleslab2/alea/tsimane/selection_files/26Jan21_merged_sites.bed

bedtools intersect -a /Genomics/ayroleslab2/alea/tsimane/selection_files/26Jan21_merged_sites.bed -b $mask -u > /Genomics/ayroleslab2/alea/tsimane/selection_files/26Jan21_merged_sites_keep.bed

###########
# merge and plot
########### 

library(data.table)
library(qqman)
library(ggplot2)
library(tidyr)

genes=read.delim('~/Dropbox/Amanda_files/future_lab/Tsimane/data/14Sep20_genes_near_50kb_windows.bed',header=F)
keep=read.delim('26Jan21_merged_sites_keep.bed',header=F)
keep$SNP<-paste(keep$V1,keep$V2,sep=':')

# results
pbs=read.delim('26Jan21_PBS_results_post_RFmix.txt')
pbs$SNP2<-paste("chr",paste(pbs$CHR,pbs$LOC,sep=':'),sep='')

iHS=read.delim('26Jan21_iHS_XP_post_RFmix.txt',header=F)
iHS$SNP<-paste("chr",paste(iHS$V2,iHS$V1,sep=':'),sep='')

iHS_filt<-subset(iHS,SNP %in% keep$SNP)
PBS_filt<-subset(pbs,SNP2 %in% keep$SNP)

# manhattan plots
ggplot(iHS_filt, aes(x=V1, y=V6)) +facet_wrap(~factor(V2),scales='free_x')+theme_bw()+geom_point(alpha=0.5)+ylab('iHS p-value')
ggplot(iHS_filt, aes(x=V1, y=V9)) +facet_wrap(~factor(V2),scales='free_x')+theme_bw()+geom_point(alpha=0.5)+ylab('XP-EHH p-value')
ggplot(PBS_filt, aes(x=LOC, y=PBS)) +facet_wrap(~factor(CHR),scales='free_x')+theme_bw()+geom_point(alpha=0.5)+ylab('PBS')

both=merge(iHS_filt,PBS_filt,by.x='SNP',by.y='SNP2',all.x=T,all.y=T)

# windows
windows<-unique(genes[,1:3])
windows$n_PBS<-NA
windows$outliers_PBS<-NA
windows$n_iHS<-NA
windows$outliers_iHS<-NA
windows$n_XP<-NA
windows$outliers_XP<-NA

thresh_pbs<-quantile(both$PBS,seq(0,1,0.01),na.rm=T)[100]
thresh_iHS<-quantile(abs(both$V5),seq(0,1,0.01),na.rm=T)[100]
thresh_xp1<-quantile(abs(both$V8),seq(0,1,0.01),na.rm=T)[100]

for (i in 1:dim(windows)[1]){
	tmp1<-windows[i,]
	tmp2<-subset(both,CHR== tmp1$V1 & LOC>= tmp1$V2 & LOC<tmp1$V3 )
	windows$n_PBS[i]<-length(which(tmp2$PBS > -100))
	windows$outliers_PBS[i]<-length(which(tmp2$PBS > thresh_pbs))

	windows$n_iHS[i]<-length(which(tmp2$V5 > -100))
	windows$n_XP[i]<-length(which(tmp2$V8 > -100))

	windows$outliers_iHS[i]<-length(which( abs(tmp2$V5) > thresh_iHS))
	windows$outliers_XP[i]<-length(which( abs(tmp2$V8) > thresh_xp1))
print(i) }

data2=subset(windows,n_PBS>0 )
thresh_pbs<-quantile(data2$outliers_PBS,seq(0,1,0.01),na.rm=T)[100]
data2=subset(windows,n_iHS>0 )
thresh_iHS<-quantile(abs(data2$outliers_iHS),seq(0,1,0.01),na.rm=T)[100]
data2=subset(windows,n_XP>0 )
thresh_xp1<-quantile(abs(data2$outliers_XP),seq(0,1,0.01),na.rm=T)[100]

# FCS
both=both[complete.cases(both),]
both<-both[order(both$PBS,decreasing = TRUE),]
both$PBS_rank<-1:dim(both)[1]

both<-both[order(abs(both$V5),decreasing = TRUE),]
both$iHS_rank<-1:dim(both)[1]

both<-both[order(abs(both$V8),decreasing = TRUE),]
both$XP_rank<-1:dim(both)[1]

both$FCS<- (-log10(both$iHS_rank/dim(both)[1]) + -log10(both$PBS_rank/dim(both)[1]) + -log10(both$XP_rank/dim(both)[1]))

ggplot(both, aes(x=V1, y=FCS)) +facet_wrap(~factor(V2),scales='free_x')+theme_bw()+geom_point(alpha=0.5)+ylab('FCS value')+xlab('Position')

windows$median_FCS<-NA
for (i in 1:dim(windows)[1]){
	tmp1<-windows[i,]
	tmp2<-subset(both,CHR== tmp1$V1 & LOC>= tmp1$V2 & LOC<tmp1$V3 )
	windows$median_FCS[i]<-median(tmp2$FCS,na.rm=T)} 
write.table(windows,'26Jan21_outliers_post_RFmix.txt',row.names=F,sep='\t')

thresh<-quantile(windows$median_FCS,seq(0,1,0.0001),na.rm=T)[9999]
candb<-subset(windows,median_FCS>=thresh)
candb$region<-paste(candb$V1,candb$V2,candb$V3,sep='_')

cand<-subset(windows,outliers_PBS>=thresh_pbs & outliers_iHS>=thresh_iHS & outliers_XP>=thresh_xp1)
cand$region<-paste(cand$V1,cand$V2,cand$V3,sep='_')

candb$criteria<-'FCS'
cand$criteria<-'3_stats'
all_cand<-rbind(candb,cand)

genes$region<-paste(genes$V1,genes$V2,genes$V3,sep='_')
candc<-merge(all_cand,genes[,c('V7','V8','region')],by='region')
write.table(candc,'26Jan21_candidate_regions.txt',row.names=F,sep='\t')

windows$region<-paste(windows$V1,windows$V2,windows$V3,sep='_')
genes_tested<-genes[which(genes$region %in% windows$region),]

