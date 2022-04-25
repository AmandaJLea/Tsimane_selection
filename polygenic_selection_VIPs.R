##################
# overlap virus interacting proteins (VIP) with Tsimane FCS measures of per-SNP selection
# then compare VIP SNPs to randomly sampled SNPs matched for various properties
##################

#!/bin/bash

vips=read.delim('29Nov20_VIPs.txt')

# remove overlapping windows
windows_nooverlap=read.delim('/Genomics/ayroleslab2/alea/ref_genomes/hg19/hg19_v0_Homo_sapiens_assembly19.50kb_windows.nooverlap_sort.txt',header=F)
windows_nooverlap$window<-paste(windows_nooverlap$V1,windows_nooverlap$V2,windows_nooverlap$V3,sep='_')

# windows in or near genes - non overlapping
windows_gene=read.delim('/Genomics/ayroleslab2/alea/ref_genomes/hg19/hg19_50kb_windows_closest_gene.txt',header=F)
windows_gene2<-subset(windows_gene,V13<10000)
windows_gene2$window<-paste(windows_gene2$V1,windows_gene2$V2,windows_gene2$V3,sep='_')

# FCS by SNP - non overlapping
windows_data=read.delim('/Genomics/ayroleslab2/alea/tsimane/selection_files/2Feb21_FCS_50kb_windows_overlap.bed',header=F)
windows_data$window=paste(windows_data$V1,windows_data$V6 ,windows_data$V7,sep='_')
windows_data<-windows_data[which(windows_data$window %in% windows_nooverlap$window),]

fcs=read.delim('/Genomics/ayroleslab2/alea/tsimane/selection_files/2Feb21_combined_selection_FCS_wGERP_wrecomb.txt',header=T)

# quantiles by window
uniq_window<-unique(fcs[,c('window','mean_recomb','conserved_bp')])

tmp<-quantile(uniq_window$mean_recomb,na.rm=T)
uniq_window$recomb_quant<-NA
uniq_window$recomb_quant[which(uniq_window$mean_recomb <= tmp[5])]<-4
uniq_window$recomb_quant[which(uniq_window$mean_recomb <= tmp[4])]<-3
uniq_window$recomb_quant[which(uniq_window$mean_recomb <= tmp[3])]<-2
uniq_window$recomb_quant[which(uniq_window$mean_recomb <= tmp[2])]<-1

tmp<-quantile(uniq_window$conserved_bp,na.rm=T)
uniq_window$gerp_quant<-NA
uniq_window$gerp_quant[which(uniq_window$conserved_bp <= tmp[5])]<-4
uniq_window$gerp_quant[which(uniq_window$conserved_bp <= tmp[4])]<-3
uniq_window$gerp_quant[which(uniq_window$conserved_bp <= tmp[3])]<-2
uniq_window$gerp_quant[which(uniq_window$conserved_bp <= tmp[2])]<-1

fcs2=unique(merge(uniq_window,fcs,by='window'))

results<-matrix(nrow=1000,ncol=dim(vips)[2])
results2<-matrix(nrow=1000,ncol=dim(vips)[2])
nsnps<-c();obs<-c()

for (z in 1:dim(vips)[2]){
genes<-unique(vips[,z]) 
windows_gene_tmp<-subset(windows_gene2,V7 %in% genes)

sig_info<-uniq_window[which(uniq_window$window %in% windows_gene_tmp$window),]
notsig_info<-uniq_window[-which(uniq_window$window %in% windows_gene_tmp$window),]

sig_fcs=unique(subset(fcs2,window %in% sig_info$window)[,c('chr','loc','fcs','recomb_quant','gerp_quant')])

sig_info2<-as.data.frame(table(sig_fcs$recomb_quant,sig_fcs$gerp_quant))

pval<-c()
median<-c()

for (k in 1:1000){
matched_fcs<-c()
for (i in 1:dim(sig_info2)[1]){
matched_windows<-(subset(notsig_info,recomb_quant==sig_info2$Var1[i] & gerp_quant==sig_info2$Var1[i] ))
matched_snps<-unique(subset(fcs2,window %in% matched_windows$window)[,c('chr','loc','fcs','recomb_quant','gerp_quant')])
matched_fcs<-c(matched_fcs,sample(matched_snps$fcs,sig_info2$Freq[i]))
} 

pval<-c(pval,wilcox.test(sig_fcs$fcs,matched_fcs)$p.value)
median<-c(median,median(matched_fcs,na.rm=T))
}

results[,z]<-pval
results2[,z]<-median
nsnps<-c(nsnps,dim(sig_fcs)[1])
obs<-c(obs,median(sig_fcs$fcs,na.rm=T))

print(z) }

vip_summary<-as.data.frame(names(vips))
vip_summary$nsnps<-nsnps
vip_summary$obs<-obs
vip_summary$null<-as.matrix(apply(results2,2,median))
vip_summary$wilcoxon<-apply(results,2,function(x) length(which(x<0.05)))
vip_summary$median<-NA

for (i in 1:dim(vips)[2]){
vip_summary$median[i]<-length(which(results2[,i]<obs[i]))}

write.table(vip_summary,'/Genomics/ayroleslab2/alea/tsimane/selection_files/2Feb21_VIP_summarized_results.txt',row.names=F,sep='\t')
write.table(results,'/Genomics/ayroleslab2/alea/tsimane/selection_files/2Feb21_VIP_null_results.txt',row.names=F,sep='\t')
write.table(results2,'/Genomics/ayroleslab2/alea/tsimane/selection_files/2Feb21_VIP_null_results2.txt',row.names=F,sep='\t')

########
# plots VIP results
########

setwd("/Users/amandalea/Dropbox/Amanda_files/future_lab/Tsimane/current_pipeline_and_data_array/Tsimane_only")
results=read.delim('2Feb21_VIP_summarized_results.txt')
results2=read.delim('2Feb21_VIP_null_results2.txt')
results$obs[which(results$Abbreviation=='HIV')]<-1.13518683

names(results2)<-results$Abbreviation
order<-results$names2[order(results$obs)]

par(mar = c(5.1, 20.1, 4.1, 2.1)) 
boxplot(results2[,order],outline=F,horizontal=T,ylim=c(0.99,1.35),yaxt = "n",xlab='Median evidence for selection')
axis(2, at=1:16, labels=results$names2[order(results$obs)],las=1)
points(results$obs[order(results$obs)],1:16,col='red',pch=20)
