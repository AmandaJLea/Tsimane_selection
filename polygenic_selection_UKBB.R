##########
# consider a window as associated with a trait if it included a SNP with a genome-wide significant association with this trait (p<10^-8)
###########

# cd /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/UKBB_pan/
library(data.table)
library(stringr)

files=read.delim('files.txt',header=F)
files2=str_split_fixed(files$V1, "-", 4)

for (i in 15:21){
x=paste("zcat ",files$V1[i],sep="")
y=fread(x)
sig=subset(y,pval_meta<10^-8 & low_confidence_MID=='FALSE')

sig$loc2<-sig$pos+1
write.table(sig[,c('chr','pos','loc2','pval_meta')],paste('trait-',files2[i,2],'-sig.bed',sep=''),row.names=F,sep='\t',col.names=F,quote=F)
print(i)}

#
windows=/Genomics/ayroleslab2/alea/ref_genomes/hg19/hg19_v0_Homo_sapiens_assembly19.50kb_windows.nooverlap_sort.txt
module load bedtools

for f in trait-301*sig.bed; do bedtools intersect -a $windows -b $f -u > windows.$f; done
for f in trait-302*sig.bed; do bedtools intersect -a $windows -b $f -u > windows.$f; done
for f in trait-300*sig.bed; do bedtools intersect -a $windows -b $f -u > windows.$f; done

##########
# combined selection scores into a Fisher’s score (FCS) equal to the sum, over the two statistics, of –log10(rank of the statistic for a given SNP/number of SNPs)
##########

setwd('~/Dropbox/Amanda_files/future_lab/Tsimane/current_pipeline_and_data_array')

keep=read.delim('26Jan21_merged_sites_keep.bed',header=F)
keep$SNP<-paste(keep$V1,keep$V2,sep=':')

# PBS results
pbs1=read.delim('Tsimane_only/2Feb21_PBS_results_post_RFmix.txt')
pbs1$SNP2<-paste("chr",pbs$SNP,sep='')
pbs_filt<-(subset(pbs1,Tsimane>0.01 & Tsimane<0.99 & CHB>0.01 & CHB<0.99 & PEL >0.01 & PEL<0.99 & SNP2 %in% keep$SNP))

# iHS/XP
iHS=read.delim('Tsimane_only/2Feb21_iHS_XP_post_RFmix.txt',header=T)
iHS$SNP<-paste("chr",paste(iHS$CHR.x,iHS$POSITION,sep=':'),sep='')
iHS_filt<-subset(iHS,SNP %in% keep$SNP)

# FCS
both=merge(iHS_filt,pbs_filt,by.x='SNP',by.y='SNP2')
both=both[complete.cases(both),]
both<-both[order(both$PBS,decreasing = TRUE),]
both$PBS_rank<-1:dim(both)[1]

both<-both[order(abs(both$IHS),decreasing = TRUE),]
both$iHS_rank<-1:dim(both)[1]

both<-both[order(abs(both$XPEHH_TSI_CHB),decreasing = TRUE),]
both$XP_rank<-1:dim(both)[1]

both$FCS<- (-log10(both$iHS_rank/dim(both)[1]) + -log10(both$PBS_rank/dim(both)[1]) + -log10(both$XP_rank/dim(both)[1]))
write.table(both,'Tsimane_only/2Feb21_FCS_by_SNP.txt',row.names=F,sep='\t',quote=F)

##########
# computed for each genomic window, associated or not with the trait, the average FCS, the proportion of conserved SNP positions based on GERP scores > 2, and the recombination rate using the combined HapMap genetic map, to account for the confounding effects of background selection.
##########

windows=/Genomics/ayroleslab2/alea/ref_genomes/hg19/hg19_v0_Homo_sapiens_assembly19.50kb_windows_sort.txt

rm /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/GERP_conserved_elements/hg19_all_elements.bed; touch /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/GERP_conserved_elements/hg19_all_elements.bed

for f in {1..22}; do awk -v CHR=$f '{OFS=""; print CHR,"\t",$1,"\t",$2}' /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/GERP_conserved_elements/hg19_chr${f}_elems.txt >> /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/GERP_conserved_elements/hg19_all_elements.bed; done

bedtools intersect -a /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/GERP_conserved_elements/hg19_all_elements.bed -b $windows -wo > /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/GERP_conserved_elements/hg19_50kb_windows_overlap.txt

rm /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/1000GP_Phase3/genetic_map_all_chr.txt; touch /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/1000GP_Phase3/genetic_map_all_chr.txt

for f in {1..22}; do awk -v CHR=$f '{OFS=""; print CHR,"\t",$1,"\t",$1+1,"\t",$2}' /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/1000GP_Phase3/genetic_map_chr${f}_combined_b37.txt | tail -n +2 >> /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/1000GP_Phase3/genetic_map_all_chr.txt; done

bedtools intersect -a /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/1000GP_Phase3/genetic_map_all_chr.txt -b $windows -wo > /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/1000GP_Phase3/hg19_50kb_windows_overlap.txt

## newly computed for Tsimane

windows=/Genomics/ayroleslab2/alea/ref_genomes/hg19/hg19_v0_Homo_sapiens_assembly19.50kb_windows.nooverlap_sort.txt

bedtools intersect -a /Genomics/ayroleslab2/alea/tsimane/selection_files/2Feb21_FCS_by_SNP_bed.txt -b $windows -wo > /Genomics/ayroleslab2/alea/tsimane/selection_files/2Feb21_FCS_50kb_windows_overlap.bed

##

library(data.table)
recomb=fread('/Genomics/ayroleslab2/alea/ref_genomes/public_datasets/1000GP_Phase3/hg19_50kb_windows_overlap.txt',header=F)
gerp=fread('/Genomics/ayroleslab2/alea/ref_genomes/public_datasets/GERP_conserved_elements/hg19_50kb_windows_overlap.txt',header=F)
fcs=fread('/Genomics/ayroleslab2/alea/tsimane/selection_files/2Feb21_FCS_50kb_windows_overlap.bed',header=F)

fcs$window<-paste(fcs$V5,fcs$V6,fcs$V7,sep='_')
gerp$window<-paste(gerp$V4,gerp$V5,gerp$V6,sep='_')
recomb$window<-paste(recomb$V5,recomb$V6,recomb$V7,sep='_')

# conserved bases by window
gerp_sum<-aggregate(gerp$V7~gerp$window,FUN=sum)
names(gerp_sum)<-c('window','conserved_bp')

# recomb by window
tmp1<-aggregate(recomb$V4~recomb$window,FUN=sum)
tmp2<-aggregate(recomb$V4~recomb$window,FUN=length)
identical(tmp1[,1],tmp2[,1])
recomb_sum<-as.data.frame(cbind(tmp1[,2],tmp2[,2]))
recomb_sum$window<-tmp1[,1]
recomb_sum$mean_recomb<- recomb_sum$V1 /50000

all=merge(recomb_sum,gerp_sum,by='window',all.x=T,all.y=T)
all=merge(all,fcs[,-c(3,5,6,7,8)],by='window',all.x=T,all.y=T)

names(all)<-c("window","recomb_sum","recomb_count","mean_recomb","conserved_bp","chr","loc","fcs") 

# remove overlapping windows
windows=read.delim('/Genomics/ayroleslab2/alea/ref_genomes/hg19/hg19_v0_Homo_sapiens_assembly19.50kb_windows.nooverlap_sort.txt',header=F)
windows$window<-paste(windows$V1,windows$V2,windows$V3,sep='_')
fcs=subset(all,window %in% windows$window)

write.table(fcs,'/Genomics/ayroleslab2/alea/tsimane/selection_files/2Feb21_combined_selection_FCS_wGERP_wrecomb.txt',row.names=F,sep='\t')

##########
# to test for polygenic selection, we generated a null distribution by randomly sampling x windows (x being the number of windows associated with a tested trait) among windows with a similar number of SNPs, proportion of GERP > 2 sites and recombination rate observed in the trait-associated windows. 
# We then calculated the average of the mean of the FCS across the x resampled windows. 
# We resampled 100,000 sets of x windows for each trait. 
##########

# cd /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/UKBB_pan
# for f in `cat files.txt` ; do cat temp.R | sed -e s/FILEINFO/${f}/g > temp.${f}.R; rm temp.${f}.sh; touch temp.${f}.sh; echo '#!/bin/bash' >> temp.${f}.sh; echo "module load R/3.5.2" >> temp.${f}.sh; echo "Rscript "temp.${f}.R >> temp.${f}.sh;done
# ls *temp.3*sh > commands1.sh

#!/bin/bash

windows=read.delim('/Genomics/ayroleslab2/alea/ref_genomes/hg19/hg19_v0_Homo_sapiens_assembly19.50kb_windows.nooverlap_sort.txt',header=F)
windows$window<-paste(windows$V1,windows$V2,windows$V3,sep='_')

fcs=read.delim('/Genomics/ayroleslab2/alea/tsimane/selection_files/2Feb21_combined_selection_FCS_wGERP_wrecomb.txt',header=T)

sig=read.delim(paste('windows.trait-',FILEINFO,'-sig.bed',sep=''),header=F)
sig$window<-paste(sig$V1,sig$V2,sig$V3,sep='_')
sig=subset(sig,window %in% windows$window)

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

sig_info<-uniq_window[which(uniq_window$window %in% sig$window),]
notsig_info<-uniq_window[-which(uniq_window$window %in% sig$window),]

fcs=unique(merge(uniq_window,fcs,by='window'))
sig_fcs=unique(subset(fcs,window %in% sig_info$window)[,c('chr','loc','fcs','recomb_quant','gerp_quant')])

sig_info2<-as.data.frame(table(sig_fcs$recomb_quant,sig_fcs$gerp_quant))

# null distribution
pval<-c()
median<-c()

for (k in 1:10000){

matched_fcs<-c()
for (i in 1:dim(sig_info2)[1]){
matched_windows<-(subset(notsig_info,recomb_quant==sig_info2$Var1[i] & gerp_quant==sig_info2$Var1[i] ))
matched_snps<-unique(subset(fcs,window %in% matched_windows$window)[,c('chr','loc','fcs','recomb_quant','gerp_quant')])
matched_fcs<-c(matched_fcs,sample(matched_snps$fcs,sig_info2$Freq[i]))
} 

print(k)
pval<-c(pval,wilcox.test(sig_fcs$fcs,matched_fcs)$p.value)
median<-c(median,median(matched_fcs,na.rm=T))
}

output<-as.data.frame(cbind(pval,median))
output$obs_median<-median(sig_fcs$fcs,na.rm=T)

write.table(output,'2Feb21_results_Tsimane_FILEINFO.txt',row.names=F,sep='\t')

##########
# To test for significance, we computed a resampling P-value by calculating the proportion of resampled windows which mean FCS was higher than that observed for the tested trait. All P-values for polygenic adaptation were then adjusted for multiple testing by the Benjamini-Hochberg method, to account for the number of traits tested, and traits with an adjusted p < 0.05 were considered as candidates for polygenic selection.
##########

# cd /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/UKBB_pan

library(stringr)
files2=read.delim('files.txt',header=F)

fcs=read.delim('/Genomics/ayroleslab2/alea/tsimane/selection_files/2Feb21_combined_selection_FCS_wGERP_wrecomb.txt',header=T)

pval<-c()
pval2<-c()
mean_null<-c()
mean_obs<-c()

for (i in c(1:42,44:49)){
x=paste("2Feb21_results_Tsimane_",files2[i,1],'.txt',sep="")
y=read.delim(x)

pval<-c(pval,length(which(y$pval<0.05)))
pval2<-c(pval2,length(which(y$median<y$obs_median)))
mean_null<-c(mean_null,median(y$median))
mean_obs<-c(mean_obs,median(y$obs_median))
}

output<-as.data.frame(cbind(pval,pval2,mean_null,mean_obs))
output$trait_id<-files2$V1[c(1:42,44:49)]

########
# plots of UKBB polygenic selection results
########

# cd /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/UKBB_pan

library(stringr)
files=read.delim('files.txt',header=F)
keep<-c(30000,30010,30020,30030,30040,30050,30060,30080,30090,30120,30130,30140,30150,30180,30190,30200,30210,30220,30240,30250,30710)
files2=subset(files,V1 %in% keep)

fcs=read.delim('/Genomics/ayroleslab2/alea/tsimane/selection_files/2Feb21_combined_selection_FCS_wGERP_wrecomb.txt',header=T)

pval<-c()
pval2<-c()
mean_null<-c()
mean_obs<-c()
output<-matrix(nrow=10000,ncol=21)

for (i in c(1:21)){
x=paste("2Feb21_results_Tsimane_",keep[i],'.txt',sep="")
y=read.delim(x)

pval2<-c(pval2,length(which(y$median<y$obs_median)))
mean_null<-c(mean_null,median(y$median))
mean_obs<-c(mean_obs,median(y$obs_median))
output[,i]<-y$median }

output_sum<-as.data.frame(cbind(pval2,mean_null,mean_obs))
output_sum$trait_id<-c('White blood cell (leukocyte) count','Red blood cell (erythrocyte) count','Haemoglobin concentration','Haematocrit percentage','Mean corpuscular volume','Mean corpuscular haemoglobin','Mean corpuscular haemoglobin concentration','Platelet count','Platelet crit','Lymphocyte count','Monocyte count','Neutrophill count','Eosinophill count','Lymphocyte percentage','Monocyte percentage','Neutrophill percentage','Eosinophil percentage','Basophill percentage','Reticulocyte percentage','Reticulocyte count','C-reactive protein')

##

pdf(file='temp.pdf')
par(mar = c(5.1, 20.1, 4.1, 2.1)) 
boxplot(output,outline=F,horizontal=T,ylim=c(1.05,1.4),yaxt = "n",xlab='Median evidence for selection')
axis(2, at=1:21, labels=output_sum$trait_id,las=1)
points(output_sum$mean_obs,1:21,col='red',pch=20)
dev.off()
