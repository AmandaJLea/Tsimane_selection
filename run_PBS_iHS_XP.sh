###########
# allele freq by pop
########### 

cd /Genomics/ayroleslab2/alea/tsimane/selection_files

path_plink=/Genomics/grid/users/alea/programs/plink_1.90
path_out=/Genomics/ayroleslab2/alea/tsimane/selection_files
# Tsimane + 1000G data merge
path_data=/Genomics/ayroleslab2/alea/tsimane/selection_files/phased_haps_Jan2021/TSI1_MEG_filt3_w1000G_ALLchr.masked.merged.vcf.gz

# FST prep
grep 'Tsim' 26Jan21_unadmixed_1KG_TSI.txt > TSI_CHB.txt
grep 'CHB' 26Jan21_unadmixed_1KG_TSI.txt >> TSI_CHB.txt

grep 'Tsim' 26Jan21_unadmixed_1KG_TSI.txt > TSI_PEL.txt
grep 'PEL' 26Jan21_unadmixed_1KG_TSI.txt >> TSI_PEL.txt

grep 'CHB' 26Jan21_unadmixed_1KG_TSI.txt > CHB_PEL.txt
grep 'PEL' 26Jan21_unadmixed_1KG_TSI.txt >> CHB_PEL.txt

# FST calc
$path_plink --vcf $path_data --out TSI_CHB_FST --within TSI_CHB.txt --fst --geno 0.1
$path_plink --vcf $path_data --out TSI_PEL_FST --within TSI_PEL.txt --fst --geno 0.1
$path_plink --vcf $path_data --out CHB_PEL_FST --within CHB_PEL.txt --fst --geno 0.1

# allele freq
$path_plink --vcf $path_data --out TSI1_MEG_filt3_w1000G_ALLchr.masked.merged --within 26Jan21_unadmixed_1KG_TSI.txt --freq

###########
# calc PBS
########### 

library(data.table)

# all Tsimane, unadmixed PEL
x<-fread("TSI_CHB_FST.fst",header=T)
y<-fread("TSI_PEL_FST.fst",header=T)
z<-fread("CHB_PEL_FST.fst",header=T)

identical(x$SNP,y$SNP)
identical(z$SNP,y$SNP)

x$T1<- -log(1-x$FST) 
x$T2<- -log(1-y$FST) 
x$T3<- -log(1-z$FST) 

x$PBS<- (x$T2 + x$T1 - x$T3)/2

tmp<-fread('/Genomics/ayroleslab2/alea/tsimane/selection_files/TSI1_MEG_filt3_w1000G_ALLchr.masked.merged.frq.strat',sep=' ')
library(reshape2)
data_wide <- dcast(tmp, SNP ~ CLST, value.var="MAF")
data_wide=merge(data_wide,unique(tmp[,c(1,2,4,5)]),by='SNP')

data_wide2_tmp1<-apply(data_wide[,c(2:4)],1,function(x) length(which(x>0.01)))
data_wide2_tmp2<-apply(data_wide[,c(2:4)],1,function(x) length(which(x<0.99)))
data_wide2=data_wide[which(data_wide2_tmp1==3 & data_wide2_tmp2==3),]

both2<-merge(x[,c('SNP','PBS')],data_wide2,by='SNP')
write.table(both2,'26Jan21_PBS_results_post_RFmix.txt',row.names=F,sep='\t',quote=F)

###########
# calc iHS & XP-EHH
########### 

# split by chr and fill AA
path_data=/Genomics/ayroleslab2/alea/tsimane/selection_files/phased_haps_Jan2021/TSI1_MEG_filt3_w1000G_ALLchr.masked.merged.vcf.gz
tabix -p vcf $path_data

export PERL5LIB=$PERL5LIB:~/programs/vcftools/src/perl/
path_out=/Genomics/ayroleslab2/alea/tsimane/selection_files/phased_haps_Jan2021
module load vcftools
module load samtools

for chr in {1..22}; 
do bcftools view -r ${chr} -O z -o $path_out/temp_${chr}.vcf.gz $path_data; 
gzip -cd $path_out/temp_${chr}.vcf.gz | fill-aa -a /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/human_ancestor_GRCh37_e59/human_ancestor_${chr}.fa.gz | bgzip -c > $path_out/TSI1_MEG_filt3_w1000G_${chr}.masked.withAA.vcf.gz;
done

rm $path_out/temp_*.vcf.gz

##

library(rehh)
library(data.table)
library(vcfR)

info=read.delim('26Jan21_unadmixed_1KG_TSI.txt',header=F)
pel_unadmixed=subset(info,V3=='PEL')
chb=subset(info,V3=='CHB')
tsi_all=subset(info,V3=='Tsimane' )

for (i in 1:22) {
	hh <- data2haplohh(hap_file = paste("TSI1_MEG_filt3_w1000G_",i,".masked.withAA.vcf.gz",sep=''),polarize_vcf = TRUE,min_perc_geno.mrk = 50)
	
	# iHS - unadmixed tsimane
	tsi_all$hap1<-paste(tsi_all$V1,1,sep='_')
	tsi_all$hap2<-paste(tsi_all$V1,1,sep='_')
	hh_tsi = subset(hh, select.hap = c(tsi_all$hap1,tsi_all$hap2), min_maf = 0.05,min_perc_geno.mrk = 50)
	res<-scan_hh(hh_tsi,threads=4,phased = TRUE)

	# XP-EHH - unadmixed tsimane vs CHB
	chb$hap1<-paste(chb$V1,1,sep='_')
	chb$hap2<-paste(chb$V1,1,sep='_')
	
	keep<-c()
	for (f in 1:dim(chb)[1]) {
	keep<-c(keep,grep(chb$V1[f],hap.names(hh))) }
	
	hh_chb= subset(hh, select.hap = keep, min_maf = 0.05)
	res_chb<-scan_hh(hh_chb,threads=4,phased = TRUE)
	
	if (i == 1) {
    wgscan_xp <- res_chb
  } else {
    wgscan_xp <- rbind(wgscan_xp, res_chb)
  }
  
  	if (i == 1) {
    wgscan_ihs <- res
  } else {
    wgscan_ihs <- rbind(wgscan_ihs, res)
  }
print(i)}

# also add data on DAF and AAF

res.ihs<-ihh2ihs(wgscan_ihs,include_freq=TRUE)
res.ihs.df2<-as.data.frame(res.ihs$ihs)
	
xpehh2 <- ies2xpehh(scan_pop1 =  wgscan_ihs,
                           scan_pop2 =  wgscan_xp,
                           popname1 = "TSI",
                           popname2 = "CHB")
                     
all=merge(res.ihs.df2,xpehh2,by='POSITION',all.x=T,all.y=T)
write.table(all,paste("26Jan21_iHS_XP_post_RFmix.txt",sep=''),row.names=F,col.names=F,sep='\t',quote=F)

