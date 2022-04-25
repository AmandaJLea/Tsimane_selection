#############
# Get the joint SFS
#############

#To run fastimcoal2 you'll first need a SFS for the population(s) included in your model. because I was trying to simulate PBS stats I used a joint SFS between my study population, Peruvians (PEL) and Han Chinese (CHB) from the 1000G project.

#I used a program called easySFS for this (https://github.com/isaacovercast/easySFS) back in 2018, but there may be a simpler tool now

# VCF order = CHB, PEL, TSI

module load python/current
path_out=/Genomics/ayroleslab2/alea/tsimane/int_data_files

path_vcf=$path_out/TSI1_MEG_filt3_w1000G_4demog_neut3.vcf
path_pop=$path_out/TSI1_MEG_filt3_w1000G_4demog_neut3.pop.txt

# get folded SFS
~/programs/easySFS.py -i $path_vcf -p $path_pop -o $path_out/demography --proj=206,60,382 --prefix CHB_PEL_TSI -a -f

~/programs/easySFS.py -i $path_vcf -p $path_pop -o $path_out/demography_proj --proj=60,60,60 --prefix CHB_PEL_TSI_proj -a -f

###############

# plot the SFS - projected down to same # of IDs

data1=read.delim('TSI_MAFpop0.txt')
plot( c(2:dim(data1)[2])-1,(t(data1[1,-1])),pch=20,xlab='Frequency class of allele',ylab='Number of sites',type='l',xlim=c(0,10),main='Folded SFS',ylim=c(0,25000))

data1=read.delim('PEL_MAFpop0.txt')
lines(c(2:dim(data1)[2])-1,(t(data1[1,-1])),col='red')

data1=read.delim('CHB_MAFpop0.txt')
lines(c(2:dim(data1)[2])-1,(t(data1[1,-1])),col='blue')

#############
# Estimate parameters of demographic model
#############

# if not, you can estimate the parameters of your demographic model using fastsimcoal2. You'll need to create a .est and .tpl file as shown in the manual and then run several (Laurent Excoffier recommends ~100) iterations as in the command below
# see also https://speciationgenomics.github.io/fastsimcoal2/

#!/bin/bash

# for f in {1..100} ; do cat estimate.sh | sed -e s/ITER/$f/g > estimate.$f.sh; done
# touch commands1.sh; for f in {1..100} ; do echo "sh "estimate.$f.sh >> commands1.sh; done

cd /Genomics/ayroleslab2/alea/tsimane/int_data_files/demography2/fastsimcoal2/

cp CHB_PEL_TSI.tpl CHB_PEL_TSI_ITER.tpl
cp CHB_PEL_TSI.est CHB_PEL_TSI_ITER.est
cp CHB_PEL_TSI_jointMAFpop1_0.obs CHB_PEL_TSI_ITER_jointMAFpop1_0.obs
cp CHB_PEL_TSI_jointMAFpop2_0.obs CHB_PEL_TSI_ITER_jointMAFpop2_0.obs
cp CHB_PEL_TSI_jointMAFpop2_1.obs CHB_PEL_TSI_ITER_jointMAFpop2_1.obs

~/programs/fsc26_linux64/fsc26 -t CHB_PEL_TSI_ITER.tpl -e CHB_PEL_TSI_ITER.est -n 100000 -m -M -L 40 –c0 -C 10
# ~/programs/fsc26_linux64/fsc26 -t CHB_PEL_TSI.tpl -e CHB_PEL_TSI.est -n 100000 -m -M -L 40 –c0 -C 10

rm CHB_PEL_TSI_ITER.est
rm CHB_PEL_TSI_ITER.tpl
rm CHB_PEL_TSI_ITER_jointMAFpop1_0.obs
rm CHB_PEL_TSI_ITER_jointMAFpop2_0.obs
rm CHB_PEL_TSI_ITER_jointMAFpop2_1.obs

# this will output, among other things, the likelihood of each set of simulated parameters. You can compare over the 100 iterations to find the best likelihood and use those parameters to simulate your genetic data.

rm best_lik_param2.txt; touch best_lik_param2.txt
for f in {1..100} ; do cat CHB_PEL_TSI_${f}/CHB_PEL_TSI_${f}.bestlhoods | tail -1 >> best_lik_param2.txt; done

#############
# simulate DNA sequence under the inferred model
#############

#Now that you have maxlikelihood parameters for your demographic model, you can simulate DNA sequence under this model to use for PBS and iHS

cd /Genomics/ayroleslab2/alea/tsimane/int_data_files/demography2/fastsimcoal2/

~/programs/fsc26_linux64/fsc26 -i CHB_PEL_TSI_maxL.par -n 100 -j -q -m -g -s0

#where dnachb_pel_xal_maxL.par looks just like your .tpl file from step 2 but with the maxlikelihood parameters substituted for the variables
#have a look at the fastsimcoal2 documentation for the best way to get the amount of SNPs/sequence data for your purposes.

#############
# convert to vcf
#############

cd /Genomics/ayroleslab2/alea/tsimane/int_data_files/demography2/fastsimcoal2/CHB_PEL_TSI_maxL

##

install.packages("remotes")
remotes::install_github("dinmatias/reconproGS")
library(reconproGS)

source('arp2vcf_mod.R')

filePath1 <- "CHB_PEL_TSI_maxL_1_1.arp"
x <- arp2vcf(filePath = filePath1, outFile="CHB_PEL_TSI_maxL_1_1.vcf" )
# names(x)[9:dim(x)[2]]<-paste('ID',9:dim(x)[2],sep='')

for (i in c(2:100)){
filePath1 <- paste('CHB_PEL_TSI_maxL_1_',i,'.arp',sep='')
x2 <- arp2vcf(filePath = filePath1, outFile=paste('/Genomics/ayroleslab2/alea/tsimane/int_data_files/demography2/fastsimcoal2/CHB_PEL_TSI_maxL/CHB_PEL_TSI_maxL_1_',i,'.vcf',sep='') )
x<-rbind(x,x2)
}

x3<-subset(x,REF %in% c(0,1,2,3) & ALT %in% c(0,1,2,3))
names(x3)[10:dim(x3)[2]]<-paste("ID",(10:dim(x3)[2])-9,sep='')
write.table(x3,'CHB_PEL_TSI_maxL_23Feb21.vcf',row.names=F,sep='\t',quote=F)

#############
# calculate Fst for downstream PBS
#############

path_plink=/Genomics/grid/users/alea/programs/plink_1.90
$path_plink --vcf CHB_PEL_TSI_maxL_23Feb21.vcf --out CHB_PEL_TSI_maxL_23Feb21 --make-bed --allow-extra-chr

grep -v 'Pop1' pop_info.txt > pop_info_23.txt
grep -v 'Pop3' pop_info.txt > pop_info_12.txt
grep -v 'Pop2' pop_info.txt > pop_info_13.txt

$path_plink --bfile CHB_PEL_TSI_maxL_23Feb21 --out CHB_PEL_TSI_maxL_23Feb21_2vs3 --within pop_info_23.txt --fst --allow-extra-chr
$path_plink --bfile CHB_PEL_TSI_maxL_23Feb21 --out CHB_PEL_TSI_maxL_23Feb21_1vs3 --within pop_info_13.txt --fst --allow-extra-chr
$path_plink --bfile CHB_PEL_TSI_maxL_23Feb21 --out CHB_PEL_TSI_maxL_23Feb21_1vs2 --within pop_info_12.txt --fst --allow-extra-chr
$path_plink --bfile CHB_PEL_TSI_maxL_23Feb21 --out CHB_PEL_TSI_maxL_23Feb21 --within pop_info.txt --freq --allow-extra-chr


