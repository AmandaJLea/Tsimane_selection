################
# get new GRM
################

awk '{OFS="\t"; print $1,$2,$3,$4,$5,1}' TOPMedImputed_clean_ALLchr_array.fam > temp.fam
mv temp.fam TOPMedImputed_clean_ALLchr_array.fam

~/programs/gemma-0.98.1-linux-static -bfile TOPMedImputed_clean_ALLchr_array -gk 2 -o TOPMedImputed_clean_ALLchr_array

# king
~/programs/king -b TOPMedImputed_clean_ALLchr_array.bed --kinship --prefix king

################
# set up for GEMMA
################

pca=read.delim('TOPMedImputed_clean_ALLchr_array.eigenvec',header=F)
pheno=read.csv('/Genomics/ayroleslab2/alea/tsimane/data_from_Ben_Heather/genetics_sample_phenotype_clean_9Dec20.csv')

pheno=merge(pca[,c(2:7)],pheno,by.x='V2',by.y='pid')
names(pheno)[1]<-'pid'

# 10,11,14,2,20,5,6
geno1=read.delim('/Genomics/ayroleslab2/alea/tsimane/int_data_files/TOPMedImputed_clean_chr10_candidates.traw')
geno2=read.delim('/Genomics/ayroleslab2/alea/tsimane/int_data_files/TOPMedImputed_clean_chr11_candidates.traw')
geno3=read.delim('/Genomics/ayroleslab2/alea/tsimane/int_data_files/TOPMedImputed_clean_chr14_candidates.traw')
geno4=read.delim('/Genomics/ayroleslab2/alea/tsimane/int_data_files/TOPMedImputed_clean_chr2_candidates.traw')
geno5=read.delim('/Genomics/ayroleslab2/alea/tsimane/int_data_files/TOPMedImputed_clean_chr20_candidates.traw')
geno6=read.delim('/Genomics/ayroleslab2/alea/tsimane/int_data_files/TOPMedImputed_clean_chr5_candidates.traw')
geno7=read.delim('/Genomics/ayroleslab2/alea/tsimane/int_data_files/TOPMedImputed_clean_chr6_candidates.traw')
genoCAND=rbind(geno1,geno2,geno3,geno4,geno5,geno6,geno7)

genoALL=read.delim('TOPMedImputed_clean_ALLchr_array.traw',nrow=1)
genoCAND2=genoCAND[,names(genoALL)]
genoALL2=rbind(genoCAND2,genoALL)

geno1_ids=read.delim('/Genomics/ayroleslab2/alea/tsimane/int_data_files/TOPMedImputed_clean_ALLchr_array.fam',header=F)

grm=read.delim('/Genomics/ayroleslab2/alea/tsimane/int_data_files/output/TOPMedImputed_clean_ALLchr_array.sXX.txt',header=F)
grm_ids<-geno1_ids
colnames(grm)<-grm_ids$V2
rownames(grm)<-grm_ids$V2

# phenotype file
pheno2<-pheno[complete.cases(pheno[,c(12,13)]),]
pheno3<-subset(pheno2, pid %in% names(grm))
write.table(pheno3[,c(17:45)],'GEMMA_phenotype.txt',row.names=F,sep='\t',col.names=F,quote=F)

# covariate file
write.table(cbind(1,pheno3[,c(2:6,12,13)]),'GEMMA_covariate.txt',row.names=F,sep='\t',col.names=F,quote=F)

# grm
grm2<-grm[as.character(pheno3$pid),as.character(pheno3$pid)]
tmp<-cbind(as.character(colnames(grm2)),as.character(rownames(grm2)),as.character(pheno3$pid))
write.table(format(grm2,digits=4),'GEMMA_GRM.txt',row.names=F,sep='\t',col.names=F,quote=F)

# genotype file
names(genoALL2)[-c(1:6)]<-as.character(geno1_ids$V2)
genoALL3<-genoALL2[,as.character(pheno3$pid)]
genoALL3<-cbind(genoALL2[,c(2,5,6)],genoALL3)

# need loc, ref, alt
write.table(genoALL3,'GEMMA_genotypes_subset.txt',row.names=F,sep='\t',col.names=F,quote=F)

################
# run GEMMA
################

for k in {1..29}; do sed -e s/TRAITNUMBER/$k/g GEMMA.sh > gemma$k.sh; done
for k in {1..29}; do echo "sh "gemma$k.sh >> commands1.sh; done

awk '{OFS="\t";print $1,$7,$8}' GEMMA_covariate.txt > GEMMA_covariate_v2.txt

#!/bin/bash

k=TRAITNUMBER

~/programs/gemma-0.98.1-linux-static -k GEMMA_GRM.txt -g GEMMA_genotypes_subset.txt -c GEMMA_covariate_v2.txt -p GEMMA_phenotype.txt -lmm 1 -n $k -o TOPMedImputed_clean_ALLchr_array.trait$k.v2 

################
# output
################

i=1
output=read.delim(paste('output/TOPMedImputed_clean_ALLchr_array.trait',i,'.v2.assoc.txt',sep=''))
output$trait<-names(pheno2[,c(17:45)])[i]
output<-unique(output)
output$qvalue<-p.adjust(output$p_wald,method='BH')

for (i in 2:29){

output2=read.delim(paste('output/TOPMedImputed_clean_ALLchr_array.trait',i,'.v2.assoc.txt',sep=''))
output2$trait<-names(pheno2[,c(17:45)])[i]
output2<-unique(output2)
output2$qvalue<-p.adjust(output2$p_wald,method='BH')

output=rbind(output,output2)
}

write.table(output,'31Jan21_GEMMA_results.txt',row.names=F,sep='\t')

