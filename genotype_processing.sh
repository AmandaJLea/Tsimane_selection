path_plink=/Genomics/grid/users/alea/programs/plink_1.90
path_king=/Genomics/grid/users/alea/programs/king

path_resources=/Genomics/ayroleslab2/alea/ref_genomes/public_datasets
path_out=/Genomics/ayroleslab2/alea/tsimane/int_data_files
path_AMR=/Genomics/ayroleslab2/alea/tsimane/1000G_compare_ids.txt
# with han chinese
path_AMR2=/Genomics/ayroleslab2/alea/tsimane/1000G_compare_ids2.txt
# PEL, GBR, YRI, han chinese (unrelated)
path_AMR3=/Genomics/ayroleslab2/alea/tsimane/1000G_compare_ids3.txt

path_array=/Genomics/ayroleslab2/alea/ref_genomes/public_datasets/ALL.chip.omni_broad_sanger_combined.20140818.snps.genotypes.vcf.gz
path_harmonizer=~/programs/GenotypeHarmonizer-1.4.20-SNAPSHOT/GenotypeHarmonizer.jar

# Tsimane data, raw
grep 'NA' /Genomics/ayroleslab2/alea/tsimane/data_from_Ben_Heather/raw_genotyped/TSI1_MEG.fam > to_remove.txt
path_tsimane=/Genomics/ayroleslab2/alea/tsimane/data_from_Ben_Heather/raw_genotyped/TSI1_MEG

# filter Tsimane - missingness, remove 1KG
$path_plink --bfile $path_tsimane --geno 0.1 --mind 0.1 --make-bed --out $path_out/TSI1_MEG_filt --chr 1-22 --remove to_remove.txt --het --allow-no-sex

# filter Tsimane - remove het, add HWE
$path_plink --bfile $path_out/TSI1_MEG_filt --geno 0.1 --mind 0.1 --make-bed --out $path_out/TSI1_MEG_filt2 --chr 1-22 --remove to_remove2.txt --het --allow-no-sex --hwe 0.00000001 --biallelic-only

# make new fam
awk '$2=(FNR FS $1)' $path_out/TSI1_MEG_filt2.fam |  awk '{print $2,$2,0,0,0,-9}' > temp.txt 
mv temp.txt $path_out/TSI1_MEG_filt2.fam

###########
# filter Tsimane
########### 

# use unrelated set from PC-Relate
# filter Tsimane - relatedness
$path_plink --bfile $path_out/TSI1_MEG_filt2 --geno 0.1 --mind 0.1 --make-bed --out $path_out/TSI1_MEG_filt3 --chr 1-22 --keep TSI1_MEG_filt2_PCRelate_unrelated.txt --hwe 0.00000001 --het --allow-no-sex

# change rsID
awk '{OFS="\t";print $1,$1":"$4,$3,$4,$5,$6}' /Genomics/ayroleslab2/alea/tsimane/int_data_files/TSI1_MEG_filt3.bim > temp.bim
mv temp.bim /Genomics/ayroleslab2/alea/tsimane/int_data_files/TSI1_MEG_filt3.bim

###########
# harmonize to 1KG WGS
###########

#!/bin/bash

# for f in {1..22} ; do cat harm.sh | sed -e s/CHROMNAME/$f/g > harm.$f.sh; done
# for f in {1..22} ; do echo "sh "harm.$f.sh >> commands1.sh; done
path_harmonizer=~/programs/GenotypeHarmonizer-1.4.20-SNAPSHOT/GenotypeHarmonizer.jar
path_resources=/Genomics/ayroleslab2/alea/ref_genomes/public_datasets

f=CHROMNAME
java -Xmx50g -jar $path_harmonizer --input /Genomics/ayroleslab2/alea/tsimane/int_data_files/TSI1_MEG_filt3 --ref $path_resources/shapeit_files/ALL.chr${f}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --output /Genomics/ayroleslab2/alea/tsimane/int_data_files/TSI1_MEG_filt3_harm_chr${f} --outputType PLINK_BED --update-reference-allele --ambiguousSnpFilter --update-id --chrFilter ${f}

###########
# filter 1KG WGS and merge with Tsimane
########### 

#!/bin/bash

# for f in {1..22} ; do cat phase.sh | sed -e s/CHROMNAME/$f/g > phase.$f.sh; done
# for f in {1..22} ; do echo "sh "phase.$f.sh >> commands1.sh; done

f=CHROMNAME
path_resources=/Genomics/ayroleslab2/alea/ref_genomes/public_datasets
path_plink=/Genomics/grid/users/alea/programs/plink_1.90
path_out=/Genomics/ayroleslab2/alea/tsimane/int_data_files
path_AMR3=/Genomics/ayroleslab2/alea/tsimane/1000G_compare_ids3.txt

# 1000G data filter
$path_plink --vcf $path_resources/shapeit_files/ALL.chr${f}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --geno 0.1 --make-bed --out $path_resources/shapeit_files/ALL.chr${f}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_AMR --keep $path_AMR3 --allow-extra-chr --chr $f --snps-only

# 1000G data merge
path_tsimane=/Genomics/ayroleslab2/alea/tsimane/int_data_files/TSI1_MEG_filt3_harm_chr${f}
path_1kg=$path_resources/shapeit_files/ALL.chr${f}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_AMR

$path_plink --bfile $path_tsimane --bmerge $path_1kg --make-bed --out $path_out/TSI1_MEG_filt3_w1000G_WGS_${f} --geno 0.025 --keep-allele-order --allow-no-sex

remove=/Genomics/ayroleslab2/alea/tsimane/int_data_files/TSI1_MEG_filt3_w1000G_WGS_${f}-merge.missnp
$path_plink --bfile $path_tsimane --make-bed --out temp1_${f} --exclude $remove
$path_plink --bfile $path_1kg --make-bed --out temp2_${f} --exclude $remove

$path_plink --bfile temp1_${f} --bmerge temp2_${f} --make-bed --out $path_out/TSI1_MEG_filt3_w1000G_WGS_${f} --geno 0.025 --allow-no-sex --keep-allele-order

###########
# phase
########### 

#!/bin/bash

# for f in {1..22} ; do cat phase.sh | sed -e s/NUMBER/$f/g > phase.$f.sh; done
chr=NUMBER

path_plink=/Genomics/grid/users/alea/programs/plink_1.90
path_resources=/Genomics/ayroleslab2/alea/ref_genomes/public_datasets
path_out=/Genomics/ayroleslab2/alea/tsimane/int_data_files

# phase w/ref panel
# version of files without ref panel doesn't have 2 after 'phased' in output
map=/Genomics/ayroleslab2/alea/ref_genomes/public_datasets/1000GP_Phase3/genetic_map_chr${chr}_combined_b37.txt
phased_haps=/Genomics/ayroleslab2/alea/tsimane/selection_files/TSI1_MEG_filt3_w1000G_WGS_chr${chr}.phased2.haps
phased_sample=/Genomics/ayroleslab2/alea/tsimane/selection_files/TSI1_MEG_filt3_w1000G_WGS_chr${chr}.phased2.sample
path_ref=/Genomics/ayroleslab2/alea/ref_genomes/public_datasets/shapeit_files/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes

/Genomics/grid/users/alea/programs/shapeit.v2.904.2.6.32-696.18.7.el6.x86_64/bin/shapeit --input-bed $path_out/TSI1_MEG_filt3_w1000G_WGS_${chr} \
        --input-map $map --output-log $path_ref.log\
        --output-max $phased_haps $phased_sample --thread 8 --force --input-ref $path_ref.hap.gz $path_ref.legend.gz $path_ref.samples

