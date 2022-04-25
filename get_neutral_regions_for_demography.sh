# filtered for missing genotype data, HWE, and relatedness; harmonized to 1KG
# not filtered for MAC or MAF
# still includes Moseten and admixed Tsimane, admixture not masked
path_out=/Genomics/ayroleslab2/alea/tsimane/int_data_files

#############
# merge files and remove admixed IDs and Moseten
#############

data1=/Genomics/ayroleslab2/alea/tsimane/int_data_files/TSI1_MEG_filt3_w1000G_WGS_1
~/programs/plink_1.90 --bfile $data1 --merge-list merge_list.txt --keep demography.txt --out $path_out/TSI1_MEG_filt3_w1000G_4demog --make-bed --allow-no-sex

#############
# get bed files of non neutral sites to remove
#############

# The (putatively) neutral regions are defined as: Outside of genic regions, conserved regions, CpG islands, 1kg mask

# region filtering
cpgs=/Genomics/ayroleslab2/alea/ref_genomes/public_datasets/CpGs_hg38_v2.bed
mask=/Genomics/ayroleslab2/alea/ref_genomes/public_datasets/20160622.allChr.pilot_mask_v2.bed
super=/Genomics/ayroleslab2/alea/ref_genomes/public_datasets/genomicSuperDups_hg38.bed
cat $cpgs $super > /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/temp.bed

awk '$4=(FNR FS)' /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/temp.bed > /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/temp_v2.bed

awk '$1 !~ /_/' /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/temp_v2.bed > /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/temp_v3.bed

# lift to hg19
~/programs/liftOver /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/temp_v3.bed ~/programs/hg38ToHg19.over.chain.gz /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/temp_v3_hg19.bed unMapped

~/programs/liftOver $mask ~/programs/hg38ToHg19.over.chain.gz /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/20160622.allChr.pilot_mask_hg19.bed unMapped

sed -e s/chr//g /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/20160622.allChr.pilot_mask_hg19.bed | awk '{print $0,NR}' | grep -v 'random' | grep -v 'Un' > /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/20160622.allChr.pilot_mask_hg19_nochr.bed

sed -e s/chr//g /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/temp_v3_hg19.bed | awk '{print $0,NR}' | grep -v 'Un' | grep -v 'random'> /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/temp_v3_hg19_nochr.bed

#############
# remove non neutral sites
#############

mask_hg19=/Genomics/ayroleslab2/alea/ref_genomes/public_datasets/20160622.allChr.pilot_mask_hg19_nochr.bed
exclude_hg19=/Genomics/ayroleslab2/alea/ref_genomes/public_datasets/temp_v3_hg19_nochr.bed

# filter regions - part 1
~/programs/plink_1.90 --bfile $path_out/TSI1_MEG_filt3_w1000G_4demog --make-bed --extract range $mask_hg19 --exclude range $exclude_hg19 --allow-extra-chr --out $path_out/TSI1_MEG_filt3_w1000G_4demog_neut1

# filter regions - part 2
conserved=/Genomics/ayroleslab2/alea/ref_genomes/public_datasets/GERP_conserved_elements/hg19_all_elements.bed
awk '{print $0,NR}' $conserved > /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/GERP_conserved_elements/hg19_all_elements_v2.bed

# reorder to CHB, PEL, TSI
~/programs/plink_1.90 --bfile $path_out/TSI1_MEG_filt3_w1000G_4demog_neut1 --make-bed --exclude range temp.txt --out $path_out/TSI1_MEG_filt3_w1000G_4demog_neut2 -geno 0 --allow-no-sex --indiv-sort f reorder.txt

# export to vcf
~/programs/plink2 --bfile $path_out/TSI1_MEG_filt3_w1000G_4demog_neut2 --recode vcf --out $path_out/TSI1_MEG_filt3_w1000G_4demog_neut3 --chr 1-22


