# Count mutations for full TGP dataset and do bootstrap counts
./count_muts.sh df_v4_dfreq_MAF_biallelic_DAF98.agesorted.txt 330978 > TGP_allcounts.txt
./boot_counts.sh df_v4_dfreq_MAF_biallelic_DAF98.agesorted.txt 330978 > bootstraps/TGP_boots_allcounts.txt

# Subset into separate files for each  population
awk '{if($14 > 0.0008) print $0}' df_v4_dfreq_MAF_biallelic_DAF98.agesorted.txt > AFR_tgppops_der.sorted.txt
awk '{if($16 > 0.001) print $0}' df_v4_dfreq_MAF_biallelic_DAF98.agesorted.txt > EAS_tgppops_der.sorted.txt
awk '{if($17 > 0.001) print $0}' df_v4_dfreq_MAF_biallelic_DAF98.agesorted.txt > EUR_tgppops_der.sorted.txt
awk '{if($18 > 0.001) print $0}' df_v4_dfreq_MAF_biallelic_DAF98.agesorted.txt > SAS_tgppops_der.sorted.txt

# Count number of rows in each to break into 100 bins for all counts
wc -l AFR_tgppops_der.sorted.txt
# 19906666 AFR_tgppops_der.sorted.txt
./count_muts.sh AFR_tgppops_der.sorted.txt 199066 > AFR_allcounts.txt
./boot_counts.sh AFR_tgppops_der.sorted.txt 199066 > bootstraps/AFR_boots_allcounts.txt

wc -l EAS_tgppops_der.sorted.txt
# 9400408 EAS_tgppops_der.sorted.txt
./count_muts.sh EAS_tgppops_der.sorted.txt 94004 > EAS_allcounts.txt
./boot_counts.sh EAS_tgppops_der.sorted.txt 94004 > bootstraps/EAS_boots_allcounts.txt

wc -l EUR_tgppops_der.sorted.txt
# 10006675 EUR_tgppops_der.sorted.txt
./count_muts.sh EUR_tgppops_der.sorted.txt 100066 > EUR_allcounts.txt
./boot_counts.sh EUR_tgppops_der.sorted.txt 100066 > bootstraps/EUR_boots_allcounts.txt

wc -l SAS_tgppops_der.sorted.txt
# 11258460 SAS_tgppops_der.sorted.txt
./count_muts.sh SAS_tgppops_der.sorted.txt 112584 > SAS_allcounts.txt
./boot_counts.sh SAS_tgppops_der.sorted.txt 112584 > bootstraps/SAS_boots_allcounts.txt

# Find row numbers for 10k and 1k generation cutoff
cut -f8 df_v4_dfreq_MAF_biallelic_DAF98.agesorted.txt | grep -En ^100[0-9]{2}
cut -f8 AFR_tgppops_der.sorted.txt | grep -En ^100[0-9]{2}
cut -f8 EUR_tgppops_der.sorted.txt | grep -En ^100[0-9]{2}
cut -f8 EAS_tgppops_der.sorted.txt | grep -En ^100[0-9]{2}
cut -f8 SAS_tgppops_der.sorted.txt | grep -En ^100[0-9]{2}

# tot, AFR, EAS, EUR, SAS
#	100kga		10kga		1kga
# 33091767	27434387	17409867
# 19900936	14291019	5410609
# 9395289	5596644		3623724
# 10001088	5817294		3217441
# 11252993	7118739		4326103

head -n 27434387 df_v4_dfreq_MAF_biallelic_DAF98.agesorted.txt > TGP_tgppops_10kga.txt
./count_muts.sh TGP_tgppops_10kga.txt 274343 > TGP_10kgacounts.txt
./boot_counts.sh TGP_tgppops_10kga.txt 274343 > bootstraps/TGP_boots_10kga.txt

head -n 14291019 AFR_tgppops_der.sorted.txt > AFR_tgppops_10kga.txt
./count_muts.sh AFR_tgppops_10kga.txt 142910 > AFR_10kgacounts.txt
./boot_counts.sh AFR_tgppops_10kga.txt 142910 > bootstraps/AFR_boots_10kga.txt

head -n 5596644 EAS_tgppops_der.sorted.txt > EAS_tgppops_10kga.txt
./count_muts.sh EAS_tgppops_10kga.txt 55966 > EAS_10kgacounts.txt
./boot_counts.sh EAS_tgppops_10kga.txt 55966 > bootstraps/EAS_boots_10kga.txt

head -n 5817294 EUR_tgppops_der.sorted.txt > EUR_tgppops_10kga.txt
./count_muts.sh EUR_tgppops_10kga.txt 58172 > EUR_10kgacounts.txt
./boot_counts.sh EUR_tgppops_10kga.txt 58172 > bootstraps/EUR_boots_10kga.txt

head -n 7118739 SAS_tgppops_der.sorted.txt > SAS_tgppops_10kga.txt
./count_muts.sh SAS_tgppops_10kga.txt 71187 > SAS_10kgacounts.txt
./boot_counts.sh SAS_tgppops_10kga.txt 71187 > bootstraps/SAS_boots_10kga.txt

##

head -n 17409867 df_v4_dfreq_MAF_biallelic_DAF98.agesorted.txt > TGP_tgppops_1kga.txt
./count_muts.sh TGP_tgppops_1kga.txt 174098 > TGP_1kgacounts.txt
./boot_counts.sh TGP_tgppops_1kga.txt 174098 > bootstraps/TGP_boots_1kga.txt

head -n 5410609 AFR_tgppops_der.sorted.txt > AFR_tgppops_1kga.txt
./count_muts.sh AFR_tgppops_1kga.txt 54106 > AFR_1kgacounts.txt
./boot_counts.sh AFR_tgppops_1kga.txt 54106 > bootstraps/AFR_boots_1kga.txt

head -n 3623724 EAS_tgppops_der.sorted.txt > EAS_tgppops_1kga.txt
./count_muts.sh EAS_tgppops_1kga.txt 36237 > EAS_1kgacounts.txt
./boot_counts.sh EAS_tgppops_1kga.txt 36237 > bootstraps/EAS_boots_1kga.txt

head -n 3217441 EUR_tgppops_der.sorted.txt > EUR_tgppops_1kga.txt
./count_muts.sh EUR_tgppops_1kga.txt 32174 > EUR_1kgacounts.txt
./boot_counts.sh EUR_tgppops_1kga.txt 32174 > bootstraps/EUR_boots_1kga.txt

head -n 4326103 SAS_tgppops_der.sorted.txt > SAS_tgppops_1kga.txt
./count_muts.sh SAS_tgppops_1kga.txt 43261 > SAS_1kgacounts.txt
./boot_counts.sh SAS_tgppops_1kga.txt 43261 > bootstraps/SAS_boots_1kga.txt
