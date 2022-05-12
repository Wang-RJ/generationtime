# Row cutoffs
# tot, AFR, EAS, EUR, SAS
#	100kga		10kga		1kga
# 33091767	27434387	17409867
# 19900936	14291019	5410609
# 9395289	5596644		3623724
# 10001088	5817294		3217441
# 11252993	7118739		4326103

# Split dataset into shared and private subsets

awk '{if(($14 > 0.0008) && 
  ($16 > 0.001) &&
  ($17 > 0.001) &&
  ($18 > 0.001))
    print $0}' df_v4_dfreq_MAF_biallelic_DAF98.agesorted.txt > SHA_tgppops_der.sorted.txt

awk '{if(($14 > 0.0008) && 
  ($16 == 0) &&
  ($17 == 0) &&
  ($18 == 0))
    print $0}' df_v4_dfreq_MAF_biallelic_DAF98.agesorted.txt > AFR_private.sorted.txt

awk '{if(($14 == 0) && 
  ($16 > 0.001) &&
  ($17 == 0) &&
  ($18 == 0))
    print $0}' df_v4_dfreq_MAF_biallelic_DAF98.agesorted.txt > EAS_private.sorted.txt
	
awk '{if(($14 == 0) && 
  ($16 == 0) &&
  ($17 > 0.001) &&
  ($18 == 0))
    print $0}' df_v4_dfreq_MAF_biallelic_DAF98.agesorted.txt > EUR_private.sorted.txt
	
awk '{if(($14 == 0) && 
  ($16 == 0) &&
  ($17 == 0) &&
  ($18 > 0.001))
    print $0}' df_v4_dfreq_MAF_biallelic_DAF98.agesorted.txt > SAS_private.sorted.txt

# Count shared mutations between populations in total dataset and at 10kga
./count_shared.sh df_v4_dfreq_MAF_biallelic_DAF98.agesorted.txt 330978 > TGP_allsharing.txt
head -n 27434387 df_v4_dfreq_MAF_biallelic_DAF98.agesorted.txt > tmp_tgppops 
./count_shared.sh tmp_tgppops 274343 > TGP_10kgasharing.txt

wc -l SHA_tgppops_der.sorted.txt
# 5043360 SHA_tgppops_der.sorted.txt
./count_muts.sh SHA_tgppops_der.sorted.txt 50433 > SHA_allcounts.txt

wc -l AFR_private.sorted.txt
# 11359486 AFR_private.sorted.txt
wc -l EAS_private.sorted.txt 
# 2923625 EAS_private.sorted.txt
wc -l EUR_private.sorted.txt 
# 1488404 EUR_private.sorted.txt
wc -l SAS_private.sorted.txt 
# 2991482 SAS_private.sorted.txt

./count_muts.sh AFR_private.sorted.txt 113594 > AFR_privatecounts.txt
./count_muts.sh EAS_private.sorted.txt 29236 > EAS_privatecounts.txt
./count_muts.sh EUR_private.sorted.txt 14884 > EUR_privatecounts.txt
./count_muts.sh SAS_private.sorted.txt 29914 > SAS_privatecounts.txt

cut -f8 AFR_private.sorted.txt | grep -En ^1100 | head
# 4749696
cut -f8 EAS_private.sorted.txt | grep -En ^1100 | head
# 2849337
cut -f8 EUR_private.sorted.txt | grep -En ^1100 | head
# 1468389
cut -f8 SAS_private.sorted.txt | grep -En ^1100 | head
# 2848666

head -n 4749696 AFR_private.sorted.txt > AFR_tgppops_priv1kga.txt
./count_muts.sh AFR_tgppops_priv1kga.txt 47496 > AFR_1kgaprivatecounts.txt
./boot_counts.sh AFR_tgppops_priv1kga.txt 47496 > bootstraps/AFR_boots_priv1kga.txt

head -n 2849337 EAS_private.sorted.txt > EAS_tgppops_priv1kga.txt
./count_muts.sh EAS_tgppops_priv1kga.txt 28493 > EAS_1kgaprivatecounts.txt
./boot_counts.sh EAS_tgppops_priv1kga.txt 284931 > bootstraps/EAS_boots_priv1kga.txt

head -n 1468389 EUR_private.sorted.txt > EUR_tgppops_priv1kga.txt
./count_muts.sh EUR_tgppops_priv1kga.txt 14683 > EUR_1kgaprivatecounts.txt
./boot_counts.sh EUR_tgppops_priv1kga.txt 14683 > bootstraps/EUR_boots_priv1kga.txt

head -n 2848666 SAS_private.sorted.txt > SAS_tgppops_priv1kga.txt
./count_muts.sh SAS_tgppops_priv1kga.txt 28486 > SAS_1kgaprivatecounts.txt
./boot_counts.sh SAS_tgppops_priv1kga.txt 28486 > bootstraps/SAS_boots_priv1kga.txt
