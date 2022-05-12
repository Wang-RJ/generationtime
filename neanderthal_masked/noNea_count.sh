# preprocess by removing variants with neanderthal masked flag, column 22
# output goes into neanderthal_masked/ folder to avoid name collision

awk '{if ($22 == "False") print $0}' df_v4_dfreq_MAF_biallelic_DAF98.agesorted.txt > df_v4_noNea.agesorted.txt

wc -l df_v4_noNea.agesorted.txt
# 20181826
./count_muts.sh df_v4_noNea.agesorted.txt 201818 > TGP_noNea_allcounts.txt

# Subset into different populations
awk '{if($14 > 0.0008) print $0}' df_v4_noNea.agesorted.txt > AFR_noNea_der.sorted.txt
awk '{if($16 > 0.001) print $0}' df_v4_noNea.agesorted.txt > EAS_noNea_der.sorted.txt
awk '{if($17 > 0.001) print $0}' df_v4_noNea.agesorted.txt > EUR_noNea_der.sorted.txt
awk '{if($18 > 0.001) print $0}' df_v4_noNea.agesorted.txt > SAS_noNea_der.sorted.txt

# Count mutations, 100 bins for entire file
wc -l AFR_noNea_der.sorted.txt
# 12100671
./count_muts.sh AFR_noNea_der.sorted.txt 121006 > AFR_noNea.allcounts.txt

wc -l EAS_noNea_der.sorted.txt
# 5578058
./count_muts.sh EAS_noNea_der.sorted.txt 55780 > EAS_noNea.allcounts.txt

wc -l EUR_noNea_der.sorted.txt
# 5937563
./count_muts.sh EUR_noNea_der.sorted.txt 59375 > EUR_noNea.allcounts.txt

wc -l SAS_noNea_der.sorted.txt
# 6705361
./count_muts.sh SAS_noNea_der.sorted.txt 67053 > SAS_noNea.allcounts.txt

wc -l EES_noNea_der.sorted.txt
# 10670856
./count_muts.sh EES_noNea_der.sorted.txt 106708 > EES_noNea.allcounts.txt

# Find allele age cut off, line number for 10k generations ago
cut -f8 df_v4_noNea.agesorted.txt | grep -En ^100[0-9]{2}
cut -f8 AFR_noNea_der.sorted.txt | grep -En ^100[0-9]{2}
cut -f8 EAS_noNea_der.sorted.txt | grep -En ^100[0-9]{2}
cut -f8 EUR_noNea_der.sorted.txt | grep -En ^100[0-9]{2}
cut -f8 SAS_noNea_der.sorted.txt | grep -En ^100[0-9]{2}

# tot, AFR, EAS, EUR, SAS
##	  10kga
# 16780964
# 8719721
# 3337918
# 3461152
# 4261044

head -n 16780964 df_v4_noNea.agesorted.txt > TGP_noNea_10kga.txt
./count_muts.sh TGP_noNea_10kga.txt 167809 > TGP_10kgacounts.txt

head -n 8719721 AFR_noNea_der.sorted.txt > AFR_noNea_10kga.txt
./count_muts.sh AFR_noNea_10kga.txt 87197 > AFR_10kgacounts.txt

head -n 3337918 EAS_noNea_der.sorted.txt > EAS_noNea_10kga.txt
./count_muts.sh EAS_noNea_10kga.txt 33379 > EAS_10kgacounts.txt

head -n 3461152 EUR_noNea_der.sorted.txt > EUR_noNea_10kga.txt
./count_muts.sh EUR_noNea_10kga.txt 34611 > EUR_10kgacounts.txt

head -n 4261044 SAS_noNea_der.sorted.txt > SAS_noNea_10kga.txt
./count_muts.sh SAS_noNea_10kga.txt 42610 > SAS_10kgacounts.txt

head -n 8079251 EES_noNea_der.sorted.txt > EES_noNea_10kga.txt
./count_muts.sh EES_noNea_10kga.txt 80792 > EES_10kgacounts.txt