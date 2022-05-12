# Sort dataset by recombination rate and split by quintiles

sort -g -k11,11 TGP_tgppops_10kga.txt > TGP_tgppops_10kga.recsorted.txt

head -n 5486877 TGP_tgppops_10kga.recsorted.txt > TGP_tgppops_10kga.recbin1.txt
tail -n +5486878 TGP_tgppops_10kga.recsorted.txt | head -n 5486877 > TGP_tgppops_10kga.recbin2.txt
tail -n +10973755 TGP_tgppops_10kga.recsorted.txt | head -n 5486877 > TGP_tgppops_10kga.recbin3.txt
tail -n +16460633 TGP_tgppops_10kga.recsorted.txt | head -n 5486877 > TGP_tgppops_10kga.recbin4.txt
tail -n +21947510 TGP_tgppops_10kga.recsorted.txt > TGP_tgppops_10kga.recbin5.txt

# Re-sort each quintile by allele age to prepare for call to count_muts
sort -g -k8,8 TGP_tgppops_10kga.recbin1.txt > TGP_tgppops_10kga.recbin1.agesorted.txt
sort -g -k8,8 TGP_tgppops_10kga.recbin2.txt > TGP_tgppops_10kga.recbin2.agesorted.txt
sort -g -k8,8 TGP_tgppops_10kga.recbin3.txt > TGP_tgppops_10kga.recbin3.agesorted.txt
sort -g -k8,8 TGP_tgppops_10kga.recbin4.txt > TGP_tgppops_10kga.recbin4.agesorted.txt
sort -g -k8,8 TGP_tgppops_10kga.recbin5.txt > TGP_tgppops_10kga.recbin5.agesorted.txt

./count_muts.sh TGP_tgppops_10kga.recbin1.agesorted.txt 54868 > TGP_10kga.recbin1.txt
./count_muts.sh TGP_tgppops_10kga.recbin2.agesorted.txt 54868 > TGP_10kga.recbin2.txt
./count_muts.sh TGP_tgppops_10kga.recbin3.agesorted.txt 54868 > TGP_10kga.recbin3.txt
./count_muts.sh TGP_tgppops_10kga.recbin4.agesorted.txt 54868 > TGP_10kga.recbin4.txt
./count_muts.sh TGP_tgppops_10kga.recbin5.agesorted.txt 54868 > TGP_10kga.recbin5.txt
