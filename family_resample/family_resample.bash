# redraw ages using lower and upper 95 CI, feeds these values to python draw95CI.py which draws from a normal
cut -f23,24 df_v4_dfreq_MAF_biallelic_DAF98.agesorted.txt | python draw95CI.py > age_redrawn.columns
paste df_v4_dfreq_MAF_biallelic_DAF98.agesorted.txt age_redrawn.columns > df_v4_redraw.txt

# redrawn mean age in columns 25-34 of df_v4_redraw.txt
# sort each redrawn data frame
for i in {25..34}
do
  paste <(cut -f1-7 df_v4_redraw.txt) <(cut -f$i df_v4_redraw.txt) <(cut -f9-34 df_v4_redraw.txt) > df_v4_redraw.$i
  sort -g -k$i,$i df_v4_redraw.$i > df_v4_redraw.$i.sorted
  rm df_v4_redraw.$i
done

# Find 10kga cutoff for each redraw and write to cut10kga_step
for i in {25..34}
do
  cut -f8 df_v4_redraw.$i.sorted | grep -En ^100[0-9]{2} | head -n1 >> cut10kga.idx
done

cut -F: '{print $1}' cut10kga.idx > cut10kga
awk '{print substr($1, 1, length($1)-2)}' cut10kga > cut10kga_step

# count_muts_zeros, set negative allele ages to zero
for i in {1..10}
do
head -n $(tail -n +$i cut10kga | head -n1) df_v4_redraw.$(($i + 24)).sorted > TGP_redraw$(($i + 24))_10kga.txt
done

for i in {1..10}
do
./count_muts_zeros.sh TGP_redraw$(($i + 24))_10kga.txt $(tail -n +$i cut10kga_step | head -n1) > TGP_10kgacounts.redraw$(($i + 24))
done
