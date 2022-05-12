# boot_counts.sh
#!/bin/bash
##
## Count the number of mutations in each spectral class among 100 bins from a bootstrapped list of dated variants.
## Same procedure as count_muts.sh, but iterates over resample of variants in each bin.
##
## Syntax:
##   ./count_muts.sh <variant file> <bin size>
## Input variant file must be ordered by bootstrap iteration and sorted by allele age, bin size in number of rows (single iteration)
##
## Output column order:
## step, A_C, ACC_ATC, A_G, A_T, C_A, CCC_CTC, C_G, CpG, C_T, G_A, G_C, G_T, T_A, T_C, TCC_TTC, TCT_TTT, T_G, bin_age
##

infile=$1
step=$2
upper=$(($step*100))
for((i=1;i<$upper;i+=$step))
do
for((j=0;j<100;j++))
do
echo -en $i'\t'
tail -n +$i $infile| head -n $step | shuf -n $step -r |
awk 'BEGIN { arr["A_C"] = 0; arr["A_G"] = 0; arr["A_T"] = 0; arr["C_A"] = 0; arr["C_G"] = 0; arr["C_T"] = 0; arr["xavg"] = 0;
             arr["T_G"] = 0; arr["T_C"] = 0; arr["T_A"] = 0; arr["G_T"] = 0; arr["G_C"] = 0; arr["G_A"] = 0; arr["CpG"] = 0;
             arr["TCC_TTC"] = 0; arr["ACC_ATC"] = 0; arr["TCT_TTT"] = 0; arr["CCC_CTC"] = 0 }
{x += $8; mut = $5"_"$6;
if($9 ~ /CG./ || $9 ~ /.CG/) mut = "CpG";
if(($9 == "TCC" && $10 == "TTC") || ($9 == "GGA" && $10 == "GAA")) mut = "TCC_TTC";
if(($9 == "ACC" && $10 == "ATC") || ($9 == "GGT" && $10 == "GAT")) mut = "ACC_ATC";
if(($9 == "TCT" && $10 == "TTT") || ($9 == "AGA" && $10 == "AAA")) mut = "TCT_TTT";
if(($9 == "CCC" && $10 == "CTC") || ($9 == "GGG" && $10 == "GAG")) mut = "CCC_CTC";
if(mut in arr) arr[mut]++ }
END {arr["xavg"] = x / NR; for (key in arr) { print key,arr[key] }}' | sort -k1,1 | cut -d' ' -f2 | tr '\n' '\t'
echo
done
done