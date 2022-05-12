# count_shared.sh
#!/bin/bash
#
## Count shared mutations between populations by bin
## Syntax:
##   ./count_muts.sh <variant file> <bin size>
## Input variant file must be sorted by allele age, bin size in number of rows
##
## bitkeys:
## 1 = AFR
## 2 = EAS
## 4 = EUR
## 8 = SAS

infile=$1
step=$2
upper=$(($step*100))
for((i=1;i<$upper;i+=$step))
do
tail -n +$i $infile| head -n $step | 
awk 'BEGIN { 
  arr["xavg"] = 0;
  for(i = 1; i <= 16; i++)
    arr[i] = 0
  }
  {
  bitkey = 0;
  if($14 > 0.0008) bitkey += 1;
  if($16 > 0.001) bitkey += 2;
  if($17 > 0.001) bitkey += 4;
  if($18 > 0.001) bitkey += 8;
  arr[bitkey]++
  x += $8
  }
END {
  arr["xavg"] = x / NR; for (key in arr) { print key,arr[key] }
  }' | sort -k1,1 | cut -d' ' -f2 | tr '\n' '\t'
echo
done
