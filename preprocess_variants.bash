# Preprocess df_v4.txt, data from 1000 Genomes Project with allele ages dated by GEVA (Albers and McVean)
# Columns in df_v4.txt ordered as (from GEVA first row, from 1000 Genomes second row):
# 	Chr, Pos, AlleleRef, AlleleAlt, AlleleAnc, DataSource, AgeMedian_Jnt,
# 	Ancestral, Derived, RecRate, ID, rsid, Ref, Alt, AFR_freq, AMR_freq, EAS_freq, EUR_freq, SAS_freq,
#   sameAlt(flag), Neanderthal_Masked(flag), AgeCI95Lower_Jnt(GEVA), AgeCI95Upper_Jnt(GEVA)

# Check whether the reference allele is the alternate allele and return the appropriate derived frequency
# Is the alt allele ($14) the same one referenced in the dating ($4)
# if derived ($6) == reference ($3), then return (1 - sum of alt freqs)
# else, return indexed entry

awk '{split($14, alt, ",");
if($4 == alt[1]) idx = 1;
else if($4 == alt[2]) idx = 2;
else if($4 == alt[3]) idx = 3;
else idx = 0;
for(i=15;i<=19;i++) {
	split($i,alt_freqs,",");
	if($6 == $3) {
		sum = 0; for(key in alt_freqs) { sum += alt_freqs[key] } output = 1 - sum }
	else
		output = alt_freqs[idx];
	printf "%f\t",output }
print ""}' df_v4.txt > df_v4_dfreq.only.txt

# Paste appropriate derived frequency columns into file 
paste <(cut -f1-11,13,14,20-22 df_v4.txt) df_v4_dfreq.only.txt > df_v4_derivedfreqs.txt

# Grab minor allele frequencies from 1000 Genomes dataset, match by rsid
# Start by creating two column text file with rsid and global MAF

zgrep -v '##' 1000GENOMES-phase_3.gvf.gz | cut -f9 | grep -o 'global_minor_allele_frequency=.*' | tr ';' '\t' | cut -f1 > 1000GENOMES_MAFonly.txt
zgrep 'global_minor_allele_frequency' 1000GENOMES-phase_3.gvf.gz | cut -f9 | grep -o 'rs[0-9]*' > 1000GENOMES_rsids.txt

paste 1000GENOMES_rsids.txt <(tr '=' '\t' < 1000GENOMES_MAFonly.txt | tr '|' '\t' | cut --complement -f1) > 1000GENOMES_MAF.txt

# Sort and join to variant file based on rsid

sort -k1,1 1000GENOMES_MAF.txt > 1000GENOMES_MAF.sorted.txt
sort -k12,12 df_v4_derivedfreqs.txt > df_v4_der.rsidsorted.txt

join -t $'\t' -1 12 -2 1 df_v4_der.rsidsorted.txt 1000GENOMES_MAF.sorted.txt > df_v4_dfreq_MAF.txt

# Paste together, reordering minor allele frequencies
paste <(cut -f2-12 df_v4_dfreq_MAF.txt) <(cut -f1 df_v4_dfreq_MAF.txt) <(cut -f13,17-21,23-25 df_v4_dfreq_MAF.txt) <(cut -f14-16 df_v4_dfreq_MAF.txt) > df_v4_dfreq_MAF.tmp
mv df_v4_dfreq_MAF.tmp df_v4_dfreq_MAF.txt

# Re-sort based on allele age
sort -g -k8,8 df_v4_dfreq_MAF.txt > df_v4_dfreq_MAF.agesorted.txt

# Filter out biallelic sites
grep -v ',' df_v4_dfreq_MAF.agesorted.txt > df_v4_dfreq_MAF_biallelic.agesorted.txt
# Filter out variants where derived frequency is > 98%
awk '{if ($19 == 0) minor_allele = $3; else minor_allele = $4;
if(minor_allele == $6) DAF = $20; else DAF = 1 - $20; if(DAF <= 0.98) print $0}' df_v4_dfreq_MAF_biallelic.agesorted.txt > df_v4_dfreq_MAF_biallelic_DAF98.agesorted.txt

# Order of columns in output: 
# Chr, Pos, Ref, Alt, Derived, Ancestral, Source, Allele_Age, triplet_derived, triplet_ancestral, rec_rate, rsid, derived,
# AFR_freq, AMR_freq, EAS_freq, EUR_freq, SAS_freq, minor_flag, MAF, allele_counts, Neanderthal_maskflag, AgeCI95Lower_Jnt, AgeCI95Upper_Jnt