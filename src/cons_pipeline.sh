#modify format to include index, which matches MRK column in baypass output files
awk '{print NR "\t" $0}' ../baypass2/vep_baypass_r0.98_d3_L3_M30_q0.99_a35_thin100.txt | tr " " "\t" > vep_baypass_r0.98_d3_L3_M30_q0.99_a35_thin100.sites

#created bed file from the list of sites used in baypass
awk '{print $2 "\t" ($3-1) "\t" $4}' vep_baypass_r0.98_d3_L3_M30_q0.99_a35_thin100.sites > vep_baypass_r0.98_d3_L3_M30_q0.99_a35_thin100.bed

#run mummer -- all genomes against Tribolium reference
#outputs a bed file for all unique exact matches for each reference
bash run_mummer.sh

#get intersection between all reference genomes and sites used in baypass
#awk command converts multi-locus intervals into a single nucleotide position per file line
multiIntersectBed -i vep_baypass_r0.98_d3_L3_M30_q0.99_a35_thin100.bed *_mummer.bed | awk '$4 == 17' | awk '{for (i=1; i<=($3 - $2); i++) print $1 "\t" ($2 +i - 1) "\t" ($2 + i )}' > conserved.bed

#similar to grep -Ff, subset the baypass sites for those found to be conserved in the steps above.
awk 'NR==FNR{A[$1,$3]=$3; next}{$3=A[$2,$3]; if($3!="") print}' conserved.bed vep_baypass_r0.98_d3_L3_M30_q0.99_a35_thin100.sites > conserved_mrk.txt


#produce baypass file with only conserved sites
python get_conserved_pij.py
