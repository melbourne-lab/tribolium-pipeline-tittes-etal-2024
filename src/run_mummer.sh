for i in *fna;
do
     #find exact matches of 18+, sort by chrom name and position
	mummer -mum -load idx -F -b -l 18 -L -n ../short-header-Tribolium_castaneum.fa ${i} | awk '{print $1 "\t" ($2 - 1) "\t" ($2 + $4 - 1)}' | grep -v ">" | sort -k1,1d -k2,2n -k 3,3n -u > ${i}_mummer.bed

done
