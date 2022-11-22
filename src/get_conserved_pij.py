#file with conserved loci
conserved = "conserved_mrk.txt"
cons_mrk = 0 #column with marker numbers

#file with all allele frequencies modeled by aux
#OLD all_pij= "../tribolium_AllPops_autosomes_auxModel_summary_yij_pij.out"
all_pij= "../baypass2/aux_model_summary_yij_pij.out"
pij_col = 1 #column with marker numbers

# OLD out_file = open("../tribolium_AllPops_autosomes_auxModel_summary_yij_pij_CONSERVED.out", "w")
out_file = open("aux_model_summary_yij_pij_CONSERVED.out", "w")
print("POP MRK M_Y SD_Y M_P SD_P M_Pstd SD_Pstd DELTA_Y ACC_Y DELTA_P ACC_P", file = out_file)

cons_dict = dict()

with open(conserved) as cons:
	for line in cons:
		mark = line.strip().split()[cons_mrk]
		if mark not in cons_dict:
			cons_dict[mark] = "" #key with empty value


with open(all_pij) as pij:
	for line in pij:
		mrk = line.strip().split()[pij_col]
		if mrk in cons_dict:
			print(line.strip(), file = out_file)
