import argparse
import os
import re
import io
import subprocess
import numpy as np
import re
from joblib import Parallel, delayed
import multiprocessing


parser = argparse.ArgumentParser(description="Produce baypass format from sorted bams with various filtering specifications")

parser.add_argument("-r", "--ref_freq", type = float,
		    help = "Float between 0.5 and 1.0 to define upper bounder filter of the global reference allele frequency.")

parser.add_argument("-q", "--qual_prob", type = float, default = 0.99,
		    help = "Float between 0 and 1.0 to filter bases by quality (defaults to 0.99)")

parser.add_argument("-a", "--phred_offset", type = int, default = 35,
		    help = "Ascii encoding offset used to calculate base probabilities (defaults to 35, almost always the right choice now days).")

parser.add_argument("-d", "--minor_count", type = int, 
	            help = "Positive integer for the minimum number of times the minor allele must be seen.")

parser.add_argument("-L", "--min_pop", type = int, 
	            help = "Positive integer for the minimum number of alleles counted in every population (sites with any failing pop are excluded).")

parser.add_argument("-M", "--max_pop", type = int,
	            help = "Positive integer for the maximum number of alleles counted in every population (sites with any failing pop are excluded).")

parser.add_argument('-l', '--positions_file', type=str, 
		    help='Optional file of reference position to pass to samtools mpileup.')

parser.add_argument('-p', '--prefix', type=str, 
		    help='prefix for naming output files (meta data and baypass)')

parser.add_argument("-n", "--num_cores", type = int, metavar = "int cores", default=1, 
		    help = "Number of cores to use. If multiple cores are not available, do not use flag, or set int cores to zero")

parser.add_argument("reference", type=str, help="The indexed reference genome each bam was aligned to")

parser.add_argument("bam_list", type=str, help="file containing path to bam file on each line")

args = parser.parse_args()
qual_prob = args.qual_prob

nucs = np.array(["A", "T", "G", "C", "N", "a", "t", "g", "c", "n"])


suffix = "r{0}_d{1}_L{2}_M{3}_q{4}_a{5}.txt".format(args.ref_freq, args.minor_count, args.min_pop, args.max_pop, args.qual_prob, args.phred_offset)
#baypass_str = "baypass_" + suffix
vep_str = "vep_baypass_" + suffix

#baypass = open(baypass_str, "w")
vep = open(vep_str, "w")

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

bam_count = file_len(args.bam_list)
#print("there are", bam_count, " bams")


def base_prob(phred, offset = 35):
	return np.array([1 - 1/np.power(10, (ord(asci) - offset)/10) for asci in phred])

def parse_seq(seq_str, qual_str, ref):
    
	seq_str = re.sub('\^.', '', seq_str.replace(".", ref).replace(",", ref).replace("$", "").upper())
    
	broken_seq = re.split('[+-]\d+', seq_str)
    
	seq_array =  np.array(list(broken_seq[0] + ''.join(
	[
	broken_seq[i][int(re.findall('\d+', seq_str)[i-1]):]
		for i in list(range(1, len(broken_seq)))
	]
	)
	))

	qs = base_prob(qual_str, offset = 35)
	seq_filtered = ''.join(seq_array[qs >= qual_prob])
	nucs = ["A", "T", "G", "C"]
	return [np.asscalar(np.char.count(seq_filtered, nuc)) for nuc in nucs]



def make_baypass(locus):
	bi_site = list(set(np.where(locus > 0)[1]))
	if len(bi_site) == 2:
		bi_pair = [" ".join((i[bi_site]).astype(str)) for i in locus]
		return " ".join(bi_pair)
	else:
		raise ValueError("This site is not bi-alleleic. Check your code!")


def parse_site(pileup):

	#get data parsed
	mp = pileup.strip().split("\t")
	chrom, pos, ref = mp[0:3] #site data
	pop_bam = mp[3:]
	
	#set up indices in groups of 3
	idx = list(range(0,len(pop_bam)-1, 3))

	#process read strings
	dp_array =  np.array([int(pop_bam[i]) for i in idx])

	dp_test = np.min(dp_array) >= args.min_pop and np.max(dp_array) <= args.max_pop

	if dp_test:


		read_list = [pop_bam[i+1] for i in idx]
		read_str = ' '.join(read_list).replace("*", "")
		ref_freq = (read_str.count(",") + read_str.count(".")) / np.sum(dp_array)
	
		freq_test = (1 - args.ref_freq) <= ref_freq <= args.ref_freq
	
		if freq_test:

			n_pop = list(range(len(read_list)))
			phred_list = [pop_bam[i+2] for i in idx]
			
			locus = np.array(Parallel(n_jobs = args.num_cores)(delayed(parse_seq)(read_list[i], phred_list[i], ref) for i in n_pop))
			
			locus_sums = locus.sum(axis = 1)
			locus_max_sums = locus_sums.max()
			locus_min_sums = locus_sums.min()
			dp_test2 = locus_min_sums >= args.min_pop and locus_max_sums <= args.max_pop

			if dp_test2:
			#locus = np.array([parse_seq(read_list[i], phred_list[i], ref) for i in n_pop])
			
			#print(locus)

				al_counts = np.sum(locus, axis = 0)
				bi_test = len(np.where(al_counts > 0)[0]) == 2
				if bi_test:

					non_zero = np.where(al_counts > 0)
					count_test = np.min(al_counts[non_zero]) >= args.minor_count
					if count_test:	
						l_baypass = make_baypass(locus)
						#print(len(l_baypass.split())
						if len(l_baypass.split()) == int(2*bam_count):
							alleles = nucs[list(set(np.where(locus > 0)[1]))]
							alt = "".join(alleles[alleles != ref].astype(str))
							print(chrom, pos, pos, ref + "/" + alt, ref_freq, l_baypass, file = vep)
							
							#l_baypass = make_baypass(locus)
							#print(l_baypass, file = baypass)
							#print("dp", np.min(dp_array), np.max(dp_array))
							#print("qual", locus_min_sums, locus_max_sums)
							#print("site", l_baypass)
							#print("site", len(l_baypass.split()))
	
							if alt == ref:
								raise ValueError("Invariant site! Abort!")



#construct mpileup command with and without positions argument
if args.positions_file:
	cmd = "samtools mpileup -f {0} -l {1} -b {2}".format(args.reference, args.positions_file, args.bam_list).split()
else:
	cmd = "samtools mpileup -f {0} -b {1}".format(args.reference, args.bam_list).split()

proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
for line in io.TextIOWrapper(proc.stdout, encoding="utf-8"):
	parse_site(line)
	

########
#legacy#
########

#def parse_read(reads, phreds, qual_min = 0):
#
#	ref_split = re.findall('\d+', reads)
#	broken_seq = re.split('[+-]\d+', reads)
#	reads = broken_seq[0] + ''.join([broken_seq[i][int(ref_split[i-1]):] for i in list(range(1, len(broken_seq)))])
#
#	depth = len(reads)
#	nuc_dict = {"A": 0, "T": 0, "G": 0, "C":0, "N":0}
#
#	qual_pos = 0
#
#	i_c = 0
#	if depth < 1:
#		return [0, 0, 0, 0]
#	else:
#		while i_c < depth:
#			base = reads[i_c]
#			if base == "^":
#				i_c += 2
#			base = reads[i_c]
#			if base in nucs:
#				if base_prob(phreds[qual_pos]) >= qual_min:
#					nuc_dict[base] += 1
#					qual_pos += 1
#				i_c +=1
#			if base not in nucs:
#				i_c +=1
#
#		nuc_count = [nuc_dict["A"], nuc_dict["T"], nuc_dict["G"], nuc_dict["C"]]
#		return nuc_count

	#print( min(dp_array) + 100)
	#print(args.max_pop)

	#locus = list(range(0, len(idx)))
	
	#print(parse_read(''.join(read_list), ''.join(phred_list), ref))	

	#for l in locus:
	#	print(read_list[l], phred_list[l])


	#locus_count = np.array([parse_read(read_list[l], phred_list[l]) for l in locus])
	
	
	#bi_allele = locus_test(locus_count, dp_min = args.min_pop, dp_max = args.max_pop, count_min = args.minor_count, freq_max = args.maj_freq)
	

	#if bi_allele is not False:

	#	alleles = nucs[list(set(np.where(locus_count > 0)[1]))]
	#	alt = "".join(alleles[alleles != ref].astype(str))
	#	ref_freq = str(np.array(bi_allele["freq"])[np.where(nucs == ref)][0])
		
	#	print(make_baypass(locus_count), file = baypass)
	#	print(chrom, pos, pos, ref + "/" + alt, ref_freq, bi_allele["counts"], file = vep, flush = True)
		

#def base_prob(ascii, offset = 35):
#	return 1 - 1/np.power(10, (ord(ascii) - offset)/10)

#def locus_test(locus, dp_min, dp_max, count_min, freq_max):
#	if locus.sum() > 0:
#		al_counts = np.sum(locus, axis = 0)
#		al_freq = al_counts / np.sum(locus)
#		non_zero = np.where(al_counts > 0)
#		pop_dp = np.sum(locus, axis = 1)
#    
#		count_test = np.min(al_counts[non_zero]) >= count_min
#		freq_test = np.max(al_freq) <= freq_max
#		dp_test = np.logical_and(pop_dp.min() >= dp_min, pop_dp.max() <= dp_max).all()
#		bi_test = len(np.where(al_counts > 0)[0]) == 2
#
#		if count_test and freq_test and dp_test and bi_test:
#			return {"freq":al_freq, "counts": al_counts}
#		else:
#			return False
#	else:
#		return False


