##filter and construct input files
python3 fais_baypass.py -r 0.98 -q 0.99 -d 3 -L 3 -M 30 -n 4 ../short-header-Tribolium_castaneum.fa bam_list > fais_baypass.log

##thin loci, 100 bp apart, keeping sites with max variation
Rscript filter_baypass.R

##run core baypass model to produce 
./g_baypass -npop 96 -gfile  baypass_r0.98_d3_L3_M30_q0.99_a35_thin100.txt \
        -poolsizefile tribolium.poolsize \
         -nthreads 4 -outprefix core_model > core_model.log

#run aux model using teh covariance matrix made in the previous command
./g_baypass -npop 96 -gfile baypass_r0.98_d3_L3_M30_q0.99_a35_thin100.txt \
	-poolsizefile tribolium.poolsize -efile tribolium.covariates \
	-omegafile core_model_mat_omega.out \
	-d0yij 8 -auxmodel -isingbeta 0.9 -nthreads 4 -outprefix aux_model > aux_model.log
