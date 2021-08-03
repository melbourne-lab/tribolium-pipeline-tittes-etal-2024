rule all:
    input:
        "figures/figure2.pdf",
        "figures/figure2.png",
        "data_out/outlier_df.csv",
        "data_out/outlier_categories.csv",
        "figures/figure1.pdf",
        "figures/figure1.png",
        "data_out/core_freq_PCA.csv",
        "figures/figure3.pdf",
        "figures/figure3.png",
        "figures/figure4.pdf",
        "figures/figure4.png",
        "data_out/trait_outlier_table.csv",
        "data_out/outlier_conserved_df.csv",
        "figures/figure5.pdf",
        "figures/figure5.png",
        "data_out/conserved_categories.csv",
        "data_out/highimpact_categories.csv",
        "data_out/allsites_categories.csv"


helpers = ["data/Z.csv", "data/conserved_sites/vep_baypass_r0.98_d3_L3_M30_q0.99_a35_thin100.sites", "data/VEP_allsites.vcf", ]

rule figure2:
    input:
        helpers,
        "src/FIGURE2_baypass_BFmc.R", 
        "src/helper_functions.R",
        "data/baypass2/aux_model_summary_betai.out",
        "data/baypass2/vep_baypass_r0.98_d3_L3_M30_q0.99_a35_thin100.txt",
    output:
        "figures/figure2.pdf",
        "figures/figure2.png",
        "data_out/outlier_df.csv",
        "data_out/outlier.tsv",
        "data_out/outlier_categories.csv"
    shell:
        "Rscript src/FIGURE2_baypass_BFmc.R"


rule figure1:
    input:
        helpers,
        "src/FIGURE1_pca_viz.R",
        "src/helper_functions.R",
        "data/baypass2/core_model_mat_omega.out",
        "data/baypass2/core_model_summary_yij_pij.out",
        "data_out/outlier_df.csv"
    output:
        "figures/figure1.pdf",
        "figures/figure1.png",
        "data_out/core_freq_PCA.csv"
    shell:
        "Rscript src/FIGURE1_pca_viz.R"

rule figure3and4:
    input:
        helpers,
        "src/helper_functions.R",
        "src/make_traitdata.R",
        "src/FIGURE_3_4_trait_association.R",
        "data_out/outlier_df.csv",
        "data_out/core_freq_PCA.csv",
        "data/baypass2/vep_baypass_r0.98_d3_L3_M30_q0.99_a35_thin100.txt",
        "data/baypass2/aux_model_summary_yij_pij.out"
    output:
        "figures/figure3.pdf",
        "figures/figure3.png",
        "figures/figure4.pdf",
        "figures/figure4.png",
        "data_out/trait_outlier_table.csv"
    shell:
        "Rscript src/FIGURE_3_4_trait_association.R"


rule figure5:
    input:
        helpers,
        "src/FIGURE5_conserved_sites.R",
        "src/helper_functions.R",
        "src/make_traitdata.R",
        "data_out/outlier_df.csv",
        "data/conserved_sites/aux_model_summary_yij_pij_CONSERVED.out",
    output:
        "data_out/outlier_conserved_df.csv",
        "figures/figure5.pdf",
        "figures/figure5.png",
        "data_out/conserved_categories.csv",
        "data_out/highimpact_categories.csv"
    shell:
        "Rscript src/FIGURE5_conserved_sites.R"


rule categories:
    input:
        "data_out/outlier_categories.csv",
        "data_out/conserved_categories.csv",
        "data_out/highimpact_categories.csv"
    output:
        "data_out/allsites_categories.csv"
    shell:
        "Rscript src/TABLE2_allsites_catogries.R"
