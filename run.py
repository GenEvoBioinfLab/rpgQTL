import os
import argparse
import pandas as pd
from tensorqtl import tensorqtl, cis, post, genotypeio
import rpgQTL

def main():

    ## parameters
    parser = argparse.ArgumentParser(description="Parameters of rpgQTL")
    parser.add_argument("expression_file", help="Input expression .bed.gz file.")
    parser.add_argument("covariates_file", help="Input covariates .txt file.")
    parser.add_argument("genotype_file", help="Input genotype .txt.gz file.")
    parser.add_argument("variant_file", help="Input variant .txt.gz file.")
    parser.add_argument("rpg_file", help="Input RPG .bed file.")
    parser.add_argument("gene_set_file", default=None, help="Input gene set file.")
    parser.add_argument("output_dir", help="Output folder")
    parser.add_argument("l_window", type=int, default=2000000, help="Boundary for cis-region.")
    parser.add_argument("s_window", type=int, default=0, help="Cis-region where any SNPs within this range are always used.")
    parser.add_argument("RType", choices=["remove", "l_window", "s_window"], default='remove', help="How to deal with SNPs not in the given RPG file.")
    parser.add_argument("Etype", choices=["default", "1M", "gene_set"], default='default', help="How to calculate independent eQTL.")
    parser.add_argument("seed", type=int, default=123456, help="Random seed")
    args = parser.parse_args()

    ## output folder and file name
    os.makedirs(args.result_dir, exist_ok=True)
    prefix = "RPG.%d." % (args.s_window)

    ## load data
    print("Loading input...")
    phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(args.expression_file)
    covariates_df = pd.read_csv(args.covariates_file, sep='\t', index_col=0).T
    genotype_df = pd.read_csv(args.genotype_file, sep="\t", index_col=0)
    variant_df = pd.read_csv(args.variant_file, sep="\t", index_col=0)
    rpg_df = pd.read_csv(args.rpg_file, sep="\t", header=None)

    ## nominal run
    rpgQTL.run_nominal(genotype_df, variant_df, phenotype_df, phenotype_pos_df, covariates_df,
        rpg_df, l_window=args.l_window, s_window=args.s_window, RType=args.RType,
        output_dir=args.output_dir, prefix=prefix)

    ## permutation run
    egenes_df = rpgQTL.run_permutation(genotype_df, variant_df, phenotype_df, phenotype_pos_df, covariates_df, 
        rpg_df, l_window=args.l_window, s_window=args.s_window, RType=args.RType, seed=args.seed)

    ## calculate q-values
    post.calculate_qvalues(egenes_df, fdr=0.05, qvalue_lambda=0.85)
    egenes_df.to_csv("%s/%s.egenes.tsv" % (args.output_dir, prefix), index=True)

    ## significant pairs
    pairs_df = post.get_significant_pairs(egenes_df, "%s/%s" % (args.output_dir, prefix))
    pairs_df.to_csv("%s/%s.sig_pairs.tsv" % (args.output_dir, prefix), index=False)

    ## independent eQTL
    if args.Etype == "default":
        ## By default, calculate the independent eQTL using FDR<=0.05 from permutation run.
        indep_df = rpgQTL.run_independent(genotype_df, variant_df, egenes_df, phenotype_df, phenotype_pos_df, covariates_df, 
            rpg_df, l_window=args.l_window, s_window=args.s_window, RType=args.RType, seed=args.seed)
        indep_df.to_csv("%s/%s.indep.tsv" % (args.output_dir, prefix), index=False)

    elif args.Etype == "1M":
        ## Define eGenes as those in the current s_window Mb runs,
        ## then use the results from 1Mb runs, keep those 'eGenes' only,
        ## call independent eQTLs from these 'eGenes' using 1Mb result
        egenes_df1 = pd.read_csv("%s/%s.egenes.tsv" % (args.output_dir, prefix), index_col=0)

        if os.path.exists("%s/RPG.1000000.egenes.tsv" % (args.output_dir)):  
            print("Etype set to 1M. eGenes result from 1M found. Loading...")
            egenes_df2 = pd.read_csv("%s/RPG.1000000.egenes.tsv" % (args.output_dir), index_col=0)
        else:
            print("Etype set to 1M. Calculating the eGenes from 1M")
            rpgQTL.run_nominal(genotype_df, variant_df, phenotype_df, phenotype_pos_df, covariates_df,
                rpg_df, l_window=args.l_window, s_window=1000000, RType=args.RType,
                output_dir=args.output_dir, prefix="RPG.1000000.")
            egenes_df2 = rpgQTL.run_permutation(genotype_df, variant_df, phenotype_df, phenotype_pos_df, covariates_df, 
                rpg_df, l_window=args.l_window, s_window=1000000, RType=args.RType, seed=args.seed)
            post.calculate_qvalues(egenes_df2, fdr=0.05, qvalue_lambda=0.85)
            egenes_df2.to_csv("%s/RPG.1000000.egenes.tsv" % (args.output_dir), index=True)
            egenes_df2 = pd.read_csv("%s/RPG.1000000.egenes.tsv" % (args.output_dir), index_col=0)       

        egenes_df = egenes_df2[egenes_df2.index.isin(egenes_df1[egenes_df1['qval'] < 0.05].index)].copy()

        ## fdr=1 to include all genes after filtering
        indep_df = rpgQTL.map_independent(genotype_df, variant_df, egenes_df, phenotype_df, phenotype_pos_df, covariates_df,
            rpg_df, l_window=args.l_window, s_window=1000000, seed=args.seed, fdr=1)
        indep_df.to_csv("%s/%s.indep1M.tsv" % (args.output_dir, prefix), index=False)

    elif args.Etype == "gene_set":
        ## Use the given gene list to call indep eqtl, regardless of FDR
        eGene_list = pd.read_csv(args.gene_set_file, header=None)[0].values
        indep_df = rpgQTL.map_independent(genotype_df, variant_df, egenes_df, phenotype_df, phenotype_pos_df, covariates_df,
            rpg_df, l_window=args.l_window, s_window=args.s_window, RType=args.RType, seed=args.seed, eGene_list=eGene_list)
        indep_df.to_csv("%s/%s.indepGeneSet.tsv" % (args.output_dir, prefix), index=False)

if __name__ == "__main__":
    main()
