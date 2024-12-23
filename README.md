# rpgQTL
Regions per gene (rpg) QTL

This is a package based on [TensorQTL](https://github.com/broadinstitute/tensorqtl) with modifications (cis-eQTL module only).  
Instead of a single parameter of the cis-window size, rpgQTL allows user to specify regions (could be discontinuous) for each gene. Only SNPs in the corresponding regions will be used for corresponding gene cis-eQTL calling.

### Installation and dependancy
See detailed tutorial [here](https://github.com/GenEvoBioinfLab/rpgQTL/blob/main/install/README.md).  

### Using rpgQTL in command line
```
python run.py [expression_file] [covariates_file] [genotype_file] [variant_file] \
    [rpg_file] [gene_set_file] [output_dir] \
    [l_window] [s_window] [RType] [EType] [seed]
```

##### input parameters
1. `expression_file`: Expression file in .bed.gz.
2. `covariates_file`: Covariates (e.g. PEER factor, sex, etc.) file in txt.
3. `genotype_file`: Genotype file in .txt.gz.
4. `variant_file`:
For expression and covatiates, you can download from GTEx v8 data from the GTEx portal.  
For genotype, see instructions in GTEx to access the controlled data.  
Also see instruction in [tensorqtl](https://github.com/broadinstitute/tensorqtl).  
See example/input for file format example.  
5. `rpg_file`: Each row correspond to one candidate genomic region for eqtl detection of one gene. Each gene could contain multiple regions (rows). The required first four columns should be: chromosome, start of the region, ends of the region, and gene name.  
6. `gene_set_file`: If given and if EType == "gene_set", only genes in this file will be used to calculate eQTLs. Each row should be one gene name. Default: None.
7. `output_dir`: Output folder.

##### rpgQTL parameters
8. `l_window` (int; default 2,000,000):  
The max boundary of cis-window. Only regions within this cis-window will be consider as valid regions. You should consider trans-eqtl calling methods if your interested SNPs are located outside of this window size.  
9. `s_window` (int; default 0):  
A fixed cis-window for all genes where SNPs in these regions will always be included, even if they are not included in the rpg_file. Setting s_window to 1,000,000 and rpg_file to empty will be equivalent to the original cis-eqtl calling using tensorqtl with a window size of 1,000,000.
10. `RType` (one of ['remove', 'l_window', 's_window'], default 'remove'):  
Different ways to handle genes that have genotype and expression data, but no candidate regions.  
  - 'remove': The genes will be completely removed from calculation.  
  - 'l_window': All SNPs within the cis-window defined by 'l_window' will be used for such genes.  
  - 's_window': All SNPs within the cis-window defined by 's_window' will be used for such genes.  
11. `EType` (one of ['default', '1M', 'gene_set'], default 'default'):
  - 'default': Call independent eQTLs using FDR<=0.05.
  - '1M': Firstly, we define eGenes as default. Secondly, using the permutation result (regardless of FDR) with s_window=1,000,000 (will do the calculation if not yet done) and keep the genes in eGenes. Finally, call independent eQTLs using this 1Mb permutation result with the defined eGenes.
  - 'gene_set': Using the given gene set as eGenes. Call independent eQTLs on this gene set, regardless of FDR.

### Using rpgQTL as package
`rpgQTL.run_nominal`: similar to `tensorqtl.cis.map_nominal`. Conduct the nominal run to get nominal p-values.  
`rpgQTL.run_permutation`: similar to `tensorqtl.cis.map_cis`. Conduct the permutation run to get beta-adjusted p-values.  
`rpgQTL.run_independent`: similar to `tensorqtl.cis.map_independent`. Conduct the forward-backward run to get independent eQTLs.  
Other parameters in rpgQTL functions are the same as those in the corresponding functions from `tensorqtl`. Notice the rpgQTL only support for the basic cis-eQTL calculations, and thus some tensorQTL parameters are not supported (e.g. all parameters related to "interaction" from tensorqtl.cis.map_nominal are not supported in rpgQTL.run_nominal).
  
See `example.ipynb` for the example script.
