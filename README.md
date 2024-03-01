# rpgQTL
Regions per gene (rpg) QTL

This is a package based on [TensorQTL](https://github.com/broadinstitute/tensorqtl) with modifications (cis-eQTL module only).  
Instead of a single parameter of the cis-window size, rpgQTL allows user to specify regions (could be discontinuous) for each gene. Only SNPs in the corresponding regions will be used for corresponding gene cis-eQTL calling.

### Install
Install directly from this repository:
```
$ git clone https://github.com/gersteinlab/rpgQTL.git
$ cd rpgQTL
$ pip install -r install/requirements.txt .
```
A detailed tutorial of setting up the environment from scratch could be found [here](https://github.com/GenEvoBioinfLab/rpgQTL/blob/main/install/README.md).

### Requirement
```
numpy
pandas
tensorqtl
```
Notice that tensorqtl depends on pandas-plink, pytorch and other packages. For furthur details, see the installing instruction in [tensorqtl](https://github.com/broadinstitute/tensorqtl) or [pytorch](https://pytorch.org/get-started/locally/).

### Example
See `example.ipynb` for the example script. The following is based on the example but applied to the general usage of the rpgQTL package.

### Input files
1. `plink_prefix_path`: Genotype VCF in PLINK format  
2. `phenotype_df, phenotype_pos_df`: expression file  
3. `covariates_df`: covariates (e.g. PEER factor, sex, etc.)  
For expression and covatiates, you can download from GTEx v8 data from the GTEx portal. For genotype, see instructions in GTEx to access the controlled data.
Also see instruction in [tensorqtl](https://github.com/broadinstitute/tensorqtl).  
4. `rpg_df`: `pandas.DataFrame` or `str`  
If `pandas.DataFrame`, this would be one large dataframe. Each row correspond to one candidate genomic region for eqtl detection of one gene. Each gene could contain multiple regions (rows). The required first four columns should be: chromosome, start of the region, ends of the region, and gene name  
If `str`, this would be the path to a directory. For each file in the directory, the file name is a gene name and the file contains the regions for only that gene. The format of each file is the same as the above pandas.DataFrame.

### rpgQTL functions
`rpgQTL.run_nominal`: similar to `tensorqtl.cis.map_nominal`. Conduct the nominal run to get nominal p-values.
`rpgQTL.run_permutation`: similar to `tensorqtl.cis.map_cis`. Conduct the permutation run to get beta-adjusted p-values.
`rpgQTL.run_independent`: similar to `tensorqtl.cis.map_independent`. Conduct the forward-backward run to get independent eQTLs.

### rpgQTL Parameters
- `l_window` (int; default 2000000):  
The max boundary of cis-window. Only regions within this cis-window will be consider as valid regions. You should consider trans-eqtl calling methods if your interested SNPs are located outside of this window size.  
- `s_window` (int; default 0):  
A fix cis-window for all genes that will always be included, even if they are not included in the rpg_df. Setting s_window to 1000000 and rpg_df to empty dataframe will be equivalent to the original cis-eqtl calling using tensorqtl with window=1000000.  
- `NonHiCType` ('remove', 'l_window' or 's_window'):  
Different ways to deal with genes that have genotype and expression data, but no candidate regions.
  - 'remove': The genes will be completely removed from calculation.  
  - 'l_window': All SNPs within the cis-window defined by 'l_window' will be used.  
  - 's_window': All SNPs within the cis-window defined by 's_window' will be used.  
- `eGene_list` (list of str; default None):
For `run_independent` only. When given, only independent eQTLs for this gene list will be calculated. By default, this is set to None, so that all genes filter by the FDR threshold are used.

Other parameters in rpgQTL functions are the same as those in the corresponding functions from `tensorqtl`. Notice the rpgQTL only support for the basic cis-eQTL calculations, and thus some tensorQTL parameters are not supported (e.g. all parameters related to "interaction" from tensorqtl.cis.map_nominal are not supported in rpgQTL.run_nominal).
