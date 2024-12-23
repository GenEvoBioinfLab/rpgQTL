import os
import time
import torch
import numpy as np
import pandas as pd
import scipy.stats as stats
from collections import OrderedDict
from tensorqtl import tensorqtl, core, cis, post, genotypeio

output_dtype_dict = {
        'num_var':np.int32,
        'beta_shape1':np.float32,
        'beta_shape2':np.float32,
        'true_df':np.float32,
        'pval_true_df':np.float64,
        'variant_id':str,
        'tss_distance':np.int32,
        'ma_samples':np.int32,
        'ma_count':np.int32,
        'maf':np.float32,
        'ref_factor':np.int32,
        'pval_nominal':np.float64,
        'slope':np.float32,
        'slope_se':np.float32,
        'pval_perm':np.float64,
        'pval_beta':np.float64}

class background:
    def __init__(self, max_prefetch=10):
        self.max_prefetch = max_prefetch
    def __call__(self,gen):
        def bg_generator(*args,**kwargs):
            return genotypeio.BackgroundGenerator(gen(*args,**kwargs), max_prefetch=self.max_prefetch)
        return bg_generator


class rpg_igc(object):
    def __init__(self, genotype_df, variant_df, phenotype_df, phenotype_pos_df,
                 rpg_df, l_window=1000000, s_window=300000, NonHiCType='s_window'):
        
        assert (genotype_df.index==variant_df.index).all()
        assert (phenotype_df.index==phenotype_df.index.unique()).all()
        self.genotype_df = genotype_df
        self.variant_df = variant_df.copy()
        self.variant_df['index'] = np.arange(variant_df.shape[0])
        n_samples = phenotype_df.shape[1]
        # check for constant phenotypes and drop
        self.phenotype_df = phenotype_df.copy()
        self.phenotype_pos_df = phenotype_pos_df.copy()        
        m = np.all(phenotype_df.values == phenotype_df.values[:,[0]], 1)
        if m.any():
            print(f'    ** dropping {np.sum(m)} constant phenotypes')
            self.phenotype_df = phenotype_df.loc[~m]
            self.phenotype_pos_df = phenotype_pos_df.loc[~m]

        self.phenotype_tss = self.phenotype_pos_df['tss'].to_dict()
        phenotype_chr = self.phenotype_pos_df['chr'].to_dict()
        variant_chrs = self.variant_df['chrom'].unique()
        phenotype_chrs = self.phenotype_pos_df['chr'].unique()
        self.chrs = [i for i in phenotype_chrs if i in variant_chrs]
        chr_variant_dfs = {c:g[['pos', 'index']] for c,g in self.variant_df.groupby('chrom')}
        
        valid_ix = []
        self.cis_ranges = {}
        for k,phenotype_id in enumerate(self.phenotype_df.index,1):
            if np.mod(k, 1000) == 0:
                print(f'\r  * checking phenotypes: {k}/{self.phenotype_df.shape[0]}', end='' if k != phenotype_df.shape[0] else None)

            tss = self.phenotype_tss[phenotype_id]
            chrom = phenotype_chr[phenotype_id]

            candidate_pos = chr_variant_dfs[chrom]['pos']
            candidate_index = chr_variant_dfs[chrom]['index']

            bool_lwindow = candidate_pos.between(tss - l_window, tss + l_window)
            bool_swindow = candidate_pos.between(tss - s_window, tss + s_window)


            if isinstance(rpg_df, pd.DataFrame):
                rpg_phenotypes = rpg_df[3].values
            elif isinstance(rpg_df, str):
                rpg_phenotypes = np.array(os.listdir(rpg_df))
            else:
                raise Exception('not valid rpg_df')

            if phenotype_id in rpg_phenotypes:
                # for genes have input candidate regions
                
                if isinstance(rpg_df, pd.DataFrame):
                    phenotype_rpg_coor = rpg_df[rpg_df[3] == phenotype_id].iloc[:,[1,2]].values
                elif isinstance(rpg_df, str):
                    phenotype_rpg_coor = pd.read_csv(rpg_df+"/"+phenotype_id, sep="\t", header=None).iloc[:,[1,2]].values
                    
                bool_rpg = pd.Series([False]*len(bool_lwindow), index=bool_lwindow.index) # init with all False
                for rpg_coor in phenotype_rpg_coor:
                    cond = candidate_pos.between(rpg_coor[0], rpg_coor[1])
                    bool_rpg = bool_rpg|cond # add True

                # use everything in swindow (TSS+-nbp, strict promoter)
                # use everything in RPG (all input candidate regions, e.g. enhancer or 'wide' promoter)
                # should be within lwindow, remove snp beyond this max boundary (cis vs trans-eQTL boundary)
                r = candidate_index[(bool_swindow|bool_rpg)&bool_lwindow].values

            else:
                # for genes with no input               
                # use swindow only
                if NonHiCType == 's_window':
                    r = candidate_index[bool_swindow].values
                # use lwindow only
                elif NonHiCType == 'l_window':
                    r = candidate_index[bool_lwindow].values
                # remove such gene
                elif NonHiCType == 'remove':
                    r = np.array([])

            if len(r) > 0:
                valid_ix.append(phenotype_id)
                self.cis_ranges[phenotype_id] = r
        
        if len(valid_ix) != self.phenotype_df.shape[0]:
            print('    ** dropping {} phenotypes without variants'.format(self.phenotype_df.shape[0] - len(valid_ix)))
            self.phenotype_df = self.phenotype_df.loc[valid_ix]
            self.phenotype_pos_df = self.phenotype_pos_df.loc[valid_ix]
            self.phenotype_tss = self.phenotype_pos_df['tss'].to_dict()
            phenotype_chr = self.phenotype_pos_df['chr'].to_dict()
        self.n_phenotypes = self.phenotype_df.shape[0]

    @background(max_prefetch=6)
    def generate_data(self, chrom=None, verbose=False):
        chr_offset = 0
        if chrom is None:
            phenotype_ids = self.phenotype_df.index
            chr_offset = 0
        else:
            phenotype_ids = self.phenotype_pos_df[self.phenotype_pos_df['chr']==chrom].index
            offset_dict = {i:j for i,j in zip(*np.unique(self.phenotype_pos_df['chr'], return_index=True))}
            if chrom in offset_dict.keys():
                chr_offset = offset_dict[chrom]
            else:
                (nulla, nullb, nullc, nulld) = (None, None, None, None)
                yield nulla,nullb,nullc,nulld

        index_dict = {j:i for i,j in enumerate(self.phenotype_df.index)}

        for k,phenotype_id in enumerate(phenotype_ids, chr_offset+1):
            if verbose:
                genotypeio.print_progress(k, self.n_phenotypes, 'phenotype')
            p = self.phenotype_df.values[index_dict[phenotype_id]]
            r = self.cis_ranges[phenotype_id]
            yield p, self.genotype_df.iloc[r,].values, r, phenotype_id


def run_nominal(genotype_df, variant_df, phenotype_df, phenotype_pos_df, covariates_df,
                    rpg_df, l_window=1000000, s_window=300000, NonHiCType='s_window',
                    output_dir='.', prefix='prefix', write_stats=True, logger=None, verbose=True):

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    if logger is None:
        logger = core.SimpleLogger()

    logger.write('cis-QTL mapping: nominal associations for all variant-phenotype pairs')
    logger.write('  * {} samples'.format(phenotype_df.shape[1]))
    logger.write('  * {} phenotypes'.format(phenotype_df.shape[0]))
    
    assert np.all(phenotype_df.columns==covariates_df.index)
    logger.write('  * {} covariates'.format(covariates_df.shape[1]))
    residualizer = core.Residualizer(torch.tensor(covariates_df.values, dtype=torch.float32).to(device))
    dof = phenotype_df.shape[1] - 2 - covariates_df.shape[1]
    logger.write('  * {} variants'.format(variant_df.shape[0]))
    
    genotype_ix = np.array([genotype_df.columns.tolist().index(i) for i in phenotype_df.columns])
    genotype_ix_t = torch.from_numpy(genotype_ix).to(device)

    igc = rpg_igc(genotype_df, variant_df, phenotype_df, phenotype_pos_df,
                                rpg_df, l_window, s_window, NonHiCType)
    best_assoc = []
    start_time = time.time()
    k = 0
    logger.write('  * Computing associations')
    
    for chrom in igc.chrs:
        logger.write('    Mapping chromosome {}'.format(chrom))
        n = 0
        for i in igc.phenotype_pos_df[igc.phenotype_pos_df['chr']==chrom].index:
            j = igc.cis_ranges[i]
            n += len(j)

        chr_res = OrderedDict()
        chr_res['phenotype_id'] = []
        chr_res['variant_id'] = []
        chr_res['tss_distance'] = np.empty(n, dtype=np.int32)
        chr_res['maf'] =          np.empty(n, dtype=np.float32)
        chr_res['ma_samples'] =   np.empty(n, dtype=np.int32)
        chr_res['ma_count'] =     np.empty(n, dtype=np.int32)
        chr_res['pval_nominal'] = np.empty(n, dtype=np.float64)
        chr_res['slope'] =        np.empty(n, dtype=np.float32)
        chr_res['slope_se'] =     np.empty(n, dtype=np.float32)

        start = 0
        for k, (phenotype, genotypes, genotype_range, phenotype_id) in enumerate(igc.generate_data(chrom=chrom, verbose=verbose), k+1):

            if phenotype is None:
                logger.write('Skipping chromosome {} with no valid phenotype'.format(chrom))
                break

            phenotype_t = torch.tensor(phenotype, dtype=torch.float).to(device)
            genotypes_t = torch.tensor(genotypes, dtype=torch.float).to(device)
            genotypes_t = genotypes_t[:,genotype_ix_t]
            core.impute_mean(genotypes_t)

            variant_ids = variant_df.iloc[genotype_range].index
            tss_distance = np.int32(variant_df.iloc[genotype_range]['pos'].values - igc.phenotype_tss[phenotype_id])

            res = cis.calculate_cis_nominal(genotypes_t, phenotype_t, residualizer=residualizer)
            tstat, slope, slope_se, maf, ma_samples, ma_count = [i.cpu().numpy() for i in res]
            n = len(variant_ids)

            if n > 0:
                chr_res['phenotype_id'].extend([phenotype_id]*n)
                chr_res['variant_id'].extend(variant_ids)
                chr_res['tss_distance'][start:start+n] = tss_distance
                chr_res['maf'][start:start+n] = maf
                chr_res['ma_samples'][start:start+n] = ma_samples
                chr_res['ma_count'][start:start+n] = ma_count
                chr_res['pval_nominal'][start:start+n] = tstat
                chr_res['slope'][start:start+n] = slope
                chr_res['slope_se'][start:start+n] = slope_se
            start += n
        
        # if no valid data, skip to next chromosome
        if phenotype is None:
            continue

        logger.write('    time elapsed: {:.2f} min'.format((time.time()-start_time)/60))

        if start < len(chr_res['maf']):
            for x in chr_res:
                chr_res[x] = chr_res[x][:start]

        if write_stats:
            chr_res_df = pd.DataFrame(chr_res)
            m = chr_res_df['pval_nominal'].notnull()
            chr_res_df.loc[m, 'pval_nominal'] = 2*stats.t.cdf(-chr_res_df.loc[m, 'pval_nominal'].abs(), dof)
            print('    * writing output')
            chr_res_df.to_parquet(os.path.join(output_dir, '{}.cis_qtl_pairs.{}.parquet'.format(prefix, chrom)))

    logger.write('done.')


def run_permutation(genotype_df, variant_df, phenotype_df, phenotype_pos_df, covariates_df,
                rpg_df, l_window=1000000, s_window=300000, NonHiCType='s_window',
                beta_approx=True, nperm=10000, logger=None, seed=None, verbose=True):

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    if logger is None:
        logger = core.SimpleLogger()

    logger.write('cis-QTL mapping: empirical p-values for phenotypes')
    logger.write('  * {} samples'.format(phenotype_df.shape[1]))
    logger.write('  * {} phenotypes'.format(phenotype_df.shape[0]))
    assert np.all(phenotype_df.columns==covariates_df.index)
    logger.write('  * {} covariates'.format(covariates_df.shape[1]))
    residualizer = core.Residualizer(torch.tensor(covariates_df.values, dtype=torch.float32).to(device))
    dof = phenotype_df.shape[1] - 2 - covariates_df.shape[1]
    logger.write('  * {} variants'.format(genotype_df.shape[0]))

    genotype_ix = np.array([genotype_df.columns.tolist().index(i) for i in phenotype_df.columns])
    genotype_ix_t = torch.from_numpy(genotype_ix).to(device)

    n_samples = phenotype_df.shape[1]
    ix = np.arange(n_samples)
    if seed is not None:
        logger.write('  * using seed {}'.format(seed))
        np.random.seed(seed)
    permutation_ix_t = torch.LongTensor(np.array([np.random.permutation(ix) for i in range(nperm)])).to(device)

    res_df = []
    igc = rpg_igc(genotype_df, variant_df, phenotype_df, phenotype_pos_df,
                                rpg_df, l_window, s_window, NonHiCType)
    if igc.n_phenotypes == 0:
        raise ValueError('No valid phenotypes found.')
    start_time = time.time()
    logger.write('  * computing permutations')

    for k, (phenotype, genotypes, genotype_range, phenotype_id) in enumerate(igc.generate_data(verbose=verbose), 1):
        if phenotype is None:
            continue
            
        genotypes_t = torch.tensor(genotypes, dtype=torch.float).to(device)
        genotypes_t = genotypes_t[:,genotype_ix_t]
        core.impute_mean(genotypes_t)

        # filter monomorphic variants
        mono_t = (genotypes_t == genotypes_t[:, [0]]).all(1)
        if mono_t.any():
            genotypes_t = genotypes_t[~mono_t]
            genotype_range = genotype_range[~mono_t.cpu()]
            logger.write('    * WARNING: excluding {} monomorphic variants'.format(mono_t.sum()))

        if genotypes_t.shape[0] == 0:
            logger.write('WARNING: skipping {} (no valid variants)'.format(phenotype_id))
            continue

        phenotype_t = torch.tensor(phenotype, dtype=torch.float).to(device)

        res = cis.calculate_cis_permutations(genotypes_t, phenotype_t, permutation_ix_t, residualizer=residualizer)
        r_nominal, std_ratio, var_ix, r2_perm, g = [i.cpu().numpy() for i in res]
        var_ix = genotype_range[var_ix]
        variant_id = variant_df.index[var_ix]
        tss_distance = variant_df['pos'].values[var_ix] - igc.phenotype_tss[phenotype_id]
        res_s = cis.prepare_cis_output(r_nominal, r2_perm, std_ratio, g, genotypes_t.shape[0], dof, variant_id, tss_distance, phenotype_id, nperm=nperm)
        if beta_approx:
            res_s[['pval_beta', 'beta_shape1', 'beta_shape2', 'true_df', 'pval_true_df']] = core.calculate_beta_approx_pval(r2_perm, r_nominal*r_nominal, dof)
        res_df.append(res_s)

    res_df = pd.concat(res_df, axis=1, sort=False).T
    res_df.index.name = 'phenotype_id'
    logger.write('  Time elapsed: {:.2f} min'.format((time.time()-start_time)/60))
    logger.write('done.')
    return res_df.astype(output_dtype_dict).infer_objects()


def run_independent(genotype_df, variant_df, cis_df, phenotype_df, phenotype_pos_df, covariates_df,
                        rpg_df, l_window=1000000, s_window=300000, NonHiCType='s_window', eGene_list=None,
                        fdr=0.05, fdr_col='qval', nperm=10000, logger=None, seed=None, verbose=True):

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    assert np.all(phenotype_df.index==phenotype_pos_df.index)
    assert np.all(covariates_df.index==phenotype_df.columns)
    if logger is None:
        logger = core.SimpleLogger()
    
    if eGene_list:
        signif_df = cis_df[cis_df['phenotype_id'].isin(eGene_list)].copy()
    else:
        signif_df = cis_df[cis_df[fdr_col]<=fdr].copy()
        
    cols = [
        'num_var', 'beta_shape1', 'beta_shape2', 'true_df', 'pval_true_df',
        'variant_id', 'tss_distance', 'ma_samples', 'ma_count', 'maf', 'ref_factor',
        'pval_nominal', 'slope', 'slope_se', 'pval_perm', 'pval_beta',
    ]
    signif_df = signif_df[cols]
    signif_threshold = signif_df['pval_beta'].max()
    # subset significant phenotypes
    ix = phenotype_df.index[phenotype_df.index.isin(signif_df.index)]
    
    logger.write('cis-QTL mapping: conditionally independent variants')
    logger.write('  * {} samples'.format(phenotype_df.shape[1]))
    logger.write('  * {}/{} significant phenotypes'.format(signif_df.shape[0], cis_df.shape[0]))
    logger.write('  * {} covariates'.format(covariates_df.shape[1]))
    logger.write('  * {} variants'.format(genotype_df.shape[0]))
    # print('Significance threshold: {}'.format(signif_threshold))
    phenotype_df = phenotype_df.loc[ix]
    phenotype_pos_df = phenotype_pos_df.loc[ix]

    genotype_ix = np.array([genotype_df.columns.tolist().index(i) for i in phenotype_df.columns])
    genotype_ix_t = torch.from_numpy(genotype_ix).to(device)
    dof = phenotype_df.shape[1] - 2 - covariates_df.shape[1]
    ix_dict = {i:k for k,i in enumerate(genotype_df.index)}

    # permutation indices
    n_samples = phenotype_df.shape[1]
    ix = np.arange(n_samples)
    if seed is not None:
        logger.write('  * using seed {}'.format(seed))
        np.random.seed(seed)
    permutation_ix_t = torch.LongTensor(np.array([np.random.permutation(ix) for i in range(nperm)])).to(device)

    res_df = []
    igc = rpg_igc(genotype_df, variant_df, phenotype_df, phenotype_pos_df,
                                rpg_df, l_window, s_window, NonHiCType)
    if igc.n_phenotypes == 0:
        raise ValueError('No valid phenotypes found.')
    logger.write('  * computing independent QTLs')
    start_time = time.time()

    for k, (phenotype, genotypes, genotype_range, phenotype_id) in enumerate(igc.generate_data(verbose=verbose), 1):
        # if no valid data, skip to next batch
        if phenotype is None:
            continue
        
        # copy genotypes to GPU
        phenotype_t = torch.tensor(phenotype, dtype=torch.float).to(device)
        genotypes_t = torch.tensor(genotypes, dtype=torch.float).to(device)
        genotypes_t = genotypes_t[:,genotype_ix_t]
        core.impute_mean(genotypes_t)

        # 1) forward pass
        forward_df = [signif_df.loc[phenotype_id]]  # initialize results with top variant
        covariates = covariates_df.values.copy()  # initialize covariates
        dosage_dict = {}
        while True:
            # add variant to covariates
            variant_id = forward_df[-1]['variant_id']
            ig = genotype_df.values[ix_dict[variant_id], genotype_ix]
            dosage_dict[variant_id] = ig
            covariates = np.hstack([covariates, ig.reshape(-1,1)]).astype(np.float32)
            dof = phenotype_df.shape[1] - 2 - covariates.shape[1]
            residualizer = core.Residualizer(torch.tensor(covariates, dtype=torch.float32).to(device))

            res = cis.calculate_cis_permutations(genotypes_t, phenotype_t, permutation_ix_t, residualizer=residualizer)
            r_nominal, std_ratio, var_ix, r2_perm, g = [i.cpu().numpy() for i in res]
            x = core.calculate_beta_approx_pval(r2_perm, r_nominal*r_nominal, dof)
            # add to list if empirical p-value passes significance threshold
            if x[0] <= signif_threshold:
                var_ix = genotype_range[var_ix]
                variant_id = variant_df.index[var_ix]
                tss_distance = variant_df['pos'].values[var_ix] - igc.phenotype_tss[phenotype_id]
                res_s = cis.prepare_cis_output(r_nominal, r2_perm, std_ratio, g, genotypes.shape[0], dof, variant_id, tss_distance, phenotype_id, nperm=nperm)
                res_s[['pval_beta', 'beta_shape1', 'beta_shape2', 'true_df', 'pval_true_df']] = x
                forward_df.append(res_s)
            else:
                break
        forward_df = pd.concat(forward_df, axis=1, sort=False).T
        dosage_df = pd.DataFrame(dosage_dict)

        # 2) backward pass
        if forward_df.shape[0]>1:
            back_df = []
            variant_set = set()
            for k,i in enumerate(forward_df['variant_id'], 1):
                covariates = np.hstack([
                    covariates_df.values,
                    dosage_df[np.setdiff1d(forward_df['variant_id'], i)].values,
                ])
                dof = phenotype_df.shape[1] - 2 - covariates.shape[1]
                residualizer = core.Residualizer(torch.tensor(covariates, dtype=torch.float32).to(device))

                res = cis.calculate_cis_permutations(genotypes_t, phenotype_t, permutation_ix_t, residualizer=residualizer)
                r_nominal, std_ratio, var_ix, r2_perm, g = [i.cpu().numpy() for i in res]
                var_ix = genotype_range[var_ix]
                variant_id = variant_df.index[var_ix]
                x = core.calculate_beta_approx_pval(r2_perm, r_nominal*r_nominal, dof)
                if x[0] <= signif_threshold and variant_id not in variant_set:
                    tss_distance = variant_df['pos'].values[var_ix] - igc.phenotype_tss[phenotype_id]
                    res_s = cis.prepare_cis_output(r_nominal, r2_perm, std_ratio, g, genotypes.shape[0], dof, variant_id, tss_distance, phenotype_id, nperm=nperm)
                    res_s[['pval_beta', 'beta_shape1', 'beta_shape2', 'true_df', 'pval_true_df']] = x
                    res_s['rank'] = k
                    back_df.append(res_s)
                    variant_set.add(variant_id)
            if len(back_df)>0:
                res_df.append(pd.concat(back_df, axis=1, sort=False).T)
        else:  # single independent variant
            forward_df['rank'] = 1
            res_df.append(forward_df)

    res_df = pd.concat(res_df, axis=0, sort=False)
    res_df.index.name = 'phenotype_id'
    logger.write('  Time elapsed: {:.2f} min'.format((time.time()-start_time)/60))
    logger.write('done.')
    return res_df.reset_index().astype(output_dtype_dict)
