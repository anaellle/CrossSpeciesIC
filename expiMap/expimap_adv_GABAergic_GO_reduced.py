print('Import packages and data')

import warnings
warnings.simplefilter(action='ignore')
import scanpy as sc
import anndata
print('imported anndata')
import numpy as np
print('imported numpy')
import gc
import pandas as pd 
import os
from datetime import date
from matplotlib.pyplot import rc_context
pd.set_option('display.max_columns', None)
import torch
import scarches as sca
import gdown
print(gc.isenabled())

sc.set_figure_params(frameon=False)
sc.set_figure_params(dpi=200)
sc.set_figure_params(figsize=(4, 4))
torch.set_printoptions(precision=3, sci_mode=False, edgeitems=7)

path_project = '/nfs/research/irene/anaelle/CrossSpeciesIC'
print(path_project)
path_scripts = os.path.join(path_project, 'expiMap')
print(path_scripts)
path_data = os.path.join(path_project, 'data')
print(path_data)
path_data_res = os.path.join(path_scripts, 'expimap_results')
print(path_data_res)
path_models = os.path.join(path_scripts, 'expimap_models')
print(path_models)

adata = sc.read_h5ad(os.path.join(path_data,'all_hhm_reduced.h5ad'))

GO = pd.read_csv(os.path.join(path_data,'GO_BP_human_gene_id_binary_matrix.csv'))

print('Prepare the data')
print('Prep anndata object')
adata.layers["raw_counts"] = adata.X.copy()

mdata = adata[adata.obs.species == 'mouse']

mdata = mdata[mdata.obs.homolog_class_label == 'GABAergic']

print('Prep GPs matrix')

GO = GO[GO.ensembl_gene_id.isin(mdata.var.index)]

mdata = mdata[:,mdata.var.human_ensembl_id.isin(GO.ensembl_gene_id)]

GO = GO.set_index('ensembl_gene_id')
GO['human_ensembl_id'] = GO.index

GOs = GO.reindex(mdata.var.index)

print('index equal :',sum(GOs.index == mdata.var.index))

del GOs[GOs.columns[-1]]

GPs = GOs.to_numpy()

print('GPs matrix:',np.shape(GPs))

#len(GPs.sum(0)) ## this get the sum for each column i.e. GP

#len(GPs.sum(1)) ## this get the sum for each row i.e. gene

#sum(GPs.sum(0)>25)

terms = GOs.columns.values

mdata.varm['I'] = GPs
mdata.uns['terms'] = terms

#sum(mdata.varm['I'].sum(1)>0)

print('Start expimap tutorial')

print('Select GPs and genes')

sc.pp.normalize_total(mdata)
sc.pp.log1p(mdata)
sc.pp.highly_variable_genes(mdata,
                            n_top_genes = 3000,
                            batch_key='species',
                            subset =True)

## We select only the GPs with more than 10 genes, we have 354 GPs

#sum(mdata.varm['I'].sum(0)>10)

select_terms = mdata.varm['I'].sum(0)>8

#select_terms

mdata.uns['terms'] = np.array(mdata.uns['terms'])[select_terms].tolist()

#len(mdata.uns['terms'])

mdata.varm['I'] = mdata.varm['I'][:,select_terms]

print('I shape : ', np.shape(mdata.varm['I']))

## We keep only the genes present in the GPs

mdata._inplace_subset_var(mdata.varm['I'].sum(1)>0)

mdata.X = mdata.layers['raw_counts'].copy()

print('First train on mouse')

intr_cvae = sca.models.EXPIMAP(
    adata=mdata,
    condition_key='species',
    hidden_layer_sizes=[300, 300, 300],
    recon_loss='nb'
)

ALPHA = 0.7

early_stopping_kwargs = {
    "early_stopping_metric": "val_unweighted_loss", # val_unweighted_loss
    "threshold": 0,
    "patience": 50,
    "reduce_lr": True,
    "lr_patience": 13,
    "lr_factor": 0.1,
}
intr_cvae.train(
    n_epochs=400,
    alpha_epoch_anneal=100,
    alpha=ALPHA,
    alpha_kl=0.6,
    weight_decay=0.,
    early_stopping_kwargs=early_stopping_kwargs,
    use_early_stopping=True,
    monitor_only_val=False,
    seed=2020,
)

print('end of first training')
      
MEAN = False

mdata.obsm['X_cvae'] = intr_cvae.get_latent(mean=MEAN, only_active=True)

sc.pp.neighbors(mdata, use_rep='X_cvae')
sc.tl.umap(mdata)

#mdata.write_h5ad(os.path.join(path_data_res,'mdata_reduced_GABAergic__'+str(date.today())+'.h5ad'),compression='gzip')

print('Second train on the human')

hdata = adata[adata.obs.species == 'human']

hdata = hdata[hdata.obs.homolog_class_label == 'GABAergic']

hdata = hdata[:,hdata.var.human_ensembl_id.isin(mdata.var_names)]

hdata.uns['terms'] = mdata.uns['terms']

q_intr_cvae = sca.models.EXPIMAP.load_query_data(hdata, intr_cvae,
                                                 unfreeze_ext=True,
                                                 n_ext=450)

q_intr_cvae.train(n_epochs=250,
                  alpha_epoch_anneal=120,
                  alpha_kl=0.1,
                  weight_decay=0.,
                  alpha_l1=0.96,
                  gamma_ext=0.7,
                  gamma_epoch_anneal=50,
                  beta=3.,
                  seed=2020,
                  use_early_stopping=False)

print('end of second train')
      
q_intr_cvae.save(os.path.join(path_models,'hhm_reduced_GO_GABAergic_adv_'+str(date.today())))

hmdata = sc.AnnData.concatenate(mdata, hdata, batch_key='batch_join', uns_merge='same')

hmdata.obsm['X_cvae'] = q_intr_cvae.get_latent(hmdata.X, hmdata.obs['species'], mean=MEAN, only_active=True)

sc.pp.neighbors(hmdata, use_rep='X_cvae')
sc.tl.umap(hmdata)

hmdata.write_h5ad(os.path.join(path_data_res,'hmdata_reduced_GABAergic_adv_'+str(date.today())+'.h5ad'),compression='gzip')