print('Import packages and paths')

import warnings
warnings.simplefilter(action='ignore')
import scanpy as sc
import anndata
import numpy as np
import gc
import pandas as pd 
import os
from datetime import date
from matplotlib.pyplot import rc_context
pd.set_option('display.max_columns', None)
import torch
import scarches as sca
import gdown
gc.isenabled()

sc.logging.print_header()

sc.set_figure_params(frameon=False)
sc.set_figure_params(dpi=200)
sc.set_figure_params(figsize=(4, 4))
torch.set_printoptions(precision=3, sci_mode=False, edgeitems=7)

path_project = '/nfs/research/irene/anaelle/CrossSpeciesIC'
path_scripts = os.path.join(path_project, 'expiMap')
path_data = os.path.join(path_project, 'data')
path_data_res = os.path.join(path_scripts, 'expimap_results')
path_models = os.path.join(path_scripts, 'expimap_models')

print('Import data')

homolog_human_mouse = sc.read_h5ad(os.path.join(path_data,'all_homolog_human_mouse.h5ad'))

homolog_human_mouse.layers["raw_counts"] = homolog_human_mouse.X.copy()

#homolog_human_mouse.obs['species'] = np.repeat(['human','mouse'], [76533,159738], axis=0)

mouse_data = homolog_human_mouse[homolog_human_mouse.obs.species == 'mouse']

mouse_data.var = mouse_data.var.set_index('human_gene_name')
mouse_data.var['human_gene_name'] = mouse_data.var.index

print('Start analysis')

print('Import reactome and add it to object')

sca.utils.add_annotations(mouse_data, 'reactome.gmt', min_genes=6, clean=True)

mouse_data._inplace_subset_var(mouse_data.varm['I'].sum(1)>0)

print('Normalisation & hvg')

sc.pp.normalize_total(mouse_data)

sc.pp.log1p(mouse_data)

sc.pp.highly_variable_genes(mouse_data,n_top_genes=2000, batch_key='species',subset=True)

select_terms = mouse_data.varm['I'].sum(0)>6

mouse_data.uns['terms'] = np.array(mouse_data.uns['terms'])[select_terms].tolist()

mouse_data.varm['I'] = mouse_data.varm['I'][:, select_terms]

mouse_data._inplace_subset_var(mouse_data.varm['I'].sum(1)>0)

mouse_data.X = mouse_data.layers["raw_counts"].copy()

print('Genes, GPs : ',mouse_data.varm['I'].shape)

print('ExpiMap Model')

print('First train')

intr_cvae = sca.models.EXPIMAP(
    adata=mouse_data,
    condition_key='species',
    hidden_layer_sizes=[256, 256, 256],
    recon_loss='nb'
)

ALPHA = 0.8

early_stopping_kwargs = {
    "early_stopping_metric": "val_unweighted_loss", # val_unweighted_loss
    "threshold": 0,
    "patience": 50,
    "reduce_lr": True,
    "lr_patience": 13,
    "lr_factor": 0.1,
}

print('start training')

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
    print_stats=True,
)

print('end training')

MEAN = False

mouse_data.obsm['X_cvae'] = intr_cvae.get_latent(mean=MEAN, only_active=True)

print('writing the results')

mouse_data.write_h5ad(os.path.join(path_data_res,'expimap_mouse_'+str(date.today())+'.h5ad'),compression='gzip')

print('Second train')

human_data = homolog_human_mouse[homolog_human_mouse.obs.species == 'human']

human_data = human_data[:,human_data.var.human_gene_name.isin(mouse_data.var_names)]

human_data.var = human_data.var.set_index('human_gene_name')
human_data.var['human_gene_name'] = human_data.var.index

human_data.uns['terms'] = mouse_data.uns['terms']

q_intr_cvae = sca.models.EXPIMAP.load_query_data(human_data, intr_cvae)

print('start training')
q_intr_cvae.train(n_epochs=500, 
                  alpha_epoch_anneal=100, 
                  weight_decay=0., 
                  alpha_kl=0.08, 
                  seed=2020, 
                  print_stats = True,
                  use_early_stopping=True
                 )

print('end training')
q_intr_cvae.save(os.path.join(path_models,'homolog_human_mouse_'+str(date.today()))

print('concatenate the results')
human_mouse = sc.AnnData.concatenate(mouse_data, human_data, batch_key='batch_join', uns_merge='same')

human_mouse.obsm['X_cvae'] = q_intr_cvae.get_latent(human_mouse.X, human_mouse.obs['species'], mean=MEAN, only_active=True)

print('write the results')
human_mouse.write_h5ad(os.path.join(path_data_res,'expimap_human_mouse_'+str(date.today())+'.h5ad'),compression='gzip')

print('end of script')