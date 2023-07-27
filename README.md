# Cross Species Integration and Comparison of single cell RNA sequencing data

This repository is the result of my 3 months internship at the EMBL-EBI. The first aim was to test if the integration tool **expiMap** was able to perform cross-species integration. Simulations not giving satisfactory results, it suggests that expiMap does not allows cross-species integration.
Moving on, the second purpose was to test **SATURN**, a promising cross-species integration tool.
To perform these integrations, I used scRNA-seq data from the primary motor cortex beteween *Homo sapiens*, *Mus musculus* and *Drosophila*, starting with only mouse and human.

## Data preprocessing

This concern the repository `data`. Each file gives an anndata object, .h5ad file, where the data are cleaned if necessary. For the mouse I used the `mouse_data_from_cellxgene_yuyao.ipynb` file.

## Expimap
### Anndata input construction

To use expiMap, the first step is to create the input. In `expiMap/object_human_mouse` there are all the scripts I made to create the final anndata object with the homologous genes from the human and the mouse. The left part of the matrix is the human counts and the right part is the mouse counts.
1. `ensembl_query.ipynb`: Get the association between the human genes and the mouse genes and create a dataframe with all the homologous genes and their information.
2. `one2one_human_mouse.ipynb`: Create the o2o part of the matrix and the `.var` dataframe associated. For each gene, the raw count of each cell from the human are concatenated with the raw count of each cell from the mouse.
3. `one2many_human_mouse.ipynb`: Create the o2m part of the matrix and the `.var` dataframe associated. For each gene, the raw count of each cell from the human are concatenated with the counts average of the multiple homologous genes associated from the mouse, for each cell.
4. `many2one_human_mouse.ipynb`: Create the m2o part of the matrix and the `.var` dataframe associated. For each gene, the count of each cell from the mouse are concatenated with the counts average of the multiple homologous genes associated from the human, for each cell.
5. `many2many_human_mouse.ipynb`: Create the m2m part of the matrix and the `.var` dataframe associated. For each group of homologous genes, the counts average of each cell from the human are concatenated with the counts average of each cell from the mouse.
6. `all_homologous_genes.ipynb`: Create the final object by concatenating the o2o, o2m, m2o and m2m part as well as their `.var` dataframe.
7. `reduced_matrix.ipynb`: Create a reduced version of the final object by keeping 20% of each cell type above a threshold.

### Training

The training files are in the `expiMap` folder, I followed both the [basic tutorial](https://docs.scarches.org/en/latest/expimap_surgery_pipeline_basic.html#Basic-tutorial-for-query-to-reference-maping-using-expiMap) and the [advanced tutorial](https://docs.scarches.org/en/latest/expimap_surgery_pipeline_advanced.html). Each type of training has a `.ipynb` file, to do most of the tutorial and visualise the results, and a `.py` file which is used to actually performed the training.
- `expimap_training.ipynb` & `expimap.py`: first attempt of expimap with the `reactome.gmt` annotation file provide in the tutorial and with the whole dataset.
- `expimap-GABAergic.ipynb` & `expimap_GABAergic.py`: second attempt with only the GABAergic cells and the `reactome` annotation.
- `expimap_advanced_GABAergic_GO.ipynb` & `expimap_adv_GABAergic_GO_reduced.py`: third attemps with only the GABAergic cells using the reduced dataset plus using the `Gene Ontology` annotation.
- `GO_matrix.ipynb`: Looking through the `GO_BP_human_gene_id_binary_matrix.csv` file.

## SATURN

The repository `SATURN` is a clone of the project [SATURN](https://github.com/snap-stanford/SATURN). All the original files are still present without modification. All information about SATURN is available in the `SATURN/README.md` file. I added scripts to use SATURN to perform a cross-species integration with the human and the mouse dataset.

### Preprocessing

1. `data_preprocessing.ipynb`: inspired by the script `Vignettes/frog_zebrafish_embryogenesis/dataloader.ipynb`, I load each anndata object and apply the modification required to run SATURN including adding the species name, filtering the cells and the genes, putting the gene name as the `.var_names`. This script also create for the human and mouse 3 other objects each being one specific cell type among *GABAergic*, *Glutamatergic* and *Non-Neuronal*. All the objects are save in the directory `Vignettes/frog_zebrafish_embryogenesis/data`.
2. `gene_names.ipynb`: get human information from ensembl, not needed.
3. `protein_embeddings/gpe_human_mouse.ipynb`: based on the `protein_embeddings/Generate Protein Embedding.ipynb` script, download the reference proteom for mouse and human and generate the embeddings files. Before generating the embeddings files, `gpe_human.sh` and `gpe_mouse.sh` need to be runned, it requires to download the [ESM repository](https://github.com/facebookresearch/esm) and `checkout commit` [839c5b82c6cd9e18baa7a88dcbed3bd4b6d48e47](https://github.com/facebookresearch/esm/commit/839c5b82c6cd9e18baa7a88dcbed3bd4b6d48e47).

### Running SATURN

All the futur scripts are based on `Vignettes/frog_zebrafish_embryogenesis/Train SATURN.ipynb` file and are situated in the directory `Vignettes/frog_zebrafish_embryogenesis/`.
1. `training_human_mouse.ipynb`: create the input for training SATURN, use the whole human and mouse datasets and visualise the results. To run the training, use the `saturn_training.sh` file.
2. `training_human_mouse_GABAergic.ipynb`: same as `training_human_mouse.ipynb` but only for *GABAergic* cells, use the `saturn_training_GABAergic.sh` file to run the training.
3. `training_human_mouse_Glutamatergic.ipynb`: same as `training_human_mouse.ipynb`, only for *Glutamatergic* cells, use the `saturn_training_Glutamatergic.sh` file to run the training.
4. `training_human_mouse_Non-Neuronal.ipynb`: same as `training_human_mouse.ipynb`, only for *Non-Neuronal* cells, use the `saturn_training_Non-Neuronal.sh` file to run the training.
5. `results.ipynb`: visualise the results of the 4 previous integration, PCA and UMAP plots.

## Data availability

| Species    | Source             | Name                | Link  |
| ---------- |--------------------| --------------------| ------------------------------------------------------------------|
| Human      | Allen brain map    | Human M1 10X        | https://portal.brain-map.org/atlases-and-data/rnaseq/human-m1-10x |
| Mouse      | CellxGene          | 10X nuclei v3 Broad | https://cellxgene.cziscience.com/collections/ae1420fe-6630-46ed-8b3d-cc6056a66467 |
| Drosophila | The Fly Cell Atlas | E-MTAB-10519        | http://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/sc_experiments/E-MTAB-10519/ |

## Packages 
### expiMap
See the tutorials to download the packages required.

### SATURN
Use the requirements.txt file and follow the indications present in the `SATURN/README.md` file.

### Packages version
- scanpy==1.9.3
- anndata==0.9.1
- umap==0.5.3
- numpy==1.24.3
- scipy==1.10.1
- pandas==2.0.1
- scikit-learn==1.2.2
- python-igraph==0.10.4
- scarches==0.5.8
- biomart==0.9.2
