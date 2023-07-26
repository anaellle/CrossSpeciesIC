# Cross Species Integration and Comparison of single cell RNA sequencing data

This repository is the result of my 3 months internship at the EMBL-EBI. The first aim was to test if the integration tool expiMap was able to perform cross-species integration. Simulations not giving satisfactory results, it suggests that expiMap does not allows cross-species integration.
Moving on, the second purpose was to test SATURN, a promising cross-species integration tool.
To perform these integrations, I used scRNA-seq data from the primary motor cortex beteween human, mouse and fly.

## Data preprocessing


## Expimap

### Anndata input construction

### Training
- reactome
- GO


## SATURN


## Data availability

| Species    | Source           | Name                | Link  |
| ---------- |------------------| --------------------| ------------------------------------------------------------------|
| Human      | Allen brain map  | Human M1 10X        | https://portal.brain-map.org/atlases-and-data/rnaseq/human-m1-10x |
| Mouse      | CellxGene        | 10X nuclei v3 Broad | https://cellxgene.cziscience.com/collections/ae1420fe-6630-46ed-8b3d-cc6056a66467 |
| Drosophila | The Fly Cell Atlas | E-MTAB-10519        | http://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/sc_experiments/E-MTAB-10519/ |

## Packages version
- scanpy==1.9.3
- anndata==0.9.1
- umap==0.5.3
- numpy==1.24.3
- scipy==1.10.1
- pandas==2.0.1
- scikit-learn==1.2.2
- statsmodels==0.14.0
- python-igraph==0.10.4
- louvain==0.8.0
- pynndescent==0.5.10
