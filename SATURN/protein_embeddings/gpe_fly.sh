#¡/bin/bash

export TORCH_HOME=/nfs/research/irene/anaelle/Scripts/SATURN/torch_home
cd /nfs/research/irene/anaelle/Scripts/esm/scripts/ 
CUDA_VISIBLE_DEVICES=2 
cmd=" python extract.py \
esm1b_t33_650M_UR50S \
/nfs/research/irene/anaelle/Scripts/SATURN/protein_embeddings/data/Drosophila_melanogaster.BDGP6.32.pep.all.clean.fa \
/nfs/research/irene/anaelle/Scripts/SATURN/protein_embeddings/data/Drosophila_melanogaster.BDGP6.32.pep.all.clean.fa_esm1b \
--include mean  --truncate"

bsub -gpu "num=2:j_exclusive=no" -q gpu -M 50G -n 2 "${cmd}"