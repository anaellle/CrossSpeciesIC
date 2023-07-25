#¡/bin/bash

cmd="python3 ../../train-saturn.py --in_data=data/human_mouse_run.csv --in_label_col=cell_type --ref_label_col=cell_type --num_macrogenes=2000 --hv_genes=8000 --score_adata --work_dir=. "

bsub -gpu "num=2:j_exclusive=no" -q gpu -M 100G -n 2 "${cmd}"