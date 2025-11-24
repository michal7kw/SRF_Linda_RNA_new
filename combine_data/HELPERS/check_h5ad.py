import scanpy as sc
from pathlib import Path

rna_file = Path("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Linda_RNA/combine_data/results_from_raw/annotation_final.h5ad")
adata_rna = sc.read_h5ad(rna_file)

print(list(adata_rna.obs))