# -*- coding: utf-8 -*-
# RNA velocity

#
#python 3.8.20
#conda create -n velocity python=3.8
#conda install scanpy==1.9.8 scvelo==0.2.4 anndata==0.8.0 numpy==1.22.0 pandas==1.5.3 matplotlib==3.6.0
#conda install loompy==3.0.6
#conda install cellrank(version1.1.0)



## 
## 
import loompy
import scanpy as sc
import numpy as np
import pandas as pd
import scvelo as scv
import h5py
import anndata
from scipy.io import mmread
import os



files=['~/BYT504821/velocyto/BYT504821.loom',
       '~/GZQ512227/velocyto/GZQ512227.loom',
       '~/HZM491546/velocyto/HZM491546.loom',
       '~/LJY0512610/velocyto/LJY0512610.loom',
       '~/MYX505319/velocyto/MYX505319.loom',
       '~/PYT505362/velocyto/PYT505362.loom',
       '~/RYH507008CG/velocyto/RYH507008CG.loom',
       '~/XJN511104/velocyto/XJN511104.loom',
       '~/ZXR505897/velocyto/ZXR505897.loom',
       '~/GMS0533679/velocyto/GMS0533679.loom',
       '~/WXW534406/velocyto/WXW534406.loom',
       '~/WY0530177/velocyto/WY0530177.loom']
output_filename='~/combined.loom'
loompy.combine(files, output_filename, key="Accession")

# 
combine_data = sc.read('~/combined.loom', cache=True)





## 
import loompy
import scanpy as sc
import numpy as np
import pandas as pd
import scvelo as scv
import h5py
import anndata
from scipy.io import mmread
import os
import matplotlib.pyplot as plt
#import omicverse as ov
#import cellrank as cr


## scanpy: anndata
# load sparse matrix
counts = mmread("~/counts.mtx").T.tocsr()   

# 
genes = pd.read_csv("~/genes.csv", header=None).squeeze("columns") 
cells = pd.read_csv("~/barcodes.csv", header=None).squeeze("columns") 

# AnnData
adata = sc.AnnData(X=counts)
print(adata)
adata.var_names = genes
adata.obs_names = cells

# 
metadata = pd.read_csv("~/metadata.csv", index_col=0)
metadata = metadata.loc[adata.obs_names.intersection(metadata.index)]
adata.obs = metadata

# 
umap_coords = pd.read_csv("~/umap_coords.csv", index_col=0)
umap_coords = umap_coords.loc[adata.obs_names]
adata.obsm["X_umap"] = umap_coords.values

harmony_coords = pd.read_csv("~/harmony_coords.csv", index_col=0)
harmony_coords = harmony_coords.loc[adata.obs_names]
adata.obsm["X_harmony"] = harmony_coords.values
# print("AnnData successfully created with Seurat data!")







## spliced/unspliced 
# read metadata
adata_loom = sc.read_loom("~/combined.loom")

# 
print(os.path.exists("E:/bd_loom/combined/home/deer/test/bd/combined.loom"))  
# 
try:
    with loompy.connect("E:/bd_loom/combined/home/deer/test/bd/combined.loom") as ds:
        print("right cell number:", ds.shape[1], "gene number:", ds.shape[0])
except Exception as e:
    print("wrong", str(e))

# 
# 
common_cells = list(set(adata.obs_names) & set(adata_loom.obs_names))   # 
# 
adata.obs_names
'''
Out:
Index(['AAACCCAAGAGATTCA-1_1', 'AAACCCAAGGTTCACT-1_1', 'AAACCCACAAGCGGAT-1_1',
       'AAACCCACACGGCGTT-1_1', 'AAACCCAGTAATTGGA-1_1', 'AAACCCATCACCGGTG-1_1',
       'AAACCCATCACTACGA-1_1', 'AAACCCATCAGTGCGC-1_1', 'AAACCCATCAGTGTGT-1_1',
       'AAACCCATCCTACAAG-1_1',
       ...
       'TTTGGTTCACTGATTG-1_13', 'TTTGGTTGTTCAACGT-1_13',
       'TTTGGTTTCGTGTGGC-1_13', 'TTTGTTGAGACGGTCA-1_13',
       'TTTGTTGAGCAGAAAG-1_13', 'TTTGTTGAGGAGGGTG-1_13',
       'TTTGTTGCATTGAGCT-1_13', 'TTTGTTGGTAAGATAC-1_13',
       'TTTGTTGGTCTGCCTT-1_13', 'TTTGTTGGTGTGCTTA-1_13'],
      dtype='object', length=104525)
'''
adata_loom.obs_names
'''
Out: 
Index(['BYT504821:AAAGGTAGTAGATGTAx', 'BYT504821:AAAGGATCACAGTACTx',
       'BYT504821:AAACGCTTCAGACAAAx', 'BYT504821:AAATGGAGTCCGCAGTx',
       'BYT504821:AACAACCAGACCATAAx', 'BYT504821:AAACGCTCACTGGATTx',
       'BYT504821:AAAGTCCGTTCCCAAAx', 'BYT504821:AAATGGAGTGCCTGACx',
       'BYT504821:AACAAGACATGAGGGTx', 'BYT504821:AACAACCTCGAGATAAx',
       ...
       'WY0530177:TTTCGATTCCACGTGGx', 'WY0530177:TTTGATCTCTTCCACGx',
       'WY0530177:TTTGGTTCATGTCGTAx', 'WY0530177:TTTGGTTGTCTGTCCTx',
       'WY0530177:TTTGGAGGTTAAACAGx', 'WY0530177:TTTGTTGGTACCTGTAx',
       'WY0530177:TTTCGATGTGTGTTTGx', 'WY0530177:TTTGGTTGTCGCATCGx',
       'WY0530177:TTTGTTGCAGCGGTTCx', 'WY0530177:TTTGGAGAGGCATGCAx'],
      dtype='object', name='CellID', length=123623)
'''
# 

# 
# 
adata.obs_names = [x.split("-")[0] for x in adata.obs_names]
adata_loom.obs_names = [x.split(":")[1][:16] for x in adata_loom.obs_names]

# 
'''
common_cells = list(set(adata.obs_names) & set(adata_loom.obs_names))
adata_filtered = adata[adata.obs_names.isin(common_cells)].copy()
adata_loom_filtered = adata_loom[adata_loom.obs_names.isin(common_cells)].copy()
# 
common_genes = list(set(adata_filtered.var_names) & set(adata_loom_filtered.var_names))
adata_filtered = adata_filtered[:, common_genes].copy()
adata_loom_filtered = adata_loom_filtered[:, common_genes].copy()
'''

# 
print("adata cell :", adata.obs_names.duplicated().sum())
print("adata_loom cell:", adata_loom.obs_names.duplicated().sum())
print("adata gene:", adata.var_names.duplicated().sum())
print("adata_loom gene:", adata_loom.var_names.duplicated().sum())

#
adata = adata[:, ~adata.var_names.duplicated()].copy()
adata_loom = adata_loom[:, ~adata_loom.var_names.duplicated()].copy()
#
adata = adata[ ~adata.obs_names.duplicated(),:].copy()
adata_loom = adata_loom[~adata_loom.obs_names.duplicated(),:].copy()

#
common_cells = list(set(adata.obs_names) & set(adata_loom.obs_names))
adata_filtered = adata[adata.obs_names.isin(common_cells)].copy()
adata_loom_filtered = adata_loom[adata_loom.obs_names.isin(common_cells)].copy()
common_genes = list(set(adata_filtered.var_names) & set(adata_loom_filtered.var_names))
adata_filtered = adata_filtered[:, common_genes].copy()
adata_loom_filtered = adata_loom_filtered[:, common_genes].copy()





## 
adata_filtered.layers["spliced"] = adata_loom_filtered.layers["spliced"]
adata_filtered.layers["unspliced"] = adata_loom_filtered.layers["unspliced"]
# 
adata_filtered.layers["spliced"] = adata_filtered.layers["spliced"].astype(np.float32)
adata_filtered.layers["unspliced"] = adata_filtered.layers["unspliced"].astype(np.float32)


# 
scv.pl.proportions(adata_filtered)
adata_filtered.obs['group'] = adata_filtered.obs['group'].astype('category')
scv.pl.proportions(adata_filtered, groupby='group')
sc.pl.highest_expr_genes(adata_filtered, n_top=20, )
# sc.pl.violin(adata_filtered, keys=['percent.mt'], groupby='celltype')
# sc.pl.violin(adata_filtered, keys=['percent.rp'], groupby='celltype')
# sc.pl.violin(adata_filtered, keys=['percent.hb'], groupby='celltype')


## RNA velocity
# 
scv.pp.filter_and_normalize(adata_filtered, min_shared_counts=30, n_top_genes=2000)
# scv.pp.log1p(adata_filtered)
scv.pp.moments(adata_filtered, n_pcs=30, n_neighbors=30)

#
scv.tl.velocity(adata_filtered, mode='stochastic')
scv.tl.velocity_graph(adata_filtered)

scv.utils.get_transition_matrix()


## 
scv.pl.velocity_embedding(adata_filtered,
                          basis="umap", 
                          color="celltype", 
                          arrow_length=3, arrow_size=2, dpi=120)
scv.pl.velocity_embedding_stream(adata_filtered, 
                                 basis="umap", 
                                 color="celltype",
                                 figsize=(8, 6), 
                                 show=False)

# save PDF
plt.savefig("~/velocity_embedding_celltype.pdf", 
            dpi=300, bbox_inches="tight", format="pdf")
plt.show()



# latent time
scv.tl.recover_dynamics(adata_filtered, n_jobs=12)
scv.tl.velocity(adata_filtered, mode="dynamical")
scv.tl.velocity_graph(adata_filtered)
scv.tl.latent_time(adata_filtered)
'''
scv.pl.scatter(
    adata_filtered,
    color='latent_time',
    color_map='viridis',
    size=50,
    basis='umap',
    title='Latent Time'
)
'''

# RNA velocity
colors = ["#27447C","#73ABCF","#C72228","#9EAAD1","#168676","#F3B169","#B88640"]
# 
scv.pl.velocity_embedding_stream(adata_filtered, basis='umap', color='celltype', 
                                 palette=colors, colorbar=True,
                                 legend_loc='right margin',
                                 legend_fontsize=10,
                                 legend_fontweight='normal',
                                 size=40,
                                 figsize=(8, 7),
                                 save="UMAP_stream.pdf")   
#print(os.getcwd())




scv.pl.velocity_embedding(adata_filtered, 
                          arrow_length=3, 
                          color = 'celltype', 
                          legend_loc = 'on data', 
                          arrow_size=1.4, dpi=150)
scv.pl.velocity_embedding_stream(adata_filtered, 
                                 basis="umap", 
                                 color="celltype", 
                                 cmap="coolwarm", 
                                 figsize=(8, 6))
scv.pl.velocity_embedding_stream(adata_filtered, 
                                 basis="umap", 
                                 color="latent_time", 
                                 cmap="coolwarm", 
                                 figsize=(8, 6))
scv.pl.velocity_embedding_grid(adata_filtered, 
                               basis='umap',
                               color="celltype", 
                               arrow_length=3,
                               arrow_size=2,
                               figsize=(8, 6))







'''
#
top_genes = adata_filtered.var['fit_likelihood'].sort_values(ascending=False).index[:300]
scv.pl.heatmap(adata_filtered, var_names=top_genes, sortby='latent_time', col_color='celltype', n_convolve=100) 


#
scv.tl.rank_velocity_genes(adata_filtered, groupby='celltype', min_corr=.3)
df = pd.DataFrame(adata_filtered.uns['rank_velocity_genes']['names'])
df.head()
df.to_csv('~/rank_velocity_genes.csv', index=False)
# 
kwargs = dict(frameon=True, size=10, linewidth=1.5, add_outline='S100A12/MMP9, TANK/NFKBIA')
scv.pl.scatter(adata_filtered, df['S100A12/MMP9'][:5], ylabel='S100A12/MMP9', **kwargs)
scv.pl.scatter(adata_filtered, df['TANK/NFKBIA'][:5], ylabel='TANK/NFKBIA', **kwargs)

scv.pl.scatter(adata_filtered, basis=df['CD4_c3_TCF7'][:5], ncols=5, frameon=False,color="celltype")



# 
scv.tl.velocity_confidence(adata_filtered)
keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(adata_filtered, c=keys, cmap='viridis', perc=[5, 95])
# 
scv.pl.scatter(
    adata_filtered, 
    color=['velocity_confidence', 'celltype'],
    basis='umap',
    legend_loc='right',
    ncols=1
)
adata_high_conf = adata_filtered[adata_filtered.obs['velocity_confidence'] > 0.5].copy()



# 
scv.tl.velocity_pseudotime(adata_filtered)
scv.pl.scatter(adata_filtered, color='velocity_pseudotime', cmap='gnuplot')
'''






# 
adata_filtered.write("~/adata_velocity.h5ad")





