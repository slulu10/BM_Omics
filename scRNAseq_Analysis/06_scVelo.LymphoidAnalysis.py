import scanpy as sc
import anndata
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix
import numpy as np
import os
import pandas as pd

prefix = 'Seurat.MetsGBM.CombinedAnalysis.Apr2022.Lymphoid'
subset = 'CD8'

# load sparse matrix:
filename = prefix+"."+subset+".counts.mtx"
X = io.mmread(filename)

# create anndata object
adata = anndata.AnnData(
    X=X.transpose().tocsr()
)

# load cell metadata:
filename = prefix+"."+subset+".metadata.csv"
cell_meta = pd.read_csv(filename)

# load gene names:
filename = prefix+"."+subset+".gene_names.csv"
with open(filename, 'r') as f:gene_names=f.read().splitlines()

# set anndata observations and index obs by barcodes, var by gene names
adata.obs = cell_meta
adata.obs.index = adata.obs['barcode']
adata.var.index = gene_names

# load dimensional reduction:
filename = prefix+"."+subset+".DM.csv"
dm = pd.read_csv(filename)
dm.index = adata.obs.index

filename = prefix+"."+subset+".pca.csv"
pca = pd.read_csv(filename)
pca.index = adata.obs.index

# set pca/dm and umap
adata.obsm['X_DM'] = dm.to_numpy()
adata.obsm['X_pca'] = pca.to_numpy()
adata.obsm['X_umap'] = np.vstack((adata.obs['UMAP_1'].to_numpy(), adata.obs['UMAP_2'].to_numpy())).T

# plot a UMAP colored by sampleID to test:
filename = "."+prefix+"."+subset+".umap"
sc.pl.umap(adata, color=['ID'], frameon=False, save=filename)

# save dataset as anndata format
filename = prefix+"."+subset+".scVelo.h5ad"
adata.write(filename)


#####Step 1: Load and integrate data ############################################################
import scvelo as scv
import cellrank as cr
import numpy as np
import anndata as ad
import matplotlib.pyplot as plt
import igraph

scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo', facecolor='white', dpi=100, frameon=False)
cr.settings.verbosity = 2

# reload dataset
filename = prefix+"."+subset+".scVelo.h5ad"
adata = sc.read_h5ad(filename)

# load loom files for spliced/unspliced matrices for each sample:
with open("Seurat.MetsGBM.CombinedAnalysis.Apr2022.Lymphoid.Looms.csv", 'r') as f:Looms=f.read().splitlines()

i=0
print(Looms[i])
ldatasub = scv.read(Looms[i], cache=True)
barcodes = [bc.split(':')[1] for bc in ldatasub.obs.index.tolist()]
sample = [bc.split(':')[0] for bc in ldatasub.obs.index.tolist()][0]
barcodes = [bc[0:len(bc)-1] + '-' + sample for bc in barcodes]
ldatasub.obs.index = barcodes
ldatasub.var_names_make_unique()
ldata = ldatasub

for i in range(1,len(Looms)):
    print(Looms[i])
    ldatasub = scv.read(Looms[i], cache=True)
    barcodes = [bc.split(':')[1] for bc in ldatasub.obs.index.tolist()]
    sample = [bc.split(':')[0] for bc in ldatasub.obs.index.tolist()][0]
    barcodes = [bc[0:len(bc)-1] + '-' + sample for bc in barcodes]
    ldatasub.obs.index = barcodes
    ldatasub.var_names_make_unique()
    ldata = ldata.concatenate([ldatasub])


# merge matrices into the original adata object
adata = scv.utils.merge(adata, ldata)

# plot umap to check
filename = "."+prefix+"."+subset+".celltype.pdf"
sc.pl.umap(adata, color='celltype', frameon=False, legend_loc='on data', title='', save=filename)

##### Part 2: Computing RNA velocity using scVelo #####################################
##Finally we can actually use scVelo to compute RNA velocity. scVelo allows the user to use the steady-state model from the 
##original 2018 publication as well as their updated dynamical model from the 2020 publication.

#First we inspect the in each of our cell clusters. Next we perform some pre-processing, and then we compute RNA velocity using the steady-state model (stochastic option).
filename = prefix+"."+subset+".celltype.pdf"
scv.pl.proportions(adata, groupby='celltype', figsize=(15, 5), save=filename)

# save dataset as anndata format
filename = prefix+"."+subset+".scVelo.h5ad"
adata.write(filename)

# pre-process
#scv.pp.filter_and_normalize(adata)
#scv.pp.moments(adata)
#scv.pp.filter_and_normalize(adata, enforce=True)#, min_shared_counts=20, n_top_genes=2000)
scv.pp.filter_genes(adata, min_shared_counts=20)
scv.pp.normalize_per_cell(adata, enforce=True)
scv.pp.filter_genes_dispersion(adata, n_top_genes=2000)
scv.pp.log1p(adata)
scv.pp.moments(adata)#, n_pcs=30, n_neighbors=30)

# compute velocity
scv.tl.velocity(adata, mode='stochastic') #stochastic model
#or
scv.tl.recover_dynamics(adata)
scv.tl.velocity(adata, mode='dynamical') #dynamical model

scv.tl.velocity_graph(adata)

# Visualize velocity fields
#We can get a broad visualiztion of RNA velocity across all genes and all cells by visualizing a vector field on top of 
#the 2D dimensional reduction.
#filename = prefix+"."+subset+".embedding.Stochastic.DM.pdf"
#scv.pl.velocity_embedding(adata, basis='DM', frameon=False, save=filename)

adata.obs['celltype'].cat.reorder_categories(['CD8-IL7R-CD69-Tm','CD8-ITGA1-ITGAL-Trm','CD8-TCF7-CD226-Tprog.exh','CD8-ISG High',
    'CD8-GZMK-CCL4-CCL5-Teff','CD8-GZMH-GZMA-CD52-Teff','CD8-VCAM1-IFNG-Tex','CD8-CXCL13-LAG3-Tex'], inplace=True)

cols = ["#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#800000","#FDBF6F","#FF7F00"]
#cols = ["#A6CEE3","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FF7F00"]


filename = prefix+"."+subset+".embedding_grid.Dynamical.DM.pdf"
scv.pl.velocity_embedding_grid(adata, basis='DM', color='celltype', save=filename, title='', scale=0.6, 
    palette=cols)

dm = adata.obsm['X_DM']
adata.obsm['X_DM2'] = np.vstack((dm[:,0], dm[:,2])).T

filename = prefix+"."+subset+".embedding_grid.Dynamical.DM2.pdf"
scv.pl.velocity_embedding_grid(adata, basis='DM2', color='celltype', save=filename, title='', scale=0.8, 
    palette=cols)


adata.obsm['X_DM2'] = np.vstack((dm[:,1], dm[:,2])).T

filename = prefix+"."+subset+".embedding_grid.Dynamical.DM3.pdf"
scv.pl.velocity_embedding_grid(adata, basis='DM2', color='celltype', save=filename, title='', scale=0.8, 
    palette=cols)


filename = prefix+"."+subset+".embedding_stream.Dynamical.png"
scv.pl.velocity_embedding_stream(adata, basis='DM', color=['celltype'], save=filename, title='',palette=cols,legend_loc="none")

filename = prefix+"."+subset+".scVelo.h5ad"
adata.write(filename)

