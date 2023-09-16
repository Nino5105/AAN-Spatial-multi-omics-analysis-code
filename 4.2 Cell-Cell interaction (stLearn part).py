# 1.Import modules

import stlearn as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc

# 2.Control_1

data = st.Read10X("./Control_1/outs/")
data.var_names_make_unique()
st.add.image(adata=data,
             imgpath="./Control_1/outs/spatial/tissue_hires_image.png",
             library_id="Control_1", visium=True)
st.pp.filter_genes(data, min_cells=3)
st.pp.normalize_total(data)
spot_mixtures = pd.read_csv('./ST_Control_1_meta.csv', index_col=0)
labels = spot_mixtures.loc[:,"spot_type"].values.astype(str)
print(labels)
print(spot_mixtures)
spot_mixtures['ID'] = spot_mixtures.index.values
data.obs['ID'] = data.obs_names.values
cell = pd.merge(spot_mixtures,data.obs['ID'],on=['ID'])['ID']
data = data[data.obs['ID'].isin(cell)]
print('Spot mixture order correct?: ',
      np.all(spot_mixtures.index.values==data.obs_names.values)) 
data.obs['cell_type'] = labels
data.obs['cell_type'] = data.obs['cell_type'].astype('category')
data.uns['cell_type'] = spot_mixtures 
st.pl.cluster_plot(data, use_label='cell_type',fname="Control_1_celltype.pdf")
lrs = st.tl.cci.load_lrs(['connectomeDB2020_lit'], species='mouse')
print(len(lrs))
Manual_lr = np.array(["Spp1_Cd44","Spp1_Itgav","Spp1_Itgb1","Spp1_Itgb5","Tnfsf12_Tnfrsf12a","Ccl8_Ccr2","Ccl2_Ccr2","Ccl6_Ccr2","Ccl12_Ccr2","Ccl5_Ccr5","Ccl8_Ccr5","H2-D1_Cd8a","H2-K1_Cd8a","H2-D1_Cd8b1","H2-K1_Cd8b1","H2-Ab1_Cd4","H2-Eb1_Cd4","Itga9_Vcam1","Itga4_Vcam1","Itgb1_Vcam1","Pros1_Axl","Gas6_Axl","Lgals9_Cd44"])
st.tl.cci.run(data, Manual_lr,
                  min_spots = 1, 
                  distance=None,
                  n_pairs=10000,
                  n_cpus=48
                  )
lr_info = data.uns['lr_summary']
print('\n', lr_info)
st.tl.cci.adj_pvals(data, correct_axis='spot',
                   pval_adj_cutoff=0.05, adj_method='fdr_bh')
data.write("Control_1_backup.h5ad")

# 3.Control_2

data = st.Read10X("./Control_2/outs/")
data.var_names_make_unique()
st.add.image(adata=data,
             imgpath="./Control_2/outs/spatial/tissue_hires_image.png",
             library_id="Control_2", visium=True)
st.pp.filter_genes(data, min_cells=3)
st.pp.normalize_total(data)
spot_mixtures = pd.read_csv('./ST_Control_2_meta.csv', index_col=0)
labels = spot_mixtures.loc[:,"spot_type"].values.astype(str)
print(labels)
print(spot_mixtures)
spot_mixtures['ID'] = spot_mixtures.index.values
data.obs['ID'] = data.obs_names.values
cell = pd.merge(spot_mixtures,data.obs['ID'],on=['ID'])['ID']
data = data[data.obs['ID'].isin(cell)]
print('Spot mixture order correct?: ',
      np.all(spot_mixtures.index.values==data.obs_names.values))
data.obs['cell_type'] = labels
data.obs['cell_type'] = data.obs['cell_type'].astype('category')
data.uns['cell_type'] = spot_mixtures
st.pl.cluster_plot(data, use_label='cell_type',fname="Control_2_celltype.pdf")
lrs = st.tl.cci.load_lrs(['connectomeDB2020_lit'], species='mouse')
print(len(lrs))
Manual_lr = np.array(["Spp1_Cd44","Spp1_Itgav","Spp1_Itgb1","Spp1_Itgb5","Tnfsf12_Tnfrsf12a","Ccl8_Ccr2","Ccl2_Ccr2","Ccl6_Ccr2","Ccl12_Ccr2","Ccl5_Ccr5","Ccl8_Ccr5","H2-D1_Cd8a","H2-K1_Cd8a","H2-D1_Cd8b1","H2-K1_Cd8b1","H2-Ab1_Cd4","H2-Eb1_Cd4","Itga9_Vcam1","Itga4_Vcam1","Itgb1_Vcam1","Pros1_Axl","Gas6_Axl","Lgals9_Cd44"])
st.tl.cci.run(data, Manual_lr,
                  min_spots = 1,
                  distance=None, 
                  n_pairs=10000,
                  n_cpus=48
                  )
lr_info = data.uns['lr_summary'] 
print('\n', lr_info)
st.tl.cci.adj_pvals(data, correct_axis='spot',
                   pval_adj_cutoff=0.05, adj_method='fdr_bh')
data.write("Control_2_backup.h5ad")

# 4.AAN_2W_1
data = st.Read10X("./AAN_2W_1/outs/")
data.var_names_make_unique()
st.add.image(adata=data,
             imgpath="./AAN_2W_1/outs/spatial/tissue_hires_image.png",
             library_id="AAN_2W_1", visium=True)
st.pp.filter_genes(data, min_cells=3)
st.pp.normalize_total(data)
spot_mixtures = pd.read_csv('./ST_AAN_2W_1_meta.csv', index_col=0)
labels = spot_mixtures.loc[:,"spot_type"].values.astype(str)
print(labels)
print(spot_mixtures)
spot_mixtures['ID'] = spot_mixtures.index.values
data.obs['ID'] = data.obs_names.values
cell = pd.merge(spot_mixtures,data.obs['ID'],on=['ID'])['ID']
data = data[data.obs['ID'].isin(cell)]
print('Spot mixture order correct?: ',
      np.all(spot_mixtures.index.values==data.obs_names.values))
data.obs['cell_type'] = labels
data.obs['cell_type'] = data.obs['cell_type'].astype('category')
data.uns['cell_type'] = spot_mixtures
st.pl.cluster_plot(data, use_label='cell_type',fname="AAN_2W_1_celltype.pdf")
lrs = st.tl.cci.load_lrs(['connectomeDB2020_lit'], species='mouse')
print(len(lrs))
Manual_lr = np.array(["Spp1_Cd44","Spp1_Itgav","Spp1_Itgb1","Spp1_Itgb5","Tnfsf12_Tnfrsf12a","Ccl8_Ccr2","Ccl2_Ccr2","Ccl6_Ccr2","Ccl12_Ccr2","Ccl5_Ccr5","Ccl8_Ccr5","H2-D1_Cd8a","H2-K1_Cd8a","H2-D1_Cd8b1","H2-K1_Cd8b1","H2-Ab1_Cd4","H2-Eb1_Cd4","Itga9_Vcam1","Itga4_Vcam1","Itgb1_Vcam1","Pros1_Axl","Gas6_Axl","Lgals9_Cd44"])
st.tl.cci.run(data, Manual_lr,
                  min_spots = 1,
                  distance=None, 
                  n_pairs=10000,
                  n_cpus=48
                  )
lr_info = data.uns['lr_summary']
print('\n', lr_info)
st.tl.cci.adj_pvals(data, correct_axis='spot',
                   pval_adj_cutoff=0.05, adj_method='fdr_bh')
data.write("AAN_2W_1_backup.h5ad")

# 5.AAN_2W_2
data = st.Read10X("./AAN_2W_2/outs/")
data.var_names_make_unique()
st.add.image(adata=data,
             imgpath="./AAN_2W_2/outs/spatial/tissue_hires_image.png",
             library_id="AAN_2W_2", visium=True)
st.pp.filter_genes(data, min_cells=3)
st.pp.normalize_total(data)
spot_mixtures = pd.read_csv('./ST_AAN_2W_2_meta.csv', index_col=0)
labels = spot_mixtures.loc[:,"spot_type"].values.astype(str)
print(labels)
print(spot_mixtures)
spot_mixtures['ID'] = spot_mixtures.index.values
data.obs['ID'] = data.obs_names.values
cell = pd.merge(spot_mixtures,data.obs['ID'],on=['ID'])['ID']
data = data[data.obs['ID'].isin(cell)]
print('Spot mixture order correct?: ',
      np.all(spot_mixtures.index.values==data.obs_names.values))
data.obs['cell_type'] = labels
data.obs['cell_type'] = data.obs['cell_type'].astype('category')
data.uns['cell_type'] = spot_mixtures
st.pl.cluster_plot(data, use_label='cell_type',fname="AAN_2W_2_celltype.pdf")
lrs = st.tl.cci.load_lrs(['connectomeDB2020_lit'], species='mouse')
print(len(lrs))
Manual_lr = np.array(["Spp1_Cd44","Spp1_Itgav","Spp1_Itgb1","Spp1_Itgb5","Tnfsf12_Tnfrsf12a","Ccl8_Ccr2","Ccl2_Ccr2","Ccl6_Ccr2","Ccl12_Ccr2","Ccl5_Ccr5","Ccl8_Ccr5","H2-D1_Cd8a","H2-K1_Cd8a","H2-D1_Cd8b1","H2-K1_Cd8b1","H2-Ab1_Cd4","H2-Eb1_Cd4","Itga9_Vcam1","Itga4_Vcam1","Itgb1_Vcam1","Pros1_Axl","Gas6_Axl","Lgals9_Cd44"])
st.tl.cci.run(data, Manual_lr,
                  min_spots = 1,
                  distance=None, 
                  n_pairs=10000,
                  n_cpus=48
                  )
lr_info = data.uns['lr_summary']
print('\n', lr_info)
st.tl.cci.adj_pvals(data, correct_axis='spot',
                   pval_adj_cutoff=0.05, adj_method='fdr_bh')
data.write("AAN_2W_2_backup.h5ad")

# 6.AAN_4W_1
data = st.Read10X("./AAN_4W_1/outs/")
data.var_names_make_unique()
st.add.image(adata=data,
             imgpath="./AAN_4W_1/outs/spatial/tissue_hires_image.png",
             library_id="AAN_4W_1", visium=True)

st.pp.filter_genes(data, min_cells=3)
st.pp.normalize_total(data)
spot_mixtures = pd.read_csv('./ST_AAN_4W_1_meta.csv', index_col=0)
labels = spot_mixtures.loc[:,"spot_type"].values.astype(str)
print(labels)
print(spot_mixtures)
spot_mixtures['ID'] = spot_mixtures.index.values
data.obs['ID'] = data.obs_names.values
cell = pd.merge(spot_mixtures,data.obs['ID'],on=['ID'])['ID']
data = data[data.obs['ID'].isin(cell)]
print('Spot mixture order correct?: ',
      np.all(spot_mixtures.index.values==data.obs_names.values))
data.obs['cell_type'] = labels
data.obs['cell_type'] = data.obs['cell_type'].astype('category')
data.uns['cell_type'] = spot_mixtures
st.pl.cluster_plot(data, use_label='cell_type',fname="AAN_4W_1_celltype.pdf")
lrs = st.tl.cci.load_lrs(['connectomeDB2020_lit'], species='mouse')
print(len(lrs))
Manual_lr = np.array(["Spp1_Cd44","Spp1_Itgav","Spp1_Itgb1","Spp1_Itgb5","Tnfsf12_Tnfrsf12a","Ccl8_Ccr2","Ccl2_Ccr2","Ccl6_Ccr2","Ccl12_Ccr2","Ccl5_Ccr5","Ccl8_Ccr5","H2-D1_Cd8a","H2-K1_Cd8a","H2-D1_Cd8b1","H2-K1_Cd8b1","H2-Ab1_Cd4","H2-Eb1_Cd4","Itga9_Vcam1","Itga4_Vcam1","Itgb1_Vcam1","Pros1_Axl","Gas6_Axl","Lgals9_Cd44"])
st.tl.cci.run(data, Manual_lr,
                  min_spots = 1,
                  distance=None, 
                  n_pairs=10000,
                  n_cpus=48
                  )
lr_info = data.uns['lr_summary']
print('\n', lr_info)
st.tl.cci.adj_pvals(data, correct_axis='spot',
                   pval_adj_cutoff=0.05, adj_method='fdr_bh')
data.write("AAN_4W_1_backup.h5ad")

# 7.AAN_4W_2
data = st.Read10X("./AAN_4W_2/outs/")
data.var_names_make_unique()
st.add.image(adata=data,
             imgpath="./AAN_4W_2/outs/spatial/tissue_hires_image.png",
             library_id="AAN_4W_2", visium=True)
st.pp.filter_genes(data, min_cells=3)
st.pp.normalize_total(data)
spot_mixtures = pd.read_csv('./ST_AAN_4W_2_meta.csv', index_col=0)
labels = spot_mixtures.loc[:,"spot_type"].values.astype(str)
print(labels)
print(spot_mixtures)
spot_mixtures['ID'] = spot_mixtures.index.values
data.obs['ID'] = data.obs_names.values
cell = pd.merge(spot_mixtures,data.obs['ID'],on=['ID'])['ID']
data = data[data.obs['ID'].isin(cell)]
print('Spot mixture order correct?: ',
      np.all(spot_mixtures.index.values==data.obs_names.values))
data.obs['cell_type'] = labels
data.obs['cell_type'] = data.obs['cell_type'].astype('category')
data.uns['cell_type'] = spot_mixtures
st.pl.cluster_plot(data, use_label='cell_type',fname="AAN_4W_2_celltype.pdf")
lrs = st.tl.cci.load_lrs(['connectomeDB2020_lit'], species='mouse')
print(len(lrs))
Manual_lr = np.array(["Spp1_Cd44","Spp1_Itgav","Spp1_Itgb1","Spp1_Itgb5","Tnfsf12_Tnfrsf12a","Ccl8_Ccr2","Ccl2_Ccr2","Ccl6_Ccr2","Ccl12_Ccr2","Ccl5_Ccr5","Ccl8_Ccr5","H2-D1_Cd8a","H2-K1_Cd8a","H2-D1_Cd8b1","H2-K1_Cd8b1","H2-Ab1_Cd4","H2-Eb1_Cd4","Itga9_Vcam1","Itga4_Vcam1","Itgb1_Vcam1","Pros1_Axl","Gas6_Axl","Lgals9_Cd44"])
st.tl.cci.run(data, Manual_lr,
                  min_spots = 1,
                  distance=None, 
                  n_pairs=10000,
                  n_cpus=48
                  )
lr_info = data.uns['lr_summary']
print('\n', lr_info)
st.tl.cci.adj_pvals(data, correct_axis='spot',
                   pval_adj_cutoff=0.05, adj_method='fdr_bh')
data.write("AAN_4W_2_backup.h5ad")

# 8.Spatial LR scores

Control_1 = sc.read_h5ad("Control_1_backup.h5ad")
Control_2 = sc.read_h5ad("Control_2_backup.h5ad")
AAN_2W_1 = sc.read_h5ad("AAN_2W_1_backup.h5ad")
AAN_2W_2 = sc.read_h5ad("AAN_2W_2_backup.h5ad")
AAN_4W_1 = sc.read_h5ad("AAN_4W_1_backup.h5ad")
AAN_4W_2 = sc.read_h5ad("AAN_4W_2_backup.h5ad")

stats = ['lr_scores', 'p_vals', 'p_adjs', '-log10(p_adjs)']
best_lr = ["Ccl8_Ccr2","Ccl6_Ccr2","Ccl5_Ccr5","Ccl8_Ccr5","H2-K1_Cd8a"]
a = best_lr[0]
fig, ax = plt.subplots(ncols=3, nrows=2,figsize=(3,2))
axes = ax.flatten()
st.pl.lr_result_plot(Control_1, use_result='lr_sig_scores', use_lr=a, show_color_bar=True, ax=axes[0],dpi=600,size=4)
axes[0].set_title(f'{a} lr_sig_scores\nin Control_1')
st.pl.lr_result_plot(AAN_2W_1, use_result='lr_sig_scores', use_lr=a, show_color_bar=True, ax=axes[1],dpi=600,size=4)
axes[1].set_title(f'{a} lr_sig_scores\nin AAN_2W_1')
st.pl.lr_result_plot(AAN_4W_1, use_result='lr_sig_scores', use_lr=a, show_color_bar=True, ax=axes[2],dpi=600,size=4)
axes[2].set_title(f'{a} lr_sig_scores\nin AAN_4W_1')
st.pl.lr_result_plot(Control_2, use_result='lr_sig_scores', use_lr=a, show_color_bar=True, ax=axes[3],dpi=600,size=4)
axes[3].set_title(f'{a} lr_sig_scores\nin Control_2')
st.pl.lr_result_plot(AAN_2W_2, use_result='lr_sig_scores', use_lr=a, show_color_bar=True, ax=axes[4],dpi=600,size=4)
axes[4].set_title(f'{a} lr_sig_scores\nin AAN_2W_2')
st.pl.lr_result_plot(AAN_4W_2, use_result='lr_sig_scores', use_lr=a, show_color_bar=True, ax=axes[5],dpi=600,size=4)
axes[5].set_title(f'{a} lr_sig_scores\nin AAN_4W_2')

a = best_lr[1]
fig, ax = plt.subplots(ncols=3, nrows=2,figsize=(3,2))
axes = ax.flatten()
st.pl.lr_result_plot(Control_1, use_result='lr_sig_scores', use_lr=a, show_color_bar=True, ax=axes[0],dpi=600,size=4)
axes[0].set_title(f'{a} lr_sig_scores\nin Control_1')
st.pl.lr_result_plot(AAN_2W_1, use_result='lr_sig_scores', use_lr=a, show_color_bar=True, ax=axes[1],dpi=600,size=4)
axes[1].set_title(f'{a} lr_sig_scores\nin AAN_2W_1')
st.pl.lr_result_plot(AAN_4W_1, use_result='lr_sig_scores', use_lr=a, show_color_bar=True, ax=axes[2],dpi=600,size=4)
axes[2].set_title(f'{a} lr_sig_scores\nin AAN_4W_1')
st.pl.lr_result_plot(Control_2, use_result='lr_sig_scores', use_lr=a, show_color_bar=True, ax=axes[3],dpi=600,size=4)
axes[3].set_title(f'{a} lr_sig_scores\nin Control_2')
st.pl.lr_result_plot(AAN_2W_2, use_result='lr_sig_scores', use_lr=a, show_color_bar=True, ax=axes[4],dpi=600,size=4)
axes[4].set_title(f'{a} lr_sig_scores\nin AAN_2W_2')
st.pl.lr_result_plot(AAN_4W_2, use_result='lr_sig_scores', use_lr=a, show_color_bar=True, ax=axes[5],dpi=600,size=4)
axes[5].set_title(f'{a} lr_sig_scores\nin AAN_4W_2')

a = best_lr[2]
fig, ax = plt.subplots(ncols=3, nrows=2,figsize=(3,2))
axes = ax.flatten()
st.pl.lr_result_plot(Control_1, use_result='lr_sig_scores', use_lr=a, show_color_bar=True, ax=axes[0],dpi=600,size=4)
axes[0].set_title(f'{a} lr_sig_scores\nin Control_1')
st.pl.lr_result_plot(AAN_2W_1, use_result='lr_sig_scores', use_lr=a, show_color_bar=True, ax=axes[1],dpi=600,size=4)
axes[1].set_title(f'{a} lr_sig_scores\nin AAN_2W_1')
st.pl.lr_result_plot(AAN_4W_1, use_result='lr_sig_scores', use_lr=a, show_color_bar=True, ax=axes[2],dpi=600,size=4)
axes[2].set_title(f'{a} lr_sig_scores\nin AAN_4W_1')
st.pl.lr_result_plot(Control_2, use_result='lr_sig_scores', use_lr=a, show_color_bar=True, ax=axes[3],dpi=600,size=4)
axes[3].set_title(f'{a} lr_sig_scores\nin Control_2')
st.pl.lr_result_plot(AAN_2W_2, use_result='lr_sig_scores', use_lr=a, show_color_bar=True, ax=axes[4],dpi=600,size=4)
axes[4].set_title(f'{a} lr_sig_scores\nin AAN_2W_2')
st.pl.lr_result_plot(AAN_4W_2, use_result='lr_sig_scores', use_lr=a, show_color_bar=True, ax=axes[5],dpi=600,size=4)
axes[5].set_title(f'{a} lr_sig_scores\nin AAN_4W_2')

a = best_lr[3]
fig, ax = plt.subplots(ncols=3, nrows=2,figsize=(3,2))
axes = ax.flatten()
st.pl.lr_result_plot(Control_1, use_result='lr_sig_scores', use_lr=a, show_color_bar=True, ax=axes[0],dpi=600,size=4)
axes[0].set_title(f'{a} lr_sig_scores\nin Control_1')
st.pl.lr_result_plot(AAN_2W_1, use_result='lr_sig_scores', use_lr=a, show_color_bar=True, ax=axes[1],dpi=600,size=4)
axes[1].set_title(f'{a} lr_sig_scores\nin AAN_2W_1')
st.pl.lr_result_plot(AAN_4W_1, use_result='lr_sig_scores', use_lr=a, show_color_bar=True, ax=axes[2],dpi=600,size=4)
axes[2].set_title(f'{a} lr_sig_scores\nin AAN_4W_1')
st.pl.lr_result_plot(Control_2, use_result='lr_sig_scores', use_lr=a, show_color_bar=True, ax=axes[3],dpi=600,size=4)
axes[3].set_title(f'{a} lr_sig_scores\nin Control_2')
st.pl.lr_result_plot(AAN_2W_2, use_result='lr_sig_scores', use_lr=a, show_color_bar=True, ax=axes[4],dpi=600,size=4)
axes[4].set_title(f'{a} lr_sig_scores\nin AAN_2W_2')
st.pl.lr_result_plot(AAN_4W_2, use_result='lr_sig_scores', use_lr=a, show_color_bar=True, ax=axes[5],dpi=600,size=4)
axes[5].set_title(f'{a} lr_sig_scores\nin AAN_4W_2')

a = best_lr[4]
fig, ax = plt.subplots(ncols=3, nrows=2,figsize=(3,2))
axes = ax.flatten()
st.pl.lr_result_plot(Control_1, use_result='lr_sig_scores', use_lr=a, show_color_bar=True, ax=axes[0],dpi=600,size=4)
axes[0].set_title(f'{a} lr_sig_scores\nin Control_1')
st.pl.lr_result_plot(AAN_2W_1, use_result='lr_sig_scores', use_lr=a, show_color_bar=True, ax=axes[1],dpi=600,size=4)
axes[1].set_title(f'{a} lr_sig_scores\nin AAN_2W_1')
st.pl.lr_result_plot(AAN_4W_1, use_result='lr_sig_scores', use_lr=a, show_color_bar=True, ax=axes[2],dpi=600,size=4)
axes[2].set_title(f'{a} lr_sig_scores\nin AAN_4W_1')
st.pl.lr_result_plot(Control_2, use_result='lr_sig_scores', use_lr=a, show_color_bar=True, ax=axes[3],dpi=600,size=4)
axes[3].set_title(f'{a} lr_sig_scores\nin Control_2')
st.pl.lr_result_plot(AAN_2W_2, use_result='lr_sig_scores', use_lr=a, show_color_bar=True, ax=axes[4],dpi=600,size=4)
axes[4].set_title(f'{a} lr_sig_scores\nin AAN_2W_2')
st.pl.lr_result_plot(AAN_4W_2, use_result='lr_sig_scores', use_lr=a, show_color_bar=True, ax=axes[5],dpi=600,size=4)
axes[5].set_title(f'{a} lr_sig_scores\nin AAN_4W_2')
st.pl.lr_summary(Control_1, n_top=30, figsize=(3,2))




