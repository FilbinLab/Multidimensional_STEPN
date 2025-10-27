# making the plot showing expression of the zftarela_fusion genes combined in celltypes
import pandas as pd
from scipy.sparse import csr_matrix
from tqdm import tqdm
import scanpy as sc
import scipy
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns

# read the metadata file from which we can get the information for each sample (like raw data location etc)
metadata = pd.read_excel('~/Multidimensional_STEPN/4.Xenium/resources/SampleIdentifier.xlsx')
# limit to just the zftarela samples
metadata = metadata[metadata["Subtype"].isin(["ZFTA-RELA"])]
# get the list of sample_names
sample_names = list(metadata['SampleName'])

# location where we saved the h5ads from all the samples: 8_neighborhood, because we were doing neighborhood analysis in python.
data_dir = '/n/data1/dfci/pedonc/filbin/lab/users/shk490/ependymoma' # where the data for epn project is stored (the data dir argument of SetUpEpendymomaGlobalVarsGeneral function in epn_functions.R)
h5ad_directory = f'{data_dir}/raw_data/h5ad_files' # inside the data_dir, where we store the h5ad files corresponding to each sample

anndata_list = []
# loop through all the sample_names, and read in their h5ad and append to the list of h5ads
for sample_name in tqdm(sample_names):
	# sample_name = 'STEPN16_Region_1'
	adata = sc.read_h5ad(f'{h5ad_directory}/{sample_name}.h5ad')
	anndata_list.append(adata)

# combine the individual anndatas into one anndata object
combined_adata = sc.concat(anndata_list)

# rename the MES-like celltype in combined_adata to MES/Hypoxia
combined_adata.obs['Metaprogram'] = combined_adata.obs['Metaprogram'].str.replace(pat = 'MES-like', repl = 'MES/Hypoxia')

# create a new gene in our combined_adata object, which has the three genes added. For this, need to add information to .X and .var 
markers = ['ZFTA_RELA_Fusion1', 'ZFTA_RELA_Fusion2', 'ZFTA_RELA_Fusion3']
zfta_fusion_combined_expression_matrix = csr_matrix(combined_adata[:, markers].X.sum(axis = 1))
# we will generate new adata object, using the new_X and new var
new_X = scipy.sparse.hstack([combined_adata.X, zfta_fusion_combined_expression_matrix])
new_var_row = pd.DataFrame(index = ['ZFTA_RELA_Fusion_combined'])
new_var = pd.concat([combined_adata.var, new_var_row])
new_adata = ad.AnnData(new_X)
new_adata.var = new_var
new_adata.obs = combined_adata.obs

# finally make the plot
order_programs = [ "Neuroepithelial-like", 'Embryonic-like', "Radial-glia-like", 'Embryonic-neuronal-like', "Neuronal-like", 
					"Ependymal-like", "MES/Hypoxia", "T-cells",  "Myeloid", "Endothelial", "Oligodendrocytes", 'Neurons']
cmap_reversed = plt.get_cmap('mako').reversed()
dp = sc.pl.dotplot(new_adata, ['ZFTA_RELA_Fusion_combined'], 'Metaprogram',
					figsize = (4,4),
					categories_order = order_programs,
					vmin = 0.245,
					vmax = 0.9,
					standard_scale = 'var', 
					return_fig = True)
dp.style(dot_edge_color='black', dot_edge_lw=0.5, cmap=cmap_reversed)
dp.savefig('/home/shk490/plot.pdf') # change target destination of plot appropriately if required





















