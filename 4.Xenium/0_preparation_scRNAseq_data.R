# source the functions from FLXenium folder, which loads useful libraries and has 
# a bunch of functions useful for general analysis of spatial data
miceadds::source.all('~/FLXenium/functions/')

# source the ependymoma project specific functions
source('~/ependymoma/xenium/scripts_revisions/resources/epn_functions.R')

# get the global vars for epn
all_vars <- SetUpEpendymomaGlobalVarsGeneral()
colors <- all_vars$colors

# Load seurat object with both malignant and normal cells
seurat_object <- qread(glue('{all_vars$preparation_dir}/data/seurat_obj_ST_normal_malig_annotated.qs'))

# Load seurat object with just mal cells (which has finer annotations of mal celltypes)
seurat_object_mal <- qread(glue('{all_vars$preparation_dir}/data/seurat_obj_malignant_annotated2.qs'))

# now, we will generate the three references (for zfta-rela, zfta-clusters, and st-yap1)
seurat_obj_zr <- CreateReferenceFromBaseSC_EPN(seurat_object, seurat_object_mal, 'ZR')
qsave(seurat_obj_zr, glue('{all_vars$preparation_dir}/data/ZR_Xenium_projection.qs'))

seurat_obj_nc <- CreateReferenceFromBaseSC_EPN(seurat_object, seurat_object_mal, 'NC')
qsave(seurat_obj_nc, glue('{all_vars$preparation_dir}/data/NC_Xenium_projection.qs'))

seurat_obj_yap <- CreateReferenceFromBaseSC_EPN(seurat_object, seurat_object_mal, 'YAP1')
qsave(seurat_obj_yap, glue('{all_vars$preparation_dir}/data/YAP1_Xenium_projection.qs'))

