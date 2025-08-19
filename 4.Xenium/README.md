# Xenium analysis


## 0. Software version requirements 

All code was run in R version 4.3.1, Python v3.10.11, and Seurat v.5.0.2. 


## 1. Preparation of single-cell object from ZFTA-RELA patients for label transfer
Script `0_preparation_scRNAseq_data.R` subsets the sc/snRNA-seq object to ZFTA-RELA patients only (includes both malignant and non-malignant cells) and applies sctransform normalization.

## 2. Pre-processing of 10X Xenium output files
Script `1_preprocessing_v2.R` preprocesses the output files from 10X Xenium by subsetting the files to cells with nCount>0, running SCTransform-based normalization, PCA and UMAP dimensionality reduction.

To run all samples in parallel, use:

```
#!/bin/bash

# Loop through each line in the input file
while IFS=, read -r SampleName SampleID; do
  # Remove leading and trailing whitespace (including carriage return characters)
  SampleName=$(echo "$SampleName" | tr -d '\r')

  # Execute your shell script with the sanitized FolderName
  sh /n/scratch/users/s/sad167/EPN/Xenium/scripts_revisions/1c_preprocessing_v2.sh $SampleName
done < /n/scratch/users/s/sad167/EPN/Xenium/scripts_revisions/SampleIdentifier.csv
```



## 3. Assign cell identities
Script `2_assign_programs_automatic.R` takes the sc/snRNA-seq object identities, and projects them onto the Xenium data. It then scores each cell for the identities used tp build the Xenium panel, and assign a final cell identity to the Xenium cells.

To run all samples in parallel, use:

```
#!/bin/bash

# Loop through each line in the input file
while IFS=, read -r SampleName SampleID; do
  # Remove leading and trailing whitespace (including carriage return characters)
  SampleName=$(echo "$SampleName" | tr -d '\r')

  # Execute your shell script with the sanitized FolderName
  sh /n/scratch/users/s/sad167/EPN/Xenium/scripts_revisions/2c_assign_programs_automatic.sh $SampleName
done < /n/scratch/users/s/sad167/EPN/Xenium/scripts_revisions/SampleIdentifier.csv
```

## 4. Niche analysis
Script `3_NicheAnalysis.R` performs spatial niche analysis.


To run all samples in parallel, use:

```
#!/bin/bash

# Loop through each line in the input file
while IFS=, read -r SampleName SampleID; do
  # Remove leading and trailing whitespace (including carriage return characters)
  SampleName=$(echo "$SampleName" | tr -d '\r')

  # Execute your shell script with the sanitized FolderName
  sh /n/scratch/users/s/sad167/EPN/Xenium/scripts_revisions/3c_NicheAnalysis.sh $SampleName
done < /n/scratch/users/s/sad167/EPN/Xenium/scripts_revisions/SampleIdentifier.csv
```

## 5. Spatial coherence score
Script `4_coherence.R` calculates spatial coherence score for each sample

To run all samples in parallel, use:

```
#!/bin/bash

# Loop through each line in the input file
while IFS=, read -r SampleName SampleID; do
  # Remove leading and trailing whitespace (including carriage return characters)
  SampleName=$(echo "$SampleName" | tr -d '\r')

  # Execute your shell script with the sanitized FolderName
  sh /n/scratch/users/s/sad167/EPN/Xenium/scripts_revisions/4c_coherence.sh $SampleName
done < /n/scratch/users/s/sad167/EPN/Xenium/scripts_revisions/SampleIdentifier.csv
```

To plot results, use script in the folder `plots/Plot_coherence_score.R`


## 6. Plot
All scripts to plot final results can be found in the `plot` folder.

## 7. CellCharter analysis 
All scripts to run the CellCharter analysis and plot results can be found in the `cellcharter` folder. 

CellCharter analysis is performed in three steps using python.
(1) Preprocessing by `scanpy`: `ST_EPN_scanpy.ipynb`. 
(2) Dimension reduction by `scVI`: `ST_EPN_scvi.ipynb`. 
(3) cellcharter analysis: `ST_EPN_cellcharter_k26.ipynb`.

Figures are generated using R.
(1) Figure-4C: Proportions of metaprograms of spatial clusters (`SpatialCluster_Metaprogram_Composition.Rmd`)
(2) Figure-4D: Heatmap of spatial clusters propotions of samples (`Sample_Spatial_Cluster_Heatmap.Rmd`)
(3) Figure-E6D: Proportions of spatial clusters of samples (`Sample_Spatial_Cluster_Heatmap.Rmd`)
(4) Figure-E6F: Pairwise correlation in metaprogram compositions between CellCharter and Meta-niche analysis (`SpatialCluster_MetaNiche_Correlation.Rmd`)
(5) Figure-E6G: Comparing entropies of sample composition between Tumor-enriched and TME-enriched spatial clusters (`Tumor_TME_Entropy_Comparison.Rmd`)



