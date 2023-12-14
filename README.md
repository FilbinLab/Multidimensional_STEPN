# Ependymoma2023

* **Samples used:** `4EP49, 3EP8, BT1743, 4EP51, 7EP35, 7EP41, 11EP22, 3EP54, 3EP67, 7EP1, 11EP8, 7EP9, 16EP8, 4EP53, MUV43R1, MUV43R2, MUV43R4, WEPN9, 4EP44, 9EP47, 9EP35, BT268, I128034, I128034R1, I128034R2, 1230717, 1599417, 2192118, RD050618, RD141019, RD251119, 4717EP17, FR0671, I050420, 4EP46, 9EP45, 11EP21, 3EP29`

* FASTQ files can be found in the following directory: `/labs/mfilbin/homes/biagi/PanCancer/Ependymoma`

* **Code to reproduce manuscript analyses**

Fig. 1: oncoplot summarizing the clinical information of the samples profiled by single-cell analysis


test

## Coculture processing
Experiment to determine the cell state shifts between EP1NS cells cultures alone (monoculture) or with rat neurons and astrocytes (coculture).
In total, 3x4 plates coculture and 1x4 plates monoculture were sequenced. 1x run of coculture was bad QC, and therefore removed.
The remaining goodQC runs are stored under the following folders:

231207 - CoCulture, 4 plates `/labs/mfilbin/Demultiplexing/231207`

231116 - CoCulture, 4 plates `/labs/mfilbin/homes/biagi/Demultiplexing/231116` 

231103 - MonoCulture, 4 plates `/labs/mfilbin/homes/biagi/Demultiplexing/231103`


1. Demultiplexing. Samples demultiplexed using standard code under Demultiplexing.

2. QC. Script '0_qc.R' was used to remove low quality cells

3. Processing. Script '2_Processing_coculture_only.R' used to create a Seurat object, score cells for metaprograms identified in frozen ZFTA-RELA patient cohort.

4. Plots. Script '3_plots_Mono_co_culture_only.R' used to generate plots for final figures