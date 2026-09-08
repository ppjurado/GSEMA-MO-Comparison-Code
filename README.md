# Code for paper:  Effect size meta-analysis framework for integrating pathway alterations across multiomics data

In this repository there are 3 R scripts, in order of execution:

-   LUAD_Enrichment_MO_Analysis.R:
    -   Analysis of LUAD dataset with DPM and GSEMA, heatmaps for
        pathways of interest and exploratory plots
-   Simulation_Testing_set_parameters.R:
    -   Simulation study for the caracterization of the relationship
        between both methods significance with intrinsic characteristic
        of pathways
-   Simulation_Testing_distribution_parametters.R:
    -   Simulation study for the caracterization of the relationship
        between both methods significance and the distribution of the
        simulated effect

After the execution of each script, r memory should be cleaned and r
session restarted, as results and data needed between scripts are
written in Ready/. Therefor, make sure to define this root folder as the
working directory

In order to run the script, LUAD expression and proteomics data must be
downloaded from LinkedOmics database:
<https://kb.linkedomics.org/download#LUAD> Precisely, needed files in
/LinkedOmicsKB are:

LUAD_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt
LUAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt
LUAD_RNAseq_gene_RSEM_coding_UQ_1500_log2_Normal.txt
LUAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt

The folder structure is:

``` text
Figures/              
├── Heatmaps/         <-- Heatmaps of relevant pathways
├── STdp/             <-- Plot results of second script 
│   ├── expr_cor/      
│   ├── expr_sd/
│   ├── prot_prop/
│   ├── size/
│   ├── test_sets/
├── STsp/             <-- Plot results of third script
LinkedOmicsKB         <-- Omics Data folder (Download omic datasets here)
Ready/                <-- For data saving and sharing between scripts
Tables/               <-- Results in .xlsx
```

Scripts were originally run with the following sesion info:

<details>

R version 4.5.2 (2025-10-31 ucrt) Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64 (build 26200)

Matrix products: default LAPACK version 3.12.1

attached base packages: [1] parallel stats graphics grDevices utils
datasets methods base

other attached packages: [1] ActivePathways_2.0.6 GSEMA_0.99.4
doRNG_1.8.6.3 rngtools_1.5.2 doParallel_1.0.17 iterators_1.0.14
foreach_1.5.2 msigdbr_25.1.1

loaded via a namespace (and not attached): [1] Rdpack_2.6.6 DBI_1.2.3
pbapply_1.7-4 GSEABase_1.70.1 rlang_1.1.7 magrittr_2.0.4\
[7] matrixStats_1.5.0 compiler_4.5.2 RSQLite_2.4.6 png_0.1-8 vctrs_0.7.1
pkgconfig_2.0.3\
[13] SpatialExperiment_1.18.1 crayon_1.5.3 fastmap_1.2.0 magick_2.9.0
XVector_0.48.0 graph_1.86.0\
[19] UCSC.utils_1.4.0 bit_4.6.0 cachem_1.1.0 beachmat_2.24.0
GenomeInfoDb_1.44.3 jsonlite_2.0.0\
[25] progress_1.2.3 blob_1.3.0 rhdf5filters_1.20.0 DelayedArray_0.34.1
Rhdf5lib_1.30.0 BiocParallel_1.42.1\
[31] irlba_2.3.7 prettyunits_1.2.0 R6_2.6.1 RColorBrewer_1.1-3
limma_3.64.3 GenomicRanges_1.60.0\
[37] numDeriv_2016.8-1.1 Rcpp_1.1.1 assertthat_0.2.1
SummarizedExperiment_1.38.1 GSVA_2.2.0 IRanges_2.42.0\
[43] Matrix_1.7-4 tidyselect_1.2.1 rstudioapi_0.18.0 abind_1.4-8
codetools_0.2-20 metafor_4.8-0\
[49] curl_7.0.0 plyr_1.8.9 lattice_0.22-7 tibble_3.3.1 S7_0.2.1
Biobase_2.68.0\
[55] KEGGREST_1.48.1 Biostrings_2.76.0 pillar_1.11.1
MatrixGenerics_1.20.0 metadat_1.4-0 stats4_4.5.2\
[61] generics_0.1.4 mathjaxr_2.0-0 ggplot2_4.0.2 S4Vectors_0.46.0
hms_1.1.4 sparseMatrixStats_1.20.0\
[67] scales_1.4.0 xtable_1.8-8 glue_1.8.0 pheatmap_1.0.13 tools_4.5.2
data.table_1.18.2.1\
[73] ScaledMatrix_1.16.0 annotate_1.86.1 babelgene_22.9 XML_3.99-0.22
rhdf5_2.52.1 grid_4.5.2\
[79] impute_1.82.0 rbibutils_2.4.1 AnnotationDbi_1.70.0
SingleCellExperiment_1.30.1 nlme_3.1-168 GenomeInfoDbData_1.2.14\
[85] BiocSingular_1.24.0 HDF5Array_1.36.0 cli_3.6.5 rsvd_1.0.5
S4Arrays_1.8.1 dplyr_1.2.0\
[91] gtable_0.3.6 digest_0.6.39 BiocGenerics_0.54.1 SparseArray_1.8.1
farver_2.1.2 rjson_0.2.23\
[97] memoise_2.0.1 lifecycle_1.0.5 h5mread_1.0.1 httr_1.4.8
statmod_1.5.1 bit64_4.6.0-1

</details>
