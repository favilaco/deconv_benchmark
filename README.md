Source code (R statistical programming language, v3.6) to reproduce the results described in the article:

> *Avila Cobos F, Alquicira-Hernandez J, Powell J, Mestdagh P and De Preter K.* **Comprehensive benchmarking of computational deconvolution of transcriptomics data.** *(bioRxiv; https://doi.org/10.1101/2020.01.10.897116)*

DATASETS
========
Here we provide an **example folder** (named "example"; see *"Folder requirements & running the deconvolution"*) that can be directly used. It contains an artificial single-cell RNA-seq dataset made of 5 artificial cell types; 200 cells per cell type and 80 genes.

The **other five external datasets** (together with the necessary metadata) can be downloaded from their respective sources:

* Baron: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84133 (Specifically, GSM2230757 to GSM2230760 for human pancreatic islands)
* GSE81547: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81547
* E-MTAB-5061: https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-5061/
* PBMCs: https://support.10xgenomics.com/single-cell-geneexpression/datasets/1.1.0/fresh_68k_pbmc_donor_a 
* kidney.HCL: https://figshare.com/articles/HCL_DGE_Data/7235471

Regarding E-MTAB-5061: cells with "not_applicable", "unclassified” and “co-expression_cell" labels were excluded and only cells coming from six healthy patients (non-diabetic) were kept.

The following line is needed for fresh installations of Linux (Debian):
`sudo apt-get install curl libcurl4-openssl-dev libssl-dev zlib1g-dev r-base-dev libxml2-dev`


R 3.6.0: REQUIRED PACKAGES AND PACKAGE DEPENDENCIES:
===================================================
Code to be run before running any deconvolution (to be run in **R >= 3.6.0**):
```
packages <- c("devtools", "BiocManager","data.table","ggplot2","tidyverse",
			  "Matrix","matrixStats",
			  "gtools",
			  "foreach","doMC","doSNOW", #for parallelism
			  "Seurat","sctransform", #sc-specific normalization
			  "nnls","FARDEEP","MASS","glmnet","ComICS","dtangle") #bulk deconvolution methods

for (i in packages){ install.packages(i, character.only = TRUE)}

# Installation using BiocManager:
# Some packages that didn't work with install.packages (e.g. may not be present in a CRAN repository chosen by the user)
packages3 = c('limma','edgeR','DESeq2','pcaMethods','BiocParallel','preprocessCore','SingleR','scater','SingleCellExperiment','Linnorm','DeconRNASeq','multtest','GSEABase','annotate','genefilter','preprocessCore','graph','MAST','Biobase') #last two are required by DWLS and MuSiC, respectively.
for (i in packages3){ BiocManager::install(i, character.only = TRUE)}

# Dependencies for CellMix: 'NMF', 'csSAM', 'GSEABase', 'annotate', 'genefilter', 'preprocessCore', 'limSolve', 'corpcor', 'graph', 'BiocInstaller'
packages2 = c('NMF','csSAM','limSolve','corpcor')
for (i in packages2){ install.packages(i, character.only = TRUE)}

# Special instructions for CellMix and DSA
install.packages("BiocInstaller", repos="http://bioconductor.org/packages/3.7/bioc/")
system('wget http://web.cbio.uct.ac.za/~renaud/CRAN/src/contrib/CellMix_1.6.2.tar.gz')
system("R CMD INSTALL CellMix_1.6.2.tar.gz")
system('wget https://github.com/zhandong/DSA/raw/master/Package/version_1.0/DSA_1.0.tar.gz')
system("R CMD INSTALL DSA_1.0.tar.gz")

# Following packages come from Github
devtools::install_github("GfellerLab/EPIC", build_vignettes=TRUE) #requires knitr
devtools::install_github("xuranw/MuSiC") 
devtools::install_bitbucket("yuanlab/dwls", ref="default")
devtools::install_github("meichendong/SCDC")
devtools::install_github("rosedu1/deconvSeq")
devtools::install_github("cozygene/bisque")

```

Users interested in the **generation of pseudo-bulk mixtures from scRNA-seq data** can use the *"Generator"* function that is located inside **helper_functions.R**


References to other methods included in our benchmark:
======================================================
While our work has a **BSD (3-clause)** license, you **may need** to obtain a license to use the individual normalization/deconvolution methods (e.g. CIBERSORT. The source code for CIBERSORT needs to be asked to the authors at https://cibersort.stanford.edu).

| method | ref |
|--------|----------|
| OLS | Chambers, J., Hastie, T. & Pregibon, D. Statistical Models in S. in Compstat (eds. Momirović, K. & Mildner, V.) 317–321 (Physica-Verlag HD, 1990). doi:10.1007/978-3-642-50096-1_48 |
| nnls | Mullen, K. M. & van Stokkum, I. H. M. nnls: The Lawson-Hanson algorithm for non-negative least squares (NNLS). R package version 1.4. https://CRAN.R-project.org/package=nnls |
| FARDEEP | Hao, Y., Yan, M., Lei, Y. L. & Xie, Y. Fast and Robust Deconvolution of Tumor Infiltrating Lymphocyte from Expression Profiles using Least Trimmed Squares. bioRxiv 358366 (2018) doi:10.1101/358366 |
| MASS: Robust linear regression (RLR) | Ripley, B. et al. MASS: Support Functions and Datasets for Venables and Ripley’s MASS. (2002) |
| DeconRNASeq | Gong, T. & Szustakowski, J. D. DeconRNASeq: a statistical framework for deconvolution of heterogeneous tissue samples based on mRNA-Seq data. Bioinforma. Oxf. Engl. 29, 1083–1085 (2013) |
| CellMix: DSA, ssKL, ssFrobenius  | Gaujoux, R. & Seoighe, C. CellMix: a comprehensive toolbox for gene expression deconvolution. Bioinformatics 29, 2211–2212 (2013) |
| DCQ | Altboum, Z. et al. Digital cell quantification identifies global immune cell dynamics during influenza infection. Mol. Syst. Biol. 10, 720 (2014) |
| glmnet: lasso, ridge, elastic net | Friedman, J., Hastie, T. & Tibshirani, R. Regularization Paths for Generalized Linear Models via Coordinate Descent. J. Stat. Softw. 33, 1–22 (2010) |
| EPIC | Racle, J., Jonge, K. de, Baumgaertner, P., Speiser, D. E. & Gfeller, D. Simultaneous enumeration of cancer and immune cell types from bulk tumor gene expression data. eLife 6, e26476 (2017) |
| dtangle | Hunt, G. J., Freytag, S., Bahlo, M. & Gagnon-Bartsch, J. A. dtangle: accurate and robust cell type deconvolution. Bioinformatics 35, 2093–2099 (2019) |
| CIBERSORT | Newman, A. M. et al. Robust enumeration of cell subsets from tissue expression profiles. Nat. Methods 12, 453–457 (2015) |
|--------|----------|
| BisqueRNA | Jew, B. et al. Accurate estimation of cell composition in bulk expression through robust integration of single-cell information. bioRxiv 669911 (2019) doi:10.1101/669911 |
| deconvSeq | Du, R., Carey, V. & Weiss, S. T. deconvSeq: deconvolution of cell mixture distribution in sequencing data. Bioinformatics doi:10.1093/bioinformatics/btz444 |
| DWLS | Tsoucas, D. et al. Accurate estimation of cell-type composition from gene expression data. Nat. Commun. 10, 1–9 (2019) |
| MuSiC | Wang, X., Park, J., Susztak, K., Zhang, N. R. & Li, M. Bulk tissue cell type deconvolution with multi-subject single-cell expression reference. Nat. Commun. 10, 380 (2019) |
| SCDC | Dong, M. et al. SCDC: Bulk Gene Expression Deconvolution by Multiple Single-Cell RNA Sequencing References. Briefings in Bioinformatics (2020), bbz166, https://doi.org/10.1093/bib/bbz166 |
|--------|----------|
| SCTransform / regularized negative binomial regression (RNBR) | Hafemeister, C. & Satija, R. Normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression. Genome Biology (2019) doi:10.1186/s13059-019-1874-1 |
| Linnorm | Yip, S. H., Wang, P., Kocher, J.-P. A., Sham, P. C. & Wang, J. Linnorm: improved statistical analysis for single cell RNA-seq expression data. Nucleic Acids Res. 45, e179–e179 (2017) |
| scran | L. Lun, A. T., Bach, K. & Marioni, J. C. Pooling across cells to normalize single-cell RNA sequencing data with many zero counts. Genome Biol. 17, 75 (2016) |
| scater | McCarthy, D. J., Campbell, K. R., Lun, A. T. L. & Wills, Q. F. Scater: pre-processing, quality control, normalization and visualization of single-cell RNA-seq data in R. Bioinformatics 33, 1179–1186 (2017) |
| Quantile normalization (QN) | Bolstad, B. M., Irizarry, R. A., Åstrand, M. & Speed, T. P. A comparison of normalization methods for high density oligonucleotide array data based on variance and bias. Bioinformatics 19, 185–193 (2003) |
| Upper quartile (UQ) | Bullard, J. H., Purdom, E., Hansen, K. D. & Dudoit, S. Evaluation of statistical methods for normalization and differential expression in mRNA-Seq experiments. BMC Bioinformatics 11, 94 (2010) |
| Trimmed mean of M-values (TMM) | Robinson, M. D. & Oshlack, A. A scaling normalization method for differential expression analysis of RNA-seq data. Genome Biol. 11, R25 (2010) |
| Transcripts per million (TPM) | Li, B., Ruotti, V., Stewart, R. M., Thomson, J. A. & Dewey, C. N. RNA-Seq gene expression estimation with read mapping uncertainty. Bioinformatics 26, 493–500 (2010) |
| LogNormalize | LogNormalize function (part of "Seurat"). R Documentation. https://www.rdocumentation.org/packages/Seurat/versions/3.1.1/topics/LogNormalize ; Butler, A., Hoffman, P., Smibert, P. et al. Integrating single-cell transcriptomic data across different conditions, technologies, and species. Nat Biotechnol 36, 411–420 (2018) doi:10.1038/nbt.4096 |
| Variance stabilization transformation (VST) & Median of ratios | Anders, S. & Huber, W. Differential expression analysis for sequence count data. Genome Biol. 11, R106 (2010) |


FOLDER REQUIREMENTS & RUNNING THE DECONVOLUTION
===============================================

a) Folder structure:
```
.
├── example
│   ├── example.rds
│   └── example_phenoData.txt
├── baron
│   ├── sc_baron.rds
│   └── baron_phenoData.txt
├── GSE81547
│   ├── sc_GSE81547.rds
│   └── GSE81547_phenoData.txt
...

├── helper_functions.R
├── Master_deconvolution.R
└── CIBERSORT.R
```

b) Minimally the following (tab-separated) columns being part of the metadata: "cellID", "cellType", "sampleID". Optionally, other columns may be present (e.g. "gender","disease").

```
# For the baron dataset, it should look like:

		     cellID  cellType sampleID
human1_lib3.final_cell_0178     delta   human1
human1_lib2.final_cell_0498     delta   human1
...
```

c) Each single-cell RNA-seq input ("sc_input") dataset is a integer matrix containing gene names as rows and cellID as columns.


d) Make the following choices:
```
	i) a specific dataset (from "example","baron","GSE81547","E-MTAB-5061","PBMCs")
	ii) data transformation (from "none","log","sqrt","vst"); with "none" meaning linear scale
	iii) type of deconvolution method (from "bulk","sc")
		iii.1) For "bulk" methods:
			iii.1.1) choose normalization method among: "column","row","mean","column_z-score","global_z-score","column_min-max","global_min-max","LogNormalize","QN","TMM","UQ", "median_ratios", "TPM"
			iii.1.2) Marker selection strategy from "all", "pos_fc", "top_50p_logFC", "bottom_50p_logFC", "top_50p_AveExpr", "bottom_50p_AveExpr", "top_n2", "random5" (see main manuscript for more details).
			iii.1.3) choose deconvolution method among: "CIBERSORT","DeconRNASeq","OLS","nnls","FARDEEP","RLR","DCQ","elastic_net","lasso","ridge","EPIC","DSA","ssKL","ssFrobenius","dtangle".

		iii.2) For "sc" methods:
			iii.2.1) choose normalization method for both the reference matrix (scC) and the pseudo-bulk matrix (scT) among: "column","row","mean","column_z-score","global_z-score","column_min-max","global_min-max","LogNormalize","QN","TMM","UQ", "median_ratios", "TPM", "SCTransform","scran","scater","Linnorm" (last 4 are single-cell-specific)
			iii.2.2.) choose deconvolution method among: "MuSiC","BisqueRNA","DWLS","deconvSeq","SCDC"

	iv) Number of cells to be used to make the pseudo-bulk mixtures (multiple of 100)
	v) Cell type to be removed from the reference matrix ("none" for the full matrix; this is dataset dependent: e.g. "alpha" from baron dataset)
	vi) Number of available cores (by default 1, can be enlarged if more resources available)
```

R example calls
===============

For bulk:
---------

```
# With the example we provided with this repository + no cell type removed:
Rscript Master_deconvolution.R example none bulk TMM all nnls 100 none 1
	#Expected output:
	#        RMSE   Pearson
	#1     0.0351    0.9866


# With the example we provided with this repository + "cell_type_1" removed:
Rscript Master_deconvolution.R example none bulk TMM all nnls 100 cell_type_1 1
	#Expected output:
	#       RMSE   Pearson
	#1    0.1038    0.9379


# With baron (or GSE81547, E-MTAB-5061, PBMCs) + no cell type removed:
Rscript Master_deconvolution.R baron none bulk TMM all nnls 100 none 1
	#Expected output:
	#       RMSE   Pearson
	#1    0.0724    0.8961


# With baron + delta cells removed:
Rscript Master_deconvolution.R baron none bulk TMM all nnls 100 delta 1
	#Expected output:
	#        RMSE   Pearson
	#1     0.0887    0.8197
```


For single-cell:
----------------

```
# With the example we provided with this repository + no cell type removed::
Rscript Master_deconvolution.R example none sc TMM TMM MuSiC 100 none 1
	#Expected output:
	#        RMSE   Pearson
	#1     0.0351    0.9866


# With the example we provided with this repository + "cell_type_1" removed:
Rscript Master_deconvolution.R example none sc TMM TMM MuSiC 100 cell_type_1 1
	#Expected output:
	#       RMSE   Pearson
	#1    0.1044    0.9376


# With baron (or GSE81547, E-MTAB-5061, PBMCs) + no cell type removed:
Rscript Master_deconvolution.R baron none sc TMM TMM MuSiC 100 none 1
	#Expected output:
	#        RMSE   Pearson
	#1     0.0488     0.953


# With baron + delta cells removed:
Rscript Master_deconvolution.R baron none sc TMM TMM MuSiC 100 delta 1
	#Expected output:
	#        RMSE   Pearson
	#1      0.073    0.8799
```

sessionInfo() files Linux & macOS
----------------------------------
Please see "sessionInfo_Linux.txt" and "sessionInfo_macOS.txt" in this repository.
