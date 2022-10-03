# CpmERCCutoff

## Overview
The `CpmERCCutoff` package implements the empirical method to derive log~2~ counts per million (CPM) cutoff to filter out lowly expressed genes using 
ERCC spike-ins as described in Goll and Bosinger et.al (2022)<doi:10.1101/2022.06.23.497396>. This package utilizes the synthetic mRNA control pairs developed by the 
External RNA Controls Consortium (ERCC) (ERCC 1 / ERCC 2) that are spiked into sample pairs at known ratios at various absolute abundances. The relationship 
between the observed and expected fold changes is then used to empirically determine an optimal log2 CPM cutoff for filtering out lowly expressed genes.
  
Filtering out lowly expressed genes prior to differential expressed gene (DEG) analysis is common practice. 
Often, 1 CPM is used as a cut off. The choice of cutoff influences downstream analysis, specifically fold change accuracy, 
differential gene expression assessments, and gene coverage. However, this 1 CPM choice is rather arbitrary. Synthetic ERCC synthetic mRNA spike in pairs 
(ERCC 1 / ERCC 2) that are spiked into sample pairs at known ratios at various absolute abundances provide a way to empirically assess fold change accuracy 
for various CPM levels.  The included 'exp_input' data contains the expected spike in ERCC data, which can be obtained from the 'ERCC Controls Analysis' manual
located on Thermo Fisher's ERCC RNA Spike-In Mix product [page](https://assets.thermofisher.com/TFS-Assets/LSG/manuals/cms_095046.txt).  The `obs_input` data frame 
included in the package is used as input for the observed expression of the 92 ERCC RNA spike-ins and stores the coverage-normalized read CPM that mapped to the 
respective ERCC sequence. Typically, prior to CPM calculation, the read count data is normalized for any systematic differences between samples, for example, 
by using TMM normalization as implemented in the edgeR package.
  
The `CpmERCCutoff` package fits a 3rd order polynomial equation to empirically determine a CPM gene expression cutoff that ensures a high correlation 
with expected fold changes based on  ERCC spike-in pairs.  The polynomial trend shows that the stronger the filtering, the more accurate the fold 
changes but fewer genes are included in the analysis. This package balances these parameters using the lower 95% bound of the bootstrap confidence interval 
as the default CPM threshold. The use of internal controls accounts for technical parameters such as read length, read coverage, and other experimental factors. It
requires ERCC1 and ERCC2 spiked into multiple sample pairs.

## Installation
install the stable version this R package from [CRAN](https://cran.r-project.org/web/packages/CpmERCCutoff/index.html) with:  
`install.packages("CpmERCCutoff")`
  
## Release Notes / What's news
v1.0.0 - Initial Release - 2022-09-13
  
## Citation
If you use `CpmERCCutoff`, please cite:  
The Vacc-SeqQC Project: Benchmarking RNA-Seq for Clinical Vaccine Studies  
Johannes B. Goll, Steven E. Bosinger, Travis L. Jensen, Hasse Walum, Tyler Grimes, Gregory D. Tharp, Muktha S Natrajan, Azra Blazevic, Rich Head, 
Casey E. Gelber, Nirav B. Patel, Patrick Sanz, Nadine G. Rouphael, Evan J. Anderson, Mark J. Mulligan, Daniel F Hoft  
bioRxiv 2022.06.23.497396; doi: https://doi.org/10.1101/2022.06.23.497396  
  
