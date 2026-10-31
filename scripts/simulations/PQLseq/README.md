# PQLseq: Efficient Mixed Model Analysis of Count Data in Large-Scale Genomic Sequencing Studies

PQLseq is an efficient tool designed for differential analysis of large-scale RNA sequencing (RNAseq) data and Bisulfite sequencing (BSseq) data in the presence of individual relatedness and population structure. PQLseq first fits a Generalized Linear Mixed Model (GLMM) with adjusted covariates, predictor of interest and random effects to account for population structure and individual relatedness, and then performs Wald tests for each gene in RNAseq or site in BSseq.

[![](https://www.r-pkg.org/badges/version/PQLseq?color=green)](https://cran.r-project.org/package=PQLseq)
[![](http://cranlogs.r-pkg.org/badges/grand-total/PQLseq?color=green)](https://cran.r-project.org/package=PQLseq)
[![](http://cranlogs.r-pkg.org/badges/last-month/PQLseq?color=green)](https://cran.r-project.org/package=PQLseq)
[![](http://cranlogs.r-pkg.org/badges/last-week/PQLseq?color=green)](https://cran.r-project.org/package=PQLseq)


:arrow_double_down: Installation
-------------------

Get the released version from CRAN:

``` r
install.packages("PQLseq")
```

How to cite `PQLseq`
-------------------
Shiquan Sun*, Jiaqiang Zhu*, Sahar Mozaffari, Carole Ober, Mengjie Chen, and Xiang Zhou#. *Heritability estimation and differential analysis of count data with generalized linear mixed models in genomic sequencing studies.* Bioinformatics 35, no. 3 (2019): 487-496.
