# R Scripts for DEGs

This folder is where I analyzed all of the data from the different datasets [REMBRANDT](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108474) and [Henry Ford Hospital](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse4290). 

I started learning how to use R and followed many different tutorials, including but not limited to:
- [Analyze your own microarray data in R/Bioconductor - BITS wiki](https://wiki.bits.vib.be/index.php/Analyze_your_own_microarray_data_in_R/Bioconductor)
- [Introduction to gene expression microarray analysis in R and Bioconductor](https://gtk-teaching.github.io/Microarrays-R/)
- [Analyzing DNA Microarray Data Using Bioconductor](https://bioconductor.org/help/course-materials/2002/JAX02/jax-A.pdf)
- [Linear Models for Microarray and RNA-Seq Data User's Guide](https://www.bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/usersguide.pdf)

This file also includes some preprocessing where I was exploring how to use PCA and is also an indicator of how much I learned. 

## housekeeping-genes.R

There are multiple ways to normalize microarray data, and one of the most common ways is Robust Multi-Average (RMA) Normalization. However, qPCR analysis usually uses housekeeping genes as a baseline, so I used common housekeeping genes for glioblastoma from [this paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4426273/) to determine if these genes could be used as a baseline. The way I determined this is by checking if the expression of these genes were statistically different between GBM samples and normal samples. Ultimately, the conclusion was that housekeeping genes could not be used and that RMA normalization is better. 

## astrocytoma-analysis.R

This file is where I analyzes all of the genes that were enriched between astrocytoma and normal samples. It allowed me to see which genes were upregulated and downregulated between these samples. 

The same thing is true for gbm-analysis.R and oligodendroglioma-analysis.R.

