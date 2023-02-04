# Leveraging Differential Gene Expression Analysis and Agent-Based Modeling to Detect Biomarkers of Brain Cancer Subtypes

By Rianna Santra

This GitHub repository holds the code that I developed for this five-month long independent project. Most of the code was developed in the months of October through December, and comments in those programs were written in December and January. The GitHub repository was made in January.

## R Scripts for DEGs

This is the first part of the project, where I used [R/RStudio](https://posit.co) to find DEGs for different subtypes of brain cancer. Many packages from [Bioconductor](https://bioconductor.org) were used in this project. Some code was made with the help of this [tutorial](https://gtk-teaching.github.io/Microarrays-R/).

### Gene Enrichment Analysis

There were more than 400 DEGs found for each subtype, so [GSEA](https://www.gsea-msigdb.org/gsea/index.jsp) was used to find the meaning of these genes. Some code was made from the help of this [tutorial](https://alexslemonade.github.io/refinebio-examples/02-microarray/pathway-analysis_microarray_02_gsea.html).

## Python Notebook Scripts

After DEGs were found, they were put through multiple models, with SVMs as the main model. The main libraries were [NumPy](https://numpy.org), [Pandas](https://pandas.pydata.org), and [Scikit-Learn](https://scikit-learn.org/stable/), also known as sklearn. The package [MatPlotLib](https://matplotlib.org) was used for plotting most of the results. 

## Simulations

Simulations were then created using [NetLogo](https://ccl.northwestern.edu/netlogo/download.shtml). The first [Tumor](https://ccl.northwestern.edu/netlogo/models/Tumor) model created in this website served as a skeleton for the models developed in this repository. The [User Manual](https://ccl.northwestern.edu/netlogo/docs/), served helpful for this project.

