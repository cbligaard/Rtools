# Rtools
A collection of R scripts for various purposes. Currently only contains functions for GSEA of gene expression data. 

Based on the script for GSEA barcode plotting from https://github.com/PeeperLab/Rtoolbox, I have implemented my own versions of the scripts for handling the output from the [fgsea](https://bioconductor.org/packages/release/bioc/html/fgsea.html) R package.

### replotFGSEA
Like the original version from PeeperLab, replotGSEA re-plots data from a GSEA in R.
The function relies on processed gene expression data which was the input to fgseaLabel (stats). It further takes the name of a gene set (gene.set.name), a list of gene sets (gene.sets), the fgsea object (fgsea), the name of the variable that is used to rank genes (class.name) and a range of enrichment scores for consistent plotting across multiple gene sets (enrichment.score.range). The last argument may be omitted. Use is exemplified in the script https://github.com/cbligaard/Rtools/R/run_fgsea.R.

![replot_FGSEA plot](images/replot_FGSEA.png?raw=true "Title")

### replot_multiFGSEA
This function plots several barcode plots as one by overlaying the enrichment plots and displaying the barcodes for each separately. It takes the same arguments as above, except for the name of the gene set to plot. Use is exemplified in the script https://github.com/cbligaard/Rtools/R/run_fgsea.R.

![replot_multiFGSEA plot](images/replot_multiFGSEA.png?raw=true "Title")
