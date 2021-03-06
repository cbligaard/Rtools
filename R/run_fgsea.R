## Script for running fgsea analysis and plotting using the replotGSEA functions
## Christina Bligaard Pedersen

# Loading libraries
library(fgsea)
library(qusage)
library(curatedBreastData)

# Sourcing functions - remember to specify path
source('replotFGSEA.R')


# Load gene sets from MSigDB
# Download the gene sets of interest from: https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp

# As an example, read C2: curated gene sets
Hs.c2 <- read.gmt('c2.all.v7.1.symbols.gmt')


# Loading data sets to test workflow
data(curatedBreastDataExprSetList) # Loads 34 breast cancer data sets

# Extracting a gene expression dataset to process
expr.mat <- processExpressionSetList(curatedBreastDataExprSetList[33])[[1]]@assayData$exprs

# Get the clinical data for these patients
data(clinicalData)
meta <- clinicalData$clinicalTable[match(colnames(expr.mat), clinicalData$clinicalTable$patient_ID),]

# Define variable to test for (the variable defining the groups of interest) - here we choose ER status (positive vs. negative)
ERstatus <- meta$ER_preTrt

# Select gene sets of interest - here three found to be upregulated in ER+ breast cancers in other studies
gene_lists <- list('DOANE_BREAST_CANCER_ESR1_UP' = Hs.c2$DOANE_BREAST_CANCER_ESR1_UP,
                   'VANTVEER_BREAST_CANCER_ESR1_UP' = Hs.c2$VANTVEER_BREAST_CANCER_ESR1_UP,
                   'YANG_BREAST_CANCER_ESR1_UP' = Hs.c2$YANG_BREAST_CANCER_ESR1_UP)


# Run label-permuting fgsea - removing really small gene sets
fgsea <- fgseaLabel(gene_lists, mat = expr.mat, labels = ERstatus, nperm = 10000, minSize = 20) # More stability with many permutations
fgsea
# No significance here - but let's look at the plots anyway



#### Plots ####

# Get the ranking the expression value for each gene used in fgsea
corRanks <- var(scale(t(expr.mat)), scale(as.numeric(ERstatus))[, 1])[,1] # These ranks are extremely close to a signal-to-noise ratio, but using a more naive calculation based on means and sum of sd's
stats <- corRanks[order(corRanks, decreasing=TRUE)] 


# Separate barcode plots per gene set
for (g in names(gene_lists)) {
  replotFGSEA(stats = stats, gene.set.name = g, gene.sets = gene_lists, fgsea = fgsea, class.name = 'ER status',
             enrichment.score.range = c(-0.25, max(fgsea$ES)))
}


# Combined barcode plot for all gene sets in a list
replot_multiFGSEA(stats = stats, gene.sets = gene_lists, fgsea = fgsea, class.name = 'ER status')


