# BCRClonoCalling10x
---------------------

Package makes clonotype calls for cells and sequences using the hamming distance

## Install 
----------

To download in R run the following commands

```R
library('devtools')
install_github('lwhitmore/BCRClonoCalling10x')
library(BCRClonoCalling10x)
```

## How to run BCRClonoCalling10x
--------------------------------

```R 
BCR.CallClono.HD(contig.list,   seq="aa", V.gene=TRUE, CDR3=TRUE, J.gene=FALSE, chain="both",
     hammingthreshold=0.7, cluster.plot=TRUE,graph.plot=FALSE, results_folder=getwd(), verbose=FALSE)
```

* contig.list = R list of VDJ sequences from cell ranger 
* seq = sequence type to use nucleotide (nt) or amino acid (aa)
* V.gene = whether to include V gene in clonotyping 
* CDR3 = whether to include CDR3 gene in clonotyping 
* J.gene = whether to include J gene in clonotyping
* chain = whether to include heavy, light or both chains
* hammingthreshold = what threshold to use for sequence similarity. All sequences above threshold are a part of the same clonotype
* cluster.plot = matrix of hamming distances 
* graph.plot = graph of cells connected by hamming distances 
* results_folder = folder to output figures and result tables 
* verbose = increase terminal output
 