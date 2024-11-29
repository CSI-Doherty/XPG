# XPG: Xenium Gene Panel Generation 

*XPG* is a novel algorithm that iteratively selects genes of interest for panel generation. We validated XPG using three single-cell RNA sequencing datasets, demonstrating that smaller panels designed by XPG can capture over 80% of the neighbourhood connections that would be identified using the full transcriptome. This approach offers a robust solution for enhancing the efficiency and accuracy of spatial transcriptomic analyses.
![image](https://github.com/user-attachments/assets/c3d2c804-93ae-46ea-a89d-196809199af8)



Installing from GitHub
```
library(devtools)   
devtools::install_github("CSI-Doherty/XPG")
```

Usage
```
output = XPG(g, sp, graph, seu, core, ct, k, cat = 'group')

Input
g      #Genes of interest
sp     #Starting panel gene list
graph  #Nearest Neighbour graph in Seurat object
seu    #Seurat object
k      #number of neighbours
ct     #celltype column in Seurat object
core   #number of available cores
cat    #average per cell type = 'group' or per cell = 'cell'

```

