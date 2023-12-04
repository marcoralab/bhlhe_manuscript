## Code for GRN construction

### Overview

Scripts in this repository provide an example of GRN network construction from single-cell/single-nucleus data corresponding to the manuscript entitled "BHLHE40/41 regulate macrophage/microglia responses associated with Alzheimer’s disease and other disorders of lipid-rich tissues" by Podlesny-Drabiniok, Novikova et al.

Currently, the manuscript can be found on bioRxiv: https://www.biorxiv.org/content/10.1101/2023.02.13.528372v1 

### Building a GRN

`get_metacells_for_a_dataset.R` file demonstrates generation of metacells from single-cell or single-nucleus RNA-seq data. The original file supplied is a count matrix, where rows are genes and colnames are cells. The subsequent functions remove lowly abundant  genes and cells with low counts. Cells with high mitochondrial gene expression are also removed (see mt-geneList.upd.csv for mitochondrial gene list). The data are then normalized and metacells are generated (number of neighbors for MetaCells and number of MetaCells are selected by the user).   The output is a N x M matrix, where N are the genes and M are MetaCells. The code was adapted from the Single Cell Analysis Boot Camp at Columbia University (instructors Lukas Vlahos and Pasquale Laise) in 2019 and accompanying functions are provided in the `pisces_2019_functions` folder. 


`calculate_threshold.sh`, `bootstraps.sh` and `consolidate.sh` were adapted from Califano lab ARACNe documentation: https://califano.c2b2.columbia.edu/aracne and are used to build a GRN using the output MetaCells file and list of desired transcription factors. The resulting file provides a mutual information value for every significant TF-gene interaction.


### References

1. Podlesny-Drabiniok A, Novikova G, Liu Y, Dunst J, Temizer R, Giannarelli C, Marro S, Kreslavsky T, Marcora E, Goate AM. BHLHE40/41 regulate macrophage/microglia responses associated with Alzheimer's disease and other disorders of lipid-rich tissues. bioRxiv [Preprint]. 2023 Feb 13:2023.02.13.528372. doi: 10.1101/2023.02.13.528372
2. Vlahos L, Obradovic A, Worley J, Tan X, Howe A, Laise P, et al. Systematic, Protein Activity-based Characterization of Single Cell State. bioRxiv. 2023. p. 2021.05.20.445002. doi:10.1101/2021.05.20.445002
3. Margolin AA, Nemenman I, Basso K, Wiggins C, Stolovitzky G, Dalla Favera R, et al. ARACNE: an algorithm for the reconstruction of gene regulatory networks in a mammalian cellular context. BMC Bioinformatics. 2006;7 Suppl 1: S7.


