# PlosComp2017
This file contains the MATLAB functions an the data used in the study.

"On the effects of alternative optima in context-specific metabolic predictions"
Semidán Robaina-Estévez, Zoran Nikoloski
Plos Computational Biology, 2017


-The file "WORKFLOW" includes all the commands followed during this study, from expression data mapping to the generation of all results presented.

-For consistency we have included the fastcore functions, as provided in the publication:

"Vlassis N, Pacheco MP, Sauter T. Fast Reconstruction of Compact Context-Specific Metabolic Network Models. PLoS Comput Biol. 2014;10(1)" 

The human liver core set from the above publication, "CLiver" is also provided in Data&GEMs.

-We also provide the iMAT function that we have used in this study. It has been obtained by adapting the original code in the COBRA toolbox, as to use iMAT with the gurobi solver. We also provide our implementation (no implementation found) of the iMAT flux variability method that was proposed and used in:

"Shlomi T, Cabili MN, Herrgård MJ, Palsson BØ, Ruppin E. Network-based prediction of human tissue-specific metabolism. Nat Biotechnol. 2008;26(9):1003–10"

-In the case of the CORDA study, we took the classification of the reactions in Recon1 (HC, MC, NC and OT) as well as the liverCORDAnew model and the metabolic test from the original publication:

"Schultz A, Qutub AA. Reconstruction of Tissue-Specific Metabolic Networks Using CORDA. PLoS Comput Biol [Internet]. Public Library of Science; 2016"

The programs depend on a working gurobi installation, as well as its MATLAB wrapper function. An academic license to gurobi as well as the software can be freely obtained from "http://www.gurobi.com/"


