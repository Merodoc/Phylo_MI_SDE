# Evolutionary Modelling Toolkit

This is working code for evolutionary modelling of quantitative traits on a phylogeny. 

2 Major Working Code Folders:

**Workflow_Code:**

Contains active testing environment code for Multiple Imputation, Phylogenetic Tree Interaction and Molecular Modelling:

**Sample_Workflow.Rmd:**
R Notebook detailing the multiple imputation workflow for phylogenetic data in practice.
Dependencies: mice, plyr, ape, PVR (R)
PVR allows for the generation of distance based eigenvectors to act as predictors in Multiple Imputation run through mice.
This workflow allows for the generation of a dataset filling missing values for given phylogenetic data, recommended continuation is to perform analysis on the imputed datasets and pool resultant analyses

**Mars_MI.R**
Active development code for MI of phylogenetic data, see Sample_Workflow for example working:
dependencies: mice, plyr, ape, PVR (R)

**Phylo_SDE.jl**
Active development code for tree analyses and evolution, current model parameters allow null model evolution of traits on a time-scaled phylogeny:
End goal of parameter estimation and historical value prediction of quantitative trait values on a phylogenetic tree, through SDE models.
dependencies: Phylo, CSV, DataFrames, Statistics (Julia)

**Hard_Spheres.jl**
Initial development code for Hard-Spheres Molecular Dynamics simulation, developing alongside SDE models as alternative evolutionary models to OU and BM

**Prac Examples**
Collection of working code for package tutorials, test code for models, example codes - not active development here for posterity

