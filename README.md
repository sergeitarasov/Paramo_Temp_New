# PARAMO pipeline: outputs and notes
1. Tree  plot with NHPP rates
	1. `NHPP_KDE.r` file to construct tree plots and edge profiles for Entire phenotype. Checked by ST.
	2. `NHPP_KDE_BR1.r` file to construct tree plots and edge profiles for body regions at level 1 (15 BRs). Checked by ST.
	3. `NHPP_KDE_BR2.r` file to construct tree plots and edge profiles for body regions at level 2 (main BRs). Not Checked by ST.
2. Edge profile plot with NHPP rates
	1. See `NHPP_KDE` family of files. Note, there is an error in edge_profiles4plotting() function.
3. Phenotype diffusion plot and animation (early version uses based on multidimensional scaling scaling)
	1. All works, see `Morphospace_animation.r`
4. Organismal plot with body regions and rates
	1. All works, see `Plot_Images.r`

# Outline: PARAMO pipeline
1. Creat Stochastic Maps for your characters
	1. get maps using RevBayes or ML approach
	2. discretize them into bins
	3. save in proper zip format
2. Amalgamte maps for a body region (BR) of interest
3. Apply PARAMO's NHPP inference
	1. get state changes over branches
	2. use Kernel density estimation (KDE) with Markov kernel
	3. apply Loess smoothing
4. Plot your results 
	1. tree with rates & edge profiles
	2. Run multidimensional scaling for phenotype diffusion
	3. create organismal vector image for organismal plot

# To DO:
1. Create 95% credible interval inference for NHPP
2. Maybe change Loess to Splines in NHPP plots
3. Fix edge_profiles4plotting() for creating edge profiles
4. Maybe change multidimensional scaling to UMAP in the Phenotype diffusion plot
