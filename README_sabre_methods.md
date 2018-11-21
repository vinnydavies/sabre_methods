Code for SABRE methods (including eSABRE). Any problems contact vinny.davies@glasgow.ac.uk

This file contains code for the following methods

eSABRE, nWAIC and biWAIC from \cite{Davies18}

Conjugate SABRE, Semi-Conjugate SABRE and Binary Mask SABRE from \cite{Davies17_COST}

Script Test_Script.R is provided to run all of the code, load simulated data and give an example of the analysis

Code requires the installation of the following R package: mvtnorm, Rlab, coda, msm, pscl, lme4, dummies, doBy, matrixStats

@article{Davies18,
	author = {Vinny Davies and William T. Harvey and Richard Reeve and Dirk Husmeier},
	title = {Improving the identification of antigenic sites in the {H1N1} Influenza virus through accounting for the experimental structure in a sparse hierarchical {B}ayesian model},
	journal = {In Submission},
	year = {2018},
	volume = { },
	number = { },
	pages = { },
}

@article{Davies17_COST,
	author = {Vinny Davies and Richard Reeve and William T. Harvey and Francois F Maree and Dirk Husmeier},
	title = {A sparse hierarchical {B}ayesian model for detecting relevant antigenic sites in virus evolution},
	journal = {Computational Statistics},
	year = {2017},
	volume = {32},
	number = {3},
	pages = {803--843},
}