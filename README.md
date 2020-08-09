# code
This repository contains all the code we write in the course of solving the hackathon problems.

The problems we are currently looking at are:

### DDT1-06 - Nucleotide analogue library (Main problem)

Development of nucleotide analogs library by performing virtual screening using molecular docking methodologies at the active site of RNA dependent RNA polymerase enzyme of SAS-COV-2
 
### DDT1-11 - Glycan docking (Main problem)

To design small molecules that can target the glycan shield of SARS-CoV-2 using any methodology. The sugar molecules (glycans) should be correctly modelled as per latest MassSpec data. Molecules binding these sugar moieties can be identified from (a) virtual screening (b) Literature search of molecules binding different sugar moieties (c) machine learning/AI approaches from a data set of molecules that bind sugar moieties. The strength of the binding must be shown by estimating binding free energy on the glycosylated residue of choice. The top 100 or top 25% compounds from your list should be validated by providing a binding site and binding free energy from free energy calculations.
 
### DDT2-04 - Generate molecules and dock (Main problem)

The sequence identity of the COVID19 protease and that of SARS-CoV is high, hence by using the known SARS-CoV protease drugs generate possible drugs using machine learning methods to generate novel drug like candidates. Use a variational autoencoder with SMILES representation to generate novel molecules from the trained continuous latent space. 

+ docking, MD of the best molecules
 
### DDT2-10 - GANs for peptide (Main problem)

The challenge is to either build GAN’S for bioactive peptide generation from scratch using python based deep learning frameworks or customize existing GAN implementations developed for new molecule generation based on SMILES . Success criterion is a python pipeline that utilizes GAN’s to generate potential bioactive peptides < 2000 kDa.
 
### DDT2-12 - Improve MOLS algorithm (Main problem)

To do :  Parallelize the code of imolsdock program to improve its computational time. We can also try to come up with new scoring functions to improve its efficiency(I am not sure how feasible this will be because they have tried out various scoring functions and decided upon the current one).  

The source codes for the programs are in FORTRAN language. The code to be parallelized are : smols subroutine in smmols.f and sminimize subroutine in smminimiz.f. They have given inconsistent instructions on parallelization in two places -  one says just parallelize and other says parallelize for CUDA enabled GPUs. Not sure which one. If anyone is interested, smdrive.f has the main program.

The next part is to dock(flexible receptor docking) a library of fda approved drugs against main protease.
 
### DDT2-14 - Liver Injury (Main problem)

Build a pipeline/model or use scripting language like Python to predict toxic effects based on the chemical structure of known DILI drugs.
 
smiles-pipeline contains code for predicting various outputs from smiles data, which will be useful for DDT2-04 and DDT2-14

Using the nr-ar.smiles dataset, we can get an F1 score of 0.4 on androgen receptivity using morgan fingerprints


