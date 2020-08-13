Contains a bunch of smiles datasets

# nr-ar.smiles

Each line of the file contains a smiles string, Compound name and a 1 or 0 depending on its androgen receptor activity. Line separator is a space.

I deleted the following five string from the file which rdkit says are not valid strings

[NH4+].[NH4+].F[Si--](F)(F)(F)(F)F	NCGC00166021-01	0

[Cl-][Pt]1([Cl-])NCCN1	NCGC00186461-01	0

[NH4+].[NH4+].[Cl-][Pt++]([Cl-])([Cl-])[Cl-]	NCGC00260349-01	0

[Na+].[Na+].F[Si--](F)(F)(F)(F)F	NCGC00255590-01	0

O.O.O.O.O=C1O[Mg]2(OC(=O)C3=CC=CC=C3O2)OC4=CC=CC=C14	NCGC00181305-01	0

Taken from https://tripod.nih.gov/tox21/challenge/data.jsp which has many more interesting datasets I haven't checked out yet


# trainDILI.csv and testDILI.csv

From 
Predicting drug-induced liver injury: The importance of data curation
Author links open overlay panelEleniKotsampasakouFlorianeMontanariGerhard F.Ecker
https://doi.org/10.1016/j.tox.2017.06.003

Merge test and merge train are from the same paper and contain many interesting feature vectors.