#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 17:21:32 2020

@author: julian

This file is for visualizing molecules in the SMILES format using rdkit
CCOC(=O)C(O)(C1=CC=C(Cl)C=C1)C2=CC=C(Cl)C=C2
"""
from rdkit import Chem, DataStructs
from rdkit.Chem import Draw, Descriptors, AllChem
import numpy as np

"""Useful commands"""  
mol = Chem.MolFromSmiles('CCCO')
smiles = Chem.MolToSmiles(mol)
img = Draw.MolToImage(mol)
smiles_list = ['CCCO',
               'CCOC(=O)C(O)(C1=CC=C(Cl)C=C1)C2=CC=C(Cl)C=C2',
               '[H][C@]12C[C@@]3([H])[C@]4([H])CC[C@@]5([H])C[C@H](CC[C@]5(C)[C@@]4([H])CC[C@]3(C)[C@@]1([H])[C@H](C)[C@]6(CC[C@H](C)CN6)O2)O[C@@H]7O[C@H](CO)[C@H](O[C@@H]8O[C@H](CO)[C@@H](O)[C@H](O[C@@H]9OC[C@@H](O)[C@H](O)[C@H]9O)[C@H]8O[C@@H]%10O[C@H](CO)[C@@H](O)[C@H](O)[C@H]%10O)[C@H](O)[C@H]7O']
mol_list = list(map(Chem.MolFromSmiles, smiles_list))
imgs = Draw.MolsToGridImage(mol_list)
pattern = Chem.MolFromSmiles('C(=O)')
for mol in mol_list:
    print(mol.HasSubstructMatch(pattern))
    
glycine = Chem.MolFromSmiles('C(C(=O)O)N')
glyprint = AllChem.GetMorganFingerprintAsBitVect(glycine, radius=2, nBits=1024)
glyp_arr = np.zeros((1,))
DataStructs.ConvertToNumpyArray(glyprint, glyp_arr)
bi={}
glyanalyse = AllChem.GetMorganFingerprintAsBitVect(glycine, radius=2, nBits=1024, bitInfo=bi)
prints = [(glycine,x,bi) for x in glyanalyse.GetOnBits()]
Draw.DrawMorganBits(prints, molsPerRow=4, legends=[str(x) for x in glyanalyse.GetOnBits()])

#Tanimoto Similarity: Gives number of common on bits divided by total combined on bits
cysteine = Chem.MolFromSmiles('C([C@@H](C(=O)O)N)S')
cyanalyse = AllChem.GetMorganFingerprintAsBitVect(cysteine, radius=2, nBits=1024)

DataStructs.TanimotoSimilarity(glyanalyse,cyanalyse)

#apparently pickling molecules is much faster than 
#reparsing SMILES strings