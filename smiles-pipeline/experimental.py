#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 17:21:32 2020

@author: julian

This file is for visualizing molecules in the SMILES format using rdkit
CCOC(=O)C(O)(C1=CC=C(Cl)C=C1)C2=CC=C(Cl)C=C2
"""
"""Useful commands
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
"""
from rdkit import Chem, DataStructs
from rdkit.Chem import Draw, Descriptors, AllChem

import numpy as np
import pandas as pd


from sklearn.linear_model import LogisticRegression
from sklearn.naive_bayes import GaussianNB
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, f1_score

data = pd.read_csv('data/nr-ar.smiles',
                   sep='\t',
                   names=['SMILES','Identifier','Activity'])
minipos = data[data['Activity']==1].head()
minineg = data[100:105]
minidata = pd.concat([minipos,minineg])

def validSMILES(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    return True

def morganArrayFromSmiles(smiles, nBits):
    """Accepts a SMILES string and a size parameter nBits and
    produces a numpy array of the required size"""
    RADIUS = 2
    mol = Chem.MolFromSmiles(smiles)
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=RADIUS, nBits=nBits)
    arr = np.zeros((1,))
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr
    
def createDataset(df):
    """takes a pandas dataframe with smiles and a target and creates the 
    X feature matrix and y target vector to train an ML models"""
    df['Valid'] = df['SMILES'].apply(validSMILES)
    cleandf = df[df['Valid']==True]
    X = [morganArrayFromSmiles(sm,2048) for sm in cleandf['SMILES']]
    X = np.array(X)
    y = np.array(cleandf['Activity'])
    assert len(X)==len(y)
    return X,y

X,y = createDataset(data)

Xtrain, Xtest, ytrain, ytest = train_test_split(X, y, random_state=1)

model = LogisticRegression(class_weight='balanced')
#model = LogisticRegression()
#model = GaussianNB()
model.fit(Xtrain, ytrain)                  
y_model = model.predict(Xtest)   
ytrain_model = model.predict(Xtrain)
print('Train accuracy: ',accuracy_score(ytrain, ytrain_model))
print('Test accuracy: ',accuracy_score(ytest, y_model))
print('F1 score: ',f1_score(ytest, y_model))

"""Comments
with 50 bits in the fingerprint
A Gaussian Naive Bayes Model gets to 84 percent accuracy 
but this is meaningless. F1 score is 0.2

with 2048 bits 
GNB
Train accuracy:  0.8306968790081232
Test accuracy:  0.791025641025641
F1 score:  0.12834224598930483
"""




























