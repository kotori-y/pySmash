# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 15:05:37 2020

@Author: Zhi-Jiang Yang, Dong-Sheng Cao
@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China
@Homepage: http://www.scbdd.com
@Mail: yzjkid9@gmail.com; oriental-cds@163.com
@Blog: https://blog.iamkotori.com

♥I love Princess Zelda forever♥
"""


import os
import pandas as pd
from rdkit import Chem
from smash import CircularLearner


TEST_DIR = os.path.dirname(os.path.abspath(__file__))
file = os.path.join(TEST_DIR, 'Carc.txt')

data = pd.read_csv(file, sep='\t')
mols = data.SMILES.map(lambda x: Chem.MolFromSmiles(x))
labels = data.Label.values


cirLearner = CircularLearner(
        minRadius=1, maxRadius=4, 
        maxFragment=True, nJobs=-1
    )


if "__main__" == __name__:
    sigPvalue, sigMatrix = cirLearner.fit(
            mols, labels, 
            accCutoff=0.7, pCutoff=0.05, Bonferroni=True
          )
    
    print("The number of significant circular fragments: ", 
          str(len(sigPvalue)), sep='')
    
    print("The shape of circular fragment matrix: ",
          str(sigMatrix.shape), sep="")

    y_pred, predMatrix = cirLearner.predict(mols)
    print(f"Predicted Label: {y_pred}")
    print(f"Predicted Matri Shape: {predMatrix.shape}")