# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 21:35:44 2020

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
from smash import PathLearner


TEST_DIR = os.path.dirname(os.path.abspath(__file__))
file = os.path.join(TEST_DIR, 'Carc.txt')

data = pd.read_csv(file, sep='\t')
mols = data.SMILES.map(lambda x: Chem.MolFromSmiles(x))
labels = data.Label.values


pathLearner = PathLearner(
        minRadius=1, maxRadius=4, 
        maxFragment=True, nJobs=-1
    )


if "__main__" == __name__:
    sigPvalue, sigMatrix = pathLearner.fit(mols, labels)
    
    print("The number of significant path fragments: ", 
          str(len(sigPvalue)), sep='')
    
    print("The shape of path fragment matrix: ",
          str(sigMatrix.shape), sep="")