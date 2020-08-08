# -*- encoding: utf-8 -*-
'''
Created on 2020/07/28 10:18:15

@Author: Zhi-Jiang Yang, Dong-Sheng Cao
@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China
@Homepage: http://www.scbdd.com
@Mail: yzjkid9@gmail.com; oriental-cds@163.com
@Blog: https://blog.iamkotori.com

♥I love Princess Zelda forever♥
'''


import pandas as pd
from rdkit import Chem
from itertools import compress
from smash import CircularLearner, PathLeanrner, FunctionGroupLearner
from sklearn.metrics import accuracy_score, recall_score, precision_score
from sklearn.model_selection import train_test_split
from imblearn.metrics import sensitivity_score, specificity_score
from sklearn.ensemble import RandomForestClassifier
from xgboost import XGBClassifier

SEEDS = (71, 14, 76, 68, 42, 67, 46, 73, 56, 92,)



def buildXgbModel(X_train, y_train, **kwagrs):
    
    model = XGBClassifier(**kwagrs)
    model.fit(X_train, y_train)
    
    return model
    
def evaluateModel(traindedModel, X_test, y_test, **kwagrs):
    
    y_pred = kwagrs.get('y_pred') or traindedModel.predict(X_test)
    
    acc = accuracy_score(y_test, y_pred)
    sen = sensitivity_score(y_test, y_pred)
    spe = specificity_score(y_test, y_pred)
    pre = precision_score(y_test, y_pred)
    
    return dict(zip(['Accurary', 'Sensitivity', 'Specificity', 'Precision'],
                    [acc, sen, spe, pre]
                ))



if '__main__' == __name__:
    # Loading data
    data = pd.read_csv('./carc.txt', sep='\t')
    
    mols = data.SMILES.map(lambda x: Chem.MolFromSmiles(x))
    labels = data.Label.values
    
    # circular = CircularLearner(minRadius=1, maxRadius=6, nJobs=4)
    # circular.fit(mols, labels, pThreshold=0.05, accuracy=0.7)
    
    # pd.set_option('display.max_colwidth', -1)
    # circular.meanPvalue.to_html('out.html', escape=False)
    
    # y_pred, pred_matrix = circular.predict(mols)
    
    
    
    # print('{:12s}: {:.3f}'.format('Accuracy', accuracy_score(y_true, y_pred)))
    # print('{:12s}: {:.3f}'.format('Sensitivity', sensitivity_score(y_true, y_pred)))
    # print('{:12s}: {:.3f}'.format('Specificity', specificity_score(y_true, y_pred)))
    # print('{:12s}: {:.3f}'.format('Precision', precision_score(y_true, y_pred)))
    
    