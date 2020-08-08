# -*- coding: utf-8 -*-
"""
Created on Sun Jun 14 17:03:20 2020

@Author: Zhi-Jiang Yang, Dong-Sheng Cao
@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China
@Homepage: http://www.scbdd.com
@Mail: yzjkid9@gmail.com; oriental-cds@163.com
@Blog: https://blog.iamkotori.com

♥I love Princess Zelda forever♥
"""


from decimal import Decimal
import multiprocessing as mp

import os
import pandas as pd
import numpy as np
import scipy as sc
import bz2
import _pickle as cPickle

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import SVG
from rdkit.Chem import rdDepictor

try:
    from .fingerprints import Circular, Path, FunctionGroup
except Exception:
    from fingerprints import Circular, Path, FunctionGroup

try:
    pd.set_option('display.max_colwidth', None)
except:
    pd.set_option('display.max_colwidth', -1)


class NotFittedError(Exception):
    pass


def Pvalue(n, m, ns, ms):
    """
    get P-value

    Parameters
    ----------
    n : int
        n compounds has m compounds with aspecific activity label.
    m : int
        n compounds has m compounds with aspecific activity label.
    ns : int
        a substructure is found in ns compounds.
    ms : int
        the amount of those compounds with the specific activity label is ms.

    Yields
    ------
    res : float
        P-value.

    """
    for i in range(ms, ns+1):
        numerator = Decimal(sc.math.factorial(float(ns)))
        denominatorA = Decimal(sc.math.factorial(
            float(i))) * Decimal(sc.math.factorial(float(ns-i)))
        denominatorB = (m/n)**float(i)
        denominatorC = (1 - (m/n))**(ns-i)
        p_value = float(numerator/denominatorA) * denominatorB * denominatorC
        yield p_value


class BaseLearner(object):

    def __init__(self):
        """
        """
        self.meanPvalue = None
        self.meanMatrix = None

    def GetMatrix(self, mols, **kwgrs):
        pass

    def fit(self, mols,
            labels, aimLabel=1,
            minNum=5, pThreshold=0.05,
            accuracy=None, Bonferroni=False,
            svg=True):
        """
        """
        matrix = self.GetMatrix(mols, svg=svg)

        bo = (matrix.sum(axis=0) >= minNum).values
        matrix = matrix.loc[:, bo]

        n = len(labels)
        m = (labels == aimLabel).sum()

        meanPvalue = {}
        for col, val in matrix.iteritems():

            ns = val.sum()
            ms = (val[labels == aimLabel] == 1).sum()

            pvalue = sum(Pvalue(n, m, ns, ms))
            if pvalue <= pThreshold:
                meanPvalue[col] = pvalue

        meanPvalue = pd.DataFrame(meanPvalue, index=['Pvalue']).T

        if Bonferroni:
            meanPvalue['Pvalue'] = meanPvalue.Pvalue.values * len(meanPvalue)
            meanPvalue = meanPvalue[meanPvalue.Pvalue <= pThreshold]
        else:
            pass

        meanMatrix = matrix.reindex(meanPvalue.index, axis=1)
        meanPvalue['Total'] = meanMatrix.sum(axis=0)
        meanPvalue['Hitted'] = meanMatrix[labels == 1].sum(axis=0)
        meanPvalue['Accuracy'] = meanPvalue.Hitted/meanPvalue.Total
        meanPvalue['Coverage'] = meanPvalue['Hitted']/m

        if accuracy is not None:
            meanPvalue = meanPvalue[meanPvalue.Accuracy >= accuracy]
            meanMatrix = matrix.reindex(meanPvalue.index, axis=1)
        else:
            pass

        self.substructure = self.substructure.reindex(meanPvalue.index)
        meanPvalue = pd.concat(
            (meanPvalue, self.substructure), axis=1, sort=False)
        self.meanPvalue, self.meanMatrix = meanPvalue, meanMatrix
        return self

    def predict(self, mols):
        if self.meanPvalue is None or self.meanMatrix is None:
            raise NotFittedError(
                "This instance is not fitted yet. Call 'fit' with appropriate arguments before using this method.")

        predMatrix = self.GetMatrix(mols, svg=True)
        cols = set(predMatrix.columns) & set(self.meanPvalue.index)
        predMatrix = predMatrix.loc[:, cols].reset_index(drop=True)
        y_pred = (predMatrix.sum(axis=1) > 0) + 0

        return y_pred.values, predMatrix

    def savePvalue(self, file):

        pd.set_option('colheader_justify', 'center')   # FOR TABLE <th>

        html_string = """<html>
        <head><title>HTML Pandas Dataframe with CSS</title></head>
        <style>
        /* includes alternating gray and white with on-hover color */

        .mystyle {
            font-size: 11pt; 
            font-family: Arial;
            border-collapse: collapse; 
            border: 1px solid silver;

        }

        .mystyle td, th {
            padding: 5px;
        }

        .mystyle tr:nth-child(even) {
            background: #E0E0E0;
        }

        .mystyle tr:hover {
            background: silver;
            cursor: pointer;
        }
        </style>
        <body>
            %s
        </body>
    </html>"""
        html = html_string%(self.meanPvalue.to_html(classes='mystyle', escape=False))
        with open(file, 'w') as f:
            f.write(html)
        f.close()
        return html

    def saveModel(self, file):

        with bz2.BZ2File(file + '.pbz2', 'w') as f: 
            cPickle.dump(self, f)
        f.close()
    


class CircularLearner(BaseLearner, Circular):

    def __init__(self,
                 maxRadius=2, minRadius=1,
                 nBits=1024, folded=False,
                 maxFragment=True, nJobs=1,
                 ):
        """
        """
        Circular.__init__(self,
                          maxRadius, minRadius,
                          nBits, folded,
                          maxFragment, nJobs)

    def GetMatrix(self, mols, **kwgrs):

        matrix = self.GetCircularMatrix(mols, **kwgrs)
        return matrix


class PathLeanrner(BaseLearner, Path):

    def __init__(self,
                 minPath=1, maxPath=7,
                 nBits=1024, folded=False,
                 nJobs=1, maxFragment=True):
        """
        """
        Path.__init__(self,
                      minPath, maxPath,
                      nBits, folded, nJobs,
                      maxFragment)

    def GetMatrix(self, mols, **kwgrs):

        matrix = self.GetPathMatrix(mols, **kwgrs)
        return matrix


class FunctionGroupLearner(BaseLearner, FunctionGroup):

    def __init__(self, nJobs=1):

        FunctionGroup.__init__(self, nJobs)

    def GetMatrix(self, mols, **kwgrs):

        matrix = self.GetFunctionGroupsMatrix(mols, **kwgrs)
        return matrix 


if '__main__' == __name__:
    from rdkit import Chem
    from itertools import compress

    data = pd.read_csv(r'tests\Carc\Carc.txt', sep='\t')
    # data = data.sample(n=100)

    mols = data.SMILES.map(lambda x: Chem.MolFromSmiles(x))
    bo = mols.notna().values

    mols = list(compress(mols, bo))
    y_true = data.Label.values[bo]

    # circular = CircularLearner(minRadius=1, maxRadius=6,
    #                            maxFragment=True, nJobs=4)

    # circular.fit(mols, y_true)
    # y_pred, predMatrix = circular.predict(mols)
    # print(y_pred)
    # print(circular.meanPvalue)

    # path = PathLeanrner(minPath=1,
    #                     maxPath=7, nJobs=4,
    #                     maxFragment=True)
    # path.fit(mols, y_true, svg=True)
    # print(path.meanPvalue)

    fg = FunctionGroupLearner(nJobs=4)

    fg.fit(mols, y_true)
    print(fg.predict(mols)[-1])
    # print('Done !!!')
