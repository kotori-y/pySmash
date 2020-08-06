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

import pandas as pd
import numpy as np
import scipy as sc

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import SVG
from rdkit.Chem import rdDepictor

try:
    from .fingerprints import Circular, Path, FunctionGroup
except Exception:
    from fingerprints import Circular, Path, FunctionGroup


class NotFittedError(Exception):
    pass
    # def __init__(self, msg):
    #     self.message = msg


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


def HighlightAtoms(mol, highlightAtoms, figsize=[400, 200], kekulize=True):
    """This function is used for showing which part
    of fragment matched the SMARTS by the id of atoms
    This function is derived from Scopy.

    :param mol: The molecule to be visualized
    :type mol: rdkit.Chem.rdchem.Mol
    :param highlightAtoms: The atoms to be highlighted
    :type highlightAtoms: tuple
    :param figsize: The resolution ratio of figure
    :type figsize: list
    :return: a figure with highlighted molecule
    :rtype: IPython.core.display.SVG

    """
    def _revised(svg_words):
        """
        """
        svg_words = svg_words.replace(
            'stroke-width:2px', 'stroke-width:1.5px').replace(
                'font-size:17px', 'font-size:15px').replace(
                    'stroke-linecap:butt', 'stroke-linecap:square').replace(
                        'fill:#FFFFFF', 'fill:none').replace(
                        'svg:', '')
        return svg_words

    mc = Chem.Mol(mol.ToBinary())

    if kekulize:
        try:
            Chem.Kekulize(mc)
        except:
            mc = Chem.Mol(mol.ToBinary())
    if not mc.GetNumConformers():
        rdDepictor.Compute2DCoords(mc)
    drawer = rdMolDraw2D.MolDraw2DSVG(*figsize)
    drawer.DrawMolecule(mc, highlightAtoms=highlightAtoms)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    # It seems that the svg renderer used doesn't quite hit the spec.
    # Here are some fixes to make it work in the notebook, although I think
    # the underlying issue needs to be resolved at the generation step
    return SVG(_revised(svg))


def ShowResult(subMatrix, subPvalue, label_field='Label',
               smiles_field='SMILES', smarts_field='SMARTS',
               pvalue_field='Val', aim_label=1, topx=50):

    if smiles_field is not None:
        subMatrix = subMatrix.set_index(smiles_field)
    if smarts_field is not None:
        subPvalue = subPvalue.set_index(smarts_field)

    labels = subMatrix[label_field].values
    subPvalue = subPvalue.sort_values(pvalue_field)

    imgs = []
    smas = subPvalue.index[:topx]
    for sma in smas:
        patt = Chem.MolFromSmarts(sma)
        smi = subMatrix[(subMatrix[sma] == 1) & (labels == aim_label)
                        ].sample(n=1).index.values[0]
        mol = Chem.MolFromSmiles(smi)
        atoms = mol.GetSubstructMatches(patt)[0]
        svg = HighlightAtoms(mol, atoms)
        imgs.append(svg.data)

    ns = subMatrix.loc[:, smas].sum(axis=0).values
    ms = subMatrix.loc[:, smas][labels == aim_label].sum(axis=0).values
    pvalues = subPvalue[pvalue_field][:topx].map(
        lambda x: '{:.3e}'.format(x)).values
    acc = list(map(lambda x: '{:.3f}'.format(x), ms/ns))

    out = pd.DataFrame({'SMARTS': smas,
                        'Total Number of Compounds': ns,
                        'Number of Positive Compounds': ms,
                        'p-value': pvalues,
                        'Accuracy': acc,
                        'Substructure': imgs})
    return out
    


class CircularLearner(Circular):

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
        self.meanPvalue = None
        self.meanMatrix = None

    def fit(self, mols,
            labels, aimLabel=1,
            minNum=5, pThreshold=0.05,
            accuracy=None, Bonferroni=False,
            svg=True):
        """
        """
        matrix = self.GetCircularMatrix(mols, svg=svg)

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

        predMatrix = self.GetCircularMatrix(mols)
        cols = set(predMatrix.columns) & set(self.meanPvalue.index)
        predMatrix = predMatrix.loc[:, cols].reset_index(drop=True)
        y_pred = (predMatrix.sum(axis=1) > 0) + 0

        return y_pred.values, predMatrix


class PathLeanrner(Path):

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

    def fit(self, mols,
            labels, aimLabel=1,
            minNum=5, pThreshold=0.05,
            accuracy=0.70, svg=True):
        """
        """
        matrix = self.GetPathMatrix(mols, svg=svg)

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

        meanMatrix = matrix.reindex(meanPvalue.keys(), axis=1)
        meanPvalue = pd.DataFrame(meanPvalue, index=['Pvalue']).T

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

        predMatrix = self.GetCircularMatrix(mols)
        cols = set(predMatrix.columns) & set(self.meanPvalue.index)
        predMatrix = predMatrix.loc[:, cols].reset_index(drop=True)
        y_pred = (predMatrix.sum(axis=1) > 0) + 0

        return y_pred.values, predMatrix


class FunctionGroupLearner(FunctionGroup):

    def __init__(self, nJobs=1):

        FunctionGroup.__init__(self, nJobs)

    def fit(self, mols, 
            labels, aimLabel=1,
            minNum=5, pThreshold=0.05,
            accuracy=0.70, svg=True):

        matrix = self.GetFunctionGroupsMatrix(mols, svg=svg)
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

        meanMatrix = matrix.reindex(meanPvalue.keys(), axis=1)
        meanPvalue = pd.DataFrame(meanPvalue, index=['Pvalue']).T

        meanPvalue['Total'] = meanMatrix.sum(axis=0)
        meanPvalue['Hitted'] = meanMatrix[labels == 1].sum(axis=0)
        meanPvalue['Accuracy'] = meanPvalue.Hitted/meanPvalue.Total
        meanPvalue['Coverage'] = meanPvalue['Hitted']/m

        if accuracy:
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

        predMatrix = self.GetFunctionGroupsMatrix(mols, svg=True)
        
        cols = set(predMatrix.columns) & set(self.meanPvalue.index)
        predMatrix = predMatrix.loc[:, cols].reset_index(drop=True)
        y_pred = (predMatrix.sum(axis=1) > 0) + 0

        return y_pred.values, predMatrix


if '__main__' == __name__:
    from rdkit import Chem
    from itertools import compress

    data = pd.read_csv(r'tests\Canc\Canc.txt', sep='\t')
    # data = data.sample(n=100)

    mols = data.SMILES.map(lambda x: Chem.MolFromSmiles(x))
    bo = mols.notna().values

    mols = list(compress(mols, bo))
    y_true = data.Label.values[bo]

    # circular = CircularLearner(minRadius=1, maxRadius=6,
    #                            maxFragment=True, nJobs=4)

    # circular.fit(mols, y_true)
    # # y_pred, predMatrix = circular.predict(mols)
    # # print(y_pred)
    # print(circular.meanPvalue)

    # path = PathLeanrner(minPath=1,
    #                     maxPath=3, nJobs=4,
    #                     maxFragment=True)
    # path.fit(mols, y_true, svg=True)
    # print(path.)

    fg = FunctionGroupLearner(nJobs=4)

    fg.fit(mols, y_true)
    print(fg.predict(mols)[0])
    print('Done !!!')