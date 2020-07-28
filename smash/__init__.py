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
    # pd.set_option('display.max_colwidth', -1)
    # out.to_html('out_html.html', escape=False)


class MeaningfulCircular(Circular):

    def __init__(self, mols,
                 maxRadius=2, minRadius=1,
                 nBits=1024, folded=False,
                 maxFragment=True, nJobs=1):
        """
        """
        Circular.__init__(self, mols,
                          maxRadius, minRadius,
                          nBits, folded,
                          maxFragment, nJobs)

    def Pvalue(self, n, m, val, aimLabel):
        """
        """
        ns = val.sum()
        ms = (val[labels] == aimLabel).sum()
        pvalue = sum(Pvalue(n, m, ns, ms))
        return pvalue

    def GetMeaningfulCircularMatrix(self, labels, aimLabel=1,
                                    minNum=5, pThreshold=0.05, accuracy=0.70):
        """
        """
        matrix = self.GetCircularMatrix()

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
            meanPvalue = meanPvalue[meanPvalue.Accuracy>=accuracy]
            meanMatrix = matrix.reindex(meanPvalue.index, axis=1)
        else:
            pass
        return meanPvalue, meanMatrix

        # pValue = mp.Manager().dict()
        # pool = mp.Pool(self.nJobs)
        # for col, val in matrix.iteritems():
        #     pool.apply_async(self.Pvalue,
        #                      args=(n, m,
        #                            col, val,
        #                            aimLabel, pValue
        #                            )
        #                      )
        # pool.close()
        # pool.join()
        # return res
if '__main__' == __name__:
    from rdkit import Chem
    import openbabel as ob
    from itertools import compress

    def obsmitosmile(smiles):
        conv = ob.OBConversion()
        conv.SetInAndOutFormats("smi", "can")
        conv.SetOptions("K", conv.OUTOPTIONS)
        mol = ob.OBMol()
        conv.ReadString(mol, smiles)
        smile = conv.WriteString(mol)
        smile = smile.replace('\t\n', '')
        return smile

    def getMol(smiles):
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            mol = Chem.MolFromSmiles(obsmitosmile(smiles)) 
        return mol

    # smis = [
    #     'C1=CC=CC(C(Br)C)=C1',
    #     'C1=CC2NC(=O)CC3C=2C(C(=O)C2C=CC=CC=23)=C1',
    #     'C1=CC=C2C(=O)C3C=CNC=3C(=O)C2=C1',
    #     'C1=NC(CCN)=CN1',
    #     'C1CCCC(CCO)C1',
    #     'C1=CC=C2N=C(O)C=CC2=C1',
    #     'C(OC)1=C(C)C=C2OC[C@]([H])3OC4C(C)=C(OC)C=CC=4C(=O)[C@@]3([H])C2=C1C',
    #     'C1=C2N=CC=NC2=C2N=CNC2=C1',
    #     'C1=C(O)C=CC(O)=C1',
    #     'C1=CC2NC(=O)CC3C=2C(C(=O)C2C=CC=CC=23)=C1',
    #     'C1=CC=C2C(=O)C3C=CNC=3C(=O)C2=C1',
    #     'C(OC)1=C(C)C=C2OC[C@]([H])3OC4C(C)=C(OC)C=CC=4C(=O)[C@@]3([H])C2=C1C']



    data = pd.read_csv(r'tests\Canc\Canc.txt', sep='\t')
    # data = data.sample(n=100)

    mols = data.SMILES.map(getMol)
    bo = mols.notna().values
    
    mols = list(compress(mols, bo))
    y_true = data.Label.values[bo]

    circular = MeaningfulCircular(mols, folded=False, maxRadius=6,
                                  minRadius=3, maxFragment=True, nJobs=4)
    meanPvalue, meanMatrix = circular.GetMeaningfulCircularMatrix(y_true)
    print(meanMatrix)

