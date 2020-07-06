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

import pandas as pd
import numpy as np
import scipy as sc

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import SVG
from rdkit.Chem import rdDepictor
from .fingerprints import Circular, Path, FunctionGroup


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
        denominatorA = Decimal(sc.math.factorial(float(i))) * Decimal(sc.math.factorial(float(ns-i)))
        denominatorB = (m/n)**float(i)
        denominatorC = (1 - (m/n))**(ns-i)
        p_value = float(numerator/denominatorA) * denominatorB * denominatorC
        yield p_value
  

def HighlightAtoms(mol,highlightAtoms,figsize=[400,200],kekulize=True):
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
        svg_words =  svg_words.replace(
            'stroke-width:2px','stroke-width:1.5px').replace(
                'font-size:17px','font-size:15px').replace(
                    'stroke-linecap:butt','stroke-linecap:square').replace(
                        'fill:#FFFFFF','fill:none').replace(
                        'svg:','')                
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
    drawer.DrawMolecule(mc,highlightAtoms=highlightAtoms)
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
        smi = subMatrix[(subMatrix[sma]==1) & (labels==aim_label)
                        ].sample(n=1).index.values[0]
        mol = Chem.MolFromSmiles(smi)
        atoms = mol.GetSubstructMatches(patt)[0]
        svg = HighlightAtoms(mol, atoms)
        imgs.append(svg.data)
    
    ns = subMatrix.loc[:, smas].sum(axis=0).values
    ms = subMatrix.loc[:, smas][labels==aim_label].sum(axis=0).values
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
    