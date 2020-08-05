# -*- coding: utf-8 -*-
"""
Created on Sun Jun 14 15:13:29 2020

@Author: Zhi-Jiang Yang, Dong-Sheng Cao
@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China
@Homepage: http://www.scbdd.com
@Mail: yzjkid9@gmail.com; oriental-cds@163.com
@Blog: https://blog.iamkotori.com

♥I love Princess Zelda forever♥
"""


from rdkit import Chem
from rdkit.Chem.rdmolops import RDKFingerprint
from rdkit.Chem.rdmolops import UnfoldedRDKFingerprintCountBased
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Draw import _getRDKitEnv

__all__ = ['GetFoldedPathFragment',
           'GetUnfoldedPathFragment',
           'GetPathFragment']


def _DisposePathFragments(mol, bitInfo, maxFragment=True, svg=False):

    def func(frags):
        sets = [set(e) for e in frags]
        us = []
        for e in sets:
            if any(e < s for s in sets):
                continue
            else:
                us.append(str(e))
        return us

    fragments = {}
    keys = list(bitInfo.keys())
    if maxFragment:
        dic = {str(set(x)): k for k, v in bitInfo.items() for x in v}
        bitInfo = [i for j in bitInfo.values() for i in j]
        bitInfo = func(bitInfo)
        for k in bitInfo:
            bondPath = list(eval(k))
            smi, svgImg = DrawRDKitEnv(mol, bondPath)
            fragments[dic[k]] = (smi, svgImg.replace('\n','')) if svg else smi
    else:
        pass

    return (keys, fragments)


def DrawRDKitEnv(mol, bondPath, molSize=(150, 150), baseRad=0.3, useSVG=True,
                 aromaticColor=(0.9, 0.9, 0.2), extraColor=(0.9, 0.9, 0.9), nonAromaticColor=None,
                 drawOptions=None, **kwargs):
    menv = _getRDKitEnv(mol, bondPath, baseRad, aromaticColor,
                        extraColor, nonAromaticColor, **kwargs)
    submol = menv.submol
    subMol = Chem.MolToSmiles(submol, canonical=True, isomericSmiles=False)

    # Drawing
    if useSVG:
        drawer = rdMolDraw2D.MolDraw2DSVG(molSize[0], molSize[1])
    else:
        drawer = rdMolDraw2D.MolDraw2DCairo(molSize[0], molSize[1])

    if drawOptions is None:
        drawopt = drawer.drawOptions()
        drawopt.continuousHighlight = False
    else:
        drawOptions.continuousHighlight = False
        drawer.SetDrawOptions(drawOptions)

    drawer.DrawMolecule(menv.submol, highlightAtoms=menv.highlightAtoms,
                        highlightAtomColors=menv.atomColors, highlightBonds=menv.highlightBonds,
                        highlightBondColors=menv.bondColors, highlightAtomRadii=menv.highlightRadii,
                        **kwargs)
    drawer.FinishDrawing()
    return subMol, drawer.GetDrawingText()


def GetFoldedPathFragment(mol,
                          minPath=1, maxPath=7,
                          nBits=2048, maxFragment=True,
                          svg=False):
    """Calculate folded path fragment.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        the aim molecule.
    minPath : int, optional
        minimum number of bonds to include in the subgraphs. The default is 1.
    maxPath : int, optional
        maximum number of bonds to include in the subgraphs. The default is 1.
    nBits : int, optional
        the number of bit of morgan. The default is 2048.

    Returns
    -------
    dict
        the key of dict is the number of bit (count from 0), 
        and value is a list of list, in each list element the bond consist of 
        substrcuture

    """
    bitInfo = {}
    fp = RDKFingerprint(mol,
                        minPath=minPath,
                        maxPath=maxPath,
                        fpSize=nBits,
                        bitInfo=bitInfo)

    fragments = _DisposePathFragments(
        mol, bitInfo, maxFragment=maxFragment, svg=svg)
    return fragments


def GetUnfoldedPathFragment(mol,
                            minPath=1, maxPath=7,
                            maxFragment=True, svg=False):
    """
    Calculate fingerprint and return the info of each bit

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        the aim molecule.
    minPath : int, optional
        minimum number of bonds to include in the subgraphs. The default is 1.
    maxPath : int, optional
        maximum number of bonds to include in the subgraphs. The default is 1.

    Returns
    -------
    dict
        the key of dict is the number of bit (count from 0), 
        and value is a list of list, in each list element the bond consist of 
        substrcuture

    """
    bitInfo = {}
    fp = UnfoldedRDKFingerprintCountBased(mol,
                                          minPath=minPath,
                                          maxPath=maxPath,
                                          bitInfo=bitInfo)

    fragments = _DisposePathFragments(mol, bitInfo, maxFragment=maxFragment, svg=svg)
    return fragments


def GetPathFragment(mol,
                    minPath=1, maxPath=7,
                    nBits=2048, folded=False,
                    maxFragment=True, svg=False):
    """


    Parameters
    ----------
    mol : TYPE
        DESCRIPTION.
    minPath : TYPE, optional
        DESCRIPTION. The default is 1.
    maxPath : TYPE, optional
        DESCRIPTION. The default is 7.
    nBits : TYPE, optional
        DESCRIPTION. The default is 2048.

    Returns
    -------
    dict
        the substruct of mol and frequency.

    """
    if folded:
        bitInfo = GetFoldedPathFragment(mol,
                                        minPath=minPath, maxPath=maxPath,
                                        nBits=nBits, maxFragment=maxFragment,
                                        svg=svg)
    else:
        bitInfo = GetUnfoldedPathFragment(mol,
                                          minPath=minPath, maxPath=maxPath,
                                          maxFragment=maxFragment, svg=svg)

    # substrcutures = []
    # head = getBeginEndAtom(mol)

    # for info in bitInfo.values():
    #     atomsToUse = [head[bond] for bond in info[0]]
    #     atomsToUse = set(sum(atomsToUse, []))
    #     smi = Chem.MolFragmentToSmiles(mol, atomsToUse, bondsToUse=info[0])
    #     substrcutures.append(smi)

    # substrcutures = list(set(substrcutures))
    # return substrcutures
    return bitInfo


if '__main__' == __name__:

    mol = Chem.MolFromSmiles('CN(C1=CC=C(C=C1)N=N/C2=CC=CC=C2)C')
    frag = GetPathFragment(mol, maxFragment=True, svg=True)
    print(frag)
